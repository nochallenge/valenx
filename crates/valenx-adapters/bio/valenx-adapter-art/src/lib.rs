//! # valenx-adapter-art
//!
//! Adapter for [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/)
//! — Weichun Huang's NIEHS Illumina-platform read simulator. ART
//! synthesises FASTQs that match the empirical error profile of a
//! given Illumina sequencing system (HiSeq 2500, HiSeq X, MiSeq v3,
//! NextSeq 500, MiniSeq) so downstream pipelines can be validated
//! against a known-truth reference at controlled coverage and read
//! length. The `art_illumina` binary is the workhorse — companion
//! tools cover 454 (`art_454`) and SOLiD (`art_SOLiD`), which this
//! adapter does not surface.
//!
//! **Phase 31 — subprocess wrapper around `art_illumina`.** The user
//! supplies a reference FASTA plus the per-system error profile and
//! coverage parameters via `[bio.art]` in `case.toml`. `prepare()`
//! resolves the reference against the case directory, validates the
//! sequencing system / coverage / fragment-size combination, and
//! composes the `art_illumina` invocation. `run()` streams the
//! generation via the shared subprocess runner; ART writes a
//! `<prefix>.fq` (single-end) or `<prefix>1.fq` + `<prefix>2.fq`
//! pair (paired-end) into the workdir.
//!
//! On `collect()` we walk the workdir for any `<prefix>*.fq`
//! (simulated-reads FASTQ) plus `<prefix>*.aln` (the per-read
//! alignment record ART writes alongside, useful for validating
//! aligner accuracy against the simulated truth).

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::fs;
use std::path::Path;
use std::time::Duration;

use semver::Version;

use valenx_core::{
    adapter_helpers::{detect_tool_version_semver, find_on_path, live_provenance},
    error::RunPhase,
    subprocess, Adapter, AdapterError, AdapterInfo, Capabilities, Case, LicenseMode, Physics,
    PreparedJob, ProbeReport, RunContext, RunReport, VersionRange,
};
use valenx_fields::{
    artifact::{Artifact, ArtifactKind},
    Results,
};

use crate::case_input::ArtInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(ArtAdapter::new())
}

pub struct ArtAdapter;

impl ArtAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for ArtAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "art";
/// ART's binary candidates. `art_illumina` is the canonical name from
/// NIEHS source builds, Bioconda, and Homebrew; the `art_454` /
/// `art_SOLiD` companions live alongside it but cover platforms this
/// adapter does not surface.
const BINARIES: &[&str] = &["art_illumina"];

impl Adapter for ArtAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "ART",
            // ART's published version line is ChocolateCherryCake
            // (the long-running "2.5.x" series since 2016). Bioconda
            // ships 2.5.8; Homebrew ships 2.5.8. Floor at 2.5.0
            // covers every reasonable install; 3.0 reserves room
            // for an eventual major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(2, 5, 0),
                max_exclusive: Version::new(3, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-3.0",
            docs_url: "https://www.niehs.nih.gov/research/resources/software/biostatistics/art/",
            homepage_url:
                "https://www.niehs.nih.gov/research/resources/software/biostatistics/art/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `art_illumina` (no args) prints a usage banner that
                // includes the version on stderr; some 2.5.x point
                // releases learned `--version`. Try both — the
                // detector reads stdout+stderr regardless.
                let found_version = detect_tool_version_semver(&binary_path, &["--version", ""]);
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    warnings: Vec::new(),
                    required_env: Vec::new(),
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "ART 2.5+ required; install via `apt install art-nextgen-simulation-tools`, \
                       `brew install art`, or `conda install -c bioconda art`"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = ArtInput::from_case_dir(&case.path)?;

        // Round-3 security fix: `output_prefix` is fed to ART's `-o`
        // and becomes the basename of every output FASTQ. A hostile
        // value like `"../../etc/cron.d/x"` would otherwise let ART
        // write outside the workdir.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_prefix,
            "[bio.art].output_prefix",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Resolve the reference path against the case directory if
        // relative. Same convention as every other Phase 17/18 bio
        // adapter — `reference = "ref.fa"` next to `case.toml`.
        let source_reference = if input.reference.is_absolute() {
            input.reference.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(
            &case.path,
            &input.reference,
        )?
        };
        if !source_reference.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.art].reference `{}` not found (resolved {})",
                    input.reference.display(),
                    source_reference.display()
                ),
            });
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "ART 2.5+ required; install via `apt install art-nextgen-simulation-tools`, \
                       `brew install art`, or `conda install -c bioconda art`"
                .into(),
        })?;

        // Compose `art_illumina -ss <ss> -i <ref> -l <len> -f <cov>
        // -o <prefix> [-p -m <mean> -s <sd>] [extras...]`. ART's
        // `-o` takes a prefix (no extension) — single-end runs land
        // as `<prefix>.fq`, paired-end runs as `<prefix>1.fq` +
        // `<prefix>2.fq`. The cwd for the subprocess is the workdir,
        // so a relative prefix lands inside.
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from("-ss"),
            OsString::from(&input.sequencing_system),
            OsString::from("-i"),
            source_reference.into_os_string(),
            OsString::from("-l"),
            OsString::from(input.read_length.to_string()),
            OsString::from("-f"),
            OsString::from(format!("{}", input.fold_coverage)),
            OsString::from("-o"),
            OsString::from(&input.output_prefix),
        ];
        if input.paired_end {
            native_command.push(OsString::from("-p"));
            native_command.push(OsString::from("-m"));
            native_command.push(OsString::from(format!("{}", input.fragment_mean)));
            native_command.push(OsString::from("-s"));
            native_command.push(OsString::from(format!("{}", input.fragment_sd)));
        }
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Whole-genome 30x simulation finishes in minutes on a
            // single core; very deep coverage on large references
            // can run for an hour or two. 4 hours is generous.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting ART", |line| {
            let mut hint = subprocess::Hint::default();
            // ART's stdout chatter (mirrored here through the line
            // handler) emits a banner at startup, "Profile" /
            // "Reference" lines mid-setup, and a "Done!" sentinel
            // at the end. Lift the obvious markers to UI ticks.
            if line.contains("Done!") || line.contains("Total CPU time") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("Reference Sequences") || line.contains("Profile") {
                hint.progress = Some((30.0, line.to_string()));
            } else if line.contains("read length") || line.contains("fold coverage") {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.to_ascii_lowercase().contains("error") {
                hint.warning = Some(line.trim().to_string());
            }
            hint
        })?;
        Ok(RunReport {
            exit_code: report.exit_code,
            wall_time: report.wall_time,
            converged: Some(true),
            residual_history: Vec::new(),
            warnings: report.warnings,
            final_phase: Some(RunPhase::Shutdown),
        })
    }

    fn collect(&self, job: &PreparedJob) -> Result<Results, AdapterError> {
        // Provenance: hash the case.toml as the canonical input —
        // the simulated FASTQ filenames depend on the chosen prefix
        // and pairing mode, so case.toml is the stable choice.
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "ART",
            "unknown",
            &job.workdir.join("case.toml"),
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk the workdir top-level. ART writes `<prefix>.fq` (or
        // `<prefix>1.fq` + `<prefix>2.fq` for paired-end) plus
        // `<prefix>.aln` / `<prefix>1.aln` + `<prefix>2.aln`. We
        // surface every `.fq` as Tabular (the canonical FASTQ
        // payload) and every `.aln` as Log (the per-read
        // alignment record).
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-art", ?e, "workdir read failed");
                return Ok(results);
            }
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if !path.is_file() {
                continue;
            }
            let ext = path
                .extension()
                .and_then(|s| s.to_str())
                .map(|s| s.to_ascii_lowercase());
            let (kind, label) = match ext.as_deref() {
                Some("fq") | Some("fastq") => {
                    (ArtifactKind::Tabular, "ART simulated reads".to_string())
                }
                Some("aln") => (ArtifactKind::Log, "ART alignment record".to_string()),
                _ => continue,
            };
            artefacts.push(Artifact {
                path,
                kind,
                checksum: None,
                label,
            });
        }
        artefacts.sort_by(|a, b| a.path.cmp(&b.path));
        results.artifacts = artefacts;
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        // Bio-specific Capability variants land in a follow-up task;
        // ribbon contributions are already enough for the registry
        // to surface the adapter.
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.art.simulate"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = ArtAdapter::new().info();
        assert_eq!(info.id, "art");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "GPL-3.0");
        assert_eq!(info.display_name, "ART");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = ArtAdapter::new().info();
        // 2.5 is the floor we test against (every distro ships
        // >= 2.5.0); 3.0 reserves room for an eventual major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(2, 5, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(3, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = ArtAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.art.simulate"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = ArtAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
