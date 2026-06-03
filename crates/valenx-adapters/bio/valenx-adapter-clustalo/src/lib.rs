//! # valenx-adapter-clustalo
//!
//! Adapter for [Clustal Omega](http://www.clustal.org/omega/) — the
//! third major workhorse of multiple-sequence alignment alongside MAFFT
//! and MUSCLE. Clustal Omega replaces the older ClustalW with a
//! HHalign-based progressive aligner that scales linearly to hundreds
//! of thousands of sequences while retaining accuracy comparable to
//! the iterative MAFFT strategies on small-to-medium inputs.
//!
//! **Phase 18.7 — subprocess wrapper around `clustalo`.** The user
//! supplies a multi-FASTA via `[bio.clustalo]` in `case.toml`.
//! `prepare()` resolves it against the case directory and composes
//! `clustalo -i <input> -o <basename>.<ext> --outfmt=<fmt>
//! --threads=<N>`, picking the extension to match the requested output
//! format so `collect()` can find the result deterministically.
//!
//! Unlike MAFFT (which writes to stdout and needs a custom run loop to
//! redirect FDs), Clustal Omega has a real `-o` flag — its output goes
//! directly to a file, leaving stdout/stderr free for the standard
//! shared subprocess runner. That makes this the simplest of the three
//! MSA adapters.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;
pub mod native;

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

use crate::case_input::ClustaloInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(ClustaloAdapter::new())
}

pub struct ClustaloAdapter;

impl ClustaloAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for ClustaloAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "clustalo";
/// Clustal Omega's binary candidates. The canonical name is
/// `clustalo` everywhere — Bioconda, Homebrew, Debian, and Fedora
/// all ship under that exact name.
const BINARIES: &[&str] = &["clustalo"];

impl Adapter for ClustaloAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "Clustal Omega",
            // Clustal Omega's 1.2 line has been the stable series since
            // 2014; current releases land in the 1.2.x band (1.2.4 is
            // the de facto distro version). Floor at 1.2.0; upper bound
            // 2.0 reserves room for the next major.
            version_range: VersionRange {
                min_inclusive: Version::new(1, 2, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-2.0",
            docs_url: "http://www.clustal.org/omega/",
            homepage_url: "http://www.clustal.org/omega/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `clustalo --version` prints "1.2.4" on stdout — bare
                // semver, no prefix, parses cleanly.
                let found_version = detect_tool_version_semver(&binary_path, &["--version", ""]);
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    warnings: Vec::new(),
                    required_env: Vec::new(),
                })
            }
            // Native Rust fallback via valenx-align progressive+iterative MSA.
            None => Ok(ProbeReport {
                ok: true,
                found_version: None,
                binary_path: None,
                warnings: vec![
                    "clustalo binary not found; using native Rust progressive+iterative MSA \
                     (valenx-align). Install Clustal Omega 1.2.0+ via apt/brew/conda for the \
                     full ClustalΩ HHalign-based progressive aligner."
                        .to_string(),
                ],
                required_env: Vec::new(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = ClustaloInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.clustalo].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Resolve the input FASTA against the case directory if
        // relative. Same convention as every other Phase 17/18 bio
        // adapter — `input = "seqs.fa"` next to `case.toml`.
        let source_input = if input.input.is_absolute() {
            input.input.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(
            &case.path,
            &input.input,
        )?
        };
        if !source_input.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.clustalo].input `{}` not found (resolved {})",
                    input.input.display(),
                    source_input.display()
                ),
            });
        }

        // Pick the output filename: `<basename>.<ext>` where ext matches
        // the requested outfmt. The native path always writes FASTA, so
        // use `.fasta` as the extension when falling back.
        let ext = extension_for_outfmt(&input.outfmt);
        let output_filename = format!("{}.{}", input.output_basename, ext);

        // Write native_params.toml for the native path.
        let native_params = native::NativeMsaParams {
            input_path: source_input
                .to_str()
                .ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "input path is not valid UTF-8: {}",
                        source_input.display()
                    ))
                })?
                .to_string(),
            output_name: output_filename.clone(),
            refine: true,
            max_iterations: 8,
        };
        native::write_params(workdir, &native_params)?;

        let native_command: Vec<OsString> = match find_on_path(BINARIES) {
            Some(binary_path) => {
                let mut cmd: Vec<OsString> = vec![
                    binary_path.into_os_string(),
                    OsString::from("-i"),
                    source_input.into_os_string(),
                    OsString::from("-o"),
                    OsString::from(&output_filename),
                    OsString::from(format!("--outfmt={}", input.outfmt)),
                    OsString::from(format!("--threads={}", input.threads)),
                ];
                for arg in &input.extra_args {
                    cmd.push(OsString::from(arg));
                }
                cmd
            }
            None => vec![OsString::from(native::NATIVE_SENTINEL)],
        };

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Small inputs (~100 sequences) finish in seconds; large
            // ones (>10k sequences) can run for hours despite Clustal
            // Omega's HHalign scaling. 4 hours covers the long tail
            // without being absurd, matching MAFFT/BLAST.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        // Native Rust path: progressive + iterative MSA, no subprocess.
        if job.native_command.first().map(|s| s.as_os_str())
            == Some(native::NATIVE_SENTINEL.as_ref())
        {
            return native::run_native(&job.workdir, ctx);
        }

        let report = subprocess::run(job, ctx, "starting Clustal Omega", |line| {
            let mut hint = subprocess::Hint::default();
            // Clustal Omega is fairly quiet by default; explicit
            // `ERROR` markers are the only thing worth lifting as
            // warnings — anything else is just startup chatter.
            if line.contains("ERROR") || line.contains("FATAL") {
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
        // Surface every file in the workdir that starts with the
        // configured output basename — covers the alignment file
        // itself plus anything Clustal Omega writes adjacent (a
        // `--guidetree-out` file, for instance, if the user passed it
        // via extra_args). Provenance hashes the case.toml since we
        // don't know the exact alignment path here without re-parsing
        // (the basename lives in [bio.clustalo].output_basename).
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "Clustal Omega",
            "unknown",
            &job.workdir.join("case.toml"),
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Recover the basename from the prepared command. In subprocess
        // mode, `prepare()` pushed `-o <basename>.<ext>` so the filename
        // sits at index 4. In native mode the command is just the
        // sentinel, so we read the output_name from native_params.toml
        // instead.
        let basename = {
            let from_cmd = job
                .native_command
                .get(4)
                .map(Path::new)
                .and_then(|p| p.file_stem())
                .and_then(|s| s.to_str())
                .map(str::to_string);
            match from_cmd {
                Some(b) if !b.is_empty() && b != "valenx:native:msa" => b,
                _ => native::read_params(&job.workdir)
                    .ok()
                    .and_then(|p| {
                        Path::new(&p.output_name)
                            .file_stem()
                            .and_then(|s| s.to_str())
                            .map(str::to_string)
                    })
                    .unwrap_or_default(),
            }
        };

        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-clustalo", ?e, "workdir read failed");
                return Ok(results);
            }
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if !path.is_file() {
                continue;
            }
            let file_name = path
                .file_name()
                .and_then(|s| s.to_str())
                .map(|s| s.to_string());
            let ext = path
                .extension()
                .and_then(|s| s.to_str())
                .map(|s| s.to_ascii_lowercase());

            let starts_with_basename = matches!(
                (&file_name, basename.as_str()),
                (Some(name), bn) if !bn.is_empty() && name.starts_with(bn)
            );

            let (kind, label) = if starts_with_basename {
                (ArtifactKind::Tabular, "Clustal Omega alignment".to_string())
            } else if ext.as_deref() == Some("log") {
                (ArtifactKind::Log, "Clustal Omega log".to_string())
            } else {
                continue;
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
        // ribbon contributions are already enough for the registry to
        // surface the adapter.
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.clustalo.align"],
        }
    }
}

/// Pick a file extension that matches Clustal Omega's `--outfmt` value.
/// The mapping covers the four formats users actually pick day-to-day
/// plus a safe default for anything Clustal Omega itself accepts but
/// that we don't have a canonical extension for. Picking sensibly here
/// keeps the workdir self-describing — a stray `alignment.aln` is
/// instantly recognisable as the Clustal-format output, while a
/// `alignment.fasta` makes downstream tools (which dispatch on
/// extension) Just Work.
fn extension_for_outfmt(outfmt: &str) -> &'static str {
    match outfmt {
        "clustal" => "aln",
        "fasta" => "fasta",
        "phylip" => "phy",
        "vienna" => "vie",
        "nexus" => "nex",
        _ => "aln",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = ClustaloAdapter::new().info();
        assert_eq!(info.id, "clustalo");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "GPL-2.0");
        assert_eq!(info.display_name, "Clustal Omega");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = ClustaloAdapter::new().info();
        // 1.2.0 is the floor we test against; 2.0 reserves room for
        // an eventual major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(1, 2, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(2, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = ClustaloAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.clustalo.align"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = ClustaloAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    #[test]
    fn probe_always_succeeds() {
        let report = ClustaloAdapter::new().probe().unwrap();
        assert!(report.ok, "probe() must return ok=true in all modes");
    }

    #[test]
    fn native_sentinel_is_stable() {
        assert_eq!(native::NATIVE_SENTINEL, "valenx:native:msa");
    }
}
