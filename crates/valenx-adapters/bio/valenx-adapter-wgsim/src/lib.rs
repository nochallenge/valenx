//! # valenx-adapter-wgsim
//!
//! Adapter for [wgsim](https://github.com/lh3/wgsim) — Heng Li's
//! classic Whole-Genome SIMulator that ships alongside samtools.
//! wgsim emits paired-end FASTQs from a reference under a uniform
//! sequencing-error model with configurable insert size, read length,
//! and per-base error rate. Unlike ART (which models per-platform
//! empirical error profiles), wgsim is deliberately simple — the
//! canonical "small + classic" simulator for fast smoke-testing of
//! mappers and variant callers when realistic error spectra are not
//! required.
//!
//! **Phase 31 — subprocess wrapper around `wgsim`.** The user
//! supplies a reference FASTA plus the desired output FASTQ names
//! and pair count via `[bio.wgsim]` in `case.toml`. `prepare()`
//! resolves the reference, validates the read-length / fragment /
//! error parameters, and composes the `wgsim` invocation. `run()`
//! streams the generation via the shared subprocess runner; wgsim
//! writes both FASTQs into the workdir at the names the user
//! specified.
//!
//! On `collect()` we surface the two FASTQs as Tabular artifacts
//! so downstream pipelines (or the UI's results pane) can find
//! them by their canonical labels.

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

use crate::case_input::WgsimInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(WgsimAdapter::new())
}

pub struct WgsimAdapter;

impl WgsimAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for WgsimAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "wgsim";
/// wgsim's binary candidates. The canonical name is `wgsim` from
/// every distro that packages samtools.
const BINARIES: &[&str] = &["wgsim"];

impl Adapter for WgsimAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "wgsim",
            // wgsim is versioned alongside samtools; the source
            // tree pins to the parent samtools 1.x line and a
            // separate 1.0+ tag has been the long-running stable
            // series. 2.0 reserves room for an eventual major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(1, 0, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "MIT",
            docs_url: "https://github.com/lh3/wgsim",
            homepage_url: "https://github.com/lh3/wgsim",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `wgsim` (no args) prints a short usage banner that
                // includes a "Version:" line on stderr; the detector
                // reads stdout+stderr regardless.
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
                hint: "wgsim 1.0+ required; install via `apt install wgsim`, \
                       `brew install wgsim`, or `conda install -c bioconda wgsim`"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = WgsimInput::from_case_dir(&case.path)?;

        // Round-10 H3: `output1` and `output2` are `PathBuf` from
        // case.toml and flowed into `workdir.join(...)` in collect()
        // with no validation. wgsim writes a single FASTQ file per
        // mate (not a directory), so basename-only is correct.
        for (field, value) in [
            ("[bio.wgsim].output1", &input.output1),
            ("[bio.wgsim].output2", &input.output2),
        ] {
            if let Some(s) = value.to_str() {
                valenx_core::adapter_helpers::validate_output_basename(s, field)
                    .map_err(|e| AdapterError::InvalidCase {
                        case_path: case.path.join("case.toml"),
                        reason: format!("{e}"),
                    })?;
            } else {
                return Err(AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!("{field}: non-UTF-8 path rejected"),
                });
            }
        }

        fs::create_dir_all(workdir)?;

        // Stage `case.toml` into the workdir so collect() can recover
        // the configured `output_basename` for prefix-filtering output
        // artifacts. Without this stage, the basename filter silently
        // degrades to "match everything".
        let staged_case_toml = workdir.join("case.toml");
        let source_case_toml = case.path.join("case.toml");
        if source_case_toml.is_file() {
            fs::copy(&source_case_toml, &staged_case_toml)
                .map_err(|e| AdapterError::Other(anyhow::anyhow!("stage case.toml: {e}")))?;
        }

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
                    "[bio.wgsim].reference `{}` not found (resolved {})",
                    input.reference.display(),
                    source_reference.display()
                ),
            });
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "wgsim 1.0+ required; install via `apt install wgsim`, \
                       `brew install wgsim`, or `conda install -c bioconda wgsim`"
                .into(),
        })?;

        // Compose `wgsim -N <num_pairs> -1 <length1> -2 <length2>
        // -d <fragment_size> -e <error_rate> <reference> <out1> <out2>
        // [extras...]`. The two output paths are interpreted relative
        // to the subprocess cwd (the workdir), so they land inside.
        //
        // Round-3 fix: extras MUST come after the positional triple
        // (`<reference> <out1> <out2>`). Pre-fix they were appended
        // between `-e <error_rate>` and the positionals, which let a
        // hostile case.toml slip an extra positional in (e.g.
        // `extra_args = ["phantom_positional"]`) and shift `<out2>`
        // onto whatever wgsim parsed `phantom_positional` as — making
        // it overwrite a different file than the user expected.
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from("-N"),
            OsString::from(input.num_pairs.to_string()),
            OsString::from("-1"),
            OsString::from(input.length1.to_string()),
            OsString::from("-2"),
            OsString::from(input.length2.to_string()),
            OsString::from("-d"),
            OsString::from(input.fragment_size.to_string()),
            OsString::from("-e"),
            OsString::from(format!("{}", input.error_rate)),
        ];
        native_command.push(source_reference.into_os_string());
        native_command.push(OsString::from(input.output1.as_os_str()));
        native_command.push(OsString::from(input.output2.as_os_str()));
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // wgsim is fast — 1M pair simulation finishes in
            // seconds; whole-genome at high coverage runs for
            // a few minutes. 1 hour is generous.
            estimated_runtime: Some(Duration::from_secs(60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting wgsim", |line| {
            let mut hint = subprocess::Hint::default();
            // wgsim prints a few setup lines and a per-chromosome
            // counter on stderr; the runner forwards both stdout
            // and stderr through this handler at log level. We
            // can lift the obvious markers to coarse UI ticks.
            if line.contains("seed:") || line.contains("Random seed") {
                hint.progress = Some((10.0, line.to_string()));
            } else if line.contains("simulation") || line.contains("haplotype") {
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
        // the FASTQ filenames are user-chosen and can vary, so
        // case.toml is the stable choice.
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "wgsim",
            "unknown",
            &job.workdir.join("case.toml"),
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);

        // Re-resolve the user-chosen output paths against the
        // workdir. wgsim writes them there because the subprocess
        // cwd is the workdir; we don't walk the directory because
        // the pair (output1, output2) is exactly what the user
        // asked for and labelling them "read 1" / "read 2"
        // matches the case.toml schema.
        let input = match WgsimInput::from_case_dir(&job.workdir) {
            Ok(i) => i,
            Err(e) => {
                tracing::warn!(target: "valenx-wgsim", ?e, "collect could not re-read case.toml");
                return Ok(results);
            }
        };
        let mut artefacts: Vec<Artifact> = Vec::new();
        let r1 = if input.output1.is_absolute() {
            input.output1.clone()
        } else {
            job.workdir.join(&input.output1)
        };
        let r2 = if input.output2.is_absolute() {
            input.output2.clone()
        } else {
            job.workdir.join(&input.output2)
        };
        artefacts.push(Artifact {
            path: r1,
            kind: ArtifactKind::Tabular,
            checksum: None,
            label: "wgsim read 1".to_string(),
        });
        artefacts.push(Artifact {
            path: r2,
            kind: ArtifactKind::Tabular,
            checksum: None,
            label: "wgsim read 2".to_string(),
        });
        results.artifacts = artefacts;
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        // Bio-specific Capability variants land in a follow-up task;
        // ribbon contributions are already enough for the registry
        // to surface the adapter.
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.wgsim.simulate"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = WgsimAdapter::new().info();
        assert_eq!(info.id, "wgsim");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "MIT");
        assert_eq!(info.display_name, "wgsim");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = WgsimAdapter::new().info();
        // 1.0 is the floor; 2.0 reserves room for an eventual
        // major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(1, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(2, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = WgsimAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.wgsim.simulate"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = WgsimAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    /// Round-10 H3 RED→GREEN: `output1` flowed into
    /// `workdir.join(...)` with no validation. Hostile
    /// `output1 = "../etc/passwd"` is now rejected.
    #[test]
    fn prepare_rejects_output1_path_traversal() {
        use valenx_test_utils::tempdir;
        let d = tempdir("wgsim-out1-trav");
        std::fs::write(d.join("ref.fa"), b">x\nACGT\n").unwrap();
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "wgsim.simulate"

[bio.wgsim]
reference     = "ref.fa"
output1       = "../etc/passwd"
output2       = "r2.fq"
num_pairs     = 1000
length1       = 100
length2       = 100
fragment_size = 500
error_rate    = 0.01
"#,
        )
        .unwrap();
        let case = Case {
            id: "trav".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = WgsimAdapter::new()
            .prepare(&case, &workdir)
            .unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("[bio.wgsim].output1"),
            "expected [bio.wgsim].output1 in error, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
