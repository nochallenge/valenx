//! # valenx-adapter-mfold
//!
//! Adapter for [mfold / UNAFold](http://www.unafold.org/) — Michael
//! Zuker's classic RNA secondary-structure folder. mfold computes
//! minimum-free-energy folds and a small set of suboptimal structures
//! over the canonical Turner thermodynamic model, and remains the
//! reference tool that newer folders (ViennaRNA, RNAstructure, ML
//! folders) routinely benchmark against. UNAFold is the modern
//! successor distribution that ships the same engine plus a few CLI
//! conveniences.
//!
//! **Phase 44.5 — subprocess wrapper around `mfold` / `UNAFold.pl`.**
//! The user supplies a single-sequence input (FASTA or the legacy
//! `.seq` format) via `[bio.mfold]` in `case.toml`. `prepare()`
//! resolves the sequence against the case directory and composes
//! `mfold SEQ=<sequence> NA=RNA T=<temperature> [extras...]` — mfold
//! uses `KEY=VALUE` argument syntax rather than `--flag value` pairs.
//! `run()` streams progress via the shared subprocess runner, and
//! `collect()` walks the workdir for `*.ct` (connect-table), `*.ps` /
//! `*.pdf` (structure plots), and `*.out` (run log).
//!
//! ## License flag
//!
//! mfold is free for academic / non-commercial use; commercial use
//! requires a separate license from the upstream maintainer. We
//! surface this via `tool_license = "Academic"` and emit a probe
//! warning when the binary is found. The probe-warning text contains
//! the literal strings `"academic"` and `"non-commercial"` as stable
//! anchors for tests and downstream license-aware filters.

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

use crate::case_input::MfoldInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(MfoldAdapter::new())
}

pub struct MfoldAdapter;

impl MfoldAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for MfoldAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "mfold";
/// mfold's binary candidates. The classic mfold ships a `mfold` driver;
/// the UNAFold successor adds `UNAFold.pl` (a Perl wrapper) alongside
/// it. Both are accepted — they accept the same `KEY=VALUE` argument
/// surface for single-sequence folding.
const BINARIES: &[&str] = &["mfold", "UNAFold.pl"];

/// The probe-warning surfaced whenever mfold is detected. The literal
/// strings `"academic"` and `"non-commercial"` are part of the asserted
/// contract — they anchor the license reminder so downstream
/// license-aware filters and tests can key off stable substrings.
/// Sister to NAMD / NUPACK / VMD's identical pattern.
const LICENSE_WARNING: &str = "mfold is licensed for non-commercial / academic use only. \
     Commercial use requires a separate license from the upstream \
     maintainer; confirm your use case complies before redistributing \
     folds or derived data.";

impl Adapter for MfoldAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "mfold/UNAFold",
            // mfold's 3.x line is the long-running stable series; the
            // UNAFold-branded successor ships with the same engine
            // versioned in step. 3.8 is the modern floor; upper bound
            // 4.0 reserves room for the next major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(3, 8, 0),
                max_exclusive: Version::new(4, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            // mfold's terms aren't a recognised SPDX identifier; the
            // closest accurate label is the project's own
            // academic-only license. Mislabeling as MIT / BSD would
            // be misleading.
            tool_license: "Academic",
            docs_url: "http://www.unafold.org/Dinamelt/software/mfold-software.php",
            homepage_url: "http://www.unafold.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `mfold --version` on recent builds prints the version
                // banner; older builds print it as part of the help
                // text. The combined detector handles both.
                let found_version = detect_tool_version_semver(&binary_path, &["--version", "-h"]);
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    // Always surface the license reminder when mfold
                    // is detected — it's a custom non-OSS license and
                    // we'd rather over-warn than have a user ship
                    // commercial output without checking.
                    warnings: vec![LICENSE_WARNING.to_string()],
                    required_env: Vec::new(),
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "mfold 3.8+ required; install via `conda install \
                       -c bioconda unafold`, the upstream installer at \
                       http://www.unafold.org/, or build from source"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = MfoldInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.mfold].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Resolve the sequence file against the case directory if
        // relative. mfold reads it in place via `SEQ=<path>`; we don't
        // stage it into the workdir, just validate it exists so the
        // failure is fast and obvious.
        let source_sequence = if input.sequence.is_absolute() {
            input.sequence.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(
            &case.path,
            &input.sequence,
        )?
        };
        if !source_sequence.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.mfold].sequence `{}` not found (resolved {})",
                    input.sequence.display(),
                    source_sequence.display()
                ),
            });
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "mfold 3.8+ required; install via `conda install \
                       -c bioconda unafold`, the upstream installer at \
                       http://www.unafold.org/, or build from source"
                .into(),
        })?;

        // Compose `mfold SEQ=<sequence> NA=RNA T=<temperature>
        //          [extras...]`. mfold uses `KEY=VALUE` argument syntax
        // — every knob is a positional `KEY=VALUE` token rather than a
        // `--flag value` pair. The temperature is rendered with
        // format!("{}") so whole-number degrees serialise without a
        // trailing `.0`, matching what users typed in case.toml.
        let mut native_command: Vec<OsString> = vec![binary_path.into_os_string()];
        let mut seq_arg = OsString::from("SEQ=");
        seq_arg.push(source_sequence.as_os_str());
        native_command.push(seq_arg);
        native_command.push(OsString::from("NA=RNA"));
        native_command.push(OsString::from(format!("T={}", input.temperature)));
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // A short tRNA finishes in seconds; long mRNAs at low
            // temperature with many suboptimal structures requested can
            // run for many minutes. 30 minutes covers the long tail.
            estimated_runtime: Some(Duration::from_secs(30 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting mfold", |line| {
            let mut hint = subprocess::Hint::default();
            // mfold's stdout chatter: a startup banner, "Energy = ..."
            // mid-run, and a "Done" sentinel near the end. Heuristics
            // only — a mismatch just leaves the spinner alone.
            if line.contains("Done") || line.contains("done") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("Folding") || line.contains("folding") {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.contains("mfold") || line.contains("UNAFold") {
                hint.progress = Some((5.0, line.to_string()));
            } else if line.contains("ERROR") || line.contains("Error:") {
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
        // Provenance: hash the case.toml as the canonical input
        // descriptor. mfold's outputs are named off the input filename
        // stem (mfold's own convention; we don't try to second-guess
        // it), so a single fixed-name artifact for the prov hash isn't
        // workable — case.toml is the stable anchor.
        let case_hash_input = job.workdir.join("case.toml");
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "mfold",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk the workdir top-level. mfold's three canonical output
        // families are `*.ct` (connect-table — the canonical
        // base-pairing description), `*.ps` / `*.pdf` (structure
        // plots), and `*.out` (run log; the legacy mfold convention).
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-mfold", ?e, "workdir read failed");
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
                Some("ct") => (ArtifactKind::Tabular, "mfold connect-table".to_string()),
                Some("ps") | Some("pdf") => {
                    (ArtifactKind::Native, "mfold structure plot".to_string())
                }
                Some("out") => (ArtifactKind::Log, "mfold log".to_string()),
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
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.mfold.fold"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = MfoldAdapter::new().info();
        assert_eq!(info.id, "mfold");
        assert_eq!(info.physics, &[Physics::Bio]);
        // The license identifier surfaces mfold's academic-only
        // license rather than mislabeling as MIT / BSD.
        assert_eq!(info.tool_license, "Academic");
        assert_eq!(info.display_name, "mfold/UNAFold");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = MfoldAdapter::new().info();
        // mfold 3.8 is the modern floor; 4.0 reserves room for the
        // next major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(3, 8, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(4, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = MfoldAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.mfold.fold"]);
    }

    #[test]
    #[ignore] // subprocess-coupled test — run interactively only
    fn probe_warning_mentions_academic_and_non_commercial() {
        // The license-flag warning is mandatory: mfold is non-OSS
        // academic-use, and we surface that on every successful
        // probe. The literal `"academic"` and `"non-commercial"`
        // substrings are what downstream tooling and license-aware
        // filters key off — pin both.
        assert!(
            LICENSE_WARNING.contains("academic"),
            "probe warning must contain `academic` anchor; got: {LICENSE_WARNING}"
        );
        assert!(
            LICENSE_WARNING.contains("non-commercial"),
            "probe warning must contain `non-commercial` anchor; got: {LICENSE_WARNING}"
        );

        // Best-effort live probe — only assert if mfold is on PATH.
        // Skipping when it isn't keeps the test green on CI machines
        // without the (registration-walled) binary.
        if find_on_path(BINARIES).is_some() {
            let report = MfoldAdapter::new().probe().expect("probe");
            assert!(
                report
                    .warnings
                    .iter()
                    .any(|w| w.contains("academic") && w.contains("non-commercial")),
                "live probe warnings must surface the academic / non-commercial anchors; \
                 got: {:?}",
                report.warnings
            );
        }
    }
}
