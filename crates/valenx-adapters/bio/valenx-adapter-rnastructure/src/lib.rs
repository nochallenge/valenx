//! # valenx-adapter-rnastructure
//!
//! Adapter for [RNAstructure](https://rna.urmc.rochester.edu/RNAstructure.html) —
//! David Mathews' lab at the University of Rochester ships a suite of
//! C++ tools for RNA secondary-structure prediction with the
//! Mathews-lab thermodynamic parameter set. `Fold` is the workhorse:
//! reads a `.seq` (or `.fa`) sequence file, writes a connectivity
//! table (`.ct`) holding the MFE plus a configurable number of
//! suboptimal structures.
//!
//! **Phase 28 — subprocess wrapper around `Fold`.** The user supplies
//! a sequence file and an output `.ct` path via `[bio.rnastructure]`
//! in `case.toml`. `prepare()` resolves both against the case
//! directory and composes the invocation:
//!
//! ```text
//! Fold <input> <output> -m <max_structures> -p <max_percent> -t <kelvin> [extras...]
//! ```
//!
//! Unlike ViennaRNA's `RNAfold`, `Fold` writes its primary output
//! (the CT file) directly to the path passed as the second positional
//! argument, so no stdout-redirect dance is needed — the shared
//! `subprocess::run` runner handles it like every other Phase 17/18
//! single-binary bio adapter (BWA, minimap2, etc.).
//!
//! ## License
//!
//! RNAstructure ships under the BSD 3-Clause license — full open-
//! source / commercial use, no academic-only restriction. We surface
//! that accurately via `tool_license = "BSD-3-Clause"` and emit no
//! probe warning.

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

use crate::case_input::RnaStructureInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(RnaStructureAdapter::new())
}

pub struct RnaStructureAdapter;

impl RnaStructureAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for RnaStructureAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "rnastructure";
/// RNAstructure's binary candidates. The Fold tool installs as the
/// canonical capital `Fold` — the upstream binary distribution, the
/// Bioconda package, and the source build all use this casing
/// (mirrors STAR, FastTree).
const BINARIES: &[&str] = &["Fold"];

impl Adapter for RnaStructureAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "RNAstructure",
            // RNAstructure's 6.4.x line is the modern stable series
            // (6.4 in 2021, point releases through 6.4.x). 6.4
            // covers the canonical CLI surface scripts target;
            // upper-bound 7.0 reserves room for an eventual major.
            version_range: VersionRange {
                min_inclusive: Version::new(6, 4, 0),
                max_exclusive: Version::new(7, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "BSD-3-Clause",
            docs_url: "https://rna.urmc.rochester.edu/RNAstructure.html",
            homepage_url: "https://rna.urmc.rochester.edu/RNAstructureWeb/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `Fold --version` prints "Fold X.Y.Z" to stdout; the
                // generic detector picks up the leading SemVer.
                let found_version = detect_tool_version_semver(&binary_path, &["--version"]);
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
                hint: "RNAstructure 6.4+ required; install via \
                       `conda install -c bioconda rnastructure` or \
                       download from https://rna.urmc.rochester.edu/RNAstructure.html"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = RnaStructureInput::from_case_dir(&case.path)?;

        // Round-10 H3: `output` is `PathBuf` and pre-fix flowed into
        // `workdir.join(&input.output)`. Validate as a basename
        // before the join so `output = "../etc/passwd"` is rejected.
        if let Some(s) = input.output.to_str() {
            valenx_core::adapter_helpers::validate_output_basename(s, "[bio.rnastructure].output")
                .map_err(|e| AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!("{e}"),
                })?;
        } else {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: "[bio.rnastructure].output: non-UTF-8 path rejected".into(),
            });
        }

        fs::create_dir_all(workdir)?;

        // Resolve the sequence input against the case directory if
        // relative.
        let source_input = if input.input.is_absolute() {
            input.input.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.input)?
        };
        if !source_input.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.rnastructure].input `{}` not found (resolved {})",
                    input.input.display(),
                    source_input.display()
                ),
            });
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "RNAstructure 6.4+ required; install via \
                       `conda install -c bioconda rnastructure` or \
                       download from https://rna.urmc.rochester.edu/RNAstructure.html"
                .into(),
        })?;

        // Compose the Fold invocation:
        //   `Fold <input> <output> -m <max_structures> -p <max_percent> -t <kelvin> [extras...]`
        // Output path is resolved against the workdir so Fold writes
        // the CT into our staging directory rather than next to the
        // input. Other Phase 17/18 single-binary adapters use the
        // same pattern.
        let output_in_workdir = workdir.join(&input.output);
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            source_input.into_os_string(),
            output_in_workdir.into_os_string(),
            OsString::from("-m"),
            OsString::from(input.max_structures.to_string()),
            OsString::from("-p"),
            OsString::from(input.max_percent.to_string()),
            OsString::from("-t"),
            OsString::from(format_temperature(input.temperature)),
        ];
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        // Stash the output filename so collect() can find the CT
        // without re-parsing the case.
        let environment: Vec<(OsString, OsString)> = vec![(
            OsString::from("VALENX_RNASTRUCTURE_OUTPUT"),
            OsString::from(input.output.as_os_str()),
        )];

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment,
            // Fold on a single moderate-length sequence runs in
            // seconds; long sequences (~1 kb+) with many suboptimals
            // can take minutes. 30 minutes covers the long tail.
            estimated_runtime: Some(Duration::from_secs(30 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting RNAstructure Fold", |line| {
            let mut hint = subprocess::Hint::default();
            // Fold's progress chatter is sparse:
            //   * "Reading the parameter files" — startup
            //   * "Energies calculated" — main step done
            //   * "Writing CT file" — wrap-up
            // Heuristics; mismatches just leave the spinner alone.
            if line.contains("Writing CT file") || line.contains("CT file written") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("Energies calculated") || line.contains("Calculating") {
                hint.progress = Some((60.0, line.to_string()));
            } else if line.contains("Reading the parameter") || line.contains("Reading sequence") {
                hint.progress = Some((20.0, line.to_string()));
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
        let output_rel = job
            .environment
            .iter()
            .find(|(k, _)| k == "VALENX_RNASTRUCTURE_OUTPUT")
            .map(|(_, v)| v.clone());
        let out_path = output_rel.map(|rel| job.workdir.join(std::path::PathBuf::from(&rel)));

        // Provenance: hash the staged CT output if present, else
        // case.toml.
        let case_hash_input = match &out_path {
            Some(p) if p.is_file() => p.clone(),
            _ => job.workdir.join("case.toml"),
        };
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "RNAstructure",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        if let Some(p) = out_path {
            if p.is_file() {
                artefacts.push(Artifact {
                    path: p,
                    kind: ArtifactKind::Native,
                    checksum: None,
                    label: "RNAstructure connectivity table".to_string(),
                });
            }
        }
        artefacts.sort_by(|a, b| a.path.cmp(&b.path));
        results.artifacts = artefacts;
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.rnastructure.fold"],
        }
    }
}

/// Format a Kelvin temperature for Fold's `-t` flag. Shortest round-
/// trip rendering — whole-number values render without trailing `.0`.
fn format_temperature(temp: f64) -> String {
    format!("{temp}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = RnaStructureAdapter::new().info();
        assert_eq!(info.id, "rnastructure");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "BSD-3-Clause");
        assert_eq!(info.display_name, "RNAstructure");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = RnaStructureAdapter::new().info();
        // RNAstructure's 6.4.x is the current stable; 7.0 reserves
        // room for an eventual major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(6, 4, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(7, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = RnaStructureAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.rnastructure.fold"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = RnaStructureAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    /// Round-10 H3 RED→GREEN: `output` flowed into
    /// `workdir.join(...)` with no validation. Hostile
    /// `output = "../etc/passwd"` is now rejected.
    #[test]
    fn prepare_rejects_output_path_traversal() {
        use valenx_test_utils::tempdir;
        let d = tempdir("rnastructure-output-trav");
        std::fs::write(d.join("seq.fa"), b">x\nACGU\n").unwrap();
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "rnastructure.fold"

[bio.rnastructure]
input          = "seq.fa"
output         = "../etc/passwd"
max_structures = 1
max_percent    = 10
temperature    = 310.15
"#,
        )
        .unwrap();
        let case = Case {
            id: "trav".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = RnaStructureAdapter::new()
            .prepare(&case, &workdir)
            .unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("[bio.rnastructure].output"),
            "expected [bio.rnastructure].output in error, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
