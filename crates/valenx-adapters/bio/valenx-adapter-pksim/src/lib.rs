//! # valenx-adapter-pksim
//!
//! Adapter for [PK-Sim](https://www.open-systems-pharmacology.org/) —
//! the Open Systems Pharmacology suite's whole-body physiologically
//! based pharmacokinetic (PBPK) simulation engine. PK-Sim consumes a
//! `.pksim5` XML project file describing compounds, formulations,
//! individuals or populations, and the simulation protocol; it
//! integrates the resulting ODE system and writes time-course
//! concentration tables plus a JSON metadata block.
//!
//! **Phase 45 — subprocess wrapper around the PK-Sim CLI.** The user
//! supplies the project file via `[bio.pksim]` in `case.toml`;
//! `prepare()` composes
//! `pksim --project <project> --output <basename> [extras...]`,
//! `run()` streams progress via the shared subprocess runner, and
//! `collect()` walks the workdir for the canonical `<basename>*.csv`,
//! `<basename>*.json`, and `*.log` outputs.

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

use crate::case_input::PkSimInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(PkSimAdapter::new())
}

pub struct PkSimAdapter;

impl PkSimAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for PkSimAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "pksim";
/// PK-Sim's binary candidates. The conda-forge / Linux release ships
/// a lowercase `pksim`; the upstream Windows build exposes the
/// `PKSim.CLI` executable. Either is sufficient.
const BINARIES: &[&str] = &["pksim", "PKSim.CLI"];

impl Adapter for PkSimAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "PK-Sim",
            // PK-Sim 11 is the modern CLI line (the OSP suite jumped
            // to 11.x in 2023); upper bound 13.0 reserves room for
            // the next two minor majors.
            version_range: VersionRange {
                min_inclusive: Version::new(11, 0, 0),
                max_exclusive: Version::new(13, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-2.0",
            docs_url: "https://docs.open-systems-pharmacology.org/",
            homepage_url: "https://www.open-systems-pharmacology.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `pksim --version` prints the OSP release banner;
                // older builds also accept `-v`. Try the long flag
                // first.
                let found_version = detect_tool_version_semver(&binary_path, &["--version", "-v"]);
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
                hint: "PK-Sim 11+ required; install via the Open \
                       Systems Pharmacology suite installer at \
                       https://www.open-systems-pharmacology.org/ \
                       or `conda install -c conda-forge pksim`"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = PkSimInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.pksim].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

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

        // Resolve the project file against the case directory via
        // `confined_join` — rejects absolute paths and `..` traversal
        // out of the case sandbox. Round-8 sibling-field sweep: the
        // stale "plain `case.path.join` is correct" comment was wrong;
        // PK-Sim reads the project file but a hostile case bundle
        // shouldn't be able to point `project` at `/etc/passwd`.
        let source_project = valenx_core::adapter_helpers::confined_join(&case.path, &input.project)?;
        if !source_project.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.pksim].project `{}` not found (resolved {})",
                    input.project.display(),
                    source_project.display()
                ),
            });
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "PK-Sim 11+ required; install via the Open \
                       Systems Pharmacology suite installer at \
                       https://www.open-systems-pharmacology.org/ \
                       or `conda install -c conda-forge pksim`"
                .into(),
        })?;

        // Compose `pksim --project <project> --output <basename>
        //          [extras...]`. PK-Sim writes `<basename>*.csv` and
        // `<basename>*.json` into the cwd (workdir).
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from("--project"),
            source_project.into_os_string(),
            OsString::from("--output"),
            OsString::from(&input.output_basename),
        ];
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Single-individual PBPK simulations finish in seconds;
            // population simulations over thousands of virtual
            // subjects can stretch to an hour or more. 4 hours
            // covers the long tail.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting PK-Sim", |line| {
            let mut hint = subprocess::Hint::default();
            // PK-Sim's CLI emits a startup banner ("PK-Sim ..."),
            // per-individual progress markers, and a completion
            // sentinel near the end of run.
            if line.contains("Simulation completed") || line.contains("simulation finished") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.starts_with("PK-Sim") {
                hint.progress = Some((5.0, line.to_string()));
            } else if line.contains("Individual ") || line.contains("Population ") {
                hint.progress = Some((50.0, line.to_string()));
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
        // descriptor. PK-Sim's project file isn't staged into the
        // workdir (read in place), so the case.toml is the closest
        // stable input fingerprint we have here.
        let case_hash_input = job.workdir.join("case.toml");
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "PK-Sim",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Read the output basename back from the case.toml so we
        // can prefix-filter `<basename>*.csv` / `<basename>*.json`
        // outputs.
        let basename = case_input::PkSimInput::from_case_dir(&job.workdir)
            .ok()
            .map(|i| i.output_basename);

        // Walk the workdir top-level. PK-Sim writes
        // `<basename>*.csv` (simulation tables), `<basename>*.json`
        // (metadata) and a `*.log` file.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-pksim", ?e, "workdir read failed");
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
            let stem = path
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("")
                .to_string();
            let stem_matches_basename = match basename.as_deref() {
                Some(b) => stem.starts_with(b),
                None => true,
            };
            match ext.as_deref() {
                Some("csv") => {
                    if !stem_matches_basename {
                        continue;
                    }
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Tabular,
                        checksum: None,
                        label: "PK-Sim simulation results".to_string(),
                    });
                }
                Some("json") => {
                    if !stem_matches_basename {
                        continue;
                    }
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Tabular,
                        checksum: None,
                        label: "PK-Sim metadata".to_string(),
                    });
                }
                Some("log") => {
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Log,
                        checksum: None,
                        label: "PK-Sim log".to_string(),
                    });
                }
                _ => continue,
            }
        }
        artefacts.sort_by(|a, b| a.path.cmp(&b.path));
        results.artifacts = artefacts;
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.pksim.simulate"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = PkSimAdapter::new().info();
        assert_eq!(info.id, "pksim");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "GPL-2.0");
        assert_eq!(info.display_name, "PK-Sim");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = PkSimAdapter::new().info();
        // PK-Sim 11 is the modern CLI line; 13.0 reserves room for
        // the next two majors.
        assert_eq!(info.version_range.min_inclusive, Version::new(11, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(13, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = PkSimAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.pksim.simulate"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = PkSimAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    #[test]
    fn prepare_rejects_traversal_project_path() {
        // Round-8 RED→GREEN: `[bio.pksim].project` now routes through
        // `confined_join`. Pre-fix a stale comment claimed plain
        // `case.path.join` was correct because PK-Sim "reads it in
        // place" — that's not a sandbox argument, just a staging one.
        //
        // We use `../etc/passwd` (relative traversal) for cross-platform
        // portability — see the bcftools test for the rationale.
        use valenx_test_utils::tempdir;
        let d = tempdir("pksim-traversal");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "pksim.simulate"

[bio.pksim]
project         = "../etc/passwd"
output_basename = "sim"
"#,
        )
        .unwrap();
        let case = Case {
            id: "pksim-traversal".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = PkSimAdapter::new().prepare(&case, &workdir).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("absolute") || msg.contains("escape") || msg.contains("`..`") || msg.contains("traversal"),
            "expected confined_join rejection on project, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
