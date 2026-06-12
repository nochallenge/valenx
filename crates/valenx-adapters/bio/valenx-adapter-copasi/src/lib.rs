//! # valenx-adapter-copasi
//!
//! Adapter for [COPASI](http://copasi.org/) — the COmplex PAthway
//! SImulator. The de-facto desktop suite for biochemical pathway and
//! ODE-based systems-biology models, descended from the Gepasi
//! lineage. The CLI binary is `CopasiSE` ("Self-Executing"): a
//! headless task runner that reads a COPASI native `.cps` model (or
//! an SBML `.xml`) and executes the simulation / scan / fitting tasks
//! defined inside.
//!
//! **Phase 32 — subprocess wrapper around `CopasiSE`.** The user
//! supplies a model file via `[bio.copasi]` in `case.toml`; an
//! optional `report = "report.csv"` adds `--save <path>` so the run
//! output lands at a known location, and `run_all = true` flips the
//! flag that asks `CopasiSE` to execute every task defined in the
//! file rather than just the primary one. `prepare()` composes
//! `CopasiSE [--save <report>] <model> [extras...]`; `run()` streams
//! through the shared subprocess runner.
//!
//! On `collect()` we surface the explicit report path when supplied,
//! and otherwise walk the workdir top-level for `.csv` / `.txt` files
//! (COPASI's tabular outputs).

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::fs;
use std::path::{Path, PathBuf};
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

use crate::case_input::CopasiInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(CopasiAdapter::new())
}

pub struct CopasiAdapter;

impl CopasiAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for CopasiAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "copasi";
/// COPASI's headless CLI binary. The capital `C-S-E` spelling is
/// canonical: `CopasiSE` = "COPASI Self-Executing".
const BINARIES: &[&str] = &["CopasiSE"];

impl Adapter for CopasiAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "COPASI",
            // COPASI 4.x is the long-running stable line; 4.40 is a
            // recent floor that ships SBML L3v2 + the task scheduler
            // every Phase 32 model relies on. 5.0 reserves room for
            // an eventual major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(4, 40, 0),
                max_exclusive: Version::new(5, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "Artistic-2.0",
            docs_url: "http://copasi.org/Support/User_Manual/",
            homepage_url: "http://copasi.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `CopasiSE --version` prints something like
                // "COPASI 4.40 (Build 287)" on stdout.
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
                hint: "COPASI 4.40+ required; download CopasiSE from \
                       http://copasi.org/Download/ or install via \
                       `conda install -c conda-forge copasi`"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = CopasiInput::from_case_dir(&case.path)?;

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

        // Resolve the model path against the case directory if
        // relative. Same convention as every other Phase 17/18 bio
        // adapter — `model = "pathway.cps"` next to `case.toml`.
        let source_model = if input.model.is_absolute() {
            input.model.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.model)?
        };
        if !source_model.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.copasi].model `{}` not found (resolved {})",
                    input.model.display(),
                    source_model.display()
                ),
            });
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "COPASI 4.40+ required; download CopasiSE from \
                       http://copasi.org/Download/ or install via \
                       `conda install -c conda-forge copasi`"
                .into(),
        })?;

        // Compose `CopasiSE [--save <report>] <model> [extras...]`.
        // `--save` redirects task output to a single file (as opposed
        // to whatever the model's report definitions imply). We keep
        // `--save` adjacent to the model so the line reads naturally
        // even when the user adds extras.
        let mut native_command: Vec<OsString> = vec![binary_path.into_os_string()];
        if let Some(report) = &input.report {
            native_command.push(OsString::from("--save"));
            native_command.push(OsString::from(report));
        }
        native_command.push(source_model.into_os_string());
        if input.run_all {
            // `--scheduled` runs every task with the "executable" flag
            // set in the file rather than just the first. Off by
            // default — typical COPASI files have one primary task.
            native_command.push(OsString::from("--scheduled"));
        }
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // ODE / scan / fitting on a typical pathway model finishes
            // in seconds-to-minutes; long parameter scans on stiff
            // systems can take an hour. 4 hours covers the long tail.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting COPASI", |line| {
            let mut hint = subprocess::Hint::default();
            // CopasiSE's stdout is sparse — task summaries and
            // success / failure markers. Lift the obvious milestones
            // to coarse UI ticks.
            if line.contains("Task finished") || line.contains("Simulation finished") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.to_ascii_lowercase().contains("error")
                || line.to_ascii_lowercase().contains("failed")
            {
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
        // Read the staged input back out so we can identify the
        // explicit report path the user pinned (if any). Failure to
        // re-read is non-fatal — collect still walks for tabular
        // outputs.
        let input = CopasiInput::from_case_dir(&job.workdir).ok();
        let explicit_report: Option<PathBuf> = input.as_ref().and_then(|i| {
            i.report.as_ref().map(|r| {
                if r.is_absolute() {
                    r.clone()
                } else {
                    job.workdir.join(r)
                }
            })
        });

        // Provenance: hash the explicit report if it landed; else
        // case.toml so the provenance block stays well-formed even
        // for partial / failed runs.
        let case_hash_input = match &explicit_report {
            Some(p) if p.is_file() => p.clone(),
            _ => job.workdir.join("case.toml"),
        };
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "COPASI",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        if let Some(p) = explicit_report {
            if p.is_file() {
                artefacts.push(Artifact {
                    path: p,
                    kind: ArtifactKind::Tabular,
                    checksum: None,
                    label: "COPASI report".to_string(),
                });
            }
        } else {
            // No pinned report — walk the workdir top-level for
            // `.csv` / `.txt` files. COPASI's report definitions
            // typically write tabular text into the run cwd; a single
            // pass is enough to surface them.
            let entries = match fs::read_dir(&job.workdir) {
                Ok(e) => e,
                Err(e) => {
                    tracing::warn!(target: "valenx-copasi", ?e, "workdir read failed");
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
                    Some("csv") | Some("txt") => {
                        (ArtifactKind::Tabular, "COPASI output".to_string())
                    }
                    _ => continue,
                };
                artefacts.push(Artifact {
                    path,
                    kind,
                    checksum: None,
                    label,
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
            ribbon_contributions: vec!["bio.copasi.simulate"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = CopasiAdapter::new().info();
        assert_eq!(info.id, "copasi");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "Artistic-2.0");
        assert_eq!(info.display_name, "COPASI");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = CopasiAdapter::new().info();
        assert_eq!(info.version_range.min_inclusive, Version::new(4, 40, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(5, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = CopasiAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.copasi.simulate"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = CopasiAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
