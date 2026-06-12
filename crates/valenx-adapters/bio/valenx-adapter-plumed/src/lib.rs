//! # valenx-adapter-plumed
//!
//! Adapter for [PLUMED](https://www.plumed.org/) — the open-source
//! enhanced-sampling and free-energy plug-in that wraps every major
//! MD engine (GROMACS, LAMMPS, AMBER, NAMD, OpenMM). PLUMED defines
//! collective variables (RMSD, dihedrals, distances, contact maps),
//! biases (metadynamics, well-tempered metad, umbrella sampling,
//! ABF), and a reweighting framework that turns biased trajectories
//! back into unbiased free-energy surfaces.
//!
//! **Phase 5.5 — subprocess wrapper around `plumed driver`.** The
//! `plumed driver` sub-command runs PLUMED standalone over a pre-
//! computed trajectory: read frames from `--mf_xtc <traj>`, evaluate
//! the collective variables defined in `--plumed <plumed.dat>`,
//! write COLVAR / bias / HILLS files into the workdir.
//!
//! The user supplies the PLUMED input file plus an XTC trajectory
//! via `[bio.plumed]` in `case.toml`; `prepare()` composes the
//! `plumed driver` invocation, `run()` streams progress via the
//! shared subprocess runner. `collect()` walks the workdir for the
//! canonical `<output_basename>*.dat` (COLVAR) and `*.bias` files.

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

use crate::case_input::PlumedInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(PlumedAdapter::new())
}

pub struct PlumedAdapter;

impl PlumedAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for PlumedAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "plumed";
/// PLUMED's binary candidate. Source builds, conda-forge, and
/// Bioconda all expose the canonical lowercase `plumed` launcher.
const BINARIES: &[&str] = &["plumed"];

impl Adapter for PlumedAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "PLUMED",
            // PLUMED 2.9 (2023) is the modern stable line — the
            // `driver` sub-command, the metadynamics / OPES bias
            // family, and the Python interface we lean on are all
            // mature there. Upper bound 3.0 reserves room for the
            // long-promised next major.
            version_range: VersionRange {
                min_inclusive: Version::new(2, 9, 0),
                max_exclusive: Version::new(3, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "LGPL-3.0",
            docs_url: "https://www.plumed.org/doc",
            homepage_url: "https://www.plumed.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `plumed --version` (and bare `plumed info --version`
                // for older lines) prints the release on stdout. The
                // generic detector tries both forms so we cover the
                // 2.9 launcher and any conda-packaged variants.
                let found_version =
                    detect_tool_version_semver(&binary_path, &["--version", "info --version"]);
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
                hint: "PLUMED 2.9+ required; install via \
                       `conda install -c conda-forge plumed`, \
                       `apt install plumed`, or build from source"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = PlumedInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.plumed].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Resolve the PLUMED input file against the case directory if
        // relative. `plumed driver` reads it via `--plumed <path>`;
        // PLUMED writes COLVAR / bias outputs into the cwd, which is
        // the workdir.
        let source_dat = if input.plumed_dat.is_absolute() {
            input.plumed_dat.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.plumed_dat)?
        };
        if !source_dat.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.plumed].plumed_dat `{}` not found (resolved {})",
                    input.plumed_dat.display(),
                    source_dat.display()
                ),
            });
        }

        // Resolve the trajectory similarly. Trajectories are commonly
        // tens of GB so we don't copy — PLUMED reads them by path.
        let source_traj = if input.trajectory.is_absolute() {
            input.trajectory.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.trajectory)?
        };
        if !source_traj.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.plumed].trajectory `{}` not found (resolved {})",
                    input.trajectory.display(),
                    source_traj.display()
                ),
            });
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "PLUMED 2.9+ required; install via \
                       `conda install -c conda-forge plumed`, \
                       `apt install plumed`, or build from source"
                .into(),
        })?;

        // Compose `plumed driver --plumed <dat> --mf_xtc <traj>
        // --kt <kt> [extras...]`. `--mf_xtc` is the GROMACS XTC
        // flag; users running DCD / TRR trajectories can override
        // by passing `--mf_dcd <path>` via `extra_args` (the bare
        // `--mf_xtc` arg here is then ignored as redundant).
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from("driver"),
            OsString::from("--plumed"),
            source_dat.into_os_string(),
            OsString::from("--mf_xtc"),
            source_traj.into_os_string(),
            OsString::from("--kt"),
            OsString::from(format!("{}", input.kt)),
        ];
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // PLUMED driver runs are single-pass over the trajectory
            // and finish in seconds for short trajectories, minutes
            // for production runs. 2 hours is a generous default that
            // covers reweighting over multi-microsecond trajectories.
            estimated_runtime: Some(Duration::from_secs(2 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting PLUMED", |line| {
            let mut hint = subprocess::Hint::default();
            // PLUMED's stderr chatter on a `driver` run: a startup
            // banner ("PLUMED: PLUMED is starting"), per-frame
            // progress when verbose, and a "PLUMED: Cycles..." or
            // "PLUMED: Finishing" sentinel at end-of-run.
            if line.contains("PLUMED: Finishing") || line.contains("Cycles total") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("PLUMED: PLUMED is starting") {
                hint.progress = Some((5.0, line.to_string()));
            } else if line.contains("PLUMED error") || line.contains("ERROR") {
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
        // descriptor. PLUMED's COLVAR filenames are baked into the
        // `plumed.dat` script itself (PRINT FILE=...), so we can't
        // assume a single fixed-name artifact for the prov hash.
        let case_hash_input = job.workdir.join("case.toml");
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "PLUMED",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk the workdir top-level. PLUMED writes COLVAR files as
        // `<basename>*.dat` and bias / HILLS surfaces as
        // `<basename>*.bias`. Restrict to outputs whose stem starts
        // with the configured `output_basename` so unrelated `.dat`
        // files (e.g. the user's input `plumed.dat` if staged) don't
        // pollute the artefact list.
        let basename = read_output_basename(&job.workdir);
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-plumed", ?e, "workdir read failed");
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
            let stem_ok = match basename.as_deref() {
                Some(b) => stem.starts_with(b),
                None => true,
            };
            if !stem_ok {
                continue;
            }
            let (kind, label) = match ext.as_deref() {
                Some("dat") => (ArtifactKind::Tabular, "PLUMED COLVAR output".to_string()),
                Some("bias") => (ArtifactKind::Tabular, "PLUMED bias".to_string()),
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
            ribbon_contributions: vec!["bio.plumed.analyze"],
        }
    }
}

/// Re-read the `[bio.plumed].output_basename` from a staged
/// `case.toml` for collect()-time output filtering. Returns `None`
/// when the case.toml is missing or unparseable — collect() then
/// accepts everything (best-effort).
fn read_output_basename(workdir: &Path) -> Option<String> {
    // Round-23 sweep: bound staged case.toml at MAX_PROJECT_FILE_BYTES.
    let text = valenx_core::io_caps::read_capped_to_string(
        &workdir.join("case.toml"),
        valenx_core::project::loader::MAX_PROJECT_FILE_BYTES as usize,
    )
    .ok()?;
    let parsed: toml::Value = toml::from_str(&text).ok()?;
    parsed
        .get("bio")?
        .get("plumed")?
        .get("output_basename")?
        .as_str()
        .map(str::to_string)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = PlumedAdapter::new().info();
        assert_eq!(info.id, "plumed");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "LGPL-3.0");
        assert_eq!(info.display_name, "PLUMED");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = PlumedAdapter::new().info();
        // PLUMED 2.9+ is the modern stable line; upper bound 3.0
        // reserves room for the next major.
        assert_eq!(info.version_range.min_inclusive, Version::new(2, 9, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(3, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = PlumedAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.plumed.analyze"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = PlumedAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
