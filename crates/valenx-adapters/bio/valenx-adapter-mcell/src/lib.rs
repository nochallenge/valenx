//! # valenx-adapter-mcell
//!
//! Adapter for [MCell](https://mcell.org/) — the Salk Institute
//! (Stiles, Bartol) spatial stochastic cell-scale simulator. MCell
//! performs Monte Carlo diffusion of individual molecules through
//! realistic 3D subcellular geometries, with reactions on surfaces
//! and in volumes. The whole model is described in MDL ("Model
//! Description Language") files (`*.mdl`) and MCell writes per-time
//! reaction-data tables (`*.dat`), visualization data (`*.dx`), and a
//! run log into the workdir.
//!
//! **Phase 32.5 — subprocess wrapper around `mcell`.** The user
//! supplies the MDL file via `[bio.mcell]` in `case.toml`;
//! `prepare()` composes `mcell [-seed N] <mdl> [extras...]`,
//! `run()` streams progress via the shared subprocess runner, and
//! `collect()` walks the workdir for the canonical `*.dat`, `*.dx`,
//! and `*.log` outputs.

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

use crate::case_input::McellInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(McellAdapter::new())
}

pub struct McellAdapter;

impl McellAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for McellAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "mcell";
/// MCell's binary candidate. Source builds, conda-forge, and the
/// upstream installer all expose the canonical lowercase `mcell`.
const BINARIES: &[&str] = &["mcell"];

impl Adapter for McellAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "MCell",
            // MCell 4.0 is the current stable line (the 4.x rewrite
            // brought parallel/spatial extensions and the modern
            // CellBlender Python bindings). Upper bound 5.0 reserves
            // room for the next major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(4, 0, 0),
                max_exclusive: Version::new(5, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-2.0",
            docs_url: "https://mcell.org/documentation/",
            homepage_url: "https://mcell.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `mcell --version` prints the release banner; older
                // builds also accept `-v`. Try the long flag first.
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
                hint: "MCell 4.0+ required; install via \
                       `conda install -c conda-forge mcell`, the \
                       upstream installer at https://mcell.org/, \
                       or build from source"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = McellInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve the MDL file against the case directory if
        // relative. MCell reads it as the sole positional argument;
        // outputs are written into the cwd (workdir).
        let source_mdl = if input.mdl.is_absolute() {
            input.mdl.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(
            &case.path,
            &input.mdl,
        )?
        };
        if !source_mdl.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.mcell].mdl `{}` not found (resolved {})",
                    input.mdl.display(),
                    source_mdl.display()
                ),
            });
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "MCell 4.0+ required; install via \
                       `conda install -c conda-forge mcell`, the \
                       upstream installer at https://mcell.org/, \
                       or build from source"
                .into(),
        })?;

        // Compose `mcell [-seed N] <mdl> [extras...]`. MCell takes
        // the MDL file as the only positional argument; the optional
        // `-seed N` flag (two separate args, NOT `-seed=N`) overrides
        // the default deterministic seed (1); other CLI options
        // (`-quiet`, `-iterations N`, `-checkpoint_infile ...`) are
        // appended via `extra_args`.
        let mut native_command: Vec<OsString> = vec![binary_path.into_os_string()];
        if let Some(seed) = input.seed {
            native_command.push(OsString::from("-seed"));
            native_command.push(OsString::from(seed.to_string()));
        }
        native_command.push(source_mdl.into_os_string());
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // A small reaction-diffusion model with a few species
            // finishes in seconds; large multi-compartment models
            // with millions of molecules and surface chemistry can
            // run for many hours. 4 hours covers the long tail.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting MCell", |line| {
            let mut hint = subprocess::Hint::default();
            // MCell's stdout chatter: a startup banner ("MCell"),
            // "Iteration:" lines per simulated step, and a "Finished
            // ... iterations" / "Exiting." sentinel near end of run.
            if line.contains("Exiting.")
                || line.contains("Finished")
                || line.contains("simulation complete")
            {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.starts_with("MCell") {
                hint.progress = Some((5.0, line.to_string()));
            } else if line.contains("Iteration:") || line.contains("iteration") {
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
        // descriptor. MCell's output filenames are baked into the
        // MDL (`OUTPUT_BUFFER_SIZE`, `REACTION_DATA_OUTPUT { ... }`,
        // `VIZ_OUTPUT { ... }`), so we can't assume a single
        // fixed-name artifact for the prov hash.
        let case_hash_input = job.workdir.join("case.toml");
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "mcell",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk the workdir top-level. MCell's three canonical output
        // families are `*.dat` (per-step reaction data tables — the
        // primary observable output), `*.dx` (DREAMM/OpenDX
        // visualization data for CellBlender / pyMCell preview), and
        // `*.log` (run log).
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-mcell", ?e, "workdir read failed");
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
                Some("dat") => (ArtifactKind::Tabular, "MCell reaction data".to_string()),
                Some("dx") => (ArtifactKind::Native, "MCell visualization data".to_string()),
                Some("log") => (ArtifactKind::Log, "MCell log".to_string()),
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
            ribbon_contributions: vec!["bio.mcell.simulate"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = McellAdapter::new().info();
        assert_eq!(info.id, "mcell");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "GPL-2.0");
        assert_eq!(info.display_name, "MCell");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = McellAdapter::new().info();
        // MCell 4.0+ is the modern stable line; upper bound 5.0
        // reserves room for the next major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(4, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(5, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = McellAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.mcell.simulate"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = McellAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
