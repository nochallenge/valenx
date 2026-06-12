//! # valenx-adapter-smoldyn
//!
//! Adapter for [Smoldyn](https://www.smoldyn.org/) — Steve Andrews'
//! spatial stochastic reaction-diffusion simulator. Particles diffuse
//! in 1D / 2D / 3D continuous space, react with each other and with
//! geometric surfaces according to user-defined chemistry, and the
//! whole simulation is described in a single plain-text configuration
//! file. Smoldyn writes per-time-step tables (`*.txt`, `*.dat`) and
//! a run log into the workdir.
//!
//! **Phase 32.5 — subprocess wrapper around `smoldyn`.** The user
//! supplies the config file via `[bio.smoldyn]` in `case.toml`;
//! `prepare()` composes the `smoldyn <config> [extras...]` invocation,
//! `run()` streams progress via the shared subprocess runner, and
//! `collect()` walks the workdir for the canonical `*.txt`, `*.dat`,
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

use crate::case_input::SmoldynInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(SmoldynAdapter::new())
}

pub struct SmoldynAdapter;

impl SmoldynAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for SmoldynAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "smoldyn";
/// Smoldyn's binary candidate. Source builds, conda-forge, and the
/// upstream installer all expose the canonical lowercase `smoldyn`.
const BINARIES: &[&str] = &["smoldyn"];

impl Adapter for SmoldynAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "Smoldyn",
            // Smoldyn 2.7x is the modern stable line (2.70 landed in
            // 2021 with the C++17 rewrite of the I/O layer). Upper
            // bound 3.0 reserves room for the next major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(2, 70, 0),
                max_exclusive: Version::new(3, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "LGPL-2.1",
            docs_url: "https://www.smoldyn.org/SmoldynManual.pdf",
            homepage_url: "https://www.smoldyn.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `smoldyn --version` prints the release banner; older
                // builds also accept `-V`. Try the long flag first.
                let found_version = detect_tool_version_semver(&binary_path, &["--version", "-V"]);
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
                hint: "Smoldyn 2.70+ required; install via \
                       `conda install -c conda-forge smoldyn`, the \
                       upstream installer at https://www.smoldyn.org/, \
                       or build from source"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = SmoldynInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve the config file against the case directory if
        // relative. Smoldyn reads it as the sole positional argument;
        // outputs are written into the cwd (workdir).
        let source_config = if input.config.is_absolute() {
            input.config.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.config)?
        };
        if !source_config.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.smoldyn].config `{}` not found (resolved {})",
                    input.config.display(),
                    source_config.display()
                ),
            });
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "Smoldyn 2.70+ required; install via \
                       `conda install -c conda-forge smoldyn`, the \
                       upstream installer at https://www.smoldyn.org/, \
                       or build from source"
                .into(),
        })?;

        // Compose `smoldyn <config> [extras...]`. Smoldyn takes the
        // config file as the only positional argument; everything
        // else (`-w` to suppress warnings, `-q` for quiet, `--define`
        // to override config variables, ...) is a flag the user
        // appends via `extra_args`.
        let mut native_command: Vec<OsString> =
            vec![binary_path.into_os_string(), source_config.into_os_string()];
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // A small reaction-diffusion model with ~1k particles
            // finishes in seconds; multi-million-particle simulations
            // with intricate surface chemistry routinely run for an
            // hour or more. 4 hours covers the long tail.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting Smoldyn", |line| {
            let mut hint = subprocess::Hint::default();
            // Smoldyn's stdout chatter: a startup banner ("Smoldyn"),
            // per-step progress when graphics are off ("Time: ..."),
            // and a "Simulation complete" sentinel near the end of run.
            if line.contains("Simulation complete") || line.contains("simulation complete") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.starts_with("Smoldyn") {
                hint.progress = Some((5.0, line.to_string()));
            } else if line.contains("Time:") || line.contains("time =") {
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
        // descriptor. Smoldyn's output filenames are baked into the
        // config file (`output_files data.txt`, `cmd N output ...`),
        // so we can't assume a single fixed-name artifact for the
        // prov hash.
        let case_hash_input = job.workdir.join("case.toml");
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "smoldyn",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk the workdir top-level. Smoldyn's three canonical
        // output families are `*.txt` (per-step output tables — the
        // most common), `*.dat` (alternative tabular data), and
        // `*.log` (run log).
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-smoldyn", ?e, "workdir read failed");
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
                Some("txt") => (ArtifactKind::Tabular, "Smoldyn output table".to_string()),
                Some("dat") => (ArtifactKind::Tabular, "Smoldyn data".to_string()),
                Some("log") => (ArtifactKind::Log, "Smoldyn log".to_string()),
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
            ribbon_contributions: vec!["bio.smoldyn.simulate"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = SmoldynAdapter::new().info();
        assert_eq!(info.id, "smoldyn");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "LGPL-2.1");
        assert_eq!(info.display_name, "Smoldyn");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = SmoldynAdapter::new().info();
        // Smoldyn 2.70+ is the modern stable line; upper bound 3.0
        // reserves room for the next major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(2, 70, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(3, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = SmoldynAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.smoldyn.simulate"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = SmoldynAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
