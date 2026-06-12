//! # valenx-adapter-cpptraj
//!
//! Adapter for [cpptraj](https://amber-md.github.io/cpptraj/) —
//! AmberTools' canonical trajectory analysis tool. cpptraj reads
//! Amber `.prmtop` / `.parm7` topologies plus
//! `.nc` / `.dcd` / `.mdcrd` trajectories, runs an analysis script
//! authored in cpptraj's domain language (`trajin`, `rms`,
//! `radgyr`, `hbond`, `volume`, `clustering`, ...), and writes
//! results into the workdir as `.dat` (per-frame tables), `.agr`
//! (XmGrace plot data), or `.gnu` (gnuplot scripts).
//!
//! **Phase 5.5 — subprocess wrapper around `cpptraj`.** The user
//! supplies a topology and an analysis script via `[bio.cpptraj]`
//! in `case.toml`; `prepare()` composes the
//! `cpptraj -p <topology> -i <script> [extras...]` invocation,
//! `run()` streams progress via the shared subprocess runner, and
//! `collect()` walks the workdir for the canonical `.dat` / `.agr`
//! / `.gnu` outputs.

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

use crate::case_input::CpptrajInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(CpptrajAdapter::new())
}

pub struct CpptrajAdapter;

impl CpptrajAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for CpptrajAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "cpptraj";
/// cpptraj's binary candidate. AmberTools, conda-forge, and
/// Bioconda all expose the canonical lowercase `cpptraj`.
const BINARIES: &[&str] = &["cpptraj"];

impl Adapter for CpptrajAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "cpptraj",
            // cpptraj 6.x is the modern stable line shipped with
            // AmberTools 23+ (2023). Upper bound 7.0 reserves room
            // for the next major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(6, 0, 0),
                max_exclusive: Version::new(7, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-3.0",
            docs_url: "https://amber-md.github.io/cpptraj/",
            homepage_url: "https://github.com/Amber-MD/cpptraj",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `cpptraj --version` prints a banner with the
                // release on stdout; older lines accept `-V` too.
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
                hint: "cpptraj 6.0+ required; install AmberTools via \
                       `conda install -c conda-forge ambertools`, \
                       `apt install ambertools`, or build from source"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = CpptrajInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve the analysis script against the case directory if
        // relative. cpptraj reads it via `-i <script>`; outputs are
        // written into the cwd (workdir).
        let source_script = if input.script.is_absolute() {
            input.script.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.script)?
        };
        if !source_script.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.cpptraj].script `{}` not found (resolved {})",
                    input.script.display(),
                    source_script.display()
                ),
            });
        }

        // Resolve the Amber topology similarly. We don't copy — the
        // topology can be tens of MB and cpptraj reads it by path.
        let source_topology = if input.topology.is_absolute() {
            input.topology.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.topology)?
        };
        if !source_topology.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.cpptraj].topology `{}` not found (resolved {})",
                    input.topology.display(),
                    source_topology.display()
                ),
            });
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "cpptraj 6.0+ required; install AmberTools via \
                       `conda install -c conda-forge ambertools`, \
                       `apt install ambertools`, or build from source"
                .into(),
        })?;

        // Compose `cpptraj -p <topology> -i <script> [extras...]`.
        // `-p` and `-i` are the canonical cpptraj flags for "topology
        // for the next trajectory" and "input script". Both must come
        // before any positional trajectory arguments the user added
        // via `extra_args` (e.g. `-y traj.nc`).
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from("-p"),
            source_topology.into_os_string(),
            OsString::from("-i"),
            source_script.into_os_string(),
        ];
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Trajectory analysis on a short MD run finishes in
            // seconds; multi-microsecond ensembles with hbond /
            // clustering routinely run for an hour or more. 4
            // hours is a generous default that covers the long
            // tail.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting cpptraj", |line| {
            let mut hint = subprocess::Hint::default();
            // cpptraj's stdout chatter: a startup banner ("CPPTRAJ:
            // Trajectory Analysis"), per-action progress ("Reading
            // 'traj.nc' as ..."), and a "TIME:" / "INTERNAL TIMING"
            // sentinel block at end-of-run.
            if line.contains("INTERNAL TIMING") || line.contains("TIME:") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("CPPTRAJ:") {
                hint.progress = Some((5.0, line.to_string()));
            } else if line.contains("Reading '") {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.contains("Error:") || line.contains("ERROR") {
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
        // descriptor. cpptraj's output filenames are baked into
        // the analysis script (`hbond out hbonds.dat`,
        // `radgyr out rg.dat`, ...), so we can't assume a single
        // fixed-name artifact for the prov hash.
        let case_hash_input = job.workdir.join("case.toml");
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "cpptraj",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk the workdir top-level. cpptraj's three canonical
        // output families are `.dat` (per-frame tables — the most
        // common), `.agr` (XmGrace plot scripts), and `.gnu`
        // (gnuplot output).
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-cpptraj", ?e, "workdir read failed");
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
                Some("dat") => (ArtifactKind::Tabular, "cpptraj analysis output".to_string()),
                Some("agr") => (ArtifactKind::Tabular, "cpptraj XmGrace output".to_string()),
                Some("gnu") => (ArtifactKind::Log, "cpptraj gnuplot output".to_string()),
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
            ribbon_contributions: vec!["bio.cpptraj.analyze"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = CpptrajAdapter::new().info();
        assert_eq!(info.id, "cpptraj");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "GPL-3.0");
        assert_eq!(info.display_name, "cpptraj");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = CpptrajAdapter::new().info();
        // cpptraj 6.x is the modern stable line shipped with
        // AmberTools 23+; upper bound 7.0 reserves room for the
        // next major.
        assert_eq!(info.version_range.min_inclusive, Version::new(6, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(7, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = CpptrajAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.cpptraj.analyze"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = CpptrajAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
