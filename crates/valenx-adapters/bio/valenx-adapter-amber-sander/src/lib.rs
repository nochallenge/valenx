//! # valenx-adapter-amber-sander
//!
//! Subprocess adapter for [AmberTools sander](https://ambermd.org/AmberTools.php)
//! — the OSS portion of AMBER's molecular-dynamics engine, sister to
//! NAMD / GROMACS / LAMMPS / OpenMM. sander itself is GPL-3.0; only
//! the GPU-accelerated `pmemd.cuda` requires the proprietary AMBER
//! license. **Phase 5.6 — straight subprocess wrapper, no academic
//! flagging needed.**
//!
//! sander reads three files at runtime:
//!   * an Amber topology (`-p`, `.prmtop` / `.parm7`),
//!   * starting coordinates (`-c`, `.inpcrd` / `.rst7`),
//!   * an mdin control deck (`-i`, `.in` / `.mdin`).
//!
//! Outputs share a single stem: `<basename>.out` (mdout text log),
//! `<basename>.rst` (restart coordinates), `<basename>.nc` (NetCDF
//! trajectory), and `<basename>.mdinfo` (periodic checkpoint).
//! `prepare()` builds
//! `sander -O -i <config> -p <topology> -c <coordinates> -o <basename>.out -r <basename>.rst -x <basename>.nc [extras...]`;
//! the `-O` flag overwrites pre-existing outputs (sander otherwise
//! exits with `Could not open <file> for writing`). `collect()` walks
//! the workdir for files whose stem starts with the configured
//! `output_basename` and surfaces the four output families.

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

use crate::case_input::SanderInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(SanderAdapter::new())
}

pub struct SanderAdapter;

impl SanderAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for SanderAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "sander";
/// sander's canonical binary name. AmberTools, conda-forge, and
/// Bioconda all expose it lowercase.
const BINARIES: &[&str] = &["sander"];

/// Hint surfaced when sander isn't on PATH. Mirrors cpptraj's
/// AmberTools install hint convention.
const INSTALL_HINT: &str = "sander 22+ required; install AmberTools via \
                            `conda install -c conda-forge ambertools`, \
                            `apt install ambertools`, or build from source";

impl Adapter for SanderAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "AmberTools sander",
            // AmberTools 22 (2022) is the floor we test against;
            // upper bound 26 reserves room for the next major
            // release line.
            version_range: VersionRange {
                min_inclusive: Version::new(22, 0, 0),
                max_exclusive: Version::new(26, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            // sander itself is GPL-3.0; only `pmemd.cuda` requires
            // the proprietary AMBER license. We label sander
            // honestly — no academic-license warning needed.
            tool_license: "GPL-3.0",
            docs_url: "https://ambermd.org/AmberTools.php",
            homepage_url: "https://ambermd.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `sander --version` is supported in modern releases;
                // `-V` worked on older lines. The combined detector
                // covers both.
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
                hint: INSTALL_HINT.into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = SanderInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.sander].output_basename",
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

        // Resolve the three input files against the case directory if
        // relative. We do not stage them into the workdir — Amber
        // topologies and trajectories can be hundreds of MB, and
        // sander reads them by path. Existence checks here surface
        // misconfiguration before launching the binary.
        let resolved_topology = if input.topology.is_absolute() {
            input.topology.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.topology)?
        };
        if !resolved_topology.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.sander].topology `{}` not found (resolved {})",
                    input.topology.display(),
                    resolved_topology.display()
                ),
            });
        }

        let resolved_coordinates = if input.coordinates.is_absolute() {
            input.coordinates.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.coordinates)?
        };
        if !resolved_coordinates.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.sander].coordinates `{}` not found (resolved {})",
                    input.coordinates.display(),
                    resolved_coordinates.display()
                ),
            });
        }

        let resolved_config = if input.config.is_absolute() {
            input.config.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.config)?
        };
        if !resolved_config.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.sander].config `{}` not found (resolved {})",
                    input.config.display(),
                    resolved_config.display()
                ),
            });
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: INSTALL_HINT.into(),
        })?;

        // Compose
        //   sander -O -i <config> -p <topology> -c <coordinates>
        //          -o <basename>.out -r <basename>.rst -x <basename>.nc
        //          [extras...]
        // `-O` overwrites pre-existing outputs; without it sander
        // refuses to clobber on re-runs. Outputs land in the workdir
        // (cwd) — bare filenames are intentional.
        let basename = &input.output_basename;
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from("-O"),
            OsString::from("-i"),
            resolved_config.into_os_string(),
            OsString::from("-p"),
            resolved_topology.into_os_string(),
            OsString::from("-c"),
            resolved_coordinates.into_os_string(),
            OsString::from("-o"),
            OsString::from(format!("{basename}.out")),
            OsString::from("-r"),
            OsString::from(format!("{basename}.rst")),
            OsString::from("-x"),
            OsString::from(format!("{basename}.nc")),
        ];
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // sander runs span seconds (single-step minimisation) to
            // multi-day production trajectories. 8 hours covers a
            // typical equilibration / short production batch; longer
            // jobs override via the executor.
            estimated_runtime: Some(Duration::from_secs(8 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting sander", |line| {
            let mut hint = subprocess::Hint::default();
            // sander's stdout / mdout chatter:
            //   * "          Setting new box info" — input loaded
            //   * "NSTEP =" — integrator stepping (frequent)
            //   * "wallclock() was called" — shutdown banner
            //   * "Error" / "ERROR" / "FATAL" — surface as warnings
            // Heuristics; mismatches just leave the spinner alone.
            if line.contains("wallclock() was called") || line.contains("Total wall time") {
                hint.progress = Some((95.0, line.trim().to_string()));
            } else if line.starts_with("NSTEP") || line.contains("NSTEP =") {
                hint.progress = Some((50.0, line.trim().to_string()));
            } else if line.contains("Setting new box info") || line.contains("Begin reading energy")
            {
                hint.progress = Some((10.0, line.trim().to_string()));
            } else if line.contains("FATAL") || line.contains("Error") || line.contains("ERROR") {
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
        // Re-parse the staged case.toml so we can filter outputs by
        // `output_basename`. Failure is non-fatal — fall back to
        // accepting every match (matches cpptraj's "best-effort
        // walk" stance).
        let basename = SanderInput::from_case_dir(&job.workdir)
            .ok()
            .map(|i| i.output_basename);

        // Provenance: hash the staged mdout if present (sander's
        // canonical run output), otherwise the case.toml. This keeps
        // the run-id deterministic across the partial / failed and
        // successful run branches.
        let case_hash_input = match basename.as_deref() {
            Some(b) => {
                let mdout = job.workdir.join(format!("{b}.out"));
                if mdout.is_file() {
                    mdout
                } else {
                    job.workdir.join("case.toml")
                }
            }
            None => job.workdir.join("case.toml"),
        };
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "sander",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);

        // Walk the workdir top-level for sander's four canonical
        // output families: mdout (`.out`), restart coordinates
        // (`.rst`), NetCDF trajectory (`.nc`), and the mdinfo
        // checkpoint (`.mdinfo`). When `output_basename` is known we
        // restrict matches to files whose stem starts with it; this
        // avoids picking up unrelated `.out` / `.rst` files the user
        // may have left in the workdir.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-sander", ?e, "workdir read failed");
                return Ok(results);
            }
        };
        let mut artefacts: Vec<Artifact> = Vec::new();
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
                Some("out") => (ArtifactKind::Log, "sander mdout".to_string()),
                Some("nc") => (ArtifactKind::Native, "sander NetCDF trajectory".to_string()),
                Some("rst") => (
                    ArtifactKind::Native,
                    "sander restart coordinates".to_string(),
                ),
                Some("mdinfo") => (ArtifactKind::Log, "sander mdinfo".to_string()),
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
        // Bio-specific Capability variants land in a follow-up task;
        // ribbon contributions are already enough for the registry to
        // surface the adapter.
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.sander.simulate"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = SanderAdapter::new().info();
        assert_eq!(info.id, "sander");
        assert_eq!(info.display_name, "AmberTools sander");
        assert_eq!(info.physics, &[Physics::Bio]);
        // sander itself is GPL-3.0; only pmemd.cuda is proprietary.
        assert_eq!(info.tool_license, "GPL-3.0");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = SanderAdapter::new().info();
        // AmberTools 22 (2022) is the floor we test against; upper
        // bound 26 reserves room for the next major release line.
        assert_eq!(info.version_range.min_inclusive, Version::new(22, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(26, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = SanderAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.sander.simulate"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = SanderAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
