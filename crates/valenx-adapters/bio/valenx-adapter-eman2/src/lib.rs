//! # valenx-adapter-eman2
//!
//! Adapter for [EMAN2](https://eman2.org/) — Steve Ludtke's
//! broad-spectrum cryo-EM image-processing package. EMAN2 is the
//! "Swiss army knife" of single-particle cryo-EM: particle picking,
//! 2D classification, initial-model building, 3D refinement (CTF
//! corrected, with simultaneous tilt-pair handling), and a sprawling
//! Python toolkit (`e2*.py`) for everything in between.
//!
//! **Phase 36 — subprocess wrapper around `e2refine_easy.py`.** The
//! user supplies a particle list, an initial 3D model, an output
//! basename, a target resolution, and a point-group symmetry via
//! `[bio.eman2]` in `case.toml`. `prepare()` resolves the inputs and
//! composes the invocation; `run()` streams via the shared subprocess
//! runner. EMAN2 prints structured progress to stdout — initialisation
//! banner, per-iteration timing, "Done." on completion — which the
//! standard handler picks up for free.
//!
//! On `collect()` we walk the workdir for EMAN2's results directory —
//! `<output_basename>_<NN>/` — surfacing the 3D reconstructions
//! (`threed_*.hdf`) and the run log (`log.txt`).

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

use crate::case_input::Eman2Input;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(Eman2Adapter::new())
}

pub struct Eman2Adapter;

impl Eman2Adapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for Eman2Adapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "eman2";
/// EMAN2's "easy" refinement entry point. The full toolkit ships
/// dozens of `e2*.py` scripts; this one is the high-level driver
/// that orchestrates the others.
const BINARIES: &[&str] = &["e2refine_easy.py"];

impl Adapter for Eman2Adapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "EMAN2",
            // EMAN2's 2.99 line is the current pre-3.0 stable release
            // (the 2.x series has carried the project for over a
            // decade). 2.99.0 is the floor we test against; upper
            // bound 3.0 reserves room for the long-rumoured 3.x line.
            version_range: VersionRange {
                min_inclusive: Version::new(2, 99, 0),
                max_exclusive: Version::new(3, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "BSD-3-Clause",
            docs_url: "https://eman2.org/",
            homepage_url: "https://eman2.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `e2refine_easy.py --version` prints the EMAN2
                // version on stdout; the combined scanner picks it
                // up.
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
                hint: "EMAN2 2.99+ required; install via `conda install -c cryoem eman2` \
                       or follow the project install guide at https://eman2.org/"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = Eman2Input::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.eman2].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Resolve the particle stack against the case directory.
        let source_particles = if input.particles.is_absolute() {
            input.particles.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(
            &case.path,
            &input.particles,
        )?
        };
        if !source_particles.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.eman2].particles `{}` not found (resolved {})",
                    input.particles.display(),
                    source_particles.display()
                ),
            });
        }

        // Resolve the initial 3D model.
        let source_model = if input.model.is_absolute() {
            input.model.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(
            &case.path,
            &input.model,
        )?
        };
        if !source_model.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.eman2].model `{}` not found (resolved {})",
                    input.model.display(),
                    source_model.display()
                ),
            });
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "EMAN2 2.99+ required; install via `conda install -c cryoem eman2` \
                       or follow the project install guide at https://eman2.org/"
                .into(),
        })?;

        // Compose `e2refine_easy.py --input <particles> --model <model>
        // --path <output_basename> --targetres <res> --sym <sym>
        // --threads <N> [extras...]`.
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from("--input"),
            source_particles.into_os_string(),
            OsString::from("--model"),
            source_model.into_os_string(),
            OsString::from("--path"),
            OsString::from(&input.output_basename),
            OsString::from("--targetres"),
            OsString::from(format!("{}", input.target_resolution)),
            OsString::from("--sym"),
            OsString::from(&input.symmetry),
            OsString::from("--threads"),
            OsString::from(input.threads.to_string()),
        ];
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // EMAN2 refinement spans minutes (test datasets) to days
            // (full single-particle reconstructions). 24 hours is a
            // generous default.
            estimated_runtime: Some(Duration::from_secs(24 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting EMAN2", |line| {
            let mut hint = subprocess::Hint::default();
            // EMAN2's progress chatter on stdout: "Iteration N",
            // "Resolution = X angstroms", "Done." sentinel at the
            // end. Lift to coarse UI ticks.
            if line.contains("Done.") || line.contains("done.") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("Iteration") || line.contains("Resolution") {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.contains("Initializing") || line.contains("Reading") {
                hint.progress = Some((15.0, line.to_string()));
            }
            if line.contains("ERROR") || line.contains("Error:") {
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
        // Provenance: hash the workdir's case.toml as the canonical
        // input descriptor. EMAN2's actual outputs land under the
        // results subdirectory and are walked separately for the
        // artifact list.
        let case_hash_input = job.workdir.join("case.toml");
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "EMAN2",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk recursively. EMAN2 turns the user's `--path` argument
        // into a `<basename>_NN/` results directory under the workdir
        // and writes `threed_*.hdf` reconstructions plus a `log.txt`
        // run log inside it.
        walk_artifacts(&job.workdir, &mut artefacts);
        artefacts.sort_by(|a, b| a.path.cmp(&b.path));
        results.artifacts = artefacts;
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.eman2.refine"],
        }
    }
}

/// Recursively walk `dir` and push any matching EMAN2 artifact under
/// `out`. Bounded depth via the natural recursion against `is_dir`.
fn walk_artifacts(dir: &Path, out: &mut Vec<Artifact>) {
    let entries = match fs::read_dir(dir) {
        Ok(e) => e,
        Err(e) => {
            tracing::warn!(target: "valenx-eman2", ?e, "workdir read failed: {}", dir.display());
            return;
        }
    };
    for entry in entries.flatten() {
        let path = entry.path();
        if path.is_dir() {
            walk_artifacts(&path, out);
            continue;
        }
        if !path.is_file() {
            continue;
        }
        let name = match path.file_name().and_then(|s| s.to_str()) {
            Some(n) => n,
            None => continue,
        };
        let lower = name.to_ascii_lowercase();
        // EMAN2's filename conventions inside the results directory:
        //   threed_NN.hdf  — 3D reconstructions, one per iteration
        //   log.txt        — per-run log file
        let (kind, label) = if lower.starts_with("threed_") && lower.ends_with(".hdf") {
            (ArtifactKind::Native, "EMAN2 reconstruction".to_string())
        } else if lower == "log.txt" {
            (ArtifactKind::Log, "EMAN2 log".to_string())
        } else {
            continue;
        };
        out.push(Artifact {
            path,
            kind,
            checksum: None,
            label,
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = Eman2Adapter::new().info();
        assert_eq!(info.id, "eman2");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "BSD-3-Clause");
        assert_eq!(info.display_name, "EMAN2");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = Eman2Adapter::new().info();
        // EMAN2 2.99 is the current pre-3.0 stable line; upper
        // bound 3.0 reserves room for the next major.
        assert_eq!(info.version_range.min_inclusive, Version::new(2, 99, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(3, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = Eman2Adapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.eman2.refine"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = Eman2Adapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
