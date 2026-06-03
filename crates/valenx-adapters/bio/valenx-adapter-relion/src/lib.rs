//! # valenx-adapter-relion
//!
//! Adapter for [RELION](https://relion.readthedocs.io/) — Sjors
//! Scheres' REgularised LIkelihood OptimisatioN suite. RELION drives
//! the Bayesian-inference 3D reconstruction pipeline at the heart of
//! modern single-particle cryo-EM: particle classification, 3D
//! refinement, CTF estimation, and post-processing. It is the de-facto
//! workhorse in cryo-EM facilities worldwide.
//!
//! **Phase 36 — subprocess wrapper around `relion_refine`.** The user
//! supplies a particle stack (STAR), an initial reference (MRC), and
//! an output basename via `[bio.relion]` in `case.toml`. `prepare()`
//! resolves the inputs, picks the single-process or MPI binary based
//! on `mpi_procs`, and composes the invocation. `run()` streams via
//! the shared subprocess runner.
//!
//! ## MPI dispatch
//!
//! `mpi_procs == 1` invokes `relion_refine` directly (the
//! single-process binary). `mpi_procs > 1` switches to
//! `mpirun -n <N> relion_refine_mpi` so MPI-aware ranks coordinate
//! across nodes / GPUs. If `mpirun` isn't on PATH for the multi-rank
//! case the adapter returns a friendly install hint pointing at
//! OpenMPI / MPICH.
//!
//! On `collect()` we walk the workdir for RELION's prefix-named
//! outputs: `<output_basename>*_class*.mrc` reconstructions,
//! `<output_basename>*_data.star` particle assignments, and
//! `<output_basename>*_model.star` model summaries.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::fs;
use std::path::Path;
use std::time::Duration;

use semver::Version;

use valenx_core::{
    adapter_helpers::{confined_join, detect_tool_version_semver, find_on_path, live_provenance},
    error::RunPhase,
    subprocess, Adapter, AdapterError, AdapterInfo, Capabilities, Case, LicenseMode, Physics,
    PreparedJob, ProbeReport, RunContext, RunReport, VersionRange,
};
use valenx_fields::{
    artifact::{Artifact, ArtifactKind},
    Results,
};

use crate::case_input::RelionInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(RelionAdapter::new())
}

pub struct RelionAdapter;

impl RelionAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for RelionAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "relion";

/// Upper bound we accept for `mpi_procs`. A typo like `99999` would
/// otherwise translate straight into `mpirun -n 99999` and OOM the
/// host before even reaching the solver — 256 is more than any
/// realistic single-node RELION refinement, and most multi-rank runs
/// land between 4 and 64. Tune up if a future supercomputer-class
/// host genuinely needs more.
const MAX_MPI_PROCS: u32 = 256;
/// Single-process RELION refine binary candidates.
const BINARIES: &[&str] = &["relion_refine"];
/// MPI-enabled RELION refine binary candidates. RELION ships these as
/// separate binaries (suffixed `_mpi`) so the launcher knows which
/// transport to use.
const MPI_BINARIES: &[&str] = &["relion_refine_mpi"];
/// MPI launcher candidates. `mpirun` is the canonical OpenMPI / MPICH
/// launcher; missing it on the multi-rank path is a hard error.
const MPI_LAUNCHERS: &[&str] = &["mpirun"];

impl Adapter for RelionAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "RELION",
            // RELION 4.0 is the current stable line (the long-running
            // 3.1 series is the predecessor). 4.0 is the floor we test
            // against; upper bound 6.0 reserves room for the next
            // major.
            version_range: VersionRange {
                min_inclusive: Version::new(4, 0, 0),
                max_exclusive: Version::new(6, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-2.0",
            docs_url: "https://relion.readthedocs.io/",
            homepage_url: "https://relion.readthedocs.io/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `relion_refine --version` prints the version string
                // on stdout; the combined scanner picks it up.
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
                hint: "RELION 4.0+ required; install via the project \
                       instructions at https://relion.readthedocs.io/ \
                       or `conda install -c conda-forge relion`"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = RelionInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.relion].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        // Bound `mpi_procs`. A typo like `99999` would otherwise be
        // handed straight to `mpirun -n 99999` and OOM the host before
        // even reaching `relion_refine_mpi`. Reject before staging.
        if input.mpi_procs > MAX_MPI_PROCS {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.relion].mpi_procs = {} exceeds the safety cap of \
                     {MAX_MPI_PROCS}; lower the value or raise `MAX_MPI_PROCS` \
                     in the adapter source if your host genuinely supports it",
                    input.mpi_procs
                ),
            });
        }

        fs::create_dir_all(workdir)?;

        // Resolve the particle-stack STAR against the case directory
        // via `confined_join` — a shared case bundle should not be
        // able to point `particles` or `reference` at arbitrary host
        // files (sandbox requirement).
        let source_particles = confined_join(&case.path, &input.particles)?;
        if !source_particles.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.relion].particles `{}` not found (resolved {})",
                    input.particles.display(),
                    source_particles.display()
                ),
            });
        }

        // Resolve the initial reference MRC.
        let source_reference = confined_join(&case.path, &input.reference)?;
        if !source_reference.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.relion].reference `{}` not found (resolved {})",
                    input.reference.display(),
                    source_reference.display()
                ),
            });
        }

        // Build the command. Single-process: `relion_refine ...`.
        // Multi-rank: `mpirun -n <N> relion_refine_mpi ...`.
        let mut native_command: Vec<OsString> = Vec::new();
        if input.mpi_procs > 1 {
            // Multi-rank: prepend `mpirun -n <N>`. If mpirun is missing
            // we surface a helpful install-hint InvalidCase rather than
            // a cryptic ToolNotInstalled — the user already has RELION
            // installed and just needs the MPI launcher.
            let mpi_launcher =
                find_on_path(MPI_LAUNCHERS).ok_or_else(|| AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.relion].mpi_procs = {} requires `mpirun` on PATH; \
                         install OpenMPI (`apt install openmpi-bin`, \
                         `brew install open-mpi`) or MPICH (`apt install mpich`) \
                         to enable multi-rank RELION runs",
                        input.mpi_procs
                    ),
                })?;
            let mpi_binary =
                find_on_path(MPI_BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
                    name: INFO_ID,
                    hint: "RELION's MPI binary `relion_refine_mpi` not found on PATH; \
                           a single-process RELION install (no MPI build) cannot serve \
                           a multi-rank job"
                        .into(),
                })?;
            native_command.push(mpi_launcher.into_os_string());
            native_command.push(OsString::from("-n"));
            native_command.push(OsString::from(input.mpi_procs.to_string()));
            native_command.push(mpi_binary.into_os_string());
        } else {
            let binary_path =
                find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
                    name: INFO_ID,
                    hint: "RELION 4.0+ required; install via the project \
                           instructions at https://relion.readthedocs.io/ \
                           or `conda install -c conda-forge relion`"
                        .into(),
                })?;
            native_command.push(binary_path.into_os_string());
        }

        // Common arg trail. RELION's CLI accepts long-form flags;
        // `--i` is the particle stack, `--ref` is the initial volume,
        // `--o` is the output basename (RELION prepends it to every
        // result file), `--angpix` is the pixel size, `--j` is the
        // OpenMP thread count per MPI rank.
        native_command.push(OsString::from("--i"));
        native_command.push(source_particles.into_os_string());
        native_command.push(OsString::from("--ref"));
        native_command.push(source_reference.into_os_string());
        native_command.push(OsString::from("--o"));
        native_command.push(OsString::from(&input.output_basename));
        native_command.push(OsString::from("--angpix"));
        native_command.push(OsString::from(format!("{}", input.angpix)));
        native_command.push(OsString::from("--j"));
        native_command.push(OsString::from(input.threads.to_string()));
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // RELION 3D refinement runs span minutes (small synthetic
            // datasets) to days (full single-particle datasets). 24
            // hours is a generous default.
            estimated_runtime: Some(Duration::from_secs(24 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting RELION", |line| {
            let mut hint = subprocess::Hint::default();
            // RELION's progress chatter: "Iteration N of M" lines,
            // "Estimated remaining time" timing reports, and a
            // "Done!" sentinel at the end. Lift to coarse UI ticks.
            if line.contains("Done!") || line.contains("done!") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("Iteration") {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.contains("Estimating") || line.contains("Initialising") {
                hint.progress = Some((15.0, line.to_string()));
            }
            if line.contains("ERROR") || line.contains("ERR ") {
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
        // input descriptor. RELION's actual outputs live under the
        // run-prefix directory and are walked separately for the
        // artifact list.
        let case_hash_input = job.workdir.join("case.toml");
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "RELION",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk recursively. RELION's `--o <prefix>` may include a
        // subdirectory component (e.g. `Refine3D/run`); the artifacts
        // land in `<workdir>/Refine3D/run_class*.mrc` etc.
        walk_artifacts(&job.workdir, &mut artefacts);
        artefacts.sort_by(|a, b| a.path.cmp(&b.path));
        results.artifacts = artefacts;
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.relion.refine"],
        }
    }
}

/// Recursively walk `dir` and push any matching RELION artifact under
/// `out`. Bounded depth via the natural recursion against `is_dir`;
/// we stop descending into anything that isn't a real directory entry.
fn walk_artifacts(dir: &Path, out: &mut Vec<Artifact>) {
    let entries = match fs::read_dir(dir) {
        Ok(e) => e,
        Err(e) => {
            tracing::warn!(target: "valenx-relion", ?e, "workdir read failed: {}", dir.display());
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

        // RELION's filename conventions:
        //   <prefix>_class001.mrc, <prefix>_class002.mrc, ...  — class averages / volumes
        //   <prefix>_data.star                                 — particle assignments
        //   <prefix>_model.star                                — model summary
        //
        // We classify by suffix-pattern rather than the user's exact
        // basename so subdirectories named after the basename are
        // still recognised.
        let lower = name.to_ascii_lowercase();
        let (kind, label) = if lower.contains("_class") && lower.ends_with(".mrc") {
            (ArtifactKind::Native, "RELION class average".to_string())
        } else if lower.ends_with("_data.star") {
            (
                ArtifactKind::Tabular,
                "RELION particle assignments".to_string(),
            )
        } else if lower.ends_with("_model.star") {
            (ArtifactKind::Log, "RELION model summary".to_string())
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
        let info = RelionAdapter::new().info();
        assert_eq!(info.id, "relion");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "GPL-2.0");
        assert_eq!(info.display_name, "RELION");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = RelionAdapter::new().info();
        // RELION 4.0 is the current stable line; upper bound 6.0
        // reserves room for the next major.
        assert_eq!(info.version_range.min_inclusive, Version::new(4, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(6, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = RelionAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.relion.refine"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = RelionAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
