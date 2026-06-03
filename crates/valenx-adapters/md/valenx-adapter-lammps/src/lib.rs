//! # valenx-adapter-lammps
//!
//! Subprocess adapter for LAMMPS. **Phase 5 — live for NVE
//! Lennard-Jones demos.** Uses the typed `LammpsInput` + input-deck
//! writer; extends cleanly to more ensembles / potentials by
//! growing the enums rather than rewriting the generator.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;
pub mod input_writer;
pub mod log_parser;

use std::ffi::OsString;
use std::fs;
use std::path::Path;
use std::time::Duration;

use semver::Version;

use valenx_core::{
    adapter_helpers::{find_on_path, first_workdir_match},
    error::RunPhase,
    subprocess, Adapter, AdapterError, AdapterInfo, Capabilities, Capability, Case, LicenseMode,
    Physics, PreparedJob, ProbeReport, RunContext, RunReport, VersionRange,
};
use valenx_fields::{
    artifact::{Artifact, ArtifactKind},
    Results,
};

use crate::case_input::{LammpsInput, Potential};
use crate::input_writer::{DEFAULT_INPUT_FILENAME, DUMP_FILENAME, THERMO_FILENAME};

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(LammpsAdapter::new())
}

pub struct LammpsAdapter;

impl LammpsAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for LammpsAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "lammps";
const BINARIES: &[&str] = &["lmp", "lmp_serial", "lmp_mpi"];

impl Adapter for LammpsAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "LAMMPS",
            version_range: VersionRange {
                min_inclusive: Version::new(2023, 8, 0),
                max_exclusive: Version::new(2030, 1, 0),
            },
            physics: &[Physics::MolecularDynamics],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-2.0-only",
            docs_url: "https://docs.lammps.org/",
            homepage_url: "https://www.lammps.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                let found_version = valenx_core::adapter_helpers::detect_tool_version_semver(
                    &binary_path,
                    &["-h", "--help", "-help"],
                );
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
                hint: "LAMMPS required; install via `conda install -c conda-forge lammps` \
                       or from source (lmp / lmp_serial on PATH)"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let (_header, input) = LammpsInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Stage external data / potential files into the workdir so
        // relative paths in the deck resolve.
        stage_external_files(&case.path, workdir, &input)?;

        let deck_path = workdir.join(DEFAULT_INPUT_FILENAME);
        input_writer::write_to_file(&input, &deck_path)?;

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "no `lmp` on PATH; install LAMMPS".into(),
        })?;

        // `lmp -in in.lammps` is the canonical invocation.
        let native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from("-in"),
            OsString::from(DEFAULT_INPUT_FILENAME),
        ];

        // Extremely rough — 10 µs per atom per step, 1k atoms, N
        // steps → bounded by N * 10 µs.
        let estimated_runtime = Some(Duration::from_millis(
            input.run_steps.saturating_mul(1).max(100),
        ));

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            estimated_runtime,
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting LAMMPS", |line| {
            let mut hint = subprocess::Hint::default();
            if let Some(pct) = lammps_progress_hint(line) {
                hint.progress = Some((pct, line.to_string()));
            }
            if line.contains("ERROR") || line.contains("WARNING") {
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
        // case_path: the canonical input deck the writer always emits
        // as `in.lammps`. mesh_path: any external `read_data` file
        // staged into the workdir (LAMMPS data files conventionally
        // use `.data` / `.lmp`).
        let case_path = job.workdir.join(DEFAULT_INPUT_FILENAME);
        let mesh_path = first_workdir_match(&job.workdir, &["data", "lmp"]);
        let prov = valenx_core::adapter_helpers::live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "LAMMPS",
            "unknown",
            &case_path,
            mesh_path.as_deref(),
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);

        // Known output files.
        let known = [
            (
                DUMP_FILENAME,
                ArtifactKind::Tabular,
                "LAMMPS trajectory dump",
            ),
            (THERMO_FILENAME, ArtifactKind::Tabular, "LAMMPS thermo log"),
            ("log.lammps", ArtifactKind::Log, "LAMMPS log"),
            (
                DEFAULT_INPUT_FILENAME,
                ArtifactKind::Other,
                "LAMMPS input deck",
            ),
        ];
        for (name, kind, label) in known {
            let path = job.workdir.join(name);
            if path.is_file() {
                results.artifacts.push(Artifact {
                    path,
                    kind,
                    checksum: None,
                    label: label.to_string(),
                });
            }
        }

        results.artifacts.sort_by(|a, b| a.path.cmp(&b.path));
        // Parse the thermo block out of log.lammps and surface every
        // (column, step) pair as a ScalarRecord. Step counts as a
        // TimeKey::Iteration so the report layer can chart energy /
        // temperature / pressure / volume vs step. Failures are
        // skipped silently — the artifact stays listed.
        let log_path = job.workdir.join("log.lammps");
        // Round-23 named finding: bound the thermo-log read at
        // MAX_LAMMPS_LOG_BYTES (256 MiB) — pre-fix a hostile or
        // runaway log would slurp into memory before the parser
        // surfaced any thermo column.
        if let Ok(text) = valenx_core::io_caps::read_capped_to_string(
            &log_path,
            valenx_core::io_caps::MAX_LAMMPS_LOG_BYTES as usize,
        ) {
            let series = log_parser::parse_log(&text);
            for record in log_parser::to_canonical_scalars(&series) {
                results.scalars.insert(record);
            }
        }
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: vec![Capability::MdClassical],
            ribbon_contributions: vec!["md.lammps.nvt", "md.lammps.npt", "md.lammps.minimize"],
        }
    }
}

/// Coarse progress hints for LAMMPS stdout banners.
fn lammps_progress_hint(line: &str) -> Option<f32> {
    if line.contains("Reading data file") || line.contains("Created") {
        Some(10.0)
    } else if line.contains("Setting up") {
        Some(25.0)
    } else if line.contains("Per MPI rank memory allocation") {
        Some(40.0)
    } else if line.contains("Loop time") {
        Some(95.0)
    } else if line.contains("Total wall time") {
        Some(99.0)
    } else {
        None
    }
}

/// Stage any external file (read_data path, EAM potential) into the
/// workdir so the deck's relative path resolves. Non-existent paths
/// are left for LAMMPS to complain about at run time.
fn stage_external_files(
    case_dir: &Path,
    workdir: &Path,
    input: &LammpsInput,
) -> Result<(), AdapterError> {
    use case_input::Initialization;

    let mut to_stage: Vec<std::path::PathBuf> = Vec::new();
    if let Initialization::ReadData { path } = &input.initialization {
        to_stage.push(path.clone());
    }
    if let Potential::Eam { path, .. } = &input.potential {
        to_stage.push(path.clone());
    }
    for rel in to_stage {
        if rel.is_absolute() {
            continue;
        }
        // Round-9 hardening: `read_data` path + EAM potential path
        // are user-supplied data that get *copied* into the workdir;
        // wrap with `confined_join` so a hostile case can't ask
        // LAMMPS to stage `../../etc/passwd`. Propagate the rejection
        // — silently continuing would mask the threat the helper
        // exists to catch.
        let src = valenx_core::adapter_helpers::confined_join(case_dir, &rel)?;
        if !src.is_file() {
            continue;
        }
        let Some(name) = rel.file_name() else {
            continue;
        };
        let dst = workdir.join(name);
        if src != dst {
            fs::copy(&src, &dst)?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn info_is_md_domain() {
        let info = LammpsAdapter::new().info();
        assert_eq!(info.id, "lammps");
        assert_eq!(info.physics, &[Physics::MolecularDynamics]);
    }

    #[test]
    fn collect_uses_live_provenance_with_real_case_hash() {
        let workdir = tempdir("lammps-collect-prov");
        let case_path = workdir.join(DEFAULT_INPUT_FILENAME);
        let case_bytes = b"units lj\nrun 0\n";
        std::fs::write(&case_path, case_bytes).expect("write in.lammps");

        let job = PreparedJob {
            workdir: workdir.clone(),
            native_command: Vec::new(),
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: false,
        };
        let results = LammpsAdapter::new().collect(&job).expect("collect");
        let prov = &results.provenance;

        assert_eq!(prov.adapter, INFO_ID);
        assert!(!prov.adapter_version.is_empty());
        assert_eq!(prov.tool, "LAMMPS");
        assert!(!prov.run_id.is_empty(), "run_id empty — stub still wired?");
        assert_eq!(
            prov.case_hash,
            valenx_core::adapter_helpers::sha256_hex_file(&case_path)
        );

        cleanup_lp(&workdir);
    }

    fn cleanup_lp(d: &std::path::Path) {
        let _ = std::fs::remove_dir_all(d);
    }

    #[test]
    fn progress_hints_monotonic() {
        let pts = [
            lammps_progress_hint("Reading data file ..."),
            lammps_progress_hint("Setting up Verlet run ..."),
            lammps_progress_hint("Per MPI rank memory allocation (min/avg/max) ="),
            lammps_progress_hint("Loop time of 0.1234 on 1 procs"),
            lammps_progress_hint("Total wall time: 0:00:01"),
        ];
        let mut last = 0.0f32;
        for (i, p) in pts.iter().enumerate() {
            let v = p.expect("banner");
            assert!(v >= last, "step {i}: {last} -> {v}");
            last = v;
        }
    }

    /// Round-9 RED→GREEN: `stage_external_files` joined `read_data` +
    /// EAM potential paths with bare `case_dir.join`. Wrap with
    /// `confined_join` so a hostile case can't ask LAMMPS to stage
    /// `../../etc/passwd` into the workdir.
    #[test]
    fn stage_external_files_rejects_read_data_path_traversing_outside_case_dir() {
        use case_input::{BoundaryCondition, Ensemble, Initialization, Potential, Units};
        let case_dir = tempdir("lammps-readdata-trav-case");
        let workdir = tempdir("lammps-readdata-trav-work");
        let input = LammpsInput {
            units: Units::Real,
            boundary: [BoundaryCondition::P; 3],
            atom_style: "atomic".into(),
            initialization: Initialization::ReadData {
                path: std::path::PathBuf::from("../../etc/passwd"),
            },
            potential: Potential::LjCut {
                epsilon: 1.0,
                sigma: 1.0,
                cutoff: 2.5,
            },
            ensemble: Ensemble::Nve,
            run_steps: 1,
            timestep: 0.001,
            initial_temperature: None,
            thermo_every: 100,
            dump_every: 0,
        };
        let err = stage_external_files(&case_dir, &workdir, &input).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("..") || msg.contains("stay within") || msg.contains("escape"),
            "expected confined_join rejection, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&case_dir);
        let _ = std::fs::remove_dir_all(&workdir);
    }
}
