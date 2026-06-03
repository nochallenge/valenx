//! # valenx-adapter-gromacs
//!
//! Adapter for GROMACS — biomolecular MD engine.
//!
//! **Phase 5 — live for `gmx mdrun` on a pre-built `.tpr`.** The
//! `gmx grompp` preprocessing step is left to the user — its
//! topology + force-field handling is outside the scope a single
//! adapter can sanely manage, and most workflows already keep the
//! `.tpr` as the "ready-to-run" artifact. `prepare()` parses
//! `[md.gromacs]` from case.toml, stages the .tpr into the workdir,
//! and builds `gmx mdrun -s <tpr> -deffnm <name> [-nt N]`. `run()`
//! spawns it. `collect()` walks for `.trr` / `.xtc` / `.edr` /
//! `.gro` / `.log` outputs.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::fs;
use std::path::Path;
use std::time::Duration;

use semver::Version;

use valenx_core::{
    adapter_helpers::{confined_join, find_on_path, first_workdir_match},
    error::RunPhase,
    subprocess, Adapter, AdapterError, AdapterInfo, Capabilities, Capability, Case, LicenseMode,
    Physics, PreparedJob, ProbeReport, RunContext, RunReport, VersionRange,
};
use valenx_fields::{
    artifact::{Artifact, ArtifactKind},
    Results,
};

use crate::case_input::GromacsInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(GromacsAdapter::new())
}

pub struct GromacsAdapter;

impl GromacsAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for GromacsAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "gromacs";
const BINARIES: &[&str] = &["gmx", "gmx_mpi"];

impl Adapter for GromacsAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "GROMACS",
            version_range: VersionRange {
                min_inclusive: Version::new(2023, 0, 0),
                max_exclusive: Version::new(2030, 0, 0),
            },
            physics: &[Physics::MolecularDynamics],
            license_mode: LicenseMode::Subprocess,
            tool_license: "LGPL-2.1-or-later",
            docs_url: "https://manual.gromacs.org/",
            homepage_url: "https://www.gromacs.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                let found_version = valenx_core::adapter_helpers::detect_tool_version_semver(
                    &binary_path,
                    &["--version", "-version"],
                );
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    warnings: vec!["gromacs adapter runs `gmx mdrun` on a pre-built .tpr — \
                         the user is responsible for `gmx grompp`"
                        .into()],
                    required_env: Vec::new(),
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "GROMACS 2023+ required; install from gromacs.org".into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = GromacsInput::from_case_dir(&case.path)?;
        fs::create_dir_all(workdir)?;

        // `confined_join` rejects absolute paths and `..` traversal —
        // a shared case bundle should not be able to point `tpr` at
        // an arbitrary host file.
        let source = confined_join(&case.path, &input.tpr)?;
        if !source.is_file() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "tpr not found at {} (resolve relative to case dir)",
                source.display()
            )));
        }
        let tpr_name = source
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("system.tpr")
            .to_string();
        fs::copy(&source, workdir.join(&tpr_name))
            .map_err(|e| AdapterError::Other(anyhow::anyhow!("stage {}: {e}", source.display())))?;

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "no `gmx` / `gmx_mpi` on PATH".into(),
        })?;

        // Build: gmx mdrun -s <tpr> -deffnm <name> [-nt N]
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from("mdrun"),
            OsString::from("-s"),
            OsString::from(&tpr_name),
            OsString::from("-deffnm"),
            OsString::from(&input.deffnm),
        ];
        if let Some(nt) = input.nt {
            native_command.push(OsString::from("-nt"));
            native_command.push(OsString::from(nt.to_string()));
        }

        // MD trajectories vary wildly; default 4 hours as the GUI
        // cancellation ceiling. Production runs use SLURM time_limit.
        let estimated_runtime = Some(Duration::from_secs(4 * 60 * 60));

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            estimated_runtime,
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting GROMACS mdrun", |line| {
            let mut hint = subprocess::Hint::default();
            if let Some(pct) = gromacs_progress_hint(line) {
                hint.progress = Some((pct, line.to_string()));
            }
            // GROMACS prints "Fatal error:" on hard failures and
            // "Note:" / "WARNING:" for advisories. Surface the
            // serious ones.
            if line.contains("Fatal error") || line.contains("Error in user input") {
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
        // case_path: prefer the staged .tpr for the canonical hash
        // (it captures topology + parameters in one binary blob).
        // Fall back to .mdp if no .tpr is present.
        let case_path = first_workdir_match(&job.workdir, &["tpr"])
            .or_else(|| first_workdir_match(&job.workdir, &["mdp"]))
            .unwrap_or_else(|| job.workdir.join("(no-tpr-found)"));
        let mesh_path = first_workdir_match(&job.workdir, &["gro", "pdb", "top"]);
        let prov = valenx_core::adapter_helpers::live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "GROMACS",
            "unknown",
            &case_path,
            mesh_path.as_deref(),
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);

        // GROMACS produces several output flavours under the
        // `-deffnm` prefix:
        //   .trr / .xtc — trajectory (binary)
        //   .edr        — energy (binary)
        //   .gro        — final coordinates
        //   .log        — text log
        //   .cpt        — checkpoint
        let classifications: &[(&str, ArtifactKind, &str)] = &[
            ("trr", ArtifactKind::Native, "GROMACS .trr trajectory"),
            (
                "xtc",
                ArtifactKind::Native,
                "GROMACS .xtc trajectory (compressed)",
            ),
            ("edr", ArtifactKind::Native, "GROMACS .edr energy file"),
            ("gro", ArtifactKind::Other, "GROMACS .gro coordinates"),
            ("pdb", ArtifactKind::Other, "PDB coordinates"),
            ("top", ArtifactKind::Other, "GROMACS topology"),
            ("tpr", ArtifactKind::Other, "GROMACS tpr (binary deck)"),
            ("mdp", ArtifactKind::Other, "GROMACS mdp parameters"),
            ("cpt", ArtifactKind::Other, "GROMACS checkpoint"),
            ("log", ArtifactKind::Log, "GROMACS log"),
        ];
        if let Ok(entries) = fs::read_dir(&job.workdir) {
            for entry in entries.flatten() {
                let path = entry.path();
                let ext = path
                    .extension()
                    .and_then(|s| s.to_str())
                    .map(|s| s.to_ascii_lowercase());
                let Some(ext) = ext else { continue };
                if let Some(&(_, kind, label)) =
                    classifications.iter().find(|(e, ..)| *e == ext.as_str())
                {
                    results.artifacts.push(Artifact {
                        path,
                        kind,
                        checksum: None,
                        label: label.to_string(),
                    });
                }
            }
        }
        results.artifacts.sort_by(|a, b| a.path.cmp(&b.path));
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: vec![Capability::MdClassical],
            ribbon_contributions: vec!["md.gromacs.grompp", "md.gromacs.mdrun"],
        }
    }
}

/// Coarse progress hints for `gmx mdrun` stdout banners.
fn gromacs_progress_hint(line: &str) -> Option<f32> {
    if line.contains("Reading file") || line.contains("Reading checkpoint") {
        Some(10.0)
    } else if line.contains("Will use") && line.contains("particle-particle") {
        Some(25.0)
    } else if line.contains("starting mdrun") || line.contains("Starting") {
        Some(40.0)
    } else if line.contains("step ") {
        // Mid-run ticks. Stay flat so the bar doesn't oscillate.
        Some(50.0)
    } else if line.contains("Writing final coordinates") {
        Some(95.0)
    } else if line.contains("Performance") || line.contains("Finished mdrun") {
        Some(98.0)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_core::adapter_helpers::sha256_hex_file;
    use valenx_test_utils::tempdir;

    #[test]
    fn info_is_well_formed() {
        let info = GromacsAdapter::new().info();
        assert_eq!(info.id, "gromacs");
    }

    #[test]
    fn collect_uses_live_provenance_with_real_case_hash() {
        let workdir = tempdir("gromacs-collect");
        let mdp_path = workdir.join("md.mdp");
        let mdp_bytes = b"integrator = md\nnsteps = 1000\n";
        std::fs::write(&mdp_path, mdp_bytes).expect("write .mdp");

        let job = PreparedJob {
            workdir: workdir.clone(),
            native_command: Vec::new(),
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: false,
        };
        let results = GromacsAdapter::new().collect(&job).expect("collect");
        let prov = &results.provenance;

        assert_eq!(prov.adapter, INFO_ID);
        assert!(!prov.adapter_version.is_empty());
        assert_eq!(prov.tool, "GROMACS");
        assert!(!prov.run_id.is_empty(), "run_id empty — stub still wired?");
        assert_eq!(prov.case_hash, sha256_hex_file(&mdp_path));

        cleanup(&workdir);
    }

    fn cleanup(d: &std::path::Path) {
        let _ = std::fs::remove_dir_all(d);
    }

    #[test]
    fn progress_hints_are_monotonic() {
        let pts = [
            gromacs_progress_hint("Reading file system.tpr, VERSION 2024"),
            gromacs_progress_hint("Will use 8 particle-particle and 2 PME ranks"),
            gromacs_progress_hint("starting mdrun 'Lysozyme'"),
            gromacs_progress_hint("step 1000  Time   2.000"),
            gromacs_progress_hint("Writing final coordinates."),
            gromacs_progress_hint("Performance:    100.0   ns/day"),
        ];
        let mut last = 0.0_f32;
        for (i, p) in pts.iter().enumerate() {
            let v = p.unwrap_or_else(|| panic!("step {i} returned None"));
            assert!(v >= last, "step {i}: {last} -> {v}");
            last = v;
        }
    }

    #[test]
    fn collect_classifies_trajectory_outputs() {
        let workdir = tempdir("gromacs-gmx-collect");
        for (name, content) in [
            ("md.trr", &b"binary"[..]),
            ("md.xtc", &b"binary"[..]),
            ("md.edr", &b"binary"[..]),
            ("md.gro", &b"# coords"[..]),
            ("md.log", &b"text log"[..]),
            ("system.tpr", &b"binary"[..]),
            ("ignored.txt", &b"unrelated"[..]),
        ] {
            std::fs::write(workdir.join(name), content).unwrap();
        }
        let job = PreparedJob {
            workdir: workdir.clone(),
            native_command: Vec::new(),
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: false,
        };
        let results = GromacsAdapter::new().collect(&job).expect("collect");
        // 6 known artifacts; ignored.txt skipped.
        assert_eq!(results.artifacts.len(), 6);
        let labels: Vec<&str> = results.artifacts.iter().map(|a| a.label.as_str()).collect();
        assert!(labels.iter().any(|l| l.contains(".trr trajectory")));
        assert!(labels.iter().any(|l| l.contains(".xtc trajectory")));
        assert!(labels.iter().any(|l| l.contains(".edr energy")));
        assert!(labels.iter().any(|l| l.contains("tpr")));
        cleanup(&workdir);
    }

    #[test]
    fn prepare_stages_tpr_and_builds_command() {
        let case_dir = tempdir("gromacs-gmx-prepare");
        std::fs::write(
            case_dir.join("case.toml"),
            "[md.gromacs]\ntpr = \"system.tpr\"\nnt = 4\ndeffnm = \"prod\"\n",
        )
        .unwrap();
        std::fs::write(case_dir.join("system.tpr"), b"binary").unwrap();
        let workdir = tempdir("gromacs-gmx-prepare-wd");
        let case = Case {
            id: "gmx-test".into(),
            path: case_dir.clone(),
        };
        let r = GromacsAdapter::new().prepare(&case, &workdir);
        if find_on_path(BINARIES).is_none() {
            assert!(matches!(r, Err(AdapterError::ToolNotInstalled { .. })));
            cleanup(&case_dir);
            cleanup(&workdir);
            return;
        }
        let job = r.expect("prepare");
        assert!(workdir.join("system.tpr").is_file());
        let cmd: Vec<String> = job
            .native_command
            .iter()
            .map(|s| s.to_string_lossy().into_owned())
            .collect();
        assert!(cmd.iter().any(|s| s == "mdrun"));
        assert!(cmd.iter().any(|s| s == "-s"));
        assert!(cmd.iter().any(|s| s == "system.tpr"));
        assert!(cmd.iter().any(|s| s == "-deffnm"));
        assert!(cmd.iter().any(|s| s == "prod"));
        assert!(cmd.iter().any(|s| s == "-nt"));
        assert!(cmd.iter().any(|s| s == "4"));
        cleanup(&case_dir);
        cleanup(&workdir);
    }

    #[test]
    fn prepare_missing_tpr_is_actionable() {
        let case_dir = tempdir("gromacs-gmx-no-tpr");
        std::fs::write(
            case_dir.join("case.toml"),
            "[md.gromacs]\ntpr = \"missing.tpr\"\n",
        )
        .unwrap();
        let workdir = tempdir("gromacs-gmx-no-tpr-wd");
        let case = Case {
            id: "gmx-test".into(),
            path: case_dir.clone(),
        };
        assert!(GromacsAdapter::new().prepare(&case, &workdir).is_err());
        cleanup(&case_dir);
        cleanup(&workdir);
    }
}
