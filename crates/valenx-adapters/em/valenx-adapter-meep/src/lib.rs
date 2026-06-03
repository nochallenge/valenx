//! # valenx-adapter-meep
//!
//! Adapter for Meep (MIT) — finite-difference time-domain photonics.
//! Python-driven via generated scripts; complements openEMS for the
//! photonics/optical end of the EM spectrum.
//!
//! **Phase 6 — live for batch script execution.** `prepare()`
//! parses `[em.meep]` from case.toml, stages the simulation
//! script (Python or legacy Scheme `.ctl`) into the workdir, and
//! builds either `python <script>` or `meep <ctl>`. With `np > 1`
//! the invocation wraps in `mpirun -np N`. `run()` spawns it.
//! `collect()` walks the workdir for `.h5` field dumps and the
//! input script for the audit trail.

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

use crate::case_input::MeepInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(MeepAdapter::new())
}

pub struct MeepAdapter;

impl MeepAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for MeepAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "meep";
/// Python-mode script interpreters (preferred — modern Meep is
/// Python-bound).
const PYTHON_BINARIES: &[&str] = &["python3", "python"];
/// Legacy Scheme `meep` interpreter for `.ctl` scripts.
const SCHEME_BINARIES: &[&str] = &["meep"];
/// Used by probe — succeeds if either interpreter is present.
const BINARIES: &[&str] = &["python3", "python", "meep"];

impl Adapter for MeepAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "Meep",
            version_range: VersionRange {
                min_inclusive: Version::new(1, 28, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Em],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-2.0-or-later",
            docs_url: "https://meep.readthedocs.io/",
            homepage_url: "https://meep.readthedocs.io/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                let found_version = valenx_core::adapter_helpers::detect_tool_version_semver(
                    &binary_path,
                    &["--version", "-V"],
                );
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    warnings: vec!["meep adapter requires the meep Python module installed in \
                         the same Python environment for Python-mode scripts"
                        .into()],
                    required_env: Vec::new(),
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "Meep required; install via `pip install meep` or conda-forge".into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = MeepInput::from_case_dir(&case.path)?;
        fs::create_dir_all(workdir)?;

        let source = confined_join(&case.path, &input.script)?;
        if !source.is_file() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "script not found at {} (resolve relative to case dir)",
                source.display()
            )));
        }
        let staged_name = source
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("sim.py")
            .to_string();
        let staged = workdir.join(&staged_name);
        fs::copy(&source, &staged)
            .map_err(|e| AdapterError::Other(anyhow::anyhow!("stage {}: {e}", source.display())))?;

        // Pick the right interpreter: python for `.py` (or
        // python=true), meep for legacy `.ctl`.
        let binaries: &[&str] = if input.python {
            PYTHON_BINARIES
        } else {
            SCHEME_BINARIES
        };
        let binary_path = find_on_path(binaries).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: if input.python {
                "no `python3` / `python` on PATH (and the meep module \
                     must be installed in that environment)"
                    .into()
            } else {
                "no `meep` (Scheme/CTL) interpreter on PATH; switch to \
                     Python by setting `python = true`"
                    .into()
            },
        })?;

        // Build the argv. mpirun wraps when np > 1.
        let mut native_command: Vec<OsString> = Vec::new();
        if let Some(np) = input.np {
            if np > 1 {
                let mpirun = find_on_path(&["mpirun", "mpiexec"]).ok_or_else(|| {
                    AdapterError::ToolNotInstalled {
                        name: INFO_ID,
                        hint: "np > 1 needs mpirun / mpiexec on PATH".into(),
                    }
                })?;
                native_command.push(mpirun.into_os_string());
                native_command.push(OsString::from("-np"));
                native_command.push(OsString::from(np.to_string()));
            }
        }
        native_command.push(binary_path.into_os_string());
        native_command.push(OsString::from(&staged_name));

        // FDTD photonics sims at moderate resolution finish in
        // minutes; large 3D sims are hours. 1-hour default ceiling.
        let estimated_runtime = Some(Duration::from_secs(60 * 60));

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            estimated_runtime,
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting Meep", |line| {
            let mut hint = subprocess::Hint::default();
            if let Some(pct) = meep_progress_hint(line) {
                hint.progress = Some((pct, line.to_string()));
            }
            // Meep's diagnostics have several flavours. Surface
            // both warning lines and Python tracebacks as warnings
            // so the GUI residual panel picks them up.
            if line.starts_with("meep: error")
                || line.contains("Traceback (most recent call last)")
                || line.contains("Error:")
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
        // Meep is Python-driven (or legacy Scheme) — case_path is the
        // first .py / .ctl in the workdir (the simulation script).
        let case_path = first_workdir_match(&job.workdir, &["py", "ctl"])
            .unwrap_or_else(|| job.workdir.join("(no-script-found)"));
        let prov = valenx_core::adapter_helpers::live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "Meep",
            "unknown",
            &case_path,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);

        // Walk the workdir for outputs. Meep's most common dump is
        // HDF5 (.h5 from `meep.simulation.output_*` calls); some
        // scripts emit text logs alongside.
        let classifications: &[(&str, ArtifactKind, &str)] = &[
            ("h5", ArtifactKind::Native, "Meep HDF5 field dump"),
            ("hdf5", ArtifactKind::Native, "Meep HDF5 field dump"),
            ("py", ArtifactKind::Other, "Meep Python script"),
            ("ctl", ArtifactKind::Other, "Meep Scheme/CTL script"),
            ("log", ArtifactKind::Log, "Meep log"),
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
            capabilities: vec![Capability::EmFdtdTimeDomain],
            ribbon_contributions: vec!["em.meep.fdtd"],
        }
    }
}

/// Coarse progress hints for Meep stdout banners. Meep's
/// per-step output uses the form `time: 1.23 / total: 10.0`
/// for time-stepping; we capture the start banner and the first
/// few step banners forward of that.
fn meep_progress_hint(line: &str) -> Option<f32> {
    if line.contains("Initializing structure") {
        Some(15.0)
    } else if line.contains("time for set_epsilon") {
        Some(30.0)
    } else if line.contains("on time step") {
        // Mid-run ticks. Stay flat at 50 % so the bar doesn't
        // flicker on every step banner.
        Some(50.0)
    } else if line.contains("run 0 finished") || line.contains("Field decay") {
        Some(95.0)
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
    fn info_is_em() {
        let info = MeepAdapter::new().info();
        assert_eq!(info.physics, &[Physics::Em]);
    }

    #[test]
    fn collect_uses_live_provenance_with_real_case_hash() {
        let workdir = tempdir("meep-collect");
        let script_path = workdir.join("sim.py");
        let script_bytes = b"import meep as mp\n# trivial sim\n";
        std::fs::write(&script_path, script_bytes).expect("write sim.py");

        let job = PreparedJob {
            workdir: workdir.clone(),
            native_command: Vec::new(),
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: false,
        };
        let results = MeepAdapter::new().collect(&job).expect("collect");
        let prov = &results.provenance;

        assert_eq!(prov.adapter, INFO_ID);
        assert!(!prov.adapter_version.is_empty());
        assert_eq!(prov.tool, "Meep");
        assert!(!prov.run_id.is_empty(), "run_id empty — stub still wired?");
        assert_eq!(prov.case_hash, sha256_hex_file(&script_path));

        cleanup(&workdir);
    }

    fn cleanup(d: &std::path::Path) {
        let _ = std::fs::remove_dir_all(d);
    }

    #[test]
    fn progress_hints_are_monotonic() {
        let pts = [
            meep_progress_hint("Initializing structure..."),
            meep_progress_hint("time for set_epsilon = 0.00s"),
            meep_progress_hint("on time step 100 (time=10), 0.0001 s/step"),
            meep_progress_hint("run 0 finished at t=20.0"),
        ];
        let mut last = 0.0_f32;
        for (i, p) in pts.iter().enumerate() {
            let v = p.unwrap_or_else(|| panic!("step {i} returned None"));
            assert!(v >= last, "step {i}: {last} -> {v}");
            last = v;
        }
    }

    #[test]
    fn collect_classifies_h5_and_script() {
        let workdir = tempdir("meep-collect-h5");
        std::fs::write(workdir.join("ez-000050.00.h5"), b"binary").unwrap();
        std::fs::write(workdir.join("sim.py"), b"import meep as mp\n").unwrap();
        std::fs::write(workdir.join("ignored.txt"), b"unrelated").unwrap();
        let job = PreparedJob {
            workdir: workdir.clone(),
            native_command: Vec::new(),
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: false,
        };
        let results = MeepAdapter::new().collect(&job).expect("collect");
        assert_eq!(results.artifacts.len(), 2);
        let labels: Vec<&str> = results.artifacts.iter().map(|a| a.label.as_str()).collect();
        assert!(labels.iter().any(|l| l.contains("HDF5")));
        assert!(labels.iter().any(|l| l.contains("Python script")));
        cleanup(&workdir);
    }

    #[test]
    fn prepare_python_path_stages_and_builds_command() {
        let case_dir = tempdir("meep-prepare");
        std::fs::write(
            case_dir.join("case.toml"),
            "[em.meep]\nscript = \"ring.py\"\n",
        )
        .unwrap();
        std::fs::write(case_dir.join("ring.py"), b"import meep\n").unwrap();
        let workdir = tempdir("meep-prepare-wd");
        let case = Case {
            id: "meep-test".into(),
            path: case_dir.clone(),
        };
        let r = MeepAdapter::new().prepare(&case, &workdir);
        if find_on_path(PYTHON_BINARIES).is_none() {
            assert!(matches!(r, Err(AdapterError::ToolNotInstalled { .. })));
            cleanup(&case_dir);
            cleanup(&workdir);
            return;
        }
        let job = r.expect("prepare");
        assert!(workdir.join("ring.py").is_file());
        let cmd: Vec<String> = job
            .native_command
            .iter()
            .map(|s| s.to_string_lossy().into_owned())
            .collect();
        assert!(cmd.iter().any(|s| s == "ring.py"));
        cleanup(&case_dir);
        cleanup(&workdir);
    }

    #[test]
    fn prepare_missing_script_is_actionable() {
        let case_dir = tempdir("meep-no-script");
        std::fs::write(
            case_dir.join("case.toml"),
            "[em.meep]\nscript = \"missing.py\"\n",
        )
        .unwrap();
        let workdir = tempdir("meep-no-script-wd");
        let case = Case {
            id: "meep-test".into(),
            path: case_dir.clone(),
        };
        assert!(MeepAdapter::new().prepare(&case, &workdir).is_err());
        cleanup(&case_dir);
        cleanup(&workdir);
    }
}
