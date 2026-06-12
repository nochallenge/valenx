//! # valenx-adapter-su2
//!
//! Adapter for SU2 — compressible CFD with strong adjoint and
//! shape-optimisation support. Complements OpenFOAM on the
//! compressible / aero-optimisation side.
//!
//! **Phase 1.5 — live for batch SU2_CFD runs.** `prepare()` parses
//! `[cfd.su2]` from case.toml, stages the cfg + mesh into the
//! workdir, and builds a `SU2_CFD <cfg>` invocation. `run()` spawns
//! it via the shared subprocess runner. `collect()` walks the
//! workdir for `.vtu` / `.vtk` outputs and parses them into the
//! canonical Field catalog.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Duration;

use semver::Version;

use valenx_core::{
    adapter_helpers::{find_on_path, first_workdir_match},
    error::RunPhase,
    io_caps::{read_capped_to_bytes, MAX_VTK_FILE_BYTES},
    subprocess, Adapter, AdapterError, AdapterInfo, Capabilities, Capability, Case, LicenseMode,
    Physics, PreparedJob, ProbeReport, RunContext, RunReport, VersionRange,
};
use valenx_fields::{
    artifact::{Artifact, ArtifactKind},
    Results,
};

use crate::case_input::Su2Input;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(Su2Adapter::new())
}

pub struct Su2Adapter;

impl Su2Adapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for Su2Adapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "su2";
const BINARIES: &[&str] = &["SU2_CFD", "SU2_SOL"];

impl Adapter for Su2Adapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "SU2",
            version_range: VersionRange {
                min_inclusive: Version::new(8, 0, 0),
                max_exclusive: Version::new(9, 0, 0),
            },
            physics: &[Physics::Cfd],
            license_mode: LicenseMode::Subprocess,
            tool_license: "LGPL-2.1-only",
            docs_url: "https://su2code.github.io/docs_v7/home/",
            homepage_url: "https://su2code.github.io/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                let found_version = valenx_core::adapter_helpers::detect_tool_version_semver(
                    &binary_path,
                    &["--version", "-v"],
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
                hint: "SU2 8.0+ required; install from su2code.github.io".into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = Su2Input::from_case_dir(&case.path)?;
        fs::create_dir_all(workdir)?;

        // Stage the SU2 .cfg into the workdir.
        let cfg_source = valenx_core::adapter_helpers::confined_join(&case.path, &input.config)?;
        if !cfg_source.is_file() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "config not found at {} (resolve relative to case dir)",
                cfg_source.display()
            )));
        }
        let cfg_name = cfg_source
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("case.cfg")
            .to_string();
        fs::copy(&cfg_source, workdir.join(&cfg_name)).map_err(|e| {
            AdapterError::Other(anyhow::anyhow!("stage {}: {e}", cfg_source.display()))
        })?;

        // Stage the mesh if specified — keeps the cfg's
        // `MESH_FILENAME=` reference resolvable.
        // Round-9 hardening: `mesh` is user-supplied data and gets
        // copied into the workdir; wrap with `confined_join` so a
        // hostile case can't ask SU2 to read `../../etc/passwd`.
        if let Some(mesh) = &input.mesh {
            let mesh_source = valenx_core::adapter_helpers::confined_join(&case.path, mesh)?;
            if !mesh_source.is_file() {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "mesh not found at {} (resolve relative to case dir)",
                    mesh_source.display()
                )));
            }
            let mesh_name = mesh_source
                .file_name()
                .and_then(|s| s.to_str())
                .unwrap_or("mesh.su2")
                .to_string();
            fs::copy(&mesh_source, workdir.join(&mesh_name)).map_err(|e| {
                AdapterError::Other(anyhow::anyhow!("stage {}: {e}", mesh_source.display()))
            })?;
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "no SU2 binary on PATH".into(),
        })?;

        let native_command: Vec<OsString> =
            vec![binary_path.into_os_string(), OsString::from(&cfg_name)];

        // Threading: SU2 reads OMP_NUM_THREADS at startup. Setting
        // it via the prepared job's env keeps the submitter shell
        // unaffected and the value reproducible per-run.
        let mut environment: Vec<(OsString, OsString)> = Vec::new();
        if let Some(n) = input.n_threads {
            environment.push((
                OsString::from("OMP_NUM_THREADS"),
                OsString::from(n.to_string()),
            ));
        }

        // Aero sims at production resolution can take hours; default
        // to 2 hours as the GUI cancellation ceiling. Users can pick
        // a smaller value for unit-test cases via the SLURM
        // time_limit knob.
        let estimated_runtime = Some(Duration::from_secs(2 * 60 * 60));

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment,
            estimated_runtime,
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting SU2_CFD", |line| {
            let mut hint = subprocess::Hint::default();
            if let Some(pct) = su2_progress_hint(line) {
                hint.progress = Some((pct, line.to_string()));
            }
            // SU2 prints `Error in ...` for fatal issues.
            if line.contains("Error in ") || line.contains("ERROR:") {
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
        // Real provenance: hash whatever inputs are recoverable from
        // the workdir. SU2 typically writes its config as
        // `<case>.cfg` next to the run; we look for the first one.
        // Mesh hash comes from the first .su2 / .cgns file we find.
        // Both fall back to empty Sha256Hex sentinels when missing,
        // matching stub_provenance's policy.
        let case_path = first_workdir_match(&job.workdir, &["cfg"])
            .unwrap_or_else(|| job.workdir.join("(no-config-found)"));
        let mesh_path = first_workdir_match(&job.workdir, &["su2", "cgns", "msh"]);
        let prov = valenx_core::adapter_helpers::live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "SU2",
            "unknown", // tool_version is per-run; could fill from probe later
            &case_path,
            mesh_path.as_deref(),
            None, // tools.lock lives at the project level, not the workdir
            0.0,  // wall_time fills in once SU2 run() is implemented
        );
        let mut results = Results::empty(INFO_ID, prov);
        // SU2 writes VTK output (default `.vtu` from WRT_FORMAT, or
        // `.vtk` from `WRT_FORMAT = LEGACY_VTK`). Walk the workdir
        // for either extension; vtk_dispatch routes each file to
        // the right reader based on its magic prefix. Also
        // registers the path as an Artifact so users can open it
        // from the Results pane.
        for path in walk_vtk_files(&job.workdir) {
            results.artifacts.push(Artifact {
                path: path.clone(),
                kind: ArtifactKind::VizData,
                checksum: None,
                label: format!(
                    "SU2 VTK output ({})",
                    path.file_name()
                        .map(|s| s.to_string_lossy().to_string())
                        .unwrap_or_default()
                ),
            });
            // Round-20 L1: cap the per-VTK read so a corrupted (or
            // adversarial) workdir with a multi-GB `.vtu` / `.vtk`
            // file can't OOM the renderer before `vtk_dispatch`
            // validates the magic bytes. Production CFD VTU
            // snapshots top out around 1 GiB for an HPC mesh; the
            // 4 GiB cap is generous while refusing the
            // `cat /dev/zero > big.vtk` DoS.
            let Ok(bytes) = read_capped_to_bytes(&path, MAX_VTK_FILE_BYTES) else {
                continue;
            };
            let stem = path.file_stem().and_then(|s| s.to_str()).unwrap_or("su2");
            if let Ok((_mesh, fields)) = valenx_fields::vtk_dispatch::load_canonical(&bytes, stem) {
                for f in fields {
                    results.fields.insert(f);
                }
            }
            // Per-file parse failures skip silently — the Artifact
            // is still listed so the user can open the file by hand.
        }
        results.artifacts.sort_by(|a, b| a.path.cmp(&b.path));
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: vec![
                Capability::CfdCompressible,
                Capability::CfdSteady,
                Capability::CfdTransient,
                Capability::CfdTurbulenceRans,
                Capability::CfdAdjointOptimization,
            ],
            ribbon_contributions: vec![
                "cfd.su2.compressible",
                "cfd.su2.adjoint",
                "cfd.su2.optimisation",
            ],
        }
    }
}

/// Recursively walk a workdir for `.vtu` (XML) or `.vtk` (legacy
/// binary) files. Used by [`Su2Adapter::collect`] to discover SU2's
/// VTK output without depending on a specific filename convention.
/// Coarse progress hints for SU2_CFD stdout banners. SU2 prints
/// recognisable section headers on startup + an iteration counter
/// during the solve. We map the start-up banners forward; mid-run
/// iteration ticks stay at 50 % so the bar doesn't flicker.
fn su2_progress_hint(line: &str) -> Option<f32> {
    if line.contains("---------------- Reading the mesh") {
        Some(10.0)
    } else if line.contains("---------------- Setting up the configuration") {
        Some(20.0)
    } else if line.contains("---------------- Numerical Grid") {
        Some(30.0)
    } else if line.contains("---------------- Begin Solver") {
        Some(40.0)
    } else if line.starts_with("ITERATION") || line.contains("MG Level") {
        Some(50.0)
    } else if line.contains("---------------- Solution") {
        Some(90.0)
    } else if line.contains("Exit Success") || line.contains("---- End SU2") {
        Some(98.0)
    } else {
        None
    }
}

fn walk_vtk_files(root: &Path) -> Vec<PathBuf> {
    let mut found: Vec<PathBuf> = Vec::new();
    let mut stack: Vec<PathBuf> = vec![root.to_path_buf()];
    while let Some(dir) = stack.pop() {
        let Ok(entries) = std::fs::read_dir(&dir) else {
            continue;
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if let Ok(ft) = entry.file_type() {
                if ft.is_dir() {
                    stack.push(path);
                    continue;
                }
            }
            let ext_match = path
                .extension()
                .and_then(|s| s.to_str())
                .map(|s| {
                    let s = s.to_ascii_lowercase();
                    s == "vtu" || s == "vtk"
                })
                .unwrap_or(false);
            if ext_match {
                found.push(path);
            }
        }
    }
    found.sort();
    found
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn info_is_cfd() {
        let info = Su2Adapter::new().info();
        assert_eq!(info.physics, &[Physics::Cfd]);
    }

    #[test]
    fn progress_hints_are_monotonic() {
        let pts = [
            su2_progress_hint("---------------- Reading the mesh ----------------"),
            su2_progress_hint("---------------- Setting up the configuration ----------------"),
            su2_progress_hint("---------------- Numerical Grid Adaption ----------------"),
            su2_progress_hint("---------------- Begin Solver ----------------"),
            su2_progress_hint("ITERATION 100"),
            su2_progress_hint("---------------- Solution Output ----------------"),
            su2_progress_hint("Exit Success (SU2_CFD)"),
        ];
        let mut last = 0.0_f32;
        for (i, p) in pts.iter().enumerate() {
            let v = p.unwrap_or_else(|| panic!("step {i} returned None"));
            assert!(v >= last, "step {i}: {last} -> {v}");
            last = v;
        }
    }

    #[test]
    fn prepare_stages_cfg_and_mesh_with_threading_env() {
        let case_dir = std::env::temp_dir().join(format!(
            "valenx-su2-prep-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&case_dir).unwrap();
        std::fs::write(
            case_dir.join("case.toml"),
            r#"
[cfd.su2]
config = "wing.cfg"
mesh = "wing.su2"
n_threads = 4
"#,
        )
        .unwrap();
        std::fs::write(case_dir.join("wing.cfg"), b"% placeholder cfg").unwrap();
        std::fs::write(case_dir.join("wing.su2"), b"% placeholder mesh").unwrap();

        let workdir = case_dir.join("workdir");
        let case = Case {
            id: "su2-test".into(),
            path: case_dir.clone(),
        };
        let r = Su2Adapter::new().prepare(&case, &workdir);
        if find_on_path(BINARIES).is_none() {
            // No SU2 binary on PATH — should error with ToolNotInstalled
            // since prepare() needs to record the binary location.
            assert!(matches!(r, Err(AdapterError::ToolNotInstalled { .. })));
            let _ = std::fs::remove_dir_all(&case_dir);
            return;
        }
        let job = r.expect("prepare");
        // Both files staged.
        assert!(workdir.join("wing.cfg").is_file());
        assert!(workdir.join("wing.su2").is_file());
        // Native command: SU2_CFD wing.cfg
        let cmd: Vec<String> = job
            .native_command
            .iter()
            .map(|s| s.to_string_lossy().into_owned())
            .collect();
        assert!(cmd.iter().any(|s| s == "wing.cfg"));
        // OMP_NUM_THREADS env present with the right value.
        let omp = job
            .environment
            .iter()
            .find(|(k, _)| k.to_string_lossy() == "OMP_NUM_THREADS")
            .expect("OMP_NUM_THREADS");
        assert_eq!(omp.1.to_string_lossy(), "4");
        let _ = std::fs::remove_dir_all(&case_dir);
    }

    #[test]
    fn prepare_missing_cfg_is_actionable_error() {
        let case_dir = std::env::temp_dir().join(format!(
            "valenx-su2-no-cfg-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&case_dir).unwrap();
        std::fs::write(
            case_dir.join("case.toml"),
            "[cfd.su2]\nconfig = \"missing.cfg\"\n",
        )
        .unwrap();
        let workdir = case_dir.join("workdir");
        let case = Case {
            id: "su2-test".into(),
            path: case_dir.clone(),
        };
        let r = Su2Adapter::new().prepare(&case, &workdir);
        assert!(r.is_err());
        let _ = std::fs::remove_dir_all(&case_dir);
    }

    #[test]
    fn walk_vtk_files_picks_up_both_extensions_recursively() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-su2-walk-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(tmp.join("sub")).unwrap();
        std::fs::write(tmp.join("flow.vtu"), b"placeholder").unwrap();
        std::fs::write(tmp.join("sub").join("history.vtk"), b"placeholder").unwrap();
        std::fs::write(tmp.join("history.csv"), b"not vtk").unwrap();

        let found = walk_vtk_files(&tmp);
        assert_eq!(found.len(), 2);
        // Lexicographic sort puts flow.vtu before sub/history.vtk
        assert!(found[0].ends_with("flow.vtu"));
        assert!(found[1].ends_with("history.vtk"));
        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn collect_populates_real_provenance_hashes_when_inputs_present() {
        // Set up a workdir with a .cfg + .su2 mesh; collect()
        // should hash both into Provenance::case_hash + mesh_hash.
        let tmp = std::env::temp_dir().join(format!(
            "valenx-su2-prov-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        std::fs::write(tmp.join("airfoil.cfg"), b"% SU2 cfg fixture").unwrap();
        std::fs::write(tmp.join("airfoil.su2"), b"NDIM= 2 ...").unwrap();
        let job = PreparedJob {
            workdir: tmp.clone(),
            native_command: vec!["true".into()],
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: false,
        };
        let results = Su2Adapter::new().collect(&job).expect("collect");
        // Real hashes are non-empty; the empty-sentinel from
        // stub_provenance must not leak through.
        assert_ne!(results.provenance.case_hash.0, "");
        assert_ne!(results.provenance.mesh_hash.0, "");
        assert_eq!(results.provenance.case_hash.0.len(), 64);
        assert_eq!(results.provenance.mesh_hash.0.len(), 64);
        // run_id is fresh per-call (UUIDv4 shape, 36 chars).
        assert_eq!(results.provenance.run_id.len(), 36);
        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn collect_lists_vtk_artifacts_and_loads_fields() {
        // Minimal end-to-end: write a synthesised legacy-binary
        // .vtk file under a fake workdir, build a PreparedJob
        // pointing at that workdir, and call collect(). The result
        // should list the file as a VizData artifact AND populate
        // the field catalog with the scalar.
        let tmp = std::env::temp_dir().join(format!(
            "valenx-su2-collect-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        let vtk = tmp.join("flow.vtk");
        let mut buf: Vec<u8> = Vec::new();
        buf.extend_from_slice(b"# vtk DataFile Version 3.0\n");
        buf.extend_from_slice(b"su2 smoke\n");
        buf.extend_from_slice(b"BINARY\n");
        buf.extend_from_slice(b"DATASET UNSTRUCTURED_GRID\n");
        buf.extend_from_slice(b"POINTS 4 float\n");
        for v in [
            0.0_f32, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
        ] {
            buf.extend_from_slice(&v.to_be_bytes());
        }
        buf.push(b'\n');
        buf.extend_from_slice(b"CELLS 1 5\n");
        for v in [4u32, 0, 1, 2, 3] {
            buf.extend_from_slice(&v.to_be_bytes());
        }
        buf.push(b'\n');
        buf.extend_from_slice(b"CELL_TYPES 1\n");
        buf.extend_from_slice(&10u32.to_be_bytes());
        buf.push(b'\n');
        buf.extend_from_slice(b"POINT_DATA 4\n");
        buf.extend_from_slice(b"SCALARS Mach float 1\n");
        buf.extend_from_slice(b"LOOKUP_TABLE default\n");
        for v in [0.5_f32, 0.6, 0.7, 0.8] {
            buf.extend_from_slice(&v.to_be_bytes());
        }
        buf.push(b'\n');
        std::fs::write(&vtk, buf).unwrap();

        let job = PreparedJob {
            workdir: tmp.clone(),
            native_command: vec!["true".into()],
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: false,
        };
        let results = Su2Adapter::new().collect(&job).expect("collect");
        assert_eq!(results.artifacts.len(), 1);
        assert!(results.artifacts[0].path.ends_with("flow.vtk"));
        assert!(
            results.fields.names().any(|n| n == "Mach"),
            "expected Mach field; got: {:?}",
            results.fields.names().collect::<Vec<_>>()
        );
        let _ = std::fs::remove_dir_all(&tmp);
    }

    /// Round-9 RED→GREEN: `[cfd.su2].mesh` used to be joined with
    /// bare `case.path.join`. Wrap with `confined_join`.
    #[test]
    fn prepare_rejects_mesh_traversing_outside_case_dir() {
        use valenx_test_utils::tempdir;
        let d = tempdir("su2-mesh-trav");
        std::fs::write(d.join("case.cfg"), b"MESH_FILENAME=mesh.su2\n").unwrap();
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "cfd"
solver  = "su2"

[cfd.su2]
config = "case.cfg"
mesh   = "../../etc/passwd"
"#,
        )
        .unwrap();
        let case = Case {
            id: "su2-mesh-trav".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = Su2Adapter::new().prepare(&case, &workdir).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("..") || msg.contains("stay within") || msg.contains("escape"),
            "expected confined_join rejection, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    /// Round-20 L1 RED→GREEN: an over-cap `.vtk` file in the
    /// SU2 workdir must be skipped silently (continue) rather than
    /// slurped into memory. The pre-fix `std::fs::read(&path)` would
    /// have allocated the full file size before `vtk_dispatch` even
    /// saw the magic bytes.
    #[test]
    fn collect_skips_oversize_vtk_file() {
        use std::io::Write;
        let tmp = std::env::temp_dir().join(format!(
            "valenx-su2-r20l1-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        // Sparse 5 GiB file — past the 4 GiB MAX_VTK_FILE_BYTES cap.
        // set_len gives us the size without writing 5 GiB of zeros.
        let vtk = tmp.join("oversize.vtu");
        let mut f = std::fs::File::create(&vtk).unwrap();
        f.set_len(MAX_VTK_FILE_BYTES + 1).unwrap();
        // Write a real byte at the start so `# vtk` magic checks have
        // something to point at if they ever run (they don't here —
        // the size check trips first).
        f.write_all(b"#").unwrap();
        drop(f);
        let job = PreparedJob {
            workdir: tmp.clone(),
            native_command: vec!["true".into()],
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: false,
        };
        let results = Su2Adapter::new().collect(&job).expect("collect ok");
        // Over-cap file is still listed as an Artifact (the discovery
        // walk runs before the read), but the field catalog is empty
        // because the bounded read returned Err and we `continue`d.
        assert_eq!(results.artifacts.len(), 1);
        assert!(results.artifacts[0].path.ends_with("oversize.vtu"));
        assert!(
            results.fields.names().next().is_none(),
            "an oversize VTK must not populate fields (it should have been skipped)"
        );
        let _ = std::fs::remove_dir_all(&tmp);
    }
}
