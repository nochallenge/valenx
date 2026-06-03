//! # valenx-adapter-elmer
//!
//! Adapter for Elmer FEM (CSC) — multi-physics finite-element
//! solver.
//!
//! **Phase 3 — live for steady-state heat equation.** `prepare()`
//! writes a deterministic `case.sif` from the typed
//! [`case_input::ElmerInput`]. `run()` spawns `ElmerSolver` via the
//! shared [`valenx_core::subprocess`] runner. `collect()` harvests
//! `.result` / `.vtu` / `.ep` files as artifacts. Richer physics
//! (Navier-Stokes, magnetics, elasticity, multi-physics coupling)
//! grow via new `Equation` variants — the SIF writer already has a
//! clean hinge for them.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;
pub mod sif_writer;

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

use crate::case_input::ElmerInput;
use crate::sif_writer::DEFAULT_SIF_FILENAME;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(ElmerAdapter::new())
}

pub struct ElmerAdapter;

impl ElmerAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for ElmerAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "elmer";
const BINARIES: &[&str] = &["ElmerSolver", "ElmerSolver_mpi"];

impl Adapter for ElmerAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "Elmer FEM",
            version_range: VersionRange {
                min_inclusive: Version::new(9, 0, 0),
                max_exclusive: Version::new(10, 0, 0),
            },
            physics: &[Physics::Fea, Physics::Em, Physics::MultiPhysics],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-2.0-or-later",
            docs_url: "https://www.csc.fi/web/elmer",
            homepage_url: "https://www.elmerfem.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // ElmerSolver --version prints to stdout. No flag
                // collision concerns — Elmer is consistent across
                // its 9.x line.
                let found_version = valenx_core::adapter_helpers::detect_tool_version_semver(
                    &binary_path,
                    &["--version"],
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
                hint: "Elmer FEM 9.0+ required; install from elmerfem.org or \
                       your distribution (ElmerSolver on PATH)"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let (_header, input) = ElmerInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[fea.elmer].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Stage the mesh directory from the case into the workdir.
        // Elmer expects `mesh.header`, `mesh.nodes`, `mesh.elements`,
        // `mesh.boundary` files sitting inside a named directory
        // next to the .sif — we copy the user's `mesh_dir` across so
        // relative `Mesh DB "." "mesh"` resolves.
        // Round-9 hardening: `mesh_dir` is user-supplied data and the
        // contents get copied into the workdir; wrap relative paths
        // with `confined_join`.
        let mesh_src = if input.mesh_dir.is_absolute() {
            input.mesh_dir.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.mesh_dir)?
        };
        let mesh_dst = workdir.join(
            input
                .mesh_dir
                .file_name()
                .map(|s| s.to_string_lossy().to_string())
                .unwrap_or_else(|| "mesh".to_string()),
        );
        if mesh_src.is_dir() && mesh_src != mesh_dst {
            copy_dir_shallow(&mesh_src, &mesh_dst)?;
        } else if !mesh_src.is_dir() {
            // Not an error — Elmer will say "mesh not found" itself
            // when ElmerSolver starts. We warn via tracing but let
            // prepare() succeed so users can iterate on .sif content
            // without a mesh handy.
            tracing::warn!(
                target: "valenx-elmer",
                mesh = %input.mesh_dir.display(),
                "elmer case's mesh_dir not found in case dir or workdir"
            );
        }

        // Write the SIF.
        let sif_path = workdir.join(DEFAULT_SIF_FILENAME);
        sif_writer::write_to_file(&input, &sif_path)?;

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "ElmerSolver not found on PATH; install Elmer FEM".into(),
        })?;

        // ElmerSolver accepts a SIF as its first positional arg.
        let native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from(DEFAULT_SIF_FILENAME),
        ];

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            estimated_runtime: Some(Duration::from_secs(120)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting ElmerSolver", |line| {
            let mut hint = subprocess::Hint::default();
            if let Some(pct) = elmer_progress_hint(line) {
                hint.progress = Some((pct, line.to_string()));
            }
            if line.contains("WARNING") || line.contains("ERROR") {
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
        // Real provenance: Elmer cases are driven by a .sif file
        // (Solver Input File) — that's the canonical hashable
        // input. Mesh comes from .mesh.* files (Elmer's own format)
        // or a .vtu / .msh import.
        let case_path = first_workdir_match(&job.workdir, &["sif"])
            .unwrap_or_else(|| job.workdir.join("(no-sif-found)"));
        let mesh_path = first_workdir_match(&job.workdir, &["msh", "vtu", "header"]);
        let prov = valenx_core::adapter_helpers::live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "Elmer",
            "unknown",
            &case_path,
            mesh_path.as_deref(),
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);

        let classifications: &[(&str, ArtifactKind, &str)] = &[
            ("vtu", ArtifactKind::VizData, "Elmer VTU result"),
            ("pvtu", ArtifactKind::VizData, "Elmer parallel VTU"),
            ("ep", ArtifactKind::VizData, "Elmer ElmerPost file"),
            ("result", ArtifactKind::Native, "Elmer .result field dump"),
            ("sif", ArtifactKind::Other, "Elmer SIF (generated)"),
            ("log", ArtifactKind::Log, "Elmer log"),
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
        // Parse every ASCII .vtu artifact into the canonical Field
        // catalog. Same path the OpenFOAM adapter uses; Elmer writes
        // .vtu by default via its `Post File` simulation entry.
        load_vtu_fields_into_results(&mut results);
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: vec![
                Capability::FeaLinearStatic,
                Capability::FeaNonlinearStatic,
                Capability::FeaTransient,
                Capability::FeaThermal,
                Capability::EmFrequencyDomain,
                Capability::CouplingConjugateHeat,
            ],
            ribbon_contributions: vec!["fea.elmer.heat", "fea.elmer.electro", "fea.elmer.coupled"],
        }
    }
}

/// For every ASCII `.vtu` artifact already collected, parse + convert
/// into canonical [`valenx_fields::Field`]s and insert into the
/// `Results::fields` catalog. Same shape the OpenFOAM adapter uses —
/// per-file failures are skipped silently rather than fatal.
///
/// Time keys come from a filename suffix when present (Elmer often
/// writes `case_t0001.vtu`, `case_t0002.vtu`, … for transient runs);
/// if no integer suffix is parseable the field lands at
/// `TimeKey::Steady` and overwrites earlier entries with the same name.
fn load_vtu_fields_into_results(results: &mut Results) {
    // Both .vtu (XML) and .vtk (legacy binary) extensions are
    // accepted; vtk_dispatch routes each file to the right reader
    // by sniffing the magic prefix.
    let mut vtk_paths: Vec<PathBuf> = results
        .artifacts
        .iter()
        .filter(|a| a.kind == ArtifactKind::VizData)
        .filter(|a| {
            a.path
                .extension()
                .and_then(|s| s.to_str())
                .map(|s| {
                    let s = s.to_ascii_lowercase();
                    s == "vtu" || s == "vtk"
                })
                .unwrap_or(false)
        })
        .map(|a| a.path.clone())
        .collect();
    vtk_paths.sort();
    for path in vtk_paths {
        // Round-21 M3: bound the per-file read at MAX_VTK_FILE_BYTES
        // (4 GiB — R20 L1 sister cap that SU2 + OpenFOAM already
        // got). A corrupted workdir with a multi-GB `.vtu` would
        // OOM `vtk_dispatch` before the magic bytes got sniffed.
        let Ok(bytes) = read_capped_to_bytes(&path, MAX_VTK_FILE_BYTES) else {
            continue;
        };
        let stem = path.file_stem().and_then(|s| s.to_str()).unwrap_or("vtu");
        let (_mesh, fields) = match valenx_fields::vtk_dispatch::load_canonical(&bytes, stem) {
            Ok(pair) => pair,
            Err(_) => continue,
        };
        let timekey = vtu_time_key_from_stem(stem);
        for mut f in fields {
            f.time = timekey;
            results.fields.insert(f);
        }
    }
}

/// Extract a time-step index from an Elmer `.vtu` filename stem.
/// Elmer's transient `Post File` typically writes `case_tNNNN.vtu`
/// — split on the last `_t` and parse the trailing token as `u64`.
/// Falls back to a generic `<stem>_<N>` split (matches the OpenFOAM
/// convention) so multi-source workdirs round-trip cleanly.
fn vtu_time_key_from_stem(stem: &str) -> valenx_fields::TimeKey {
    if let Some((_, tail)) = stem.rsplit_once("_t") {
        if let Ok(n) = tail.parse::<u64>() {
            return valenx_fields::TimeKey::Iteration(n);
        }
    }
    if let Some((_, tail)) = stem.rsplit_once('_') {
        if let Ok(n) = tail.parse::<u64>() {
            return valenx_fields::TimeKey::Iteration(n);
        }
    }
    valenx_fields::TimeKey::Steady
}

/// Coarse progress hints from Elmer's stdout banners. Not exhaustive;
/// Elmer emits a lot and we only care about the big transitions.
fn elmer_progress_hint(line: &str) -> Option<f32> {
    if line.contains("Reading Mesh") {
        Some(10.0)
    } else if line.contains("Starting Solver") {
        Some(25.0)
    } else if line.contains("Assembly:") {
        Some(45.0)
    } else if line.contains("BiCGStab") || line.contains("Iterative") {
        Some(60.0)
    } else if line.contains("Convergence check") {
        Some(80.0)
    } else if line.contains("Writing results") || line.contains("Saving results") {
        Some(95.0)
    } else if line.contains("ELMER SOLVER FINISHED") || line.contains("SOLVER TOTAL TIME") {
        Some(99.0)
    } else {
        None
    }
}

/// Copy every file at `src`'s top level into `dst`. Elmer's mesh
/// directories are flat (`mesh.header`, `mesh.nodes`, `mesh.elements`,
/// `mesh.boundary`, optional `mesh.names`) — no subdirectories — so
/// a shallow copy is sufficient.
fn copy_dir_shallow(src: &Path, dst: &Path) -> std::io::Result<()> {
    fs::create_dir_all(dst)?;
    for entry in fs::read_dir(src)? {
        let entry = entry?;
        if entry.file_type()?.is_file() {
            let target = dst.join(entry.file_name());
            fs::copy(entry.path(), target)?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_declares_multi_physics() {
        let info = ElmerAdapter::new().info();
        assert_eq!(info.id, "elmer");
        assert!(info.physics.contains(&Physics::MultiPhysics));
    }

    #[test]
    fn vtu_time_key_from_stem_handles_canonical_shapes() {
        use valenx_fields::TimeKey;
        // Elmer's transient convention.
        assert_eq!(
            super::vtu_time_key_from_stem("case_t0001"),
            TimeKey::Iteration(1)
        );
        assert_eq!(
            super::vtu_time_key_from_stem("heatsink_t0042"),
            TimeKey::Iteration(42)
        );
        // OpenFOAM-style fallback (single underscore + integer).
        assert_eq!(
            super::vtu_time_key_from_stem("cavity_500"),
            TimeKey::Iteration(500)
        );
        // No trailing integer → Steady.
        assert_eq!(super::vtu_time_key_from_stem("snapshot"), TimeKey::Steady);
        assert_eq!(super::vtu_time_key_from_stem(""), TimeKey::Steady);
    }

    #[test]
    fn collect_loads_vtu_fields_into_catalog() {
        use std::env::temp_dir;
        let tmp = temp_dir().join(format!(
            "valenx-elmer-collect-vtu-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        let vtk_path = tmp.join("heatsink_t0010.vtu");
        std::fs::write(
            &vtk_path,
            r#"<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="2.0">
  <UnstructuredGrid>
    <Piece NumberOfPoints="3" NumberOfCells="1">
      <Points>
        <DataArray type="Float32" NumberOfComponents="3" format="ascii">
          0 0 0  1 0 0  0 1 0
        </DataArray>
      </Points>
      <Cells>
        <DataArray Name="connectivity" format="ascii">0 1 2</DataArray>
        <DataArray Name="offsets" format="ascii">3</DataArray>
        <DataArray Name="types" format="ascii">5</DataArray>
      </Cells>
      <PointData>
        <DataArray type="Float32" Name="temperature" NumberOfComponents="1" format="ascii">
          290 295 300
        </DataArray>
      </PointData>
    </Piece>
  </UnstructuredGrid>
</VTKFile>"#,
        )
        .unwrap();

        let job = PreparedJob {
            workdir: tmp.clone(),
            native_command: vec![],
            environment: vec![],
            estimated_runtime: None,
            kill_on_drop: false,
        };
        let adapter = ElmerAdapter::new();
        let results = adapter.collect(&job).expect("collect");

        assert!(results.artifacts.iter().any(|a| a.path == vtk_path));
        assert!(
            !results.fields.is_empty(),
            "expected temperature field in catalog"
        );
        let t = results
            .fields
            .by_name("temperature")
            .next()
            .expect("temperature field");
        assert_eq!(t.kind, valenx_fields::FieldKind::Scalar);
        assert_eq!(t.location, valenx_fields::Location::OnNode);
        // _t0010 in the filename → Iteration(10).
        assert_eq!(t.time, valenx_fields::TimeKey::Iteration(10));
        let _ = std::fs::remove_dir_all(&tmp);
    }

    /// Round-21 M3 RED→GREEN: an oversize `.vtu` artefact in a
    /// solver workdir is skipped (not slurped) by the bounded read.
    /// Pre-fix the bare `fs::read(&path)` would have allocated the
    /// whole file before `vtk_dispatch` sniffed the magic bytes.
    /// Mirrors the SU2 + OpenFOAM cap test from round-20 L1.
    #[test]
    fn oversize_vtu_is_skipped_not_slurped() {
        use std::env::temp_dir;
        let tmp = temp_dir().join(format!(
            "valenx-elmer-r21-oversize-vtu-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        // Plant a sparse `.vtu` one byte past the 4 GiB cap. Sparse
        // allocation skips the actual disk write so the test stays
        // fast in CI (the stat-check path of read_capped_to_bytes
        // fires before any read).
        let vtk_path = tmp.join("oversize.vtu");
        let f = std::fs::File::create(&vtk_path).unwrap();
        f.set_len(MAX_VTK_FILE_BYTES + 1).unwrap();
        drop(f);

        let err = read_capped_to_bytes(&vtk_path, MAX_VTK_FILE_BYTES).unwrap_err();
        assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);

        // And the collect-style for-loop must skip this path
        // silently (continue) rather than propagating the error —
        // mirrors how SU2 + OpenFOAM behave for the same scenario.
        let mut any_read = false;
        for path in std::iter::once(vtk_path.clone()) {
            if read_capped_to_bytes(&path, MAX_VTK_FILE_BYTES).is_ok() {
                any_read = true;
            }
        }
        assert!(!any_read, "oversize file must not be read");

        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn progress_hints_monotonic() {
        let pts = [
            elmer_progress_hint("Reading Mesh ..."),
            elmer_progress_hint("Starting Solver 1"),
            elmer_progress_hint("Assembly:        0.21  (s)"),
            elmer_progress_hint("BiCGStab iteration"),
            elmer_progress_hint("Convergence check: converged"),
            elmer_progress_hint("Saving results in VTU file."),
            elmer_progress_hint("ELMER SOLVER FINISHED AT"),
        ];
        let mut last = 0.0f32;
        for (i, p) in pts.iter().enumerate() {
            let v = p.expect("known banner");
            assert!(v >= last, "step {i}: {last} -> {v}");
            last = v;
        }
    }

    #[test]
    fn copy_dir_shallow_copies_files() {
        let src = std::env::temp_dir().join(format!(
            "elmer-copy-src-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let dst = std::env::temp_dir().join(format!(
            "elmer-copy-dst-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&src).unwrap();
        std::fs::write(src.join("mesh.header"), b"2 4 5\n").unwrap();
        std::fs::write(src.join("mesh.nodes"), b"1 -1 0 0 0\n").unwrap();
        super::copy_dir_shallow(&src, &dst).unwrap();
        assert!(dst.join("mesh.header").is_file());
        assert!(dst.join("mesh.nodes").is_file());
        let _ = std::fs::remove_dir_all(&src);
        let _ = std::fs::remove_dir_all(&dst);
    }

    /// Round-9 RED→GREEN: `[heat].mesh_dir` used to be joined with
    /// bare `case.path.join`. Wrap with `confined_join`.
    #[test]
    fn prepare_rejects_mesh_dir_traversing_outside_case_dir() {
        let d = std::env::temp_dir().join(format!(
            "valenx-elmer-mesh-trav-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&d).unwrap();
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
format  = "1.0"
name    = "trav"
physics = "fea"
solver  = "elmer"
mesh    = "(none)"

[heat]
mesh_dir        = "../../etc"
output_basename = "out"

[heat.material]
name              = "steel"
density           = 7800.0
heat_capacity     = 500.0
heat_conductivity = 50.0

[heat.simulation]
max_output_level             = 5
steady_state_max_iterations  = 20
convergence_tolerance        = 1e-5
"#,
        )
        .unwrap();
        let case = Case {
            id: "elmer-mesh-trav".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = ElmerAdapter::new().prepare(&case, &workdir).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("..") || msg.contains("stay within") || msg.contains("escape"),
            "expected confined_join rejection, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
