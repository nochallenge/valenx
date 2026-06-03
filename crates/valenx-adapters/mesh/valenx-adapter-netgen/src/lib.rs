//! # valenx-adapter-netgen
//!
//! Adapter for the Netgen unstructured mesher (Joachim Schoeberl).
//! Strong on curved-geometry tetrahedral meshing; complements gmsh.
//!
//! **Phase 2 — live for batch CSG / BREP meshing.** `prepare()`
//! parses `[meshing.netgen]` from case.toml and stages the
//! geometry source into the workdir. `run()` invokes
//! `netgen -batchmode -geofile=<src> -meshfile=<out>` (with an
//! optional `-meshsize=` from the case). `collect()` discovers
//! the produced `.vol` (or `.vol.gz`) and any input source files.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;
pub mod vol_parser;

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

use crate::case_input::NetgenInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(NetgenAdapter::new())
}

pub struct NetgenAdapter;

impl NetgenAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for NetgenAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "netgen";
const BINARIES: &[&str] = &["netgen", "ngscxx"];

impl Adapter for NetgenAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "Netgen",
            version_range: VersionRange {
                min_inclusive: Version::new(6, 2, 0),
                max_exclusive: Version::new(7, 0, 0),
            },
            physics: &[Physics::Meshing],
            license_mode: LicenseMode::Subprocess,
            tool_license: "LGPL-2.1-or-later",
            docs_url: "https://docu.ngsolve.org/nightly/",
            homepage_url: "https://ngsolve.org/",
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
                hint: "Netgen 6.2+ required; install via NGSolve distribution".into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = NetgenInput::from_case_dir(&case.path)?;
        fs::create_dir_all(workdir)?;

        // Stage the geometry source into the workdir under its
        // original name so netgen's relative-path lookups work.
        let source = valenx_core::adapter_helpers::confined_join(
            &case.path,
            &input.geometry_file,
        )?;
        if !source.is_file() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "geometry_file not found at {} (resolve relative to case dir)",
                source.display()
            )));
        }
        let staged_name = source
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("geom.geo")
            .to_string();
        let staged = workdir.join(&staged_name);
        fs::copy(&source, &staged)
            .map_err(|e| AdapterError::Other(anyhow::anyhow!("stage {}: {e}", source.display())))?;

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "no Netgen binary on PATH".into(),
        })?;

        // netgen invocation:
        //   netgen -batchmode -geofile=<src> -meshfile=<out>
        //          [-meshsize=<h>]
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from("-batchmode"),
            OsString::from(format!("-geofile={staged_name}")),
            OsString::from(format!("-meshfile={out}", out = input.output)),
        ];
        if let Some(ms) = input.mesh_size {
            native_command.push(OsString::from(format!("-meshsize={ms}")));
        }

        // Sub-second meshes are common; complex curved geometry can
        // take minutes. 5-min default is generous.
        let estimated_runtime = Some(Duration::from_secs(5 * 60));

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            estimated_runtime,
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting Netgen", |line| {
            let mut hint = subprocess::Hint::default();
            if let Some(pct) = netgen_progress_hint(line) {
                hint.progress = Some((pct, line.to_string()));
            }
            // Netgen surfaces "Error" / "ERROR" on its own diagnostics.
            if line.starts_with("Error") || line.contains("ERROR") {
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
        // case_path: Netgen's native CSG inputs (.geo for 3D, .geo2d
        // for 2D) or imported BREP (.step / .stp / .iges / .igs /
        // .brep). mesh_path: produced .vol / .vol.gz when the
        // mesher has run.
        let case_path = first_workdir_match(
            &job.workdir,
            &["geo", "geo2d", "step", "stp", "iges", "igs", "brep"],
        )
        .unwrap_or_else(|| job.workdir.join("(no-source-found)"));
        let mesh_path = first_workdir_match(&job.workdir, &["vol", "gz"]);
        let prov = valenx_core::adapter_helpers::live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "Netgen",
            "unknown",
            &case_path,
            mesh_path.as_deref(),
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);

        // Walk workdir for produced + source artifacts. Netgen
        // emits .vol (text) or .vol.gz (gzipped) as the primary
        // output; everything else is the input source we staged.
        let classifications: &[(&str, ArtifactKind, &str)] = &[
            ("vol", ArtifactKind::Native, "Netgen .vol mesh"),
            ("gz", ArtifactKind::Native, "Netgen .vol.gz (gzipped)"),
            ("geo", ArtifactKind::Other, "Netgen CSG source"),
            ("geo2d", ArtifactKind::Other, "Netgen 2D CSG source"),
            ("step", ArtifactKind::Other, "STEP geometry"),
            ("stp", ArtifactKind::Other, "STEP geometry"),
            ("iges", ArtifactKind::Other, "IGES geometry"),
            ("igs", ArtifactKind::Other, "IGES geometry"),
            ("brep", ArtifactKind::Other, "BREP geometry"),
        ];
        let mut vol_path: Option<std::path::PathBuf> = None;
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
                    if ext == "vol" && vol_path.is_none() {
                        vol_path = Some(path.clone());
                    }
                    results.artifacts.push(Artifact {
                        path,
                        kind,
                        checksum: None,
                        label: label.to_string(),
                    });
                }
            }
        }

        // If a `.vol` is present, parse it into a canonical Mesh and
        // serialise as `mesh.canonical.json` next to it. The app's
        // post-run hook auto-loads that file into the viewport, which
        // closes the meshing UX loop the same way it does for gmsh.
        if let Some(p) = &vol_path {
            match vol_parser::parse_file(p, &format!("netgen-{}", job.workdir.display())) {
                Ok(mesh) if !mesh.nodes.is_empty() => {
                    let canonical_path = job.workdir.join("mesh.canonical.json");
                    if let Ok(bytes) = serde_json::to_vec_pretty(&mesh) {
                        if valenx_core::io_caps::atomic_write_bytes(&canonical_path, &bytes).is_ok()
                        {
                            results.artifacts.push(Artifact {
                                path: canonical_path,
                                kind: ArtifactKind::VizData,
                                checksum: None,
                                label: format!(
                                    "canonical mesh · {} nodes · {} elements",
                                    mesh.stats.node_count, mesh.stats.element_count,
                                ),
                            });
                        }
                    }
                }
                Ok(_) => {
                    // Parsed but empty (placeholder file or malformed
                    // body) — don't emit a canonical mesh stub.
                }
                Err(e) => {
                    tracing::warn!(target: "valenx-netgen", ?e, ?p, "vol parse failed");
                }
            }
        }
        results.artifacts.sort_by(|a, b| a.path.cmp(&b.path));
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: vec![
                Capability::Meshing3D,
                Capability::MeshingUnstructured,
                Capability::MeshingPrismLayers,
            ],
            ribbon_contributions: vec!["mesh.netgen.generate"],
        }
    }
}

/// Coarse progress hints for Netgen stdout banners. Based on the
/// messages netgen 6.2+ prints in `-batchmode`. Non-monotonic
/// banners are deliberately skipped so the GUI bar only moves
/// forward.
fn netgen_progress_hint(line: &str) -> Option<f32> {
    if line.contains("Start Findpoints") {
        Some(15.0)
    } else if line.contains("Surface meshing") {
        Some(35.0)
    } else if line.contains("Volume meshing") {
        Some(60.0)
    } else if line.contains("Optimization") || line.contains("Optimize") {
        Some(85.0)
    } else if line.contains("Mesh successful")
        || line.contains("Statistics:")
        || line.contains("Save mesh")
    {
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
    fn info_is_well_formed() {
        let info = NetgenAdapter::new().info();
        assert_eq!(info.id, "netgen");
    }

    #[test]
    fn collect_uses_live_provenance_with_real_case_hash() {
        let workdir = tempdir("netgen-collect");
        let case_path = workdir.join("model.geo");
        let case_bytes = b"algebraic3d\nsolid box = orthobrick (0,0,0; 1,1,1);\n";
        std::fs::write(&case_path, case_bytes).expect("write .geo");

        let job = PreparedJob {
            workdir: workdir.clone(),
            native_command: Vec::new(),
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: false,
        };
        let results = NetgenAdapter::new().collect(&job).expect("collect");
        let prov = &results.provenance;

        assert_eq!(prov.adapter, INFO_ID);
        assert!(!prov.adapter_version.is_empty());
        assert_eq!(prov.tool, "Netgen");
        assert!(!prov.run_id.is_empty(), "run_id empty — stub still wired?");
        assert_eq!(prov.case_hash, sha256_hex_file(&case_path));

        cleanup(&workdir);
    }

    fn cleanup(d: &std::path::Path) {
        let _ = std::fs::remove_dir_all(d);
    }

    #[test]
    fn progress_hints_are_monotonic() {
        let pts = [
            netgen_progress_hint("Start Findpoints"),
            netgen_progress_hint("Surface meshing..."),
            netgen_progress_hint("Volume meshing..."),
            netgen_progress_hint("Optimization step 1"),
            netgen_progress_hint("Mesh successful"),
        ];
        let mut last = 0.0_f32;
        for (i, p) in pts.iter().enumerate() {
            let v = p.unwrap_or_else(|| panic!("step {i} returned None"));
            assert!(v >= last, "step {i}: {last} -> {v}");
            last = v;
        }
    }

    #[test]
    fn collect_emits_canonical_json_for_parseable_vol() {
        // A real (one-tet) Netgen .vol should round-trip through
        // vol_parser into mesh.canonical.json so the post-run hook
        // can drop it into the viewport. The placeholder test
        // (`collect_classifies_vol_and_geo_outputs`) covers the
        // negative case where the .vol body is unparseable.
        let workdir = tempdir("netgen-canonical");
        let vol = "\
mesh3d
dimension
3

points
4
0 0 0
1 0 0
0 1 0
0 0 1

volumeelements
1
1 4 1 2 3 4

endmesh
";
        std::fs::write(workdir.join("mesh.vol"), vol).unwrap();
        let job = PreparedJob {
            workdir: workdir.clone(),
            native_command: Vec::new(),
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: false,
        };
        let results = NetgenAdapter::new().collect(&job).expect("collect");
        // Both .vol and mesh.canonical.json were emitted.
        let labels: Vec<&str> = results.artifacts.iter().map(|a| a.label.as_str()).collect();
        assert!(
            labels.iter().any(|l| l.contains("canonical mesh")),
            "expected canonical-mesh artifact; got: {labels:?}"
        );
        let canonical = workdir.join("mesh.canonical.json");
        assert!(canonical.is_file(), "mesh.canonical.json must exist");
        // The serialised JSON round-trips back into Mesh and has the
        // right node count.
        let bytes = std::fs::read(&canonical).unwrap();
        let mesh: valenx_mesh::Mesh = serde_json::from_slice(&bytes).unwrap();
        assert_eq!(mesh.nodes.len(), 4);
        assert_eq!(mesh.element_blocks.len(), 1);
        cleanup(&workdir);
    }

    #[test]
    fn collect_classifies_vol_and_geo_outputs() {
        let workdir = tempdir("netgen-collect");
        for (name, content) in [
            ("mesh.vol", &b"NETGEN\n"[..]),
            ("source.geo", &b"algebraic3d\n"[..]),
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
        let results = NetgenAdapter::new().collect(&job).expect("collect");
        // .vol + .geo classified; .txt skipped.
        assert_eq!(results.artifacts.len(), 2);
        let labels: Vec<&str> = results.artifacts.iter().map(|a| a.label.as_str()).collect();
        assert!(labels.iter().any(|l| l.contains(".vol mesh")));
        assert!(labels.iter().any(|l| l.contains("CSG source")));
        cleanup(&workdir);
    }

    #[test]
    fn prepare_stages_geometry_and_builds_command() {
        if find_on_path(BINARIES).is_none() {
            // No netgen binary — verify the not-installed path.
            let case_dir = tempdir("netgen-prepare-no-bin");
            std::fs::write(
                case_dir.join("case.toml"),
                "[meshing.netgen]\ngeometry_file = \"shape.geo\"\n",
            )
            .unwrap();
            std::fs::write(case_dir.join("shape.geo"), b"algebraic3d\n").unwrap();
            let workdir = tempdir("netgen-prepare-no-bin-wd");
            let case = Case {
                id: "netgen-test".into(),
                path: case_dir.clone(),
            };
            let r = NetgenAdapter::new().prepare(&case, &workdir);
            assert!(matches!(r, Err(AdapterError::ToolNotInstalled { .. })));
            cleanup(&case_dir);
            cleanup(&workdir);
            return;
        }
        // Binary present — exercise the happy path.
        let case_dir = tempdir("netgen-prepare");
        std::fs::write(
            case_dir.join("case.toml"),
            "[meshing.netgen]\ngeometry_file = \"shape.geo\"\nmesh_size = 0.5\noutput = \"out.vol\"\n",
        )
        .unwrap();
        std::fs::write(case_dir.join("shape.geo"), b"algebraic3d\n").unwrap();
        let workdir = tempdir("netgen-prepare-wd");
        let case = Case {
            id: "netgen-test".into(),
            path: case_dir.clone(),
        };
        let job = NetgenAdapter::new()
            .prepare(&case, &workdir)
            .expect("prepare");
        // Geometry staged.
        assert!(workdir.join("shape.geo").is_file());
        // Command structure: netgen -batchmode -geofile=shape.geo
        // -meshfile=out.vol -meshsize=0.5
        let cmd: Vec<String> = job
            .native_command
            .iter()
            .map(|s| s.to_string_lossy().into_owned())
            .collect();
        assert!(cmd.iter().any(|s| s == "-batchmode"));
        assert!(cmd.iter().any(|s| s == "-geofile=shape.geo"));
        assert!(cmd.iter().any(|s| s == "-meshfile=out.vol"));
        assert!(cmd.iter().any(|s| s == "-meshsize=0.5"));
        cleanup(&case_dir);
        cleanup(&workdir);
    }
}
