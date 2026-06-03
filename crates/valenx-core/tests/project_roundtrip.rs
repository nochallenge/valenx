//! Integration test: load the repository's `tests/fixtures/
//! minimal.valenx` project, verify its shape, round-trip it through
//! save/load in a temp dir, and check the result matches.
//!
//! Keeping this test inside the `valenx-core` crate (rather than the
//! workspace-root `tests/`) means it runs under `cargo test -p
//! valenx-core` without needing a separate test crate.

use std::fs;
use std::path::PathBuf;

use valenx_core::project::{CaseDef, LoadedProject, ProjectLoadError};
use valenx_core::workflow::{PortType, Workflow, WorkflowEdge, WorkflowNode};

/// Absolute path to the fixture at `tests/fixtures/minimal.valenx`.
fn fixture_path() -> PathBuf {
    // CARGO_MANIFEST_DIR = .../crates/valenx-core at test time.
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .expect("workspace root")
        .join("tests")
        .join("fixtures")
        .join("minimal.valenx")
}

#[test]
fn loads_minimal_fixture() {
    let project = LoadedProject::load(fixture_path()).expect("load fixture");

    // Header
    assert_eq!(project.project.project.format, "1.0");
    assert_eq!(project.project.project.name, "minimal");

    // Units defaults applied for anything missing; the fixture sets
    // them explicitly.
    assert_eq!(project.project.units.length, "m");
    assert_eq!(project.project.units.temperature, "K");

    // Geometry entry resolves as project-relative.
    assert_eq!(project.project.geometry.entries.len(), 1);
    let g = &project.project.geometry.entries[0];
    assert_eq!(g.id, "box");
    assert_eq!(g.format, "stl");

    // Cases ordered list. Steady (simpleFoam), transient (pimpleFoam),
    // FEA (calculix.static), heat (elmer.heat), gmsh-mesh, and
    // netgen-mesh demo cases all ship in the fixture so cross-crate
    // tests can exercise CFD + FEA + heat + both meshers from a
    // single project load.
    assert_eq!(
        project.case_names(),
        &[
            "box-mesh".to_string(),
            "cfd-steady".to_string(),
            "cfd-transient".to_string(),
            "fea-cantilever".to_string(),
            "heat-cube".to_string(),
            "netgen-cylinder".to_string(),
        ]
    );
    let cfd = project.cases.get("cfd-steady").expect("cfd-steady loaded");
    assert_eq!(cfd.case.physics, "cfd");
    assert_eq!(cfd.case.solver, "openfoam.simpleFoam");
    assert!(cfd.has_section("flow"));
    assert!(cfd.has_section("boundaries"));
    assert!(cfd.has_section("solve"));

    let transient = project
        .cases
        .get("cfd-transient")
        .expect("cfd-transient loaded");
    assert_eq!(transient.case.physics, "cfd");
    assert_eq!(transient.case.solver, "openfoam.pimpleFoam");
    // Transient cases carry a [solve.transient] subsection — the loader
    // doesn't validate its shape (that's the OpenFOAM adapter's job),
    // it just preserves the section so the adapter can read it.
    assert!(transient.has_section("solve"));

    let fea = project
        .cases
        .get("fea-cantilever")
        .expect("fea-cantilever loaded");
    assert_eq!(fea.case.physics, "fea");
    assert_eq!(fea.case.solver, "calculix.static");
    // The FEA structural section uses [structural] (not [flow]) —
    // verify the project loader preserves it as opaquely as it does
    // the CFD blocks.
    assert!(fea.has_section("structural"));

    // Elmer heat-conduction demo. Solver tag is `elmer.heat`; the
    // [heat] section carries the equation block the SIF writer
    // consumes (material, boundaries, simulation knobs).
    let heat = project.cases.get("heat-cube").expect("heat-cube loaded");
    assert_eq!(heat.case.physics, "fea");
    assert_eq!(heat.case.solver, "elmer.heat");
    assert!(heat.has_section("heat"));

    // gmsh box-mesh demo. Lives at the front of the order so a quick
    // "open the project and click the first case" demo path doesn't
    // require running a full CFD solve before the user can see what
    // the meshing UX feels like.
    let mesh_box = project.cases.get("box-mesh").expect("box-mesh loaded");
    assert_eq!(mesh_box.case.physics, "meshing");
    assert_eq!(mesh_box.case.solver, "gmsh.delaunay");
    assert!(mesh_box.has_section("mesh"));

    // Netgen CSG cylinder demo. The sibling `cylinder.geo` is staged
    // by the adapter at run time; the project loader doesn't read
    // it but the case.toml round-trip must preserve the
    // [meshing.netgen] block.
    let mesh_cyl = project
        .cases
        .get("netgen-cylinder")
        .expect("netgen-cylinder loaded");
    assert_eq!(mesh_cyl.case.physics, "meshing");
    assert_eq!(mesh_cyl.case.solver, "netgen.csg");
    assert!(mesh_cyl.has_section("meshing"));
}

#[test]
fn tools_lock_parses() {
    let project = LoadedProject::load(fixture_path()).expect("load fixture");
    let lock = project.tools_lock.as_ref().expect("tools.lock present");
    assert_eq!(lock.format, "1.0");
    let openfoam = lock.get("openfoam").expect("openfoam entry");
    assert_eq!(openfoam.version, "v2406");
    let gmsh = lock.get("gmsh").expect("gmsh entry");
    assert_eq!(gmsh.version, "4.12.2");
}

#[test]
fn roundtrip_save_then_reload() {
    let source = LoadedProject::load(fixture_path()).expect("load fixture");

    // Copy the fixture into a temp dir, point a new LoadedProject at
    // the copy, save it, reload. This exercises the writer without
    // mutating the checked-in fixture.
    let tmp = tempdir_for_test("valenx-core-roundtrip");
    copy_dir_recursive(&source.root, &tmp).expect("copy fixture to temp");

    // Load from the copy, save it back, then reload.
    let copy = LoadedProject::load(&tmp).expect("load copy");
    copy.save().expect("save copy");
    let reloaded = LoadedProject::load(&tmp).expect("reload copy");

    // Key invariants preserved across the roundtrip.
    assert_eq!(reloaded.project.project.name, source.project.project.name);
    assert_eq!(
        reloaded.case_names(),
        source.case_names(),
        "case order preserved"
    );
    let src_case: &CaseDef = source.cases.get("cfd-steady").unwrap();
    let dst_case: &CaseDef = reloaded.cases.get("cfd-steady").unwrap();
    assert_eq!(src_case.case.solver, dst_case.case.solver);
    assert_eq!(src_case.case.physics, dst_case.case.physics);

    cleanup(&tmp);
}

#[test]
fn rejects_missing_directory() {
    let err = LoadedProject::load("/does/not/exist/hopefully")
        .expect_err("loading a missing dir must fail");
    match err {
        ProjectLoadError::Io { .. } | ProjectLoadError::NotADirectory { .. } => {}
        other => panic!("expected Io or NotADirectory, got {other:?}"),
    }
}

#[test]
fn workflow_validate_and_order() {
    // A simple DAG that matches the fixture's conceptual pipeline.
    let wf = Workflow {
        nodes: vec![
            WorkflowNode {
                id: "geometry".into(),
                adapter_id: "valenx-stl".into(),
                inputs: Default::default(),
                outputs: [("geom".to_string(), PortType::Geometry)].into(),
                config: Default::default(),
            },
            WorkflowNode {
                id: "mesh".into(),
                adapter_id: "gmsh".into(),
                inputs: [("geom".to_string(), PortType::Geometry)].into(),
                outputs: [("mesh".to_string(), PortType::Mesh)].into(),
                config: Default::default(),
            },
            WorkflowNode {
                id: "solve".into(),
                adapter_id: "openfoam".into(),
                inputs: [("mesh".to_string(), PortType::Mesh)].into(),
                outputs: [("results".to_string(), PortType::Results)].into(),
                config: Default::default(),
            },
        ],
        edges: vec![
            WorkflowEdge {
                from: "geometry".into(),
                from_port: "geom".into(),
                to: "mesh".into(),
                to_port: "geom".into(),
            },
            WorkflowEdge {
                from: "mesh".into(),
                from_port: "mesh".into(),
                to: "solve".into(),
                to_port: "mesh".into(),
            },
        ],
    };
    wf.validate().expect("workflow validates");
    let order = wf.topo_order().expect("topo order");
    let i_geom = order.iter().position(|s| s == "geometry").unwrap();
    let i_mesh = order.iter().position(|s| s == "mesh").unwrap();
    let i_solve = order.iter().position(|s| s == "solve").unwrap();
    assert!(i_geom < i_mesh && i_mesh < i_solve);
}

// ---------------------------------------------------------------------------
// Lightweight helpers — we avoid pulling `tempfile` as a dev-dep for
// one test; the OS temp dir plus a predictable suffix is enough.
// ---------------------------------------------------------------------------

fn tempdir_for_test(name: &str) -> PathBuf {
    let mut dir = std::env::temp_dir();
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    dir.push(format!("{name}-{nanos}"));
    fs::create_dir_all(&dir).expect("create temp dir");
    dir
}

fn cleanup(dir: &std::path::Path) {
    // Best-effort cleanup. Failures do not fail the test; the OS
    // will sweep the temp dir eventually.
    let _ = fs::remove_dir_all(dir);
}

fn copy_dir_recursive(src: &std::path::Path, dst: &std::path::Path) -> std::io::Result<()> {
    fs::create_dir_all(dst)?;
    for entry in fs::read_dir(src)? {
        let entry = entry?;
        let ty = entry.file_type()?;
        let src_path = entry.path();
        let dst_path = dst.join(entry.file_name());
        if ty.is_dir() {
            copy_dir_recursive(&src_path, &dst_path)?;
        } else if ty.is_file() {
            fs::copy(&src_path, &dst_path)?;
        }
    }
    Ok(())
}
