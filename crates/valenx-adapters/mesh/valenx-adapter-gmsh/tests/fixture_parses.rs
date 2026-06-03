//! Lock-in test: the bundled `tests/fixtures/minimal.valenx/cases/
//! box-mesh/case.toml` must parse cleanly through the gmsh
//! adapter's `mesh_input` module.

use std::path::PathBuf;

use valenx_adapter_gmsh::mesh_input::{Domain, MeshSpec};

fn fixture_case_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .ancestors()
        .nth(4)
        .expect("workspace root")
        .join("tests")
        .join("fixtures")
        .join("minimal.valenx")
        .join("cases")
        .join("box-mesh")
}

#[test]
fn box_mesh_fixture_parses_through_gmsh_mesh_input() {
    let dir = fixture_case_dir();
    assert!(
        dir.is_dir(),
        "fixture missing — expected {} to exist",
        dir.display()
    );

    let (header, spec) = MeshSpec::from_case_dir(&dir).expect("parse box-mesh");

    assert_eq!(header.physics, "meshing");
    assert_eq!(header.solver, "gmsh.delaunay");

    // The fixture sets up a `box` domain at the origin with side 1.
    match spec.domain {
        Domain::Box { origin, size } => {
            for axis in 0..3 {
                assert!(
                    (origin[axis] - 0.0).abs() < 1e-9,
                    "origin[{axis}] = {}",
                    origin[axis]
                );
                assert!(
                    (size[axis] - 1.0).abs() < 1e-9,
                    "size[{axis}] = {}",
                    size[axis]
                );
            }
        }
        other => panic!("expected Box domain; got {other:?}"),
    }

    // Characteristic length must round-trip. The fixture's
    // `characteristic_length = 0.1` is consumed as `char_length_max`
    // per the parser's alias chain.
    assert!(
        (spec.sizes.char_length_max - 0.1).abs() < 1e-9,
        "char_length_max = {}",
        spec.sizes.char_length_max
    );
}
