//! Lock-in test: the bundled `tests/fixtures/minimal.valenx/cases/
//! netgen-cylinder/case.toml` must parse cleanly through the Netgen
//! adapter's `case_input` module.
//!
//! Pairs with `valenx-adapter-elmer/tests/fixture_parses.rs` —
//! every fixture case shipped in this rev now has a dedicated
//! integration test that locks the schema contract in place.

use std::path::PathBuf;

use valenx_adapter_netgen::case_input::NetgenInput;

fn fixture_case_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .ancestors()
        .nth(4)
        .expect("workspace root")
        .join("tests")
        .join("fixtures")
        .join("minimal.valenx")
        .join("cases")
        .join("netgen-cylinder")
}

#[test]
fn netgen_cylinder_fixture_parses_through_netgen_case_input() {
    let dir = fixture_case_dir();
    assert!(
        dir.is_dir(),
        "fixture missing — expected {} to exist",
        dir.display()
    );

    let input = NetgenInput::from_case_dir(&dir).expect("parse netgen-cylinder");

    // Geometry source resolves to the sibling `.geo` we ship.
    assert_eq!(input.geometry_file, PathBuf::from("cylinder.geo"));
    // Mesh size is the documented 0.1 mm characteristic edge length.
    assert!((input.mesh_size.unwrap_or(0.0) - 0.1).abs() < 1e-9);
    // Output stays canonical.
    assert_eq!(input.output, "mesh.vol");

    // Sibling `.geo` should also exist — netgen needs to read it.
    let geo = dir.join("cylinder.geo");
    assert!(
        geo.is_file(),
        "expected sibling `cylinder.geo` next to case.toml at {}",
        geo.display()
    );
}
