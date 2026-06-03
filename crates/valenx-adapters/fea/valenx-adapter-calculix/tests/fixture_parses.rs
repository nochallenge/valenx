//! Lock-in test: the bundled `tests/fixtures/minimal.valenx/cases/
//! fea-cantilever/case.toml` must parse cleanly through the
//! CalculiX adapter's `case_input` module.
//!
//! Pairs with the elmer / netgen / openfoam fixture-parse tests —
//! every fixture case shipped today has a dedicated integration
//! test that locks the schema contract in place.

use std::path::PathBuf;

use valenx_adapter_calculix::case_input::LinearStaticInput;

fn fixture_case_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .ancestors()
        .nth(4)
        .expect("workspace root")
        .join("tests")
        .join("fixtures")
        .join("minimal.valenx")
        .join("cases")
        .join("fea-cantilever")
}

#[test]
fn fea_cantilever_fixture_parses_through_calculix_case_input() {
    let dir = fixture_case_dir();
    assert!(
        dir.is_dir(),
        "fixture missing — expected {} to exist",
        dir.display()
    );

    let (header, _input) = LinearStaticInput::from_case_dir(&dir).expect("parse fea-cantilever");

    assert_eq!(header.physics, "fea");
    assert_eq!(header.solver, "calculix.static");
}
