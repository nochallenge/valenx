//! Lock-in test: the bundled `tests/fixtures/minimal.valenx/cases/
//! cfd-steady` and `cfd-transient` cases must parse cleanly through
//! the OpenFOAM adapter's `case_input` module. Pre-fix, no test
//! exercised the bundled fixtures against the actual adapter — a
//! schema tightening could have silently dropped fields.

use std::path::PathBuf;

use valenx_adapter_openfoam::case_input::SimpleFoamInput;

fn fixture_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .ancestors()
        .nth(4)
        .expect("workspace root")
        .join("tests")
        .join("fixtures")
        .join("minimal.valenx")
        .join("cases")
}

#[test]
fn cfd_steady_fixture_parses_through_openfoam_case_input() {
    let dir = fixture_root().join("cfd-steady");
    assert!(
        dir.is_dir(),
        "fixture missing — expected {} to exist",
        dir.display()
    );
    let (header, _input) = SimpleFoamInput::from_case_dir(&dir).expect("parse cfd-steady");
    assert_eq!(header.physics, "cfd");
    assert_eq!(header.solver, "openfoam.simpleFoam");
}

#[test]
fn cfd_transient_fixture_parses_through_openfoam_case_input() {
    let dir = fixture_root().join("cfd-transient");
    assert!(dir.is_dir(), "fixture missing — {}", dir.display());
    let (header, _input) = SimpleFoamInput::from_case_dir(&dir).expect("parse cfd-transient");
    assert_eq!(header.physics, "cfd");
    assert_eq!(header.solver, "openfoam.pimpleFoam");
}
