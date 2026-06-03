//! Phase 101 — STEP AP203 (ISO 10303-203, Configuration Controlled
//! 3D Design) writer.
//!
//! ## What OCCT does
//!
//! OCCT's `STEPControl_Writer` with `STEPControl_AsIs` mode + the
//! `write.step.schema` parameter set to `"AP203"` emits a STEP file
//! conforming to the ISO 10303-203 schema. AP203 is the narrowest of
//! the still-relevant STEP application protocols — pure CAD geometry
//! (faces, edges, vertices, B-spline surfaces, NURBS curves) with
//! product structure but no PMI / GD&T. It's the format every CAD
//! package since 1994 can read.
//!
//! ## v1 status
//!
//! **Honest implementation** — delegates to
//! [`valenx_step_iges::step::write`], which uses `truck-stepio` to
//! emit STEP text. truck-stepio's schema header is AP203/AP214-flavoured
//! ("CONFIG_CONTROL_DESIGN" / "AUTOMOTIVE_DESIGN"); the geometry it
//! writes is a strict subset that any AP203-only reader will accept,
//! so this entry point is round-trip-compatible with the
//! valenx-step-iges importer and with mainstream CAD packages.

use std::path::Path;

use valenx_cad::Solid;

use crate::error::OcctExchangeError;

/// Write `solid` to `path` as ISO 10303-203 STEP text.
///
/// # Errors
///
/// - [`OcctExchangeError::BadInput`] if `path` has no extension or an
///   extension other than `.step` / `.stp`.
/// - [`OcctExchangeError::Backend`] when `valenx-step-iges` refuses the
///   solid (e.g. mesh-backed solids, empty boundary).
/// - [`OcctExchangeError::Io`] for filesystem failures.
///
/// # Example
///
/// ```no_run
/// use std::path::PathBuf;
/// use valenx_occt_exchange::step_ap203_writer;
/// let cube = valenx_cad::box_solid(10.0, 10.0, 10.0).unwrap();
/// step_ap203_writer(&cube, &PathBuf::from("cube_ap203.step")).unwrap();
/// ```
pub fn step_ap203_writer(solid: &Solid, path: &Path) -> Result<(), OcctExchangeError> {
    validate_step_extension(path)?;
    valenx_step_iges::step::write(solid, path)
        .map_err(|e| OcctExchangeError::Backend(format!("step::write: {e}")))
}

/// Reject paths that aren't STEP — keeps the dispatcher honest and
/// surfaces typos immediately rather than letting them flow through
/// to the backend as a confusing `Unsupported` error.
fn validate_step_extension(path: &Path) -> Result<(), OcctExchangeError> {
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .map(str::to_ascii_lowercase);
    match ext.as_deref() {
        Some("step") | Some("stp") => Ok(()),
        Some(other) => Err(OcctExchangeError::bad_input(
            "path",
            format!("extension must be .step or .stp; got .{other}"),
        )),
        None => Err(OcctExchangeError::bad_input(
            "path",
            "missing extension; expected .step or .stp",
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn rejects_non_step_extension() {
        let cube = valenx_cad::box_solid(1.0, 1.0, 1.0).unwrap();
        let err = step_ap203_writer(&cube, &PathBuf::from("a.obj")).unwrap_err();
        assert_eq!(err.code(), "occt_exchange.bad_input");
    }

    #[test]
    fn rejects_missing_extension() {
        let cube = valenx_cad::box_solid(1.0, 1.0, 1.0).unwrap();
        let err = step_ap203_writer(&cube, &PathBuf::from("noext")).unwrap_err();
        assert_eq!(err.code(), "occt_exchange.bad_input");
    }

    #[test]
    fn accepts_step_and_stp() {
        assert!(validate_step_extension(&PathBuf::from("a.step")).is_ok());
        assert!(validate_step_extension(&PathBuf::from("a.STP")).is_ok());
        assert!(validate_step_extension(&PathBuf::from("a.Stp")).is_ok());
    }
}
