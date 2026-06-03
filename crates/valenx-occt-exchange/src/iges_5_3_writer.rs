//! Phase 109 — IGES 5.3 full-spec writer (surfaces + curves +
//! annotations).
//!
//! ## What OCCT does
//!
//! `IGESControl_Writer::AddShape` + `Model().AddEntity` emits any
//! `TopoDS_Shape` as the appropriate IGES entity type:
//!
//! - Solid → entity 186 (Manifold Solid B-Rep Object).
//! - Shell → entity 514 (Shell), built of entities 510 (Face),
//!   508 (Loop), 504 (Edge).
//! - Surface → entity 128 (Rational B-Spline Surface) or trimmed
//!   variant 144 (Trimmed Surface) wrapping a 128.
//! - Curve → 126 (Rational B-Spline Curve) or 110 (Line) / 100 (Arc)
//!   for analytic curves.
//! - Annotations (text) → 212 (General Note) + 214 (Leader Arrow).
//!
//! ## v1 status
//!
//! **Honest implementation** for the wireframe subset only —
//! delegates to [`valenx_step_iges::iges::write`], which emits
//! entities 110 (Line), 100 (Arc), 116 (Point), and 124
//! (Transformation Matrix) in the standard 80-column IGES grid.
//! Full BRep / surface output is deferred to Phase 109.5 (which
//! will reuse the trimmed-surface writer from
//! [`crate::iges_trimmed_surface_writer()`] plus the wireframe path
//! as fallback for entities truck can't yet expose).

use std::path::Path;

use valenx_cad::Solid;

use crate::error::OcctExchangeError;

/// Write `solid` to `path` as an IGES 5.3 file.
///
/// # Errors
///
/// - [`OcctExchangeError::BadInput`] if the extension isn't `.iges`
///   or `.igs`.
/// - [`OcctExchangeError::Backend`] for backend failures.
/// - [`OcctExchangeError::Io`] for filesystem failures.
pub fn iges_5_3_writer(solid: &Solid, path: &Path) -> Result<(), OcctExchangeError> {
    validate_iges_extension(path)?;
    valenx_step_iges::iges::write(solid, path)
        .map_err(|e| OcctExchangeError::Backend(format!("iges::write: {e}")))
}

fn validate_iges_extension(path: &Path) -> Result<(), OcctExchangeError> {
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .map(str::to_ascii_lowercase);
    match ext.as_deref() {
        Some("iges") | Some("igs") => Ok(()),
        Some(other) => Err(OcctExchangeError::bad_input(
            "path",
            format!("extension must be .iges or .igs; got .{other}"),
        )),
        None => Err(OcctExchangeError::bad_input(
            "path",
            "missing extension; expected .iges or .igs",
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn rejects_wrong_extension() {
        let cube = valenx_cad::box_solid(1.0, 1.0, 1.0).unwrap();
        let err = iges_5_3_writer(&cube, &PathBuf::from("a.step")).unwrap_err();
        assert_eq!(err.code(), "occt_exchange.bad_input");
    }
}
