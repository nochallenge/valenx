//! Phase 110 — IGES 5.3 full-spec reader including BREP entities
//! (186, 502, 510).
//!
//! ## What OCCT does
//!
//! `IGESControl_Reader::ReadFile` + `TransferRoots` walks the IGES
//! Directory + Parameter Data sections and assembles a
//! `TopoDS_Shape`. Coverage spans the full 1986 spec:
//!
//! - 100 (Circular Arc), 110 (Line), 116 (Point), 124 (Transformation
//!   Matrix) — the wireframe core.
//! - 126 (Rational B-Spline Curve), 128 (Rational B-Spline Surface),
//!   142 (Curve on Parametric Surface), 144 (Trimmed Surface) — the
//!   surface core.
//! - 186 (Manifold Solid B-Rep Object), 502 (Vertex List), 504
//!   (Edge List), 508 (Loop), 510 (Face), 514 (Shell) — the BRep
//!   topology core.
//! - 212 (General Note), 314 (Color Definition) — annotation.
//!
//! ## v1 status
//!
//! **Honest implementation** for the wireframe subset only —
//! delegates to [`valenx_step_iges::iges::read`]. The full BRep
//! entity family (186 / 502 / 504 / 508 / 510 / 514) needs a
//! dedicated parser; Phase 110.5 will layer that on top of the
//! existing reader using the same hand-rolled approach the rest of
//! `valenx-step-iges::iges` uses. For trimmed surfaces, callers
//! should use [`crate::iges_trimmed_surface_reader()`] today.

use std::path::Path;

use valenx_cad::Solid;

use crate::error::OcctExchangeError;

/// Read an IGES 5.3 file from `path` and return its first solid.
///
/// # Errors
///
/// - [`OcctExchangeError::BadInput`] if the extension isn't `.iges`
///   or `.igs`.
/// - [`OcctExchangeError::Parse`] for malformed IGES text.
/// - [`OcctExchangeError::Io`] for filesystem failures.
pub fn iges_5_3_reader(path: &Path) -> Result<Solid, OcctExchangeError> {
    validate_iges_extension(path)?;
    valenx_step_iges::iges::read(path).map_err(map_iges_err)
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

fn map_iges_err(err: valenx_step_iges::StepIgesError) -> OcctExchangeError {
    use valenx_step_iges::StepIgesError;
    match err {
        StepIgesError::Io(e) => OcctExchangeError::Io(e),
        StepIgesError::ParseError(msg) => OcctExchangeError::parse("iges file", msg),
        other => OcctExchangeError::Backend(format!("iges::read: {other}")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn rejects_wrong_extension() {
        let err = iges_5_3_reader(&PathBuf::from("a.step")).unwrap_err();
        assert_eq!(err.code(), "occt_exchange.bad_input");
    }
}
