//! Phase 102 — STEP AP203 (ISO 10303-203) reader.
//!
//! ## What OCCT does
//!
//! `STEPControl_Reader::ReadFile` + `TransferRoots` ingests a STEP
//! file, validates it against the requested schema (AP203 here),
//! and resolves it to a `TopoDS_Shape`. AP203's geometric coverage
//! (faces / edges / vertices / NURBS) round-trips losslessly.
//!
//! ## v1 status
//!
//! **Honest implementation** — delegates to
//! [`valenx_step_iges::step::read`]. The truck-stepio reader accepts
//! both AP203 and AP214 files (the schemas differ in PDM metadata,
//! not geometry) so this entry point handles AP203 input correctly
//! today. The returned [`Solid`] is mesh-backed (truck-stepio
//! tessellates incoming BRep shells) — see the valenx-step-iges
//! module header for the BRep-import deferral note.

use std::path::Path;

use valenx_cad::Solid;

use crate::error::OcctExchangeError;

/// Read a STEP AP203 file from `path` and return its first solid.
///
/// # Errors
///
/// - [`OcctExchangeError::BadInput`] if `path` has no extension or an
///   extension other than `.step` / `.stp`.
/// - [`OcctExchangeError::Parse`] if the file is malformed or has no
///   shells.
/// - [`OcctExchangeError::Backend`] for other backend failures.
/// - [`OcctExchangeError::Io`] for filesystem failures.
pub fn step_ap203_reader(path: &Path) -> Result<Solid, OcctExchangeError> {
    validate_step_extension(path)?;
    valenx_step_iges::step::read(path).map_err(map_step_err)
}

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

/// Map `StepIgesError` into the right `OcctExchangeError` variant —
/// `ParseError` is a `Parse` (the *file* is bad), everything else is
/// `Backend` (or an `Io` if it was originally an io::Error inside the
/// backend's Io variant).
fn map_step_err(err: valenx_step_iges::StepIgesError) -> OcctExchangeError {
    use valenx_step_iges::StepIgesError;
    match err {
        StepIgesError::Io(e) => OcctExchangeError::Io(e),
        StepIgesError::ParseError(msg) => OcctExchangeError::parse("step file", msg),
        other => OcctExchangeError::Backend(format!("step::read: {other}")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn rejects_non_step_extension() {
        let err = step_ap203_reader(&PathBuf::from("a.obj")).unwrap_err();
        assert_eq!(err.code(), "occt_exchange.bad_input");
    }

    #[test]
    fn rejects_missing_extension() {
        let err = step_ap203_reader(&PathBuf::from("noext")).unwrap_err();
        assert_eq!(err.code(), "occt_exchange.bad_input");
    }
}
