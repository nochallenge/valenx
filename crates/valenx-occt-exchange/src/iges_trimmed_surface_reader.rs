//! Phase 112 — IGES Type 144 trimmed-surface reader.
//!
//! ## What OCCT does
//!
//! `IGESToBRep_BRepEntity::TransferTrimmedSurface` resolves the
//! entity 144 + child 128 / 142 chain into a `TopoDS_Face` with a
//! trimming domain. The `tolerance` parameter on the directory entry
//! controls how aggressively the trim curves get healed during BRep
//! construction.
//!
//! ## v1 status
//!
//! **Honest implementation** — delegates to
//! [`valenx_step_iges::iges_trimmed::read`], the Phase 20 production
//! reader. Round-trip with [`crate::iges_trimmed_surface_writer()`]
//! preserves surface degree, knots, control points, weights, and
//! trim chains exactly.

use std::path::Path;

use valenx_step_iges::iges_trimmed::TrimmedSurfaceFile;

use crate::error::OcctExchangeError;

/// Read an IGES file at `path` and recover its trimmed-surface
/// hierarchy.
///
/// # Errors
///
/// - [`OcctExchangeError::BadInput`] if the extension isn't `.iges`
///   or `.igs`.
/// - [`OcctExchangeError::Parse`] for malformed IGES records.
/// - [`OcctExchangeError::Io`] for filesystem failures.
pub fn iges_trimmed_surface_reader(path: &Path) -> Result<TrimmedSurfaceFile, OcctExchangeError> {
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .map(str::to_ascii_lowercase);
    match ext.as_deref() {
        Some("iges") | Some("igs") => {}
        Some(other) => {
            return Err(OcctExchangeError::bad_input(
                "path",
                format!("extension must be .iges or .igs; got .{other}"),
            ));
        }
        None => {
            return Err(OcctExchangeError::bad_input(
                "path",
                "missing extension; expected .iges or .igs",
            ));
        }
    }
    valenx_step_iges::iges_trimmed::read(path).map_err(map_iges_err)
}

fn map_iges_err(err: valenx_step_iges::StepIgesError) -> OcctExchangeError {
    use valenx_step_iges::StepIgesError;
    match err {
        StepIgesError::Io(e) => OcctExchangeError::Io(e),
        StepIgesError::ParseError(msg) => OcctExchangeError::parse("iges file", msg),
        other => OcctExchangeError::Backend(format!("iges_trimmed::read: {other}")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn rejects_wrong_extension() {
        let err = iges_trimmed_surface_reader(&PathBuf::from("a.step")).unwrap_err();
        assert_eq!(err.code(), "occt_exchange.bad_input");
    }
}
