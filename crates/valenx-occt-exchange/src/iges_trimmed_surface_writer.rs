//! Phase 111 — IGES Type 144 trimmed-surface writer.
//!
//! ## What OCCT does
//!
//! `IGESToBRep_BRepEntity::TransferTrimmedSurface` writes a
//! `TopoDS_Face` with a trimming loop as IGES entity 144 wrapping
//! an entity 128 (Rational B-Spline Surface) and a chain of entity
//! 142 (Curve on Parametric Surface) entries that define the
//! 2D trimming domain plus their 3D counterparts. The PD record
//! layout follows the IGES 5.3 spec exactly.
//!
//! ## v1 status
//!
//! **Honest implementation** — delegates to
//! [`valenx_step_iges::iges_trimmed::write`], the Phase 20
//! production writer. Accepts the same
//! [`valenx_step_iges::iges_trimmed::TrimmedSurfaceFile`] structure
//! and emits the matching IGES text. Round-trippable with
//! [`crate::iges_trimmed_surface_reader()`].

use std::path::Path;

use valenx_step_iges::iges_trimmed::TrimmedSurfaceFile;

use crate::error::OcctExchangeError;

/// Write `file` to `path` as an IGES 5.3 file containing one or
/// more trimmed-surface (entity 144) entries.
///
/// # Errors
///
/// - [`OcctExchangeError::BadInput`] if the extension isn't `.iges`
///   or `.igs`.
/// - [`OcctExchangeError::Backend`] for backend failures.
/// - [`OcctExchangeError::Io`] for filesystem failures.
pub fn iges_trimmed_surface_writer(
    file: &TrimmedSurfaceFile,
    path: &Path,
) -> Result<(), OcctExchangeError> {
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
    valenx_step_iges::iges_trimmed::write(file, path)
        .map_err(|e| OcctExchangeError::Backend(format!("iges_trimmed::write: {e}")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;
    use valenx_step_iges::iges_trimmed::TrimmedSurfaceFile;

    #[test]
    fn rejects_wrong_extension() {
        let file = TrimmedSurfaceFile::default();
        let err = iges_trimmed_surface_writer(&file, &PathBuf::from("a.step")).unwrap_err();
        assert_eq!(err.code(), "occt_exchange.bad_input");
    }
}
