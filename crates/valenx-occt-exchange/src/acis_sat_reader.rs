//! Phase 115 — ACIS .sat (Standard ACIS Text) reader.
//!
//! ## What OCCT does
//!
//! See [`crate::acis_sat_writer()`] for the format overview. The
//! reader-side equivalent in `OcctSpatial` is `SatReader::Read`,
//! which parses the record-oriented text into a `TopoDS_Shape` via
//! Spatial's BRep-to-OCCT bridge. The reverse path from ACIS
//! topology to OCCT topology is roughly one-to-one for analytic
//! surfaces and curves; NURBS conversion is exact since both
//! kernels use the same standard mathematical representation.
//!
//! ## v1 status
//!
//! **Stub** — see Phase 114.5 note in the writer. The reader will
//! follow the same path: parse the record stream into an
//! intermediate AST, then bridge into [`valenx_cad::Solid`] via the
//! mesh-backed fallback. Lossy but functional for legacy ACIS file
//! visualisation.

use std::path::Path;

use valenx_cad::Solid;

use crate::error::OcctExchangeError;

/// Read an ACIS .sat file from `path`.
///
/// # Errors
///
/// Always [`OcctExchangeError::NotYetImplemented`] in v1.
pub fn acis_sat_reader(_path: &Path) -> Result<Solid, OcctExchangeError> {
    Err(OcctExchangeError::not_yet("acis_sat_reader"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn acis_sat_reader_is_stub() {
        let err = acis_sat_reader(&PathBuf::from("a.sat")).unwrap_err();
        assert_eq!(err.code(), "occt_exchange.not_yet_implemented");
    }
}
