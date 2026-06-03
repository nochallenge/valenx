//! Phase 114 — ACIS .sat (Standard ACIS Text) writer.
//!
//! ## What OCCT does
//!
//! OCCT itself does not natively write `.sat` — ACIS is the
//! competing geometric kernel sold by Spatial Corp / Dassault.
//! OCCT users typically interop with ACIS via the third-party
//! `SAT3D` plugin or via STEP. For the purposes of feature parity
//! with the wider "OCCT data exchange" ecosystem, this entry point
//! covers what `OcctSpatial::SatWriter` would do in a fully wired
//! installation:
//!
//! ASCII text starting with a five-line header (`<version> <records>
//! <fileinfo> <units> ...`) followed by space-separated record lines
//! of the form `<index> <type> <props> ;`. Each record describes
//! one ACIS topology entity (body, lump, shell, face, loop, edge,
//! vertex, point, surface, curve). The format is proprietary and
//! Spatial actively guards the spec.
//!
//! ## v1 status
//!
//! **Stub** — emitting a `.sat` requires a reverse-engineered
//! subset of the Spatial spec. Phase 114.5 will deliver this by
//! tessellating the BRep to a triangle mesh (we already have that
//! via `valenx-mesh`) then writing a degraded `.sat` whose only
//! topology entities are the mesh facets-as-faces. That's enough to
//! ingest in legacy ACIS viewers but loses surface continuity;
//! callers wanting full round-trip should keep using STEP.

use std::path::Path;

use valenx_cad::Solid;

use crate::error::OcctExchangeError;

/// Write `solid` to `path` as ACIS .sat ASCII.
///
/// # Errors
///
/// Always [`OcctExchangeError::NotYetImplemented`] in v1.
pub fn acis_sat_writer(_solid: &Solid, _path: &Path) -> Result<(), OcctExchangeError> {
    Err(OcctExchangeError::not_yet("acis_sat_writer"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn acis_sat_writer_is_stub() {
        let cube = valenx_cad::box_solid(1.0, 1.0, 1.0).unwrap();
        let err = acis_sat_writer(&cube, &PathBuf::from("a.sat")).unwrap_err();
        assert_eq!(err.code(), "occt_exchange.not_yet_implemented");
    }
}
