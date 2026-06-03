//! Phase 118 — Siemens JT (Jupiter Tessellation) writer.
//!
//! ## What OCCT does
//!
//! JT is Siemens' lightweight visualisation + collaboration format,
//! standardised as ISO 14306. The format ships geometry as a hybrid
//! payload: optional precise BRep (XT BREP segment, derived from the
//! Parasolid kernel) plus mandatory tessellated LODs for fast
//! viewport display. A JT file is a binary container with a tree of
//! segments (TOC + LOD pyramid + PMI + metadata).
//!
//! OCCT itself does not ship a writer; the Siemens / Mentor
//! Graphics JT Open Toolkit is the canonical implementation. OEM
//! integrations (PLM systems, MBOM tools) commonly pair OCCT with
//! the JT Open Toolkit through a thin bridge.
//!
//! ## v1 status
//!
//! **Stub** — implementing the JT segment hierarchy correctly is a
//! multi-week project against a 600-page spec. Phase 118.5 will
//! deliver a degraded "tessellated LOD only" writer that emits a
//! single-LOD JT containing the triangle mesh of `solid` (via
//! `valenx-mesh`'s STL backend, retargeted). That's enough for the
//! 80% case (viewport display in Teamcenter / JT2Go) without the
//! BRep round-trip.

use std::path::Path;

use valenx_cad::Solid;

use crate::error::OcctExchangeError;

/// Write `solid` to `path` as a JT (.jt) file.
///
/// # Errors
///
/// Always [`OcctExchangeError::NotYetImplemented`] in v1.
pub fn jt_writer(_solid: &Solid, _path: &Path) -> Result<(), OcctExchangeError> {
    Err(OcctExchangeError::not_yet("jt_writer"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn jt_writer_is_stub() {
        let cube = valenx_cad::box_solid(1.0, 1.0, 1.0).unwrap();
        let err = jt_writer(&cube, &PathBuf::from("a.jt")).unwrap_err();
        assert_eq!(err.code(), "occt_exchange.not_yet_implemented");
    }
}
