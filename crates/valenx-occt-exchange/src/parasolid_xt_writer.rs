//! Phase 116 — Parasolid X_T transmit-format writer.
//!
//! ## What OCCT does
//!
//! Parasolid is the geometric kernel behind Siemens NX, Solid Edge,
//! SOLIDWORKS, Onshape, and IronCAD. OCCT does not natively write
//! `.x_t` — the transmit format is owned by Siemens, who license
//! Parasolid integration to commercial CAD vendors. OCCT users
//! typically interop with Parasolid via STEP. The closest
//! comparable functionality is `OcctSiemens::ParasolidWriter` in a
//! fully licensed integration.
//!
//! On the wire `.x_t` is a text-based schema-versioned format with
//! a four-line header followed by tagged record blocks for each
//! topology entity (`BODY`, `SHELL`, `FACE`, `LOOP`, `EDGE`,
//! `VERTEX`, `POINT`, `SURFACE`, `CURVE`). It supports all major
//! surface types Parasolid emits: B-splines, analytic surfaces
//! (plane, cylinder, cone, sphere, torus), blends, swept surfaces,
//! and procedural surfaces. Binary `.x_b` is the same content with
//! length-prefixed records.
//!
//! ## v1 status
//!
//! **Stub** — emitting `.x_t` requires the Parasolid spec.
//! Phase 116.5 will deliver this by going through a STEP roundtrip
//! plus a STEP-to-Parasolid-text rewriter that handles the analytic
//! surface family and falls back to NURBS for everything else.
//! Mainstream CAD packages will accept the result though the round-
//! trip will be lossy for non-NURBS surfaces.

use std::path::Path;

use valenx_cad::Solid;

use crate::error::OcctExchangeError;

/// Write `solid` to `path` as Parasolid X_T text.
///
/// # Errors
///
/// Always [`OcctExchangeError::NotYetImplemented`] in v1.
pub fn parasolid_xt_writer(_solid: &Solid, _path: &Path) -> Result<(), OcctExchangeError> {
    Err(OcctExchangeError::not_yet("parasolid_xt_writer"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn parasolid_xt_writer_is_stub() {
        let cube = valenx_cad::box_solid(1.0, 1.0, 1.0).unwrap();
        let err = parasolid_xt_writer(&cube, &PathBuf::from("a.x_t")).unwrap_err();
        assert_eq!(err.code(), "occt_exchange.not_yet_implemented");
    }
}
