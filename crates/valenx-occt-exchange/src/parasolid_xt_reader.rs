//! Phase 117 — Parasolid X_T transmit-format reader.
//!
//! ## What OCCT does
//!
//! See [`crate::parasolid_xt_writer()`] for the format overview. The
//! reader-side equivalent (`OcctSiemens::ParasolidReader::Read` in
//! a licensed integration) parses the tagged record stream into a
//! `TopoDS_Shape` via Siemens' OCCT bridge. Conversion of the
//! analytic-surface family is lossless; procedural surfaces convert
//! to a NURBS approximation at a controlled tolerance.
//!
//! ## v1 status
//!
//! **Stub** — see Phase 116.5 plan. Reader will follow the writer's
//! lossy round-trip in reverse.

use std::path::Path;

use valenx_cad::Solid;

use crate::error::OcctExchangeError;

/// Read a Parasolid X_T file from `path`.
///
/// # Errors
///
/// Always [`OcctExchangeError::NotYetImplemented`] in v1.
pub fn parasolid_xt_reader(_path: &Path) -> Result<Solid, OcctExchangeError> {
    Err(OcctExchangeError::not_yet("parasolid_xt_reader"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn parasolid_xt_reader_is_stub() {
        let err = parasolid_xt_reader(&PathBuf::from("a.x_t")).unwrap_err();
        assert_eq!(err.code(), "occt_exchange.not_yet_implemented");
    }
}
