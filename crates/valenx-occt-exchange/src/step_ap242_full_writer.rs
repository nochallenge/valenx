//! Phase 105 — STEP AP242 (Managed Model-based 3D Engineering) full
//! writer with assembly + GD&T + kinematics metadata.
//!
//! ## What OCCT does
//!
//! `STEPCAFControl_Writer` with the AP242 schema flag emits a STEP
//! file conforming to the most ambitious of the still-relevant STEP
//! application protocols. AP242 carries everything AP203/AP214 does
//! plus:
//!
//! - Product Manufacturing Information (PMI): geometric dimensions
//!   and tolerances bound to specific BRep faces.
//! - Kinematics: joints, motion paths, link assemblies.
//! - Composite materials: layered laminate stack-ups.
//! - Visual rendering hints: presentation styles for downstream PDM
//!   viewers.
//!
//! ## v1 status
//!
//! **Honest hybrid implementation** — writes geometry via
//! [`valenx_step_iges::step::write`] (truck-stepio emits AP203/AP214-
//! compatible BRep that AP242 readers accept as the geometric
//! payload), then appends the caller's
//! [`valenx_step_iges::ap242::Ap242Metadata`] block via
//! [`valenx_step_iges::ap242::append_metadata`]. The metadata is
//! written as STEP-syntax comments rather than the full AP242 entity
//! graph, so an AP242 reader will round-trip the geometry but only
//! recover the metadata via the matching scan-based reader
//! ([`valenx_step_iges::ap242::parse_metadata`]). Phase 105.5 will
//! replace the comment-based metadata with the real
//! `PROPERTY_DEFINITION_REPRESENTATION` / `GEOMETRIC_TOLERANCE`
//! entity graph.

use std::path::Path;

use valenx_cad::Solid;
use valenx_step_iges::ap242::Ap242Metadata;

use crate::error::OcctExchangeError;

/// Write `solid` + `metadata` to `path` as an AP242-flavoured STEP
/// file. The metadata block follows the geometry in the file body.
///
/// # Errors
///
/// - [`OcctExchangeError::BadInput`] if the extension is wrong.
/// - [`OcctExchangeError::Backend`] when truck-stepio refuses the
///   solid or the metadata append fails.
/// - [`OcctExchangeError::Io`] for filesystem failures.
pub fn step_ap242_full_writer(
    solid: &Solid,
    metadata: &Ap242Metadata,
    path: &Path,
) -> Result<(), OcctExchangeError> {
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .map(str::to_ascii_lowercase);
    match ext.as_deref() {
        Some("step") | Some("stp") => {}
        Some(other) => {
            return Err(OcctExchangeError::bad_input(
                "path",
                format!("extension must be .step or .stp; got .{other}"),
            ));
        }
        None => {
            return Err(OcctExchangeError::bad_input(
                "path",
                "missing extension; expected .step or .stp",
            ));
        }
    }
    valenx_step_iges::step::write(solid, path)
        .map_err(|e| OcctExchangeError::Backend(format!("step::write: {e}")))?;
    if !metadata.is_empty() {
        valenx_step_iges::ap242::append_metadata(path, metadata)
            .map_err(|e| OcctExchangeError::Backend(format!("ap242::append_metadata: {e}")))?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn rejects_wrong_extension() {
        let cube = valenx_cad::box_solid(1.0, 1.0, 1.0).unwrap();
        let md = Ap242Metadata::default();
        let err = step_ap242_full_writer(&cube, &md, &PathBuf::from("a.iges")).unwrap_err();
        assert_eq!(err.code(), "occt_exchange.bad_input");
    }
}
