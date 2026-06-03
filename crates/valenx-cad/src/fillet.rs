//! Edge filleting — currently a typed-error stub.
//!
//! truck 0.6 does NOT expose a per-edge fillet algorithm; the
//! upstream maintainers' roadmap puts that work behind the
//! still-WIP `truck-shapeops` rewrite. Rather than ship a fake
//! function that returns the unmodified solid, [`fillet_edges`]
//! always returns [`CadError::NotImplemented`]. The toolbox UI shows
//! the error string verbatim so users know to fall back to the
//! "Open in FreeCAD" escape hatch for filleting workflows.
//!
//! When truck eventually ships a fillet API we'll route to it from
//! this module without changing the public signature.

use crate::solid::{CadError, Solid};

/// Fillet all edges of `solid` with the given (uniform) radius.
///
/// Currently always returns [`CadError::NotImplemented`] — see the
/// module docs. The function signature is stable so callers (the
/// Mesh Toolbox in `valenx-app`, integration tests) don't have to
/// guess at when the capability will arrive.
pub fn fillet_edges(_solid: &Solid, radius: f64) -> Result<Solid, CadError> {
    if !radius.is_finite() {
        return Err(CadError::InvalidParam(format!(
            "fillet radius must be finite, got {radius}"
        )));
    }
    if radius <= 0.0 {
        return Err(CadError::InvalidParam(format!(
            "fillet radius must be strictly positive, got {radius}"
        )));
    }
    Err(CadError::NotImplemented {
        op: "fillet_edges",
        reason: "truck 0.6 does not expose a per-edge fillet operation \
                 — use 'Open in FreeCAD' for filleting workflows until \
                 the upstream API ships"
            .to_string(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::box_solid;

    #[test]
    fn fillet_rejects_invalid_radius() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        assert!(matches!(
            fillet_edges(&cube, -0.1),
            Err(CadError::InvalidParam(_))
        ));
        assert!(matches!(
            fillet_edges(&cube, 0.0),
            Err(CadError::InvalidParam(_))
        ));
        assert!(matches!(
            fillet_edges(&cube, f64::NAN),
            Err(CadError::InvalidParam(_))
        ));
    }

    #[test]
    fn fillet_returns_not_implemented_for_valid_radius() {
        // A positive, finite radius is structurally valid — the
        // honest "not yet" response is NotImplemented, NOT a fake
        // success or InvalidParam.
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        match fillet_edges(&cube, 0.1) {
            Err(CadError::NotImplemented { op, reason }) => {
                assert_eq!(op, "fillet_edges");
                assert!(reason.contains("truck"));
            }
            other => panic!("expected NotImplemented, got {other:?}"),
        }
    }
}
