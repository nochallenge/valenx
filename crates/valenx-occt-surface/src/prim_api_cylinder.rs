//! Phase 79 — `BRepPrimAPI_MakeCylinder` (axis-aligned + axis-variant
//! convenience).
//!
//! ## What OCCT does
//!
//! `BRepPrimAPI_MakeCylinder` constructs a closed cylindrical solid.
//! Three constructor overloads:
//!
//! 1. `MakeCylinder(R, H)` — base disk on XY plane, axis along +Z.
//! 2. `MakeCylinder(Ax2, R, H)` — base centred on the user's `gp_Ax2`
//!    coordinate system (origin + Z direction + X reference).
//! 3. `MakeCylinder(R, H, Angle)` — sector cylinder (open along
//!    `[0, Angle]` around the axis).
//!
//! All three return a `TopoDS_Solid` whose `IsClosed()` is true for
//! `Angle = 2π`.
//!
//! ## v1 status
//!
//! **Honest implementation** for the axis-aligned case via
//! `valenx_cad::cylinder`, and (Phase 79.5) for the **`gp_Ax2`
//! axis-variant** case — [`prim_api_cylinder_on_axis`] builds the
//! canonical +Z cylinder and re-orients it onto an arbitrary
//! `(origin, axis)` frame with a rigid rotate-then-translate. The
//! sector-cylinder (`MakeCylinder(R, H, Angle)`) overload remains a
//! follow-up.

use valenx_cad::Solid;

use crate::error::OcctSurfaceError;

/// Closed cylinder centred on the origin with axis along +Z.
///
/// # Errors
///
/// [`OcctSurfaceError::TruckLimit`] when the underlying
/// `valenx-cad` builder rejects the parameters (zero/negative
/// dimensions, non-finite inputs).
///
/// # Example
///
/// ```
/// use valenx_occt_surface::prim_api_cylinder;
/// let cyl = prim_api_cylinder(1.0, 2.0).unwrap();
/// assert!(cyl.faces() > 0);
/// ```
pub fn prim_api_cylinder(radius: f64, height: f64) -> Result<Solid, OcctSurfaceError> {
    valenx_cad::cylinder(radius, height)
        .map_err(|e| OcctSurfaceError::TruckLimit(format!("cylinder: {e:?}")))
}

/// Closed cylinder on an arbitrary `gp_Ax2` frame — OCCT's
/// `MakeCylinder(Ax2, R, H)` overload (Phase 79.5).
///
/// `axis_origin` is the centre of the base disk; `axis_dir` is the
/// direction the cylinder extrudes along (need not be unit length).
/// The base disk is perpendicular to `axis_dir`.
///
/// Implementation: build the canonical +Z cylinder via
/// [`prim_api_cylinder`], rotate it so +Z aligns with the
/// (normalised) `axis_dir`, then translate the base centre from the
/// origin to `axis_origin`. Both transforms are rigid, so the result
/// is still a closed BRep `Solid` with identical topology.
///
/// # Errors
///
/// - [`OcctSurfaceError::BadInput`] when `axis_dir` is the zero
///   vector or contains a non-finite component.
/// - [`OcctSurfaceError::TruckLimit`] when the underlying
///   `valenx-cad` cylinder builder rejects `radius` / `height`.
///
/// # Example
///
/// ```
/// use valenx_occt_surface::prim_api_cylinder::prim_api_cylinder_on_axis;
/// // Cylinder lying along +X, base at (5, 0, 0).
/// let cyl = prim_api_cylinder_on_axis(1.0, 3.0, [5.0, 0.0, 0.0], [1.0, 0.0, 0.0]).unwrap();
/// assert!(cyl.faces() > 0);
/// ```
pub fn prim_api_cylinder_on_axis(
    radius: f64,
    height: f64,
    axis_origin: [f64; 3],
    axis_dir: [f64; 3],
) -> Result<Solid, OcctSurfaceError> {
    let dir_len = (axis_dir[0].powi(2) + axis_dir[1].powi(2) + axis_dir[2].powi(2)).sqrt();
    if !dir_len.is_finite() || dir_len < f64::EPSILON {
        return Err(OcctSurfaceError::bad_input(
            "axis_dir",
            "axis direction must be a non-zero finite vector",
        ));
    }
    if axis_origin.iter().any(|c| !c.is_finite()) {
        return Err(OcctSurfaceError::bad_input(
            "axis_origin",
            "origin must be finite",
        ));
    }

    // Canonical +Z cylinder, base disk centred on the origin.
    let cyl = prim_api_cylinder(radius, height)?;

    // Unit target axis.
    let target = [
        axis_dir[0] / dir_len,
        axis_dir[1] / dir_len,
        axis_dir[2] / dir_len,
    ];

    // Rotation that carries +Z onto `target`. The rotation axis is
    // z × target; the angle is acos(z · target). When `target` is
    // already ±Z the cross product degenerates — handle both.
    let dot = target[2].clamp(-1.0, 1.0); // z · target
    let oriented = if dot > 1.0 - 1e-12 {
        // Already +Z — no rotation needed.
        cyl
    } else if dot < -1.0 + 1e-12 {
        // Antiparallel: 180° about any axis perpendicular to Z (use X).
        cyl.rotated((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), std::f64::consts::PI)
            .map_err(|e| OcctSurfaceError::TruckLimit(format!("cylinder rotate: {e}")))?
    } else {
        // General case: axis = z × target, angle = acos(dot).
        let axis = (-target[1], target[0], 0.0); // (0,0,1) × target
        let axis_len = (axis.0 * axis.0 + axis.1 * axis.1 + axis.2 * axis.2).sqrt();
        let unit_axis = (axis.0 / axis_len, axis.1 / axis_len, axis.2 / axis_len);
        cyl.rotated((0.0, 0.0, 0.0), unit_axis, dot.acos())
            .map_err(|e| OcctSurfaceError::TruckLimit(format!("cylinder rotate: {e}")))?
    };

    oriented
        .translated(axis_origin[0], axis_origin[1], axis_origin[2])
        .map_err(|e| OcctSurfaceError::TruckLimit(format!("cylinder translate: {e}")))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cylinder_builds() {
        let cyl = prim_api_cylinder(1.0, 2.0).unwrap();
        assert!(cyl.faces() > 0);
        assert!(cyl.vertices() > 0);
    }

    #[test]
    fn cylinder_rejects_zero_radius() {
        let err = prim_api_cylinder(0.0, 1.0).unwrap_err();
        assert_eq!(err.code(), "occt_surface.truck_limit");
    }

    #[test]
    fn cylinder_on_axis_builds_for_x_axis() {
        let cyl = prim_api_cylinder_on_axis(1.0, 3.0, [5.0, 0.0, 0.0], [1.0, 0.0, 0.0]).unwrap();
        // Re-orientation is a rigid transform — topology is preserved,
        // so face/vertex counts match the canonical +Z cylinder.
        let canon = prim_api_cylinder(1.0, 3.0).unwrap();
        assert_eq!(cyl.faces(), canon.faces());
        assert_eq!(cyl.vertices(), canon.vertices());
    }

    #[test]
    fn cylinder_on_axis_builds_for_z_and_antiparallel() {
        // +Z passthrough.
        let up = prim_api_cylinder_on_axis(1.0, 2.0, [0.0; 3], [0.0, 0.0, 1.0]).unwrap();
        assert!(up.faces() > 0);
        // Antiparallel (-Z).
        let down = prim_api_cylinder_on_axis(1.0, 2.0, [0.0; 3], [0.0, 0.0, -1.0]).unwrap();
        assert!(down.faces() > 0);
    }

    #[test]
    fn cylinder_on_axis_rejects_zero_axis() {
        let err = prim_api_cylinder_on_axis(1.0, 2.0, [0.0; 3], [0.0, 0.0, 0.0]).unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn cylinder_on_axis_rejects_bad_radius() {
        let err = prim_api_cylinder_on_axis(-1.0, 2.0, [0.0; 3], [0.0, 0.0, 1.0]).unwrap_err();
        assert_eq!(err.code(), "occt_surface.truck_limit");
    }
}
