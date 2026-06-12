//! Phase 85 — `BRepPrimAPI_MakeRevol` (revolve a profile around an
//! axis).
//!
//! ## What OCCT does
//!
//! `BRepPrimAPI_MakeRevol(profile, axis, angle = 2π)` revolves a
//! `TopoDS_Shape` (typically a `TopoDS_Wire` or `TopoDS_Face`) around
//! the supplied `gp_Ax1`. The result is a `TopoDS_Solid` for a full
//! 2π revolution of a closed planar face, or a `TopoDS_Shell` for
//! partial revolutions or open profiles.
//!
//! Used heavily for symmetric solids that aren't covered by the
//! cylinder/cone/sphere/torus primitives — bottles, lampshades,
//! mushroom-head fasteners.
//!
//! ## v1 status
//!
//! **Honest implementation** (Phase 85.5) — builds a closed planar
//! face from the XY-plane `profile_xy` polygon via
//! [`truck_modeling::builder::try_attach_plane`], then revolves it
//! with [`truck_modeling::builder::rsweep`] around the supplied axis.
//! A full 2π revolution of a closed face yields a `Solid`; partial
//! revolutions also produce a `Solid` (the two end caps are the
//! revolved face at angle 0 and angle θ). The axis is normalised
//! before being handed to truck (which requires a unit axis).

use truck_modeling::{builder, Point3, Rad, Solid as TruckSolid, Vector3, Wire};
use valenx_cad::Solid;

use crate::error::OcctSurfaceError;

/// Revolve a planar polygon profile around an axis.
///
/// `profile_xy` lies in the XY plane and is traversed in order; the
/// polygon is auto-closed (a final edge links the last point back to
/// the first), so callers pass the *distinct* corners only.
/// `axis_origin` and `axis_dir` define the rotation axis in world
/// space; `angle_rad` is the sweep amount (use `2π` for a full
/// revolution).
///
/// For the revolution to produce a valid (non-self-intersecting)
/// BRep the profile must lie entirely on one side of the axis — the
/// same constraint OCCT imposes. truck does not validate this; a
/// profile straddling the axis will produce a degenerate solid.
///
/// # Errors
///
/// - [`OcctSurfaceError::BadInput`] for fewer than 3 profile points,
///   a non-finite profile point, a zero/non-finite axis direction,
///   or a zero/non-finite sweep angle.
/// - [`OcctSurfaceError::TruckLimit`] when truck rejects the wire
///   (e.g. the profile is not planar-attachable — self-intersecting
///   polygon).
///
/// # Example
///
/// ```
/// use valenx_occt_surface::prim_api_revol;
/// // Revolve an off-axis rectangle around the Y axis → a washer-ish
/// // solid of revolution.
/// let solid = prim_api_revol(
///     &[(1.0, 0.0), (2.0, 0.0), (2.0, 1.0), (1.0, 1.0)],
///     [0.0, 0.0, 0.0],
///     [0.0, 1.0, 0.0],
///     std::f64::consts::TAU,
/// )
/// .unwrap();
/// assert!(solid.faces() > 0);
/// ```
pub fn prim_api_revol(
    profile_xy: &[(f64, f64)],
    axis_origin: [f64; 3],
    axis_dir: [f64; 3],
    angle_rad: f64,
) -> Result<Solid, OcctSurfaceError> {
    if profile_xy.len() < 3 {
        return Err(OcctSurfaceError::bad_input(
            "profile_xy",
            format!("need at least 3 points, got {}", profile_xy.len()),
        ));
    }
    for (x, y) in profile_xy {
        if !x.is_finite() || !y.is_finite() {
            return Err(OcctSurfaceError::bad_input(
                "profile_xy",
                format!("profile point ({x}, {y}) contains a non-finite value"),
            ));
        }
    }
    if axis_origin.iter().any(|c| !c.is_finite()) {
        return Err(OcctSurfaceError::bad_input(
            "axis_origin",
            "axis origin must be finite",
        ));
    }
    let dir_len = (axis_dir[0].powi(2) + axis_dir[1].powi(2) + axis_dir[2].powi(2)).sqrt();
    if !dir_len.is_finite() || dir_len < f64::EPSILON {
        return Err(OcctSurfaceError::bad_input(
            "axis_dir",
            "axis direction must be a non-zero finite vector",
        ));
    }
    if !angle_rad.is_finite() || angle_rad.abs() < f64::EPSILON {
        return Err(OcctSurfaceError::bad_input(
            "angle_rad",
            format!("must be a non-zero finite angle, got {angle_rad}"),
        ));
    }

    // Build the closed profile wire in the XY plane.
    let verts: Vec<_> = profile_xy
        .iter()
        .map(|(x, y)| builder::vertex(Point3::new(*x, *y, 0.0)))
        .collect();
    let mut edges = Vec::with_capacity(verts.len());
    for i in 0..verts.len() {
        let next = (i + 1) % verts.len();
        edges.push(builder::line(&verts[i], &verts[next]));
    }
    let wire: Wire = edges.into();

    // Attach a planar face — truck rejects non-attachable (e.g.
    // self-intersecting) wires here.
    let face = builder::try_attach_plane(&[wire])
        .map_err(|e| OcctSurfaceError::TruckLimit(format!("revol profile: {e:?}")))?;

    // truck's rsweep requires a unit axis.
    let axis = Vector3::new(
        axis_dir[0] / dir_len,
        axis_dir[1] / dir_len,
        axis_dir[2] / dir_len,
    );
    let origin = Point3::new(axis_origin[0], axis_origin[1], axis_origin[2]);

    let solid: TruckSolid = builder::rsweep(&face, origin, axis, Rad(angle_rad));
    Ok(Solid::from_truck(solid))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn revol_rejects_short_profile() {
        let err =
            prim_api_revol(&[(0.0, 0.0), (1.0, 0.0)], [0.0; 3], [0.0, 0.0, 1.0], 1.5).unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn revol_rejects_zero_axis() {
        let err = prim_api_revol(
            &[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)],
            [0.0; 3],
            [0.0; 3],
            1.5,
        )
        .unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn revol_rejects_zero_angle() {
        let err = prim_api_revol(
            &[(1.0, 0.0), (2.0, 0.0), (2.0, 1.0)],
            [0.0; 3],
            [0.0, 1.0, 0.0],
            0.0,
        )
        .unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn revol_full_revolution_builds_solid() {
        // An off-axis rectangle revolved 2π around Y is a tube solid.
        let solid = prim_api_revol(
            &[(1.0, 0.0), (2.0, 0.0), (2.0, 1.0), (1.0, 1.0)],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            std::f64::consts::TAU,
        )
        .unwrap();
        assert!(solid.faces() > 0);
        assert!(solid.vertices() > 0);
    }

    #[test]
    fn revol_partial_revolution_builds_solid() {
        let solid = prim_api_revol(
            &[(1.0, 0.0), (2.0, 0.0), (2.0, 1.0), (1.0, 1.0)],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            std::f64::consts::FRAC_PI_2,
        )
        .unwrap();
        assert!(solid.faces() > 0);
    }
}
