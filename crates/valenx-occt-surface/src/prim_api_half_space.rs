//! Phase 87 — `BRepPrimAPI_MakeHalfSpace`.
//!
//! ## What OCCT does
//!
//! `BRepPrimAPI_MakeHalfSpace(face, ref_point)` builds a topological
//! half-space — an unbounded `TopoDS_Solid` bounded by the supplied
//! `TopoDS_Face`. The reference point disambiguates which side of
//! the face is "inside" the half-space.
//!
//! Half-spaces are unusual: they're infinite and so can't be
//! tessellated directly, but they're indispensable for boolean
//! operations that want to "cut off everything below a plane" without
//! constructing a giant bounding box first. The classic usage:
//!
//! ```text
//! Solid result = BRepAlgoAPI_Cut(part, half_space)
//! ```
//!
//! cuts the part by an arbitrary plane in one call.
//!
//! ## v1 status
//!
//! **Honest implementation** (Phase 87.5). truck's `Solid` enum has
//! no unbounded-solid variant, so a true *infinite* half-space is not
//! representable. Instead — exactly as OCCT's own boolean code does
//! internally when it needs a bounded approximation — this builds a
//! very large **oriented box** (`half_extent` controls the size; the
//! default [`HALF_SPACE_EXTENT`] is large enough to swallow any
//! real-world model). One face of the box lies on the bounding
//! plane; the box extends toward the half-space interior. The result
//! is a normal closed BRep `Solid`, so it composes with the existing
//! `valenx_cad::difference` / `intersection` boolean operators with
//! no special-casing.

use truck_modeling::{builder, Point3, Vector3, Wire};
use valenx_cad::Solid;

use crate::error::OcctSurfaceError;

/// Edge length used for the bounded half-space approximation when the
/// caller does not supply one. 1e9 model units — large enough that
/// any practical part sits well inside it, small enough that the
/// `f64` geometry kernel stays numerically well-conditioned (a
/// literal `f64::MAX / 2` would overflow truck's curve evaluation).
pub const HALF_SPACE_EXTENT: f64 = 1.0e9;

/// Build a (bounded-approximation) half-space bounded by the plane
/// through `plane_point` with outward normal `plane_normal`.
///
/// The resulting solid occupies the side of the plane the normal
/// points **away from** — i.e. `plane_normal` is the outward normal
/// of the bounding face, the same convention as OCCT, where the
/// reference point lies on the interior side. Pass the negated
/// normal to flip which side is solid.
///
/// `half_extent` is the box half-edge length; pass
/// [`HALF_SPACE_EXTENT`] for the default. The box is centred so its
/// `+normal` face lies exactly on the plane.
///
/// # Errors
///
/// [`OcctSurfaceError::BadInput`] when `plane_normal` is the zero
/// vector, any input is non-finite, or `half_extent` is not strictly
/// positive.
///
/// # Example
///
/// ```
/// use valenx_occt_surface::prim_api_half_space::{prim_api_half_space, HALF_SPACE_EXTENT};
/// // Half-space below the XY plane (normal points +Z, solid is -Z).
/// let hs = prim_api_half_space([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], HALF_SPACE_EXTENT)
///     .unwrap();
/// assert_eq!(hs.faces(), 6);
/// ```
pub fn prim_api_half_space(
    plane_point: [f64; 3],
    plane_normal: [f64; 3],
    half_extent: f64,
) -> Result<Solid, OcctSurfaceError> {
    if plane_point.iter().any(|c| !c.is_finite()) {
        return Err(OcctSurfaceError::bad_input(
            "plane_point",
            "plane point must be finite",
        ));
    }
    let n_len =
        (plane_normal[0].powi(2) + plane_normal[1].powi(2) + plane_normal[2].powi(2)).sqrt();
    if !n_len.is_finite() || n_len < f64::EPSILON {
        return Err(OcctSurfaceError::bad_input(
            "plane_normal",
            "plane normal must be a non-zero finite vector",
        ));
    }
    if !half_extent.is_finite() || half_extent <= 0.0 {
        return Err(OcctSurfaceError::bad_input(
            "half_extent",
            format!("must be strictly positive and finite, got {half_extent}"),
        ));
    }

    // Unit outward normal.
    let n = [
        plane_normal[0] / n_len,
        plane_normal[1] / n_len,
        plane_normal[2] / n_len,
    ];

    // Build an orthonormal in-plane basis (u, v) ⟂ n.
    let (u, v) = in_plane_basis(n);

    // The bounding face is a square of side 2*half_extent centred on
    // `plane_point`, in the (u, v) plane. We sweep it by
    // -2*half_extent * n so the box body lies on the -normal
    // (interior) side, with the original face flush with the plane.
    let e = half_extent;
    let corners = [
        offset(plane_point, u, -e, v, -e),
        offset(plane_point, u, e, v, -e),
        offset(plane_point, u, e, v, e),
        offset(plane_point, u, -e, v, e),
    ];
    let verts: Vec<_> = corners
        .iter()
        .map(|c| builder::vertex(Point3::new(c[0], c[1], c[2])))
        .collect();
    let mut edges = Vec::with_capacity(4);
    for i in 0..4 {
        edges.push(builder::line(&verts[i], &verts[(i + 1) % 4]));
    }
    let wire: Wire = edges.into();
    let face = builder::try_attach_plane(&[wire])
        .map_err(|err| OcctSurfaceError::TruckLimit(format!("half-space face: {err:?}")))?;

    // Sweep into the interior (-normal) by the full box depth.
    let depth = Vector3::new(-n[0] * 2.0 * e, -n[1] * 2.0 * e, -n[2] * 2.0 * e);
    let solid = builder::tsweep(&face, depth);
    Ok(Solid::from_truck(solid))
}

/// Two orthonormal vectors spanning the plane perpendicular to unit
/// vector `n`. Picks the world axis least aligned with `n` as a seed
/// so the cross products never degenerate.
fn in_plane_basis(n: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    // Seed = world axis with smallest |component of n|.
    let ax = n[0].abs();
    let ay = n[1].abs();
    let az = n[2].abs();
    let seed = if ax <= ay && ax <= az {
        [1.0, 0.0, 0.0]
    } else if ay <= az {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };
    let u = normalize(cross(n, seed));
    let v = cross(n, u); // already unit: n ⟂ u, both unit
    (u, v)
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn normalize(a: [f64; 3]) -> [f64; 3] {
    let l = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
    [a[0] / l, a[1] / l, a[2] / l]
}

fn offset(p: [f64; 3], u: [f64; 3], su: f64, v: [f64; 3], sv: f64) -> [f64; 3] {
    [
        p[0] + u[0] * su + v[0] * sv,
        p[1] + u[1] * su + v[1] * sv,
        p[2] + u[2] * su + v[2] * sv,
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn half_space_below_xy_plane_is_a_box() {
        let hs = prim_api_half_space([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], HALF_SPACE_EXTENT)
            .unwrap();
        // A swept square is a closed box: 6 faces, 12 edges, 8 verts.
        assert_eq!(hs.faces(), 6);
        assert_eq!(hs.edges(), 12);
        assert_eq!(hs.vertices(), 8);
    }

    #[test]
    fn half_space_arbitrary_normal_builds() {
        let hs = prim_api_half_space([1.0, 2.0, 3.0], [1.0, 1.0, 1.0], 100.0).unwrap();
        assert_eq!(hs.faces(), 6);
    }

    #[test]
    fn half_space_rejects_zero_normal() {
        let err = prim_api_half_space([0.0; 3], [0.0; 3], 1.0).unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn half_space_rejects_bad_extent() {
        let err = prim_api_half_space([0.0; 3], [0.0, 0.0, 1.0], 0.0).unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
        let err = prim_api_half_space([0.0; 3], [0.0, 0.0, 1.0], -5.0).unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }
}
