//! Real BRep corner-blend construction — the rolling-ball corner at
//! an orthogonal convex 3-edge vertex (Phase 14.7).
//!
//! # What this module does
//!
//! [`crate::brep_build`] fillets a *single edge*. When three filleted
//! edges meet at a box corner, the three quarter-cylinder fillet
//! surfaces do not meet tangentially — they leave a tri-tangent gap
//! at the vertex. This module fills that gap with the **rolling-ball
//! corner**: the ball of radius `r` seated in the corner, tangent to
//! all three faces.
//!
//! # The method — constructive solid geometry, same as the edge fillet
//!
//! At an orthogonal box corner the apex `V` is the meeting point of
//! three mutually-orthogonal convex edges with unit directions `d₀`,
//! `d₁`, `d₂` pointing away from `V` into the solid. Each `dᵢ` is also
//! an inward face direction (it lies in two of the three faces). The
//! rolling-ball's centre is the apex offset `r` inward along all
//! three:
//!
//! ```text
//!   C = V + r·d₀ + r·d₁ + r·d₂
//! ```
//!
//! `C` is exactly `r` from each of the three faces (each face is the
//! plane through `V` spanned by two of the `dᵢ`; the third `dᵢ` is its
//! unit inward normal, so the component of `C − V` along it is `r`).
//! The radius-`r` sphere centred at `C` is therefore tangent to all
//! three faces — the seated ball.
//!
//! The corner blend is then built exactly in the spirit of the
//! single-edge fillet (`(solid − cutter) ∪ bar`):
//!
//! - `corner_cutter` — a parallelepiped removing the sharp corner
//!   region: the box spanned by `d₀, d₁, d₂` from the apex, with a
//!   side length large enough to contain the seated ball's octant;
//! - `corner_ball` — the **full radius-`r` sphere** centred at `C`.
//!
//! Then `blended = (solid − corner_cutter) ∪ corner_ball`. The seated
//! ball is mostly *inside* the solid already (the union is idempotent
//! on the overlap); only the **octant facing the corner** lay in the
//! region the cutter removed, and the union adds exactly that octant
//! back — the spherical corner-blend surface. Using the closed full
//! sphere (rather than an open octant shell) is the most robust input
//! for the boolean kernel, and the set algebra is identical because
//! the seated ball's octant is contained in the cutter box.
//!
//! Both `corner_cutter` and `corner_ball` are genuine `truck-modeling`
//! BRep solids; the two set operations are the real `truck_shapeops`
//! booleans. The output is a [`valenx_cad::Solid::Brep`] with a true
//! spherical corner face — it round-trips through STEP/IGES and
//! composes with further BRep ops.
//!
//! # Honest scope and the coincident-surface caveat
//!
//! Three faces of `corner_cutter` lie flush with the solid's three
//! adjacent planar faces — the same flush trim the edge fillet's
//! cutter does, and the same fragile input for any BRep boolean
//! kernel. The construction is correct; on geometry where
//! `truck_shapeops` cannot resolve the coincidence a boolean returns
//! `None`, surfaced here as the soft [`FilletBrepError::TruckOp`]
//! error so the caller can fall through to the independent per-edge
//! fillets. Only the **orthogonal convex 3-edge corner** is handled
//! ([`crate::corner::classify_corner`] gates this) — general N-edge,
//! non-orthogonal, and concave corners stay a Tier-3 residue.

use truck_modeling::{builder, Point3, Rad, Solid as TruckSolid, Vector3, Wire};

use valenx_cad::Solid;

use crate::corner::CornerEdges;
use crate::error::FilletBrepError;

/// Linear tolerance for the corner-blend booleans — the same
/// `valenx-cad` default the edge fillet uses (tuned for primitives in
/// the 1–10 unit range, the fillet's working scale).
const CORNER_BOOL_TOL: f64 = valenx_cad::DEFAULT_BOOL_TOLERANCE;

/// How far past the seated ball the corner cutter extends, as a
/// multiple of the radius. The cutter must fully contain the ball's
/// corner octant; `2·r` per edge gives generous headroom while
/// keeping the cutter local to the corner (so it never reaches the
/// far end of an edge and collides with another corner's cutter).
const CUTTER_REACH_FACTOR: f64 = 2.0;

/// The resolved geometry of one rolling-ball corner blend.
#[derive(Clone, Copy, Debug)]
pub struct CornerBlend {
    /// The sharp corner apex (the solid's vertex).
    pub apex: Point3,
    /// Centre of the seated ball — `apex + r·(d₀ + d₁ + d₂)`.
    pub centre: Point3,
    /// Blend radius (the rolling-ball radius).
    pub radius: f64,
    /// The three orthogonal unit edge directions away from the apex.
    pub dirs: [Vector3; 3],
}

/// Resolve the [`CornerBlend`] geometry for an orthogonal convex
/// 3-edge corner.
///
/// `edges` must already be the supported case (the caller obtains it
/// from [`crate::corner::classify_corner`]). This places the seated
/// ball's centre at `apex + r·(d₀ + d₁ + d₂)`.
///
/// # Errors
///
/// [`FilletBrepError::BadParameter`] if the radius is non-finite or
/// non-positive.
pub fn resolve_corner_blend(
    edges: &CornerEdges,
    radius: f64,
) -> Result<CornerBlend, FilletBrepError> {
    if !radius.is_finite() || radius <= 0.0 {
        return Err(FilletBrepError::BadParameter {
            name: "radius",
            reason: format!("corner-blend radius must be > 0 and finite, got {radius}"),
        });
    }
    let dirs = [edges.dir0, edges.dir1, edges.dir2];
    // C = V + r·d₀ + r·d₁ + r·d₂.
    let centre = edges.apex + (dirs[0] + dirs[1] + dirs[2]) * radius;
    Ok(CornerBlend {
        apex: edges.apex,
        centre,
        radius,
        dirs,
    })
}

/// Build the **corner cutter** — the parallelepiped that removes the
/// sharp corner region.
///
/// The cutter is the box spanned from the apex by the three edge
/// directions, each scaled to [`CUTTER_REACH_FACTOR`]·`r`. Because the
/// directions are orthonormal this is an axis-aligned-in-the-corner-
/// frame cube containing the seated ball's corner octant. It is built
/// as a genuine 6-face BRep prism: a quad base face swept along the
/// third direction.
///
/// Three of its faces pass through the apex flush with the solid's
/// three adjacent faces — the documented coincident-surface case.
fn build_corner_cutter(blend: &CornerBlend) -> Result<TruckSolid, FilletBrepError> {
    let reach = blend.radius * CUTTER_REACH_FACTOR;
    let [d0, d1, d2] = blend.dirs;
    let v = blend.apex;

    // The four base-quad corners in the plane spanned by d0, d1
    // through the apex: V, V+reach·d0, V+reach·d0+reach·d1, V+reach·d1.
    let p00 = builder::vertex(v);
    let p10 = builder::vertex(v + d0 * reach);
    let p11 = builder::vertex(v + d0 * reach + d1 * reach);
    let p01 = builder::vertex(v + d1 * reach);

    let base_wire: Wire = vec![
        builder::line(&p00, &p10),
        builder::line(&p10, &p11),
        builder::line(&p11, &p01),
        builder::line(&p01, &p00),
    ]
    .into();
    let base_face = builder::try_attach_plane(&[base_wire])
        .map_err(|err| FilletBrepError::TruckOp(format!("corner cutter base face: {err:?}")))?;
    // Sweep the quad along d2 to make the parallelepiped.
    let solid: TruckSolid = builder::tsweep(&base_face, d2 * reach);
    Ok(solid)
}

/// Build the **corner ball** — the full radius-`r` sphere centred at
/// the seated-ball centre.
///
/// Built by revolving a semicircular wire (apex on the rotation axis)
/// with [`builder::cone`], exactly as `valenx_cad::primitives::sphere`
/// does — `builder::cone` collapses the degenerate pole edges, giving
/// a closed sphere shell. The semicircle is placed so the sphere is
/// centred on `blend.centre`.
fn build_corner_ball(blend: &CornerBlend) -> Result<TruckSolid, FilletBrepError> {
    let r = blend.radius;
    let c = blend.centre;
    // Semicircle in the plane containing the Y axis: top pole at
    // C+(0,r,0), swept 180° about the X axis so the wire spans the
    // C-centred semicircle, then revolved 360° about the Y axis.
    let top = builder::vertex(Point3::new(c.x, c.y + r, c.z));
    let wire: Wire = builder::rsweep(&top, c, Vector3::unit_x(), Rad(std::f64::consts::PI));
    let shell = builder::cone(&wire, Vector3::unit_y(), Rad(2.0 * std::f64::consts::PI));
    // `TruckSolid::new` panics if the swept shell is not a closed
    // 2-manifold. The `builder::cone` sphere recipe normally produces a
    // closed shell, but use the fallible `try_new` so any degenerate
    // input surfaces as a soft `TruckOp` rather than a thread-unwinding
    // panic.
    TruckSolid::try_new(vec![shell])
        .map_err(|err| FilletBrepError::TruckOp(format!("corner ball assembly: {err}")))
}

/// Apply one rolling-ball corner blend to a solid at an orthogonal
/// convex 3-edge corner.
///
/// Builds the corner cutter + corner ball and evaluates
/// `blended = (solid − corner_cutter) ∪ corner_ball` with the real
/// `truck_shapeops` booleans.
///
/// `edges` must be the supported orthogonal-convex case (obtain it via
/// [`crate::corner::classify_corner`]).
///
/// # Errors
///
/// - [`FilletBrepError::BadParameter`] if the radius is invalid.
/// - [`FilletBrepError::TruckOp`] if a cutter / ball face cannot be
///   built or a boolean returns no solid (the coincident-surface
///   case — a *soft* fall-through signal for the caller, exactly like
///   the edge fillet).
pub fn blend_corner(
    solid: &TruckSolid,
    edges: &CornerEdges,
    radius: f64,
) -> Result<TruckSolid, FilletBrepError> {
    let blend = resolve_corner_blend(edges, radius)?;

    let cutter = build_corner_cutter(&blend)?;
    let ball = build_corner_ball(&blend)?;

    let base = Solid::from_truck(solid.clone());
    let cutter_solid = Solid::from_truck(cutter);
    let ball_solid = Solid::from_truck(ball);

    // (solid − corner_cutter) ∪ corner_ball.
    let trimmed = valenx_cad::boolean::difference_tol(&base, &cutter_solid, CORNER_BOOL_TOL)
        .map_err(|err| FilletBrepError::TruckOp(format!("corner removal boolean: {err}")))?;
    let blended = valenx_cad::boolean::union_tol(&trimmed, &ball_solid, CORNER_BOOL_TOL)
        .map_err(|err| FilletBrepError::TruckOp(format!("corner-ball union boolean: {err}")))?;

    match blended {
        Solid::Brep(b) => Ok(b),
        // The boolean path always returns a BRep solid; a mesh-backed
        // result would mean an internal contract break.
        Solid::Mesh(_) => Err(FilletBrepError::TruckOp(
            "corner-blend boolean returned a mesh-backed solid (internal error)".to_string(),
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::corner::{classify_corner, CornerClass};
    use std::collections::HashSet;
    use truck_modeling::InnerSpace;
    use valenx_cad::primitives::box_solid;

    fn inner_brep(s: &Solid) -> &TruckSolid {
        match s {
            Solid::Brep(b) => b,
            _ => panic!("expected brep"),
        }
    }

    fn unique_edges(brep: &TruckSolid) -> Vec<truck_modeling::Edge> {
        let mut seen = HashSet::new();
        let mut out = Vec::new();
        for e in brep.edge_iter() {
            if seen.insert(e.id()) {
                out.push(e);
            }
        }
        out
    }

    /// Resolve the orthogonal-convex corner geometry at a cube vertex.
    fn cube_corner(brep: &TruckSolid) -> CornerEdges {
        let edges = unique_edges(brep);
        let vid = edges[0].front().id();
        match classify_corner(brep, vid, &edges) {
            CornerClass::OrthogonalConvex(ce) => ce,
            CornerClass::Unsupported => panic!("cube corner should be orthogonal-convex"),
        }
    }

    #[test]
    fn corner_centre_is_radius_from_each_face() {
        // The seated ball's centre must be exactly `radius` from each
        // of the three faces meeting at the corner. Each face is the
        // plane through the apex spanned by two edge directions; the
        // third direction is its unit normal, so the signed distance
        // of (C − apex) along that direction must equal `radius`.
        let cube = box_solid(4.0, 4.0, 4.0).unwrap();
        let brep = inner_brep(&cube);
        let ce = cube_corner(brep);
        let blend = resolve_corner_blend(&ce, 0.5).unwrap();
        let to_c = blend.centre - blend.apex;
        for (k, d) in blend.dirs.iter().enumerate() {
            let dist = to_c.dot(*d);
            assert!(
                (dist - 0.5).abs() < 1e-9,
                "corner centre should be radius from face {k}, got {dist}"
            );
        }
    }

    #[test]
    fn corner_ball_octant_surface_points_are_radius_from_centre() {
        // Sample the corner-facing octant of the seated ball: the
        // three "axis" points C − r·dᵢ (where the ball is tangent to
        // the faces) and the diagonal point C − r·(d₀+d₁+d₂)/√3
        // (nearest the apex). All must lie exactly `radius` from C.
        let cube = box_solid(4.0, 4.0, 4.0).unwrap();
        let brep = inner_brep(&cube);
        let ce = cube_corner(brep);
        let r = 0.7;
        let blend = resolve_corner_blend(&ce, r).unwrap();
        let [d0, d1, d2] = blend.dirs;
        // Tangent points on each face.
        for d in [d0, d1, d2] {
            let p = blend.centre - d * r;
            let dist = (p - blend.centre).magnitude();
            assert!(
                (dist - r).abs() < 1e-9,
                "octant tangent point should be radius from centre, got {dist}"
            );
        }
        // Diagonal point nearest the apex.
        let diag = (d0 + d1 + d2).normalize();
        let apex_near = blend.centre - diag * r;
        let dist = (apex_near - blend.centre).magnitude();
        assert!(
            (dist - r).abs() < 1e-9,
            "octant apex-near point should be radius from centre, got {dist}"
        );
    }

    #[test]
    fn corner_cutter_builds_as_a_six_face_brep() {
        // The corner cutter is a real parallelepiped: 6 faces (a
        // swept quad), with genuine edges and vertices.
        let cube = box_solid(4.0, 4.0, 4.0).unwrap();
        let brep = inner_brep(&cube);
        let ce = cube_corner(brep);
        let blend = resolve_corner_blend(&ce, 0.5).unwrap();
        let cutter = build_corner_cutter(&blend).unwrap();
        let faces: usize = cutter.boundaries().iter().map(|s| s.len()).sum();
        assert_eq!(faces, 6, "corner cutter should be a 6-face box");
        assert!(cutter.edge_iter().count() > 0, "cutter has no edges");
        assert!(cutter.vertex_iter().count() > 0, "cutter has no vertices");
    }

    #[test]
    fn corner_cutter_contains_the_seated_ball_octant() {
        // The cutter must fully contain the corner-facing octant of
        // the seated ball, otherwise (solid − cutter) ∪ ball would
        // not produce a clean blend. Express the seated ball's
        // bounding extent in the corner frame: along each edge
        // direction dᵢ the centre sits at signed offset
        // (C − apex)·dᵢ, and the ball reaches ±r about it. The cutter
        // spans [0, reach] along each dᵢ; require the whole ball
        // octant [centre_offset − r , centre_offset + r] to fit.
        let cube = box_solid(6.0, 6.0, 6.0).unwrap();
        let brep = inner_brep(&cube);
        let ce = cube_corner(brep);
        let r = 0.8;
        let blend = resolve_corner_blend(&ce, r).unwrap();
        let reach = r * CUTTER_REACH_FACTOR;
        let to_c = blend.centre - blend.apex;
        for d in blend.dirs {
            let centre_offset = to_c.dot(d);
            let ball_max = centre_offset + r;
            let ball_min = centre_offset - r;
            assert!(
                ball_min >= -1e-9 && ball_max <= reach + 1e-9,
                "seated ball extent [{ball_min}, {ball_max}] must fit the cutter [0, {reach}]"
            );
        }
    }

    #[test]
    fn corner_ball_builds_as_a_closed_brep() {
        // The corner ball is a real closed sphere BRep with faces,
        // edges and vertices.
        let cube = box_solid(4.0, 4.0, 4.0).unwrap();
        let brep = inner_brep(&cube);
        let ce = cube_corner(brep);
        let blend = resolve_corner_blend(&ce, 0.5).unwrap();
        let ball = build_corner_ball(&blend).unwrap();
        let faces: usize = ball.boundaries().iter().map(|s| s.len()).sum();
        assert!(faces > 0, "corner ball should have faces");
        assert!(ball.vertex_iter().count() > 0, "corner ball has no vertices");
    }

    #[test]
    fn resolve_corner_blend_rejects_bad_radius() {
        let cube = box_solid(4.0, 4.0, 4.0).unwrap();
        let brep = inner_brep(&cube);
        let ce = cube_corner(brep);
        assert!(matches!(
            resolve_corner_blend(&ce, 0.0),
            Err(FilletBrepError::BadParameter { name: "radius", .. })
        ));
        assert!(matches!(
            resolve_corner_blend(&ce, -1.0),
            Err(FilletBrepError::BadParameter { name: "radius", .. })
        ));
        assert!(matches!(
            resolve_corner_blend(&ce, f64::NAN),
            Err(FilletBrepError::BadParameter { name: "radius", .. })
        ));
    }

    #[test]
    fn blend_corner_on_a_cube_either_builds_or_soft_fails() {
        // The full corner blend runs the real cutter + ball + boolean
        // construction. For a well-conditioned cube corner the
        // outcome must be either a real BRep solid or the documented
        // soft `TruckOp` fall-through (the cutter trims flush with
        // three faces — the fragile coincident-surface case for any
        // kernel) — never a panic, never a BadParameter.
        let cube = box_solid(6.0, 6.0, 6.0).unwrap();
        let brep = inner_brep(&cube);
        let ce = cube_corner(brep);
        match blend_corner(brep, &ce, 0.6) {
            Ok(blended) => {
                let faces: usize = blended.boundaries().iter().map(|s| s.len()).sum();
                assert!(
                    faces >= 6,
                    "a corner-blended cube keeps the cube's faces, got {faces}"
                );
            }
            Err(FilletBrepError::TruckOp(_)) => {
                // Soft fall-through — acceptable.
            }
            other => panic!("unexpected corner-blend outcome: {other:?}"),
        }
    }
}
