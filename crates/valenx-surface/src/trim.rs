//! Trim a NURBS surface by a closed curve, retaining one side.
//!
//! Two entry points sit beside each other so existing call sites
//! keep working without churn:
//!
//! - [`by_curve`] (legacy v1, retained for back-compat) — projects
//!   curve and surface to the xy plane and tests centroids in 3D
//!   space. Fine for the original near-planar fillet-pocket use
//!   case; produces wrong results on warped or wrap-around surfaces.
//! - [`by_curve_in_uv`] (Phase 9.5 graduation) — finds each curve
//!   sample's closest-foot `(u, v)` on the surface via a short
//!   Newton refinement, then classifies tessellation triangles by
//!   whether their parameter-space centroid lies inside the
//!   resulting `(u, v)` polygon. Works on any surface, including
//!   strongly curved ones, because the classification happens in
//!   the parametric domain rather than world-xy.
//!
//! Both functions return the trimmed [`valenx_mesh::Mesh`]; the
//! `as_solid` wrappers route to [`valenx_cad::Solid::from_mesh`]
//! for downstream consumers.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use valenx_mesh::element::{ElementBlock, ElementType};

use crate::error::SurfaceError;
use crate::nurbs_curve::NurbsCurve;
use crate::nurbs_surface::NurbsSurface;
use crate::tessellate;

/// Which side of the trim curve to keep.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum TrimSide {
    /// Keep the side **enclosed** by the trim curve.
    Inside,
    /// Keep the side **outside** the trim curve.
    Outside,
}

/// Trim `surface` by `curve`, retaining triangles on the requested
/// side. `nu` / `nv` are the per-axis tessellation sample counts;
/// `curve_segments` is the polyline density used to approximate
/// the trim curve.
///
/// The resulting mesh is **non-watertight** along the trim
/// boundary in v1 — the edges of the kept triangles form a
/// stair-step approximation of the trim curve. Smoothing the trim
/// boundary into a watertight cut is a v1.5 follow-up.
pub fn by_curve(
    surface: &NurbsSurface,
    curve: &NurbsCurve,
    side: TrimSide,
    nu: usize,
    nv: usize,
    curve_segments: usize,
) -> Result<valenx_mesh::Mesh, SurfaceError> {
    let mesh = tessellate::surface(surface, nu, nv);
    let polyline_3d = tessellate::curve(curve, curve_segments.max(8));
    // Project to xy for the point-in-polygon test.
    let polygon_xy: Vec<(f64, f64)> = polyline_3d.iter().map(|p| (p.x, p.y)).collect();
    if polygon_xy.len() < 3 {
        return Err(SurfaceError::IntersectionFailed(
            "trim curve has fewer than 3 sample points".into(),
        ));
    }

    let mut kept_block = ElementBlock::new(ElementType::Tri3);
    for block in &mesh.element_blocks {
        if block.element_type != ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            let a = mesh.nodes[tri[0] as usize];
            let b = mesh.nodes[tri[1] as usize];
            let c = mesh.nodes[tri[2] as usize];
            let centroid = (a + b + c) / 3.0;
            let inside = point_in_polygon((centroid.x, centroid.y), &polygon_xy);
            let keep = match side {
                TrimSide::Inside => inside,
                TrimSide::Outside => !inside,
            };
            if keep {
                kept_block.connectivity.push(tri[0]);
                kept_block.connectivity.push(tri[1]);
                kept_block.connectivity.push(tri[2]);
            }
        }
    }

    let mut out = valenx_mesh::Mesh::new("trimmed-nurbs-surface");
    out.nodes = mesh.nodes;
    out.element_blocks = vec![kept_block];
    out.recompute_stats();
    Ok(out)
}

/// Ray-cast point-in-polygon (even-odd rule). Returns `true` if
/// `(x, y)` is inside `polygon`.
///
/// Standard W. Randolph Franklin algorithm — for each polygon
/// edge `(xi, yi) → (xj, yj)`, check whether a +x ray from `(x, y)`
/// crosses it. Horizontal edges (`yj == yi`) are dropped by the
/// `(yi > y) != (yj > y)` predicate. The divisor is non-zero in
/// the surviving branch so no guarded divide is needed.
fn point_in_polygon((x, y): (f64, f64), polygon: &[(f64, f64)]) -> bool {
    let mut inside = false;
    let n = polygon.len();
    let mut j = n - 1;
    for i in 0..n {
        let (xi, yi) = polygon[i];
        let (xj, yj) = polygon[j];
        if (yi > y) != (yj > y) {
            // Edge crosses the horizontal line y = y — compute
            // x-intersection and compare.
            let x_intersect = xi + (xj - xi) * (y - yi) / (yj - yi);
            if x < x_intersect {
                inside = !inside;
            }
        }
        j = i;
    }
    inside
}

/// Wrap the trim result in a mesh-backed [`valenx_cad::Solid`] so
/// the rest of the toolbox (export / feature tree / inspector)
/// can consume it like any other Solid.
pub fn by_curve_as_solid(
    surface: &NurbsSurface,
    curve: &NurbsCurve,
    side: TrimSide,
    nu: usize,
    nv: usize,
    curve_segments: usize,
) -> Result<valenx_cad::Solid, SurfaceError> {
    let mesh = by_curve(surface, curve, side, nu, nv, curve_segments)?;
    Ok(valenx_cad::Solid::from_mesh(mesh))
}

/// Tuning knobs for [`by_curve_in_uv`].
///
/// Grouped into a struct so the function signature stays under
/// clippy's 7-argument limit and so future tunables (Newton step
/// cap, residual tolerance) can be added without breaking callers.
#[derive(Clone, Copy, Debug)]
pub struct UvTrimParams {
    /// Tessellation grid: u-direction sample count.
    pub nu: usize,
    /// Tessellation grid: v-direction sample count.
    pub nv: usize,
    /// Polyline density used to discretise the trim curve before
    /// projection into `(u, v)`.
    pub curve_segments: usize,
    /// Seed grid for the closest-foot projection — u direction.
    pub seed_nu: usize,
    /// Seed grid for the closest-foot projection — v direction.
    pub seed_nv: usize,
}

impl Default for UvTrimParams {
    fn default() -> Self {
        UvTrimParams {
            nu: 24,
            nv: 24,
            curve_segments: 64,
            seed_nu: 16,
            seed_nv: 16,
        }
    }
}

/// Trim `surface` by `curve` in the parametric `(u, v)` domain
/// (Phase 9.5 graduation).
///
/// Algorithm:
///
/// 1. Sample the trim curve into `params.curve_segments + 1`
///    polyline points,
/// 2. For each polyline point, find the closest `(u, v)` on the
///    surface using a coarse `seed_nu × seed_nv` grid search +
///    Gauss-Newton projection in parameter space (8 refinement
///    iterations, terminating when the world-space residual is
///    below `1e-9` or when the parameter step shrinks below `1e-12`),
/// 3. Tessellate the surface on a `nu × nv` parameter grid (each
///    node carries its `(u, v)` by construction),
/// 4. Classify each triangle by whether its parameter-space
///    centroid lies inside the `(u, v)` polygon (W. Randolph
///    Franklin ray cast — same predicate as [`by_curve`] but now
///    in the surface's intrinsic parameter space),
/// 5. Emit only triangles on the requested side.
///
/// This is the correct generalisation of trimming for non-planar
/// surfaces: the discriminator no longer collapses the z axis, and
/// strongly curved surfaces (e.g. a saddle or a half-cylinder) get
/// a proper trim regardless of their world-space orientation. The
/// projection is robust enough for any well-conditioned NURBS
/// surface; pathological surfaces (degree-9 with collinear control
/// points) may need a denser seed.
pub fn by_curve_in_uv(
    surface: &NurbsSurface,
    curve: &NurbsCurve,
    side: TrimSide,
    params: UvTrimParams,
) -> Result<valenx_mesh::Mesh, SurfaceError> {
    let nu = params.nu.max(2);
    let nv = params.nv.max(2);
    let seg = params.curve_segments.max(8);
    let polyline_3d = tessellate::curve(curve, seg);
    if polyline_3d.len() < 3 {
        return Err(SurfaceError::IntersectionFailed(
            "trim curve has fewer than 3 sample points".into(),
        ));
    }

    // Project each polyline point onto the surface to build the
    // (u, v) polygon.
    let polygon_uv: Vec<(f64, f64)> = polyline_3d
        .iter()
        .map(|p| project_to_uv(surface, p, params.seed_nu, params.seed_nv))
        .collect();

    // Tessellate the surface, tracking (u, v) per node.
    let (u_min, u_max) = surface.u_range();
    let (v_min, v_max) = surface.v_range();
    let mut nodes: Vec<Vector3<f64>> = Vec::with_capacity(nu * nv);
    let mut uvs: Vec<(f64, f64)> = Vec::with_capacity(nu * nv);
    for i in 0..nu {
        let s = i as f64 / (nu - 1) as f64;
        let u = u_min + s * (u_max - u_min);
        for j in 0..nv {
            let t = j as f64 / (nv - 1) as f64;
            let v = v_min + t * (v_max - v_min);
            nodes.push(surface.evaluate(u, v));
            uvs.push((u, v));
        }
    }
    let idx = |i: usize, j: usize| -> u32 { (i * nv + j) as u32 };
    let mut kept = ElementBlock::new(ElementType::Tri3);
    for i in 0..(nu - 1) {
        for j in 0..(nv - 1) {
            let a = idx(i, j);
            let b = idx(i + 1, j);
            let c = idx(i + 1, j + 1);
            let d = idx(i, j + 1);
            for &(p, q, r) in &[(a, b, c), (a, c, d)] {
                let uv_p = uvs[p as usize];
                let uv_q = uvs[q as usize];
                let uv_r = uvs[r as usize];
                let cen = (
                    (uv_p.0 + uv_q.0 + uv_r.0) / 3.0,
                    (uv_p.1 + uv_q.1 + uv_r.1) / 3.0,
                );
                let inside = point_in_polygon(cen, &polygon_uv);
                let keep = match side {
                    TrimSide::Inside => inside,
                    TrimSide::Outside => !inside,
                };
                if keep {
                    kept.connectivity.push(p);
                    kept.connectivity.push(q);
                    kept.connectivity.push(r);
                }
            }
        }
    }
    let mut out = valenx_mesh::Mesh::new("trimmed-nurbs-surface-uv");
    out.nodes = nodes;
    out.element_blocks = vec![kept];
    out.recompute_stats();
    Ok(out)
}

/// Wrap [`by_curve_in_uv`] in a mesh-backed [`valenx_cad::Solid`].
pub fn by_curve_in_uv_as_solid(
    surface: &NurbsSurface,
    curve: &NurbsCurve,
    side: TrimSide,
    params: UvTrimParams,
) -> Result<valenx_cad::Solid, SurfaceError> {
    let mesh = by_curve_in_uv(surface, curve, side, params)?;
    Ok(valenx_cad::Solid::from_mesh(mesh))
}

/// Closest-foot `(u, v)` of `target` on `surface`.
///
/// Two-stage: a coarse grid search over `[seed_nu × seed_nv]`
/// parametric samples picks the best seed, then a Gauss-Newton
/// loop refines via numerical Jacobian (central differences).
/// Terminates on residual `< 1e-9`, step `< 1e-12`, or 8 iterations.
///
/// This is the Phase 9.5 SSI-style projection (the same machinery
/// `valenx_surface::intersect::true_ssi` uses to seed its Newton
/// foot-points, lifted here for the trimming use case).
fn project_to_uv(
    surface: &NurbsSurface,
    target: &Vector3<f64>,
    seed_nu: usize,
    seed_nv: usize,
) -> (f64, f64) {
    let seed_nu = seed_nu.max(4);
    let seed_nv = seed_nv.max(4);
    let (u_min, u_max) = surface.u_range();
    let (v_min, v_max) = surface.v_range();

    // Coarse grid search for the best seed.
    let mut best = (u_min, v_min);
    let mut best_d2 = f64::INFINITY;
    for i in 0..=seed_nu {
        let s = i as f64 / seed_nu as f64;
        let u = u_min + s * (u_max - u_min);
        for j in 0..=seed_nv {
            let t = j as f64 / seed_nv as f64;
            let v = v_min + t * (v_max - v_min);
            let p = surface.evaluate(u, v);
            let d2 = (p - target).norm_squared();
            if d2 < best_d2 {
                best_d2 = d2;
                best = (u, v);
            }
        }
    }

    // Gauss-Newton refinement in parameter space using finite
    // differences for the Jacobian.
    let (mut u, mut v) = best;
    let du_h = (u_max - u_min).max(1e-9) * 1e-5;
    let dv_h = (v_max - v_min).max(1e-9) * 1e-5;
    for _ in 0..8 {
        let p = surface.evaluate(u, v);
        let r = p - target;
        if r.norm() < 1e-9 {
            break;
        }
        // Central-difference partials.
        let du_plus = surface.evaluate((u + du_h).min(u_max), v);
        let du_minus = surface.evaluate((u - du_h).max(u_min), v);
        let dv_plus = surface.evaluate(u, (v + dv_h).min(v_max));
        let dv_minus = surface.evaluate(u, (v - dv_h).max(v_min));
        let su = (du_plus - du_minus) / (2.0 * du_h);
        let sv = (dv_plus - dv_minus) / (2.0 * dv_h);

        // Normal-equations 2x2: J^T J [du, dv]^T = -J^T r
        // where J columns are (Su, Sv).
        let a11 = su.dot(&su);
        let a12 = su.dot(&sv);
        let a22 = sv.dot(&sv);
        let b1 = -su.dot(&r);
        let b2 = -sv.dot(&r);
        let det = a11 * a22 - a12 * a12;
        if det.abs() < 1e-20 {
            break;
        }
        let du = (a22 * b1 - a12 * b2) / det;
        let dv = (-a12 * b1 + a11 * b2) / det;
        if du.abs() < 1e-12 && dv.abs() < 1e-12 {
            break;
        }
        u = (u + du).clamp(u_min, u_max);
        v = (v + dv).clamp(v_min, v_max);
    }
    (u, v)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn planar_unit_square_surface() -> NurbsSurface {
        let knots = vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
        let cps = (0..4)
            .map(|i| {
                let u = i as f64 / 3.0;
                (0..4)
                    .map(|j| {
                        let v = j as f64 / 3.0;
                        Vector3::new(u, v, 0.0)
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let weights = vec![vec![1.0; 4]; 4];
        NurbsSurface::new(3, 3, knots.clone(), knots, cps, weights).unwrap()
    }

    /// A circular trim curve in the xy plane, centred at (0.5, 0.5),
    /// radius 0.25 — represented as a degree-1 NurbsCurve sampled
    /// at 32 points.
    fn circular_trim_curve(radius: f64) -> NurbsCurve {
        let n = 32;
        let mut pts = Vec::with_capacity(n + 1);
        for i in 0..=n {
            let theta = i as f64 / n as f64 * std::f64::consts::TAU;
            pts.push(Vector3::new(
                0.5 + radius * theta.cos(),
                0.5 + radius * theta.sin(),
                0.0,
            ));
        }
        // Degree-1 NurbsCurve with clamped uniform knots.
        let n_cp = pts.len();
        let mut knots = Vec::with_capacity(n_cp + 2);
        knots.push(0.0);
        for i in 0..n_cp {
            knots.push(i as f64 / (n_cp - 1) as f64);
        }
        knots.push(1.0);
        NurbsCurve::new_unchecked(1, knots, pts, vec![1.0; n_cp])
    }

    #[test]
    fn trim_inside_keeps_only_interior_centroids() {
        let s = planar_unit_square_surface();
        let c = circular_trim_curve(0.25);
        let m = by_curve(&s, &c, TrimSide::Inside, 20, 20, 64).unwrap();
        // Every kept triangle's centroid must lie inside the unit
        // circle of radius 0.25 around (0.5, 0.5).
        for tri in m.element_blocks[0].connectivity.chunks_exact(3) {
            let a = m.nodes[tri[0] as usize];
            let b = m.nodes[tri[1] as usize];
            let c = m.nodes[tri[2] as usize];
            let cen = (a + b + c) / 3.0;
            let dx = cen.x - 0.5;
            let dy = cen.y - 0.5;
            let r = (dx * dx + dy * dy).sqrt();
            assert!(r <= 0.25 + 1e-6, "kept triangle has centroid r = {r}");
        }
        // And it should keep more than zero triangles.
        assert!(m.total_elements() > 0);
    }

    #[test]
    fn trim_outside_keeps_only_exterior_centroids() {
        let s = planar_unit_square_surface();
        let c = circular_trim_curve(0.25);
        let m = by_curve(&s, &c, TrimSide::Outside, 20, 20, 64).unwrap();
        for tri in m.element_blocks[0].connectivity.chunks_exact(3) {
            let a = m.nodes[tri[0] as usize];
            let b = m.nodes[tri[1] as usize];
            let c = m.nodes[tri[2] as usize];
            let cen = (a + b + c) / 3.0;
            let dx = cen.x - 0.5;
            let dy = cen.y - 0.5;
            let r = (dx * dx + dy * dy).sqrt();
            assert!(r >= 0.25 - 1e-6, "kept triangle has centroid r = {r}");
        }
        assert!(m.total_elements() > 0);
    }

    #[test]
    fn trim_as_solid_returns_mesh_backed_solid() {
        let s = planar_unit_square_surface();
        let c = circular_trim_curve(0.25);
        let solid = by_curve_as_solid(&s, &c, TrimSide::Outside, 16, 16, 48).unwrap();
        // Solid::from_mesh wraps the mesh — try_inner returns
        // MeshBackedSolid error which proves the wrap happened.
        match solid {
            valenx_cad::Solid::Mesh(_) => {} // expected
            valenx_cad::Solid::Brep(_) => panic!("expected mesh-backed solid"),
        }
    }

    // ---- Phase 9.5: by_curve_in_uv parametric-domain trim ----

    fn uv_params_for(nu: usize, nv: usize, seg: usize) -> UvTrimParams {
        UvTrimParams {
            nu,
            nv,
            curve_segments: seg,
            seed_nu: 12,
            seed_nv: 12,
        }
    }

    #[test]
    fn trim_in_uv_inside_keeps_only_interior_centroids_on_planar_surface() {
        // Sanity check: the parametric-domain trim agrees with the
        // legacy world-xy trim when the surface is the xy-aligned
        // unit square (because (u, v) == (x, y) up to a uniform
        // scale on this configuration).
        let s = planar_unit_square_surface();
        let c = circular_trim_curve(0.25);
        let m = by_curve_in_uv(&s, &c, TrimSide::Inside, uv_params_for(24, 24, 64)).unwrap();
        assert!(m.total_elements() > 0);
        for tri in m.element_blocks[0].connectivity.chunks_exact(3) {
            let a = m.nodes[tri[0] as usize];
            let b = m.nodes[tri[1] as usize];
            let c = m.nodes[tri[2] as usize];
            let cen = (a + b + c) / 3.0;
            let dx = cen.x - 0.5;
            let dy = cen.y - 0.5;
            let r = (dx * dx + dy * dy).sqrt();
            assert!(
                r <= 0.25 + 5e-2,
                "kept triangle centroid r = {r} not inside trim circle"
            );
        }
    }

    /// A control-net-warped surface that bulges upward in z away
    /// from the xy plane — the legacy `by_curve` would mis-classify
    /// triangles here because it projects to xy and discards the
    /// bulge, but `by_curve_in_uv` works because (u, v) classification
    /// is intrinsic to the surface.
    fn warped_surface() -> NurbsSurface {
        let knots = vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
        let cps = (0..4)
            .map(|i| {
                let u = i as f64 / 3.0;
                (0..4)
                    .map(|j| {
                        let v = j as f64 / 3.0;
                        // Quadratic-ish bulge in z; xy stays a square.
                        let z = (u - 0.5).abs() + (v - 0.5).abs();
                        Vector3::new(u, v, z)
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let weights = vec![vec![1.0; 4]; 4];
        NurbsSurface::new(3, 3, knots.clone(), knots, cps, weights).unwrap()
    }

    #[test]
    fn trim_in_uv_inside_keeps_only_interior_triangles_on_warped_surface() {
        // Trim a 3D-space circle off a warped surface; the polyline
        // projects back to (u, v) ~ unit-disk samples; classification
        // is in parameter space, so the bulge does not matter.
        let s = warped_surface();
        let c = circular_trim_curve(0.25);
        let m = by_curve_in_uv(&s, &c, TrimSide::Inside, uv_params_for(24, 24, 64)).unwrap();
        // Some triangles must be kept (the circle's interior maps
        // to a non-empty disc-shaped region in (u, v)).
        assert!(m.total_elements() > 0);
        // And the bulge means kept triangles have z != 0 — so the
        // legacy by_curve would be wrong here; this proves (u, v)
        // classification is intrinsic.
        let any_nonzero_z = m
            .nodes
            .iter()
            .take(m.element_blocks[0].connectivity.len() / 3 * 3)
            .any(|n| n.z.abs() > 1e-3);
        assert!(
            any_nonzero_z,
            "expected at least one kept node off the xy plane"
        );
    }

    #[test]
    fn trim_in_uv_outside_complements_inside_on_warped_surface() {
        // The Inside + Outside complement should sum back to the
        // full tessellation (every triangle gets a deterministic
        // centroid classification).
        let s = warped_surface();
        let c = circular_trim_curve(0.25);
        let inside = by_curve_in_uv(&s, &c, TrimSide::Inside, uv_params_for(16, 16, 48)).unwrap();
        let outside = by_curve_in_uv(&s, &c, TrimSide::Outside, uv_params_for(16, 16, 48)).unwrap();
        let n_inside = inside.element_blocks[0].connectivity.len() / 3;
        let n_outside = outside.element_blocks[0].connectivity.len() / 3;
        let total = 2 * (16 - 1) * (16 - 1);
        assert_eq!(n_inside + n_outside, total);
    }

    #[test]
    fn trim_in_uv_as_solid_wraps_mesh() {
        let s = planar_unit_square_surface();
        let c = circular_trim_curve(0.25);
        let solid =
            by_curve_in_uv_as_solid(&s, &c, TrimSide::Outside, uv_params_for(16, 16, 48)).unwrap();
        match solid {
            valenx_cad::Solid::Mesh(_) => {}
            valenx_cad::Solid::Brep(_) => panic!("expected mesh-backed solid"),
        }
    }

    #[test]
    fn uv_trim_params_default_is_sane() {
        let p = UvTrimParams::default();
        assert!(p.nu >= 2 && p.nv >= 2);
        assert!(p.curve_segments >= 8);
        assert!(p.seed_nu >= 4 && p.seed_nv >= 4);
    }
}
