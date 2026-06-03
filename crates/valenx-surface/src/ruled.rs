//! Ruled surfaces (Phase 19E).
//!
//! A ruled surface is the locus swept by a moving straight line. The
//! three constructors supported in v1 are:
//!
//! - [`between_curves`]: linear interpolation between two arbitrary
//!   NURBS curves `c1` and `c2`. The result is a degree-`(1, d)`
//!   NURBS surface where `d = max(c1.degree, c2.degree)`. The two
//!   curves must share a parameter range `[0, 1]` (we normalise
//!   internally); the input curves are degree-elevated to the same
//!   degree if they differ, then both are used as the v-direction
//!   control rows of the output. Knots in the v direction are
//!   inherited from the elevated curve(s); the u-direction knot
//!   vector is the canonical degree-1 `[0, 0, 1, 1]`.
//!
//! - [`extrude_along_vector`]: special case where `c2 = c1 + vector`
//!   (i.e. translation extrusion). Returns the same shape as
//!   `between_curves` but skips the degree-elevation step.
//!
//! - [`cone_from_apex`]: degenerate ruled surface — every ruling
//!   ends at the apex point. The `apex_row` is a row of copies of
//!   the apex (so the resulting surface evaluated at `u = 1` returns
//!   the apex everywhere in `v`).

use nalgebra::Vector3;

use crate::error::SurfaceError;
use crate::nurbs_curve::NurbsCurve;
use crate::nurbs_surface::NurbsSurface;

/// Build a ruled surface between two NURBS curves.
///
/// Each ruling at parameter `v` is the line segment from `c1(v)` to
/// `c2(v)`. The u parameter parameterises the ruling (`u = 0` →
/// point on `c1`, `u = 1` → point on `c2`).
pub fn between_curves(c1: &NurbsCurve, c2: &NurbsCurve) -> Result<NurbsSurface, SurfaceError> {
    // 1. Degree-match.
    let target_degree = c1.degree.max(c2.degree);
    let elevated_1 = if c1.degree < target_degree {
        c1.elevate_degree(target_degree - c1.degree)?
    } else {
        c1.clone()
    };
    let elevated_2 = if c2.degree < target_degree {
        c2.elevate_degree(target_degree - c2.degree)?
    } else {
        c2.clone()
    };

    // 2. Knot-match in the v direction. v1 requires matching knot
    //    vectors after elevation; if they differ, we conservatively
    //    insert the union of both knot sets into both curves so they
    //    end up with the same knots and consequently the same CP
    //    count.
    let merged_v_knots = merge_knot_vectors(&elevated_1.knots, &elevated_2.knots);
    let row_1 = refine_curve_to_knot_set(&elevated_1, &merged_v_knots)?;
    let row_2 = refine_curve_to_knot_set(&elevated_2, &merged_v_knots)?;

    // 3. Build a 2 x nv control grid: row 0 = row_1, row 1 = row_2.
    let nv = row_1.control_points.len();
    if row_2.control_points.len() != nv {
        return Err(SurfaceError::BadKnotVector {
            reason: "ruled surface: post-refinement CP counts disagree".into(),
        });
    }
    let cps: Vec<Vec<Vector3<f64>>> =
        vec![row_1.control_points.clone(), row_2.control_points.clone()];
    let weights: Vec<Vec<f64>> = vec![row_1.weights.clone(), row_2.weights.clone()];
    let u_knots = vec![0.0, 0.0, 1.0, 1.0];
    NurbsSurface::new(1, target_degree, u_knots, row_1.knots.clone(), cps, weights)
}

/// Extrude `curve` along `vector` to produce a ruled surface.
pub fn extrude_along_vector(
    curve: &NurbsCurve,
    vector: Vector3<f64>,
) -> Result<NurbsSurface, SurfaceError> {
    let c2_cps: Vec<Vector3<f64>> = curve.control_points.iter().map(|p| p + vector).collect();
    let c2 = NurbsCurve::new(
        curve.degree,
        curve.knots.clone(),
        c2_cps,
        curve.weights.clone(),
    )?;
    between_curves(curve, &c2)
}

/// Build a ruled surface that connects `curve` to a single `apex`
/// point. Every ruling terminates at the apex.
pub fn cone_from_apex(
    curve: &NurbsCurve,
    apex: Vector3<f64>,
) -> Result<NurbsSurface, SurfaceError> {
    let n_cps = curve.control_points.len();
    let apex_cps = vec![apex; n_cps];
    let apex_weights = curve.weights.clone();
    let c2 = NurbsCurve::new(curve.degree, curve.knots.clone(), apex_cps, apex_weights)?;
    between_curves(curve, &c2)
}

// ===== Helpers =====

fn merge_knot_vectors(a: &[f64], b: &[f64]) -> Vec<f64> {
    // Union of both knot vectors preserving multiplicities. The
    // result must be non-decreasing.
    let mut out: Vec<f64> = Vec::with_capacity(a.len() + b.len());
    let mut i = 0;
    let mut j = 0;
    while i < a.len() && j < b.len() {
        if (a[i] - b[j]).abs() < 1.0e-12 {
            out.push(a[i]);
            i += 1;
            j += 1;
        } else if a[i] < b[j] {
            out.push(a[i]);
            i += 1;
        } else {
            out.push(b[j]);
            j += 1;
        }
    }
    out.extend_from_slice(&a[i..]);
    out.extend_from_slice(&b[j..]);
    out
}

fn refine_curve_to_knot_set(
    c: &NurbsCurve,
    target_knots: &[f64],
) -> Result<NurbsCurve, SurfaceError> {
    if c.knots == target_knots {
        return Ok(c.clone());
    }
    // Insert any knots present in target_knots beyond what c has.
    let mut current = c.clone();
    // Build a frequency map for both vectors.
    let mut cur_idx = 0;
    let mut tgt_idx = 0;
    let mut to_insert: Vec<f64> = Vec::new();
    while tgt_idx < target_knots.len() {
        if cur_idx < current.knots.len()
            && (current.knots[cur_idx] - target_knots[tgt_idx]).abs() < 1.0e-12
        {
            cur_idx += 1;
            tgt_idx += 1;
        } else if cur_idx < current.knots.len() && current.knots[cur_idx] < target_knots[tgt_idx] {
            cur_idx += 1;
        } else {
            to_insert.push(target_knots[tgt_idx]);
            tgt_idx += 1;
        }
    }
    for u in to_insert {
        if u <= current.parameter_range().0 || u >= current.parameter_range().1 {
            // Skip boundary clamp knots — those are baked in.
            continue;
        }
        current = current.insert_knot(u)?;
    }
    Ok(current)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn line_curve(a: Vector3<f64>, b: Vector3<f64>) -> NurbsCurve {
        let knots = vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
        let cps = vec![a, a + (b - a) / 3.0, a + 2.0 * (b - a) / 3.0, b];
        NurbsCurve::new(3, knots, cps, vec![1.0; 4]).unwrap()
    }

    #[test]
    fn ruled_between_two_parallel_lines_is_planar() {
        let c1 = line_curve(Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 0.0, 0.0));
        let c2 = line_curve(Vector3::new(0.0, 1.0, 0.0), Vector3::new(1.0, 1.0, 0.0));
        let s = between_curves(&c1, &c2).unwrap();
        // The surface should be the unit square in z=0.
        for &(u, v) in &[(0.0_f64, 0.0_f64), (0.5, 0.3), (1.0, 1.0), (0.7, 0.8)] {
            let p = s.evaluate(u, v);
            assert!(p.z.abs() < 1.0e-9, "p.z = {}", p.z);
            // x = v (since v parameterises along the line), y = u
            // (since u interpolates between c1.y=0 and c2.y=1).
            assert!((p.x - v).abs() < 1.0e-9, "p.x = {}, v = {}", p.x, v);
            assert!((p.y - u).abs() < 1.0e-9, "p.y = {}, u = {}", p.y, u);
        }
    }

    #[test]
    fn extrude_along_vector_works() {
        let c = line_curve(Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 0.0, 0.0));
        let s = extrude_along_vector(&c, Vector3::new(0.0, 0.0, 2.0)).unwrap();
        // At u=0 along the curve, z should be 0; at u=1, z = 2.
        let p0 = s.evaluate(0.0, 0.5);
        let p1 = s.evaluate(1.0, 0.5);
        assert!(p0.z.abs() < 1.0e-9);
        assert!((p1.z - 2.0).abs() < 1.0e-9);
    }

    #[test]
    fn cone_from_apex_collapses_at_apex() {
        let c = line_curve(Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 0.0, 0.0));
        let apex = Vector3::new(0.5, 0.5, 1.0);
        let s = cone_from_apex(&c, apex).unwrap();
        // At u=1 (apex side), the surface evaluates to the apex for
        // any v.
        for &v in &[0.0_f64, 0.25, 0.5, 0.75, 1.0] {
            let p = s.evaluate(1.0, v);
            assert!((p - apex).norm() < 1.0e-9, "apex check: {p:?}");
        }
    }
}
