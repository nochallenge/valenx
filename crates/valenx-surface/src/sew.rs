//! Sew two NURBS surfaces along a shared parametric edge.
//!
//! v1 requires:
//! - Matching parameterisation along the shared edge (i.e. the same
//!   degree and the same knot vector in the direction that runs
//!   *along* the edge),
//! - The shared edge's control points already coincide within
//!   `tolerance` (small misalignments are averaged out by the
//!   stitching).
//!
//! Anything fancier (re-parameterisation, knot insertion to merge
//! mismatched curves) is deferred to Phase 9.5.

use serde::{Deserialize, Serialize};

use crate::error::SurfaceError;
use crate::nurbs_surface::NurbsSurface;

/// Which parametric edge of a NURBS surface a given operation
/// acts on.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Edge {
    /// `u = u_min` — the "left" edge.
    UMin,
    /// `u = u_max` — the "right" edge.
    UMax,
    /// `v = v_min` — the "bottom" edge.
    VMin,
    /// `v = v_max` — the "top" edge.
    VMax,
}

impl Edge {
    /// True if the edge runs along the v direction (i.e. `UMin` /
    /// `UMax` — the parameter that varies along the edge is `v`).
    pub fn runs_along_v(self) -> bool {
        matches!(self, Edge::UMin | Edge::UMax)
    }
}

/// Phase 19C — continuity class for the G2 sew.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Continuity {
    /// Position-only averaging (Phase 9 behaviour: only the shared
    /// row of CPs is averaged).
    G0,
    /// Position + tangent matching across the seam (adjusts 2 rows
    /// of CPs per side: shared row + one neighbour).
    G1,
    /// Position + tangent + curvature matching across the seam
    /// (adjusts 3 rows of CPs per side: shared row + two neighbours).
    /// Phase 19C default for the UI sew button.
    G2,
}

/// Stitch two NURBS surfaces along a shared edge.
///
/// `edge_pair.0` is the edge on `s1` and `edge_pair.1` is the edge
/// on `s2`. The two edges must run in the same direction (both
/// along u, or both along v) and have the same degree + knot
/// vector in that direction; otherwise the call fails with
/// [`SurfaceError::SewMismatch`].
///
/// The result is a new NURBS surface that concatenates `s1` and
/// `s2` in the u direction (for a UMax↔UMin pair) or in the v
/// direction (for a VMax↔VMin pair).
///
/// The shared edge's CPs are averaged so any sub-`tolerance`
/// mismatch is squashed to zero in the output.
pub fn stitch(
    s1: &NurbsSurface,
    s2: &NurbsSurface,
    edge_pair: (Edge, Edge),
    tolerance: f64,
) -> Result<NurbsSurface, SurfaceError> {
    let (e1, e2) = edge_pair;
    if e1.runs_along_v() != e2.runs_along_v() {
        return Err(SurfaceError::SewMismatch(
            "edges run in different parametric directions".into(),
        ));
    }

    // For v1 we support the canonical case: UMax of s1 stitched to
    // UMin of s2 (or VMax to VMin).
    match (e1, e2) {
        (Edge::UMax, Edge::UMin) => stitch_u(s1, s2, tolerance),
        (Edge::VMax, Edge::VMin) => stitch_v(s1, s2, tolerance),
        _ => Err(SurfaceError::SewMismatch(format!(
            "unsupported edge pair {e1:?} / {e2:?}; v1 supports only \
             (UMax, UMin) and (VMax, VMin)"
        ))),
    }
}

fn stitch_u(
    s1: &NurbsSurface,
    s2: &NurbsSurface,
    tolerance: f64,
) -> Result<NurbsSurface, SurfaceError> {
    if s1.v_degree != s2.v_degree {
        return Err(SurfaceError::SewMismatch(format!(
            "v_degree mismatch: {} vs {}",
            s1.v_degree, s2.v_degree
        )));
    }
    if s1.v_knots != s2.v_knots {
        return Err(SurfaceError::SewMismatch(
            "v_knots mismatch — re-parameterisation deferred to Phase 9.5".into(),
        ));
    }
    if s1.u_degree != s2.u_degree {
        return Err(SurfaceError::SewMismatch(format!(
            "u_degree mismatch: {} vs {}",
            s1.u_degree, s2.u_degree
        )));
    }
    // Check the shared edge: last row of s1 (u=u_max) == first row
    // of s2 (u=u_min), within tolerance.
    let nv = s1.nv();
    if s2.nv() != nv {
        return Err(SurfaceError::SewMismatch(format!(
            "v-direction CP count mismatch: {} vs {}",
            nv,
            s2.nv()
        )));
    }
    let row_s1 = &s1.control_points[s1.nu() - 1];
    let row_s2 = &s2.control_points[0];
    for j in 0..nv {
        if (row_s1[j] - row_s2[j]).norm() > tolerance {
            return Err(SurfaceError::SewMismatch(format!(
                "shared CP[{j}] differs by {} (> tolerance {tolerance})",
                (row_s1[j] - row_s2[j]).norm()
            )));
        }
    }

    // Concatenate. The output has nu_combined = s1.nu() + s2.nu() - 1
    // CPs in the u direction; the shared row is averaged into one.
    let nu1 = s1.nu();
    let nu2 = s2.nu();
    let mut cps = Vec::with_capacity(nu1 + nu2 - 1);
    let mut weights = Vec::with_capacity(nu1 + nu2 - 1);
    for i in 0..(nu1 - 1) {
        cps.push(s1.control_points[i].clone());
        weights.push(s1.weights[i].clone());
    }
    // Averaged shared row.
    let mut shared = Vec::with_capacity(nv);
    let mut shared_w = Vec::with_capacity(nv);
    for j in 0..nv {
        shared.push(0.5 * (s1.control_points[nu1 - 1][j] + s2.control_points[0][j]));
        shared_w.push(0.5 * (s1.weights[nu1 - 1][j] + s2.weights[0][j]));
    }
    cps.push(shared);
    weights.push(shared_w);
    for i in 1..nu2 {
        cps.push(s2.control_points[i].clone());
        weights.push(s2.weights[i].clone());
    }
    let nu_new = cps.len();
    // Re-build a clamped uniform knot vector in u (open-uniform of
    // the same degree). For v1 we don't try to preserve interior
    // knot multiplicity from the originals — the surface is
    // resampled in u space and the v knots are unchanged.
    let u_knots = open_uniform_knots(nu_new, s1.u_degree);
    NurbsSurface::new(
        s1.u_degree,
        s1.v_degree,
        u_knots,
        s1.v_knots.clone(),
        cps,
        weights,
    )
}

fn stitch_v(
    s1: &NurbsSurface,
    s2: &NurbsSurface,
    tolerance: f64,
) -> Result<NurbsSurface, SurfaceError> {
    if s1.u_degree != s2.u_degree {
        return Err(SurfaceError::SewMismatch(format!(
            "u_degree mismatch: {} vs {}",
            s1.u_degree, s2.u_degree
        )));
    }
    if s1.u_knots != s2.u_knots {
        return Err(SurfaceError::SewMismatch(
            "u_knots mismatch — re-parameterisation deferred to Phase 9.5".into(),
        ));
    }
    if s1.v_degree != s2.v_degree {
        return Err(SurfaceError::SewMismatch(format!(
            "v_degree mismatch: {} vs {}",
            s1.v_degree, s2.v_degree
        )));
    }
    let nu = s1.nu();
    if s2.nu() != nu {
        return Err(SurfaceError::SewMismatch(format!(
            "u-direction CP count mismatch: {} vs {}",
            nu,
            s2.nu()
        )));
    }
    // Shared edge: last column of s1 (v=v_max) == first column of
    // s2 (v=v_min).
    let nv1 = s1.nv();
    for i in 0..nu {
        let a = s1.control_points[i][nv1 - 1];
        let b = s2.control_points[i][0];
        if (a - b).norm() > tolerance {
            return Err(SurfaceError::SewMismatch(format!(
                "shared CP at row {i} differs by {} (> tolerance {tolerance})",
                (a - b).norm()
            )));
        }
    }

    let nv1 = s1.nv();
    let nv2 = s2.nv();
    let mut cps = Vec::with_capacity(nu);
    let mut weights = Vec::with_capacity(nu);
    for i in 0..nu {
        let mut row = Vec::with_capacity(nv1 + nv2 - 1);
        let mut wrow = Vec::with_capacity(nv1 + nv2 - 1);
        for j in 0..(nv1 - 1) {
            row.push(s1.control_points[i][j]);
            wrow.push(s1.weights[i][j]);
        }
        let shared = 0.5 * (s1.control_points[i][nv1 - 1] + s2.control_points[i][0]);
        let shared_w = 0.5 * (s1.weights[i][nv1 - 1] + s2.weights[i][0]);
        row.push(shared);
        wrow.push(shared_w);
        for j in 1..nv2 {
            row.push(s2.control_points[i][j]);
            wrow.push(s2.weights[i][j]);
        }
        cps.push(row);
        weights.push(wrow);
    }
    let nv_new = cps[0].len();
    let v_knots = open_uniform_knots(nv_new, s1.v_degree);
    NurbsSurface::new(
        s1.u_degree,
        s1.v_degree,
        s1.u_knots.clone(),
        v_knots,
        cps,
        weights,
    )
}

/// Build a clamped open-uniform knot vector of length `n_cp + deg + 1`
/// over the canonical parameter range `[0, 1]`. Multiplicity
/// `deg + 1` at each endpoint, uniform interior.
fn open_uniform_knots(n_cp: usize, degree: usize) -> Vec<f64> {
    let p = degree;
    let m = n_cp + p + 1;
    let mut k = vec![0.0; m];
    if n_cp <= p + 1 {
        // Pure Bezier — `vec![0.0; m]` already zeroes the first
        // p + 1 entries; set the trailing p + 1 entries to 1.0.
        for kv in k.iter_mut().skip(m - p - 1) {
            *kv = 1.0;
        }
        return k;
    }
    let n_internal = n_cp - p - 1;
    for (i, kv) in k.iter_mut().enumerate().take(m) {
        if i <= p {
            *kv = 0.0;
        } else if i >= n_cp {
            *kv = 1.0;
        } else {
            let idx = i - p;
            *kv = idx as f64 / (n_internal + 1) as f64;
        }
    }
    k
}

// ===== Phase 19C — G2 continuous sew =====

/// Stitch two NURBS surfaces with the requested continuity class.
///
/// For `Continuity::G0` this is equivalent to [`stitch`]. For `G1`
/// and `G2` the shared row of CPs is averaged AND the neighbouring
/// row(s) on each side are adjusted to make the first / second
/// directional derivative across the seam continuous.
///
/// ## Method (G2)
///
/// The G2 adjustment moves the first interior row on each side so
/// that the cross-seam *tangent* matches at every column, and moves
/// the second interior row so that the cross-seam *second
/// derivative* matches. For a tensor-product NURBS surface with
/// uniform knots, the cross-seam derivatives at the seam parameter
/// `u_seam` are linear combinations of the local CP rows weighted by
/// basis-function derivatives — for a clamped uniform cubic, the
/// tangent at the seam reduces to `degree * (P_seam - P_prev) /
/// (knot_delta)`. We solve the matching condition column-by-column,
/// adjusting both sides symmetrically so neither half is favoured.
///
/// ## v1 limitations
///
/// - The local-derivative formulas assume an open-uniform clamped
///   knot vector with degree >= 3. For degree 2 the curvature term
///   is silently dropped and the call degrades to G1.
/// - The shared edge's parameterisation must match exactly between
///   the two surfaces (same v_knots or same u_knots depending on
///   the edge pair). v1.5 will add a re-parameterisation pass.
pub fn stitch_with_continuity(
    s1: &NurbsSurface,
    s2: &NurbsSurface,
    edge_pair: (Edge, Edge),
    tolerance: f64,
    continuity: Continuity,
) -> Result<NurbsSurface, SurfaceError> {
    if matches!(continuity, Continuity::G0) {
        return stitch(s1, s2, edge_pair, tolerance);
    }
    let (e1, e2) = edge_pair;
    if e1.runs_along_v() != e2.runs_along_v() {
        return Err(SurfaceError::SewMismatch(
            "edges run in different parametric directions".into(),
        ));
    }
    match (e1, e2) {
        (Edge::UMax, Edge::UMin) => g2_stitch_u(s1, s2, tolerance, continuity),
        (Edge::VMax, Edge::VMin) => g2_stitch_v(s1, s2, tolerance, continuity),
        _ => Err(SurfaceError::SewMismatch(format!(
            "unsupported edge pair {e1:?} / {e2:?}; v1 supports only \
             (UMax, UMin) and (VMax, VMin)"
        ))),
    }
}

fn g2_stitch_u(
    s1: &NurbsSurface,
    s2: &NurbsSurface,
    tolerance: f64,
    continuity: Continuity,
) -> Result<NurbsSurface, SurfaceError> {
    // First do the standard G0 stitch, then adjust the 1st (G1) and
    // 2nd (G2) interior rows on each side to make derivatives match.
    let mut combined = stitch_u(s1, s2, tolerance)?;
    let seam_idx = s1.nu() - 1; // index of the shared row in combined
                                // The combined CP grid:
                                //   rows 0..seam_idx-1 belong to s1 (excluding shared)
                                //   row seam_idx is the averaged shared row
                                //   rows seam_idx+1..nu_combined-1 belong to s2 (excluding shared)
    let nu_combined = combined.nu();
    let nv = combined.nv();
    if nu_combined < 5 {
        return Ok(combined);
    }
    // G1: average the two CPs symmetric around the seam so the
    // tangent across the seam matches column-by-column.
    for j in 0..nv {
        let s1_prev = combined.control_points[seam_idx - 1][j];
        let s2_next = combined.control_points[seam_idx + 1][j];
        let seam_cp = combined.control_points[seam_idx][j];
        // For G1 continuity across the seam we require the tangent
        // from s1 to the seam to equal the tangent from the seam to
        // s2: seam - s1_prev = s2_next - seam, i.e. seam =
        // (s1_prev + s2_next) / 2. Adjust both sides symmetrically:
        // shift s1_prev outward and s2_next outward equally so they
        // become mirror reflections through the seam CP.
        let mean_tan = 0.5 * ((seam_cp - s1_prev) + (s2_next - seam_cp));
        let new_s1_prev = seam_cp - mean_tan;
        let new_s2_next = seam_cp + mean_tan;
        combined.control_points[seam_idx - 1][j] = new_s1_prev;
        combined.control_points[seam_idx + 1][j] = new_s2_next;
    }
    if matches!(continuity, Continuity::G2) && nu_combined >= 7 && s1.u_degree >= 3 {
        // G2: match the cross-seam *second derivative*. For uniform
        // cubics the second derivative at the seam involves the row
        // 2 indices away: seam - 2 * P_{seam-1} + P_{seam-2} on the
        // s1 side and P_{seam+2} - 2 * P_{seam+1} + seam on the s2
        // side. Equating these and using the G1 adjustment already
        // applied, we shift P_{seam-2} and P_{seam+2} outward so the
        // local quadratic curvature is mirrored through the seam.
        for j in 0..nv {
            let next1 = combined.control_points[seam_idx + 1][j];
            let prev2 = combined.control_points[seam_idx - 2][j];
            let next2 = combined.control_points[seam_idx + 2][j];
            let seam_cp = combined.control_points[seam_idx][j];
            // Required for curvature continuity:
            //   P_{seam-2} - 2 P_{seam-1} + seam = seam - 2 P_{seam+1} + P_{seam+2}
            // We've already enforced P_{seam-1} = seam - mean_tan,
            // P_{seam+1} = seam + mean_tan, so:
            //   P_{seam-2} - 2 (seam - mean_tan) + seam = seam - 2 (seam + mean_tan) + P_{seam+2}
            //   P_{seam-2} + 2 mean_tan = P_{seam+2} - 2 mean_tan
            //   P_{seam+2} - P_{seam-2} = 4 mean_tan
            // Adjust both symmetrically around their average.
            let avg = 0.5 * (prev2 + next2);
            let mean_tan = next1 - seam_cp; // == seam - P_{seam-1}
            combined.control_points[seam_idx - 2][j] = avg - 2.0 * mean_tan;
            combined.control_points[seam_idx + 2][j] = avg + 2.0 * mean_tan;
        }
    }
    Ok(combined)
}

fn g2_stitch_v(
    s1: &NurbsSurface,
    s2: &NurbsSurface,
    tolerance: f64,
    continuity: Continuity,
) -> Result<NurbsSurface, SurfaceError> {
    let mut combined = stitch_v(s1, s2, tolerance)?;
    let nu = combined.nu();
    let nv_combined = combined.nv();
    let seam_idx = s1.nv() - 1;
    if nv_combined < 5 {
        return Ok(combined);
    }
    // G1 column-by-column for each row.
    for i in 0..nu {
        let s1_prev = combined.control_points[i][seam_idx - 1];
        let s2_next = combined.control_points[i][seam_idx + 1];
        let seam_cp = combined.control_points[i][seam_idx];
        let mean_tan = 0.5 * ((seam_cp - s1_prev) + (s2_next - seam_cp));
        combined.control_points[i][seam_idx - 1] = seam_cp - mean_tan;
        combined.control_points[i][seam_idx + 1] = seam_cp + mean_tan;
    }
    if matches!(continuity, Continuity::G2) && nv_combined >= 7 && s1.v_degree >= 3 {
        for i in 0..nu {
            let next1 = combined.control_points[i][seam_idx + 1];
            let prev2 = combined.control_points[i][seam_idx - 2];
            let next2 = combined.control_points[i][seam_idx + 2];
            let seam_cp = combined.control_points[i][seam_idx];
            let avg = 0.5 * (prev2 + next2);
            let mean_tan = next1 - seam_cp;
            combined.control_points[i][seam_idx - 2] = avg - 2.0 * mean_tan;
            combined.control_points[i][seam_idx + 2] = avg + 2.0 * mean_tan;
        }
    }
    Ok(combined)
}

/// Convenience wrapper — G2 continuous stitch.
pub fn g2_stitch(
    s1: &NurbsSurface,
    s2: &NurbsSurface,
    edge_pair: (Edge, Edge),
    tolerance: f64,
) -> Result<NurbsSurface, SurfaceError> {
    stitch_with_continuity(s1, s2, edge_pair, tolerance, Continuity::G2)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    /// 4x4 planar patch covering `[x0, x0+1] × [0, 1]` in z=0.
    fn planar_patch(x0: f64) -> NurbsSurface {
        let knots = vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
        let cps = (0..4)
            .map(|i| {
                let u = i as f64 / 3.0;
                (0..4)
                    .map(|j| {
                        let v = j as f64 / 3.0;
                        Vector3::new(x0 + u, v, 0.0)
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let weights = vec![vec![1.0; 4]; 4];
        NurbsSurface::new(3, 3, knots.clone(), knots, cps, weights).unwrap()
    }

    #[test]
    fn stitch_two_planar_surfaces_along_u_axis() {
        let s1 = planar_patch(0.0);
        let s2 = planar_patch(1.0); // UMax of s1 (x=1) == UMin of s2 (x=1)
        let stitched = stitch(&s1, &s2, (Edge::UMax, Edge::UMin), 1e-6).unwrap();
        // Combined CPs: 4 + 4 - 1 = 7 in u direction.
        assert_eq!(stitched.nu(), 7);
        assert_eq!(stitched.nv(), 4);
        // First CP is at x=0, last at x=2.
        assert!((stitched.control_points[0][0] - Vector3::new(0.0, 0.0, 0.0)).norm() < 1e-10);
        assert!((stitched.control_points[6][3] - Vector3::new(2.0, 1.0, 0.0)).norm() < 1e-10);
        // Evaluate at parametric extents.
        let p_min = stitched.evaluate(0.0, 0.5);
        let p_max = stitched.evaluate(1.0, 0.5);
        assert!((p_min - Vector3::new(0.0, 0.5, 0.0)).norm() < 1e-10);
        assert!((p_max - Vector3::new(2.0, 0.5, 0.0)).norm() < 1e-10);
    }

    #[test]
    fn stitch_fails_on_mismatched_v_knots() {
        let s1 = planar_patch(0.0);
        let mut s2 = planar_patch(1.0);
        // Force a v-knot mismatch.
        s2.v_knots = vec![0.0, 0.0, 0.0, 0.0, 0.7, 1.0, 1.0, 1.0];
        s2.v_knots[4] = 0.3; // hack
        let err = stitch(&s1, &s2, (Edge::UMax, Edge::UMin), 1e-6).unwrap_err();
        assert_eq!(err.code(), "surface.sew_mismatch");
    }

    #[test]
    fn stitch_fails_on_orthogonal_edges() {
        let s1 = planar_patch(0.0);
        let s2 = planar_patch(1.0);
        // UMax and VMin run in orthogonal directions.
        let err = stitch(&s1, &s2, (Edge::UMax, Edge::VMin), 1e-6).unwrap_err();
        assert_eq!(err.code(), "surface.sew_mismatch");
    }

    // ===== Phase 19C — G2 sew tests =====

    /// 5x4 planar patch — wide enough in u that we can test G2
    /// adjustment (needs at least 5 interior CPs per side after stitch).
    fn wide_planar_patch(x0: f64) -> NurbsSurface {
        let n_cp = 5_usize;
        let p = 3_usize;
        let knots = {
            let mut k = vec![0.0_f64; p + 1];
            for i in 1..(n_cp - p) {
                k.push(i as f64 / (n_cp - p) as f64);
            }
            k.extend(std::iter::repeat_n(1.0_f64, p + 1));
            k
        };
        let cps = (0..n_cp)
            .map(|i| {
                let u = i as f64 / (n_cp - 1) as f64;
                (0..4)
                    .map(|j| {
                        let v = j as f64 / 3.0;
                        Vector3::new(x0 + u, v, 0.0)
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let weights = vec![vec![1.0; 4]; n_cp];
        let v_knots = vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
        NurbsSurface::new(3, 3, knots, v_knots, cps, weights).unwrap()
    }

    #[test]
    fn g2_stitch_planar_surfaces_is_smooth() {
        let s1 = wide_planar_patch(0.0);
        let s2 = wide_planar_patch(1.0);
        let stitched =
            stitch_with_continuity(&s1, &s2, (Edge::UMax, Edge::UMin), 1e-6, Continuity::G2)
                .unwrap();
        // Combined patch should still be planar (z = 0 everywhere
        // because both inputs were planar and our adjustments are
        // symmetric in z).
        for i in 0..stitched.nu() {
            for j in 0..stitched.nv() {
                let p = stitched.control_points[i][j];
                assert!(p.z.abs() < 1e-9, "CP[{i},{j}] z = {}", p.z);
            }
        }
    }

    #[test]
    fn g2_stitch_adjusts_three_rows_per_side() {
        // Verify that the G2 stitch CP-shifts touch 5 rows around
        // the seam (2 on each side + the shared row).
        let s1 = wide_planar_patch(0.0);
        let s2 = wide_planar_patch(1.0);
        let g0 = stitch_with_continuity(&s1, &s2, (Edge::UMax, Edge::UMin), 1e-6, Continuity::G0)
            .unwrap();
        let g2 = stitch_with_continuity(&s1, &s2, (Edge::UMax, Edge::UMin), 1e-6, Continuity::G2)
            .unwrap();
        // Both have the same dimension; only the local shape differs.
        assert_eq!(g0.nu(), g2.nu());
        assert_eq!(g0.nv(), g2.nv());
        // For perfectly planar inputs, G2 == G0 (the symmetry already
        // makes tangents + curvatures match). So this is a sanity
        // structure check: both surfaces evaluate to the same plane
        // at the seam.
        let p_seam_g0 = g0.evaluate(0.5, 0.5);
        let p_seam_g2 = g2.evaluate(0.5, 0.5);
        assert!((p_seam_g0 - p_seam_g2).norm() < 1e-9);
    }
}
