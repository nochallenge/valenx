//! Knot insertion, knot removal, and degree elevation for NURBS
//! curves and surfaces (Phase 19A).
//!
//! ## Algorithms
//!
//! - **Boehm's algorithm** for single knot insertion (Piegl & Tiller,
//!   *The NURBS Book*, Algorithm A5.1). Insert a new knot `u_bar` into
//!   the knot vector; the curve geometry is unchanged but the control
//!   polygon is refined — `degree` new CPs replace `degree`
//!   neighbours. We work in homogeneous (wP, w) coordinates so the
//!   rational basis stays valid.
//! - **Tiller's algorithm** for knot removal (Piegl & Tiller A5.8).
//!   Try to remove an interior knot; accept the removal only if the
//!   resulting curve stays within `tolerance` of the original.
//! - **Degree elevation** by 1 (Piegl & Tiller A5.9 simplified): the
//!   Prautzsch closed-form (special case for Bezier segments) applied
//!   per Bezier segment after Bezier-decomposition. v1 implementation
//!   elevates a clamped curve segment-by-segment using the
//!   binomial-combination identity:
//!
//!   ```text
//!   Q_i = sum_{j=max(0, i-1)}^{min(p, i)} (C(p, j) * C(1, i-j) / C(p+1, i)) * P_j
//!   ```
//!
//!   (the formula for elevating a Bezier of degree p by 1 to a Bezier
//!   of degree p+1; we apply this per-segment after splitting at
//!   internal knots, then stitch back into a single open-uniform
//!   knot vector).
//!
//! All operations preserve the rational basis by carrying weights
//! through the homogeneous projection.

use nalgebra::Vector4;

use crate::error::SurfaceError;
use crate::nurbs_curve::NurbsCurve;
use crate::nurbs_surface::NurbsSurface;

// ===== Curve-level operations =====

impl NurbsCurve {
    /// Insert knot `u_bar` into this curve via Boehm's algorithm.
    ///
    /// Geometry is preserved: `new.evaluate(u) == self.evaluate(u)` to
    /// floating-point precision. Insertion multiplies the local CP
    /// neighbourhood by knot-distance ratios and produces `degree`
    /// new CPs that replace `degree - 1` originals.
    ///
    /// Returns an error if `u_bar` is outside the parameter range or
    /// would push knot multiplicity above `degree + 1` (where the
    /// curve becomes degenerate).
    pub fn insert_knot(&self, u_bar: f64) -> Result<Self, SurfaceError> {
        let p = self.degree;
        let (u_min, u_max) = self.parameter_range();
        if !(u_min..=u_max).contains(&u_bar) {
            return Err(SurfaceError::EvaluationOutOfRange { u: u_bar });
        }

        // Find the span k such that knots[k] <= u_bar < knots[k+1].
        // For u_bar == u_max, find_knot_span clamps to n-1; that's
        // already correct for the insertion algorithm.
        let k = self.find_knot_span(u_bar);

        // Current multiplicity of u_bar in the knot vector.
        let mut s = 0_usize;
        for &kn in &self.knots {
            if (kn - u_bar).abs() < 1.0e-12 {
                s += 1;
            }
        }
        if s + 1 > p + 1 {
            return Err(SurfaceError::BadKnotVector {
                reason: format!(
                    "inserting u_bar={u_bar} would exceed max multiplicity p+1={}",
                    p + 1
                ),
            });
        }

        // Work in homogeneous coordinates (wP, w) so rational weights
        // carry through correctly.
        let n_old = self.control_points.len();
        let mut q: Vec<Vector4<f64>> = Vec::with_capacity(n_old + 1);

        // Copy the unaffected leading CPs.
        for i in 0..=(k - p) {
            q.push(homog(&self.control_points[i], self.weights[i]));
        }
        // Compute the p new CPs in the affected range.
        for i in (k - p + 1)..=(k - s) {
            let alpha = (u_bar - self.knots[i]) / (self.knots[i + p - s] - self.knots[i]);
            let prev = homog(&self.control_points[i - 1], self.weights[i - 1]);
            let curr = homog(&self.control_points[i], self.weights[i]);
            q.push((1.0 - alpha) * prev + alpha * curr);
        }
        // Copy the unaffected trailing CPs.
        for i in (k - s)..n_old {
            q.push(homog(&self.control_points[i], self.weights[i]));
        }

        // Build new knot vector with u_bar inserted before knots[k+1].
        let mut new_knots = Vec::with_capacity(self.knots.len() + 1);
        new_knots.extend_from_slice(&self.knots[..=k]);
        new_knots.push(u_bar);
        new_knots.extend_from_slice(&self.knots[k + 1..]);

        let (cps, weights) = dehomog_vec(&q);
        NurbsCurve::new(p, new_knots, cps, weights)
    }

    /// Try to remove one instance of the knot at index `r` (a knot
    /// of multiplicity `s` at parameter value `u_r = knots[r]`). The
    /// removal is accepted only if the maximum deviation of the new
    /// curve from the old curve at sample points is `<= tolerance`.
    ///
    /// Returns the new curve on success, or the original error from
    /// validation if removal would corrupt the curve.
    ///
    /// v1 simplification: we accept the removal unconditionally when
    /// the local-error proxy `alfi * |P[i+1] - 0.5*(P[i]+P[i+2])|` is
    /// `<= tolerance`. The full Tiller-Wang test (sample-based
    /// deviation) is a v1.5 upgrade.
    pub fn remove_knot(&self, u: f64, tolerance: f64) -> Result<Self, SurfaceError> {
        let p = self.degree;
        let (u_min, u_max) = self.parameter_range();
        if !(u_min..=u_max).contains(&u) {
            return Err(SurfaceError::EvaluationOutOfRange { u });
        }
        // Find the highest knot index r such that knots[r] == u.
        let mut r: Option<usize> = None;
        let mut s: usize = 0;
        for (i, &kn) in self.knots.iter().enumerate() {
            if (kn - u).abs() < 1.0e-12 {
                r = Some(i);
                s += 1;
            }
        }
        let r = r.ok_or_else(|| SurfaceError::BadKnotVector {
            reason: format!("knot {u} not found in knot vector"),
        })?;
        if s == 0 {
            return Err(SurfaceError::BadKnotVector {
                reason: format!("knot {u} has multiplicity 0"),
            });
        }
        // Removing past multiplicity p+1 makes the curve degenerate.
        if s > p {
            return Err(SurfaceError::BadKnotVector {
                reason: format!("cannot remove knot at multiplicity {s} (> degree {p})"),
            });
        }

        // Piegl-Tiller A5.8 simplified: compute the new CP triangle
        // and accept if the linearised error is within tolerance.
        let n = self.control_points.len() - 1; // last CP index
        let m = self.knots.len() - 1; // last knot index

        // First / last affected CP indices.
        let first = r - p;
        let last = r - s;
        if last + 1 > n {
            return Err(SurfaceError::BadKnotVector {
                reason: "removal indices out of bounds".into(),
            });
        }

        // Work in homogeneous coords.
        let mut q: Vec<Vector4<f64>> = self
            .control_points
            .iter()
            .zip(&self.weights)
            .map(|(p, w)| homog(p, *w))
            .collect();

        // Local triangle of new CPs as we shift the row.
        // We do a single pass (t=0 in Piegl-Tiller language: remove
        // one instance of the knot).
        let mut temp = vec![Vector4::zeros(); 2 * p + 1];
        let mut i = first;
        let mut j = last;
        let mut ii: i64 = 0;
        let mut jj: i64 = (j as i64) - (i as i64);
        temp[ii as usize] = q[(first as i64 - 1).max(0) as usize];
        temp[jj as usize + 1] = q[last + 1];
        while (j as i64) - (i as i64) > 0 {
            let denom_i = self.knots[i + p + 1] - self.knots[i];
            let denom_j = self.knots[j + p + 1] - self.knots[j];
            let alfi = if denom_i.abs() < 1.0e-12 {
                0.0
            } else {
                (u - self.knots[i]) / denom_i
            };
            let alfj = if denom_j.abs() < 1.0e-12 {
                0.0
            } else {
                (u - self.knots[j]) / denom_j
            };
            temp[ii as usize + 1] = if alfi.abs() < 1.0e-12 {
                temp[ii as usize]
            } else {
                (q[i] - (1.0 - alfi) * temp[ii as usize]) / alfi
            };
            temp[jj as usize] = if (1.0 - alfj).abs() < 1.0e-12 {
                temp[jj as usize + 1]
            } else {
                (q[j] - alfj * temp[jj as usize + 1]) / (1.0 - alfj)
            };
            i += 1;
            j = j.saturating_sub(1);
            ii += 1;
            jj -= 1;
        }
        // Linearised error proxy: deviation between the two computed
        // candidate CPs at the meeting index.
        let err = if (j as i64) - (i as i64) < 0 {
            (temp[ii as usize] - temp[jj as usize + 1]).norm()
        } else {
            let alfi = (u - self.knots[i]) / (self.knots[i + p + 1] - self.knots[i]).max(1.0e-12);
            (q[i] - (alfi * temp[jj as usize + 1] + (1.0 - alfi) * temp[ii as usize])).norm()
        };
        if err > tolerance {
            return Err(SurfaceError::BadKnotVector {
                reason: format!(
                    "knot removal at u={u} would exceed tolerance {tolerance} (err={err:.2e})"
                ),
            });
        }

        // Apply: rewrite the affected CP slice. The compute loop stores
        // the new left-side CP for index `i` at `temp[ii + 1]` (with
        // `ii = i - first`) and the new right-side CP for index `j` at
        // `temp[jj] = temp[j - first]`, so the writes mirror those.
        i = first;
        j = last;
        while (j as i64) - (i as i64) > 0 {
            q[i] = temp[i - first + 1];
            q[j] = temp[j - first];
            i += 1;
            j = j.saturating_sub(1);
        }

        // Drop the now-redundant CP. After the loop `i` is the meeting
        // index — the single stale CP for an even-width affected span
        // (`i == j`), and the duplicated meeting CP for an odd span
        // (`i == j + 1`). Either way index `i` is the one to remove.
        let _ = (n, m); // silence unused warnings on some configurations
        let cp_remove_idx = i;
        q.remove(cp_remove_idx.min(q.len() - 1));
        let mut new_knots = self.knots.clone();
        new_knots.remove(r);
        let (cps, weights) = dehomog_vec(&q);
        NurbsCurve::new(p, new_knots, cps, weights)
    }

    /// Elevate degree by `by` (must be >= 1). v1 implementation
    /// elevates one level at a time by decomposing into Bezier
    /// segments via knot insertion to multiplicity `degree`, elevating
    /// each segment using the closed-form Bezier-elevation formula,
    /// then collecting the result back into a single clamped knot
    /// vector.
    pub fn elevate_degree(&self, by: usize) -> Result<Self, SurfaceError> {
        if by == 0 {
            return Ok(self.clone());
        }
        let mut current = self.clone();
        for _ in 0..by {
            current = current.elevate_degree_by_one()?;
        }
        Ok(current)
    }

    fn elevate_degree_by_one(&self) -> Result<Self, SurfaceError> {
        let p = self.degree;
        let new_p = p + 1;

        // 1. Decompose to Bezier segments by inserting every internal
        //    knot to multiplicity p. The number of Bezier segments
        //    equals `n_segments = (#distinct interior knots) + 1`.
        let segments = bezier_decompose_curve(self)?;

        // 2. Elevate each Bezier segment from degree p to degree p+1.
        let mut elevated_segments: Vec<Vec<Vector4<f64>>> = Vec::with_capacity(segments.len());
        for seg in &segments {
            let q = elevate_bezier_segment(seg, p);
            elevated_segments.push(q);
        }

        // 3. Stitch back. The shared CP at the join of two segments
        //    is identical in both, so we drop the duplicate.
        let mut joined: Vec<Vector4<f64>> = Vec::new();
        for (k, seg) in elevated_segments.iter().enumerate() {
            if k == 0 {
                joined.extend_from_slice(seg);
            } else {
                joined.extend_from_slice(&seg[1..]);
            }
        }

        // 4. Rebuild a clamped knot vector. The interior knots are
        //    the original distinct interior knots, each with
        //    multiplicity new_p (since each Bezier segment has
        //    multiplicity new_p at its join).
        let distinct_interior = distinct_interior_knots(&self.knots, p);
        let n_new_cps = joined.len();
        let mut new_knots = Vec::with_capacity(n_new_cps + new_p + 1);
        for _ in 0..=new_p {
            new_knots.push(*self.knots.first().unwrap());
        }
        for &u in &distinct_interior {
            for _ in 0..new_p {
                new_knots.push(u);
            }
        }
        for _ in 0..=new_p {
            new_knots.push(*self.knots.last().unwrap());
        }
        // The above may produce more knots than expected if joined
        // already collapsed; trim or pad to `n_new_cps + new_p + 1`.
        let target = n_new_cps + new_p + 1;
        if new_knots.len() > target {
            // Trim from the interior (after the leading clamp).
            let excess = new_knots.len() - target;
            new_knots.drain(new_p + 1..new_p + 1 + excess);
        } else if new_knots.len() < target {
            // Pad by repeating the last knot.
            let last = *new_knots.last().unwrap();
            while new_knots.len() < target {
                new_knots.push(last);
            }
        }

        let (cps, weights) = dehomog_vec(&joined);
        NurbsCurve::new(new_p, new_knots, cps, weights)
    }
}

// ===== Surface-level operations =====

impl NurbsSurface {
    /// Insert knot `u_bar` into the u-knot vector. Geometry preserved.
    pub fn insert_knot_u(&self, u_bar: f64) -> Result<Self, SurfaceError> {
        // Treat each v-isoparametric strip as a NURBS curve in u and
        // insert the knot there, then reassemble.
        let nv = self.nv();
        let mut new_rows: Vec<Vec<Vector4<f64>>> = Vec::with_capacity(self.nu() + 1);
        let mut new_u_knots: Option<Vec<f64>> = None;

        for j in 0..nv {
            let cps_u: Vec<nalgebra::Vector3<f64>> =
                (0..self.nu()).map(|i| self.control_points[i][j]).collect();
            let ws_u: Vec<f64> = (0..self.nu()).map(|i| self.weights[i][j]).collect();
            let curve = NurbsCurve::new(self.u_degree, self.u_knots.clone(), cps_u, ws_u)?;
            let inserted = curve.insert_knot(u_bar)?;
            if new_u_knots.is_none() {
                new_u_knots = Some(inserted.knots.clone());
                for _ in 0..inserted.control_points.len() {
                    new_rows.push(Vec::with_capacity(nv));
                }
            }
            for (i, (cp, w)) in inserted
                .control_points
                .iter()
                .zip(inserted.weights.iter())
                .enumerate()
            {
                new_rows[i].push(homog(cp, *w));
            }
        }
        let new_u_knots = new_u_knots.expect("at least one column processed");
        // Build cps/weights grids.
        let new_nu = new_rows.len();
        let mut cps = Vec::with_capacity(new_nu);
        let mut weights = Vec::with_capacity(new_nu);
        for row in &new_rows {
            let mut r_cps = Vec::with_capacity(nv);
            let mut r_ws = Vec::with_capacity(nv);
            for h in row {
                let (p, w) = dehomog(h);
                r_cps.push(p);
                r_ws.push(w);
            }
            cps.push(r_cps);
            weights.push(r_ws);
        }
        NurbsSurface::new(
            self.u_degree,
            self.v_degree,
            new_u_knots,
            self.v_knots.clone(),
            cps,
            weights,
        )
    }

    /// Insert knot `v_bar` into the v-knot vector. Geometry preserved.
    pub fn insert_knot_v(&self, v_bar: f64) -> Result<Self, SurfaceError> {
        // Symmetric to insert_knot_u: treat each u-isoparametric row
        // as a NURBS curve in v.
        let nu = self.nu();
        let mut new_cols: Vec<Vec<Vector4<f64>>> = Vec::with_capacity(self.nv() + 1);
        let mut new_v_knots: Option<Vec<f64>> = None;

        for i in 0..nu {
            let cps_v = self.control_points[i].clone();
            let ws_v = self.weights[i].clone();
            let curve = NurbsCurve::new(self.v_degree, self.v_knots.clone(), cps_v, ws_v)?;
            let inserted = curve.insert_knot(v_bar)?;
            if new_v_knots.is_none() {
                new_v_knots = Some(inserted.knots.clone());
                for _ in 0..inserted.control_points.len() {
                    new_cols.push(Vec::with_capacity(nu));
                }
            }
            for (j, (cp, w)) in inserted
                .control_points
                .iter()
                .zip(inserted.weights.iter())
                .enumerate()
            {
                new_cols[j].push(homog(cp, *w));
            }
        }
        let new_v_knots = new_v_knots.expect("at least one row processed");
        let new_nv = new_cols.len();
        let mut cps = vec![Vec::with_capacity(new_nv); nu];
        let mut weights = vec![Vec::with_capacity(new_nv); nu];
        for col in &new_cols {
            for (i, h) in col.iter().enumerate() {
                let (p, w) = dehomog(h);
                cps[i].push(p);
                weights[i].push(w);
            }
        }
        NurbsSurface::new(
            self.u_degree,
            self.v_degree,
            self.u_knots.clone(),
            new_v_knots,
            cps,
            weights,
        )
    }

    /// Elevate u-degree by `by`. Symmetric to curve elevation.
    pub fn elevate_degree_u(&self, by: usize) -> Result<Self, SurfaceError> {
        if by == 0 {
            return Ok(self.clone());
        }
        let nv = self.nv();
        // Elevate every v-isoparametric strip independently.
        let mut elevated_curves: Vec<NurbsCurve> = Vec::with_capacity(nv);
        for j in 0..nv {
            let cps_u: Vec<nalgebra::Vector3<f64>> =
                (0..self.nu()).map(|i| self.control_points[i][j]).collect();
            let ws_u: Vec<f64> = (0..self.nu()).map(|i| self.weights[i][j]).collect();
            let curve = NurbsCurve::new(self.u_degree, self.u_knots.clone(), cps_u, ws_u)?;
            elevated_curves.push(curve.elevate_degree(by)?);
        }
        let new_u_knots = elevated_curves[0].knots.clone();
        let new_nu = elevated_curves[0].control_points.len();
        let mut cps = vec![Vec::with_capacity(nv); new_nu];
        let mut weights = vec![Vec::with_capacity(nv); new_nu];
        for curve in &elevated_curves {
            for (i, (cp, w)) in curve
                .control_points
                .iter()
                .zip(curve.weights.iter())
                .enumerate()
            {
                cps[i].push(*cp);
                weights[i].push(*w);
            }
        }
        NurbsSurface::new(
            self.u_degree + by,
            self.v_degree,
            new_u_knots,
            self.v_knots.clone(),
            cps,
            weights,
        )
    }

    /// Elevate v-degree by `by`. Symmetric to `elevate_degree_u`.
    pub fn elevate_degree_v(&self, by: usize) -> Result<Self, SurfaceError> {
        if by == 0 {
            return Ok(self.clone());
        }
        let nu = self.nu();
        let mut elevated_curves: Vec<NurbsCurve> = Vec::with_capacity(nu);
        for i in 0..nu {
            let curve = NurbsCurve::new(
                self.v_degree,
                self.v_knots.clone(),
                self.control_points[i].clone(),
                self.weights[i].clone(),
            )?;
            elevated_curves.push(curve.elevate_degree(by)?);
        }
        let new_v_knots = elevated_curves[0].knots.clone();
        let new_nv = elevated_curves[0].control_points.len();
        let mut cps = Vec::with_capacity(nu);
        let mut weights = Vec::with_capacity(nu);
        for curve in &elevated_curves {
            cps.push(curve.control_points.clone());
            weights.push(curve.weights.clone());
        }
        let _ = new_nv; // silence
        NurbsSurface::new(
            self.u_degree,
            self.v_degree + by,
            self.u_knots.clone(),
            new_v_knots,
            cps,
            weights,
        )
    }
}

// ===== Helpers =====

fn homog(p: &nalgebra::Vector3<f64>, w: f64) -> Vector4<f64> {
    Vector4::new(p.x * w, p.y * w, p.z * w, w)
}

fn dehomog(h: &Vector4<f64>) -> (nalgebra::Vector3<f64>, f64) {
    let w = h.w;
    if w.abs() < 1.0e-30 {
        (nalgebra::Vector3::new(h.x, h.y, h.z), 1.0)
    } else {
        (nalgebra::Vector3::new(h.x / w, h.y / w, h.z / w), w)
    }
}

fn dehomog_vec(hs: &[Vector4<f64>]) -> (Vec<nalgebra::Vector3<f64>>, Vec<f64>) {
    let mut cps = Vec::with_capacity(hs.len());
    let mut ws = Vec::with_capacity(hs.len());
    for h in hs {
        let (p, w) = dehomog(h);
        cps.push(p);
        ws.push(w);
    }
    (cps, ws)
}

/// Return the distinct interior knot values of a clamped knot vector
/// of degree `p`. Interior = strictly between the first multiplicity
/// (degree+1) and the last multiplicity (degree+1).
fn distinct_interior_knots(knots: &[f64], p: usize) -> Vec<f64> {
    if knots.len() < 2 * (p + 1) {
        return Vec::new();
    }
    let mut out = Vec::new();
    let mut last: Option<f64> = None;
    for &k in &knots[p + 1..knots.len() - p - 1] {
        if last.map(|l| (k - l).abs() > 1.0e-12).unwrap_or(true) {
            out.push(k);
            last = Some(k);
        }
    }
    out
}

/// Decompose a clamped NURBS curve into a sequence of Bezier segments
/// (in homogeneous coords). Each segment has `degree + 1` CPs.
fn bezier_decompose_curve(c: &NurbsCurve) -> Result<Vec<Vec<Vector4<f64>>>, SurfaceError> {
    let p = c.degree;
    // Insert each distinct interior knot up to multiplicity p.
    let mut current = c.clone();
    let interior = distinct_interior_knots(&c.knots, p);
    for &u in &interior {
        // Count current multiplicity.
        let mut mult = 0;
        for &kn in &current.knots {
            if (kn - u).abs() < 1.0e-12 {
                mult += 1;
            }
        }
        for _ in mult..p {
            current = current.insert_knot(u)?;
        }
    }
    // Now slice into (p+1)-CP groups, with consecutive groups sharing
    // the boundary CP.
    let n_segments = interior.len() + 1;
    let mut segments = Vec::with_capacity(n_segments);
    let stride = p; // each segment starts `p` after the previous start
    for seg_idx in 0..n_segments {
        let start = seg_idx * stride;
        let mut seg = Vec::with_capacity(p + 1);
        for k in 0..=p {
            let idx = start + k;
            if idx >= current.control_points.len() {
                return Err(SurfaceError::BadKnotVector {
                    reason: format!(
                        "bezier decomposition: CP index {idx} out of range ({})",
                        current.control_points.len()
                    ),
                });
            }
            seg.push(homog(&current.control_points[idx], current.weights[idx]));
        }
        segments.push(seg);
    }
    Ok(segments)
}

/// Elevate a single Bezier segment (in homogeneous coords) from
/// degree `p` to degree `p + 1`. Returns the `p + 2` new homogeneous
/// CPs.
fn elevate_bezier_segment(seg: &[Vector4<f64>], p: usize) -> Vec<Vector4<f64>> {
    // Q_0 = P_0
    // Q_i = (i / (p+1)) * P_{i-1} + (1 - i / (p+1)) * P_i,  1 <= i <= p
    // Q_{p+1} = P_p
    let new_p = p + 1;
    let mut q = Vec::with_capacity(new_p + 1);
    q.push(seg[0]);
    for i in 1..=p {
        let alpha = (i as f64) / ((p + 1) as f64);
        q.push(alpha * seg[i - 1] + (1.0 - alpha) * seg[i]);
    }
    q.push(seg[p]);
    q
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    fn cubic_bezier() -> NurbsCurve {
        let knots = vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
        let cps = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 2.0, 0.0),
            Vector3::new(2.0, 2.0, 0.0),
            Vector3::new(3.0, 0.0, 0.0),
        ];
        let weights = vec![1.0; 4];
        NurbsCurve::new(3, knots, cps, weights).unwrap()
    }

    #[test]
    fn insert_knot_preserves_geometry_on_cubic_bezier() {
        let c = cubic_bezier();
        let inserted = c.insert_knot(0.5).unwrap();
        // The curve evaluated at the same parameters must match.
        for &u in &[0.0_f64, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0] {
            let a = c.evaluate(u);
            let b = inserted.evaluate(u);
            assert!((a - b).norm() < 1.0e-8, "u={u}: {a:?} vs {b:?}");
        }
        // Knot vector grew by 1.
        assert_eq!(inserted.knots.len(), c.knots.len() + 1);
        // CP count grew by 1.
        assert_eq!(inserted.control_points.len(), c.control_points.len() + 1);
    }

    #[test]
    fn insert_knot_rejects_out_of_range() {
        let c = cubic_bezier();
        assert!(c.insert_knot(-0.1).is_err());
        assert!(c.insert_knot(1.1).is_err());
    }

    #[test]
    fn insert_then_remove_round_trips() {
        let c = cubic_bezier();
        let inserted = c.insert_knot(0.4).unwrap();
        let removed = inserted.remove_knot(0.4, 1.0e-6).unwrap();
        // Should match original CPs within tolerance.
        assert_eq!(removed.control_points.len(), c.control_points.len());
        for (a, b) in removed.control_points.iter().zip(&c.control_points) {
            assert!((a - b).norm() < 1.0e-6, "{a:?} vs {b:?}");
        }
    }

    #[test]
    fn elevate_degree_preserves_geometry_on_cubic_bezier() {
        let c = cubic_bezier();
        let elevated = c.elevate_degree(1).unwrap();
        assert_eq!(elevated.degree, 4);
        // Sample-match.
        for &u in &[0.0_f64, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0] {
            let a = c.evaluate(u);
            let b = elevated.evaluate(u);
            assert!((a - b).norm() < 1.0e-8, "u={u}: {a:?} vs {b:?}");
        }
    }

    fn planar_xy_surface() -> NurbsSurface {
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

    #[test]
    fn surface_insert_knot_u_preserves_geometry() {
        let s = planar_xy_surface();
        let s2 = s.insert_knot_u(0.4).unwrap();
        for &(u, v) in &[(0.1_f64, 0.3_f64), (0.5, 0.5), (0.9, 0.2)] {
            let a = s.evaluate(u, v);
            let b = s2.evaluate(u, v);
            assert!((a - b).norm() < 1.0e-8, "({u},{v}): {a:?} vs {b:?}");
        }
    }

    #[test]
    fn surface_insert_knot_v_preserves_geometry() {
        let s = planar_xy_surface();
        let s2 = s.insert_knot_v(0.6).unwrap();
        for &(u, v) in &[(0.1_f64, 0.3_f64), (0.5, 0.5), (0.9, 0.2)] {
            let a = s.evaluate(u, v);
            let b = s2.evaluate(u, v);
            assert!((a - b).norm() < 1.0e-8, "({u},{v}): {a:?} vs {b:?}");
        }
    }

    #[test]
    fn surface_elevate_degree_u_preserves_geometry() {
        let s = planar_xy_surface();
        let s2 = s.elevate_degree_u(1).unwrap();
        assert_eq!(s2.u_degree, 4);
        for &(u, v) in &[(0.1_f64, 0.3_f64), (0.5, 0.5), (0.9, 0.2)] {
            let a = s.evaluate(u, v);
            let b = s2.evaluate(u, v);
            assert!((a - b).norm() < 1.0e-8, "({u},{v}): {a:?} vs {b:?}");
        }
    }

    #[test]
    fn surface_elevate_degree_v_preserves_geometry() {
        let s = planar_xy_surface();
        let s2 = s.elevate_degree_v(1).unwrap();
        assert_eq!(s2.v_degree, 4);
        for &(u, v) in &[(0.1_f64, 0.3_f64), (0.5, 0.5), (0.9, 0.2)] {
            let a = s.evaluate(u, v);
            let b = s2.evaluate(u, v);
            assert!((a - b).norm() < 1.0e-8, "({u},{v}): {a:?} vs {b:?}");
        }
    }
}
