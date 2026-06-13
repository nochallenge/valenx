//! NURBS surface — bidirectional tensor product of two NURBS basis
//! sets.
//!
//! Stored as:
//! - `u_degree`, `v_degree`: polynomial degrees in u + v directions,
//! - `u_knots`, `v_knots`: knot vectors,
//! - `control_points[i][j]`: outer index is u, inner is v
//!   (i.e. `control_points[i]` is the i-th *row* of CPs in v),
//! - `weights[i][j]`: same indexing as `control_points`.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::error::SurfaceError;
use crate::nurbs_curve::{basis_functions, find_knot_span};

/// A 3D NURBS surface.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct NurbsSurface {
    /// Degree in the u parameter.
    pub u_degree: usize,
    /// Degree in the v parameter.
    pub v_degree: usize,
    /// Knot vector in u — non-decreasing, length `nu + u_degree + 1`.
    pub u_knots: Vec<f64>,
    /// Knot vector in v — non-decreasing, length `nv + v_degree + 1`.
    pub v_knots: Vec<f64>,
    /// 2D grid of control points indexed as `[i (u)][j (v)]`.
    pub control_points: Vec<Vec<Vector3<f64>>>,
    /// 2D grid of weights with the same indexing as `control_points`.
    pub weights: Vec<Vec<f64>>,
}

impl NurbsSurface {
    /// Construct a validated NURBS surface.
    ///
    /// Validates degrees, knot-vector lengths, knot monotonicity,
    /// and the rectangularity of the control-point grid.
    pub fn new(
        u_degree: usize,
        v_degree: usize,
        u_knots: Vec<f64>,
        v_knots: Vec<f64>,
        control_points: Vec<Vec<Vector3<f64>>>,
        weights: Vec<Vec<f64>>,
    ) -> Result<Self, SurfaceError> {
        let surface = Self::new_unchecked(
            u_degree,
            v_degree,
            u_knots,
            v_knots,
            control_points,
            weights,
        );
        surface.validate()?;
        Ok(surface)
    }

    /// Check the NURBS-surface invariant on an already-built surface:
    /// valid u/v degrees, a non-empty rectangular control-point grid with
    /// enough CPs per direction, knot-vector lengths `n + degree + 1`,
    /// non-decreasing knots, and a weights grid matching the CP grid.
    /// [`Self::new`] runs this at construction; a surface obtained another
    /// way (deserialised, or via [`Self::new_unchecked`]) can re-check
    /// itself with this before its indexing evaluation methods are used.
    pub fn validate(&self) -> Result<(), SurfaceError> {
        if !(1..=9).contains(&self.u_degree) {
            return Err(SurfaceError::BadDegree(self.u_degree));
        }
        if !(1..=9).contains(&self.v_degree) {
            return Err(SurfaceError::BadDegree(self.v_degree));
        }

        let nu = self.control_points.len();
        if nu == 0 {
            return Err(SurfaceError::BadKnotVector {
                reason: "empty control_points grid".into(),
            });
        }
        let nv = self.control_points[0].len();
        for row in &self.control_points {
            if row.len() != nv {
                return Err(SurfaceError::BadKnotVector {
                    reason: "control_points grid is not rectangular".into(),
                });
            }
        }
        if nu < self.u_degree + 1 {
            return Err(SurfaceError::BadKnotVector {
                reason: format!(
                    "need at least {} u-direction CPs for degree {}, got {}",
                    self.u_degree + 1,
                    self.u_degree,
                    nu
                ),
            });
        }
        if nv < self.v_degree + 1 {
            return Err(SurfaceError::BadKnotVector {
                reason: format!(
                    "need at least {} v-direction CPs for degree {}, got {}",
                    self.v_degree + 1,
                    self.v_degree,
                    nv
                ),
            });
        }

        let u_expected = nu + self.u_degree + 1;
        if self.u_knots.len() != u_expected {
            return Err(SurfaceError::BadKnotVector {
                reason: format!("u_knots: expected {u_expected}, got {}", self.u_knots.len()),
            });
        }
        let v_expected = nv + self.v_degree + 1;
        if self.v_knots.len() != v_expected {
            return Err(SurfaceError::BadKnotVector {
                reason: format!("v_knots: expected {v_expected}, got {}", self.v_knots.len()),
            });
        }
        for w in self.u_knots.windows(2) {
            if w[1] < w[0] {
                return Err(SurfaceError::BadKnotVector {
                    reason: "u_knots must be non-decreasing".into(),
                });
            }
        }
        for w in self.v_knots.windows(2) {
            if w[1] < w[0] {
                return Err(SurfaceError::BadKnotVector {
                    reason: "v_knots must be non-decreasing".into(),
                });
            }
        }

        if self.weights.len() != nu || self.weights.iter().any(|r| r.len() != nv) {
            return Err(SurfaceError::BadKnotVector {
                reason: "weights grid shape ≠ control_points grid shape".into(),
            });
        }

        Ok(())
    }

    /// Skip validation — caller asserts well-formedness.
    pub fn new_unchecked(
        u_degree: usize,
        v_degree: usize,
        u_knots: Vec<f64>,
        v_knots: Vec<f64>,
        control_points: Vec<Vec<Vector3<f64>>>,
        weights: Vec<Vec<f64>>,
    ) -> Self {
        Self {
            u_degree,
            v_degree,
            u_knots,
            v_knots,
            control_points,
            weights,
        }
    }

    /// Number of control points in u.
    pub fn nu(&self) -> usize {
        self.control_points.len()
    }

    /// Number of control points in v.
    pub fn nv(&self) -> usize {
        self.control_points[0].len()
    }

    /// Valid u-parameter range.
    pub fn u_range(&self) -> (f64, f64) {
        (self.u_knots[self.u_degree], self.u_knots[self.nu()])
    }

    /// Valid v-parameter range.
    pub fn v_range(&self) -> (f64, f64) {
        (self.v_knots[self.v_degree], self.v_knots[self.nv()])
    }

    /// Evaluate the surface at `(u, v)` using the tensor product of
    /// u and v basis functions.
    ///
    /// `S(u, v) = Σ_i Σ_j N_i^p(u) N_j^q(v) w_ij P_ij / Σ_i Σ_j N_i^p(u) N_j^q(v) w_ij`
    pub fn evaluate(&self, u: f64, v: f64) -> Vector3<f64> {
        let nu = self.nu();
        let nv = self.nv();
        let span_u = find_knot_span(u, &self.u_knots, self.u_degree, nu);
        let span_v = find_knot_span(v, &self.v_knots, self.v_degree, nv);
        let basis_u = basis_functions(span_u, u, self.u_degree, &self.u_knots);
        let basis_v = basis_functions(span_v, v, self.v_degree, &self.v_knots);

        let mut num = Vector3::zeros();
        let mut den = 0.0_f64;
        for (i, bu) in basis_u.iter().enumerate() {
            let u_idx = span_u - self.u_degree + i;
            for (j, bv) in basis_v.iter().enumerate() {
                let v_idx = span_v - self.v_degree + j;
                let w = self.weights[u_idx][v_idx];
                let b = bu * bv;
                let wb = w * b;
                num += self.control_points[u_idx][v_idx] * wb;
                den += wb;
            }
        }
        if den.abs() < 1e-30 {
            num
        } else {
            num / den
        }
    }

    /// The **first partial derivative** `∂S/∂u` at `(u, v)` — the surface
    /// tangent along the u iso-curve (the direction of travel as `u` increases
    /// at fixed `v`).
    ///
    /// Central finite difference over the valid [`u_range`](Self::u_range),
    /// with the stencil clamped (one-sided) at the domain edges so it never
    /// evaluates out of range; 6–7 digits of accuracy, the surface analogue of
    /// the curve's finite-difference derivative.
    pub fn partial_u(&self, u: f64, v: f64) -> Vector3<f64> {
        let (u_min, u_max) = self.u_range();
        let h = ((u_max - u_min) * 1e-4).max(1e-9);
        let u_lo = (u - h).max(u_min);
        let u_hi = (u + h).min(u_max);
        let d = u_hi - u_lo;
        if d.abs() < 1e-30 {
            return Vector3::zeros();
        }
        (self.evaluate(u_hi, v) - self.evaluate(u_lo, v)) / d
    }

    /// The **first partial derivative** `∂S/∂v` at `(u, v)` — the surface
    /// tangent along the v iso-curve. Central finite difference over the valid
    /// [`v_range`](Self::v_range), clamped at the edges (see
    /// [`partial_u`](Self::partial_u)).
    pub fn partial_v(&self, u: f64, v: f64) -> Vector3<f64> {
        let (v_min, v_max) = self.v_range();
        let k = ((v_max - v_min) * 1e-4).max(1e-9);
        let v_lo = (v - k).max(v_min);
        let v_hi = (v + k).min(v_max);
        let d = v_hi - v_lo;
        if d.abs() < 1e-30 {
            return Vector3::zeros();
        }
        (self.evaluate(u, v_hi) - self.evaluate(u, v_lo)) / d
    }

    /// The **unit surface normal** at `(u, v)` — the normalised cross product
    /// `(∂S/∂u × ∂S/∂v)/|∂S/∂u × ∂S/∂v|` of the two tangent vectors
    /// ([`partial_u`](Self::partial_u), [`partial_v`](Self::partial_v)).
    ///
    /// Points to the side given by the right-hand rule of the (u, v)
    /// parameterisation — the outward face for a counter-clockwise patch. This
    /// is the vector shading, offsetting and ray intersection need. Returns the
    /// zero vector at a degenerate point where the tangents are parallel or
    /// vanish (`|∂S/∂u × ∂S/∂v| ≈ 0`, e.g. a pole of the parameterisation),
    /// where the normal is undefined.
    pub fn normal(&self, u: f64, v: f64) -> Vector3<f64> {
        let cross = self.partial_u(u, v).cross(&self.partial_v(u, v));
        let mag = cross.norm();
        if mag < 1e-12 {
            Vector3::zeros()
        } else {
            cross / mag
        }
    }

    /// The **surface area** of the whole patch — `∫∫ |∂S/∂u × ∂S/∂v| du dv` over
    /// the full [`u_range`](Self::u_range) × [`v_range`](Self::v_range). A
    /// convenience for [`area_between`](Self::area_between) over the entire
    /// domain.
    pub fn area(&self) -> f64 {
        let (u_min, u_max) = self.u_range();
        let (v_min, v_max) = self.v_range();
        self.area_between(u_min, u_max, v_min, v_max)
    }

    /// The **area of a sub-patch** `[u0, u1] × [v0, v1]` —
    /// `∫∫ |∂S/∂u × ∂S/∂v| du dv`, the true geometric area of that region of
    /// the surface (the cross-product magnitude is the local area-scaling of
    /// the parameterisation, so the result is parameterisation-independent).
    ///
    /// Computed by 2-D composite **Simpson's rule** on the area element
    /// `|∂S/∂u × ∂S/∂v|`, with the domain first split at every interior knot in
    /// each direction so no panel straddles a knot (where the tangents lose
    /// smoothness). Parameters are clamped to the valid ranges; a reversed or
    /// degenerate box (`u1 ≤ u0` or `v1 ≤ v0`) returns `0.0`.
    pub fn area_between(&self, u0: f64, u1: f64, v0: f64, v1: f64) -> f64 {
        const PANELS: usize = 16; // even, per knot-span cell, each direction
        let (u_min, u_max) = self.u_range();
        let (v_min, v_max) = self.v_range();
        let ua = u0.clamp(u_min, u_max);
        let ub = u1.clamp(u_min, u_max);
        let va = v0.clamp(v_min, v_max);
        let vb = v1.clamp(v_min, v_max);
        if ub <= ua || vb <= va {
            return 0.0;
        }
        let u_breaks = knot_breakpoints(&self.u_knots, ua, ub);
        let v_breaks = knot_breakpoints(&self.v_knots, va, vb);
        let mut total = 0.0;
        for us in u_breaks.windows(2) {
            for vs in v_breaks.windows(2) {
                total += self.area_cell(us[0], us[1], vs[0], vs[1], PANELS);
            }
        }
        total
    }

    /// 2-D composite-Simpson area element integral over one smooth cell.
    fn area_cell(&self, ua: f64, ub: f64, va: f64, vb: f64, n: usize) -> f64 {
        let hu = (ub - ua) / n as f64;
        let hv = (vb - va) / n as f64;
        // Composite-Simpson weight for node index `i` of `n` (even) panels.
        let weight = |i: usize| -> f64 {
            if i == 0 || i == n {
                1.0
            } else if i % 2 == 1 {
                4.0
            } else {
                2.0
            }
        };
        let mut sum = 0.0;
        for i in 0..=n {
            let u = ua + hu * i as f64;
            let wu = weight(i);
            for j in 0..=n {
                let v = va + hv * j as f64;
                let area_element = self.partial_u(u, v).cross(&self.partial_v(u, v)).norm();
                sum += wu * weight(j) * area_element;
            }
        }
        (hu / 3.0) * (hv / 3.0) * sum
    }
}

/// Breakpoints `[a, …interior knots in (a, b)…, b]` so each Simpson panel lies
/// within one smooth knot span.
fn knot_breakpoints(knots: &[f64], a: f64, b: f64) -> Vec<f64> {
    let eps = 1e-12;
    let mut breaks = vec![a];
    let mut prev = a;
    for &k in knots {
        if k > a + eps && k < b - eps && k > prev + eps {
            breaks.push(k);
            prev = k;
        }
    }
    breaks.push(b);
    breaks
}

#[cfg(test)]
mod tests {
    use super::*;

    fn open_uniform_knots(n_cp: usize, degree: usize) -> Vec<f64> {
        // Clamped uniform knot vector — multiplicity (degree+1) at
        // each endpoint, uniform interior.
        let p = degree;
        let m = n_cp + p + 1;
        let mut k = vec![0.0; m];
        let n_internal = n_cp - p - 1;
        for (i, kv) in k.iter_mut().enumerate().take(m) {
            if i <= p {
                *kv = 0.0;
            } else if i >= n_cp {
                *kv = 1.0;
            } else {
                let idx = i - p; // 1-based interior index
                *kv = idx as f64 / (n_internal + 1) as f64;
            }
        }
        k
    }

    #[test]
    fn rejects_bad_degrees() {
        let s = NurbsSurface::new(
            0,
            3,
            vec![0.0; 5],
            vec![0.0; 8],
            vec![vec![Vector3::zeros(); 4]; 4],
            vec![vec![1.0; 4]; 4],
        );
        assert_eq!(s.unwrap_err().code(), "surface.bad_degree");
    }

    #[test]
    fn rejects_non_rectangular_grid() {
        let s = NurbsSurface::new(
            3,
            3,
            open_uniform_knots(4, 3),
            open_uniform_knots(4, 3),
            vec![
                vec![Vector3::zeros(); 4],
                vec![Vector3::zeros(); 3], // wrong length
                vec![Vector3::zeros(); 4],
                vec![Vector3::zeros(); 4],
            ],
            vec![vec![1.0; 4]; 4],
        );
        assert_eq!(s.unwrap_err().code(), "surface.bad_knot_vector");
    }

    #[test]
    fn accepts_well_formed_4x4_cubic() {
        let s = NurbsSurface::new(
            3,
            3,
            open_uniform_knots(4, 3),
            open_uniform_knots(4, 3),
            vec![vec![Vector3::zeros(); 4]; 4],
            vec![vec![1.0; 4]; 4],
        );
        assert!(s.is_ok());
        let s = s.unwrap();
        assert_eq!(s.nu(), 4);
        assert_eq!(s.nv(), 4);
        assert_eq!(s.u_range(), (0.0, 1.0));
        assert_eq!(s.v_range(), (0.0, 1.0));
    }

    // ===== Phase 9B — surface evaluation tests =====

    /// Build a 4x4 planar bezier surface lying in z=0 with the CPs
    /// at the corners of a unit square (and 1/3, 2/3 sample points
    /// on the inner CPs for a perfect plane).
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

    #[test]
    fn planar_surface_at_corners_returns_corner_cps() {
        let s = planar_unit_square_surface();
        let p00 = s.evaluate(0.0, 0.0);
        let p10 = s.evaluate(1.0, 0.0);
        let p01 = s.evaluate(0.0, 1.0);
        let p11 = s.evaluate(1.0, 1.0);
        assert!((p00 - Vector3::new(0.0, 0.0, 0.0)).norm() < 1e-10);
        assert!((p10 - Vector3::new(1.0, 0.0, 0.0)).norm() < 1e-10);
        assert!((p01 - Vector3::new(0.0, 1.0, 0.0)).norm() < 1e-10);
        assert!((p11 - Vector3::new(1.0, 1.0, 0.0)).norm() < 1e-10);
    }

    #[test]
    fn planar_surface_centroid_is_centre_of_square() {
        let s = planar_unit_square_surface();
        let mid = s.evaluate(0.5, 0.5);
        assert!(
            (mid - Vector3::new(0.5, 0.5, 0.0)).norm() < 1e-10,
            "midpoint = {mid:?}"
        );
    }

    #[test]
    fn planar_surface_arbitrary_point_is_on_plane() {
        let s = planar_unit_square_surface();
        // Every point on the surface has z=0 and (x, y) inside the
        // unit square, with x ≈ u and y ≈ v for a perfectly planar
        // tensor-product Bezier.
        for &(u, v) in &[(0.1_f64, 0.7_f64), (0.3, 0.4), (0.9, 0.2)] {
            let p = s.evaluate(u, v);
            assert!((p.x - u).abs() < 1e-10, "p.x={} vs u={}", p.x, u);
            assert!((p.y - v).abs() < 1e-10, "p.y={} vs v={}", p.y, v);
            assert!(p.z.abs() < 1e-10);
        }
    }

    // ===== differential geometry: partials / normal / area =====

    /// A quarter cylinder of radius `r`, height `h`: the rational-quadratic
    /// quarter circle in xy (CPs (r,0),(r,r),(0,r), weights 1,√2/2,1) extruded
    /// linearly along z from 0 to `h`.
    fn quarter_cylinder(r: f64, h: f64) -> NurbsSurface {
        let w = std::f64::consts::FRAC_1_SQRT_2;
        // control_points[i (u)][j (v)] — u is the circle, v the z extrusion.
        let cps = vec![
            vec![Vector3::new(r, 0.0, 0.0), Vector3::new(r, 0.0, h)],
            vec![Vector3::new(r, r, 0.0), Vector3::new(r, r, h)],
            vec![Vector3::new(0.0, r, 0.0), Vector3::new(0.0, r, h)],
        ];
        let weights = vec![vec![1.0, 1.0], vec![w, w], vec![1.0, 1.0]];
        NurbsSurface::new(
            2,
            1,
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0],
            cps,
            weights,
        )
        .unwrap()
    }

    #[test]
    fn planar_surface_normal_is_constant_and_area_is_exact() {
        // The unit-square patch S(u,v) = (u, v, 0): tangents are the axes, the
        // normal is +z everywhere, and the area is exactly 1.
        let s = planar_unit_square_surface();
        for &(u, v) in &[(0.2_f64, 0.3_f64), (0.5, 0.5), (0.8, 0.1)] {
            let n = s.normal(u, v);
            assert!(
                (n - Vector3::new(0.0, 0.0, 1.0)).norm() < 1e-6,
                "normal {n:?} at ({u},{v})"
            );
            assert!((n.norm() - 1.0).abs() < 1e-9);
        }
        assert!((s.partial_u(0.5, 0.5) - Vector3::new(1.0, 0.0, 0.0)).norm() < 1e-6);
        assert!((s.partial_v(0.5, 0.5) - Vector3::new(0.0, 1.0, 0.0)).norm() < 1e-6);
        assert!((s.area() - 1.0).abs() < 1e-6, "planar area {}", s.area());
    }

    #[test]
    fn quarter_cylinder_normal_is_radial_and_area_is_lateral() {
        // GROUND TRUTH: a quarter cylinder of radius r and height h has a radial
        // outward unit normal (⊥ the z axis) and lateral area = (π/2)·r·h — the
        // quarter circumference times the height, independent of the NURBS
        // parameterization.
        let (r, h) = (2.0_f64, 3.0_f64);
        let s = quarter_cylinder(r, h);
        for &(u, v) in &[(0.2_f64, 0.3_f64), (0.5, 0.5), (0.8, 0.7)] {
            let p = s.evaluate(u, v);
            // The surface lies on the cylinder x²+y²=r².
            assert!(
                ((p.x * p.x + p.y * p.y) - r * r).abs() < 1e-9,
                "off cylinder at ({u},{v})"
            );
            let n = s.normal(u, v);
            assert!((n.norm() - 1.0).abs() < 1e-9, "unit normal at ({u},{v})");
            assert!(n.z.abs() < 1e-6, "normal ⊥ axis, n.z={}", n.z);
            // Radial outward: parallel to the point's xy projection.
            let radial = Vector3::new(p.x, p.y, 0.0).normalize();
            assert!(
                (n - radial).norm() < 1e-3,
                "normal {n:?} vs radial {radial:?}"
            );
        }
        let expected = std::f64::consts::FRAC_PI_2 * r * h; // (π/2)·r·h
        assert!(
            (s.area() - expected).abs() / expected < 1e-4,
            "cylinder lateral area {} != (π/2)·r·h = {expected}",
            s.area()
        );
    }

    #[test]
    fn surface_area_between_is_additive_and_clamped() {
        let s = planar_unit_square_surface();
        let whole = s.area(); // 1.0
        let left = s.area_between(0.0, 0.5, 0.0, 1.0);
        let right = s.area_between(0.5, 1.0, 0.0, 1.0);
        assert!(
            (left + right - whole).abs() < 1e-6,
            "{left} + {right} != {whole}"
        );
        assert!((left - 0.5).abs() < 1e-6, "left-half area {left}");
        // Degenerate and reversed boxes return 0.
        assert_eq!(s.area_between(0.5, 0.5, 0.0, 1.0), 0.0);
        assert_eq!(s.area_between(1.0, 0.0, 0.0, 1.0), 0.0);
        assert_eq!(s.area_between(0.0, 1.0, 1.0, 0.0), 0.0);
    }
}
