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
}
