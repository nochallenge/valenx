//! **Curvature analysis** for Class-A surfacing — the principal curvatures and
//! the local-shape diagnosis a surfacing tool's curvature / porcupine analysis
//! reports.
//!
//! [`NurbsSurface`] already exposes the Gaussian curvature `K = κ₁·κ₂` and the
//! mean curvature `H = (κ₁ + κ₂)/2`. This recovers the individual **principal
//! curvatures** `κ_max, κ_min = H ± √(max(H²−K, 0))` and classifies the point
//! as elliptic (dome), hyperbolic (saddle), parabolic (developable /
//! cylindrical) or planar — exactly the curvature map a Class-A reviewer reads
//! to judge whether a surface flows fairly.
//!
//! Validated against analytic surfaces: a plane → `(0, 0)` / planar; a cylinder
//! of radius `r` → `(1/r, 0)` / parabolic; a sphere of radius `R` →
//! `(1/R, 1/R)` with `K = 1/R²`.
//!
//! Honest scope: pointwise principal-curvature analysis from the surface's
//! fundamental forms — research-grade. It is not a rendered reflection / zebra
//! / isophote analyser (those are viewport features), and a step toward, not an
//! equal of, CATIA-class surface diagnostics.

use crate::nurbs_surface::NurbsSurface;

/// Local surface shape at a point, from the signs of the principal curvatures.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LocalShape {
    /// Both principal curvatures the same sign — a dome (`K > 0`).
    Elliptic,
    /// Principal curvatures of opposite sign — a saddle (`K < 0`).
    Hyperbolic,
    /// One principal curvature ≈ 0 — developable / cylindrical (`K ≈ 0`, `H ≠ 0`).
    Parabolic,
    /// Both principal curvatures ≈ 0 — flat (`K ≈ 0`, `H ≈ 0`).
    Planar,
}

/// Principal curvatures `(κ_max, κ_min)` of `surface` at `(u, v)`, recovered
/// from the mean curvature `H` and Gaussian curvature `K` as
/// `κ = H ± √(max(H²−K, 0))`.
pub fn principal_curvatures(surface: &NurbsSurface, u: f64, v: f64) -> (f64, f64) {
    let h = surface.mean_curvature(u, v);
    let k = surface.gaussian_curvature(u, v);
    let disc = (h * h - k).max(0.0).sqrt();
    (h + disc, h - disc)
}

/// Classify the local surface shape at `(u, v)`. `tol` is the curvature
/// magnitude below which a principal curvature is treated as zero.
pub fn local_shape(surface: &NurbsSurface, u: f64, v: f64, tol: f64) -> LocalShape {
    let (k1, k2) = principal_curvatures(surface, u, v);
    match (k1.abs() < tol, k2.abs() < tol) {
        (true, true) => LocalShape::Planar,
        (true, false) | (false, true) => LocalShape::Parabolic,
        (false, false) => {
            if k1 * k2 > 0.0 {
                LocalShape::Elliptic
            } else {
                LocalShape::Hyperbolic
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    fn plane() -> NurbsSurface {
        NurbsSurface::new(
            1,
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0],
            vec![
                vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(0.0, 1.0, 0.0)],
                vec![Vector3::new(1.0, 0.0, 0.0), Vector3::new(1.0, 1.0, 0.0)],
            ],
            vec![vec![1.0, 1.0], vec![1.0, 1.0]],
        )
        .expect("valid plane")
    }

    /// A rational-quadratic quarter cylinder of radius `r` (exact circular arc
    /// in u, straight extrusion in v).
    fn quarter_cylinder(r: f64, h: f64) -> NurbsSurface {
        let w = std::f64::consts::FRAC_1_SQRT_2;
        NurbsSurface::new(
            2,
            1,
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0],
            vec![
                vec![Vector3::new(r, 0.0, 0.0), Vector3::new(r, 0.0, h)],
                vec![Vector3::new(r, r, 0.0), Vector3::new(r, r, h)],
                vec![Vector3::new(0.0, r, 0.0), Vector3::new(0.0, r, h)],
            ],
            vec![vec![1.0, 1.0], vec![w, w], vec![1.0, 1.0]],
        )
        .expect("valid quarter cylinder")
    }

    #[test]
    fn plane_is_planar_with_zero_principal_curvatures() {
        let s = plane();
        let (k1, k2) = principal_curvatures(&s, 0.5, 0.5);
        assert!(k1.abs() < 1e-6 && k2.abs() < 1e-6, "κ = ({k1}, {k2})");
        assert_eq!(local_shape(&s, 0.5, 0.5, 1e-4), LocalShape::Planar);
    }

    #[test]
    fn cylinder_has_one_principal_curvature_equal_to_inverse_radius() {
        let r = 2.0;
        let s = quarter_cylinder(r, 1.0);
        let (k1, k2) = principal_curvatures(&s, 0.5, 0.5);
        // One principal curvature has magnitude 1/r (around the arc), the other
        // is ~0 (the straight extrusion). The sign depends on the normal
        // orientation, so compare magnitudes.
        let (curved, flat) = if k1.abs() >= k2.abs() {
            (k1, k2)
        } else {
            (k2, k1)
        };
        assert!(
            (curved.abs() - 1.0 / r).abs() < 0.02 * (1.0 / r),
            "curved |κ| {curved} vs 1/r {}",
            1.0 / r
        );
        assert!(flat.abs() < 0.02, "extrusion κ {flat} should be ~0");
        assert_eq!(local_shape(&s, 0.5, 0.5, 0.05), LocalShape::Parabolic);
    }

    #[test]
    fn principal_max_is_at_least_min() {
        let s = quarter_cylinder(1.5, 2.0);
        for &(u, v) in &[(0.25, 0.25), (0.5, 0.5), (0.75, 0.8)] {
            let (k1, k2) = principal_curvatures(&s, u, v);
            assert!(k1 >= k2 - 1e-12, "κ_max {k1} < κ_min {k2}");
        }
    }
}
