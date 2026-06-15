//! Rectangular 0/45/90 rosette reduction and principal-strain analysis.
//!
//! A *rectangular* (also called a 45-degree) rosette bonds three
//! single-element strain gauges to a surface at 0, 45, and 90 degrees
//! measured from the chosen x-axis. Each gauge reports the normal
//! strain along its own axis. Inverting the strain-transformation law
//! recovers the full in-plane (Cartesian) strain state.
//!
//! # Strain-transformation law
//!
//! The normal strain along a line at angle `theta` from the x-axis is
//!
//! ```text
//! eps(theta) = eps_x cos^2(theta) + eps_y sin^2(theta)
//!            + gamma_xy sin(theta) cos(theta)
//! ```
//!
//! Substituting the three gauge angles gives the closed-form inverse
//! used by [`reduce`]:
//!
//! ```text
//! eps_x    = eps_0
//! eps_y    = eps_90
//! gamma_xy = 2 eps_45 - eps_0 - eps_90
//! ```
//!
//! # Principal strains
//!
//! The in-plane strain tensor is the symmetric matrix
//!
//! ```text
//! [ eps_x      gamma_xy/2 ]
//! [ gamma_xy/2 eps_y      ]
//! ```
//!
//! whose eigenvalues are the principal strains
//!
//! ```text
//! eps_1,2 = (eps_x + eps_y)/2 +/- sqrt( ((eps_x - eps_y)/2)^2
//!                                       + (gamma_xy/2)^2 )
//! ```
//!
//! and whose eigenvectors point along the principal directions, the
//! larger-eigenvalue direction making an angle
//! `theta_p = 0.5 atan2(gamma_xy, eps_x - eps_y)` with the x-axis.

use nalgebra::Matrix2;
use serde::{Deserialize, Serialize};

/// The three normal-strain readings of a rectangular 0/45/90 rosette.
///
/// Strains are stored as dimensionless engineering strain (for example
/// `0.001` for 1000 microstrain). The sign convention is tension
/// positive.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RosetteReadings {
    /// Gauge bonded along the x-axis (0 degrees).
    pub eps_0: f64,
    /// Gauge bonded at 45 degrees from the x-axis.
    pub eps_45: f64,
    /// Gauge bonded along the y-axis (90 degrees).
    pub eps_90: f64,
}

impl RosetteReadings {
    /// Construct a reading triple.
    pub fn new(eps_0: f64, eps_45: f64, eps_90: f64) -> Self {
        Self {
            eps_0,
            eps_45,
            eps_90,
        }
    }
}

/// The recovered in-plane (Cartesian) strain state.
///
/// `gamma_xy` is the *engineering* shear strain (twice the tensor shear
/// component), matching the convention used in every introductory
/// strength-of-materials text.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CartesianStrain {
    /// Direct strain along x.
    pub eps_x: f64,
    /// Direct strain along y.
    pub eps_y: f64,
    /// Engineering shear strain in the xy-plane.
    pub gamma_xy: f64,
}

impl CartesianStrain {
    /// Construct a Cartesian strain state directly.
    pub fn new(eps_x: f64, eps_y: f64, gamma_xy: f64) -> Self {
        Self {
            eps_x,
            eps_y,
            gamma_xy,
        }
    }

    /// The symmetric 2x2 strain tensor
    /// `[[eps_x, gamma_xy/2], [gamma_xy/2, eps_y]]`.
    ///
    /// Note the off-diagonal entry is the *tensor* shear
    /// (`gamma_xy / 2`), not the engineering shear, so that the matrix
    /// transforms as a true second-order tensor and its eigenvalues are
    /// the principal strains.
    pub fn tensor(&self) -> Matrix2<f64> {
        let t = 0.5 * self.gamma_xy;
        Matrix2::new(self.eps_x, t, t, self.eps_y)
    }

    /// Normal strain along a line at `theta` radians from the x-axis,
    /// via the forward strain-transformation law.
    ///
    /// This is the inverse of [`reduce`]: feeding the recovered state at
    /// 0, 45, and 90 degrees reproduces the original gauge readings.
    pub fn normal_strain_at(&self, theta: f64) -> f64 {
        let (s, c) = theta.sin_cos();
        self.eps_x * c * c + self.eps_y * s * s + self.gamma_xy * s * c
    }
}

/// Principal-strain solution: the two principal strains and the
/// orientation of the major principal axis.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PrincipalStrain {
    /// Major (algebraically larger) principal strain.
    pub eps_1: f64,
    /// Minor (algebraically smaller) principal strain.
    pub eps_2: f64,
    /// Angle (radians) from the x-axis to the major principal
    /// direction, in the range `(-pi/2, pi/2]`.
    pub theta_p: f64,
}

impl PrincipalStrain {
    /// Maximum in-plane shear strain, `gamma_max = eps_1 - eps_2`.
    ///
    /// This is the engineering shear strain on the plane of maximum
    /// distortion, oriented 45 degrees from the principal axes.
    pub fn max_shear(&self) -> f64 {
        self.eps_1 - self.eps_2
    }
}

/// Reduce three rectangular 0/45/90 gauge readings to the Cartesian
/// strain state.
///
/// Applies the closed-form inverse of the strain-transformation law:
///
/// ```text
/// eps_x    = eps_0
/// eps_y    = eps_90
/// gamma_xy = 2 eps_45 - eps_0 - eps_90
/// ```
///
/// # Examples
///
/// ```
/// use valenx_strainrosette::{reduce, RosetteReadings};
///
/// let s = reduce(RosetteReadings::new(0.0010, 0.0006, 0.0002));
/// assert!((s.eps_x - 0.0010).abs() < 1e-12);
/// assert!((s.eps_y - 0.0002).abs() < 1e-12);
/// // gamma_xy = 2 * 0.0006 - 0.0010 - 0.0002 = 0.0
/// assert!(s.gamma_xy.abs() < 1e-12);
/// ```
pub fn reduce(r: RosetteReadings) -> CartesianStrain {
    CartesianStrain {
        eps_x: r.eps_0,
        eps_y: r.eps_90,
        gamma_xy: 2.0 * r.eps_45 - r.eps_0 - r.eps_90,
    }
}

/// Compute principal strains and the major-axis orientation from a
/// Cartesian strain state.
///
/// The principal strains are the closed-form eigenvalues of the strain
/// tensor; the angle uses the four-quadrant `atan2` form so the
/// reported `theta_p` always names the *major* (larger-eigenvalue)
/// direction.
///
/// # Examples
///
/// ```
/// use valenx_strainrosette::{principal_strains, CartesianStrain};
///
/// // Pure uniaxial strain along x: principals are eps_x and 0,
/// // aligned with the axes.
/// let p = principal_strains(CartesianStrain::new(0.002, 0.0, 0.0));
/// assert!((p.eps_1 - 0.002).abs() < 1e-12);
/// assert!(p.eps_2.abs() < 1e-12);
/// assert!(p.theta_p.abs() < 1e-12);
/// ```
pub fn principal_strains(s: CartesianStrain) -> PrincipalStrain {
    let mean = 0.5 * (s.eps_x + s.eps_y);
    let half_diff = 0.5 * (s.eps_x - s.eps_y);
    let tensor_shear = 0.5 * s.gamma_xy;
    let radius = (half_diff * half_diff + tensor_shear * tensor_shear).sqrt();

    let eps_1 = mean + radius;
    let eps_2 = mean - radius;

    // Angle to the MAJOR principal axis. atan2(gamma_xy, eps_x - eps_y)
    // gives 2 * theta in (-pi, pi]; halving lands theta_p in
    // (-pi/2, pi/2]. When the state is isotropic (radius == 0) the
    // direction is undefined, so report 0.
    let theta_p = if radius == 0.0 {
        0.0
    } else {
        0.5 * s.gamma_xy.atan2(s.eps_x - s.eps_y)
    };

    PrincipalStrain {
        eps_1,
        eps_2,
        theta_p,
    }
}

/// One-shot convenience: reduce gauge readings and immediately solve the
/// principal strains.
///
/// Equivalent to `principal_strains(reduce(r))`.
pub fn principal_from_readings(r: RosetteReadings) -> PrincipalStrain {
    principal_strains(reduce(r))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

    const EPS: f64 = 1e-12;

    #[test]
    fn reduction_recovers_eps_x_and_eps_y_directly() {
        // eps_x must equal the 0-degree gauge, eps_y the 90-degree gauge.
        let r = RosetteReadings::new(7.0e-4, -3.0e-4, 2.0e-4);
        let s = reduce(r);
        assert!((s.eps_x - 7.0e-4).abs() < EPS, "eps_x = {x}", x = s.eps_x);
        assert!((s.eps_y - 2.0e-4).abs() < EPS, "eps_y = {y}", y = s.eps_y);
    }

    #[test]
    fn gamma_formula_matches_textbook() {
        // gamma_xy = 2*eps_45 - eps_0 - eps_90.
        let r = RosetteReadings::new(1.0e-3, 4.0e-4, -2.0e-4);
        let s = reduce(r);
        let expected = 2.0 * 4.0e-4 - 1.0e-3 - (-2.0e-4);
        assert!(
            (s.gamma_xy - expected).abs() < EPS,
            "gamma_xy = {g}, expected = {expected}",
            g = s.gamma_xy
        );
    }

    #[test]
    fn forward_transform_round_trips_through_the_three_gauges() {
        // Feeding the recovered state back through the forward law at
        // 0/45/90 must reproduce the original readings exactly.
        let r = RosetteReadings::new(9.0e-4, 1.0e-4, -5.0e-4);
        let s = reduce(r);
        assert!((s.normal_strain_at(0.0) - r.eps_0).abs() < EPS);
        assert!((s.normal_strain_at(FRAC_PI_4) - r.eps_45).abs() < EPS);
        assert!((s.normal_strain_at(FRAC_PI_2) - r.eps_90).abs() < EPS);
    }

    #[test]
    fn principal_strains_match_closed_form_invariants() {
        // The sum of principals equals the trace (eps_x + eps_y); the
        // difference equals 2*radius. Check both invariants on a known
        // state.
        let s = CartesianStrain::new(6.0e-4, 2.0e-4, 8.0e-4);
        let p = principal_strains(s);

        let trace = s.eps_x + s.eps_y;
        assert!(
            ((p.eps_1 + p.eps_2) - trace).abs() < EPS,
            "sum {sum} != trace {trace}",
            sum = p.eps_1 + p.eps_2
        );

        let half_diff = 0.5 * (s.eps_x - s.eps_y);
        let tensor_shear = 0.5 * s.gamma_xy;
        let radius = (half_diff * half_diff + tensor_shear * tensor_shear).sqrt();
        assert!((p.eps_1 - (0.5 * trace + radius)).abs() < EPS);
        assert!((p.eps_2 - (0.5 * trace - radius)).abs() < EPS);
        assert!(p.eps_1 >= p.eps_2);
    }

    #[test]
    fn principal_strains_equal_nalgebra_eigenvalues() {
        // Cross-validate the closed-form principals against nalgebra's
        // symmetric eigensolver on the strain tensor.
        let s = CartesianStrain::new(5.0e-4, -1.0e-4, 9.0e-4);
        let p = principal_strains(s);

        let eig = s.tensor().symmetric_eigen();
        let mut vals = [eig.eigenvalues[0], eig.eigenvalues[1]];
        vals.sort_by(|a, b| b.partial_cmp(a).unwrap());

        assert!(
            (p.eps_1 - vals[0]).abs() < 1e-10,
            "eps_1 {e1} vs eig {v}",
            e1 = p.eps_1,
            v = vals[0]
        );
        assert!(
            (p.eps_2 - vals[1]).abs() < 1e-10,
            "eps_2 {e2} vs eig {v}",
            e2 = p.eps_2,
            v = vals[1]
        );
    }

    #[test]
    fn principal_angle_diagonalises_the_tensor() {
        // Rotating the strain tensor by -theta_p must zero the shear and
        // place eps_1 on the leading diagonal. This is the defining
        // property of the principal direction.
        let s = CartesianStrain::new(6.0e-4, 2.0e-4, 8.0e-4);
        let p = principal_strains(s);

        // Normal strain along theta_p should equal eps_1.
        let along = s.normal_strain_at(p.theta_p);
        assert!(
            (along - p.eps_1).abs() < 1e-12,
            "normal at theta_p = {along}, eps_1 = {e1}",
            e1 = p.eps_1
        );
        // Normal strain perpendicular should equal eps_2.
        let perp = s.normal_strain_at(p.theta_p + FRAC_PI_2);
        assert!((perp - p.eps_2).abs() < 1e-12, "perp = {perp}");
    }

    #[test]
    fn principal_angle_is_in_canonical_range() {
        // theta_p must always land in (-pi/2, pi/2].
        let cases = [
            CartesianStrain::new(3.0e-4, 9.0e-4, 5.0e-4),
            CartesianStrain::new(9.0e-4, 3.0e-4, -5.0e-4),
            CartesianStrain::new(-2.0e-4, 4.0e-4, 7.0e-4),
        ];
        for s in cases {
            let p = principal_strains(s);
            assert!(
                p.theta_p > -FRAC_PI_2 - EPS && p.theta_p <= FRAC_PI_2 + EPS,
                "theta_p out of range: {t}",
                t = p.theta_p
            );
        }
    }

    #[test]
    fn uniaxial_state_recovers_applied_strain_and_axis() {
        // A pure uniaxial strain of 0.0012 along x: gauge 0 reads it
        // fully, gauge 90 reads -nu * eps (but for the kinematic check
        // we use eps_y = 0), gauge 45 reads half. Principal 1 must equal
        // the applied strain, principal 2 zero, and the major axis lie
        // along x (theta_p == 0).
        let applied = 1.2e-3;
        let r = RosetteReadings::new(applied, 0.5 * applied, 0.0);
        let s = reduce(r);
        assert!((s.eps_x - applied).abs() < EPS);
        assert!(s.eps_y.abs() < EPS);
        // gamma_xy = 2*(applied/2) - applied - 0 = 0 for this aligned case.
        assert!(s.gamma_xy.abs() < EPS, "gamma_xy = {g}", g = s.gamma_xy);

        let p = principal_strains(s);
        assert!((p.eps_1 - applied).abs() < EPS, "eps_1 = {e}", e = p.eps_1);
        assert!(p.eps_2.abs() < EPS, "eps_2 = {e}", e = p.eps_2);
        assert!(p.theta_p.abs() < EPS, "theta_p = {t}", t = p.theta_p);
    }

    #[test]
    fn pure_shear_state_has_equal_and_opposite_principals_at_45_deg() {
        // eps_x = eps_y = 0, gamma_xy = g => principals +/- g/2 with the
        // major axis at 45 degrees. A classic textbook Mohr's-circle
        // case.
        let g = 8.0e-4;
        let s = CartesianStrain::new(0.0, 0.0, g);
        let p = principal_strains(s);
        assert!((p.eps_1 - 0.5 * g).abs() < EPS, "eps_1 = {e}", e = p.eps_1);
        assert!((p.eps_2 + 0.5 * g).abs() < EPS, "eps_2 = {e}", e = p.eps_2);
        assert!(
            (p.theta_p - FRAC_PI_4).abs() < EPS,
            "theta_p = {t}",
            t = p.theta_p
        );
        assert!(
            (p.max_shear() - g).abs() < EPS,
            "gamma_max = {m}",
            m = p.max_shear()
        );
    }

    #[test]
    fn isotropic_state_has_undefined_axis_reported_as_zero() {
        // eps_x = eps_y, gamma_xy = 0 => the Mohr circle degenerates to a
        // point; every direction is principal, so theta_p is reported as 0.
        let s = CartesianStrain::new(3.0e-4, 3.0e-4, 0.0);
        let p = principal_strains(s);
        assert!((p.eps_1 - 3.0e-4).abs() < EPS);
        assert!((p.eps_2 - 3.0e-4).abs() < EPS);
        assert!(p.theta_p.abs() < EPS);
        assert!(p.max_shear().abs() < EPS);
    }

    #[test]
    fn principal_angle_full_atan2_quadrant_is_used() {
        // With eps_x < eps_y and positive shear the naive
        // 0.5*atan(gamma/(eps_x-eps_y)) form would mis-place the axis;
        // the atan2 form keeps theta_p naming the major direction. Verify
        // by direct diagonalisation.
        let s = CartesianStrain::new(1.0e-4, 9.0e-4, 6.0e-4);
        let p = principal_strains(s);
        let along = s.normal_strain_at(p.theta_p);
        assert!(
            (along - p.eps_1).abs() < 1e-12,
            "major axis check failed: {along} vs {e1}",
            e1 = p.eps_1
        );
        // And theta_p stays within the canonical half-open interval.
        assert!(p.theta_p > -FRAC_PI_2 && p.theta_p <= FRAC_PI_2);
        // Sanity: PI constant referenced so the range bound is explicit.
        assert!(2.0 * p.theta_p > -PI && 2.0 * p.theta_p <= PI);
    }
}
