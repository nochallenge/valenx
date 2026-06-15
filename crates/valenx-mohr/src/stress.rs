//! Plane-stress state and the closed-form Mohr's-circle transformation.
//!
//! Sign and angle conventions follow the standard mechanics-of-materials
//! presentation (e.g. Gere & Goodno, Hibbeler):
//!
//! - Normal stresses are positive in tension.
//! - The shear stress `txy` is positive when it acts in the `+y`
//!   direction on the face whose outward normal points in `+x`.
//! - The plane angle `theta` is measured counter-clockwise from the `x`
//!   axis to the outward normal of the cut plane.
//!
//! For a plane whose normal is rotated by `theta`, the transformed
//! normal and shear stresses are the standard double-angle relations
//!
//! ```text
//! sx'  = (sx+sy)/2 + (sx-sy)/2 * cos(2 theta) + txy * sin(2 theta)
//! txy' = -(sx-sy)/2 * sin(2 theta) + txy * cos(2 theta)
//! ```

use serde::{Deserialize, Serialize};

use crate::error::MohrError;

/// An in-plane stress state described by its two normal components and
/// the in-plane shear component.
///
/// Construct one with [`StressState::new`], which validates that every
/// component is finite. The transformation methods then expose the
/// classic Mohr's-circle results.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct StressState {
    /// Normal stress on the `x` face (positive in tension).
    pub sx: f64,
    /// Normal stress on the `y` face (positive in tension).
    pub sy: f64,
    /// In-plane shear stress `txy`.
    pub txy: f64,
}

/// The two principal stresses of a plane-stress state.
///
/// By convention `s1 >= s2`. These are the maximum and minimum normal
/// stresses over all plane orientations, and they coincide with the
/// eigenvalues of the symmetric 2x2 stress tensor
/// `[[sx, txy], [txy, sy]]`.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PrincipalStresses {
    /// Maximum principal stress (algebraically largest).
    pub s1: f64,
    /// Minimum principal stress (algebraically smallest).
    pub s2: f64,
}

/// The normal and shear stress acting on a single inclined plane.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PlaneStress {
    /// Normal stress on the plane (positive in tension).
    pub normal: f64,
    /// Shear stress on the plane.
    pub shear: f64,
}

impl StressState {
    /// Construct a plane-stress state, validating that every component
    /// is finite.
    ///
    /// # Errors
    ///
    /// Returns [`MohrError::NonFinite`] if any of `sx`, `sy`, or `txy`
    /// is `NaN` or infinite.
    pub fn new(sx: f64, sy: f64, txy: f64) -> Result<Self, MohrError> {
        check_finite("sx", sx)?;
        check_finite("sy", sy)?;
        check_finite("txy", txy)?;
        Ok(Self { sx, sy, txy })
    }

    /// The average normal stress `(sx + sy) / 2`.
    ///
    /// This is the location of the centre of Mohr's circle on the
    /// normal-stress axis, and it is invariant under rotation of the
    /// plane (it equals `(s1 + s2) / 2`).
    pub fn mean_normal(&self) -> f64 {
        0.5 * (self.sx + self.sy)
    }

    /// The radius of Mohr's circle,
    /// `sqrt(((sx - sy) / 2)^2 + txy^2)`.
    ///
    /// This equals the maximum in-plane shear stress and half the
    /// difference of the principal stresses.
    pub fn radius(&self) -> f64 {
        let h = 0.5 * (self.sx - self.sy);
        (h * h + self.txy * self.txy).sqrt()
    }

    /// The two principal stresses `s1 >= s2`.
    ///
    /// ```text
    /// s1,2 = (sx+sy)/2 +/- sqrt(((sx-sy)/2)^2 + txy^2)
    /// ```
    pub fn principal_stresses(&self) -> PrincipalStresses {
        let c = self.mean_normal();
        let r = self.radius();
        PrincipalStresses {
            s1: c + r,
            s2: c - r,
        }
    }

    /// The maximum in-plane shear stress, `(s1 - s2) / 2`.
    ///
    /// Numerically identical to [`StressState::radius`]; provided under
    /// its mechanical name for readability.
    pub fn max_shear(&self) -> f64 {
        self.radius()
    }

    /// The principal angle `theta_p` (radians), measured
    /// counter-clockwise from the `x` axis to the normal of the plane
    /// carrying the maximum principal stress `s1`.
    ///
    /// Defined by `tan(2 theta_p) = 2 txy / (sx - sy)`. The
    /// quadrant-correct two-argument arctangent is used and the result
    /// is halved, giving a value in `(-pi/2, pi/2]`. For a hydrostatic
    /// state (`sx == sy` and `txy == 0`) every direction is principal;
    /// `0.0` is returned.
    pub fn principal_angle(&self) -> f64 {
        let num = 2.0 * self.txy;
        let den = self.sx - self.sy;
        if num == 0.0 && den == 0.0 {
            return 0.0;
        }
        0.5 * num.atan2(den)
    }

    /// The normal and shear stress on the plane whose outward normal is
    /// rotated `theta` radians counter-clockwise from the `x` axis.
    ///
    /// At `theta == 0` this returns the original face stresses
    /// (`normal == sx`, `shear == txy`).
    ///
    /// # Errors
    ///
    /// Returns [`MohrError::NonFiniteAngle`] if `theta` is `NaN` or
    /// infinite.
    pub fn stress_on_plane(&self, theta: f64) -> Result<PlaneStress, MohrError> {
        if !theta.is_finite() {
            return Err(MohrError::NonFiniteAngle {
                name: "theta",
                value: theta,
            });
        }
        let c = self.mean_normal();
        let h = 0.5 * (self.sx - self.sy);
        let (s2t, c2t) = (2.0 * theta).sin_cos();
        Ok(PlaneStress {
            normal: c + h * c2t + self.txy * s2t,
            shear: -h * s2t + self.txy * c2t,
        })
    }
}

/// Reject a non-finite stress component at the constructor boundary.
fn check_finite(name: &'static str, value: f64) -> Result<(), MohrError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(MohrError::NonFinite { name, value })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tolerance for analytic float comparisons.
    const EPS: f64 = 1e-9;

    #[test]
    fn rejects_non_finite_components() {
        assert!(StressState::new(f64::NAN, 0.0, 0.0).is_err());
        assert!(StressState::new(0.0, f64::INFINITY, 0.0).is_err());
        assert!(StressState::new(0.0, 0.0, f64::NEG_INFINITY).is_err());
        let err = StressState::new(0.0, f64::NAN, 0.0).unwrap_err();
        assert_eq!(err.code(), "mohr.non_finite");
    }

    #[test]
    fn rejects_non_finite_angle() {
        let s = StressState::new(10.0, 0.0, 0.0).unwrap();
        assert!(s.stress_on_plane(f64::NAN).is_err());
        assert!(s.stress_on_plane(f64::INFINITY).is_err());
        let err = s.stress_on_plane(f64::NAN).unwrap_err();
        assert_eq!(err.code(), "mohr.non_finite_angle");
    }

    #[test]
    fn principals_known_value() {
        // Classic worked example (Hibbeler): sx=-20, sy=90, txy=60.
        // s1 = 35 + sqrt(55^2 + 60^2) = 35 + 81.394... = 116.39 MPa
        // s2 = 35 - 81.394... = -46.39 MPa
        let s = StressState::new(-20.0, 90.0, 60.0).unwrap();
        let p = s.principal_stresses();
        let expected = (55.0_f64 * 55.0 + 60.0 * 60.0).sqrt();
        assert!((p.s1 - (35.0 + expected)).abs() < EPS, "s1 was {}", p.s1);
        assert!((p.s2 - (35.0 - expected)).abs() < EPS, "s2 was {}", p.s2);
        // Ordering guarantee.
        assert!(p.s1 >= p.s2);
    }

    #[test]
    fn sum_of_normals_is_invariant() {
        // s1 + s2 == sx + sy for arbitrary states (trace invariance).
        let cases = [
            (10.0, 4.0, 3.0),
            (-20.0, 90.0, 60.0),
            (0.0, 0.0, 7.5),
            (-5.0, -12.0, -2.0),
        ];
        for (sx, sy, txy) in cases {
            let s = StressState::new(sx, sy, txy).unwrap();
            let p = s.principal_stresses();
            assert!(
                (p.s1 + p.s2 - (sx + sy)).abs() < EPS,
                "invariant failed for ({sx}, {sy}, {txy}): s1+s2 = {}",
                p.s1 + p.s2
            );
        }
    }

    #[test]
    fn principals_are_eigenvalues_of_the_tensor() {
        // The symmetric tensor [[sx, txy], [txy, sy]] has eigenvalues
        // (tr/2) +/- sqrt((tr/2)^2 - det). Compare to principal_stresses.
        let cases = [(10.0, 4.0, 3.0), (-20.0, 90.0, 60.0), (7.0, -3.0, -4.0)];
        for (sx, sy, txy) in cases {
            let s = StressState::new(sx, sy, txy).unwrap();
            let p = s.principal_stresses();
            let tr = sx + sy;
            let det = sx * sy - txy * txy;
            let disc = ((tr * 0.5) * (tr * 0.5) - det).sqrt();
            let lam_max = tr * 0.5 + disc;
            let lam_min = tr * 0.5 - disc;
            assert!(
                (p.s1 - lam_max).abs() < EPS,
                "({sx},{sy},{txy}) s1={} eig_max={lam_max}",
                p.s1
            );
            assert!(
                (p.s2 - lam_min).abs() < EPS,
                "({sx},{sy},{txy}) s2={} eig_min={lam_min}",
                p.s2
            );
            // Eigenvalue residual: det(A - lambda I) == 0 to tolerance.
            for lam in [p.s1, p.s2] {
                let residual = (sx - lam) * (sy - lam) - txy * txy;
                assert!(
                    residual.abs() < 1e-6,
                    "residual {residual} for lambda {lam}"
                );
            }
        }
    }

    #[test]
    fn max_shear_equals_half_principal_difference() {
        let s = StressState::new(-20.0, 90.0, 60.0).unwrap();
        let p = s.principal_stresses();
        assert!((s.max_shear() - 0.5 * (p.s1 - p.s2)).abs() < EPS);
    }

    #[test]
    fn shear_vanishes_on_principal_plane() {
        // Rotating to the principal angle must zero the shear, and the
        // normal stress there must equal a principal stress.
        let s = StressState::new(-20.0, 90.0, 60.0).unwrap();
        let theta_p = s.principal_angle();
        let on = s.stress_on_plane(theta_p).unwrap();
        assert!(
            on.shear.abs() < 1e-6,
            "shear on principal plane = {}",
            on.shear
        );
        let p = s.principal_stresses();
        let hits_principal = (on.normal - p.s1).abs() < 1e-6 || (on.normal - p.s2).abs() < 1e-6;
        assert!(hits_principal, "normal {} is neither principal", on.normal);
        // The perpendicular plane is also shear-free.
        let on_perp = s
            .stress_on_plane(theta_p + std::f64::consts::FRAC_PI_2)
            .unwrap();
        assert!(on_perp.shear.abs() < 1e-6);
    }

    #[test]
    fn transform_at_zero_returns_original_face() {
        let s = StressState::new(40.0, -10.0, 15.0).unwrap();
        let on = s.stress_on_plane(0.0).unwrap();
        assert!((on.normal - 40.0).abs() < EPS, "normal {}", on.normal);
        assert!((on.shear - 15.0).abs() < EPS, "shear {}", on.shear);
    }

    #[test]
    fn transform_at_ninety_swaps_normals_and_flips_shear() {
        // A 90-degree rotation maps the x face onto the y face: the
        // normal becomes sy and the shear flips sign.
        let s = StressState::new(40.0, -10.0, 15.0).unwrap();
        let on = s.stress_on_plane(std::f64::consts::FRAC_PI_2).unwrap();
        assert!((on.normal - (-10.0)).abs() < 1e-6, "normal {}", on.normal);
        assert!((on.shear - (-15.0)).abs() < 1e-6, "shear {}", on.shear);
    }

    #[test]
    fn transformed_normals_sum_to_invariant_at_any_angle() {
        // sx'(theta) + sx'(theta + 90deg) == sx + sy for any theta.
        let s = StressState::new(40.0, -10.0, 15.0).unwrap();
        let theta = 0.37_f64;
        let a = s.stress_on_plane(theta).unwrap();
        let b = s
            .stress_on_plane(theta + std::f64::consts::FRAC_PI_2)
            .unwrap();
        assert!((a.normal + b.normal - (s.sx + s.sy)).abs() < 1e-9);
    }

    #[test]
    fn point_on_plane_lies_on_mohrs_circle() {
        // (normal - centre)^2 + shear^2 == radius^2 for every angle.
        let s = StressState::new(40.0, -10.0, 15.0).unwrap();
        let c = s.mean_normal();
        let r = s.radius();
        for k in 0..12 {
            let theta = k as f64 * std::f64::consts::PI / 12.0;
            let on = s.stress_on_plane(theta).unwrap();
            let lhs = (on.normal - c) * (on.normal - c) + on.shear * on.shear;
            assert!(
                (lhs - r * r).abs() < 1e-6,
                "off-circle at theta {theta}: lhs {lhs}, r^2 {}",
                r * r
            );
        }
    }

    #[test]
    fn pure_shear_principals_are_plus_minus_tau() {
        // sx=sy=0, txy=t  =>  s1=+t, s2=-t, max shear=t, principal at 45deg.
        let t = 25.0;
        let s = StressState::new(0.0, 0.0, t).unwrap();
        let p = s.principal_stresses();
        assert!((p.s1 - t).abs() < EPS);
        assert!((p.s2 + t).abs() < EPS);
        assert!((s.max_shear() - t).abs() < EPS);
        assert!((s.principal_angle() - std::f64::consts::FRAC_PI_4).abs() < EPS);
    }

    #[test]
    fn uniaxial_principals_match_inputs() {
        // Pure uniaxial tension: principals are (sx, 0), max shear sx/2.
        let s = StressState::new(50.0, 0.0, 0.0).unwrap();
        let p = s.principal_stresses();
        assert!((p.s1 - 50.0).abs() < EPS);
        assert!(p.s2.abs() < EPS);
        assert!((s.max_shear() - 25.0).abs() < EPS);
        // Already aligned with principal axes => angle 0.
        assert!(s.principal_angle().abs() < EPS);
    }

    #[test]
    fn hydrostatic_state_has_zero_radius_and_defined_angle() {
        // sx=sy, txy=0: circle degenerates to a point; angle defaults 0.
        let s = StressState::new(30.0, 30.0, 0.0).unwrap();
        assert!(s.radius().abs() < EPS);
        let p = s.principal_stresses();
        assert!((p.s1 - 30.0).abs() < EPS);
        assert!((p.s2 - 30.0).abs() < EPS);
        assert!(s.principal_angle().abs() < EPS);
        // Every plane is shear-free.
        let on = s.stress_on_plane(1.234).unwrap();
        assert!(on.shear.abs() < EPS);
        assert!((on.normal - 30.0).abs() < EPS);
    }

    #[test]
    fn serde_round_trip() {
        let s = StressState::new(40.0, -10.0, 15.0).unwrap();
        let json = serde_json::to_string(&s).unwrap();
        let back: StressState = serde_json::from_str(&json).unwrap();
        assert_eq!(s, back);

        let p = s.principal_stresses();
        let pj = serde_json::to_string(&p).unwrap();
        let pback: PrincipalStresses = serde_json::from_str(&pj).unwrap();
        assert_eq!(p, pback);
    }
}
