//! Square-thread power / lead screw model.
//!
//! All torque expressions follow the standard textbook development for
//! a square-thread power screw (Shigley / Budynas, *Mechanical
//! Engineering Design*). Acme and trapezoidal threads are approximated
//! by the square-thread relations here; a true Acme treatment divides
//! the friction term by `cos(alpha_n)` for the thread half-angle, which
//! this crate deliberately omits to keep the closed forms transparent.
//!
//! # Symbols
//!
//! `F`   axial load,
//! `dm`  mean (pitch-line) diameter,
//! `p`   single-thread pitch,
//! `l`   lead `= p * starts`,
//! `mu`  coefficient of friction between screw and nut,
//! `lambda`  lead (helix) angle, `tan(lambda) = l / (pi * dm)`.
//!
//! Collar friction is not modelled; the torques returned are the
//! thread-only components.

use serde::{Deserialize, Serialize};

use crate::error::ScrewThreadError;

/// Validated geometry + friction state for a square-thread power
/// screw.
///
/// Build instances with [`ScrewThread::new`], which enforces the
/// physical admissibility of every field. The accessors and analysis
/// methods then operate on a guaranteed-valid screw, so they are
/// infallible apart from [`ScrewThread::raise_torque`], which can still
/// hit the `pi*dm - mu*l` singularity for extreme leads.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct ScrewThread {
    /// Mean (pitch-line) diameter `dm`. Same length unit as the pitch.
    mean_diameter: f64,
    /// Single-thread pitch `p` (axial distance between adjacent
    /// threads). Same length unit as the mean diameter.
    pitch: f64,
    /// Number of independent thread starts (1 for a single-start
    /// screw, 2 for double-start, ...).
    starts: u32,
    /// Coefficient of friction `mu` between the screw and the nut
    /// (dimensionless, non-negative).
    friction: f64,
}

impl ScrewThread {
    /// Build a validated power screw.
    ///
    /// # Parameters
    ///
    /// `mean_diameter` mean (pitch-line) diameter, must be `> 0`;
    /// `pitch` single-thread pitch, must be `> 0`;
    /// `starts` number of thread starts, must be `>= 1`;
    /// `friction` Coulomb friction coefficient, must be `>= 0` and
    /// finite.
    ///
    /// # Errors
    ///
    /// Returns [`ScrewThreadError::InvalidParameter`] if any argument is
    /// non-finite or outside the admissible range described above.
    pub fn new(
        mean_diameter: f64,
        pitch: f64,
        starts: u32,
        friction: f64,
    ) -> Result<Self, ScrewThreadError> {
        if !mean_diameter.is_finite() || mean_diameter <= 0.0 {
            return Err(ScrewThreadError::InvalidParameter {
                name: "mean_diameter",
                reason: "must be a finite value > 0",
                value: mean_diameter,
            });
        }
        if !pitch.is_finite() || pitch <= 0.0 {
            return Err(ScrewThreadError::InvalidParameter {
                name: "pitch",
                reason: "must be a finite value > 0",
                value: pitch,
            });
        }
        if starts == 0 {
            return Err(ScrewThreadError::InvalidParameter {
                name: "starts",
                reason: "must be >= 1",
                value: starts as f64,
            });
        }
        if !friction.is_finite() || friction < 0.0 {
            return Err(ScrewThreadError::InvalidParameter {
                name: "friction",
                reason: "must be a finite value >= 0",
                value: friction,
            });
        }
        Ok(Self {
            mean_diameter,
            pitch,
            starts,
            friction,
        })
    }

    /// Mean (pitch-line) diameter `dm`.
    pub fn mean_diameter(&self) -> f64 {
        self.mean_diameter
    }

    /// Single-thread pitch `p`.
    pub fn pitch(&self) -> f64 {
        self.pitch
    }

    /// Number of thread starts.
    pub fn starts(&self) -> u32 {
        self.starts
    }

    /// Coefficient of friction `mu`.
    pub fn friction(&self) -> f64 {
        self.friction
    }

    /// Lead `l = p * starts` — the axial advance per full turn.
    pub fn lead(&self) -> f64 {
        self.pitch * self.starts as f64
    }

    /// Lead (helix) angle `lambda = atan(l / (pi * dm))`, in radians.
    ///
    /// This is the angle the thread helix makes with a plane normal to
    /// the screw axis.
    pub fn lead_angle(&self) -> f64 {
        (self.lead() / (std::f64::consts::PI * self.mean_diameter)).atan()
    }

    /// Friction angle `phi = atan(mu)`, in radians.
    ///
    /// The screw is self-locking exactly when `phi >= lambda` (see
    /// [`ScrewThread::is_self_locking`]).
    pub fn friction_angle(&self) -> f64 {
        self.friction.atan()
    }

    /// Raising torque `T_R` required to drive `axial_load` *against*
    /// the load direction (for example, lifting):
    ///
    /// `T_R = (F * dm / 2) * (l + pi * mu * dm) / (pi * dm - mu * l)`
    ///
    /// # Parameters
    ///
    /// `axial_load` axial force `F` resisted by the screw (same force
    /// unit as the returned torque's force component; with `dm` in mm
    /// and `F` in N the torque is in N·mm).
    ///
    /// # Errors
    ///
    /// Returns [`ScrewThreadError::SingularDenominator`] when
    /// `pi*dm - mu*l` is at or below a small epsilon, i.e. the lead is
    /// so large that the raise-torque expression is singular.
    pub fn raise_torque(&self, axial_load: f64) -> Result<f64, ScrewThreadError> {
        let dm = self.mean_diameter;
        let l = self.lead();
        let mu = self.friction;
        let pi = std::f64::consts::PI;

        let denominator = pi * dm - mu * l;
        if denominator <= f64::EPSILON {
            return Err(ScrewThreadError::SingularDenominator {
                denominator,
                mean_diameter: dm,
                lead: l,
                friction: mu,
            });
        }
        Ok((axial_load * dm / 2.0) * (l + pi * mu * dm) / denominator)
    }

    /// Lowering torque `T_L` required to drive `axial_load` *with* the
    /// load direction (for example, paying a weight down):
    ///
    /// `T_L = (F * dm / 2) * (pi * mu * dm - l) / (pi * dm + mu * l)`
    ///
    /// A positive result means torque must be applied to lower the
    /// load (the screw is self-locking); a negative result means the
    /// load would back-drive the screw and the returned magnitude is
    /// the holding (braking) torque required to keep it static.
    ///
    /// The denominator `pi*dm + mu*l` is strictly positive for any
    /// valid screw, so this method is infallible.
    pub fn lower_torque(&self, axial_load: f64) -> f64 {
        let dm = self.mean_diameter;
        let l = self.lead();
        let mu = self.friction;
        let pi = std::f64::consts::PI;

        let denominator = pi * dm + mu * l;
        (axial_load * dm / 2.0) * (pi * mu * dm - l) / denominator
    }

    /// Mechanical efficiency of the raising motion (thread friction
    /// only), in the open interval `(0, 1)` for any `mu > 0`.
    ///
    /// Defined as the ratio of the frictionless raising torque
    /// `T_0 = F * l / (2*pi)` to the actual raising torque `T_R`:
    ///
    /// `e = (F * l) / (2 * pi * T_R)`.
    ///
    /// Equivalently `e = tan(lambda) * (1 - mu * tan(lambda)) /
    /// (tan(lambda) + mu)`. With `mu = 0` the efficiency is `1`
    /// (ideal screw); as `mu` grows the efficiency falls toward `0`.
    ///
    /// # Errors
    ///
    /// Propagates [`ScrewThreadError::SingularDenominator`] from
    /// [`ScrewThread::raise_torque`].
    pub fn efficiency(&self) -> Result<f64, ScrewThreadError> {
        // Use a unit load: the load cancels in the ratio, so the
        // result is independent of F (for any F > 0).
        let f = 1.0;
        let t_r = self.raise_torque(f)?;
        let t_0 = f * self.lead() / (2.0 * std::f64::consts::PI);
        Ok(t_0 / t_r)
    }

    /// Whether the screw is self-locking — i.e. an axial load alone
    /// cannot back-drive it, so it stays put with no holding torque.
    ///
    /// The criterion is `friction angle >= lead angle`
    /// (`phi >= lambda`), which is algebraically equivalent to
    /// `mu >= tan(lambda) = l / (pi * dm)` and to a non-negative
    /// lowering torque `T_L >= 0`.
    pub fn is_self_locking(&self) -> bool {
        self.friction_angle() >= self.lead_angle()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tolerance for floating-point ground-truth comparisons.
    const EPS: f64 = 1e-9;

    fn sample() -> ScrewThread {
        // Square-thread power screw, single start.
        // dm = 25 mm, p = 5 mm, starts = 1, mu = 0.08.
        ScrewThread::new(25.0, 5.0, 1, 0.08).unwrap()
    }

    #[test]
    fn lead_is_pitch_times_starts() {
        let single = ScrewThread::new(25.0, 5.0, 1, 0.1).unwrap();
        assert!((single.lead() - 5.0).abs() < EPS);

        let double = ScrewThread::new(25.0, 5.0, 2, 0.1).unwrap();
        assert!((double.lead() - 10.0).abs() < EPS);

        let triple = ScrewThread::new(25.0, 5.0, 3, 0.1).unwrap();
        assert!((triple.lead() - 15.0).abs() < EPS);
    }

    #[test]
    fn lead_angle_matches_atan_definition() {
        let s = sample();
        let expected = (s.lead() / (std::f64::consts::PI * s.mean_diameter())).atan();
        assert!((s.lead_angle() - expected).abs() < EPS);

        // Closed-form check: l = 5, dm = 25 -> tan(lambda) = 5/(25*pi).
        let tan_lambda = s.lead_angle().tan();
        assert!((tan_lambda - 5.0 / (25.0 * std::f64::consts::PI)).abs() < EPS);
    }

    #[test]
    fn friction_angle_matches_atan_mu() {
        let s = sample();
        assert!((s.friction_angle() - 0.08_f64.atan()).abs() < EPS);
        assert!((s.friction_angle().tan() - 0.08).abs() < EPS);
    }

    #[test]
    fn raise_torque_matches_closed_form() {
        let s = sample();
        let f = 6500.0; // N
        let dm = 25.0;
        let l = 5.0;
        let mu = 0.08;
        let pi = std::f64::consts::PI;
        let expected = (f * dm / 2.0) * (l + pi * mu * dm) / (pi * dm - mu * l);
        let got = s.raise_torque(f).unwrap();
        assert!((got - expected).abs() < 1e-6);
        // Spot value (N·mm) computed independently.
        // num = 5 + pi*0.08*25 = 5 + 6.283185307 = 11.283185307
        // den = pi*25 - 0.08*5 = 78.53981634 - 0.4 = 78.13981634
        // T = (6500*25/2) * 11.283185307 / 78.13981634 ≈ 11732.29 N·mm
        assert!((got - 11_732.287_7).abs() < 1e-2);
    }

    #[test]
    fn lower_torque_matches_closed_form() {
        let s = sample();
        let f = 6500.0;
        let dm = 25.0;
        let l = 5.0;
        let mu = 0.08;
        let pi = std::f64::consts::PI;
        let expected = (f * dm / 2.0) * (pi * mu * dm - l) / (pi * dm + mu * l);
        let got = s.lower_torque(f);
        assert!((got - expected).abs() < 1e-6);
    }

    #[test]
    fn raise_exceeds_lower_torque_for_typical_screw() {
        let s = sample();
        let f = 1000.0;
        let raise = s.raise_torque(f).unwrap();
        let lower = s.lower_torque(f);
        assert!(raise > lower);
        // For a self-locking screw the lower torque is positive but
        // strictly smaller than the raise torque.
        assert!(lower > 0.0);
    }

    #[test]
    fn efficiency_in_open_unit_interval() {
        let s = sample();
        let e = s.efficiency().unwrap();
        assert!(e > 0.0, "efficiency {e} must be > 0");
        assert!(e < 1.0, "efficiency {e} must be < 1");
    }

    #[test]
    fn efficiency_equals_trig_identity() {
        let s = sample();
        let lambda = s.lead_angle();
        let mu = s.friction();
        let tan_lambda = lambda.tan();
        let identity = tan_lambda * (1.0 - mu * tan_lambda) / (tan_lambda + mu);
        let e = s.efficiency().unwrap();
        assert!((e - identity).abs() < 1e-9);
    }

    #[test]
    fn frictionless_screw_is_perfectly_efficient() {
        let s = ScrewThread::new(25.0, 5.0, 1, 0.0).unwrap();
        let e = s.efficiency().unwrap();
        assert!((e - 1.0).abs() < EPS);
    }

    #[test]
    fn efficiency_decreases_with_friction() {
        let lo = ScrewThread::new(25.0, 5.0, 1, 0.05).unwrap();
        let hi = ScrewThread::new(25.0, 5.0, 1, 0.20).unwrap();
        assert!(lo.efficiency().unwrap() > hi.efficiency().unwrap());
    }

    #[test]
    fn self_locking_when_friction_angle_exceeds_lead_angle() {
        // mu = 0.08, lead angle is small -> self-locking.
        let locking = sample();
        assert!(locking.friction_angle() > locking.lead_angle());
        assert!(locking.is_self_locking());
        // Self-locking <=> lowering torque is non-negative.
        assert!(locking.lower_torque(1000.0) >= 0.0);
    }

    #[test]
    fn not_self_locking_for_steep_low_friction_screw() {
        // Large lead (multi-start, coarse pitch), tiny friction ->
        // lead angle exceeds friction angle -> back-drives.
        let overhauling = ScrewThread::new(20.0, 12.0, 4, 0.02).unwrap();
        assert!(overhauling.lead_angle() > overhauling.friction_angle());
        assert!(!overhauling.is_self_locking());
        // Back-driving screw: lowering torque is negative (the load
        // drives the screw; a brake torque of |T_L| is needed).
        assert!(overhauling.lower_torque(1000.0) < 0.0);
    }

    #[test]
    fn self_locking_boundary_is_consistent_with_tan_lead_angle() {
        // Construct a screw whose friction exactly equals tan(lambda).
        let dm = 30.0;
        let p = 6.0;
        let starts = 1;
        let l = p * starts as f64;
        let tan_lambda = l / (std::f64::consts::PI * dm);
        let s = ScrewThread::new(dm, p, starts, tan_lambda).unwrap();
        // At the boundary mu == tan(lambda) => phi == lambda => the
        // lowering torque is (within rounding) zero.
        assert!(s.lower_torque(1000.0).abs() < 1e-6);
        // friction_angle and lead_angle coincide at the boundary.
        assert!((s.friction_angle() - s.lead_angle()).abs() < EPS);
        assert!(s.is_self_locking());
    }

    #[test]
    fn raise_torque_scales_linearly_with_load() {
        let s = sample();
        let t1 = s.raise_torque(1000.0).unwrap();
        let t2 = s.raise_torque(2000.0).unwrap();
        assert!((t2 - 2.0 * t1).abs() < 1e-6);
    }

    #[test]
    fn rejects_non_positive_diameter() {
        let err = ScrewThread::new(0.0, 5.0, 1, 0.1).unwrap_err();
        assert_eq!(err.code(), "screwthread.invalid_parameter");
    }

    #[test]
    fn rejects_non_positive_pitch() {
        let err = ScrewThread::new(25.0, -1.0, 1, 0.1).unwrap_err();
        assert_eq!(err.code(), "screwthread.invalid_parameter");
    }

    #[test]
    fn rejects_zero_starts() {
        let err = ScrewThread::new(25.0, 5.0, 0, 0.1).unwrap_err();
        assert_eq!(err.code(), "screwthread.invalid_parameter");
    }

    #[test]
    fn rejects_negative_friction() {
        let err = ScrewThread::new(25.0, 5.0, 1, -0.1).unwrap_err();
        assert_eq!(err.code(), "screwthread.invalid_parameter");
    }

    #[test]
    fn rejects_non_finite_parameter() {
        assert!(ScrewThread::new(f64::NAN, 5.0, 1, 0.1).is_err());
        assert!(ScrewThread::new(25.0, f64::INFINITY, 1, 0.1).is_err());
        assert!(ScrewThread::new(25.0, 5.0, 1, f64::NAN).is_err());
    }

    #[test]
    fn singular_denominator_is_reported() {
        // Force pi*dm - mu*l <= 0 by making the lead huge and friction
        // large relative to the mean diameter.
        // dm = 1, p = 100, starts = 1 -> l = 100; mu = 1.0
        // denom = pi*1 - 1*100 = 3.14159 - 100 < 0.
        let s = ScrewThread::new(1.0, 100.0, 1, 1.0).unwrap();
        let err = s.raise_torque(10.0).unwrap_err();
        assert_eq!(err.code(), "screwthread.singular_denominator");
        // Efficiency, which leans on raise_torque, propagates it.
        assert!(s.efficiency().is_err());
    }
}
