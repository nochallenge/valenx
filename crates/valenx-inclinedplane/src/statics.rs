//! Frictional statics of a block resting on a ramp.
//!
//! A [`Ramp`] bundles a slope angle, a Coulomb friction coefficient,
//! and a block weight. [`Ramp::forces`] returns the full
//! [`RampForces`] breakdown (normal force, slope force, friction,
//! raise/lower effort, self-locking flag) in one shot, while the
//! individual accessors compute each quantity on its own.
//!
//! All forces are expressed in the same unit as the supplied weight
//! `weight` (newtons, pounds-force, ...); the angle is in radians.

use serde::{Deserialize, Serialize};

use crate::error::InclinedPlaneError;
use crate::geometry::IdealRamp;

/// A block of given weight on a ramp with Coulomb friction.
///
/// The slope angle is reused from [`IdealRamp`] (and therefore shares
/// its `(0, pi/2)` validation). The friction coefficient must be
/// finite and non-negative; the weight must be finite and strictly
/// positive.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Ramp {
    geometry: IdealRamp,
    mu: f64,
    weight: f64,
}

/// The complete static force breakdown for a [`Ramp`].
///
/// Each field is in the same force unit as the originating
/// [`Ramp::weight`]. `effort_to_lower` may be negative — see its field
/// docs and [`RampForces::is_self_locking`].
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RampForces {
    /// Normal reaction `N = W cos(theta)`.
    pub normal: f64,
    /// Down-slope gravity component `W sin(theta)`.
    pub slope_force: f64,
    /// Maximum available Coulomb friction `mu * N`.
    pub friction: f64,
    /// Slope-parallel effort to raise the load,
    /// `W (sin(theta) + mu cos(theta))`.
    pub effort_to_raise: f64,
    /// Slope-parallel effort to lower / restrain the load,
    /// `W (sin(theta) - mu cos(theta))`.
    ///
    /// A positive value is the force needed (acting up the slope) to
    /// keep the block from sliding down; a non-positive value means
    /// friction alone holds the block and the ramp is self-locking.
    pub effort_to_lower: f64,
    /// `true` when the ramp is self-locking (the block will not slide
    /// under gravity alone).
    pub is_self_locking: bool,
}

impl Ramp {
    /// Build a ramp from a slope angle (radians), friction coefficient,
    /// and block weight.
    ///
    /// # Errors
    ///
    /// Returns [`InclinedPlaneError::BadParameter`] if the angle is not
    /// in `(0, pi/2)`, if `mu` is non-finite or negative, or if
    /// `weight` is non-finite or not strictly positive.
    pub fn new(angle_rad: f64, mu: f64, weight: f64) -> Result<Self, InclinedPlaneError> {
        let geometry = IdealRamp::from_angle(angle_rad)?;
        Self::from_geometry(geometry, mu, weight)
    }

    /// Build a ramp from an existing [`IdealRamp`] plus friction and
    /// weight.
    ///
    /// # Errors
    ///
    /// Returns [`InclinedPlaneError::BadParameter`] if `mu` is
    /// non-finite or negative, or if `weight` is non-finite or not
    /// strictly positive.
    pub fn from_geometry(
        geometry: IdealRamp,
        mu: f64,
        weight: f64,
    ) -> Result<Self, InclinedPlaneError> {
        if !mu.is_finite() || mu < 0.0 {
            return Err(InclinedPlaneError::bad_parameter(
                "mu",
                format!("friction coefficient must be finite and >= 0, got {mu}"),
            ));
        }
        if !weight.is_finite() || weight <= 0.0 {
            return Err(InclinedPlaneError::bad_parameter(
                "weight",
                format!("must be finite and > 0, got {weight}"),
            ));
        }
        Ok(Self {
            geometry,
            mu,
            weight,
        })
    }

    /// The underlying frictionless geometry.
    pub fn geometry(&self) -> IdealRamp {
        self.geometry
    }

    /// The Coulomb friction coefficient `mu`.
    pub fn mu(&self) -> f64 {
        self.mu
    }

    /// The block weight `W`.
    pub fn weight(&self) -> f64 {
        self.weight
    }

    /// Normal reaction `N = W cos(theta)`.
    pub fn normal_force(&self) -> f64 {
        self.weight * self.geometry.angle_rad().cos()
    }

    /// Down-slope gravity component `W sin(theta)` — the force that
    /// tries to slide the block down the ramp.
    pub fn slope_force(&self) -> f64 {
        self.weight * self.geometry.angle_rad().sin()
    }

    /// Maximum available Coulomb friction force `mu * N`.
    pub fn friction_force(&self) -> f64 {
        self.mu * self.normal_force()
    }

    /// Slope-parallel effort to raise the load steadily up the ramp,
    /// `F_up = W (sin(theta) + mu cos(theta))`.
    ///
    /// Friction opposes the upward motion, so it adds to the
    /// gravitational term.
    pub fn effort_to_raise(&self) -> f64 {
        let theta = self.geometry.angle_rad();
        self.weight * (theta.sin() + self.mu * theta.cos())
    }

    /// Slope-parallel effort to lower / restrain the load,
    /// `F_down = W (sin(theta) - mu cos(theta))`.
    ///
    /// When this is positive it is the holding force (directed up the
    /// slope) required to stop the block sliding down; when it is zero
    /// or negative, friction alone is enough and the ramp is
    /// self-locking (see [`Ramp::is_self_locking`]).
    pub fn effort_to_lower(&self) -> f64 {
        let theta = self.geometry.angle_rad();
        self.weight * (theta.sin() - self.mu * theta.cos())
    }

    /// The friction angle `phi = atan(mu)` in radians.
    ///
    /// This is the steepest slope on which the block can rest without
    /// sliding: at `theta = phi` the down-slope gravity component
    /// exactly equals the maximum friction.
    pub fn friction_angle_rad(&self) -> f64 {
        self.mu.atan()
    }

    /// Whether the ramp is self-locking: the block will not slide
    /// under gravity alone.
    ///
    /// This holds exactly when the friction angle is at least the slope
    /// angle, `atan(mu) >= theta`, equivalently `mu >= tan(theta)`,
    /// equivalently the down-slope force does not exceed the available
    /// friction (`W sin(theta) <= mu W cos(theta)`).
    pub fn is_self_locking(&self) -> bool {
        self.friction_angle_rad() >= self.geometry.angle_rad()
    }

    /// The **actual mechanical advantage** when raising the load,
    /// `AMA = W / F_up = 1 / (sin(theta) + mu cos(theta))`.
    ///
    /// The load-to-effort ratio actually achieved against friction. It
    /// equals the ideal [`IdealRamp::ideal_mechanical_advantage`]
    /// (`1 / sin(theta)`) only in the frictionless limit `mu = 0`; any
    /// friction lowers it. Independent of the block weight, since both `W`
    /// and `F_up` scale with it.
    pub fn actual_mechanical_advantage(&self) -> f64 {
        self.weight / self.effort_to_raise()
    }

    /// The **efficiency** of the ramp when raising the load,
    /// `eta = sin(theta) / (sin(theta) + mu cos(theta))` — the ratio of the
    /// frictionless effort to the actual effort, equivalently the actual
    /// mechanical advantage over the ideal one. Lies in `(0, 1]`.
    ///
    /// `eta = 1` exactly when frictionless (`mu = 0`); any friction makes
    /// it strictly less. A self-locking ramp ([`Ramp::is_self_locking`])
    /// always has `eta <= 1/2` — the classic result that a machine which
    /// holds its own load under friction wastes at least half its input
    /// work.
    pub fn efficiency(&self) -> f64 {
        self.slope_force() / self.effort_to_raise()
    }

    /// Compute every static quantity in one pass.
    pub fn forces(&self) -> RampForces {
        RampForces {
            normal: self.normal_force(),
            slope_force: self.slope_force(),
            friction: self.friction_force(),
            effort_to_raise: self.effort_to_raise(),
            effort_to_lower: self.effort_to_lower(),
            is_self_locking: self.is_self_locking(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    /// A 30-degree ramp, mu = 0.2, weight = 100. Hand-computed
    /// reference values used across several tests.
    fn ramp_30deg() -> Ramp {
        Ramp::new(std::f64::consts::FRAC_PI_6, 0.2, 100.0).unwrap()
    }

    #[test]
    fn rejects_bad_parameters() {
        // Bad angle (vertical).
        assert!(Ramp::new(std::f64::consts::FRAC_PI_2, 0.2, 100.0).is_err());
        // Negative friction.
        assert!(Ramp::new(std::f64::consts::FRAC_PI_6, -0.1, 100.0).is_err());
        // Non-finite friction.
        assert!(Ramp::new(std::f64::consts::FRAC_PI_6, f64::NAN, 100.0).is_err());
        // Zero / negative weight.
        assert!(Ramp::new(std::f64::consts::FRAC_PI_6, 0.2, 0.0).is_err());
        assert!(Ramp::new(std::f64::consts::FRAC_PI_6, 0.2, -5.0).is_err());
    }

    #[test]
    fn normal_force_is_weight_times_cos() {
        let ramp = ramp_30deg();
        // N = 100 * cos(30) = 100 * (sqrt(3)/2) = 86.6025403784...
        let expected = 100.0 * (3.0_f64).sqrt() / 2.0;
        assert!((ramp.normal_force() - expected).abs() < EPS);
    }

    #[test]
    fn slope_force_is_weight_times_sin() {
        let ramp = ramp_30deg();
        // W sin(30) = 100 * 0.5 = 50 exactly.
        assert!((ramp.slope_force() - 50.0).abs() < 1e-11);
    }

    #[test]
    fn friction_is_mu_times_normal() {
        let ramp = ramp_30deg();
        // mu N = 0.2 * 86.60254... = 17.320508...
        let expected = 0.2 * 100.0 * (3.0_f64).sqrt() / 2.0;
        assert!((ramp.friction_force() - expected).abs() < EPS);
    }

    #[test]
    fn effort_to_raise_adds_friction() {
        let ramp = ramp_30deg();
        // W(sin30 + mu cos30) = 100*(0.5 + 0.2*0.8660254) = 67.320508...
        let expected = 100.0 * (0.5 + 0.2 * (3.0_f64).sqrt() / 2.0);
        assert!((ramp.effort_to_raise() - expected).abs() < EPS);
        // Raising must cost more than the frictionless slope force.
        assert!(ramp.effort_to_raise() > ramp.slope_force());
    }

    #[test]
    fn effort_to_lower_subtracts_friction() {
        let ramp = ramp_30deg();
        // W(sin30 - mu cos30) = 100*(0.5 - 0.2*0.8660254) = 32.679492...
        let expected = 100.0 * (0.5 - 0.2 * (3.0_f64).sqrt() / 2.0);
        assert!((ramp.effort_to_lower() - expected).abs() < EPS);
        // Still positive here -> not self-locking.
        assert!(ramp.effort_to_lower() > 0.0);
        assert!(!ramp.is_self_locking());
    }

    #[test]
    fn raise_minus_lower_is_twice_friction() {
        // F_up - F_down = 2 mu W cos(theta) = 2 * friction.
        let ramp = ramp_30deg();
        let delta = ramp.effort_to_raise() - ramp.effort_to_lower();
        assert!((delta - 2.0 * ramp.friction_force()).abs() < EPS);
    }

    #[test]
    fn self_locking_when_friction_angle_exceeds_slope() {
        // 20-degree ramp, mu = 0.5 -> friction angle atan(0.5)=26.57deg
        // > 20deg, so it locks and needs no holding force.
        let theta = 20.0_f64.to_radians();
        let ramp = Ramp::new(theta, 0.5, 100.0).unwrap();
        assert!(ramp.is_self_locking());
        // effort_to_lower is the holding force; self-locking => <= 0.
        assert!(ramp.effort_to_lower() <= 0.0);
        assert!(ramp.friction_force() >= ramp.slope_force());
    }

    #[test]
    fn not_self_locking_when_slope_exceeds_friction_angle() {
        // 40-degree ramp, mu = 0.5 -> friction angle 26.57deg < 40deg.
        let theta = 40.0_f64.to_radians();
        let ramp = Ramp::new(theta, 0.5, 100.0).unwrap();
        assert!(!ramp.is_self_locking());
        assert!(ramp.effort_to_lower() > 0.0);
        assert!(ramp.slope_force() > ramp.friction_force());
    }

    #[test]
    fn self_locking_boundary_at_friction_angle() {
        // Set the slope exactly equal to the friction angle: the
        // criterion is `>=`, so the boundary is self-locking and the
        // holding force is (within rounding) zero.
        let mu = 0.3_f64;
        let theta = mu.atan();
        let ramp = Ramp::new(theta, mu, 100.0).unwrap();
        assert!(ramp.is_self_locking());
        assert!(ramp.effort_to_lower().abs() < 1e-9);
        // slope force and friction coincide at the boundary.
        assert!((ramp.slope_force() - ramp.friction_force()).abs() < 1e-9);
    }

    #[test]
    fn frictionless_ramp_never_locks_and_efforts_equal_slope_force() {
        // mu = 0: friction angle is 0 < theta, so never self-locking;
        // raise == lower == slope force == W sin(theta).
        let ramp = Ramp::new(std::f64::consts::FRAC_PI_6, 0.0, 100.0).unwrap();
        assert!(!ramp.is_self_locking());
        assert!((ramp.effort_to_raise() - ramp.slope_force()).abs() < EPS);
        assert!((ramp.effort_to_lower() - ramp.slope_force()).abs() < EPS);
        // Frictionless MA matches the ideal geometry.
        assert!((ramp.geometry().ideal_mechanical_advantage() - 2.0).abs() < 1e-12);
    }

    #[test]
    fn forces_struct_matches_individual_accessors() {
        let ramp = ramp_30deg();
        let f = ramp.forces();
        assert!((f.normal - ramp.normal_force()).abs() < EPS);
        assert!((f.slope_force - ramp.slope_force()).abs() < EPS);
        assert!((f.friction - ramp.friction_force()).abs() < EPS);
        assert!((f.effort_to_raise - ramp.effort_to_raise()).abs() < EPS);
        assert!((f.effort_to_lower - ramp.effort_to_lower()).abs() < EPS);
        assert_eq!(f.is_self_locking, ramp.is_self_locking());
    }

    #[test]
    fn friction_angle_is_atan_mu() {
        let ramp = ramp_30deg();
        assert!((ramp.friction_angle_rad() - 0.2_f64.atan()).abs() < EPS);
    }

    #[test]
    fn efficiency_and_ama_hand_values() {
        let ramp = ramp_30deg();
        // denom = sin30 + mu cos30 = 0.5 + 0.2*sqrt(3)/2 = 0.6732050808...
        let denom = 0.5 + 0.2 * (3.0_f64).sqrt() / 2.0;
        // eta = sin30 / denom = 0.7427813527...
        assert!((ramp.efficiency() - 0.5 / denom).abs() < EPS);
        // AMA = W / F_up = 1 / denom = 1.4855627054...
        assert!((ramp.actual_mechanical_advantage() - 1.0 / denom).abs() < EPS);
    }

    #[test]
    fn frictionless_efficiency_is_one_and_ama_equals_ideal() {
        // mu = 0: no friction wasted, so efficiency is exactly 1 and the
        // actual MA collapses onto the ideal geometric MA (1/sin30 = 2).
        let ramp = Ramp::new(std::f64::consts::FRAC_PI_6, 0.0, 100.0).unwrap();
        assert!((ramp.efficiency() - 1.0).abs() < EPS);
        let ima = ramp.geometry().ideal_mechanical_advantage();
        assert!((ramp.actual_mechanical_advantage() - ima).abs() < EPS);
        assert!((ramp.actual_mechanical_advantage() - 2.0).abs() < EPS);
    }

    #[test]
    fn efficiency_equals_actual_over_ideal_ma() {
        // The defining identity linking the two: eta = AMA / IMA.
        let ramp = ramp_30deg();
        let ima = ramp.geometry().ideal_mechanical_advantage();
        let ratio = ramp.actual_mechanical_advantage() / ima;
        assert!((ramp.efficiency() - ratio).abs() < EPS);
    }

    #[test]
    fn efficiency_and_ama_independent_of_weight() {
        // Both are pure ratios, so scaling the load leaves them fixed.
        let light = Ramp::new(std::f64::consts::FRAC_PI_6, 0.2, 5.0).unwrap();
        let heavy = Ramp::new(std::f64::consts::FRAC_PI_6, 0.2, 5000.0).unwrap();
        assert!((light.efficiency() - heavy.efficiency()).abs() < EPS);
        assert!(
            (light.actual_mechanical_advantage() - heavy.actual_mechanical_advantage()).abs() < EPS
        );
    }

    #[test]
    fn self_locking_efficiency_at_most_one_half() {
        // Classic result: a self-locking machine has eta <= 1/2, while a
        // freely-sliding ramp exceeds 1/2.
        let locking = Ramp::new(20.0_f64.to_radians(), 0.5, 100.0).unwrap();
        assert!(locking.is_self_locking());
        assert!(
            locking.efficiency() <= 0.5,
            "eta = {}",
            locking.efficiency()
        );
        let free = Ramp::new(40.0_f64.to_radians(), 0.5, 100.0).unwrap();
        assert!(!free.is_self_locking());
        assert!(free.efficiency() > 0.5, "eta = {}", free.efficiency());
    }

    #[test]
    fn more_friction_lowers_efficiency_and_ama() {
        let theta = std::f64::consts::FRAC_PI_6;
        let low = Ramp::new(theta, 0.1, 100.0).unwrap();
        let high = Ramp::new(theta, 0.4, 100.0).unwrap();
        assert!(high.efficiency() < low.efficiency());
        assert!(high.actual_mechanical_advantage() < low.actual_mechanical_advantage());
    }

    #[test]
    fn serde_round_trips() {
        let ramp = ramp_30deg();
        let json = serde_json::to_string(&ramp).unwrap();
        let back: Ramp = serde_json::from_str(&json).unwrap();
        assert_eq!(ramp, back);
    }
}
