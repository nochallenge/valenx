//! Back-drive (self-locking) analysis.
//!
//! A power screw is **self-locking** when an axial load on the nut
//! cannot, by itself, spin the screw — i.e. the load stays put with no
//! holding torque. The textbook criterion compares the thread's lead
//! angle `lambda` against the friction angle `phi = atan(mu)`:
//!
//! ```text
//! self-locking            <=>   mu  >=  tan(lambda)     (lambda <= phi)
//! back-drivable           <=>   mu  <   tan(lambda)     (lambda >  phi)
//! ```
//!
//! Equivalently, since `tan(lambda) = lead / (pi * d_m)`, the screw is
//! back-drivable exactly when the lead is steep enough that
//! `lead / (pi * d_m) > mu`.
//!
//! The boundary `mu == tan(lambda)` is the marginal case; this module
//! classifies it as self-locking (`>=`), matching the conservative
//! convention used in machine-design texts (Shigley).

use crate::error::{require_non_negative_friction, LeadScrewError};
use crate::screw::LeadScrew;

/// Result of a back-drive analysis for a screw / friction pairing.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BackDrive {
    /// Lead angle of the thread (radians).
    pub lead_angle_rad: f64,
    /// Friction angle `phi = atan(mu)` (radians).
    pub friction_angle_rad: f64,
    /// `true` when the screw is self-locking (an axial load alone cannot
    /// back-drive it); `false` when the load can spin the screw.
    pub self_locking: bool,
}

impl BackDrive {
    /// `true` when an axial load **can** spin the screw — the logical
    /// negation of [`self_locking`](Self::self_locking).
    pub fn back_drivable(&self) -> bool {
        !self.self_locking
    }

    /// Margin = `friction_angle - lead_angle` (radians).
    ///
    /// Positive when self-locking (with headroom), zero at the marginal
    /// boundary, negative when back-drivable. A handy single number for
    /// "how far from the locking threshold are we".
    pub fn locking_margin_rad(&self) -> f64 {
        self.friction_angle_rad - self.lead_angle_rad
    }
}

impl LeadScrew {
    /// Classify whether this screw self-locks under a nut/thread
    /// friction coefficient `mu`.
    ///
    /// Implements `self_locking <=> mu >= tan(lambda)`, where
    /// `lambda` is the screw's [`lead_angle_rad`](LeadScrew::lead_angle_rad).
    ///
    /// # Errors
    ///
    /// Returns [`LeadScrewError::NegativeFriction`] /
    /// [`LeadScrewError::NotFinite`] if `mu` is negative or non-finite.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_leadscrew::LeadScrew;
    /// // Shallow lead (1 mm) on a fat screw (10 mm) -> small lead angle,
    /// // self-locking under typical steel-on-steel mu ~ 0.15.
    /// let fine = LeadScrew::new(1.0, 10.0).unwrap();
    /// assert!(fine.back_drive(0.15).unwrap().self_locking);
    ///
    /// // A steep multi-start lead (12 mm) on the same screw back-drives.
    /// let steep = LeadScrew::new(12.0, 10.0).unwrap();
    /// assert!(steep.back_drive(0.15).unwrap().back_drivable());
    /// ```
    pub fn back_drive(&self, mu: f64) -> Result<BackDrive, LeadScrewError> {
        let mu = require_non_negative_friction("mu", mu)?;
        let lead_angle_rad = self.lead_angle_rad();
        let friction_angle_rad = mu.atan();
        // Compare directly on tangents to avoid any atan round-trip:
        // self-locking  <=>  mu >= tan(lambda).
        let self_locking = mu >= lead_angle_rad.tan();
        Ok(BackDrive {
            lead_angle_rad,
            friction_angle_rad,
            self_locking,
        })
    }

    /// The friction coefficient at the self-locking boundary for this
    /// screw: `mu_crit = tan(lambda)`.
    ///
    /// Any `mu >= mu_crit` self-locks; any `mu < mu_crit` back-drives.
    /// Because `tan(lambda) = lead / (pi * d_m)`, this is also the
    /// minimum friction needed to hold a load without a brake.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_leadscrew::LeadScrew;
    /// use std::f64::consts::PI;
    /// let screw = LeadScrew::new(2.0, 8.0).unwrap();
    /// // tan(lambda) = lead / (pi * d_m) = 2 / (pi*8).
    /// let mu_crit = screw.critical_friction();
    /// assert!((mu_crit - 2.0 / (PI * 8.0)).abs() < 1e-12);
    /// ```
    pub fn critical_friction(&self) -> f64 {
        self.lead_angle_rad().tan()
    }

    /// Ideal **screw efficiency** raising a load under a nut/thread
    /// friction coefficient `mu`:
    ///
    /// `eta = tan(lambda) / tan(lambda + phi)`,   `phi = atan(mu)`.
    ///
    /// This is the square-thread power-screw efficiency (Shigley): the
    /// fraction of the input torque-work delivered as useful axial work
    /// when raising a load, computed from the geometry (through the lead
    /// angle `lambda`) and the friction alone — i.e. the `eta` a caller
    /// would otherwise have to supply to [`LeadScrew::thrust_n`]. It is
    /// `1` when frictionless (`mu = 0`) and falls as friction rises.
    ///
    /// A self-locking screw is necessarily inefficient: whenever
    /// `mu >= tan(lambda)` (the [`LeadScrew::back_drive`] self-locking
    /// condition) the efficiency is below `0.5`. The result is clamped to
    /// `[0, 1]`, so a deeply self-locking screw (where `lambda + phi`
    /// approaches `90°`) reports `0` rather than a non-physical value.
    ///
    /// # Errors
    ///
    /// Returns [`LeadScrewError::NegativeFriction`] /
    /// [`LeadScrewError::NotFinite`] if `mu` is negative or non-finite.
    pub fn screw_efficiency(&self, mu: f64) -> Result<f64, LeadScrewError> {
        let mu = require_non_negative_friction("mu", mu)?;
        let lambda = self.lead_angle_rad();
        let phi = mu.atan();
        let eta = lambda.tan() / (lambda + phi).tan();
        Ok(eta.clamp(0.0, 1.0))
    }
}

#[cfg(test)]
mod tests {
    use crate::LeadScrew;

    const EPS: f64 = 1e-9;

    #[test]
    fn efficiency_is_unity_when_frictionless() {
        // mu = 0 => phi = 0 => eta = tan(lambda)/tan(lambda) = 1.
        let screw = LeadScrew::new(2.0, 8.0).unwrap();
        let eta = screw.screw_efficiency(0.0).unwrap();
        assert!((eta - 1.0).abs() < EPS, "eta = {eta}");
    }

    #[test]
    fn efficiency_matches_closed_form() {
        // eta = tan(lambda) / tan(lambda + atan(mu)), recomputed.
        let screw = LeadScrew::new(2.0, 8.0).unwrap();
        let mu = 0.15_f64;
        let lambda = screw.lead_angle_rad();
        let expected = lambda.tan() / (lambda + mu.atan()).tan();
        let eta = screw.screw_efficiency(mu).unwrap();
        assert!((eta - expected).abs() < 1e-12, "eta = {eta} vs {expected}");
    }

    #[test]
    fn self_locking_implies_efficiency_below_half() {
        // The classic theorem: a self-locking screw has eta < 0.5.
        // A fine 1 mm lead on a fat 10 mm screw self-locks at mu = 0.15.
        let screw = LeadScrew::new(1.0, 10.0).unwrap();
        let mu = 0.15;
        assert!(screw.back_drive(mu).unwrap().self_locking);
        assert!(screw.screw_efficiency(mu).unwrap() < 0.5);
        // And exactly at the self-locking boundary mu = tan(lambda).
        let mu_crit = screw.critical_friction();
        assert!(screw.back_drive(mu_crit).unwrap().self_locking);
        assert!(screw.screw_efficiency(mu_crit).unwrap() < 0.5);
    }

    #[test]
    fn efficiency_decreases_with_friction_and_stays_in_unit_interval() {
        let screw = LeadScrew::new(4.0, 8.0).unwrap();
        let mut prev = f64::INFINITY;
        for &mu in &[0.0, 0.05, 0.1, 0.2, 0.4] {
            let eta = screw.screw_efficiency(mu).unwrap();
            assert!(eta < prev, "eta should fall with friction: {eta} !< {prev}");
            assert!((0.0..=1.0).contains(&eta), "eta out of [0,1]: {eta}");
            prev = eta;
        }
    }

    #[test]
    fn efficiency_rejects_bad_friction() {
        let screw = LeadScrew::new(2.0, 8.0).unwrap();
        assert!(screw.screw_efficiency(-0.1).is_err());
        assert!(screw.screw_efficiency(f64::NAN).is_err());
    }
}
