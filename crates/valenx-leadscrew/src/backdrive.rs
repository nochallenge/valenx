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
}
