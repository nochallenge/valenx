//! The [`LeadScrew`] geometry type and its kinematic / static queries.
//!
//! A lead screw converts the rotation of a threaded shaft into linear
//! translation of a nut. Two geometric numbers fully describe the ideal
//! kinematics:
//!
//! - **lead** `l` — the axial distance the nut advances per **one full
//!   revolution** of the screw (mm/rev). For a single-start thread this
//!   equals the pitch; for an `s`-start thread it is `s * pitch`.
//! - **pitch diameter** `d_m` — the mean diameter at which the thread
//!   forces act (mm). Only needed for the lead-angle / back-drive math.
//!
//! All public methods are exact evaluations of the textbook constant-
//! lead screw relations; none allocate and none can panic for an input
//! that the validating constructor accepted.

use crate::error::{require_positive, LeadScrewError};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// An idealized constant-lead power screw (lead screw or ball screw).
///
/// Construct with [`LeadScrew::new`], which validates that both the
/// lead and the pitch diameter are finite and strictly positive. The
/// fields are public for serialization but the invariants are only
/// guaranteed for values produced through the constructor.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct LeadScrew {
    /// Lead — axial advance per revolution (mm/rev). Strictly positive.
    pub lead_mm: f64,
    /// Pitch (mean) diameter of the thread (mm). Strictly positive.
    pub pitch_diameter_mm: f64,
}

impl LeadScrew {
    /// Build a validated lead screw.
    ///
    /// # Errors
    ///
    /// Returns [`LeadScrewError::NotPositive`] or
    /// [`LeadScrewError::NotFinite`] if either `lead_mm` or
    /// `pitch_diameter_mm` is non-positive or non-finite.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_leadscrew::LeadScrew;
    /// // An 8 mm-lead, 8 mm pitch-diameter screw (a common "T8x8" leadscrew).
    /// let screw = LeadScrew::new(8.0, 8.0).unwrap();
    /// assert!(LeadScrew::new(0.0, 8.0).is_err());
    /// ```
    pub fn new(lead_mm: f64, pitch_diameter_mm: f64) -> Result<Self, LeadScrewError> {
        let lead_mm = require_positive("lead_mm", lead_mm)?;
        let pitch_diameter_mm = require_positive("pitch_diameter_mm", pitch_diameter_mm)?;
        Ok(Self {
            lead_mm,
            pitch_diameter_mm,
        })
    }

    /// Linear speed of the nut for a given screw speed.
    ///
    /// `v = lead * rpm`. With `lead` in mm/rev and `rpm` in rev/min the
    /// result is **mm/min**.
    ///
    /// # Errors
    ///
    /// Returns an error if `rpm` is non-finite or non-positive.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_leadscrew::LeadScrew;
    /// let screw = LeadScrew::new(2.0, 8.0).unwrap();
    /// // 2 mm/rev * 300 rev/min = 600 mm/min.
    /// let v = screw.linear_speed_mm_per_min(300.0).unwrap();
    /// assert!((v - 600.0).abs() < 1e-9);
    /// ```
    pub fn linear_speed_mm_per_min(&self, rpm: f64) -> Result<f64, LeadScrewError> {
        let rpm = require_positive("rpm", rpm)?;
        Ok(self.lead_mm * rpm)
    }

    /// Linear speed of the nut in **mm/s** for a given screw speed in
    /// rev/min.
    ///
    /// Identical to [`linear_speed_mm_per_min`](Self::linear_speed_mm_per_min)
    /// divided by 60.
    ///
    /// # Errors
    ///
    /// Returns an error if `rpm` is non-finite or non-positive.
    pub fn linear_speed_mm_per_s(&self, rpm: f64) -> Result<f64, LeadScrewError> {
        Ok(self.linear_speed_mm_per_min(rpm)? / 60.0)
    }

    /// Axial thrust produced by an applied screw torque.
    ///
    /// `F = 2 * pi * eta * T / lead`.
    ///
    /// Consistent units are the caller's responsibility. With torque
    /// `T` in **N·mm** and `lead` in **mm/rev**, the thrust comes out in
    /// **newtons**, because one revolution does `2 * pi * T` (N·mm) of
    /// input work and advances the load by `lead` (mm), so the ideal
    /// force is `2 * pi * T / lead` newtons, scaled by the efficiency
    /// `eta` (fraction of input work delivered as axial work).
    ///
    /// # Errors
    ///
    /// Returns [`LeadScrewError::NotPositive`] /
    /// [`LeadScrewError::NotFinite`] for a bad `torque`, and
    /// [`LeadScrewError::EfficiencyOutOfRange`] for an `eta` outside
    /// `(0, 1]`.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_leadscrew::LeadScrew;
    /// let screw = LeadScrew::new(2.0, 8.0).unwrap();
    /// // Ideal (eta = 1): F = 2*pi*1.0 / 2.0 = pi newtons.
    /// let f = screw.thrust_n(1.0, 1.0).unwrap();
    /// assert!((f - std::f64::consts::PI).abs() < 1e-9);
    /// ```
    pub fn thrust_n(&self, torque_n_mm: f64, eta: f64) -> Result<f64, LeadScrewError> {
        let torque_n_mm = require_positive("torque_n_mm", torque_n_mm)?;
        let eta = crate::error::require_efficiency("eta", eta)?;
        Ok(2.0 * PI * eta * torque_n_mm / self.lead_mm)
    }

    /// Screw torque required to hold or drive a target axial thrust.
    ///
    /// The inverse of [`thrust_n`](Self::thrust_n):
    /// `T = F * lead / (2 * pi * eta)`. Same unit convention — with
    /// `force` in newtons and `lead` in mm the torque is in N·mm.
    ///
    /// # Errors
    ///
    /// Returns an error for a non-positive / non-finite `force` or an
    /// `eta` outside `(0, 1]`.
    pub fn torque_for_thrust_n_mm(&self, force_n: f64, eta: f64) -> Result<f64, LeadScrewError> {
        let force_n = require_positive("force_n", force_n)?;
        let eta = crate::error::require_efficiency("eta", eta)?;
        Ok(force_n * self.lead_mm / (2.0 * PI * eta))
    }

    /// Linear resolution — distance travelled per microstep.
    ///
    /// `resolution = lead / microsteps_per_rev`. A finer microstepping
    /// (more steps per revolution) yields a smaller — i.e. finer —
    /// linear increment.
    ///
    /// # Errors
    ///
    /// Returns [`LeadScrewError::ZeroMicrosteps`] if
    /// `microsteps_per_rev` is `0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_leadscrew::LeadScrew;
    /// let screw = LeadScrew::new(8.0, 8.0).unwrap();
    /// // A 200-step motor at 16x microstepping = 3200 microsteps/rev.
    /// // 8 mm / 3200 = 0.0025 mm per microstep.
    /// let r = screw.resolution_mm(3200).unwrap();
    /// assert!((r - 0.0025).abs() < 1e-12);
    /// ```
    pub fn resolution_mm(&self, microsteps_per_rev: u32) -> Result<f64, LeadScrewError> {
        if microsteps_per_rev == 0 {
            return Err(LeadScrewError::ZeroMicrosteps(0));
        }
        Ok(self.lead_mm / f64::from(microsteps_per_rev))
    }

    /// Lead angle (helix angle of the thread at the pitch diameter), in
    /// **radians**.
    ///
    /// `lambda = atan(lead / (pi * d_m))`. This is the angle of the
    /// thread helix measured from the plane perpendicular to the screw
    /// axis; it governs both the mechanical advantage and whether the
    /// screw can back-drive.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_leadscrew::LeadScrew;
    /// let screw = LeadScrew::new(2.0, 8.0).unwrap();
    /// let lambda = screw.lead_angle_rad();
    /// // atan(2 / (pi*8)) ~= 0.07943 rad.
    /// assert!((lambda - (2.0f64 / (std::f64::consts::PI * 8.0)).atan()).abs() < 1e-12);
    /// ```
    pub fn lead_angle_rad(&self) -> f64 {
        (self.lead_mm / (PI * self.pitch_diameter_mm)).atan()
    }

    /// Lead angle in **degrees** — convenience wrapper around
    /// [`lead_angle_rad`](Self::lead_angle_rad).
    pub fn lead_angle_deg(&self) -> f64 {
        self.lead_angle_rad().to_degrees()
    }
}
