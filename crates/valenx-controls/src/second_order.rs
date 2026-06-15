//! Standard second-order system and its step-response metrics.
//!
//! A canonical second-order system has the closed-loop transfer function
//!
//! ```text
//!                 wn^2
//! G(s) = -----------------------
//!         s^2 + 2 zeta wn s + wn^2
//! ```
//!
//! parameterised by the **undamped natural frequency** `wn > 0` (rad/s)
//! and the **damping ratio** `zeta >= 0` (dimensionless). Almost every
//! textbook performance figure for a unit-step input is a closed form in
//! these two numbers; this module computes them.
//!
//! ## Models (all standard, Ogata / Franklin closed forms)
//!
//! For the **underdamped** case `0 <= zeta < 1` the damped frequency is
//! `wd = wn * sqrt(1 - zeta^2)` and:
//!
//! - **Percent overshoot** `Mp = exp(-pi * zeta / sqrt(1 - zeta^2))`,
//!   reported as a fraction (multiply by 100 for a percentage). It
//!   depends on `zeta` *only*, not on `wn`.
//! - **Peak time** `tp = pi / wd` — when the first (largest) overshoot
//!   occurs.
//! - **Settling time** `ts ~ 4 / (zeta * wn)` for the 2% band (a `3 /
//!   (zeta*wn)` 5% variant is also provided). Governed by the real part
//!   `zeta * wn` of the poles.
//! - **Rise time** `tr = (pi - acos(zeta)) / wd`, the 0-100% rise for an
//!   underdamped system.
//!
//! For `zeta >= 1` (critically damped / overdamped) there is no
//! oscillation: the overshoot is exactly zero, and the peak / rise-time
//! formulae above are singular (`wd = 0`), so those accessors return a
//! [`DomainError`](crate::error::ControlsError::DomainError).

use serde::{Deserialize, Serialize};

use crate::error::{ControlsError, Result};

/// A canonical second-order system parameterised by its natural
/// frequency and damping ratio.
///
/// Construct with [`SecondOrder::new`], which validates that `wn` is a
/// finite positive number and `zeta` is a finite non-negative number.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SecondOrder {
    /// Undamped natural frequency `wn` (rad/s), strictly positive.
    pub wn: f64,
    /// Damping ratio `zeta` (dimensionless), non-negative. `< 1` is
    /// underdamped, `== 1` critically damped, `> 1` overdamped.
    pub zeta: f64,
}

/// The qualitative damping regime of a second-order system.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum DampingRegime {
    /// `zeta == 0`: poles on the imaginary axis, sustained oscillation.
    Undamped,
    /// `0 < zeta < 1`: complex-conjugate poles, decaying oscillation.
    Underdamped,
    /// `zeta == 1`: a repeated real pole, fastest non-oscillatory.
    CriticallyDamped,
    /// `zeta > 1`: two distinct real poles, no oscillation.
    Overdamped,
}

/// A bundle of unit-step performance metrics for an underdamped system.
///
/// Produced by [`SecondOrder::step_metrics`]. Every field is in SI
/// (seconds, or a dimensionless fraction for the overshoot).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct StepMetrics {
    /// Peak (fractional) overshoot `Mp = exp(-pi*zeta/sqrt(1-zeta^2))`.
    /// Multiply by 100 for a percentage. Zero for `zeta >= 1`.
    pub overshoot: f64,
    /// Peak time `tp = pi / wd` (s) — when the first overshoot occurs.
    pub peak_time: f64,
    /// 2% settling time `ts = 4 / (zeta*wn)` (s).
    pub settling_time: f64,
    /// 0-100% rise time `tr = (pi - acos(zeta)) / wd` (s).
    pub rise_time: f64,
    /// Damped natural frequency `wd = wn*sqrt(1-zeta^2)` (rad/s).
    pub damped_frequency: f64,
}

impl SecondOrder {
    /// Construct a second-order system from its natural frequency and
    /// damping ratio.
    ///
    /// # Errors
    ///
    /// Returns [`ControlsError::InvalidParameter`] if `wn` is not finite
    /// or not strictly positive, or if `zeta` is not finite or is
    /// negative.
    pub fn new(wn: f64, zeta: f64) -> Result<Self> {
        if !wn.is_finite() || wn <= 0.0 {
            return Err(ControlsError::invalid(
                "wn",
                "natural frequency must be finite and > 0",
            ));
        }
        if !zeta.is_finite() || zeta < 0.0 {
            return Err(ControlsError::invalid(
                "zeta",
                "damping ratio must be finite and >= 0",
            ));
        }
        Ok(Self { wn, zeta })
    }

    /// The qualitative [`DampingRegime`] of this system.
    pub fn regime(&self) -> DampingRegime {
        if self.zeta == 0.0 {
            DampingRegime::Undamped
        } else if self.zeta < 1.0 {
            DampingRegime::Underdamped
        } else if self.zeta == 1.0 {
            DampingRegime::CriticallyDamped
        } else {
            DampingRegime::Overdamped
        }
    }

    /// `true` when `zeta < 1` (the oscillatory case).
    pub fn is_underdamped(&self) -> bool {
        self.zeta < 1.0
    }

    /// Damped natural frequency `wd = wn * sqrt(1 - zeta^2)` (rad/s).
    ///
    /// Real and positive only for `zeta < 1`; it is zero at `zeta == 1`
    /// and undefined (the square root of a negative) for `zeta > 1`.
    ///
    /// # Errors
    ///
    /// Returns [`ControlsError::DomainError`] for `zeta >= 1`, where the
    /// damped frequency is not a positive real and the dependent
    /// peak/rise-time formulae are singular.
    pub fn damped_frequency(&self) -> Result<f64> {
        if self.zeta >= 1.0 {
            return Err(ControlsError::domain(
                "damped frequency is not positive-real for zeta >= 1",
            ));
        }
        Ok(self.wn * (1.0 - self.zeta * self.zeta).sqrt())
    }

    /// Fractional peak overshoot of the unit-step response.
    ///
    /// `Mp = exp(-pi * zeta / sqrt(1 - zeta^2))` for `0 <= zeta < 1`;
    /// **exactly `0.0`** for `zeta >= 1` (a critically- or over-damped
    /// system never overshoots). Multiply by 100 for a percentage.
    ///
    /// This is total and infallible: the `zeta >= 1` case has the
    /// physically correct value of zero rather than a domain error,
    /// because "no overshoot" is a meaningful answer there.
    pub fn overshoot(&self) -> f64 {
        if self.zeta >= 1.0 {
            return 0.0;
        }
        if self.zeta == 0.0 {
            // Undamped: sustained oscillation, the step response peaks at
            // exactly twice the final value -> 100% overshoot.
            return 1.0;
        }
        let ratio = self.zeta / (1.0 - self.zeta * self.zeta).sqrt();
        (-std::f64::consts::PI * ratio).exp()
    }

    /// Peak overshoot expressed as a percentage (`100 * overshoot`).
    pub fn percent_overshoot(&self) -> f64 {
        100.0 * self.overshoot()
    }

    /// Peak time `tp = pi / wd` (s) — when the first overshoot occurs.
    ///
    /// # Errors
    ///
    /// Returns [`ControlsError::DomainError`] for `zeta >= 1`: a
    /// non-oscillatory system has no overshoot peak, and the formula is
    /// singular (`wd = 0`).
    pub fn peak_time(&self) -> Result<f64> {
        let wd = self.damped_frequency()?;
        Ok(std::f64::consts::PI / wd)
    }

    /// 2% settling time `ts ~ 4 / (zeta * wn)` (s).
    ///
    /// The standard engineering estimate for the time after which the
    /// response stays within ±2% of its final value. Governed by the
    /// real part `zeta * wn` of the dominant poles.
    ///
    /// # Errors
    ///
    /// Returns [`ControlsError::DomainError`] for `zeta == 0`, where the
    /// undamped response never decays and the settling time is infinite.
    pub fn settling_time_2pct(&self) -> Result<f64> {
        if self.zeta == 0.0 {
            return Err(ControlsError::domain(
                "settling time is infinite for an undamped (zeta = 0) system",
            ));
        }
        Ok(4.0 / (self.zeta * self.wn))
    }

    /// 5% settling time `ts ~ 3 / (zeta * wn)` (s).
    ///
    /// The same estimate for the looser ±5% band.
    ///
    /// # Errors
    ///
    /// Returns [`ControlsError::DomainError`] for `zeta == 0` (see
    /// [`settling_time_2pct`](Self::settling_time_2pct)).
    pub fn settling_time_5pct(&self) -> Result<f64> {
        if self.zeta == 0.0 {
            return Err(ControlsError::domain(
                "settling time is infinite for an undamped (zeta = 0) system",
            ));
        }
        Ok(3.0 / (self.zeta * self.wn))
    }

    /// 0-100% rise time `tr = (pi - acos(zeta)) / wd` (s).
    ///
    /// The textbook rise time for an underdamped system — the time for
    /// the step response to climb from 0% to 100% of its final value the
    /// first time.
    ///
    /// # Errors
    ///
    /// Returns [`ControlsError::DomainError`] for `zeta >= 1`, where `wd
    /// = 0` and the underdamped rise-time formula does not apply.
    pub fn rise_time(&self) -> Result<f64> {
        let wd = self.damped_frequency()?;
        // acos(zeta) is the pole angle from the negative real axis; the
        // numerator (pi - acos(zeta)) is the phase swept to first reach
        // the final value.
        Ok((std::f64::consts::PI - self.zeta.acos()) / wd)
    }

    /// All underdamped step-response metrics in one [`StepMetrics`]
    /// bundle.
    ///
    /// # Errors
    ///
    /// Returns [`ControlsError::DomainError`] for `zeta >= 1` (the
    /// oscillatory metrics are undefined) or `zeta == 0` (infinite
    /// settling time). For those regimes query the individual,
    /// well-defined accessors instead — e.g. [`overshoot`](Self::overshoot),
    /// which is `0.0` for `zeta >= 1`.
    pub fn step_metrics(&self) -> Result<StepMetrics> {
        let damped_frequency = self.damped_frequency()?;
        Ok(StepMetrics {
            overshoot: self.overshoot(),
            peak_time: self.peak_time()?,
            settling_time: self.settling_time_2pct()?,
            rise_time: self.rise_time()?,
            damped_frequency,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tolerance for floating-point ground-truth comparisons.
    const EPS: f64 = 1e-9;

    #[test]
    fn rejects_non_physical_parameters() {
        assert!(SecondOrder::new(0.0, 0.5).is_err());
        assert!(SecondOrder::new(-1.0, 0.5).is_err());
        assert!(SecondOrder::new(f64::NAN, 0.5).is_err());
        assert!(SecondOrder::new(f64::INFINITY, 0.5).is_err());
        assert!(SecondOrder::new(1.0, -0.1).is_err());
        assert!(SecondOrder::new(1.0, f64::NAN).is_err());
        // A normal underdamped system is accepted.
        assert!(SecondOrder::new(1.0, 0.5).is_ok());
    }

    #[test]
    fn regime_classification() {
        assert_eq!(
            SecondOrder::new(2.0, 0.0).unwrap().regime(),
            DampingRegime::Undamped
        );
        assert_eq!(
            SecondOrder::new(2.0, 0.5).unwrap().regime(),
            DampingRegime::Underdamped
        );
        assert_eq!(
            SecondOrder::new(2.0, 1.0).unwrap().regime(),
            DampingRegime::CriticallyDamped
        );
        assert_eq!(
            SecondOrder::new(2.0, 2.0).unwrap().regime(),
            DampingRegime::Overdamped
        );
    }

    #[test]
    fn overshoot_is_zero_for_critical_and_overdamped() {
        // VALIDATE: zeta >= 1 -> 0 overshoot.
        let critical = SecondOrder::new(3.0, 1.0).unwrap();
        let over = SecondOrder::new(3.0, 2.5).unwrap();
        assert!(
            critical.overshoot().abs() < EPS,
            "Mp = {}",
            critical.overshoot()
        );
        assert!(over.overshoot().abs() < EPS, "Mp = {}", over.overshoot());
        // The oscillatory metrics are domain errors there.
        assert!(critical.peak_time().is_err());
        assert!(over.step_metrics().is_err());
    }

    #[test]
    fn overshoot_increases_as_damping_falls() {
        // VALIDATE: lower zeta -> more overshoot. wn is irrelevant to Mp.
        let zetas = [0.9_f64, 0.7, 0.5, 0.3, 0.1];
        let mps: Vec<f64> = zetas
            .iter()
            .map(|&z| SecondOrder::new(5.0, z).unwrap().overshoot())
            .collect();
        for pair in mps.windows(2) {
            assert!(
                pair[1] > pair[0],
                "overshoot should rise as zeta falls: {} then {}",
                pair[0],
                pair[1]
            );
        }
        // Overshoot depends only on zeta, not wn.
        let a = SecondOrder::new(1.0, 0.4).unwrap().overshoot();
        let b = SecondOrder::new(100.0, 0.4).unwrap().overshoot();
        assert!((a - b).abs() < EPS, "{a} vs {b}");
    }

    #[test]
    fn overshoot_half_damping_is_about_sixteen_percent() {
        // VALIDATE: Mp(0.5) ~ 16%. Closed form:
        // exp(-pi*0.5/sqrt(0.75)) = 0.1630335348...
        let mp = SecondOrder::new(4.0, 0.5).unwrap().overshoot();
        assert!((mp - 0.163_033_534_8).abs() < 1e-9, "Mp(0.5) = {mp}");
        // Sanity: percent form is ~16.3.
        let pct = SecondOrder::new(4.0, 0.5).unwrap().percent_overshoot();
        assert!((pct - 16.303_353_5).abs() < 1e-6, "{pct}%");
    }

    #[test]
    fn undamped_overshoots_one_hundred_percent() {
        // zeta = 0: the step response peaks at 2x final value -> Mp = 1.
        let mp = SecondOrder::new(2.0, 0.0).unwrap().overshoot();
        assert!((mp - 1.0).abs() < EPS, "Mp = {mp}");
    }

    #[test]
    fn settling_time_falls_with_zeta_wn_product() {
        // VALIDATE: settling falls as zeta*wn grows.
        let slow = SecondOrder::new(1.0, 0.2)
            .unwrap()
            .settling_time_2pct()
            .unwrap();
        let mid = SecondOrder::new(2.0, 0.4)
            .unwrap()
            .settling_time_2pct()
            .unwrap();
        let fast = SecondOrder::new(5.0, 0.7)
            .unwrap()
            .settling_time_2pct()
            .unwrap();
        assert!(slow > mid && mid > fast, "{slow} > {mid} > {fast}");

        // Exact closed form: ts = 4 / (zeta*wn).
        let s = SecondOrder::new(2.0, 0.5).unwrap();
        assert!((s.settling_time_2pct().unwrap() - 4.0 / 1.0).abs() < EPS);
        assert!((s.settling_time_5pct().unwrap() - 3.0 / 1.0).abs() < EPS);
        // Undamped never settles.
        assert!(SecondOrder::new(2.0, 0.0)
            .unwrap()
            .settling_time_2pct()
            .is_err());
    }

    #[test]
    fn peak_time_matches_pi_over_wd() {
        // wn = 2, zeta = 0.5 -> wd = 2*sqrt(0.75) = 1.7320508; tp = pi/wd.
        let s = SecondOrder::new(2.0, 0.5).unwrap();
        let wd = s.damped_frequency().unwrap();
        assert!((wd - 1.732_050_808).abs() < 1e-6, "wd = {wd}");
        let tp = s.peak_time().unwrap();
        assert!((tp - std::f64::consts::PI / wd).abs() < EPS, "tp = {tp}");
        // Critically damped -> no peak.
        assert!(SecondOrder::new(2.0, 1.0).unwrap().peak_time().is_err());
    }

    #[test]
    fn rise_time_closed_form() {
        // tr = (pi - acos(zeta)) / wd. For zeta = 0.5, acos(0.5) = pi/3.
        let s = SecondOrder::new(2.0, 0.5).unwrap();
        let wd = s.damped_frequency().unwrap();
        let expected = (std::f64::consts::PI - (0.5_f64).acos()) / wd;
        assert!((s.rise_time().unwrap() - expected).abs() < EPS);
        // Rise time is positive and smaller than peak time for an
        // underdamped system.
        assert!(s.rise_time().unwrap() > 0.0);
        assert!(s.rise_time().unwrap() < s.peak_time().unwrap());
    }

    #[test]
    fn step_metrics_bundle_is_consistent_with_accessors() {
        let s = SecondOrder::new(3.0, 0.6).unwrap();
        let m = s.step_metrics().unwrap();
        assert!((m.overshoot - s.overshoot()).abs() < EPS);
        assert!((m.peak_time - s.peak_time().unwrap()).abs() < EPS);
        assert!((m.settling_time - s.settling_time_2pct().unwrap()).abs() < EPS);
        assert!((m.rise_time - s.rise_time().unwrap()).abs() < EPS);
        assert!((m.damped_frequency - s.damped_frequency().unwrap()).abs() < EPS);
    }

    #[test]
    fn step_metrics_round_trips_through_json() {
        let m = SecondOrder::new(3.0, 0.6).unwrap().step_metrics().unwrap();
        let json = serde_json::to_string(&m).unwrap();
        let back: StepMetrics = serde_json::from_str(&json).unwrap();
        assert_eq!(m, back);
    }
}
