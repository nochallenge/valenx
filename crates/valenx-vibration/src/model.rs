//! The single-degree-of-freedom (SDOF) mass-spring-damper model.
//!
//! A point mass `m` (kg) is restrained by a linear spring of stiffness
//! `k` (N/m) and a viscous damper of coefficient `c` (N·s/m). Its free
//! motion `x(t)` obeys the linear ODE
//!
//! ```text
//! m x'' + c x' + k x = 0
//! ```
//!
//! From the three physical constants this module derives the standard
//! modal descriptors:
//!
//! - **Undamped natural frequency** `wn = sqrt(k/m)` (rad/s).
//! - **Damping ratio** `zeta = c / (2*sqrt(k*m)) = c / c_crit`, the
//!   damping as a fraction of the critical value
//!   `c_crit = 2*sqrt(k*m)`.
//! - **Damped natural frequency** `wd = wn*sqrt(1 - zeta^2)` (rad/s) —
//!   the frequency at which an *underdamped* system actually oscillates.
//!   Only real for `zeta < 1`.
//!
//! These follow the standard treatment in Rao, *Mechanical Vibrations*,
//! and Thomson, *Theory of Vibration with Applications*.

use crate::error::VibrationError;
use serde::{Deserialize, Serialize};

/// The damping regime of an SDOF system, set by its damping ratio
/// `zeta`.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum DampingRegime {
    /// `zeta = 0`: undamped — the system oscillates forever at `wn`.
    Undamped,
    /// `0 < zeta < 1`: underdamped — decaying oscillation at `wd`.
    Underdamped,
    /// `zeta = 1`: critically damped — the fastest non-oscillating
    /// return to equilibrium (a repeated real root).
    CriticallyDamped,
    /// `zeta > 1`: overdamped — a slow, non-oscillating return (two
    /// distinct real roots).
    Overdamped,
}

/// A linear single-degree-of-freedom mass-spring-damper system.
///
/// Construct with [`SdofSystem::new`], which validates that `m > 0`,
/// `k > 0` and `c >= 0`. All derived modal quantities are then exact
/// closed-form functions of these three numbers.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SdofSystem {
    mass_kg: f64,
    stiffness_n_per_m: f64,
    damping_n_s_per_m: f64,
}

impl SdofSystem {
    /// Build a system from mass `m` (kg), stiffness `k` (N/m) and
    /// viscous damping coefficient `c` (N·s/m).
    ///
    /// # Errors
    ///
    /// Returns [`VibrationError::BadParameter`] if `m` or `k` is not
    /// strictly positive, if `c` is negative, or if any input is not a
    /// finite number.
    pub fn new(
        mass_kg: f64,
        stiffness_n_per_m: f64,
        damping_n_s_per_m: f64,
    ) -> Result<Self, VibrationError> {
        require_positive("mass_kg", mass_kg)?;
        require_positive("stiffness_n_per_m", stiffness_n_per_m)?;
        require_non_negative("damping_n_s_per_m", damping_n_s_per_m)?;
        Ok(Self {
            mass_kg,
            stiffness_n_per_m,
            damping_n_s_per_m,
        })
    }

    /// Convenience constructor from `wn` (rad/s) and `zeta` with unit
    /// mass (`m = 1` kg).
    ///
    /// Picks `k = wn^2` and `c = 2*zeta*wn` so the resulting system has
    /// exactly the requested natural frequency and damping ratio. Handy
    /// for tests and modal "what-if" studies where only the ratios
    /// matter.
    ///
    /// # Errors
    ///
    /// Returns [`VibrationError::BadParameter`] if `wn <= 0`, if
    /// `zeta < 0`, or if either is non-finite.
    pub fn from_modal(natural_freq_rad_s: f64, damping_ratio: f64) -> Result<Self, VibrationError> {
        require_positive("natural_freq_rad_s", natural_freq_rad_s)?;
        require_non_negative("damping_ratio", damping_ratio)?;
        let k = natural_freq_rad_s * natural_freq_rad_s;
        let c = 2.0 * damping_ratio * natural_freq_rad_s;
        Self::new(1.0, k, c)
    }

    /// Mass `m` (kg).
    pub fn mass_kg(&self) -> f64 {
        self.mass_kg
    }

    /// Spring stiffness `k` (N/m).
    pub fn stiffness_n_per_m(&self) -> f64 {
        self.stiffness_n_per_m
    }

    /// Viscous damping coefficient `c` (N·s/m).
    pub fn damping_n_s_per_m(&self) -> f64 {
        self.damping_n_s_per_m
    }

    /// Undamped natural frequency `wn = sqrt(k/m)` (rad/s).
    pub fn natural_freq_rad_s(&self) -> f64 {
        (self.stiffness_n_per_m / self.mass_kg).sqrt()
    }

    /// Undamped natural frequency in hertz, `f_n = wn / (2*pi)`.
    pub fn natural_freq_hz(&self) -> f64 {
        self.natural_freq_rad_s() / std::f64::consts::TAU
    }

    /// Undamped natural period `T_n = 2*pi / wn` (s).
    pub fn natural_period_s(&self) -> f64 {
        std::f64::consts::TAU / self.natural_freq_rad_s()
    }

    /// Critical damping coefficient `c_crit = 2*sqrt(k*m)` (N·s/m).
    ///
    /// This is the smallest `c` for which the system no longer
    /// oscillates; the damping ratio is `zeta = c / c_crit`.
    pub fn critical_damping(&self) -> f64 {
        2.0 * (self.stiffness_n_per_m * self.mass_kg).sqrt()
    }

    /// Damping ratio `zeta = c / (2*sqrt(k*m))` (dimensionless).
    pub fn damping_ratio(&self) -> f64 {
        self.damping_n_s_per_m / self.critical_damping()
    }

    /// Classify the [`DampingRegime`] from the damping ratio.
    ///
    /// Because `zeta` is a computed float, the boundaries (`zeta = 0`
    /// for undamped, `zeta = 1` for critical) are compared with a small
    /// relative tolerance rather than exact equality.
    pub fn regime(&self) -> DampingRegime {
        let zeta = self.damping_ratio();
        // Absolute tolerance on the dimensionless ratio.
        const TOL: f64 = 1e-12;
        if zeta <= TOL {
            DampingRegime::Undamped
        } else if (zeta - 1.0).abs() <= 1e-9 {
            DampingRegime::CriticallyDamped
        } else if zeta < 1.0 {
            DampingRegime::Underdamped
        } else {
            DampingRegime::Overdamped
        }
    }

    /// Damped natural frequency `wd = wn*sqrt(1 - zeta^2)` (rad/s).
    ///
    /// This is the actual oscillation frequency of the decaying free
    /// response.
    ///
    /// # Errors
    ///
    /// Returns [`VibrationError::NotApplicable`] unless the system is
    /// underdamped (`zeta < 1`); a critically- or over-damped system
    /// does not oscillate, so `wd` is not defined.
    pub fn damped_freq_rad_s(&self) -> Result<f64, VibrationError> {
        let zeta = self.damping_ratio();
        if zeta >= 1.0 {
            return Err(VibrationError::NotApplicable(format!(
                "damped frequency requires zeta < 1 (underdamped); got zeta = {zeta}"
            )));
        }
        Ok(self.natural_freq_rad_s() * (1.0 - zeta * zeta).sqrt())
    }

    /// Damped natural period `T_d = 2*pi / wd` (s).
    ///
    /// # Errors
    ///
    /// Propagates [`VibrationError::NotApplicable`] from
    /// [`damped_freq_rad_s`](Self::damped_freq_rad_s) for non-underdamped
    /// systems.
    pub fn damped_period_s(&self) -> Result<f64, VibrationError> {
        Ok(std::f64::consts::TAU / self.damped_freq_rad_s()?)
    }
}

/// Validate that a named parameter is finite and strictly positive.
fn require_positive(name: &'static str, value: f64) -> Result<(), VibrationError> {
    if !value.is_finite() {
        return Err(VibrationError::BadParameter {
            name,
            reason: format!("must be a finite number, got {value}"),
        });
    }
    if value <= 0.0 {
        return Err(VibrationError::BadParameter {
            name,
            reason: format!("must be strictly positive, got {value}"),
        });
    }
    Ok(())
}

/// Validate that a named parameter is finite and non-negative.
fn require_non_negative(name: &'static str, value: f64) -> Result<(), VibrationError> {
    if !value.is_finite() {
        return Err(VibrationError::BadParameter {
            name,
            reason: format!("must be a finite number, got {value}"),
        });
    }
    if value < 0.0 {
        return Err(VibrationError::BadParameter {
            name,
            reason: format!("must be non-negative, got {value}"),
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tolerance for floating-point comparisons in this module's tests.
    const EPS: f64 = 1e-9;

    #[test]
    fn natural_frequency_is_sqrt_k_over_m() {
        // k = 400 N/m, m = 1 kg  =>  wn = sqrt(400) = 20 rad/s.
        let sys = SdofSystem::new(1.0, 400.0, 0.0).expect("valid");
        assert!((sys.natural_freq_rad_s() - 20.0).abs() < EPS);

        // k = 100, m = 4  =>  wn = sqrt(25) = 5 rad/s.
        let sys = SdofSystem::new(4.0, 100.0, 0.0).expect("valid");
        assert!((sys.natural_freq_rad_s() - 5.0).abs() < EPS);
    }

    #[test]
    fn natural_freq_hz_and_period_are_consistent() {
        let sys = SdofSystem::new(1.0, 400.0, 0.0).expect("valid");
        // f = wn / 2pi; T = 1/f = 2pi/wn.
        let expected_hz = 20.0 / std::f64::consts::TAU;
        assert!((sys.natural_freq_hz() - expected_hz).abs() < EPS);
        assert!((sys.natural_period_s() - 1.0 / expected_hz).abs() < EPS);
        // T = 2pi/wn directly.
        assert!((sys.natural_period_s() - std::f64::consts::TAU / 20.0).abs() < EPS);
    }

    #[test]
    fn critical_damping_and_ratio() {
        // m = 2, k = 8  =>  c_crit = 2*sqrt(16) = 8.
        let sys = SdofSystem::new(2.0, 8.0, 0.0).expect("valid");
        assert!((sys.critical_damping() - 8.0).abs() < EPS);

        // c = 4 of c_crit = 8  =>  zeta = 0.5.
        let sys = SdofSystem::new(2.0, 8.0, 4.0).expect("valid");
        assert!((sys.damping_ratio() - 0.5).abs() < EPS);
    }

    #[test]
    fn zeta_equals_one_is_critically_damped() {
        // Set c exactly to c_crit = 2*sqrt(k*m).
        let m = 3.0_f64;
        let k = 27.0_f64;
        let c_crit = 2.0 * (k * m).sqrt();
        let sys = SdofSystem::new(m, k, c_crit).expect("valid");
        assert!((sys.damping_ratio() - 1.0).abs() < EPS);
        assert_eq!(sys.regime(), DampingRegime::CriticallyDamped);
        // wd is not defined for a critically-damped system.
        assert!(sys.damped_freq_rad_s().is_err());
    }

    #[test]
    fn regime_classification_across_zeta() {
        let undamped = SdofSystem::from_modal(10.0, 0.0).expect("valid");
        assert_eq!(undamped.regime(), DampingRegime::Undamped);

        let under = SdofSystem::from_modal(10.0, 0.25).expect("valid");
        assert_eq!(under.regime(), DampingRegime::Underdamped);

        let over = SdofSystem::from_modal(10.0, 2.0).expect("valid");
        assert_eq!(over.regime(), DampingRegime::Overdamped);
    }

    #[test]
    fn damped_frequency_for_underdamped() {
        // wn = 10, zeta = 0.6  =>  wd = 10*sqrt(1 - 0.36) = 10*0.8 = 8.
        let sys = SdofSystem::from_modal(10.0, 0.6).expect("valid");
        let wd = sys.damped_freq_rad_s().expect("underdamped");
        assert!((wd - 8.0).abs() < EPS);
        // Damped period matches 2pi/wd.
        assert!((sys.damped_period_s().expect("ud") - std::f64::consts::TAU / 8.0).abs() < EPS);
    }

    #[test]
    fn damped_frequency_below_natural_for_light_damping() {
        // For any 0 < zeta < 1, wd < wn but close for small zeta.
        let sys = SdofSystem::from_modal(50.0, 0.02).expect("valid");
        let wd = sys.damped_freq_rad_s().expect("ud");
        assert!(wd < sys.natural_freq_rad_s());
        // sqrt(1 - 0.0004) ~ 0.9998, so within 0.05% of wn.
        assert!((wd / sys.natural_freq_rad_s() - 1.0).abs() < 5e-4);
    }

    #[test]
    fn overdamped_has_no_damped_frequency() {
        let sys = SdofSystem::from_modal(10.0, 1.5).expect("valid");
        let err = sys.damped_freq_rad_s().expect_err("overdamped");
        assert_eq!(err.code(), "vibration.not_applicable");
    }

    #[test]
    fn from_modal_round_trips_wn_and_zeta() {
        let sys = SdofSystem::from_modal(7.5, 0.3).expect("valid");
        assert!((sys.natural_freq_rad_s() - 7.5).abs() < EPS);
        assert!((sys.damping_ratio() - 0.3).abs() < EPS);
    }

    #[test]
    fn rejects_non_positive_mass_and_stiffness() {
        assert!(SdofSystem::new(0.0, 1.0, 0.0).is_err());
        assert!(SdofSystem::new(-1.0, 1.0, 0.0).is_err());
        assert!(SdofSystem::new(1.0, 0.0, 0.0).is_err());
        assert!(SdofSystem::new(1.0, -5.0, 0.0).is_err());
    }

    #[test]
    fn rejects_negative_damping_and_non_finite() {
        assert!(SdofSystem::new(1.0, 1.0, -0.1).is_err());
        assert!(SdofSystem::new(f64::NAN, 1.0, 0.0).is_err());
        assert!(SdofSystem::new(1.0, f64::INFINITY, 0.0).is_err());
    }

    #[test]
    fn bad_parameter_reports_name_and_category() {
        let err = SdofSystem::new(-2.0, 1.0, 0.0).expect_err("bad mass");
        assert_eq!(err.code(), "vibration.bad_parameter");
        assert_eq!(err.category(), crate::error::ErrorCategory::Input);
        // The offending parameter name appears in the message.
        assert!(format!("{err}").contains("mass_kg"));
    }
}
