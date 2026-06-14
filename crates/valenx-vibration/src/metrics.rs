//! Steady-state resonance amplification and the logarithmic decrement.
//!
//! Two classic SDOF results that do not need a time-domain solve:
//!
//! ## Harmonic magnification factor
//!
//! Drive the system with a harmonic force `F0 cos(w t)`. The steady-state
//! displacement amplitude `X`, normalised by the static deflection
//! `X_static = F0/k`, is the **magnification factor** (Rao §3.4):
//!
//! ```text
//! M(r, zeta) = 1 / sqrt( (1 - r^2)^2 + (2*zeta*r)^2 )
//! ```
//!
//! where `r = w / wn` is the frequency ratio. For light damping `M`
//! peaks sharply just below `r = 1` (resonance); the exact peak is at
//! `r_peak = sqrt(1 - 2*zeta^2)` (real only for `zeta < 1/sqrt(2)`),
//! where `M_peak = 1 / (2*zeta*sqrt(1 - zeta^2))`. At resonance proper
//! (`r = 1`), `M = 1/(2*zeta)`.
//!
//! ## Logarithmic decrement
//!
//! For an *underdamped* free decay, the ratio of two displacement peaks
//! `n` cycles apart gives the **logarithmic decrement**
//!
//! ```text
//! delta = (1/n) * ln(x_i / x_(i+n))
//! ```
//!
//! which relates to the damping ratio by
//!
//! ```text
//! delta = 2*pi*zeta / sqrt(1 - zeta^2)   <=>   zeta = delta / sqrt(4*pi^2 + delta^2).
//! ```
//!
//! These are exact inverses, so a measured decay recovers `zeta`.

use crate::error::VibrationError;
use crate::model::SdofSystem;

/// Two-pi, reused below for the decrement relations.
const TWO_PI: f64 = std::f64::consts::TAU;

/// Steady-state magnification factor `M(r, zeta)` for a harmonic force.
///
/// `frequency_ratio` is `r = w / wn` (the forcing frequency over the
/// natural frequency). The result is the dimensionless ratio of dynamic
/// to static displacement amplitude,
/// `M = 1 / sqrt((1 - r^2)^2 + (2*zeta*r)^2)`.
///
/// # Errors
///
/// Returns [`VibrationError::BadParameter`] if `frequency_ratio` is
/// negative or non-finite. (At `r = 0` and `zeta = 0` the value is
/// exactly `1`; the undamped resonance `r = 1, zeta = 0` is `+inf`,
/// which is returned as-is rather than erroring.)
pub fn magnification_factor(
    system: &SdofSystem,
    frequency_ratio: f64,
) -> Result<f64, VibrationError> {
    if !frequency_ratio.is_finite() || frequency_ratio < 0.0 {
        return Err(VibrationError::BadParameter {
            name: "frequency_ratio",
            reason: format!("must be finite and non-negative, got {frequency_ratio}"),
        });
    }
    let zeta = system.damping_ratio();
    let r = frequency_ratio;
    let r2 = r * r;
    let denom_sq = (1.0 - r2) * (1.0 - r2) + (2.0 * zeta * r) * (2.0 * zeta * r);
    Ok(1.0 / denom_sq.sqrt())
}

/// The frequency ratio `r_peak = sqrt(1 - 2*zeta^2)` at which the
/// magnification factor is maximal.
///
/// # Errors
///
/// Returns [`VibrationError::NotApplicable`] when `zeta >= 1/sqrt(2)`
/// (`~0.707`): for such heavy damping the response has no interior peak
/// — `M` decreases monotonically from `r = 0`, so there is no resonance.
pub fn resonant_frequency_ratio(system: &SdofSystem) -> Result<f64, VibrationError> {
    let zeta = system.damping_ratio();
    let radicand = 1.0 - 2.0 * zeta * zeta;
    if radicand <= 0.0 {
        return Err(VibrationError::NotApplicable(format!(
            "no resonant peak for zeta >= 1/sqrt(2); got zeta = {zeta}"
        )));
    }
    Ok(radicand.sqrt())
}

/// The peak magnification factor `M_peak = 1 / (2*zeta*sqrt(1 - zeta^2))`.
///
/// This is the value of [`magnification_factor`] at
/// [`resonant_frequency_ratio`]. It is sometimes called the dynamic
/// amplification at resonance; for a lightly-damped system it is close
/// to the quality factor `Q = 1/(2*zeta)`.
///
/// # Errors
///
/// Returns [`VibrationError::NotApplicable`] when `zeta >= 1/sqrt(2)`
/// (no interior peak) or when `zeta = 0` (the peak is infinite).
pub fn peak_magnification(system: &SdofSystem) -> Result<f64, VibrationError> {
    let zeta = system.damping_ratio();
    if zeta <= 0.0 {
        return Err(VibrationError::NotApplicable(
            "peak magnification is infinite for an undamped system (zeta = 0)".to_string(),
        ));
    }
    // Reuse the peak-existence check.
    resonant_frequency_ratio(system)?;
    Ok(1.0 / (2.0 * zeta * (1.0 - zeta * zeta).sqrt()))
}

/// The magnification factor exactly at resonance `r = 1`,
/// `M(1) = 1/(2*zeta)`.
///
/// This is the half-power-point reference, not the true maximum (which
/// sits at [`resonant_frequency_ratio`] for `zeta > 0`).
///
/// # Errors
///
/// Returns [`VibrationError::NotApplicable`] for an undamped system,
/// where the response at `r = 1` is unbounded.
pub fn magnification_at_resonance(system: &SdofSystem) -> Result<f64, VibrationError> {
    let zeta = system.damping_ratio();
    if zeta <= 0.0 {
        return Err(VibrationError::NotApplicable(
            "response at r = 1 is unbounded for an undamped system".to_string(),
        ));
    }
    Ok(1.0 / (2.0 * zeta))
}

/// Logarithmic decrement `delta = 2*pi*zeta / sqrt(1 - zeta^2)` for a
/// given system.
///
/// # Errors
///
/// Returns [`VibrationError::NotApplicable`] unless the system is
/// underdamped (`zeta < 1`); the decrement is only defined for a
/// decaying oscillation.
pub fn logarithmic_decrement(system: &SdofSystem) -> Result<f64, VibrationError> {
    let zeta = system.damping_ratio();
    if zeta >= 1.0 {
        return Err(VibrationError::NotApplicable(format!(
            "logarithmic decrement requires zeta < 1 (underdamped); got zeta = {zeta}"
        )));
    }
    Ok(TWO_PI * zeta / (1.0 - zeta * zeta).sqrt())
}

/// Recover the damping ratio from a logarithmic decrement,
/// `zeta = delta / sqrt(4*pi^2 + delta^2)`.
///
/// This is the exact inverse of [`logarithmic_decrement`]: round-tripping
/// any `0 <= zeta < 1` through `delta` and back returns the original
/// ratio. The result is always in `[0, 1)`.
///
/// # Errors
///
/// Returns [`VibrationError::BadParameter`] if `delta` is negative or
/// non-finite.
pub fn damping_ratio_from_decrement(delta: f64) -> Result<f64, VibrationError> {
    if !delta.is_finite() || delta < 0.0 {
        return Err(VibrationError::BadParameter {
            name: "delta",
            reason: format!("logarithmic decrement must be finite and non-negative, got {delta}"),
        });
    }
    Ok(delta / (4.0 * std::f64::consts::PI * std::f64::consts::PI + delta * delta).sqrt())
}

/// Estimate the logarithmic decrement from two measured peak amplitudes
/// `n` cycles apart: `delta = (1/n) * ln(earlier / later)`.
///
/// `earlier_peak` and `later_peak` are the (positive) displacement
/// amplitudes at peaks separated by `cycles` full damped periods.
///
/// # Errors
///
/// - [`VibrationError::BadParameter`] if `cycles == 0`.
/// - [`VibrationError::InvalidDecay`] if either amplitude is not
///   strictly positive and finite, or if `later_peak > earlier_peak`
///   (which would imply growth, not decay).
pub fn decrement_from_peaks(
    earlier_peak: f64,
    later_peak: f64,
    cycles: u32,
) -> Result<f64, VibrationError> {
    if cycles == 0 {
        return Err(VibrationError::BadParameter {
            name: "cycles",
            reason: "number of cycles between peaks must be >= 1".to_string(),
        });
    }
    if !(earlier_peak.is_finite() && later_peak.is_finite())
        || earlier_peak <= 0.0
        || later_peak <= 0.0
    {
        return Err(VibrationError::InvalidDecay(format!(
            "peak amplitudes must be finite and strictly positive, got {earlier_peak} and {later_peak}"
        )));
    }
    if later_peak > earlier_peak {
        return Err(VibrationError::InvalidDecay(format!(
            "later peak {later_peak} exceeds earlier peak {earlier_peak}: not a decay"
        )));
    }
    Ok((earlier_peak / later_peak).ln() / f64::from(cycles))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tolerance for analytic comparisons.
    const EPS: f64 = 1e-9;

    #[test]
    fn static_magnification_is_one_at_zero_frequency() {
        // r = 0  =>  M = 1 for any damping.
        for zeta in [0.0_f64, 0.1, 0.5, 0.9] {
            let sys = SdofSystem::from_modal(10.0, zeta).expect("valid");
            assert!((magnification_factor(&sys, 0.0).expect("ok") - 1.0).abs() < EPS);
        }
    }

    #[test]
    fn magnification_at_resonance_is_half_over_zeta() {
        // r = 1  =>  M = 1/(2 zeta).
        let sys = SdofSystem::from_modal(10.0, 0.05).expect("valid");
        let m = magnification_factor(&sys, 1.0).expect("ok");
        assert!((m - 1.0 / (2.0 * 0.05)).abs() < EPS);
        // And the dedicated helper agrees.
        assert!((magnification_at_resonance(&sys).expect("ok") - 10.0).abs() < EPS);
    }

    #[test]
    fn resonance_peak_is_near_wn_for_light_damping() {
        // For small zeta, the peak ratio r_peak = sqrt(1 - 2 zeta^2) ~ 1,
        // and the peak amplitude is large (sharp resonance).
        let zeta = 0.02;
        let sys = SdofSystem::from_modal(30.0, zeta).expect("valid");
        let r_peak = resonant_frequency_ratio(&sys).expect("peak exists");
        assert!((r_peak - (1.0 - 2.0 * zeta * zeta).sqrt()).abs() < EPS);
        // Within 0.1% of r = 1 for this light damping.
        assert!((r_peak - 1.0).abs() < 1e-3);

        // The peak is genuinely a maximum: M(r_peak) >= M at r slightly
        // away on either side, and is large.
        let m_peak = magnification_factor(&sys, r_peak).expect("ok");
        let m_lo = magnification_factor(&sys, r_peak - 0.01).expect("ok");
        let m_hi = magnification_factor(&sys, r_peak + 0.01).expect("ok");
        assert!(m_peak >= m_lo && m_peak >= m_hi);
        assert!(m_peak > 20.0, "light damping should amplify strongly");
    }

    #[test]
    fn peak_magnification_matches_closed_form() {
        let zeta = 0.2;
        let sys = SdofSystem::from_modal(10.0, zeta).expect("valid");
        let expected = 1.0 / (2.0 * zeta * (1.0 - zeta * zeta).sqrt());
        assert!((peak_magnification(&sys).expect("ok") - expected).abs() < EPS);
        // It equals M evaluated at the peak ratio.
        let r_peak = resonant_frequency_ratio(&sys).expect("ok");
        let m_at_peak = magnification_factor(&sys, r_peak).expect("ok");
        assert!((peak_magnification(&sys).expect("ok") - m_at_peak).abs() < 1e-9);
    }

    #[test]
    fn no_resonant_peak_for_heavy_damping() {
        // zeta = 1/sqrt(2) is the threshold; above it, no interior peak.
        let sys = SdofSystem::from_modal(10.0, 0.8).expect("valid");
        assert!(resonant_frequency_ratio(&sys).is_err());
        assert!(peak_magnification(&sys).is_err());
    }

    #[test]
    fn undamped_resonance_helpers_error() {
        let sys = SdofSystem::from_modal(10.0, 0.0).expect("valid");
        assert!(magnification_at_resonance(&sys).is_err());
        assert!(peak_magnification(&sys).is_err());
    }

    #[test]
    fn log_decrement_relates_to_zeta() {
        // delta = 2 pi zeta / sqrt(1 - zeta^2).
        let zeta = 0.1;
        let sys = SdofSystem::from_modal(10.0, zeta).expect("valid");
        let expected = TWO_PI * zeta / (1.0 - zeta * zeta).sqrt();
        assert!((logarithmic_decrement(&sys).expect("ok") - expected).abs() < EPS);
    }

    #[test]
    fn decrement_and_zeta_are_exact_inverses() {
        // Round-trip several damping ratios through delta and back.
        for zeta in [0.01_f64, 0.05, 0.2, 0.5, 0.7, 0.99] {
            let sys = SdofSystem::from_modal(10.0, zeta).expect("valid");
            let delta = logarithmic_decrement(&sys).expect("ok");
            let recovered = damping_ratio_from_decrement(delta).expect("ok");
            assert!((recovered - zeta).abs() < 1e-12, "zeta = {zeta}");
        }
    }

    #[test]
    fn small_zeta_decrement_approx_two_pi_zeta() {
        // For very light damping, delta ~ 2 pi zeta (since sqrt(1-z^2)~1).
        let zeta = 1e-3;
        let sys = SdofSystem::from_modal(10.0, zeta).expect("valid");
        let delta = logarithmic_decrement(&sys).expect("ok");
        assert!((delta - TWO_PI * zeta).abs() < 1e-7);
    }

    #[test]
    fn decrement_undefined_for_critical_and_over_damped() {
        let crit = SdofSystem::from_modal(10.0, 1.0).expect("valid");
        assert!(logarithmic_decrement(&crit).is_err());
        let over = SdofSystem::from_modal(10.0, 1.7).expect("valid");
        assert!(logarithmic_decrement(&over).is_err());
    }

    #[test]
    fn decrement_from_peaks_matches_definition() {
        // Two peaks a factor e apart, one cycle  =>  delta = ln(e) = 1.
        let delta = decrement_from_peaks(std::f64::consts::E, 1.0, 1).expect("ok");
        assert!((delta - 1.0).abs() < EPS);

        // Same ratio over 2 cycles  =>  delta halves.
        let delta2 = decrement_from_peaks(std::f64::consts::E, 1.0, 2).expect("ok");
        assert!((delta2 - 0.5).abs() < EPS);
    }

    #[test]
    fn measured_decay_recovers_zeta_end_to_end() {
        // Generate exact peak amplitudes from a known system, then
        // recover its damping ratio purely from those two numbers —
        // the end-to-end log-decrement identification check.
        let zeta = 0.08;
        let sys = SdofSystem::from_modal(25.0, zeta).expect("valid");
        let true_delta = logarithmic_decrement(&sys).expect("ok");

        // Peaks n cycles apart shrink by e^(-n*delta).
        let n = 5;
        let earlier = 1.0;
        let later = earlier * (-(f64::from(n)) * true_delta).exp();

        let measured_delta = decrement_from_peaks(earlier, later, n).expect("ok");
        let recovered_zeta = damping_ratio_from_decrement(measured_delta).expect("ok");
        assert!((measured_delta - true_delta).abs() < 1e-9);
        assert!((recovered_zeta - zeta).abs() < 1e-9);
    }

    #[test]
    fn decrement_from_peaks_validates_inputs() {
        // Zero cycles rejected.
        assert!(decrement_from_peaks(2.0, 1.0, 0).is_err());
        // Non-positive amplitude rejected.
        assert!(decrement_from_peaks(0.0, 1.0, 1).is_err());
        assert!(decrement_from_peaks(2.0, -1.0, 1).is_err());
        // Growth (later > earlier) rejected.
        let err = decrement_from_peaks(1.0, 2.0, 1).expect_err("growth");
        assert_eq!(err.code(), "vibration.invalid_decay");
    }

    #[test]
    fn negative_frequency_ratio_rejected() {
        let sys = SdofSystem::from_modal(10.0, 0.1).expect("valid");
        assert!(magnification_factor(&sys, -0.5).is_err());
        assert!(magnification_factor(&sys, f64::NAN).is_err());
    }

    #[test]
    fn damping_ratio_from_decrement_validates() {
        assert!(damping_ratio_from_decrement(-0.1).is_err());
        assert!(damping_ratio_from_decrement(f64::INFINITY).is_err());
        // Zero decrement => zero damping.
        assert!((damping_ratio_from_decrement(0.0).expect("ok")).abs() < EPS);
    }
}
