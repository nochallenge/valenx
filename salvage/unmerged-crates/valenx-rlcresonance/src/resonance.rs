//! Resonant frequency, quality factor and bandwidth.
//!
//! All values are in SI base units: inductance `L` in henries (H),
//! capacitance `C` in farads (F), resistance `R` in ohms (Ω), and every
//! frequency in hertz (Hz). The governing relations are
//!
//! ```text
//! f0 = 1 / (2 pi sqrt(L C))            resonant frequency
//! Q  = (1 / R) sqrt(L / C)             series loaded quality factor
//! BW = f0 / Q = R / (2 pi L)           -3 dB bandwidth
//! ```
//!
//! and the round-trip identity `Q * BW = f0` holds exactly (it is just
//! the definition `BW = f0 / Q` rearranged), which the test-suite checks
//! against hand-worked textbook numbers.

use crate::error::{require_non_negative, require_positive, Result, RlcError};

/// One full turn in radians, `2 pi`.
const TWO_PI: f64 = std::f64::consts::TAU;

/// Undamped resonant (natural) frequency of an LC pair, in hertz.
///
/// `f0 = 1 / (2 pi sqrt(L C))`. This is the frequency at which the
/// inductive and capacitive reactances cancel; it is identical for the
/// series and parallel ideal topologies and is independent of `R`.
///
/// # Errors
///
/// Returns [`RlcError::NonPositive`] if `inductance_h` or
/// `capacitance_f` is not finite and `> 0`.
///
/// # Example
///
/// ```
/// use valenx_rlcresonance::resonant_frequency_hz;
/// // L = 1 mH, C = 1 uF -> f0 = 1/(2 pi sqrt(1e-9)) ~= 5032.92 Hz.
/// let f0 = resonant_frequency_hz(1e-3, 1e-6).unwrap();
/// assert!((f0 - 5032.921210448704).abs() < 1e-6);
/// ```
pub fn resonant_frequency_hz(inductance_h: f64, capacitance_f: f64) -> Result<f64> {
    let l = require_positive("inductance", inductance_h)?;
    let c = require_positive("capacitance", capacitance_f)?;
    Ok(1.0 / (TWO_PI * (l * c).sqrt()))
}

/// Undamped resonant angular frequency, in radians per second.
///
/// `omega0 = 1 / sqrt(L C) = 2 pi f0`.
///
/// # Errors
///
/// Returns [`RlcError::NonPositive`] if `inductance_h` or
/// `capacitance_f` is not finite and `> 0`.
pub fn resonant_angular_frequency_rad_s(inductance_h: f64, capacitance_f: f64) -> Result<f64> {
    let l = require_positive("inductance", inductance_h)?;
    let c = require_positive("capacitance", capacitance_f)?;
    Ok(1.0 / (l * c).sqrt())
}

/// Loaded quality factor of a **series** RLC circuit (dimensionless).
///
/// `Q = (1 / R) sqrt(L / C) = omega0 L / R`. A higher `Q` means a
/// sharper resonance and a narrower bandwidth. As `R -> 0` the tank
/// becomes lossless and `Q` diverges; that singular case is reported as
/// [`RlcError::Singular`] rather than returning `+Inf`.
///
/// # Errors
///
/// Returns [`RlcError::NonPositive`] for a non-positive `inductance_h`
/// or `capacitance_f`, [`RlcError::Negative`] for a negative or
/// non-finite `resistance_ohm`, and [`RlcError::Singular`] when
/// `resistance_ohm` is exactly `0` (lossless, `Q` undefined).
///
/// # Example
///
/// ```
/// use valenx_rlcresonance::series_quality_factor;
/// // R = 10, L = 1 mH, C = 1 uF -> Q = (1/10) sqrt(1e-3/1e-6) ~= 3.1623.
/// let q = series_quality_factor(10.0, 1e-3, 1e-6).unwrap();
/// assert!((q - 3.1622776601683795).abs() < 1e-9);
/// ```
pub fn series_quality_factor(
    resistance_ohm: f64,
    inductance_h: f64,
    capacitance_f: f64,
) -> Result<f64> {
    let r = require_non_negative("resistance", resistance_ohm)?;
    let l = require_positive("inductance", inductance_h)?;
    let c = require_positive("capacitance", capacitance_f)?;
    if r == 0.0 {
        return Err(RlcError::Singular(
            "lossless series tank (R = 0): quality factor diverges",
        ));
    }
    Ok((1.0 / r) * (l / c).sqrt())
}

/// Loaded quality factor of a **parallel** RLC circuit (dimensionless).
///
/// For the standard parallel tank with a resistance across the `L`/`C`
/// pair, `Q = R sqrt(C / L) = omega0 R C`. It is the reciprocal of the
/// series expression (with `R` playing the dual role), so a *large*
/// shunt `R` gives a *high* `Q`. As `R -> infinity` (or `omega0 -> 0`)
/// the result is unbounded, but any finite, positive `R` is well
/// defined here.
///
/// # Errors
///
/// Returns [`RlcError::NonPositive`] for a non-positive `inductance_h`,
/// `capacitance_f` or `resistance_ohm`. (Unlike the series case a
/// parallel `R = 0` is a short across the tank, which has `Q = 0` but is
/// a degenerate topology, so a strictly positive `R` is required.)
///
/// # Example
///
/// ```
/// use valenx_rlcresonance::parallel_quality_factor;
/// // R = 10k, L = 1 mH, C = 1 uF -> Q = 10000 sqrt(1e-6/1e-3) ~= 316.23.
/// let q = parallel_quality_factor(10_000.0, 1e-3, 1e-6).unwrap();
/// assert!((q - 316.22776601683796).abs() < 1e-9);
/// ```
pub fn parallel_quality_factor(
    resistance_ohm: f64,
    inductance_h: f64,
    capacitance_f: f64,
) -> Result<f64> {
    let r = require_positive("resistance", resistance_ohm)?;
    let l = require_positive("inductance", inductance_h)?;
    let c = require_positive("capacitance", capacitance_f)?;
    Ok(r * (c / l).sqrt())
}

/// `-3 dB` bandwidth of a series RLC circuit, in hertz.
///
/// `BW = f0 / Q = R / (2 pi L)`. This is the width of the band over
/// which the response stays within `1/sqrt(2)` of its peak. It is
/// computed directly from `R` and `L` (no `C` dependence), so it stays
/// finite even in the lossless `R = 0` limit, where it correctly equals
/// zero.
///
/// # Errors
///
/// Returns [`RlcError::NonPositive`] for a non-positive `inductance_h`,
/// and [`RlcError::Negative`] for a negative or non-finite
/// `resistance_ohm`.
///
/// # Example
///
/// ```
/// use valenx_rlcresonance::series_bandwidth_hz;
/// // R = 10, L = 1 mH -> BW = 10/(2 pi 1e-3) ~= 1591.55 Hz.
/// let bw = series_bandwidth_hz(10.0, 1e-3).unwrap();
/// assert!((bw - 1591.5494309189533).abs() < 1e-6);
/// ```
pub fn series_bandwidth_hz(resistance_ohm: f64, inductance_h: f64) -> Result<f64> {
    let r = require_non_negative("resistance", resistance_ohm)?;
    let l = require_positive("inductance", inductance_h)?;
    Ok(r / (TWO_PI * l))
}

/// `-3 dB` bandwidth from a resonant frequency and quality factor, in hertz.
///
/// `BW = f0 / Q`. The dual of [`quality_from_f0_bandwidth`]; together
/// they encode the round-trip identity `Q * BW = f0` that the tests pin
/// to textbook numbers.
///
/// # Errors
///
/// Returns [`RlcError::NonPositive`] if either `f0_hz` or
/// `quality_factor` is not finite and `> 0`.
///
/// # Example
///
/// ```
/// use valenx_rlcresonance::bandwidth_from_f0_quality;
/// let bw = bandwidth_from_f0_quality(5032.921210448704, 3.1622776601683795).unwrap();
/// assert!((bw - 1591.5494309189533).abs() < 1e-6);
/// ```
pub fn bandwidth_from_f0_quality(f0_hz: f64, quality_factor: f64) -> Result<f64> {
    let f0 = require_positive("f0", f0_hz)?;
    let q = require_positive("quality_factor", quality_factor)?;
    Ok(f0 / q)
}

/// Quality factor from a resonant frequency and `-3 dB` bandwidth.
///
/// `Q = f0 / BW`. The dual of [`bandwidth_from_f0_quality`]; rearranges
/// the same `Q * BW = f0` identity.
///
/// # Errors
///
/// Returns [`RlcError::NonPositive`] if either `f0_hz` or
/// `bandwidth_hz` is not finite and `> 0`.
///
/// # Example
///
/// ```
/// use valenx_rlcresonance::quality_from_f0_bandwidth;
/// let q = quality_from_f0_bandwidth(5032.921210448704, 1591.5494309189533).unwrap();
/// assert!((q - 3.1622776601683795).abs() < 1e-9);
/// ```
pub fn quality_from_f0_bandwidth(f0_hz: f64, bandwidth_hz: f64) -> Result<f64> {
    let f0 = require_positive("f0", f0_hz)?;
    let bw = require_positive("bandwidth", bandwidth_hz)?;
    Ok(f0 / bw)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::RlcError;

    /// Absolute tolerance for tight closed-form comparisons.
    const EPS: f64 = 1e-9;
    /// Looser tolerance for derived quantities carrying a `2 pi` factor.
    const EPS_HZ: f64 = 1e-6;

    #[test]
    fn f0_matches_textbook_1mh_1uf() {
        // L = 1 mH, C = 1 uF. sqrt(LC) = sqrt(1e-9) = 3.1622776601e-5 s.
        // f0 = 1 / (2 pi * 3.1622776601e-5) = 5032.9212104 Hz.
        let f0 = resonant_frequency_hz(1e-3, 1e-6).unwrap();
        assert!((f0 - 5032.921210448704).abs() < EPS_HZ, "f0 = {f0}");
    }

    #[test]
    fn f0_unit_lc_is_one_over_two_pi() {
        // L = C = 1 -> sqrt(LC) = 1 -> f0 = 1/(2 pi).
        let f0 = resonant_frequency_hz(1.0, 1.0).unwrap();
        assert!((f0 - 1.0 / std::f64::consts::TAU).abs() < EPS, "f0 = {f0}");
    }

    #[test]
    fn omega0_is_two_pi_times_f0() {
        let l = 2.7e-3;
        let c = 4.7e-9;
        let f0 = resonant_frequency_hz(l, c).unwrap();
        let w0 = resonant_angular_frequency_rad_s(l, c).unwrap();
        assert!((w0 - std::f64::consts::TAU * f0).abs() < EPS, "w0 = {w0}");
        // And omega0 = 1/sqrt(LC) directly.
        assert!((w0 - 1.0 / (l * c).sqrt()).abs() < EPS, "w0 = {w0}");
    }

    #[test]
    fn series_q_matches_textbook() {
        // R = 10, L = 1 mH, C = 1 uF.
        // sqrt(L/C) = sqrt(1e-3/1e-6) = sqrt(1000) = 31.6227766.
        // Q = (1/10) * 31.6227766 = 3.16227766.
        let q = series_quality_factor(10.0, 1e-3, 1e-6).unwrap();
        assert!((q - 3.1622776601683795).abs() < EPS, "q = {q}");
    }

    #[test]
    fn series_q_equals_omega0_l_over_r() {
        // Q = omega0 * L / R is an equivalent closed form; check it agrees.
        let (r, l, c) = (4.7, 3.3e-3, 2.2e-7);
        let q = series_quality_factor(r, l, c).unwrap();
        let w0 = resonant_angular_frequency_rad_s(l, c).unwrap();
        assert!((q - w0 * l / r).abs() < EPS, "q = {q}");
    }

    #[test]
    fn parallel_q_is_reciprocal_form_of_series() {
        // For the same R, L, C the parallel Q = R sqrt(C/L) is the
        // reciprocal (in R) of the series Q = (1/R) sqrt(L/C):
        // series_Q * parallel_Q = (1/R) sqrt(L/C) * R sqrt(C/L) = 1.
        let (r, l, c) = (12.0, 1e-3, 1e-6);
        let qs = series_quality_factor(r, l, c).unwrap();
        let qp = parallel_quality_factor(r, l, c).unwrap();
        assert!((qs * qp - 1.0).abs() < EPS, "qs*qp = {}", qs * qp);
    }

    #[test]
    fn parallel_q_matches_textbook() {
        // R = 10k, L = 1 mH, C = 1 uF.
        // sqrt(C/L) = sqrt(1e-6/1e-3) = sqrt(1e-3) = 0.0316227766.
        // Q = 10000 * 0.0316227766 = 316.227766.
        let q = parallel_quality_factor(10_000.0, 1e-3, 1e-6).unwrap();
        assert!((q - 316.22776601683796).abs() < EPS, "q = {q}");
    }

    #[test]
    fn series_bandwidth_matches_textbook() {
        // R = 10, L = 1 mH -> BW = 10/(2 pi * 1e-3) = 1591.549431 Hz.
        let bw = series_bandwidth_hz(10.0, 1e-3).unwrap();
        assert!((bw - 1591.5494309189533).abs() < EPS_HZ, "bw = {bw}");
    }

    #[test]
    fn lossless_bandwidth_is_zero() {
        // R = 0 is the ideal lossless limit; BW = R/(2 pi L) = 0 exactly.
        let bw = series_bandwidth_hz(0.0, 1e-3).unwrap();
        assert!(bw.abs() < EPS, "bw = {bw}");
    }

    #[test]
    fn ground_truth_q_times_bw_equals_f0() {
        // THE headline ground-truth identity: Q * BW = f0, exactly,
        // for the series circuit, computed from independent formulas.
        let (r, l, c) = (10.0, 1e-3, 1e-6);
        let f0 = resonant_frequency_hz(l, c).unwrap();
        let q = series_quality_factor(r, l, c).unwrap();
        let bw = series_bandwidth_hz(r, l).unwrap();
        assert!((q * bw - f0).abs() < EPS_HZ, "q*bw = {} f0 = {f0}", q * bw);
    }

    #[test]
    fn ground_truth_identity_holds_for_a_second_circuit() {
        // Re-check Q*BW = f0 with unrelated component values so the
        // identity is not an artefact of the round decade values above.
        let (r, l, c) = (47.0, 6.8e-3, 3.3e-9);
        let f0 = resonant_frequency_hz(l, c).unwrap();
        let q = series_quality_factor(r, l, c).unwrap();
        let bw = series_bandwidth_hz(r, l).unwrap();
        assert!((q * bw - f0).abs() < EPS_HZ, "q*bw = {} f0 = {f0}", q * bw);
    }

    #[test]
    fn bandwidth_and_quality_helpers_round_trip() {
        // BW = f0/Q and Q = f0/BW are inverses: feed one into the other.
        let f0 = 5032.921210448704;
        let q = 3.1622776601683795;
        let bw = bandwidth_from_f0_quality(f0, q).unwrap();
        let q_back = quality_from_f0_bandwidth(f0, bw).unwrap();
        assert!((q_back - q).abs() < EPS, "q_back = {q_back}");
        // And the explicit identity Q * BW = f0.
        assert!((q * bw - f0).abs() < EPS_HZ, "q*bw = {}", q * bw);
    }

    #[test]
    fn bandwidth_helper_agrees_with_series_formula() {
        // bandwidth_from_f0_quality(f0, Q) must equal series_bandwidth_hz.
        let (r, l, c) = (22.0, 1.5e-3, 1e-7);
        let f0 = resonant_frequency_hz(l, c).unwrap();
        let q = series_quality_factor(r, l, c).unwrap();
        let bw_a = bandwidth_from_f0_quality(f0, q).unwrap();
        let bw_b = series_bandwidth_hz(r, l).unwrap();
        assert!((bw_a - bw_b).abs() < EPS_HZ, "bw_a = {bw_a} bw_b = {bw_b}");
    }

    #[test]
    fn rejects_non_positive_inductance_and_capacitance() {
        assert_eq!(
            resonant_frequency_hz(0.0, 1e-6).unwrap_err().code(),
            "rlc.non_positive"
        );
        assert_eq!(
            resonant_frequency_hz(1e-3, -1.0).unwrap_err().code(),
            "rlc.non_positive"
        );
        assert!(matches!(
            resonant_frequency_hz(f64::NAN, 1e-6),
            Err(RlcError::NonPositive {
                name: "inductance",
                ..
            })
        ));
        assert!(matches!(
            resonant_frequency_hz(1e-3, f64::INFINITY),
            Err(RlcError::NonPositive {
                name: "capacitance",
                ..
            })
        ));
    }

    #[test]
    fn rejects_negative_resistance_but_accepts_zero_in_bandwidth() {
        assert_eq!(
            series_bandwidth_hz(-1.0, 1e-3).unwrap_err().code(),
            "rlc.negative"
        );
        assert!(matches!(
            series_bandwidth_hz(f64::NAN, 1e-3),
            Err(RlcError::Negative {
                name: "resistance",
                ..
            })
        ));
        // Zero resistance is the valid lossless limit for bandwidth.
        assert!(series_bandwidth_hz(0.0, 1e-3).is_ok());
    }

    #[test]
    fn series_q_zero_resistance_is_singular() {
        let err = series_quality_factor(0.0, 1e-3, 1e-6).unwrap_err();
        assert_eq!(err.code(), "rlc.singular");
        assert!(matches!(err, RlcError::Singular(_)));
    }

    #[test]
    fn parallel_q_requires_positive_resistance() {
        assert_eq!(
            parallel_quality_factor(0.0, 1e-3, 1e-6).unwrap_err().code(),
            "rlc.non_positive"
        );
        assert!(matches!(
            parallel_quality_factor(-5.0, 1e-3, 1e-6),
            Err(RlcError::NonPositive {
                name: "resistance",
                ..
            })
        ));
    }

    #[test]
    fn helpers_reject_non_positive_inputs() {
        assert_eq!(
            bandwidth_from_f0_quality(0.0, 3.0).unwrap_err().code(),
            "rlc.non_positive"
        );
        assert_eq!(
            quality_from_f0_bandwidth(1000.0, 0.0).unwrap_err().code(),
            "rlc.non_positive"
        );
    }
}
