//! First-order **RC** filters: a resistor and a capacitor in series.
//!
//! ## Model
//!
//! A series resistor `R` (ohms) and capacitor `C` (farads) form a
//! single-pole filter whose corner (cutoff) frequency is the textbook
//!
//! ```text
//! fc = 1 / (2 * pi * R * C)        [Hz]
//! ```
//!
//! Taking the output **across the capacitor** gives a **low-pass**
//! response; taking it **across the resistor** gives a **high-pass**
//! response. Writing the normalised frequency `x = f / fc`, the two
//! transfer functions are
//!
//! ```text
//! low-pass:   H_lp(f) = 1 / (1 + j x)
//! high-pass:  H_hp(f) = j x / (1 + j x)
//! ```
//!
//! with magnitudes
//!
//! ```text
//! |H_lp(f)| = 1 / sqrt(1 + x^2)
//! |H_hp(f)| = x / sqrt(1 + x^2)
//! ```
//!
//! At `f = fc` (`x = 1`) both magnitudes equal `1 / sqrt(2) ~ 0.7071`,
//! i.e. exactly `-3 dB` — which is why `fc` is the *cutoff* or
//! *half-power* frequency.

use crate::error::{check_component, check_frequency, Result};
use crate::response::Response;
use core::f64::consts::PI;
use serde::{Deserialize, Serialize};

/// Which terminal the output is taken from, selecting the response shape
/// of an [`RcFilter`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum RcKind {
    /// Output across the capacitor: a **low-pass** response (passes DC,
    /// attenuates high frequencies).
    LowPass,
    /// Output across the resistor: a **high-pass** response (blocks DC,
    /// passes high frequencies).
    HighPass,
}

/// A first-order RC filter defined by its resistance and capacitance.
///
/// Construct one with [`RcFilter::new`], which validates that both
/// component values are strictly-positive and finite, then query its
/// [`RcFilter::cutoff_hz`], [`RcFilter::magnitude`],
/// [`RcFilter::phase_rad`], or full [`RcFilter::response`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct RcFilter {
    resistance_ohm: f64,
    capacitance_f: f64,
    kind: RcKind,
}

impl RcFilter {
    /// Build an RC filter from a resistance `R` (ohms), a capacitance
    /// `C` (farads), and the output terminal [`RcKind`].
    ///
    /// # Errors
    ///
    /// Returns [`crate::FilterError::InvalidComponent`] if `R` or `C`
    /// is not a strictly-positive finite number.
    pub fn new(resistance_ohm: f64, capacitance_f: f64, kind: RcKind) -> Result<Self> {
        let resistance_ohm = check_component("R", resistance_ohm)?;
        let capacitance_f = check_component("C", capacitance_f)?;
        Ok(Self {
            resistance_ohm,
            capacitance_f,
            kind,
        })
    }

    /// Convenience constructor for a **low-pass** RC filter.
    ///
    /// # Errors
    ///
    /// See [`RcFilter::new`].
    pub fn low_pass(resistance_ohm: f64, capacitance_f: f64) -> Result<Self> {
        Self::new(resistance_ohm, capacitance_f, RcKind::LowPass)
    }

    /// Convenience constructor for a **high-pass** RC filter.
    ///
    /// # Errors
    ///
    /// See [`RcFilter::new`].
    pub fn high_pass(resistance_ohm: f64, capacitance_f: f64) -> Result<Self> {
        Self::new(resistance_ohm, capacitance_f, RcKind::HighPass)
    }

    /// The configured resistance `R`, in ohms.
    #[must_use]
    pub fn resistance_ohm(&self) -> f64 {
        self.resistance_ohm
    }

    /// The configured capacitance `C`, in farads.
    #[must_use]
    pub fn capacitance_f(&self) -> f64 {
        self.capacitance_f
    }

    /// Which terminal the output is taken from.
    #[must_use]
    pub fn kind(&self) -> RcKind {
        self.kind
    }

    /// The `-3 dB` cutoff frequency `fc = 1 / (2 * pi * R * C)`, in hertz.
    ///
    /// This is the same corner frequency for both the low-pass and the
    /// high-pass configuration — only the *side* that is attenuated
    /// differs.
    #[must_use]
    pub fn cutoff_hz(&self) -> f64 {
        1.0 / (2.0 * PI * self.resistance_ohm * self.capacitance_f)
    }

    /// The time constant `tau = R * C`, in seconds.
    ///
    /// Related to the cutoff by `fc = 1 / (2 * pi * tau)`.
    #[must_use]
    pub fn time_constant_s(&self) -> f64 {
        self.resistance_ohm * self.capacitance_f
    }

    /// The linear magnitude `|H(f)|` at frequency `f` (hertz).
    ///
    /// For the low-pass kind this is `1 / sqrt(1 + x^2)`; for the
    /// high-pass kind it is `x / sqrt(1 + x^2)`, where `x = f / fc`.
    ///
    /// # Errors
    ///
    /// Returns [`crate::FilterError::InvalidFrequency`] if `f` is
    /// negative or non-finite.
    pub fn magnitude(&self, freq_hz: f64) -> Result<f64> {
        let freq_hz = check_frequency(freq_hz)?;
        let x = freq_hz / self.cutoff_hz();
        let denom = (1.0 + x * x).sqrt();
        Ok(match self.kind {
            RcKind::LowPass => 1.0 / denom,
            RcKind::HighPass => x / denom,
        })
    }

    /// The phase shift `arg H(f)` at frequency `f` (hertz), in radians.
    ///
    /// For the low-pass kind the phase is `-atan(x)` (output lags,
    /// `0` at DC falling toward `-pi/2`); for the high-pass kind it is
    /// `atan(1 / x)` (output leads, `+pi/2` at DC falling toward `0`),
    /// where `x = f / fc`.
    ///
    /// At DC (`f = 0`) the high-pass phase is taken as its `+pi/2`
    /// limit.
    ///
    /// # Errors
    ///
    /// Returns [`crate::FilterError::InvalidFrequency`] if `f` is
    /// negative or non-finite.
    pub fn phase_rad(&self, freq_hz: f64) -> Result<f64> {
        let freq_hz = check_frequency(freq_hz)?;
        let x = freq_hz / self.cutoff_hz();
        Ok(match self.kind {
            RcKind::LowPass => -x.atan(),
            RcKind::HighPass => {
                if x == 0.0 {
                    PI / 2.0
                } else {
                    (1.0 / x).atan()
                }
            }
        })
    }

    /// The full [`Response`] (magnitude and phase) at frequency `f`.
    ///
    /// # Errors
    ///
    /// Returns [`crate::FilterError::InvalidFrequency`] if `f` is
    /// negative or non-finite.
    pub fn response(&self, freq_hz: f64) -> Result<Response> {
        Ok(Response::new(
            self.magnitude(freq_hz)?,
            self.phase_rad(freq_hz)?,
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::FilterError;

    /// Loose tolerance for floating-point comparisons.
    const EPS: f64 = 1e-9;

    /// `1 / sqrt(2)`, the half-power magnitude, written out for clarity.
    const INV_SQRT2: f64 = core::f64::consts::FRAC_1_SQRT_2;

    #[test]
    fn cutoff_matches_textbook_formula() {
        // R = 1 kOhm, C = 1 uF -> fc = 1/(2*pi*1e3*1e-6) ~ 159.1549 Hz.
        let f = RcFilter::low_pass(1_000.0, 1e-6).unwrap();
        let expected = 1.0 / (2.0 * PI * 1_000.0 * 1e-6);
        assert!((f.cutoff_hz() - expected).abs() < EPS);
        assert!((f.cutoff_hz() - 159.154_943_091_895_34).abs() < 1e-6);
    }

    #[test]
    fn time_constant_relates_to_cutoff() {
        let f = RcFilter::low_pass(2_200.0, 4.7e-9).unwrap();
        let tau = f.time_constant_s();
        assert!((tau - 2_200.0 * 4.7e-9).abs() < EPS);
        // fc = 1 / (2*pi*tau).
        assert!((f.cutoff_hz() - 1.0 / (2.0 * PI * tau)).abs() < EPS);
    }

    #[test]
    fn gain_is_one_over_sqrt2_at_cutoff_lowpass() {
        let f = RcFilter::low_pass(1_000.0, 1e-6).unwrap();
        let mag = f.magnitude(f.cutoff_hz()).unwrap();
        assert!((mag - INV_SQRT2).abs() < EPS);
        // ... which is exactly -3 dB.
        let db = Response::new(mag, 0.0).magnitude_db();
        assert!((db - (-3.010_299_956_639_812)).abs() < 1e-9);
    }

    #[test]
    fn gain_is_one_over_sqrt2_at_cutoff_highpass() {
        let f = RcFilter::high_pass(1_000.0, 1e-6).unwrap();
        let mag = f.magnitude(f.cutoff_hz()).unwrap();
        assert!((mag - INV_SQRT2).abs() < EPS);
    }

    #[test]
    fn lowpass_passes_dc_and_attenuates_high_freq() {
        let f = RcFilter::low_pass(1_000.0, 1e-6).unwrap();
        let fc = f.cutoff_hz();
        // DC -> unity gain.
        assert!((f.magnitude(0.0).unwrap() - 1.0).abs() < EPS);
        // A decade above cutoff: |H| = 1/sqrt(1+100) ~ 0.0995 (~ -20 dB).
        let mag_decade = f.magnitude(10.0 * fc).unwrap();
        assert!((mag_decade - 1.0 / (101.0_f64).sqrt()).abs() < EPS);
        // Strictly monotone decreasing from DC outwards.
        assert!(f.magnitude(0.5 * fc).unwrap() > f.magnitude(fc).unwrap());
        assert!(f.magnitude(fc).unwrap() > f.magnitude(2.0 * fc).unwrap());
        assert!(mag_decade < 0.1);
    }

    #[test]
    fn highpass_blocks_dc_and_passes_high_freq() {
        let f = RcFilter::high_pass(1_000.0, 1e-6).unwrap();
        let fc = f.cutoff_hz();
        // DC -> blocked (zero gain).
        assert!(f.magnitude(0.0).unwrap().abs() < EPS);
        // Far above cutoff -> approaches unity.
        let mag_high = f.magnitude(1_000.0 * fc).unwrap();
        assert!(mag_high > 0.999 && mag_high <= 1.0);
        // Monotone increasing toward the passband.
        assert!(f.magnitude(fc).unwrap() > f.magnitude(0.5 * fc).unwrap());
    }

    #[test]
    fn lowpass_and_highpass_are_power_complementary() {
        // |H_lp|^2 + |H_hp|^2 = 1 at every frequency (same R, C).
        let lp = RcFilter::low_pass(1_000.0, 1e-6).unwrap();
        let hp = RcFilter::high_pass(1_000.0, 1e-6).unwrap();
        for &freq in &[0.0, 50.0, 159.0, 1_000.0, 12_345.0] {
            let s = lp.magnitude(freq).unwrap().powi(2) + hp.magnitude(freq).unwrap().powi(2);
            assert!((s - 1.0).abs() < EPS, "failed at {freq} Hz");
        }
    }

    #[test]
    fn lowpass_phase_is_minus_45_deg_at_cutoff() {
        let f = RcFilter::low_pass(1_000.0, 1e-6).unwrap();
        let phase = f.phase_rad(f.cutoff_hz()).unwrap();
        assert!((phase - (-PI / 4.0)).abs() < EPS);
        // DC -> 0 rad; far above -> approaches -pi/2.
        assert!(f.phase_rad(0.0).unwrap().abs() < EPS);
        assert!((f.phase_rad(1e6 * f.cutoff_hz()).unwrap() + PI / 2.0).abs() < 1e-5);
    }

    #[test]
    fn highpass_phase_is_plus_45_deg_at_cutoff() {
        let f = RcFilter::high_pass(1_000.0, 1e-6).unwrap();
        let phase = f.phase_rad(f.cutoff_hz()).unwrap();
        assert!((phase - PI / 4.0).abs() < EPS);
        // DC -> +pi/2 (limit); far above -> approaches 0.
        assert!((f.phase_rad(0.0).unwrap() - PI / 2.0).abs() < EPS);
        assert!(f.phase_rad(1e6 * f.cutoff_hz()).unwrap().abs() < 1e-5);
    }

    #[test]
    fn response_bundles_magnitude_and_phase() {
        let f = RcFilter::low_pass(1_000.0, 1e-6).unwrap();
        let fc = f.cutoff_hz();
        let r = f.response(fc).unwrap();
        assert!((r.magnitude - f.magnitude(fc).unwrap()).abs() < EPS);
        assert!((r.phase_rad - f.phase_rad(fc).unwrap()).abs() < EPS);
        assert!((r.phase_deg() - (-45.0)).abs() < 1e-7);
    }

    #[test]
    fn rejects_non_physical_components() {
        assert!(matches!(
            RcFilter::low_pass(0.0, 1e-6),
            Err(FilterError::InvalidComponent { field: "R", .. })
        ));
        assert!(matches!(
            RcFilter::low_pass(1_000.0, -1.0),
            Err(FilterError::InvalidComponent { field: "C", .. })
        ));
        assert!(matches!(
            RcFilter::low_pass(f64::NAN, 1e-6),
            Err(FilterError::InvalidComponent { field: "R", .. })
        ));
        assert!(matches!(
            RcFilter::high_pass(f64::INFINITY, 1e-6),
            Err(FilterError::InvalidComponent { field: "R", .. })
        ));
    }

    #[test]
    fn rejects_non_physical_frequency() {
        let f = RcFilter::low_pass(1_000.0, 1e-6).unwrap();
        assert!(matches!(
            f.magnitude(-1.0),
            Err(FilterError::InvalidFrequency { .. })
        ));
        assert!(matches!(
            f.phase_rad(f64::NAN),
            Err(FilterError::InvalidFrequency { .. })
        ));
    }
}
