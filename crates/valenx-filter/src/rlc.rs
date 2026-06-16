//! Series **RLC** resonant circuit: resistor, inductor, and capacitor.
//!
//! ## Model
//!
//! A resistor `R` (ohms), inductor `L` (henries), and capacitor `C`
//! (farads) in series resonate at the textbook
//!
//! ```text
//! f0 = 1 / (2 * pi * sqrt(L * C))        [Hz]
//! ```
//!
//! the frequency at which the inductive and capacitive reactances
//! cancel. The sharpness of the resonance is the **quality factor**
//!
//! ```text
//! Q = (1 / R) * sqrt(L / C)              [dimensionless]
//! ```
//!
//! and the `-3 dB` **bandwidth** between the two half-power frequencies
//! either side of `f0` is
//!
//! ```text
//! BW = f0 / Q                            [Hz]
//! ```
//!
//! A larger `Q` therefore means a **narrower** bandwidth — a more
//! sharply-peaked, more selective resonance. (`Q = (1/R)*sqrt(L/C)` is
//! the standard *series*-RLC expression, where a smaller `R` raises
//! `Q`.)

use crate::error::{check_component, Result};
use core::f64::consts::PI;
use serde::{Deserialize, Serialize};

/// A series RLC circuit defined by its three component values.
///
/// Construct one with [`RlcCircuit::new`], which validates that every
/// component value is strictly-positive and finite, then read off its
/// [`RlcCircuit::resonant_hz`], [`RlcCircuit::quality_factor`], and
/// [`RlcCircuit::bandwidth_hz`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct RlcCircuit {
    resistance_ohm: f64,
    inductance_h: f64,
    capacitance_f: f64,
}

impl RlcCircuit {
    /// Build a series RLC circuit from `R` (ohms), `L` (henries), and
    /// `C` (farads).
    ///
    /// # Errors
    ///
    /// Returns [`crate::FilterError::InvalidComponent`] if any of `R`,
    /// `L`, or `C` is not a strictly-positive finite number.
    pub fn new(resistance_ohm: f64, inductance_h: f64, capacitance_f: f64) -> Result<Self> {
        let resistance_ohm = check_component("R", resistance_ohm)?;
        let inductance_h = check_component("L", inductance_h)?;
        let capacitance_f = check_component("C", capacitance_f)?;
        Ok(Self {
            resistance_ohm,
            inductance_h,
            capacitance_f,
        })
    }

    /// The configured resistance `R`, in ohms.
    #[must_use]
    pub fn resistance_ohm(&self) -> f64 {
        self.resistance_ohm
    }

    /// The configured inductance `L`, in henries.
    #[must_use]
    pub fn inductance_h(&self) -> f64 {
        self.inductance_h
    }

    /// The configured capacitance `C`, in farads.
    #[must_use]
    pub fn capacitance_f(&self) -> f64 {
        self.capacitance_f
    }

    /// The resonant frequency `f0 = 1 / (2 * pi * sqrt(L * C))`, in hertz.
    #[must_use]
    pub fn resonant_hz(&self) -> f64 {
        1.0 / (2.0 * PI * (self.inductance_h * self.capacitance_f).sqrt())
    }

    /// The resonant **angular** frequency
    /// `omega0 = 1 / sqrt(L * C) = 2 * pi * f0`, in radians per second.
    #[must_use]
    pub fn resonant_omega(&self) -> f64 {
        1.0 / (self.inductance_h * self.capacitance_f).sqrt()
    }

    /// The quality factor `Q = (1 / R) * sqrt(L / C)` (dimensionless).
    ///
    /// This is the standard *series*-RLC expression: a smaller series
    /// resistance gives a higher `Q` and hence a sharper resonance.
    #[must_use]
    pub fn quality_factor(&self) -> f64 {
        (1.0 / self.resistance_ohm) * (self.inductance_h / self.capacitance_f).sqrt()
    }

    /// The `-3 dB` bandwidth `BW = f0 / Q`, in hertz.
    ///
    /// For a series RLC this also equals `R / (2 * pi * L)`, which the
    /// implementation does **not** assume — it is computed directly as
    /// `f0 / Q` so the `BW = f0 / Q` relationship is exact by
    /// construction.
    #[must_use]
    pub fn bandwidth_hz(&self) -> f64 {
        self.resonant_hz() / self.quality_factor()
    }

    /// The lower `-3 dB` half-power frequency `f0 - BW / 2`, in hertz.
    ///
    /// This is the standard symmetric-in-`f` narrow-band approximation
    /// (accurate for `Q >> 1`); the geometric mean of the two returned
    /// edge frequencies is therefore close to, but for low `Q` not
    /// exactly, `f0`.
    #[must_use]
    pub fn lower_cutoff_hz(&self) -> f64 {
        self.resonant_hz() - self.bandwidth_hz() / 2.0
    }

    /// The upper `-3 dB` half-power frequency `f0 + BW / 2`, in hertz.
    ///
    /// See [`RlcCircuit::lower_cutoff_hz`] for the approximation used.
    #[must_use]
    pub fn upper_cutoff_hz(&self) -> f64 {
        self.resonant_hz() + self.bandwidth_hz() / 2.0
    }

    /// The **exact** lower `-3 dB` half-power frequency, in hertz.
    ///
    /// Unlike [`RlcCircuit::lower_cutoff_hz`] (the narrow-band
    /// `f0 - BW/2` approximation), this solves the half-power condition
    /// `|omega L - 1/(omega C)| = R` exactly, giving
    ///
    /// ```text
    /// f_low = f0 * (sqrt(1 + 1/(2Q)^2) - 1/(2Q))
    /// ```
    ///
    /// The exact pair is symmetric in the **geometric** (not arithmetic)
    /// sense — their product is exactly `f0^2`, so `sqrt(f_low*f_high)`
    /// is the resonant frequency — yet their difference is still exactly
    /// the bandwidth `f0/Q`. As `Q -> infinity` this converges to the
    /// narrow-band value.
    #[must_use]
    pub fn lower_cutoff_exact_hz(&self) -> f64 {
        let inv_2q = 1.0 / (2.0 * self.quality_factor());
        self.resonant_hz() * ((1.0 + inv_2q * inv_2q).sqrt() - inv_2q)
    }

    /// The **exact** upper `-3 dB` half-power frequency, in hertz.
    ///
    /// The exact counterpart of [`RlcCircuit::upper_cutoff_hz`]:
    ///
    /// ```text
    /// f_high = f0 * (sqrt(1 + 1/(2Q)^2) + 1/(2Q))
    /// ```
    ///
    /// See [`RlcCircuit::lower_cutoff_exact_hz`] for the geometric-symmetry
    /// and bandwidth identities the exact pair satisfies.
    #[must_use]
    pub fn upper_cutoff_exact_hz(&self) -> f64 {
        let inv_2q = 1.0 / (2.0 * self.quality_factor());
        self.resonant_hz() * ((1.0 + inv_2q * inv_2q).sqrt() + inv_2q)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::FilterError;

    /// Loose tolerance for floating-point comparisons.
    const EPS: f64 = 1e-9;

    #[test]
    fn resonant_frequency_matches_textbook_formula() {
        // L = 10 mH, C = 1 uF -> f0 = 1/(2*pi*sqrt(1e-2*1e-6)).
        let c = RlcCircuit::new(10.0, 10e-3, 1e-6).unwrap();
        let expected = 1.0 / (2.0 * PI * (10e-3 * 1e-6_f64).sqrt());
        assert!((c.resonant_hz() - expected).abs() < EPS);
        // ~ 1591.549 Hz.
        assert!((c.resonant_hz() - 1_591.549_430_918_953_4).abs() < 1e-6);
    }

    #[test]
    fn angular_frequency_is_two_pi_times_resonant() {
        let c = RlcCircuit::new(10.0, 10e-3, 1e-6).unwrap();
        assert!((c.resonant_omega() - 2.0 * PI * c.resonant_hz()).abs() < EPS);
        assert!((c.resonant_omega() - 1.0 / (10e-3 * 1e-6_f64).sqrt()).abs() < EPS);
    }

    #[test]
    fn quality_factor_matches_textbook_formula() {
        // Q = (1/R)*sqrt(L/C) = (1/10)*sqrt(1e-2/1e-6) = 0.1*100 = 10.
        let c = RlcCircuit::new(10.0, 10e-3, 1e-6).unwrap();
        assert!((c.quality_factor() - 10.0).abs() < EPS);
    }

    #[test]
    fn bandwidth_equals_f0_over_q() {
        let c = RlcCircuit::new(10.0, 10e-3, 1e-6).unwrap();
        let bw = c.bandwidth_hz();
        assert!((bw - c.resonant_hz() / c.quality_factor()).abs() < EPS);
        // For series RLC, BW should also equal R/(2*pi*L) — independent check.
        let bw_alt = 10.0 / (2.0 * PI * 10e-3);
        assert!((bw - bw_alt).abs() < 1e-7);
    }

    #[test]
    fn higher_q_gives_narrower_bandwidth() {
        // Same L and C (so same f0); lower R raises Q.
        let low_q = RlcCircuit::new(100.0, 10e-3, 1e-6).unwrap(); // Q = 1
        let high_q = RlcCircuit::new(10.0, 10e-3, 1e-6).unwrap(); // Q = 10

        // Same resonant frequency.
        assert!((low_q.resonant_hz() - high_q.resonant_hz()).abs() < EPS);
        // Higher Q.
        assert!(high_q.quality_factor() > low_q.quality_factor());
        // ... and therefore strictly narrower bandwidth.
        assert!(high_q.bandwidth_hz() < low_q.bandwidth_hz());
        // Concretely BW scales as 1/Q: Q 1 -> 10 means BW shrinks 10x.
        assert!((low_q.bandwidth_hz() / high_q.bandwidth_hz() - 10.0).abs() < 1e-6);
    }

    #[test]
    fn cutoffs_straddle_resonance_by_half_bandwidth() {
        let c = RlcCircuit::new(10.0, 10e-3, 1e-6).unwrap();
        let f0 = c.resonant_hz();
        let bw = c.bandwidth_hz();
        assert!((c.upper_cutoff_hz() - c.lower_cutoff_hz() - bw).abs() < 1e-7);
        // Midpoint of the two edges is f0.
        let mid = 0.5 * (c.lower_cutoff_hz() + c.upper_cutoff_hz());
        assert!((mid - f0).abs() < 1e-7);
    }

    #[test]
    fn exact_cutoffs_match_closed_form() {
        let c = RlcCircuit::new(10.0, 10e-3, 1e-6).unwrap(); // Q = 10
        let f0 = c.resonant_hz();
        let inv_2q = 1.0 / (2.0 * c.quality_factor());
        let root = (1.0 + inv_2q * inv_2q).sqrt();
        assert!((c.lower_cutoff_exact_hz() - f0 * (root - inv_2q)).abs() < EPS);
        assert!((c.upper_cutoff_exact_hz() - f0 * (root + inv_2q)).abs() < EPS);
    }

    #[test]
    fn exact_cutoffs_bandwidth_is_exact_at_any_q() {
        // The exact edges' difference equals f0/Q exactly, at every Q —
        // unlike where the narrow-band approximation's center drifts.
        for &r in &[1.0, 10.0, 100.0, 500.0] {
            let c = RlcCircuit::new(r, 10e-3, 1e-6).unwrap();
            let bw = c.upper_cutoff_exact_hz() - c.lower_cutoff_exact_hz();
            assert!((bw - c.bandwidth_hz()).abs() < 1e-6, "R={r}");
        }
    }

    #[test]
    fn exact_cutoffs_geometric_mean_is_f0() {
        // f0 is the GEOMETRIC mean of the exact half-power edges (the band
        // is symmetric in log frequency): f_low*f_high = f0^2 exactly,
        // unlike the arithmetic-midpoint approximation.
        for &r in &[1.0, 50.0, 100.0] {
            let c = RlcCircuit::new(r, 10e-3, 1e-6).unwrap();
            let f0 = c.resonant_hz();
            let prod = c.lower_cutoff_exact_hz() * c.upper_cutoff_exact_hz();
            assert!((prod - f0 * f0).abs() / (f0 * f0) < 1e-12, "R={r}");
            assert!((prod.sqrt() - f0).abs() < 1e-6, "R={r}");
        }
    }

    #[test]
    fn exact_converges_to_approx_at_high_q() {
        // As Q grows the exact edges approach the narrow-band f0 +/- BW/2.
        let c = RlcCircuit::new(0.1, 10e-3, 1e-6).unwrap(); // Q = 1000
        let f0 = c.resonant_hz();
        assert!((c.lower_cutoff_exact_hz() - c.lower_cutoff_hz()).abs() / f0 < 1e-6);
        assert!((c.upper_cutoff_exact_hz() - c.upper_cutoff_hz()).abs() / f0 < 1e-6);
    }

    #[test]
    fn exact_edges_are_golden_ratio_at_q_one() {
        // At Q = 1 the exact edges are f0/phi and f0*phi (golden-ratio
        // conjugates), well away from the approximation's 0.5 f0 / 1.5 f0.
        let c = RlcCircuit::new(100.0, 10e-3, 1e-6).unwrap(); // Q = 1
        let f0 = c.resonant_hz();
        assert!((c.lower_cutoff_exact_hz() / f0 - 0.618_033_988_75).abs() < 1e-6);
        assert!((c.upper_cutoff_exact_hz() / f0 - 1.618_033_988_75).abs() < 1e-6);
        // The exact lower edge is strictly above the narrow-band value.
        assert!(c.lower_cutoff_exact_hz() > c.lower_cutoff_hz());
    }

    #[test]
    fn rejects_non_physical_components() {
        assert!(matches!(
            RlcCircuit::new(0.0, 10e-3, 1e-6),
            Err(FilterError::InvalidComponent { field: "R", .. })
        ));
        assert!(matches!(
            RlcCircuit::new(10.0, -1.0, 1e-6),
            Err(FilterError::InvalidComponent { field: "L", .. })
        ));
        assert!(matches!(
            RlcCircuit::new(10.0, 10e-3, f64::NAN),
            Err(FilterError::InvalidComponent { field: "C", .. })
        ));
        assert!(matches!(
            RlcCircuit::new(10.0, f64::INFINITY, 1e-6),
            Err(FilterError::InvalidComponent { field: "L", .. })
        ));
    }
}
