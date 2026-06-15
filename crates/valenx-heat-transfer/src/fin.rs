//! Extended-surface (straight fin) heat rate and efficiency.
//!
//! ## Model
//!
//! For a straight fin of uniform cross-section — cross-sectional area
//! `A_c`, perimeter `P`, conductivity `k`, length `L`, base
//! over-temperature `θ_b = T_base − T_∞` and uniform film coefficient
//! `h` — the 1D fin equation `d²θ/dx² − m²θ = 0` (with
//! `m = √(h·P / (k·A_c))`) has the classic closed-form solution. With
//! an **adiabatic (insulated) tip** the base heat rate is
//!
//! ```text
//! q_fin = √(h · P · k · A_c) · θ_b · tanh(m · L)
//! ```
//!
//! and the **fin efficiency** — actual heat dissipated divided by the
//! ideal heat that would be dissipated if the entire fin sat at the
//! base temperature, `q_max = h · A_fin · θ_b` with
//! `A_fin = P · L` — is
//!
//! ```text
//! η = tanh(m · L) / (m · L)        ∈ (0, 1].
//! ```
//!
//! Because `tanh(x)/x ∈ (0, 1]` for `x > 0` and `→ 1` as `x → 0`, the
//! efficiency is bounded in `(0, 1]`: a vanishingly short or perfectly
//! conductive fin approaches `η = 1`, while a long or poorly conducting
//! fin has `η → 0`.
//!
//! A **convective-tip** correction is also provided using the standard
//! "corrected length" approximation `L_c = L + A_c / P`, evaluating the
//! adiabatic solution at `L_c`.
//!
//! References: Incropera §3.6, Cengel §3-6. Results are the idealised
//! constant-property, uniform-`h`, 1D solutions.

use serde::{Deserialize, Serialize};

use crate::error::{require_finite, require_positive, Result};

/// Tip boundary condition for a straight fin.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum TipCondition {
    /// Insulated (adiabatic) tip — the canonical
    /// `q = √(hPkA)·θ_b·tanh(mL)` closed form.
    Adiabatic,
    /// Convecting tip, approximated by the corrected length
    /// `L_c = L + A_c / P`.
    ConvectiveCorrected,
}

/// A straight fin of uniform cross-section.
///
/// SI units throughout: `length_m` in m, `perimeter_m` in m,
/// `area_c_m2` (cross-section) in m², `k_w_per_mk` in W/(m·K),
/// `h_w_per_m2k` in W/(m²·K).
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Fin {
    /// Fin length `L` from base to tip (m).
    pub length_m: f64,
    /// Cross-sectional perimeter `P` (m).
    pub perimeter_m: f64,
    /// Cross-sectional area `A_c` (m²).
    pub area_c_m2: f64,
    /// Thermal conductivity `k` (W/(m·K)).
    pub k_w_per_mk: f64,
    /// Convection coefficient `h` (W/(m²·K)).
    pub h_w_per_m2k: f64,
}

impl Fin {
    /// Build a validated fin.
    ///
    /// # Errors
    ///
    /// Returns [`HeatTransferError::BadParameter`](crate::HeatTransferError::BadParameter)
    /// if any of the five geometric/material parameters is not finite
    /// and strictly positive.
    pub fn new(
        length_m: f64,
        perimeter_m: f64,
        area_c_m2: f64,
        k_w_per_mk: f64,
        h_w_per_m2k: f64,
    ) -> Result<Self> {
        Ok(Self {
            length_m: require_positive("length_m", length_m)?,
            perimeter_m: require_positive("perimeter_m", perimeter_m)?,
            area_c_m2: require_positive("area_c_m2", area_c_m2)?,
            k_w_per_mk: require_positive("k_w_per_mk", k_w_per_mk)?,
            h_w_per_m2k: require_positive("h_w_per_m2k", h_w_per_m2k)?,
        })
    }

    /// Convenience constructor for a thin **rectangular** fin of width
    /// `w` and thickness `t`: `P = 2(w + t)`, `A_c = w · t`.
    ///
    /// # Errors
    ///
    /// Returns an error if any argument is not finite and strictly
    /// positive.
    pub fn rectangular(
        length_m: f64,
        width_m: f64,
        thickness_m: f64,
        k_w_per_mk: f64,
        h_w_per_m2k: f64,
    ) -> Result<Self> {
        let w = require_positive("width_m", width_m)?;
        let t = require_positive("thickness_m", thickness_m)?;
        Fin::new(length_m, 2.0 * (w + t), w * t, k_w_per_mk, h_w_per_m2k)
    }

    /// Fin parameter `m = √(h·P / (k·A_c))` (1/m).
    pub fn m(&self) -> f64 {
        (self.h_w_per_m2k * self.perimeter_m / (self.k_w_per_mk * self.area_c_m2)).sqrt()
    }

    /// Effective length used for the heat-rate evaluation given the tip
    /// condition: `L` (adiabatic) or `L_c = L + A_c / P` (convective).
    pub fn effective_length(&self, tip: TipCondition) -> f64 {
        match tip {
            TipCondition::Adiabatic => self.length_m,
            TipCondition::ConvectiveCorrected => self.length_m + self.area_c_m2 / self.perimeter_m,
        }
    }

    /// Base heat rate `q = √(h·P·k·A_c) · θ_b · tanh(m·L_eff)` (W).
    ///
    /// `theta_b` is the base over-temperature `T_base − T_∞`; a positive
    /// value yields a positive (outgoing) heat rate.
    ///
    /// # Errors
    ///
    /// Returns an error if `theta_b` is non-finite.
    pub fn heat_rate(&self, theta_b: f64, tip: TipCondition) -> Result<f64> {
        let theta_b = require_finite("theta_b", theta_b)?;
        let m = self.m();
        let big_m = (self.h_w_per_m2k * self.perimeter_m * self.k_w_per_mk * self.area_c_m2).sqrt();
        Ok(big_m * theta_b * (m * self.effective_length(tip)).tanh())
    }

    /// Fin efficiency `η = tanh(m·L_eff) / (m·L_eff) ∈ (0, 1]`.
    ///
    /// The ratio of the actual heat dissipated to the ideal heat that
    /// an isothermal (base-temperature) fin would dissipate. Always in
    /// `(0, 1]`, approaching `1` as `m·L_eff → 0`.
    pub fn efficiency(&self, tip: TipCondition) -> f64 {
        let ml = self.m() * self.effective_length(tip);
        if ml == 0.0 {
            // Limit tanh(x)/x → 1; unreachable for validated (positive)
            // fins but kept for total-function safety.
            1.0
        } else {
            ml.tanh() / ml
        }
    }

    /// Ideal (maximum) heat rate `q_max = h · A_fin · θ_b` (W) for an
    /// isothermal fin, where the fin surface area
    /// `A_fin = P · L_eff`.
    ///
    /// # Errors
    ///
    /// Returns an error if `theta_b` is non-finite.
    pub fn max_heat_rate(&self, theta_b: f64, tip: TipCondition) -> Result<f64> {
        let theta_b = require_finite("theta_b", theta_b)?;
        let a_fin = self.perimeter_m * self.effective_length(tip);
        Ok(self.h_w_per_m2k * a_fin * theta_b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    /// Reusable textbook-ish fin: aluminium pin, moderate convection.
    fn sample() -> Fin {
        // L = 0.05 m, P = 0.02 m, A_c = 2.5e-5 m^2, k = 180, h = 25.
        Fin::new(0.05, 0.02, 2.5e-5, 180.0, 25.0).unwrap()
    }

    #[test]
    fn m_matches_definition() {
        let f = sample();
        // m = sqrt(25 * 0.02 / (180 * 2.5e-5)) = sqrt(0.5 / 0.0045)
        //   = sqrt(111.111...) ≈ 10.5409 1/m.
        let expected = (25.0 * 0.02 / (180.0 * 2.5e-5_f64)).sqrt();
        assert!((f.m() - expected).abs() < 1e-9);
        assert!((f.m() - 10.540925533894598).abs() < 1e-9);
    }

    #[test]
    fn heat_rate_matches_canonical_formula() {
        let f = sample();
        let theta_b = 80.0;
        // q = sqrt(h*P*k*A_c) * θ_b * tanh(m*L).
        let big_m = (25.0 * 0.02 * 180.0 * 2.5e-5_f64).sqrt();
        let m = (25.0 * 0.02 / (180.0 * 2.5e-5_f64)).sqrt();
        let expected = big_m * theta_b * (m * 0.05).tanh();
        let got = f.heat_rate(theta_b, TipCondition::Adiabatic).unwrap();
        assert!((got - expected).abs() < 1e-9);
        assert!(got > 0.0);
    }

    #[test]
    fn efficiency_is_in_unit_interval() {
        let f = sample();
        let eta = f.efficiency(TipCondition::Adiabatic);
        assert!(eta > 0.0 && eta <= 1.0, "eta out of (0,1]: {eta}");
        // η = tanh(mL)/(mL).
        let ml = f.m() * 0.05;
        assert!((eta - ml.tanh() / ml).abs() < 1e-12);
    }

    #[test]
    fn efficiency_approaches_one_for_short_fin() {
        // A very short / very conductive fin: mL → 0 so η → 1.
        let f = Fin::new(1e-4, 0.02, 2.5e-5, 1.0e5, 5.0).unwrap();
        let eta = f.efficiency(TipCondition::Adiabatic);
        assert!(eta > 0.999 && eta <= 1.0, "eta = {eta}");
    }

    #[test]
    fn efficiency_drops_for_long_fin() {
        // A long, low-k fin dissipates far below the isothermal ideal.
        let short = Fin::new(0.02, 0.02, 2.5e-5, 180.0, 25.0).unwrap();
        let long = Fin::new(0.20, 0.02, 2.5e-5, 180.0, 25.0).unwrap();
        let eta_short = short.efficiency(TipCondition::Adiabatic);
        let eta_long = long.efficiency(TipCondition::Adiabatic);
        assert!(eta_long < eta_short);
        assert!(eta_long < 0.5);
    }

    #[test]
    fn efficiency_definition_matches_q_over_qmax() {
        // η must equal q_fin / q_max independently of the closed form.
        let f = sample();
        let theta_b = 80.0;
        let q = f.heat_rate(theta_b, TipCondition::Adiabatic).unwrap();
        let q_max = f.max_heat_rate(theta_b, TipCondition::Adiabatic).unwrap();
        let eta_ratio = q / q_max;
        let eta_formula = f.efficiency(TipCondition::Adiabatic);
        assert!((eta_ratio - eta_formula).abs() < 1e-9);
    }

    #[test]
    fn convective_tip_raises_effective_length() {
        let f = sample();
        let l_adi = f.effective_length(TipCondition::Adiabatic);
        let l_corr = f.effective_length(TipCondition::ConvectiveCorrected);
        assert!(l_corr > l_adi);
        // L_c = L + A_c/P = 0.05 + 2.5e-5/0.02 = 0.05125.
        assert!((l_corr - 0.05125).abs() < EPS);
    }

    #[test]
    fn convective_tip_dissipates_more_than_adiabatic() {
        // A convecting tip adds surface, so it sheds more heat.
        let f = sample();
        let q_adi = f.heat_rate(80.0, TipCondition::Adiabatic).unwrap();
        let q_conv = f
            .heat_rate(80.0, TipCondition::ConvectiveCorrected)
            .unwrap();
        assert!(q_conv > q_adi);
    }

    #[test]
    fn rectangular_helper_sets_geometry() {
        // w = 0.04, t = 0.002 -> P = 2*(0.042) = 0.084, A_c = 8e-5.
        let f = Fin::rectangular(0.05, 0.04, 0.002, 180.0, 25.0).unwrap();
        assert!((f.perimeter_m - 0.084).abs() < EPS);
        assert!((f.area_c_m2 - 8e-5).abs() < 1e-12);
    }

    #[test]
    fn rejects_non_positive_inputs() {
        assert!(Fin::new(0.0, 0.02, 2.5e-5, 180.0, 25.0).is_err());
        assert!(Fin::new(0.05, 0.02, 2.5e-5, 180.0, -1.0).is_err());
        assert!(Fin::rectangular(0.05, 0.0, 0.002, 180.0, 25.0).is_err());
    }

    #[test]
    fn heat_rate_rejects_nan_theta() {
        let f = sample();
        assert!(f.heat_rate(f64::NAN, TipCondition::Adiabatic).is_err());
    }
}
