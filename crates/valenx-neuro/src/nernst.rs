//! Nernst equilibrium potential — the membrane voltage at which an ion's
//! electrical and diffusive fluxes exactly balance.
//!
//! `E = (R·T)/(z·F) · ln([ion]_out / [ion]_in)`
//!
//! with the gas constant `R`, absolute temperature `T`, ionic valence `z`, and
//! the Faraday constant `F`. This is the reversal potential each ionic current
//! drives the membrane toward; with physiological gradients it reproduces the
//! familiar resting set-points `E_K ≈ −90 mV`, `E_Na ≈ +60 mV`,
//! `E_Ca ≈ +130 mV`.
//!
//! Exact closed form — no fit. The two concentrations may be in any unit so
//! long as they share it, since only their ratio enters.

/// Universal gas constant, J·mol⁻¹·K⁻¹ (CODATA 2018).
pub const GAS_CONSTANT: f64 = 8.314_462_618;
/// Faraday constant, C·mol⁻¹ (CODATA 2018).
pub const FARADAY: f64 = 96_485.332_12;
/// Mammalian body temperature in kelvin (37 °C).
pub const BODY_TEMPERATURE_K: f64 = 310.15;

/// Nernst equilibrium potential in **millivolts** for an ion of valence `z`
/// at absolute temperature `temp_k` (K), from its outside (`c_out`) and inside
/// (`c_in`) concentrations.
///
/// The concentrations may be in any unit provided both share it — only the
/// ratio `c_out / c_in` matters. `z` must be non-zero and the concentrations
/// strictly positive; otherwise the result is non-finite (`±∞` or `NaN`),
/// which the caller is expected to guard.
pub fn nernst_potential_mv(temp_k: f64, z: f64, c_out: f64, c_in: f64) -> f64 {
    // (R·T)/(z·F) is in volts; the ×1e3 converts to millivolts.
    1.0e3 * (GAS_CONSTANT * temp_k) / (z * FARADAY) * (c_out / c_in).ln()
}

/// The thermal voltage `R·T/F` in **millivolts** — the Nernst slope for a
/// monovalent ion (≈ 26.7 mV at body temperature, equivalently ≈ 61.5 mV per
/// tenfold concentration ratio).
pub fn thermal_voltage_mv(temp_k: f64) -> f64 {
    1.0e3 * GAS_CONSTANT * temp_k / FARADAY
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn thermal_voltage_is_about_26_7_mv_at_body_temp() {
        let vt = thermal_voltage_mv(BODY_TEMPERATURE_K);
        assert!((vt - 26.73).abs() < 0.05, "R·T/F should be ≈ 26.7 mV, got {vt}");
        // Per-decade slope is R·T/F · ln(10) ≈ 61.5 mV.
        let decade = vt * 10.0_f64.ln();
        assert!((decade - 61.5).abs() < 0.3, "decade slope ≈ 61.5 mV, got {decade}");
    }

    #[test]
    fn potassium_reversal_is_near_minus_95_mv() {
        // [K]o = 4 mM, [K]i = 140 mM, z = +1 → the textbook resting E_K.
        let e_k = nernst_potential_mv(BODY_TEMPERATURE_K, 1.0, 4.0, 140.0);
        assert!((e_k - (-95.0)).abs() < 2.0, "E_K should be ≈ −95 mV, got {e_k}");
    }

    #[test]
    fn sodium_reversal_is_near_plus_60_mv() {
        // [Na]o = 145 mM, [Na]i = 15 mM, z = +1.
        let e_na = nernst_potential_mv(BODY_TEMPERATURE_K, 1.0, 145.0, 15.0);
        assert!((e_na - 60.6).abs() < 2.0, "E_Na should be ≈ +61 mV, got {e_na}");
    }

    #[test]
    fn divalent_ion_has_half_the_potential() {
        // Same gradient, valence 2 ⇒ exactly half the monovalent potential.
        let e1 = nernst_potential_mv(BODY_TEMPERATURE_K, 1.0, 100.0, 1.0);
        let e2 = nernst_potential_mv(BODY_TEMPERATURE_K, 2.0, 100.0, 1.0);
        assert!((e2 - 0.5 * e1).abs() < 1e-9, "z=2 should halve E");
    }

    #[test]
    fn equal_concentrations_give_zero() {
        assert!(nernst_potential_mv(BODY_TEMPERATURE_K, 1.0, 50.0, 50.0).abs() < 1e-12);
    }

    #[test]
    fn flipping_the_gradient_flips_the_sign() {
        let a = nernst_potential_mv(BODY_TEMPERATURE_K, 1.0, 10.0, 1.0);
        let b = nernst_potential_mv(BODY_TEMPERATURE_K, 1.0, 1.0, 10.0);
        assert!((a + b).abs() < 1e-9, "swapping out/in should negate E");
    }

    #[test]
    fn warmer_temperature_scales_the_slope_linearly() {
        // E ∝ T at fixed gradient: doubling absolute temperature doubles E.
        let cold = nernst_potential_mv(150.0, 1.0, 10.0, 1.0);
        let hot = nernst_potential_mv(300.0, 1.0, 10.0, 1.0);
        assert!((hot - 2.0 * cold).abs() < 1e-9, "E should scale with T");
    }
}
