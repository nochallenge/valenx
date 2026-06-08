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

/// The **inverse Nernst relation** — the trans-membrane concentration ratio
/// `c_out/c_in = exp(z·E_rev / V_T)` implied by a measured reversal (equilibrium)
/// potential `reversal_potential_mv` `E_rev` (mV) for an ion of valence `valence` `z`
/// at absolute temperature `temp_k`, with `V_T` the thermal voltage
/// [`thermal_voltage_mv`]. It inverts [`nernst_potential_mv`]: given the equilibrium
/// potential read off an I–V curve, it recovers the out/in concentration gradient that
/// sets it. `E_rev = 0` gives a ratio of `1` (equal concentrations); a negative
/// `E_rev` (for `z > 0`) gives a ratio below `1` (the cation is more concentrated
/// inside). Like [`nernst_potential_mv`] it is total — non-physical input
/// (`temp_k ≤ 0`) yields a non-finite result the caller is expected to guard.
pub fn nernst_concentration_ratio(temp_k: f64, valence: f64, reversal_potential_mv: f64) -> f64 {
    (valence * reversal_potential_mv / thermal_voltage_mv(temp_k)).exp()
}

/// The thermal voltage `R·T/F` in **millivolts** — the Nernst slope for a
/// monovalent ion (≈ 26.7 mV at body temperature, equivalently ≈ 61.5 mV per
/// tenfold concentration ratio).
pub fn thermal_voltage_mv(temp_k: f64) -> f64 {
    1.0e3 * GAS_CONSTANT * temp_k / FARADAY
}

/// The **per-decade Nernst slope** `S = (R·T/F)·ln(10)` in millivolts — the change in a
/// monovalent ion's equilibrium potential per *tenfold* change in its concentration
/// ratio. It is `ln(10) ≈ 2.303` times the [`thermal_voltage_mv`] (the per-e-fold slope),
/// the classic **≈ 59 mV/decade at 25 °C** (≈ 61.5 mV at body temperature) quoted for
/// ion-selective and pH electrodes, and the decade form of the Nernst equation
/// `E = (S/z)·log₁₀(c_out/c_in)`. Like its siblings it is a bare closed form; a
/// non-positive `temp_k` yields a non-physical result the caller is expected to guard.
pub fn nernst_slope_per_decade_mv(temp_k: f64) -> f64 {
    thermal_voltage_mv(temp_k) * std::f64::consts::LN_10
}

/// The **Ussing flux ratio** `M_in/M_out = (c_out/c_in)·exp(−z·V/V_T)` — the ratio
/// of an ion's unidirectional **influx to efflux** across the membrane under
/// combined diffusion and electrical drift (Ussing, 1949), at membrane potential
/// `v_membrane_mv` `V` (mV, inside − outside), valence `z`, outside/inside
/// concentrations `c_out`/`c_in`, and absolute temperature `temp_k`; `V_T` is the
/// thermal voltage [`thermal_voltage_mv`].
///
/// It is the classic test for **passive** (purely electrodiffusive) transport: if
/// an ion's measured unidirectional fluxes obey this ratio it crosses the membrane
/// down its electrochemical gradient alone, while a systematic deviation betrays
/// active transport or carrier coupling. It is the kinetic generalisation of the
/// [`nernst_potential_mv`] equilibrium: at `V = E_Nernst` the two unidirectional
/// fluxes balance and the ratio is exactly `1` (no net flux); below the reversal
/// potential the ratio exceeds `1` (net influx for a cation), above it falls below
/// `1` (net efflux). With `V = 0` it reduces to the bare concentration ratio
/// `c_out/c_in` (pure diffusion), and with equal concentrations to `exp(−z·V/V_T)`
/// (pure drift). Like [`nernst_potential_mv`] it is total: non-physical input
/// (`z = 0`, non-positive concentration) yields a non-finite result the caller is
/// expected to guard.
pub fn ussing_flux_ratio(temp_k: f64, z: f64, c_out: f64, c_in: f64, v_membrane_mv: f64) -> f64 {
    (c_out / c_in) * (-z * v_membrane_mv / thermal_voltage_mv(temp_k)).exp()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nernst_concentration_ratio_inverts_the_nernst_potential() {
        let t = BODY_TEMPERATURE_K;

        // Round-trip: the inverse recovers c_out/c_in exactly from nernst_potential_mv.
        for &(z, c_out, c_in) in &[(1.0_f64, 5.0_f64, 140.0_f64), (1.0, 145.0, 12.0), (2.0, 2.0, 0.1)]
        {
            let e = nernst_potential_mv(t, z, c_out, c_in);
            let ratio = nernst_concentration_ratio(t, z, e);
            assert!(
                (ratio - c_out / c_in).abs() / (c_out / c_in) < 1e-9,
                "inverse Nernst round-trip"
            );
        }

        // Threads thermal_voltage_mv: ratio = exp(z·E/V_T).
        let (z, e) = (2.0_f64, -45.0_f64);
        assert!(
            (nernst_concentration_ratio(t, z, e) - (z * e / thermal_voltage_mv(t)).exp()).abs()
                < 1e-12,
            "ratio = exp(z·E/V_T)"
        );

        // E_rev = 0 → equal concentrations (ratio 1).
        assert!((nernst_concentration_ratio(t, 1.0, 0.0) - 1.0).abs() < 1e-12, "E=0 → 1");

        // Worked: a K⁺-like ion at E = −90 mV (z = 1) has c_out/c_in ≈ 0.034 — about 29×
        // more concentrated inside.
        let r = nernst_concentration_ratio(t, 1.0, -90.0);
        assert!((r - 0.0345).abs() < 1e-2, "K⁺-like ratio ≈ 0.034, got {r}");

        // Monotonic increasing in E_rev for a cation.
        assert!(
            nernst_concentration_ratio(t, 1.0, -30.0) < nernst_concentration_ratio(t, 1.0, 30.0),
            "ratio rises with E_rev"
        );
    }

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
    fn nernst_slope_per_decade_mv_is_the_decade_nernst_slope() {
        use std::f64::consts::LN_10;

        // Threads thermal_voltage_mv: S = V_T · ln(10).
        for &t in &[298.15_f64, 310.15, 273.15] {
            assert!(
                (nernst_slope_per_decade_mv(t) - thermal_voltage_mv(t) * LN_10).abs()
                    <= 1e-12 * nernst_slope_per_decade_mv(t),
                "S = V_T·ln(10)"
            );
        }

        // Threads nernst_potential_mv via the decade form E = (S/z)·log10(c_out/c_in)
        // (ln(10)·log10(x) = ln(x), a different code path than R·T/zF·ln).
        for &(t, c_out, c_in) in &[(310.15_f64, 145.0_f64, 15.0_f64), (298.15, 4.0, 140.0)] {
            let from_slope = nernst_slope_per_decade_mv(t) * (c_out / c_in).log10();
            assert!(
                (nernst_potential_mv(t, 1.0, c_out, c_in) - from_slope).abs() <= 1e-9 * from_slope.abs(),
                "E = S·log10(ratio) for z=1"
            );
            // z = 2 carries the 1/z.
            let from_slope_z2 = nernst_slope_per_decade_mv(t) * (c_out / c_in).log10() / 2.0;
            assert!(
                (nernst_potential_mv(t, 2.0, c_out, c_in) - from_slope_z2).abs()
                    <= 1e-9 * from_slope_z2.abs(),
                "E = (S/z)·log10(ratio)"
            );
        }

        // Textbook: ≈ 59 mV/decade at 25 °C, ≈ 61.5 mV at body temperature.
        assert!((nernst_slope_per_decade_mv(298.15) - 59.16).abs() < 0.3, "≈ 59 mV/decade at 25 °C");
        assert!((nernst_slope_per_decade_mv(310.15) - 61.5).abs() < 0.3, "≈ 61.5 mV at body temp");

        // Linear in absolute temperature.
        assert!(
            (nernst_slope_per_decade_mv(2.0 * 298.15) - 2.0 * nernst_slope_per_decade_mv(298.15))
                .abs()
                <= 1e-12 * nernst_slope_per_decade_mv(2.0 * 298.15),
            "linear in T"
        );
    }

    #[test]
    fn warmer_temperature_scales_the_slope_linearly() {
        // E ∝ T at fixed gradient: doubling absolute temperature doubles E.
        let cold = nernst_potential_mv(150.0, 1.0, 10.0, 1.0);
        let hot = nernst_potential_mv(300.0, 1.0, 10.0, 1.0);
        assert!((hot - 2.0 * cold).abs() < 1e-9, "E should scale with T");
    }

    #[test]
    fn ussing_flux_ratio_is_unity_at_equilibrium() {
        let t = BODY_TEMPERATURE_K;
        // STRONG cross-check: at the Nernst reversal potential the unidirectional
        // influx and efflux balance, so the flux ratio is exactly 1 — for several
        // ions/valences. Ties #215 to nernst_potential_mv AND thermal_voltage_mv.
        for &(z, c_out, c_in) in &[
            (1.0_f64, 4.0_f64, 140.0_f64), // K⁺
            (1.0, 145.0, 15.0),            // Na⁺
            (2.0, 2.0, 1.0e-4),            // Ca²⁺
            (-1.0, 110.0, 10.0),           // Cl⁻
        ] {
            let e = nernst_potential_mv(t, z, c_out, c_in);
            let r = ussing_flux_ratio(t, z, c_out, c_in, e);
            assert!((r - 1.0).abs() < 1e-12, "flux ratio = 1 at E for z={z}: got {r}");
        }
        // At V = 0 the ratio is the bare concentration ratio (pure diffusion, no drift).
        let r0 = ussing_flux_ratio(t, 1.0, 4.0, 140.0, 0.0);
        assert!((r0 - 4.0 / 140.0).abs() < 1e-12, "V=0 → c_out/c_in, got {r0}");
        // Equal concentrations → pure drift exp(−z·V/V_T): 1 at V=0, e at V=−V_T (z=1).
        assert!((ussing_flux_ratio(t, 1.0, 50.0, 50.0, 0.0) - 1.0).abs() < 1e-12, "equal conc, V=0 → 1");
        let vt = thermal_voltage_mv(t);
        assert!(
            (ussing_flux_ratio(t, 1.0, 50.0, 50.0, -vt) - 1.0_f64.exp()).abs() < 1e-9,
            "equal conc, V=−V_T → e"
        );
        // Monotonic in V for a cation (z=1): more negative V → larger influx ratio.
        let (k, m, n) = (
            ussing_flux_ratio(t, 1.0, 4.0, 140.0, -120.0),
            ussing_flux_ratio(t, 1.0, 4.0, 140.0, -95.0),
            ussing_flux_ratio(t, 1.0, 4.0, 140.0, -60.0),
        );
        assert!(k > m && m > n, "cation flux ratio rises as V drops: {k} {m} {n}");
        // Direction sanity: below E_K the ratio exceeds 1 (net K⁺ influx), above it < 1.
        let e_k = nernst_potential_mv(t, 1.0, 4.0, 140.0);
        assert!(ussing_flux_ratio(t, 1.0, 4.0, 140.0, e_k - 20.0) > 1.0, "below E_K → influx");
        assert!(ussing_flux_ratio(t, 1.0, 4.0, 140.0, e_k + 20.0) < 1.0, "above E_K → efflux");
    }
}
