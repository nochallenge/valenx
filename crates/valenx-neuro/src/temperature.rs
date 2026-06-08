//! Q10 temperature scaling — how the *rate* of a temperature-dependent process
//! (ion-channel gating, diffusion, reaction kinetics) changes with temperature.
//!
//! A process with temperature coefficient `Q10` runs `Q10×` faster for every
//! 10 °C of warming:
//!
//! ```text
//! rate(T) = rate_ref · Q10^((T − T_ref)/10)
//! ```
//!
//! Hodgkin–Huxley gating kinetics are usually assigned `Q10 ≈ 3`, so the squid
//! rate constants (measured at 6.3 °C) run roughly 30× faster at mammalian body
//! temperature — which is why an uncorrected HH model produces a sluggish spike
//! and must be temperature-scaled before it reproduces a fast mammalian one.

/// A typical `Q10` for ion-channel gating kinetics (≈ 3).
pub const TYPICAL_GATING_Q10: f64 = 3.0;

/// Scale a rate constant `rate_ref` (measured at `ref_temp_c`, °C) to
/// temperature `temp_c` (°C) with temperature coefficient `q10`:
/// `rate = rate_ref · q10^((temp_c − ref_temp_c)/10)`. At the reference
/// temperature the rate is unchanged; each +10 °C multiplies it by `q10`.
pub fn q10_scale(rate_ref: f64, q10: f64, temp_c: f64, ref_temp_c: f64) -> f64 {
    rate_ref * q10.powf((temp_c - ref_temp_c) / 10.0)
}

/// The **temperature at which a `Q10`-scaled process runs at a target rate** (°C) — the
/// inverse of [`q10_scale`] for temperature. Given a `rate` relative to its reference
/// `rate_ref` (measured at `ref_temp_c`, °C) and the temperature coefficient `q10`,
///
/// ```text
/// T = ref_temp_c + 10 · ln(rate / rate_ref) / ln(q10)
/// ```
///
/// recovers the temperature that produces that rate — the thermal-acclimation question
/// (at what temperature does a channel-gating or enzyme process reach a target rate?).
/// At `rate = rate_ref` it returns `ref_temp_c`, and for `q10 > 1` a faster rate maps to a
/// higher temperature. Like its siblings [`q10_scale`] and [`q10_from_rates`] it is a bare
/// closed form: a non-positive `rate` / `rate_ref` (where `ln` is undefined) or `q10 ≤ 0`
/// or `q10 = 1` (where `ln(q10)` is `0`) yields a non-finite result the caller is expected
/// to guard.
pub fn temperature_for_q10_rate(rate: f64, rate_ref: f64, q10: f64, ref_temp_c: f64) -> f64 {
    ref_temp_c + 10.0 * (rate / rate_ref).ln() / q10.ln()
}

/// Recover the **temperature coefficient `Q10`** from two rate measurements — the
/// inverse of [`q10_scale`]. Given a rate `rate_cold` at `temp_cold_c` (°C) and
/// `rate_hot` at `temp_hot_c` (°C),
/// `Q10 = (rate_hot/rate_cold)^(10 / (temp_hot_c − temp_cold_c))`: the factor by
/// which the process speeds up per 10 °C. This is how a `Q10` is *measured* — fit
/// from rates at two temperatures — whereas [`q10_scale`] *applies* a known `Q10`.
/// A temperature-independent process (equal rates) gives `Q10 = 1`; faster-when-
/// warmer kinetics give `Q10 > 1` (biological gating is typically `≈ 2–3`), and a
/// process that slows on warming gives `Q10 < 1`. The result is a property of the
/// pair, unchanged if both the rate and the temperature labels are swapped. Like
/// [`q10_scale`] it is a bare closed form: a zero temperature span or non-positive
/// rate yields a non-finite result the caller is expected to guard.
pub fn q10_from_rates(rate_cold: f64, rate_hot: f64, temp_cold_c: f64, temp_hot_c: f64) -> f64 {
    (rate_hot / rate_cold).powf(10.0 / (temp_hot_c - temp_cold_c))
}

/// The **Arrhenius activation energy** `Eₐ` (J/mol) implied by a temperature
/// coefficient `q10` taken around absolute temperature `temp_k` (K). The Arrhenius
/// rate `k = A·exp(−Eₐ/RT)` predicts a *temperature-dependent* coefficient
/// `Q10 = exp(10·Eₐ / (R·T·(T+10)))`; inverting it at the reference temperature
/// gives
///
/// ```text
/// Eₐ = R · T · (T + 10) · ln(Q10) / 10
/// ```
///
/// where `R` is the molar gas constant (`crate::nernst::GAS_CONSTANT`). This is the
/// standard bridge between the *empirical* `Q10` (how much faster per 10 °C — see
/// [`q10_from_rates`]) and the *mechanistic* activation energy of the underlying
/// reaction, letting a measured channel-gating or reaction `Q10` be reported as an
/// `Eₐ` (gating energies cluster near ~50–100 kJ/mol). A temperature-independent
/// process (`Q10 = 1`) gives `Eₐ = 0`. Like its siblings [`q10_scale`] and
/// [`q10_from_rates`] this is a bare closed form, so a non-positive `q10` (where
/// `ln` is undefined) or `temp_k` yields a non-finite result the caller is expected
/// to guard.
pub fn arrhenius_activation_energy_from_q10(q10: f64, temp_k: f64) -> f64 {
    crate::nernst::GAS_CONSTANT * temp_k * (temp_k + 10.0) * q10.ln() / 10.0
}

/// The **temperature coefficient `Q10`** implied by an Arrhenius activation energy
/// `activation_energy` `Eₐ` (J/mol) taken around absolute temperature `temp_k` `T`
/// (K) — the inverse of [`arrhenius_activation_energy_from_q10`]. Inverting
/// `Eₐ = R·T·(T+10)·ln(Q10)/10` gives
///
/// ```text
/// Q10 = exp(10·Eₐ / (R·T·(T+10)))
/// ```
///
/// with `R` the molar gas constant (`crate::nernst::GAS_CONSTANT`). This is the
/// *mechanistic → empirical* direction (the reverse of
/// [`arrhenius_activation_energy_from_q10`]): given a reaction's activation energy,
/// predict how much faster it runs per 10 °C. A zero activation energy gives
/// `Q10 = 1` (temperature-independent); a larger `Eₐ` is more temperature-sensitive
/// (larger `Q10`). Like its siblings [`q10_scale`] / [`q10_from_rates`] it is a bare
/// closed form, so a non-finite or non-positive `temp_k` yields a non-finite result
/// the caller is expected to guard.
pub fn q10_from_activation_energy(activation_energy: f64, temp_k: f64) -> f64 {
    (10.0 * activation_energy / (crate::nernst::GAS_CONSTANT * temp_k * (temp_k + 10.0))).exp()
}

/// The **Arrhenius rate ratio** `k(T)/k(T_ref) = exp(Ea/R · (1/T_ref − 1/T))` — the
/// exact temperature dependence of a reaction or channel-gating rate with activation
/// energy `activation_energy` `Ea` (J/mol), referenced to absolute temperature
/// `temp_k_ref` and evaluated at `temp_k`, where `R` is the molar gas constant
/// [`crate::nernst::GAS_CONSTANT`]. This is the textbook Arrhenius law as a ratio,
/// without the `Q₁₀` linearisation: a warmer temperature (`temp_k > temp_k_ref`) gives
/// a ratio above `1` for `Ea > 0`. Over a 10 °C step it equals the `Q₁₀`
/// [`q10_from_activation_energy`]. Like the other rate functions it is total —
/// non-physical input (non-positive temperature) yields a non-finite result the caller
/// is expected to guard.
pub fn arrhenius_rate_ratio(activation_energy: f64, temp_k_ref: f64, temp_k: f64) -> f64 {
    (activation_energy / crate::nernst::GAS_CONSTANT * (1.0 / temp_k_ref - 1.0 / temp_k)).exp()
}

/// The **temperature for a target Arrhenius rate ratio** `T = 1 / (1/T_ref − R·ln(ratio)/Ea)`
/// (K) — the inverse of [`arrhenius_rate_ratio`]: the absolute temperature at which a
/// reaction or channel-gating rate with activation energy `activation_energy` `Ea` (J/mol),
/// referenced to `temp_k_ref`, runs `rate_ratio` times its reference rate (`R` is
/// [`crate::nernst::GAS_CONSTANT`]). A `rate_ratio` of `1` returns `temp_k_ref`; for `Ea > 0`
/// a faster target (`> 1`) needs a warmer temperature. Total — non-physical input yields a
/// non-finite result the caller is expected to guard.
pub fn temperature_for_arrhenius_rate_ratio(
    rate_ratio: f64,
    activation_energy: f64,
    temp_k_ref: f64,
) -> f64 {
    1.0 / (1.0 / temp_k_ref - crate::nernst::GAS_CONSTANT * rate_ratio.ln() / activation_energy)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn temperature_for_arrhenius_rate_ratio_inverts_the_rate_ratio() {
        let (ea, t_ref) = (50_000.0_f64, 298.15_f64);

        // Round-trips arrhenius_rate_ratio for several temperatures.
        for &t in &[280.0, 298.15, 310.0, 330.0] {
            let ratio = arrhenius_rate_ratio(ea, t_ref, t);
            assert!(
                (temperature_for_arrhenius_rate_ratio(ratio, ea, t_ref) - t).abs() <= 1e-9 * t,
                "round-trip T = {t}"
            );
        }

        // Identity: rate_ratio = 1 → T_ref (ln 1 = 0).
        assert!(
            (temperature_for_arrhenius_rate_ratio(1.0, ea, t_ref) - t_ref).abs() <= 1e-9 * t_ref,
            "ratio 1 → T_ref"
        );

        // Monotonic (Ea > 0): a faster target needs a warmer T, a slower one a cooler T.
        assert!(temperature_for_arrhenius_rate_ratio(2.0, ea, t_ref) > t_ref, "2× → warmer");
        assert!(temperature_for_arrhenius_rate_ratio(0.5, ea, t_ref) < t_ref, "0.5× → cooler");

        // Worked closed form.
        let expected = 1.0 / (1.0 / t_ref - crate::nernst::GAS_CONSTANT * 2.0_f64.ln() / ea);
        assert!(
            (temperature_for_arrhenius_rate_ratio(2.0, ea, t_ref) - expected).abs()
                <= 1e-9 * expected,
            "worked value"
        );
    }

    #[test]
    fn arrhenius_rate_ratio_threads_the_q10_conversions() {
        use crate::nernst::GAS_CONSTANT as R;
        let t = 300.0;

        // No temperature change → ratio 1.
        assert!((arrhenius_rate_ratio(85_000.0, t, t) - 1.0).abs() < 1e-12, "ratio(T,T) = 1");

        // Definition: k(T)/k(T_ref) = exp(Ea/R·(1/T_ref − 1/T)).
        let ea = 85_000.0;
        let expected = (ea / R * (1.0 / t - 1.0 / 310.0)).exp();
        assert!(
            (arrhenius_rate_ratio(ea, t, 310.0) - expected).abs() <= 1e-12 * expected,
            "Arrhenius definition"
        );

        // Threads q10_from_activation_energy: the rate ratio over a 10 °C step IS the Q10.
        for &e in &[40_000.0_f64, 85_000.0, 120_000.0] {
            let q10 = q10_from_activation_energy(e, t);
            assert!(
                (arrhenius_rate_ratio(e, t, t + 10.0) - q10).abs() <= 1e-12 * q10,
                "ratio(T, T+10) = Q10"
            );
        }

        // Round-trips arrhenius_activation_energy_from_q10: Ea(q10) then ratio over 10° = q10.
        for &q10 in &[2.0_f64, 3.0, 4.5] {
            let e = arrhenius_activation_energy_from_q10(q10, t);
            assert!(
                (arrhenius_rate_ratio(e, t, t + 10.0) - q10).abs() <= 1e-9 * q10,
                "round-trip via activation energy"
            );
        }

        // Monotonic increasing in T for Ea > 0 (warming speeds the reaction).
        assert!(
            arrhenius_rate_ratio(ea, t, 320.0) > arrhenius_rate_ratio(ea, t, 305.0),
            "warmer → faster"
        );
    }

    #[test]
    fn q10_from_activation_energy_inverts_arrhenius() {
        // Eₐ = 0 ⇒ temperature-independent ⇒ Q10 = 1 exactly.
        assert!((q10_from_activation_energy(0.0, 300.0) - 1.0).abs() < 1e-12, "Eₐ=0 → Q10=1");

        // STRONG round-trip cross-check threading arrhenius_activation_energy_from_q10
        // (#233): Q10 → Eₐ → Q10 and Eₐ → Q10 → Eₐ both recover exactly — two
        // independent closed forms (exp of the group vs ln of it).
        for &(q10, t) in &[(3.0_f64, 279.45_f64), (2.0, 293.15), (2.5, 310.15), (1.5, 300.0)] {
            let ea = arrhenius_activation_energy_from_q10(q10, t);
            assert!(
                (q10_from_activation_energy(ea, t) - q10).abs() < 1e-9,
                "Q10→Eₐ→Q10 at q10={q10}"
            );
        }
        for &(ea, t) in &[(74000.0_f64, 279.45_f64), (50000.0, 300.0), (100000.0, 310.0)] {
            let q10 = q10_from_activation_energy(ea, t);
            assert!(
                (arrhenius_activation_energy_from_q10(q10, t) - ea).abs() / ea < 1e-12,
                "Eₐ→Q10→Eₐ at ea={ea}"
            );
        }

        // Physical anchor: HH gating Eₐ ≈ 74 kJ/mol at 6.3 °C ⇒ Q10 ≈ 3.
        let q10_hh = q10_from_activation_energy(74_000.0, 6.3 + 273.15);
        assert!((q10_hh - 3.0).abs() < 0.05, "HH gating Q10 ≈ 3, got {q10_hh}");

        // Monotonic increasing in Eₐ (a higher activation energy is more T-sensitive).
        let t = 300.0;
        assert!(
            q10_from_activation_energy(40_000.0, t) < q10_from_activation_energy(80_000.0, t),
            "Q10 rises with Eₐ"
        );
    }

    #[test]
    fn q10_scaling_matches_the_definition() {
        // No change at the reference temperature.
        assert!((q10_scale(5.0, 3.0, 6.3, 6.3) - 5.0).abs() < 1e-12);
        // +10 °C multiplies by Q10, +20 °C by Q10².
        assert!((q10_scale(1.0, 3.0, 16.3, 6.3) - 3.0).abs() < 1e-12);
        assert!((q10_scale(1.0, 3.0, 26.3, 6.3) - 9.0).abs() < 1e-12);
        // −10 °C divides by Q10.
        assert!((q10_scale(1.0, 3.0, -3.7, 6.3) - 1.0 / 3.0).abs() < 1e-12);
        // Linear in the reference rate.
        assert!((q10_scale(2.0, 3.0, 16.3, 6.3) - 2.0 * q10_scale(1.0, 3.0, 16.3, 6.3)).abs() < 1e-12);
    }

    #[test]
    fn temperature_for_q10_rate_inverts_q10_scale() {
        // Round-trip: recover T from the rate q10_scale produces (the exact inverse).
        for &(r0, q, t, tref) in &[
            (1.0_f64, 2.0_f64, 30.0_f64, 20.0_f64),
            (0.5, 3.0, 15.0, 25.0),
            (2.0, 2.5, 40.0, 10.0),
        ] {
            let recovered = temperature_for_q10_rate(q10_scale(r0, q, t, tref), r0, q, tref);
            assert!((recovered - t).abs() <= 1e-12 * t.abs(), "T = T_ref + 10·ln(r/r0)/ln(q10)");
        }

        // Worked: rate = 2·3² = 18 at Q10 = 3 is two decades above ref → +20 °C.
        assert!(
            (temperature_for_q10_rate(18.0, 2.0, 3.0, 10.0) - 30.0).abs() < 1e-12,
            "2·3² → ref + 20"
        );

        // Identity: rate == rate_ref → the reference temperature (ln 1 = 0).
        assert!(
            (temperature_for_q10_rate(5.0, 5.0, 2.0, 22.0) - 22.0).abs() < 1e-12,
            "rate = ref → T_ref"
        );

        // Monotonic: with Q10 > 1, a higher rate maps to a higher temperature.
        assert!(
            temperature_for_q10_rate(4.0, 1.0, 2.0, 20.0)
                > temperature_for_q10_rate(2.0, 1.0, 2.0, 20.0),
            "faster → warmer"
        );

        // Bare closed form (like q10_scale): non-physical input is non-finite, not 0.
        assert!(temperature_for_q10_rate(-1.0, 1.0, 2.0, 20.0).is_nan(), "ln of negative");
        assert!(!temperature_for_q10_rate(4.0, 1.0, 1.0, 20.0).is_finite(), "q10 = 1 → /0");
        assert!(temperature_for_q10_rate(4.0, -1.0, 2.0, 20.0).is_nan(), "negative rate_ref");
    }

    #[test]
    fn q10_from_rates_inverts_q10_scaling() {
        // STRONG round-trip cross-check: scale a rate by a known Q10 over a span,
        // then fit Q10 back from the two rates — it recovers exactly. Ties #221 to
        // q10_scale (apply ↔ fit): the impl is (k2/k1)^(10/ΔT); the check composes
        // the independent q10_scale forward map.
        for &(q10, t_ref, t) in &[(3.0_f64, 6.3_f64, 37.0_f64), (2.0, 20.0, 30.0), (2.5, 0.0, 25.0)] {
            let k1 = 1.7; // arbitrary reference rate
            let k2 = q10_scale(k1, q10, t, t_ref);
            assert!((q10_from_rates(k1, k2, t_ref, t) - q10).abs() < 1e-9, "round-trip Q10={q10}");
        }
        // A 10 °C span makes Q10 the bare rate ratio (the definition).
        assert!((q10_from_rates(2.0, 6.0, 20.0, 30.0) - 3.0).abs() < 1e-12, "ΔT=10 → Q10 = k2/k1");
        // Equal rates → temperature-independent → Q10 = 1.
        assert!((q10_from_rates(4.0, 4.0, 10.0, 25.0) - 1.0).abs() < 1e-12, "equal rates → Q10 = 1");
        // Faster-when-warmer → Q10 > 1; slower-when-warmer → Q10 < 1.
        assert!(q10_from_rates(1.0, 5.0, 10.0, 30.0) > 1.0, "speeds up → Q10 > 1");
        assert!(q10_from_rates(5.0, 1.0, 10.0, 30.0) < 1.0, "slows down → Q10 < 1");
        // A property of the pair: swapping both the rate and the temperature labels
        // leaves Q10 unchanged.
        assert!(
            (q10_from_rates(2.0, 6.0, 20.0, 30.0) - q10_from_rates(6.0, 2.0, 30.0, 20.0)).abs() < 1e-12,
            "symmetric in the (rate, temp) labelling"
        );
    }

    #[test]
    fn squid_to_mammalian_gating_is_about_thirty_times_faster() {
        // Hodgkin–Huxley squid kinetics (6.3 °C) corrected to body temperature
        // (37 °C) with the typical gating Q10 ≈ 3 ⇒ ~30× faster.
        let factor = q10_scale(1.0, TYPICAL_GATING_Q10, 37.0, 6.3);
        assert!((25.0..=35.0).contains(&factor), "squid→mammal factor {factor}");
        // Warming always speeds the process up relative to the reference.
        assert!(factor > 1.0);
    }

    #[test]
    fn arrhenius_activation_energy_matches_q10_and_threads_q10_from_rates() {
        use crate::nernst::GAS_CONSTANT as R;
        // A temperature-independent process (Q10 = 1) has zero activation energy.
        assert!(arrhenius_activation_energy_from_q10(1.0, 283.15).abs() < 1e-9, "Q10=1 → Ea=0");
        // Monotone in Q10: a steeper temperature dependence ⇒ larger Ea (> 0).
        let ea2 = arrhenius_activation_energy_from_q10(2.0, 283.15);
        let ea3 = arrhenius_activation_energy_from_q10(3.0, 283.15);
        assert!(ea3 > ea2 && ea2 > 0.0, "monotone: Ea(3)={ea3} > Ea(2)={ea2} > 0");
        // Round-trip back to Q10 via the Arrhenius law Q10 = exp(10·Ea/(R·T·(T+10))).
        for &(q10, t_k) in &[(2.0_f64, 290.0_f64), (3.0, 279.45), (2.5, 310.15)] {
            let ea = arrhenius_activation_energy_from_q10(q10, t_k);
            let q10_back = (10.0 * ea / (R * t_k * (t_k + 10.0))).exp();
            assert!((q10_back - q10).abs() / q10 < 1e-9, "round-trip Q10 {q10} → {q10_back}");
        }
        // STRONG non-tautological cross-check threading the independent q10_from_rates:
        // pick an Ea, build two Arrhenius rates 10 K apart, FIT a Q10 from them with
        // q10_from_rates, then recover Ea — closing Ea → rates → Q10 → Ea. The 10 °C
        // span makes q10_from_rates exact, so recovery is to round-off.
        for &(ea, temp_cold_c) in &[(74_000.0_f64, 6.3_f64), (50_000.0, 20.0), (95_000.0, 0.0)] {
            let t_cold_k = temp_cold_c + 273.15;
            let k_cold = (-ea / (R * t_cold_k)).exp(); // Arrhenius rate, A = 1
            let k_hot = (-ea / (R * (t_cold_k + 10.0))).exp();
            let q10 = q10_from_rates(k_cold, k_hot, temp_cold_c, temp_cold_c + 10.0);
            let ea_back = arrhenius_activation_energy_from_q10(q10, t_cold_k);
            assert!((ea_back - ea).abs() / ea < 1e-9, "Ea→rates→Q10→Ea: {ea} vs {ea_back}");
        }
        // Physical anchor: HH gating Q10 ≈ 3 around 6.3 °C ⇒ Ea ≈ 74 kJ/mol.
        let ea_gating = arrhenius_activation_energy_from_q10(3.0, 6.3 + 273.15);
        assert!((60_000.0..=90_000.0).contains(&ea_gating), "gating Ea ≈ 74 kJ/mol, got {ea_gating}");
    }
}
