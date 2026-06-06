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

#[cfg(test)]
mod tests {
    use super::*;

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
}
