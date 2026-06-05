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
    fn squid_to_mammalian_gating_is_about_thirty_times_faster() {
        // Hodgkin–Huxley squid kinetics (6.3 °C) corrected to body temperature
        // (37 °C) with the typical gating Q10 ≈ 3 ⇒ ~30× faster.
        let factor = q10_scale(1.0, TYPICAL_GATING_Q10, 37.0, 6.3);
        assert!((25.0..=35.0).contains(&factor), "squid→mammal factor {factor}");
        // Warming always speeds the process up relative to the reference.
        assert!(factor > 1.0);
    }
}
