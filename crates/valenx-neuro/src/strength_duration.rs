//! Strength–duration relationship: stimulation threshold vs pulse width.
//!
//! For a space-clamped membrane patch the rectangular-pulse threshold follows
//! Lapicque's hyperbola  `I_th(w) = I_rh · (1 + chronaxie/w)`: long pulses
//! approach the **rheobase** `I_rh`, short pulses need ~constant **charge**.
//! The **chronaxie** is the width at which threshold = 2·rheobase — the
//! standard single-number summary of a fiber's excitability. We recover both
//! from a Hodgkin–Huxley patch by bisection and check the Lapicque shape.

use crate::membrane::{HhMembrane, ImplicitCable};

/// Does a single HH patch fire (overshoot 0 mV) for a rectangular pulse of
/// amplitude `amp` (µA/cm²) and `width` ms?
fn fires(amp: f64, width: f64) -> bool {
    let mut c = ImplicitCable::uniform(1, HhMembrane::at_rest(), 100.0, 238.0, 35.4);
    let peak = c.stimulate_block(amp, 1.0, width, width + 20.0, 0.01);
    peak[0] > 0.0
}

/// Threshold amplitude (µA/cm²) for a rectangular pulse of `width` ms, by
/// bisection — the smallest amplitude that fires an action potential.
pub fn threshold_amplitude(width: f64) -> f64 {
    let mut lo = 0.0_f64;
    let mut hi = 1.0e6_f64;
    for _ in 0..40 {
        let mid = 0.5 * (lo + hi);
        if fires(mid, width) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    hi
}

/// Rheobase: threshold amplitude for a long (50 ms) pulse (µA/cm²).
pub fn rheobase() -> f64 {
    threshold_amplitude(50.0)
}

/// Chronaxie: the pulse width (ms) at which threshold = 2 × rheobase, by
/// bisection on width (threshold falls monotonically with width).
pub fn chronaxie() -> f64 {
    let target = 2.0 * rheobase();
    let mut lo = 0.001_f64;
    let mut hi = 50.0_f64;
    for _ in 0..30 {
        let mid = 0.5 * (lo + hi);
        if threshold_amplitude(mid) > target {
            lo = mid; // threshold still too high → need a wider pulse
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chronaxie_doubles_rheobase_and_is_plausible() {
        let rh = rheobase();
        let cx = chronaxie();
        // HH membrane τ = Rm·Cm ≈ 3.3 ms → chronaxie ≈ 0.5·τ, order ~1–2 ms.
        assert!((0.1..20.0).contains(&cx), "chronaxie in a plausible ms range; got {cx}");
        let i_at_cx = threshold_amplitude(cx);
        assert!(
            (i_at_cx / rh - 2.0).abs() < 0.2,
            "by definition I(chronaxie) ≈ 2·rheobase: {i_at_cx:.1} vs 2×{rh:.1}"
        );
    }

    #[test]
    fn short_pulses_need_constant_charge() {
        // Lapicque/Weiss: at short widths the membrane integrates charge, so the
        // threshold CHARGE Q = I·w is ~constant (independent of width).
        let q_005 = threshold_amplitude(0.05) * 0.05;
        let q_010 = threshold_amplitude(0.10) * 0.10;
        let q_025 = threshold_amplitude(0.25) * 0.25;
        assert!((q_005 / q_010 - 1.0).abs() < 0.1, "charge ~constant: {q_005:.2} vs {q_010:.2}");
        assert!((q_025 / q_010 - 1.0).abs() < 0.1, "charge ~constant: {q_025:.2} vs {q_010:.2}");
    }

    #[test]
    fn rheobase_is_finite_and_positive() {
        let rh = rheobase();
        assert!(rh > 0.0 && rh.is_finite() && rh < 1.0e4, "rheobase plausible; got {rh}");
    }

    #[test]
    fn threshold_rises_as_pulse_shortens() {
        // Lapicque: shorter pulses need more current.
        let i_long = threshold_amplitude(5.0);
        let i_short = threshold_amplitude(0.1);
        assert!(i_short > i_long, "short pulse needs more current: {i_short} vs {i_long}");
    }
}
