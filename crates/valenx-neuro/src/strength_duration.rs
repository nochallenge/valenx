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

/// The analytic Lapicque threshold amplitude `I_th(w) = I_rh·(1 + chronaxie/w)`
/// (µA/cm²) at pulse width `width_ms`, from a measured `rheobase` and
/// `chronaxie_ms`. This is the closed-form strength–duration hyperbola the
/// bisection [`threshold_amplitude`] traces out — fast to evaluate at any width
/// once the two summary parameters are known. `None` for a non-positive width.
pub fn lapicque_threshold(rheobase: f64, chronaxie_ms: f64, width_ms: f64) -> Option<f64> {
    if width_ms > 0.0 {
        Some(rheobase * (1.0 + chronaxie_ms / width_ms))
    } else {
        None
    }
}

/// The Weiss threshold charge `Q(w) = I_th·w = I_rh·(w + chronaxie)`
/// (µA·ms/cm²) at pulse width `width_ms` — the charge–duration line. It rises
/// linearly with width at slope `rheobase`, so the minimum charge
/// `Q_min = I_rh·chronaxie` is the `w → 0` intercept: short pulses are
/// charge-limited, long pulses current-limited. `None` for a non-positive
/// width.
pub fn weiss_threshold_charge(rheobase: f64, chronaxie_ms: f64, width_ms: f64) -> Option<f64> {
    if width_ms > 0.0 {
        Some(rheobase * (width_ms + chronaxie_ms))
    } else {
        None
    }
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

    #[test]
    fn lapicque_and_weiss_models_match_the_textbook_relations() {
        let i_rh = 10.0;
        let cx = 2.0;
        // At width = chronaxie the threshold is exactly 2·rheobase (the definition).
        assert!((lapicque_threshold(i_rh, cx, cx).unwrap() - 2.0 * i_rh).abs() < 1e-12);
        // A long pulse approaches the rheobase from above.
        assert!((lapicque_threshold(i_rh, cx, 1000.0).unwrap() - i_rh).abs() < 0.05);
        // Weiss charge is linear in width: slope = rheobase, intercept = I_rh·cx.
        let q_min = weiss_threshold_charge(i_rh, cx, 1e-9).unwrap(); // w → 0
        assert!((q_min - i_rh * cx).abs() < 1e-6, "Q_min {q_min}");
        let q1 = weiss_threshold_charge(i_rh, cx, 1.0).unwrap();
        let q2 = weiss_threshold_charge(i_rh, cx, 2.0).unwrap();
        assert!((q2 - q1 - i_rh).abs() < 1e-12, "slope = rheobase");
        // The two models agree: Q(w) = I_th(w)·w.
        let w = 0.5;
        let q = weiss_threshold_charge(i_rh, cx, w).unwrap();
        assert!((q - lapicque_threshold(i_rh, cx, w).unwrap() * w).abs() < 1e-12);
        // A non-positive width is undefined in both.
        assert!(lapicque_threshold(i_rh, cx, 0.0).is_none());
        assert!(weiss_threshold_charge(i_rh, cx, -1.0).is_none());
    }
}
