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

/// The **minimum stimulating charge** `Q_min = I_rh·chronaxie` (µA·ms/cm²) — the
/// `w → 0` intercept of the Weiss charge–duration line [`weiss_threshold_charge`],
/// from the `rheobase` `I_rh` and `chronaxie_ms`. It is the *charge-axis* asymptote
/// of the strength–duration curve: as the pulse shrinks the threshold *current*
/// diverges as `1/w`, but the *charge* `I_th·w` falls to this finite floor — the
/// least charge that can ever excite the membrane. It complements the rheobase (the
/// current-axis, long-pulse asymptote): short pulses are charge-limited at `Q_min`,
/// long pulses current-limited at `I_rh`. Linear in both factors.
pub fn minimum_stimulating_charge(rheobase: f64, chronaxie_ms: f64) -> f64 {
    rheobase * chronaxie_ms
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

    #[test]
    fn minimum_stimulating_charge_is_the_weiss_intercept() {
        let (rh, cx) = (10.0, 0.5); // 10 µA/cm² rheobase, 0.5 ms chronaxie
        let q_min = minimum_stimulating_charge(rh, cx);
        // Worked point: Q_min = I_rh·chronaxie = 5 µA·ms/cm².
        assert!((q_min - 5.0).abs() < 1e-12, "Q_min = rheobase·chronaxie = 5, got {q_min}");
        // STRONG cross-check: Q_min is the w → 0 intercept of the Weiss line
        // Q(w) = I_rh·(w + chronaxie), so Q_min = Q(w) − I_rh·w for ANY width, and the
        // short-pulse floor lies below the charge required at any finite width.
        for w in [0.1_f64, 0.5, 2.0, 10.0] {
            let q = weiss_threshold_charge(rh, cx, w).unwrap();
            assert!((q_min - (q - rh * w)).abs() < 1e-12, "Q_min = Weiss(w) − I_rh·w at w={w}");
            assert!(q_min < q, "Q_min < Q(w) at w={w}");
        }
        // Linear in both factors.
        assert!((minimum_stimulating_charge(2.0 * rh, cx) - 2.0 * q_min).abs() < 1e-12, "∝ rheobase");
        assert!((minimum_stimulating_charge(rh, 2.0 * cx) - 2.0 * q_min).abs() < 1e-12, "∝ chronaxie");
        // The model's own rheobase()·chronaxie() gives a finite positive Q_min.
        let q_model = minimum_stimulating_charge(rheobase(), chronaxie());
        assert!(q_model > 0.0 && q_model.is_finite(), "model Q_min plausible: {q_model}");
    }
}
