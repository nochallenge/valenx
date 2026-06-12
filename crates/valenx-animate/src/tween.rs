//! Easing curves.

use serde::{Deserialize, Serialize};

/// Easing mode applied when interpolating *into* a keyframe.
#[derive(Copy, Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub enum TweenMode {
    /// Linear lerp.
    #[default]
    Linear,
    /// Accelerate from rest (quadratic).
    EaseIn,
    /// Decelerate into rest (quadratic).
    EaseOut,
    /// Smooth-step: ease in then ease out.
    EaseInOut,
    /// Hermite cubic with zero tangent at both ends. This is the
    /// special case of [`TweenMode::Hermite`] with `m0 == m1 == 0`
    /// and, for that tangent pair, is mathematically identical to
    /// [`TweenMode::EaseInOut`]'s smooth-step (`3t² − 2t³`).
    Cubic,
    /// **Real parameterised Hermite cubic** (Phase 29.5). Interpolates
    /// `0 → 1` over `t ∈ [0, 1]` using the standard Hermite basis
    /// with caller-supplied **end tangents** `m0` (slope leaving the
    /// start) and `m1` (slope arriving at the end). This is the genuine
    /// Hermite tween — `EaseIn`-like with `m0 = 0, m1 > 0`,
    /// `EaseOut`-like with `m0 > 0, m1 = 0`, an overshoot with a large
    /// `m1`, etc. With `m0 = m1 = 0` it reduces to [`TweenMode::Cubic`].
    Hermite {
        /// Tangent (slope) at `t = 0`.
        m0: f64,
        /// Tangent (slope) at `t = 1`.
        m1: f64,
    },
}

impl TweenMode {
    /// Map normalized `t` in `[0, 1]` through the easing curve.
    pub fn apply(self, t: f64) -> f64 {
        let t = t.clamp(0.0, 1.0);
        match self {
            TweenMode::Linear => t,
            TweenMode::EaseIn => t * t,
            TweenMode::EaseOut => 1.0 - (1.0 - t).powi(2),
            TweenMode::EaseInOut => t * t * (3.0 - 2.0 * t),
            // Zero-tangent Hermite ≡ smooth-step.
            TweenMode::Cubic => hermite(t, 0.0, 1.0, 0.0, 0.0),
            // Real Hermite with caller-supplied end tangents.
            TweenMode::Hermite { m0, m1 } => hermite(t, 0.0, 1.0, m0, m1),
        }
    }

    /// Short label for UI dropdowns.
    pub fn label(self) -> &'static str {
        match self {
            TweenMode::Linear => "Linear",
            TweenMode::EaseIn => "EaseIn",
            TweenMode::EaseOut => "EaseOut",
            TweenMode::EaseInOut => "EaseInOut",
            TweenMode::Cubic => "Cubic",
            TweenMode::Hermite { .. } => "Hermite",
        }
    }
}

/// Cubic Hermite interpolation between `p0` and `p1` over `t ∈ [0, 1]`
/// with end tangents `m0` (at `p0`) and `m1` (at `p1`).
///
/// Uses the standard Hermite basis:
/// `h00 = 2t³−3t²+1`, `h10 = t³−2t²+t`,
/// `h01 = −2t³+3t²`, `h11 = t³−t²`, giving
/// `p(t) = h00·p0 + h10·m0 + h01·p1 + h11·m1`.
pub fn hermite(t: f64, p0: f64, p1: f64, m0: f64, m1: f64) -> f64 {
    let t2 = t * t;
    let t3 = t2 * t;
    let h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
    let h10 = t3 - 2.0 * t2 + t;
    let h01 = -2.0 * t3 + 3.0 * t2;
    let h11 = t3 - t2;
    h00 * p0 + h10 * m0 + h01 * p1 + h11 * m1
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn linear_is_identity() {
        assert!((TweenMode::Linear.apply(0.0) - 0.0).abs() < 1e-12);
        assert!((TweenMode::Linear.apply(0.5) - 0.5).abs() < 1e-12);
        assert!((TweenMode::Linear.apply(1.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn ease_in_slow_start() {
        let v = TweenMode::EaseIn.apply(0.5);
        assert!(v < 0.5); // slower than linear at midpoint
    }

    #[test]
    fn ease_out_fast_start() {
        let v = TweenMode::EaseOut.apply(0.5);
        assert!(v > 0.5);
    }

    #[test]
    fn smooth_step_hits_endpoints() {
        assert!((TweenMode::EaseInOut.apply(0.0) - 0.0).abs() < 1e-12);
        assert!((TweenMode::EaseInOut.apply(1.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn clamping_works() {
        assert!((TweenMode::Linear.apply(-1.0) - 0.0).abs() < 1e-12);
        assert!((TweenMode::Linear.apply(2.0) - 1.0).abs() < 1e-12);
    }

    // --- Phase 29.5 Hermite tests ---

    #[test]
    fn cubic_is_zero_tangent_hermite() {
        // Cubic must equal a Hermite with both tangents zero — and
        // that pair coincides with smooth-step.
        for &t in &[0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0] {
            let cubic = TweenMode::Cubic.apply(t);
            let hzero = TweenMode::Hermite { m0: 0.0, m1: 0.0 }.apply(t);
            let smooth = TweenMode::EaseInOut.apply(t);
            assert!(
                (cubic - hzero).abs() < 1e-12,
                "cubic != zero-tangent hermite"
            );
            assert!(
                (cubic - smooth).abs() < 1e-12,
                "zero-tangent hermite != smoothstep"
            );
        }
    }

    #[test]
    fn hermite_hits_endpoints_regardless_of_tangents() {
        // Whatever the tangents, p(0) == 0 and p(1) == 1.
        for &(m0, m1) in &[(0.0, 0.0), (1.0, 1.0), (3.0, -2.0), (-5.0, 5.0)] {
            let h = TweenMode::Hermite { m0, m1 };
            assert!(h.apply(0.0).abs() < 1e-12, "Hermite p(0) should be 0");
            assert!(
                (h.apply(1.0) - 1.0).abs() < 1e-12,
                "Hermite p(1) should be 1"
            );
        }
    }

    #[test]
    fn hermite_respects_start_tangent() {
        // A large positive m0 makes the curve climb fast near t=0 —
        // a real tangent effect smooth-step (m0=0) cannot produce.
        let steep = TweenMode::Hermite { m0: 3.0, m1: 0.0 };
        let flat = TweenMode::Hermite { m0: 0.0, m1: 0.0 };
        let small_t = 0.05;
        assert!(
            steep.apply(small_t) > flat.apply(small_t) + 0.05,
            "a larger start tangent must lift the curve earlier"
        );
    }

    #[test]
    fn hermite_can_overshoot() {
        // A large *start* tangent shoots the curve up past 1 in the
        // interior; the (small) arrival tangent then pulls it back down
        // to exactly 1 at t=1 — a genuine interior overshoot. (A large
        // *arrival* tangent, by contrast, only steepens the approach at
        // t=1 and does not lift the curve above 1 within [0, 1].)
        let overshoot = TweenMode::Hermite { m0: 4.0, m1: 0.0 };
        let peak = (1..20)
            .map(|i| overshoot.apply(i as f64 / 20.0))
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(
            peak > 1.0,
            "a large start tangent should overshoot past the endpoint, got peak {peak}"
        );
        // ...and it still lands exactly on the endpoint.
        assert!((overshoot.apply(1.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn hermite_helper_matches_linear_endpoints() {
        // The bare helper between arbitrary p0/p1.
        assert!((hermite(0.0, 2.0, 9.0, 0.0, 0.0) - 2.0).abs() < 1e-12);
        assert!((hermite(1.0, 2.0, 9.0, 0.0, 0.0) - 9.0).abs() < 1e-12);
    }

    #[test]
    fn label_covers_hermite() {
        assert_eq!(TweenMode::Hermite { m0: 1.0, m1: 0.0 }.label(), "Hermite");
    }
}
