//! **Abstract, probabilistic** engagement outcomes.
//!
//! ## Dual-use boundary (read this)
//!
//! This module is the one place a mission sim touches "combat", and it stays
//! deliberately **abstract**. There are exactly two mechanisms, and **neither**
//! contains any lethality, targeting, fire-control, or kill-chain model:
//!
//! 1. **A single probabilistic outcome** ([`resolve_pk`]). Probability-of-kill
//!    (`Pk`) is an **INPUT parameter** in `[0, 1]` — supplied by the scenario,
//!    *not* computed from any weapon/target physics. The outcome is one seeded
//!    Bernoulli draw against that input: `hit` with probability `Pk`, `miss`
//!    otherwise. `Pk = 1` always hits, `Pk = 0` never (the benchmark pin). What
//!    a "hit" *means* is just an abstract state change (an entity's `alive` flag)
//!    — there is no damage, penetration, or weapon model of any kind.
//!
//! 2. **Aggregate attrition** ([`lanchester_square_step`], [`lanchester_run`]).
//!    Lanchester's classic *square law* for two aggregate forces `A`, `B`:
//!    `dA/dt = −b·B`, `dB/dt = −a·A`, where `a`, `b` are abstract attrition-rate
//!    coefficients (effectiveness inputs, not lethality models). This is a
//!    century-old operations-research ODE used across logistics and policy
//!    wargaming; it has the conserved quantity `a·A² − b·B²` (the "square law")
//!    and a closed-form `cosh`/`sinh` solution, both of which the tests pin.
//!
//! Nothing here selects a target, computes where to aim, models a warhead, or
//! optimizes for harm. It is force-on-force *bookkeeping* at the level of
//! probabilities and aggregate counts.

use crate::error::{require_non_negative, require_probability, MissionError};
use crate::rng::SplitMix64;

/// The outcome of a single abstract engagement.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EngagementOutcome {
    /// The probabilistic draw succeeded (an abstract "hit").
    Hit,
    /// The probabilistic draw failed (an abstract "miss").
    Miss,
}

/// Resolve one abstract engagement: a seeded Bernoulli draw against the **input**
/// probability-of-kill `pk`.
///
/// `pk` is a parameter, not a computed lethality. The draw uses the supplied
/// deterministic [`SplitMix64`], so the same RNG state reproduces the same
/// outcome. `pk = 1.0` is always [`EngagementOutcome::Hit`]; `pk = 0.0` is always
/// [`EngagementOutcome::Miss`].
///
/// # Errors
///
/// [`MissionError::OutOfProbabilityRange`] if `pk` is outside `[0, 1]`.
pub fn resolve_pk(rng: &mut SplitMix64, pk: f64) -> Result<EngagementOutcome, MissionError> {
    let pk = require_probability("pk", pk)?;
    Ok(if rng.bernoulli(pk) {
        EngagementOutcome::Hit
    } else {
        EngagementOutcome::Miss
    })
}

/// The aggregate strengths of the two sides in a Lanchester model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ForceState {
    /// Strength of side A (non-negative, abstract units).
    pub a: f64,
    /// Strength of side B (non-negative, abstract units).
    pub b: f64,
}

/// One fixed-step RK4 update of the Lanchester **square law**
/// `dA/dt = −b·B`, `dB/dt = −a·A`, clamped so neither force goes negative.
///
/// `a` and `b` are abstract per-capita attrition-rate coefficients (effectiveness
/// inputs). RK4 on this *linear* system is extremely accurate; the conserved
/// quantity `a·A² − b·B²` is held to integration tolerance, which the tests pin.
///
/// # Errors
///
/// [`MissionError::Negative`] if any of `a`, `b`, the strengths, or `dt` is
/// negative / not finite.
pub fn lanchester_square_step(
    state: ForceState,
    a_rate: f64,
    b_rate: f64,
    dt: f64,
) -> Result<ForceState, MissionError> {
    let a_rate = require_non_negative("a_rate", a_rate)?;
    let b_rate = require_non_negative("b_rate", b_rate)?;
    require_non_negative("force A", state.a)?;
    require_non_negative("force B", state.b)?;
    let dt = require_non_negative("dt", dt)?;

    // Derivatives of (A, B): dA = -b_rate * B, dB = -a_rate * A.
    let deriv = |a: f64, b: f64| (-b_rate * b, -a_rate * a);

    let (ka1, kb1) = deriv(state.a, state.b);
    let (ka2, kb2) = deriv(state.a + 0.5 * dt * ka1, state.b + 0.5 * dt * kb1);
    let (ka3, kb3) = deriv(state.a + 0.5 * dt * ka2, state.b + 0.5 * dt * kb2);
    let (ka4, kb4) = deriv(state.a + dt * ka3, state.b + dt * kb3);

    let a = state.a + (dt / 6.0) * (ka1 + 2.0 * ka2 + 2.0 * ka3 + ka4);
    let b = state.b + (dt / 6.0) * (kb1 + 2.0 * kb2 + 2.0 * kb3 + kb4);

    // A force cannot drop below zero; once it hits zero it is annihilated.
    Ok(ForceState {
        a: a.max(0.0),
        b: b.max(0.0),
    })
}

/// Integrate the Lanchester square law from `initial` over total time
/// `duration` in `steps` equal RK4 sub-steps.
///
/// Returns the final [`ForceState`]. With `steps == 0` or `duration == 0` the
/// initial state is returned unchanged (a clean no-op, never a divide-by-zero).
///
/// # Errors
///
/// [`MissionError::Negative`] if `duration` (or any rate / strength) is negative
/// / not finite.
pub fn lanchester_run(
    initial: ForceState,
    a_rate: f64,
    b_rate: f64,
    duration: f64,
    steps: usize,
) -> Result<ForceState, MissionError> {
    let duration = require_non_negative("duration", duration)?;
    if steps == 0 || duration == 0.0 {
        // Validate the rates/strengths even on the no-op path so bad inputs
        // still fail loud.
        let _ = lanchester_square_step(initial, a_rate, b_rate, 0.0)?;
        return Ok(initial);
    }
    let dt = duration / steps as f64;
    let mut state = initial;
    for _ in 0..steps {
        state = lanchester_square_step(state, a_rate, b_rate, dt)?;
    }
    Ok(state)
}

/// The Lanchester square-law invariant `a·A² − b·B²`, conserved by the exact
/// dynamics. Exposed so callers (and tests) can check conservation.
#[must_use]
pub fn square_law_invariant(state: ForceState, a_rate: f64, b_rate: f64) -> f64 {
    a_rate * state.a * state.a - b_rate * state.b * state.b
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-6 * b.abs().max(1.0)
    }

    // ---- BENCHMARK PIN: Pk = 1 always hits, Pk = 0 never (seeded) -----------
    #[test]
    fn pk_one_always_hits_pk_zero_never() {
        let mut rng = SplitMix64::new(12345);
        for _ in 0..10_000 {
            assert_eq!(resolve_pk(&mut rng, 1.0).unwrap(), EngagementOutcome::Hit);
            assert_eq!(resolve_pk(&mut rng, 0.0).unwrap(), EngagementOutcome::Miss);
        }
    }

    #[test]
    fn pk_draw_is_reproducible_and_about_the_right_frequency() {
        // Same seed -> identical outcome sequence.
        let mut a = SplitMix64::new(99);
        let mut b = SplitMix64::new(99);
        for _ in 0..1000 {
            assert_eq!(
                resolve_pk(&mut a, 0.5).unwrap(),
                resolve_pk(&mut b, 0.5).unwrap()
            );
        }
        // Frequency of hits ≈ Pk.
        let mut rng = SplitMix64::new(7);
        let n = 100_000;
        let pk = 0.42;
        let hits = (0..n)
            .filter(|_| resolve_pk(&mut rng, pk).unwrap() == EngagementOutcome::Hit)
            .count();
        assert!((hits as f64 / n as f64 - pk).abs() < 0.01);
    }

    #[test]
    fn pk_out_of_range_is_rejected() {
        let mut rng = SplitMix64::new(1);
        assert!(resolve_pk(&mut rng, 1.5).is_err());
        assert!(resolve_pk(&mut rng, -0.1).is_err());
        assert!(resolve_pk(&mut rng, f64::NAN).is_err());
    }

    // ---- BENCHMARK PIN: square law conserves a*A^2 - b*B^2 ------------------
    #[test]
    fn square_law_conserves_the_invariant() {
        let initial = ForceState { a: 100.0, b: 80.0 };
        let (a_rate, b_rate) = (0.05, 0.03);
        let inv0 = square_law_invariant(initial, a_rate, b_rate);
        // March many steps; the invariant must hold throughout.
        let mut state = initial;
        let dt = 0.01;
        for _ in 0..2000 {
            state = lanchester_square_step(state, a_rate, b_rate, dt).unwrap();
            // Stop checking once a force is annihilated (the clamp breaks the
            // smooth invariant by design).
            if state.a <= 0.0 || state.b <= 0.0 {
                break;
            }
            let inv = square_law_invariant(state, a_rate, b_rate);
            assert!(close(inv, inv0), "invariant drifted: {inv} vs {inv0}");
        }
    }

    // ---- BENCHMARK PIN: matches the closed-form solution --------------------
    #[test]
    fn matches_the_closed_form_solution() {
        // For dA/dt = -b B, dB/dt = -a A the exact solution is
        //   A(t) = A0 cosh(k t) - B0 sqrt(b/a) sinh(k t)
        //   B(t) = B0 cosh(k t) - A0 sqrt(a/b) sinh(k t)
        // with k = sqrt(a b). Verify the RK4 integrator reproduces it before
        // either force hits zero.
        let (a0, b0) = (100.0_f64, 90.0_f64);
        let (a_rate, b_rate) = (0.04_f64, 0.05_f64);
        let k = (a_rate * b_rate).sqrt();
        let t = 1.5;
        let exact_a = a0 * (k * t).cosh() - b0 * (b_rate / a_rate).sqrt() * (k * t).sinh();
        let exact_b = b0 * (k * t).cosh() - a0 * (a_rate / b_rate).sqrt() * (k * t).sinh();
        assert!(
            exact_a > 0.0 && exact_b > 0.0,
            "test epoch must precede annihilation"
        );

        let got = lanchester_run(ForceState { a: a0, b: b0 }, a_rate, b_rate, t, 5000).unwrap();
        assert!(close(got.a, exact_a), "A: {} vs {}", got.a, exact_a);
        assert!(close(got.b, exact_b), "B: {} vs {}", got.b, exact_b);
    }

    #[test]
    fn stronger_side_wins_and_weaker_is_annihilated() {
        // Equal rates, A starts much larger -> B should be wiped out, A survives
        // with strength sqrt(A0^2 - B0^2) by the square law.
        //
        // The closed-form survivor count is sqrt(A0^2 - B0^2) = sqrt(100^2 -
        // 60^2) = 80. We approach it from above: the non-negativity clamp that
        // stops B at 0 injects a little "energy" (the smooth invariant no longer
        // holds across annihilation), so A lands slightly *over* 80 and the
        // overshoot shrinks as the step size shrinks. Pin it to a coarse tol that
        // reflects this honest discretisation artefact rather than pretending the
        // clamp is exact, and check the convergence direction.
        let initial = ForceState { a: 100.0, b: 60.0 };
        let rate = 0.1;
        let coarse = lanchester_run(initial, rate, rate, 100.0, 2_000).unwrap();
        let fine = lanchester_run(initial, rate, rate, 100.0, 200_000).unwrap();
        assert!(
            coarse.b <= 1e-3 && fine.b <= 1e-3,
            "weaker force should be annihilated"
        );
        assert!(
            coarse.a >= 80.0 && fine.a >= 80.0,
            "survivors approach 80 from above"
        );
        // Finer steps -> closer to the analytic 80, and within ~1% even coarse.
        assert!(
            (fine.a - 80.0) < (coarse.a - 80.0),
            "refining converges toward 80"
        );
        assert!(
            (fine.a - 80.0).abs() < 0.8,
            "fine survivors {} ~ 80",
            fine.a
        );
    }

    #[test]
    fn zero_duration_or_zero_steps_is_a_noop() {
        let initial = ForceState { a: 10.0, b: 7.0 };
        assert_eq!(
            lanchester_run(initial, 0.1, 0.1, 0.0, 100).unwrap(),
            initial
        );
        assert_eq!(lanchester_run(initial, 0.1, 0.1, 5.0, 0).unwrap(), initial);
    }

    #[test]
    fn negative_inputs_are_rejected() {
        let initial = ForceState { a: 10.0, b: 7.0 };
        assert!(lanchester_square_step(initial, -0.1, 0.1, 1.0).is_err());
        assert!(lanchester_square_step(ForceState { a: -1.0, b: 7.0 }, 0.1, 0.1, 1.0).is_err());
        assert!(lanchester_run(initial, 0.1, 0.1, -5.0, 10).is_err());
    }
}
