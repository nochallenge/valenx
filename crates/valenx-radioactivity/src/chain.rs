//! Two-member parent -> daughter decay chain (the Bateman solution).
//!
//! A parent nuclide P decays into a daughter D, which is itself
//! radioactive and decays onward (to a stable or unmodelled grand-
//! daughter). With decay constants `lambda_p` and `lambda_d` the
//! populations obey the coupled first-order system
//!
//! ```text
//! dNp/dt = -lambda_p * Np
//! dNd/dt =  lambda_p * Np - lambda_d * Nd
//! ```
//!
//! whose closed-form solution is the **two-member Bateman equation**.
//! With an initial daughter population `Nd0` it is
//!
//! ```text
//! Np(t) = Np0 * exp(-lambda_p * t)
//! Nd(t) = Np0 * lambda_p / (lambda_d - lambda_p)
//!             * (exp(-lambda_p * t) - exp(-lambda_d * t))
//!         + Nd0 * exp(-lambda_d * t)
//! ```
//!
//! Starting from a pure parent (`Nd0 = 0`) the daughter population starts
//! at zero, **grows in** as the parent feeds it, reaches a single maximum
//! and then **decays away** — the characteristic grow-in-then-decay
//! transient. The daughter activity peaks at
//!
//! ```text
//! t_max = ln(lambda_d / lambda_p) / (lambda_d - lambda_p)
//! ```
//!
//! When the parent is the longer-lived member (`lambda_p < lambda_d`) the
//! system settles into **transient equilibrium**, where the daughter and
//! parent activities track each other at the fixed ratio
//! `A_d / A_p -> lambda_d / (lambda_d - lambda_p)`. In the extreme
//! `lambda_p << lambda_d` this ratio tends to one — **secular
//! equilibrium**, `A_d ~ A_p`.

use serde::{Deserialize, Serialize};

use crate::decay::Nuclide;
use crate::error::{RadioactivityError, Result};

/// A parent -> daughter chain of two radioactive nuclides.
///
/// Holds the parent and daughter [`Nuclide`]s. All time-valued outputs
/// share the time unit of the two decay constants (which must match each
/// other). Construct with [`DecayChain::new`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct DecayChain {
    parent: Nuclide,
    daughter: Nuclide,
}

impl DecayChain {
    /// Builds a parent -> daughter chain.
    ///
    /// # Errors
    ///
    /// [`RadioactivityError::DegenerateChain`] if the parent and daughter
    /// decay constants are equal to within `1e-12` (relative): the
    /// Bateman daughter term carries a `1 / (lambda_d - lambda_p)` factor
    /// whose equal-constant limit, though finite, is not special-cased
    /// here. Use distinct constants.
    pub fn new(parent: Nuclide, daughter: Nuclide) -> Result<Self> {
        let lp = parent.decay_constant();
        let ld = daughter.decay_constant();
        // Relative closeness guard; both constants are strictly positive
        // by `Nuclide`'s own invariant, so the denominator is safe.
        if (ld - lp).abs() <= 1e-12 * lp.max(ld) {
            return Err(RadioactivityError::DegenerateChain { lambda: lp });
        }
        Ok(Self { parent, daughter })
    }

    /// The parent nuclide.
    #[inline]
    pub fn parent(&self) -> Nuclide {
        self.parent
    }

    /// The daughter nuclide.
    #[inline]
    pub fn daughter(&self) -> Nuclide {
        self.daughter
    }

    /// Parent population at time `t`, `Np(t) = Np0 * exp(-lambda_p * t)`.
    ///
    /// The parent is unaffected by the daughter, so this is just the
    /// single-nuclide law — see [`Nuclide::remaining`].
    ///
    /// # Errors
    ///
    /// [`RadioactivityError::NonPositive`] if `np0` is not strictly
    /// positive, or if `t` is negative or non-finite.
    pub fn parent_population(&self, np0: f64, t: f64) -> Result<f64> {
        self.parent.remaining(np0, t)
    }

    /// Daughter population at time `t` from the two-member Bateman
    /// equation, given initial parent `np0` and initial daughter `nd0`.
    ///
    /// Pass `nd0 = 0.0` for the pure-parent grow-in case.
    ///
    /// # Errors
    ///
    /// [`RadioactivityError::NonPositive`] if `np0` is not strictly
    /// positive, if `nd0` is negative or non-finite, or if `t` is
    /// negative or non-finite.
    pub fn daughter_population(&self, np0: f64, nd0: f64, t: f64) -> Result<f64> {
        let np0 = RadioactivityError::require_positive("np0", np0)?;
        let nd0 = RadioactivityError::require_non_negative("nd0", nd0)?;
        let t = RadioactivityError::require_non_negative("t", t)?;

        let lp = self.parent.decay_constant();
        let ld = self.daughter.decay_constant();

        let grow_in = np0 * lp / (ld - lp) * ((-lp * t).exp() - (-ld * t).exp());
        let carried = nd0 * (-ld * t).exp();
        Ok(grow_in + carried)
    }

    /// Parent activity at time `t`, `A_p(t) = lambda_p * Np(t)`.
    ///
    /// # Errors
    ///
    /// As [`parent_population`](Self::parent_population).
    pub fn parent_activity(&self, np0: f64, t: f64) -> Result<f64> {
        let np = self.parent_population(np0, t)?;
        Ok(self.parent.decay_constant() * np)
    }

    /// Daughter activity at time `t`, `A_d(t) = lambda_d * Nd(t)`.
    ///
    /// # Errors
    ///
    /// As [`daughter_population`](Self::daughter_population).
    pub fn daughter_activity(&self, np0: f64, nd0: f64, t: f64) -> Result<f64> {
        let nd = self.daughter_population(np0, nd0, t)?;
        Ok(self.daughter.decay_constant() * nd)
    }

    /// Time at which the daughter population (and, for the pure-parent
    /// case, its activity) reaches its maximum,
    /// `t_max = ln(lambda_d / lambda_p) / (lambda_d - lambda_p)`.
    ///
    /// This is the instant where `dNd/dt = 0`, i.e. where the daughter's
    /// formation rate `lambda_p * Np` equals its decay rate
    /// `lambda_d * Nd`. It is independent of `Np0`.
    ///
    /// Returns a positive time whenever the constants differ; the
    /// expression is symmetric and well-defined for both `lambda_d >
    /// lambda_p` and `lambda_d < lambda_p`.
    pub fn daughter_peak_time(&self) -> f64 {
        let lp = self.parent.decay_constant();
        let ld = self.daughter.decay_constant();
        (ld / lp).ln() / (ld - lp)
    }

    /// Limiting daughter-to-parent activity ratio under **transient
    /// equilibrium**, `lambda_d / (lambda_d - lambda_p)`.
    ///
    /// Physically meaningful only when the parent outlives the daughter
    /// (`lambda_p < lambda_d`), where the daughter activity asymptotically
    /// settles to this constant multiple of the parent activity. As
    /// `lambda_p / lambda_d -> 0` the ratio tends to one, recovering
    /// **secular equilibrium** (`A_d ~ A_p`).
    ///
    /// For `lambda_p > lambda_d` (a shorter-lived parent) no equilibrium
    /// forms; the returned value is the bare algebraic expression and the
    /// caller should not interpret it as a physical ratio.
    pub fn transient_activity_ratio(&self) -> f64 {
        let lp = self.parent.decay_constant();
        let ld = self.daughter.decay_constant();
        ld / (ld - lp)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    /// Build a chain from two half-lives (in the same time unit).
    fn chain_from_half_lives(tp: f64, td: f64) -> DecayChain {
        DecayChain::new(
            Nuclide::from_half_life(tp).unwrap(),
            Nuclide::from_half_life(td).unwrap(),
        )
        .unwrap()
    }

    #[test]
    fn degenerate_chain_is_rejected() {
        let n = Nuclide::from_decay_constant(0.5).unwrap();
        let err = DecayChain::new(n, n).expect_err("equal constants must fail");
        assert_eq!(err.code(), "radioactivity.degenerate_chain");
    }

    #[test]
    fn daughter_starts_at_zero_for_pure_parent() {
        // Nd(0) = 0 when nd0 = 0.
        let chain = chain_from_half_lives(10.0, 2.0);
        let nd0_pop = chain.daughter_population(1.0e6, 0.0, 0.0).unwrap();
        assert!(close(nd0_pop, 0.0, 1e-9), "Nd(0) = {nd0_pop}");
    }

    #[test]
    fn daughter_grows_then_decays_transient() {
        // GROUND TRUTH (transient equilibrium, longer-lived parent):
        // the daughter activity rises from 0, hits a single peak, then
        // falls. Sample before, at, and after the analytic peak time.
        let chain = chain_from_half_lives(100.0, 10.0); // lambda_p < lambda_d
        let np0 = 1.0e9;
        let t_peak = chain.daughter_peak_time();
        assert!(t_peak > 0.0, "t_peak = {t_peak}");

        let before = chain.daughter_activity(np0, 0.0, 0.5 * t_peak).unwrap();
        let at = chain.daughter_activity(np0, 0.0, t_peak).unwrap();
        let after = chain.daughter_activity(np0, 0.0, 1.5 * t_peak).unwrap();

        // Strictly rising up to the peak, strictly falling past it.
        assert!(before < at, "before {before} should be < peak {at}");
        assert!(after < at, "after {after} should be < peak {at}");
        assert!(before > 0.0, "daughter activity should be positive");
    }

    #[test]
    fn daughter_peak_time_is_a_stationary_point() {
        // At t_max, dNd/dt = 0, i.e. lambda_p * Np = lambda_d * Nd, i.e.
        // the parent and daughter ACTIVITIES are momentarily equal.
        let chain = chain_from_half_lives(50.0, 7.0);
        let np0 = 5.0e8;
        let t_peak = chain.daughter_peak_time();
        let ap = chain.parent_activity(np0, t_peak).unwrap();
        let ad = chain.daughter_activity(np0, 0.0, t_peak).unwrap();
        assert!(close(ap, ad, 1e-3 * ap), "A_p = {ap}, A_d = {ad} at peak");
    }

    #[test]
    fn peak_time_matches_closed_form() {
        // GROUND TRUTH: t_max = ln(ld/lp)/(ld-lp).
        let chain = chain_from_half_lives(80.0, 6.0);
        let lp = chain.parent().decay_constant();
        let ld = chain.daughter().decay_constant();
        let expected = (ld / lp).ln() / (ld - lp);
        assert!(close(chain.daughter_peak_time(), expected, 1e-12));
    }

    #[test]
    fn transient_equilibrium_activity_ratio() {
        // GROUND TRUTH: long after the peak, with lambda_p < lambda_d,
        // A_d / A_p -> lambda_d / (lambda_d - lambda_p).
        let chain = chain_from_half_lives(1000.0, 20.0);
        let lp = chain.parent().decay_constant();
        let ld = chain.daughter().decay_constant();
        let predicted = ld / (ld - lp);
        assert!(close(chain.transient_activity_ratio(), predicted, 1e-15));

        // The realised ratio at any finite t has the exact closed form
        //   A_d / A_p = predicted * (1 - exp(-(ld - lp) * t)),
        // which is GROUND TRUTH for the approach to transient equilibrium.
        // Pin it tightly at a moderate time, then confirm it converges to
        // `predicted` far out (residual ~ exp(-(ld-lp) t) -> 0).
        let np0 = 1.0e12;
        let t_mid = 10.0 * chain.daughter().half_life();
        let exact = predicted * (1.0 - (-(ld - lp) * t_mid).exp());
        let realised = chain.daughter_activity(np0, 0.0, t_mid).unwrap()
            / chain.parent_activity(np0, t_mid).unwrap();
        assert!(close(realised, exact, 1e-9), "{realised} vs exact {exact}");

        // Many half-lives out, the residual is negligible and the ratio
        // matches the asymptote to a tight tolerance.
        let t_far = 30.0 * chain.daughter().half_life();
        let far = chain.daughter_activity(np0, 0.0, t_far).unwrap()
            / chain.parent_activity(np0, t_far).unwrap();
        assert!(
            close(far, predicted, 1e-6),
            "far ratio {far} vs {predicted}"
        );
    }

    #[test]
    fn secular_equilibrium_ratio_tends_to_one() {
        // GROUND TRUTH: when lambda_p << lambda_d the daughter activity
        // equals the parent activity (secular equilibrium); the limiting
        // ratio differs from 1 by ~ lambda_p / lambda_d.
        let chain = chain_from_half_lives(1.0e6, 1.0); // parent vastly longer
        let lp = chain.parent().decay_constant();
        let ld = chain.daughter().decay_constant();
        let ratio = chain.transient_activity_ratio();
        assert!(close(ratio, 1.0, 1e-5), "secular ratio = {ratio}");
        // The first-order correction is lambda_p / lambda_d.
        assert!(close(ratio - 1.0, lp / ld, 1e-9));

        // And the realised activities equalise out at large t: many
        // daughter half-lives out the transient residual exp(-(ld-lp) t)
        // is negligible, so A_d / A_p -> 1.
        let np0 = 6.022e23;
        let t = 20.0 * chain.daughter().half_life();
        let ap = chain.parent_activity(np0, t).unwrap();
        let ad = chain.daughter_activity(np0, 0.0, t).unwrap();
        assert!(close(ad / ap, 1.0, 1e-3), "A_d / A_p = {}", ad / ap);
    }

    #[test]
    fn bateman_matches_hand_evaluation() {
        // GROUND TRUTH: independent evaluation of the closed form at a
        // fixed point, with lambda_p = 0.1, lambda_d = 0.4, Np0 = 1000,
        // Nd0 = 0, t = 5.
        let chain = DecayChain::new(
            Nuclide::from_decay_constant(0.1).unwrap(),
            Nuclide::from_decay_constant(0.4).unwrap(),
        )
        .unwrap();
        let np0 = 1000.0;
        let t = 5.0;
        // Nd = Np0 * lp/(ld-lp) * (e^{-lp t} - e^{-ld t}).
        let lp: f64 = 0.1;
        let ld: f64 = 0.4;
        let expected = np0 * lp / (ld - lp) * ((-lp * t).exp() - (-ld * t).exp());
        let got = chain.daughter_population(np0, 0.0, t).unwrap();
        assert!(
            close(got, expected, 1e-9),
            "Nd = {got}, expected {expected}"
        );
    }

    #[test]
    fn carried_initial_daughter_decays_exponentially() {
        // With np0 negligible relative to a big nd0 we still must honour
        // the nd0 * exp(-ld t) carried term; check it decays at lambda_d.
        let chain = chain_from_half_lives(1.0e9, 3.0);
        let nd0 = 1.0e6;
        // Tiny parent feed so the carried term dominates the comparison.
        let at_zero = chain.daughter_population(1.0, nd0, 0.0).unwrap();
        let at_one_half_life = chain
            .daughter_population(1.0, nd0, chain.daughter().half_life())
            .unwrap();
        assert!(close(at_zero, nd0, 1.0), "Nd(0) ~ nd0, got {at_zero}");
        // After one daughter half-life the carried part has halved; the
        // parent feed adds a tiny positive amount, so allow a small band.
        assert!(
            at_one_half_life > nd0 / 2.0 && at_one_half_life < nd0 / 2.0 + 100.0,
            "Nd(t_half_d) = {at_one_half_life}"
        );
    }

    #[test]
    fn parent_population_is_pure_exponential() {
        // The parent ignores the daughter: Np(t) must match the lone law.
        let chain = chain_from_half_lives(14.0, 3.0);
        let lone = chain.parent();
        let np0 = 2.0e5;
        for &t in &[0.0, 7.0, 14.0, 28.0] {
            let from_chain = chain.parent_population(np0, t).unwrap();
            let from_lone = lone.remaining(np0, t).unwrap();
            assert!(close(from_chain, from_lone, 1e-9), "t = {t}");
        }
    }

    #[test]
    fn population_inputs_are_validated() {
        let chain = chain_from_half_lives(10.0, 2.0);
        assert!(chain.daughter_population(0.0, 0.0, 1.0).is_err()); // np0 <= 0
        assert!(chain.daughter_population(1.0, -1.0, 1.0).is_err()); // nd0 < 0
        assert!(chain.daughter_population(1.0, 0.0, -1.0).is_err()); // t < 0
    }
}
