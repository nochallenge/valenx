//! **Monte-Carlo engagement analysis** — pure statistics over the existing
//! abstract engagement model.
//!
//! ## What this is (and is not)
//!
//! This module runs the **same** [`crate::Scenario`] many times with the
//! crate's seeded [`crate::SplitMix64`] PRNG and aggregates the per-run
//! [`crate::OutcomeMetrics`] into an [`OutcomeStats`]. It adds **no** new
//! mechanism of any kind: every run is the existing deterministic scenario, and
//! the only stochastic element remains the abstract probability-of-kill (`Pk`)
//! Bernoulli draw that already lives in [`crate::engagement::resolve_pk`].
//!
//! There is **no** lethality model, **no** targeting / fire-control, and **no**
//! kill-chain logic here — this is purely the *statistics of an abstract input
//! probability*, the same dual-use posture as the rest of the crate. It answers
//! operations-research questions ("how often does blue prevail across the random
//! draws, and with what spread of survivors?") by sampling, exactly as a
//! think-tank constructive sim reports a probability of mission success.
//!
//! ## Determinism
//!
//! Reproducibility is preserved end-to-end. A single base seed deterministically
//! produces the **per-run** seeds via one [`crate::SplitMix64`]; the same base
//! seed therefore replays a bit-for-bit identical Monte-Carlo ensemble (and hence
//! identical statistics) on every run and machine. The PRNG is **not** used for
//! any security purpose.

use crate::error::MissionError;
use crate::scenario::Scenario;
use crate::SplitMix64;

/// Which side prevailed in a single run.
///
/// "Prevail" is the abstract operations-research outcome: the side with strictly
/// more survivors at the stop time. Equal survivor counts are a [`Self::Draw`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Prevailed {
    /// Blue had strictly more survivors than red.
    Blue,
    /// Red had strictly more survivors than blue.
    Red,
    /// Both sides ended with the same number of survivors.
    Draw,
}

impl Prevailed {
    /// Classify a single run from its per-side survivor counts: the side with
    /// strictly more survivors prevails; equal counts are a [`Self::Draw`].
    #[must_use]
    pub fn from_survivors(survivors_blue: usize, survivors_red: usize) -> Self {
        match survivors_blue.cmp(&survivors_red) {
            std::cmp::Ordering::Greater => Self::Blue,
            std::cmp::Ordering::Less => Self::Red,
            std::cmp::Ordering::Equal => Self::Draw,
        }
    }
}

/// The mean, sample standard deviation, and 95% confidence interval of the mean
/// for one sampled quantity across the Monte-Carlo ensemble.
///
/// The CI is the standard large-sample interval `mean ± 1.96 · std / √n` (the
/// standard error of the mean). With a single sample the std and CI half-width
/// are `0`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SummaryStat {
    /// Sample mean.
    pub mean: f64,
    /// Sample standard deviation (Bessel-corrected, `n - 1`; `0` for `n == 1`).
    pub std: f64,
    /// Lower bound of the 95% confidence interval for the mean.
    pub ci95_lo: f64,
    /// Upper bound of the 95% confidence interval for the mean.
    pub ci95_hi: f64,
}

impl SummaryStat {
    /// Compute mean / std / 95% CI of the mean from a slice of samples.
    ///
    /// `samples` must be non-empty (the caller guarantees `n >= 1`). The std is
    /// the Bessel-corrected sample standard deviation (`0` when there is a single
    /// sample); the CI half-width is `1.96 · std / √n`.
    #[must_use]
    fn from_samples(samples: &[f64]) -> Self {
        let n = samples.len();
        debug_assert!(n >= 1, "SummaryStat needs at least one sample");
        let mean = samples.iter().sum::<f64>() / n as f64;
        let std = if n >= 2 {
            let var = samples.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / (n as f64 - 1.0);
            var.sqrt()
        } else {
            0.0
        };
        // 95% normal critical value; the standard error is std / sqrt(n).
        let half = 1.96 * std / (n as f64).sqrt();
        Self {
            mean,
            std,
            ci95_lo: mean - half,
            ci95_hi: mean + half,
        }
    }
}

/// Aggregate statistics over an `n`-run Monte-Carlo ensemble of the **same**
/// abstract scenario.
///
/// Outcome probabilities are simple sample frequencies; the survivor / exchange
/// summaries are [`SummaryStat`]s; and [`Self::blue_survivor_histogram`] is the
/// binned distribution of blue survivors (one bin per possible integer count,
/// `0..=max_blue`), so it always sums back to `n`.
#[derive(Debug, Clone, PartialEq)]
pub struct OutcomeStats {
    /// Number of runs sampled (`>= 1`).
    pub runs: usize,
    /// Fraction of runs in which blue prevailed (strictly more blue survivors).
    pub p_blue_prevails: f64,
    /// Fraction of runs in which red prevailed (strictly more red survivors).
    pub p_red_prevails: f64,
    /// Fraction of runs that ended in a draw (equal survivors).
    pub p_draw: f64,
    /// Mean / std / 95% CI of blue survivors per run.
    pub blue_survivors: SummaryStat,
    /// Mean / std / 95% CI of red survivors per run.
    pub red_survivors: SummaryStat,
    /// Mean / std / 95% CI of the exchange ratio per run.
    ///
    /// The exchange ratio is `red_losses / blue_losses` where a side's losses are
    /// `initial - survivors`. To stay finite when blue takes no losses, the
    /// per-run ratio is `red_losses / max(blue_losses, 1)` (a conventional
    /// guard); it is `0` for a run with no red losses. With no blue force at all
    /// the ratio is reported as `0` for that run.
    pub exchange_ratio: SummaryStat,
    /// Histogram of blue survivors: `blue_survivor_histogram[k]` is the number of
    /// runs that ended with exactly `k` blue survivors, for `k` in
    /// `0..=max_blue`. Always sums to [`Self::runs`].
    pub blue_survivor_histogram: Vec<usize>,
}

impl OutcomeStats {
    /// The most-frequent blue-survivor count (the histogram mode), or `0` for an
    /// all-zero histogram. Convenience for a readout.
    #[must_use]
    pub fn modal_blue_survivors(&self) -> usize {
        self.blue_survivor_histogram
            .iter()
            .enumerate()
            .max_by_key(|&(_, &c)| c)
            .map(|(k, _)| k)
            .unwrap_or(0)
    }
}

/// Run `scenario` `runs` times with deterministic per-run seeds and aggregate the
/// per-run outcomes into an [`OutcomeStats`].
///
/// Each run is an independent replay of the **same** abstract engagement model:
/// a fresh [`Scenario`] identical to `scenario` but re-seeded from a per-run seed.
/// The per-run seeds are generated deterministically from the scenario's own seed
/// by a single [`SplitMix64`], so the whole ensemble — and therefore every
/// statistic — is bit-for-bit reproducible for a given base seed. No new
/// lethality / targeting is introduced; this is pure sampling over the existing
/// abstract `Pk` draws.
///
/// "Blue prevails" means strictly more blue survivors than red survivors at the
/// stop time (equal counts are a draw). The exchange-ratio convention is
/// documented on [`OutcomeStats::exchange_ratio`].
///
/// # Errors
///
/// [`MissionError::NonPositive`] if `runs == 0` (a Monte-Carlo with no runs has
/// no defined statistics — fail loud rather than return `NaN`s). Any error from
/// an individual [`Scenario::run`] is propagated verbatim.
pub fn monte_carlo(scenario: &Scenario, runs: usize) -> Result<OutcomeStats, MissionError> {
    if runs == 0 {
        return Err(MissionError::NonPositive {
            quantity: "monte_carlo runs",
            value: 0.0,
        });
    }

    // Initial per-side strengths (constant across runs — the scenario entities
    // are the same every time), used for losses / exchange ratio / histogram size.
    let initial_blue = scenario
        .entities()
        .iter()
        .filter(|e| e.side == crate::entity::Side::Blue)
        .count();
    let initial_red = scenario
        .entities()
        .iter()
        .filter(|e| e.side == crate::entity::Side::Red)
        .count();

    // One PRNG drives the per-run seeds, so the ensemble is reproducible.
    let mut seeder = SplitMix64::new(scenario.base_seed());

    let mut blue_wins = 0usize;
    let mut red_wins = 0usize;
    let mut draws = 0usize;
    let mut blue_samples: Vec<f64> = Vec::with_capacity(runs);
    let mut red_samples: Vec<f64> = Vec::with_capacity(runs);
    let mut exch_samples: Vec<f64> = Vec::with_capacity(runs);
    // One bin per possible blue-survivor count, 0..=initial_blue.
    let mut histogram = vec![0usize; initial_blue + 1];

    for _ in 0..runs {
        let run_seed = seeder.next_u64();
        let run = scenario.with_seed(run_seed);
        let res = run.run()?;
        let sb = res.metrics.survivors_blue;
        let sr = res.metrics.survivors_red;

        match Prevailed::from_survivors(sb, sr) {
            Prevailed::Blue => blue_wins += 1,
            Prevailed::Red => red_wins += 1,
            Prevailed::Draw => draws += 1,
        }

        blue_samples.push(sb as f64);
        red_samples.push(sr as f64);

        // Exchange ratio: red losses per blue loss (guarded to stay finite).
        let blue_losses = initial_blue.saturating_sub(sb);
        let red_losses = initial_red.saturating_sub(sr);
        let exch = red_losses as f64 / (blue_losses.max(1)) as f64;
        exch_samples.push(exch);

        // Histogram bin (sb is in 0..=initial_blue by construction).
        if let Some(slot) = histogram.get_mut(sb) {
            *slot += 1;
        }
    }

    debug_assert_eq!(
        histogram.iter().sum::<usize>(),
        runs,
        "histogram bins must sum to the run count"
    );

    Ok(OutcomeStats {
        runs,
        p_blue_prevails: blue_wins as f64 / runs as f64,
        p_red_prevails: red_wins as f64 / runs as f64,
        p_draw: draws as f64 / runs as f64,
        blue_survivors: SummaryStat::from_samples(&blue_samples),
        red_survivors: SummaryStat::from_samples(&red_samples),
        exchange_ratio: SummaryStat::from_samples(&exch_samples),
        blue_survivor_histogram: histogram,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::entity::{Entity, Mover, Side};
    use nalgebra::Vector3;

    /// A co-located opposing pair: one blue with engagement range + a given `Pk`,
    /// one red sitting inside that range. This is a **single-shot** duel — the
    /// tick step exceeds the stop time, so the scenario evaluates exactly one
    /// engagement tick (at `t = 0`). Per run blue therefore annihilates red with
    /// probability exactly `Pk` (one Bernoulli draw), which makes the Monte-Carlo
    /// statistics analytic. Blue carries no sensor / is never targeted, so it
    /// always survives.
    fn duel(pk: f64, seed: u64) -> Scenario {
        let blue =
            Entity::new(Vector3::zeros(), Side::Blue, Mover::Static, 0.0, 100.0, pk).unwrap();
        let red = Entity::new(
            Vector3::new(10.0, 0.0, 0.0),
            Side::Red,
            Mover::Static,
            0.0,
            0.0,
            0.0,
        )
        .unwrap();
        // tick_dt (2.0) > stop (1.0) -> a single tick at t = 0 -> one Pk draw.
        Scenario::new(vec![blue, red], 1.0, 2.0, seed).unwrap()
    }

    #[test]
    fn duel_is_single_shot_one_draw_per_run() {
        // Pin the helper's "exactly one engagement draw" property: with Pk in
        // (0,1) a single run records exactly one Engagement event (so the
        // per-run kill probability is exactly Pk, the basis of the analytic
        // tests below).
        use crate::scenario::Event;
        let res = duel(0.5, 1).run().unwrap();
        let engagements = res
            .timeline
            .iter()
            .filter(|e| matches!(e.event, Event::Engagement { .. }))
            .count();
        assert_eq!(
            engagements, 1,
            "single-shot duel must have exactly one draw"
        );
    }

    #[test]
    fn runs_zero_is_rejected_not_panicked() {
        let scn = duel(0.5, 1);
        let err = monte_carlo(&scn, 0).unwrap_err();
        assert!(
            matches!(err, MissionError::NonPositive { .. }),
            "runs == 0 must fail loud, got {err:?}"
        );
    }

    #[test]
    fn fixed_seed_is_reproducible() {
        // Same base seed -> bit-identical statistics (the per-run seed stream is
        // deterministic).
        let scn = duel(0.5, 0xABCD_1234);
        let a = monte_carlo(&scn, 500).unwrap();
        let b = monte_carlo(&scn, 500).unwrap();
        assert_eq!(a, b, "fixed seed must reproduce the whole OutcomeStats");
    }

    #[test]
    fn histogram_bins_sum_to_n() {
        let scn = duel(0.5, 7);
        let n = 333;
        let stats = monte_carlo(&scn, n).unwrap();
        assert_eq!(
            stats.blue_survivor_histogram.iter().sum::<usize>(),
            n,
            "histogram bin counts must sum to the number of runs"
        );
        // One blue entity -> survivors are 0 or 1 -> 2 bins.
        assert_eq!(stats.blue_survivor_histogram.len(), 2);
    }

    #[test]
    fn pk_one_is_deterministic_blue_always_prevails() {
        // Pk = 1: blue always kills red on the first tick, so blue prevails in
        // EVERY run regardless of seed. Blue (sensor-less, never targeted) always
        // survives -> mean blue survivors is exactly 1 with a zero-width CI.
        let scn = duel(1.0, 42);
        let stats = monte_carlo(&scn, 200).unwrap();
        assert_eq!(stats.p_blue_prevails, 1.0, "Pk=1 -> blue always prevails");
        assert_eq!(stats.p_red_prevails, 0.0);
        assert_eq!(stats.p_draw, 0.0);
        assert!((stats.blue_survivors.mean - 1.0).abs() < 1e-12);
        assert_eq!(stats.blue_survivors.std, 0.0);
        assert!((stats.red_survivors.mean - 0.0).abs() < 1e-12);
        // Every run is one blue survivor -> all mass in bin 1.
        assert_eq!(stats.blue_survivor_histogram[1], 200);
        assert_eq!(stats.blue_survivor_histogram[0], 0);
    }

    #[test]
    fn pk_zero_is_deterministic_red_never_dies() {
        // Pk = 0: blue never kills red, so both sides keep their single entity
        // every run -> every run is a draw (1 vs 1).
        let scn = duel(0.0, 99);
        let stats = monte_carlo(&scn, 150).unwrap();
        assert_eq!(stats.p_draw, 1.0, "Pk=0 -> 1v1 every run is a draw");
        assert_eq!(stats.p_blue_prevails, 0.0);
        assert_eq!(stats.p_red_prevails, 0.0);
        assert!((stats.red_survivors.mean - 1.0).abs() < 1e-12);
        // No blue losses ever -> exchange ratio is 0 every run.
        assert_eq!(stats.exchange_ratio.mean, 0.0);
    }

    #[test]
    fn blue_prevail_frequency_tracks_pk_within_ci() {
        // Analytic check: in this 1v1 duel blue prevails iff it kills red, which
        // happens with probability exactly Pk per run. So P(blue prevails) is a
        // binomial frequency with true mean Pk; the sample frequency must sit
        // within a few standard errors of Pk.
        let pk = 0.3;
        let n = 20_000;
        let scn = duel(pk, 0x5EED);
        let stats = monte_carlo(&scn, n).unwrap();
        let se = (pk * (1.0 - pk) / n as f64).sqrt();
        assert!(
            (stats.p_blue_prevails - pk).abs() < 4.0 * se,
            "P(blue prevails) = {} should be ~Pk = {pk} (se = {se})",
            stats.p_blue_prevails
        );
        // Red survivors here are Bernoulli(1 - Pk): the analytic mean is 1 - Pk
        // and must lie inside the reported 95% CI for n this large.
        let expected_red_mean = 1.0 - pk;
        assert!(
            stats.red_survivors.ci95_lo <= expected_red_mean
                && expected_red_mean <= stats.red_survivors.ci95_hi,
            "expected red mean {expected_red_mean} must lie in the 95% CI [{}, {}]",
            stats.red_survivors.ci95_lo,
            stats.red_survivors.ci95_hi
        );
        // P(blue) + P(red) + P(draw) is a partition of unity.
        let total = stats.p_blue_prevails + stats.p_red_prevails + stats.p_draw;
        assert!(
            (total - 1.0).abs() < 1e-12,
            "outcome probabilities sum to 1"
        );
    }

    #[test]
    fn summary_stat_single_sample_has_zero_spread() {
        let s = SummaryStat::from_samples(&[3.0]);
        assert_eq!(s.mean, 3.0);
        assert_eq!(s.std, 0.0);
        assert_eq!(s.ci95_lo, 3.0);
        assert_eq!(s.ci95_hi, 3.0);
    }
}
