//! Demography: population-size schedules, drift and effective size.
//!
//! Real populations do not hold a constant census size. A
//! [`DemographicSchedule`] is a piecewise-constant size history played
//! *forward* in time: generation by generation it tells the simulator
//! how many individuals the next generation should hold.
//!
//! Convenience constructors cover the textbook scenarios:
//!
//! - [`DemographicSchedule::constant`] — a flat size.
//! - [`DemographicSchedule::bottleneck`] — a sharp contraction for a
//!   fixed window, then recovery.
//! - [`DemographicSchedule::exponential_growth`] — a geometric
//!   expansion (or decline) at a fixed per-generation rate.
//! - [`DemographicSchedule::piecewise`] — an arbitrary list of
//!   `(generation, size)` change-points.
//!
//! This module also exposes the closed-form **drift** quantities that
//! make popgen popgen: the per-generation loss of heterozygosity
//! `1/(2N)`, the variance of an allele frequency after one generation
//! of binomial sampling, and the expected time to fixation/loss.

use crate::error::{PopgenError, Result};
use serde::{Deserialize, Serialize};

/// A forward-in-time piecewise-constant population-size history.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DemographicSchedule {
    /// `(start_generation, size)` change-points, sorted ascending by
    /// generation. The size of the first segment applies from
    /// generation 0; each later segment's size applies from its start
    /// generation onward.
    changes: Vec<(usize, usize)>,
}

impl DemographicSchedule {
    /// A constant population of `n` individuals for all generations.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] if `n == 0`.
    pub fn constant(n: usize) -> Result<Self> {
        if n == 0 {
            return Err(PopgenError::invalid("n", "size must be positive"));
        }
        Ok(DemographicSchedule {
            changes: vec![(0, n)],
        })
    }

    /// A bottleneck: size `n0` until generation `start`, size
    /// `bottleneck_n` for `duration` generations, then back to `n0`.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] on any zero size or zero duration.
    pub fn bottleneck(
        n0: usize,
        start: usize,
        duration: usize,
        bottleneck_n: usize,
    ) -> Result<Self> {
        if n0 == 0 || bottleneck_n == 0 {
            return Err(PopgenError::invalid("n", "sizes must be positive"));
        }
        if duration == 0 {
            return Err(PopgenError::invalid("duration", "must be positive"));
        }
        Ok(DemographicSchedule {
            changes: vec![(0, n0), (start, bottleneck_n), (start + duration, n0)],
        })
    }

    /// A geometric size change: from `n0` at generation 0, multiplied
    /// by `(1 + rate)` each generation for `generations` steps, sampled
    /// into change-points. A `rate > 0` is expansion, `rate < 0`
    /// decline.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] if `n0 == 0`, `rate <= -1`, or
    /// `generations == 0`.
    pub fn exponential_growth(n0: usize, rate: f64, generations: usize) -> Result<Self> {
        if n0 == 0 {
            return Err(PopgenError::invalid("n0", "size must be positive"));
        }
        if rate <= -1.0 {
            return Err(PopgenError::invalid("rate", "growth rate must exceed -1"));
        }
        if generations == 0 {
            return Err(PopgenError::invalid("generations", "must be positive"));
        }
        let mut changes = Vec::with_capacity(generations);
        let mut size = n0 as f64;
        for g in 0..generations {
            changes.push((g, size.round().max(1.0) as usize));
            size *= 1.0 + rate;
        }
        Ok(DemographicSchedule { changes })
    }

    /// An arbitrary piecewise schedule from `(generation, size)`
    /// change-points. The list is sorted; the earliest change-point's
    /// size covers generation 0 onward.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] on an empty list or a zero size.
    pub fn piecewise(mut changes: Vec<(usize, usize)>) -> Result<Self> {
        if changes.is_empty() {
            return Err(PopgenError::invalid("changes", "schedule is empty"));
        }
        if changes.iter().any(|&(_, n)| n == 0) {
            return Err(PopgenError::invalid("size", "every size must be positive"));
        }
        changes.sort_by_key(|&(g, _)| g);
        Ok(DemographicSchedule { changes })
    }

    /// The population size at `generation`.
    pub fn size_at(&self, generation: usize) -> usize {
        let mut size = self.changes[0].1;
        for &(g, n) in &self.changes {
            if g <= generation {
                size = n;
            } else {
                break;
            }
        }
        size
    }

    /// The size of the founding generation (generation 0).
    pub fn initial_size(&self) -> usize {
        self.changes[0].1
    }

    /// The change-points.
    pub fn change_points(&self) -> &[(usize, usize)] {
        &self.changes
    }
}

/// Closed-form genetic-drift quantities for a diploid population of
/// census size `n` (so `2n` chromosomes).
///
/// The *effective* population size is supplied separately because it
/// usually differs from the census size; if unknown, pass the census
/// size and the formulae reduce to the idealised Wright-Fisher case.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Drift {
    /// Effective (diploid) population size used in the drift formulae.
    pub effective_size: f64,
}

impl Drift {
    /// A drift model with the given diploid effective size `ne`.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] if `ne <= 0`.
    pub fn new(ne: f64) -> Result<Self> {
        if ne <= 0.0 {
            return Err(PopgenError::invalid("effective_size", "must be positive"));
        }
        Ok(Drift { effective_size: ne })
    }

    /// Per-generation fraction of heterozygosity lost to drift:
    /// `1 / (2 * Ne)`.
    pub fn heterozygosity_decay(&self) -> f64 {
        1.0 / (2.0 * self.effective_size)
    }

    /// Expected heterozygosity after `t` generations of drift starting
    /// from `h0`: `H_t = H_0 * (1 - 1/(2 Ne))^t`.
    pub fn expected_heterozygosity(&self, h0: f64, t: usize) -> f64 {
        h0 * (1.0 - self.heterozygosity_decay()).powi(t as i32)
    }

    /// Variance of a derived-allele frequency after **one** generation
    /// of Wright-Fisher binomial sampling, starting from frequency `p`:
    /// `Var = p (1 - p) / (2 Ne)`.
    pub fn allele_frequency_variance(&self, p: f64) -> f64 {
        let p = p.clamp(0.0, 1.0);
        p * (1.0 - p) / (2.0 * self.effective_size)
    }

    /// Probability that a derived allele currently at frequency `p`
    /// eventually reaches fixation under pure drift — which is simply
    /// `p` itself (the classic neutral-fixation result).
    pub fn fixation_probability(&self, p: f64) -> f64 {
        p.clamp(0.0, 1.0)
    }

    /// Expected time (in generations) to fixation of a new neutral
    /// mutation (initial frequency `1/(2Ne)`): approximately `4 Ne`
    /// generations.
    pub fn expected_fixation_time(&self) -> f64 {
        4.0 * self.effective_size
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constant_schedule_is_flat() {
        let s = DemographicSchedule::constant(100).unwrap();
        assert_eq!(s.size_at(0), 100);
        assert_eq!(s.size_at(999), 100);
        assert_eq!(s.initial_size(), 100);
    }

    #[test]
    fn bottleneck_dips_then_recovers() {
        let s = DemographicSchedule::bottleneck(1000, 10, 5, 50).unwrap();
        assert_eq!(s.size_at(0), 1000);
        assert_eq!(s.size_at(9), 1000);
        assert_eq!(s.size_at(10), 50); // bottleneck begins
        assert_eq!(s.size_at(14), 50);
        assert_eq!(s.size_at(15), 1000); // recovered
    }

    #[test]
    fn exponential_growth_increases() {
        let s = DemographicSchedule::exponential_growth(100, 0.1, 20).unwrap();
        assert_eq!(s.size_at(0), 100);
        assert!(s.size_at(19) > s.size_at(0));
        // After 10 generations: 100 * 1.1^10 ~ 259.
        assert!((s.size_at(10) as i64 - 259).abs() <= 2);
    }

    #[test]
    fn piecewise_sorts_change_points() {
        let s = DemographicSchedule::piecewise(vec![(20, 30), (0, 10), (10, 20)]).unwrap();
        assert_eq!(s.size_at(0), 10);
        assert_eq!(s.size_at(10), 20);
        assert_eq!(s.size_at(25), 30);
    }

    #[test]
    fn schedule_rejects_bad_input() {
        assert!(DemographicSchedule::constant(0).is_err());
        assert!(DemographicSchedule::bottleneck(100, 0, 0, 10).is_err());
        assert!(DemographicSchedule::exponential_growth(100, -1.5, 5).is_err());
        assert!(DemographicSchedule::piecewise(vec![]).is_err());
        assert!(DemographicSchedule::piecewise(vec![(0, 0)]).is_err());
    }

    #[test]
    fn drift_heterozygosity_decays() {
        let d = Drift::new(50.0).unwrap();
        // 1/(2*50) = 0.01.
        assert!((d.heterozygosity_decay() - 0.01).abs() < 1e-12);
        let h10 = d.expected_heterozygosity(0.5, 10);
        assert!(h10 < 0.5 && h10 > 0.4);
        // Monotone decline.
        assert!(d.expected_heterozygosity(0.5, 20) < h10);
    }

    #[test]
    fn drift_variance_and_fixation() {
        let d = Drift::new(100.0).unwrap();
        // Var = p(1-p)/(2Ne) = 0.25/200 = 0.00125.
        assert!((d.allele_frequency_variance(0.5) - 0.00125).abs() < 1e-12);
        // Neutral fixation probability equals the frequency.
        assert!((d.fixation_probability(0.3) - 0.3).abs() < 1e-12);
        // Fixation time ~ 4 Ne.
        assert!((d.expected_fixation_time() - 400.0).abs() < 1e-9);
    }

    #[test]
    fn drift_rejects_nonpositive_ne() {
        assert!(Drift::new(0.0).is_err());
        assert!(Drift::new(-5.0).is_err());
    }
}
