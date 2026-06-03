//! A small deterministic random-number generator.
//!
//! Every simulator in this crate — Wright-Fisher forward, Kingman
//! coalescent, the ARG, Approximate Bayesian Computation — needs
//! pseudo-random numbers. Rather than pull in the `rand` crate this
//! module ships a tiny self-contained PCG-XSH-RR 64/32 generator
//! (O'Neill 2014). It is the same algorithm `valenx-phylo` uses, kept
//! local here so `valenx-popgen` has no extra dependency.
//!
//! It is **deterministic and reproducible** — same seed, same stream —
//! which is exactly what reproducible population-genetics simulation
//! wants. It is *not* cryptographically secure.
//!
//! Beyond the uniform draws this generator exposes the four discrete /
//! continuous distributions a forward simulator needs every generation:
//! [`Rng::exponential`] (coalescent waiting times), [`Rng::poisson`]
//! (mutation / recombination event counts), [`Rng::binomial`] (Wright-
//! Fisher allele resampling) and [`Rng::normal`] (quantitative-genetics
//! environmental noise).

/// A PCG-XSH-RR 64->32 pseudo-random generator.
#[derive(Debug, Clone)]
pub struct Rng {
    state: u64,
    inc: u64,
}

impl Rng {
    /// PCG's default multiplier constant.
    const MULT: u64 = 6_364_136_223_846_793_005;

    /// Creates a generator from a 64-bit seed. The stream selector is
    /// derived from the seed so distinct seeds give distinct streams.
    pub fn new(seed: u64) -> Self {
        let mut rng = Rng {
            state: 0,
            inc: (seed << 1) | 1,
        };
        // Standard PCG seeding ritual.
        rng.next_u32();
        rng.state = rng.state.wrapping_add(seed);
        rng.next_u32();
        rng
    }

    /// Draws the next 32-bit value.
    pub fn next_u32(&mut self) -> u32 {
        let old = self.state;
        self.state = old.wrapping_mul(Self::MULT).wrapping_add(self.inc);
        // XSH-RR output permutation.
        let xorshifted = (((old >> 18) ^ old) >> 27) as u32;
        let rot = (old >> 59) as u32;
        xorshifted.rotate_right(rot)
    }

    /// Draws the next 64-bit value (two `u32` draws).
    pub fn next_u64(&mut self) -> u64 {
        let hi = self.next_u32() as u64;
        let lo = self.next_u32() as u64;
        (hi << 32) | lo
    }

    /// Draws a uniform `f64` in the half-open interval `[0, 1)`.
    pub fn uniform(&mut self) -> f64 {
        // 53-bit mantissa precision.
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Draws a uniform `f64` in `[lo, hi)`. Collapses to `lo` if
    /// `hi <= lo`.
    pub fn uniform_range(&mut self, lo: f64, hi: f64) -> f64 {
        if hi <= lo {
            lo
        } else {
            lo + self.uniform() * (hi - lo)
        }
    }

    /// Draws a uniform integer in `[0, n)`. Returns 0 if `n == 0`.
    /// Uses rejection sampling so the result is unbiased.
    pub fn below(&mut self, n: usize) -> usize {
        if n == 0 {
            return 0;
        }
        let n = n as u64;
        let zone = u64::MAX - (u64::MAX % n);
        loop {
            let v = self.next_u64();
            if v < zone {
                return (v % n) as usize;
            }
        }
    }

    /// Draws a Bernoulli outcome: `true` with probability `p`.
    ///
    /// `p` is clamped to `[0, 1]`.
    pub fn bernoulli(&mut self, p: f64) -> bool {
        self.uniform() < p.clamp(0.0, 1.0)
    }

    /// Draws an exponential variate with the given `rate` (`lambda`).
    ///
    /// Mean `1/rate`. A non-positive rate yields `f64::INFINITY`.
    pub fn exponential(&mut self, rate: f64) -> f64 {
        if rate <= 0.0 {
            return f64::INFINITY;
        }
        // Inverse-CDF: -ln(1-u)/lambda. `1 - uniform()` is in (0, 1].
        -(1.0 - self.uniform()).ln() / rate
    }

    /// Draws a standard normal variate (mean 0, variance 1) via the
    /// Box-Muller transform.
    pub fn normal(&mut self) -> f64 {
        let u1 = (1.0 - self.uniform()).max(1e-300);
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }

    /// Draws a normal variate with the given `mean` and standard
    /// deviation `sd`. A negative `sd` is treated as its absolute value.
    pub fn normal_with(&mut self, mean: f64, sd: f64) -> f64 {
        mean + self.normal() * sd.abs()
    }

    /// Draws a Poisson variate with mean `lambda`.
    ///
    /// For small `lambda` (< 30) uses Knuth's multiplicative algorithm;
    /// for large `lambda` uses a normal approximation with a continuity
    /// correction, which is accurate enough for mutation-count draws and
    /// avoids unbounded loops. A non-positive `lambda` returns 0.
    pub fn poisson(&mut self, lambda: f64) -> u64 {
        if lambda <= 0.0 {
            return 0;
        }
        if lambda < 30.0 {
            // Knuth: multiply uniforms until the product drops below
            // e^-lambda.
            let limit = (-lambda).exp();
            let mut k = 0u64;
            let mut product = 1.0;
            loop {
                product *= self.uniform();
                if product <= limit {
                    return k;
                }
                k += 1;
                if k > 1_000_000 {
                    return k; // pathological guard
                }
            }
        } else {
            // Normal approximation: N(lambda, lambda).
            let v = self.normal_with(lambda, lambda.sqrt()) + 0.5;
            if v < 0.0 {
                0
            } else {
                v as u64
            }
        }
    }

    /// Draws a binomial variate: the number of successes in `n`
    /// independent trials each succeeding with probability `p`.
    ///
    /// For small `n` (<= 64) counts Bernoulli trials directly; for
    /// large `n` uses an inversion through individual trials capped so
    /// the loop is bounded. `p` is clamped to `[0, 1]`.
    pub fn binomial(&mut self, n: u64, p: f64) -> u64 {
        let p = p.clamp(0.0, 1.0);
        if n == 0 || p == 0.0 {
            return 0;
        }
        if p == 1.0 {
            return n;
        }
        if n <= 64 {
            (0..n).filter(|_| self.uniform() < p).count() as u64
        } else {
            // Geometric-skip sampling: jump from one success to the
            // next, which is fast when p is small, and exact.
            let mut successes = 0u64;
            let mut index = 0u64;
            let log1mp = (1.0 - p).ln();
            loop {
                // Distance to the next success ~ Geometric(p).
                let gap = (self.uniform().max(1e-300).ln() / log1mp).floor() as u64;
                index += gap + 1;
                if index > n {
                    return successes;
                }
                successes += 1;
            }
        }
    }

    /// Picks an index in `[0, weights.len())` proportional to
    /// `weights`. Returns 0 on an empty or all-zero weight slice.
    pub fn weighted_index(&mut self, weights: &[f64]) -> usize {
        let total: f64 = weights.iter().map(|w| w.max(0.0)).sum();
        if total <= 0.0 {
            return 0;
        }
        let mut target = self.uniform() * total;
        for (i, &w) in weights.iter().enumerate() {
            target -= w.max(0.0);
            if target <= 0.0 {
                return i;
            }
        }
        weights.len() - 1
    }

    /// Fisher-Yates shuffle of a slice in place.
    pub fn shuffle<T>(&mut self, slice: &mut [T]) {
        for i in (1..slice.len()).rev() {
            let j = self.below(i + 1);
            slice.swap(i, j);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn is_deterministic_for_a_seed() {
        let mut a = Rng::new(42);
        let mut b = Rng::new(42);
        for _ in 0..100 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }

    #[test]
    fn distinct_seeds_diverge() {
        let mut a = Rng::new(1);
        let mut b = Rng::new(2);
        let same = (0..50).filter(|_| a.next_u32() == b.next_u32()).count();
        assert!(same < 5, "streams suspiciously correlated: {same}/50");
    }

    #[test]
    fn uniform_stays_in_unit_interval() {
        let mut rng = Rng::new(7);
        for _ in 0..10_000 {
            let u = rng.uniform();
            assert!((0.0..1.0).contains(&u));
        }
    }

    #[test]
    fn exponential_mean_is_about_right() {
        let mut rng = Rng::new(99);
        let n = 50_000;
        let mean: f64 =
            (0..n).map(|_| rng.exponential(2.0)).sum::<f64>() / n as f64;
        assert!((mean - 0.5).abs() < 0.05, "mean = {mean}");
    }

    #[test]
    fn poisson_mean_is_about_right() {
        let mut rng = Rng::new(13);
        let n = 60_000;
        // Small-lambda branch.
        let mean: f64 =
            (0..n).map(|_| rng.poisson(3.0) as f64).sum::<f64>() / n as f64;
        assert!((mean - 3.0).abs() < 0.06, "small mean = {mean}");
        // Large-lambda branch.
        let mean2: f64 =
            (0..n).map(|_| rng.poisson(100.0) as f64).sum::<f64>() / n as f64;
        assert!((mean2 - 100.0).abs() < 1.0, "large mean = {mean2}");
    }

    #[test]
    fn binomial_mean_is_about_right() {
        let mut rng = Rng::new(21);
        let n = 40_000;
        // Small-n branch.
        let mean: f64 =
            (0..n).map(|_| rng.binomial(20, 0.3) as f64).sum::<f64>() / n as f64;
        assert!((mean - 6.0).abs() < 0.1, "small mean = {mean}");
        // Large-n branch.
        let mean2: f64 = (0..n)
            .map(|_| rng.binomial(1000, 0.1) as f64)
            .sum::<f64>()
            / n as f64;
        assert!((mean2 - 100.0).abs() < 1.0, "large mean = {mean2}");
    }

    #[test]
    fn binomial_edges() {
        let mut rng = Rng::new(5);
        assert_eq!(rng.binomial(0, 0.5), 0);
        assert_eq!(rng.binomial(10, 0.0), 0);
        assert_eq!(rng.binomial(10, 1.0), 10);
    }

    #[test]
    fn normal_mean_and_spread() {
        let mut rng = Rng::new(31);
        let n = 50_000;
        let xs: Vec<f64> = (0..n).map(|_| rng.normal_with(5.0, 2.0)).collect();
        let mean: f64 = xs.iter().sum::<f64>() / n as f64;
        let var: f64 =
            xs.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n as f64;
        assert!((mean - 5.0).abs() < 0.1, "mean = {mean}");
        assert!((var - 4.0).abs() < 0.3, "var = {var}");
    }

    #[test]
    fn weighted_index_respects_weights() {
        let mut rng = Rng::new(5);
        let weights = [0.0, 10.0, 0.0];
        for _ in 0..200 {
            assert_eq!(rng.weighted_index(&weights), 1);
        }
    }

    #[test]
    fn shuffle_is_a_permutation() {
        let mut rng = Rng::new(8);
        let mut v: Vec<usize> = (0..50).collect();
        rng.shuffle(&mut v);
        v.sort_unstable();
        assert_eq!(v, (0..50).collect::<Vec<_>>());
    }
}
