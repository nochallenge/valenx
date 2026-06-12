//! A small deterministic random-number generator.
//!
//! Phylogenetic simulation (coalescent, birth-death, sequence
//! evolution) and bootstrap resampling all need pseudo-random numbers.
//! Rather than pull in the `rand` crate, this module ships a tiny
//! self-contained PCG-XSH-RR 64/32 generator (O'Neill 2014).
//!
//! It is **deterministic and reproducible** — same seed, same stream —
//! which is exactly what reproducible scientific simulation wants. It
//! is *not* cryptographically secure; do not use it for anything that
//! needs unpredictability.

/// A PCG-XSH-RR 64→32 pseudo-random generator.
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

    /// Draws a uniform `f64` in `[lo, hi)`.
    ///
    /// If `hi <= lo` the result collapses to `lo`.
    pub fn uniform_range(&mut self, lo: f64, hi: f64) -> f64 {
        if hi <= lo {
            lo
        } else {
            lo + self.uniform() * (hi - lo)
        }
    }

    /// Draws a uniform integer in `[0, n)`.
    ///
    /// Returns 0 if `n == 0`. Uses rejection sampling so the result is
    /// unbiased.
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

    /// Draws an exponential variate with the given `rate` (`λ`).
    ///
    /// Mean `1/rate`. `rate` must be positive; a non-positive rate
    /// yields `f64::INFINITY`.
    pub fn exponential(&mut self, rate: f64) -> f64 {
        if rate <= 0.0 {
            return f64::INFINITY;
        }
        // Inverse-CDF: -ln(1-u)/λ. `1 - uniform()` is in (0, 1].
        -(1.0 - self.uniform()).ln() / rate
    }

    /// Draws a standard normal variate (mean 0, variance 1) via the
    /// Box-Muller transform.
    pub fn normal(&mut self) -> f64 {
        let u1 = (1.0 - self.uniform()).max(1e-300);
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }

    /// Draws a Gamma variate with the given `shape` (`k`) and `scale`
    /// (`θ`), using Marsaglia & Tsang's method.
    ///
    /// Mean `k·θ`. `shape` and `scale` must be positive.
    pub fn gamma(&mut self, shape: f64, scale: f64) -> f64 {
        if shape <= 0.0 || scale <= 0.0 {
            return 0.0;
        }
        if shape < 1.0 {
            // Boost: Γ(k) = Γ(k+1) · U^(1/k).
            let g = self.gamma(shape + 1.0, 1.0);
            return g * self.uniform().max(1e-300).powf(1.0 / shape) * scale;
        }
        let d = shape - 1.0 / 3.0;
        let c = 1.0 / (9.0 * d).sqrt();
        loop {
            let x = self.normal();
            let v = (1.0 + c * x).powi(3);
            if v <= 0.0 {
                continue;
            }
            let u = self.uniform();
            if u < 1.0 - 0.0331 * x.powi(4) || u.ln() < 0.5 * x * x + d * (1.0 - v + v.ln()) {
                return d * v * scale;
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
    fn below_is_in_range_and_covers() {
        let mut rng = Rng::new(11);
        let mut seen = [false; 6];
        for _ in 0..1000 {
            let v = rng.below(6);
            assert!(v < 6);
            seen[v] = true;
        }
        assert!(seen.iter().all(|&s| s), "did not cover 0..6");
    }

    #[test]
    fn exponential_mean_is_about_right() {
        let mut rng = Rng::new(99);
        let n = 50_000;
        let mean: f64 = (0..n).map(|_| rng.exponential(2.0)).sum::<f64>() / n as f64;
        // Expected mean 1/λ = 0.5; allow Monte-Carlo slack.
        assert!((mean - 0.5).abs() < 0.05, "mean = {mean}");
    }

    #[test]
    fn gamma_mean_is_about_right() {
        let mut rng = Rng::new(123);
        let n = 50_000;
        let mean: f64 = (0..n).map(|_| rng.gamma(2.0, 1.5)).sum::<f64>() / n as f64;
        // Expected mean k·θ = 3.0.
        assert!((mean - 3.0).abs() < 0.1, "mean = {mean}");
    }

    #[test]
    fn weighted_index_respects_weights() {
        let mut rng = Rng::new(5);
        let weights = [0.0, 10.0, 0.0];
        for _ in 0..200 {
            // All weight on index 1.
            assert_eq!(rng.weighted_index(&weights), 1);
        }
    }
}
