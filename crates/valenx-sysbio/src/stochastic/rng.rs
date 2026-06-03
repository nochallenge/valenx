//! Deterministic pseudo-random number generator.
//!
//! The three stochastic simulators (SSA, tau-leaping, next-reaction)
//! and the global-sampling sensitivity analysis all need random
//! numbers. To keep the dependency tree minimal — and, more
//! importantly, to keep every simulation **bit-for-bit reproducible**
//! from its seed — the crate ships its own small PRNG rather than
//! pulling in `rand`.
//!
//! [`Rng`] is `splitmix64`: a single 64-bit state advanced by a
//! well-tested mixing function (Steele, Lea & Flood, 2014). It is not
//! cryptographic, but its statistical quality is more than adequate
//! for Monte-Carlo kinetics, and it has a long period and excellent
//! avalanche behaviour. Every reaction-trajectory result in this
//! crate is therefore exactly reproducible given the seed.

/// A small, fast, deterministic PRNG (`splitmix64`).
#[derive(Debug, Clone)]
pub struct Rng {
    state: u64,
}

impl Rng {
    /// Create a generator from a 64-bit seed. Seed `0` is remapped so
    /// the generator never starts from a degenerate state.
    pub fn new(seed: u64) -> Self {
        Rng {
            state: if seed == 0 { 0x9E3779B97F4A7C15 } else { seed },
        }
    }

    /// Next raw 64-bit value.
    pub fn next_u64(&mut self) -> u64 {
        // splitmix64 step.
        self.state = self.state.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }

    /// A uniform `f64` in the half-open interval `[0, 1)`.
    pub fn uniform(&mut self) -> f64 {
        // 53 random mantissa bits → exactly representable.
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// A uniform `f64` in `[lo, hi)`.
    pub fn uniform_range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.uniform()
    }

    /// An exponentially distributed sample with the given rate
    /// `lambda` (mean `1/lambda`). The inverse-CDF transform — the
    /// waiting-time draw at the heart of the Gillespie SSA.
    pub fn exponential(&mut self, lambda: f64) -> f64 {
        if lambda <= 0.0 {
            return f64::INFINITY;
        }
        // 1 - u avoids ln(0) when u == 0.
        -(1.0 - self.uniform()).ln() / lambda
    }

    /// A Poisson-distributed non-negative integer with mean `mean`,
    /// using Knuth's multiplicative algorithm for small means and a
    /// normal approximation (rounded, clamped) for large means. Drives
    /// the tau-leaping reaction-firing counts.
    pub fn poisson(&mut self, mean: f64) -> u64 {
        if mean <= 0.0 {
            return 0;
        }
        if mean < 30.0 {
            let l = (-mean).exp();
            let mut k = 0u64;
            let mut p = 1.0;
            loop {
                k += 1;
                p *= self.uniform();
                if p < l {
                    return k - 1;
                }
                if k > 10_000 {
                    return k; // defensive ceiling
                }
            }
        } else {
            // Gaussian approximation: N(mean, mean).
            let g = self.normal() * mean.sqrt() + mean;
            g.round().max(0.0) as u64
        }
    }

    /// A standard-normal sample via the Box-Muller transform.
    pub fn normal(&mut self) -> f64 {
        let u1 = (1.0 - self.uniform()).max(1e-300);
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deterministic_from_seed() {
        let mut a = Rng::new(42);
        let mut b = Rng::new(42);
        for _ in 0..1000 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }

    #[test]
    fn uniform_in_unit_interval() {
        let mut r = Rng::new(7);
        for _ in 0..10_000 {
            let u = r.uniform();
            assert!((0.0..1.0).contains(&u));
        }
    }

    #[test]
    fn uniform_mean_is_about_half() {
        let mut r = Rng::new(99);
        let n = 100_000;
        let mean: f64 = (0..n).map(|_| r.uniform()).sum::<f64>() / n as f64;
        assert!((mean - 0.5).abs() < 0.01, "mean {mean}");
    }

    #[test]
    fn exponential_mean_matches_rate() {
        let mut r = Rng::new(123);
        let lambda = 2.0;
        let n = 200_000;
        let mean: f64 = (0..n).map(|_| r.exponential(lambda)).sum::<f64>() / n as f64;
        // Expected mean is 1/lambda = 0.5.
        assert!((mean - 0.5).abs() < 0.01, "mean {mean}");
    }

    #[test]
    fn poisson_mean_matches_parameter() {
        let mut r = Rng::new(555);
        let n = 100_000;
        let mean: f64 =
            (0..n).map(|_| r.poisson(4.0) as f64).sum::<f64>() / n as f64;
        assert!((mean - 4.0).abs() < 0.05, "mean {mean}");
    }

    #[test]
    fn poisson_large_mean_uses_approximation() {
        let mut r = Rng::new(8);
        let n = 50_000;
        let mean: f64 =
            (0..n).map(|_| r.poisson(100.0) as f64).sum::<f64>() / n as f64;
        assert!((mean - 100.0).abs() < 1.0, "mean {mean}");
    }

    #[test]
    fn normal_is_roughly_standard() {
        let mut r = Rng::new(2024);
        let n = 100_000;
        let xs: Vec<f64> = (0..n).map(|_| r.normal()).collect();
        let mean: f64 = xs.iter().sum::<f64>() / n as f64;
        let var: f64 =
            xs.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n as f64;
        assert!(mean.abs() < 0.02, "mean {mean}");
        assert!((var - 1.0).abs() < 0.05, "var {var}");
    }
}
