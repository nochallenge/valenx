//! A small deterministic random-number generator.
//!
//! The stochastic parts of an MD engine — the Langevin / Brownian
//! integrator, the Andersen and velocity-rescale thermostats, the
//! Maxwell-Boltzmann velocity initialiser — all need pseudo-random
//! numbers. Rather than pull in the `rand` crate this module ships a
//! tiny self-contained PCG-XSH-RR 64/32 generator (O'Neill 2014). It
//! is the same algorithm `valenx-popgen` and `valenx-phylo` use, kept
//! local here so `valenx-md` has no extra dependency.
//!
//! It is **deterministic and reproducible** — same seed, same stream —
//! which is exactly what a reproducible simulation wants. It is *not*
//! cryptographically secure.

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
        rng.next_u32();
        rng.state = rng.state.wrapping_add(seed);
        rng.next_u32();
        rng
    }

    /// Draws the next 32-bit value.
    pub fn next_u32(&mut self) -> u32 {
        let old = self.state;
        self.state = old.wrapping_mul(Self::MULT).wrapping_add(self.inc);
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
    /// `p` is clamped to `[0, 1]`.
    pub fn bernoulli(&mut self, p: f64) -> bool {
        self.uniform() < p.clamp(0.0, 1.0)
    }

    /// Draws a standard normal variate (mean 0, variance 1) via the
    /// Box-Muller transform.
    pub fn normal(&mut self) -> f64 {
        let u1 = (1.0 - self.uniform()).max(1e-300);
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }

    /// Draws a normal variate with the given `mean` and standard
    /// deviation `sd`. A negative `sd` is treated as its absolute
    /// value.
    pub fn normal_with(&mut self, mean: f64, sd: f64) -> f64 {
        mean + self.normal() * sd.abs()
    }

    /// Draws a chi-squared variate with `dof` degrees of freedom.
    ///
    /// Used by the velocity-rescale (Bussi) thermostat. Implemented as
    /// the sum of `dof` squared standard normals — exact, and `dof` is
    /// small in practice (it equals the system's degree-of-freedom
    /// count, but the thermostat only ever needs `sum_of_noises^2`
    /// plus one extra normal, so this direct form is used for the
    /// small auxiliary draws and a normal approximation for large
    /// `dof`).
    pub fn chi_squared(&mut self, dof: usize) -> f64 {
        if dof == 0 {
            return 0.0;
        }
        if dof <= 64 {
            (0..dof).map(|_| self.normal().powi(2)).sum()
        } else {
            // Wilson-Hilferty: chi^2(k) ~ k * (1 - 2/(9k) + Z*sqrt(2/(9k)))^3.
            let k = dof as f64;
            let z = self.normal();
            let t = 1.0 - 2.0 / (9.0 * k) + z * (2.0 / (9.0 * k)).sqrt();
            (k * t * t * t).max(0.0)
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
    fn normal_mean_and_spread() {
        let mut rng = Rng::new(31);
        let n = 50_000;
        let xs: Vec<f64> = (0..n).map(|_| rng.normal_with(5.0, 2.0)).collect();
        let mean: f64 = xs.iter().sum::<f64>() / n as f64;
        let var: f64 = xs.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n as f64;
        assert!((mean - 5.0).abs() < 0.1, "mean = {mean}");
        assert!((var - 4.0).abs() < 0.3, "var = {var}");
    }

    #[test]
    fn chi_squared_mean_is_about_right() {
        let mut rng = Rng::new(13);
        let n = 40_000;
        // Small-dof branch: E[chi^2(k)] = k.
        let mean: f64 = (0..n).map(|_| rng.chi_squared(10)).sum::<f64>() / n as f64;
        assert!((mean - 10.0).abs() < 0.3, "small mean = {mean}");
        // Large-dof branch.
        let mean2: f64 = (0..n).map(|_| rng.chi_squared(200)).sum::<f64>() / n as f64;
        assert!((mean2 - 200.0).abs() < 5.0, "large mean = {mean2}");
    }
}
