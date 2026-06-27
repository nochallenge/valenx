//! Deterministic, seeded PRNG used for every stochastic draw in this crate.
//!
//! A constructive simulation is only useful for analysis if it is **exactly
//! reproducible** — the same seed must replay the same timeline on every run and
//! machine, or you cannot regression-test it or compare two scenario variants.
//! This crate therefore takes **no `rand` dependency** and carries its own tiny
//! [`SplitMix64`] generator (the same deterministic finaliser-based generator
//! used in `valenx-sensors`, `valenx-uq`, and `valenx-photogrammetry`). It is
//! **not** intended for any cryptographic / security purpose.

/// Minimal SplitMix64 PRNG.
///
/// SplitMix64 is a well-known, high-quality finaliser-based generator
/// (Steele/Lea/Flood; the seeding RNG of the xoshiro family). Seeding it with
/// the same value reproduces the same stream exactly on every run and platform.
#[derive(Debug, Clone)]
pub struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    /// Seed the generator.
    #[must_use]
    pub fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    /// Next raw 64-bit output.
    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform `f64` in `[0, 1)` from the top 53 bits.
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// A Bernoulli trial: `true` with probability `p`, `false` otherwise.
    ///
    /// `p <= 0` always returns `false` and `p >= 1` always returns `true`
    /// (so a probability-of-kill of exactly `0` *never* kills and exactly `1`
    /// *always* kills, deterministically — the boundary the tests pin). For an
    /// interior `p` the draw is `u < p` with `u` uniform in `[0, 1)`; because
    /// the uniform can equal `0` but never `1`, the boundaries are exact.
    pub fn bernoulli(&mut self, p: f64) -> bool {
        if p <= 0.0 {
            false
        } else if p >= 1.0 {
            true
        } else {
            self.next_f64() < p
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn same_seed_same_stream() {
        let mut a = SplitMix64::new(0xABCD_1234);
        let mut b = SplitMix64::new(0xABCD_1234);
        for _ in 0..1000 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }

    #[test]
    fn uniform_in_unit_interval() {
        let mut r = SplitMix64::new(7);
        for _ in 0..10_000 {
            let x = r.next_f64();
            assert!((0.0..1.0).contains(&x));
        }
    }

    #[test]
    fn bernoulli_boundaries_are_deterministic() {
        let mut r = SplitMix64::new(42);
        for _ in 0..1000 {
            assert!(!r.bernoulli(0.0), "p=0 must never fire");
            assert!(r.bernoulli(1.0), "p=1 must always fire");
            assert!(!r.bernoulli(-0.5), "p<0 clamps to never");
            assert!(r.bernoulli(1.5), "p>1 clamps to always");
        }
    }

    #[test]
    fn bernoulli_frequency_is_about_right() {
        let mut r = SplitMix64::new(0x5EED);
        let n = 200_000;
        let p = 0.3;
        let hits = (0..n).filter(|_| r.bernoulli(p)).count();
        let freq = hits as f64 / n as f64;
        // Monte-Carlo error is O(1/√n); a generous tolerance.
        assert!((freq - p).abs() < 0.01, "freq = {freq}");
    }
}
