//! A small deterministic random-number generator for stochastic
//! structure sampling.
//!
//! Boltzmann sampling needs uniform `[0, 1)` draws, and the sampler
//! must be *reproducible* — the same seed must always yield the same
//! set of structures. This is a hand-rolled `splitmix64` generator:
//! tiny, fast, fully deterministic, with no external dependency. It is
//! not cryptographic and is not meant to be — reproducibility is the
//! goal.

/// A deterministic `splitmix64` pseudo-random generator.
#[derive(Clone, Debug)]
pub struct Rng {
    state: u64,
}

impl Rng {
    /// Creates a generator seeded with `seed`.
    pub fn new(seed: u64) -> Self {
        // Avoid the all-zero state degenerating.
        Rng {
            state: seed ^ 0x9E37_79B9_7F4A_7C15,
        }
    }

    /// Returns the next raw 64-bit value (the `splitmix64` step).
    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Returns a uniform `f64` in `[0, 1)`.
    pub fn next_f64(&mut self) -> f64 {
        // 53 bits of mantissa precision.
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deterministic_for_a_seed() {
        let a: Vec<u64> = {
            let mut r = Rng::new(42);
            (0..8).map(|_| r.next_u64()).collect()
        };
        let b: Vec<u64> = {
            let mut r = Rng::new(42);
            (0..8).map(|_| r.next_u64()).collect()
        };
        assert_eq!(a, b);
    }

    #[test]
    fn different_seeds_differ() {
        let mut r1 = Rng::new(1);
        let mut r2 = Rng::new(2);
        assert_ne!(r1.next_u64(), r2.next_u64());
    }

    #[test]
    fn floats_in_unit_interval() {
        let mut r = Rng::new(7);
        for _ in 0..1000 {
            let x = r.next_f64();
            assert!((0.0..1.0).contains(&x), "draw {x} out of [0,1)");
        }
    }
}
