//! A tiny deterministic pseudo-random generator.
//!
//! The read simulators ([`crate::simulate`]) and the subsampling
//! utilities ([`crate::util::subsample`]) need reproducible randomness
//! without pulling the `rand` crate (a heavy multi-crate dependency the
//! Round-6 budget keeps off the tree). [`Rng`] is a `SplitMix64`
//! generator — a well-known, statistically decent, 64-bit-state
//! splittable PRNG. It is **not** cryptographic; it is exactly what a
//! reproducible scientific simulation wants: same seed in, same stream
//! out, on every platform.

/// A seeded `SplitMix64` pseudo-random generator.
///
/// `SplitMix64` is the seeding generator recommended alongside
/// `xoshiro`; on its own it is a fine general-purpose stream. Every
/// method is deterministic in the seed.
#[derive(Clone, Debug)]
pub struct Rng {
    state: u64,
}

impl Rng {
    /// Builds a generator from a 64-bit seed. The same seed always
    /// yields the same stream.
    pub fn new(seed: u64) -> Self {
        // Avoid the all-zero fixed point degenerating the first draws.
        Rng {
            state: seed ^ 0x9E37_79B9_7F4A_7C15,
        }
    }

    /// Draws the next raw 64-bit value (the `SplitMix64` step).
    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Draws a `f64` uniformly in the half-open interval `[0, 1)`.
    pub fn next_f64(&mut self) -> f64 {
        // Top 53 bits → mantissa precision, exactly the standard idiom.
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Draws a `usize` uniformly in `[0, n)`. Returns `0` when `n == 0`.
    pub fn below(&mut self, n: usize) -> usize {
        if n == 0 {
            return 0;
        }
        (self.next_u64() % n as u64) as usize
    }

    /// Draws an `i64` uniformly in the inclusive range `[lo, hi]`.
    /// Returns `lo` when `hi <= lo`.
    pub fn range_i64(&mut self, lo: i64, hi: i64) -> i64 {
        if hi <= lo {
            return lo;
        }
        let span = (hi - lo + 1) as u64;
        lo + (self.next_u64() % span) as i64
    }

    /// Returns `true` with probability `p` (clamped to `[0, 1]`).
    pub fn chance(&mut self, p: f64) -> bool {
        self.next_f64() < p.clamp(0.0, 1.0)
    }

    /// Draws a standard-normal `f64` (mean 0, variance 1) via the
    /// Box-Muller transform. Each call consumes two uniforms.
    pub fn next_gaussian(&mut self) -> f64 {
        // Guard the log against an exact zero.
        let u1 = self.next_f64().max(1e-12);
        let u2 = self.next_f64();
        (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }

    /// Draws a normally distributed `f64` with the given `mean` and
    /// `std_dev`. A negative `std_dev` is treated as its absolute value.
    pub fn next_normal(&mut self, mean: f64, std_dev: f64) -> f64 {
        mean + self.next_gaussian() * std_dev.abs()
    }

    /// Performs an in-place Fisher-Yates shuffle of `slice`.
    pub fn shuffle<T>(&mut self, slice: &mut [T]) {
        let n = slice.len();
        for i in (1..n).rev() {
            let j = self.below(i + 1);
            slice.swap(i, j);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn same_seed_same_stream() {
        let mut a = Rng::new(42);
        let mut b = Rng::new(42);
        for _ in 0..32 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }

    #[test]
    fn different_seeds_diverge() {
        let mut a = Rng::new(1);
        let mut b = Rng::new(2);
        assert_ne!(a.next_u64(), b.next_u64());
    }

    #[test]
    fn uniforms_in_unit_interval() {
        let mut r = Rng::new(7);
        for _ in 0..1000 {
            let x = r.next_f64();
            assert!((0.0..1.0).contains(&x), "{x} out of range");
        }
    }

    #[test]
    fn below_is_bounded() {
        let mut r = Rng::new(7);
        for _ in 0..1000 {
            assert!(r.below(10) < 10);
        }
        assert_eq!(r.below(0), 0);
    }

    #[test]
    fn range_is_inclusive_and_bounded() {
        let mut r = Rng::new(99);
        for _ in 0..1000 {
            let v = r.range_i64(-5, 5);
            assert!((-5..=5).contains(&v));
        }
        assert_eq!(r.range_i64(3, 3), 3);
    }

    #[test]
    fn gaussian_is_roughly_centred() {
        let mut r = Rng::new(123);
        let n = 5000;
        let mean: f64 = (0..n).map(|_| r.next_gaussian()).sum::<f64>() / n as f64;
        assert!(mean.abs() < 0.1, "mean drifted: {mean}");
    }

    #[test]
    fn shuffle_preserves_multiset() {
        let mut r = Rng::new(5);
        let mut v: Vec<i32> = (0..50).collect();
        r.shuffle(&mut v);
        v.sort_unstable();
        assert_eq!(v, (0..50).collect::<Vec<_>>());
    }
}
