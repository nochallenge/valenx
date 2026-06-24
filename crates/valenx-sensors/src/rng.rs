//! Deterministic, seeded PRNG used for all sensor noise in this crate.
//!
//! Simulated sensor readings are only useful for regression testing if they are
//! reproducible, so the crate takes **no `rand` dependency** and instead carries
//! its own tiny [`SplitMix64`] generator — the same deterministic, seeded
//! generator used in `valenx-uq` and `valenx-photogrammetry`. Normals are drawn
//! with the Box–Muller transform; one Box–Muller evaluation yields *two*
//! independent standard normals, so a single cached value avoids wasting half of
//! every transform.

/// Minimal SplitMix64 PRNG.
///
/// SplitMix64 is a well-known, high-quality finaliser-based generator
/// (Steele/Lea/Flood; the seeding RNG of the xoshiro family). Seeding it with
/// the same value reproduces the same stream exactly on every run and platform,
/// so two sensors built with the same seed emit identical noise. It is **not**
/// intended for any cryptographic / security purpose.
#[derive(Debug, Clone)]
pub struct SplitMix64 {
    state: u64,
    /// A standard-normal value produced by Box–Muller but not yet consumed.
    /// Box–Muller emits normals in pairs; this caches the second one.
    cached_normal: Option<f64>,
}

impl SplitMix64 {
    /// Seed the generator.
    #[must_use]
    pub fn new(seed: u64) -> Self {
        Self {
            state: seed,
            cached_normal: None,
        }
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

    /// A single standard-normal `N(0, 1)` sample via Box–Muller, caching the
    /// paired second value for the next call.
    pub fn next_standard_normal(&mut self) -> f64 {
        if let Some(z) = self.cached_normal.take() {
            return z;
        }
        let (z0, z1) = self.box_muller_pair();
        self.cached_normal = Some(z1);
        z0
    }

    /// A `N(mean, std)` sample. A non-positive `std` collapses to `mean`
    /// (degenerate but never `NaN`); callers that require `std > 0` validate that
    /// themselves at construction time.
    pub fn next_normal(&mut self, mean: f64, std: f64) -> f64 {
        if std > 0.0 {
            mean + std * self.next_standard_normal()
        } else {
            mean
        }
    }

    /// A pair of independent standard-normal samples via the Box–Muller
    /// transform.
    fn box_muller_pair(&mut self) -> (f64, f64) {
        // Guard the log against u1 == 0 so the radius is always finite.
        let u1 = self.next_f64().max(f64::MIN_POSITIVE);
        let u2 = self.next_f64();
        let r = (-2.0 * u1.ln()).sqrt();
        let theta = 2.0 * std::f64::consts::PI * u2;
        (r * theta.cos(), r * theta.sin())
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
    fn normals_have_roughly_right_moments() {
        let mut r = SplitMix64::new(0x5EED);
        let n = 200_000;
        let (mut sum, mut sum_sq) = (0.0_f64, 0.0_f64);
        for _ in 0..n {
            let z = r.next_standard_normal();
            sum += z;
            sum_sq += z * z;
        }
        let mean = sum / n as f64;
        let var = sum_sq / n as f64 - mean * mean;
        // Monte-Carlo error is O(1/√n); generous tolerances.
        assert!(mean.abs() < 0.02, "mean = {mean}");
        assert!((var - 1.0).abs() < 0.05, "var = {var}");
    }

    #[test]
    fn zero_std_collapses_to_mean() {
        let mut r = SplitMix64::new(1);
        assert_eq!(r.next_normal(3.5, 0.0), 3.5);
        assert_eq!(r.next_normal(3.5, -1.0), 3.5);
    }
}
