//! Deterministic, seeded PRNG used for all sampling in this crate.
//!
//! UQ is only useful if its results are reproducible, so the crate takes **no
//! `rand` dependency** and instead carries its own tiny [`SplitMix64`]
//! generator (the same one used in `valenx-photogrammetry`). Normals are drawn
//! with the Box–Muller transform; one Box–Muller evaluation yields *two*
//! independent standard normals, so a single cached value avoids wasting half
//! of every transform.

/// Minimal SplitMix64 PRNG.
///
/// SplitMix64 is a well-known, high-quality finaliser-based generator
/// (Steele/Lea/Flood; the seeding RNG of the xoshiro family). Seeding it with
/// the same value reproduces the same stream exactly on every run and
/// platform. It is **not** intended for any cryptographic / security purpose.
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

    /// Uniform `f64` in `[lo, hi)`. If `hi <= lo` the result collapses to
    /// `lo` (degenerate but never `NaN`); callers that need `lo < hi` validate
    /// that themselves (see [`crate::Distribution`]).
    pub fn next_range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_f64()
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
