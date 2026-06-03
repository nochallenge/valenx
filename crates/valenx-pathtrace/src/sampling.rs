//! Pseudo-random numbers and the importance-sampling primitives the
//! Monte-Carlo integrator draws from.
//!
//! # The RNG
//!
//! [`Rng`] is a **PCG32** generator (O'Neill 2014) — a 64-bit LCG state
//! advanced through a permutation step. It is small (one `u64`), fast,
//! statistically excellent for graphics, and — crucially —
//! deterministically seedable per pixel, so a render is exactly
//! reproducible. The standard library's `HashMap` RNG is not seedable
//! and `rand` would be a heavyweight dependency for what is a dozen
//! lines of well-known code.
//!
//! # Importance sampling
//!
//! Path tracing converges fastest when the directions it samples are
//! distributed *like the function being integrated*. The diffuse
//! reflectance integral has a `cos θ` weight, so [`cosine_hemisphere`]
//! draws directions with a `pdf = cos θ / π` — the diffuse BRDF's own
//! shape — and the `cos θ` and the `1/π` then cancel out of the
//! estimator, leaving the noise-free `albedo` factor (this is why a
//! Lambert bounce in [`crate::tracer`] multiplies throughput by just
//! the albedo).

/// A small, fast, seedable pseudo-random generator (PCG32).
///
/// One generator per pixel (seeded from the pixel coordinate + the
/// sample index) makes a render deterministic and trivially
/// parallelisable — no shared mutable state.
#[derive(Clone, Debug)]
pub struct Rng {
    /// The 64-bit LCG state.
    state: u64,
    /// The LCG increment — must be odd; selects the stream.
    inc: u64,
}

impl Rng {
    /// Seed a generator from two integers.
    ///
    /// `seq` selects the random *stream* (two RNGs with the same
    /// `state` but different `seq` are decorrelated); `seed` is the
    /// initial state. The path tracer passes the pixel index as `seq`
    /// so neighbouring pixels never share a sequence.
    pub fn new(seed: u64, seq: u64) -> Rng {
        // PCG seeding ritual: set the increment (forced odd), run one
        // step, fold in the seed, run another step.
        let mut rng = Rng {
            state: 0,
            inc: (seq << 1) | 1,
        };
        rng.next_u32();
        rng.state = rng.state.wrapping_add(seed);
        rng.next_u32();
        rng
    }

    /// Draw the next raw 32-bit value (the PCG32 output function).
    #[inline]
    pub fn next_u32(&mut self) -> u32 {
        let old = self.state;
        // LCG step with the PCG multiplier.
        self.state = old
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(self.inc);
        // Output permutation: xorshift-high then rotate.
        let xorshifted = (((old >> 18) ^ old) >> 27) as u32;
        let rot = (old >> 59) as u32;
        xorshifted.rotate_right(rot)
    }

    /// A uniform `f32` in the half-open interval `[0, 1)`.
    ///
    /// Built from the top 24 bits of a raw `u32` (the `f32` mantissa
    /// width) so every representable value is equally likely and the
    /// result never reaches exactly 1.0.
    #[inline]
    pub fn next_f32(&mut self) -> f32 {
        // 24-bit mantissa → divide by 2^24.
        (self.next_u32() >> 8) as f32 * (1.0 / 16_777_216.0)
    }
}

use crate::math::{ortho_basis, Vec3};

/// Draw a direction on the hemisphere about the unit normal `n`,
/// distributed with a **cosine-weighted** density `pdf = cos θ / π`.
///
/// Returns the world-space direction. Cosine-weighted sampling is the
/// importance-sampling match for a Lambertian (diffuse) BRDF: because
/// the sample density already carries the `cos θ` of the rendering
/// equation and the `1/π` of the Lambert BRDF, both cancel in the
/// Monte-Carlo estimator and a diffuse bounce reduces to a plain
/// `albedo` multiply with no extra variance.
///
/// # Method
///
/// Malley's method: sample a point uniformly on the unit disk
/// (concentric mapping of the two input randoms `u1`, `u2`) and project
/// it up onto the hemisphere. The projection automatically produces the
/// `cos θ` distribution. The disk point's `(x, y)` become the tangent
/// components and `z = √(1 − x² − y²)` the normal component, all rotated
/// into world space by the orthonormal basis of `n`.
#[inline]
pub fn cosine_hemisphere(n: Vec3, u1: f32, u2: f32) -> Vec3 {
    // Concentric disk mapping (Shirley) — lower distortion than the
    // naive `r = √u`, `θ = 2πu` polar map.
    let a = 2.0 * u1 - 1.0;
    let b = 2.0 * u2 - 1.0;
    let (r, phi) = if a == 0.0 && b == 0.0 {
        (0.0, 0.0)
    } else if a * a > b * b {
        (a, std::f32::consts::FRAC_PI_4 * (b / a))
    } else {
        (
            b,
            std::f32::consts::FRAC_PI_2 - std::f32::consts::FRAC_PI_4 * (a / b),
        )
    };
    let dx = r * phi.cos();
    let dy = r * phi.sin();
    // Project the disk point onto the hemisphere: z carries the cosine.
    let dz = (1.0 - dx * dx - dy * dy).max(0.0).sqrt();
    // Rotate the local (dx, dy, dz) into the world frame of `n`.
    let (tangent, bitangent) = ortho_basis(n);
    tangent
        .scale(dx)
        .add(bitangent.scale(dy))
        .add(n.scale(dz))
        .normalized()
        .unwrap_or(n)
}

/// The probability density of [`cosine_hemisphere`] for a direction at
/// angle cosine `cos_theta` from the normal: `cos θ / π`.
///
/// Exposed so the integrator (and tests) can reason about the estimator
/// weight explicitly.
#[inline]
pub fn cosine_hemisphere_pdf(cos_theta: f32) -> f32 {
    cos_theta.max(0.0) * std::f32::consts::FRAC_1_PI
}

/// Sample a point uniformly inside the unit disk, returned as
/// `(x, y)` — used to jitter camera rays for anti-aliasing and (scaled)
/// for depth-of-field lens sampling.
#[inline]
pub fn uniform_disk(u1: f32, u2: f32) -> (f32, f32) {
    let r = u1.max(0.0).sqrt();
    let theta = std::f32::consts::TAU * u2;
    (r * theta.cos(), r * theta.sin())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::vec3;

    #[test]
    fn rng_is_deterministic_for_a_fixed_seed() {
        let mut a = Rng::new(42, 7);
        let mut b = Rng::new(42, 7);
        for _ in 0..100 {
            assert_eq!(a.next_u32(), b.next_u32(), "same seed must replay");
        }
    }

    #[test]
    fn rng_streams_decorrelate() {
        // Different `seq` values give different sequences.
        let mut a = Rng::new(1, 1);
        let mut b = Rng::new(1, 2);
        let mut same = 0;
        for _ in 0..64 {
            if a.next_u32() == b.next_u32() {
                same += 1;
            }
        }
        assert!(same < 4, "distinct streams should rarely coincide");
    }

    #[test]
    fn next_f32_stays_in_unit_interval() {
        let mut rng = Rng::new(123, 4);
        for _ in 0..10_000 {
            let x = rng.next_f32();
            assert!((0.0..1.0).contains(&x), "f32 sample {x} out of [0,1)");
        }
    }

    #[test]
    fn next_f32_is_roughly_uniform() {
        // A coarse bucket test: 8 bins, 80k samples, each bin within
        // 15% of the expected 1/8.
        let mut rng = Rng::new(999, 1);
        let mut bins = [0u32; 8];
        let n = 80_000;
        for _ in 0..n {
            let b = (rng.next_f32() * 8.0) as usize;
            bins[b.min(7)] += 1;
        }
        let expected = n as f32 / 8.0;
        for (i, &c) in bins.iter().enumerate() {
            let rel = (c as f32 - expected).abs() / expected;
            assert!(rel < 0.15, "bin {i} count {c} skewed (rel {rel})");
        }
    }

    #[test]
    fn cosine_hemisphere_samples_lie_in_the_normal_hemisphere() {
        // Every cosine-weighted sample must be a unit vector on the
        // *positive* side of the normal.
        let n = vec3(0.0, 0.0, 1.0);
        let mut rng = Rng::new(7, 7);
        for _ in 0..5000 {
            let d = cosine_hemisphere(n, rng.next_f32(), rng.next_f32());
            assert!((d.length() - 1.0).abs() < 1e-4, "sample not unit");
            assert!(d.dot(n) >= -1e-4, "sample below the hemisphere");
        }
    }

    #[test]
    fn cosine_hemisphere_mean_cosine_is_two_thirds() {
        // The expected cosine of a cosine-weighted hemisphere sample is
        // ∫cosθ·(cosθ/π)dω = 2/3. A Monte-Carlo average of many
        // samples must land near it — this confirms the *distribution*,
        // not just the support.
        let n = vec3(0.0, 0.0, 1.0);
        let mut rng = Rng::new(2024, 11);
        let mut sum = 0.0f64;
        let count = 200_000;
        for _ in 0..count {
            let d = cosine_hemisphere(n, rng.next_f32(), rng.next_f32());
            sum += d.dot(n).max(0.0) as f64;
        }
        let mean = sum / count as f64;
        assert!(
            (mean - 2.0 / 3.0).abs() < 0.01,
            "mean cosine {mean} should be ~0.667"
        );
    }

    #[test]
    fn cosine_pdf_matches_the_definition() {
        // pdf(θ) = cosθ/π; straight overhead → 1/π, grazing → 0.
        assert!((cosine_hemisphere_pdf(1.0) - std::f32::consts::FRAC_1_PI).abs() < 1e-6);
        assert_eq!(cosine_hemisphere_pdf(0.0), 0.0);
        assert_eq!(cosine_hemisphere_pdf(-0.5), 0.0, "below horizon → 0");
    }

    #[test]
    fn uniform_disk_samples_stay_inside_the_disk() {
        let mut rng = Rng::new(55, 3);
        for _ in 0..5000 {
            let (x, y) = uniform_disk(rng.next_f32(), rng.next_f32());
            assert!(x * x + y * y <= 1.0 + 1e-5, "disk sample outside r=1");
        }
    }
}
