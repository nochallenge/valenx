//! ORB-style oriented binary descriptor.
//!
//! For each keypoint this module computes:
//!
//! 1. **Orientation** by the *intensity centroid* method (Rosten/Rublee):
//!    over a circular patch of radius [`PATCH_RADIUS`] around the corner,
//!    the first-order moments `m10 = Σ x·I` and `m01 = Σ y·I` define a
//!    centroid direction `θ = atan2(m01, m10)`. The corner "points" from
//!    its centre toward the brighter side; this gives a repeatable angle
//!    under in-plane rotation.
//!
//! 2. A **256-bit steered-BRIEF** descriptor: 256 fixed pairs of sample
//!    offsets `(a_i, b_i)` are rotated by `θ` and compared,
//!    `bit_i = I(a_i^θ) > I(b_i^θ)`. Steering the pattern by the patch
//!    orientation is the "oriented BRIEF" (`steered-BRIEF`) construction
//!    at the heart of ORB.
//!
//! ## Honesty note — this is *steered* BRIEF, not the full learned rBRIEF
//!
//! ORB's full descriptor (`rBRIEF`) adds a second step on top of steered
//! BRIEF: a *greedy learning* pass that selects, from a large pool of
//! candidate tests, the 256 binary tests that are individually
//! high-variance and mutually *decorrelated* (trained on a patch corpus).
//! That learning step recovers the variance that naive steering loses and
//! is what makes rBRIEF maximally discriminative.
//!
//! **This crate ships oriented-FAST + steered-BRIEF, with a deterministic
//! Gaussian-sampled test pattern — it does NOT perform the learned test
//! selection / decorrelation.** The pattern is generated once from a fixed
//! seed (a tiny in-crate SplitMix64 PRNG, so no `rand` dependency and full
//! cross-run/cross-machine reproducibility), giving an isotropic Gaussian
//! cloud of pairs within the patch. It is a faithful steered-BRIEF, and a
//! reasonable ORB approximation, but it is not bit-compatible with
//! OpenCV's learned `ORB_pattern` and is somewhat less discriminative than
//! true rBRIEF. Stage-2 matching (Hamming distance) works the same way
//! regardless.

use crate::image::GrayImage;
use crate::keypoint::Keypoint;

/// Radius, in pixels, of the circular patch used for both the
/// intensity-centroid orientation and the BRIEF sampling cloud.
pub const PATCH_RADIUS: i32 = 15;

/// Number of binary tests, i.e. descriptor length in bits.
pub const DESCRIPTOR_BITS: usize = 256;

/// Descriptor length in bytes (`DESCRIPTOR_BITS / 8`).
pub const DESCRIPTOR_BYTES: usize = DESCRIPTOR_BITS / 8;

/// One BRIEF test: two integer sample offsets relative to the keypoint,
/// before steering.
type TestPair = ((i32, i32), (i32, i32));

/// Deterministic 256-pair Gaussian-sampled BRIEF pattern, built lazily on
/// first use. Generated from a fixed PRNG seed so the descriptor is
/// identical on every run and platform.
fn brief_pattern() -> &'static [TestPair; DESCRIPTOR_BITS] {
    use std::sync::OnceLock;
    static PATTERN: OnceLock<[TestPair; DESCRIPTOR_BITS]> = OnceLock::new();
    PATTERN.get_or_init(build_pattern)
}

/// Build the fixed BRIEF test pattern: 256 pairs of points drawn from an
/// isotropic Gaussian (σ = patch/5, the ORB setting) centred on the
/// keypoint and clamped to the patch disc. Uses a small SplitMix64 PRNG
/// so the result is fully deterministic with no external dependency.
fn build_pattern() -> [TestPair; DESCRIPTOR_BITS] {
    let mut rng = SplitMix64::new(0x5EED_C0DE_1234_ABCD);
    // ORB samples from a Gaussian with σ = patchSize/5 ≈ (2R+1)/5.
    let sigma = (2.0 * PATCH_RADIUS as f64 + 1.0) / 5.0;
    let limit = PATCH_RADIUS; // keep samples inside the patch disc radius

    let mut pattern: [TestPair; DESCRIPTOR_BITS] = [((0, 0), (0, 0)); DESCRIPTOR_BITS];
    for slot in pattern.iter_mut() {
        let a = gaussian_point(&mut rng, sigma, limit);
        let b = gaussian_point(&mut rng, sigma, limit);
        *slot = (a, b);
    }
    pattern
}

/// Draw one Gaussian-distributed integer offset within the patch disc.
/// Rejection-resamples until the point lies inside the radius-`limit`
/// circle so steering can never push a tap outside the intended patch.
fn gaussian_point(rng: &mut SplitMix64, sigma: f64, limit: i32) -> (i32, i32) {
    let lf = limit as f64;
    loop {
        let (g0, g1) = rng.next_gaussian_pair();
        let x = (g0 * sigma).round();
        let y = (g1 * sigma).round();
        if x * x + y * y <= lf * lf {
            return (x as i32, y as i32);
        }
        // If the draw landed outside the disc, shrink and retry on the next
        // loop iteration with a fresh sample.
    }
}

/// Compute the intensity-centroid orientation (radians, `(-π, π]`) of the
/// patch centred at integer `(cx, cy)`.
///
/// Sample taps that would fall outside the image are clamped to the
/// nearest edge pixel, so corners near (but not on) the border still get a
/// stable angle without any out-of-bounds access.
fn centroid_orientation(img: &GrayImage, cx: i32, cy: i32) -> f32 {
    let mut m01: i64 = 0;
    let mut m10: i64 = 0;
    let r = PATCH_RADIUS;
    let r2 = r * r;
    for dy in -r..=r {
        for dx in -r..=r {
            if dx * dx + dy * dy > r2 {
                continue;
            }
            let v = i64::from(sample_clamped(img, cx + dx, cy + dy));
            m10 += i64::from(dx) * v;
            m01 += i64::from(dy) * v;
        }
    }
    (m01 as f32).atan2(m10 as f32)
}

/// Compute the steered-BRIEF descriptor for the keypoint at integer
/// `(cx, cy)` with orientation `theta` (radians).
fn steered_brief(img: &GrayImage, cx: i32, cy: i32, theta: f32) -> [u8; DESCRIPTOR_BYTES] {
    let (s, c) = theta.sin_cos();
    let mut desc = [0u8; DESCRIPTOR_BYTES];

    for (i, &((ax, ay), (bx, by))) in brief_pattern().iter().enumerate() {
        // Rotate each test point by theta, then sample (edge-clamped).
        let (rax, ray) = rotate(ax, ay, s, c);
        let (rbx, rby) = rotate(bx, by, s, c);
        let ia = sample_clamped(img, cx + rax, cy + ray);
        let ib = sample_clamped(img, cx + rbx, cy + rby);
        if ia > ib {
            desc[i >> 3] |= 1u8 << (i & 7);
        }
    }
    desc
}

/// Rotate an integer offset `(x, y)` by the angle whose sine/cosine are
/// `(s, c)`, returning the rounded integer offset.
#[inline]
fn rotate(x: i32, y: i32, s: f32, c: f32) -> (i32, i32) {
    let xf = x as f32;
    let yf = y as f32;
    let rx = c * xf - s * yf;
    let ry = s * xf + c * yf;
    (rx.round() as i32, ry.round() as i32)
}

/// Sample `(x, y)`, clamping out-of-range coordinates to the nearest edge
/// pixel. Guarantees no out-of-bounds access for any keypoint, including
/// those a few pixels from the border whose steered taps reach outside.
#[inline]
fn sample_clamped(img: &GrayImage, x: i32, y: i32) -> u8 {
    let cx = x.clamp(0, img.width as i32 - 1) as usize;
    let cy = y.clamp(0, img.height as i32 - 1) as usize;
    img.at(cx, cy)
}

/// Assign an intensity-centroid orientation and a 256-bit steered-BRIEF
/// descriptor to a single keypoint, returning the (oriented) keypoint and
/// its descriptor bytes.
///
/// The input keypoint's position is rounded to the nearest integer pixel
/// for sampling; the returned keypoint carries the computed orientation in
/// [`Keypoint::angle`].
#[must_use]
pub fn describe_keypoint(img: &GrayImage, kp: &Keypoint) -> (Keypoint, [u8; DESCRIPTOR_BYTES]) {
    let cx = kp.x.round() as i32;
    let cy = kp.y.round() as i32;
    let theta = centroid_orientation(img, cx, cy);
    let desc = steered_brief(img, cx, cy, theta);
    let mut oriented = *kp;
    oriented.angle = theta;
    (oriented, desc)
}

/// Hamming distance between two descriptors: the number of differing bits.
/// The natural similarity metric for binary descriptors (Stage-2 matching
/// uses it). Range `0..=256`.
#[must_use]
pub fn hamming_distance(a: &[u8; DESCRIPTOR_BYTES], b: &[u8; DESCRIPTOR_BYTES]) -> u32 {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x ^ y).count_ones())
        .sum()
}

/// Minimal SplitMix64 PRNG — used only to generate the fixed BRIEF test
/// pattern deterministically, so the crate needs no `rand` dependency.
///
/// SplitMix64 is a well-known, high-quality finalizer-based generator
/// (Steele/Lea/Flood, the seeding RNG of the xoshiro family). It is *not*
/// used for any security purpose.
struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    /// Seed the generator.
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    /// Next raw 64-bit output.
    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform `f64` in `[0, 1)` from the top 53 bits.
    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// A pair of independent standard-normal samples via the Box–Muller
    /// transform.
    fn next_gaussian_pair(&mut self) -> (f64, f64) {
        // Guard the log against u1 == 0.
        let u1 = (self.next_f64()).max(f64::MIN_POSITIVE);
        let u2 = self.next_f64();
        let r = (-2.0 * u1.ln()).sqrt();
        let theta = 2.0 * std::f64::consts::PI * u2;
        (r * theta.cos(), r * theta.sin())
    }
}
