//! The render target — a high-dynamic-range float framebuffer and the
//! tone-mapping step that turns it into a displayable low-dynamic-range
//! image.
//!
//! A path tracer accumulates radiance, which is unbounded (a light
//! source is far brighter than 1.0). The renderer therefore writes into
//! an [`HdrFramebuffer`] of raw `f32` linear-RGB radiance, and only at
//! the very end maps that to an 8-bit [`LdrImage`] with a tone curve +
//! sRGB gamma — exactly the HDR-then-display split a real renderer
//! uses.

use crate::math::Vec3;

/// Maximum pixel count [`HdrFramebuffer::try_new`] will allocate.
/// 8K × 8K = 67,108,864 pixels — generous for any real path-trace
/// target and small enough that the accumulator (3 × f32 per pixel
/// = 12 bytes) stays under a gigabyte. The cap guards against the
/// `width * height` overflow class: pre-fix `HdrFramebuffer::new(65536, 65536)`
/// would silently wrap the multiplication, then allocate the wrapped
/// value's worth of `Vec3` — typically a few MB — but `add_sample`
/// indexed by the real `width * height` would write WAY out of bounds.
pub const MAX_FRAMEBUFFER_PIXELS: usize = 67_108_864;

/// Errors raised by [`HdrFramebuffer::try_new`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FramebufferError {
    /// `width * height` either overflowed `usize` or exceeded
    /// [`MAX_FRAMEBUFFER_PIXELS`].
    TooLarge {
        /// Image width that was requested.
        width: u32,
        /// Image height that was requested.
        height: u32,
    },
}

impl std::fmt::Display for FramebufferError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FramebufferError::TooLarge { width, height } => write!(
                f,
                "framebuffer dimensions {width}×{height} overflow or exceed the \
                 {MAX_FRAMEBUFFER_PIXELS}-pixel cap"
            ),
        }
    }
}

impl std::error::Error for FramebufferError {}

/// A high-dynamic-range render target — `width · height` linear-RGB
/// radiance samples, accumulated as the renderer draws samples.
#[derive(Clone, Debug)]
pub struct HdrFramebuffer {
    /// Image width in pixels.
    pub width: u32,
    /// Image height in pixels.
    pub height: u32,
    /// Accumulated radiance, row-major, `width · height` entries. Each
    /// entry is the **sum** of every sample drawn for that pixel;
    /// divide by [`Self::sample_count`] for the mean.
    pub accum: Vec<Vec3>,
    /// Number of samples accumulated into every pixel so far.
    pub sample_count: u32,
}

impl HdrFramebuffer {
    /// A black framebuffer of the given size, zero samples accumulated.
    ///
    /// # Panics
    ///
    /// Panics when `width * height` overflows `usize` or exceeds
    /// [`MAX_FRAMEBUFFER_PIXELS`]. Round-6 fix made the panicking
    /// path explicit (pre-fix it silently wrapped the multiplication
    /// and corrupted the index space). Callers that want recoverable
    /// failure should use [`HdrFramebuffer::try_new`].
    pub fn new(width: u32, height: u32) -> HdrFramebuffer {
        Self::try_new(width, height).expect("HdrFramebuffer::new: dimensions overflow")
    }

    /// Fallible constructor — returns
    /// [`FramebufferError::TooLarge`] instead of panicking when
    /// `width * height` overflows or exceeds the pixel cap.
    pub fn try_new(width: u32, height: u32) -> Result<HdrFramebuffer, FramebufferError> {
        let n = (width as usize)
            .checked_mul(height as usize)
            .ok_or(FramebufferError::TooLarge { width, height })?;
        if n > MAX_FRAMEBUFFER_PIXELS {
            return Err(FramebufferError::TooLarge { width, height });
        }
        Ok(HdrFramebuffer {
            width,
            height,
            accum: vec![Vec3::ZERO; n],
            sample_count: 0,
        })
    }

    /// Add one sample's radiance into pixel `(x, y)`.
    ///
    /// A non-finite sample (a stray NaN / inf — possible from a
    /// degenerate BRDF evaluation) is dropped so one bad path cannot
    /// poison the whole pixel. Out-of-range coordinates are ignored.
    #[inline]
    pub fn add_sample(&mut self, x: u32, y: u32, radiance: Vec3) {
        if x >= self.width || y >= self.height {
            return;
        }
        if !radiance.is_finite() {
            return;
        }
        let i = (y as usize) * (self.width as usize) + (x as usize);
        self.accum[i] = self.accum[i].add(radiance);
    }

    /// Record that one more sample-per-pixel has been completed.
    ///
    /// The renderer calls this once after a full image pass so
    /// [`Self::mean`] divides by the right count.
    pub fn finish_sample(&mut self) {
        self.sample_count += 1;
    }

    /// The mean radiance of pixel `(x, y)` — the accumulated sum
    /// divided by the sample count. Returns black before any sample
    /// has been finished.
    #[inline]
    pub fn mean(&self, x: u32, y: u32) -> Vec3 {
        if x >= self.width || y >= self.height || self.sample_count == 0 {
            return Vec3::ZERO;
        }
        let i = (y as usize) * (self.width as usize) + (x as usize);
        self.accum[i].scale(1.0 / self.sample_count as f32)
    }

    /// Tone-map the HDR buffer to a displayable 8-bit [`LdrImage`].
    ///
    /// See [`tonemap_pixel`] for the curve.
    pub fn to_ldr(&self, exposure: f32) -> LdrImage {
        let mut pixels = Vec::with_capacity(self.accum.len() * 3);
        for y in 0..self.height {
            for x in 0..self.width {
                let rgb = tonemap_pixel(self.mean(x, y), exposure);
                pixels.push(to_u8(rgb[0]));
                pixels.push(to_u8(rgb[1]));
                pixels.push(to_u8(rgb[2]));
            }
        }
        LdrImage {
            width: self.width,
            height: self.height,
            pixels,
        }
    }
}

/// A tone-mapped, gamma-encoded 8-bit image — the displayable result.
#[derive(Clone, Debug, PartialEq)]
pub struct LdrImage {
    /// Image width in pixels.
    pub width: u32,
    /// Image height in pixels.
    pub height: u32,
    /// Row-major RGB8 pixels, `width · height · 3` bytes.
    pub pixels: Vec<u8>,
}

impl LdrImage {
    /// The RGB triple at pixel `(x, y)`. Out-of-range coordinates
    /// return black.
    pub fn pixel(&self, x: u32, y: u32) -> [u8; 3] {
        if x >= self.width || y >= self.height {
            return [0, 0, 0];
        }
        let i = ((y as usize) * (self.width as usize) + (x as usize)) * 3;
        [self.pixels[i], self.pixels[i + 1], self.pixels[i + 2]]
    }
}

/// Tone-map one linear-RGB radiance value to a `[0, 1]` display triple.
///
/// The pipeline is:
///
/// 1. **Exposure** — multiply the radiance by `exposure` (a manual
///    "shutter / ISO" knob).
/// 2. **Tone curve** — the ACES filmic approximation (Narkowicz 2015):
///    a rational curve that compresses highlights gracefully and keeps
///    the shadows from crushing, the de-facto game-engine tone map.
/// 3. **sRGB gamma** — the IEC 61966-2-1 transfer function, so the
///    `[0, 1]` linear value is encoded for a standard display.
#[inline]
pub fn tonemap_pixel(radiance: Vec3, exposure: f32) -> [f32; 3] {
    let exposed = radiance.scale(exposure);
    let mut out = [0.0f32; 3];
    for (k, c) in [exposed.x, exposed.y, exposed.z].into_iter().enumerate() {
        let mapped = aces_filmic(c.max(0.0));
        out[k] = linear_to_srgb(mapped.clamp(0.0, 1.0));
    }
    out
}

/// The ACES filmic tone curve (Krzysztof Narkowicz's fitted
/// approximation): `(x·(a·x + b)) / (x·(c·x + d) + e)`.
///
/// Maps an unbounded non-negative linear value into `[0, 1)`, rolling
/// highlights off smoothly.
#[inline]
fn aces_filmic(x: f32) -> f32 {
    const A: f32 = 2.51;
    const B: f32 = 0.03;
    const C: f32 = 2.43;
    const D: f32 = 0.59;
    const E: f32 = 0.14;
    ((x * (A * x + B)) / (x * (C * x + D) + E)).clamp(0.0, 1.0)
}

/// Linear → sRGB gamma encode for one channel.
#[inline]
fn linear_to_srgb(c: f32) -> f32 {
    if c <= 0.003_130_8 {
        12.92 * c
    } else {
        1.055 * c.powf(1.0 / 2.4) - 0.055
    }
}

/// Quantise a `[0, 1]` channel to an 8-bit value, rounding to nearest.
#[inline]
fn to_u8(c: f32) -> u8 {
    (c.clamp(0.0, 1.0) * 255.0 + 0.5) as u8
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::vec3;

    #[test]
    fn fresh_framebuffer_is_black() {
        let fb = HdrFramebuffer::new(8, 8);
        assert_eq!(fb.mean(3, 3), Vec3::ZERO);
        assert_eq!(fb.sample_count, 0);
    }

    #[test]
    fn accumulating_samples_averages_correctly() {
        let mut fb = HdrFramebuffer::new(2, 2);
        // Two samples into pixel (0, 0): radiance 2 and radiance 4.
        fb.add_sample(0, 0, vec3(2.0, 2.0, 2.0));
        fb.finish_sample();
        fb.add_sample(0, 0, vec3(4.0, 4.0, 4.0));
        fb.finish_sample();
        // Mean = (2 + 4) / 2 = 3.
        let m = fb.mean(0, 0);
        assert!((m.x - 3.0).abs() < 1e-6, "mean should be 3, got {}", m.x);
    }

    #[test]
    fn non_finite_samples_are_dropped() {
        let mut fb = HdrFramebuffer::new(1, 1);
        fb.add_sample(0, 0, vec3(1.0, 1.0, 1.0));
        // A NaN sample must not poison the accumulator.
        fb.add_sample(0, 0, vec3(f32::NAN, 0.0, 0.0));
        fb.finish_sample();
        let m = fb.mean(0, 0);
        assert!(m.is_finite(), "NaN sample leaked into the buffer");
        assert!((m.x - 1.0).abs() < 1e-6);
    }

    #[test]
    fn out_of_range_samples_are_ignored() {
        let mut fb = HdrFramebuffer::new(4, 4);
        // Should not panic.
        fb.add_sample(99, 99, vec3(1.0, 1.0, 1.0));
        fb.finish_sample();
        assert_eq!(fb.mean(0, 0), Vec3::ZERO);
    }

    #[test]
    fn tonemap_is_monotone_and_bounded() {
        // Brighter input → brighter (or equal) output, always in [0,1].
        let mut prev = -1.0f32;
        for i in 0..50 {
            let r = i as f32 * 0.5;
            let out = tonemap_pixel(vec3(r, r, r), 1.0);
            assert!((0.0..=1.0).contains(&out[0]), "tonemap out of range");
            assert!(out[0] >= prev - 1e-4, "tonemap must be monotone");
            prev = out[0];
        }
    }

    #[test]
    fn tonemap_black_is_black_and_huge_is_near_white() {
        let black = tonemap_pixel(Vec3::ZERO, 1.0);
        assert!(black[0] < 1e-3, "black radiance should map to black");
        let bright = tonemap_pixel(vec3(1000.0, 1000.0, 1000.0), 1.0);
        assert!(bright[0] > 0.95, "a huge radiance should map near white");
    }

    #[test]
    fn exposure_scales_brightness() {
        let r = vec3(0.5, 0.5, 0.5);
        let dim = tonemap_pixel(r, 0.5);
        let bright = tonemap_pixel(r, 2.0);
        assert!(bright[0] > dim[0], "more exposure → brighter pixel");
    }

    #[test]
    fn to_ldr_produces_the_right_buffer_size() {
        let mut fb = HdrFramebuffer::new(5, 3);
        fb.finish_sample();
        let img = fb.to_ldr(1.0);
        assert_eq!(img.width, 5);
        assert_eq!(img.height, 3);
        assert_eq!(img.pixels.len(), 5 * 3 * 3, "RGB8 = w·h·3 bytes");
    }

    #[test]
    fn ldr_pixel_accessor_round_trips() {
        let mut fb = HdrFramebuffer::new(2, 2);
        // A bright white pixel at (1, 1).
        fb.add_sample(1, 1, vec3(50.0, 50.0, 50.0));
        fb.finish_sample();
        let img = fb.to_ldr(1.0);
        let white = img.pixel(1, 1);
        assert!(white[0] > 200, "bright pixel should be near-white");
        let black = img.pixel(0, 0);
        assert!(black[0] < 20, "untouched pixel should be near-black");
    }

    #[test]
    fn try_new_rejects_dimensions_past_max_pixels() {
        // Round-6 RED→GREEN: pre-fix `HdrFramebuffer::new(65536, 65536)`
        // silently wrapped the `width * height` multiplication and
        // produced a tiny `accum` vector that `add_sample` then
        // indexed wildly out of bounds. The cap rejects the
        // configuration up front.
        let err = HdrFramebuffer::try_new(65536, 65536).unwrap_err();
        assert!(matches!(err, FramebufferError::TooLarge { .. }));
        // The cap edge: at exactly the cap (8192×8192 = 67_108_864),
        // try_new succeeds (we have ~768 MiB of room, which is past
        // the cap-of-the-product check by design).
        assert!(
            HdrFramebuffer::try_new(8192, 8192).is_ok(),
            "8192² is exactly at MAX_FRAMEBUFFER_PIXELS"
        );
        // 8193² exceeds the cap.
        assert!(HdrFramebuffer::try_new(8193, 8193).is_err());
    }

    #[test]
    fn try_new_rejects_overflow_via_u32_max() {
        // `u32::MAX as usize * u32::MAX as usize` overflows the
        // usize space on every supported target (LP64 or LLP64).
        let err = HdrFramebuffer::try_new(u32::MAX, u32::MAX).unwrap_err();
        assert!(matches!(err, FramebufferError::TooLarge { .. }));
    }
}
