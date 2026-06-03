//! Edge-avoiding **à-trous wavelet denoiser** (Dammertz et al. 2010) —
//! a classical, **non-ML** Monte-Carlo image denoiser.
//!
//! # The problem
//!
//! A path-traced image at a low sample count is *noisy* — every pixel
//! is a Monte-Carlo estimate whose variance falls only as `1/spp`, so
//! a clean image costs hundreds of samples. A **denoiser** removes that
//! noise as a post-process so a far cheaper render looks finished.
//!
//! # The algorithm — "Edge-Avoiding À-Trous Wavelet Transform for
//! Fast Global Illumination Filtering"
//!
//! This is the Dammertz/Sewtz/Hanika/Lensch (HPG 2010) filter. A
//! plain blur would remove the noise *and* the image's real edges. The
//! à-trous filter keeps the edges by guiding the blur with **per-pixel
//! feature buffers** the renderer produces for free alongside the
//! colour:
//!
//! - **albedo** — the surface base colour at the pixel (no lighting);
//! - **normal** — the surface shading normal;
//! - **depth** — the camera-space hit distance.
//!
//! These feature ("G-buffer") channels are essentially noise-free —
//! they come from the *first* hit, not the Monte-Carlo light
//! transport. The filter blurs the noisy colour with a 5×5 B-spline
//! kernel, but **weights each tap** by how similar its albedo, normal
//! and depth are to the centre pixel. Across a real edge (a different
//! surface, a silhouette) the feature weights collapse toward zero, so
//! the blur does not cross it — the edge is preserved while the
//! flat-shaded interior is smoothed.
//!
//! "À-trous" ("with holes") means the kernel is applied at growing
//! **dilations**: iteration `i` spaces its 5×5 taps `2ⁱ` pixels apart.
//! A handful of iterations therefore filters a wide neighbourhood at
//! `O(n)` cost per iteration instead of the `O(n·r²)` of one big
//! kernel — the defining trick of the wavelet transform.
//!
//! # Honest scope
//!
//! A real, faithful implementation of the 2010 paper. It is a
//! **single-frame, spatial** denoiser: it does not do the temporal
//! reprojection / accumulation a real-time renderer (SVGF) layers on
//! top, nor does it consume a per-pixel variance estimate to drive the
//! colour weight adaptively (the paper's optional refinement). Those
//! are documented, additive extensions. What ships is the genuine
//! edge-avoiding à-trous filter: a noisy constant image denoises to
//! the constant, and an albedo / normal edge survives the filter — the
//! two properties the tests check.

use crate::framebuffer::HdrFramebuffer;
use crate::math::Vec3;

/// The per-pixel **guide buffers** an à-trous denoise needs alongside
/// the noisy colour — the (near-)noise-free feature channels.
///
/// All three are `width · height`, row-major, matching the colour
/// [`HdrFramebuffer`]'s layout. A renderer fills them from the
/// **primary hit** (the first surface each camera ray strikes), so
/// they carry no Monte-Carlo noise and make ideal edge guides.
#[derive(Clone, Debug)]
pub struct GuideBuffers {
    /// Image width in pixels.
    pub width: u32,
    /// Image height in pixels.
    pub height: u32,
    /// Surface albedo (base colour, unlit) at each pixel's primary
    /// hit. An albedo discontinuity marks a material boundary the
    /// filter must not blur across.
    pub albedo: Vec<Vec3>,
    /// Shading normal at each pixel's primary hit, world space. A
    /// normal discontinuity marks a geometric edge / silhouette.
    pub normal: Vec<Vec3>,
    /// Camera-space depth (hit distance) at each pixel. A depth jump
    /// marks an occlusion boundary.
    pub depth: Vec<f32>,
}

impl GuideBuffers {
    /// An all-zero guide set of the given size — a renderer fills it
    /// in as it shoots primary rays.
    pub fn new(width: u32, height: u32) -> GuideBuffers {
        let n = (width as usize) * (height as usize);
        GuideBuffers {
            width,
            height,
            albedo: vec![Vec3::ZERO; n],
            normal: vec![Vec3::ZERO; n],
            depth: vec![0.0; n],
        }
    }

    /// Record the primary-hit features for pixel `(x, y)`.
    ///
    /// A renderer calls this once per pixel while generating primary
    /// rays. Out-of-range coordinates are ignored.
    pub fn set(&mut self, x: u32, y: u32, albedo: Vec3, normal: Vec3, depth: f32) {
        if x >= self.width || y >= self.height {
            return;
        }
        let i = (y as usize) * (self.width as usize) + (x as usize);
        self.albedo[i] = albedo;
        self.normal[i] = normal;
        self.depth[i] = depth;
    }

    #[inline]
    fn index(&self, x: u32, y: u32) -> usize {
        (y as usize) * (self.width as usize) + (x as usize)
    }
}

/// Tunable weights of the edge-avoiding à-trous filter — the σ
/// (sigma) parameters of the paper.
///
/// Each σ controls how sharply the filter rejects a tap whose feature
/// differs from the centre. A **smaller** σ → a more selective filter
/// (stops at fainter edges, removes less noise); a **larger** σ → a
/// more aggressive blur (smoother, but feature edges blur sooner).
#[derive(Clone, Copy, Debug)]
pub struct AtrousParams {
    /// Number of à-trous iterations. Each doubles the kernel dilation,
    /// so `n` iterations filter a `~2ⁿ`-pixel neighbourhood. 5 is the
    /// paper's default and plenty for a typical render.
    pub iterations: u32,
    /// Colour-similarity σ. Rejects a tap whose *colour* is far from
    /// the centre — the term that actually removes the noise; set
    /// generous (a large value) so noise is smoothed, since the
    /// feature terms guard the real edges.
    pub sigma_color: f32,
    /// Normal-similarity σ. A small value keeps the filter from
    /// blurring across a geometric edge / silhouette.
    pub sigma_normal: f32,
    /// Depth-similarity σ. Keeps the filter from blurring across an
    /// occlusion boundary (a near object in front of a far one).
    pub sigma_depth: f32,
    /// Albedo-similarity σ. Keeps the filter from blurring across a
    /// material boundary (two surfaces of different base colour).
    pub sigma_albedo: f32,
}

impl Default for AtrousParams {
    /// The paper's defaults, tuned for an HDR linear-radiance image:
    /// 5 iterations, a loose colour σ, tight feature σ's.
    fn default() -> Self {
        AtrousParams {
            iterations: 5,
            sigma_color: 8.0,
            sigma_normal: 0.2,
            sigma_depth: 1.0,
            sigma_albedo: 0.25,
        }
    }
}

/// The 5-tap **B-spline** row `(1, 4, 6, 4, 1)/16`; the 2-D à-trous
/// kernel is its outer product. This is the wavelet smoothing filter
/// the paper uses.
const KERNEL_1D: [f32; 5] = [1.0 / 16.0, 4.0 / 16.0, 6.0 / 16.0, 4.0 / 16.0, 1.0 / 16.0];

/// Denoise a noisy [`HdrFramebuffer`] with the edge-avoiding à-trous
/// wavelet filter, guided by the [`GuideBuffers`].
///
/// Returns a fresh framebuffer holding the denoised radiance (the
/// input is not modified). The guide buffers must match the colour
/// buffer's dimensions; on a mismatch the input is returned
/// unchanged.
///
/// # Method
///
/// For each of `params.iterations` iterations the current image is
/// convolved with the 5×5 B-spline kernel at dilation `step = 2ⁱ`.
/// Every tap `q` around centre `p` is weighted by
///
/// ```text
///   w(p,q) = h(q) · w_color · w_normal · w_depth · w_albedo
/// ```
///
/// where `h(q)` is the B-spline kernel coefficient and each `w_*` is
/// `exp(−‖feature_p − feature_q‖² / σ²)` — a Gaussian on the feature
/// difference. The output is the weight-normalised tap sum. Because a
/// real edge makes one or more feature weights vanish, the blur stops
/// at the edge while smoothing the noisy flat regions.
pub fn denoise_atrous(
    color: &HdrFramebuffer,
    guides: &GuideBuffers,
    params: &AtrousParams,
) -> HdrFramebuffer {
    let w = color.width;
    let h = color.height;
    if guides.width != w || guides.height != h {
        // Dimension mismatch — cannot denoise; pass the input through.
        return color.clone();
    }
    if w == 0 || h == 0 {
        return color.clone();
    }

    // Work on the per-pixel *mean* radiance (the displayed image),
    // not the running accumulator sum.
    let n = (w as usize) * (h as usize);
    let mut current: Vec<Vec3> = Vec::with_capacity(n);
    for y in 0..h {
        for x in 0..w {
            current.push(color.mean(x, y));
        }
    }

    let iterations = params.iterations.max(1);
    let mut scratch: Vec<Vec3> = vec![Vec3::ZERO; n];

    for it in 0..iterations {
        // À-trous dilation: taps spaced 2^it pixels apart.
        let step = 1i32 << it;
        atrous_iteration(&current, &mut scratch, guides, params, w, h, step);
        std::mem::swap(&mut current, &mut scratch);
    }

    // Pack the denoised means back into a framebuffer with one
    // "sample" so `mean(x, y)` returns the denoised value directly.
    let mut out = HdrFramebuffer::new(w, h);
    out.accum.copy_from_slice(&current);
    out.sample_count = 1;
    out
}

/// One à-trous iteration: convolve `src` into `dst` with the 5×5
/// B-spline kernel at dilation `step`, edge-weighted by the guides.
fn atrous_iteration(
    src: &[Vec3],
    dst: &mut [Vec3],
    guides: &GuideBuffers,
    params: &AtrousParams,
    w: u32,
    h: u32,
    step: i32,
) {
    // Guard the σ's against a zero divisor.
    let inv_color = 1.0 / (params.sigma_color.max(1e-4)).powi(2);
    let inv_normal = 1.0 / (params.sigma_normal.max(1e-4)).powi(2);
    let inv_depth = 1.0 / (params.sigma_depth.max(1e-4)).powi(2);
    let inv_albedo = 1.0 / (params.sigma_albedo.max(1e-4)).powi(2);

    for y in 0..h {
        for x in 0..w {
            let p = guides.index(x, y);
            let c_p = src[p];
            let n_p = guides.normal[p];
            let d_p = guides.depth[p];
            let a_p = guides.albedo[p];

            let mut sum = Vec3::ZERO;
            let mut weight_sum = 0.0f32;

            // 5×5 B-spline neighbourhood at the current dilation.
            for (ky, &kh_y) in KERNEL_1D.iter().enumerate() {
                let oy = (ky as i32 - 2) * step;
                let sy = y as i32 + oy;
                if sy < 0 || sy >= h as i32 {
                    continue;
                }
                for (kx, &kh_x) in KERNEL_1D.iter().enumerate() {
                    let ox = (kx as i32 - 2) * step;
                    let sx = x as i32 + ox;
                    if sx < 0 || sx >= w as i32 {
                        continue;
                    }
                    let q = guides.index(sx as u32, sy as u32);
                    let kernel = kh_x * kh_y;

                    // --- edge-stopping feature weights ---
                    // Colour: removes the noise. A large difference is
                    // tolerated (σ_color is loose) so the filter still
                    // smooths; the feature terms below guard real edges.
                    let dc = src[q].sub(c_p);
                    let w_color = (-dc.length_sq() * inv_color).exp();

                    // Normal: a Gaussian on the normal difference —
                    // collapses across a geometric silhouette.
                    let dn = guides.normal[q].sub(n_p);
                    let w_normal = (-dn.length_sq() * inv_normal).exp();

                    // Depth: scaled by the dilation step so a fixed
                    // depth gradient is judged consistently across
                    // iterations (the paper's per-iteration depth term).
                    let dd = (guides.depth[q] - d_p) / step as f32;
                    let w_depth = (-(dd * dd) * inv_depth).exp();

                    // Albedo: collapses across a material boundary.
                    let da = guides.albedo[q].sub(a_p);
                    let w_albedo = (-da.length_sq() * inv_albedo).exp();

                    let weight = kernel * w_color * w_normal * w_depth * w_albedo;
                    sum = sum.add(src[q].scale(weight));
                    weight_sum += weight;
                }
            }

            // Normalise by the accumulated weight. A fully-isolated
            // pixel (every neighbour rejected) keeps its own value.
            dst[p] = if weight_sum > 1e-8 {
                sum.scale(1.0 / weight_sum)
            } else {
                c_p
            };
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::vec3;

    /// Build a framebuffer of the given per-pixel mean radiances.
    fn framebuffer_from(w: u32, h: u32, pixels: &[Vec3]) -> HdrFramebuffer {
        let mut fb = HdrFramebuffer::new(w, h);
        fb.accum.copy_from_slice(pixels);
        fb.sample_count = 1;
        fb
    }

    /// A flat guide set — every pixel the same surface. Edges therefore
    /// come only from the colour buffer.
    fn flat_guides(w: u32, h: u32) -> GuideBuffers {
        let mut g = GuideBuffers::new(w, h);
        for y in 0..h {
            for x in 0..w {
                g.set(x, y, vec3(0.5, 0.5, 0.5), vec3(0.0, 0.0, 1.0), 5.0);
            }
        }
        g
    }

    /// A noisy *constant* image must denoise back to (near) the
    /// constant — the headline correctness check: a denoiser may not
    /// shift the mean of a flat region, only remove its variance.
    #[test]
    fn noisy_constant_image_denoises_to_the_constant() {
        let w = 16;
        let h = 16;
        let mean = 0.5f32;
        // A constant image perturbed by deterministic pseudo-noise.
        let mut pixels = Vec::with_capacity((w * h) as usize);
        let mut seed = 0x1234_5678u32;
        for _ in 0..(w * h) {
            // Cheap LCG noise in [-0.3, 0.3].
            seed = seed.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            let r = (seed >> 8) as f32 / 16_777_216.0;
            let noise = (r - 0.5) * 0.6;
            pixels.push(vec3(mean + noise, mean + noise, mean + noise));
        }
        let fb = framebuffer_from(w, h, &pixels);
        let guides = flat_guides(w, h);
        let params = AtrousParams::default();
        let denoised = denoise_atrous(&fb, &guides, &params);

        // The denoised image should be far flatter — measure variance
        // before and after over the interior (avoid the borders).
        let variance = |buf: &HdrFramebuffer| -> f32 {
            let mut sum = 0.0f64;
            let mut sum_sq = 0.0f64;
            let mut count = 0u32;
            for y in 3..h - 3 {
                for x in 3..w - 3 {
                    let v = buf.mean(x, y).x as f64;
                    sum += v;
                    sum_sq += v * v;
                    count += 1;
                }
            }
            let m = sum / count as f64;
            ((sum_sq / count as f64) - m * m) as f32
        };
        let before = variance(&fb);
        let after = variance(&denoised);
        assert!(
            after < before * 0.25,
            "denoise should cut the variance: before {before}, after {after}"
        );
        // And it must not have shifted the mean of the flat field.
        let mut total = 0.0f64;
        let mut n = 0u32;
        for y in 3..h - 3 {
            for x in 3..w - 3 {
                total += denoised.mean(x, y).x as f64;
                n += 1;
            }
        }
        let denoised_mean = (total / n as f64) as f32;
        assert!(
            (denoised_mean - mean).abs() < 0.03,
            "denoised mean {denoised_mean} should stay at the constant {mean}"
        );
    }

    /// An **albedo edge** must survive the filter. The colour buffer
    /// has a sharp left/right step that coincides with an albedo step
    /// in the guide buffer; the denoiser must NOT blur across it.
    #[test]
    fn an_albedo_edge_is_preserved() {
        let w = 24;
        let h = 8;
        // Colour: left half dark (0.2), right half bright (0.9).
        let mut pixels = Vec::with_capacity((w * h) as usize);
        for _y in 0..h {
            for x in 0..w {
                let v = if x < w / 2 { 0.2 } else { 0.9 };
                pixels.push(vec3(v, v, v));
            }
        }
        let fb = framebuffer_from(w, h, &pixels);
        // Guides: albedo steps at the same column; normal & depth flat.
        let mut g = GuideBuffers::new(w, h);
        for y in 0..h {
            for x in 0..w {
                let a = if x < w / 2 {
                    vec3(0.1, 0.1, 0.1)
                } else {
                    vec3(0.9, 0.9, 0.9)
                };
                g.set(x, y, a, vec3(0.0, 0.0, 1.0), 5.0);
            }
        }
        let denoised = denoise_atrous(&fb, &g, &AtrousParams::default());

        // The pixels just left and right of the seam must keep their
        // distinct values — the step must not have washed out.
        let left = denoised.mean(w / 2 - 1, h / 2).x;
        let right = denoised.mean(w / 2, h / 2).x;
        assert!(
            (right - left) > 0.55,
            "albedo edge should be preserved: left {left}, right {right}"
        );
        // The dark side stays dark, the bright side stays bright.
        assert!(left < 0.35, "dark side bled bright: {left}");
        assert!(right > 0.75, "bright side bled dark: {right}");
    }

    /// A **normal edge** (a geometric silhouette) is likewise
    /// preserved — a colour step that coincides with a normal step in
    /// the guide buffer is not blurred across.
    #[test]
    fn a_normal_edge_is_preserved() {
        let w = 24;
        let h = 8;
        let mut pixels = Vec::with_capacity((w * h) as usize);
        for _y in 0..h {
            for x in 0..w {
                let v = if x < w / 2 { 0.15 } else { 0.85 };
                pixels.push(vec3(v, v, v));
            }
        }
        let fb = framebuffer_from(w, h, &pixels);
        // Guides: normal flips between two faces at the seam.
        let mut g = GuideBuffers::new(w, h);
        for y in 0..h {
            for x in 0..w {
                let nrm = if x < w / 2 {
                    vec3(0.0, 0.0, 1.0)
                } else {
                    vec3(1.0, 0.0, 0.0)
                };
                g.set(x, y, vec3(0.5, 0.5, 0.5), nrm, 5.0);
            }
        }
        let denoised = denoise_atrous(&fb, &g, &AtrousParams::default());
        let left = denoised.mean(w / 2 - 1, h / 2).x;
        let right = denoised.mean(w / 2, h / 2).x;
        assert!(
            (right - left) > 0.5,
            "normal edge should be preserved: left {left}, right {right}"
        );
    }

    /// Inside a flat-feature region the filter genuinely smooths: two
    /// noisy neighbours end up closer together after the denoise.
    #[test]
    fn flat_region_is_smoothed() {
        let w = 12;
        let h = 12;
        let mut pixels = vec![vec3(0.5, 0.5, 0.5); (w * h) as usize];
        // Two strongly-different adjacent interior pixels, same feature.
        pixels[(6 * w + 5) as usize] = vec3(0.9, 0.9, 0.9);
        pixels[(6 * w + 6) as usize] = vec3(0.1, 0.1, 0.1);
        let fb = framebuffer_from(w, h, &pixels);
        let guides = flat_guides(w, h);
        let denoised = denoise_atrous(&fb, &guides, &AtrousParams::default());
        let a = denoised.mean(5, 6).x;
        let b = denoised.mean(6, 6).x;
        // After smoothing they should be much closer than 0.8 apart.
        assert!(
            (a - b).abs() < 0.4,
            "flat-region neighbours should be smoothed together: {a} vs {b}"
        );
    }

    /// A dimension mismatch between colour and guides returns the
    /// input unchanged rather than panicking.
    #[test]
    fn dimension_mismatch_passes_the_image_through() {
        let fb = framebuffer_from(8, 8, &vec![vec3(0.3, 0.3, 0.3); 64]);
        let guides = GuideBuffers::new(4, 4); // wrong size
        let out = denoise_atrous(&fb, &guides, &AtrousParams::default());
        assert_eq!(out.width, 8);
        assert!((out.mean(2, 2).x - 0.3).abs() < 1e-6);
    }

    /// More iterations widen the effective filter — a single bright
    /// outlier spreads its energy further. We confirm the denoise of a
    /// lone spike is smoother with more iterations.
    #[test]
    fn more_iterations_widen_the_filter() {
        let w = 32;
        let h = 32;
        let mut pixels = vec![vec3(0.0, 0.0, 0.0); (w * h) as usize];
        pixels[(16 * w + 16) as usize] = vec3(10.0, 10.0, 10.0);
        let fb = framebuffer_from(w, h, &pixels);
        let guides = flat_guides(w, h);

        let one = denoise_atrous(
            &fb,
            &guides,
            &AtrousParams {
                iterations: 1,
                ..AtrousParams::default()
            },
        );
        let many = denoise_atrous(
            &fb,
            &guides,
            &AtrousParams {
                iterations: 5,
                ..AtrousParams::default()
            },
        );
        // A pixel well away from the spike picks up more energy with
        // more iterations (the filter reached further).
        let far_one = one.mean(16, 22).x;
        let far_many = many.mean(16, 22).x;
        assert!(
            far_many > far_one,
            "more iterations should spread the spike further: {far_one} → {far_many}"
        );
        // Energy is conserved-ish: neither denoise invents huge energy.
        let total = |buf: &HdrFramebuffer| -> f32 {
            let mut s = 0.0;
            for y in 0..h {
                for x in 0..w {
                    s += buf.mean(x, y).x;
                }
            }
            s
        };
        // The B-spline kernel is normalised, so the total stays close
        // to the input's single spike of 10.
        assert!(total(&many) < 20.0, "denoise should not amplify energy");
    }
}
