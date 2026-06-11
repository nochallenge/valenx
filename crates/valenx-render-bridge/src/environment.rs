//! HDR environment lighting + image-based lighting (IBL).
//!
//! ## What this is
//!
//! An environment map is a high-dynamic-range image wrapped around the
//! whole scene that acts as a light source — the standard way to limb
//! a model with realistic, photographically-captured illumination
//! ("studio HDR", "outdoor sky", and so on). This module ships:
//!
//! - [`EnvironmentMap`] — an equirectangular (latitude-longitude)
//!   floating-point RGB image plus an intensity multiplier and a
//!   yaw rotation.
//! - [`EnvironmentMap::from_radiance_hdr`] — a real decoder for the
//!   **Radiance `.hdr` / `.pic` RGBE** format (the de-facto HDR
//!   environment-map container). The format is small and fully
//!   specified, so no external image-codec dependency is pulled in —
//!   the render-bridge stays pure-data.
//! - [`EnvironmentMap::sample_direction`] — bilinear equirectangular
//!   lookup of the radiance arriving from a world direction (used for
//!   specular reflections / the background).
//! - [`EnvironmentMap::diffuse_irradiance`] — the **IBL diffuse
//!   convolution**: the cosine-weighted hemisphere integral of the
//!   map about a surface normal, i.e. the irradiance a Lambertian
//!   surface with that normal receives from the environment.
//! - [`IrradianceMap`] — a small precomputed irradiance map so a
//!   renderer can look diffuse lighting up per-pixel instead of
//!   re-integrating the hemisphere every shading call.
//!
//! ## Honest scope
//!
//! This is real IBL maths — the diffuse convolution and the RGBE
//! decode are exactly what a production renderer does. What is *not*
//! here:
//!
//! - **Specular prefiltering** (the split-sum roughness-mip chain) is
//!   not built — only the diffuse term is convolved. A renderer can
//!   still do sharp mirror reflections via
//!   [`EnvironmentMap::sample_direction`]; a roughness-aware
//!   prefilter is a bounded follow-up.
//! - **Importance sampling** — the convolution is a uniform
//!   hemisphere quadrature; for a very high-contrast map (a tiny
//!   bright sun) importance sampling would converge faster. The
//!   uniform quadrature is correct, just not the fastest.
//! - The render-bridge stays pure-data: this module computes the
//!   lighting maths and decodes the file; uploading the result to a
//!   GPU is the desktop shell's job.

use crate::error::RenderError;

/// Maximum width or height (in pixels) accepted from a Radiance `.hdr`
/// header before any pixel buffer is allocated.
///
/// The resolution line (`-Y h +X w`) is attacker-controlled: a ~35-byte
/// header such as `-Y 1000000000 +X 1000000000` passes the byte-length
/// cap that [`crate::scene::EnvironmentRef::load`] applies to the file
/// read, yet asks `from_radiance_hdr` to allocate ~10^18 pixels — an
/// instant out-of-memory abort — and the unchecked `width * height`
/// `usize` multiply wraps in release builds. 16 384 (16K per side) is
/// generous for any real environment map (16384² ≈ 268M pixels), and is
/// validated together with a `checked_mul` pixel ceiling *before* the
/// first allocation.
const MAX_HDR_DIM: usize = 16_384;

/// An equirectangular HDR environment map.
///
/// Pixels are stored row-major, `width · height` RGB triples of
/// linear-light floats. The equirectangular convention used
/// throughout this module: the image `u ∈ [0, 1)` maps to azimuth
/// `φ ∈ [0, 2π)` and `v ∈ [0, 1)` maps to polar angle `θ ∈ [0, π]`
/// (`v = 0` is straight up, `+Y`).
#[derive(Clone, Debug, PartialEq)]
pub struct EnvironmentMap {
    /// Image width in pixels.
    pub width: usize,
    /// Image height in pixels.
    pub height: usize,
    /// `width · height` linear-RGB float triples, row-major.
    pub pixels: Vec<[f32; 3]>,
    /// Scalar multiplier applied to every sample — exposure control.
    pub intensity: f32,
    /// Yaw rotation of the map about the world `+Y` axis, in radians.
    /// Lets the user spin the environment without re-decoding.
    pub yaw: f32,
}

impl EnvironmentMap {
    /// Build a map from a raw float-RGB buffer.
    ///
    /// # Errors
    ///
    /// [`RenderError::BadParameter`] if the dimensions are zero or the
    /// pixel count does not equal `width · height`.
    pub fn from_pixels(
        width: usize,
        height: usize,
        pixels: Vec<[f32; 3]>,
    ) -> Result<Self, RenderError> {
        if width == 0 || height == 0 {
            return Err(RenderError::BadParameter {
                name: "dimensions",
                reason: "environment map dimensions must be non-zero".into(),
            });
        }
        if pixels.len() != width * height {
            return Err(RenderError::BadParameter {
                name: "pixels",
                reason: format!(
                    "expected {} pixels for {width}×{height}, got {}",
                    width * height,
                    pixels.len()
                ),
            });
        }
        Ok(Self {
            width,
            height,
            pixels,
            intensity: 1.0,
            yaw: 0.0,
        })
    }

    /// A uniform (constant-colour) environment — a flat-lit "studio"
    /// with no detail. Useful as a default and for tests.
    pub fn uniform(color: [f32; 3]) -> Self {
        Self {
            width: 1,
            height: 1,
            pixels: vec![color],
            intensity: 1.0,
            yaw: 0.0,
        }
    }

    /// Decode a Radiance `.hdr` / `.pic` RGBE image.
    ///
    /// Supports the standard `32-bit_rle_rgbe` format: the ASCII
    /// header (terminated by a blank line), the resolution line
    /// (`-Y h +X w`), and both the new-style adaptive run-length
    /// encoding and the flat / old-style scanline encodings. Each
    /// RGBE pixel `(r, g, b, e)` decodes to linear RGB
    /// `c · 2^(e-128) / 256`.
    ///
    /// # Errors
    ///
    /// [`RenderError::BadParameter`] for a malformed magic line,
    /// missing resolution line, an unsupported (non-`-Y +X`)
    /// orientation, or a truncated pixel stream.
    pub fn from_radiance_hdr(bytes: &[u8]) -> Result<Self, RenderError> {
        let bad = |reason: &str| RenderError::BadParameter {
            name: "radiance_hdr",
            reason: reason.to_string(),
        };

        // --- header ---
        // First line must be a Radiance magic (`#?RADIANCE` or
        // `#?RGBE`).
        let mut cursor = 0usize;
        let first = read_line(bytes, &mut cursor).ok_or_else(|| bad("empty file"))?;
        if !first.starts_with("#?") {
            return Err(bad("missing Radiance `#?` magic line"));
        }
        // Header lines until a blank line.
        loop {
            let line = read_line(bytes, &mut cursor)
                .ok_or_else(|| bad("header not terminated by a blank line"))?;
            if line.is_empty() {
                break;
            }
            // FORMAT must be the RGBE variant if present.
            if let Some(fmt) = line.strip_prefix("FORMAT=") {
                let fmt = fmt.trim();
                if fmt != "32-bit_rle_rgbe" && fmt != "32-bit_rle_xyze" {
                    return Err(bad("only 32-bit RLE RGBE/XYZE HDR is supported"));
                }
            }
        }
        // Resolution line: `-Y <height> +X <width>`.
        let res = read_line(bytes, &mut cursor)
            .ok_or_else(|| bad("missing resolution line"))?;
        let parts: Vec<&str> = res.split_whitespace().collect();
        if parts.len() != 4 || parts[0] != "-Y" || parts[2] != "+X" {
            return Err(bad(
                "only the standard `-Y h +X w` HDR orientation is supported",
            ));
        }
        let height: usize = parts[1]
            .parse()
            .map_err(|_| bad("bad height in resolution line"))?;
        let width: usize = parts[3]
            .parse()
            .map_err(|_| bad("bad width in resolution line"))?;
        if width == 0 || height == 0 {
            return Err(bad("zero-sized HDR image"));
        }
        // Validate the attacker-controlled dimensions BEFORE allocating
        // anything: a malicious header can declare astronomically large
        // dims (the file itself is only a few dozen bytes), which would
        // OOM-abort the process on `Vec::with_capacity` and overflow the
        // unchecked `width * height` `usize` multiply in release.
        if width > MAX_HDR_DIM || height > MAX_HDR_DIM {
            return Err(bad(
                "HDR image dimensions exceed the maximum supported size",
            ));
        }
        let pixel_count = width
            .checked_mul(height)
            .ok_or_else(|| bad("HDR image dimensions overflow"))?;

        // --- pixel data ---
        let mut pixels = Vec::with_capacity(pixel_count);
        let body = &bytes[cursor..];
        let mut bi = 0usize;
        for _ in 0..height {
            let mut rgbe_row = vec![[0u8; 4]; width];
            decode_scanline(body, &mut bi, width, &mut rgbe_row)
                .map_err(|e| bad(&e))?;
            for rgbe in &rgbe_row {
                pixels.push(rgbe_to_linear(*rgbe));
            }
        }
        Self::from_pixels(width, height, pixels)
    }

    /// Bilinearly sample the radiance arriving from world `direction`
    /// (need not be unit length).
    ///
    /// Returns linear RGB scaled by [`intensity`](Self::intensity).
    /// The map's [`yaw`](Self::yaw) rotation about `+Y` is applied
    /// before the lookup. A degenerate (zero) direction returns the
    /// map's first pixel.
    pub fn sample_direction(&self, direction: [f32; 3]) -> [f32; 3] {
        let dir = match normalize3(direction) {
            Some(d) => d,
            None => return scale3(self.pixels[0], self.intensity),
        };
        let (u, v) = direction_to_equirect(dir, self.yaw);
        let raw = self.sample_uv(u, v);
        scale3(raw, self.intensity)
    }

    /// Bilinear lookup at equirectangular coordinates `(u, v)`,
    /// `u` wrapping and `v` clamping. Does *not* apply the intensity
    /// multiplier.
    pub fn sample_uv(&self, u: f32, v: f32) -> [f32; 3] {
        if self.width == 1 && self.height == 1 {
            return self.pixels[0];
        }
        // Pixel-centre sampling: texel (i, j) covers [i, i+1).
        let fx = u.rem_euclid(1.0) * self.width as f32 - 0.5;
        let fy = v.clamp(0.0, 1.0) * self.height as f32 - 0.5;
        let x0 = fx.floor() as isize;
        let y0 = fy.floor() as isize;
        let tx = fx - x0 as f32;
        let ty = fy - y0 as f32;
        let c00 = self.texel(x0, y0);
        let c10 = self.texel(x0 + 1, y0);
        let c01 = self.texel(x0, y0 + 1);
        let c11 = self.texel(x0 + 1, y0 + 1);
        let top = lerp3(c00, c10, tx);
        let bot = lerp3(c01, c11, tx);
        lerp3(top, bot, ty)
    }

    /// Fetch texel `(x, y)` with `x` wrapping (azimuth is periodic)
    /// and `y` clamping (the poles).
    fn texel(&self, x: isize, y: isize) -> [f32; 3] {
        let xi = x.rem_euclid(self.width as isize) as usize;
        let yi = y.clamp(0, self.height as isize - 1) as usize;
        self.pixels[yi * self.width + xi]
    }

    /// Diffuse irradiance arriving at a Lambertian surface whose
    /// outward unit normal is `normal`.
    ///
    /// This is the cosine-weighted hemisphere integral of the
    /// environment radiance:
    ///
    /// ```text
    /// E(n) = ∫_Ω L(ω) · max(0, n·ω) dω
    /// ```
    ///
    /// approximated by a uniform `samples × samples` quadrature over
    /// the hemisphere about `normal`. The result already includes the
    /// map's [`intensity`](Self::intensity); divide by π and multiply
    /// by albedo to get outgoing diffuse radiance.
    ///
    /// `samples` controls accuracy — 16–32 is plenty for a smooth
    /// map. The integral is normalised so a uniform white environment
    /// of radiance `L` returns exactly `π·L` (the analytic value of
    /// the cosine-weighted hemisphere integral).
    pub fn diffuse_irradiance(&self, normal: [f32; 3], samples: usize) -> [f32; 3] {
        let n = match normalize3(normal) {
            Some(v) => v,
            None => return [0.0; 3],
        };
        let samples = samples.max(2);
        // Orthonormal tangent frame about the normal.
        let (tangent, bitangent) = ortho_basis(n);
        let mut acc = [0.0f64; 3];
        let mut weight = 0.0f64;
        // Stratified hemisphere sampling: `samples` steps in the polar
        // angle θ and `samples` in the azimuth ψ.
        for i in 0..samples {
            // θ in (0, π/2); sample at stratum centres.
            let theta = std::f64::consts::FRAC_PI_2 * (i as f64 + 0.5) / samples as f64;
            let (st, ct) = (theta.sin(), theta.cos());
            for j in 0..samples {
                let psi = std::f64::consts::TAU * (j as f64 + 0.5) / samples as f64;
                // Direction in the local hemisphere frame.
                let local = [
                    st * psi.cos(),
                    st * psi.sin(),
                    ct, // +z is the normal
                ];
                // To world coordinates.
                let world = [
                    tangent[0] as f64 * local[0]
                        + bitangent[0] as f64 * local[1]
                        + n[0] as f64 * local[2],
                    tangent[1] as f64 * local[0]
                        + bitangent[1] as f64 * local[1]
                        + n[1] as f64 * local[2],
                    tangent[2] as f64 * local[0]
                        + bitangent[2] as f64 * local[1]
                        + n[2] as f64 * local[2],
                ];
                let radiance = self.sample_direction([
                    world[0] as f32,
                    world[1] as f32,
                    world[2] as f32,
                ]);
                // Cosine-weighted: the measure of the hemisphere
                // quadrature stratum is `sinθ·cosθ·dθ·dψ`. Accumulate
                // radiance · cosθ · sinθ.
                let w = ct * st;
                acc[0] += radiance[0] as f64 * w;
                acc[1] += radiance[1] as f64 * w;
                acc[2] += radiance[2] as f64 * w;
                weight += w;
            }
        }
        // The exact cosine-weighted hemisphere integral of a constant
        // radiance L is π·L. The quadrature sum of `cosθ·sinθ` over
        // the same grid is `weight`; normalising by `weight` and
        // scaling by π reproduces that analytic value, so a uniform
        // environment integrates to exactly π·L.
        if weight <= 0.0 {
            return [0.0; 3];
        }
        let k = std::f64::consts::PI / weight;
        [
            (acc[0] * k) as f32,
            (acc[1] * k) as f32,
            (acc[2] * k) as f32,
        ]
    }

    /// Precompute an [`IrradianceMap`] of the given resolution.
    ///
    /// Each texel of the output stores [`Self::diffuse_irradiance`]
    /// for the normal pointing in that texel's equirectangular
    /// direction — the renderer then looks diffuse lighting up by
    /// direction instead of re-integrating the hemisphere per shading
    /// call.
    ///
    /// `out_width`/`out_height` are small (32×16 is typical — the
    /// irradiance signal is very low-frequency). `samples` is the
    /// per-texel hemisphere quadrature resolution.
    ///
    /// # Errors
    ///
    /// [`RenderError::BadParameter`] for zero output dimensions.
    pub fn prefilter_irradiance(
        &self,
        out_width: usize,
        out_height: usize,
        samples: usize,
    ) -> Result<IrradianceMap, RenderError> {
        if out_width == 0 || out_height == 0 {
            return Err(RenderError::BadParameter {
                name: "irradiance_dimensions",
                reason: "irradiance map dimensions must be non-zero".into(),
            });
        }
        let mut pixels = Vec::with_capacity(out_width * out_height);
        for j in 0..out_height {
            let v = (j as f32 + 0.5) / out_height as f32;
            for i in 0..out_width {
                let u = (i as f32 + 0.5) / out_width as f32;
                let dir = equirect_to_direction(u, v);
                pixels.push(self.diffuse_irradiance(dir, samples));
            }
        }
        Ok(IrradianceMap {
            width: out_width,
            height: out_height,
            pixels,
        })
    }
}

/// A precomputed diffuse-irradiance map — the convolved, low-frequency
/// result of [`EnvironmentMap::prefilter_irradiance`].
#[derive(Clone, Debug, PartialEq)]
pub struct IrradianceMap {
    /// Width in texels.
    pub width: usize,
    /// Height in texels.
    pub height: usize,
    /// `width · height` irradiance triples, row-major.
    pub pixels: Vec<[f32; 3]>,
}

impl IrradianceMap {
    /// Look up the diffuse irradiance for a surface `normal` by
    /// nearest-texel equirectangular lookup. The irradiance signal is
    /// very smooth, so nearest sampling on a small map is adequate.
    pub fn irradiance(&self, normal: [f32; 3]) -> [f32; 3] {
        let n = match normalize3(normal) {
            Some(v) => v,
            None => return [0.0; 3],
        };
        let (u, v) = direction_to_equirect(n, 0.0);
        let x = ((u.rem_euclid(1.0) * self.width as f32) as usize).min(self.width - 1);
        let y = ((v.clamp(0.0, 1.0) * self.height as f32) as usize).min(self.height - 1);
        self.pixels[y * self.width + x]
    }
}

/// One roughness level of a [`PrefilteredEnvironment`] — the
/// environment map convolved with the GGX specular lobe for a single
/// roughness value.
#[derive(Clone, Debug, PartialEq)]
pub struct PrefilteredLevel {
    /// Perceptual roughness this level was convolved at, `[0, 1]`.
    pub roughness: f32,
    /// Equirectangular width of this level (coarser for rougher
    /// levels — a rough convolution is low-frequency).
    pub width: usize,
    /// Equirectangular height of this level.
    pub height: usize,
    /// `width · height` linear-RGB radiance triples, row-major.
    pub pixels: Vec<[f32; 3]>,
}

/// A **prefiltered specular environment** — the `roughness`-indexed
/// mip chain of the split-sum image-based-lighting approximation
/// (Karis 2013, the UE4 IBL technique).
///
/// The specular IBL integral
/// `∫ Lᵢ(l)·BRDF(l,v)·(n·l) dl` is split into two factors that are
/// each precomputed once:
///
/// 1. the **prefiltered environment** — this type — the environment
///    map convolved with the GGX lobe at each roughness level, so a
///    renderer samples it by `(reflection direction, roughness)`;
/// 2. the **BRDF integration LUT** — a `(n·v, roughness) → (scale,
///    bias)` table — see `valenx_render_bridge::pbr`'s LUT.
///
/// The specular IBL is then reconstructed as
/// `prefiltered(r, roughness) · (F₀·scale + bias)`.
#[derive(Clone, Debug, PartialEq)]
pub struct PrefilteredEnvironment {
    /// The roughness levels, ascending in roughness. `levels[0]` is the
    /// sharp (roughness ≈ 0) mirror level.
    pub levels: Vec<PrefilteredLevel>,
}

impl PrefilteredEnvironment {
    /// Sample the prefiltered radiance arriving along world
    /// `reflection` at perceptual `roughness`.
    ///
    /// The roughness selects (and linearly blends between) the two
    /// bracketing mip levels; within a level the lookup is a bilinear
    /// equirectangular sample. This is exactly the GPU texture fetch a
    /// real-time renderer issues against a prefiltered cubemap.
    pub fn sample(&self, reflection: [f32; 3], roughness: f32) -> [f32; 3] {
        if self.levels.is_empty() {
            return [0.0; 3];
        }
        let dir = match normalize3(reflection) {
            Some(d) => d,
            None => return self.sample_level(0, 1.0, 1.0),
        };
        let (u, v) = direction_to_equirect(dir, 0.0);
        let r = roughness.clamp(0.0, 1.0);
        // Map roughness to a fractional level index. Levels are stored
        // ascending in roughness; assume they are evenly spaced.
        let last = self.levels.len() - 1;
        let scaled = r * last as f32;
        let lo = (scaled.floor() as usize).min(last);
        let hi = (lo + 1).min(last);
        let frac = scaled - lo as f32;
        let a = self.sample_level(lo, u, v);
        let b = self.sample_level(hi, u, v);
        [
            a[0] + (b[0] - a[0]) * frac,
            a[1] + (b[1] - a[1]) * frac,
            a[2] + (b[2] - a[2]) * frac,
        ]
    }

    /// Bilinear equirectangular lookup within level `idx`.
    fn sample_level(&self, idx: usize, u: f32, v: f32) -> [f32; 3] {
        let lvl = &self.levels[idx];
        if lvl.width == 0 || lvl.height == 0 {
            return [0.0; 3];
        }
        if lvl.width == 1 && lvl.height == 1 {
            return lvl.pixels[0];
        }
        let fx = u.rem_euclid(1.0) * lvl.width as f32 - 0.5;
        let fy = v.clamp(0.0, 1.0) * lvl.height as f32 - 0.5;
        let x0 = fx.floor() as isize;
        let y0 = fy.floor() as isize;
        let tx = fx - x0 as f32;
        let ty = fy - y0 as f32;
        let texel = |x: isize, y: isize| -> [f32; 3] {
            let xi = x.rem_euclid(lvl.width as isize) as usize;
            let yi = y.clamp(0, lvl.height as isize - 1) as usize;
            lvl.pixels[yi * lvl.width + xi]
        };
        let top = lerp3(texel(x0, y0), texel(x0 + 1, y0), tx);
        let bot = lerp3(texel(x0, y0 + 1), texel(x0 + 1, y0 + 1), tx);
        lerp3(top, bot, ty)
    }
}

impl EnvironmentMap {
    /// Prefilter the environment into a [`PrefilteredEnvironment`] — the
    /// **specular split-sum mip chain** (Phase 30.7).
    ///
    /// For each of `levels` roughness values (evenly spaced over
    /// `[0, 1]`) the environment is convolved with the GGX specular
    /// lobe by **GGX importance sampling**: `samples` half-vectors are
    /// drawn from the GGX distribution about the level's reflection
    /// direction (using the standard `n = v = r` simplification), each
    /// reflected to a light direction, and the environment radiance
    /// there is accumulated weighted by `n·l`. Roughness 0 reproduces
    /// a sharp mirror; roughness 1 a broad, near-diffuse blur.
    ///
    /// Coarser (rougher) levels are stored at a lower resolution — a
    /// heavily-convolved map is low-frequency, so a small grid loses
    /// nothing and keeps the precompute cheap.
    ///
    /// `base_width`/`base_height` size the sharp (roughness-0) level;
    /// `levels` is the number of roughness steps (5–6 is typical);
    /// `samples` is the GGX importance-sample count per texel (32–64
    /// is plenty for a smooth map).
    ///
    /// # Errors
    ///
    /// [`RenderError::BadParameter`] for zero dimensions or fewer than
    /// 2 levels (a single level cannot span the roughness range).
    pub fn prefilter_specular(
        &self,
        base_width: usize,
        base_height: usize,
        levels: usize,
        samples: usize,
    ) -> Result<PrefilteredEnvironment, RenderError> {
        if base_width == 0 || base_height == 0 {
            return Err(RenderError::BadParameter {
                name: "prefilter_dimensions",
                reason: "prefiltered environment dimensions must be non-zero".into(),
            });
        }
        if levels < 2 {
            return Err(RenderError::BadParameter {
                name: "levels",
                reason: "need at least 2 roughness levels".into(),
            });
        }
        let samples = samples.max(1);
        let mut out_levels = Vec::with_capacity(levels);
        for li in 0..levels {
            let roughness = li as f32 / (levels - 1) as f32;
            // Halve the resolution every other level — rough levels
            // are low-frequency. Keep at least a 4×2 grid.
            let shrink = 1usize << (li / 2);
            let w = (base_width / shrink).max(4);
            let h = (base_height / shrink).max(2);
            let mut pixels = Vec::with_capacity(w * h);
            for j in 0..h {
                let v = (j as f32 + 0.5) / h as f32;
                for i in 0..w {
                    let u = (i as f32 + 0.5) / w as f32;
                    let dir = equirect_to_direction(u, v);
                    pixels.push(self.prefilter_direction(dir, roughness, samples));
                }
            }
            out_levels.push(PrefilteredLevel {
                roughness,
                width: w,
                height: h,
                pixels,
            });
        }
        Ok(PrefilteredEnvironment { levels: out_levels })
    }

    /// Convolve the environment with the GGX lobe for one reflection
    /// `direction` at one `roughness` — the per-texel kernel of
    /// [`Self::prefilter_specular`].
    ///
    /// Uses the split-sum simplification `n = v = r`: the half-vectors
    /// are GGX-importance-sampled about `direction`, each reflected to
    /// a light direction `l`, and the environment radiance there is
    /// accumulated weighted by `n·l` (the standard weighting that
    /// biases toward the lobe centre and is normalised out).
    fn prefilter_direction(&self, direction: [f32; 3], roughness: f32, samples: usize) -> [f32; 3] {
        let n = match normalize3(direction) {
            Some(d) => d,
            None => return [0.0; 3],
        };
        // Roughness 0: a perfect mirror — just sample the environment.
        if roughness <= 1e-4 {
            return self.sample_direction(n);
        }
        let alpha = roughness * roughness;
        let (tangent, bitangent) = ortho_basis(n);
        let mut acc = [0.0f64; 3];
        let mut total_weight = 0.0f64;
        // The split-sum assumption: the view and the normal both equal
        // the reflection direction.
        let view = n;
        for s in 0..samples {
            // Hammersley low-discrepancy 2-D sample.
            let (u1, u2) = hammersley(s, samples);
            // GGX importance-sampled half-vector in the tangent frame.
            let phi = std::f32::consts::TAU * u1;
            let cos_theta = (((1.0 - u2) / (1.0 + (alpha * alpha - 1.0) * u2)).max(0.0)).sqrt();
            let sin_theta = (1.0 - cos_theta * cos_theta).max(0.0).sqrt();
            let h_local = [
                sin_theta * phi.cos(),
                sin_theta * phi.sin(),
                cos_theta,
            ];
            // Half-vector to world space.
            let h = [
                tangent[0] * h_local[0] + bitangent[0] * h_local[1] + n[0] * h_local[2],
                tangent[1] * h_local[0] + bitangent[1] * h_local[1] + n[1] * h_local[2],
                tangent[2] * h_local[0] + bitangent[2] * h_local[1] + n[2] * h_local[2],
            ];
            // Reflect the view about the half-vector to get the light
            // direction:  l = 2·(v·h)·h − v.
            let vh = view[0] * h[0] + view[1] * h[1] + view[2] * h[2];
            let l = [
                2.0 * vh * h[0] - view[0],
                2.0 * vh * h[1] - view[1],
                2.0 * vh * h[2] - view[2],
            ];
            let n_dot_l = n[0] * l[0] + n[1] * l[1] + n[2] * l[2];
            if n_dot_l > 0.0 {
                let radiance = self.sample_direction(l);
                acc[0] += radiance[0] as f64 * n_dot_l as f64;
                acc[1] += radiance[1] as f64 * n_dot_l as f64;
                acc[2] += radiance[2] as f64 * n_dot_l as f64;
                total_weight += n_dot_l as f64;
            }
        }
        if total_weight <= 0.0 {
            // Degenerate — fall back to the sharp sample.
            return self.sample_direction(n);
        }
        [
            (acc[0] / total_weight) as f32,
            (acc[1] / total_weight) as f32,
            (acc[2] / total_weight) as f32,
        ]
    }
}

/// The Hammersley low-discrepancy sequence point `i` of `n` — a 2-D
/// `(u1, u2)` quasi-random sample. `u1` is `i/n`; `u2` is the
/// radical-inverse (Van der Corput) of `i` in base 2. Low-discrepancy
/// sampling converges far faster than plain pseudo-random for the GGX
/// prefilter integral.
fn hammersley(i: usize, n: usize) -> (f32, f32) {
    // Van der Corput radical inverse, base 2: reverse the bits of `i`
    // (swap pairs of every bit-width down to single bits).
    let mut bits = (i as u32).rotate_right(16);
    bits = ((bits & 0x55555555) << 1) | ((bits & 0xAAAAAAAA) >> 1);
    bits = ((bits & 0x33333333) << 2) | ((bits & 0xCCCCCCCC) >> 2);
    bits = ((bits & 0x0F0F0F0F) << 4) | ((bits & 0xF0F0F0F0) >> 4);
    bits = ((bits & 0x00FF00FF) << 8) | ((bits & 0xFF00FF00) >> 8);
    let radical_inverse = bits as f32 * 2.328_306_4e-10; // / 2^32
    let u1 = if n == 0 { 0.0 } else { i as f32 / n as f32 };
    (u1, radical_inverse)
}

// --- equirectangular <-> direction conversions ---

/// World direction → equirectangular `(u, v)`. `+Y` is up; `yaw`
/// rotates the lookup about `+Y`.
fn direction_to_equirect(dir: [f32; 3], yaw: f32) -> (f32, f32) {
    // Apply the inverse yaw so a positive `yaw` spins the map.
    let (s, c) = (-yaw).sin_cos();
    let x = c * dir[0] - s * dir[2];
    let z = s * dir[0] + c * dir[2];
    let y = dir[1].clamp(-1.0, 1.0);
    // Azimuth φ ∈ [0, 2π).
    let phi = z.atan2(x);
    let u = (phi / std::f32::consts::TAU).rem_euclid(1.0);
    // Polar θ ∈ [0, π], v = 0 at +Y.
    let theta = y.acos();
    let v = theta / std::f32::consts::PI;
    (u, v)
}

/// Equirectangular `(u, v)` → unit world direction (`+Y` up).
fn equirect_to_direction(u: f32, v: f32) -> [f32; 3] {
    let phi = u * std::f32::consts::TAU;
    let theta = v * std::f32::consts::PI;
    let st = theta.sin();
    [st * phi.cos(), theta.cos(), st * phi.sin()]
}

// --- RGBE decode ---

/// Decode an RGBE quad to linear RGB: `c · 2^(e-128) / 256`.
fn rgbe_to_linear(rgbe: [u8; 4]) -> [f32; 3] {
    let e = rgbe[3];
    if e == 0 {
        return [0.0; 3];
    }
    // ldexp(1, e-136) == 2^(e-128) / 256.
    let f = libm_ldexp(1.0, e as i32 - 136);
    [
        rgbe[0] as f32 * f,
        rgbe[1] as f32 * f,
        rgbe[2] as f32 * f,
    ]
}

/// `mantissa · 2^exp` without pulling in `libm` — a small loop over
/// the exponent is fine here (HDR exponents are bounded).
fn libm_ldexp(mantissa: f32, exp: i32) -> f32 {
    // f32::powi is in std; this is exact for the bounded exponents an
    // RGBE byte produces.
    mantissa * 2.0f32.powi(exp)
}

/// Decode one HDR scanline of `width` RGBE pixels into `out`.
///
/// Handles the new-style adaptive RLE (signalled by a `2 2 hi lo`
/// header where `hi·256+lo == width`) and falls back to the
/// flat / old-RLE scanline otherwise.
fn decode_scanline(
    body: &[u8],
    bi: &mut usize,
    width: usize,
    out: &mut [[u8; 4]],
) -> Result<(), String> {
    let need = |bi: usize, n: usize| -> Result<(), String> {
        if bi + n > body.len() {
            Err("truncated HDR pixel stream".to_string())
        } else {
            Ok(())
        }
    };
    // New-style adaptive RLE applies only for 8..=0x7fff-wide rows.
    if (8..=0x7fff).contains(&width) {
        need(*bi, 4)?;
        let header = [
            body[*bi],
            body[*bi + 1],
            body[*bi + 2],
            body[*bi + 3],
        ];
        if header[0] == 2
            && header[1] == 2
            && ((header[2] as usize) << 8 | header[3] as usize) == width
        {
            *bi += 4;
            // Four component planes, each RLE-encoded separately.
            for comp in 0..4 {
                let mut x = 0usize;
                while x < width {
                    need(*bi, 1)?;
                    let count = body[*bi];
                    *bi += 1;
                    if count > 128 {
                        // A run of `count-128` copies of one value.
                        let run = (count as usize) - 128;
                        need(*bi, 1)?;
                        let value = body[*bi];
                        *bi += 1;
                        for _ in 0..run {
                            if x >= width {
                                return Err("HDR RLE run overflow".into());
                            }
                            out[x][comp] = value;
                            x += 1;
                        }
                    } else {
                        // `count` literal bytes.
                        let lit = count as usize;
                        need(*bi, lit)?;
                        for _ in 0..lit {
                            if x >= width {
                                return Err("HDR literal overflow".into());
                            }
                            out[x][comp] = body[*bi];
                            *bi += 1;
                            x += 1;
                        }
                    }
                }
            }
            return Ok(());
        }
    }
    // Flat / old-style scanline: `width` raw RGBE quads, with the
    // old-RLE escape (an RGBE of `1 1 1 n` repeats the previous
    // pixel `n` times, possibly shifted for long runs).
    let mut x = 0usize;
    let mut shift = 0u32;
    while x < width {
        need(*bi, 4)?;
        let q = [
            body[*bi],
            body[*bi + 1],
            body[*bi + 2],
            body[*bi + 3],
        ];
        *bi += 4;
        if q[0] == 1 && q[1] == 1 && q[2] == 1 {
            // Old-style RLE repeat. `shift` grows by 8 per consecutive
            // escape; a crafted run of zero-length escapes (`1 1 1 0`,
            // which never advances `x`) would drive it to 64 and make
            // the bare `<< shift` panic in debug builds ("attempt to
            // shift left with overflow"). `checked_shl` yields the same
            // benign `run = 0` the release build already produces, so a
            // hostile scanline errors out instead of panicking. Real
            // files never reach this: a valid run is bounded by the
            // scanline width (<= 0x7fff), needing at most two escapes.
            let run = (q[3] as usize).checked_shl(shift).unwrap_or(0);
            let prev = if x > 0 { out[x - 1] } else { [0u8; 4] };
            for _ in 0..run {
                if x >= width {
                    return Err("HDR old-RLE overflow".into());
                }
                out[x] = prev;
                x += 1;
            }
            // `saturating_add` so a hostile all-escape stream can't wrap
            // `shift` (a u32) either — together with the `checked_shl`
            // above this makes the escape branch panic-free on any input.
            shift = shift.saturating_add(8);
        } else {
            out[x] = q;
            x += 1;
            shift = 0;
        }
    }
    Ok(())
}

/// Read one `\n`-terminated line from `bytes` starting at `*cursor`,
/// advancing the cursor past the newline. Returns the line without
/// its terminator, or `None` at end-of-input.
fn read_line(bytes: &[u8], cursor: &mut usize) -> Option<String> {
    if *cursor >= bytes.len() {
        return None;
    }
    let start = *cursor;
    let mut end = start;
    while end < bytes.len() && bytes[end] != b'\n' {
        end += 1;
    }
    let line = String::from_utf8_lossy(&bytes[start..end])
        .trim_end_matches('\r')
        .to_string();
    *cursor = (end + 1).min(bytes.len());
    Some(line)
}

// --- small vector helpers ---

fn normalize3(v: [f32; 3]) -> Option<[f32; 3]> {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len < 1e-12 {
        None
    } else {
        Some([v[0] / len, v[1] / len, v[2] / len])
    }
}

fn scale3(v: [f32; 3], s: f32) -> [f32; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

fn lerp3(a: [f32; 3], b: [f32; 3], t: f32) -> [f32; 3] {
    [
        a[0] + (b[0] - a[0]) * t,
        a[1] + (b[1] - a[1]) * t,
        a[2] + (b[2] - a[2]) * t,
    ]
}

/// An orthonormal `(tangent, bitangent)` pair spanning the plane
/// perpendicular to the unit vector `n`.
fn ortho_basis(n: [f32; 3]) -> ([f32; 3], [f32; 3]) {
    // Pick the world axis least aligned with `n` as a seed.
    let seed = if n[0].abs() <= n[1].abs() && n[0].abs() <= n[2].abs() {
        [1.0, 0.0, 0.0]
    } else if n[1].abs() <= n[2].abs() {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };
    let t = normalize3(cross3(seed, n)).unwrap_or([1.0, 0.0, 0.0]);
    let b = cross3(n, t);
    (t, b)
}

fn cross3(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Encode a linear-RGB pixel back to an RGBE quad — the inverse of
    /// [`rgbe_to_linear`], used to build synthetic HDR test data.
    fn linear_to_rgbe(c: [f32; 3]) -> [u8; 4] {
        let max = c[0].max(c[1]).max(c[2]);
        if max < 1e-32 {
            return [0, 0, 0, 0];
        }
        let (frac, exp) = frexp(max);
        let scale = frac * 256.0 / max;
        [
            (c[0] * scale) as u8,
            (c[1] * scale) as u8,
            (c[2] * scale) as u8,
            (exp + 128) as u8,
        ]
    }

    /// Decompose `x` into `(mantissa ∈ [0.5, 1), exponent)`.
    fn frexp(x: f32) -> (f32, i32) {
        if x == 0.0 {
            return (0.0, 0);
        }
        let mut exp = 0;
        let mut m = x;
        while m >= 1.0 {
            m *= 0.5;
            exp += 1;
        }
        while m < 0.5 {
            m *= 2.0;
            exp -= 1;
        }
        (m, exp)
    }

    /// Build a minimal flat-scanline HDR byte stream of one colour.
    fn synthetic_hdr(width: usize, height: usize, color: [f32; 3]) -> Vec<u8> {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(b"#?RADIANCE\n");
        bytes.extend_from_slice(b"FORMAT=32-bit_rle_rgbe\n");
        bytes.extend_from_slice(b"\n");
        bytes.extend_from_slice(format!("-Y {height} +X {width}\n").as_bytes());
        let rgbe = linear_to_rgbe(color);
        // Force the flat path: a width < 8 never takes the new RLE.
        for _ in 0..(width * height) {
            bytes.extend_from_slice(&rgbe);
        }
        bytes
    }

    #[test]
    fn from_pixels_rejects_bad_dimensions() {
        assert!(EnvironmentMap::from_pixels(0, 4, vec![]).is_err());
        assert!(EnvironmentMap::from_pixels(2, 2, vec![[0.0; 3]; 3]).is_err());
        assert!(EnvironmentMap::from_pixels(2, 2, vec![[0.0; 3]; 4]).is_ok());
    }

    #[test]
    fn rgbe_round_trips_through_linear() {
        // A mid-grey pixel encodes and decodes back to ~itself.
        let original = [0.5f32, 0.25, 0.75];
        let rgbe = linear_to_rgbe(original);
        let back = rgbe_to_linear(rgbe);
        for k in 0..3 {
            assert!(
                (back[k] - original[k]).abs() < 0.01,
                "channel {k}: {} vs {}",
                back[k],
                original[k]
            );
        }
    }

    #[test]
    fn rgbe_zero_exponent_is_black() {
        assert_eq!(rgbe_to_linear([200, 200, 200, 0]), [0.0; 3]);
    }

    #[test]
    fn radiance_hdr_decodes_a_flat_image() {
        let bytes = synthetic_hdr(4, 3, [0.6, 0.3, 0.1]);
        let map = EnvironmentMap::from_radiance_hdr(&bytes).unwrap();
        assert_eq!(map.width, 4);
        assert_eq!(map.height, 3);
        assert_eq!(map.pixels.len(), 12);
        // Every pixel decodes to roughly the encoded colour.
        for p in &map.pixels {
            assert!((p[0] - 0.6).abs() < 0.02, "got {p:?}");
            assert!((p[1] - 0.3).abs() < 0.02);
            assert!((p[2] - 0.1).abs() < 0.02);
        }
    }

    #[test]
    fn radiance_hdr_rejects_missing_magic() {
        let bad = b"not an hdr file\n\n-Y 1 +X 1\n";
        assert!(EnvironmentMap::from_radiance_hdr(bad).is_err());
    }

    #[test]
    fn radiance_hdr_rejects_bad_resolution_line() {
        let bad = b"#?RADIANCE\n\n+X 4 -Y 3\n";
        assert!(EnvironmentMap::from_radiance_hdr(bad).is_err());
    }

    #[test]
    fn radiance_hdr_rejects_oversized_dimensions_before_allocating() {
        // A minimal valid header that declares dimensions one past the
        // cap, then EOF (no pixel body). Pre-fix, `from_radiance_hdr`
        // ran `Vec::with_capacity(width * height)` immediately — for the
        // real-world attack header `-Y 1000000000 +X 1000000000` that is
        // a ~10^18-element request that OOM-aborts the process (and the
        // multiply wraps in release). `MAX_HDR_DIM + 1` per side proves
        // the dimension guard rejects the header *before* any allocation
        // while keeping the test fast and memory-safe.
        let over = MAX_HDR_DIM + 1;
        let header = format!("#?RADIANCE\nFORMAT=32-bit_rle_rgbe\n\n-Y {over} +X {over}\n");
        let err = EnvironmentMap::from_radiance_hdr(header.as_bytes())
            .expect_err("oversized HDR dimensions must be rejected");
        assert!(
            matches!(err, RenderError::BadParameter { .. }),
            "expected BadParameter, got {err:?}"
        );

        // A height past the cap with a sane width is rejected too.
        let header2 = format!("#?RADIANCE\nFORMAT=32-bit_rle_rgbe\n\n-Y {over} +X 4\n");
        assert!(EnvironmentMap::from_radiance_hdr(header2.as_bytes()).is_err());
    }

    #[test]
    fn radiance_hdr_accepts_dimensions_at_the_cap() {
        // A header exactly at the cap must still parse the header and
        // proceed (here it then errors on the truncated pixel stream,
        // NOT on the dimension check — the dims themselves are legal).
        // We use a 1-row image at max width to avoid allocating the full
        // 16384² buffer: width is at the cap, height is 1.
        let header = format!(
            "#?RADIANCE\nFORMAT=32-bit_rle_rgbe\n\n-Y 1 +X {MAX_HDR_DIM}\n"
        );
        let err = EnvironmentMap::from_radiance_hdr(header.as_bytes())
            .expect_err("truncated body should error");
        let msg = err.to_string();
        // It must fail on the missing pixels, not reject the dimension.
        assert!(
            !msg.contains("exceed the maximum"),
            "a dimension at the cap must be accepted, got: {msg}"
        );
    }

    #[test]
    fn radiance_hdr_old_rle_shift_does_not_overflow() {
        // Regression: the old-style RLE escape `1 1 1 n` repeats the
        // previous pixel `n << shift` times and bumps `shift` by 8 after
        // every consecutive escape. A crafted run of zero-length escapes
        // (`01 01 01 00`) never advances `x`, so the width guard never
        // fires and `shift` climbs 8, 16, 24 … 64. Pre-fix the 9th escape
        // evaluated `0usize << 64`, panicking in debug builds with
        // "attempt to shift left with overflow". The decoder must instead
        // fail gracefully (truncated row) without panicking.
        let mut bytes = Vec::new();
        bytes.extend_from_slice(b"#?RADIANCE\n");
        bytes.extend_from_slice(b"FORMAT=32-bit_rle_rgbe\n\n");
        bytes.extend_from_slice(b"-Y 1 +X 4\n"); // width 4 (< 8) → old-style scanline
        // Twelve zero-length escapes drive `shift` well past usize::BITS
        // (and past 32 on a 32-bit target) without ever filling the row.
        for _ in 0..12 {
            bytes.extend_from_slice(&[1, 1, 1, 0]);
        }
        // Must return an error (the 4-pixel row is never filled), not panic.
        assert!(EnvironmentMap::from_radiance_hdr(&bytes).is_err());
    }

    #[test]
    fn radiance_hdr_decodes_new_style_rle() {
        // A 16-wide row (>= 8 → eligible for adaptive RLE). Encode each
        // of the 4 component planes as a single run of one value.
        let mut bytes = Vec::new();
        bytes.extend_from_slice(b"#?RADIANCE\n");
        bytes.extend_from_slice(b"FORMAT=32-bit_rle_rgbe\n\n");
        bytes.extend_from_slice(b"-Y 1 +X 16\n");
        // New-RLE scanline header: 2, 2, width-hi, width-lo.
        bytes.extend_from_slice(&[2, 2, 0, 16]);
        // Component planes r, g, b, e — each "run of 16 copies of v".
        for v in [120u8, 90, 60, 130] {
            bytes.push(128 + 16); // run length 16
            bytes.push(v);
        }
        let map = EnvironmentMap::from_radiance_hdr(&bytes).unwrap();
        assert_eq!(map.width, 16);
        assert_eq!(map.height, 1);
        let expect = rgbe_to_linear([120, 90, 60, 130]);
        for p in &map.pixels {
            assert!((p[0] - expect[0]).abs() < 1e-6, "got {p:?}");
        }
    }

    #[test]
    fn sample_direction_of_uniform_map_is_constant() {
        let map = EnvironmentMap::uniform([0.4, 0.4, 0.4]);
        let a = map.sample_direction([1.0, 0.0, 0.0]);
        let b = map.sample_direction([0.0, -1.0, 0.3]);
        assert_eq!(a, b);
        assert!((a[0] - 0.4).abs() < 1e-6);
    }

    #[test]
    fn intensity_scales_the_sample() {
        let mut map = EnvironmentMap::uniform([0.5, 0.5, 0.5]);
        map.intensity = 3.0;
        let s = map.sample_direction([0.0, 1.0, 0.0]);
        assert!((s[0] - 1.5).abs() < 1e-6);
    }

    #[test]
    fn diffuse_irradiance_of_white_environment_is_pi() {
        // The cosine-weighted hemisphere integral of a uniform white
        // environment of radiance 1 is exactly π — the canonical IBL
        // sanity check.
        let map = EnvironmentMap::uniform([1.0, 1.0, 1.0]);
        let e = map.diffuse_irradiance([0.0, 1.0, 0.0], 24);
        for (k, &channel) in e.iter().enumerate() {
            assert!(
                (channel - std::f32::consts::PI).abs() < 1e-3,
                "channel {k} irradiance {channel} should be π"
            );
        }
    }

    #[test]
    fn diffuse_irradiance_is_higher_facing_a_bright_hemisphere() {
        // A map that is bright in the upper half (+Y) and dark in the
        // lower half. A normal pointing up must collect far more
        // irradiance than one pointing down.
        let w = 8;
        let h = 8;
        let mut pixels = vec![[0.0f32; 3]; w * h];
        for j in 0..h {
            // v = 0 is +Y (up); the top half is bright.
            let bright = j < h / 2;
            for i in 0..w {
                pixels[j * w + i] = if bright {
                    [2.0, 2.0, 2.0]
                } else {
                    [0.0, 0.0, 0.0]
                };
            }
        }
        let map = EnvironmentMap::from_pixels(w, h, pixels).unwrap();
        let up = map.diffuse_irradiance([0.0, 1.0, 0.0], 24);
        let down = map.diffuse_irradiance([0.0, -1.0, 0.0], 24);
        assert!(
            up[0] > down[0] + 1.0,
            "up-facing irradiance {} should dwarf down-facing {}",
            up[0],
            down[0]
        );
        // The up normal sees only the bright hemisphere → close to π·2.
        assert!(up[0] > std::f32::consts::PI, "up irradiance {}", up[0]);
    }

    #[test]
    fn prefilter_irradiance_builds_a_lookup_map() {
        let map = EnvironmentMap::uniform([0.5, 0.5, 0.5]);
        let irr = map.prefilter_irradiance(16, 8, 16).unwrap();
        assert_eq!(irr.width, 16);
        assert_eq!(irr.height, 8);
        assert_eq!(irr.pixels.len(), 128);
        // Uniform environment → uniform irradiance (= π · 0.5).
        let look = irr.irradiance([0.3, 0.7, -0.2]);
        let expect = std::f32::consts::PI * 0.5;
        assert!((look[0] - expect).abs() < 1e-2, "got {look:?}");
    }

    #[test]
    fn prefilter_rejects_zero_dimensions() {
        let map = EnvironmentMap::uniform([1.0; 3]);
        assert!(map.prefilter_irradiance(0, 8, 8).is_err());
    }

    #[test]
    fn equirect_round_trips_a_direction() {
        // direction → (u, v) → direction should be (near-)identity.
        for dir in [
            [0.0f32, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            normalize3([1.0, 1.0, 1.0]).unwrap(),
        ] {
            let (u, v) = direction_to_equirect(dir, 0.0);
            let back = equirect_to_direction(u, v);
            for k in 0..3 {
                assert!(
                    (back[k] - dir[k]).abs() < 1e-4,
                    "axis {k}: {} vs {}",
                    back[k],
                    dir[k]
                );
            }
        }
    }

    #[test]
    fn prefilter_specular_rejects_bad_parameters() {
        let map = EnvironmentMap::uniform([1.0; 3]);
        // Zero dimensions.
        assert!(map.prefilter_specular(0, 8, 4, 16).is_err());
        // Fewer than 2 levels.
        assert!(map.prefilter_specular(16, 8, 1, 16).is_err());
        // A valid call succeeds.
        assert!(map.prefilter_specular(16, 8, 4, 16).is_ok());
    }

    #[test]
    fn prefilter_specular_of_uniform_environment_stays_uniform() {
        // Convolving a flat-colour environment with any lobe leaves it
        // the same flat colour — the GGX prefilter must not invent
        // energy or colour.
        let map = EnvironmentMap::uniform([0.5, 0.5, 0.5]);
        let pre = map.prefilter_specular(16, 8, 5, 32).unwrap();
        assert_eq!(pre.levels.len(), 5);
        // Sample at several roughnesses — all should read ~0.5.
        for &r in &[0.0f32, 0.25, 0.5, 0.75, 1.0] {
            let s = pre.sample([1.0, 0.0, 0.0], r);
            for (k, &c) in s.iter().enumerate() {
                assert!(
                    (c - 0.5).abs() < 1e-2,
                    "prefiltered channel {k} at roughness {r} = {c}, expected ~0.5"
                );
            }
        }
    }

    #[test]
    fn prefilter_roughness_zero_is_a_sharp_mirror() {
        // At roughness 0 the prefilter must reproduce the raw
        // environment (a perfect mirror reflection). Build a map with
        // a distinct bright spot and confirm the sharp level still
        // sees it strongly.
        let w = 16;
        let h = 8;
        let mut pixels = vec![[0.1f32; 3]; w * h];
        // One very bright texel.
        pixels[2 * w + 4] = [20.0, 20.0, 20.0];
        let map = EnvironmentMap::from_pixels(w, h, pixels).unwrap();
        let pre = map.prefilter_specular(16, 8, 5, 32).unwrap();
        // The sharp level (roughness 0) should be much higher-contrast
        // than the roughest level: take the brightest texel of each.
        let sharp_max = pre.levels[0]
            .pixels
            .iter()
            .map(|p| p[0])
            .fold(0.0f32, f32::max);
        let rough_max = pre
            .levels
            .last()
            .unwrap()
            .pixels
            .iter()
            .map(|p| p[0])
            .fold(0.0f32, f32::max);
        assert!(
            sharp_max > rough_max,
            "sharp level peak {sharp_max} should exceed the blurred peak {rough_max}"
        );
    }

    #[test]
    fn prefilter_rougher_levels_are_smoother() {
        // A high-contrast environment: convolving it with a wider GGX
        // lobe (higher roughness) must reduce the contrast. Measure
        // contrast as (max − min) of the red channel per level.
        let w = 32;
        let h = 16;
        let mut pixels = vec![[0.0f32; 3]; w * h];
        // A bright stripe down one column.
        for j in 0..h {
            pixels[j * w + 8] = [10.0, 10.0, 10.0];
        }
        let map = EnvironmentMap::from_pixels(w, h, pixels).unwrap();
        let pre = map.prefilter_specular(32, 16, 5, 48).unwrap();
        let contrast = |lvl: &PrefilteredLevel| -> f32 {
            let mx = lvl.pixels.iter().map(|p| p[0]).fold(0.0f32, f32::max);
            let mn = lvl.pixels.iter().map(|p| p[0]).fold(f32::MAX, f32::min);
            mx - mn
        };
        let sharp_contrast = contrast(&pre.levels[0]);
        let rough_contrast = contrast(pre.levels.last().unwrap());
        assert!(
            rough_contrast < sharp_contrast,
            "rough level contrast {rough_contrast} should be below sharp {sharp_contrast}"
        );
    }

    #[test]
    fn hammersley_samples_lie_in_unit_square() {
        // Every Hammersley point must fall inside [0, 1)².
        for i in 0..64 {
            let (u1, u2) = hammersley(i, 64);
            assert!((0.0..1.0).contains(&u1), "u1 {u1} out of range");
            assert!((0.0..1.0).contains(&u2), "u2 {u2} out of range");
        }
        // The first point is the origin (i = 0).
        let (u1, u2) = hammersley(0, 64);
        assert!(u1 == 0.0 && u2 == 0.0);
    }

    #[test]
    fn yaw_rotates_the_lookup() {
        // A map bright on one azimuth column only; a 180° yaw must
        // move which direction sees the bright column.
        let w = 4;
        let h = 2;
        let mut pixels = vec![[0.0f32; 3]; w * h];
        // Column 0 bright.
        for j in 0..h {
            pixels[j * w] = [5.0, 5.0, 5.0];
        }
        let mut map = EnvironmentMap::from_pixels(w, h, pixels).unwrap();
        let dir = [1.0, 0.0, 0.0];
        let before = map.sample_direction(dir);
        map.yaw = std::f32::consts::PI;
        let after = map.sample_direction(dir);
        assert!(
            (before[0] - after[0]).abs() > 1.0,
            "yaw should change the sampled radiance: {before:?} vs {after:?}"
        );
    }
}
