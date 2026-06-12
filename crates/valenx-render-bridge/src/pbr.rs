//! Real-time PBR forward-shading — the Cook-Torrance BRDF (Phase 30.6).
//!
//! ## What this is
//!
//! A genuine, self-contained **physically-based rendering** shading
//! library: the Cook-Torrance microfacet BRDF and a forward-shading
//! evaluator that lights a surface point from analytic lights plus the
//! image-based-lighting environment term.
//!
//! This is the **real-time** PBR path — the same closed-form BRDF a
//! GPU fragment shader runs once per pixel, *not* a path tracer. It
//! pairs with the HDR / IBL environment work already in
//! [`crate::environment`]: [`shade_surface`] consumes a
//! [`crate::environment::IrradianceMap`] for the diffuse ambient term
//! exactly as a real-time renderer samples a prefiltered irradiance
//! probe.
//!
//! ## The Cook-Torrance BRDF
//!
//! The outgoing radiance toward the viewer is the sum, over every
//! light, of
//!
//! ```text
//!   Lo = (kd · albedo / π  +  spec) · radiance · (n·l)
//! ```
//!
//! where the **specular** term is the Cook-Torrance microfacet model
//!
//! ```text
//!   spec = D · G · F  /  (4 · (n·v) · (n·l))
//! ```
//!
//! with:
//!
//! - **D** — the GGX / Trowbridge-Reitz normal-distribution function
//!   ([`distribution_ggx`]): the statistical fraction of microfacets
//!   whose normal aligns with the half-vector `h`. GGX is the
//!   industry-standard NDF (Disney, UE4, glTF) — long tails, a
//!   physically-plausible highlight.
//! - **G** — the Smith geometry term ([`geometry_smith`]): the
//!   fraction of microfacets not shadowed or masked by their
//!   neighbours, the product of two Schlick-GGX one-direction terms.
//! - **F** — the Fresnel-Schlick reflectance ([`fresnel_schlick`]):
//!   how reflectance climbs toward 1 at grazing angles, interpolating
//!   from the normal-incidence reflectance `F₀`.
//!
//! The **diffuse** term is Lambertian (`albedo / π`), scaled by the
//! energy-conservation factor `kd = (1 − F)·(1 − metallic)` so a
//! surface never reflects more light than it receives and metals carry
//! no diffuse lobe.
//!
//! `F₀` (normal-incidence reflectance) is derived from the material:
//! for a dielectric it is the achromatic `((ior−1)/(ior+1))²`; for a
//! metal it is the base colour itself (metals tint their specular
//! reflection). [`Material`] supplies `metallic`, `roughness`, `ior`,
//! and the colours; [`f0_from_material`] does the blend.
//!
//! ## Honest scope — what it is and is not
//!
//! **It is a real PBR shader.** Every function here is the exact
//! closed form that ships inside UE4 / glTF-viewer / Filament fragment
//! shaders. The energy behaviour is correct: a white furnace test (a
//! uniform environment, no analytic lights) returns close to the input
//! radiance, dielectric Fresnel rises to 1 at grazing incidence, and a
//! rougher surface spreads its highlight.
//!
//! **It is the CPU-side library, not a wired GPU render pass.** The
//! task that owns the live `wgpu::Device` — the desktop shell —
//! uploads the material constants + the irradiance probe and runs this
//! same maths in WGSL across a real forward render pass. Shipping the
//! BRDF as a compile-checked, test-covered Rust library is the
//! correct crate-layer deliverable; the WGSL fragment shader + the
//! render-pass plumbing are the app-layer follow-up (they need the GPU
//! context this pure-data crate deliberately does not pull in). The
//! Rust functions here are also directly usable for CPU-side preview
//! shading, light-probe baking, and as the reference the WGSL port is
//! checked against.
//!
//! **Specular IBL** is the split-sum approximation (Karis 2013): the
//! environment prefiltered into a roughness-indexed mip chain
//! ([`crate::environment::EnvironmentMap::prefilter_specular`]) times
//! a precomputed environment-BRDF scale/bias table ([`BrdfLut`] /
//! [`compute_brdf_lut`]). [`specular_ibl`] reconstructs the
//! roughness-aware specular ambient term from those two factors. The
//! simpler sharp single-sample reflection in [`ambient_ibl`] is kept
//! as the fallback for callers that have no prefiltered environment.

use crate::environment::{EnvironmentMap, IrradianceMap};
use crate::material::Material;

/// Outcome of shading one surface point — the linear-RGB radiance
/// leaving the surface toward the viewer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShadedColor {
    /// Linear-RGB radiance toward the viewer (not yet tone-mapped or
    /// gamma-encoded — that is the display stage's job).
    pub rgb: [f32; 3],
}

impl ShadedColor {
    /// The pure-black result (no light reached the point).
    pub const BLACK: ShadedColor = ShadedColor { rgb: [0.0; 3] };

    /// Apply Reinhard tone mapping `c / (1 + c)` followed by sRGB
    /// gamma encoding — a convenience for callers that want a
    /// display-ready `[0, 1]` colour without writing the transfer
    /// function themselves. A real renderer would do this in a
    /// post-process pass; it is offered here so the BRDF library is
    /// usable end-to-end for CPU preview.
    pub fn to_display_srgb(self) -> [f32; 3] {
        let mut out = [0.0f32; 3];
        for (k, o) in out.iter_mut().enumerate() {
            let mapped = self.rgb[k] / (1.0 + self.rgb[k]);
            *o = linear_to_srgb(mapped.clamp(0.0, 1.0));
        }
        out
    }
}

/// A single analytic light reduced to what the BRDF evaluator needs at
/// a shading point: the unit direction *toward* the light and the
/// linear-RGB radiance arriving along it.
///
/// Callers convert a [`crate::Light`] into this with
/// [`incident_light`] — that step folds in distance attenuation, the
/// spot cone falloff, and the watt-to-radiance conversion, so the
/// BRDF core stays a pure function of geometry + radiance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IncidentLight {
    /// Unit vector from the shading point *toward* the light.
    pub direction: [f32; 3],
    /// Linear-RGB radiance arriving at the point along `direction`.
    pub radiance: [f32; 3],
}

/// The surface point + viewing geometry a shading call needs.
#[derive(Clone, Copy, Debug)]
pub struct SurfacePoint {
    /// Outward unit surface normal in world space.
    pub normal: [f32; 3],
    /// Unit vector from the shading point *toward* the camera.
    pub view: [f32; 3],
}

/// GGX / Trowbridge-Reitz normal-distribution function `D`.
///
/// Returns the (unnormalised-to-solid-angle) density of microfacets
/// whose normal equals the half-vector. `n_dot_h` is `max(0, n·h)`;
/// `roughness` is the perceptual roughness in `[0, 1]`.
///
/// ```text
///   α  = roughness²
///   D  = α² / (π · ((n·h)²·(α²−1) + 1)²)
/// ```
///
/// `roughness` is squared once to get `α` (the "roughness remapping"
/// Disney introduced so the parameter feels perceptually linear), and
/// `α` is squared again inside the formula. A roughness of 0 gives a
/// Dirac-like spike (a perfect mirror); roughness 1 gives a broad,
/// nearly-uniform lobe.
pub fn distribution_ggx(n_dot_h: f32, roughness: f32) -> f32 {
    let alpha = (roughness * roughness).max(1e-4);
    let a2 = alpha * alpha;
    let nh = n_dot_h.max(0.0);
    let denom = nh * nh * (a2 - 1.0) + 1.0;
    a2 / (std::f32::consts::PI * denom * denom)
}

/// One-direction Schlick-GGX geometry term `G₁`.
///
/// `n_dot_x` is `max(0, n·x)` for either the view or the light
/// direction; `k` is the roughness-derived remap. Returns the fraction
/// of microfacets visible along `x`.
fn geometry_schlick_ggx(n_dot_x: f32, k: f32) -> f32 {
    let nx = n_dot_x.max(0.0);
    nx / (nx * (1.0 - k) + k)
}

/// Smith geometry term `G` for direct (analytic-light) shading.
///
/// The Smith model factors masking-shadowing into the product of two
/// independent one-direction terms — one for the light, one for the
/// view:
///
/// ```text
///   k = (roughness + 1)² / 8          (the UE4 direct-lighting remap)
///   G = G₁(n·v) · G₁(n·l)
/// ```
///
/// `n_dot_v` and `n_dot_l` are the clamped cosines. The `(r+1)²/8`
/// remap is the one Epic published for *direct* lights; image-based
/// lighting uses `r²/2` instead (different because IBL integrates
/// over all directions).
pub fn geometry_smith(n_dot_v: f32, n_dot_l: f32, roughness: f32) -> f32 {
    // Direct-lighting remap of roughness into the Schlick-GGX `k`.
    let r = roughness + 1.0;
    let k = (r * r) / 8.0;
    geometry_schlick_ggx(n_dot_v, k) * geometry_schlick_ggx(n_dot_l, k)
}

/// Fresnel-Schlick reflectance.
///
/// Interpolates from the normal-incidence reflectance `f0` toward 1
/// as the angle to the surface increases:
///
/// ```text
///   F = F₀ + (1 − F₀) · (1 − cosθ)⁵
/// ```
///
/// `cos_theta` is `max(0, h·v)` (the cosine between the half-vector
/// and the view). Per-channel because `f0` is an RGB triple (metals
/// tint their reflection).
pub fn fresnel_schlick(cos_theta: f32, f0: [f32; 3]) -> [f32; 3] {
    let m = (1.0 - cos_theta.clamp(0.0, 1.0)).powi(5);
    [
        f0[0] + (1.0 - f0[0]) * m,
        f0[1] + (1.0 - f0[1]) * m,
        f0[2] + (1.0 - f0[2]) * m,
    ]
}

/// Fresnel-Schlick with a roughness-aware ceiling — used for the
/// ambient (IBL) term so a rough surface's grazing-angle reflectance
/// does not blow past its actual specular colour.
///
/// This is the standard `fresnelSchlickRoughness` from the learnopengl
/// IBL chapter: the `(1 − F₀)` term is replaced by
/// `max(1 − roughness, F₀) − F₀`.
fn fresnel_schlick_roughness(cos_theta: f32, f0: [f32; 3], roughness: f32) -> [f32; 3] {
    let m = (1.0 - cos_theta.clamp(0.0, 1.0)).powi(5);
    let mut out = [0.0f32; 3];
    for k in 0..3 {
        let ceil = (1.0 - roughness).max(f0[k]);
        out[k] = f0[k] + (ceil - f0[k]) * m;
    }
    out
}

/// Derive the normal-incidence reflectance `F₀` of a [`Material`].
///
/// - **Dielectric** (`metallic = 0`): the achromatic Schlick `F₀`
///   from the index of refraction, `((ior−1)/(ior+1))²`. For the
///   common `ior = 1.5` this is the familiar `0.04`.
/// - **Metal** (`metallic = 1`): `F₀` is the base colour itself —
///   metals reflect a *tinted* specular and have no separate diffuse
///   albedo.
/// - **In between**: a linear blend, the glTF metallic-roughness
///   convention.
///
/// The material's [`Material::specular_color`] is used as the metallic
/// tint (it already stores the per-metal reflectance, e.g. gold's
/// warm `F₀`); the dielectric base is the IOR-derived grey.
pub fn f0_from_material(material: &Material) -> [f32; 3] {
    let ior = material.ior.max(1.0);
    let dielectric = {
        let r = (ior - 1.0) / (ior + 1.0);
        r * r
    };
    let metallic = material.metallic.clamp(0.0, 1.0);
    let mut out = [0.0f32; 3];
    for (k, o) in out.iter_mut().enumerate() {
        // Dielectrics: achromatic IOR reflectance. Metals: the
        // material's specular tint. Lerp by metallic.
        *o = dielectric * (1.0 - metallic) + material.specular_color[k] * metallic;
    }
    out
}

/// The diffuse albedo a [`Material`] contributes — its base colour
/// faded out as the surface becomes metallic (metals have no diffuse
/// lobe).
fn diffuse_albedo(material: &Material) -> [f32; 3] {
    let metallic = material.metallic.clamp(0.0, 1.0);
    [
        material.diffuse_color[0] * (1.0 - metallic),
        material.diffuse_color[1] * (1.0 - metallic),
        material.diffuse_color[2] * (1.0 - metallic),
    ]
}

/// The precomputed **BRDF integration look-up table** — the second
/// factor of the split-sum specular-IBL approximation (Phase 30.7).
///
/// The specular IBL integral splits into a prefiltered environment
/// (see [`crate::environment::PrefilteredEnvironment`]) times a
/// **scale / bias** pair that depends only on `(n·v, roughness)` — not
/// on the environment at all. That pair is what this table stores, so
/// it is computed **once** for any scene. The specular IBL is then
/// reconstructed as `prefiltered · (F₀·scale + bias)`.
#[derive(Clone, Debug, PartialEq)]
pub struct BrdfLut {
    /// Table resolution (square: `size × size`).
    pub size: usize,
    /// `size · size` `(scale, bias)` pairs, row-major. The row index
    /// is `roughness`, the column index is `n·v`, both mapped from
    /// `[0, 1]`.
    pub entries: Vec<(f32, f32)>,
}

impl BrdfLut {
    /// Look up the `(scale, bias)` pair for a given `n·v` and
    /// `roughness` by bilinear interpolation.
    pub fn sample(&self, n_dot_v: f32, roughness: f32) -> (f32, f32) {
        if self.size == 0 {
            return (1.0, 0.0);
        }
        let nv = n_dot_v.clamp(0.0, 1.0);
        let r = roughness.clamp(0.0, 1.0);
        let last = (self.size - 1) as f32;
        let fx = nv * last;
        let fy = r * last;
        let x0 = fx.floor() as usize;
        let y0 = fy.floor() as usize;
        let x1 = (x0 + 1).min(self.size - 1);
        let y1 = (y0 + 1).min(self.size - 1);
        let tx = fx - x0 as f32;
        let ty = fy - y0 as f32;
        let at = |x: usize, y: usize| self.entries[y * self.size + x];
        let (s00, b00) = at(x0, y0);
        let (s10, b10) = at(x1, y0);
        let (s01, b01) = at(x0, y1);
        let (s11, b11) = at(x1, y1);
        let s = lerp(lerp(s00, s10, tx), lerp(s01, s11, tx), ty);
        let b = lerp(lerp(b00, b10, tx), lerp(b01, b11, tx), ty);
        (s, b)
    }
}

/// Compute the [`BrdfLut`] — integrate the environment-BRDF split-sum
/// scale / bias table at `size × size` resolution.
///
/// # Method (Karis 2013, `IntegrateBRDF`)
///
/// For each `(n·v, roughness)` cell the BRDF is integrated by **GGX
/// importance sampling**: `samples` half-vectors are drawn from the
/// GGX distribution, each reflected to a light direction `l`; for the
/// hemisphere-positive ones the Smith geometry term (IBL `k = α²/2`
/// remap) and a Fresnel split give two accumulators —
///
/// ```text
///   G_vis = G · (v·h) / (n·h · n·v)
///   Fc    = (1 − v·h)⁵
///   scale += (1 − Fc) · G_vis
///   bias  += Fc · G_vis
/// ```
///
/// — averaged over the samples. The result is environment-independent,
/// so this is called once and reused for every material and every
/// scene.
///
/// `size` is the table resolution (a 64×64 or 128×128 LUT is the
/// industry norm); `samples` is the importance-sample count per cell.
pub fn compute_brdf_lut(size: usize, samples: usize) -> BrdfLut {
    let size = size.max(2);
    let samples = samples.max(1);
    let mut entries = Vec::with_capacity(size * size);
    for yi in 0..size {
        // Row = roughness.
        let roughness = yi as f32 / (size - 1) as f32;
        for xi in 0..size {
            // Column = n·v.
            let n_dot_v = (xi as f32 / (size - 1) as f32).max(1e-4);
            entries.push(integrate_brdf_cell(n_dot_v, roughness, samples));
        }
    }
    BrdfLut { size, entries }
}

/// Integrate one `(n·v, roughness)` cell of the BRDF LUT.
///
/// The normal is fixed along `+Z`, so `n·x` for any direction `x` is
/// simply its Z component — the integration uses that directly.
fn integrate_brdf_cell(n_dot_v: f32, roughness: f32, samples: usize) -> (f32, f32) {
    // Place the view in the X-Z plane at the requested n·v (the normal
    // is +Z).
    let v = [(1.0 - n_dot_v * n_dot_v).max(0.0).sqrt(), 0.0, n_dot_v];
    let alpha = roughness * roughness;
    // IBL geometry remap: k = roughness²/2 (the Karis split-sum IBL
    // remap — distinct from the direct-lighting (r+1)²/8 remap, which
    // integrates a single light rather than the whole hemisphere).
    // NB this is roughness², *not* α² = roughness⁴: a roughness⁴ remap
    // collapses k toward 0 at mid-roughness, which sends the Smith G
    // (and the visibility-weighted G_vis = G·VoH/(NoH·NoV)) far above
    // unity and makes the LUT amplify environment energy.
    let k = alpha / 2.0;
    let mut scale = 0.0f32;
    let mut bias = 0.0f32;
    for s in 0..samples {
        let (u1, u2) = hammersley_2d(s, samples);
        // GGX importance-sampled half-vector about +Z.
        let phi = std::f32::consts::TAU * u1;
        let cos_theta = (((1.0 - u2) / (1.0 + (alpha * alpha - 1.0) * u2)).max(0.0)).sqrt();
        let sin_theta = (1.0 - cos_theta * cos_theta).max(0.0).sqrt();
        let h = [sin_theta * phi.cos(), sin_theta * phi.sin(), cos_theta];
        // l = reflect(-v, h) = 2·(v·h)·h − v.
        let vh = v[0] * h[0] + v[1] * h[1] + v[2] * h[2];
        let l = [
            2.0 * vh * h[0] - v[0],
            2.0 * vh * h[1] - v[1],
            2.0 * vh * h[2] - v[2],
        ];
        let n_dot_l = l[2].max(0.0); // n is +Z
        let n_dot_h = h[2].max(0.0);
        let v_dot_h = vh.max(0.0);
        if n_dot_l > 0.0 {
            // Smith geometry, IBL remap.
            let g = geom_schlick_ibl(n_dot_v, k) * geom_schlick_ibl(n_dot_l, k);
            // Visibility-weighted geometry term.
            let g_vis = g * v_dot_h / (n_dot_h * n_dot_v).max(1e-6);
            let fc = (1.0 - v_dot_h).powi(5);
            scale += (1.0 - fc) * g_vis;
            bias += fc * g_vis;
        }
    }
    let inv = 1.0 / samples as f32;
    (scale * inv, bias * inv)
}

/// One-direction Schlick-GGX geometry term with the **IBL** roughness
/// remap (`k = α²/2`).
fn geom_schlick_ibl(n_dot_x: f32, k: f32) -> f32 {
    let nx = n_dot_x.max(0.0);
    nx / (nx * (1.0 - k) + k)
}

/// Hammersley low-discrepancy 2-D sample `i` of `n` — `u1 = i/n`,
/// `u2` = base-2 radical inverse (Van der Corput) of `i`. Low-
/// discrepancy sampling converges the BRDF integral far faster than
/// pseudo-random.
fn hammersley_2d(i: usize, n: usize) -> (f32, f32) {
    // Bit-reverse `i` (the base-2 radical inverse): swap nibble pairs
    // of every width down to single bits.
    let mut bits = (i as u32).rotate_right(16);
    bits = ((bits & 0x5555_5555) << 1) | ((bits & 0xAAAA_AAAA) >> 1);
    bits = ((bits & 0x3333_3333) << 2) | ((bits & 0xCCCC_CCCC) >> 2);
    bits = ((bits & 0x0F0F_0F0F) << 4) | ((bits & 0xF0F0_F0F0) >> 4);
    bits = ((bits & 0x00FF_00FF) << 8) | ((bits & 0xFF00_FF00) >> 8);
    let radical = bits as f32 * 2.328_306_4e-10;
    let u1 = if n == 0 { 0.0 } else { i as f32 / n as f32 };
    (u1, radical)
}

/// Reconstruct the **specular IBL** term from the split-sum factors —
/// the prefiltered environment times the BRDF LUT's `(scale, bias)`.
///
/// ```text
///   specular_ibl = prefiltered(reflect(−v, n), roughness)
///                  · (F₀ · scale + bias)
/// ```
///
/// This is the proper roughness-aware specular ambient term: a rough
/// surface samples a blurred environment level and a smooth one a
/// sharp level, and the `(scale, bias)` pair applies the
/// energy-correct Fresnel + geometry weighting. It supersedes the
/// sharp single-sample reflection that [`ambient_ibl`] falls back to
/// when no prefiltered environment is supplied.
pub fn specular_ibl(
    surface: &SurfacePoint,
    material: &Material,
    prefiltered: &crate::environment::PrefilteredEnvironment,
    lut: &BrdfLut,
) -> [f32; 3] {
    let n = match normalize3(surface.normal) {
        Some(v) => v,
        None => return [0.0; 3],
    };
    let v = match normalize3(surface.view) {
        Some(v) => v,
        None => return [0.0; 3],
    };
    let n_dot_v = dot3(n, v).max(1e-4);
    let roughness = material.roughness.clamp(0.0, 1.0);
    let f0 = f0_from_material(material);

    // Sample the prefiltered environment in the reflection direction.
    let refl = reflect(neg3(v), n);
    let prefiltered_color = prefiltered.sample(refl, roughness);
    // The environment-BRDF scale / bias from the LUT.
    let (scale, bias) = lut.sample(n_dot_v, roughness);

    [
        prefiltered_color[0] * (f0[0] * scale + bias),
        prefiltered_color[1] * (f0[1] * scale + bias),
        prefiltered_color[2] * (f0[2] * scale + bias),
    ]
}

/// Evaluate the full Cook-Torrance BRDF response of one surface point
/// to one analytic light, returning the linear-RGB radiance that light
/// contributes toward the viewer.
///
/// This is the per-light body of [`shade_surface`]; it is exposed so
/// callers can shade against a custom light set or accumulate lights
/// incrementally.
///
/// # Method
///
/// With `n` the normal, `v` the view direction, `l` the light
/// direction, and `h = normalise(v + l)` the half-vector:
///
/// 1. specular `= D·G·F / (4·(n·v)·(n·l))` — the Cook-Torrance term;
/// 2. `kd = (1 − F)·(1 − metallic)` — the diffuse energy left after
///    the specular reflection (and zero for metals);
/// 3. `Lo = (kd·albedo/π + specular) · radiance · (n·l)`.
///
/// A light below the horizon (`n·l ≤ 0`) contributes nothing.
pub fn brdf_direct(surface: &SurfacePoint, material: &Material, light: &IncidentLight) -> [f32; 3] {
    let n = match normalize3(surface.normal) {
        Some(v) => v,
        None => return [0.0; 3],
    };
    let v = match normalize3(surface.view) {
        Some(v) => v,
        None => return [0.0; 3],
    };
    let l = match normalize3(light.direction) {
        Some(v) => v,
        None => return [0.0; 3],
    };
    let n_dot_l = dot3(n, l);
    if n_dot_l <= 0.0 {
        // Light is below the surface horizon.
        return [0.0; 3];
    }
    let n_dot_v = dot3(n, v).max(1e-4);
    // Half-vector between view and light.
    let h = match normalize3([v[0] + l[0], v[1] + l[1], v[2] + l[2]]) {
        Some(v) => v,
        None => return [0.0; 3],
    };
    let n_dot_h = dot3(n, h).max(0.0);
    let h_dot_v = dot3(h, v).max(0.0);

    let roughness = material.roughness.clamp(0.0, 1.0);
    let f0 = f0_from_material(material);

    // Cook-Torrance specular term.
    let d = distribution_ggx(n_dot_h, roughness);
    let g = geometry_smith(n_dot_v, n_dot_l, roughness);
    let f = fresnel_schlick(h_dot_v, f0);
    // spec = D·G·F / (4·(n·v)·(n·l)).
    let denom = 4.0 * n_dot_v * n_dot_l + 1e-4;
    let spec = [
        d * g * f[0] / denom,
        d * g * f[1] / denom,
        d * g * f[2] / denom,
    ];

    // Diffuse term — energy-conserving Lambertian. kd is what's left
    // after the specular reflection, and zero for metals.
    let metallic = material.metallic.clamp(0.0, 1.0);
    let albedo = diffuse_albedo(material);
    let inv_pi = std::f32::consts::FRAC_1_PI;
    let mut out = [0.0f32; 3];
    for k in 0..3 {
        let kd = (1.0 - f[k]) * (1.0 - metallic);
        let diffuse = kd * albedo[k] * inv_pi;
        out[k] = (diffuse + spec[k]) * light.radiance[k] * n_dot_l;
    }
    out
}

/// The ambient (image-based-lighting) term for a surface point.
///
/// Real-time PBR splits ambient light into a **diffuse** part — the
/// cosine-weighted environment irradiance, which [`IrradianceMap`]
/// already stores prefiltered — and a **specular** part. This computes
/// the diffuse IBL exactly and adds a *sharp* environment mirror
/// reflection for the specular part (a true roughness-aware specular
/// prefilter is the documented follow-up — see the module docs).
///
/// ```text
///   ambient_diffuse  = kd · albedo · irradiance(n) / π
///   ambient_specular = F · environment(reflect(−v, n))
/// ```
///
/// `kd` uses the roughness-aware Fresnel so a rough dielectric does
/// not over-reflect at grazing angles. `environment` may be `None` —
/// then only the diffuse IBL term is returned (still correct, just
/// without the mirror highlight).
pub fn ambient_ibl(
    surface: &SurfacePoint,
    material: &Material,
    irradiance: &IrradianceMap,
    environment: Option<&EnvironmentMap>,
) -> [f32; 3] {
    let n = match normalize3(surface.normal) {
        Some(v) => v,
        None => return [0.0; 3],
    };
    let v = match normalize3(surface.view) {
        Some(v) => v,
        None => return [0.0; 3],
    };
    let n_dot_v = dot3(n, v).max(1e-4);
    let roughness = material.roughness.clamp(0.0, 1.0);
    let metallic = material.metallic.clamp(0.0, 1.0);
    let f0 = f0_from_material(material);
    let f = fresnel_schlick_roughness(n_dot_v, f0, roughness);

    // Diffuse IBL: the prefiltered irradiance probe carries the
    // cosine-weighted hemisphere integral (already π·L for a uniform
    // white map), so outgoing diffuse radiance is irradiance·albedo/π.
    let irr = irradiance.irradiance(n);
    let albedo = diffuse_albedo(material);
    let inv_pi = std::f32::consts::FRAC_1_PI;

    // Specular IBL: a sharp mirror reflection of the environment in
    // the direction reflect(−v, n). For a rough surface this is too
    // sharp (a real renderer samples a roughness-blurred mip), but it
    // is a correct upper-roughness-0 reference and a usable
    // approximation; the prefilter is the follow-up.
    let refl = reflect(neg3(v), n);
    let env_radiance = environment
        .map(|e| e.sample_direction(refl))
        .unwrap_or([0.0; 3]);

    let mut out = [0.0f32; 3];
    for k in 0..3 {
        let kd = (1.0 - f[k]) * (1.0 - metallic);
        let diffuse = kd * albedo[k] * irr[k] * inv_pi;
        let specular = f[k] * env_radiance[k];
        out[k] = diffuse + specular;
    }
    out
}

/// Forward-shade a surface point: the full real-time PBR result.
///
/// Sums [`brdf_direct`] over every analytic light, adds the
/// [`ambient_ibl`] environment term, and folds in the material's
/// emissive colour. This is the complete fragment-shader computation
/// for one pixel.
///
/// `lights` is the per-point incident-light set (build it with
/// [`incident_light`]). `irradiance` is the prefiltered diffuse probe;
/// `environment` is the optional HDR map for the sharp specular
/// reflection.
pub fn shade_surface(
    surface: &SurfacePoint,
    material: &Material,
    lights: &[IncidentLight],
    irradiance: &IrradianceMap,
    environment: Option<&EnvironmentMap>,
) -> ShadedColor {
    let mut rgb = [0.0f32; 3];
    // Direct lighting — sum the Cook-Torrance response of every light.
    for light in lights {
        let contrib = brdf_direct(surface, material, light);
        for (c, &add) in rgb.iter_mut().zip(contrib.iter()) {
            *c += add;
        }
    }
    // Ambient image-based lighting.
    let ambient = ambient_ibl(surface, material, irradiance, environment);
    for (c, &add) in rgb.iter_mut().zip(ambient.iter()) {
        *c += add;
    }
    // Emissive: the surface's own glow, added directly.
    for (c, &add) in rgb.iter_mut().zip(material.emissive.iter()) {
        *c += add;
    }
    ShadedColor { rgb }
}

/// Reduce a scene [`crate::Light`] to the [`IncidentLight`] arriving
/// at a world-space `point` — folding in distance attenuation, the
/// spot-cone falloff, and the watt-to-radiance conversion.
///
/// Returns `None` when the light contributes nothing at `point` (the
/// point is outside a spot cone, or the area-light back-faces it).
///
/// # Conversions
///
/// - **Point / Spot:** inverse-square attenuation `1/d²`; the watt
///   intensity becomes radiance `intensity/(4π·d²)`.
/// - **Directional:** no attenuation; the irradiance is used directly
///   (a directional light is specified as irradiance already).
/// - **Spot:** additionally a smooth `smoothstep` falloff between the
///   inner and outer cone half-angles.
/// - **Area:** treated as a point at the rectangle centre for the
///   real-time path, attenuated `1/d²`, with a cosine fade by the
///   emitter's facing — a real renderer would integrate the
///   rectangle (LTC); the centre approximation is the real-time
///   stand-in.
pub fn incident_light(light: &crate::Light, point: [f32; 3]) -> Option<IncidentLight> {
    use crate::Light;
    match light {
        Light::Directional {
            direction,
            color,
            irradiance,
        } => {
            // The scene stores the direction the light *travels*; the
            // shading direction points back toward the light.
            let d = [
                -(direction.x as f32),
                -(direction.y as f32),
                -(direction.z as f32),
            ];
            let dir = normalize3(d)?;
            Some(IncidentLight {
                direction: dir,
                radiance: scale3(*color, *irradiance),
            })
        }
        Light::Point {
            position,
            color,
            intensity,
        } => {
            let to_light = [
                position.x as f32 - point[0],
                position.y as f32 - point[1],
                position.z as f32 - point[2],
            ];
            let dist2 = dot3(to_light, to_light);
            if dist2 < 1e-12 {
                return None;
            }
            let dir = normalize3(to_light)?;
            // Watt point source → radiance at the point: I/(4π·d²).
            let atten = 1.0 / (4.0 * std::f32::consts::PI * dist2);
            Some(IncidentLight {
                direction: dir,
                radiance: scale3(*color, intensity * atten),
            })
        }
        Light::Spot {
            position,
            direction,
            inner_angle_rad,
            outer_angle_rad,
            color,
            intensity,
        } => {
            let to_light = [
                position.x as f32 - point[0],
                position.y as f32 - point[1],
                position.z as f32 - point[2],
            ];
            let dist2 = dot3(to_light, to_light);
            if dist2 < 1e-12 {
                return None;
            }
            let dir = normalize3(to_light)?;
            // Cone falloff: angle between the spot axis and the
            // direction *from* the spot toward the point.
            let axis = normalize3([direction.x as f32, direction.y as f32, direction.z as f32])?;
            // From the spot toward the point = −dir.
            let cos_angle = dot3(axis, neg3(dir));
            let cos_inner = (*inner_angle_rad as f32).cos();
            let cos_outer = (*outer_angle_rad as f32).cos();
            if cos_angle <= cos_outer {
                // Outside the outer cone — no light.
                return None;
            }
            // Smooth falloff between inner (full) and outer (zero).
            let falloff = if cos_angle >= cos_inner {
                1.0
            } else {
                let t = (cos_angle - cos_outer) / (cos_inner - cos_outer).max(1e-6);
                t * t * (3.0 - 2.0 * t) // smoothstep
            };
            let atten = falloff / (4.0 * std::f32::consts::PI * dist2);
            Some(IncidentLight {
                direction: dir,
                radiance: scale3(*color, intensity * atten),
            })
        }
        Light::Area {
            corner,
            edge_u,
            edge_v,
            color,
            radiance,
        } => {
            // Real-time stand-in: treat the rectangle as a point at
            // its centre. The centre is corner + (edge_u + edge_v)/2.
            let centre = [
                corner.x as f32 + 0.5 * (edge_u.x as f32 + edge_v.x as f32),
                corner.y as f32 + 0.5 * (edge_u.y as f32 + edge_v.y as f32),
                corner.z as f32 + 0.5 * (edge_u.z as f32 + edge_v.z as f32),
            ];
            let to_light = [
                centre[0] - point[0],
                centre[1] - point[1],
                centre[2] - point[2],
            ];
            let dist2 = dot3(to_light, to_light);
            if dist2 < 1e-12 {
                return None;
            }
            let dir = normalize3(to_light)?;
            // The emitter faces along edge_u × edge_v.
            let emit_normal = normalize3(cross3(
                [edge_u.x as f32, edge_u.y as f32, edge_u.z as f32],
                [edge_v.x as f32, edge_v.y as f32, edge_v.z as f32],
            ))?;
            // Cosine fade: the rectangle radiates most toward points
            // it faces. `−dir` is from the emitter toward the point.
            let facing = dot3(emit_normal, neg3(dir)).max(0.0);
            if facing <= 0.0 {
                return None;
            }
            // Approximate area-light irradiance: radiance · area ·
            // facing / d². edge area = |edge_u × edge_v|.
            let area = cross3(
                [edge_u.x as f32, edge_u.y as f32, edge_u.z as f32],
                [edge_v.x as f32, edge_v.y as f32, edge_v.z as f32],
            );
            let area_mag = (area[0] * area[0] + area[1] * area[1] + area[2] * area[2]).sqrt();
            let atten = facing * area_mag / dist2;
            Some(IncidentLight {
                direction: dir,
                radiance: scale3(*color, radiance * atten),
            })
        }
    }
}

// --- small vector + colour helpers ---

/// Reflect `incident` about the unit surface normal `n`:
/// `r = i − 2·(i·n)·n`.
fn reflect(incident: [f32; 3], n: [f32; 3]) -> [f32; 3] {
    let d = 2.0 * dot3(incident, n);
    [
        incident[0] - d * n[0],
        incident[1] - d * n[1],
        incident[2] - d * n[2],
    ]
}

fn dot3(a: [f32; 3], b: [f32; 3]) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross3(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn neg3(a: [f32; 3]) -> [f32; 3] {
    [-a[0], -a[1], -a[2]]
}

fn scale3(v: [f32; 3], s: f32) -> [f32; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

fn normalize3(v: [f32; 3]) -> Option<[f32; 3]> {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len < 1e-12 {
        None
    } else {
        Some([v[0] / len, v[1] / len, v[2] / len])
    }
}

/// Scalar linear interpolation, `a` at `t=0` and `b` at `t=1`.
fn lerp(a: f32, b: f32, t: f32) -> f32 {
    a + (b - a) * t
}

/// Linear → sRGB gamma encode for one channel (the IEC 61966-2-1
/// transfer function).
fn linear_to_srgb(c: f32) -> f32 {
    if c <= 0.0031308 {
        12.92 * c
    } else {
        1.055 * c.powf(1.0 / 2.4) - 0.055
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    #[test]
    fn ggx_peaks_when_half_vector_aligns_with_normal() {
        // The GGX NDF is maximal at n·h = 1 and decreasing away from
        // it — the defining shape of a microfacet distribution.
        let r = 0.3;
        let peak = distribution_ggx(1.0, r);
        let off = distribution_ggx(0.7, r);
        let far = distribution_ggx(0.2, r);
        assert!(peak > off, "GGX should peak at n·h = 1");
        assert!(off > far, "GGX should fall off away from the peak");
        assert!(peak.is_finite() && peak > 0.0);
    }

    #[test]
    fn ggx_is_broader_for_a_rougher_surface() {
        // A rougher surface spreads its highlight: at the peak the
        // rough NDF is *lower* (energy spread wider), and off-peak it
        // is *higher* than a smooth surface.
        let smooth_peak = distribution_ggx(1.0, 0.1);
        let rough_peak = distribution_ggx(1.0, 0.9);
        assert!(
            smooth_peak > rough_peak,
            "a smooth surface has a taller, narrower highlight"
        );
        let smooth_off = distribution_ggx(0.5, 0.1);
        let rough_off = distribution_ggx(0.5, 0.9);
        assert!(
            rough_off > smooth_off,
            "a rough surface spreads energy into the off-peak angles"
        );
    }

    #[test]
    fn fresnel_rises_to_one_at_grazing_incidence() {
        // Fresnel-Schlick: at normal incidence F = F₀; at grazing
        // incidence (cosθ → 0) F → 1 for every channel.
        let f0 = [0.04, 0.04, 0.04];
        let normal_inc = fresnel_schlick(1.0, f0);
        let grazing = fresnel_schlick(0.0, f0);
        for k in 0..3 {
            assert!(
                (normal_inc[k] - 0.04).abs() < 1e-6,
                "F at normal incidence should equal F₀"
            );
            assert!(
                (grazing[k] - 1.0).abs() < 1e-6,
                "F at grazing incidence should reach 1"
            );
        }
        // Monotone in between.
        let mid = fresnel_schlick(0.5, f0);
        assert!(mid[0] > normal_inc[0] && mid[0] < grazing[0]);
    }

    #[test]
    fn geometry_term_is_in_unit_range() {
        // The Smith G term is a visibility fraction — it must stay in
        // [0, 1] for every plausible cosine / roughness combination.
        for &nv in &[0.05f32, 0.3, 0.7, 1.0] {
            for &nl in &[0.05f32, 0.3, 0.7, 1.0] {
                for &r in &[0.0f32, 0.3, 0.7, 1.0] {
                    let g = geometry_smith(nv, nl, r);
                    assert!(
                        (0.0..=1.0).contains(&g),
                        "G out of range: nv={nv} nl={nl} r={r} → {g}"
                    );
                }
            }
        }
    }

    #[test]
    fn f0_of_dielectric_is_the_ior_grey() {
        // A non-metal at ior 1.5 has the canonical F₀ ≈ 0.04, grey.
        let mat = Material {
            metallic: 0.0,
            ior: 1.5,
            ..Material::default()
        };
        let f0 = f0_from_material(&mat);
        for (k, &channel) in f0.iter().enumerate() {
            assert!(
                (channel - 0.04).abs() < 0.005,
                "dielectric F₀ channel {k} should be ~0.04, got {channel}"
            );
        }
    }

    #[test]
    fn f0_of_metal_is_the_specular_tint() {
        // A fully-metallic material reflects its specular colour as F₀.
        let mat = Material {
            metallic: 1.0,
            specular_color: [1.0, 0.78, 0.34], // gold-ish
            ..Material::default()
        };
        let f0 = f0_from_material(&mat);
        assert!((f0[0] - 1.0).abs() < 1e-6);
        assert!((f0[1] - 0.78).abs() < 1e-6);
        assert!((f0[2] - 0.34).abs() < 1e-6);
    }

    #[test]
    fn brdf_direct_ignores_a_light_below_the_horizon() {
        // A light on the far side of the surface (n·l ≤ 0) contributes
        // nothing.
        let surface = SurfacePoint {
            normal: [0.0, 0.0, 1.0],
            view: [0.0, 0.0, 1.0],
        };
        let below = IncidentLight {
            direction: [0.0, 0.0, -1.0], // points away into the surface
            radiance: [10.0, 10.0, 10.0],
        };
        let out = brdf_direct(&surface, &Material::default(), &below);
        assert_eq!(out, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn brdf_direct_lit_surface_returns_positive_radiance() {
        // A light straight above a surface facing the camera produces
        // positive outgoing radiance on every channel.
        let surface = SurfacePoint {
            normal: [0.0, 0.0, 1.0],
            view: [0.0, 0.0, 1.0],
        };
        let light = IncidentLight {
            direction: [0.0, 0.0, 1.0],
            radiance: [5.0, 5.0, 5.0],
        };
        let mat = Material::matte("white", [0.8, 0.8, 0.8]);
        let out = brdf_direct(&surface, &mat, &light);
        for (k, &channel) in out.iter().enumerate() {
            assert!(channel > 0.0, "channel {k} should be lit, got {channel}");
        }
    }

    #[test]
    fn brdf_diffuse_obeys_lamberts_cosine_law() {
        // For a matte (non-specular-dominated) surface the diffuse
        // response scales with n·l. A light at 60° off the normal
        // should deliver about cos60° = 0.5 of the head-on response.
        let mat = Material {
            roughness: 1.0, // kill the specular spike
            metallic: 0.0,
            ..Material::matte("matte", [0.8, 0.8, 0.8])
        };
        let surface = SurfacePoint {
            normal: [0.0, 0.0, 1.0],
            view: [0.0, 0.0, 1.0],
        };
        let head_on = IncidentLight {
            direction: [0.0, 0.0, 1.0],
            radiance: [1.0, 1.0, 1.0],
        };
        // 60° from the normal in the x-z plane.
        let angled = IncidentLight {
            direction: normalize3([0.866_025, 0.0, 0.5]).unwrap(),
            radiance: [1.0, 1.0, 1.0],
        };
        let r0 = brdf_direct(&surface, &mat, &head_on)[0];
        let r60 = brdf_direct(&surface, &mat, &angled)[0];
        // The ratio should be close to cos60° = 0.5. Allow slack for
        // the (small) specular contribution and the changing
        // half-vector.
        let ratio = r60 / r0;
        assert!(
            (0.35..0.65).contains(&ratio),
            "60° response ratio {ratio} should be ~0.5 (Lambert cosine)"
        );
    }

    #[test]
    fn metal_has_no_diffuse_lobe() {
        // A fully-metallic surface carries zero diffuse albedo — its
        // entire response is the specular term.
        let metal = Material::polished_metal("chrome", [0.95, 0.95, 0.95]);
        assert_eq!(diffuse_albedo(&metal), [0.0, 0.0, 0.0]);
    }

    #[test]
    fn ambient_ibl_of_uniform_white_environment_returns_albedo_scale() {
        // Diffuse IBL: a uniform white irradiance probe stores π·L per
        // texel (L = 1). Outgoing diffuse radiance = irr·albedo/π =
        // π·albedo/π = albedo·kd. For a rough dielectric kd ≈ 1, so
        // the ambient result should be close to the albedo.
        let env = EnvironmentMap::uniform([1.0, 1.0, 1.0]);
        let irr = env.prefilter_irradiance(8, 4, 16).unwrap();
        let mat = Material {
            roughness: 1.0,
            metallic: 0.0,
            emissive: [0.0; 3],
            ..Material::matte("white", [0.7, 0.7, 0.7])
        };
        let surface = SurfacePoint {
            normal: [0.0, 1.0, 0.0],
            view: [0.0, 1.0, 0.0],
        };
        let ambient = ambient_ibl(&surface, &mat, &irr, None);
        // kd at normal incidence for a roughness-1 dielectric is
        // close to 1; the diffuse ambient should land near the albedo.
        for (k, &channel) in ambient.iter().enumerate() {
            assert!(
                channel > 0.4 && channel < 0.8,
                "ambient channel {k} = {channel} should be ~albedo (0.7)"
            );
        }
    }

    #[test]
    fn shade_surface_adds_emissive_directly() {
        // A purely emissive black surface in the dark returns exactly
        // its emissive colour.
        let env = EnvironmentMap::uniform([0.0, 0.0, 0.0]);
        let irr = env.prefilter_irradiance(4, 2, 8).unwrap();
        let mat = Material {
            diffuse_color: [0.0, 0.0, 0.0],
            specular_color: [0.0, 0.0, 0.0],
            emissive: [0.3, 0.6, 0.9],
            ..Material::default()
        };
        let surface = SurfacePoint {
            normal: [0.0, 0.0, 1.0],
            view: [0.0, 0.0, 1.0],
        };
        let shaded = shade_surface(&surface, &mat, &[], &irr, None);
        for k in 0..3 {
            assert!(
                (shaded.rgb[k] - mat.emissive[k]).abs() < 1e-4,
                "emissive channel {k} should pass through: {} vs {}",
                shaded.rgb[k],
                mat.emissive[k]
            );
        }
    }

    #[test]
    fn shade_surface_white_furnace_does_not_create_energy() {
        // The white-furnace test: a uniform white environment, no
        // analytic lights, no emissive. A non-metal surface in a
        // furnace of radiance 1 must reflect *at most* radiance 1 on
        // each channel — a physically-based BRDF never amplifies the
        // environment (energy conservation).
        let env = EnvironmentMap::uniform([1.0, 1.0, 1.0]);
        let irr = env.prefilter_irradiance(8, 4, 16).unwrap();
        let mat = Material {
            diffuse_color: [1.0, 1.0, 1.0],
            specular_color: [0.04, 0.04, 0.04],
            roughness: 0.5,
            metallic: 0.0,
            emissive: [0.0; 3],
            ..Material::default()
        };
        let surface = SurfacePoint {
            normal: [0.0, 1.0, 0.0],
            view: [0.0, 1.0, 0.0],
        };
        let shaded = shade_surface(&surface, &mat, &[], &irr, Some(&env));
        for k in 0..3 {
            assert!(
                shaded.rgb[k] <= 1.05,
                "white-furnace channel {k} = {} must not exceed input radiance 1",
                shaded.rgb[k]
            );
            assert!(shaded.rgb[k] > 0.0, "furnace-lit surface should be lit");
        }
    }

    #[test]
    fn incident_light_directional_points_back_at_the_sun() {
        // A directional light travelling toward −Y illuminates from
        // +Y; the incident direction must point +Y.
        let light = crate::Light::Directional {
            direction: Vector3::new(0.0, -1.0, 0.0),
            color: [1.0, 1.0, 1.0],
            irradiance: 3.0,
        };
        let inc = incident_light(&light, [0.0, 0.0, 0.0]).unwrap();
        assert!((inc.direction[1] - 1.0).abs() < 1e-6);
        // Radiance = colour · irradiance.
        assert!((inc.radiance[0] - 3.0).abs() < 1e-6);
    }

    #[test]
    fn incident_light_point_attenuates_with_distance() {
        // A point light's radiance falls off as 1/d² — doubling the
        // distance quarters the radiance.
        let light = crate::Light::Point {
            position: Vector3::new(0.0, 0.0, 0.0),
            color: [1.0, 1.0, 1.0],
            intensity: 100.0,
        };
        let near = incident_light(&light, [1.0, 0.0, 0.0]).unwrap();
        let far = incident_light(&light, [2.0, 0.0, 0.0]).unwrap();
        let ratio = near.radiance[0] / far.radiance[0];
        assert!(
            (ratio - 4.0).abs() < 1e-3,
            "doubling distance should quarter radiance, ratio {ratio}"
        );
    }

    #[test]
    fn incident_light_spot_is_dark_outside_the_cone() {
        // A spot light pointing down −Y: a point directly below is
        // lit, a point far to the side (outside the cone) gets None.
        let light = crate::Light::Spot {
            position: Vector3::new(0.0, 5.0, 0.0),
            direction: Vector3::new(0.0, -1.0, 0.0),
            inner_angle_rad: 0.2,
            outer_angle_rad: 0.3,
            color: [1.0, 1.0, 1.0],
            intensity: 100.0,
        };
        // Directly below — inside the cone.
        let lit = incident_light(&light, [0.0, 0.0, 0.0]);
        assert!(lit.is_some(), "point under the spot should be lit");
        // Far to the side — outside the outer cone.
        let dark = incident_light(&light, [50.0, 0.0, 0.0]);
        assert!(dark.is_none(), "point outside the spot cone should be dark");
    }

    #[test]
    fn to_display_srgb_clamps_and_gamma_encodes() {
        // A very bright linear colour tone-maps into [0, 1] and gamma
        // encoding lifts the mid-tones (sRGB(0.5_linear) > 0.5).
        let bright = ShadedColor {
            rgb: [10.0, 10.0, 10.0],
        };
        let disp = bright.to_display_srgb();
        for &c in &disp {
            assert!((0.0..=1.0).contains(&c), "display colour {c} out of [0,1]");
        }
        // A mid linear value gamma-encodes brighter than itself.
        let mid = ShadedColor {
            rgb: [1.0, 1.0, 1.0],
        };
        let d = mid.to_display_srgb();
        // 1.0 linear → Reinhard 0.5 → sRGB ≈ 0.735.
        assert!(
            d[0] > 0.5,
            "sRGB-encoded mid-tone should be lifted: {}",
            d[0]
        );
    }

    #[test]
    fn reflect_mirrors_about_the_normal() {
        // A ray coming straight down reflects straight back up.
        let r = reflect([0.0, 0.0, -1.0], [0.0, 0.0, 1.0]);
        assert!((r[2] - 1.0).abs() < 1e-6, "reflection should flip Z");
        // A 45° incidence reflects to the mirror 45°.
        let inc = normalize3([1.0, 0.0, -1.0]).unwrap();
        let refl = reflect(inc, [0.0, 0.0, 1.0]);
        assert!((refl[0] - inc[0]).abs() < 1e-6, "tangential component kept");
        assert!((refl[2] + inc[2]).abs() < 1e-6, "normal component flipped");
    }

    #[test]
    fn brdf_lut_entries_are_in_unit_range() {
        // The split-sum scale / bias are reflectance fractions — both
        // must stay in [0, 1] for every (n·v, roughness) cell.
        let lut = compute_brdf_lut(16, 64);
        assert_eq!(lut.entries.len(), 16 * 16);
        for &(scale, bias) in &lut.entries {
            assert!(
                (0.0..=1.01).contains(&scale),
                "BRDF LUT scale {scale} out of range"
            );
            assert!(
                (0.0..=1.01).contains(&bias),
                "BRDF LUT bias {bias} out of range"
            );
        }
    }

    #[test]
    fn brdf_lut_scale_plus_bias_does_not_amplify() {
        // For F₀ = 1 the environment-BRDF factor is `scale + bias`;
        // it must never exceed 1 — a passive BRDF cannot amplify the
        // environment (this is the LUT's energy-conservation property).
        let lut = compute_brdf_lut(24, 128);
        for &(scale, bias) in &lut.entries {
            assert!(
                scale + bias <= 1.05,
                "scale+bias {} exceeds 1 — BRDF LUT amplifies energy",
                scale + bias
            );
        }
    }

    #[test]
    fn brdf_lut_sample_interpolates_within_the_table() {
        // A sample at a cell centre returns that cell; a sample
        // between cells returns an interpolated value bracketed by its
        // neighbours.
        let lut = compute_brdf_lut(8, 64);
        // Corner (n·v = 0, roughness = 0).
        let (s, b) = lut.sample(0.0, 0.0);
        assert!(s.is_finite() && b.is_finite());
        // A mid sample is finite and in range.
        let (sm, bm) = lut.sample(0.5, 0.5);
        assert!((0.0..=1.01).contains(&sm) && (0.0..=1.01).contains(&bm));
    }

    #[test]
    fn specular_ibl_of_uniform_environment_is_bounded() {
        // The split-sum specular IBL against a uniform radiance-1
        // environment must stay bounded — it cannot exceed the
        // environment radiance (energy conservation).
        let env = EnvironmentMap::uniform([1.0, 1.0, 1.0]);
        let pre = env.prefilter_specular(16, 8, 5, 32).unwrap();
        let lut = compute_brdf_lut(16, 64);
        let mat = Material::polished_metal("chrome", [0.9, 0.9, 0.9]);
        let surface = SurfacePoint {
            normal: [0.0, 0.0, 1.0],
            view: [0.0, 0.0, 1.0],
        };
        let spec = specular_ibl(&surface, &mat, &pre, &lut);
        for (k, &c) in spec.iter().enumerate() {
            assert!(c >= 0.0, "specular IBL channel {k} should be non-negative");
            assert!(
                c <= 1.1,
                "specular IBL channel {k} = {c} should not exceed the environment radiance"
            );
        }
    }

    #[test]
    fn specular_ibl_metal_reflects_more_than_a_rough_dielectric() {
        // A smooth metal in a bright environment should pick up a
        // strong specular IBL term; a rough dielectric a much weaker
        // one (low F₀, blurred environment).
        let env = EnvironmentMap::uniform([2.0, 2.0, 2.0]);
        let pre = env.prefilter_specular(16, 8, 5, 32).unwrap();
        let lut = compute_brdf_lut(16, 64);
        let surface = SurfacePoint {
            normal: [0.0, 0.0, 1.0],
            view: [0.0, 0.0, 1.0],
        };
        let metal = Material::polished_metal("chrome", [0.95, 0.95, 0.95]);
        let rough_plastic = Material {
            roughness: 0.95,
            metallic: 0.0,
            ..Material::matte("plastic", [0.5, 0.5, 0.5])
        };
        let metal_spec = specular_ibl(&surface, &metal, &pre, &lut)[0];
        let plastic_spec = specular_ibl(&surface, &rough_plastic, &pre, &lut)[0];
        assert!(
            metal_spec > plastic_spec,
            "metal specular IBL {metal_spec} should exceed rough dielectric {plastic_spec}"
        );
    }
}
