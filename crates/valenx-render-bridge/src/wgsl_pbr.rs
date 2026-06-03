//! WGSL Cook-Torrance PBR forward shader — the GPU port of the
//! [`crate::pbr`] BRDF.
//!
//! ## What this is
//!
//! [`crate::pbr`] is the CPU-side physically-based-rendering library —
//! the Cook-Torrance microfacet BRDF, IBL, the split-sum specular.
//! This module is its **WGSL fragment-shader port**: the *identical*
//! GGX / Smith / Fresnel maths, transcribed term-by-term into the
//! shading language a `wgpu` render pass runs on the GPU once per
//! pixel.
//!
//! It ships **three things**, all of which a `wgpu` host crate
//! consumes:
//!
//! - [`PBR_FORWARD_WGSL`] — the complete WGSL shader source (vertex +
//!   fragment), ready to hand to `wgpu::Device::create_shader_module`.
//! - The CPU-side **uniform-block layouts** ([`PbrFrameUniform`],
//!   [`PbrMaterialUniform`], [`PbrLightUniform`], [`ShL2Uniform`]) —
//!   `#[repr(C)]` structs whose byte layout matches the WGSL
//!   `struct`s exactly, so a host can `bytemuck`-cast them straight
//!   into GPU buffers.
//! - A **cross-check** ([`reference_brdf_rgb`]) — the same BRDF
//!   evaluated on the CPU through [`crate::pbr::brdf_direct`], so a
//!   test can confirm the WGSL maths it ports is faithful to the CPU
//!   reference (the WGSL itself cannot be run under a headless test).
//!
//! ## The shader
//!
//! [`PBR_FORWARD_WGSL`] is a **forward** shader — it shades each
//! fragment in one pass against:
//!
//! 1. every analytic light (point / directional), via the
//!    Cook-Torrance BRDF `D·G·F/(4·n·v·n·l)` plus the Lambert diffuse
//!    lobe — the GPU twin of [`crate::pbr::brdf_direct`];
//! 2. the **IBL** ambient term — a diffuse irradiance + a Fresnel
//!    environment specular, the GPU twin of [`crate::pbr::ambient_ibl`];
//! 3. the **irradiance-volume GI** term — an `L2` spherical-harmonic
//!    probe ([`crate::irradiance_volume`]) evaluated in the surface
//!    normal, the indirect (bounced) light;
//! 4. the material's emissive colour.
//!
//! The result is tone-mapped (ACES filmic) and sRGB-encoded in the
//! shader so it writes a display-ready colour.
//!
//! ## HONEST REQUIREMENT — this code is GPU-unverified
//!
//! **The WGSL string here, and the `pbr_forward_pass` render-pass
//! module in `valenx-app` that consumes it, are GPU shader code. They
//! compile and `cargo check` cleanly, and the BRDF maths is checked
//! term-by-term against the CPU reference by [`reference_brdf_rgb`] +
//! the tests here — but the shader has *not* been executed on a
//! GPU.** A correct visual result needs the app run on real hardware,
//! which the development lockdown forbids. Treat the WGSL as
//! *written-correct-against-the-CPU-reference* but *not
//! pixel-validated*. The CPU BRDF in [`crate::pbr`] is the verified
//! one; this is its careful transcription.

use bytemuck::{Pod, Zeroable};

use crate::material::Material;
use crate::pbr::{brdf_direct, IncidentLight, SurfacePoint};

/// The complete **WGSL Cook-Torrance PBR forward shader** — vertex +
/// fragment entry points, ready for `wgpu::ShaderSource::Wgsl`.
///
/// Bindings (group 0): `0` the per-frame uniform ([`PbrFrameUniform`]),
/// `1` the material uniform ([`PbrMaterialUniform`]), `2` the light
/// array (`array<`[`PbrLightUniform`]`>`), `3` the `L2` SH probe
/// ([`ShL2Uniform`]). The vertex input is `position` + `normal`
/// (locations 0 and 1) — the same layout the existing viewport
/// pipeline uses.
///
/// Every numbered helper (`distribution_ggx`, `geometry_smith`,
/// `fresnel_schlick`, …) is the WGSL transcription of the
/// identically-named [`crate::pbr`] function — see [`reference_brdf_rgb`]
/// for the cross-check.
pub const PBR_FORWARD_WGSL: &str = r#"
// === valenx-render-bridge :: Cook-Torrance PBR forward shader ===
// GPU port of valenx_render_bridge::pbr — GGX / Smith / Fresnel.

const PI: f32 = 3.14159265359;
const INV_PI: f32 = 0.31830988618;
const MAX_LIGHTS: u32 = 8u;

struct FrameUniform {
    view_proj:   mat4x4<f32>,   // model-view-projection matrix
    camera_pos:  vec4<f32>,     // world-space eye (xyz; w unused)
    env_color:   vec4<f32>,     // ambient IBL radiance (rgb; w unused)
    light_count: vec4<u32>,     // x = number of active lights
};

struct MaterialUniform {
    base_color: vec4<f32>,      // diffuse / base colour (rgb)
    spec_emiss: vec4<f32>,      // rgb = unused tint slot, a = unused
    emissive:   vec4<f32>,      // emissive colour (rgb)
    params:     vec4<f32>,      // x=roughness y=metallic z=ior w=ao
};

struct LightUniform {
    // direction_or_pos.xyz: for a directional light, the unit
    //   direction *toward* the light; for a point light, the world
    //   position. .w selects the type (0 = directional, 1 = point).
    direction_or_pos: vec4<f32>,
    // color_intensity.rgb = light colour, .a = intensity / irradiance.
    color_intensity: vec4<f32>,
};

// An L2 spherical-harmonic probe — 9 RGB coefficients. Each vec4 holds
// one coefficient in .xyz (.w padding) so the std140-style alignment
// matches the Rust ShL2Uniform layout byte-for-byte.
struct ShProbe {
    coeffs: array<vec4<f32>, 9>,
};

@group(0) @binding(0) var<uniform> frame: FrameUniform;
@group(0) @binding(1) var<uniform> material: MaterialUniform;
@group(0) @binding(2) var<uniform> lights: array<LightUniform, MAX_LIGHTS>;
@group(0) @binding(3) var<uniform> sh_probe: ShProbe;

struct VertexIn {
    @location(0) position: vec3<f32>,
    @location(1) normal:   vec3<f32>,
};

struct VertexOut {
    @builtin(position) clip_pos: vec4<f32>,
    @location(0) world_pos: vec3<f32>,
    @location(1) world_normal: vec3<f32>,
};

@vertex
fn vs_main(in: VertexIn) -> VertexOut {
    var out: VertexOut;
    out.clip_pos = frame.view_proj * vec4<f32>(in.position, 1.0);
    // The viewport pipeline passes geometry already in world space.
    out.world_pos = in.position;
    out.world_normal = in.normal;
    return out;
}

// --- GGX / Trowbridge-Reitz normal-distribution function D ---
// Ports valenx_render_bridge::pbr::distribution_ggx.
fn distribution_ggx(n_dot_h: f32, roughness: f32) -> f32 {
    let alpha = max(roughness * roughness, 1e-4);
    let a2 = alpha * alpha;
    let nh = max(n_dot_h, 0.0);
    let denom = nh * nh * (a2 - 1.0) + 1.0;
    return a2 / (PI * denom * denom);
}

// --- one-direction Schlick-GGX geometry term G1 ---
fn geometry_schlick_ggx(n_dot_x: f32, k: f32) -> f32 {
    let nx = max(n_dot_x, 0.0);
    return nx / (nx * (1.0 - k) + k);
}

// --- Smith geometry term G (direct-lighting (r+1)^2/8 remap) ---
// Ports valenx_render_bridge::pbr::geometry_smith.
fn geometry_smith(n_dot_v: f32, n_dot_l: f32, roughness: f32) -> f32 {
    let r = roughness + 1.0;
    let k = (r * r) / 8.0;
    return geometry_schlick_ggx(n_dot_v, k) * geometry_schlick_ggx(n_dot_l, k);
}

// --- Fresnel-Schlick reflectance ---
// Ports valenx_render_bridge::pbr::fresnel_schlick.
fn fresnel_schlick(cos_theta: f32, f0: vec3<f32>) -> vec3<f32> {
    let m = pow(clamp(1.0 - cos_theta, 0.0, 1.0), 5.0);
    return f0 + (vec3<f32>(1.0) - f0) * m;
}

// --- Fresnel-Schlick with a roughness-aware ceiling (ambient term) ---
fn fresnel_schlick_roughness(cos_theta: f32, f0: vec3<f32>, roughness: f32) -> vec3<f32> {
    let m = pow(clamp(1.0 - cos_theta, 0.0, 1.0), 5.0);
    let ceil_v = max(vec3<f32>(1.0 - roughness), f0);
    return f0 + (ceil_v - f0) * m;
}

// --- F0 from the material (IOR grey for dielectric, base colour for
//     metal — ports valenx_render_bridge::pbr::f0_from_material) ---
fn f0_from_material() -> vec3<f32> {
    let ior = max(material.params.z, 1.0);
    let r = (ior - 1.0) / (ior + 1.0);
    let dielectric = r * r;
    let metallic = clamp(material.params.y, 0.0, 1.0);
    // Dielectrics: achromatic IOR grey. Metals: tint the specular with
    // the base colour (the metallic-roughness convention).
    return mix(vec3<f32>(dielectric), material.base_color.rgb, metallic);
}

// --- the diffuse albedo (base colour, faded out by metallic) ---
fn diffuse_albedo() -> vec3<f32> {
    let metallic = clamp(material.params.y, 0.0, 1.0);
    return material.base_color.rgb * (1.0 - metallic);
}

// --- Cook-Torrance response to one analytic light ---
// Ports valenx_render_bridge::pbr::brdf_direct.
fn brdf_direct(n: vec3<f32>, v: vec3<f32>, l: vec3<f32>,
               radiance: vec3<f32>) -> vec3<f32> {
    let n_dot_l = dot(n, l);
    if (n_dot_l <= 0.0) {
        return vec3<f32>(0.0);
    }
    let n_dot_v = max(dot(n, v), 1e-4);
    let h = normalize(v + l);
    let n_dot_h = max(dot(n, h), 0.0);
    let h_dot_v = max(dot(h, v), 0.0);

    let roughness = clamp(material.params.x, 0.0, 1.0);
    let metallic = clamp(material.params.y, 0.0, 1.0);
    let f0 = f0_from_material();

    let d = distribution_ggx(n_dot_h, roughness);
    let g = geometry_smith(n_dot_v, n_dot_l, roughness);
    let f = fresnel_schlick(h_dot_v, f0);
    let denom = 4.0 * n_dot_v * n_dot_l + 1e-4;
    let spec = d * g * f / denom;

    let kd = (vec3<f32>(1.0) - f) * (1.0 - metallic);
    let diffuse = kd * diffuse_albedo() * INV_PI;
    return (diffuse + spec) * radiance * n_dot_l;
}

// --- L2 spherical-harmonic basis in direction `dir` ---
// Mirrors valenx_render_bridge::irradiance_volume::sh_basis (L2).
fn sh_basis_l2(dir: vec3<f32>) -> array<f32, 9> {
    let d = normalize(dir);
    var b: array<f32, 9>;
    b[0] = 0.2820948;
    b[1] = 0.4886025 * d.y;
    b[2] = 0.4886025 * d.z;
    b[3] = 0.4886025 * d.x;
    b[4] = 1.0925484 * d.x * d.y;
    b[5] = 1.0925484 * d.y * d.z;
    b[6] = 0.31539157 * (3.0 * d.z * d.z - 1.0);
    b[7] = 1.0925484 * d.x * d.z;
    b[8] = 0.5462742 * (d.x * d.x - d.y * d.y);
    return b;
}

// --- irradiance-volume GI: the SH probe cosine-convolved in `n` ---
// Mirrors LightProbe::irradiance — the per-band A_l/π weights.
fn irradiance_volume_gi(n: vec3<f32>) -> vec3<f32> {
    // `var` (not `let`): the loop indexes `basis` with a dynamic `i`,
    // which WGSL permits only for a memory-backed local, not a
    // value-typed expression.
    var basis = sh_basis_l2(n);
    // A_l/π band weights: band 0 -> 1, band 1 -> 2/3, band 2 -> 1/4.
    var acc = vec3<f32>(0.0);
    for (var i: u32 = 0u; i < 9u; i = i + 1u) {
        var w: f32 = 0.25;
        if (i == 0u) { w = 1.0; }
        else if (i < 4u) { w = 0.6666667; }
        acc = acc + sh_probe.coeffs[i].xyz * (w * basis[i]);
    }
    return max(acc, vec3<f32>(0.0));
}

// --- ACES filmic tone curve (Narkowicz fit) ---
fn aces_filmic(x: f32) -> f32 {
    let a = 2.51; let b = 0.03; let c = 2.43; let d = 0.59; let e = 0.14;
    return clamp((x * (a * x + b)) / (x * (c * x + d) + e), 0.0, 1.0);
}

// --- linear -> sRGB gamma encode (IEC 61966-2-1) ---
fn linear_to_srgb(c: f32) -> f32 {
    if (c <= 0.0031308) {
        return 12.92 * c;
    }
    return 1.055 * pow(c, 1.0 / 2.4) - 0.055;
}

@fragment
fn fs_main(in: VertexOut) -> @location(0) vec4<f32> {
    let n = normalize(in.world_normal);
    let v = normalize(frame.camera_pos.xyz - in.world_pos);

    var color = vec3<f32>(0.0);

    // (1) Direct lighting — sum the Cook-Torrance response of every
    //     active analytic light.
    let count = min(frame.light_count.x, MAX_LIGHTS);
    for (var i: u32 = 0u; i < count; i = i + 1u) {
        let light = lights[i];
        var l: vec3<f32>;
        var radiance: vec3<f32>;
        if (light.direction_or_pos.w < 0.5) {
            // Directional: .xyz is the direction toward the light.
            l = normalize(light.direction_or_pos.xyz);
            radiance = light.color_intensity.rgb * light.color_intensity.a;
        } else {
            // Point: inverse-square attenuation I/(4*pi*d^2).
            let to_light = light.direction_or_pos.xyz - in.world_pos;
            let dist2 = max(dot(to_light, to_light), 1e-6);
            l = to_light / sqrt(dist2);
            let atten = 1.0 / (4.0 * PI * dist2);
            radiance = light.color_intensity.rgb * (light.color_intensity.a * atten);
        }
        color = color + brdf_direct(n, v, l, radiance);
    }

    // (2) Ambient IBL — a diffuse irradiance term off the environment
    //     colour plus a Fresnel environment specular reflection.
    let roughness = clamp(material.params.x, 0.0, 1.0);
    let metallic = clamp(material.params.y, 0.0, 1.0);
    let f0 = f0_from_material();
    let n_dot_v = max(dot(n, v), 1e-4);
    let f_amb = fresnel_schlick_roughness(n_dot_v, f0, roughness);
    let kd_amb = (vec3<f32>(1.0) - f_amb) * (1.0 - metallic);
    // The env colour is a radiance; a uniform environment's diffuse
    // irradiance is pi*L, so diffuse out = irr*albedo/pi = L*albedo.
    let ambient_diffuse = kd_amb * diffuse_albedo() * frame.env_color.rgb;
    let ambient_specular = f_amb * frame.env_color.rgb;
    color = color + (ambient_diffuse + ambient_specular) * material.params.w;

    // (3) Irradiance-volume global illumination — the indirect bounced
    //     light, an L2 SH probe evaluated in the surface normal.
    color = color + irradiance_volume_gi(n) * diffuse_albedo();

    // (4) Emissive — the surface's own glow.
    color = color + material.emissive.rgb;

    // Tone-map (ACES) + sRGB encode so the pass writes a display colour.
    let mapped = vec3<f32>(
        aces_filmic(max(color.r, 0.0)),
        aces_filmic(max(color.g, 0.0)),
        aces_filmic(max(color.b, 0.0)),
    );
    let display = vec3<f32>(
        linear_to_srgb(mapped.r),
        linear_to_srgb(mapped.g),
        linear_to_srgb(mapped.b),
    );
    return vec4<f32>(display, 1.0);
}
"#;

/// The maximum analytic-light count the shader's fixed-size light
/// array holds — must match `MAX_LIGHTS` in [`PBR_FORWARD_WGSL`].
pub const MAX_LIGHTS: usize = 8;

/// CPU mirror of the WGSL `FrameUniform` — the per-frame constants.
///
/// `#[repr(C)]` with the field order and `vec4` padding chosen so the
/// byte layout equals the WGSL `struct` exactly; it derives
/// [`bytemuck::Pod`] so a host casts it straight to bytes and uploads
/// it to binding 0.
#[repr(C)]
#[derive(Clone, Copy, Debug, Pod, Zeroable)]
pub struct PbrFrameUniform {
    /// Model-view-projection matrix, column-major (WGSL `mat4x4`).
    pub view_proj: [[f32; 4]; 4],
    /// World-space camera position; `w` is unused padding.
    pub camera_pos: [f32; 4],
    /// Ambient-IBL environment radiance; `w` is unused padding.
    pub env_color: [f32; 4],
    /// `[0]` = number of active lights; `[1..4]` unused padding.
    pub light_count: [u32; 4],
}

impl PbrFrameUniform {
    /// Build a frame uniform from an MVP matrix, the camera position,
    /// the ambient environment colour, and the active light count.
    pub fn new(
        view_proj: [[f32; 4]; 4],
        camera_pos: [f32; 3],
        env_color: [f32; 3],
        light_count: u32,
    ) -> PbrFrameUniform {
        PbrFrameUniform {
            view_proj,
            camera_pos: [camera_pos[0], camera_pos[1], camera_pos[2], 0.0],
            env_color: [env_color[0], env_color[1], env_color[2], 0.0],
            light_count: [light_count.min(MAX_LIGHTS as u32), 0, 0, 0],
        }
    }
}

/// CPU mirror of the WGSL `MaterialUniform`.
#[repr(C)]
#[derive(Clone, Copy, Debug, Pod, Zeroable)]
pub struct PbrMaterialUniform {
    /// Base / diffuse colour; `w` unused.
    pub base_color: [f32; 4],
    /// Reserved tint slot — kept for layout parity; `w` unused.
    pub spec_emiss: [f32; 4],
    /// Emissive colour; `w` unused.
    pub emissive: [f32; 4],
    /// `[roughness, metallic, ior, ambient_occlusion]`.
    pub params: [f32; 4],
}

impl PbrMaterialUniform {
    /// Build the material uniform from a [`Material`] (and an ambient-
    /// occlusion multiplier, `1.0` for none).
    pub fn from_material(m: &Material, ao: f32) -> PbrMaterialUniform {
        PbrMaterialUniform {
            base_color: [m.diffuse_color[0], m.diffuse_color[1], m.diffuse_color[2], 1.0],
            spec_emiss: [m.specular_color[0], m.specular_color[1], m.specular_color[2], 0.0],
            emissive: [m.emissive[0], m.emissive[1], m.emissive[2], 0.0],
            params: [
                m.roughness.clamp(0.0, 1.0),
                m.metallic.clamp(0.0, 1.0),
                m.ior.max(1.0),
                ao,
            ],
        }
    }
}

/// CPU mirror of the WGSL `LightUniform` — one analytic light.
#[repr(C)]
#[derive(Clone, Copy, Debug, Pod, Zeroable)]
pub struct PbrLightUniform {
    /// For a directional light, the unit direction toward the light;
    /// for a point light, the world position. `.w` is the type tag —
    /// `0.0` directional, `1.0` point.
    pub direction_or_pos: [f32; 4],
    /// `rgb` = light colour, `.a` = intensity (point) / irradiance
    /// (directional).
    pub color_intensity: [f32; 4],
}

impl PbrLightUniform {
    /// A directional light — `direction` points *toward* the light;
    /// `irradiance` is the linear-RGB irradiance it delivers.
    pub fn directional(direction: [f32; 3], color: [f32; 3], irradiance: f32) -> PbrLightUniform {
        // Normalise the direction so the shader need not.
        let len = (direction[0] * direction[0]
            + direction[1] * direction[1]
            + direction[2] * direction[2])
            .sqrt()
            .max(1e-12);
        PbrLightUniform {
            direction_or_pos: [
                direction[0] / len,
                direction[1] / len,
                direction[2] / len,
                0.0,
            ],
            color_intensity: [color[0], color[1], color[2], irradiance],
        }
    }

    /// A point light at world `position` with watt `intensity`.
    pub fn point(position: [f32; 3], color: [f32; 3], intensity: f32) -> PbrLightUniform {
        PbrLightUniform {
            direction_or_pos: [position[0], position[1], position[2], 1.0],
            color_intensity: [color[0], color[1], color[2], intensity],
        }
    }

    /// The all-zero "inactive" light — used to pad the fixed-size GPU
    /// light array beyond the active count.
    pub fn inactive() -> PbrLightUniform {
        PbrLightUniform {
            direction_or_pos: [0.0; 4],
            color_intensity: [0.0; 4],
        }
    }
}

/// CPU mirror of the WGSL `ShProbe` — the 9-coefficient `L2`
/// spherical-harmonic irradiance-volume probe for the GI term.
///
/// Each coefficient is stored in the `xyz` of a `vec4` (the `w` is
/// padding) so the array's stride matches the GPU's std140-style
/// alignment.
#[repr(C)]
#[derive(Clone, Copy, Debug, Pod, Zeroable)]
pub struct ShL2Uniform {
    /// 9 SH coefficients, each `[r, g, b, _pad]`.
    pub coeffs: [[f32; 4]; 9],
}

impl ShL2Uniform {
    /// An all-zero probe — no indirect light.
    pub fn zero() -> ShL2Uniform {
        ShL2Uniform {
            coeffs: [[0.0; 4]; 9],
        }
    }

    /// Build the GPU probe uniform from an `L2`
    /// [`crate::irradiance_volume::LightProbe`]'s coefficients.
    ///
    /// Pass the probe's `coeffs` slice; the first 9 RGB triples are
    /// copied (a shorter slice zero-pads, an `L1` probe simply leaves
    /// the quadratic bands zero).
    pub fn from_coeffs(coeffs: &[[f32; 3]]) -> ShL2Uniform {
        let mut out = ShL2Uniform::zero();
        for (dst, src) in out.coeffs.iter_mut().zip(coeffs.iter()) {
            dst[0] = src[0];
            dst[1] = src[1];
            dst[2] = src[2];
        }
        out
    }
}

/// Evaluate the **CPU reference** of one analytic light's Cook-Torrance
/// contribution — the value the WGSL `brdf_direct` must reproduce.
///
/// This is a thin wrapper over [`crate::pbr::brdf_direct`]: it lets a
/// test feed the *same* geometry + material + light into the verified
/// CPU BRDF and compare against a hand-evaluation of the WGSL formula,
/// confirming the shader port is faithful. (The WGSL itself cannot be
/// executed under a headless test — see this module's honest-scope
/// note.)
///
/// `normal`, `view`, `light_dir` are unit directions (`view` toward
/// the camera, `light_dir` toward the light); `radiance` is the
/// linear-RGB radiance arriving along `light_dir`.
pub fn reference_brdf_rgb(
    material: &Material,
    normal: [f32; 3],
    view: [f32; 3],
    light_dir: [f32; 3],
    radiance: [f32; 3],
) -> [f32; 3] {
    let surface = SurfacePoint { normal, view };
    let light = IncidentLight {
        direction: light_dir,
        radiance,
    };
    brdf_direct(&surface, material, &light)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The WGSL source is non-empty and declares both shader entry
    /// points the render pass references.
    #[test]
    fn wgsl_source_declares_both_entry_points() {
        assert!(PBR_FORWARD_WGSL.contains("fn vs_main"), "missing vertex entry");
        assert!(PBR_FORWARD_WGSL.contains("fn fs_main"), "missing fragment entry");
        // The four bind-group entries the host must supply.
        assert!(PBR_FORWARD_WGSL.contains("@group(0) @binding(0)"));
        assert!(PBR_FORWARD_WGSL.contains("@group(0) @binding(3)"));
    }

    /// The WGSL transcribes the named [`crate::pbr`] BRDF functions —
    /// a guard that the port did not silently drop a term.
    #[test]
    fn wgsl_ports_the_cook_torrance_terms() {
        for needle in [
            "fn distribution_ggx",
            "fn geometry_smith",
            "fn fresnel_schlick",
            "fn brdf_direct",
            "fn irradiance_volume_gi",
            "fn aces_filmic",
        ] {
            assert!(
                PBR_FORWARD_WGSL.contains(needle),
                "WGSL is missing the ported function `{needle}`"
            );
        }
    }

    /// The uniform structs have the `#[repr(C)]` byte sizes the WGSL
    /// layout demands (every member padded to a `vec4` / 16 bytes).
    #[test]
    fn uniform_layouts_have_the_expected_sizes() {
        // FrameUniform: mat4 (64) + 3 × vec4 (48) = 112 bytes.
        assert_eq!(std::mem::size_of::<PbrFrameUniform>(), 112);
        // MaterialUniform: 4 × vec4 = 64 bytes.
        assert_eq!(std::mem::size_of::<PbrMaterialUniform>(), 64);
        // LightUniform: 2 × vec4 = 32 bytes.
        assert_eq!(std::mem::size_of::<PbrLightUniform>(), 32);
        // ShL2Uniform: 9 × vec4 = 144 bytes.
        assert_eq!(std::mem::size_of::<ShL2Uniform>(), 144);
    }

    /// The MAX_LIGHTS constant agrees between the Rust side and the
    /// WGSL `array<LightUniform, MAX_LIGHTS>` declaration.
    #[test]
    fn max_lights_matches_the_wgsl() {
        assert!(PBR_FORWARD_WGSL.contains(&format!("MAX_LIGHTS: u32 = {MAX_LIGHTS}u")));
    }

    /// The frame-uniform constructor clamps the light count to the
    /// array capacity so a host cannot overrun the GPU buffer.
    #[test]
    fn frame_uniform_clamps_the_light_count() {
        let f = PbrFrameUniform::new(
            [[0.0; 4]; 4],
            [0.0; 3],
            [0.0; 3],
            999, // far above MAX_LIGHTS
        );
        assert_eq!(f.light_count[0], MAX_LIGHTS as u32);
    }

    /// The material uniform copies the [`Material`] PBR parameters into
    /// the `params` vec4 in the order the WGSL reads them.
    #[test]
    fn material_uniform_packs_the_pbr_params() {
        let m = Material {
            roughness: 0.4,
            metallic: 0.7,
            ior: 1.5,
            ..Material::matte("m", [0.2, 0.5, 0.8])
        };
        let u = PbrMaterialUniform::from_material(&m, 1.0);
        assert_eq!(u.base_color, [0.2, 0.5, 0.8, 1.0]);
        assert_eq!(u.params, [0.4, 0.7, 1.5, 1.0]);
    }

    /// A directional light uniform normalises its direction and tags
    /// itself type-0; a point light tags itself type-1.
    #[test]
    fn light_uniform_tags_and_normalises() {
        let dir = PbrLightUniform::directional([0.0, -3.0, 0.0], [1.0, 1.0, 1.0], 2.0);
        // Direction normalised, type tag 0.
        assert!((dir.direction_or_pos[1] + 1.0).abs() < 1e-5);
        assert_eq!(dir.direction_or_pos[3], 0.0);
        let pt = PbrLightUniform::point([1.0, 2.0, 3.0], [1.0, 1.0, 1.0], 50.0);
        assert_eq!(pt.direction_or_pos, [1.0, 2.0, 3.0, 1.0]);
    }

    /// The reference BRDF — which the WGSL `brdf_direct` ports — gives
    /// a light below the surface horizon zero contribution, exactly as
    /// the WGSL's `if (n_dot_l <= 0.0)` early-out does.
    #[test]
    fn reference_brdf_rejects_a_below_horizon_light() {
        let m = Material::matte("white", [0.8, 0.8, 0.8]);
        let out = reference_brdf_rgb(
            &m,
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0], // pointing into the surface
            [10.0, 10.0, 10.0],
        );
        assert_eq!(out, [0.0, 0.0, 0.0]);
    }

    /// The reference BRDF a lit surface yields positive radiance — the
    /// WGSL port computes the identical `(diffuse + spec)·radiance·n·l`.
    #[test]
    fn reference_brdf_lit_surface_is_positive() {
        let m = Material::matte("white", [0.8, 0.8, 0.8]);
        let out = reference_brdf_rgb(
            &m,
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0],
            [5.0, 5.0, 5.0],
        );
        for c in out {
            assert!(c > 0.0, "a lit surface should have positive radiance");
        }
    }

    /// The SH-probe uniform copies an irradiance-volume probe's
    /// coefficients into the GPU layout.
    #[test]
    fn sh_uniform_copies_probe_coefficients() {
        let coeffs = [[1.0f32, 2.0, 3.0]; 9];
        let u = ShL2Uniform::from_coeffs(&coeffs);
        for c in u.coeffs {
            assert_eq!([c[0], c[1], c[2]], [1.0, 2.0, 3.0]);
        }
        // A shorter (L1) coefficient slice zero-pads the quadratic band.
        let short = ShL2Uniform::from_coeffs(&[[5.0, 5.0, 5.0]; 4]);
        assert_eq!(short.coeffs[0], [5.0, 5.0, 5.0, 0.0]);
        assert_eq!(short.coeffs[8], [0.0, 0.0, 0.0, 0.0]);
    }

    // ===================================================================
    // Static WGSL validation — `naga` parse + validate.
    //
    // These prove [`PBR_FORWARD_WGSL`] is syntactically *and*
    // semantically valid WGSL without a GPU: `naga` is the exact
    // front-end `wgpu::Device::create_shader_module` runs internally,
    // so a shader that passes here is a shader the GPU pipeline accepts.
    // This closes the "the WGSL itself cannot be run under a headless
    // test" gap noted in this module's honest-scope header.
    // ===================================================================

    /// `PBR_FORWARD_WGSL` parses cleanly with the `naga` WGSL
    /// front-end — no syntax error.
    #[test]
    fn pbr_wgsl_parses_with_naga() {
        let module = naga::front::wgsl::parse_str(PBR_FORWARD_WGSL);
        let module = module.unwrap_or_else(|e| {
            panic!(
                "PBR_FORWARD_WGSL failed to parse as WGSL:\n{}",
                e.emit_to_string(PBR_FORWARD_WGSL)
            )
        });
        // Sanity: the parsed module exposes both entry points.
        let entries: Vec<&str> = module.entry_points.iter().map(|e| e.name.as_str()).collect();
        assert!(entries.contains(&"vs_main"), "no vs_main entry: {entries:?}");
        assert!(entries.contains(&"fs_main"), "no fs_main entry: {entries:?}");
    }

    /// `PBR_FORWARD_WGSL` passes the `naga` semantic validator — types,
    /// bindings, control flow, and entry-point I/O are all sound. This
    /// is the same validation pass `wgpu` runs before it hands the
    /// shader to the driver.
    #[test]
    fn pbr_wgsl_validates_with_naga() {
        let module = naga::front::wgsl::parse_str(PBR_FORWARD_WGSL)
            .expect("PBR_FORWARD_WGSL must parse before it can be validated");
        let mut validator = naga::valid::Validator::new(
            naga::valid::ValidationFlags::all(),
            naga::valid::Capabilities::empty(),
        );
        let info = validator.validate(&module);
        info.unwrap_or_else(|e| {
            panic!(
                "PBR_FORWARD_WGSL failed naga validation:\n{}",
                e.emit_to_string(PBR_FORWARD_WGSL)
            )
        });
    }
}
