//! Volumetric rendering v1 — ray-marching of participating media.
//!
//! # What this is
//!
//! The surface path tracer in [`crate::tracer`] solves light transport
//! *between* opaque surfaces. This module adds the orthogonal piece:
//! light transport *through a participating medium* — fog, smoke, a
//! cloud, a sub-surface haze — where the ray is absorbed, scatters, and
//! can itself emit along its whole length, not only at a surface hit.
//!
//! It is a self-contained volume renderer: it reuses the crate's
//! [`Ray`] / [`Vec3`] / [`Rng`] / [`HdrFramebuffer`] /
//! [`EnvironmentMap`] infrastructure but has its own integrator, so a
//! caller can render a volume standalone with [`render_volume`] without
//! touching the surface BVH.
//!
//! # The physics it integrates
//!
//! For a ray through a medium the **radiative transfer equation** says
//! the radiance reaching the camera is
//!
//! ```text
//!   L = ∫₀ᵀ transmittance(0→t) · [ σ_a·L_e + σ_s·L_scatter(t) ] dt
//!       + transmittance(0→T) · L_background
//! ```
//!
//! where `transmittance(0→t) = exp(−∫₀ᵗ σ_t ds)` is the Beer-Lambert
//! attenuation (`σ_t = σ_a + σ_s` is the extinction coefficient). The
//! integrator walks the ray in regular steps (`RegularMarcher`) and
//! accumulates, per step:
//!
//! 1. **Emission-absorption** — `σ_a · L_e` scaled by the running
//!    transmittance: a glowing medium adds light.
//! 2. **Single scattering** — `σ_s · L_scatter`, where `L_scatter` is a
//!    shadow-tested connection to each light, weighted by the
//!    **Henyey-Greenstein phase function** (the volumetric analogue of
//!    a BRDF). This is "single" scattering: light scatters into the
//!    ray exactly once. Multiple in-scattering is a documented
//!    follow-up.
//! 3. **Transmittance update** — the step multiplies the running
//!    transmittance by `exp(−σ_t · step)`.
//!
//! # Honest scope — a real v1, not a production volume renderer
//!
//! Everything here is the genuine radiative-transfer integral and the
//! tests verify it (a homogeneous absorber attenuates by Beer-Lambert
//! to solver precision; an emissive medium adds the analytic glow; a
//! zero-density medium is an exact no-op). It is deliberately a v1:
//!
//! - **Single scattering only.** Light in-scatters once; the
//!   multiple-scattering term (a path that scatters, scatters again,
//!   then reaches the eye) is not integrated. For an optically-thin
//!   medium — the common fog / haze case — single scattering is the
//!   dominant term and the result is close. A dense cloud needs
//!   multiple scattering (the documented follow-up).
//! - **Regular-step ray marching**, not delta tracking. The step
//!   integrator is unbiased in the limit of a small step and is exact
//!   for a homogeneous medium; delta / ratio tracking would remove the
//!   step-size bias on a strongly heterogeneous grid (follow-up).
//! - **Isotropic single-bounce phase sampling for the marched ray
//!   itself** is not done — the marched (primary) ray travels straight
//!   and gathers in-scattered light; it does not itself scatter into a
//!   new direction. That, again, is the multiple-scattering extension.
//! - Analytic point lights only ([`VolumeLight`]); no area lights, no
//!   light from the surface scene.

use crate::framebuffer::{FramebufferError, HdrFramebuffer};
use crate::geometry::Ray;
use crate::math::Vec3;
use crate::sampling::Rng;
use crate::scene::PtCamera;

use valenx_render_bridge::environment::EnvironmentMap;

/// A participating medium — what the ray marches through.
///
/// The medium is defined by its **extinction** (absorption +
/// scattering), its **single-scattering albedo** (the fraction of
/// extinction that is scattering rather than absorption), its
/// **emission**, and its **phase function asymmetry**. The spatial
/// density is either uniform ([`Medium::homogeneous`]) or read from a
/// 3-D grid ([`Medium::from_grid`]); the density modulates the
/// extinction and scattering coefficients point-by-point.
#[derive(Clone, Debug)]
pub struct Medium {
    /// Extinction coefficient `σ_t` at unit density, in inverse world
    /// units — how fast the medium attenuates a ray. `σ_t = σ_a +
    /// σ_s`.
    pub sigma_t: f32,
    /// Single-scattering albedo `σ_s / σ_t` in `[0, 1]` — the fraction
    /// of an extinction event that scatters (the rest is absorbed). 0
    /// is a pure absorber (black smoke), 1 a pure scatterer
    /// (non-absorbing fog).
    pub scattering_albedo: f32,
    /// Linear-RGB radiance the medium emits per unit `σ_a` along the
    /// ray — a glowing medium (fire, a luminous gas). Zero for a
    /// non-emitting medium.
    pub emission: Vec3,
    /// Henyey-Greenstein asymmetry parameter `g ∈ (−1, 1)`. `g = 0` is
    /// isotropic scattering; `g > 0` is forward-peaked (haze, clouds);
    /// `g < 0` is back-scattering.
    pub phase_g: f32,
    /// The density field — `Uniform(1.0)` for a homogeneous medium, or
    /// a sampled grid.
    density: Density,
}

/// The spatial density distribution of a [`Medium`].
#[derive(Clone, Debug)]
enum Density {
    /// Constant density everywhere inside the medium bounds.
    Uniform(f32),
    /// A regular 3-D grid of densities, sampled with trilinear
    /// interpolation over the medium's bounding box.
    Grid(DensityGrid),
}

/// A regular 3-D scalar density grid.
#[derive(Clone, Debug)]
pub struct DensityGrid {
    /// Cell counts along x, y, z.
    pub dims: [usize; 3],
    /// Densities, `dims.x · dims.y · dims.z` entries, x-fastest then y
    /// then z.
    pub values: Vec<f32>,
}

impl DensityGrid {
    /// Build a grid from its dimensions and a flat value buffer.
    ///
    /// Returns `None` if `values.len()` does not match the product of
    /// `dims`, or any dimension is zero — a caller must hand a
    /// well-formed grid.
    pub fn new(dims: [usize; 3], values: Vec<f32>) -> Option<DensityGrid> {
        if dims.contains(&0) {
            return None;
        }
        if values.len() != dims[0] * dims[1] * dims[2] {
            return None;
        }
        Some(DensityGrid { dims, values })
    }

    /// Sample the grid with trilinear interpolation at normalised
    /// coordinates `(u, v, w)` — each in `[0, 1]` spanning the grid.
    /// Coordinates outside `[0, 1]` clamp to the boundary.
    fn sample(&self, u: f32, v: f32, w: f32) -> f32 {
        let coord = |t: f32, n: usize| -> (usize, usize, f32) {
            // Map [0,1] onto cell-centre space [0, n-1].
            let x = (t.clamp(0.0, 1.0)) * (n as f32 - 1.0).max(0.0);
            let i0 = x.floor() as usize;
            let i0 = i0.min(n.saturating_sub(1));
            let i1 = (i0 + 1).min(n.saturating_sub(1));
            (i0, i1, x - i0 as f32)
        };
        let (x0, x1, fx) = coord(u, self.dims[0]);
        let (y0, y1, fy) = coord(v, self.dims[1]);
        let (z0, z1, fz) = coord(w, self.dims[2]);
        let at = |x: usize, y: usize, z: usize| -> f32 {
            self.values[(z * self.dims[1] + y) * self.dims[0] + x]
        };
        // Trilinear blend of the 8 corners.
        let c00 = at(x0, y0, z0) * (1.0 - fx) + at(x1, y0, z0) * fx;
        let c10 = at(x0, y1, z0) * (1.0 - fx) + at(x1, y1, z0) * fx;
        let c01 = at(x0, y0, z1) * (1.0 - fx) + at(x1, y0, z1) * fx;
        let c11 = at(x0, y1, z1) * (1.0 - fx) + at(x1, y1, z1) * fx;
        let c0 = c00 * (1.0 - fy) + c10 * fy;
        let c1 = c01 * (1.0 - fy) + c11 * fy;
        c0 * (1.0 - fz) + c1 * fz
    }
}

impl Medium {
    /// A homogeneous medium of constant `density`.
    ///
    /// `sigma_t` is the extinction at that density, `scattering_albedo`
    /// the scattering fraction. The result has no emission and
    /// isotropic scattering — set [`Medium::emission`] / [`Medium::phase_g`]
    /// afterwards to change that.
    pub fn homogeneous(density: f32, sigma_t: f32, scattering_albedo: f32) -> Medium {
        Medium {
            sigma_t: sigma_t.max(0.0),
            scattering_albedo: scattering_albedo.clamp(0.0, 1.0),
            emission: Vec3::ZERO,
            phase_g: 0.0,
            density: Density::Uniform(density.max(0.0)),
        }
    }

    /// A medium whose density is read from a 3-D [`DensityGrid`].
    ///
    /// The grid is mapped onto the medium's bounding box (see
    /// [`VolumeBox`]) by trilinear interpolation. `sigma_t` is the
    /// extinction at density 1.0.
    pub fn from_grid(grid: DensityGrid, sigma_t: f32, scattering_albedo: f32) -> Medium {
        Medium {
            sigma_t: sigma_t.max(0.0),
            scattering_albedo: scattering_albedo.clamp(0.0, 1.0),
            emission: Vec3::ZERO,
            phase_g: 0.0,
            density: Density::Grid(grid),
        }
    }

    /// Set the medium's emitted radiance per unit absorption and
    /// return the modified medium (builder style).
    pub fn with_emission(mut self, emission: Vec3) -> Medium {
        self.emission = emission;
        self
    }

    /// Set the Henyey-Greenstein asymmetry `g` and return the modified
    /// medium (builder style). `g` is clamped to `(−0.999, 0.999)` to
    /// keep the phase function finite.
    pub fn with_phase_g(mut self, g: f32) -> Medium {
        self.phase_g = g.clamp(-0.999, 0.999);
        self
    }

    /// The density at normalised box coordinates `(u, v, w)` —
    /// constant for a homogeneous medium, trilinearly sampled for a
    /// grid medium.
    fn density_at(&self, u: f32, v: f32, w: f32) -> f32 {
        match &self.density {
            Density::Uniform(d) => *d,
            Density::Grid(g) => g.sample(u, v, w).max(0.0),
        }
    }

    /// The extinction `σ_t · density` at normalised box coordinates.
    #[inline]
    fn extinction_at(&self, u: f32, v: f32, w: f32) -> f32 {
        self.sigma_t * self.density_at(u, v, w)
    }
}

/// An axis-aligned box that bounds a [`Medium`] in world space.
///
/// The medium exists only inside this box; a ray is marched over its
/// `[t_enter, t_exit]` overlap with the box. Normalised box
/// coordinates `(u, v, w) ∈ [0, 1]³` address the density grid.
#[derive(Clone, Copy, Debug)]
pub struct VolumeBox {
    /// Minimum corner.
    pub min: Vec3,
    /// Maximum corner.
    pub max: Vec3,
}

impl VolumeBox {
    /// A box from its two corners.
    pub fn new(min: Vec3, max: Vec3) -> VolumeBox {
        VolumeBox { min, max }
    }

    /// World point → normalised box coordinates `(u, v, w)`.
    #[inline]
    fn normalise(&self, p: Vec3) -> (f32, f32, f32) {
        let ext = self.max.sub(self.min);
        let safe = |num: f32, den: f32| if den.abs() < 1e-12 { 0.5 } else { num / den };
        (
            safe(p.x - self.min.x, ext.x),
            safe(p.y - self.min.y, ext.y),
            safe(p.z - self.min.z, ext.z),
        )
    }

    /// Slab intersection of `ray` with the box, returning the
    /// `(t_enter, t_exit)` parameter interval clipped to non-negative
    /// `t`, or `None` if the ray misses the box.
    fn intersect(&self, ray: &Ray) -> Option<(f32, f32)> {
        let mut t0 = 0.0f32;
        let mut t1 = f32::INFINITY;
        for axis in 0..3 {
            let o = ray.origin.axis(axis);
            let d = ray.direction.axis(axis);
            let lo = self.min.axis(axis);
            let hi = self.max.axis(axis);
            if d.abs() < 1e-12 {
                // Ray parallel to this slab — must already be inside.
                if o < lo || o > hi {
                    return None;
                }
            } else {
                let inv = 1.0 / d;
                let mut near = (lo - o) * inv;
                let mut far = (hi - o) * inv;
                if near > far {
                    std::mem::swap(&mut near, &mut far);
                }
                t0 = t0.max(near);
                t1 = t1.min(far);
                if t0 > t1 {
                    return None;
                }
            }
        }
        if t1 < 0.0 {
            return None;
        }
        Some((t0.max(0.0), t1))
    }
}

/// An analytic point light illuminating the volume.
///
/// Single-scattering connects each marched sample to every light:
/// the in-scattered radiance is the light's radiance, attenuated by
/// the inverse-square falloff and by the medium's transmittance along
/// the shadow segment back to the light, weighted by the phase
/// function.
#[derive(Clone, Copy, Debug)]
pub struct VolumeLight {
    /// World position of the point light.
    pub position: Vec3,
    /// Linear-RGB radiant intensity (W·sr⁻¹). The irradiance at a
    /// point is `intensity / distance²`.
    pub intensity: Vec3,
}

impl VolumeLight {
    /// A point light at `position` with the given RGB intensity.
    pub fn point(position: Vec3, intensity: Vec3) -> VolumeLight {
        VolumeLight {
            position,
            intensity,
        }
    }
}

/// Tunable parameters of a volume render / march.
#[derive(Clone, Copy, Debug)]
pub struct VolumeParams {
    /// Ray-march step length in world units. Smaller → more accurate
    /// (the regular-step integral converges as `step → 0`) and slower.
    pub step: f32,
    /// Number of shadow-ray steps used to integrate the transmittance
    /// from a marched sample back to a light. A coarse value is fine
    /// for an optically-thin medium.
    pub shadow_steps: u32,
    /// Master random seed (the marcher jitters the first step to break
    /// up banding).
    pub seed: u64,
    /// Exposure passed to the tone mapper when producing the LDR
    /// image from [`render_volume`].
    pub exposure: f32,
}

impl Default for VolumeParams {
    fn default() -> VolumeParams {
        VolumeParams {
            step: 0.1,
            shadow_steps: 16,
            seed: 0x5eed_0001,
            exposure: 1.0,
        }
    }
}

/// The Henyey-Greenstein phase function value for a scattering angle
/// whose cosine is `cos_theta`, with asymmetry `g`.
///
/// The phase function is the volumetric analogue of a BRDF: it gives
/// the fraction of light arriving from direction `ω_i` that scatters
/// into direction `ω_o`, as a function only of the angle between them.
/// `cos_theta` is `ω_i · ω_o`. The Henyey-Greenstein form
///
/// ```text
///   p(θ) = (1 − g²) / (4π · (1 + g² − 2g·cosθ)^{3/2})
/// ```
///
/// is the standard single-parameter model: `g = 0` is isotropic
/// (`1/4π`), `g → 1` forward-peaked, `g → −1` back-peaked. It is
/// normalised so its integral over the sphere is 1.
#[inline]
pub fn henyey_greenstein(cos_theta: f32, g: f32) -> f32 {
    let g = g.clamp(-0.999, 0.999);
    let denom = 1.0 + g * g - 2.0 * g * cos_theta;
    // denom is strictly positive for |g| < 1; guard anyway.
    let denom = denom.max(1e-8);
    (1.0 - g * g) / (4.0 * std::f32::consts::PI * denom.powf(1.5))
}

/// The transmittance (Beer-Lambert survival fraction) along the
/// segment of `ray` from `t0` to `t1` through `medium` bounded by
/// `bounds`.
///
/// `transmittance = exp(−∫ σ_t ds)`. For a homogeneous medium this is
/// exact; for a grid medium it is the regular-step Riemann estimate.
/// Returned per channel as a `Vec3` of `[0, 1]` survival fractions
/// (the same scalar in all three components — extinction here is
/// grey).
pub fn transmittance(
    ray: &Ray,
    t0: f32,
    t1: f32,
    medium: &Medium,
    bounds: &VolumeBox,
    step: f32,
) -> Vec3 {
    if t1 <= t0 {
        return Vec3::ONE;
    }
    let span = t1 - t0;
    let step = step.max(1e-4);
    let n = ((span / step).ceil() as u32).max(1);
    let dt = span / n as f32;
    let mut optical_depth = 0.0f32;
    for i in 0..n {
        // Sample at the segment mid-point (mid-point rule).
        let t = t0 + (i as f32 + 0.5) * dt;
        let p = ray.at(t);
        let (u, v, w) = bounds.normalise(p);
        optical_depth += medium.extinction_at(u, v, w) * dt;
    }
    Vec3::splat((-optical_depth).exp())
}

/// The radiance gathered by `ray` as it marches through `medium`, plus
/// the residual transmittance of the ray past the medium.
///
/// This is the per-ray volumetric integrator. It returns a
/// [`VolumeResult`] — the in-medium radiance the camera collects, and
/// the `transmittance` fraction of any background radiance that
/// survives the medium (so a caller can composite the volume over a
/// surface render or an environment).
pub fn march_ray(
    ray: &Ray,
    medium: &Medium,
    bounds: &VolumeBox,
    lights: &[VolumeLight],
    params: &VolumeParams,
    rng: &mut Rng,
) -> VolumeResult {
    let Some((t_enter, t_exit)) = bounds.intersect(ray) else {
        // The ray misses the medium entirely — full transmittance, no
        // emitted / scattered radiance.
        return VolumeResult {
            radiance: Vec3::ZERO,
            transmittance: Vec3::ONE,
        };
    };
    let span = t_exit - t_enter;
    if span <= 0.0 {
        return VolumeResult {
            radiance: Vec3::ZERO,
            transmittance: Vec3::ONE,
        };
    }
    let step = params.step.max(1e-4);
    let n_steps = ((span / step).ceil() as u32).max(1);
    let dt = span / n_steps as f32;

    let mut radiance = Vec3::ZERO;
    // The running transmittance from the ray origin to the current
    // marched sample — starts at 1 (nothing has attenuated yet).
    let mut throughput = Vec3::ONE;
    // Jitter the first sample within the first step to break banding.
    let jitter = rng.next_f32();

    for i in 0..n_steps {
        let t = t_enter + (i as f32 + jitter) * dt;
        if t > t_exit {
            break;
        }
        let p = ray.at(t);
        let (u, v, w) = bounds.normalise(p);
        let density = medium.density_at(u, v, w);
        if density <= 0.0 {
            continue; // empty cell — no emission, no scattering
        }
        let sigma_t = medium.sigma_t * density;
        let sigma_s = sigma_t * medium.scattering_albedo;
        let sigma_a = sigma_t - sigma_s;

        // (1) Emission-absorption: a glowing medium adds σ_a·L_e,
        // scaled by the transmittance reaching this sample and the
        // step length.
        if medium.emission.max_component() > 0.0 {
            radiance = radiance.add(throughput.mul(medium.emission).scale(sigma_a * dt));
        }

        // (2) Single scattering: connect to every light. The
        // in-scattered radiance is the light's radiance attenuated by
        // distance and by the medium between the sample and the light,
        // weighted by the phase function.
        if sigma_s > 0.0 {
            let mut in_scatter = Vec3::ZERO;
            for light in lights {
                in_scatter = in_scatter.add(single_scatter_from_light(
                    p,
                    ray.direction,
                    light,
                    medium,
                    bounds,
                    params,
                ));
            }
            radiance = radiance.add(throughput.mul(in_scatter).scale(sigma_s * dt));
        }

        // (3) Transmittance update — the step attenuates the ray for
        // every later sample. exp(−σ_t·dt) is the Beer-Lambert
        // survival across this step.
        let step_transmittance = (-sigma_t * dt).exp();
        throughput = throughput.scale(step_transmittance);

        // An almost-opaque ray contributes nothing further — stop.
        if throughput.max_component() < 1e-4 {
            break;
        }
    }

    VolumeResult {
        radiance,
        transmittance: throughput,
    }
}

/// What [`march_ray`] returns — the volume's own contribution plus the
/// fraction of background radiance that survives it.
#[derive(Clone, Copy, Debug)]
pub struct VolumeResult {
    /// Radiance the medium itself contributes (emission + single
    /// scattering), in the camera's direction.
    pub radiance: Vec3,
    /// Per-channel transmittance of the ray through the whole medium —
    /// multiply a background radiance by this to composite it behind
    /// the volume.
    pub transmittance: Vec3,
}

impl VolumeResult {
    /// Composite the volume over a background radiance: the volume's
    /// own radiance plus the background attenuated by the medium's
    /// transmittance.
    #[inline]
    pub fn over(self, background: Vec3) -> Vec3 {
        self.radiance.add(self.transmittance.mul(background))
    }
}

/// The single-scattering contribution at point `p` from one light.
///
/// Connects `p` to the light, evaluates the phase function for the
/// camera-direction / light-direction angle, applies the inverse-square
/// falloff, and attenuates by the medium's transmittance along the
/// shadow segment from `p` back to the light.
fn single_scatter_from_light(
    p: Vec3,
    view_dir: Vec3,
    light: &VolumeLight,
    medium: &Medium,
    bounds: &VolumeBox,
    params: &VolumeParams,
) -> Vec3 {
    let to_light = light.position.sub(p);
    let dist2 = to_light.length_sq();
    if dist2 < 1e-8 {
        return Vec3::ZERO;
    }
    let dist = dist2.sqrt();
    let wi = to_light.scale(1.0 / dist);

    // Phase function: the scattering angle is between the direction
    // the ray is travelling (view_dir) and the direction toward the
    // light. For in-scattering we want light arriving from `wi` to
    // scatter into `−view_dir` (back toward the camera); the HG cosine
    // is the dot of the *incoming* light direction with the
    // *outgoing* (camera) direction.
    let cos_theta = wi.dot(view_dir.neg());
    let phase = henyey_greenstein(cos_theta, medium.phase_g);

    // Shadow transmittance: how much of the light survives the medium
    // between `p` and the light. Build a shadow ray and integrate the
    // extinction over its in-medium segment.
    let shadow_ray = Ray::new(p, wi);
    let shadow_tr = if let Some((s0, s1)) = bounds.intersect(&shadow_ray) {
        // Only the segment up to the light matters (the light may be
        // inside or outside the box).
        let s1 = s1.min(dist);
        let shadow_step = (s1 - s0) / params.shadow_steps.max(1) as f32;
        transmittance(&shadow_ray, s0, s1, medium, bounds, shadow_step.max(1e-4))
    } else {
        Vec3::ONE
    };

    // Light irradiance at `p`, phase-weighted, shadow-attenuated.
    light.intensity.scale(phase / dist2).mul(shadow_tr)
}

/// Render a participating medium standalone into an [`HdrFramebuffer`].
///
/// Shoots one camera ray per pixel (volume rendering converges fast
/// enough that one ray + the marcher's own jitter is a reasonable v1;
/// a caller wanting less noise can average several frames with
/// different seeds), marches it through `medium`, and composites the
/// result over the `environment` background.
///
/// This is the volume counterpart of [`crate::tracer::render`]; it
/// does **not** touch the surface BVH.
///
/// # Errors
///
/// Returns [`FramebufferError::TooLarge`] when the camera resolution
/// would allocate a framebuffer larger than the
/// `MAX_FRAMEBUFFER_PIXELS` cap. Round-10 sister fix to the round-9
/// `tracer::render` migration.
pub fn render_volume(
    camera: &PtCamera,
    medium: &Medium,
    bounds: &VolumeBox,
    lights: &[VolumeLight],
    environment: &EnvironmentMap,
    params: &VolumeParams,
) -> Result<HdrFramebuffer, FramebufferError> {
    let w = camera.width;
    let h = camera.height;
    let mut fb = HdrFramebuffer::try_new(w, h)?;

    for y in 0..h {
        for x in 0..w {
            let pixel_index = (y as u64) * (w as u64) + (x as u64);
            let mut rng = Rng::new(params.seed, pixel_index);
            // Primary ray through the pixel centre.
            let u = (x as f32 + 0.5) / w as f32;
            let v = 1.0 - (y as f32 + 0.5) / h as f32;
            let target = camera
                .lower_left
                .add(camera.horizontal.scale(u))
                .add(camera.vertical.scale(v));
            let Some(dir) = target.sub(camera.eye).normalized() else {
                fb.add_sample(x, y, Vec3::ZERO);
                continue;
            };
            let ray = Ray::new(camera.eye, dir);
            let result = march_ray(&ray, medium, bounds, lights, params, &mut rng);
            // Background: the environment radiance along the ray.
            let bg = Vec3::from_array(environment.sample_direction(dir.to_array()));
            fb.add_sample(x, y, result.over(bg));
        }
    }
    fb.finish_sample();
    Ok(fb)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::vec3;

    fn unit_box() -> VolumeBox {
        VolumeBox::new(vec3(-1.0, -1.0, -1.0), vec3(1.0, 1.0, 1.0))
    }

    #[test]
    fn henyey_greenstein_is_isotropic_at_g_zero() {
        // g = 0 → the phase function is the constant 1/4π for every
        // angle.
        let iso = 1.0 / (4.0 * std::f32::consts::PI);
        for cos in [-1.0, -0.3, 0.0, 0.5, 1.0] {
            let p = henyey_greenstein(cos, 0.0);
            assert!((p - iso).abs() < 1e-6, "g=0 should be isotropic: {p}");
        }
    }

    #[test]
    fn henyey_greenstein_forward_peaks_for_positive_g() {
        // g > 0 → forward scattering (cosθ = 1) is stronger than
        // back-scattering (cosθ = −1).
        let fwd = henyey_greenstein(1.0, 0.7);
        let back = henyey_greenstein(-1.0, 0.7);
        assert!(fwd > back, "g>0 must forward-peak: fwd {fwd} back {back}");
    }

    #[test]
    fn henyey_greenstein_integrates_to_one() {
        // The phase function is a normalised pdf over the sphere — a
        // numeric integral over solid angle must be ≈ 1.
        let g = 0.4f32;
        let n = 2000;
        let mut acc = 0.0f64;
        for i in 0..n {
            // Uniform-cosine quadrature over the sphere.
            let cos = -1.0 + 2.0 * (i as f32 + 0.5) / n as f32;
            // dω = 2π · d(cosθ); d(cosθ) = 2/n.
            acc +=
                henyey_greenstein(cos, g) as f64 * (2.0 * std::f64::consts::PI) * (2.0 / n as f64);
        }
        assert!(
            (acc - 1.0).abs() < 0.02,
            "phase integral {acc} should be ~1"
        );
    }

    #[test]
    fn beer_lambert_attenuation_is_exact_for_a_homogeneous_absorber() {
        // The headline correctness test: a homogeneous pure absorber
        // (scattering albedo 0) attenuates a ray's transmittance by
        // exactly exp(−σ_t · path-length).
        let sigma_t = 1.5f32;
        let medium = Medium::homogeneous(1.0, sigma_t, 0.0);
        let bounds = unit_box();
        // A ray straight through the box along +X: it traverses from
        // x=−1 to x=+1, a path length of 2.
        let ray = Ray::new(vec3(-5.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0));
        let (t0, t1) = bounds.intersect(&ray).expect("ray crosses the box");
        let tr = transmittance(&ray, t0, t1, &medium, &bounds, 0.01);
        let expected = (-sigma_t * 2.0).exp();
        assert!(
            (tr.x - expected).abs() < 1e-3,
            "Beer-Lambert: transmittance {} should be exp(−σ_t·2) = {}",
            tr.x,
            expected
        );
    }

    #[test]
    fn march_through_a_homogeneous_absorber_attenuates_the_background() {
        // A pure absorber (no scattering, no emission) contributes no
        // radiance of its own; its only effect is to attenuate the
        // background by Beer-Lambert.
        let sigma_t = 0.8f32;
        let medium = Medium::homogeneous(1.0, sigma_t, 0.0);
        let bounds = unit_box();
        let ray = Ray::new(vec3(-5.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0));
        let params = VolumeParams {
            step: 0.01,
            ..VolumeParams::default()
        };
        let mut rng = Rng::new(1, 1);
        let result = march_ray(&ray, &medium, &bounds, &[], &params, &mut rng);
        // No emission, no scattering, no lights → zero own radiance.
        assert!(
            result.radiance.max_component() < 1e-5,
            "a pure absorber emits nothing: {:?}",
            result.radiance
        );
        // The background transmittance is exp(−σ_t · 2).
        let expected = (-sigma_t * 2.0).exp();
        assert!(
            (result.transmittance.x - expected).abs() < 1e-2,
            "absorber transmittance {} should be ~{}",
            result.transmittance.x,
            expected
        );
        // Composited over a white background, the pixel is dimmed.
        let composited = result.over(Vec3::ONE);
        assert!(
            (composited.x - expected).abs() < 1e-2,
            "white background should be attenuated to ~{expected}"
        );
    }

    #[test]
    fn an_emissive_medium_adds_glow() {
        // An emissive medium adds radiance even with no external
        // light — the emission-absorption term.
        let dark = Medium::homogeneous(1.0, 0.5, 0.0);
        let glowing = Medium::homogeneous(1.0, 0.5, 0.0).with_emission(vec3(2.0, 2.0, 2.0));
        let bounds = unit_box();
        let ray = Ray::new(vec3(-5.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0));
        let params = VolumeParams {
            step: 0.01,
            ..VolumeParams::default()
        };
        let env = EnvironmentMap::uniform([0.0, 0.0, 0.0]);
        let mut rng = Rng::new(2, 2);
        let dark_r = march_ray(&ray, &dark, &bounds, &[], &params, &mut rng);
        let mut rng = Rng::new(2, 2);
        let glow_r = march_ray(&ray, &glowing, &bounds, &[], &params, &mut rng);
        assert!(
            dark_r.radiance.max_component() < 1e-5,
            "non-emissive medium has no glow"
        );
        assert!(
            glow_r.radiance.x > 0.1,
            "emissive medium must add radiance, got {}",
            glow_r.radiance.x
        );
        let _ = env;
    }

    #[test]
    fn a_zero_density_medium_is_a_no_op() {
        // A medium with density 0 must be perfectly transparent: full
        // transmittance, zero own radiance — an exact no-op.
        let medium = Medium::homogeneous(0.0, 5.0, 1.0).with_emission(vec3(9.0, 9.0, 9.0));
        let bounds = unit_box();
        let ray = Ray::new(vec3(-5.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0));
        let lights = [VolumeLight::point(
            vec3(0.0, 5.0, 0.0),
            vec3(50.0, 50.0, 50.0),
        )];
        let params = VolumeParams::default();
        let mut rng = Rng::new(3, 3);
        let result = march_ray(&ray, &medium, &bounds, &lights, &params, &mut rng);
        assert!(
            result.radiance.max_component() < 1e-6,
            "zero-density medium emits / scatters nothing: {:?}",
            result.radiance
        );
        assert!(
            (result.transmittance.x - 1.0).abs() < 1e-6,
            "zero-density medium is perfectly transparent, got {}",
            result.transmittance.x
        );
        // Compositing over a background returns the background exactly.
        let bg = vec3(0.3, 0.6, 0.9);
        let out = result.over(bg);
        assert!(
            (out.sub(bg)).length() < 1e-6,
            "no-op medium changes nothing"
        );
    }

    #[test]
    fn a_ray_missing_the_box_is_a_no_op() {
        // A ray that never enters the medium box → no radiance, full
        // transmittance.
        let medium = Medium::homogeneous(1.0, 3.0, 0.5).with_emission(vec3(5.0, 5.0, 5.0));
        let bounds = unit_box();
        // A ray well above the box, travelling parallel to it.
        let ray = Ray::new(vec3(-5.0, 10.0, 0.0), vec3(1.0, 0.0, 0.0));
        let params = VolumeParams::default();
        let mut rng = Rng::new(4, 4);
        let result = march_ray(&ray, &medium, &bounds, &[], &params, &mut rng);
        assert!(result.radiance.max_component() < 1e-9);
        assert!((result.transmittance.x - 1.0).abs() < 1e-9);
    }

    #[test]
    fn single_scattering_brightens_a_lit_medium() {
        // A scattering medium with a light present must be brighter
        // than the same medium with no light — single scattering adds
        // in-scattered radiance.
        let medium = Medium::homogeneous(1.0, 1.0, 1.0); // pure scatterer
        let bounds = unit_box();
        let ray = Ray::new(vec3(0.0, 0.0, -5.0), vec3(0.0, 0.0, 1.0));
        let params = VolumeParams {
            step: 0.02,
            ..VolumeParams::default()
        };
        let lit = {
            let lights = [VolumeLight::point(
                vec3(0.0, 5.0, 0.0),
                vec3(40.0, 40.0, 40.0),
            )];
            let mut rng = Rng::new(5, 5);
            march_ray(&ray, &medium, &bounds, &lights, &params, &mut rng).radiance
        };
        let unlit = {
            let mut rng = Rng::new(5, 5);
            march_ray(&ray, &medium, &bounds, &[], &params, &mut rng).radiance
        };
        assert!(unlit.max_component() < 1e-5, "no light → no in-scatter");
        assert!(
            lit.max_component() > 0.01,
            "a lit scattering medium must glow: {lit:?}"
        );
    }

    #[test]
    fn density_grid_trilinear_sampling_recovers_corner_values() {
        // A 2×2×2 grid: trilinear sampling at the 8 corners must return
        // exactly the 8 stored values.
        let values: Vec<f32> = (0..8).map(|i| i as f32).collect();
        let grid = DensityGrid::new([2, 2, 2], values.clone()).unwrap();
        // Corner (0,0,0) → value index 0; corner (1,1,1) → index 7.
        assert!((grid.sample(0.0, 0.0, 0.0) - 0.0).abs() < 1e-6);
        assert!((grid.sample(1.0, 1.0, 1.0) - 7.0).abs() < 1e-6);
        assert!((grid.sample(1.0, 0.0, 0.0) - 1.0).abs() < 1e-6);
        // The centre is the mean of all 8 corners.
        let mean: f32 = values.iter().sum::<f32>() / 8.0;
        assert!((grid.sample(0.5, 0.5, 0.5) - mean).abs() < 1e-5);
    }

    #[test]
    fn density_grid_rejects_a_mismatched_buffer() {
        // A value buffer whose length does not match the dims is
        // rejected — the constructor never builds an invalid grid.
        assert!(DensityGrid::new([2, 2, 2], vec![0.0; 7]).is_none());
        assert!(DensityGrid::new([0, 2, 2], vec![]).is_none());
        assert!(DensityGrid::new([2, 2, 2], vec![1.0; 8]).is_some());
    }

    #[test]
    fn grid_medium_with_zero_density_grid_is_transparent() {
        // A grid medium whose grid is all-zero behaves like a
        // zero-density homogeneous medium — fully transparent.
        let grid = DensityGrid::new([3, 3, 3], vec![0.0; 27]).unwrap();
        let medium = Medium::from_grid(grid, 4.0, 1.0);
        let bounds = unit_box();
        let ray = Ray::new(vec3(-5.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0));
        let params = VolumeParams::default();
        let mut rng = Rng::new(6, 6);
        let result = march_ray(&ray, &medium, &bounds, &[], &params, &mut rng);
        assert!((result.transmittance.x - 1.0).abs() < 1e-5);
        assert!(result.radiance.max_component() < 1e-6);
    }

    #[test]
    fn render_volume_produces_a_framebuffer() {
        // The standalone volume-render entry point produces a
        // framebuffer of the camera's resolution; a lit emissive
        // medium leaves a non-black centre pixel.
        let camera = PtCamera::look_at(
            vec3(0.0, 0.0, 5.0),
            Vec3::ZERO,
            vec3(0.0, 1.0, 0.0),
            50f32.to_radians(),
            16,
            16,
        );
        let medium = Medium::homogeneous(1.0, 1.0, 0.0).with_emission(vec3(3.0, 3.0, 3.0));
        let bounds = unit_box();
        let env = EnvironmentMap::uniform([0.0, 0.0, 0.0]);
        let params = VolumeParams {
            step: 0.05,
            ..VolumeParams::default()
        };
        let fb = render_volume(&camera, &medium, &bounds, &[], &env, &params)
            .expect("render small framebuffer");
        assert_eq!(fb.width, 16);
        assert_eq!(fb.height, 16);
        // The centre pixel looks through the emissive medium → bright.
        let centre = fb.mean(8, 8);
        assert!(
            centre.max_component() > 0.05,
            "centre pixel should see the glowing medium, got {centre:?}"
        );
        // A corner pixel that misses the box sees the black env → dark.
        let corner = fb.mean(0, 0);
        assert!(
            corner.max_component() < 0.05,
            "a ray missing the box should stay dark, got {corner:?}"
        );
    }

    #[test]
    fn denser_medium_attenuates_more() {
        // Monotonicity: a higher density attenuates the background
        // more (a lower transmittance).
        let bounds = unit_box();
        let ray = Ray::new(vec3(-5.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0));
        let params = VolumeParams {
            step: 0.01,
            ..VolumeParams::default()
        };
        let thin = Medium::homogeneous(0.5, 1.0, 0.0);
        let thick = Medium::homogeneous(2.0, 1.0, 0.0);
        let mut rng = Rng::new(7, 7);
        let thin_tr = march_ray(&ray, &thin, &bounds, &[], &params, &mut rng).transmittance;
        let mut rng = Rng::new(7, 7);
        let thick_tr = march_ray(&ray, &thick, &bounds, &[], &params, &mut rng).transmittance;
        assert!(
            thick_tr.x < thin_tr.x,
            "denser medium must attenuate more: thin {} thick {}",
            thin_tr.x,
            thick_tr.x
        );
    }

    /// Round-10 M5 RED→GREEN: pre-fix `render_volume` called
    /// `HdrFramebuffer::new`, which panicked on oversized cameras.
    /// Migrated to `Result<HdrFramebuffer, FramebufferError>` — the
    /// sister fix to the round-9 `tracer::render` migration.
    #[test]
    fn render_volume_returns_too_large_for_oversized_camera() {
        let camera = PtCamera::look_at(
            vec3(0.0, 0.0, 5.0),
            Vec3::ZERO,
            vec3(0.0, 1.0, 0.0),
            50f32.to_radians(),
            100_000,
            100_000,
        );
        let medium = Medium::homogeneous(1.0, 1.0, 0.0);
        let bounds = unit_box();
        let env = EnvironmentMap::uniform([0.0, 0.0, 0.0]);
        let params = VolumeParams::default();
        let err = render_volume(&camera, &medium, &bounds, &[], &env, &params)
            .expect_err("oversized camera must be rejected");
        assert!(matches!(err, FramebufferError::TooLarge { .. }));
    }
}
