//! The Monte-Carlo path-tracing integrator — the renderer's core.
//!
//! # What it computes
//!
//! For every pixel the tracer estimates the **rendering equation**
//!
//! ```text
//!   Lo(x, ωo) = Le(x, ωo) + ∫_Ω fr(x, ωi, ωo)·Li(x, ωi)·(n·ωi) dωi
//! ```
//!
//! by Monte-Carlo integration: it shoots a primary ray through the
//! pixel, and at each surface hit it (a) gathers the emitted radiance,
//! (b) does **next-event estimation** — a direct connection to a
//! sampled point on an emitter — and (c) samples one BRDF direction and
//! recurses. Averaging many such paths per pixel converges, without
//! bias, to the true integral. More samples → less noise; the answer
//! it converges to is correct.
//!
//! # The BRDF
//!
//! Surfaces use the Cook-Torrance microfacet model from
//! [`valenx_render_bridge::pbr`] — GGX normal distribution, Smith
//! geometry term, Fresnel-Schlick — plus a Lambert diffuse lobe.
//! Sampling picks the diffuse or the specular lobe probabilistically
//! (by their relative weight) and importance-samples within it:
//! cosine-weighted for diffuse ([`crate::sampling::cosine_hemisphere`]),
//! GGX-half-vector for specular. The Monte-Carlo estimator weight is
//! `fr·cos / pdf` for the chosen lobe.
//!
//! # Honest scope
//!
//! This is a **real, unbiased v1** — global illumination, importance
//! sampling, next-event estimation, Russian roulette, HDR environment
//! lighting all work and the image converges to the correct answer.
//! It is *not* Cycles-fast and does not ship the production extras:
//! no multiple-importance-sampling weighting between the BRDF and the
//! light sample (so a pure-NEE estimate of a tiny bright light is a
//! little noisier than an MIS renderer), no bidirectional / Metropolis
//! transport, no transmission / refraction lobe (the BRDF is
//! reflection-only), no spectral rendering, no volumetrics, no
//! adaptive sampling, and the render loop is single-threaded. Each of
//! those is a documented, bounded follow-up; none changes the
//! correctness of what is here.

use valenx_render_bridge::pbr::{
    distribution_ggx, f0_from_material, fresnel_schlick, geometry_smith,
};

use crate::framebuffer::{FramebufferError, HdrFramebuffer};
use crate::geometry::{Hit, Ray};
use crate::math::{ortho_basis, vec3, Vec3};
use crate::sampling::{cosine_hemisphere, uniform_disk, Rng};
use crate::scene::{PtMaterial, Scene};

/// Tunable parameters of a render.
#[derive(Clone, Copy, Debug)]
pub struct RenderParams {
    /// Number of independent paths averaged per pixel. More samples →
    /// less Monte-Carlo noise; the variance falls as `1/samples`.
    pub samples_per_pixel: u32,
    /// Hard ceiling on the number of bounces along a path. Russian
    /// roulette usually terminates a path well before this; the cap
    /// only bounds the worst case.
    pub max_depth: u32,
    /// Master random seed — change it for a different noise pattern,
    /// keep it for a bit-exact reproduction.
    pub seed: u64,
    /// Exposure passed to the tone mapper when producing the LDR image.
    pub exposure: f32,
}

impl Default for RenderParams {
    /// A modest default: 64 samples, 8 bounces — enough to see the
    /// global-illumination effect converge on a small image.
    fn default() -> Self {
        RenderParams {
            samples_per_pixel: 64,
            max_depth: 8,
            seed: 0x5eed,
            exposure: 1.0,
        }
    }
}

/// Ray self-intersection guard. A bounce / shadow ray starts this far
/// from the surface along the normal so floating-point error does not
/// make it immediately re-hit the surface it left.
const RAY_EPSILON: f32 = 1e-3;

/// The bounce index at which Russian-roulette path termination begins.
/// The first few bounces always run (they carry most of the energy);
/// from this depth on, a path survives with a probability tied to its
/// throughput, and survivors are scaled up to keep the estimator
/// unbiased.
const RUSSIAN_ROULETTE_START: u32 = 3;

/// Render `scene` into a fresh [`HdrFramebuffer`].
///
/// Shoots `params.samples_per_pixel` jittered primary rays through every
/// pixel, traces each as a full path, and accumulates the result. The
/// returned framebuffer holds raw HDR radiance — call
/// [`HdrFramebuffer::to_ldr`] to get a displayable image.
///
/// The render is deterministic for a fixed `params.seed`: each pixel's
/// RNG stream is seeded from its linear index, so the result does not
/// depend on iteration order (and would be identical if the loop were
/// parallelised).
///
/// # Errors
///
/// Returns [`FramebufferError::TooLarge`] when the scene's camera
/// resolution would allocate a framebuffer larger than the
/// `MAX_FRAMEBUFFER_PIXELS` cap. Round-9 migrated this from a panic
/// path inside `HdrFramebuffer::new` to a fallible result so a
/// hostile scene file can't crash the host process.
pub fn render(scene: &Scene, params: &RenderParams) -> Result<HdrFramebuffer, FramebufferError> {
    let w = scene.camera.width;
    let h = scene.camera.height;
    let mut fb = HdrFramebuffer::try_new(w, h)?;
    let spp = params.samples_per_pixel.max(1);

    for s in 0..spp {
        for y in 0..h {
            for x in 0..w {
                // Per-pixel, per-sample RNG: the stream id is the
                // pixel's linear index, the seed folds in the sample
                // index and the master seed. This makes the render
                // order-independent and exactly reproducible.
                let pixel_index = (y as u64) * (w as u64) + (x as u64);
                let mut rng = Rng::new(
                    params.seed ^ (s as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15),
                    pixel_index,
                );
                let radiance = sample_pixel(scene, params, x, y, &mut rng);
                fb.add_sample(x, y, radiance);
            }
        }
        fb.finish_sample();
    }
    Ok(fb)
}

/// Trace one primary path for pixel `(x, y)` and return its radiance.
///
/// The pixel is jittered by a random sub-pixel offset (box-filter
/// anti-aliasing), the primary ray is generated from the camera frame,
/// and [`trace_path`] integrates the path.
fn sample_pixel(
    scene: &Scene,
    params: &RenderParams,
    x: u32,
    y: u32,
    rng: &mut Rng,
) -> Vec3 {
    let cam = &scene.camera;
    // Jittered normalised image coordinates in [0, 1].
    let jx = rng.next_f32();
    let jy = rng.next_f32();
    let u = (x as f32 + jx) / cam.width as f32;
    // Image-plane v: flip so pixel row 0 is the top of the image.
    let v = 1.0 - (y as f32 + jy) / cam.height as f32;
    let target = cam
        .lower_left
        .add(cam.horizontal.scale(u))
        .add(cam.vertical.scale(v));
    let dir = match target.sub(cam.eye).normalized() {
        Some(d) => d,
        None => return Vec3::ZERO,
    };
    let ray = Ray::new(cam.eye, dir);
    trace_path(scene, params, ray, rng)
}

/// Integrate one full light-transport path, returning the radiance it
/// carries back to the camera.
///
/// The loop maintains a running `throughput` (the product of every
/// BRDF/pdf weight so far) and accumulates `radiance`. At each surface:
///
/// 1. add emitted radiance (the very first hit, plus — to keep the
///    estimate of an emitter hit by a BRDF bounce unbiased — emitter
///    hits along the path);
/// 2. add the **next-event** contribution — a shadow-tested connection
///    to a sampled emitter point, weighted by the BRDF;
/// 3. **sample the BRDF** for the next direction, fold its weight into
///    `throughput`, and continue;
/// 4. apply **Russian roulette** from [`RUSSIAN_ROULETTE_START`] on.
///
/// A ray that escapes all geometry samples the HDR environment and the
/// path ends.
fn trace_path(
    scene: &Scene,
    params: &RenderParams,
    mut ray: Ray,
    rng: &mut Rng,
) -> Vec3 {
    let mut radiance = Vec3::ZERO;
    let mut throughput = Vec3::ONE;
    // `specular_bounce` tracks whether the *previous* bounce was a
    // perfect/near-perfect mirror; after such a bounce next-event
    // estimation was impossible, so an emitter hit must be counted in
    // full. After a diffuse bounce, NEE already accounted for direct
    // light, so an accidental emitter hit is *not* re-added (that would
    // double-count). The primary ray counts as a "specular" predecessor.
    let mut prev_was_specular = true;

    for depth in 0..params.max_depth {
        let hit = scene
            .bvh
            .intersect(&scene.triangles, &ray, RAY_EPSILON, f32::INFINITY);
        let Some(hit) = hit else {
            // Escaped the scene — gather the environment radiance.
            let env = sample_environment(scene, ray.direction);
            radiance = radiance.add(throughput.mul(env));
            break;
        };

        let material = &scene.materials[hit.material];

        // (1) Emitted radiance. Counted on the primary ray, and on a
        // BRDF-bounce hit only when NEE could not have found this
        // emitter (i.e. the predecessor bounce was specular).
        if material.is_emitter() && prev_was_specular {
            radiance = radiance.add(throughput.mul(material.emission));
        }

        // (2) Next-event estimation — direct light from a sampled
        // emitter point, shadow-tested.
        let direct = next_event_estimation(scene, &hit, material, ray.direction, rng);
        radiance = radiance.add(throughput.mul(direct));

        // (3) Sample a new direction from the BRDF.
        let Some(bsdf) = sample_bsdf(material, &hit, ray.direction, rng) else {
            break; // degenerate BRDF sample — terminate the path
        };
        throughput = throughput.mul(bsdf.weight);
        prev_was_specular = bsdf.is_specular;

        // (4) Russian roulette: from a few bounces in, kill low-energy
        // paths probabilistically and boost the survivors so the
        // estimator stays unbiased.
        if depth >= RUSSIAN_ROULETTE_START {
            let survive = throughput.max_component().clamp(0.02, 0.95);
            if rng.next_f32() > survive {
                break;
            }
            throughput = throughput.scale(1.0 / survive);
        }

        // Continue the path from the hit point along the sampled
        // direction, offset off the surface to dodge self-intersection.
        let origin = offset_origin(hit.position, hit.geo_normal, bsdf.direction);
        ray = Ray::new(origin, bsdf.direction);
    }

    radiance
}

/// Sample the HDR environment radiance arriving along `direction`.
fn sample_environment(scene: &Scene, direction: Vec3) -> Vec3 {
    let d = [direction.x, direction.y, direction.z];
    Vec3::from_array(scene.environment.sample_direction(d))
}

/// The result of sampling a surface's BSDF for a new path direction.
struct BsdfSample {
    /// The sampled outgoing (next-bounce) direction, unit length.
    direction: Vec3,
    /// The Monte-Carlo estimator weight `fr·cos / pdf` for this
    /// sample — what the path throughput is multiplied by.
    weight: Vec3,
    /// True if this was a (near-)specular bounce — used to decide
    /// whether a subsequent emitter hit should be counted (see
    /// [`trace_path`]).
    is_specular: bool,
}

/// Sample the Cook-Torrance + Lambert BSDF at a hit point.
///
/// Probabilistically chooses the diffuse or the specular lobe by their
/// relative weight, then importance-samples within it. Returns the
/// sampled direction and the throughput multiplier `fr·cos / pdf`.
///
/// Returns `None` for a degenerate sample (a direction that fell below
/// the surface, a zero pdf) so the caller can terminate the path.
fn sample_bsdf(
    material: &PtMaterial,
    hit: &Hit,
    incoming: Vec3,
    rng: &mut Rng,
) -> Option<BsdfSample> {
    let n = hit.normal;
    // The view direction points *back* toward where the ray came from.
    let v = incoming.neg();
    let n_dot_v = n.dot(v);
    if n_dot_v <= 0.0 {
        // The shading normal is already viewer-facing (the intersection
        // flipped it), so this should not happen; guard anyway.
        return None;
    }

    let mat = &material.pbr;
    let roughness = mat.roughness.clamp(0.03, 1.0);
    let metallic = mat.metallic.clamp(0.0, 1.0);
    let albedo = vec3(
        mat.diffuse_color[0],
        mat.diffuse_color[1],
        mat.diffuse_color[2],
    )
    .scale(1.0 - metallic);
    let f0 = Vec3::from_array(f0_from_material(mat));

    // Lobe-selection probability: weight the specular lobe by its
    // average Fresnel-0 reflectance and the diffuse lobe by its average
    // albedo, then normalise. A pure metal picks specular always; a
    // matte dielectric picks diffuse almost always.
    let diffuse_weight = albedo.max_component();
    let specular_weight = f0.max_component().max(0.04);
    let p_specular = (specular_weight / (specular_weight + diffuse_weight)).clamp(0.05, 0.95);

    if rng.next_f32() < p_specular {
        // --- specular (GGX) lobe ---
        let (dir, half) =
            sample_ggx_direction(n, v, roughness, rng.next_f32(), rng.next_f32())?;
        let n_dot_l = n.dot(dir);
        if n_dot_l <= 0.0 {
            return None; // sample went below the surface
        }
        let n_dot_h = n.dot(half).max(1e-4);
        let v_dot_h = v.dot(half).max(1e-4);

        // Cook-Torrance terms (reuse the shared real-time BRDF code).
        let d = distribution_ggx(n_dot_h, roughness);
        let g = geometry_smith(n_dot_v, n_dot_l, roughness);
        let f = Vec3::from_array(fresnel_schlick(v_dot_h, f0.to_array()));

        // The GGX half-vector sampling pdf, converted from a density
        // over half-vectors to one over light directions:
        //   pdf(l) = D·(n·h) / (4·(v·h))
        let pdf_specular = d * n_dot_h / (4.0 * v_dot_h);
        if pdf_specular <= 1e-8 {
            return None;
        }
        // Cook-Torrance BRDF value: f_spec = D·G·F / (4·(n·v)·(n·l)).
        let denom = 4.0 * n_dot_v * n_dot_l;
        let brdf = f.scale(d * g / denom);
        // Estimator weight = brdf·cos / (pdf · p_lobe).
        let weight = brdf.scale(n_dot_l / (pdf_specular * p_specular));
        Some(BsdfSample {
            direction: dir,
            weight,
            // A low-roughness GGX bounce behaves like a mirror for the
            // emitter-counting decision.
            is_specular: roughness < 0.12,
        })
    } else {
        // --- diffuse (Lambert) lobe ---
        // Cosine-weighted sampling: pdf = cosθ/π and the Lambert BRDF
        // is albedo/π, so brdf·cos/pdf collapses to just `albedo`.
        let dir = cosine_hemisphere(n, rng.next_f32(), rng.next_f32());
        let n_dot_l = n.dot(dir);
        if n_dot_l <= 0.0 {
            return None;
        }
        // The diffuse lobe is scaled by `(1 − F)` energy conservation
        // at normal incidence (so a glossy dielectric does not reflect
        // diffuse + specular > 1).
        let kd = Vec3::ONE.sub(Vec3::from_array(fresnel_schlick(n_dot_v, f0.to_array())));
        let weight = albedo.mul(kd).scale(1.0 / (1.0 - p_specular));
        Some(BsdfSample {
            direction: dir,
            weight,
            is_specular: false,
        })
    }
}

/// Sample a reflected direction from the GGX normal distribution.
///
/// Draws a microfacet half-vector `h` from the GGX distribution about
/// the shading normal, then reflects the view direction about it to get
/// the light direction `l`. Returns `(l, h)`, or `None` if the
/// reflected direction would point into the surface.
fn sample_ggx_direction(
    n: Vec3,
    v: Vec3,
    roughness: f32,
    u1: f32,
    u2: f32,
) -> Option<(Vec3, Vec3)> {
    let alpha = roughness * roughness;
    // GGX half-vector in the local tangent frame: the standard
    // inverse-CDF sampling of the GGX distribution.
    let phi = std::f32::consts::TAU * u1;
    let cos_theta = (((1.0 - u2) / (1.0 + (alpha * alpha - 1.0) * u2)).max(0.0)).sqrt();
    let sin_theta = (1.0 - cos_theta * cos_theta).max(0.0).sqrt();
    let h_local = vec3(sin_theta * phi.cos(), sin_theta * phi.sin(), cos_theta);
    // Rotate the half-vector into world space.
    let (tangent, bitangent) = ortho_basis(n);
    let h = tangent
        .scale(h_local.x)
        .add(bitangent.scale(h_local.y))
        .add(n.scale(h_local.z))
        .normalized()?;
    // Reflect the view about the half-vector → the light direction.
    let l = v.neg().reflect(h);
    if n.dot(l) <= 0.0 {
        None
    } else {
        Some((l, h))
    }
}

/// Next-event estimation — estimate the direct light at a hit point by
/// sampling a point on an emitter and shadow-testing the connection.
///
/// Picks an emitter triangle through the scene's **light tree** (a
/// power × geometric-importance hierarchy — far better than a uniform
/// emitter pick on many-light scenes; see [`crate::light_tree`]),
/// samples a uniform point on it, builds the connection direction,
/// evaluates the surface BRDF for that direction, and — if no geometry
/// blocks the segment — returns the emitter radiance scaled by the
/// BRDF, the two cosines, and the geometric `1/d²` falloff. The whole
/// thing is divided by the area-measure pdf of the light sample,
/// converting it to the solid-angle measure the rendering equation
/// integrates in.
///
/// Returns black if the scene has no emitters or the connection is
/// occluded / back-facing.
fn next_event_estimation(
    scene: &Scene,
    hit: &Hit,
    material: &PtMaterial,
    incoming: Vec3,
    rng: &mut Rng,
) -> Vec3 {
    if scene.emitters.is_empty() {
        return Vec3::ZERO;
    }
    // An emitter does not light itself.
    if material.is_emitter() {
        return Vec3::ZERO;
    }
    // Pick an emitter triangle through the light-tree importance
    // hierarchy.
    let light_sample = match scene.light_tree.sample(hit.position, hit.normal, rng) {
        Some(s) => s,
        None => return Vec3::ZERO,
    };
    let emitter_idx = light_sample.triangle_index as usize;
    let tri = &scene.triangles[emitter_idx];
    let emitter_mat = &scene.materials[tri.material];

    // Uniform point on the triangle via the √-remap barycentric sample.
    let r1 = rng.next_f32();
    let r2 = rng.next_f32();
    let su = r1.sqrt();
    let bary = (1.0 - su, su * (1.0 - r2), su * r2);
    let light_point = tri
        .v0
        .scale(bary.0)
        .add(tri.v1.scale(bary.1))
        .add(tri.v2.scale(bary.2));

    // Direction and distance from the shading point to the light point.
    let to_light = light_point.sub(hit.position);
    let dist2 = to_light.length_sq();
    if dist2 < 1e-8 {
        return Vec3::ZERO;
    }
    let dist = dist2.sqrt();
    let wi = to_light.scale(1.0 / dist);

    let n = hit.normal;
    let n_dot_l = n.dot(wi);
    if n_dot_l <= 0.0 {
        return Vec3::ZERO; // light is below the surface horizon
    }
    // The emitter's geometric normal; the light only radiates from its
    // front face toward the shading point.
    let light_n = tri.geometric_normal();
    let cos_light = light_n.dot(wi.neg());
    if cos_light <= 0.0 {
        return Vec3::ZERO; // back of the emitter faces us
    }

    // Shadow ray: is anything between the shading point and the light?
    let shadow_origin = offset_origin(hit.position, hit.geo_normal, wi);
    let shadow = Ray::new(shadow_origin, wi);
    // Stop the shadow ray just short of the light so the emitter
    // triangle itself does not register as an occluder.
    if scene
        .bvh
        .occluded(&scene.triangles, &shadow, RAY_EPSILON, dist - 2.0 * RAY_EPSILON)
    {
        return Vec3::ZERO;
    }

    // Evaluate the surface BRDF for the connection direction.
    let brdf = evaluate_brdf(material, hit, incoming, wi);

    // Area-measure pdf of the light sample: 1 / (triangle area), times
    // the light-tree selection pdf for the chosen emitter.
    let area = 0.5 * tri.double_area();
    if area < 1e-12 {
        return Vec3::ZERO;
    }
    let pdf_area = light_sample.selection_pdf / area;

    // Convert the area-measure estimate to the solid-angle integral:
    //   contribution = Le · brdf · (n·l) · cos_light / (d² · pdf_area)
    emitter_mat
        .emission
        .mul(brdf)
        .scale(n_dot_l * cos_light / (dist2 * pdf_area))
}

/// Evaluate the Cook-Torrance + Lambert BRDF value `fr(wo, wi)` for an
/// explicit incoming light direction `wi` — used by next-event
/// estimation (which needs the BRDF *value*, not a sample).
fn evaluate_brdf(material: &PtMaterial, hit: &Hit, incoming: Vec3, wi: Vec3) -> Vec3 {
    let n = hit.normal;
    let v = incoming.neg();
    let n_dot_v = n.dot(v).max(1e-4);
    let n_dot_l = n.dot(wi);
    if n_dot_l <= 0.0 {
        return Vec3::ZERO;
    }
    let mat = &material.pbr;
    let roughness = mat.roughness.clamp(0.03, 1.0);
    let metallic = mat.metallic.clamp(0.0, 1.0);
    let albedo = vec3(
        mat.diffuse_color[0],
        mat.diffuse_color[1],
        mat.diffuse_color[2],
    )
    .scale(1.0 - metallic);
    let f0 = Vec3::from_array(f0_from_material(mat));

    // Half-vector between view and light.
    let half = match v.add(wi).normalized() {
        Some(h) => h,
        None => return Vec3::ZERO,
    };
    let n_dot_h = n.dot(half).max(0.0);
    let v_dot_h = v.dot(half).max(0.0);

    // Cook-Torrance specular term.
    let d = distribution_ggx(n_dot_h, roughness);
    let g = geometry_smith(n_dot_v, n_dot_l, roughness);
    let f = Vec3::from_array(fresnel_schlick(v_dot_h, f0.to_array()));
    let spec = f.scale(d * g / (4.0 * n_dot_v * n_dot_l));

    // Energy-conserving Lambert diffuse term: albedo·(1−F)/π.
    let kd = Vec3::ONE.sub(f);
    let diffuse = albedo
        .mul(kd)
        .scale(std::f32::consts::FRAC_1_PI);

    diffuse.add(spec)
}

/// Offset a ray origin a tiny step off a surface along the geometric
/// normal so the new ray cannot immediately re-hit the surface it left.
///
/// The offset is taken toward whichever side the outgoing `direction`
/// goes, so a reflected ray is pushed out and a (hypothetical)
/// transmitted ray would be pushed in.
#[inline]
fn offset_origin(position: Vec3, geo_normal: Vec3, direction: Vec3) -> Vec3 {
    let side = if geo_normal.dot(direction) >= 0.0 {
        1.0
    } else {
        -1.0
    };
    position.add(geo_normal.scale(side * RAY_EPSILON))
}

/// Trace a single primary ray and return the radiance it gathers — a
/// thin convenience wrapper over the internal path integrator for
/// callers / tests that want one ray rather than a whole framebuffer.
pub fn trace_single_ray(scene: &Scene, params: &RenderParams, ray: Ray, rng: &mut Rng) -> Vec3 {
    trace_path(scene, params, ray, rng)
}

/// Sample a point uniformly on the unit disk — re-exported so a caller
/// can build a depth-of-field lens sampler on top of the tracer. The
/// core renderer uses a pinhole camera; lens sampling is a documented
/// follow-up.
pub fn lens_sample(rng: &mut Rng) -> (f32, f32) {
    uniform_disk(rng.next_f32(), rng.next_f32())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scene::{PtCamera, SceneBuilder};
    use valenx_render_bridge::environment::EnvironmentMap;

    /// A camera looking down −Z at the origin.
    fn front_camera(w: u32, h: u32) -> PtCamera {
        PtCamera::look_at(
            vec3(0.0, 0.0, 4.0),
            Vec3::ZERO,
            vec3(0.0, 1.0, 0.0),
            50f32.to_radians(),
            w,
            h,
        )
    }

    /// A large quad in the z = 0 plane facing the camera, of the given
    /// material — fills the view.
    fn add_backdrop(b: &mut SceneBuilder, material: usize) {
        b.add_quad(
            vec3(-10.0, -10.0, 0.0),
            vec3(10.0, -10.0, 0.0),
            vec3(10.0, 10.0, 0.0),
            vec3(-10.0, 10.0, 0.0),
            material,
        );
    }

    #[test]
    fn white_furnace_converges_to_the_environment_colour() {
        // The white-furnace test: a fully-reflective white diffuse
        // surface in a uniform environment must converge to *exactly*
        // the environment radiance — a correct, energy-conserving path
        // tracer neither loses nor invents energy. Albedo 1 means every
        // bounce preserves throughput, so the path eventually escapes
        // to the environment carrying radiance 1.
        let env_radiance = 0.6f32;
        let mut b = SceneBuilder::new(front_camera(16, 16))
            .environment(EnvironmentMap::uniform([env_radiance; 3]));
        let white = b.add_material(PtMaterial::diffuse([1.0, 1.0, 1.0]));
        add_backdrop(&mut b, white);
        let scene = b.build();
        let params = RenderParams {
            samples_per_pixel: 256,
            max_depth: 12,
            seed: 1,
            exposure: 1.0,
        };
        let fb = render(&scene, &params).expect("render small framebuffer");
        // A central pixel (looking straight at the white wall) should
        // read the environment radiance back, within Monte-Carlo noise.
        let c = fb.mean(8, 8);
        assert!(
            (c.x - env_radiance).abs() < 0.05,
            "white-furnace pixel {} should converge to env {}",
            c.x,
            env_radiance
        );
    }

    #[test]
    fn black_surface_in_a_bright_environment_stays_black() {
        // A perfectly absorbing surface reflects nothing — the pixel
        // must be black however bright the environment is. A zero
        // *albedo* alone is not a perfect absorber in a physically
        // based model: a dielectric still has a Fresnel specular lobe
        // (~4 % at the default ior 1.5 — real black plastic is shiny).
        // To express a true absorber the index contrast is removed too
        // (ior = 1 ⇒ F₀ = 0), so this honestly tests that the
        // integrator never invents energy.
        let mut b = SceneBuilder::new(front_camera(8, 8))
            .environment(EnvironmentMap::uniform([5.0, 5.0, 5.0]));
        let mut absorber = PtMaterial::diffuse([0.0, 0.0, 0.0]);
        absorber.pbr.ior = 1.0;
        let black = b.add_material(absorber);
        add_backdrop(&mut b, black);
        let scene = b.build();
        let params = RenderParams {
            samples_per_pixel: 32,
            max_depth: 6,
            seed: 2,
            exposure: 1.0,
        };
        let fb = render(&scene, &params).expect("render small framebuffer");
        let c = fb.mean(4, 4);
        assert!(c.max_component() < 0.02, "absorbing surface should be black");
    }

    #[test]
    fn direct_lighting_matches_the_analytic_value() {
        // A known-direct-lighting case. A small overhead area light of
        // area A and radiance Le illuminates a point on a flat
        // Lambertian floor of albedo ρ directly below it, at distance
        // d. With the light small enough to treat as a point, the
        // irradiance is E = Le·A·cosθ_floor·cosθ_light / d² and the
        // outgoing radiance is L = ρ/π · E. Here the light is directly
        // overhead so both cosines are 1.
        //
        // We disable global illumination noise by using a black
        // environment and reading the direct (next-event) term.
        let albedo = 0.5f32;
        let le = 8.0f32; // emitter radiance
        let light_size = 0.2f32; // small square → near-point
        let d = 3.0f32; // light height above the floor

        let mut b = SceneBuilder::new(PtCamera::look_at(
            vec3(0.0, d - 0.5, 0.001), // just above, looking down
            vec3(0.0, 0.0, 0.0),
            vec3(0.0, 0.0, -1.0),
            40f32.to_radians(),
            8,
            8,
        ))
        .environment(EnvironmentMap::uniform([0.0, 0.0, 0.0]));

        // Floor: a large quad in the y = 0 plane, normal +Y.
        let floor_mat = b.add_material(PtMaterial::diffuse([albedo; 3]));
        b.add_quad(
            vec3(-20.0, 0.0, -20.0),
            vec3(20.0, 0.0, -20.0),
            vec3(20.0, 0.0, 20.0),
            vec3(-20.0, 0.0, 20.0),
            floor_mat,
        );
        // Overhead light: a small quad at y = d, facing down (−Y) so it
        // illuminates the floor beneath it. The winding a→b→c→d below
        // has its geometric normal (b−a)×(c−a) pointing −Y.
        let light_mat = b.add_material(PtMaterial::emissive([le, le, le]));
        let hs = light_size * 0.5;
        b.add_quad(
            vec3(-hs, d, -hs),
            vec3(hs, d, -hs),
            vec3(hs, d, hs),
            vec3(-hs, d, hs),
            light_mat,
        );
        let scene = b.build();
        assert_eq!(scene.emitter_count(), 2, "two light triangles");

        let params = RenderParams {
            samples_per_pixel: 400,
            max_depth: 2, // direct light only (1 bounce)
            seed: 7,
            exposure: 1.0,
        };
        let fb = render(&scene, &params).expect("render small framebuffer");
        let measured = fb.mean(4, 4).x;

        // Analytic point-light approximation of the floor radiance
        // directly under the light.
        let area = light_size * light_size;
        let irradiance = le * area * 1.0 * 1.0 / (d * d); // both cos = 1
        let analytic = albedo * std::f32::consts::FRAC_1_PI * irradiance;

        // The light is small but not a true point; allow a generous
        // tolerance for the finite-area + Monte-Carlo error.
        let rel = (measured - analytic).abs() / analytic;
        assert!(
            rel < 0.25,
            "measured floor radiance {measured} vs analytic {analytic} (rel {rel})"
        );
    }

    #[test]
    fn a_lit_diffuse_surface_is_brighter_than_an_unlit_one() {
        // Qualitative: with an emitter present, the floor directly
        // under it must be brighter than the same floor with no light.
        let build = |with_light: bool| -> f32 {
            let mut b = SceneBuilder::new(PtCamera::look_at(
                vec3(0.0, 2.5, 0.001),
                Vec3::ZERO,
                vec3(0.0, 0.0, -1.0),
                40f32.to_radians(),
                8,
                8,
            ))
            .environment(EnvironmentMap::uniform([0.0, 0.0, 0.0]));
            let floor = b.add_material(PtMaterial::diffuse([0.7; 3]));
            b.add_quad(
                vec3(-20.0, 0.0, -20.0),
                vec3(20.0, 0.0, -20.0),
                vec3(20.0, 0.0, 20.0),
                vec3(-20.0, 0.0, 20.0),
                floor,
            );
            if with_light {
                // Overhead light wound so its normal points −Y (down)
                // toward the floor it lights.
                let light = b.add_material(PtMaterial::emissive([6.0; 3]));
                b.add_quad(
                    vec3(-0.5, 3.0, -0.5),
                    vec3(0.5, 3.0, -0.5),
                    vec3(0.5, 3.0, 0.5),
                    vec3(-0.5, 3.0, 0.5),
                    light,
                );
            }
            let scene = b.build();
            let params = RenderParams {
                samples_per_pixel: 64,
                max_depth: 3,
                seed: 3,
                exposure: 1.0,
            };
            render(&scene, &params).expect("render small framebuffer").mean(4, 4).x
        };
        let lit = build(true);
        let dark = build(false);
        assert!(dark < 0.01, "no light → dark floor, got {dark}");
        assert!(lit > 0.05, "lit floor should be visibly bright, got {lit}");
    }

    #[test]
    fn render_is_deterministic_for_a_fixed_seed() {
        // Same scene, same seed → bit-identical framebuffer.
        let make = || {
            let mut b = SceneBuilder::new(front_camera(8, 8))
                .environment(EnvironmentMap::uniform([0.5; 3]));
            let m = b.add_material(PtMaterial::diffuse([0.8, 0.4, 0.2]));
            add_backdrop(&mut b, m);
            b.build()
        };
        let params = RenderParams {
            samples_per_pixel: 8,
            max_depth: 4,
            seed: 12345,
            exposure: 1.0,
        };
        let a = render(&make(), &params).expect("render small framebuffer");
        let b = render(&make(), &params).expect("render small framebuffer");
        for i in 0..a.accum.len() {
            assert_eq!(a.accum[i], b.accum[i], "render must be reproducible");
        }
    }

    #[test]
    fn a_ray_into_the_void_returns_the_environment() {
        // With no geometry, every primary ray escapes and the pixel is
        // exactly the environment colour.
        let b = SceneBuilder::new(front_camera(4, 4))
            .environment(EnvironmentMap::uniform([0.3, 0.6, 0.9]));
        let scene = b.build();
        let params = RenderParams {
            samples_per_pixel: 16,
            max_depth: 4,
            seed: 9,
            exposure: 1.0,
        };
        let fb = render(&scene, &params).expect("render small framebuffer");
        let c = fb.mean(2, 2);
        assert!((c.x - 0.3).abs() < 1e-4, "empty scene → env colour");
        assert!((c.y - 0.6).abs() < 1e-4);
        assert!((c.z - 0.9).abs() < 1e-4);
    }

    #[test]
    fn indirect_light_bounces_colour_onto_a_neighbour() {
        // The signature global-illumination effect: a bright red wall
        // next to a white wall tints the white wall red via one diffuse
        // bounce. We compare the white wall's red channel with the
        // bounce path enabled (max_depth 4) vs. direct-only
        // (max_depth 1) and require the indirect bounce to add red.
        let render_white_wall = |depth: u32| -> Vec3 {
            let mut b = SceneBuilder::new(PtCamera::look_at(
                vec3(2.0, 0.0, 3.0),
                vec3(-1.0, 0.0, 0.0),
                vec3(0.0, 1.0, 0.0),
                45f32.to_radians(),
                12,
                12,
            ))
            .environment(EnvironmentMap::uniform([0.0, 0.0, 0.0]));
            // White wall at x = -2, facing +X (toward the camera).
            let white = b.add_material(PtMaterial::diffuse([0.9, 0.9, 0.9]));
            b.add_quad(
                vec3(-2.0, -5.0, -5.0),
                vec3(-2.0, -5.0, 5.0),
                vec3(-2.0, 5.0, 5.0),
                vec3(-2.0, 5.0, -5.0),
                white,
            );
            // Bright red wall at x = +2 facing −X (toward the white
            // wall). It is also a mild emitter so there is light to
            // bounce. The winding gives a geometric normal of −X.
            let mut red = PtMaterial::diffuse([0.9, 0.05, 0.05]);
            red.emission = vec3(2.0, 0.1, 0.1);
            let red_id = b.add_material(red);
            b.add_quad(
                vec3(2.0, -5.0, -5.0),
                vec3(2.0, -5.0, 5.0),
                vec3(2.0, 5.0, 5.0),
                vec3(2.0, 5.0, -5.0),
                red_id,
            );
            let scene = b.build();
            let params = RenderParams {
                samples_per_pixel: 200,
                max_depth: depth,
                seed: 4,
                exposure: 1.0,
            };
            render(&scene, &params).expect("render small framebuffer").mean(6, 6)
        };
        // Direct-only: the white wall gets light straight from the red
        // emitter. With more bounces, additional red indirect light is
        // gathered, so the red channel must rise.
        let direct = render_white_wall(1);
        let global = render_white_wall(5);
        assert!(
            global.x > direct.x + 0.01,
            "indirect bounces should add red light: direct {} → global {}",
            direct.x,
            global.x
        );
        // And the bounced light is red-dominated.
        assert!(
            global.x > global.z,
            "the bounced colour should be red-tinted"
        );
    }

    #[test]
    fn trace_single_ray_hits_an_emitter_directly() {
        // A ray fired straight at an emitter returns its emission.
        let mut b = SceneBuilder::new(front_camera(4, 4))
            .environment(EnvironmentMap::uniform([0.0, 0.0, 0.0]));
        let light = b.add_material(PtMaterial::emissive([3.0, 4.0, 5.0]));
        add_backdrop(&mut b, light);
        let scene = b.build();
        let params = RenderParams::default();
        let mut rng = Rng::new(1, 1);
        let ray = Ray::new(vec3(0.0, 0.0, 4.0), vec3(0.0, 0.0, -1.0));
        let r = trace_single_ray(&scene, &params, ray, &mut rng);
        assert!((r.x - 3.0).abs() < 1e-3, "should see the emitter's red");
        assert!((r.z - 5.0).abs() < 1e-3, "should see the emitter's blue");
    }

    /// Round-9 RED→GREEN: `render` used to call `HdrFramebuffer::new`
    /// which would panic when `width * height` exceeded the
    /// `MAX_FRAMEBUFFER_PIXELS` cap. Migrate to the fallible
    /// `try_new` variant and propagate `FramebufferError::TooLarge`
    /// so a hostile scene file (huge camera resolution) can't crash
    /// the host process.
    #[test]
    fn render_returns_too_large_for_giant_camera_instead_of_panicking() {
        let camera = PtCamera::look_at(
            vec3(0.0, 0.0, 4.0),
            Vec3::ZERO,
            vec3(0.0, 1.0, 0.0),
            50f32.to_radians(),
            100_000,
            100_000,
        );
        let mut b = SceneBuilder::new(camera);
        let white = b.add_material(PtMaterial::diffuse([0.8, 0.8, 0.8]));
        add_backdrop(&mut b, white);
        let scene = b.build();
        let params = RenderParams {
            samples_per_pixel: 1,
            max_depth: 1,
            seed: 1,
            exposure: 1.0,
        };
        let err = render(&scene, &params).expect_err("100k x 100k must reject");
        assert!(matches!(err, FramebufferError::TooLarge { .. }));
    }
}
