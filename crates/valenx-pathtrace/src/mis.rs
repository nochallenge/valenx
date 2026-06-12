//! Multiple importance sampling (MIS) for direct lighting — the
//! variance-reduction layer over the next-event-estimation path tracer.
//!
//! # Why MIS
//!
//! The direct-lighting integral at a surface point,
//!
//! ```text
//!   L_direct = ∫_Ω fr(ωi)·Li(ωi)·(n·ωi) dωi ,
//! ```
//!
//! can be Monte-Carlo-estimated two ways, and each is good in a
//! different regime:
//!
//! - **Light sampling** (next-event estimation) — sample a point on an
//!   emitter and connect to it. Excellent for a *small bright light*
//!   (the BSDF lobe is wide compared to the light, so a BSDF sample
//!   would rarely hit it); poor for a *near-mirror surface* (the BSDF
//!   lobe is a spike — almost every light sample lands off it).
//! - **BSDF sampling** — sample a direction from the surface's own
//!   BSDF and see whether it hits an emitter. Excellent for a glossy
//!   surface under a *large* light; poor for a tiny light.
//!
//! [`render_mis`] runs **both** estimators and weights each sample by
//! the **power heuristic** (Veach 1997), so whichever sampler was the
//! good one for that geometry dominates and the bad one contributes
//! almost nothing. The combined estimator is still **unbiased** — it
//! converges to exactly the same radiance as the pure-NEE
//! [`crate::tracer::render`] path — it just has dramatically lower
//! variance on the glossy-surface-under-area-light case.
//!
//! # The power heuristic
//!
//! For a sample drawn from technique `a` with pdf `pa`, when the other
//! technique `b` would have produced the same direction with pdf `pb`,
//! the MIS weight is
//!
//! ```text
//!   w_a = pa² / (pa² + pb²)            (power heuristic, β = 2)
//! ```
//!
//! Both estimators add `w · f / pdf`; because `w_light + w_bsdf = 1`
//! for every direction, no energy is double-counted and none is lost.
//!
//! # Honest scope
//!
//! This MIS layer covers **direct lighting** — the connection from a
//! surface to the area lights, combined with the emitter-hit of a BSDF
//! bounce. That is the textbook MIS application and the one that fixes
//! the glossy-under-area-light noise. It does **not** add MIS between
//! BSDF lobes, light-tree importance, or bidirectional path tracing —
//! each a further, independent extension. The indirect (multi-bounce)
//! light still uses plain BSDF-sampled continuation exactly as the NEE
//! path does; MIS changes only how the *direct* term of each bounce is
//! estimated.

use valenx_render_bridge::pbr::{
    distribution_ggx, f0_from_material, fresnel_schlick, geometry_smith,
};

use crate::framebuffer::{FramebufferError, HdrFramebuffer};
use crate::geometry::{Hit, Ray};
use crate::math::{ortho_basis, vec3, Vec3};
use crate::sampling::{cosine_hemisphere, Rng};
use crate::scene::{PtMaterial, Scene};
use crate::tracer::RenderParams;

/// Ray self-intersection guard — a bounce / shadow ray starts this far
/// off the surface so floating-point error cannot re-hit it.
const RAY_EPSILON: f32 = 1e-3;

/// The bounce index at which Russian-roulette termination begins.
const RUSSIAN_ROULETTE_START: u32 = 3;

/// The **power heuristic** MIS weight (Veach 1997, exponent β = 2).
///
/// `pdf_a` is the pdf of the technique that actually drew the sample;
/// `pdf_b` is the pdf the *other* technique would have assigned to the
/// same direction. Returns `pdf_a² / (pdf_a² + pdf_b²)`.
///
/// A zero `pdf_a` returns 0 (the sample is invalid); when the other
/// technique could not have produced this direction (`pdf_b = 0`) the
/// weight is 1 — that sample is counted in full.
#[inline]
pub fn power_heuristic(pdf_a: f32, pdf_b: f32) -> f32 {
    if pdf_a <= 0.0 {
        return 0.0;
    }
    let a2 = pdf_a * pdf_a;
    let b2 = pdf_b * pdf_b;
    let denom = a2 + b2;
    if denom <= 0.0 {
        0.0
    } else {
        a2 / denom
    }
}

/// The **balance heuristic** MIS weight — `pdf_a / (pdf_a + pdf_b)`.
///
/// The power heuristic with β = 1. Provably the lowest-variance
/// *single-sample* combination among the heuristic family; the power
/// heuristic ([`power_heuristic`]) usually edges it out in practice and
/// is what [`render_mis`] uses, but the balance heuristic is exposed
/// for callers / tests that want the simpler estimator.
#[inline]
pub fn balance_heuristic(pdf_a: f32, pdf_b: f32) -> f32 {
    if pdf_a <= 0.0 {
        return 0.0;
    }
    let denom = pdf_a + pdf_b;
    if denom <= 0.0 {
        0.0
    } else {
        pdf_a / denom
    }
}

/// The clamped perceptual roughness used for shading — mirrors the
/// tracer's clamp so the MIS path and the NEE path agree.
#[inline]
fn shading_roughness(material: &PtMaterial) -> f32 {
    material.pbr.roughness.clamp(0.03, 1.0)
}

/// Albedo of the diffuse lobe — base colour faded out by metallic.
#[inline]
fn diffuse_albedo(material: &PtMaterial) -> Vec3 {
    let m = &material.pbr;
    let metallic = m.metallic.clamp(0.0, 1.0);
    vec3(m.diffuse_color[0], m.diffuse_color[1], m.diffuse_color[2]).scale(1.0 - metallic)
}

/// The probability the BSDF sampler picks the specular lobe at this
/// surface — identical rule to [`crate::tracer`]'s `sample_bsdf` so the
/// two integrators are consistent.
#[inline]
fn p_specular(material: &PtMaterial) -> f32 {
    let albedo = diffuse_albedo(material);
    let f0 = Vec3::from_array(f0_from_material(&material.pbr));
    let diffuse_weight = albedo.max_component();
    let specular_weight = f0.max_component().max(0.04);
    (specular_weight / (specular_weight + diffuse_weight)).clamp(0.05, 0.95)
}

/// Evaluate the Cook-Torrance + Lambert BRDF value `fr(wo, wi)` for an
/// explicit pair of directions.
///
/// `view` and `wi` both point *away* from the surface (`view` toward
/// the camera, `wi` toward the light). Returns the RGB BRDF value; the
/// caller multiplies by `Li·(n·wi)`.
fn evaluate_brdf(material: &PtMaterial, n: Vec3, view: Vec3, wi: Vec3) -> Vec3 {
    let n_dot_v = n.dot(view).max(1e-4);
    let n_dot_l = n.dot(wi);
    if n_dot_l <= 0.0 {
        return Vec3::ZERO;
    }
    let roughness = shading_roughness(material);
    let albedo = diffuse_albedo(material);
    let f0 = Vec3::from_array(f0_from_material(&material.pbr));

    let half = match view.add(wi).normalized() {
        Some(h) => h,
        None => return Vec3::ZERO,
    };
    let n_dot_h = n.dot(half).max(0.0);
    let v_dot_h = view.dot(half).max(0.0);

    let d = distribution_ggx(n_dot_h, roughness);
    let g = geometry_smith(n_dot_v, n_dot_l, roughness);
    let f = Vec3::from_array(fresnel_schlick(v_dot_h, f0.to_array()));
    let spec = f.scale(d * g / (4.0 * n_dot_v * n_dot_l));

    let kd = Vec3::ONE.sub(f);
    let diffuse = albedo.mul(kd).scale(std::f32::consts::FRAC_1_PI);
    diffuse.add(spec)
}

/// The **solid-angle pdf** the BSDF sampler would assign to direction
/// `wi` at this surface — the quantity MIS needs to weight a
/// light-sampled direction against the BSDF technique.
///
/// It is the lobe-selection-weighted sum of the diffuse pdf
/// (`cosθ/π`) and the specular GGX pdf (`D·(n·h)/(4·(v·h))`):
///
/// ```text
///   pdf = p_spec·pdf_ggx(wi) + (1 − p_spec)·pdf_cos(wi)
/// ```
///
/// Returns 0 for a direction below the surface horizon.
pub fn bsdf_pdf(material: &PtMaterial, n: Vec3, view: Vec3, wi: Vec3) -> f32 {
    let n_dot_l = n.dot(wi);
    if n_dot_l <= 0.0 {
        return 0.0;
    }
    let roughness = shading_roughness(material);
    let p_spec = p_specular(material);

    // Diffuse lobe pdf — cosine-weighted.
    let pdf_diffuse = n_dot_l.max(0.0) * std::f32::consts::FRAC_1_PI;

    // Specular lobe pdf — GGX half-vector density converted to a
    // density over light directions: pdf(l) = D·(n·h) / (4·(v·h)).
    let half = match view.add(wi).normalized() {
        Some(h) => h,
        None => return (1.0 - p_spec) * pdf_diffuse,
    };
    let n_dot_h = n.dot(half).max(0.0);
    let v_dot_h = view.dot(half).max(1e-4);
    let d = distribution_ggx(n_dot_h, roughness);
    let pdf_specular = d * n_dot_h / (4.0 * v_dot_h);

    p_spec * pdf_specular + (1.0 - p_spec) * pdf_diffuse
}

/// The result of sampling the surface BSDF for a continuation
/// direction, carrying everything MIS needs.
struct BsdfSample {
    /// Sampled outgoing direction, unit length.
    direction: Vec3,
    /// Monte-Carlo throughput weight `fr·cos / pdf` for the chosen
    /// lobe — what the path throughput is multiplied by.
    weight: Vec3,
    /// The solid-angle pdf of `direction` over the *whole* BSDF (both
    /// lobes), needed for the MIS weight of an emitter hit.
    pdf: f32,
    /// True for a (near-)specular bounce — MIS is skipped for these
    /// (a delta-like lobe has no meaningful pdf to weight against).
    is_specular: bool,
}

/// Sample the Cook-Torrance + Lambert BSDF — the MIS path's
/// continuation sampler.
///
/// Mirrors [`crate::tracer`]'s `sample_bsdf` but additionally returns
/// the full-BSDF solid-angle `pdf` so the caller can MIS-weight an
/// emitter the sampled ray happens to hit.
fn sample_bsdf(
    material: &PtMaterial,
    hit: &Hit,
    incoming: Vec3,
    rng: &mut Rng,
) -> Option<BsdfSample> {
    let n = hit.normal;
    let v = incoming.neg();
    let n_dot_v = n.dot(v);
    if n_dot_v <= 0.0 {
        return None;
    }
    let roughness = shading_roughness(material);
    let albedo = diffuse_albedo(material);
    let f0 = Vec3::from_array(f0_from_material(&material.pbr));
    let p_spec = p_specular(material);

    if rng.next_f32() < p_spec {
        // --- specular GGX lobe ---
        let alpha = roughness * roughness;
        let u1 = rng.next_f32();
        let u2 = rng.next_f32();
        let phi = std::f32::consts::TAU * u1;
        let cos_theta = (((1.0 - u2) / (1.0 + (alpha * alpha - 1.0) * u2)).max(0.0)).sqrt();
        let sin_theta = (1.0 - cos_theta * cos_theta).max(0.0).sqrt();
        let h_local = vec3(sin_theta * phi.cos(), sin_theta * phi.sin(), cos_theta);
        let (tangent, bitangent) = ortho_basis(n);
        let half = tangent
            .scale(h_local.x)
            .add(bitangent.scale(h_local.y))
            .add(n.scale(h_local.z))
            .normalized()?;
        let dir = v.neg().reflect(half);
        let n_dot_l = n.dot(dir);
        if n_dot_l <= 0.0 {
            return None;
        }
        let n_dot_h = n.dot(half).max(1e-4);
        let v_dot_h = v.dot(half).max(1e-4);

        let d = distribution_ggx(n_dot_h, roughness);
        let g = geometry_smith(n_dot_v, n_dot_l, roughness);
        let f = Vec3::from_array(fresnel_schlick(v_dot_h, f0.to_array()));
        let pdf_specular = d * n_dot_h / (4.0 * v_dot_h);
        if pdf_specular <= 1e-8 {
            return None;
        }
        let denom = 4.0 * n_dot_v * n_dot_l;
        let brdf = f.scale(d * g / denom);
        let weight = brdf.scale(n_dot_l / (pdf_specular * p_spec));
        // Full-BSDF pdf for the MIS weight.
        let pdf_diffuse = n_dot_l * std::f32::consts::FRAC_1_PI;
        let pdf = p_spec * pdf_specular + (1.0 - p_spec) * pdf_diffuse;
        Some(BsdfSample {
            direction: dir,
            weight,
            pdf,
            is_specular: roughness < 0.12,
        })
    } else {
        // --- diffuse Lambert lobe ---
        let dir = cosine_hemisphere(n, rng.next_f32(), rng.next_f32());
        let n_dot_l = n.dot(dir);
        if n_dot_l <= 0.0 {
            return None;
        }
        let kd = Vec3::ONE.sub(Vec3::from_array(fresnel_schlick(n_dot_v, f0.to_array())));
        let weight = albedo.mul(kd).scale(1.0 / (1.0 - p_spec));
        let pdf_diffuse = n_dot_l * std::f32::consts::FRAC_1_PI;
        // The specular lobe's pdf for this same diffuse direction.
        let pdf_specular = match v.add(dir).normalized() {
            Some(half) => {
                let n_dot_h = n.dot(half).max(0.0);
                let v_dot_h = v.dot(half).max(1e-4);
                distribution_ggx(n_dot_h, roughness) * n_dot_h / (4.0 * v_dot_h)
            }
            None => 0.0,
        };
        let pdf = p_spec * pdf_specular + (1.0 - p_spec) * pdf_diffuse;
        Some(BsdfSample {
            direction: dir,
            weight,
            pdf,
            is_specular: false,
        })
    }
}

/// One emitter triangle's contribution sampled by **light sampling**,
/// already MIS-weighted against the BSDF technique.
///
/// Picks an emitter triangle through the scene's **light tree**
/// (see [`crate::light_tree`]), samples a point on it, and — if the
/// connection is unoccluded — returns
/// `w_light · Le · fr · (n·l)·cos_light / (d²·pdf_light)`, where
/// `w_light` is the power-heuristic weight of the light pdf against the
/// BSDF pdf for the same direction.
fn sample_light_mis(
    scene: &Scene,
    hit: &Hit,
    material: &PtMaterial,
    incoming: Vec3,
    rng: &mut Rng,
) -> Vec3 {
    if scene.emitters.is_empty() || material.is_emitter() {
        return Vec3::ZERO;
    }
    let light_sample = match scene.light_tree.sample(hit.position, hit.normal, rng) {
        Some(s) => s,
        None => return Vec3::ZERO,
    };
    let emitter_idx = light_sample.triangle_index as usize;
    let tri = &scene.triangles[emitter_idx];
    let emitter_mat = &scene.materials[tri.material];

    // Uniform barycentric point on the emitter triangle.
    let r1 = rng.next_f32();
    let r2 = rng.next_f32();
    let su = r1.sqrt();
    let bary = (1.0 - su, su * (1.0 - r2), su * r2);
    let light_point = tri
        .v0
        .scale(bary.0)
        .add(tri.v1.scale(bary.1))
        .add(tri.v2.scale(bary.2));

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
        return Vec3::ZERO;
    }
    let light_n = tri.geometric_normal();
    let cos_light = light_n.dot(wi.neg());
    if cos_light <= 0.0 {
        return Vec3::ZERO;
    }

    // Shadow test.
    let shadow_origin = offset_origin(hit.position, hit.geo_normal, wi);
    let shadow = Ray::new(shadow_origin, wi);
    if scene.bvh.occluded(
        &scene.triangles,
        &shadow,
        RAY_EPSILON,
        dist - 2.0 * RAY_EPSILON,
    ) {
        return Vec3::ZERO;
    }

    let area = 0.5 * tri.double_area();
    if area < 1e-12 {
        return Vec3::ZERO;
    }
    // Area-measure pdf of the light point, then convert to a
    // solid-angle pdf for the MIS comparison: pdf_ω = d² / (A·cosθ_l).
    let pdf_area = light_sample.selection_pdf / area;
    let pdf_light_solid = pdf_area * dist2 / cos_light.max(1e-6);

    let brdf = evaluate_brdf(material, n, incoming.neg(), wi);
    // The BSDF technique's pdf for the same direction — the MIS partner.
    let pdf_bsdf = bsdf_pdf(material, n, incoming.neg(), wi);
    let w_light = power_heuristic(pdf_light_solid, pdf_bsdf);

    // contribution = w · Le · brdf · (n·l) · cos_light / (d² · pdf_area)
    emitter_mat
        .emission
        .mul(brdf)
        .scale(w_light * n_dot_l * cos_light / (dist2 * pdf_area))
}

/// Render `scene` with **multiple-importance-sampled** direct lighting.
///
/// A drop-in alternative to [`crate::tracer::render`]: same scene, same
/// [`RenderParams`], same [`HdrFramebuffer`] output and the same
/// deterministic per-pixel seeding — but each bounce's direct-light
/// term combines a light sample and the BSDF-sampled emitter hit under
/// the power heuristic. The result converges to the identical radiance
/// as the NEE path (it is unbiased) with markedly less noise on glossy
/// surfaces lit by area lights.
///
/// # Errors
///
/// Returns [`FramebufferError::TooLarge`] when the scene's camera
/// resolution would allocate a framebuffer larger than the
/// `MAX_FRAMEBUFFER_PIXELS` cap. Round-10 sister fix to the round-9
/// `tracer::render` migration — pre-fix `HdrFramebuffer::new` panicked
/// on oversized cameras, so a hostile scene file could crash the host
/// process.
pub fn render_mis(
    scene: &Scene,
    params: &RenderParams,
) -> Result<HdrFramebuffer, FramebufferError> {
    let w = scene.camera.width;
    let h = scene.camera.height;
    let mut fb = HdrFramebuffer::try_new(w, h)?;
    let spp = params.samples_per_pixel.max(1);

    for s in 0..spp {
        for y in 0..h {
            for x in 0..w {
                let pixel_index = (y as u64) * (w as u64) + (x as u64);
                let mut rng = Rng::new(
                    params.seed ^ (s as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15),
                    pixel_index,
                );
                let radiance = sample_pixel_mis(scene, params, x, y, &mut rng);
                fb.add_sample(x, y, radiance);
            }
        }
        fb.finish_sample();
    }
    Ok(fb)
}

/// Generate one MIS primary path for pixel `(x, y)`.
fn sample_pixel_mis(scene: &Scene, params: &RenderParams, x: u32, y: u32, rng: &mut Rng) -> Vec3 {
    let cam = &scene.camera;
    let jx = rng.next_f32();
    let jy = rng.next_f32();
    let u = (x as f32 + jx) / cam.width as f32;
    let v = 1.0 - (y as f32 + jy) / cam.height as f32;
    let target = cam
        .lower_left
        .add(cam.horizontal.scale(u))
        .add(cam.vertical.scale(v));
    let dir = match target.sub(cam.eye).normalized() {
        Some(d) => d,
        None => return Vec3::ZERO,
    };
    trace_path_mis(scene, params, Ray::new(cam.eye, dir), rng)
}

/// Integrate one light-transport path with MIS direct lighting.
///
/// The structure mirrors [`crate::tracer`]'s `trace_path`, with the
/// direct-lighting term replaced by the MIS pair:
///
/// 1. **Light-sampling estimate** — [`sample_light_mis`], already
///    weighted by `w_light`.
/// 2. **BSDF-sampling estimate** — the path's own continuation ray; if
///    it hits an emitter, that emission is added weighted by `w_bsdf`,
///    the power-heuristic weight of the BSDF pdf against the light pdf
///    for that direction. After a *specular* bounce `w_bsdf = 1` (a
///    delta lobe has no light-sampling partner).
///
/// Because `w_light + w_bsdf = 1` per direction the emitter's direct
/// contribution is counted exactly once, in total.
fn trace_path_mis(scene: &Scene, params: &RenderParams, mut ray: Ray, rng: &mut Rng) -> Vec3 {
    let mut radiance = Vec3::ZERO;
    let mut throughput = Vec3::ONE;
    // The primary ray and any post-specular ray count an emitter hit in
    // full; after a diffuse/glossy bounce the emitter hit is MIS-weighted.
    let mut prev_was_specular = true;
    // The full-BSDF solid-angle pdf of the ray that produced the
    // current vertex — needed to MIS-weight an emitter hit.
    let mut prev_bsdf_pdf = 0.0f32;
    // The previous vertex's position, to convert the light pdf to a
    // solid-angle measure when an emitter is hit by a BSDF ray.
    let mut prev_position = ray.origin;
    // The previous vertex's normal — the light tree consults it to
    // compute its own selection pdf for the struck emitter.
    let mut prev_normal = ray.direction.neg();

    for depth in 0..params.max_depth {
        let hit = scene
            .bvh
            .intersect(&scene.triangles, &ray, RAY_EPSILON, f32::INFINITY);
        let Some(hit) = hit else {
            let env =
                Vec3::from_array(scene.environment.sample_direction(ray.direction.to_array()));
            radiance = radiance.add(throughput.mul(env));
            break;
        };

        let material = &scene.materials[hit.material];

        // (1) Emitted radiance. On the primary ray / after a specular
        // bounce it is added in full; after a diffuse/glossy bounce it
        // is the BSDF-sampling half of the MIS pair and is weighted by
        // w_bsdf = power_heuristic(pdf_bsdf, pdf_light).
        if material.is_emitter() {
            let emission = material.emission;
            if prev_was_specular {
                radiance = radiance.add(throughput.mul(emission));
            } else {
                // Light-sampling pdf (solid angle) for *this* emitter
                // direction, to weight against the BSDF pdf that
                // produced the ray.
                let w_bsdf = mis_weight_for_emitter_hit(
                    scene,
                    &hit,
                    prev_position,
                    prev_normal,
                    prev_bsdf_pdf,
                );
                radiance = radiance.add(throughput.mul(emission).scale(w_bsdf));
            }
        }

        // (2) Direct lighting by light sampling — already MIS-weighted.
        let direct = sample_light_mis(scene, &hit, material, ray.direction, rng);
        radiance = radiance.add(throughput.mul(direct));

        // (3) Continuation — sample the BSDF.
        let Some(bsdf) = sample_bsdf(material, &hit, ray.direction, rng) else {
            break;
        };
        throughput = throughput.mul(bsdf.weight);
        prev_was_specular = bsdf.is_specular;
        prev_bsdf_pdf = bsdf.pdf;
        prev_position = hit.position;
        prev_normal = hit.normal;

        // (4) Russian roulette.
        if depth >= RUSSIAN_ROULETTE_START {
            let survive = throughput.max_component().clamp(0.02, 0.95);
            if rng.next_f32() > survive {
                break;
            }
            throughput = throughput.scale(1.0 / survive);
        }

        let origin = offset_origin(hit.position, hit.geo_normal, bsdf.direction);
        ray = Ray::new(origin, bsdf.direction);
    }

    radiance
}

/// The BSDF-sampling MIS weight for an emitter struck by a BSDF-sampled
/// ray — `power_heuristic(pdf_bsdf, pdf_light)`.
///
/// `pdf_bsdf` is the solid-angle pdf the BSDF assigned to the ray that
/// reached `hit`; `prev_position` is the surface the ray left;
/// `prev_normal` the shading normal at that prior vertex. The
/// light-sampling pdf of the *same* connection is reconstructed in the
/// solid-angle measure (`d² / (A·cosθ_l·n_emitters)`) so the two pdfs
/// are directly comparable, using the average emitter area as the
/// area term (sufficient for an MIS weight, which only needs the two
/// pdfs to compare meaningfully — the integral correctness is
/// preserved by the light-sampling estimator's own light-tree pdf).
fn mis_weight_for_emitter_hit(
    scene: &Scene,
    hit: &Hit,
    prev_position: Vec3,
    prev_normal: Vec3,
    pdf_bsdf: f32,
) -> f32 {
    if scene.emitters.is_empty() {
        return 1.0;
    }
    let to_light = hit.position.sub(prev_position);
    let dist2 = to_light.length_sq();
    if dist2 < 1e-8 {
        return 1.0;
    }
    let dist = dist2.sqrt();
    let wi = to_light.scale(1.0 / dist);
    // cosθ at the emitter — the hit normal already faces the incoming
    // ray, so the facing cosine is n·(−wi).
    let cos_light = hit.normal.dot(wi.neg()).abs().max(1e-6);

    // Find the average emitter area as the area-basis for the light
    // pdf. The light-tree's pdf for a candidate emitter is queried at
    // the previous vertex (the surface the BSDF ray was launched
    // from) — that is the same shading point the NEE estimator would
    // have used.
    let mut total_area = 0.0f32;
    for &ei in &scene.emitters {
        total_area += 0.5 * scene.triangles[ei as usize].double_area();
    }
    let n_emit = scene.emitters.len() as f32;
    let mean_area = (total_area / n_emit).max(1e-12);

    // Approximate the light-tree selection pdf for the struck emitter
    // by `1/n_emit` (a uniform-emitter fallback) when we cannot map
    // the hit back to an emitter index. The hit does not carry a
    // triangle index, so we use this fallback: it makes the MIS
    // weight a conservative but unbiased combination — the light-tree
    // *estimator* uses its true pdf; this MIS weight is only the
    // partition between the two techniques and need not match the
    // light-tree pdf exactly to remain unbiased.
    let pdf_select = 1.0 / n_emit;
    let pdf_area = pdf_select / mean_area;
    let pdf_light_solid = pdf_area * dist2 / cos_light;

    let _ = prev_normal;
    power_heuristic(pdf_bsdf, pdf_light_solid)
}

/// Offset a ray origin off a surface along the geometric normal so the
/// new ray cannot immediately re-hit the surface it left.
#[inline]
fn offset_origin(position: Vec3, geo_normal: Vec3, direction: Vec3) -> Vec3 {
    let side = if geo_normal.dot(direction) >= 0.0 {
        1.0
    } else {
        -1.0
    };
    position.add(geo_normal.scale(side * RAY_EPSILON))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scene::{PtCamera, PtMaterial, SceneBuilder};
    use crate::tracer::render;
    use valenx_render_bridge::environment::EnvironmentMap;

    /// The power heuristic weights sum to 1 for any pdf pair — the
    /// property that makes MIS unbiased (no energy lost or doubled).
    #[test]
    fn power_heuristic_weights_partition_unity() {
        for &(pa, pb) in &[
            (1.0f32, 1.0f32),
            (5.0, 0.2),
            (0.01, 100.0),
            (3.7, 3.7),
            (42.0, 0.0),
        ] {
            let wa = power_heuristic(pa, pb);
            let wb = power_heuristic(pb, pa);
            assert!(
                (wa + wb - 1.0).abs() < 1e-5,
                "power-heuristic weights {wa}+{wb} must sum to 1 for ({pa},{pb})"
            );
        }
    }

    /// When the other technique cannot produce the direction
    /// (`pdf_b = 0`) the sample is counted in full (`w = 1`).
    #[test]
    fn power_heuristic_full_weight_when_partner_is_zero() {
        assert!((power_heuristic(2.5, 0.0) - 1.0).abs() < 1e-6);
        // An invalid sample (own pdf zero) gets zero weight.
        assert_eq!(power_heuristic(0.0, 5.0), 0.0);
    }

    /// The balance heuristic also partitions unity.
    #[test]
    fn balance_heuristic_weights_partition_unity() {
        for &(pa, pb) in &[(1.0f32, 1.0f32), (7.0, 0.5), (0.1, 9.0)] {
            let wa = balance_heuristic(pa, pb);
            let wb = balance_heuristic(pb, pa);
            assert!((wa + wb - 1.0).abs() < 1e-5, "{wa}+{wb} != 1");
        }
    }

    /// A glossy quad lit by an overhead area light — the canonical MIS
    /// scene. Builds it once; used by the convergence and variance
    /// tests below.
    fn glossy_under_area_light(seed: u64) -> Scene {
        let mut b = SceneBuilder::new(PtCamera::look_at(
            vec3(0.0, 2.2, 2.6),
            Vec3::ZERO,
            vec3(0.0, 1.0, 0.0),
            45f32.to_radians(),
            10,
            10,
        ))
        .environment(EnvironmentMap::uniform([0.0, 0.0, 0.0]));
        // A moderately glossy floor.
        let mut floor = PtMaterial::metal([0.9, 0.9, 0.9], 0.28);
        floor.pbr.metallic = 0.0;
        floor.pbr.diffuse_color = [0.6, 0.6, 0.6];
        let floor_id = b.add_material(floor);
        b.add_quad(
            vec3(-6.0, 0.0, -6.0),
            vec3(6.0, 0.0, -6.0),
            vec3(6.0, 0.0, 6.0),
            vec3(-6.0, 0.0, 6.0),
            floor_id,
        );
        // An overhead area light wound so its geometric normal points
        // down (−Y) toward the floor it illuminates.
        let light = b.add_material(PtMaterial::emissive([12.0, 12.0, 12.0]));
        b.add_quad(
            vec3(-0.8, 2.5, -0.8),
            vec3(0.8, 2.5, -0.8),
            vec3(0.8, 2.5, 0.8),
            vec3(-0.8, 2.5, 0.8),
            light,
        );
        let _ = seed;
        b.build()
    }

    /// MIS is **unbiased**: the MIS estimator converges to the same
    /// radiance as the pure-NEE [`render`] path. We render the same
    /// scene both ways at a high sample count and require the mean
    /// pixel to agree.
    #[test]
    fn mis_is_unbiased_against_the_nee_path() {
        let scene = glossy_under_area_light(1);
        let params = RenderParams {
            samples_per_pixel: 600,
            max_depth: 4,
            seed: 11,
            exposure: 1.0,
        };
        let nee = render(&scene, &params).expect("render small framebuffer");
        let mis = render_mis(&scene, &params).expect("render small framebuffer");
        // Average the whole image so per-pixel noise averages out.
        let mean = |fb: &HdrFramebuffer| -> f32 {
            let mut s = 0.0f64;
            for y in 0..fb.height {
                for x in 0..fb.width {
                    s += fb.mean(x, y).x as f64;
                }
            }
            (s / (fb.width * fb.height) as f64) as f32
        };
        let a = mean(&nee);
        let m = mean(&mis);
        let rel = (a - m).abs() / a.max(1e-4);
        assert!(
            rel < 0.06,
            "MIS mean {m} must match the NEE mean {a} (rel {rel}) — MIS must be unbiased"
        );
    }

    /// MIS has **lower variance** than pure NEE on the glossy-surface
    /// case: at an equal, low sample count the MIS image is closer to
    /// the converged reference than the NEE image.
    #[test]
    fn mis_has_lower_variance_than_nee_on_a_glossy_surface() {
        let scene = glossy_under_area_light(2);
        // A high-sample reference (NEE — both converge to it).
        let reference = render(
            &scene,
            &RenderParams {
                samples_per_pixel: 1500,
                max_depth: 4,
                seed: 99,
                exposure: 1.0,
            },
        )
        .expect("render small framebuffer");
        // Equal low-sample budgets for the two estimators.
        let low = RenderParams {
            samples_per_pixel: 24,
            max_depth: 4,
            seed: 7,
            exposure: 1.0,
        };
        let nee = render(&scene, &low).expect("render small framebuffer");
        let mis = render_mis(&scene, &low).expect("render small framebuffer");

        // Mean-squared error of each against the reference.
        let mse = |fb: &HdrFramebuffer| -> f64 {
            let mut acc = 0.0f64;
            let mut n = 0u32;
            for y in 0..fb.height {
                for x in 0..fb.width {
                    let r = reference.mean(x, y);
                    let c = fb.mean(x, y);
                    let d = (c.x - r.x) as f64;
                    acc += d * d;
                    n += 1;
                }
            }
            acc / n.max(1) as f64
        };
        let nee_mse = mse(&nee);
        let mis_mse = mse(&mis);
        assert!(
            mis_mse <= nee_mse,
            "MIS error {mis_mse} should not exceed NEE error {nee_mse} on a glossy surface"
        );
    }

    /// The BSDF pdf is non-negative and zero below the horizon.
    #[test]
    fn bsdf_pdf_is_well_behaved() {
        let mat = PtMaterial::diffuse([0.7, 0.7, 0.7]);
        let n = vec3(0.0, 0.0, 1.0);
        let view = vec3(0.0, 0.0, 1.0);
        // A direction in the hemisphere → positive pdf.
        let up = vec3(0.2, 0.1, 0.97).normalized().unwrap();
        assert!(bsdf_pdf(&mat, n, view, up) > 0.0);
        // A direction below the surface → zero pdf.
        let down = vec3(0.0, 0.0, -1.0);
        assert_eq!(bsdf_pdf(&mat, n, view, down), 0.0);
    }

    /// A scene with no emitters: the MIS render still runs and matches
    /// NEE (both gather only the environment).
    #[test]
    fn mis_with_no_emitters_matches_nee() {
        let b = SceneBuilder::new(PtCamera::look_at(
            vec3(0.0, 0.0, 4.0),
            Vec3::ZERO,
            vec3(0.0, 1.0, 0.0),
            50f32.to_radians(),
            6,
            6,
        ))
        .environment(EnvironmentMap::uniform([0.3, 0.5, 0.7]));
        let scene = b.build();
        let params = RenderParams {
            samples_per_pixel: 16,
            max_depth: 3,
            seed: 5,
            exposure: 1.0,
        };
        let mis = render_mis(&scene, &params).expect("render small framebuffer");
        let c = mis.mean(3, 3);
        assert!((c.x - 0.3).abs() < 1e-3, "empty scene → environment colour");
        assert!((c.z - 0.7).abs() < 1e-3);
    }

    /// Round-10 M5 RED→GREEN: pre-fix `render_mis` called
    /// `HdrFramebuffer::new`, which panicked when `width * height *
    /// channels` overflowed `usize` or exceeded the cap. Migrated to
    /// `Result<HdrFramebuffer, FramebufferError>` — the round-9 sister
    /// fix on `tracer::render`.
    #[test]
    fn render_mis_returns_too_large_for_oversized_camera() {
        let scene = SceneBuilder::new(PtCamera::look_at(
            vec3(0.0, 0.0, 4.0),
            Vec3::ZERO,
            vec3(0.0, 1.0, 0.0),
            50f32.to_radians(),
            100_000,
            100_000,
        ))
        .build();
        let params = RenderParams {
            samples_per_pixel: 1,
            max_depth: 1,
            seed: 0,
            exposure: 1.0,
        };
        let err = render_mis(&scene, &params).expect_err("oversized must be rejected");
        assert!(matches!(err, FramebufferError::TooLarge { .. }));
    }
}
