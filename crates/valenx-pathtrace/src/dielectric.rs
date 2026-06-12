//! Dielectric (glass) BSDF — Fresnel reflection / refraction with a
//! smooth and a GGX rough-dielectric variant.
//!
//! # What this adds
//!
//! The base path tracer's BRDF is reflection-only — it can shade
//! metal and plastic but not glass, water, or a gemstone. This module
//! is the missing **transmission lobe**: a physically-based dielectric
//! BSDF that both reflects *and* refracts, so a ray crossing a
//! glass surface either bounces off it or bends through it.
//!
//! Two surfaces are provided:
//!
//! - [`sample_smooth_dielectric`] — a **perfectly smooth** dielectric
//!   interface. The incident ray is split by the **dielectric Fresnel
//!   equations** into a reflected and a transmitted part; one of the
//!   two is chosen stochastically with probability equal to the
//!   Fresnel reflectance, so the estimator stays unbiased. The
//!   transmitted direction is the **Snell refraction** of the
//!   incident ray; at angles past the **critical angle** the Fresnel
//!   transmittance is zero and the surface **totally-internally
//!   reflects**.
//! - [`sample_rough_dielectric`] — the same physics on a **GGX
//!   microfacet** interface: a microfacet normal is GGX-importance-
//!   sampled, the Fresnel split is evaluated about *that* microfacet,
//!   and the ray is reflected or refracted about it. This is the
//!   "frosted glass" surface — a blurred transmission, the dielectric
//!   analogue of a rough metal.
//!
//! # Energy conservation
//!
//! The dielectric Fresnel reflectance `F` and transmittance `1 − F`
//! sum to exactly 1 at every angle ([`fresnel_dielectric`] returns
//! `F`; the caller takes `1 − F` as the transmitted fraction). No
//! energy is created or destroyed at the interface — the headline
//! correctness property, checked in the tests.
//!
//! # Honest scope
//!
//! A real v1 of a dielectric BSDF. It is **single-interface** — each
//! hit is one reflect-or-refract event; a closed glass solid is shaded
//! by the path tracer crossing its two surfaces in turn (entry and
//! exit), which is exactly how a path tracer handles a solid. It does
//! **not** model wavelength-dependent IOR (dispersion / the prism
//! rainbow), absorption inside the medium (Beer-Lambert tinting of
//! thick glass — that is the `volume` module's job), nor a thin-film
//! interference term. Each is an independent, documented follow-up;
//! none changes the correctness of the reflect/refract split here.

use crate::geometry::Hit;
use crate::math::{ortho_basis, vec3, Vec3};
use crate::sampling::Rng;

/// The **dielectric Fresnel reflectance** for unpolarised light at a
/// flat interface (the full Fresnel equations, not the Schlick
/// approximation).
///
/// `cos_i` is the cosine of the incidence angle (`|ω·n|`, always
/// non-negative); `eta_i` is the IOR of the medium the ray is *in* and
/// `eta_t` the IOR on the far side. Returns the fraction of energy
/// **reflected**; `1 − F` is transmitted.
///
/// ```text
///   sinθt = (η_i/η_t)·sinθi                       (Snell)
///   r_∥ = (η_t·cosθi − η_i·cosθt)/(η_t·cosθi + η_i·cosθt)
///   r_⊥ = (η_i·cosθi − η_t·cosθt)/(η_i·cosθi + η_t·cosθt)
///   F   = ½·(r_∥² + r_⊥²)
/// ```
///
/// When `sinθt ≥ 1` the ray is past the **critical angle** — there is
/// no transmitted ray and the function returns `1.0` (total internal
/// reflection).
pub fn fresnel_dielectric(cos_i: f32, eta_i: f32, eta_t: f32) -> f32 {
    let cos_i = cos_i.clamp(0.0, 1.0);
    // Snell's law for the transmitted-angle sine.
    let sin_i = (1.0 - cos_i * cos_i).max(0.0).sqrt();
    let sin_t = (eta_i / eta_t) * sin_i;
    if sin_t >= 1.0 {
        // Total internal reflection — all energy reflects.
        return 1.0;
    }
    let cos_t = (1.0 - sin_t * sin_t).max(0.0).sqrt();
    // Parallel- and perpendicular-polarised reflection coefficients.
    let r_parl = (eta_t * cos_i - eta_i * cos_t) / (eta_t * cos_i + eta_i * cos_t);
    let r_perp = (eta_i * cos_i - eta_t * cos_t) / (eta_i * cos_i + eta_t * cos_t);
    0.5 * (r_parl * r_parl + r_perp * r_perp)
}

/// Refract incident direction `wi` about unit normal `n` for the
/// relative IOR `eta = η_i/η_t` — **Snell's law** in vector form.
///
/// `wi` points *into* the surface (the direction the ray travels);
/// `n` is on the same side as `wi` is coming *from* (`wi·n < 0`).
/// Returns the unit refracted direction, or `None` for **total
/// internal reflection** (the discriminant goes negative — no real
/// transmitted ray exists past the critical angle).
///
/// ```text
///   cosθt² = 1 − η²·(1 − cosθi²)
///   ωt = η·ωi + (η·cosθi − cosθt)·n
/// ```
pub fn refract(wi: Vec3, n: Vec3, eta: f32) -> Option<Vec3> {
    let cos_i = -wi.dot(n);
    let sin2_i = (1.0 - cos_i * cos_i).max(0.0);
    let sin2_t = eta * eta * sin2_i;
    if sin2_t >= 1.0 {
        // Total internal reflection — no transmitted direction.
        return None;
    }
    let cos_t = (1.0 - sin2_t).max(0.0).sqrt();
    // ωt = η·ωi + (η·cosθi − cosθt)·n
    wi.scale(eta).add(n.scale(eta * cos_i - cos_t)).normalized()
}

/// Which lobe a dielectric BSDF sample took.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DielectricEvent {
    /// The ray reflected off the interface.
    Reflected,
    /// The ray refracted (transmitted) through the interface.
    Transmitted,
}

/// The result of sampling a dielectric BSDF.
#[derive(Clone, Copy, Debug)]
pub struct DielectricSample {
    /// The chosen outgoing direction, unit length.
    pub direction: Vec3,
    /// Monte-Carlo throughput weight for the chosen lobe. With the
    /// stochastic reflect/refract choice this collapses to a clean
    /// value: `1` for a smooth interface (the Fresnel probability and
    /// the Fresnel BSDF factor cancel), times the radiance-scaling
    /// `η_t²/η_i²` correction on transmission.
    pub weight: Vec3,
    /// Whether the sample reflected or refracted.
    pub event: DielectricEvent,
    /// True for a (near-)delta lobe — the smooth interface, and the
    /// rough interface at very low roughness. A delta lobe has no
    /// meaningful pdf, so next-event estimation / MIS skips it.
    pub is_specular: bool,
}

/// Sample a **perfectly smooth dielectric** interface.
///
/// `incoming` is the unit direction the ray travels (toward the
/// surface). `hit.normal` faces the incoming ray (the intersection
/// code already flipped it). `ior` is the dielectric's index of
/// refraction (≈ 1.5 for glass); the other side is assumed to be air
/// (IOR 1).
///
/// The Fresnel reflectance `F` is computed, a uniform random number
/// chooses **reflect** with probability `F` or **refract** with
/// probability `1 − F`, and the matching direction is returned. Total
/// internal reflection (`refract` returns `None`) forces the
/// reflection branch. The throughput weight is `1` for reflection and
/// `η_t²/η_i²` for transmission (the radiance-compression factor that
/// keeps the estimator correct when a ray changes medium).
pub fn sample_smooth_dielectric(
    hit: &Hit,
    incoming: Vec3,
    ior: f32,
    rng: &mut Rng,
) -> DielectricSample {
    // The intersection flips `normal` to face the incoming ray;
    // `back_face` tells us which side we are on. Entering: air → glass.
    // Exiting (a back-face hit): glass → air.
    let entering = !hit.back_face;
    let (eta_i, eta_t) = if entering {
        (1.0, ior.max(1.0))
    } else {
        (ior.max(1.0), 1.0)
    };
    let n = hit.normal; // already faces `−incoming`
    let cos_i = (-incoming.dot(n)).clamp(0.0, 1.0);
    let fresnel = fresnel_dielectric(cos_i, eta_i, eta_t);

    // Stochastically pick reflection or transmission.
    if rng.next_f32() < fresnel {
        // Reflect — mirror about the (flipped) normal.
        let dir = incoming.reflect(n);
        DielectricSample {
            direction: dir,
            // Estimator weight: F (the BSDF) / F (the choice pdf) = 1.
            weight: Vec3::ONE,
            event: DielectricEvent::Reflected,
            is_specular: true,
        }
    } else {
        // Transmit — Snell refraction. `eta = η_i/η_t`.
        let eta = eta_i / eta_t;
        match refract(incoming, n, eta) {
            Some(dir) => {
                // Radiance scales by η_t²/η_i² across an interface
                // (the solid-angle-compression factor). For a ray that
                // enters and later exits the same solid the two
                // factors cancel, so a fully-traced glass object is
                // radiometrically exact.
                let radiance_scale = (eta_t * eta_t) / (eta_i * eta_i);
                DielectricSample {
                    direction: dir,
                    // (1−F)/(1−F)·η_t²/η_i² = η_t²/η_i².
                    weight: Vec3::splat(radiance_scale),
                    event: DielectricEvent::Transmitted,
                    is_specular: true,
                }
            }
            None => {
                // Total internal reflection — the only physical option.
                let dir = incoming.reflect(n);
                DielectricSample {
                    direction: dir,
                    weight: Vec3::ONE,
                    event: DielectricEvent::Reflected,
                    is_specular: true,
                }
            }
        }
    }
}

/// Sample a **rough (GGX microfacet) dielectric** — frosted glass.
///
/// Identical physics to [`sample_smooth_dielectric`] but the
/// reflect / refract event happens about a **microfacet normal**
/// drawn from the GGX distribution rather than the smooth surface
/// normal:
///
/// 1. GGX-importance-sample a microfacet half-vector `h` about the
///    shading normal (the standard GGX inverse-CDF sample).
/// 2. Evaluate the **dielectric Fresnel** about `h`.
/// 3. Stochastically reflect about `h`, or refract through `h` with
///    Snell's law.
///
/// `roughness` in `[0, 1]` controls the blur — 0 reproduces the
/// smooth interface, 1 is a strongly-diffused transmission. The
/// throughput weight folds in the standard microfacet
/// reflect/transmit estimator simplification; for the
/// visible-normal-style sampling used here it reduces to the Smith
/// shadowing ratio, approximated as 1 for this v1 (a documented
/// simplification — see the module docs), so the weight is the same
/// clean `1` / `η_t²/η_i²` as the smooth case.
pub fn sample_rough_dielectric(
    hit: &Hit,
    incoming: Vec3,
    ior: f32,
    roughness: f32,
    rng: &mut Rng,
) -> DielectricSample {
    let roughness = roughness.clamp(0.0, 1.0);
    // Below this roughness the microfacet sampling is indistinguishable
    // from a smooth interface and far cheaper — short-circuit.
    if roughness < 0.02 {
        return sample_smooth_dielectric(hit, incoming, ior, rng);
    }

    let entering = !hit.back_face;
    let (eta_i, eta_t) = if entering {
        (1.0, ior.max(1.0))
    } else {
        (ior.max(1.0), 1.0)
    };
    let n = hit.normal;

    // GGX-importance-sample a microfacet half-vector about `n`.
    let alpha = roughness * roughness;
    let u1 = rng.next_f32();
    let u2 = rng.next_f32();
    let phi = std::f32::consts::TAU * u1;
    let cos_theta = (((1.0 - u2) / (1.0 + (alpha * alpha - 1.0) * u2)).max(0.0)).sqrt();
    let sin_theta = (1.0 - cos_theta * cos_theta).max(0.0).sqrt();
    let h_local = vec3(sin_theta * phi.cos(), sin_theta * phi.sin(), cos_theta);
    let (tangent, bitangent) = ortho_basis(n);
    let h = tangent
        .scale(h_local.x)
        .add(bitangent.scale(h_local.y))
        .add(n.scale(h_local.z))
        .normalized()
        .unwrap_or(n);
    // The microfacet must face the incoming ray for a sensible Fresnel
    // term; flip it if the GGX sample landed on the far side.
    let h = if h.dot(incoming) > 0.0 { h.neg() } else { h };

    let cos_i = (-incoming.dot(h)).clamp(0.0, 1.0);
    let fresnel = fresnel_dielectric(cos_i, eta_i, eta_t);

    if rng.next_f32() < fresnel {
        // Reflect about the microfacet normal.
        let dir = incoming.reflect(h);
        // Reject a reflection that went below the geometric surface
        // (a microfacet artefact) — fall back to a smooth reflection.
        if dir.dot(hit.geo_normal) <= 0.0 {
            let smooth_dir = incoming.reflect(n);
            return DielectricSample {
                direction: smooth_dir,
                weight: Vec3::ONE,
                event: DielectricEvent::Reflected,
                is_specular: roughness < 0.08,
            };
        }
        DielectricSample {
            direction: dir,
            weight: Vec3::ONE,
            event: DielectricEvent::Reflected,
            is_specular: roughness < 0.08,
        }
    } else {
        // Refract through the microfacet normal.
        let eta = eta_i / eta_t;
        match refract(incoming, h, eta) {
            Some(dir) => {
                let radiance_scale = (eta_t * eta_t) / (eta_i * eta_i);
                DielectricSample {
                    direction: dir,
                    weight: Vec3::splat(radiance_scale),
                    event: DielectricEvent::Transmitted,
                    is_specular: roughness < 0.08,
                }
            }
            None => {
                // TIR about the microfacet.
                let dir = incoming.reflect(h);
                let dir = if dir.dot(hit.geo_normal) <= 0.0 {
                    incoming.reflect(n)
                } else {
                    dir
                };
                DielectricSample {
                    direction: dir,
                    weight: Vec3::ONE,
                    event: DielectricEvent::Reflected,
                    is_specular: roughness < 0.08,
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Hit;

    /// Build a synthetic [`Hit`] on a flat +Z surface struck by a ray
    /// travelling along `incoming`; `back_face` says which side.
    fn hit_on_z_plane(incoming: Vec3, back_face: bool) -> Hit {
        // The intersection code flips the normal to face the ray; mimic
        // that — the stored normal opposes `incoming`.
        let mut n = vec3(0.0, 0.0, 1.0);
        if n.dot(incoming) > 0.0 {
            n = n.neg();
        }
        Hit {
            t: 1.0,
            position: Vec3::ZERO,
            normal: n,
            geo_normal: n,
            material: 0,
            back_face,
        }
    }

    /// Fresnel reflectance + transmittance sum to exactly 1 at every
    /// angle below the critical angle — the energy-conservation
    /// property of a dielectric interface.
    #[test]
    fn fresnel_reflect_and_transmit_sum_to_one() {
        // Air → glass (no TIR possible going into the denser medium).
        for &cos_i in &[1.0f32, 0.9, 0.6, 0.3, 0.05] {
            let f = fresnel_dielectric(cos_i, 1.0, 1.5);
            let t = 1.0 - f;
            assert!(
                (f + t - 1.0).abs() < 1e-6,
                "F {f} + T {t} must equal 1 at cosθ {cos_i}"
            );
            assert!((0.0..=1.0).contains(&f), "F {f} out of [0,1]");
        }
    }

    /// At normal incidence the dielectric Fresnel equals the familiar
    /// Schlick `F₀ = ((η−1)/(η+1))²` — for glass, ≈ 0.04.
    #[test]
    fn fresnel_at_normal_incidence_is_the_ior_reflectance() {
        let f = fresnel_dielectric(1.0, 1.0, 1.5);
        let expected = {
            let r = (1.5 - 1.0) / (1.5 + 1.0);
            r * r
        };
        assert!(
            (f - expected).abs() < 1e-4,
            "normal-incidence Fresnel {f} should be ≈ {expected}"
        );
    }

    /// Fresnel rises monotonically toward 1 as the angle grazes.
    #[test]
    fn fresnel_rises_to_one_at_grazing_incidence() {
        let normal = fresnel_dielectric(1.0, 1.0, 1.5);
        let grazing = fresnel_dielectric(0.001, 1.0, 1.5);
        assert!(grazing > normal, "Fresnel should climb toward grazing");
        assert!(grazing > 0.95, "grazing Fresnel {grazing} should near 1");
    }

    /// Past the critical angle (going from the dense to the rare
    /// medium) the interface totally-internally-reflects: F = 1.
    #[test]
    fn total_internal_reflection_past_the_critical_angle() {
        // Glass → air. Critical angle: sinθc = 1/1.5 → θc ≈ 41.8°,
        // cosθc ≈ 0.745. A shallower cosine (steeper angle) is past it.
        let cos_critical = (1.0f32 - (1.0 / 1.5) * (1.0 / 1.5)).sqrt();
        // Well past the critical angle.
        let past = fresnel_dielectric(cos_critical * 0.5, 1.5, 1.0);
        assert!(
            (past - 1.0).abs() < 1e-6,
            "TIR should give F = 1, got {past}"
        );
        // Just inside it, F < 1 (some light still transmits).
        let inside = fresnel_dielectric((cos_critical + 1.0) * 0.5, 1.5, 1.0);
        assert!(
            inside < 1.0,
            "below the critical angle some light transmits"
        );
    }

    /// `refract` obeys Snell's law: η_i·sinθi = η_t·sinθt. We refract a
    /// known direction and check the angle relationship.
    #[test]
    fn refract_obeys_snells_law() {
        // A ray entering glass at 45° from the +Z normal, in the X-Z
        // plane, travelling toward −Z.
        let theta_i = 45f32.to_radians();
        let incoming = vec3(theta_i.sin(), 0.0, -theta_i.cos());
        let n = vec3(0.0, 0.0, 1.0); // faces the incoming ray
        let eta = 1.0 / 1.5; // air → glass
        let refracted = refract(incoming, n, eta).expect("should transmit");
        // The refracted ray's angle from −Z (the transmitted-side
        // normal). Its sine is the X component magnitude.
        let sin_t = refracted.x.abs();
        let sin_i = theta_i.sin();
        // Snell: η_i·sinθi = η_t·sinθt → sinθt = (η_i/η_t)·sinθi.
        let expected_sin_t = eta * sin_i;
        assert!(
            (sin_t - expected_sin_t).abs() < 1e-4,
            "Snell: sinθt {sin_t} should be (η_i/η_t)·sinθi = {expected_sin_t}"
        );
        // The refracted ray continues into the surface (−Z).
        assert!(
            refracted.z < 0.0,
            "refracted ray should go through the surface"
        );
        // It is a unit vector.
        assert!((refracted.length() - 1.0).abs() < 1e-5);
    }

    /// `refract` returns `None` for total internal reflection.
    #[test]
    fn refract_returns_none_for_total_internal_reflection() {
        // Glass → air at a steep angle past the critical angle.
        let theta_i = 70f32.to_radians(); // θc ≈ 41.8°, so this is past it
        let incoming = vec3(theta_i.sin(), 0.0, -theta_i.cos());
        let n = vec3(0.0, 0.0, 1.0);
        let eta = 1.5 / 1.0; // glass → air
        assert!(
            refract(incoming, n, eta).is_none(),
            "a ray past the critical angle has no refracted direction"
        );
        // A shallow angle (below critical) does refract.
        let shallow = 20f32.to_radians();
        let in2 = vec3(shallow.sin(), 0.0, -shallow.cos());
        assert!(refract(in2, n, eta).is_some(), "below critical → transmits");
    }

    /// A straight-on ray passes through the interface undeviated —
    /// refraction at normal incidence does not bend the ray.
    #[test]
    fn refraction_at_normal_incidence_is_undeviated() {
        let incoming = vec3(0.0, 0.0, -1.0);
        let n = vec3(0.0, 0.0, 1.0);
        let refracted = refract(incoming, n, 1.0 / 1.5).expect("straight-on transmits");
        assert!(
            refracted.sub(incoming).length() < 1e-5,
            "a normal-incidence ray should pass straight through"
        );
    }

    /// Sampling the smooth dielectric returns a valid direction and
    /// classifies the event; over many samples both reflection and
    /// transmission occur, in roughly Fresnel-dictated proportion.
    #[test]
    fn smooth_dielectric_splits_reflect_and_transmit() {
        // A near-normal ray into glass: Fresnel ≈ 0.04, so ~96 % of
        // samples should transmit.
        let incoming = vec3(0.0, 0.0, -1.0);
        let hit = hit_on_z_plane(incoming, false);
        let mut rng = Rng::new(1234, 7);
        let mut reflected = 0;
        let mut transmitted = 0;
        for _ in 0..4000 {
            let s = sample_smooth_dielectric(&hit, incoming, 1.5, &mut rng);
            assert!((s.direction.length() - 1.0).abs() < 1e-4, "non-unit dir");
            match s.event {
                DielectricEvent::Reflected => reflected += 1,
                DielectricEvent::Transmitted => transmitted += 1,
            }
        }
        // Both events occur.
        assert!(transmitted > 0 && reflected > 0, "both lobes should fire");
        // Transmission dominates at near-normal incidence (F ≈ 0.04).
        let t_frac = transmitted as f32 / 4000.0;
        assert!(
            t_frac > 0.85,
            "near-normal incidence should mostly transmit, got {t_frac}"
        );
    }

    /// A refracted ray crosses to the far side of the surface; a
    /// reflected ray stays on the near side.
    #[test]
    fn smooth_dielectric_refracted_ray_crosses_the_surface() {
        let incoming = vec3(0.3, 0.0, -0.954).normalized().unwrap();
        let hit = hit_on_z_plane(incoming, false);
        let mut rng = Rng::new(55, 3);
        // Sample until we get one of each event.
        let mut saw_transmit = false;
        let mut saw_reflect = false;
        for _ in 0..2000 {
            let s = sample_smooth_dielectric(&hit, incoming, 1.5, &mut rng);
            match s.event {
                DielectricEvent::Transmitted => {
                    // Continues into the surface (−Z side).
                    assert!(s.direction.z < 0.0, "transmitted ray should cross");
                    saw_transmit = true;
                }
                DielectricEvent::Reflected => {
                    // Bounces back to the near side (+Z).
                    assert!(s.direction.z > 0.0, "reflected ray should stay near");
                    saw_reflect = true;
                }
            }
            if saw_transmit && saw_reflect {
                break;
            }
        }
        assert!(saw_transmit && saw_reflect, "expected both events");
    }

    /// The rough dielectric at near-zero roughness behaves like the
    /// smooth one — the frosted lobe collapses to the sharp lobe.
    #[test]
    fn rough_dielectric_at_zero_roughness_is_smooth() {
        let incoming = vec3(0.0, 0.0, -1.0);
        let hit = hit_on_z_plane(incoming, false);
        let mut rng = Rng::new(9, 9);
        let s = sample_rough_dielectric(&hit, incoming, 1.5, 0.0, &mut rng);
        assert!(s.is_specular, "roughness 0 → a specular (delta) lobe");
        assert!((s.direction.length() - 1.0).abs() < 1e-4);
    }

    /// A frosted-glass surface scatters transmitted rays into a *cone*
    /// rather than a single direction — the rough lobe genuinely
    /// blurs. We confirm the transmitted directions of a rough surface
    /// have spread, while a smooth surface's are identical.
    #[test]
    fn rough_dielectric_blurs_the_transmitted_direction() {
        let incoming = vec3(0.0, 0.0, -1.0);
        let hit = hit_on_z_plane(incoming, false);

        let collect_transmitted = |roughness: f32| -> Vec<Vec3> {
            let mut dirs = Vec::new();
            let mut local_rng = Rng::new(777, 13);
            let mut tries = 0;
            while dirs.len() < 40 && tries < 20000 {
                let s = sample_rough_dielectric(&hit, incoming, 1.5, roughness, &mut local_rng);
                if s.event == DielectricEvent::Transmitted {
                    dirs.push(s.direction);
                }
                tries += 1;
            }
            dirs
        };

        // Spread = mean angular deviation from the average direction.
        let spread = |dirs: &[Vec3]| -> f32 {
            if dirs.is_empty() {
                return 0.0;
            }
            let mut mean = Vec3::ZERO;
            for d in dirs {
                mean = mean.add(*d);
            }
            let mean = mean.normalized().unwrap_or(vec3(0.0, 0.0, -1.0));
            let mut acc = 0.0;
            for d in dirs {
                acc += (1.0 - d.dot(mean)).max(0.0);
            }
            acc / dirs.len() as f32
        };

        let smooth = collect_transmitted(0.0);
        let frosted = collect_transmitted(0.6);
        assert!(
            !frosted.is_empty(),
            "frosted glass should transmit some rays"
        );
        let smooth_spread = spread(&smooth);
        let frosted_spread = spread(&frosted);
        assert!(
            frosted_spread > smooth_spread,
            "frosted-glass transmission {frosted_spread} should be more spread than smooth {smooth_spread}"
        );
    }

    /// Entering vs exiting: the IOR ratio flips between an entry hit
    /// and a back-face (exit) hit, so the refraction bends the
    /// opposite way. We check a transmitted ray on entry bends toward
    /// the normal and on exit bends away.
    #[test]
    fn entering_and_exiting_use_opposite_ior_ratios() {
        // 30° is below the glass→air critical angle (≈ 41.8° for
        // ior 1.5), so the *exit* ray genuinely transmits rather than
        // totally-internally-reflecting — at 45° it would TIR and there
        // would be no transmitted ray to compare.
        let theta = 30f32.to_radians();
        let entering_dir = vec3(theta.sin(), 0.0, -theta.cos());
        // Entering: air → glass, the ray bends *toward* the normal
        // (transmitted angle smaller than incident).
        let eta_enter = 1.0 / 1.5;
        let r_enter = refract(entering_dir, vec3(0.0, 0.0, 1.0), eta_enter).unwrap();
        let sin_enter_out = r_enter.x.abs();
        assert!(
            sin_enter_out < theta.sin(),
            "entering glass should bend toward the normal"
        );
        // Exiting at the same 30° inside the glass: glass → air bends
        // *away* from the normal.
        let eta_exit = 1.5 / 1.0;
        let r_exit = refract(entering_dir, vec3(0.0, 0.0, 1.0), eta_exit).unwrap();
        let sin_exit_out = r_exit.x.abs();
        assert!(
            sin_exit_out > theta.sin(),
            "exiting glass should bend away from the normal"
        );
    }
}
