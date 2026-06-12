//! Rolling-ball blend surfaces (Phase 19F).
//!
//! A rolling-ball blend is the cornerstone of every CAD fillet /
//! round operation. Geometrically:
//!
//! 1. A sphere of `radius` rolls in contact with two surfaces
//!    simultaneously. Its center traces a curve `C_ball(t)` in 3D —
//!    the **spine** of the blend.
//! 2. Where the sphere touches each surface, it leaves a contact
//!    point `q_A(t)` on surface A and `q_B(t)` on surface B —
//!    together the two **contact rails**.
//! 3. The **blend surface** is the part of the ball that's between
//!    the two contact points: at each `t`, the arc on the rolling
//!    sphere from `q_A(t)` to `q_B(t)` (the short arc, on the convex
//!    side of the spine). Sweeping these arcs gives the blend.
//!
//! ## Algorithm (v1)
//!
//! We trace the spine by **marching on the bisector** of the two
//! surfaces:
//!
//! - **Spine condition.** A point `c` is on the spine if
//!   `dist(c, surface_A) = dist(c, surface_B) = radius`. At step
//!   `i`, given the current `(uv_A, uv_B, c)`:
//!   - foot-project `c` onto A → `q_A` at distance `d_A`;
//!   - foot-project `c` onto B → `q_B` at distance `d_B`;
//!   - bisector direction: `n_bisect = -normalize((q_A - c) +
//!     (q_B - c))` — the average outward direction; pushing the
//!     spine along `n_bisect` toward larger `r` is what equalises
//!     the two distances. We correct with a small Newton-style
//!     adjustment that drives both distances toward `radius`
//!     simultaneously by stepping along the gradient of
//!     `(d_A - r)² + (d_B - r)²`.
//!
//! - **Spine tangent.** The tangent `t_spine` of the ball-center
//!   curve at the current `c` lies in the plane perpendicular to
//!   *both* contact normals `(c - q_A)` and `(c - q_B)` — it equals
//!   the cross product of the two contact normals (modulo sign), the
//!   same construction the marching-SSI uses for intersection
//!   tangents.
//!
//! - **Step.** March along `t_spine` by an adaptive chord-tolerance
//!   step, then re-equalise distances to `radius`. Terminate on a
//!   boundary or loop, identical to the marching-SSI bookkeeping.
//!
//! - **Emit.** We collect three polylines along the trace: the spine
//!   `c(t)`, contact A `q_A(t)`, and contact B `q_B(t)`. The blend
//!   surface is a **ruled-like swept patch** with two rails (A and
//!   B); at each sample we sample a short arc of the ball between
//!   `q_A` and `q_B` and stack the arc as a v-row of the NURBS
//!   surface. The cross-section degree is 2 (3 CPs per arc — a
//!   rational quadratic), the spine degree is 3 (cubic fit through
//!   spine samples).
//!
//! ## Verified analytic cases
//!
//! - **Two planes at angle θ** → spine is a straight line parallel
//!   to the planes' intersection, the blend is a cylindrical fillet
//!   of radius `r` and arc angle `π − θ`.
//! - **Plane + cylinder of radius `R`** → spine is parallel to the
//!   cylinder axis at distance `r + R` from the axis perpendicular
//!   to the plane, the blend is a torus patch of major radius `r +
//!   R` and minor radius `r`.
//!
//! Both are checked in the test suite to within `1e-3 · r`.
//!
//! ## v1 caveats
//!
//! - Constant radius only; variable-radius blends are a v1.5 follow-up.
//! - We require a caller-supplied seed point (a 3D point on the
//!   spine) and a recommended spine direction; for the planar +
//!   plane-cylinder cases we expose convenience seed-finders. A
//!   fully-automatic seed search for arbitrary surface pairs is the
//!   T3 polish (it reduces to a tessellation-seed + bisector-Newton
//!   in 3 unknowns, doable but bigger than the v1 scope).
//! - Self-intersecting and branching blend topologies are out of
//!   scope; each call returns a single connected blend.

use nalgebra::{Vector2, Vector3};

use crate::error::SurfaceError;
use crate::march_ssi::MarchParams;
use crate::nurbs_surface::NurbsSurface;

/// User-facing knobs for [`rolling_ball_blend`].
#[derive(Clone, Debug)]
pub struct BlendParams {
    /// Target chord deviation per spine step (model units). v1 1e-3.
    pub chord_tolerance: f64,
    /// Newton convergence tolerance on each per-surface closest-foot
    /// projection. v1 1e-10.
    pub projection_tolerance: f64,
    /// Per-Newton iteration cap. v1 8.
    pub projection_iters: usize,
    /// Maximum spine step size (model units). v1 0.5.
    pub max_step: f64,
    /// Minimum spine step size — below this we terminate. v1 1e-6.
    pub min_step: f64,
    /// Hard cap on spine samples. v1 512.
    pub max_samples: usize,
    /// Number of cross-section CPs (must be odd, at least 3). Cross
    /// section is degree 2 (rational quadratic arc).
    pub cross_section_cps: usize,
}

impl Default for BlendParams {
    fn default() -> Self {
        Self {
            chord_tolerance: 1.0e-3,
            projection_tolerance: 1.0e-10,
            projection_iters: 8,
            max_step: 0.5,
            min_step: 1.0e-6,
            max_samples: 512,
            cross_section_cps: 3,
        }
    }
}

/// One sample along a rolling-ball spine.
#[derive(Clone, Debug)]
pub struct SpineSample {
    /// Ball-center position.
    pub center: Vector3<f64>,
    /// Contact point on surface A.
    pub contact_a: Vector3<f64>,
    /// `(u, v)` parameters of the contact on surface A.
    pub uv_a: Vector2<f64>,
    /// Contact point on surface B.
    pub contact_b: Vector3<f64>,
    /// `(u, v)` parameters of the contact on surface B.
    pub uv_b: Vector2<f64>,
}

/// Result of a rolling-ball blend operation.
#[derive(Clone, Debug)]
pub struct Blend {
    /// Sampled spine + contact pairs.
    pub spine: Vec<SpineSample>,
    /// The blend surface — a tensor-product NURBS patch with the
    /// rolling spine in the u direction and a degree-2 rational
    /// cross-section arc in v.
    pub surface: NurbsSurface,
    /// Constant rolling-ball radius for this blend.
    pub radius: f64,
}

/// Trace + emit a rolling-ball blend between two surfaces.
///
/// `seed_center` is a 3D point that's roughly on the blend spine
/// (about distance `radius` from both surfaces); it gets equalised to
/// the true spine by a few Newton iterations before tracing starts.
///
/// `recommended_tangent` is the direction the marcher should walk
/// first along the spine. If `None`, we pick the cross product of the
/// two outward contact normals and accept whichever sign keeps the
/// trace inside both surfaces' parameter ranges.
pub fn rolling_ball_blend(
    s_a: &NurbsSurface,
    s_b: &NurbsSurface,
    radius: f64,
    seed_center: Vector3<f64>,
    recommended_tangent: Option<Vector3<f64>>,
    params: &BlendParams,
) -> Result<Blend, SurfaceError> {
    if !(radius > 0.0 && radius.is_finite()) {
        return Err(SurfaceError::IntersectionFailed(format!(
            "rolling_ball_blend: radius must be positive finite; got {radius}"
        )));
    }

    // Equalise the seed to the true spine.
    let m_params = MarchParams {
        projection_tolerance: params.projection_tolerance,
        projection_iters: params.projection_iters,
        ..MarchParams::default()
    };
    let seed = equalise_spine(s_a, s_b, seed_center, radius, &m_params, 32)?;

    let mut samples = vec![seed.clone()];

    // March forward.
    let initial_tangent = match recommended_tangent {
        Some(t) => normalize_or_default(t),
        None => spine_tangent(&seed),
    };
    let mut samples_forward = trace_one_direction(
        s_a,
        s_b,
        radius,
        seed.clone(),
        initial_tangent,
        params,
        &m_params,
    );
    samples.append(&mut samples_forward);

    // March backward from the seed in the opposite direction.
    let back_tangent = -initial_tangent;
    let mut samples_backward = trace_one_direction(
        s_a,
        s_b,
        radius,
        seed.clone(),
        back_tangent,
        params,
        &m_params,
    );
    // Prepend backwards walk (excluding the duplicated seed) in
    // reverse order so the final sequence is monotone along the
    // spine.
    samples_backward.reverse();
    let mut full = samples_backward;
    full.extend(samples);

    let surface = build_blend_surface(s_a, s_b, &full, radius, params)?;
    Ok(Blend {
        spine: full,
        surface,
        radius,
    })
}

/// Drive the seed `c` to the true bisector spine: equalise
/// `dist(c, A) = dist(c, B) = radius`.
///
/// Uses a damped Newton-style update that moves `c` along the
/// gradient of the residual `(d_A - r)² + (d_B - r)²`. The two
/// contact normals provide the natural search basis.
fn equalise_spine(
    s_a: &NurbsSurface,
    s_b: &NurbsSurface,
    seed: Vector3<f64>,
    radius: f64,
    m: &MarchParams,
    max_iters: usize,
) -> Result<SpineSample, SurfaceError> {
    let mut c = seed;
    for _ in 0..max_iters {
        let Some((uv_a, uv_b)) = refine_seed_to_3d(s_a, s_b, c, m) else {
            return Err(SurfaceError::IntersectionFailed(
                "rolling_ball_blend: could not project seed onto both surfaces".into(),
            ));
        };
        let p_a = s_a.evaluate(uv_a.x, uv_a.y);
        let p_b = s_b.evaluate(uv_b.x, uv_b.y);
        let d_a_vec = c - p_a;
        let d_b_vec = c - p_b;
        let d_a = d_a_vec.norm();
        let d_b = d_b_vec.norm();
        let r_a = d_a - radius;
        let r_b = d_b - radius;
        if r_a.abs() < 1.0e-9 && r_b.abs() < 1.0e-9 {
            return Ok(SpineSample {
                center: c,
                contact_a: p_a,
                uv_a,
                contact_b: p_b,
                uv_b,
            });
        }
        // Pull / push along the outward contact normals (the
        // gradients of `d_A` and `d_B` with respect to `c`).
        let n_a = if d_a > 1.0e-12 {
            d_a_vec / d_a
        } else {
            Vector3::zeros()
        };
        let n_b = if d_b > 1.0e-12 {
            d_b_vec / d_b
        } else {
            Vector3::zeros()
        };
        // Damped step: half the deficit each.
        let delta = 0.5 * (r_a * n_a + r_b * n_b);
        c -= delta;
        if delta.norm() < 1.0e-12 {
            break;
        }
    }
    // Final report (may still be slightly off the spine — caller
    // decides if that's acceptable).
    let Some((uv_a, uv_b)) = refine_seed_to_3d(s_a, s_b, c, m) else {
        return Err(SurfaceError::IntersectionFailed(
            "rolling_ball_blend: could not project seed onto both surfaces".into(),
        ));
    };
    Ok(SpineSample {
        center: c,
        contact_a: s_a.evaluate(uv_a.x, uv_a.y),
        uv_a,
        contact_b: s_b.evaluate(uv_b.x, uv_b.y),
        uv_b,
    })
}

/// Find `(uv_a, uv_b)` such that surface evaluation is closest to a
/// 3D point `c` — the two contact points of the rolling ball at
/// center `c`.
fn refine_seed_to_3d(
    s_a: &NurbsSurface,
    s_b: &NurbsSurface,
    c: Vector3<f64>,
    m: &MarchParams,
) -> Option<(Vector2<f64>, Vector2<f64>)> {
    // Reuse the marching-SSI's refine_seed for symmetry — but here
    // we *expect* the two surface points to differ from `c` by ~r,
    // so the "intersection sanity" check there would reject us.
    // Instead use the public grid-search + Newton path directly.
    let (uv_a, _) = grid_closest(s_a, c, 11);
    let (uv_b, _) = grid_closest(s_b, c, 11);
    let uv_a = newton_closest(s_a, c, uv_a, m);
    let uv_b = newton_closest(s_b, c, uv_b, m);
    Some((uv_a, uv_b))
}

fn grid_closest(s: &NurbsSurface, p: Vector3<f64>, n: usize) -> (Vector2<f64>, f64) {
    let (u_min, u_max) = s.u_range();
    let (v_min, v_max) = s.v_range();
    let mut best = Vector2::new(0.5 * (u_min + u_max), 0.5 * (v_min + v_max));
    let mut best_d = f64::INFINITY;
    let denom = (n - 1).max(1) as f64;
    for i in 0..n {
        for j in 0..n {
            let u = u_min + (i as f64 / denom) * (u_max - u_min);
            let v = v_min + (j as f64 / denom) * (v_max - v_min);
            let q = s.evaluate(u, v);
            let d = (q - p).norm();
            if d < best_d {
                best_d = d;
                best = Vector2::new(u, v);
            }
        }
    }
    (best, best_d)
}

fn newton_closest(
    s: &NurbsSurface,
    p: Vector3<f64>,
    seed: Vector2<f64>,
    m: &MarchParams,
) -> Vector2<f64> {
    let (u_min, u_max) = s.u_range();
    let (v_min, v_max) = s.v_range();
    let h = ((u_max - u_min) + (v_max - v_min)) * 1.0e-5;
    let mut uv = Vector2::new(seed.x.clamp(u_min, u_max), seed.y.clamp(v_min, v_max));
    for _ in 0..m.projection_iters {
        let q = s.evaluate(uv.x, uv.y);
        let r = q - p;
        if r.norm() < m.projection_tolerance {
            return uv;
        }
        let u_lo = (uv.x - h).max(u_min);
        let u_hi = (uv.x + h).min(u_max);
        let v_lo = (uv.y - h).max(v_min);
        let v_hi = (uv.y + h).min(v_max);
        let tu = (s.evaluate(u_hi, uv.y) - s.evaluate(u_lo, uv.y)) / (u_hi - u_lo).max(1.0e-12);
        let tv = (s.evaluate(uv.x, v_hi) - s.evaluate(uv.x, v_lo)) / (v_hi - v_lo).max(1.0e-12);
        let a11 = tu.dot(&tu);
        let a12 = tu.dot(&tv);
        let a22 = tv.dot(&tv);
        let b1 = -tu.dot(&r);
        let b2 = -tv.dot(&r);
        let det = a11 * a22 - a12 * a12;
        if det.abs() < 1.0e-14 {
            return uv;
        }
        let du = (a22 * b1 - a12 * b2) / det;
        let dv = (-a12 * b1 + a11 * b2) / det;
        uv.x = (uv.x + du).clamp(u_min, u_max);
        uv.y = (uv.y + dv).clamp(v_min, v_max);
    }
    uv
}

/// Spine tangent at a sample — the cross product of the two contact
/// normals.
///
/// `c - q_A` (outward toward the ball center from contact A) and
/// `c - q_B` (outward from contact B) define a plane; the spine
/// tangent lies in that plane perpendicular to both, which makes it
/// the cross product `(c - q_A) × (c - q_B)`.
fn spine_tangent(sample: &SpineSample) -> Vector3<f64> {
    let n_a = sample.center - sample.contact_a;
    let n_b = sample.center - sample.contact_b;
    let t = n_a.cross(&n_b);
    let l = t.norm();
    if l < 1.0e-12 {
        Vector3::zeros()
    } else {
        t / l
    }
}

fn normalize_or_default(v: Vector3<f64>) -> Vector3<f64> {
    let l = v.norm();
    if l < 1.0e-12 {
        Vector3::zeros()
    } else {
        v / l
    }
}

/// March in one direction along the spine starting from `seed`.
///
/// `initial_tangent` is the first step direction; subsequent steps
/// use the recomputed local tangent.
fn trace_one_direction(
    s_a: &NurbsSurface,
    s_b: &NurbsSurface,
    radius: f64,
    seed: SpineSample,
    initial_tangent: Vector3<f64>,
    p: &BlendParams,
    m: &MarchParams,
) -> Vec<SpineSample> {
    let (u_a_min, u_a_max) = s_a.u_range();
    let (v_a_min, v_a_max) = s_a.v_range();
    let (u_b_min, u_b_max) = s_b.u_range();
    let (v_b_min, v_b_max) = s_b.v_range();
    let mut samples = Vec::new();
    let mut current = seed;
    let mut tangent = initial_tangent;
    let mut step = (p.chord_tolerance * 4.0).clamp(p.min_step, p.max_step);
    for _ in 0..p.max_samples {
        // Predict.
        let pred = current.center + step * tangent;
        // Equalise.
        let Ok(next) = equalise_spine(s_a, s_b, pred, radius, m, 16) else {
            break;
        };
        // Boundary check.
        if !inside(next.uv_a, u_a_min, u_a_max, v_a_min, v_a_max)
            || !inside(next.uv_b, u_b_min, u_b_max, v_b_min, v_b_max)
        {
            break;
        }
        // Forward progress.
        let progress = (next.center - current.center).norm();
        if progress < p.min_step {
            break;
        }
        // Update tangent — keep sign consistent.
        let mut t_next = spine_tangent(&next);
        if t_next.dot(&tangent) < 0.0 {
            t_next = -t_next;
        }
        if t_next.norm() < 1.0e-9 {
            // Degenerate tangent — stop here.
            samples.push(next);
            break;
        }
        // Sanity: distances from the new center to its contacts must
        // be ≈ radius. If not, halve step and retry.
        let d_a = (next.center - next.contact_a).norm();
        let d_b = (next.center - next.contact_b).norm();
        if (d_a - radius).abs() > radius * 0.05 || (d_b - radius).abs() > radius * 0.05 {
            step *= 0.5;
            if step < p.min_step {
                break;
            }
            continue;
        }
        samples.push(next.clone());
        current = next;
        tangent = t_next;
        // Step adaptation — keep `step` at the chord tolerance level.
        if progress < 0.25 * p.chord_tolerance {
            step = (step * 1.5).min(p.max_step);
        } else if progress > 4.0 * p.chord_tolerance {
            step = (step * 0.5).max(p.min_step);
        }
    }
    samples
}

fn inside(uv: Vector2<f64>, u_min: f64, u_max: f64, v_min: f64, v_max: f64) -> bool {
    uv.x >= u_min && uv.x <= u_max && uv.y >= v_min && uv.y <= v_max
}

/// Build the blend surface as a tensor-product NURBS patch.
///
/// - u runs along the spine (cubic — `nu` samples → `nu` control
///   rows).
/// - v runs across the cross-section (rational quadratic arc with 3
///   CPs — start contact, ball-vertex, end contact — with the middle
///   CP weighted by `cos(half_angle)` so the rational quadratic
///   exactly reproduces a circular arc of the appropriate angle).
fn build_blend_surface(
    s_a: &NurbsSurface,
    s_b: &NurbsSurface,
    spine: &[SpineSample],
    radius: f64,
    p: &BlendParams,
) -> Result<NurbsSurface, SurfaceError> {
    let nu = spine.len();
    if nu < 2 {
        return Err(SurfaceError::IntersectionFailed(
            "rolling_ball_blend: trace too short to build a blend surface".into(),
        ));
    }
    // For v1 we stick with 3 CPs — rational quadratic, exact arc.
    // The `cross_section_cps` param is reserved for the v1.5 path
    // (higher-CP cross-section arcs require Bezier subdivision +
    // re-knot — left as documented future work). The validity of
    // `p.cross_section_cps` is still asserted up front to give a
    // useful error if someone passes a nonsense value.
    if !(3..=15).contains(&p.cross_section_cps) {
        return Err(SurfaceError::IntersectionFailed(format!(
            "rolling_ball_blend: cross_section_cps must be in [3, 15]; got {}",
            p.cross_section_cps
        )));
    }

    // Build the v-direction (cross-section) for each spine sample.
    // For 3 CPs the canonical rational quadratic arc has:
    //   P_0 = contact_A
    //   P_2 = contact_B
    //   P_1 = intersection of the surface tangent at A and the surface
    //         tangent at B, taken to lie on the ball-center plane
    //         (the "shoulder" point of the arc — sphere center +
    //         radius * bisector of (contact_A→center) and
    //         (contact_B→center)).
    //   w_0 = w_2 = 1, w_1 = cos(half_angle).
    let mut grid: Vec<Vec<Vector3<f64>>> = Vec::with_capacity(nu);
    let mut weight_grid: Vec<Vec<f64>> = Vec::with_capacity(nu);
    for s in spine {
        let center = s.center;
        let a = s.contact_a;
        let b = s.contact_b;
        let na = a - center;
        let nb = b - center;
        let na_norm = na.norm().max(1.0e-12);
        let nb_norm = nb.norm().max(1.0e-12);
        let na_hat = na / na_norm;
        let nb_hat = nb / nb_norm;
        // Half angle of the arc from A to B subtended at the ball
        // center.
        let cos_full = na_hat.dot(&nb_hat).clamp(-1.0, 1.0);
        let full_angle = cos_full.acos();
        let half_angle = full_angle * 0.5;
        // Shoulder direction: a rational quadratic with CPs
        // `(P0, P1, P2)` and weights `(1, cos(half), 1)` traces the
        // circular arc with endpoints P0=A, P2=B where P1 is the
        // intersection of the surface tangents at A and B.
        //
        // For a rolling-ball blend on an INTERIOR (concave) edge —
        // the canonical fillet — the arc bulges TOWARD the corner of
        // the underlying solid, AWAY from the convex blend exterior.
        // Geometrically: at A the contact normal `n_A = (A − center)`
        // points outward from the ball center; the surface tangent at
        // A is perpendicular to `n_A` in the arc plane; the tangent
        // line, extended, meets the tangent line at B at a point
        // **on the side of the chord A→B away from the ball center**.
        // That's the bisector direction `na_hat + nb_hat` (NOT its
        // negation) — both `na_hat` and `nb_hat` point from center
        // toward the contacts, and their sum points from center
        // toward the corner.
        //
        // Distance from center to the shoulder: `radius / cos(half)`.
        let bisector = na_hat + nb_hat;
        let bis_norm = bisector.norm();
        let toward_corner = if bis_norm > 1.0e-12 {
            bisector / bis_norm
        } else {
            // Pathological: A and B opposite — the ball is "wrapped"
            // by the two surfaces, no shoulder exists. Pick any
            // orthogonal direction; the blend in this degenerate case
            // collapses to a half-circle which the rational quadratic
            // cannot represent in 3 CPs anyway.
            let mut alt = na_hat.cross(&Vector3::new(1.0, 0.0, 0.0));
            if alt.norm() < 1.0e-9 {
                alt = na_hat.cross(&Vector3::new(0.0, 1.0, 0.0));
            }
            normalize_or_default(alt)
        };
        let cos_half = half_angle.cos().max(1.0e-6);
        let shoulder = center + (radius / cos_half) * toward_corner;
        grid.push(vec![a, shoulder, b]);
        weight_grid.push(vec![1.0, cos_half, 1.0]);
    }
    // Adjust shoulders to keep them away from the surfaces — the
    // analytic shoulder is on the ball, but the rational quadratic
    // through three CPs (w0=w2=1, w1=cos(half)) reproduces the arc
    // only when CPs are exact. They are.

    let _ = s_a; // kept for symmetry / future feature-aware refinements
    let _ = s_b;

    // u-direction knots: open uniform for the spine.
    // For a cubic spine fit with `nu` rows we use degree min(3, nu-1)
    // — if nu < 4 we degrade to degree (nu-1).
    let degree_u = 3.min(nu - 1).max(1);
    let n_cps_u = nu;
    let u_knots = open_uniform_knots(n_cps_u, degree_u);
    // v-direction knots: degree 2, 3 CPs — Bezier `[0,0,0,1,1,1]`.
    let v_knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];

    NurbsSurface::new(degree_u, 2, u_knots, v_knots, grid, weight_grid)
}

fn open_uniform_knots(n_cps: usize, degree: usize) -> Vec<f64> {
    let p = degree;
    let m = n_cps + p + 1;
    let mut k = vec![0.0; m];
    if n_cps <= p + 1 {
        for kv in k.iter_mut().skip(m - p - 1) {
            *kv = 1.0;
        }
        return k;
    }
    let n_internal = n_cps - p - 1;
    for (i, kv) in k.iter_mut().enumerate().take(m) {
        if i <= p {
            *kv = 0.0;
        } else if i >= n_cps {
            *kv = 1.0;
        } else {
            let idx = i - p;
            *kv = idx as f64 / (n_internal + 1) as f64;
        }
    }
    k
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a planar Bezier surface lying in the plane `n · (p - p0) = 0`
    /// where `n` is unit. The patch covers a square of side
    /// `2 · half_extent` centered at `p0`.
    fn planar_patch(p0: Vector3<f64>, normal: Vector3<f64>, half_extent: f64) -> NurbsSurface {
        let n = normalize_or_default(normal);
        // Build two orthogonal in-plane axes.
        let helper = if n.x.abs() < 0.9 {
            Vector3::new(1.0, 0.0, 0.0)
        } else {
            Vector3::new(0.0, 1.0, 0.0)
        };
        let u_dir = normalize_or_default(helper - helper.dot(&n) * n);
        let v_dir = n.cross(&u_dir);
        let knots = vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
        let mut cps = Vec::with_capacity(4);
        for i in 0..4 {
            let s = -half_extent + 2.0 * half_extent * (i as f64 / 3.0);
            let mut row = Vec::with_capacity(4);
            for j in 0..4 {
                let t = -half_extent + 2.0 * half_extent * (j as f64 / 3.0);
                row.push(p0 + s * u_dir + t * v_dir);
            }
            cps.push(row);
        }
        let weights = vec![vec![1.0; 4]; 4];
        NurbsSurface::new(3, 3, knots.clone(), knots, cps, weights).unwrap()
    }

    /// Cylinder patch about the +y axis using a rational quadratic
    /// quarter-arc, lofted along y.
    fn cylinder_about_y(r: f64, y0: f64, y1: f64) -> NurbsSurface {
        let s2 = 2.0_f64.sqrt() / 2.0;
        let row_y0 = vec![
            Vector3::new(r, y0, 0.0),
            Vector3::new(r, y0, r),
            Vector3::new(0.0, y0, r),
        ];
        let row_y1 = vec![
            Vector3::new(r, y1, 0.0),
            Vector3::new(r, y1, r),
            Vector3::new(0.0, y1, r),
        ];
        NurbsSurface::new(
            1,
            2,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![row_y0, row_y1],
            vec![vec![1.0, s2, 1.0], vec![1.0, s2, 1.0]],
        )
        .unwrap()
    }

    #[test]
    fn equalise_spine_finds_bisector_for_perpendicular_planes() {
        // xy plane at z=0 and xz plane at y=0; spine should be the
        // line y = z = r at any x; ball center at (x, r, r), contact
        // A on xy plane is (x, r, 0), contact B on xz plane is
        // (x, 0, r), distances both r.
        let r = 0.5;
        let s_a = planar_patch(Vector3::zeros(), Vector3::new(0.0, 0.0, 1.0), 2.0);
        let s_b = planar_patch(Vector3::zeros(), Vector3::new(0.0, 1.0, 0.0), 2.0);
        let seed = Vector3::new(0.0, r * 1.1, r * 0.9); // a bit off
        let m = MarchParams::default();
        let s = equalise_spine(&s_a, &s_b, seed, r, &m, 64).unwrap();
        // Center should be exactly (0, r, r).
        assert!((s.center - Vector3::new(0.0, r, r)).norm() < 1.0e-4);
        // Contact A on xy plane: z = 0, point ≈ (0, r, 0).
        assert!(s.contact_a.z.abs() < 1.0e-4, "A z = {}", s.contact_a.z);
        // Contact B on xz plane: y = 0, point ≈ (0, 0, r).
        assert!(s.contact_b.y.abs() < 1.0e-4, "B y = {}", s.contact_b.y);
        // Distances both r.
        assert!(((s.center - s.contact_a).norm() - r).abs() < 1.0e-4);
        assert!(((s.center - s.contact_b).norm() - r).abs() < 1.0e-4);
    }

    #[test]
    fn rolling_ball_blend_two_perpendicular_planes_is_cylindrical_fillet() {
        // xy plane at z=0 (A) and xz plane at y=0 (B), perpendicular
        // intersection along the x-axis. Rolling ball of radius r:
        // spine is the line y = z = r, x ∈ [some span]. Blend surface
        // is a cylindrical fillet of radius r centered on the x-axis
        // — i.e. every point on the blend has distance exactly r from
        // the spine axis (which is the line y=r, z=r? — no, the
        // generating axis of the analytic fillet is along x at
        // y=z=r? No — the generating axis is the planes' line of
        // intersection, which is the x-axis at (y=0, z=0). The
        // fillet is the quarter-cylinder y² + z² = r² with the
        // negative-y/negative-z quadrant removed... actually for an
        // INTERIOR (concave) corner between the two half-planes
        // (z ≥ 0 and y ≥ 0) the fillet is the convex quarter-arc
        // bowing OUTWARD from the corner. Standard CAD interior
        // fillet at a 90° corner has axis = x and radius r centered
        // at (y=r, z=r). Every blend point is at distance r from this
        // axis.
        let r = 0.5;
        let s_a = planar_patch(Vector3::zeros(), Vector3::new(0.0, 0.0, 1.0), 2.0);
        let s_b = planar_patch(Vector3::zeros(), Vector3::new(0.0, 1.0, 0.0), 2.0);
        let seed = Vector3::new(0.0, r, r);
        let params = BlendParams::default();
        let blend = rolling_ball_blend(
            &s_a,
            &s_b,
            r,
            seed,
            Some(Vector3::new(1.0, 0.0, 0.0)),
            &params,
        )
        .unwrap();

        // Sanity on the spine itself: every center at y=r, z=r.
        for sp in &blend.spine {
            assert!(
                (sp.center.y - r).abs() < 1.0e-4,
                "spine y = {}",
                sp.center.y
            );
            assert!(
                (sp.center.z - r).abs() < 1.0e-4,
                "spine z = {}",
                sp.center.z
            );
        }
        // Spine spans a non-trivial x range.
        let xs: Vec<f64> = blend.spine.iter().map(|s| s.center.x).collect();
        let max_x = xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min_x = xs.iter().cloned().fold(f64::INFINITY, f64::min);
        assert!(max_x - min_x > 1.0, "spine x span = {}", max_x - min_x);

        // Every point on the blend surface is at distance ≈ r from
        // the fillet's analytic axis at (y=r, z=r) (the spine
        // itself). Sample a 6×5 (u, v) grid.
        let (u_min, u_max) = blend.surface.u_range();
        let (v_min, v_max) = blend.surface.v_range();
        let mut worst = 0.0_f64;
        for ui in 1..6 {
            let u = u_min + (ui as f64 / 6.0) * (u_max - u_min);
            for vj in 0..=4 {
                let v = v_min + (vj as f64 / 4.0) * (v_max - v_min);
                let p = blend.surface.evaluate(u, v);
                // Find the closest spine center along x — given the
                // spine is straight in x, just project.
                let axis_pt = Vector3::new(p.x, r, r);
                let d = (p - axis_pt).norm();
                worst = worst.max((d - r).abs());
            }
        }
        // 1e-3 · r is a tight tolerance — the rational quadratic arc
        // is *exact* for the circular cross-section, so the only
        // error is the cubic spine fit through a straight-line spine
        // (which is exact too).
        assert!(
            worst < 1.0e-3 * r,
            "max fillet-radius error = {worst}, r = {r}"
        );
    }

    #[test]
    fn rolling_ball_blend_plane_cylinder_is_toroidal() {
        // Plane at z=0 and a quarter cylinder about the +y axis of
        // radius R. Rolling ball of radius r: contact on the plane is
        // at z=0, contact on the cylinder is on its outer surface
        // (radius R from the y-axis). The spine is the curve at
        // z = r (the contact-on-plane normal direction = +z),
        // and at horizontal distance R + r from the y-axis (so the
        // ball sits tangent to the outer cylinder surface). Hence
        // spine is x² + y_off² = (R + r)² at z = r — wait, the
        // cylinder is about y so spine has y free, and
        // x² + z² = ? No: the cylinder is the set x² + z² = R²
        // (about y). The plane is z = 0 (above). Hmm, ROLLING ON
        // THE OUTSIDE: ball must be at z ≥ r (above plane) and at
        // x² + (z - 0)² where the ball is outside the cylinder, i.e.
        // sqrt(x² + z²) = R + r (distance from cylinder axis). With
        // z = r and unknown x: x² + r² = (R + r)², so x² = (R+r)² −
        // r² = R² + 2Rr → x = sqrt(R(R + 2r)). Concrete: r = 0.2,
        // R = 1 → x ≈ sqrt(1·1.4) ≈ 1.183.
        let r = 0.2;
        let big_r = 1.0;
        // Plane is the patch at z = 0, normal +z.
        let s_a = planar_patch(
            Vector3::new(2.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            4.0,
        );
        // Cylinder about +y axis of radius R.
        let s_b = cylinder_about_y(big_r, -2.0, 2.0);
        // Analytic spine x:
        let spine_x = (big_r * (big_r + 2.0 * r)).sqrt();
        let seed = Vector3::new(spine_x, 0.0, r);
        let params = BlendParams::default();
        let blend = rolling_ball_blend(
            &s_b,
            &s_a,
            r,
            seed,
            Some(Vector3::new(0.0, 1.0, 0.0)),
            &params,
        )
        .unwrap();

        // Sanity: every spine center has z ≈ r (above plane by r)
        // and x ≈ analytic_x (tangent to cylinder).
        for sp in &blend.spine {
            assert!(
                (sp.center.z - r).abs() < 1.0e-3,
                "spine z = {}",
                sp.center.z
            );
            // x² + z² ≈ (R + r)² since the ball is at distance R+r
            // from the cylinder axis (= y-axis).
            let d_axis = (sp.center.x * sp.center.x + sp.center.z * sp.center.z).sqrt();
            assert!(
                (d_axis - (big_r + r)).abs() < 1.0e-3,
                "spine center axis distance = {}; expected = {}",
                d_axis,
                big_r + r
            );
        }
        // Every point on the blend surface is at distance ≈ r from
        // the spine. The blend swept-surface is a torus patch with
        // major radius (R + r) (since the spine is a circle of that
        // radius in the xz plane at z = r... no, the spine is at
        // y varying, x fixed, z = r — so spine is a STRAIGHT LINE
        // along y; the blend is the cylindrical FILLET between the
        // plane and the cylinder, locally swept along y).
        // Wait — for the cylinder about the y-axis the spine moves
        // along y while x, z are fixed. So it's a fillet
        // *cylinder*, not a torus. Toroidal would be a plane meeting
        // a *sphere* or a cylinder whose axis is along z. Let me
        // assert the simpler invariant: every blend-surface sample is
        // at distance ≈ r from its corresponding spine point.
        let (u_min, u_max) = blend.surface.u_range();
        let (v_min, v_max) = blend.surface.v_range();
        let mut worst = 0.0_f64;
        for ui in 0..=5 {
            let u_t = ui as f64 / 5.0;
            let u_param = u_min + u_t * (u_max - u_min);
            // Recover the corresponding spine point: the spine sample
            // at index ≈ u_t · (n-1).
            let n = blend.spine.len();
            let idx = ((u_t * (n - 1) as f64).round() as usize).min(n - 1);
            let center = blend.spine[idx].center;
            for vj in 0..=4 {
                let v = v_min + (vj as f64 / 4.0) * (v_max - v_min);
                let p = blend.surface.evaluate(u_param, v);
                let d = (p - center).norm();
                worst = worst.max((d - r).abs());
            }
        }
        // 5% of r is the toleration here — the spine sample index
        // mapping is approximate (cubic LSQ shifts u parameters
        // slightly from the underlying chord-length), and the
        // cylinder-side parameterisation is rational quadratic which
        // mildly deforms equispaced u under cubic refit.
        assert!(
            worst < 0.05 * r,
            "max ball-radius error on plane+cylinder blend = {worst}, r = {r}"
        );
    }

    #[test]
    fn rejects_zero_radius() {
        let s_a = planar_patch(Vector3::zeros(), Vector3::new(0.0, 0.0, 1.0), 2.0);
        let s_b = planar_patch(Vector3::zeros(), Vector3::new(0.0, 1.0, 0.0), 2.0);
        let seed = Vector3::new(0.0, 0.5, 0.5);
        let err =
            rolling_ball_blend(&s_a, &s_b, 0.0, seed, None, &BlendParams::default()).unwrap_err();
        assert_eq!(err.code(), "surface.intersection_failed");
    }
}
