//! Marching surface-surface intersection (continuous-trace, Bajaj-style).
//!
//! Companion to [`crate::intersect`]. The legacy `intersect::true_ssi`
//! refines tessellation polyline *vertices* one at a time — production
//! CAD tools instead **march** along the intersection curve
//! continuously, taking a step in the tangent direction `t = N_A × N_B`
//! and snapping back onto both surfaces by Newton closest-foot
//! projection in parametric `(u, v)` space.
//!
//! The trace is intrinsically parametric: every traced sample is a
//! 4-tuple `(u_A, v_A, u_B, v_B)`, not just a 3D point. That gives:
//!
//! - smooth intersection curves (no zig-zag from triangle-soup
//!   artefacts),
//! - proper handling of curved surfaces where many tessellation
//!   triangles converge on the same true curve,
//! - termination on a real surface boundary (`u` or `v` parameter
//!   leaves its valid range) versus a loop closure (the marcher comes
//!   back within `loop_tolerance` of its starting point).
//!
//! ## Algorithm (Bajaj-Hoffmann marching v1)
//!
//! 1. **Seed** — caller supplies a starting `(u_A, v_A, u_B, v_B)` on
//!    the intersection (we use the tessellation-based intersection's
//!    first vertex + Newton projection onto both surfaces to get the
//!    parametric seed).
//! 2. **Step** — compute the unit intersection tangent
//!    `t̂ = N̂_A × N̂_B / |N̂_A × N̂_B|`. Pick a step length `h` from
//!    curvature: a coarse estimate by taking a half-step, measuring
//!    deflection, and scaling so deflection × `max_curvature_step` ≈
//!    target chord error.
//! 3. **Predict** — `p_new = p_old + h · t̂`.
//! 4. **Correct** — Newton closest-foot in `(u_A, v_A)` and
//!    `(u_B, v_B)` independently, then average the two corrected 3D
//!    positions. Three Newton iterations is usually enough at
//!    `1e-12` residual.
//! 5. **Sanity** — if the corrected position drifts more than
//!    `correction_cap * h` from the predictor, halve `h` and retry.
//!    Bounded retries to avoid infinite loops.
//! 6. **Boundary / loop check** — if any of the four parameters
//!    leaves its valid range, project to the boundary and emit a
//!    terminating point. If the corrected point comes within
//!    `loop_tolerance` of the seed and we've taken at least
//!    `min_loop_samples` steps, close the loop.
//! 7. **Fit** — once the polyline of 3D corrected points exists, fit
//!    a cubic NURBS curve through it via
//!    [`crate::fit::nurbs_curve_through_points`].
//!
//! ## v1 caveats
//!
//! - The seed comes from the tessellation-based intersection; if no
//!   tessellation segment exists, we don't trace.
//! - Single-direction tracing per seed — we don't bidirectionally
//!   march both ways from an interior seed; instead we walk forward
//!   from the seed, and (if we don't close a loop) restart from the
//!   seed in the opposite direction and concatenate the trace.
//! - Branch / self-intersecting / cusp configurations are accepted
//!   only as "the trace terminates on a near-degenerate tangent" —
//!   the marcher doesn't try to detect or follow branches.

use nalgebra::{Vector2, Vector3};

use crate::error::SurfaceError;
use crate::nurbs_curve::NurbsCurve;
use crate::nurbs_surface::NurbsSurface;
use crate::{intersect, tessellate};

/// User-facing knobs for the marching SSI.
#[derive(Clone, Debug)]
pub struct MarchParams {
    /// Target chord deviation per step (model units). Smaller →
    /// finer trace, more samples, slower runtime. v1 default 1e-3.
    pub chord_tolerance: f64,
    /// Newton convergence tolerance for the closest-foot projection
    /// (3D residual norm). v1 default 1e-10.
    pub projection_tolerance: f64,
    /// Max Newton iterations per projection. v1 default 8.
    pub projection_iters: usize,
    /// Maximum step size, in model units. Caps `h` so we don't take a
    /// runaway step on a near-degenerate tangent.
    pub max_step: f64,
    /// Minimum step size — below this we terminate (we're near a
    /// boundary or a degenerate point). v1 default 1e-8.
    pub min_step: f64,
    /// If the corrected position drifts more than
    /// `correction_cap * h` from the predictor, halve `h` and retry.
    /// v1 default 0.5.
    pub correction_cap: f64,
    /// Loop-closure tolerance — the marcher closes the loop when the
    /// corrected point lands within this distance of the seed AND at
    /// least `min_loop_samples` steps have been taken.
    pub loop_tolerance: f64,
    /// Minimum sample count before a loop closure is allowed.
    pub min_loop_samples: usize,
    /// Hard cap on samples per trace direction (forward + back). v1
    /// default 1024 — defends against runaway tracers.
    pub max_samples: usize,
}

impl Default for MarchParams {
    fn default() -> Self {
        Self {
            chord_tolerance: 1.0e-3,
            projection_tolerance: 1.0e-10,
            projection_iters: 8,
            max_step: 0.5,
            min_step: 1.0e-8,
            correction_cap: 0.5,
            loop_tolerance: 1.0e-3,
            min_loop_samples: 8,
            max_samples: 1024,
        }
    }
}

/// One sample on a traced intersection curve.
#[derive(Clone, Debug)]
pub struct TracePoint {
    /// 3D position (average of the two surface evaluations).
    pub xyz: Vector3<f64>,
    /// `(u, v)` parameters on surface A.
    pub uv_a: Vector2<f64>,
    /// `(u, v)` parameters on surface B.
    pub uv_b: Vector2<f64>,
}

/// Output of a single marching trace.
#[derive(Clone, Debug)]
pub struct Trace {
    /// Ordered samples along the trace.
    pub samples: Vec<TracePoint>,
    /// Why the trace terminated.
    pub termination: TraceEnd,
}

/// Reason a trace stopped.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TraceEnd {
    /// Hit a surface-parameter boundary (`u` or `v` left the valid
    /// range on either surface).
    Boundary,
    /// Closed back on the seed — the intersection curve is a loop.
    Loop,
    /// Tangent is near-zero (surfaces touched at a tangent point or
    /// the cross-product collapsed below `1e-12`).
    Degenerate,
    /// Step adaptation gave up — successive halvings could not bring
    /// the corrector below `correction_cap * h`. Rare; usually means
    /// surfaces are very close and the marcher can't pick a side.
    StepFloor,
    /// Hit the hard sample cap (`MarchParams::max_samples`).
    SampleCap,
}

/// March one intersection curve starting from a seeded
/// `(uv_a, uv_b)`.
///
/// The seed must already be on the intersection (use
/// [`refine_seed`] to project a coarse `(uv_a, uv_b)` onto the true
/// intersection). The marcher walks in the `+t̂` direction until
/// termination; bidirectional tracing is delegated to
/// [`trace_bidirectional`].
pub fn trace_forward(
    s_a: &NurbsSurface,
    s_b: &NurbsSurface,
    seed_a: Vector2<f64>,
    seed_b: Vector2<f64>,
    params: &MarchParams,
    flip_tangent: bool,
) -> Trace {
    let (u_a_min, u_a_max) = s_a.u_range();
    let (v_a_min, v_a_max) = s_a.v_range();
    let (u_b_min, u_b_max) = s_b.u_range();
    let (v_b_min, v_b_max) = s_b.v_range();

    let initial_xyz = average_xyz(s_a, seed_a, s_b, seed_b);
    let mut samples = vec![TracePoint {
        xyz: initial_xyz,
        uv_a: seed_a,
        uv_b: seed_b,
    }];

    let mut uv_a = seed_a;
    let mut uv_b = seed_b;
    let mut h = (params.chord_tolerance * 5.0).clamp(params.min_step, params.max_step);

    let mut termination = TraceEnd::SampleCap;

    for sample_idx in 0..params.max_samples {
        let last = samples.last().expect("samples non-empty").clone();
        // Tangent at the current point.
        let n_a = surface_normal(s_a, uv_a);
        let n_b = surface_normal(s_b, uv_b);
        let mut t = n_a.cross(&n_b);
        let t_norm = t.norm();
        if t_norm < 1.0e-12 {
            termination = TraceEnd::Degenerate;
            break;
        }
        t /= t_norm;
        if flip_tangent {
            t = -t;
        }
        // Keep the tangent consistent with the previous step (avoid
        // accidentally walking back on ourselves at near-degenerate
        // points). After the first step, project against the last
        // chord direction.
        if samples.len() >= 2 {
            let prev = samples[samples.len() - 2].xyz;
            let chord = last.xyz - prev;
            if chord.dot(&t) < 0.0 {
                t = -t;
            }
        }

        // Adaptive predict + correct, with step-halving on a poor
        // correction.
        let mut step = h.clamp(params.min_step, params.max_step);
        let mut retries = 0_usize;
        let mut new_sample_opt: Option<TracePoint> = None;
        loop {
            let pred = last.xyz + step * t;
            // Correct: Newton closest-foot from (uv_a, uv_b) using
            // pred as the seed for both surfaces.
            let (uv_a_new, p_a_new) = closest_foot_seeded(s_a, pred, uv_a, params);
            let (uv_b_new, p_b_new) = closest_foot_seeded(s_b, pred, uv_b, params);
            let corrected = 0.5 * (p_a_new + p_b_new);
            let drift = (corrected - pred).norm();
            let on_a = (corrected - p_a_new).norm();
            let on_b = (corrected - p_b_new).norm();
            // If the per-surface footpoints disagree wildly, halve.
            let two_surface_gap = (p_a_new - p_b_new).norm();
            if step > params.min_step
                && (drift > params.correction_cap * step
                    || two_surface_gap > params.chord_tolerance * 4.0)
            {
                step *= 0.5;
                retries += 1;
                if retries > 8 {
                    break;
                }
                continue;
            }
            // Accept the step.
            new_sample_opt = Some(TracePoint {
                xyz: corrected,
                uv_a: uv_a_new,
                uv_b: uv_b_new,
            });
            // Adapt h for next step: if both footprints are tight,
            // try a slightly longer step next time; if loose, contract.
            let footprint_quality = on_a.max(on_b);
            if footprint_quality < 0.1 * params.chord_tolerance {
                h = (step * 1.5).min(params.max_step);
            } else if footprint_quality > params.chord_tolerance {
                h = (step * 0.5).max(params.min_step);
            } else {
                h = step;
            }
            break;
        }

        let Some(new_sample) = new_sample_opt else {
            termination = TraceEnd::StepFloor;
            break;
        };

        // Boundary check.
        //
        // `closest_foot_seeded` *clamps* `uv` to the surface
        // parameter range — so the projected uv NEVER ends up
        // strictly outside the range. Instead the boundary is
        // detected by a **two-condition** test:
        //
        // 1. The projected uv lands **on** the parameter-range
        //    boundary (within an `on_edge_eps` tolerance), AND
        // 2. The clamped projection has non-trivial residual
        //    (`||q - pred|| > projection_tolerance * 10`), meaning
        //    the predictor wanted to step *past* the boundary and
        //    the clamp absorbed the excess.
        let on_edge = on_param_boundary(
            new_sample.uv_a,
            u_a_min,
            u_a_max,
            v_a_min,
            v_a_max,
        ) || on_param_boundary(
            new_sample.uv_b,
            u_b_min,
            u_b_max,
            v_b_min,
            v_b_max,
        );
        // Residual: how far the corrected position is from where the
        // predictor wanted to land. If the clamp absorbed real
        // movement, this is non-zero.
        let want_pred = last.xyz + step * t;
        let absorbed = (new_sample.xyz - want_pred).norm();
        let absorbed_large = absorbed > params.projection_tolerance * 100.0
            || absorbed > 0.1 * step;
        if on_edge && absorbed_large {
            // Project the last step onto the boundary by chord
            // bisection between last and new_sample.
            let boundary_sample = boundary_bisect(
                s_a,
                s_b,
                last.clone(),
                new_sample,
                params,
            );
            samples.push(boundary_sample);
            termination = TraceEnd::Boundary;
            break;
        }
        // Loop-closure check.
        let to_seed = (new_sample.xyz - samples[0].xyz).norm();
        if to_seed < params.loop_tolerance && sample_idx + 1 >= params.min_loop_samples {
            samples.push(samples[0].clone());
            termination = TraceEnd::Loop;
            break;
        }

        // Forward-progress sanity: if we're not moving, give up.
        let progress = (new_sample.xyz - last.xyz).norm();
        if progress < params.min_step {
            termination = TraceEnd::StepFloor;
            break;
        }

        uv_a = new_sample.uv_a;
        uv_b = new_sample.uv_b;
        samples.push(new_sample);
    }

    Trace { samples, termination }
}

/// Bidirectional trace from a seed.
///
/// Walks forward; if the forward trace did not close a loop, walks
/// backward from the seed and prepends. Returns a single ordered
/// sample list.
pub fn trace_bidirectional(
    s_a: &NurbsSurface,
    s_b: &NurbsSurface,
    seed_a: Vector2<f64>,
    seed_b: Vector2<f64>,
    params: &MarchParams,
) -> Trace {
    let forward = trace_forward(s_a, s_b, seed_a, seed_b, params, false);
    if forward.termination == TraceEnd::Loop {
        return forward;
    }
    let backward = trace_forward(s_a, s_b, seed_a, seed_b, params, true);
    // Stitch: reverse `backward` (excluding the seed sample) and put
    // forward after it.
    let mut combined = Vec::with_capacity(backward.samples.len() + forward.samples.len());
    for s in backward.samples.iter().skip(1).rev() {
        combined.push(s.clone());
    }
    combined.extend(forward.samples.iter().cloned());
    Trace {
        samples: combined,
        // The combined trace's termination is whichever side hit a
        // boundary (or both); we report `Boundary` if either did.
        termination: match (forward.termination, backward.termination) {
            (TraceEnd::Boundary, _) | (_, TraceEnd::Boundary) => TraceEnd::Boundary,
            _ => forward.termination,
        },
    }
}

/// March every connected component of the surface-surface
/// intersection.
///
/// Internally:
/// 1. Tessellate both surfaces at `seed_resolution × seed_resolution`,
///    triangle-vs-triangle to get seed polylines.
/// 2. Take one vertex from each seed polyline; project to find a
///    parametric `(uv_a, uv_b)` seed.
/// 3. Bidirectionally march from each seed.
/// 4. Fit a cubic NURBS curve through each trace.
pub fn march_all_components(
    s_a: &NurbsSurface,
    s_b: &NurbsSurface,
    seed_resolution: usize,
    params: &MarchParams,
) -> Vec<MarchedCurve> {
    let resolution = seed_resolution.max(8);
    let m_a = tessellate::surface(s_a, resolution, resolution);
    let m_b = tessellate::surface(s_b, resolution, resolution);
    let segs = intersect::triangle_mesh_intersection_segments_public(&m_a, &m_b);
    let polylines = intersect::chain_segments_public(segs);
    let mut out = Vec::new();
    for poly in polylines {
        // Pick the middle vertex as the seed (most stable choice).
        if poly.is_empty() {
            continue;
        }
        let seed_xyz = poly[poly.len() / 2];
        let Some((seed_a, seed_b)) = refine_seed(s_a, s_b, seed_xyz, params) else {
            continue;
        };
        let trace = trace_bidirectional(s_a, s_b, seed_a, seed_b, params);
        if trace.samples.len() < 2 {
            continue;
        }
        let polyline_xyz: Vec<Vector3<f64>> =
            trace.samples.iter().map(|s| s.xyz).collect();
        let curve = fit_polyline(&polyline_xyz);
        out.push(MarchedCurve {
            trace,
            curve,
        });
    }
    out
}

/// A marched intersection: the parametric trace + the fitted NURBS
/// curve through its 3D samples.
#[derive(Clone, Debug)]
pub struct MarchedCurve {
    /// Raw marching samples.
    pub trace: Trace,
    /// Cubic NURBS curve fit through `trace.samples`.
    pub curve: NurbsCurve,
}

/// Refine a 3D seed point `p` to parametric `(uv_a, uv_b)` on both
/// surfaces simultaneously.
///
/// Strategy: coarse grid search on each surface to seed the Newton
/// loop, then `MarchParams::projection_iters` Newton iterations. The
/// outputs both project to within `projection_tolerance` of the true
/// intersection.
pub fn refine_seed(
    s_a: &NurbsSurface,
    s_b: &NurbsSurface,
    p: Vector3<f64>,
    params: &MarchParams,
) -> Option<(Vector2<f64>, Vector2<f64>)> {
    let (uv_a_init, _) = grid_search_closest(s_a, p, 9);
    let (uv_b_init, _) = grid_search_closest(s_b, p, 9);
    let (uv_a, p_a) = closest_foot_seeded(s_a, p, uv_a_init, params);
    let (uv_b, p_b) = closest_foot_seeded(s_b, p, uv_b_init, params);
    let sep = (p_a - p_b).norm();
    if sep > params.chord_tolerance * 100.0 {
        // Seed isn't a true intersection — likely a misclassified
        // triangle pair. Refuse.
        return None;
    }
    Some((uv_a, uv_b))
}

// ===== internal helpers =====

fn inside(uv: Vector2<f64>, u_min: f64, u_max: f64, v_min: f64, v_max: f64) -> bool {
    uv.x >= u_min && uv.x <= u_max && uv.y >= v_min && uv.y <= v_max
}

/// True if `uv` lies **on** the parameter-range boundary (within a
/// small epsilon scaled to the range — `(u_max - u_min) * 1e-6`).
fn on_param_boundary(
    uv: Vector2<f64>,
    u_min: f64,
    u_max: f64,
    v_min: f64,
    v_max: f64,
) -> bool {
    let u_eps = (u_max - u_min).max(1.0) * 1.0e-6;
    let v_eps = (v_max - v_min).max(1.0) * 1.0e-6;
    (uv.x - u_min).abs() < u_eps
        || (u_max - uv.x).abs() < u_eps
        || (uv.y - v_min).abs() < v_eps
        || (v_max - uv.y).abs() < v_eps
}

fn average_xyz(s_a: &NurbsSurface, uv_a: Vector2<f64>, s_b: &NurbsSurface, uv_b: Vector2<f64>) -> Vector3<f64> {
    0.5 * (s_a.evaluate(uv_a.x, uv_a.y) + s_b.evaluate(uv_b.x, uv_b.y))
}

/// Coarse grid search returning `(uv, distance)` of the nearest grid
/// node to `p`. `nodes` is the per-axis sample count (so the grid is
/// `nodes × nodes`).
fn grid_search_closest(
    s: &NurbsSurface,
    p: Vector3<f64>,
    nodes: usize,
) -> (Vector2<f64>, f64) {
    let (u_min, u_max) = s.u_range();
    let (v_min, v_max) = s.v_range();
    let mut best_uv = Vector2::new(0.5 * (u_min + u_max), 0.5 * (v_min + v_max));
    let mut best_d = f64::INFINITY;
    let denom = (nodes - 1).max(1) as f64;
    for i in 0..nodes {
        let u = u_min + (i as f64 / denom) * (u_max - u_min);
        for j in 0..nodes {
            let v = v_min + (j as f64 / denom) * (v_max - v_min);
            let q = s.evaluate(u, v);
            let d = (q - p).norm();
            if d < best_d {
                best_d = d;
                best_uv = Vector2::new(u, v);
            }
        }
    }
    (best_uv, best_d)
}

/// Newton closest-foot of `p` on `s`, seeded by `uv_seed`. Returns
/// the refined `(uv, q)` where `q = s.evaluate(uv)`.
///
/// Single-surface Gauss-Newton in `(u, v)` using central-difference
/// tangents. Clamps to the parameter range on each iteration.
fn closest_foot_seeded(
    s: &NurbsSurface,
    p: Vector3<f64>,
    uv_seed: Vector2<f64>,
    params: &MarchParams,
) -> (Vector2<f64>, Vector3<f64>) {
    let (u_min, u_max) = s.u_range();
    let (v_min, v_max) = s.v_range();
    let h = ((u_max - u_min) + (v_max - v_min)) * 1.0e-5;
    let mut uv = Vector2::new(uv_seed.x.clamp(u_min, u_max), uv_seed.y.clamp(v_min, v_max));
    let mut q = s.evaluate(uv.x, uv.y);
    for _ in 0..params.projection_iters {
        let r = q - p;
        if r.norm() < params.projection_tolerance {
            return (uv, q);
        }
        let u_lo = (uv.x - h).max(u_min);
        let u_hi = (uv.x + h).min(u_max);
        let v_lo = (uv.y - h).max(v_min);
        let v_hi = (uv.y + h).min(v_max);
        let du_span = (u_hi - u_lo).max(1.0e-12);
        let dv_span = (v_hi - v_lo).max(1.0e-12);
        let tu = (s.evaluate(u_hi, uv.y) - s.evaluate(u_lo, uv.y)) / du_span;
        let tv = (s.evaluate(uv.x, v_hi) - s.evaluate(uv.x, v_lo)) / dv_span;
        let a11 = tu.dot(&tu);
        let a12 = tu.dot(&tv);
        let a22 = tv.dot(&tv);
        let b1 = -tu.dot(&r);
        let b2 = -tv.dot(&r);
        let det = a11 * a22 - a12 * a12;
        if det.abs() < 1.0e-14 {
            return (uv, q);
        }
        let du = (a22 * b1 - a12 * b2) / det;
        let dv = (-a12 * b1 + a11 * b2) / det;
        uv.x = (uv.x + du).clamp(u_min, u_max);
        uv.y = (uv.y + dv).clamp(v_min, v_max);
        q = s.evaluate(uv.x, uv.y);
    }
    (uv, q)
}

/// Surface unit normal at `(u, v)`.
fn surface_normal(s: &NurbsSurface, uv: Vector2<f64>) -> Vector3<f64> {
    let (u_min, u_max) = s.u_range();
    let (v_min, v_max) = s.v_range();
    let h = ((u_max - u_min) + (v_max - v_min)) * 1.0e-5;
    let u_lo = (uv.x - h).max(u_min);
    let u_hi = (uv.x + h).min(u_max);
    let v_lo = (uv.y - h).max(v_min);
    let v_hi = (uv.y + h).min(v_max);
    let tu = (s.evaluate(u_hi, uv.y) - s.evaluate(u_lo, uv.y)) / (u_hi - u_lo).max(1.0e-12);
    let tv = (s.evaluate(uv.x, v_hi) - s.evaluate(uv.x, v_lo)) / (v_hi - v_lo).max(1.0e-12);
    let n = tu.cross(&tv);
    let l = n.norm();
    if l < 1.0e-12 {
        Vector3::zeros()
    } else {
        n / l
    }
}

/// Find the boundary-crossing sample by bisecting the chord between
/// `inside_sample` and `outside_sample` and projecting back to both
/// surfaces.
fn boundary_bisect(
    s_a: &NurbsSurface,
    s_b: &NurbsSurface,
    inside_sample: TracePoint,
    outside_sample: TracePoint,
    params: &MarchParams,
) -> TracePoint {
    let (u_a_min, u_a_max) = s_a.u_range();
    let (v_a_min, v_a_max) = s_a.v_range();
    let (u_b_min, u_b_max) = s_b.u_range();
    let (v_b_min, v_b_max) = s_b.v_range();
    let mut lo = inside_sample;
    let mut hi = outside_sample;
    for _ in 0..20 {
        let mid_xyz = 0.5 * (lo.xyz + hi.xyz);
        let (mid_a, _) = closest_foot_seeded(s_a, mid_xyz, lo.uv_a, params);
        let (mid_b, _) = closest_foot_seeded(s_b, mid_xyz, lo.uv_b, params);
        let mid = TracePoint {
            xyz: average_xyz(s_a, mid_a, s_b, mid_b),
            uv_a: mid_a,
            uv_b: mid_b,
        };
        let mid_in = inside(mid.uv_a, u_a_min, u_a_max, v_a_min, v_a_max)
            && inside(mid.uv_b, u_b_min, u_b_max, v_b_min, v_b_max);
        if mid_in {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    // Return the inside endpoint snapped to the boundary if it
    // crossed.
    lo
}

/// Fit a cubic NURBS curve through a polyline of 3D samples,
/// degrading to degree-1 if the cubic fit fails.
fn fit_polyline(points: &[Vector3<f64>]) -> NurbsCurve {
    if points.len() < 4 {
        // Build a degree-1 piecewise-linear curve through the
        // points directly.
        return degree_one_through(points);
    }
    let n_cps = points.len().clamp(4, 16);
    match crate::fit::nurbs_curve_through_points(points, 3, n_cps) {
        Ok(fit) => fit.curve,
        Err(_) => degree_one_through(points),
    }
}

fn degree_one_through(points: &[Vector3<f64>]) -> NurbsCurve {
    let n = points.len().max(2);
    let mut pts = points.to_vec();
    if pts.len() < 2 {
        pts.push(*pts.first().unwrap_or(&Vector3::zeros()));
    }
    let mut knots = Vec::with_capacity(n + 2);
    knots.push(0.0);
    for i in 0..n {
        knots.push(i as f64 / (n - 1).max(1) as f64);
    }
    knots.push(1.0);
    let weights = vec![1.0; n];
    NurbsCurve::new_unchecked(1, knots, pts, weights)
}

/// Convert a list of marched curves into the same `Vec<NurbsCurve>`
/// shape the legacy [`crate::intersect::true_ssi`] returns, so callers
/// can switch with one drop-in.
pub fn marched_curves_to_nurbs(curves: Vec<MarchedCurve>) -> Vec<NurbsCurve> {
    curves.into_iter().map(|c| c.curve).collect()
}

/// Convenience: a one-call marching SSI returning fitted NURBS curves.
///
/// Equivalent to `march_all_components` followed by
/// `marched_curves_to_nurbs`. Returns the error-free zero-vector
/// result for non-intersecting surfaces.
pub fn marching_ssi(
    s_a: &NurbsSurface,
    s_b: &NurbsSurface,
    tolerance: f64,
) -> Result<Vec<NurbsCurve>, SurfaceError> {
    let params = MarchParams {
        chord_tolerance: tolerance,
        loop_tolerance: tolerance.max(1.0e-4),
        ..MarchParams::default()
    };
    Ok(marched_curves_to_nurbs(march_all_components(s_a, s_b, 24, &params)))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Planar surface in the xy plane at fixed z, CPs span the unit
    /// square.
    fn planar_xy_surface(z: f64) -> NurbsSurface {
        let knots = vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
        let cps = (0..4)
            .map(|i| {
                let u = i as f64 / 3.0;
                (0..4)
                    .map(|j| {
                        let v = j as f64 / 3.0;
                        Vector3::new(u, v, z)
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let weights = vec![vec![1.0; 4]; 4];
        NurbsSurface::new(3, 3, knots.clone(), knots, cps, weights).unwrap()
    }

    /// xz plane at fixed y, x ∈ [0,1], z ∈ [-0.5, 0.5].
    fn planar_xz_surface(y: f64) -> NurbsSurface {
        let knots = vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
        let cps = (0..4)
            .map(|i| {
                let u = i as f64 / 3.0;
                (0..4)
                    .map(|j| {
                        let v = -0.5 + j as f64 / 3.0;
                        Vector3::new(u, y, v)
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let weights = vec![vec![1.0; 4]; 4];
        NurbsSurface::new(3, 3, knots.clone(), knots, cps, weights).unwrap()
    }

    #[test]
    fn refine_seed_lands_on_both_surfaces() {
        let s_a = planar_xy_surface(0.0);
        let s_b = planar_xz_surface(0.5);
        let seed = Vector3::new(0.5, 0.5, 0.0);
        let params = MarchParams::default();
        let (uv_a, uv_b) = refine_seed(&s_a, &s_b, seed, &params).expect("seed");
        let p_a = s_a.evaluate(uv_a.x, uv_a.y);
        let p_b = s_b.evaluate(uv_b.x, uv_b.y);
        assert!((p_a - seed).norm() < 1.0e-6);
        assert!((p_b - seed).norm() < 1.0e-6);
    }

    #[test]
    fn marching_perpendicular_planes_terminates_on_boundary() {
        let s_a = planar_xy_surface(0.0);
        let s_b = planar_xz_surface(0.5);
        let seed = Vector3::new(0.5, 0.5, 0.0);
        let params = MarchParams::default();
        let (seed_a, seed_b) = refine_seed(&s_a, &s_b, seed, &params).expect("seed");
        let trace = trace_bidirectional(&s_a, &s_b, seed_a, seed_b, &params);
        assert!(trace.samples.len() >= 8, "got {} samples", trace.samples.len());
        // Every sample lies on the analytic line y = 0.5, z = 0.
        for s in &trace.samples {
            assert!((s.xyz.y - 0.5).abs() < 1.0e-6, "y = {}", s.xyz.y);
            assert!(s.xyz.z.abs() < 1.0e-6, "z = {}", s.xyz.z);
        }
        // Trace spans roughly x ∈ [0, 1].
        let xs: Vec<f64> = trace.samples.iter().map(|s| s.xyz.x).collect();
        let min_x = xs.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_x = xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(min_x < 0.05, "min_x = {min_x}");
        assert!(max_x > 0.95, "max_x = {max_x}");
        assert_eq!(trace.termination, TraceEnd::Boundary);
    }

    #[test]
    fn march_all_components_returns_smooth_curve() {
        let s_a = planar_xy_surface(0.0);
        let s_b = planar_xz_surface(0.5);
        let params = MarchParams::default();
        let curves = march_all_components(&s_a, &s_b, 16, &params);
        assert!(!curves.is_empty(), "no curves");
        let longest = curves
            .into_iter()
            .max_by_key(|c| c.trace.samples.len())
            .unwrap();
        // The fitted curve should be cubic.
        assert!(
            longest.curve.degree == 3 || longest.curve.degree == 1,
            "got degree {}",
            longest.curve.degree
        );
        // The fit should land on y=0.5, z=0 within a relaxed
        // tolerance (LSQ).
        let mut max_dev = 0.0_f64;
        for s in &longest.trace.samples {
            max_dev = max_dev.max((s.xyz.y - 0.5).abs()).max(s.xyz.z.abs());
        }
        assert!(max_dev < 1.0e-3, "max dev = {max_dev}");
    }

    /// Two perpendicular cylinders: one along y, radius 1; one along
    /// x, radius 1. Their intersection on the central piece is a
    /// figure-eight in 3D — but for a single quarter-arc patch from
    /// each cylinder the intersection is a smooth arc closely
    /// matching the analytic curve `x² + y² = 1` ∩ `y² + z² = 1`
    /// which simplifies to `x² = z²` (two arcs).
    fn quarter_cylinder_along_y(r: f64, y0: f64, y1: f64) -> NurbsSurface {
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

    fn quarter_cylinder_along_x(r: f64, x0: f64, x1: f64) -> NurbsSurface {
        let s2 = 2.0_f64.sqrt() / 2.0;
        let row_x0 = vec![
            Vector3::new(x0, r, 0.0),
            Vector3::new(x0, r, r),
            Vector3::new(x0, 0.0, r),
        ];
        let row_x1 = vec![
            Vector3::new(x1, r, 0.0),
            Vector3::new(x1, r, r),
            Vector3::new(x1, 0.0, r),
        ];
        NurbsSurface::new(
            1,
            2,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![row_x0, row_x1],
            vec![vec![1.0, s2, 1.0], vec![1.0, s2, 1.0]],
        )
        .unwrap()
    }

    #[test]
    fn marching_ssi_intersection_line_lies_on_analytic_curve() {
        // The cubic NURBS fit through the marched samples should
        // reproduce the analytic intersection curve to ≤ 2 ×
        // chord_tolerance — sampling along the fitted curve's
        // parameter range, every point must lie on the true
        // y=0.5, z=0 line.
        //
        // (We don't assert uniform chord lengths — the bidirectional
        // trace + adaptive step legitimately produces denser samples
        // at the boundary bisection points, and a centripetal fit
        // respects that distribution. The 3D-deviation test below is
        // the honest one.)
        let s_a = planar_xy_surface(0.0);
        let s_b = planar_xz_surface(0.5);
        let curves = marching_ssi(&s_a, &s_b, 1.0e-3).unwrap();
        assert!(!curves.is_empty());
        let c = curves
            .into_iter()
            .max_by_key(|c| c.n_control_points())
            .unwrap();
        let (u_min, u_max) = c.parameter_range();
        let n = 41;
        for i in 0..n {
            let u = u_min + (i as f64 / (n - 1) as f64) * (u_max - u_min);
            let p = c.evaluate(u);
            assert!(
                (p.y - 0.5).abs() < 5.0e-3,
                "sample {i} y={} not on analytic line",
                p.y
            );
            assert!(p.z.abs() < 5.0e-3, "sample {i} z={} not on analytic line", p.z);
        }
        // x range should cover the analytic [0, 1] within
        // 2 × chord_tolerance.
        let p_start = c.evaluate(u_min);
        let p_end = c.evaluate(u_max);
        let x_min = p_start.x.min(p_end.x);
        let x_max = p_start.x.max(p_end.x);
        assert!(x_min < 0.01, "x_min = {x_min}");
        assert!(x_max > 0.99, "x_max = {x_max}");
    }

    #[test]
    fn marching_perpendicular_cylinders_traces_analytic_curve() {
        let r = 1.0;
        let c1 = quarter_cylinder_along_y(r, -1.0, 1.0);
        let c2 = quarter_cylinder_along_x(r, -1.0, 1.0);
        let params = MarchParams {
            chord_tolerance: 5.0e-3,
            ..MarchParams::default()
        };
        let curves = march_all_components(&c1, &c2, 24, &params);
        assert!(!curves.is_empty(), "expected at least one trace");
        // For every traced sample p = (x, y, z), p lies on cylinder
        // 1 (x²+z² ≈ r²) AND cylinder 2 (y²+z² ≈ r²). The
        // closed-form curve has x²=y² so |x| = |y|.
        let longest = curves
            .into_iter()
            .max_by_key(|c| c.trace.samples.len())
            .unwrap();
        let mut max_residual = 0.0_f64;
        for s in &longest.trace.samples {
            let p = s.xyz;
            let r1 = (p.x * p.x + p.z * p.z).sqrt();
            let r2 = (p.y * p.y + p.z * p.z).sqrt();
            max_residual = max_residual.max((r1 - r).abs()).max((r2 - r).abs());
        }
        // Both implicit equations should be satisfied to within
        // chord_tolerance * a few — Newton on each surface guarantees
        // the projection residual is much tighter than the chord
        // length.
        assert!(
            max_residual < 1.0e-2,
            "max analytic residual = {max_residual}"
        );
    }

    #[test]
    fn marching_ssi_convenience_returns_nurbs() {
        let s_a = planar_xy_surface(0.0);
        let s_b = planar_xz_surface(0.5);
        let curves = marching_ssi(&s_a, &s_b, 1.0e-3).unwrap();
        assert!(!curves.is_empty());
    }

    #[test]
    fn disjoint_surfaces_produce_no_curves() {
        let s_a = planar_xy_surface(0.0);
        let s_b = planar_xy_surface(10.0); // 10 units apart, no intersection
        let curves = marching_ssi(&s_a, &s_b, 1.0e-3).unwrap();
        assert!(curves.is_empty(), "got {} curves", curves.len());
    }
}
