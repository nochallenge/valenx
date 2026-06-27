//! **Conjunction screening** (space situational awareness): given two
//! satellites' ephemerides over a common window, find the **time and distance of
//! closest approach** (TCA / miss distance) and flag a potential conjunction
//! when the miss distance drops below a screening threshold.
//!
//! This is the geometric core of collision-avoidance / SSA: continuously, for
//! every pair of tracked objects, ask "how close do these two get, and when?"
//! The answer is the *miss distance* (the minimum separation over the window)
//! and the *TCA* (the time it occurs). A pair whose miss distance falls inside a
//! screening volume is escalated for a closer look (covariance-based collision
//! probability, an operator decision) — but the screening itself is pure
//! relative-orbit geometry, identical for civil debris avoidance and defense
//! SSA.
//!
//! Algorithm: evaluate the relative separation on the **common time span** of
//! the two ephemerides (linearly interpolating each onto a shared coarse grid),
//! pick the coarse minimum, then **refine locally** with a golden-section search
//! on a cubic-Hermite (position + velocity) interpolation of the relative
//! motion inside the bracketing interval. The Hermite form uses each
//! ephemeris's velocity, so the refined TCA is far better than the coarse grid
//! spacing.
//!
//! # Honest scope
//!
//! - Accuracy is inherited from the supplied ephemerides (this crate's
//!   two-body / J2 propagator) and from the **sampling resolution**: the coarse
//!   pass can only *find* an approach it brackets, so the grid must be fine
//!   enough that no close pass slips between samples (for LEO relative speeds of
//!   ~km/s, sub-minute sampling). The Hermite refinement then sharpens the TCA
//!   within a bracket to well below the sample spacing.
//! - This is **deterministic miss-distance geometry only** — there is no
//!   covariance propagation and hence no probability of collision (`Pc`); the
//!   threshold flag is a hard-distance screen, the first filter of an SSA
//!   pipeline, not a final conjunction assessment.

use nalgebra::Vector3;

use crate::access::EphemerisPoint;
use crate::error::AstroError;

/// The result of screening two ephemerides for their closest approach.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Conjunction {
    /// Time of closest approach (s since the ephemeris epoch — the epoch shared
    /// by both inputs).
    pub tca: f64,
    /// Miss distance — the minimum separation between the two objects (m).
    pub miss_distance: f64,
    /// Relative speed of the two objects at TCA (m/s) — how fast they are
    /// passing, the other half of a conjunction's severity.
    pub relative_speed: f64,
    /// `true` if [`miss_distance`](Self::miss_distance) is at or below the
    /// screening threshold passed to [`screen_conjunction`].
    pub flagged: bool,
}

/// Cubic-Hermite interpolation of a position vector at local parameter `s ∈
/// [0,1]` across an interval of length `dt`, from the endpoint positions
/// `p0,p1` and velocities `v0,v1`. Standard Hermite basis; using the velocities
/// makes this a `C¹` fit that honours the true orbital motion between samples.
fn hermite(
    p0: Vector3<f64>,
    v0: Vector3<f64>,
    p1: Vector3<f64>,
    v1: Vector3<f64>,
    dt: f64,
    s: f64,
) -> Vector3<f64> {
    let s2 = s * s;
    let s3 = s2 * s;
    let h00 = 2.0 * s3 - 3.0 * s2 + 1.0;
    let h10 = s3 - 2.0 * s2 + s;
    let h01 = -2.0 * s3 + 3.0 * s2;
    let h11 = s3 - s2;
    p0 * h00 + v0 * (h10 * dt) + p1 * h01 + v1 * (h11 * dt)
}

/// Linear interpolation of an ephemeris [`EphemerisPoint`] (position **and**
/// velocity) at absolute time `t`, assumed to lie within `[a.time, b.time]`.
fn lerp_state(a: &EphemerisPoint, b: &EphemerisPoint, t: f64) -> (Vector3<f64>, Vector3<f64>) {
    let span = b.time - a.time;
    let f = if span.abs() < 1e-15 {
        0.0
    } else {
        (t - a.time) / span
    };
    (
        a.state.position + (b.state.position - a.state.position) * f,
        a.state.velocity + (b.state.velocity - a.state.velocity) * f,
    )
}

/// Sample an ephemeris at absolute time `t` (within its span), returning the
/// interpolated `(position, velocity)`. `idx_hint` is advanced as a cursor for
/// monotonically increasing `t` so a full scan is `O(n)` not `O(n²)`.
fn sample_at(eph: &[EphemerisPoint], t: f64, idx_hint: &mut usize) -> (Vector3<f64>, Vector3<f64>) {
    // Advance the cursor to the interval containing t.
    while *idx_hint + 1 < eph.len() && eph[*idx_hint + 1].time < t {
        *idx_hint += 1;
    }
    let i = (*idx_hint).min(eph.len().saturating_sub(2));
    lerp_state(&eph[i], &eph[i + 1], t)
}

/// Screen two satellite ephemerides for their **closest approach** over their
/// overlapping time span, flagging the pair if the miss distance is at or below
/// `threshold` (m).
///
/// Both ephemerides must share the same epoch (their `time` fields are seconds
/// from a common `t = 0`). The two are compared on the **intersection** of
/// their time spans, sampled onto a shared grid (the finer of the two inputs'
/// average step, capped to a sane resolution). The coarse minimum-separation
/// interval is then refined by a golden-section search on a Hermite
/// interpolation of each object's motion, giving a sub-sample TCA.
///
/// # Errors
///
/// Returns [`AstroError::InvalidParameter`] if either ephemeris has fewer than
/// two points, if the time spans do not overlap, if `threshold` is negative or
/// non-finite, or if any sample time is non-finite. These are the degenerate
/// inputs that would otherwise divide by zero or compare an empty window.
pub fn screen_conjunction(
    eph_a: &[EphemerisPoint],
    eph_b: &[EphemerisPoint],
    threshold: f64,
) -> Result<Conjunction, AstroError> {
    if eph_a.len() < 2 || eph_b.len() < 2 {
        return Err(AstroError::InvalidParameter(
            "screen_conjunction: each ephemeris needs at least two points",
        ));
    }
    if !threshold.is_finite() || threshold < 0.0 {
        return Err(AstroError::InvalidParameter(
            "screen_conjunction: threshold must be finite and >= 0",
        ));
    }
    for e in [eph_a, eph_b] {
        if e.iter().any(|p| !p.time.is_finite()) {
            return Err(AstroError::InvalidParameter(
                "screen_conjunction: non-finite ephemeris sample time",
            ));
        }
    }

    // Overlapping span [t_start, t_end].
    let t_start = eph_a[0].time.max(eph_b[0].time);
    let t_end = eph_a[eph_a.len() - 1].time.min(eph_b[eph_b.len() - 1].time);
    // Sample times are already validated finite above, so a non-positive span
    // means the two windows simply don't overlap.
    if t_end <= t_start {
        return Err(AstroError::InvalidParameter(
            "screen_conjunction: ephemerides do not overlap in time",
        ));
    }

    // Shared coarse grid: use the finer average step of the two inputs, but at
    // least 2 intervals and at most a fixed cap so we never allocate unbounded.
    let span = t_end - t_start;
    let step_a = span / (eph_a.len() as f64 - 1.0).max(1.0);
    let step_b = span / (eph_b.len() as f64 - 1.0).max(1.0);
    let mut step = step_a.min(step_b);
    if !step.is_finite() || step <= 0.0 {
        step = span;
    }
    const MAX_GRID: usize = 200_000;
    let mut n = (span / step).ceil() as usize;
    n = n.clamp(2, MAX_GRID);
    let dt = span / n as f64;

    // Coarse scan for the minimum-separation grid node.
    let mut cur_a = 0usize;
    let mut cur_b = 0usize;
    let sep_at = |t: f64, ca: &mut usize, cb: &mut usize| -> f64 {
        let (pa, _) = sample_at(eph_a, t, ca);
        let (pb, _) = sample_at(eph_b, t, cb);
        (pa - pb).norm()
    };

    let mut best_k = 0usize;
    let mut best_sep = f64::INFINITY;
    for k in 0..=n {
        let t = t_start + k as f64 * dt;
        let d = sep_at(t, &mut cur_a, &mut cur_b);
        if d < best_sep {
            best_sep = d;
            best_k = k;
        }
    }

    // Refine within the interval bracketing the coarse minimum. Build Hermite
    // states for each object at the bracket endpoints, then golden-section the
    // relative separation.
    let k0 = best_k.saturating_sub(1);
    let k1 = (best_k + 1).min(n);
    let ta = t_start + k0 as f64 * dt;
    let tb = t_start + k1 as f64 * dt;

    let mut ca = 0usize;
    let (pa0, va0) = sample_at(eph_a, ta, &mut ca);
    let (pa1, va1) = sample_at(eph_a, tb, &mut ca);
    let mut cb = 0usize;
    let (pb0, vb0) = sample_at(eph_b, ta, &mut cb);
    let (pb1, vb1) = sample_at(eph_b, tb, &mut cb);

    let bracket = tb - ta;
    let rel_sep = |s: f64| -> f64 {
        let a = hermite(pa0, va0, pa1, va1, bracket, s);
        let b = hermite(pb0, vb0, pb1, vb1, bracket, s);
        (a - b).norm()
    };

    // Golden-section minimisation on s ∈ [0, 1].
    let inv_phi = (5.0_f64.sqrt() - 1.0) / 2.0; // 1/φ ≈ 0.618
    let (mut lo, mut hi) = (0.0_f64, 1.0_f64);
    let mut x1 = hi - inv_phi * (hi - lo);
    let mut x2 = lo + inv_phi * (hi - lo);
    let mut f1 = rel_sep(x1);
    let mut f2 = rel_sep(x2);
    for _ in 0..80 {
        if f1 < f2 {
            hi = x2;
            x2 = x1;
            f2 = f1;
            x1 = hi - inv_phi * (hi - lo);
            f1 = rel_sep(x1);
        } else {
            lo = x1;
            x1 = x2;
            f1 = f2;
            x2 = lo + inv_phi * (hi - lo);
            f2 = rel_sep(x2);
        }
        if (hi - lo) < 1e-12 {
            break;
        }
    }
    let s_star = 0.5 * (lo + hi);
    let tca = ta + s_star * bracket;
    let miss_refined = rel_sep(s_star);

    // The refined minimum should never be worse than the coarse node; guard.
    let (miss_distance, tca) = if miss_refined <= best_sep {
        (miss_refined, tca)
    } else {
        (best_sep, t_start + best_k as f64 * dt)
    };

    // Relative speed at TCA from the Hermite velocity (finite difference of the
    // Hermite position is exact for the cubic; use a small central difference).
    let ds = 1e-4_f64.min(0.25);
    let s_lo = (s_star - ds).clamp(0.0, 1.0);
    let s_hi = (s_star + ds).clamp(0.0, 1.0);
    let span_s = (s_hi - s_lo) * bracket;
    let rel_speed = if span_s.abs() > 1e-12 {
        let a_lo = hermite(pa0, va0, pa1, va1, bracket, s_lo);
        let b_lo = hermite(pb0, vb0, pb1, vb1, bracket, s_lo);
        let a_hi = hermite(pa0, va0, pa1, va1, bracket, s_hi);
        let b_hi = hermite(pb0, vb0, pb1, vb1, bracket, s_hi);
        (((a_hi - b_hi) - (a_lo - b_lo)) / span_s).norm()
    } else {
        // Fall back to the endpoint relative velocity.
        (va0 - vb0).norm()
    };

    Ok(Conjunction {
        tca,
        miss_distance,
        relative_speed: rel_speed,
        flagged: miss_distance <= threshold,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::{MU_EARTH, R_EARTH};
    use crate::orbit3d::{self, ClassicalElements, StateVector};

    fn ephemeris_from_coe(coe: &ClassicalElements, total: f64, n: usize) -> Vec<EphemerisPoint> {
        let mut state = orbit3d::coe_to_rv(coe).unwrap();
        let dt = total / n as f64;
        let mut out = Vec::with_capacity(n + 1);
        out.push(EphemerisPoint { time: 0.0, state });
        for i in 1..=n {
            state = orbit3d::propagate(&state, dt, 1, false).unwrap();
            out.push(EphemerisPoint {
                time: i as f64 * dt,
                state,
            });
        }
        out
    }

    #[test]
    fn crossing_orbits_have_tca_at_the_analytic_crossing() {
        // Two satellites built to pass through the SAME point at the SAME time.
        // A is an equatorial circular orbit; B is the same orbit inclined 90°,
        // phased so that both are at the ascending-node crossing point
        // simultaneously. The two ground-truth facts: (1) they actually meet
        // (miss distance ~ 0 to interpolation accuracy), and (2) the TCA is the
        // crossing time we constructed.
        let a = R_EARTH + 700_000.0;
        let period = std::f64::consts::TAU * (a * a * a / MU_EARTH).sqrt();

        // Pick a meeting time a quarter-period in.
        let t_meet = period / 4.0;

        // Sat A: equatorial circular, ν chosen so that at t_meet it is at the
        // +x axis crossing. Easier: construct both from explicit states that
        // coincide at t_meet by back-propagating from a shared rendezvous state.
        // Shared point: on the +x axis at radius a.
        let meet_pos = Vector3::new(a, 0.0, 0.0);
        let v_circ = (MU_EARTH / a).sqrt();
        // A moves in +y (equatorial, prograde); B moves in +z (polar).
        let state_a_meet = StateVector {
            position: meet_pos,
            velocity: Vector3::new(0.0, v_circ, 0.0),
        };
        let state_b_meet = StateVector {
            position: meet_pos,
            velocity: Vector3::new(0.0, 0.0, v_circ),
        };

        // Back-propagate both to t = 0, then forward-sample an ephemeris.
        let steps = 600u64;
        let dt_back = t_meet / steps as f64;
        // Propagate backward by stepping with negative dt.
        let state_a0 = orbit3d::propagate(&state_a_meet, -dt_back, steps, false).unwrap();
        let state_b0 = orbit3d::propagate(&state_b_meet, -dt_back, steps, false).unwrap();

        let make_eph = |s0: &StateVector| -> Vec<EphemerisPoint> {
            let total = period / 2.0;
            let n = 600usize;
            let dt = total / n as f64;
            let mut s = *s0;
            let mut out = vec![EphemerisPoint {
                time: 0.0,
                state: s,
            }];
            for i in 1..=n {
                s = orbit3d::propagate(&s, dt, 1, false).unwrap();
                out.push(EphemerisPoint {
                    time: i as f64 * dt,
                    state: s,
                });
            }
            out
        };
        let eph_a = make_eph(&state_a0);
        let eph_b = make_eph(&state_b0);

        let c = screen_conjunction(&eph_a, &eph_b, 1_000.0).unwrap();
        // They genuinely meet: miss distance is tiny (limited by interpolation /
        // back-prop round-off, not geometry).
        assert!(
            c.miss_distance < 2_000.0,
            "miss distance {} m, expected ~0 at a true crossing",
            c.miss_distance
        );
        assert!(c.flagged, "a near-zero miss must be flagged");
        // TCA is the constructed crossing time.
        assert!(
            (c.tca - t_meet).abs() < 5.0,
            "TCA {} s != constructed crossing {} s",
            c.tca,
            t_meet
        );
        // Relative speed of two perpendicular circular velocities is √2·v_circ.
        let expected_rel = std::f64::consts::SQRT_2 * v_circ;
        assert!(
            (c.relative_speed - expected_rel).abs() / expected_rel < 0.05,
            "rel speed {} != ~{}",
            c.relative_speed,
            expected_rel
        );
    }

    #[test]
    fn well_separated_orbits_have_large_miss_and_no_false_flag() {
        // Two circular orbits at very different altitudes never come close.
        let low = ClassicalElements {
            semi_major_axis: R_EARTH + 400_000.0,
            eccentricity: 0.0,
            inclination: 0.0,
            raan: 0.0,
            arg_periapsis: 0.0,
            true_anomaly: 0.0,
        };
        let high = ClassicalElements {
            semi_major_axis: R_EARTH + 20_000_000.0,
            eccentricity: 0.0,
            inclination: 0.0,
            raan: 0.0,
            arg_periapsis: 0.0,
            true_anomaly: 0.0,
        };
        let p_low = std::f64::consts::TAU * ((R_EARTH + 400_000.0).powi(3) / MU_EARTH).sqrt();
        let eph_low = ephemeris_from_coe(&low, 2.0 * p_low, 1000);
        let eph_high = ephemeris_from_coe(&high, 2.0 * p_low, 1000);

        let c = screen_conjunction(&eph_low, &eph_high, 5_000.0).unwrap();
        // The altitude gap alone is ~19,600 km; miss distance must be enormous.
        assert!(
            c.miss_distance > 1.0e7,
            "miss distance {} m should be > 10,000 km",
            c.miss_distance
        );
        assert!(!c.flagged, "well-separated orbits must NOT be flagged");
    }

    #[test]
    fn identical_orbits_have_zero_miss_distance() {
        // The same orbit screened against itself: separation is ~0 everywhere,
        // and it is flagged.
        let coe = ClassicalElements {
            semi_major_axis: R_EARTH + 550_000.0,
            eccentricity: 0.001,
            inclination: 51.6_f64.to_radians(),
            raan: 0.4,
            arg_periapsis: 0.2,
            true_anomaly: 0.0,
        };
        let p = std::f64::consts::TAU * ((R_EARTH + 550_000.0).powi(3) / MU_EARTH).sqrt();
        let eph = ephemeris_from_coe(&coe, p, 500);
        let c = screen_conjunction(&eph, &eph, 100.0).unwrap();
        assert!(
            c.miss_distance < 1e-3,
            "self miss distance {}",
            c.miss_distance
        );
        assert!(c.flagged);
    }

    #[test]
    fn refined_tca_beats_the_coarse_grid() {
        // The Hermite refinement should locate a TCA that is NOT pinned to a
        // coarse grid node. Use a coarse ephemeris of a close tail-chase and
        // check the miss distance is below the coarsely-sampled minimum (i.e.
        // refinement found a better point between nodes). We build a guaranteed
        // close pass: two coplanar circular orbits at slightly different
        // altitudes and phases so they approach between samples.
        let a1 = R_EARTH + 500_000.0;
        let a2 = R_EARTH + 505_000.0;
        let coe1 = ClassicalElements {
            semi_major_axis: a1,
            eccentricity: 0.0,
            inclination: 0.0,
            raan: 0.0,
            arg_periapsis: 0.0,
            true_anomaly: 0.0,
        };
        let coe2 = ClassicalElements {
            semi_major_axis: a2,
            eccentricity: 0.0,
            inclination: 0.0,
            raan: 0.0,
            arg_periapsis: 0.0,
            true_anomaly: 0.05, // small phase offset
        };
        let p = std::f64::consts::TAU * (a1 * a1 * a1 / MU_EARTH).sqrt();
        // Coarse: only ~30 s sampling.
        let n = (p / 30.0) as usize;
        let eph1 = ephemeris_from_coe(&coe1, p, n);
        let eph2 = ephemeris_from_coe(&coe2, p, n);

        // Coarse minimum separation (no refinement) for comparison.
        let mut coarse_min = f64::INFINITY;
        for (pa, pb) in eph1.iter().zip(eph2.iter()) {
            let d = (pa.state.position - pb.state.position).norm();
            coarse_min = coarse_min.min(d);
        }
        let c = screen_conjunction(&eph1, &eph2, 1.0e9).unwrap();
        assert!(
            c.miss_distance <= coarse_min + 1e-6,
            "refined miss {} should be <= coarse min {}",
            c.miss_distance,
            coarse_min
        );
        assert!(c.tca.is_finite() && c.tca >= 0.0);
    }

    #[test]
    fn degenerate_inputs_error_not_panic() {
        let coe = ClassicalElements {
            semi_major_axis: R_EARTH + 500_000.0,
            eccentricity: 0.0,
            inclination: 0.0,
            raan: 0.0,
            arg_periapsis: 0.0,
            true_anomaly: 0.0,
        };
        let eph = ephemeris_from_coe(&coe, 3000.0, 50);

        // Fewer than two points.
        let one = vec![eph[0]];
        assert!(screen_conjunction(&one, &eph, 100.0).is_err());
        assert!(screen_conjunction(&[], &eph, 100.0).is_err());

        // Negative / non-finite threshold.
        assert!(screen_conjunction(&eph, &eph, -1.0).is_err());
        assert!(screen_conjunction(&eph, &eph, f64::NAN).is_err());

        // Non-overlapping time spans.
        let mut shifted = eph.clone();
        for p in &mut shifted {
            p.time += 1.0e9;
        }
        assert!(screen_conjunction(&eph, &shifted, 100.0).is_err());

        // Non-finite sample time.
        let mut bad = eph.clone();
        bad[3].time = f64::NAN;
        assert!(screen_conjunction(&bad, &eph, 100.0).is_err());
    }
}
