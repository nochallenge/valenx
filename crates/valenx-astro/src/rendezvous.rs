//! Relative orbital motion and rendezvous via the **Clohessy–Wiltshire
//! (Hill's) equations** — the linearised dynamics of a chaser relative to
//! a target on a circular reference orbit, expressed in the target's
//! **LVLH** frame (x = radial / "up", y = along-track / velocity
//! direction, z = cross-track / orbit-normal).
//!
//! For a target in a circular orbit of semi-major axis `a`, the **mean
//! motion** is `n = √(μ/a³)`, and the relative state evolves by the
//! closed-form **state-transition matrix** `Φ(t)`:
//!
//! ```text
//!   [ r(t) ]   [ Φrr(t)  Φrv(t) ] [ r₀ ]
//!   [ v(t) ] = [ Φvr(t)  Φvv(t) ] [ v₀ ]
//! ```
//!
//! whose four 3×3 blocks are the textbook sine/cosine expressions in
//! `nt`. From these falls out the classic result that a chaser placed at
//! a radial offset with the bounded-orbit velocity (`ẏ₀ = −2n·x₀`) traces
//! a **2:1 along-track-by-radial relative ellipse** ("football"), and the
//! **two-impulse rendezvous** that drives the relative position to a
//! target point in a chosen transfer time `T`:
//!
//! ```text
//!   Δv₁ = Φrv(T)⁻¹·(r_target − Φrr(T)·r₀) − v₀
//!   Δv₂ = v_target − ( Φvr(T)·r₀ + Φvv(T)·(v₀+Δv₁) )
//! ```
//!
//! `Φrv(T)` is **singular whenever `sin(nT) = 0`** — i.e. at every whole
//! number of half-periods, `nT = kπ` (`k = 1, 2, 3, …`): the cross-track
//! entry `sin(nT)/n` vanishes there, so the block cannot be inverted. The
//! in-plane (x–y) sub-block is *additionally* degenerate at the full
//! periods `nT = 2kπ`. The guard rejects the whole `nT = kπ` family, where
//! the boundary-value problem has no unique impulsive solution, and
//! [`two_impulse_rendezvous`] returns [`AstroError::NonConvergent`] rather
//! than inverting a singular matrix into garbage.
//!
//! All the relations are exact closed forms (built on `nalgebra`), so the
//! tests pin the STM, the 2:1 ellipse and the two-impulse
//! self-consistency (the solved burns, propagated by the same STM, null
//! the relative position to well under a millimetre) directly.

use nalgebra::{Matrix3, Vector3};
use serde::{Deserialize, Serialize};

use crate::constants::MU_EARTH;
use crate::error::{AstroError, Result};

/// Below this `|sin(nT)|` the position-to-velocity block `Φrv` is treated
/// as singular (its cross-track entry is `sin(nT)/n`, and the in-plane
/// block degenerates at the full period). Chosen well above floating
/// round-off so a near-singular transfer is rejected before the inverse
/// blows up.
const SIN_SINGULAR_TOL: f64 = 1e-6;

/// Mean motion `n = √(μ/a³)` (rad/s) of a circular reference orbit of
/// semi-major axis `a` (m), about Earth.
///
/// # Errors
///
/// Returns [`AstroError::InvalidParameter`] if `a` is non-finite or
/// non-positive (which would make `μ/a³` or its square root a silent
/// `NaN`/`Inf`).
pub fn mean_motion(semi_major_axis: f64) -> Result<f64> {
    if !semi_major_axis.is_finite() || semi_major_axis <= 0.0 {
        return Err(AstroError::InvalidParameter(
            "target semi-major axis must be finite and > 0",
        ));
    }
    Ok((MU_EARTH / semi_major_axis.powi(3)).sqrt())
}

/// The **synodic period** `T_syn = 2π / |n₁ − n₂|` (s) between two circular orbits
/// of semi-major axes `a1` and `a2` (m) — the time between successive conjunctions
/// (the two bodies lining up again), where `nₖ = √(μ/aₖ³)` is each orbit's
/// [`mean_motion`]. The faster (inner) orbit gains a full lap of phase on the slower
/// one over exactly this interval, so it is the *beat period* of the two orbital
/// frequencies — equivalently `1/T_syn = |1/T₁ − 1/T₂|` in terms of the orbital
/// periods `Tₖ = 2π/nₖ`. For two nearby orbits the beat is slow (`T_syn → ∞` as
/// `a1 → a2`, returned as `+∞` from the `2π/0` limit); for widely separated orbits
/// it approaches the shorter period.
///
/// # Errors
///
/// Returns [`AstroError::InvalidParameter`] (via [`mean_motion`]) if either
/// semi-major axis is non-finite or non-positive.
pub fn synodic_period(a1: f64, a2: f64) -> Result<f64> {
    Ok(std::f64::consts::TAU / (mean_motion(a1)? - mean_motion(a2)?).abs())
}

/// The four 3×3 blocks `(Φrr, Φrv, Φvr, Φvv)` of the Clohessy–Wiltshire
/// state-transition matrix at elapsed time `t` (s) for mean motion `n`
/// (rad/s).
///
/// # Errors
///
/// Returns [`AstroError::InvalidParameter`] if `n` is non-finite or
/// non-positive, or `t` is non-finite.
#[allow(clippy::type_complexity)]
pub fn cw_state_transition(
    n: f64,
    t: f64,
) -> Result<(Matrix3<f64>, Matrix3<f64>, Matrix3<f64>, Matrix3<f64>)> {
    if !n.is_finite() || n <= 0.0 {
        return Err(AstroError::InvalidParameter("mean motion must be > 0"));
    }
    if !t.is_finite() {
        return Err(AstroError::InvalidParameter("time must be finite"));
    }
    let s = (n * t).sin();
    let c = (n * t).cos();
    let nt = n * t;

    // Φrr: position response to initial position.
    let phi_rr = Matrix3::new(
        4.0 - 3.0 * c,
        0.0,
        0.0,
        6.0 * (s - nt),
        1.0,
        0.0,
        0.0,
        0.0,
        c,
    );
    // Φrv: position response to initial velocity.
    let phi_rv = Matrix3::new(
        s / n,
        2.0 * (1.0 - c) / n,
        0.0,
        -2.0 * (1.0 - c) / n,
        (4.0 * s - 3.0 * nt) / n,
        0.0,
        0.0,
        0.0,
        s / n,
    );
    // Φvr: velocity response to initial position.
    let phi_vr = Matrix3::new(
        3.0 * n * s,
        0.0,
        0.0,
        6.0 * n * (c - 1.0),
        0.0,
        0.0,
        0.0,
        0.0,
        -n * s,
    );
    // Φvv: velocity response to initial velocity.
    let phi_vv = Matrix3::new(c, 2.0 * s, 0.0, -2.0 * s, 4.0 * c - 3.0, 0.0, 0.0, 0.0, c);
    Ok((phi_rr, phi_rv, phi_vr, phi_vv))
}

/// Propagate a relative LVLH state `(r₀, v₀)` (m, m/s) forward by time `t`
/// (s) under the Clohessy–Wiltshire dynamics for mean motion `n`.
///
/// Returns `(r(t), v(t))`.
///
/// # Errors
///
/// As [`cw_state_transition`].
pub fn cw_propagate(
    n: f64,
    t: f64,
    r0: Vector3<f64>,
    v0: Vector3<f64>,
) -> Result<(Vector3<f64>, Vector3<f64>)> {
    if !r0.iter().all(|x| x.is_finite()) || !v0.iter().all(|x| x.is_finite()) {
        return Err(AstroError::InvalidParameter(
            "relative state must be finite",
        ));
    }
    let (rr, rv, vr, vv) = cw_state_transition(n, t)?;
    let r = rr * r0 + rv * v0;
    let v = vr * r0 + vv * v0;
    Ok((r, v))
}

/// The Δv plan of a two-impulse Clohessy–Wiltshire rendezvous.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct RendezvousPlan {
    /// First burn (m/s, LVLH): applied at `t = 0` to put the chaser on a
    /// trajectory that arrives at `r_target` at `t = T`.
    pub delta_v1: Vector3<f64>,
    /// Second burn (m/s, LVLH): applied at arrival to match the target
    /// relative velocity (nulls the approach for `v_target = 0`).
    pub delta_v2: Vector3<f64>,
    /// Total Δv magnitude `|Δv₁| + |Δv₂|` (m/s).
    pub total_delta_v: f64,
    /// Transfer time `T` (s).
    pub transfer_time: f64,
}

/// Solve the **two-impulse rendezvous** that moves the chaser from
/// relative state `(r0, v0)` to relative position `r_target` with relative
/// velocity `v_target` over a transfer time `transfer_time` (s), under the
/// Clohessy–Wiltshire dynamics for mean motion `n`.
///
/// # Errors
///
/// - [`AstroError::InvalidParameter`] for a non-finite or non-positive `n`
///   / `transfer_time`, or a non-finite state vector.
/// - [`AstroError::NonConvergent`] when `Φrv(T)` is singular —
///   `|sin(nT)| < SIN_SINGULAR_TOL` — i.e. the transfer time is a
///   whole number of half-periods (`nT = kπ`), where the boundary-value
///   problem has no unique impulsive solution.
pub fn two_impulse_rendezvous(
    n: f64,
    transfer_time: f64,
    r0: Vector3<f64>,
    v0: Vector3<f64>,
    r_target: Vector3<f64>,
    v_target: Vector3<f64>,
) -> Result<RendezvousPlan> {
    if !n.is_finite() || n <= 0.0 {
        return Err(AstroError::InvalidParameter("mean motion must be > 0"));
    }
    if !transfer_time.is_finite() || transfer_time <= 0.0 {
        return Err(AstroError::InvalidParameter(
            "transfer_time must be finite and > 0",
        ));
    }
    for vvec in [&r0, &v0, &r_target, &v_target] {
        if !vvec.iter().all(|x| x.is_finite()) {
            return Err(AstroError::InvalidParameter(
                "rendezvous state vectors must be finite",
            ));
        }
    }

    // Φrv is singular when sin(nT) = 0; reject before inverting.
    if (n * transfer_time).sin().abs() < SIN_SINGULAR_TOL {
        return Err(AstroError::NonConvergent(
            "CW Phi_rv singular: transfer time is a whole number of half-periods (sin(nT)=0)",
        ));
    }

    let (phi_rr, phi_rv, phi_vr, phi_vv) = cw_state_transition(n, transfer_time)?;
    let phi_rv_inv = phi_rv.try_inverse().ok_or(AstroError::NonConvergent(
        "CW Phi_rv is not invertible at this transfer time",
    ))?;

    // Required post-burn velocity to arrive at r_target: solve
    // r_target = Φrr·r0 + Φrv·v0_plus  =>  v0_plus = Φrv⁻¹·(r_target − Φrr·r0).
    let v0_plus = phi_rv_inv * (r_target - phi_rr * r0);
    let delta_v1 = v0_plus - v0;

    // Arrival relative velocity before the second burn.
    let v_arrival = phi_vr * r0 + phi_vv * v0_plus;
    let delta_v2 = v_target - v_arrival;

    Ok(RendezvousPlan {
        delta_v1,
        delta_v2,
        total_delta_v: delta_v1.norm() + delta_v2.norm(),
        transfer_time,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::R_EARTH;
    use std::f64::consts::TAU;

    /// Mean motion of a ~400 km circular orbit.
    fn n_leo() -> f64 {
        mean_motion(R_EARTH + 400_000.0).expect("valid")
    }

    #[test]
    fn mean_motion_matches_circular_period() {
        let a = R_EARTH + 400_000.0;
        let n = mean_motion(a).expect("valid");
        // Period 2π/n must equal the Keplerian circular period.
        let period = TAU / n;
        let expected = TAU * (a.powi(3) / MU_EARTH).sqrt();
        assert!((period - expected).abs() / expected < 1e-12);
        // ~400 km LEO period is ~92.6 min ≈ 5554 s.
        assert!((period - 5_554.0).abs() < 10.0, "period {period} s");
    }

    #[test]
    fn synodic_period_is_the_beat_of_two_mean_motions() {
        let (a1, a2) = (7.0e6, 7.5e6);
        let n1 = mean_motion(a1).expect("valid");
        let n2 = mean_motion(a2).expect("valid");
        let t_syn = synodic_period(a1, a2).expect("valid");
        // Definition: T_syn = 2π/|n1 − n2|.
        assert!(
            (t_syn - TAU / (n1 - n2).abs()).abs() / t_syn < 1e-12,
            "T_syn = 2π/|Δn|"
        );
        // Beat-of-periods identity 1/T_syn = |1/P1 − 1/P2| (Pk = 2π/nk) — an
        // independent derivation threading mean_motion.
        let (p1, p2) = (TAU / n1, TAU / n2);
        assert!(
            (1.0 / t_syn - (1.0 / p1 - 1.0 / p2).abs()).abs() * t_syn < 1e-12,
            "1/T_syn = |1/P1 − 1/P2|"
        );
        // For nearby orbits the beat is slow — longer than either orbital period.
        assert!(
            t_syn > p1 && t_syn > p2,
            "close orbits → T_syn exceeds both periods"
        );
        // Symmetric in its arguments.
        assert!(
            (synodic_period(a2, a1).expect("valid") - t_syn).abs() / t_syn < 1e-12,
            "symmetric"
        );
        // Identical orbits never re-conjunct → infinite synodic period (2π/0).
        assert!(
            synodic_period(a1, a1).expect("valid").is_infinite(),
            "equal orbits → ∞"
        );
        // A non-positive semi-major axis is rejected (propagated from mean_motion).
        assert!(synodic_period(-1.0, a2).is_err(), "a ≤ 0 → Err");
    }

    #[test]
    fn stm_is_identity_at_t0() {
        let n = n_leo();
        let (rr, rv, vr, vv) = cw_state_transition(n, 0.0).expect("valid");
        // At t=0: Φrr = I, Φvv = I, Φrv = 0, Φvr = 0.
        assert!((rr - Matrix3::identity()).norm() < 1e-12);
        assert!((vv - Matrix3::identity()).norm() < 1e-12);
        assert!(rv.norm() < 1e-12);
        assert!(vr.norm() < 1e-12);
    }

    #[test]
    fn stm_reduces_to_identity_after_full_period_in_position() {
        // After one full period the in-plane position returns (Φrr's
        // bounded terms are periodic); the along-track secular drift only
        // appears via initial velocity, so a zero-velocity radial state
        // recurs. Check Φrr at nt=2π: 4−3cos=1, 6(sin−2π) drift term, c=1.
        let n = n_leo();
        let t = TAU / n;
        let (rr, _, _, _) = cw_state_transition(n, t).expect("valid");
        // Diagonal radial/cross-track terms return to 1.
        assert!((rr[(0, 0)] - 1.0).abs() < 1e-9);
        assert!((rr[(2, 2)] - 1.0).abs() < 1e-9);
    }

    #[test]
    fn pure_along_track_offset_has_no_radial_motion() {
        // A chaser purely ahead/behind with zero relative velocity stays
        // at zero radial offset and drifts along-track (the classic
        // station-keeping drift) — a clean CW signature.
        let n = n_leo();
        let r0 = Vector3::new(0.0, 1_000.0, 0.0);
        let v0 = Vector3::zeros();
        let period = TAU / n;
        for k in 1..=8 {
            let t = period * (k as f64) / 8.0;
            let (r, _) = cw_propagate(n, t, r0, v0).expect("valid");
            assert!(r.x.abs() < 1e-6, "radial offset {} at t={t}", r.x);
        }
    }

    #[test]
    fn bounded_relative_orbit_is_a_2_to_1_football() {
        // ORACLE: a radial offset x0 with the bounded-orbit velocity
        // ẏ0 = −2n·x0 traces a centred relative ellipse whose along-track
        // semi-axis is exactly twice the radial — the 2:1 "football".
        let n = n_leo();
        let rho = 500.0;
        let r0 = Vector3::new(rho, 0.0, 0.0);
        let v0 = Vector3::new(0.0, -2.0 * n * rho, 0.0);
        let period = TAU / n;

        let (mut xmin, mut xmax) = (f64::INFINITY, f64::NEG_INFINITY);
        let (mut ymin, mut ymax) = (f64::INFINITY, f64::NEG_INFINITY);
        for i in 0..=4_000 {
            let t = period * (i as f64) / 4_000.0;
            let (r, _) = cw_propagate(n, t, r0, v0).expect("valid");
            xmin = xmin.min(r.x);
            xmax = xmax.max(r.x);
            ymin = ymin.min(r.y);
            ymax = ymax.max(r.y);
        }
        let radial_extent = xmax - xmin;
        let along_extent = ymax - ymin;
        // Radial peak-to-peak is 2ρ; along-track is 4ρ -> ratio 2.000.
        assert!(
            (radial_extent - 2.0 * rho).abs() < 1.0,
            "radial {radial_extent}"
        );
        assert!(
            (along_extent / radial_extent - 2.0).abs() < 1e-3,
            "ratio {}",
            along_extent / radial_extent
        );
        // Bounded: returns to the start after one period (no secular drift).
        let (r_end, v_end) = cw_propagate(n, period, r0, v0).expect("valid");
        assert!(
            (r_end - r0).norm() < 1e-3,
            "closure {}",
            (r_end - r0).norm()
        );
        assert!((v_end - v0).norm() < 1e-6);
    }

    #[test]
    fn two_impulse_rendezvous_nulls_relative_position() {
        // ORACLE / self-consistency: solve the two-impulse transfer, then
        // propagate (r0, v0+Δv1) with the SAME STM. The relative position
        // at T must equal r_target to < 1 mm, and the arrival velocity
        // after Δv2 must equal v_target to machine precision.
        let n = n_leo();
        let r0 = Vector3::new(-2_000.0, 5_000.0, 100.0);
        let v0 = Vector3::new(1.0, -2.0, 0.5);
        let r_target = Vector3::zeros(); // dock at the origin
        let v_target = Vector3::zeros(); // null the approach
        let period = TAU / n;
        let t = 0.25 * period; // quarter-period transfer (sin(nT)=sin(π/2)=1)

        let plan = two_impulse_rendezvous(n, t, r0, v0, r_target, v_target).expect("solvable");

        // Propagate with the first burn applied.
        let (r_arr, v_arr_minus) = cw_propagate(n, t, r0, v0 + plan.delta_v1).expect("valid");
        assert!(
            (r_arr - r_target).norm() < 1e-3,
            "miss = {} m",
            (r_arr - r_target).norm()
        );
        // The second burn nulls the arrival velocity to the target.
        let v_after = v_arr_minus + plan.delta_v2;
        assert!(
            (v_after - v_target).norm() < 1e-9,
            "residual v = {}",
            (v_after - v_target).norm()
        );
        assert!(plan.total_delta_v > 0.0);
        assert!((plan.transfer_time - t).abs() < 1e-12);
    }

    #[test]
    fn two_impulse_singular_at_whole_half_periods() {
        // Φrv is singular at nT = kπ (sin = 0): a transfer time of an
        // integer number of half-periods must return NonConvergent, not a
        // garbage inverse.
        let n = n_leo();
        let period = TAU / n;
        let r0 = Vector3::new(100.0, 200.0, 0.0);
        let v0 = Vector3::zeros();
        let z = Vector3::zeros();
        for half_periods in 1..=4 {
            let t = 0.5 * period * (half_periods as f64); // nT = kπ
            let r = two_impulse_rendezvous(n, t, r0, v0, z, z);
            assert!(
                matches!(r, Err(AstroError::NonConvergent(_))),
                "nT={half_periods}π must be singular, got {r:?}",
            );
        }
        // A non-degenerate time in between is solvable.
        assert!(two_impulse_rendezvous(n, 0.3 * period, r0, v0, z, z).is_ok());
    }

    #[test]
    fn rejects_non_physical_inputs() {
        assert!(mean_motion(0.0).is_err());
        assert!(mean_motion(-1.0).is_err());
        assert!(cw_state_transition(0.0, 10.0).is_err());
        assert!(cw_state_transition(1e-3, f64::NAN).is_err());
        let n = n_leo();
        let z = Vector3::zeros();
        assert!(two_impulse_rendezvous(n, 0.0, z, z, z, z).is_err()); // T = 0
        assert!(two_impulse_rendezvous(-1.0, 100.0, z, z, z, z).is_err()); // n < 0
        let bad = Vector3::new(f64::NAN, 0.0, 0.0);
        assert!(two_impulse_rendezvous(n, 100.0, bad, z, z, z).is_err());
    }

    #[test]
    fn rendezvous_plan_round_trips_through_json() {
        // Pin the serde derives on the nalgebra-backed plan (this relies on
        // nalgebra's `serde-serialize` feature, enabled at the workspace
        // level): the Vector3 burns must serialize and deserialize back. The
        // components are computed f64s, so JSON's shortest-decimal form can
        // land on an adjacent f64 on re-parse; assert the re-parsed values
        // match the originals within a tight relative tolerance.
        let n = n_leo();
        let period = TAU / n;
        let r0 = Vector3::new(-2_000.0, 5_000.0, 100.0);
        let v0 = Vector3::new(1.0, -2.0, 0.5);
        let z = Vector3::zeros();
        let plan = two_impulse_rendezvous(n, 0.25 * period, r0, v0, z, z).expect("solvable");

        let s = serde_json::to_string(&plan).expect("serialize plan");
        let back: RendezvousPlan = serde_json::from_str(&s).expect("deserialize plan");
        assert!((back.delta_v1 - plan.delta_v1).norm() <= plan.delta_v1.norm() * 1e-12);
        assert!((back.delta_v2 - plan.delta_v2).norm() <= plan.delta_v2.norm() * 1e-12);
        assert!((back.total_delta_v - plan.total_delta_v).abs() <= plan.total_delta_v * 1e-12);
        assert!((back.transfer_time - plan.transfer_time).abs() <= plan.transfer_time * 1e-12);
    }
}
