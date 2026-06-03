//! Lambert's problem: the two-body transfer connecting two position
//! vectors in a specified time of flight.
//!
//! Given start and end positions `r1`, `r2` (ECI, m) and a time of
//! flight `tof` (s), solve for the velocities `v1` (at `r1`) and `v2`
//! (at `r2`) of the unique single-revolution conic that links them.
//! This is the targeting primitive behind rendezvous, intercept and
//! interplanetary trajectory design — pair it with [`crate::orbit3d`]
//! to get the departure/arrival Δv against known orbits.
//!
//! Implementation: the universal-variable formulation (Bate–Mueller–
//! White / Vallado) with bisection on the universal variable `ψ`.
//! Single-revolution, prograde or retrograde; the (near-)180° degenerate
//! geometry (transfer plane undefined) returns `None`.

use nalgebra::Vector3;

use crate::constants::MU_EARTH;
use crate::orbit3d::StateVector;

/// Stumpff-style coefficients `c2(ψ)`, `c3(ψ)`.
fn stumpff(psi: f64) -> (f64, f64) {
    if psi > 1e-6 {
        let s = psi.sqrt();
        ((1.0 - s.cos()) / psi, (s - s.sin()) / psi.powf(1.5))
    } else if psi < -1e-6 {
        let s = (-psi).sqrt();
        ((1.0 - s.cosh()) / psi, (s.sinh() - s) / (-psi).powf(1.5))
    } else {
        (0.5, 1.0 / 6.0)
    }
}

/// Solve Lambert's problem about Earth (`μ = μ⊕`).
///
/// Returns `(v1, v2)` or `None` if the geometry is degenerate or the
/// iteration fails to converge.
pub fn lambert(
    r1: Vector3<f64>,
    r2: Vector3<f64>,
    tof: f64,
    prograde: bool,
) -> Option<(Vector3<f64>, Vector3<f64>)> {
    lambert_mu(r1, r2, tof, prograde, MU_EARTH)
}

/// Solve Lambert's problem for an arbitrary central-body `μ`.
pub fn lambert_mu(
    r1v: Vector3<f64>,
    r2v: Vector3<f64>,
    tof: f64,
    prograde: bool,
    mu: f64,
) -> Option<(Vector3<f64>, Vector3<f64>)> {
    if tof <= 0.0 {
        return None;
    }
    let r1 = r1v.norm();
    let r2 = r2v.norm();

    let cos_dnu = (r1v.dot(&r2v) / (r1 * r2)).clamp(-1.0, 1.0);
    let cross = r1v.cross(&r2v);

    // Transfer direction: for a prograde transfer, a negative z angular
    // momentum means the short way sweeps > 180°, so flip the A sign.
    let long_way = if prograde {
        cross.z < 0.0
    } else {
        cross.z > 0.0
    };
    let tm = if long_way { -1.0 } else { 1.0 };

    let a = tm * (r1 * r2 * (1.0 + cos_dnu)).sqrt();
    if a.abs() < 1e-9 {
        // ~180° transfer: the plane is undefined.
        return None;
    }

    let sqrt_mu = mu.sqrt();
    let mut psi = 0.0_f64;
    let (mut c2, mut c3) = (0.5_f64, 1.0 / 6.0_f64);
    let mut psi_up = 4.0 * std::f64::consts::PI * std::f64::consts::PI;
    let mut psi_low = -4.0 * std::f64::consts::PI;

    let mut y = 0.0_f64;
    let mut converged = false;
    for _ in 0..300 {
        y = r1 + r2 + a * (psi * c3 - 1.0) / c2.sqrt();
        // For a>0 the bisection must keep y >= 0; if it dips negative the
        // current ψ is too low, so raise the lower bracket.
        if a > 0.0 && y < 0.0 {
            psi_low = psi;
            psi = 0.5 * (psi_up + psi_low);
            let (nc2, nc3) = stumpff(psi);
            c2 = nc2;
            c3 = nc3;
            continue;
        }
        let chi = (y / c2).sqrt();
        let tof_calc = (chi.powi(3) * c3 + a * y.sqrt()) / sqrt_mu;

        if (tof_calc - tof).abs() < 1e-5 {
            converged = true;
            break;
        }
        if tof_calc <= tof {
            psi_low = psi;
        } else {
            psi_up = psi;
        }
        psi = 0.5 * (psi_up + psi_low);
        let (nc2, nc3) = stumpff(psi);
        c2 = nc2;
        c3 = nc3;
    }

    if !converged || !y.is_finite() || y < 0.0 {
        return None;
    }

    // Lagrange coefficients.
    let f = 1.0 - y / r1;
    let g = a * (y / mu).sqrt();
    let g_dot = 1.0 - y / r2;
    if g.abs() < 1e-12 {
        return None;
    }
    let v1 = (r2v - f * r1v) / g;
    let v2 = (g_dot * r2v - r1v) / g;
    Some((v1, v2))
}

/// Convenience: solve Lambert and return the departure / arrival states.
pub fn lambert_states(
    r1: Vector3<f64>,
    r2: Vector3<f64>,
    tof: f64,
    prograde: bool,
) -> Option<(StateVector, StateVector)> {
    let (v1, v2) = lambert(r1, r2, tof, prograde)?;
    Some((
        StateVector {
            position: r1,
            velocity: v1,
        },
        StateVector {
            position: r2,
            velocity: v2,
        },
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::R_EARTH;
    use crate::orbit3d::{propagate, StateVector};

    /// Build a true transfer by propagating a known state, then check
    /// Lambert recovers the departure velocity from the endpoints + TOF.
    fn round_trip(r1: Vector3<f64>, v1_true: Vector3<f64>, tof: f64) {
        // High-accuracy two-body propagation to the arrival point.
        let s0 = StateVector {
            position: r1,
            velocity: v1_true,
        };
        let dt = 0.5;
        let steps = (tof / dt).round() as u64;
        let s1 = propagate(&s0, dt, steps, false).expect("valid step count");

        let (v1, v2) = lambert(r1, s1.position, tof, true).expect("lambert solves");
        assert!(
            (v1 - v1_true).norm() < 5.0,
            "v1 recover off by {} m/s",
            (v1 - v1_true).norm()
        );
        // Arrival velocity should match the propagated arrival velocity.
        assert!(
            (v2 - s1.velocity).norm() < 5.0,
            "v2 recover off by {} m/s",
            (v2 - s1.velocity).norm()
        );
    }

    #[test]
    fn round_trip_equatorial_short_arc() {
        let r1 = Vector3::new(R_EARTH + 600_000.0, 0.0, 0.0);
        let v1 = Vector3::new(0.0, 8_200.0, 0.0); // prograde, elliptical
        round_trip(r1, v1, 1_200.0);
    }

    #[test]
    fn round_trip_inclined_arc() {
        let r1 = Vector3::new(R_EARTH + 500_000.0, 0.0, 0.0);
        let v1 = Vector3::new(100.0, 7_000.0, 2_500.0); // out-of-plane component
        round_trip(r1, v1, 1_500.0);
    }

    #[test]
    fn degenerate_180_degree_returns_none() {
        let r1 = Vector3::new(R_EARTH + 500_000.0, 0.0, 0.0);
        let r2 = Vector3::new(-(R_EARTH + 500_000.0), 0.0, 0.0);
        assert!(lambert(r1, r2, 3_000.0, true).is_none());
    }

    #[test]
    fn nonpositive_tof_returns_none() {
        let r1 = Vector3::new(R_EARTH + 500_000.0, 0.0, 0.0);
        let r2 = Vector3::new(0.0, R_EARTH + 500_000.0, 0.0);
        assert!(lambert(r1, r2, 0.0, true).is_none());
    }
}
