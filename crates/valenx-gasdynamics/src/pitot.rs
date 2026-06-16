//! Rayleigh supersonic Pitot-tube relation.
//!
//! When a Pitot probe faces a supersonic stream a detached bow shock stands
//! ahead of it; the probe senses the stagnation pressure `p02` *behind* the
//! normal portion of that shock, not the free-stream stagnation pressure `p01`
//! (which is unrecoverable across the shock). Referenced to the upstream
//! *static* pressure `p1`, the reading follows the **Rayleigh Pitot formula**
//! (NACA Report 1135, eq. 100):
//!
//! ```text
//! p02/p1 = [ (γ+1)² M1² / (4γ M1² − 2(γ−1)) ] ^ (γ/(γ−1))
//!          · [ (1 − γ + 2γ M1²) / (γ+1) ]
//! ```
//!
//! valid for `M1 ≥ 1`. It is the supersonic analogue of the subsonic Pitot
//! reading `p0/p` (the isentropic stagnation ratio, recovered when no shock
//! forms); the two agree at `M1 = 1`, where the shock is infinitely weak.
//! [`mach_from_pitot_ratio`] inverts the relation — the practical direction,
//! recovering the free-stream Mach number from a measured `p02/p1`.
//!
//! Equivalently `p02/p1 = (p02/p01) · (p01/p1)`, the product of the
//! normal-shock stagnation-pressure ratio and the upstream isentropic
//! stagnation ratio — the cross-check pinned by the tests.
//!
//! Reference: Anderson, *Modern Compressible Flow*; NACA Report 1135. Same
//! perfect-gas scope and caveats as the rest of the crate.

use crate::error::{check_gamma, GasError, Result};

/// The Rayleigh Pitot ratio `p02/p1` — stagnation pressure behind the bow
/// shock over upstream static pressure — for a supersonic Mach `m1` (`>= 1`).
///
/// # Errors
///
/// [`GasError::BadGamma`] for invalid `gamma`; [`GasError::BadMach`] if `m1` is
/// not finite or `< 1` (the relation is defined only for supersonic flow).
pub fn rayleigh_pitot_ratio(m1: f64, gamma: f64) -> Result<f64> {
    let gamma = check_gamma(gamma)?;
    if !m1.is_finite() || m1 < 1.0 {
        return Err(GasError::bad_mach(
            m1,
            "rayleigh_pitot",
            "must be finite and >= 1 (supersonic)",
        ));
    }
    let m2 = m1 * m1;
    let term1 = ((gamma + 1.0) * (gamma + 1.0) * m2 / (4.0 * gamma * m2 - 2.0 * (gamma - 1.0)))
        .powf(gamma / (gamma - 1.0));
    let term2 = (1.0 - gamma + 2.0 * gamma * m2) / (gamma + 1.0);
    Ok(term1 * term2)
}

/// Invert the Rayleigh Pitot relation: recover the supersonic free-stream Mach
/// number from a measured stagnation-behind-shock over upstream-static ratio
/// `p02_over_p1`. The ratio increases monotonically with `M1` for `M1 >= 1`,
/// so the root is unique; it is found by bracketing then bisection.
///
/// # Errors
///
/// [`GasError::BadGamma`] for invalid `gamma`; [`GasError::BadMach`] if the
/// ratio is not finite or is below the `M1 = 1` value (no supersonic solution).
pub fn mach_from_pitot_ratio(p02_over_p1: f64, gamma: f64) -> Result<f64> {
    let gamma = check_gamma(gamma)?;
    let min_ratio = rayleigh_pitot_ratio(1.0, gamma)?;
    if !p02_over_p1.is_finite() || p02_over_p1 < min_ratio {
        return Err(GasError::bad_mach(
            p02_over_p1,
            "mach_from_pitot_ratio",
            "p02/p1 must be finite and >= its M=1 value (supersonic regime)",
        ));
    }
    // Grow the upper bound until it straddles the target, then bisect the
    // monotone branch.
    let mut lo = 1.0;
    let mut hi = 2.0;
    while rayleigh_pitot_ratio(hi, gamma)? < p02_over_p1 {
        hi *= 2.0;
        if hi > 1.0e6 {
            break;
        }
    }
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if rayleigh_pitot_ratio(mid, gamma)? < p02_over_p1 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Ok(0.5 * (lo + hi))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::isentropic::pressure_ratio as isentropic_pressure_ratio;
    use crate::normal_shock::normal_shock;

    const G: f64 = 1.4;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn pitot_ratio_matches_naca_1135_at_m2() {
        // NACA 1135 supersonic-Pitot column, gamma = 1.4, M = 2: p02/p1 = 5.6404.
        let r = rayleigh_pitot_ratio(2.0, G).unwrap();
        assert!(close(r, 5.6404, 1e-3), "p02/p1 = {r}");
    }

    #[test]
    fn agrees_with_isentropic_at_sonic() {
        // At M = 1 the bow shock is infinitely weak, so the supersonic Pitot
        // reading collapses onto the subsonic isentropic stagnation ratio
        // p0/p = 1.2^3.5 = 1.892929 (gamma = 1.4).
        let r = rayleigh_pitot_ratio(1.0, G).unwrap();
        assert!(close(r, 1.2_f64.powf(3.5), 1e-12), "p02/p1(1) = {r}");
        assert!(close(r, isentropic_pressure_ratio(1.0, G).unwrap(), 1e-12));
    }

    #[test]
    fn equals_shock_stagnation_times_upstream_isentropic() {
        // Independent cross-check across modules: p02/p1 = (p02/p01)·(p01/p1),
        // the normal-shock stagnation-pressure ratio times the upstream
        // isentropic stagnation ratio. Holds at every supersonic Mach.
        for &m in &[1.5, 2.0, 3.0, 5.0] {
            let pitot = rayleigh_pitot_ratio(m, G).unwrap();
            let composed = normal_shock(m, G).unwrap().stagnation_pressure_ratio
                * isentropic_pressure_ratio(m, G).unwrap();
            assert!(close(pitot, composed, 1e-9), "M={m}: {pitot} vs {composed}");
        }
    }

    #[test]
    fn ratio_increases_monotonically_with_mach() {
        assert!(rayleigh_pitot_ratio(1.5, G).unwrap() < rayleigh_pitot_ratio(2.0, G).unwrap());
        assert!(rayleigh_pitot_ratio(2.0, G).unwrap() < rayleigh_pitot_ratio(4.0, G).unwrap());
    }

    #[test]
    fn inverse_round_trips() {
        for &m in &[1.2, 2.0, 2.5, 4.0] {
            let ratio = rayleigh_pitot_ratio(m, G).unwrap();
            let recovered = mach_from_pitot_ratio(ratio, G).unwrap();
            assert!(close(recovered, m, 1e-6), "M={m} -> {recovered}");
        }
    }

    #[test]
    fn rejects_bad_inputs() {
        // Subsonic Mach has no bow shock.
        assert!(rayleigh_pitot_ratio(0.8, G).is_err());
        assert!(rayleigh_pitot_ratio(f64::NAN, G).is_err());
        assert!(rayleigh_pitot_ratio(2.0, 1.0).is_err()); // gamma must be > 1
                                                          // A ratio below the M=1 value has no supersonic root.
        let min = rayleigh_pitot_ratio(1.0, G).unwrap();
        assert!(mach_from_pitot_ratio(min - 0.1, G).is_err());
    }
}
