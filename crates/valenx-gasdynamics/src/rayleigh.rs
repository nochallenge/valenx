//! Rayleigh flow — frictionless flow with heat addition in a constant-area duct.
//!
//! With no wall friction and constant area, adding heat to a compressible flow
//! drives its Mach number **toward unity** (subsonic flow accelerates, supersonic
//! flow decelerates) until it thermally chokes at `M = 1`, where the stagnation
//! temperature reaches its maximum and no further heat can be added without
//! shifting the whole operating point. The Rayleigh relations express each
//! property as a ratio to its sonic (`*`) reference value, as a function of `M`
//! and the specific-heat ratio `gamma`:
//!
//! ```text
//! p/p*      = (1 + γ) / (1 + γ M²)
//! T/T*      = M² (1 + γ)² / (1 + γ M²)²
//! ρ/ρ*      = (1 + γ M²) / ((1 + γ) M²)            (= V*/V)
//! V/V*      = (1 + γ) M² / (1 + γ M²)
//! T0/T0*    = (γ+1) M² (2 + (γ−1) M²) / (1 + γ M²)²
//! p0/p0*    = (1 + γ)/(1 + γ M²) · [ (2 + (γ−1) M²)/(γ + 1) ] ^ (γ/(γ−1))
//! ```
//!
//! The stagnation-temperature ratio `T0/T0*` is the **heat-addition parameter**:
//! it is `1` at `M = 1` (its maximum — the thermal-choking limit), and less than
//! one on both the subsonic and supersonic branches. Heating raises `T0/T0*`
//! toward unity from either side; the stagnation-pressure ratio `p0/p0*`
//! correspondingly bottoms out at `1` at `M = 1` (the heat-addition
//! stagnation-pressure-loss signature) and exceeds one elsewhere.
//!
//! Reference: Anderson, *Modern Compressible Flow*; Shapiro; NACA Report 1135.
//! Same calorically-perfect-gas scope and caveats as the rest of the crate.

use crate::error::{check_gamma, check_mach_pos, Result};

/// `1 + γ M²`, the recurring grouping in the Rayleigh relations.
fn denom(m: f64, gamma: f64) -> f64 {
    1.0 + gamma * m * m
}

/// Static-pressure ratio `p/p*`.
///
/// # Errors
/// [`crate::GasError`] for non-finite/`<= 0` Mach or invalid `gamma`.
pub fn pressure_ratio(m: f64, gamma: f64) -> Result<f64> {
    let gamma = check_gamma(gamma)?;
    let m = check_mach_pos(m, "rayleigh")?;
    Ok((1.0 + gamma) / denom(m, gamma))
}

/// Static-temperature ratio `T/T*`.
///
/// # Errors
/// As [`pressure_ratio`].
pub fn temperature_ratio(m: f64, gamma: f64) -> Result<f64> {
    let gamma = check_gamma(gamma)?;
    let m = check_mach_pos(m, "rayleigh")?;
    let r = m * (1.0 + gamma) / denom(m, gamma);
    Ok(r * r)
}

/// Static-density ratio `ρ/ρ*` (`= V*/V`).
///
/// # Errors
/// As [`pressure_ratio`].
pub fn density_ratio(m: f64, gamma: f64) -> Result<f64> {
    let gamma = check_gamma(gamma)?;
    let m = check_mach_pos(m, "rayleigh")?;
    Ok(denom(m, gamma) / ((1.0 + gamma) * m * m))
}

/// Velocity ratio `V/V*` (`= ρ*/ρ`).
///
/// # Errors
/// As [`pressure_ratio`].
pub fn velocity_ratio(m: f64, gamma: f64) -> Result<f64> {
    let gamma = check_gamma(gamma)?;
    let m = check_mach_pos(m, "rayleigh")?;
    Ok((1.0 + gamma) * m * m / denom(m, gamma))
}

/// Stagnation-temperature ratio `T0/T0*` — the heat-addition parameter.
///
/// Equals `1` at `M = 1` (its maximum; the thermal-choking limit) and is less
/// than one on both the subsonic and supersonic branches.
///
/// # Errors
/// As [`pressure_ratio`].
pub fn stagnation_temperature_ratio(m: f64, gamma: f64) -> Result<f64> {
    let gamma = check_gamma(gamma)?;
    let m = check_mach_pos(m, "rayleigh")?;
    let d = denom(m, gamma);
    Ok((gamma + 1.0) * m * m * (2.0 + (gamma - 1.0) * m * m) / (d * d))
}

/// Stagnation-pressure ratio `p0/p0*` (`>= 1`; the heat-addition
/// stagnation-pressure-loss signature, minimised at `1` at `M = 1`).
///
/// # Errors
/// As [`pressure_ratio`].
pub fn stagnation_pressure_ratio(m: f64, gamma: f64) -> Result<f64> {
    let gamma = check_gamma(gamma)?;
    let m = check_mach_pos(m, "rayleigh")?;
    let exponent = gamma / (gamma - 1.0);
    Ok((1.0 + gamma) / denom(m, gamma)
        * ((2.0 + (gamma - 1.0) * m * m) / (gamma + 1.0)).powf(exponent))
}

/// The full Rayleigh state at Mach `m`: every sonic-referenced ratio.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RayleighState {
    /// `T/T*`.
    pub temperature_ratio: f64,
    /// `p/p*`.
    pub pressure_ratio: f64,
    /// `ρ/ρ*`.
    pub density_ratio: f64,
    /// `V/V*`.
    pub velocity_ratio: f64,
    /// `T0/T0*` (the heat-addition parameter).
    pub stagnation_temperature_ratio: f64,
    /// `p0/p0*`.
    pub stagnation_pressure_ratio: f64,
}

/// Bundle all Rayleigh relations at Mach `m`.
///
/// # Errors
/// As [`pressure_ratio`].
pub fn rayleigh_state(m: f64, gamma: f64) -> Result<RayleighState> {
    Ok(RayleighState {
        temperature_ratio: temperature_ratio(m, gamma)?,
        pressure_ratio: pressure_ratio(m, gamma)?,
        density_ratio: density_ratio(m, gamma)?,
        velocity_ratio: velocity_ratio(m, gamma)?,
        stagnation_temperature_ratio: stagnation_temperature_ratio(m, gamma)?,
        stagnation_pressure_ratio: stagnation_pressure_ratio(m, gamma)?,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const G: f64 = 1.4;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn supersonic_m2_matches_naca_1135() {
        // NACA 1135 / Anderson App. B Rayleigh table, gamma = 1.4, M = 2:
        let s = rayleigh_state(2.0, G).unwrap();
        assert!(
            close(s.pressure_ratio, 0.3636, 1e-3),
            "p/p* {}",
            s.pressure_ratio
        );
        assert!(
            close(s.temperature_ratio, 0.5289, 1e-3),
            "T/T* {}",
            s.temperature_ratio
        );
        assert!(
            close(s.stagnation_temperature_ratio, 0.7934, 1e-3),
            "T0/T0* {}",
            s.stagnation_temperature_ratio
        );
        assert!(
            close(s.stagnation_pressure_ratio, 1.5031, 1e-3),
            "p0/p0* {}",
            s.stagnation_pressure_ratio
        );
        // V/V* = 1.4545; rho/rho* = V*/V is its reciprocal.
        assert!(
            close(s.velocity_ratio, 1.4545, 1e-3),
            "V/V* {}",
            s.velocity_ratio
        );
        assert!(close(s.density_ratio, 1.0 / s.velocity_ratio, 1e-12));
    }

    #[test]
    fn subsonic_m05_matches_closed_form() {
        // gamma = 1.4, M = 0.5 closed forms:
        //   p/p*   = 2.4 / 1.35              = 1.77778
        //   T/T*   = (0.5*2.4/1.35)^2        = 0.79012
        //   T0/T0* = 2.4*0.25*2.1 / 1.35^2   = 0.69136
        let s = rayleigh_state(0.5, G).unwrap();
        assert!(
            close(s.pressure_ratio, 2.4 / 1.35, 1e-9),
            "p/p* {}",
            s.pressure_ratio
        );
        let t = (0.5 * 2.4 / 1.35_f64).powi(2);
        assert!(
            close(s.temperature_ratio, t, 1e-9),
            "T/T* {}",
            s.temperature_ratio
        );
        assert!(
            close(s.stagnation_temperature_ratio, 0.69136, 1e-4),
            "T0/T0* {}",
            s.stagnation_temperature_ratio
        );
    }

    #[test]
    fn sonic_point_is_the_reference() {
        // At M = 1 every ratio is exactly 1.
        let s = rayleigh_state(1.0, G).unwrap();
        assert!(close(s.temperature_ratio, 1.0, 1e-12));
        assert!(close(s.pressure_ratio, 1.0, 1e-12));
        assert!(close(s.density_ratio, 1.0, 1e-12));
        assert!(close(s.velocity_ratio, 1.0, 1e-12));
        assert!(close(s.stagnation_temperature_ratio, 1.0, 1e-12));
        assert!(close(s.stagnation_pressure_ratio, 1.0, 1e-12));
    }

    #[test]
    fn stagnation_temperature_peaks_at_sonic() {
        // T0/T0* is the heat-addition parameter: its maximum (1) is at M = 1,
        // approached from below on both the subsonic and supersonic branches —
        // this is the thermal-choking limit.
        for &m in &[0.3, 0.5, 0.8, 1.2, 2.0, 3.0] {
            let r = stagnation_temperature_ratio(m, G).unwrap();
            assert!(r < 1.0, "T0/T0*({m}) = {r} should be < 1 off-sonic");
        }
        // Monotone increase toward the sonic maximum from each side.
        assert!(
            stagnation_temperature_ratio(0.5, G).unwrap()
                < stagnation_temperature_ratio(0.9, G).unwrap()
        );
        assert!(
            stagnation_temperature_ratio(3.0, G).unwrap()
                < stagnation_temperature_ratio(1.5, G).unwrap()
        );
    }

    #[test]
    fn stagnation_pressure_bottoms_at_sonic() {
        // p0/p0* >= 1 with the minimum (1) at M = 1.
        assert!(stagnation_pressure_ratio(0.5, G).unwrap() > 1.0);
        assert!(stagnation_pressure_ratio(2.0, G).unwrap() > 1.0);
        assert!(
            stagnation_pressure_ratio(0.9, G).unwrap() < stagnation_pressure_ratio(0.5, G).unwrap()
        );
    }

    #[test]
    fn rejects_bad_inputs() {
        assert!(pressure_ratio(0.0, G).is_err());
        assert!(pressure_ratio(-1.0, G).is_err());
        assert!(temperature_ratio(2.0, 1.0).is_err()); // gamma must be > 1
    }
}
