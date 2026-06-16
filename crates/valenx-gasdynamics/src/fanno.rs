//! Fanno flow — adiabatic flow with friction in a constant-area duct.
//!
//! With no heat transfer and constant area, wall friction drives the Mach number
//! of a compressible flow **toward unity** (subsonic flow accelerates, supersonic
//! flow decelerates) until it chokes at `M = 1`. The Fanno relations express each
//! property as a ratio to its sonic (`*`) reference value, as a function of `M`
//! and the specific-heat ratio `gamma`:
//!
//! ```text
//! T/T*      = (γ+1) / (2 + (γ−1) M²)
//! V/V*      = M · sqrt[ (γ+1) / (2 + (γ−1) M²) ]
//! ρ/ρ*      = (V*/V) = (1/M) · sqrt[ (2 + (γ−1) M²) / (γ+1) ]
//! p/p*      = (1/M) · sqrt[ (γ+1) / (2 + (γ−1) M²) ]
//! p0/p0*    = (1/M) · [ (2 + (γ−1) M²) / (γ+1) ] ^ ((γ+1)/(2(γ−1)))
//! 4f·L*/D   = (1 − M²)/(γ M²) + (γ+1)/(2γ) · ln[ (γ+1) M² / (2 + (γ−1) M²) ]
//! ```
//!
//! `4f·L*/D` is the Fanning-friction duct parameter needed to drive the flow from
//! `M` to the sonic point; it is `0` at `M = 1` and positive on both the subsonic
//! and supersonic branches.
//!
//! Reference: Anderson, *Modern Compressible Flow*; Shapiro; NACA Report 1135.
//! Same calorically-perfect-gas scope and caveats as the rest of the crate.

use crate::error::{check_gamma, check_mach_pos, Result};

/// `2 + (γ−1) M²`, the recurring grouping in the Fanno relations.
fn denom(m: f64, gamma: f64) -> f64 {
    2.0 + (gamma - 1.0) * m * m
}

/// Static-temperature ratio `T/T*`.
///
/// # Errors
/// [`crate::GasError`] for non-finite/`<= 0` Mach or invalid `gamma`.
pub fn temperature_ratio(m: f64, gamma: f64) -> Result<f64> {
    let gamma = check_gamma(gamma)?;
    let m = check_mach_pos(m, "fanno")?;
    Ok((gamma + 1.0) / denom(m, gamma))
}

/// Velocity ratio `V/V*` (equals `ρ*/ρ` by mass conservation).
///
/// # Errors
/// As [`temperature_ratio`].
pub fn velocity_ratio(m: f64, gamma: f64) -> Result<f64> {
    let gamma = check_gamma(gamma)?;
    let m = check_mach_pos(m, "fanno")?;
    Ok(m * ((gamma + 1.0) / denom(m, gamma)).sqrt())
}

/// Static-density ratio `ρ/ρ*` (`= V*/V`).
///
/// # Errors
/// As [`temperature_ratio`].
pub fn density_ratio(m: f64, gamma: f64) -> Result<f64> {
    Ok(1.0 / velocity_ratio(m, gamma)?)
}

/// Static-pressure ratio `p/p*`.
///
/// # Errors
/// As [`temperature_ratio`].
pub fn pressure_ratio(m: f64, gamma: f64) -> Result<f64> {
    let gamma = check_gamma(gamma)?;
    let m = check_mach_pos(m, "fanno")?;
    Ok((1.0 / m) * ((gamma + 1.0) / denom(m, gamma)).sqrt())
}

/// Stagnation-pressure ratio `p0/p0*` (`>= 1`; the friction entropy-loss
/// signature).
///
/// # Errors
/// As [`temperature_ratio`].
pub fn stagnation_pressure_ratio(m: f64, gamma: f64) -> Result<f64> {
    let gamma = check_gamma(gamma)?;
    let m = check_mach_pos(m, "fanno")?;
    let exponent = (gamma + 1.0) / (2.0 * (gamma - 1.0));
    Ok((1.0 / m) * (denom(m, gamma) / (gamma + 1.0)).powf(exponent))
}

/// The Fanning-friction duct parameter `4f·L*/D` to drive the flow from `M` to
/// the sonic point — `0` at `M = 1`, positive otherwise.
///
/// # Errors
/// As [`temperature_ratio`].
pub fn friction_length(m: f64, gamma: f64) -> Result<f64> {
    let gamma = check_gamma(gamma)?;
    let m = check_mach_pos(m, "fanno")?;
    let m2 = m * m;
    let term1 = (1.0 - m2) / (gamma * m2);
    let term2 = (gamma + 1.0) / (2.0 * gamma) * ((gamma + 1.0) * m2 / denom(m, gamma)).ln();
    Ok(term1 + term2)
}

/// The full Fanno state at Mach `m`: every sonic-referenced ratio plus the
/// friction length.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FannoState {
    /// `T/T*`.
    pub temperature_ratio: f64,
    /// `p/p*`.
    pub pressure_ratio: f64,
    /// `ρ/ρ*`.
    pub density_ratio: f64,
    /// `V/V*`.
    pub velocity_ratio: f64,
    /// `p0/p0*`.
    pub stagnation_pressure_ratio: f64,
    /// `4f·L*/D`.
    pub friction_length: f64,
}

/// Bundle all Fanno relations at Mach `m`.
///
/// # Errors
/// As [`temperature_ratio`].
pub fn fanno_state(m: f64, gamma: f64) -> Result<FannoState> {
    Ok(FannoState {
        temperature_ratio: temperature_ratio(m, gamma)?,
        pressure_ratio: pressure_ratio(m, gamma)?,
        density_ratio: density_ratio(m, gamma)?,
        velocity_ratio: velocity_ratio(m, gamma)?,
        stagnation_pressure_ratio: stagnation_pressure_ratio(m, gamma)?,
        friction_length: friction_length(m, gamma)?,
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
        let s = fanno_state(2.0, G).unwrap();
        assert!(
            close(s.temperature_ratio, 0.6667, 1e-3),
            "T/T* {}",
            s.temperature_ratio
        );
        assert!(
            close(s.pressure_ratio, 0.4082, 1e-3),
            "p/p* {}",
            s.pressure_ratio
        );
        assert!(
            close(s.stagnation_pressure_ratio, 1.6875, 1e-3),
            "p0/p0* {}",
            s.stagnation_pressure_ratio
        );
        assert!(
            close(s.friction_length, 0.3050, 1e-3),
            "4fL*/D {}",
            s.friction_length
        );
        assert!(
            close(s.velocity_ratio, 1.6330, 1e-3),
            "V/V* {}",
            s.velocity_ratio
        );
        // ρ/ρ* = V*/V.
        assert!(close(s.density_ratio, 1.0 / s.velocity_ratio, 1e-12));
    }

    #[test]
    fn subsonic_m05_matches_textbook() {
        let s = fanno_state(0.5, G).unwrap();
        assert!(
            close(s.temperature_ratio, 1.1429, 1e-3),
            "T/T* {}",
            s.temperature_ratio
        );
        assert!(
            close(s.pressure_ratio, 2.1381, 1e-3),
            "p/p* {}",
            s.pressure_ratio
        );
        assert!(
            close(s.stagnation_pressure_ratio, 1.3398, 1e-3),
            "p0/p0* {}",
            s.stagnation_pressure_ratio
        );
        assert!(
            close(s.friction_length, 1.0691, 1e-3),
            "4fL*/D {}",
            s.friction_length
        );
    }

    #[test]
    fn sonic_point_is_the_reference() {
        // At M = 1 every ratio is 1 and the friction length is 0.
        let s = fanno_state(1.0, G).unwrap();
        assert!(close(s.temperature_ratio, 1.0, 1e-12));
        assert!(close(s.pressure_ratio, 1.0, 1e-12));
        assert!(close(s.density_ratio, 1.0, 1e-12));
        assert!(close(s.velocity_ratio, 1.0, 1e-12));
        assert!(close(s.stagnation_pressure_ratio, 1.0, 1e-12));
        assert!(close(s.friction_length, 0.0, 1e-12));
    }

    #[test]
    fn friction_length_is_positive_on_both_branches_and_vanishes_at_sonic() {
        // Friction drives both subsonic and supersonic flow toward M = 1.
        assert!(friction_length(0.3, G).unwrap() > 0.0);
        assert!(friction_length(3.0, G).unwrap() > 0.0);
        // Monotone toward 0 approaching M = 1 from below.
        assert!(friction_length(0.5, G).unwrap() > friction_length(0.9, G).unwrap());
        assert!(friction_length(0.9, G).unwrap() > friction_length(0.99, G).unwrap());
    }

    #[test]
    fn rejects_bad_inputs() {
        assert!(temperature_ratio(0.0, G).is_err());
        assert!(temperature_ratio(-1.0, G).is_err());
        assert!(friction_length(2.0, 1.0).is_err()); // gamma must be > 1
    }
}
