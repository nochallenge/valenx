//! Darcy-Weisbach pressure-drop calculation.
//!
//! Standard form:
//!
//! ```text
//! ΔP = f * (L / D) * (ρ * v² / 2)
//! ```
//!
//! with `f` Darcy friction factor (dimensionless), `L` length, `D`
//! hydraulic diameter, `ρ` fluid density (kg/m³), `v` velocity (m/s).
//! Returns `ΔP` in Pa.

/// Darcy-Weisbach loss across a length of duct.
///
/// All inputs SI: `d` and `length` in m, `velocity` in m/s,
/// `friction_factor` dimensionless. Air at 20°C has ρ ≈ 1.204 kg/m³.
pub fn darcy_weisbach(d: f64, length: f64, velocity: f64, friction_factor: f64) -> f64 {
    if d <= 0.0 {
        return f64::INFINITY;
    }
    const RHO_AIR_20C: f64 = 1.204;
    friction_factor * (length / d) * 0.5 * RHO_AIR_20C * velocity * velocity
}

/// Same as [`darcy_weisbach`] with a caller-supplied fluid density.
pub fn darcy_weisbach_rho(d: f64, length: f64, velocity: f64, f: f64, rho: f64) -> f64 {
    if d <= 0.0 {
        return f64::INFINITY;
    }
    f * (length / d) * 0.5 * rho * velocity * velocity
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pressure_drop_is_positive_for_positive_inputs() {
        let dp = darcy_weisbach(0.3, 10.0, 5.0, 0.02);
        assert!(dp > 0.0);
    }

    #[test]
    fn pressure_drop_returns_inf_on_zero_diameter() {
        let dp = darcy_weisbach(0.0, 10.0, 5.0, 0.02);
        assert!(dp.is_infinite());
    }

    #[test]
    fn caller_density_lower_gives_lower_dp() {
        let dp_dense = darcy_weisbach_rho(0.3, 10.0, 5.0, 0.02, 1.5);
        let dp_thin = darcy_weisbach_rho(0.3, 10.0, 5.0, 0.02, 0.5);
        assert!(dp_thin < dp_dense);
    }
}
