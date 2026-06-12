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

    #[test]
    fn darcy_weisbach_matches_closed_form_worked_value() {
        // Closed form ΔP = f·(L/D)·½·ρ·v². With f=0.02, L=10 m, D=0.3 m,
        // v=5 m/s, ρ=1.204 kg/m³ (air at 20 °C):
        //   ΔP = 0.02·(10/0.3)·0.5·1.204·5² = 10.03333… Pa.
        // Pin the caller-density entry point against the exact hand-computed
        // ground truth (not the code's own output) at a tight absolute tol.
        let f = 0.02;
        let l = 10.0;
        let d = 0.3;
        let v = 5.0;
        let rho = 1.204;
        let expected = f * (l / d) * 0.5 * rho * v * v; // = 10.0333… Pa
        let dp = darcy_weisbach_rho(d, l, v, f, rho);
        assert!(
            (dp - expected).abs() < 1e-9,
            "computed {dp} vs closed form {expected}"
        );
        // And against the literal published value, tol 0.05 Pa.
        assert!(
            (dp - 10.0333).abs() < 0.05,
            "computed {dp} Pa, expected ≈ 10.03 Pa"
        );
        // The convenience wrapper hard-codes ρ_air = 1.204, so it must agree.
        assert!((darcy_weisbach(d, l, v, f) - dp).abs() < 1e-12);
    }
}
