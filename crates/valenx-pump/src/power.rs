//! **Hydraulic and shaft power.**
//!
//! The useful power a pump adds to the fluid — the *hydraulic* (or
//! *water*) power — is the weight flow rate times the head:
//!
//! ```text
//!   P_hyd = rho * g * Q * H            [W]
//! ```
//!
//! with density `rho` in kg/m³, the standard gravity `g`
//! ([`crate::G`]), flow `Q` in m³/s and head `H` in metres of fluid.
//! Dividing by a pump efficiency `eta` (0, 1] gives the *shaft* power the
//! driver must deliver, since the hydraulic output is only that fraction
//! of the mechanical input:
//!
//! ```text
//!   P_shaft = P_hyd / eta             [W]
//! ```

use crate::error::{require_finite, require_non_negative, require_positive, PumpError};
use crate::G;

/// Hydraulic (water) power added to the fluid, in watts: `ρ·g·Q·H`.
///
/// `density_kg_m3` and `flow_m3s` must be finite and non-negative;
/// `head_m` must be finite (a negative head — a pump operated as a brake
/// against the flow — yields a negative power and is permitted).
///
/// # Errors
///
/// Returns [`PumpError::BadParameter`] if any argument is out of domain.
///
/// # Examples
///
/// ```
/// use valenx_pump::power::hydraulic_power_w;
///
/// // Water (1000 kg/m³) at 0.05 m³/s against 30 m of head.
/// let p = hydraulic_power_w(1000.0, 0.05, 30.0).unwrap();
/// assert!((p - 1000.0 * 9.806_65 * 0.05 * 30.0).abs() < 1e-9);
/// ```
pub fn hydraulic_power_w(density_kg_m3: f64, flow_m3s: f64, head_m: f64) -> Result<f64, PumpError> {
    let density_kg_m3 = require_non_negative("density_kg_m3", density_kg_m3)?;
    let flow_m3s = require_non_negative("flow_m3s", flow_m3s)?;
    let head_m = require_finite("head_m", head_m)?;
    Ok(density_kg_m3 * G * flow_m3s * head_m)
}

/// Shaft (brake) power the driver must supply, in watts: `P_hyd / η`.
///
/// `efficiency` is the pump efficiency `η`, which must be finite and lie
/// in the half-open interval `(0, 1]`. The other arguments follow
/// [`hydraulic_power_w`].
///
/// # Errors
///
/// Returns [`PumpError::BadParameter`] if `efficiency` is outside
/// `(0, 1]` or any other argument is out of domain.
///
/// # Examples
///
/// ```
/// use valenx_pump::power::shaft_power_w;
///
/// // The 14.71 kW hydraulic point above at 75% efficiency.
/// let p = shaft_power_w(1000.0, 0.05, 30.0, 0.75).unwrap();
/// let hyd = 1000.0 * 9.806_65 * 0.05 * 30.0;
/// assert!((p - hyd / 0.75).abs() < 1e-9);
/// ```
pub fn shaft_power_w(
    density_kg_m3: f64,
    flow_m3s: f64,
    head_m: f64,
    efficiency: f64,
) -> Result<f64, PumpError> {
    let efficiency = require_positive("efficiency", efficiency)?;
    if efficiency > 1.0 {
        return Err(PumpError::bad(
            "efficiency",
            format!("must be <= 1, got {efficiency}"),
        ));
    }
    let hyd = hydraulic_power_w(density_kg_m3, flow_m3s, head_m)?;
    Ok(hyd / efficiency)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn hydraulic_power_equals_rho_g_q_h() {
        // 1000 * 9.80665 * 0.05 * 30 = 14709.975 W.
        let p = hydraulic_power_w(1000.0, 0.05, 30.0).unwrap();
        assert!((p - 14_709.975).abs() < 1e-6);
    }

    #[test]
    fn power_scales_linearly_in_flow_and_head() {
        let base = hydraulic_power_w(998.0, 0.04, 25.0).unwrap();
        // Doubling flow doubles power.
        let twice_q = hydraulic_power_w(998.0, 0.08, 25.0).unwrap();
        assert!((twice_q / base - 2.0).abs() < EPS);
        // Doubling head doubles power.
        let twice_h = hydraulic_power_w(998.0, 0.04, 50.0).unwrap();
        assert!((twice_h / base - 2.0).abs() < EPS);
    }

    #[test]
    fn denser_fluid_needs_more_power() {
        // Mercury vs water at the same Q and H: power scales with density.
        let water = hydraulic_power_w(1000.0, 0.01, 10.0).unwrap();
        let mercury = hydraulic_power_w(13_534.0, 0.01, 10.0).unwrap();
        assert!((mercury / water - 13.534).abs() < 1e-6);
    }

    #[test]
    fn zero_flow_or_zero_head_is_zero_power() {
        assert!((hydraulic_power_w(1000.0, 0.0, 30.0).unwrap()).abs() < EPS);
        assert!((hydraulic_power_w(1000.0, 0.05, 0.0).unwrap()).abs() < EPS);
    }

    #[test]
    fn shaft_power_is_hydraulic_over_efficiency() {
        let hyd = hydraulic_power_w(1000.0, 0.05, 30.0).unwrap();
        let shaft = shaft_power_w(1000.0, 0.05, 30.0, 0.75).unwrap();
        assert!((shaft - hyd / 0.75).abs() < 1e-6);
        // Shaft power must exceed the useful hydraulic power.
        assert!(shaft > hyd);
    }

    #[test]
    fn unit_efficiency_gives_hydraulic_power() {
        let hyd = hydraulic_power_w(1000.0, 0.03, 18.0).unwrap();
        let shaft = shaft_power_w(1000.0, 0.03, 18.0, 1.0).unwrap();
        assert!((shaft - hyd).abs() < 1e-9);
    }

    #[test]
    fn rejects_efficiency_out_of_range() {
        assert_eq!(
            shaft_power_w(1000.0, 0.05, 30.0, 0.0).unwrap_err().code(),
            "pump.bad_parameter"
        );
        assert_eq!(
            shaft_power_w(1000.0, 0.05, 30.0, 1.5).unwrap_err().code(),
            "pump.bad_parameter"
        );
    }

    #[test]
    fn rejects_negative_density() {
        assert!(hydraulic_power_w(-1.0, 0.05, 30.0).is_err());
    }
}
