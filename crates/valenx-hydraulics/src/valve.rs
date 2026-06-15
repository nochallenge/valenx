//! Control-valve flow from the flow-coefficient relation.
//!
//! A valve's throughput is summarised by a single flow coefficient
//! `Cv`. For incompressible liquid flow the volumetric flow is
//!
//! ```text
//! Q = Cv * sqrt(dP / SG)
//! ```
//!
//! where `dP` is the pressure drop across the valve and `SG` is the
//! liquid's specific gravity (water = 1, dimensionless). The result
//! shares the unit system of `Cv` and `dP`: feed `Cv` and `dP` in a
//! consistent set of units and `Q` comes out in the matching flow
//! unit. The square-root dependence means flow scales with the
//! square root of the pressure drop — quadrupling `dP` doubles `Q`.

use crate::error::HydraulicsError;

/// Liquid volumetric flow through a valve of flow coefficient `cv`
/// under pressure drop `delta_p` for a fluid of specific gravity
/// `specific_gravity`, via `Q = Cv sqrt(dP / SG)`.
///
/// The returned flow is in the unit implied by the caller's choice of
/// `cv` and `delta_p` (see the module docs); the function does not
/// impose a particular unit system. A zero pressure drop yields zero
/// flow.
///
/// # Errors
///
/// - [`HydraulicsError::NonPositive`] if `cv` or `specific_gravity` is
///   not strictly positive (or not finite).
/// - [`HydraulicsError::Negative`] if `delta_p` is negative (or not
///   finite). Zero is allowed and gives zero flow.
pub fn valve_flow(cv: f64, delta_p: f64, specific_gravity: f64) -> Result<f64, HydraulicsError> {
    let cv = HydraulicsError::require_positive("cv", cv)?;
    let sg = HydraulicsError::require_positive("specific_gravity", specific_gravity)?;
    let dp = HydraulicsError::require_non_negative("delta_p", delta_p)?;
    Ok(cv * (dp / sg).sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tolerance for floating-point comparisons.
    const EPS: f64 = 1e-9;

    #[test]
    fn matches_closed_form() {
        // Cv = 12, dP = 4, SG = 1  ->  12 * sqrt(4) = 24.
        let q = valve_flow(12.0, 4.0, 1.0).unwrap();
        assert!((q - 24.0).abs() < EPS);
    }

    #[test]
    fn specific_gravity_scaling() {
        // SG = 4 divides under the root: Cv * sqrt(dP/4) halves the
        // SG = 1 result for the same dP.
        let q_water = valve_flow(10.0, 9.0, 1.0).unwrap();
        let q_heavy = valve_flow(10.0, 9.0, 4.0).unwrap();
        assert!((q_water - 30.0).abs() < EPS);
        assert!((q_heavy - 15.0).abs() < EPS);
        assert!((q_heavy - q_water / 2.0).abs() < EPS);
    }

    #[test]
    fn flow_scales_with_sqrt_of_pressure_drop() {
        // Quadrupling dP must exactly double Q.
        let q1 = valve_flow(7.0, 5.0, 1.3).unwrap();
        let q4 = valve_flow(7.0, 20.0, 1.3).unwrap();
        assert!((q4 - 2.0 * q1).abs() < 1e-9);
    }

    #[test]
    fn nine_fold_pressure_drop_triples_flow() {
        let q1 = valve_flow(3.0, 2.0, 0.9).unwrap();
        let q9 = valve_flow(3.0, 18.0, 0.9).unwrap();
        assert!((q9 - 3.0 * q1).abs() < 1e-9);
    }

    #[test]
    fn flow_is_linear_in_cv() {
        let q1 = valve_flow(5.0, 8.0, 1.0).unwrap();
        let q2 = valve_flow(10.0, 8.0, 1.0).unwrap();
        assert!((q2 - 2.0 * q1).abs() < EPS);
    }

    #[test]
    fn zero_pressure_drop_gives_zero_flow() {
        let q = valve_flow(15.0, 0.0, 1.0).unwrap();
        assert!(q.abs() < EPS);
    }

    #[test]
    fn rejects_bad_inputs() {
        assert_eq!(
            valve_flow(0.0, 4.0, 1.0).unwrap_err().code(),
            "hydraulics.non_positive"
        );
        assert_eq!(
            valve_flow(1.0, 4.0, 0.0).unwrap_err().code(),
            "hydraulics.non_positive"
        );
        assert_eq!(
            valve_flow(1.0, -4.0, 1.0).unwrap_err().code(),
            "hydraulics.negative"
        );
    }
}
