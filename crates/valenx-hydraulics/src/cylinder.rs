//! Linear hydraulic actuators (cylinders).
//!
//! A double-acting cylinder is described by its bore (piston)
//! diameter and rod diameter. Pressure acting on a face develops a
//! force `F = p A`; the bore side exposes the full piston area while
//! the rod side exposes the annulus left around the rod, so the two
//! strokes are asymmetric.
//!
//! All lengths are metres, areas square metres, pressures pascals,
//! forces newtons, velocities metres per second, and flows cubic
//! metres per second. The relations are unit-coherent in SI, so any
//! self-consistent unit system works if the caller keeps it
//! consistent.

use crate::error::HydraulicsError;
use serde::{Deserialize, Serialize};

/// Which way a [`Cylinder`] is being driven.
///
/// The variant selects the face that pressure acts on, and therefore
/// the area used in `F = p A` and `Q = A v`.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Stroke {
    /// Extension: fluid enters the bore (cap) side, pushing the rod
    /// out. The full piston area is pressurised.
    Extend,
    /// Retraction: fluid enters the rod (annulus) side, pulling the
    /// rod in. Only the annulus area around the rod is pressurised.
    Retract,
}

impl Stroke {
    /// Short UI / log label.
    pub fn label(self) -> &'static str {
        match self {
            Stroke::Extend => "extend",
            Stroke::Retract => "retract",
        }
    }
}

/// A double-acting hydraulic cylinder defined by bore and rod
/// diameters.
///
/// Construct with [`Cylinder::new`], which validates that both
/// diameters are positive and that the rod is strictly thinner than
/// the bore (otherwise the rod side would have no area).
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Cylinder {
    /// Bore (piston) diameter, metres.
    bore_diameter_m: f64,
    /// Rod diameter, metres.
    rod_diameter_m: f64,
}

impl Cylinder {
    /// Create a cylinder from its bore and rod diameters (metres).
    ///
    /// # Errors
    ///
    /// Returns [`HydraulicsError::NonPositive`] if either diameter is
    /// not strictly positive (or not finite), and
    /// [`HydraulicsError::Geometry`] if the rod is at least as thick
    /// as the bore.
    pub fn new(bore_diameter_m: f64, rod_diameter_m: f64) -> Result<Self, HydraulicsError> {
        let bore = HydraulicsError::require_positive("bore_diameter_m", bore_diameter_m)?;
        let rod = HydraulicsError::require_positive("rod_diameter_m", rod_diameter_m)?;
        if rod >= bore {
            return Err(HydraulicsError::geometry(format!(
                "rod diameter {rod} must be < bore diameter {bore}"
            )));
        }
        Ok(Cylinder {
            bore_diameter_m: bore,
            rod_diameter_m: rod,
        })
    }

    /// Bore (piston) diameter, metres.
    pub fn bore_diameter_m(&self) -> f64 {
        self.bore_diameter_m
    }

    /// Rod diameter, metres.
    pub fn rod_diameter_m(&self) -> f64 {
        self.rod_diameter_m
    }

    /// Full piston (bore-side) area `A_bore = pi D^2 / 4`, square
    /// metres.
    pub fn bore_area_m2(&self) -> f64 {
        circle_area(self.bore_diameter_m)
    }

    /// Rod-side annulus area `A_rod = A_bore - pi d^2 / 4`, square
    /// metres.
    ///
    /// This is the bore area minus the cross-section the rod occupies.
    /// It is always strictly smaller than [`Cylinder::bore_area_m2`].
    pub fn rod_area_m2(&self) -> f64 {
        self.bore_area_m2() - circle_area(self.rod_diameter_m)
    }

    /// Effective pressurised area for `stroke`, square metres.
    ///
    /// Extending uses the bore area; retracting uses the rod annulus
    /// area.
    pub fn effective_area_m2(&self, stroke: Stroke) -> f64 {
        match stroke {
            Stroke::Extend => self.bore_area_m2(),
            Stroke::Retract => self.rod_area_m2(),
        }
    }

    /// Force developed at gauge pressure `pressure_pa` for `stroke`,
    /// newtons, via `F = p A`.
    ///
    /// Because the extend area exceeds the retract area, at equal
    /// pressure [`Stroke::Extend`] always returns the larger force.
    ///
    /// # Errors
    ///
    /// [`HydraulicsError::NonPositive`] if `pressure_pa` is not
    /// strictly positive (or not finite).
    pub fn force_n(&self, pressure_pa: f64, stroke: Stroke) -> Result<f64, HydraulicsError> {
        let p = HydraulicsError::require_positive("pressure_pa", pressure_pa)?;
        Ok(p * self.effective_area_m2(stroke))
    }

    /// Volumetric flow required to drive `stroke` at piston speed
    /// `velocity_m_s`, cubic metres per second, via `Q = A v`.
    ///
    /// # Errors
    ///
    /// [`HydraulicsError::NonPositive`] if `velocity_m_s` is not
    /// strictly positive (or not finite).
    pub fn flow_m3_s(&self, velocity_m_s: f64, stroke: Stroke) -> Result<f64, HydraulicsError> {
        let v = HydraulicsError::require_positive("velocity_m_s", velocity_m_s)?;
        Ok(self.effective_area_m2(stroke) * v)
    }

    /// Piston speed produced by a commanded flow `flow_m3_s` on
    /// `stroke`, metres per second, via `v = Q / A`.
    ///
    /// This is the exact inverse of [`Cylinder::flow_m3_s`].
    ///
    /// # Errors
    ///
    /// [`HydraulicsError::NonPositive`] if `flow_m3_s` is not strictly
    /// positive (or not finite).
    pub fn speed_m_s(&self, flow_m3_s: f64, stroke: Stroke) -> Result<f64, HydraulicsError> {
        let q = HydraulicsError::require_positive("flow_m3_s", flow_m3_s)?;
        Ok(q / self.effective_area_m2(stroke))
    }
}

/// Area of a circle of `diameter`, `pi d^2 / 4`, in the squared unit
/// of `diameter`.
fn circle_area(diameter: f64) -> f64 {
    std::f64::consts::PI * diameter * diameter / 4.0
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tolerance for floating-point comparisons.
    const EPS: f64 = 1e-9;

    fn cyl() -> Cylinder {
        // 100 mm bore, 50 mm rod.
        Cylinder::new(0.100, 0.050).unwrap()
    }

    #[test]
    fn bore_area_matches_pi_d2_over_4() {
        let c = cyl();
        let expected = std::f64::consts::PI * 0.100 * 0.100 / 4.0;
        assert!((c.bore_area_m2() - expected).abs() < EPS);
    }

    #[test]
    fn rod_area_is_bore_minus_rod_circle() {
        let c = cyl();
        let bore = std::f64::consts::PI * 0.100 * 0.100 / 4.0;
        let rod = std::f64::consts::PI * 0.050 * 0.050 / 4.0;
        let expected = bore - rod;
        assert!((c.rod_area_m2() - expected).abs() < EPS);
        // With a 50 mm rod in a 100 mm bore the rod circle is exactly
        // 1/4 of the bore, so the annulus is 3/4 of the bore.
        assert!((c.rod_area_m2() - 0.75 * bore).abs() < EPS);
    }

    #[test]
    fn force_equals_pressure_times_area_bore_side() {
        let c = cyl();
        let p = 21.0e6; // 210 bar
        let f = c.force_n(p, Stroke::Extend).unwrap();
        assert!((f - p * c.bore_area_m2()).abs() < 1e-3);
    }

    #[test]
    fn force_equals_pressure_times_area_rod_side() {
        let c = cyl();
        let p = 21.0e6;
        let f = c.force_n(p, Stroke::Retract).unwrap();
        assert!((f - p * c.rod_area_m2()).abs() < 1e-3);
    }

    #[test]
    fn extend_force_exceeds_retract_force_at_equal_pressure() {
        let c = cyl();
        let p = 15.0e6;
        let extend = c.force_n(p, Stroke::Extend).unwrap();
        let retract = c.force_n(p, Stroke::Retract).unwrap();
        assert!(extend > retract);
        // Ratio of forces equals ratio of areas.
        let area_ratio = c.bore_area_m2() / c.rod_area_m2();
        assert!((extend / retract - area_ratio).abs() < EPS);
    }

    #[test]
    fn flow_equals_area_times_velocity() {
        let c = cyl();
        let v = 0.20; // m/s
        let q = c.flow_m3_s(v, Stroke::Extend).unwrap();
        assert!((q - c.bore_area_m2() * v).abs() < EPS);
    }

    #[test]
    fn flow_doubles_when_velocity_doubles() {
        let c = cyl();
        let q1 = c.flow_m3_s(0.10, Stroke::Extend).unwrap();
        let q2 = c.flow_m3_s(0.20, Stroke::Extend).unwrap();
        assert!((q2 - 2.0 * q1).abs() < EPS);
    }

    #[test]
    fn speed_is_inverse_of_flow() {
        let c = cyl();
        let v = 0.37;
        for stroke in [Stroke::Extend, Stroke::Retract] {
            let q = c.flow_m3_s(v, stroke).unwrap();
            let back = c.speed_m_s(q, stroke).unwrap();
            assert!((back - v).abs() < EPS);
        }
    }

    #[test]
    fn retract_flow_is_smaller_than_extend_flow_for_same_speed() {
        let c = cyl();
        let v = 0.25;
        let extend = c.flow_m3_s(v, Stroke::Extend).unwrap();
        let retract = c.flow_m3_s(v, Stroke::Retract).unwrap();
        assert!(retract < extend);
    }

    #[test]
    fn rejects_rod_not_thinner_than_bore() {
        let err = Cylinder::new(0.05, 0.05).unwrap_err();
        assert_eq!(err.code(), "hydraulics.geometry");
        let err = Cylinder::new(0.05, 0.06).unwrap_err();
        assert_eq!(err.code(), "hydraulics.geometry");
    }

    #[test]
    fn rejects_non_positive_diameters() {
        assert_eq!(
            Cylinder::new(0.0, 0.01).unwrap_err().code(),
            "hydraulics.non_positive"
        );
        assert_eq!(
            Cylinder::new(0.1, -0.01).unwrap_err().code(),
            "hydraulics.non_positive"
        );
    }

    #[test]
    fn force_and_flow_reject_non_positive_inputs() {
        let c = cyl();
        assert!(c.force_n(0.0, Stroke::Extend).is_err());
        assert!(c.force_n(-1.0, Stroke::Extend).is_err());
        assert!(c.flow_m3_s(0.0, Stroke::Extend).is_err());
        assert!(c.speed_m_s(-2.0, Stroke::Retract).is_err());
    }

    #[test]
    fn stroke_label_round_trips() {
        assert_eq!(Stroke::Extend.label(), "extend");
        assert_eq!(Stroke::Retract.label(), "retract");
    }
}
