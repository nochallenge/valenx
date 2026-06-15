//! Belt geometry, transported load, and throughput.
//!
//! # Model
//!
//! A belt conveyor carries bulk material as a moving prism: a (roughly)
//! constant material cross-section `A` (m^2) translated along the belt
//! at speed `v` (m/s). The volume swept per unit time is therefore
//! `Q_v = A * v` (m^3/s) and, for bulk density `rho` (kg/m^3), the mass
//! flow rate is
//!
//! ```text
//! mdot = rho * A * v          [kg/s]
//! ```
//!
//! This is the textbook continuity relation for a moving stream of
//! constant section. The load area `A` itself is a property of the belt
//! width, the surcharge/fill geometry and the material; this crate lets
//! you either pass `A` directly or estimate it from a flat belt with a
//! prescribed fractional fill.
//!
//! # Honest scope
//!
//! The fill estimate here is a deliberately simple rectangular-fill
//! approximation (`A = fill_fraction * width^2`), useful for teaching
//! and order-of-magnitude work. Real conveyor sizing uses troughed
//! idler angles, material surcharge angles and edge-distance rules from
//! CEMA / ISO 5048; those refinements are out of scope.

use crate::error::ConveyorError;
use serde::{Deserialize, Serialize};

/// Standard gravitational acceleration in m/s^2 (CODATA conventional
/// value), shared by the throughput and power models.
pub const G: f64 = 9.80665;

/// A belt conveyor's transported material stream, expressed as a bulk
/// density, a material cross-sectional load area and a belt speed.
///
/// Construct with [`Belt::new`] (validates every field) or, when only
/// the belt width and a fractional fill are known, with
/// [`Belt::from_flat_fill`].
///
/// All quantities are SI: density in kg/m^3, area in m^2, speed in m/s.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Belt {
    /// Bulk density of the conveyed material, kg/m^3.
    pub density: f64,
    /// Material cross-sectional load area carried by the belt, m^2.
    pub load_area: f64,
    /// Belt travel speed, m/s.
    pub speed: f64,
}

impl Belt {
    /// Build a belt from an explicit load area.
    ///
    /// # Errors
    ///
    /// Returns [`ConveyorError`] if any of `density`, `load_area` or
    /// `speed` is non-finite or not strictly positive.
    pub fn new(density: f64, load_area: f64, speed: f64) -> Result<Self, ConveyorError> {
        let density = ConveyorError::require_positive("density", density)?;
        let load_area = ConveyorError::require_positive("load_area", load_area)?;
        let speed = ConveyorError::require_positive("speed", speed)?;
        Ok(Self {
            density,
            load_area,
            speed,
        })
    }

    /// Build a belt whose load area is estimated from a flat belt of the
    /// given `width` (m) filled to `fill_fraction` of a square `width^2`
    /// reference section.
    ///
    /// The load area is `A = fill_fraction * width^2`. `fill_fraction`
    /// must lie in `(0, 1]`.
    ///
    /// # Errors
    ///
    /// Returns [`ConveyorError`] if `density`, `width` or `speed` is
    /// non-finite or not strictly positive, or if `fill_fraction` is not
    /// finite or outside the inclusive range `[f64::MIN_POSITIVE, 1.0]`.
    pub fn from_flat_fill(
        density: f64,
        width: f64,
        fill_fraction: f64,
        speed: f64,
    ) -> Result<Self, ConveyorError> {
        let width = ConveyorError::require_positive("width", width)?;
        let fill_fraction =
            ConveyorError::require_range("fill_fraction", fill_fraction, f64::MIN_POSITIVE, 1.0)?;
        let load_area = fill_fraction * width * width;
        Belt::new(density, load_area, speed)
    }

    /// Volumetric throughput `Q_v = A * v`, in m^3/s.
    pub fn volumetric_flow(&self) -> f64 {
        self.load_area * self.speed
    }

    /// Mass flow rate `mdot = rho * A * v`, in kg/s.
    ///
    /// This is the headline conveyor relation: the mass of material
    /// crossing any fixed station per second.
    pub fn mass_flow(&self) -> f64 {
        self.density * self.load_area * self.speed
    }

    /// Rated capacity in tonnes per hour (t/h), the customary conveyor
    /// design unit.
    ///
    /// Equal to `mass_flow()` converted: `kg/s * 3600 s/h / 1000 kg/t`.
    pub fn capacity_tph(&self) -> f64 {
        self.mass_flow() * 3.6
    }
}

/// Standalone mass-flow helper `mdot = rho * A * v` (kg/s) for callers
/// that do not want to build a [`Belt`].
///
/// # Errors
///
/// Returns [`ConveyorError`] if any argument is non-finite or not
/// strictly positive.
pub fn mass_flow_rate(density: f64, area: f64, velocity: f64) -> Result<f64, ConveyorError> {
    let density = ConveyorError::require_positive("density", density)?;
    let area = ConveyorError::require_positive("area", area)?;
    let velocity = ConveyorError::require_positive("velocity", velocity)?;
    Ok(density * area * velocity)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn mass_flow_matches_rho_a_v() {
        // rho = 1500 kg/m^3, A = 0.20 m^2, v = 2.0 m/s
        // mdot = 1500 * 0.20 * 2.0 = 600 kg/s.
        let belt = Belt::new(1500.0, 0.20, 2.0).expect("valid belt");
        assert!((belt.mass_flow() - 600.0).abs() < EPS);
    }

    #[test]
    fn volumetric_flow_matches_a_v() {
        // A = 0.20 m^2, v = 2.0 m/s -> Q_v = 0.40 m^3/s.
        let belt = Belt::new(1500.0, 0.20, 2.0).expect("valid belt");
        assert!((belt.volumetric_flow() - 0.40).abs() < EPS);
    }

    #[test]
    fn standalone_mass_flow_matches_method() {
        let belt = Belt::new(1200.0, 0.15, 1.5).expect("valid belt");
        let standalone = mass_flow_rate(1200.0, 0.15, 1.5).expect("valid");
        assert!((standalone - belt.mass_flow()).abs() < EPS);
        // 1200 * 0.15 * 1.5 = 270 kg/s.
        assert!((standalone - 270.0).abs() < EPS);
    }

    #[test]
    fn capacity_is_mass_flow_times_3point6() {
        // 600 kg/s -> 600 * 3.6 = 2160 t/h.
        let belt = Belt::new(1500.0, 0.20, 2.0).expect("valid belt");
        assert!((belt.capacity_tph() - 2160.0).abs() < 1e-6);
    }

    #[test]
    fn capacity_scales_linearly_with_speed() {
        // Doubling the belt speed doubles the rated capacity.
        let slow = Belt::new(1500.0, 0.20, 2.0).expect("valid");
        let fast = Belt::new(1500.0, 0.20, 4.0).expect("valid");
        assert!((fast.capacity_tph() - 2.0 * slow.capacity_tph()).abs() < 1e-6);
    }

    #[test]
    fn capacity_scales_linearly_with_area() {
        // Tripling the load area triples the rated capacity.
        let small = Belt::new(1500.0, 0.10, 2.0).expect("valid");
        let big = Belt::new(1500.0, 0.30, 2.0).expect("valid");
        assert!((big.capacity_tph() - 3.0 * small.capacity_tph()).abs() < 1e-6);
    }

    #[test]
    fn flat_fill_area_is_fraction_times_width_squared() {
        // width = 1.0 m, fill = 0.25 -> A = 0.25 m^2.
        let belt = Belt::from_flat_fill(1000.0, 1.0, 0.25, 2.0).expect("valid");
        assert!((belt.load_area - 0.25).abs() < EPS);
        // mdot = 1000 * 0.25 * 2.0 = 500 kg/s.
        assert!((belt.mass_flow() - 500.0).abs() < EPS);
    }

    #[test]
    fn constructors_reject_bad_inputs() {
        assert!(Belt::new(0.0, 0.2, 2.0).is_err());
        assert!(Belt::new(1500.0, -0.2, 2.0).is_err());
        assert!(Belt::new(1500.0, 0.2, f64::NAN).is_err());
        assert!(Belt::from_flat_fill(1000.0, 1.0, 0.0, 2.0).is_err());
        assert!(Belt::from_flat_fill(1000.0, 1.0, 1.5, 2.0).is_err());
        assert!(mass_flow_rate(1.0, 1.0, -1.0).is_err());
    }

    #[test]
    fn gravity_constant_is_standard() {
        assert!((G - 9.80665).abs() < EPS);
    }
}
