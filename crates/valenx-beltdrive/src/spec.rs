//! A serialisable belt-drive description and its derived quantities.
//!
//! [`DriveSpec`] bundles the inputs of an open flat-belt drive so a
//! whole configuration can be (de)serialised (via `serde`) and analysed
//! in one call. The derived [`DriveAnalysis`] collects the headline
//! results — speed ratio, belt speed, wrap angles, capstan tension
//! ratio, centrifugal tension, and slipping-limited power — each from
//! the corresponding free function in the topic modules, so the numbers
//! match those APIs exactly.

use crate::error::BeltError;
use crate::{friction, geometry, power};
use serde::{Deserialize, Serialize};

/// Inputs describing an open flat-belt drive on parallel shafts.
///
/// Lengths are in metres, the rotational speed in revolutions per
/// second, tension in newtons, and `linear_density` in kg/m. `mu` is
/// the dimensionless belt/pulley coefficient of friction.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DriveSpec {
    /// Driver (input) pulley pitch diameter, m.
    pub driver_diameter: f64,
    /// Driven (output) pulley pitch diameter, m.
    pub driven_diameter: f64,
    /// Driver rotational speed, rev/s.
    pub driver_rev_per_sec: f64,
    /// Shaft-to-shaft centre distance, m.
    pub center_distance: f64,
    /// Belt/pulley coefficient of friction (dimensionless).
    pub mu: f64,
    /// Belt linear mass density, kg/m.
    pub linear_density: f64,
    /// Maximum allowable tight-side tension, N.
    pub t1_max: f64,
}

/// Headline results derived from a [`DriveSpec`].
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DriveAnalysis {
    /// Transmission ratio `i = D_driven / D_driver`.
    pub speed_ratio: f64,
    /// Belt linear speed at the driver pulley, m/s.
    pub belt_speed: f64,
    /// Driven-pulley rotational speed, rev/s.
    pub driven_rev_per_sec: f64,
    /// Wrap angle on the small pulley, radians.
    pub wrap_small: f64,
    /// Wrap angle on the large pulley, radians.
    pub wrap_large: f64,
    /// Capstan tension ratio `T1/T2 = exp(mu*theta)` on the small
    /// pulley (the limiting pulley, since it has the smaller wrap).
    pub tension_ratio: f64,
    /// Centrifugal tension `Tc = m*v^2`, N.
    pub centrifugal_tension: f64,
    /// Slipping-limited transmitted power, W.
    pub max_power: f64,
}

impl DriveSpec {
    /// Compute the derived [`DriveAnalysis`] for this specification.
    ///
    /// The small pulley's wrap angle is used for the capstan limit
    /// because it is the smaller of the two and therefore governs
    /// slipping.
    ///
    /// # Errors
    ///
    /// Propagates any [`BeltError`] raised by the underlying geometry,
    /// friction, and power calculations (out-of-domain inputs, a centre
    /// distance too small for the pulley sizes, or a centrifugal
    /// tension that exceeds `t1_max`).
    pub fn analyze(&self) -> Result<DriveAnalysis, BeltError> {
        let speed_ratio = geometry::speed_ratio(self.driver_diameter, self.driven_diameter)?;
        let belt_speed = geometry::belt_speed(self.driver_diameter, self.driver_rev_per_sec)?;
        let driven_rev_per_sec = geometry::driven_speed(
            self.driver_diameter,
            self.driven_diameter,
            self.driver_rev_per_sec,
        )?;
        let (wrap_small, wrap_large) = geometry::wrap_angles_open(
            self.driver_diameter / 2.0,
            self.driven_diameter / 2.0,
            self.center_distance,
        )?;
        let tension_ratio = friction::tension_ratio(self.mu, wrap_small)?;
        let centrifugal_tension = power::centrifugal_tension(self.linear_density, belt_speed)?;
        let max_power = power::max_power(
            self.t1_max,
            self.linear_density,
            belt_speed,
            self.mu,
            wrap_small,
        )?;
        Ok(DriveAnalysis {
            speed_ratio,
            belt_speed,
            driven_rev_per_sec,
            wrap_small,
            wrap_large,
            tension_ratio,
            centrifugal_tension,
            max_power,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    fn sample() -> DriveSpec {
        DriveSpec {
            driver_diameter: 0.1,
            driven_diameter: 0.2,
            driver_rev_per_sec: 25.0,
            center_distance: 0.5,
            mu: 0.3,
            linear_density: 0.4,
            t1_max: 1200.0,
        }
    }

    #[test]
    fn analysis_fields_match_the_free_functions() {
        let spec = sample();
        let a = spec.analyze().unwrap();

        let i = geometry::speed_ratio(0.1, 0.2).unwrap();
        let v = geometry::belt_speed(0.1, 25.0).unwrap();
        let (ts, tl) = geometry::wrap_angles_open(0.05, 0.1, 0.5).unwrap();
        let k = friction::tension_ratio(0.3, ts).unwrap();
        let tc = power::centrifugal_tension(0.4, v).unwrap();
        let pmax = power::max_power(1200.0, 0.4, v, 0.3, ts).unwrap();

        assert!((a.speed_ratio - i).abs() < EPS);
        assert!((a.belt_speed - v).abs() < EPS);
        assert!((a.wrap_small - ts).abs() < EPS);
        assert!((a.wrap_large - tl).abs() < EPS);
        assert!((a.tension_ratio - k).abs() < EPS);
        assert!((a.centrifugal_tension - tc).abs() < EPS);
        assert!((a.max_power - pmax).abs() < EPS);
    }

    #[test]
    fn spec_round_trips_through_json() {
        let spec = sample();
        let json = serde_json::to_string(&spec).unwrap();
        let back: DriveSpec = serde_json::from_str(&json).unwrap();
        assert_eq!(spec, back);
    }

    #[test]
    fn analysis_round_trips_through_json() {
        let a = sample().analyze().unwrap();
        let json = serde_json::to_string(&a).unwrap();
        let back: DriveAnalysis = serde_json::from_str(&json).unwrap();
        assert_eq!(a, back);
    }

    #[test]
    fn analyze_propagates_degenerate_geometry() {
        let mut spec = sample();
        // Centre distance smaller than the radius difference (0.05 m).
        spec.center_distance = 0.04;
        let err = spec.analyze().unwrap_err();
        assert_eq!(err.code(), "beltdrive.degenerate_geometry");
    }
}
