//! A simple planetary (epicyclic) gear set with a **fixed ring**.
//!
//! ## Model
//!
//! A planetary set has three coaxial members: a central **sun** gear, an
//! outer internal-tooth **ring** (annulus), and a **carrier** that holds
//! the orbiting planet gears. Their speeds obey the Willis fundamental
//! train equation:
//!
//! ```text
//! (w_sun - w_carrier) / (w_ring - w_carrier) = -N_ring / N_sun
//! ```
//!
//! The planet tooth count cancels out of the speed relationship (the
//! planets are idlers), so it does not appear. Geometry still requires
//! `N_ring = N_sun + 2 * N_planet`, hence `N_ring > N_sun`.
//!
//! This crate models the most common configuration: **ring held fixed,
//! sun is the input, carrier is the output**. Setting `w_ring = 0` and
//! solving for the input/output speed ratio gives
//!
//! ```text
//! ratio = w_sun / w_carrier = 1 + N_ring / N_sun
//! ```
//!
//! a pure reduction (always `> 1`). As with a parallel stage, speed
//! divides by the ratio and torque multiplies by it (times efficiency):
//!
//! ```text
//! carrier_speed  = sun_speed  / ratio
//! carrier_torque = sun_torque * ratio * efficiency
//! ```

use serde::{Deserialize, Serialize};

use crate::error::{check_efficiency, GearboxError};

/// A fixed-ring planetary set: sun input, carrier output.
///
/// Build with [`PlanetarySet::new`] (lossless) or
/// [`PlanetarySet::with_efficiency`]. Geometry is validated:
/// `ring > sun >= 1` is required, so the constructors never produce a
/// physically impossible set and never panic.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct PlanetarySet {
    sun_teeth: u32,
    ring_teeth: u32,
    efficiency: f64,
}

impl PlanetarySet {
    /// Build an **ideal (lossless)** fixed-ring planetary set.
    ///
    /// Efficiency defaults to `1.0`.
    ///
    /// # Errors
    ///
    /// [`GearboxError::DegeneratePlanetary`] if `sun_teeth == 0` or
    /// `ring_teeth <= sun_teeth` (a real ring must be larger than the
    /// sun to leave room for the planets).
    pub fn new(sun_teeth: u32, ring_teeth: u32) -> Result<Self, GearboxError> {
        Self::with_efficiency(sun_teeth, ring_teeth, 1.0)
    }

    /// Build a fixed-ring planetary set with an explicit `efficiency`
    /// in `(0, 1]`.
    ///
    /// # Errors
    ///
    /// - [`GearboxError::DegeneratePlanetary`] if `sun_teeth == 0` or
    ///   `ring_teeth <= sun_teeth`.
    /// - [`GearboxError::EfficiencyOutOfRange`] if `efficiency` is not
    ///   in `(0, 1]`.
    /// - [`GearboxError::NotFinite`] if `efficiency` is `NaN`/`±∞`.
    pub fn with_efficiency(
        sun_teeth: u32,
        ring_teeth: u32,
        efficiency: f64,
    ) -> Result<Self, GearboxError> {
        if sun_teeth == 0 {
            return Err(GearboxError::DegeneratePlanetary {
                reason: "sun must have >= 1 tooth".to_string(),
            });
        }
        if ring_teeth <= sun_teeth {
            return Err(GearboxError::DegeneratePlanetary {
                reason: format!("ring ({ring_teeth}) must be larger than sun ({sun_teeth})"),
            });
        }
        let efficiency = check_efficiency(efficiency)?;
        Ok(Self {
            sun_teeth,
            ring_teeth,
            efficiency,
        })
    }

    /// Sun (central) gear tooth count.
    pub fn sun_teeth(&self) -> u32 {
        self.sun_teeth
    }

    /// Ring (annulus) internal tooth count.
    pub fn ring_teeth(&self) -> u32 {
        self.ring_teeth
    }

    /// Mechanical efficiency in `(0, 1]`.
    pub fn efficiency(&self) -> f64 {
        self.efficiency
    }

    /// Per-planet tooth count implied by the meshing constraint
    /// `N_ring = N_sun + 2 * N_planet`, if it comes out to a whole
    /// number of teeth.
    ///
    /// Returns `Some(n)` when `ring - sun` is even (so the planets are
    /// integer-toothed and the sun and ring are coaxial), otherwise
    /// `None`. The speed ratio does not depend on this value — it is
    /// reported only as a geometry aid.
    pub fn planet_teeth(&self) -> Option<u32> {
        let diff = self.ring_teeth - self.sun_teeth;
        if diff % 2 == 0 {
            Some(diff / 2)
        } else {
            None
        }
    }

    /// Reduction ratio (sun speed / carrier speed) for the fixed-ring,
    /// sun-input, carrier-output case: `1 + ring / sun` (Willis).
    ///
    /// Always greater than `1`.
    pub fn ratio(&self) -> f64 {
        1.0 + self.ring_teeth as f64 / self.sun_teeth as f64
    }

    /// Carrier (output) speed for a given sun (input) speed:
    /// `sun_speed / ratio`.
    pub fn carrier_speed(&self, sun_speed: f64) -> f64 {
        sun_speed / self.ratio()
    }

    /// Carrier (output) torque for a given sun (input) torque:
    /// `sun_torque * ratio * efficiency`.
    pub fn carrier_torque(&self, sun_torque: f64) -> f64 {
        sun_torque * self.ratio() * self.efficiency
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    #[test]
    fn willis_fixed_ring_ratio() {
        // Sun 24, ring 72 -> 1 + 72/24 = 4.0.
        let p = PlanetarySet::new(24, 72).unwrap();
        assert!((p.ratio() - 4.0).abs() < EPS);
    }

    #[test]
    fn ratio_equals_one_plus_ring_over_sun() {
        // Sun 20, ring 80 -> 1 + 4 = 5.0.
        let p = PlanetarySet::new(20, 80).unwrap();
        let expected = 1.0 + 80.0 / 20.0;
        assert!((p.ratio() - expected).abs() < EPS);
        assert!((p.ratio() - 5.0).abs() < EPS);
    }

    #[test]
    fn ratio_always_above_one() {
        // Even the smallest legal ring (sun + 1) gives a reduction.
        let p = PlanetarySet::new(30, 31).unwrap();
        assert!(p.ratio() > 1.0);
        assert!((p.ratio() - (1.0 + 31.0 / 30.0)).abs() < EPS);
    }

    #[test]
    fn cross_check_against_willis_general_form() {
        // Derive carrier speed directly from the general Willis
        // equation with w_ring = 0 and compare to the closed form.
        // (w_sun - w_c)/(0 - w_c) = -N_ring/N_sun
        //   => w_c = w_sun * N_sun / (N_sun + N_ring)
        let sun = 18u32;
        let ring = 90u32;
        let w_sun = 3000.0;
        let p = PlanetarySet::new(sun, ring).unwrap();
        let w_c_general = w_sun * sun as f64 / (sun as f64 + ring as f64);
        assert!((p.carrier_speed(w_sun) - w_c_general).abs() < 1e-9);
    }

    #[test]
    fn speed_divides_torque_multiplies() {
        // Sun 24, ring 72 -> ratio 4. Lossless.
        let p = PlanetarySet::new(24, 72).unwrap();
        assert!((p.carrier_speed(2000.0) - 500.0).abs() < EPS);
        assert!((p.carrier_torque(10.0) - 40.0).abs() < EPS);
        // Ideal power conservation.
        let p_in = 2000.0 * 10.0;
        let p_out = p.carrier_speed(2000.0) * p.carrier_torque(10.0);
        assert!((p_out - p_in).abs() < 1e-9);
    }

    #[test]
    fn efficiency_reduces_carrier_torque() {
        let ideal = PlanetarySet::new(24, 72).unwrap();
        let lossy = PlanetarySet::with_efficiency(24, 72, 0.9).unwrap();
        // Speed unaffected by efficiency.
        assert!((lossy.carrier_speed(2000.0) - ideal.carrier_speed(2000.0)).abs() < EPS);
        // Torque scaled: 10 * 4 * 0.9 = 36.
        assert!((lossy.carrier_torque(10.0) - 36.0).abs() < EPS);
        assert!(lossy.carrier_torque(10.0) < ideal.carrier_torque(10.0));
    }

    #[test]
    fn planet_teeth_even_diff() {
        // ring - sun = 48, even -> 24-tooth planets.
        let p = PlanetarySet::new(24, 72).unwrap();
        assert_eq!(p.planet_teeth(), Some(24));
    }

    #[test]
    fn planet_teeth_odd_diff_is_none() {
        // ring - sun = 47, odd -> no integer planet count.
        let p = PlanetarySet::new(24, 71).unwrap();
        assert_eq!(p.planet_teeth(), None);
    }

    #[test]
    fn zero_sun_rejected() {
        let err = PlanetarySet::new(0, 50).unwrap_err();
        assert_eq!(err.code(), "gearbox.degenerate_planetary");
        assert!(matches!(err, GearboxError::DegeneratePlanetary { .. }));
    }

    #[test]
    fn ring_not_larger_than_sun_rejected() {
        // Equal counts.
        assert!(PlanetarySet::new(40, 40).is_err());
        // Ring smaller than sun.
        let err = PlanetarySet::new(40, 30).unwrap_err();
        assert_eq!(err.code(), "gearbox.degenerate_planetary");
    }

    #[test]
    fn efficiency_out_of_range_rejected() {
        assert!(PlanetarySet::with_efficiency(24, 72, 1.5).is_err());
        assert!(PlanetarySet::with_efficiency(24, 72, 0.0).is_err());
    }

    #[test]
    fn serde_round_trip() {
        let p = PlanetarySet::with_efficiency(24, 72, 0.96).unwrap();
        let json = serde_json::to_string(&p).unwrap();
        let back: PlanetarySet = serde_json::from_str(&json).unwrap();
        assert_eq!(p, back);
    }
}
