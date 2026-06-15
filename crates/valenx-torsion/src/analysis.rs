//! High-level bundled analysis of one shaft under one load case.
//!
//! [`TorsionCase`] groups a [`Shaft`] with an applied torque, a length, a
//! shear modulus and an angular speed, and [`TorsionCase::analyse`]
//! evaluates every closed-form quantity at once into a [`TorsionResult`].
//! This is the convenient entry point for a UI or CLI that wants the
//! whole answer in a single call.

use serde::{Deserialize, Serialize};

use crate::error::{require_positive, TorsionError};
use crate::response::{angle_of_twist, max_shear_stress, power, torsional_rigidity};
use crate::shaft::Shaft;

/// A complete torsion load case: a shaft plus the loads applied to it.
///
/// All scalar fields must be finite and strictly positive; this is
/// checked when the case is [`analyse`](TorsionCase::analyse)d. Keep the
/// units consistent across fields (e.g. SI throughout).
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TorsionCase {
    /// The shaft cross-section.
    pub shaft: Shaft,
    /// Applied torque `T` (`> 0`).
    pub torque: f64,
    /// Shaft length `L` over which the twist accumulates (`> 0`).
    pub length: f64,
    /// Shear modulus `G` of the material (`> 0`).
    pub shear_modulus: f64,
    /// Angular speed `omega` for the power calculation (`> 0`).
    pub angular_speed: f64,
}

impl TorsionCase {
    /// Build a load case, validating every scalar input up front.
    ///
    /// # Errors
    ///
    /// Returns [`TorsionError::NonPositive`] for the first non-positive or
    /// non-finite scalar encountered.
    pub fn new(
        shaft: Shaft,
        torque: f64,
        length: f64,
        shear_modulus: f64,
        angular_speed: f64,
    ) -> Result<Self, TorsionError> {
        let torque = require_positive("torque", torque)?;
        let length = require_positive("length", length)?;
        let shear_modulus = require_positive("shear_modulus", shear_modulus)?;
        let angular_speed = require_positive("angular_speed", angular_speed)?;
        Ok(Self {
            shaft,
            torque,
            length,
            shear_modulus,
            angular_speed,
        })
    }

    /// Evaluate all closed-form torsion quantities for this case.
    ///
    /// # Errors
    ///
    /// Propagates any [`TorsionError`] from the underlying formulas (in
    /// practice none, since [`TorsionCase::new`] already validated the
    /// inputs).
    pub fn analyse(&self) -> Result<TorsionResult, TorsionError> {
        Ok(TorsionResult {
            polar_moment: self.shaft.polar_moment(),
            max_shear_stress: max_shear_stress(&self.shaft, self.torque)?,
            angle_of_twist: angle_of_twist(
                &self.shaft,
                self.torque,
                self.length,
                self.shear_modulus,
            )?,
            torsional_rigidity: torsional_rigidity(&self.shaft, self.shear_modulus)?,
            power: power(self.torque, self.angular_speed)?,
        })
    }
}

/// The bundled results of a [`TorsionCase`] analysis.
///
/// Every field is in the unit system implied by the inputs (SI gives `J`
/// in m⁴, stress in Pa, twist in rad, rigidity in N·m², power in W).
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TorsionResult {
    /// Polar second moment of area `J` of the section.
    pub polar_moment: f64,
    /// Maximum (outer-surface) shear stress `tau_max = T r / J`.
    pub max_shear_stress: f64,
    /// Angle of twist `theta = T L / (G J)`.
    pub angle_of_twist: f64,
    /// Torsional rigidity `G J`.
    pub torsional_rigidity: f64,
    /// Transmitted power `P = T omega`.
    pub power: f64,
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    const EPS: f64 = 1e-9;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < EPS * b.abs().max(1.0)
    }

    #[test]
    fn analysis_reproduces_each_closed_form() {
        let d = 30.0_f64;
        let shaft = Shaft::solid(d).unwrap();
        let case = TorsionCase::new(shaft, 250.0, 1_500.0, 79_300.0, 12.0).unwrap();
        let r = case.analyse().unwrap();

        let j = PI * d.powi(4) / 32.0;
        assert!(close(r.polar_moment, j));
        assert!(close(r.max_shear_stress, 250.0 * (d / 2.0) / j));
        assert!(close(r.angle_of_twist, 250.0 * 1_500.0 / (79_300.0 * j)));
        assert!(close(r.torsional_rigidity, 79_300.0 * j));
        assert!(close(r.power, 250.0 * 12.0));
    }

    #[test]
    fn new_rejects_non_positive_scalars() {
        let shaft = Shaft::solid(10.0).unwrap();
        assert!(matches!(
            TorsionCase::new(shaft, -1.0, 1.0, 1.0, 1.0),
            Err(TorsionError::NonPositive { name: "torque", .. })
        ));
        assert!(matches!(
            TorsionCase::new(shaft, 1.0, 1.0, 1.0, 0.0),
            Err(TorsionError::NonPositive {
                name: "angular_speed",
                ..
            })
        ));
    }

    #[test]
    fn result_round_trips_through_json() {
        let shaft = Shaft::hollow(40.0, 20.0).unwrap();
        let case = TorsionCase::new(shaft, 333.0, 800.0, 81_000.0, 9.5).unwrap();
        let result = case.analyse().unwrap();

        let json = serde_json::to_string(&result).unwrap();
        let back: TorsionResult = serde_json::from_str(&json).unwrap();

        // serde_json's decimal text form is not guaranteed bit-exact for
        // every f64, so compare each field within a tight tolerance rather
        // than with `==`.
        assert!(close(back.polar_moment, result.polar_moment));
        assert!(close(back.max_shear_stress, result.max_shear_stress));
        assert!(close(back.angle_of_twist, result.angle_of_twist));
        assert!(close(back.torsional_rigidity, result.torsional_rigidity));
        assert!(close(back.power, result.power));
    }

    #[test]
    fn case_round_trips_through_json() {
        let shaft = Shaft::solid(15.0).unwrap();
        let case = TorsionCase::new(shaft, 120.0, 600.0, 79_300.0, 7.0).unwrap();

        let json = serde_json::to_string(&case).unwrap();
        let back: TorsionCase = serde_json::from_str(&json).unwrap();

        // Compare the reconstructed scalars within tolerance; the decimal
        // text form is not guaranteed to be bit-exact for every f64.
        assert_eq!(back.shaft, case.shaft);
        assert!(close(back.torque, case.torque));
        assert!(close(back.length, case.length));
        assert!(close(back.shear_modulus, case.shear_modulus));
        assert!(close(back.angular_speed, case.angular_speed));
    }
}
