//! Biot number — ratio of internal conduction resistance to surface
//! convection resistance.
//!
//! The Biot number
//!
//! ```text
//! Bi = h L / k
//! ```
//!
//! compares the conductive resistance **inside** a solid body to the
//! convective resistance at its **surface**. Here `h` is the convective
//! heat-transfer coefficient, `L` a characteristic length (often the
//! volume-to-surface-area ratio), and `k` the conductivity of the
//! **solid**. A small Biot number means the body conducts heat
//! internally much faster than it exchanges heat at the surface, so its
//! temperature is nearly uniform — the assumption behind the
//! lumped-capacitance model.
//!
//! The conventional rule of thumb is that lumped-capacitance analysis is
//! acceptable when `Bi < 0.1`; see [`Biot::lumped_capacitance`].
//!
//! The formula `h L / k` is identical to that of the Nusselt number (see
//! [`crate::nusselt`]) but the conductivity here is the solid's, not the
//! fluid's, so the two groups describe different physics.

use crate::error::{require_non_negative, require_positive, DimensionlessError};
use serde::{Deserialize, Serialize};

/// Conventional upper bound on the Biot number for the
/// lumped-capacitance approximation to be acceptable. The classical
/// rule-of-thumb value is 0.1.
pub const LUMPED_CAPACITANCE_LIMIT: f64 = 0.1;

/// Biot number `Bi = h L / k`, dimensionless.
///
/// Construct with [`Biot::new`]. The inner value is always finite and
/// non-negative.
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct Biot(f64);

impl Biot {
    /// Build a Biot number from `Bi = h L / k`.
    ///
    /// - `heat_transfer_coefficient` (`h`) must be finite and
    ///   non-negative (it may be zero).
    /// - `length` (`L`) must be strictly positive.
    /// - `thermal_conductivity` (`k`) must be strictly positive (it is
    ///   the denominator) and refers to the **solid** conductivity.
    ///
    /// # Errors
    ///
    /// Returns [`DimensionlessError`] if any input is non-finite or
    /// violates the domain above.
    pub fn new(
        heat_transfer_coefficient: f64,
        length: f64,
        thermal_conductivity: f64,
    ) -> Result<Self, DimensionlessError> {
        let h = require_non_negative("heat_transfer_coefficient", heat_transfer_coefficient)?;
        let length = require_positive("length", length)?;
        let k = require_positive("thermal_conductivity", thermal_conductivity)?;
        Ok(Biot(h * length / k))
    }

    /// The raw dimensionless value.
    pub fn value(&self) -> f64 {
        self.0
    }

    /// Classify whether the lumped-capacitance (spatially uniform
    /// temperature) approximation is reasonable, using the conventional
    /// threshold [`LUMPED_CAPACITANCE_LIMIT`] (`Bi < 0.1`). This is a
    /// rule of thumb, not a hard physical guarantee.
    pub fn lumped_capacitance(&self) -> LumpedCapacitance {
        if self.0 < LUMPED_CAPACITANCE_LIMIT {
            LumpedCapacitance::Valid
        } else {
            LumpedCapacitance::Invalid
        }
    }
}

/// Verdict of the lumped-capacitance validity test for a [`Biot`]
/// number.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum LumpedCapacitance {
    /// `Bi < 0.1`: internal temperature gradients are small enough that
    /// the lumped-capacitance model is a reasonable approximation.
    Valid,
    /// `Bi >= 0.1`: internal gradients matter; a spatially resolved
    /// (distributed) model should be used instead.
    Invalid,
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn biot_matches_definition() {
        // h=10, L=0.01, k=200  ->  Bi = 10*0.01/200 = 5e-4
        let bi = Biot::new(10.0, 0.01, 200.0).unwrap();
        assert!((bi.value() - 5.0e-4).abs() < 1e-12);
    }

    #[test]
    fn unit_inputs_give_unit_value() {
        let bi = Biot::new(1.0, 1.0, 1.0).unwrap();
        assert!((bi.value() - 1.0).abs() < EPS);
    }

    #[test]
    fn small_biot_is_lumped_valid() {
        // A small, highly conductive metal piece -> very low Bi.
        let bi = Biot::new(10.0, 0.01, 200.0).unwrap();
        assert_eq!(bi.lumped_capacitance(), LumpedCapacitance::Valid);
    }

    #[test]
    fn large_biot_is_lumped_invalid() {
        // A poorly conducting body with strong convection -> high Bi.
        let bi = Biot::new(500.0, 0.1, 0.2).unwrap();
        assert!(bi.value() > LUMPED_CAPACITANCE_LIMIT);
        assert_eq!(bi.lumped_capacitance(), LumpedCapacitance::Invalid);
    }

    #[test]
    fn at_threshold_is_invalid() {
        // Exactly 0.1 is treated as invalid (test is `< 0.1`).
        let bi = Biot::new(0.1, 1.0, 1.0).unwrap();
        assert!((bi.value() - LUMPED_CAPACITANCE_LIMIT).abs() < EPS);
        assert_eq!(bi.lumped_capacitance(), LumpedCapacitance::Invalid);
    }

    #[test]
    fn just_below_threshold_is_valid() {
        let bi = Biot::new(0.099, 1.0, 1.0).unwrap();
        assert_eq!(bi.lumped_capacitance(), LumpedCapacitance::Valid);
    }

    #[test]
    fn rejects_non_positive_conductivity() {
        let err = Biot::new(1.0, 1.0, 0.0).unwrap_err();
        assert_eq!(err.parameter(), "thermal_conductivity");
        assert_eq!(err.code(), "dimensionless.out-of-domain");
    }

    #[test]
    fn rejects_non_finite_length() {
        let err = Biot::new(1.0, f64::NAN, 1.0).unwrap_err();
        assert_eq!(err.code(), "dimensionless.not-finite");
        assert_eq!(err.parameter(), "length");
    }
}
