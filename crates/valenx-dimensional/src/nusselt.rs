//! Nusselt number — dimensionless convective heat-transfer coefficient.
//!
//! The Nusselt number
//!
//! ```text
//! Nu = h L / k
//! ```
//!
//! is the ratio of convective to conductive heat transfer across a
//! boundary layer. Here `h` is the convective heat-transfer
//! coefficient, `L` the characteristic length, and `k` the thermal
//! conductivity of the fluid. A value of `Nu = 1` corresponds to pure
//! conduction across the layer; larger values mean convection enhances
//! the transfer. It is dimensionless in any consistent unit system.
//!
//! Note that the Nusselt and Biot numbers share the algebraic form
//! `h L / k` but use different conductivities: the Nusselt number uses
//! the **fluid** conductivity, whereas the Biot number (see
//! [`crate::biot`]) uses the **solid** conductivity. They are physically
//! distinct despite the identical formula.

use crate::error::{require_non_negative, require_positive, DimensionlessError};
use serde::{Deserialize, Serialize};

/// Nusselt number `Nu = h L / k`, dimensionless.
///
/// Construct with [`Nusselt::new`]. The inner value is always finite and
/// non-negative.
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct Nusselt(f64);

impl Nusselt {
    /// Build a Nusselt number from `Nu = h L / k`.
    ///
    /// - `heat_transfer_coefficient` (`h`) must be finite and
    ///   non-negative (it may be zero).
    /// - `length` (`L`) must be strictly positive.
    /// - `thermal_conductivity` (`k`) must be strictly positive (it is
    ///   the denominator) and refers to the **fluid** conductivity.
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
        Ok(Nusselt(h * length / k))
    }

    /// The raw dimensionless value.
    pub fn value(&self) -> f64 {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn nusselt_matches_definition() {
        // h=100, L=0.5, k=0.6  ->  Nu = 100*0.5/0.6 = 83.333...
        let nu = Nusselt::new(100.0, 0.5, 0.6).unwrap();
        assert!((nu.value() - (100.0 * 0.5 / 0.6)).abs() < 1e-9);
    }

    #[test]
    fn pure_conduction_gives_unit_value() {
        // h L = k  ->  Nu = 1.
        let nu = Nusselt::new(0.6, 1.0, 0.6).unwrap();
        assert!((nu.value() - 1.0).abs() < EPS);
    }

    #[test]
    fn zero_coefficient_gives_zero() {
        let nu = Nusselt::new(0.0, 1.0, 1.0).unwrap();
        assert!(nu.value().abs() < EPS);
    }

    #[test]
    fn unit_inputs_give_unit_value() {
        let nu = Nusselt::new(1.0, 1.0, 1.0).unwrap();
        assert!((nu.value() - 1.0).abs() < EPS);
    }

    #[test]
    fn rejects_non_positive_conductivity() {
        let err = Nusselt::new(1.0, 1.0, -0.5).unwrap_err();
        assert_eq!(err.parameter(), "thermal_conductivity");
        assert_eq!(err.code(), "dimensionless.out-of-domain");
    }

    #[test]
    fn rejects_non_positive_length() {
        let err = Nusselt::new(1.0, 0.0, 1.0).unwrap_err();
        assert_eq!(err.parameter(), "length");
    }

    #[test]
    fn rejects_negative_coefficient() {
        let err = Nusselt::new(-1.0, 1.0, 1.0).unwrap_err();
        assert_eq!(err.parameter(), "heat_transfer_coefficient");
    }

    #[test]
    fn serde_round_trip() {
        let nu = Nusselt::new(100.0, 0.5, 0.6).unwrap();
        let json = serde_json::to_string(&nu).unwrap();
        let back: Nusselt = serde_json::from_str(&json).unwrap();
        assert!((nu.value() - back.value()).abs() < EPS);
    }
}
