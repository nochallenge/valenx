//! Hydrocarbon fuel description and molar-mass constants.
//!
//! A fuel here is a single pure hydrocarbon `CxHy` (methane, propane,
//! octane, ...). The molar masses use standard IUPAC atomic weights;
//! the air model is the textbook simplification of dry air as
//! `O2 + 3.76 N2` (i.e. 21% / 79% by mole), with everything except
//! oxygen lumped into the nitrogen carrier.

use serde::{Deserialize, Serialize};

use crate::error::CombustionError;

/// Molar mass of atomic carbon, g/mol (IUPAC standard atomic weight).
pub const M_C: f64 = 12.011;
/// Molar mass of atomic hydrogen, g/mol.
pub const M_H: f64 = 1.008;
/// Molar mass of molecular oxygen O2, g/mol.
pub const M_O2: f64 = 31.998;
/// Molar mass of molecular nitrogen N2, g/mol.
pub const M_N2: f64 = 28.013;
/// Molar mass of carbon dioxide CO2, g/mol.
pub const M_CO2: f64 = 44.009;
/// Molar mass of water H2O, g/mol.
pub const M_H2O: f64 = 18.015;

/// Moles of N2 carried per mole of O2 in the textbook dry-air model
/// (79.0 / 21.0). Multiply O2 demand by this to get the N2 that rides
/// along inert.
pub const N2_PER_O2: f64 = 3.76;

/// Apparent molar mass of one mole of "air" expressed as
/// `O2 + 3.76 N2`, divided by the (1 + 3.76) total moles — i.e. the
/// mean molar mass of the air mixture, g/mol. Used so that air masses
/// can also be reasoned about per total mole of air if desired.
pub const M_AIR_PER_MOLE: f64 = (M_O2 + N2_PER_O2 * M_N2) / (1.0 + N2_PER_O2);

/// A pure hydrocarbon fuel `CxHy`.
///
/// Construct via [`Fuel::new`], which validates the atom counts, or use
/// one of the named constructors ([`Fuel::methane`], etc.).
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Fuel {
    /// Carbon atom count `x` (>= 1).
    pub carbon: u32,
    /// Hydrogen atom count `y` (>= 1).
    pub hydrogen: u32,
}

impl Fuel {
    /// Build a fuel `CxHy`, validating that `carbon >= 1` and
    /// `hydrogen >= 1`.
    ///
    /// # Errors
    ///
    /// Returns [`CombustionError::InvalidFuel`] when either atom count is
    /// below its physical minimum.
    pub fn new(carbon: u32, hydrogen: u32) -> Result<Self, CombustionError> {
        if carbon < 1 {
            return Err(CombustionError::InvalidFuel {
                carbon,
                hydrogen,
                reason: "carbon atom count must be at least 1",
            });
        }
        if hydrogen < 1 {
            return Err(CombustionError::InvalidFuel {
                carbon,
                hydrogen,
                reason: "hydrogen atom count must be at least 1",
            });
        }
        Ok(Self { carbon, hydrogen })
    }

    /// Methane, CH4.
    pub fn methane() -> Self {
        Self {
            carbon: 1,
            hydrogen: 4,
        }
    }

    /// Propane, C3H8.
    pub fn propane() -> Self {
        Self {
            carbon: 3,
            hydrogen: 8,
        }
    }

    /// iso-Octane proxy, C8H18 (the gasoline surrogate).
    pub fn octane() -> Self {
        Self {
            carbon: 8,
            hydrogen: 18,
        }
    }

    /// Molar mass of the fuel, g/mol: `x * M_C + y * M_H`.
    pub fn molar_mass(&self) -> f64 {
        self.carbon as f64 * M_C + self.hydrogen as f64 * M_H
    }

    /// Moles of O2 required for stoichiometric (complete) combustion of
    /// one mole of fuel: `x + y / 4`.
    ///
    /// From `CxHy + (x + y/4) O2 -> x CO2 + (y/2) H2O`.
    pub fn stoich_o2_moles(&self) -> f64 {
        self.carbon as f64 + self.hydrogen as f64 / 4.0
    }
}
