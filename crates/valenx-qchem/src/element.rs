//! Periodic-table data for the elements `valenx-qchem` supports.
//!
//! The quantum-chemistry core targets the small-molecule regime, so the
//! built-in tables cover the first ten elements (H–Ne). [`Element`] is a
//! `u8` atomic number behind a thin API: symbol round-trip, atomic mass
//! (for the centre of mass), nuclear charge `Z` (for the integrals and
//! the nuclear-repulsion energy) and a covalent radius (for the
//! geometry-based bond-order guess).
//!
//! ## v1 coverage
//!
//! The element data goes up to argon (Z = 18) so a heavier atom still
//! has a symbol and a mass, but the built-in basis sets only define
//! shells through neon — a basis lookup for an out-of-range element
//! returns [`QchemError::BasisNotFound`](crate::error::QchemError).

use crate::error::{QchemError, Result};
use serde::{Deserialize, Serialize};

/// A chemical element, stored as its atomic number `Z`.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub struct Element(pub u8);

/// `(symbol, atomic mass in unified atomic mass units)` for `Z = 1..18`.
///
/// Masses are the standard atomic weights (CIAAW 2021, rounded).
const ELEMENT_DATA: [(&str, f64); 18] = [
    ("H", 1.008),
    ("He", 4.002_602),
    ("Li", 6.94),
    ("Be", 9.012_183),
    ("B", 10.81),
    ("C", 12.011),
    ("N", 14.007),
    ("O", 15.999),
    ("F", 18.998_403),
    ("Ne", 20.1797),
    ("Na", 22.989_77),
    ("Mg", 24.305),
    ("Al", 26.981_538),
    ("Si", 28.085),
    ("P", 30.973_762),
    ("S", 32.06),
    ("Cl", 35.45),
    ("Ar", 39.95),
];

/// Covalent radii in ångström (`Z = 1..18`), Cordero 2008 values.
///
/// Used only by the geometry-based bond-order guess; two atoms are
/// considered bonded when their separation is below the sum of their
/// covalent radii times a tolerance factor.
const COVALENT_RADII_ANGSTROM: [f64; 18] = [
    0.31, 0.28, // H, He
    1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58, // Li..Ne
    1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06, // Na..Ar
];

impl Element {
    /// Build an [`Element`] from an atomic symbol (case-insensitive,
    /// e.g. `"C"`, `"he"`, `"NA"`).
    ///
    /// # Errors
    ///
    /// Returns [`QchemError::Parse`] when the symbol is unknown.
    pub fn from_symbol(symbol: &str) -> Result<Self> {
        let want = symbol.trim();
        for (i, (sym, _)) in ELEMENT_DATA.iter().enumerate() {
            if sym.eq_ignore_ascii_case(want) {
                return Ok(Element((i + 1) as u8));
            }
        }
        Err(QchemError::parse(
            "xyz",
            format!("unknown element symbol `{symbol}`"),
        ))
    }

    /// Build an [`Element`] from an atomic number `Z` (`1..=18`).
    ///
    /// # Errors
    ///
    /// Returns [`QchemError::InvalidInput`] when `z` is zero or beyond
    /// the supported range.
    pub fn from_z(z: u8) -> Result<Self> {
        if z == 0 || z as usize > ELEMENT_DATA.len() {
            return Err(QchemError::invalid(format!(
                "atomic number {z} out of supported range 1..={}",
                ELEMENT_DATA.len()
            )));
        }
        Ok(Element(z))
    }

    /// The atomic number `Z` — also the nuclear charge in atomic units.
    #[inline]
    pub fn atomic_number(self) -> u8 {
        self.0
    }

    /// Nuclear charge `Z` as an `f64`, for the integral and
    /// nuclear-repulsion routines.
    #[inline]
    pub fn nuclear_charge(self) -> f64 {
        f64::from(self.0)
    }

    /// The element's atomic symbol (`"C"`, `"Ne"`, …).
    pub fn symbol(self) -> &'static str {
        ELEMENT_DATA[(self.0 - 1) as usize].0
    }

    /// Standard atomic mass in unified atomic mass units.
    pub fn atomic_mass(self) -> f64 {
        ELEMENT_DATA[(self.0 - 1) as usize].1
    }

    /// Covalent radius in ångström (Cordero 2008).
    pub fn covalent_radius_angstrom(self) -> f64 {
        COVALENT_RADII_ANGSTROM[(self.0 - 1) as usize]
    }

    /// Number of electrons in the neutral atom — equal to `Z`.
    #[inline]
    pub fn neutral_electron_count(self) -> u32 {
        u32::from(self.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn symbol_round_trip() {
        for z in 1..=18u8 {
            let e = Element::from_z(z).unwrap();
            let back = Element::from_symbol(e.symbol()).unwrap();
            assert_eq!(e, back);
        }
    }

    #[test]
    fn case_insensitive_symbol() {
        assert_eq!(Element::from_symbol("c").unwrap(), Element(6));
        assert_eq!(Element::from_symbol("NE").unwrap(), Element(10));
        assert_eq!(Element::from_symbol(" O ").unwrap(), Element(8));
    }

    #[test]
    fn unknown_symbol_errors() {
        assert!(Element::from_symbol("Xx").is_err());
        assert!(Element::from_z(0).is_err());
        assert!(Element::from_z(99).is_err());
    }

    #[test]
    fn carbon_facts() {
        let c = Element::from_symbol("C").unwrap();
        assert_eq!(c.atomic_number(), 6);
        assert_eq!(c.nuclear_charge(), 6.0);
        assert!((c.atomic_mass() - 12.011).abs() < 1.0e-6);
        assert_eq!(c.neutral_electron_count(), 6);
    }
}
