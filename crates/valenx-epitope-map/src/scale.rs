//! Per-residue propensity scales.

use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::error::EpitopeError;

/// The 20 standard amino acids in canonical column order.
pub const ALPHABET: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";

/// Column index of `residue` (either case) in `0..20`, or `None`.
pub fn aa_index(residue: u8) -> Option<usize> {
    let up = residue.to_ascii_uppercase();
    ALPHABET.iter().position(|&a| a == up)
}

/// A per-residue propensity scale (higher = more epitope-prone, for a
/// hydrophilicity scale).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PropensityScale {
    /// Values in [`ALPHABET`] order.
    values: [f64; 20],
}

impl PropensityScale {
    /// Build from a map that must provide a finite value for every standard
    /// amino acid.
    pub fn from_map(map: &HashMap<char, f64>) -> Result<Self, EpitopeError> {
        let mut values = [0.0_f64; 20];
        for (i, &b) in ALPHABET.iter().enumerate() {
            let aa = b as char;
            let v = *map
                .get(&aa)
                .ok_or(EpitopeError::ScaleIncomplete { residue: aa })?;
            if !v.is_finite() {
                return Err(EpitopeError::NonFinite {
                    what: "scale value",
                });
            }
            values[i] = v;
        }
        Ok(Self { values })
    }

    /// The propensity of `residue`, or `None` if non-standard.
    pub fn value(&self, residue: u8) -> Option<f64> {
        aa_index(residue).map(|i| self.values[i])
    }
}

/// The built-in default: a hydrophilicity scale equal to the **negated
/// Kyte-Doolittle** hydropathy (positive = hydrophilic/surface-exposed). See the
/// crate docs for why this stands in for the original Parker/Hopp-Woods scales.
pub fn hydrophilicity_kd() -> PropensityScale {
    PropensityScale {
        // negated Kyte-Doolittle, in ALPHABET order A C D E F G H I K L M N P Q R S T V W Y
        values: [
            -1.8, -2.5, 3.5, 3.5, -2.8, 0.4, 3.2, -4.5, 3.9, -3.8, -1.9, 3.5, 1.6, 3.5, 4.5, 0.8,
            0.7, -4.2, 0.9, 1.3,
        ],
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hydrophilic_residues_are_positive() {
        let s = hydrophilicity_kd();
        // charged/polar residues are hydrophilic (positive)
        for r in [b'D', b'E', b'K', b'R', b'N', b'Q', b'H'] {
            assert!(s.value(r).unwrap() > 0.0, "{}", r as char);
        }
        // aliphatic residues are hydrophobic (negative)
        for r in [b'I', b'L', b'V', b'F'] {
            assert!(s.value(r).unwrap() < 0.0, "{}", r as char);
        }
    }

    #[test]
    fn from_map_requires_all_twenty() {
        let mut m = HashMap::new();
        m.insert('A', 1.0);
        assert_eq!(
            PropensityScale::from_map(&m).unwrap_err().code(),
            "scale_incomplete"
        );
    }

    #[test]
    fn from_map_round_trips() {
        let mut m = HashMap::new();
        for &b in ALPHABET {
            m.insert(b as char, 1.0);
        }
        let s = PropensityScale::from_map(&m).unwrap();
        assert!((s.value(b'A').unwrap() - 1.0).abs() < 1e-12);
    }
}
