//! Canonical structure type: atomic / residue / chain hierarchy for
//! proteins, nucleic acids, and small molecules.
//!
//! The hierarchy mirrors the PDB and mmCIF conventions:
//!
//! - [`Structure`] holds an id and an ordered list of [`Chain`]s.
//! - [`Chain`] holds a one-character chain id and ordered [`Residue`]s.
//! - [`Residue`] holds a 3-letter residue code (`ALA`, `GLY`, `DA`, …)
//!   and ordered [`Atom`]s.
//! - [`Atom`] holds the atom name (`CA`, `N`, `O`, …), element symbol,
//!   3D position, and B-factor.
//!
//! Format readers + writers live in [`crate::format::pdb`] and
//! `crate::format::mmcif` (mmCIF lands in a follow-up).

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// Top-level structure — atoms grouped by residue, residues by chain.
#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq)]
pub struct Structure {
    /// Source identifier (e.g. PDB id, mmCIF entry id, file basename).
    pub id: String,
    /// Chains in source order.
    pub chains: Vec<Chain>,
}

/// One chain (e.g. a single polypeptide or nucleic-acid strand).
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Chain {
    /// Single-character chain id from the PDB record (`A`, `B`, …).
    pub id: char,
    /// Residues in N-terminus → C-terminus (or 5' → 3') order.
    pub residues: Vec<Residue>,
}

/// One residue (amino acid, nucleotide, or ligand).
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Residue {
    /// 3-letter residue code from PDB (`ALA`, `GLY`, `HIS`, `DA`, `RG`).
    pub name: String,
    /// 1-based residue number from the source file.
    pub seq_id: i32,
    /// Atoms in the residue (source order).
    pub atoms: Vec<Atom>,
}

/// One atom — name + element + 3D position + B-factor.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Atom {
    /// Atom name (`CA`, `N`, `C`, `O`, `OG1`, …).
    pub name: String,
    /// Element symbol (`C`, `N`, `O`, `H`, `S`, …).
    pub element: String,
    /// Cartesian coordinates in Ångströms.
    pub position: Vector3<f64>,
    /// B-factor (temperature factor) — keeps PDB → cryoEM round-trips clean.
    #[serde(default)]
    pub b_factor: f64,
}

impl Structure {
    /// Total atom count across all chains and residues.
    pub fn atom_count(&self) -> usize {
        self.chains
            .iter()
            .flat_map(|c| c.residues.iter())
            .map(|r| r.atoms.len())
            .sum()
    }

    /// Total residue count across all chains.
    pub fn residue_count(&self) -> usize {
        self.chains.iter().map(|c| c.residues.len()).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_structure_has_zero_counts() {
        let s = Structure::default();
        assert_eq!(s.atom_count(), 0);
        assert_eq!(s.residue_count(), 0);
    }

    #[test]
    fn nested_count_walks_every_chain() {
        let s = Structure {
            id: "fake".into(),
            chains: vec![
                Chain {
                    id: 'A',
                    residues: vec![Residue {
                        name: "ALA".into(),
                        seq_id: 1,
                        atoms: vec![Atom {
                            name: "CA".into(),
                            element: "C".into(),
                            position: Vector3::new(0.0, 0.0, 0.0),
                            b_factor: 0.0,
                        }],
                    }],
                },
                Chain {
                    id: 'B',
                    residues: vec![Residue {
                        name: "GLY".into(),
                        seq_id: 1,
                        atoms: vec![],
                    }],
                },
            ],
        };
        assert_eq!(s.atom_count(), 1);
        assert_eq!(s.residue_count(), 2);
    }
}
