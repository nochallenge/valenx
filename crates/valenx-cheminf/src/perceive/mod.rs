//! Molecular perception — deriving the chemistry implied by the graph.
//!
//! Five perceivers, each usable on its own:
//!
//! - [`formula`] — molecular formula + exact / average molecular
//!   weight;
//! - [`rings`] — the smallest set of smallest rings (SSSR);
//! - [`aromatic`] — Hückel 4n+2 aromaticity;
//! - [`hydrogen`] — add / remove explicit hydrogens, recompute the
//!   implicit count;
//! - [`stereo`] — CIP tetrahedral R/S and double-bond E/Z assignment.
//!
//! [`perceive_all`] runs the ring → aromaticity → stereo pipeline in
//! the right order on a molecule.

pub mod aromatic;
pub mod formula;
pub mod hydrogen;
pub mod rings;
pub mod stereo;

pub use aromatic::perceive_aromaticity;
pub use formula::{average_molecular_weight, molecular_formula, monoisotopic_mass};
pub use hydrogen::{add_explicit_hydrogens, remove_explicit_hydrogens};
pub use rings::{sssr, Ring, RingInfo};
pub use stereo::{assign_cip_labels, assign_double_bond_stereo};

use crate::molecule::Molecule;

/// Run the full perception pipeline on `mol`:
///
/// 1. recompute implicit hydrogens (so valences are current),
/// 2. SSSR ring perception,
/// 3. Hückel aromaticity (needs the rings),
/// 4. CIP R/S and double-bond E/Z stereo assignment.
///
/// Idempotent — running it twice gives the same molecule.
pub fn perceive_all(mol: &mut Molecule) {
    hydrogen::recompute_implicit_hydrogens(mol);
    let ring_info = rings::sssr(mol);
    aromatic::perceive_aromaticity(mol, &ring_info);
    stereo::assign_cip_labels(mol);
    stereo::assign_double_bond_stereo(mol);
}

/// A bundled perception result for callers that want the ring info kept
/// around (descriptors, scaffolds) rather than recomputing it.
#[derive(Clone, Debug)]
pub struct Perception {
    /// The SSSR ring analysis.
    pub rings: RingInfo,
}

/// Perceive and return the [`Perception`] bundle (mutates `mol` exactly
/// like [`perceive_all`]).
pub fn perceive(mol: &mut Molecule) -> Perception {
    hydrogen::recompute_implicit_hydrogens(mol);
    let rings = rings::sssr(mol);
    aromatic::perceive_aromaticity(mol, &rings);
    stereo::assign_cip_labels(mol);
    stereo::assign_double_bond_stereo(mol);
    // aromaticity may have changed ring membership? no — SSSR is graph-
    // only, recompute is unnecessary; return the rings we built.
    Perception { rings }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::parse_smiles;

    #[test]
    fn perceive_all_is_idempotent() {
        let mut m = parse_smiles("c1ccccc1").unwrap();
        perceive_all(&mut m);
        let snapshot = m.clone();
        perceive_all(&mut m);
        assert_eq!(m, snapshot, "perception must be idempotent");
    }

    #[test]
    fn perceive_bundle_has_rings() {
        let mut m = parse_smiles("c1ccccc1").unwrap();
        let p = perceive(&mut m);
        assert_eq!(p.rings.ring_count(), 1);
        assert!(m.atoms.iter().all(|a| a.aromatic));
    }
}
