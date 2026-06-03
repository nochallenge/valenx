//! SMILES notation — parsing and (canonical) writing.
//!
//! - [`parse_smiles`] turns a SMILES string into a [`Molecule`],
//!   filling valence-based implicit hydrogens.
//! - [`write_smiles`] emits a plain DFS SMILES string.
//! - [`write_canonical_smiles`] emits a deterministic canonical string
//!   (Morgan extended-connectivity ranking) so isomorphic graphs map to
//!   the same text — the basis for dedup keys and structure lookup.
//!
//! See the submodule docs for the supported grammar and the v1
//! simplifications.

pub mod parser;
pub mod writer;

pub use parser::{fill_implicit_hydrogens, parse_smiles};
pub use writer::{canonical_ranks, write_canonical_smiles, write_smiles};

use crate::molecule::Molecule;

/// Parse a SMILES string and immediately re-emit its canonical form —
/// a convenience for "normalise this SMILES" callers.
pub fn canonicalize(smiles: &str) -> crate::error::Result<String> {
    let mol = parse_smiles(smiles)?;
    Ok(write_canonical_smiles(&mol))
}

/// Round-trip check: does `smiles` parse, write and re-parse to a graph
/// with the same atom and bond counts? A cheap structural sanity test
/// used by the test-suite and exposed for callers validating input.
pub fn round_trips(smiles: &str) -> bool {
    let Ok(m) = parse_smiles(smiles) else {
        return false;
    };
    let w = write_smiles(&m);
    matches!(parse_smiles(&w), Ok(m2)
        if m2.atom_count() == m.atom_count() && m2.bond_count() == m.bond_count())
}

/// Re-export so `Molecule` is reachable as `smiles::Molecule` in docs.
pub type Mol = Molecule;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonicalize_round_trips() {
        let c = canonicalize("CCO").unwrap();
        let c2 = canonicalize(&c).unwrap();
        assert_eq!(c, c2);
    }

    #[test]
    fn round_trips_flag() {
        assert!(round_trips("c1ccccc1"));
        assert!(round_trips("CC(=O)O"));
        assert!(!round_trips("C1CC")); // unclosed ring
    }
}
