//! Chemical reactions — model, transforms, standardisation, tautomers.
//!
//! - [`model`] — the [`Reaction`] data model and reaction-SMILES I/O;
//! - [`transform`] — a SMIRKS-class [`ReactionTransform`] applied to a
//!   molecule, and combinatorial [`enumerate_library`];
//! - [`standardize`](mod@standardize) — salt stripping,
//!   neutralisation, the canonical parent form;
//! - [`tautomer`] — tautomer enumeration and canonical tautomer.

pub mod model;
pub mod standardize;
pub mod tautomer;
pub mod transform;

pub use model::{parse_reaction_smiles, write_reaction_smiles, Reaction};
pub use standardize::{neutralize, standardize, strip_salts};
pub use tautomer::{canonical_tautomer, enumerate_tautomers, tautomer_count, tautomer_score};
pub use transform::{enumerate_library, ReactionTransform};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol_from_smiles;

    #[test]
    fn module_surface_smoke() {
        // a reaction parses, a transform applies, a molecule standardises
        let rxn = parse_reaction_smiles("C=C>>CC").unwrap();
        assert_eq!(rxn.reactants.len(), 1);

        let t = ReactionTransform::parse("[C:1]=[C:2]>>[C:1]-[C:2]").unwrap();
        let prod = t.apply(&mol_from_smiles("C=C").unwrap());
        assert_eq!(prod.len(), 1);

        let parent = strip_salts(&mol_from_smiles("CCO.Cl").unwrap());
        assert_eq!(parent.component_count(), 1);
    }
}
