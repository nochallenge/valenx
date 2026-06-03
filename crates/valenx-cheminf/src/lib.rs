//! # valenx-cheminf — cheminformatics core
//!
//! A native-Rust replacement for
//! the classical core of RDKit, Open Babel, the non-neural portion of
//! DeepChem and datamol — pure graph and numerical algorithms, no
//! neural-network weights and no external processes.
//!
//! ## What it does
//!
//! - **Model & notation** — the [`molecule`] graph ([`Molecule`] of
//!   [`Atom`]s and [`Bond`]s); a
//!   [`smiles`] parser + canonical writer; a [`smarts`] query parser
//!   with VF2 substructure matching; a [`molfile`] MOL/SDF reader and
//!   writer; an [`inchi`] layered canonical string.
//! - **Perception** ([`perceive`]) — molecular formula and exact /
//!   average weight, smallest-set-of-smallest-rings, Hückel 4n+2
//!   aromaticity, valence-based hydrogen handling, CIP R/S and E/Z.
//! - **Coordinates & charges** — 2D depiction and distance-geometry 3D
//!   conformers in [`coords`], a force-field cleanup in
//!   [`forcefield`], Gasteiger PEOE charges in [`charge`].
//! - **Fingerprints** ([`fingerprint`]) — Morgan / ECFP, MACCS-class
//!   keys, topological path fingerprints, Tanimoto / Dice similarity.
//! - **Descriptors & scaffolds** — [`descriptors`] (logP, TPSA, HBD /
//!   HBA, rotatable bonds, Lipinski / Veber) and [`scaffold`]
//!   (Bemis-Murcko framework, maximum common substructure).
//! - **Reactions** ([`reaction`]) — a reaction model, a SMIRKS-class
//!   transform, combinatorial enumeration, standardization, tautomers.
//! - **Analysis** ([`analyze`]) — a pharmacophore feature model, a QED
//!   drug-likeness score, and a batch [`MoleculeReport`](analyze::MoleculeReport).
//!
//! ## Errors
//!
//! Every fallible function returns
//! [`Result<_, CheminfError>`](error::CheminfError). The error type
//! carries stable [`code`](error::CheminfError::code) and
//! [`category`](error::CheminfError::category) accessors for telemetry.
//!
//! ## v1 scope
//!
//! This is a real working v1, not production parity with RDKit. Each
//! module documents its own simplifications inline; the notable ones:
//! SMILES stereo is plain `@` / `@@` and `/` `\` (no extended `@TH` /
//! `@AL` / `@OH` classes); the InChI-class string is a layered
//! canonical identifier in the *spirit* of InChI, not the official
//! IUPAC algorithm; MACCS keys are a ~160-key SMARTS-defined set, not
//! the exact 166 of the MDL definition; QED uses the published
//! Bickerton desirability functions.
//!
//! ## Commercial-depth upgrades
//!
//! Three commercial-depth modules upgrade the most-used cheminformatics
//! primitives beyond the original v1 reduced-term implementations:
//!
//! - [`forcefield_mmff94`] — production MMFF94 with tabulated per-
//!   atom-type parameters from Halgren 1996 (bond / angle / stretch-
//!   bend / torsion / buffered-14-7 vdW / Coulomb), an MMFF94 atom
//!   typer over a representative subset of the published ~95 types,
//!   and an analytic gradient + steepest-descent minimiser. Replaces
//!   the original reduced FF as the default conformer-cleanup path.
//! - [`coords::etkdg`] — Experimental-Torsion-knowledge Distance
//!   Geometry (Riniker & Landrum 2015): a Gaussian-mixture
//!   torsion-preference library biases the initial dihedrals,
//!   followed by MMFF94 cleanup and multi-conformer generation with
//!   RMSD pruning.
//! - [`reaction::tautomer`] — 1,3 + 1,5 (vinylogous) shifts including
//!   aromatic ring-chain (lactam ↔ lactim), with a published-class
//!   InChI-style canonical-tautomer scoring rubric.

#![forbid(unsafe_code)]
// Surface future `&str` byte-offset slicing in clippy review — this
// crate parses untrusted text (SMILES/InChI/MOL), where
// non-char-boundary slices panic. WARN (not deny): most existing slices
// are safe ASCII; this only flags NEW ones.
#![allow(clippy::string_slice, reason = "parsers slice ASCII fixed-format records at byte offsets from find() or constant ASCII prefixes, always valid char boundaries")]

pub mod analyze;
pub mod charge;
pub mod coords;
pub mod descriptors;
pub mod element;
pub mod error;
pub mod fingerprint;
pub mod forcefield;
pub mod forcefield_mmff94;
pub mod inchi;
pub mod molecule;
pub mod molfile;
pub mod perceive;
pub mod reaction;
pub mod scaffold;
pub mod smarts;
pub mod smiles;

// --- Convenience re-exports of the most-used types --------------------

pub use error::{CheminfError, ErrorCategory, Result};
pub use molecule::{Atom, Bond, BondOrder, BondStereo, Chirality, Molecule};
pub use smiles::{parse_smiles, write_canonical_smiles, write_smiles};

/// Parse a SMILES string, run full perception (rings, aromaticity,
/// implicit hydrogens already filled by the parser) and return a
/// ready-to-analyse [`Molecule`]. The one-call entry point most
/// callers want.
pub fn mol_from_smiles(smiles: &str) -> Result<Molecule> {
    let mut mol = parse_smiles(smiles)?;
    perceive::perceive_all(&mut mol);
    Ok(mol)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn end_to_end_smoke() {
        let m = mol_from_smiles("c1ccccc1").expect("benzene parses");
        assert_eq!(m.atom_count(), 6);
        assert_eq!(m.heavy_atom_count(), 6);
    }

    #[test]
    fn crate_error_reexported() {
        let e = CheminfError::not_yet("x");
        assert_eq!(e.category_enum(), ErrorCategory::Capability);
    }
}
