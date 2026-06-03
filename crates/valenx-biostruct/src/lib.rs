//! # valenx-biostruct — macromolecular structure analysis
//!
//! Round 6 Block 8 of the Valenx roadmap. A native-Rust replacement
//! for Biopython.PDB and the structure-handling portions of ChimeraX,
//! PyMOL, VMD, 3DNA / X3DNA, Curves+, DSSR and FoldSeek — pure
//! parsing and geometry, no neural-network weights and no external
//! processes.
//!
//! ## What it does
//!
//! - **Model & I/O** — the [`structure`] hierarchy ([`Structure`] >
//!   [`Model`] > [`Chain`] > [`Residue`] > [`Atom`]); PDB and mmCIF
//!   readers and writers in [`io`]; a [`select`]ion mini-language.
//! - **Geometry** ([`geometry`]) — inter-atomic distances with a
//!   grid-accelerated neighbour search, bond / dihedral / backbone /
//!   sidechain torsions, Ramachandran classification, residue contact
//!   maps, steric-clash detection, radius of gyration / principal
//!   axes, and Shrake-Rupley SASA.
//! - **Secondary structure** ([`dssp`]) — full **Kabsch-Sander 1983**
//!   `H G I E B T S` assignment with the published H > G > I
//!   tie-breaking, parallel / antiparallel β-bridge perception with
//!   strand-extension ladders, n-turn interior labelling for `T`,
//!   bend (`S`) curvature cut-off, helix / sheet enumeration with
//!   sheet topology, and per-state residue counts for benchmarking.
//! - **Superposition & comparison** — [`superpose`] (Kabsch and
//!   quaternion optimal superposition, RMSD) and [`compare`] (the
//!   sequence-anchored iterative aligner [`compare::align`] as a v1
//!   path; **the TM-align-class iterative-DP aligner**
//!   [`compare::tmalign`] as the production default — sequence-
//!   independent, with SS / fragment / diagonal seeding and a CE-style
//!   aligned-fragment-pair variant; TM-score; FoldSeek-class
//!   structural alphabet).
//! - **Nucleic acids** ([`nucleic`]) — Watson-Crick / wobble
//!   base-pair detection, 3DNA-class base-pair and step parameters,
//!   groove widths, a single straight TLS helical axis
//!   ([`nucleic::fit_helical_axis`]), and a Curves+-class **curved**
//!   spline axis with per-bp curvature analysis
//!   ([`nucleic::fit_curved_axis`]) — natural cubic spline through
//!   the per-bp screw-axis foot-points.
//! - **Assemblies, edits, validation** — [`assembly`] (biological-
//!   assembly generation from BIOMT operators), [`edit`] (B-factor
//!   analysis, interface detection, point-mutation modelling),
//!   [`validate`] (a [`validate::ValidationReport`]), and [`analyze`]
//!   (a per-structure [`analyze::StructureReport`] bundling chain
//!   summaries, composition, secondary-structure content and Rg, plus
//!   multi-structure batch utilities).
//!
//! ## Errors
//!
//! Every fallible function returns
//! [`Result<_, BiostructError>`](error::BiostructError). The error
//! type carries stable [`code`](error::BiostructError::code) and
//! [`category`](error::BiostructError::category) accessors for
//! telemetry.
//!
//! ## Scope
//!
//! Most modules ship a real working implementation of the published
//! algorithm. Two production-depth structure-alignment / DNA-helical
//! upgrades are now the defaults:
//!
//! - The structure aligner ships **both** the v1 sequence-anchored
//!   iterative aligner ([`compare::align::align_chains`]) **and** a
//!   TM-align-class DP-refined sequence-independent aligner
//!   ([`compare::tmalign::align_chains`]) — pick the latter for
//!   sequence-divergent / structurally-similar pairs. A CE-style
//!   AFP combinatorial variant ([`compare::tmalign::align_chains_ce`])
//!   is the second method.
//! - The helical-axis fit ships **both** the single straight axis and
//!   a Curves+-class curved spline axis with per-bp curvature
//!   ([`nucleic::fit_curved_axis`]) — pick the curved one for bent
//!   duplexes.
//!
//! Documented remaining simplifications:
//! - The FoldSeek-class structural alphabet is a hand-designed
//!   quantiser in the *spirit* of the learnt 3Di alphabet, not the
//!   trained VQ-VAE.
//! - The 3DNA-class base reference frame is least-squares-fitted from
//!   ring atoms rather than the Olson-convention standard-base
//!   template.
//! - Mutation modelling places one idealised sidechain template
//!   without rotamer search or energy relaxation.
//! - The DSSP H-bond model + state set is the published Kabsch-Sander
//!   1983 algorithm; the **only** documented divergence from a
//!   reference DSSP run is the use of geometric chain-break detection
//!   (Cα–Cα distance jump) rather than `SEQRES` gap parsing.

#![forbid(unsafe_code)]

pub mod analyze;
pub mod assembly;
pub mod compare;
pub mod dssp;
pub mod edit;
pub mod error;
pub mod geometry;
pub mod io;
pub mod nucleic;
pub mod select;
pub mod structure;
pub mod superpose;
pub mod validate;

// --- Convenience re-exports of the most-used types --------------------

pub use analyze::{analyze_batch, StructureReport};
pub use compare::{
    align_chains_ce, align_chains_tm, StructureAlignment, TmAlignment, TmSeedKind,
};
pub use error::{BiostructError, ErrorCategory, Result};
pub use io::{read_structure, write_mmcif, write_pdb};
pub use nucleic::{fit_curved_axis, fit_helical_axis, CurvedHelicalAxis, HelicalAxis};
pub use select::Selection;
pub use structure::{Atom, Chain, Model, Residue, ResidueKind, Structure};
pub use superpose::{kabsch, rmsd, Superposition};
pub use validate::{validate_structure, ValidationReport};

/// Parse a coordinate file (PDB or mmCIF, auto-detected) and validate
/// the resulting hierarchy. The one-call entry point most callers
/// want.
pub fn load_structure(text: &str, id: &str) -> Result<Structure> {
    let structure = read_structure(text, id)?;
    structure.validate()?;
    Ok(structure)
}

#[cfg(test)]
mod tests {
    use super::*;

    const MINI_PDB: &str = "\
HEADER    TEST                                    01-JAN-00   1ABC
ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00 20.00           N
ATOM      2  CA  ALA A   1      11.639   6.071  -5.147  1.00 21.00           C
ATOM      3  C   ALA A   1      13.085   6.564  -5.180  1.00 19.50           C
ATOM      4  O   ALA A   1      13.339   7.742  -5.430  1.00 18.00           O
ATOM      5  CB  ALA A   1      10.768   6.882  -4.205  1.00 22.00           C
END
";

    #[test]
    fn end_to_end_load() {
        let s = load_structure(MINI_PDB, "fallback").expect("structure loads");
        assert_eq!(s.id, "1ABC");
        assert_eq!(s.atom_count(), 5);
        assert_eq!(
            s.first_model().chain("A").unwrap().observed_sequence(),
            "A"
        );
    }

    #[test]
    fn crate_error_reexported() {
        let e = BiostructError::not_yet("x");
        assert_eq!(e.category_enum(), ErrorCategory::Capability);
    }

    #[test]
    fn selection_reexported() {
        let sel = Selection::parse("name CA").unwrap();
        let s = load_structure(MINI_PDB, "x").unwrap();
        assert_eq!(sel.select(s.first_model()).len(), 1);
    }

    #[test]
    fn validation_reexported() {
        let s = load_structure(MINI_PDB, "x").unwrap();
        let report = validate_structure(&s, 0.4).unwrap();
        assert_eq!(report.structure_id, "1ABC");
    }

    #[test]
    fn structure_report_reexported() {
        // The batch per-structure analysis bundle.
        let s = load_structure(MINI_PDB, "x").unwrap();
        let report = StructureReport::analyze(&s, 0.4).unwrap();
        assert_eq!(report.id, "1ABC");
        assert_eq!(report.chains.len(), 1);
        let batch = analyze_batch(std::slice::from_ref(&s), 0.4).unwrap();
        assert_eq!(batch.len(), 1);
    }
}
