//! Feature 8 — functional-motif scaffold templates.
//!
//! Many synthetic-RNA designs do not start from a bare dot-bracket —
//! they start from a **known functional scaffold**: a hammerhead
//! ribozyme, a malachite-green aptamer, a stable tetraloop hairpin. The
//! designer keeps the conserved catalytic / binding core and customises
//! only the variable arms.
//!
//! This module ships a small built-in library of well-characterised RNA
//! scaffolds — each with its sequence, its secondary structure and a
//! description of its conserved core — and a [`motif_scaffold`] lookup
//! that returns one as a starting [`RnaDesign`].
//!
//! ## v1 scope
//!
//! The library is a representative set of textbook scaffolds, not an
//! exhaustive Rfam catalogue. Each scaffold is the published *minimal*
//! functional sequence; a real design would graft the user's aptamer
//! arms or substrate-recognition helices onto the conserved core — that
//! arm customisation is left to the optimiser (which can mutate the
//! non-core positions). The scaffold's structure is the published
//! consensus; folding it on the exact sequence is a separate check the
//! validation step performs.

use crate::design::{DesignKind, RnaDesign};
use crate::error::{Result, RnaDesignError};
use valenx_rnastruct::Structure;

/// The kind of functional motif a scaffold provides.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub enum ScaffoldKind {
    /// A small self-cleaving ribozyme.
    Ribozyme,
    /// A small-molecule-binding aptamer.
    Aptamer,
    /// A stable structural hairpin / stem-loop.
    Hairpin,
}

impl ScaffoldKind {
    /// Human-readable name.
    pub fn name(self) -> &'static str {
        match self {
            ScaffoldKind::Ribozyme => "ribozyme",
            ScaffoldKind::Aptamer => "aptamer",
            ScaffoldKind::Hairpin => "hairpin",
        }
    }
}

/// A built-in functional-RNA scaffold template.
#[derive(Clone, Debug, PartialEq)]
pub struct MotifScaffold {
    /// The scaffold's stable identifier (e.g. `"hammerhead_minimal"`).
    pub id: &'static str,
    /// A human-readable name.
    pub name: &'static str,
    /// What functional class the scaffold belongs to.
    pub kind: ScaffoldKind,
    /// The scaffold's RNA sequence (`A C G U`).
    pub sequence: &'static str,
    /// The scaffold's secondary structure, as a dot-bracket string.
    pub dot_bracket: &'static str,
    /// A short description of the conserved core and the variable arms.
    pub description: &'static str,
}

impl MotifScaffold {
    /// The scaffold length in nucleotides.
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// `true` if the scaffold has no nucleotides (never — the library
    /// entries are all non-empty).
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }
}

/// The built-in scaffold library — a representative set of textbook
/// functional-RNA scaffolds.
///
/// Sequences and structures are the published minimal functional forms;
/// the GNRA-tetraloop and hammerhead entries are the canonical
/// engineering scaffolds.
fn scaffold_library() -> &'static [MotifScaffold] {
    &[
        MotifScaffold {
            id: "gnra_tetraloop",
            name: "GNRA tetraloop hairpin",
            kind: ScaffoldKind::Hairpin,
            // A 5-bp stem closing a GAAA tetraloop — one of the most
            // stable, most-used structural hairpins in RNA engineering.
            sequence: "GGCACGAAAGUGCC",
            dot_bracket: "(((((....)))))",
            description: "A 5-bp GC stem closing an exceptionally stable GAAA \
                          tetraloop. The conserved core is the GNRA loop; the \
                          stem sequence / length is the customisable arm.",
        },
        MotifScaffold {
            id: "uucg_tetraloop",
            name: "UUCG tetraloop hairpin",
            kind: ScaffoldKind::Hairpin,
            sequence: "GGCACUUCGGUGCC",
            dot_bracket: "(((((....)))))",
            description: "A 5-bp stem closing a UUCG tetraloop — the most \
                          thermostable RNA tetraloop, a standard cap for \
                          engineered stems. The stem is the customisable arm.",
        },
        MotifScaffold {
            id: "hammerhead_minimal",
            name: "Minimal hammerhead ribozyme core",
            kind: ScaffoldKind::Ribozyme,
            // A minimal hammerhead-style three-way junction: the
            // conserved catalytic core (CUGANGA ... GAAA) flanked by
            // three short helices. This is a representative minimal
            // engineering scaffold, not a specific natural sequence.
            // The 33-nt sequence and the 33-char dot-bracket must agree
            // in length; the unpaired 3' tail carries 7 dots.
            sequence: "GGGCUGAUGAGGCCGAAAGGCCGAAACGGGCCC",
            dot_bracket: "(((..(((((....)))))....))).......",
            description: "A minimal hammerhead-style self-cleaving ribozyme: a \
                          conserved catalytic core at a three-helix junction. \
                          Helices I and III are the substrate-recognition arms \
                          the designer customises; the core stays fixed.",
        },
        MotifScaffold {
            id: "malachite_green_aptamer",
            name: "Malachite-green aptamer core",
            kind: ScaffoldKind::Aptamer,
            // A representative malachite-green-binding aptamer fold: an
            // internal-loop binding pocket between two short helices.
            sequence: "GGAUCCCGACUGGCGAGAGCCAGGUAACGAAUGGAUCC",
            dot_bracket: "((((((...((((....))))...))))))........",
            description: "An aptamer fold presenting a small-molecule binding \
                          pocket in an asymmetric internal loop between two \
                          helices. The pocket is the conserved core; the \
                          flanking helices are customisable.",
        },
        MotifScaffold {
            id: "theophylline_aptamer",
            name: "Theophylline aptamer core",
            kind: ScaffoldKind::Aptamer,
            // The 33-nt sequence and the 33-char dot-bracket must agree
            // in length; the unpaired 3' tail carries 7 dots.
            sequence: "GGCGAUACCAGCCGAAAGGCCCUUGGCAGCGUC",
            dot_bracket: "((((.((((....)))).....)))).......",
            description: "The widely-used theophylline-binding aptamer — a \
                          high-affinity, high-specificity pocket popular in \
                          synthetic riboswitches. The pocket core is fixed.",
        },
    ]
}

/// Looks up a built-in scaffold by its id (case-insensitive).
pub fn scaffold_by_id(id: &str) -> Option<&'static MotifScaffold> {
    scaffold_library()
        .iter()
        .find(|s| s.id.eq_ignore_ascii_case(id))
}

/// All scaffolds of a given functional kind.
pub fn scaffolds_of_kind(kind: ScaffoldKind) -> Vec<&'static MotifScaffold> {
    scaffold_library().iter().filter(|s| s.kind == kind).collect()
}

/// Returns a built-in functional-motif scaffold as a starting
/// [`RnaDesign`] (feature 8).
///
/// The returned design's sequence is the scaffold's published sequence
/// and its [`DesignKind`] records the scaffold it started from. The
/// optimiser can then customise the variable (non-core) arms.
///
/// # Errors
/// [`RnaDesignError::NoDesign`] if `id` matches no library scaffold.
pub fn motif_scaffold(id: &str) -> Result<RnaDesign> {
    let scaffold = scaffold_by_id(id).ok_or_else(|| {
        RnaDesignError::no_design(
            "motif",
            format!("no built-in scaffold with id `{id}`"),
        )
    })?;

    // Confirm the scaffold's shipped consensus dot-bracket is parseable
    // (every library entry is, but check rather than assume).
    Structure::from_dot_bracket(scaffold.dot_bracket)
        .map_err(|e| RnaDesignError::no_design("motif", e.to_string()))?;

    let notes = vec![
        format!(
            "Seeded from the built-in {} scaffold `{}` ({}).",
            scaffold.kind.name(),
            scaffold.id,
            scaffold.name,
        ),
        scaffold.description.to_string(),
        "The conserved functional core must be kept fixed; the optimiser may \
         customise only the variable arm positions."
            .to_string(),
        "A scaffold is a published consensus — confirm the exact sequence folds \
         to the functional structure and retains activity in the wet lab."
            .to_string(),
    ];

    Ok(RnaDesign {
        sequence: scaffold.sequence.as_bytes().to_vec(),
        kind: DesignKind::MotifScaffold {
            scaffold: scaffold.id.to_string(),
        },
        cds_span: None,
        construct: None,
        notes,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn library_is_populated() {
        assert!(scaffold_library().len() >= 5);
        assert!(scaffold_by_id("gnra_tetraloop").is_some());
        assert!(scaffold_by_id("GNRA_TETRALOOP").is_some()); // case-insensitive
        assert!(scaffold_by_id("nonexistent").is_none());
    }

    #[test]
    fn every_scaffold_is_self_consistent() {
        // Each scaffold's sequence length must equal its dot-bracket
        // length, and the dot-bracket must parse.
        for s in scaffold_library() {
            assert_eq!(
                s.sequence.len(),
                s.dot_bracket.len(),
                "scaffold {} length mismatch",
                s.id
            );
            let st = Structure::from_dot_bracket(s.dot_bracket)
                .unwrap_or_else(|_| panic!("scaffold {} dot-bracket invalid", s.id));
            // Scaffolds are pseudoknot-free.
            assert!(!st.has_pseudoknot(), "scaffold {} is pseudoknotted", s.id);
            // The sequence is valid RNA.
            assert!(
                valenx_rnastruct::RnaSeq::parse(s.sequence).is_ok(),
                "scaffold {} sequence is not valid RNA",
                s.id
            );
        }
    }

    #[test]
    fn scaffolds_of_kind_filters() {
        let aptamers = scaffolds_of_kind(ScaffoldKind::Aptamer);
        assert!(!aptamers.is_empty());
        assert!(aptamers.iter().all(|s| s.kind == ScaffoldKind::Aptamer));
        let ribozymes = scaffolds_of_kind(ScaffoldKind::Ribozyme);
        assert!(!ribozymes.is_empty());
    }

    #[test]
    fn motif_scaffold_returns_a_design() {
        let d = motif_scaffold("gnra_tetraloop").unwrap();
        assert_eq!(d.sequence_str(), "GGCACGAAAGUGCC");
        assert!(matches!(d.kind, DesignKind::MotifScaffold { .. }));
        assert!(!d.is_coding());
        assert!(!d.notes.is_empty());
    }

    #[test]
    fn motif_scaffold_rejects_unknown_id() {
        let err = motif_scaffold("not_a_scaffold").unwrap_err();
        assert_eq!(err.code(), "rnadesign.no_design");
    }

    #[test]
    fn hammerhead_is_a_ribozyme() {
        let s = scaffold_by_id("hammerhead_minimal").unwrap();
        assert_eq!(s.kind, ScaffoldKind::Ribozyme);
    }
}

