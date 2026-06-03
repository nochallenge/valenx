//! Feature 1 — the CRISPR nuclease database.
//!
//! A guide-design workflow needs, for each nuclease, the parameters
//! that decide where it can cut and what the cut looks like:
//!
//! - the **PAM** motif and which side of the protospacer it sits on;
//! - the **guide length** (protospacer length);
//! - the **cut-site offset** — where, relative to the protospacer, the
//!   double-strand (or single-strand) break falls;
//! - the **end chemistry** — blunt (Cas9-class) or 5′-overhang
//!   staggered (Cas12a-class).
//!
//! This module ships a small table of the widely-used nucleases —
//! SpCas9, SpCas9-NG, SaCas9, Cas12a / Cpf1, Cas12f, Cas13 and xCas9
//! — with *representative published* parameters. It is not an
//! exhaustive catalogue of every engineered variant; the values are
//! the ones a design tool needs to scan a locus and report a cut.
//!
//! The PAM is exposed both as a [`Nuclease`] field and, via
//! [`Nuclease::pam_spec`], as a [`valenx_genomics`] `PamSpec` so the
//! guide-design module can hand it straight to that crate's scanner —
//! the PAM logic is reused, never re-implemented here.

use serde::{Deserialize, Serialize};
use valenx_genomics::crispr::guide::{PamSide, PamSpec};

/// The mechanistic class of a CRISPR effector.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum NucleaseClass {
    /// A type-II RNA-guided DNA nuclease (Cas9 family) — a blunt
    /// double-strand break ~3 bp 5′ of the PAM.
    Cas9,
    /// A type-V RNA-guided DNA nuclease (Cas12a / Cpf1, Cas12f) — a
    /// staggered double-strand break with a 5′ overhang, distal to a
    /// 5′ PAM.
    Cas12,
    /// A type-VI RNA-guided **RNA** nuclease (Cas13) — cleaves RNA,
    /// not DNA; used for transcript knockdown rather than genome
    /// editing.
    Cas13,
}

impl NucleaseClass {
    /// Human-readable name.
    pub fn name(self) -> &'static str {
        match self {
            NucleaseClass::Cas9 => "Cas9 (type II)",
            NucleaseClass::Cas12 => "Cas12 (type V)",
            NucleaseClass::Cas13 => "Cas13 (type VI)",
        }
    }

    /// `true` when this class edits DNA (Cas9 / Cas12). Cas13 targets
    /// RNA and returns `false`.
    pub fn edits_dna(self) -> bool {
        matches!(self, NucleaseClass::Cas9 | NucleaseClass::Cas12)
    }
}

/// The chemistry of the break a nuclease leaves.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum CutEnds {
    /// A blunt double-strand break (Cas9-class).
    Blunt,
    /// A staggered double-strand break with a 5′ overhang of the given
    /// length in nucleotides (Cas12a-class — typically 4–5 nt).
    Staggered {
        /// 5′-overhang length in nucleotides.
        overhang: usize,
    },
    /// RNA cleavage — not a DNA double-strand break (Cas13).
    RnaCleavage,
}

/// A stable identifier for a catalogued nuclease.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum NucleaseId {
    /// *Streptococcus pyogenes* Cas9 — the canonical `NGG`-PAM editor.
    SpCas9,
    /// SpCas9-NG — an engineered SpCas9 with a relaxed `NG` PAM.
    SpCas9Ng,
    /// *Staphylococcus aureus* Cas9 — compact, `NNGRRT` PAM.
    SaCas9,
    /// *Acidaminococcus / Lachnospiraceae* Cas12a (Cpf1) — `TTTV`
    /// 5′ PAM, staggered cut.
    Cas12a,
    /// Cas12f (Cas14) — an ultra-compact type-V nuclease, `TTTR` PAM.
    Cas12f,
    /// Cas13 — an RNA-targeting effector (transcript knockdown).
    Cas13,
    /// xCas9 3.7 — an evolved SpCas9 with a broadened `NG` / `GAA` /
    /// `GAT` PAM range.
    XCas9,
}

impl NucleaseId {
    /// Every catalogued nuclease, in a stable order.
    pub fn all() -> [NucleaseId; 7] {
        [
            NucleaseId::SpCas9,
            NucleaseId::SpCas9Ng,
            NucleaseId::SaCas9,
            NucleaseId::Cas12a,
            NucleaseId::Cas12f,
            NucleaseId::Cas13,
            NucleaseId::XCas9,
        ]
    }
}

/// A catalogued CRISPR nuclease with the parameters a design workflow
/// needs.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Nuclease {
    /// Stable identifier.
    pub id: NucleaseId,
    /// Display name.
    pub name: String,
    /// Mechanistic class.
    pub class: NucleaseClass,
    /// The PAM motif in IUPAC codes (e.g. `"NGG"`, `"TTTV"`).
    pub pam: String,
    /// Which side of the protospacer the PAM sits on.
    pub pam_three_prime: bool,
    /// The protospacer (guide) length in nucleotides.
    pub guide_len: usize,
    /// Signed cut-site offset, in base pairs, of the break relative to
    /// the **PAM-proximal** end of the protospacer.
    ///
    /// For a 3′-PAM Cas9 the blunt cut is ~3 bp *into* the protospacer
    /// from the PAM, so the offset is `-3` (counted from the PAM-side
    /// protospacer end, negative = toward the PAM-distal end). For a
    /// 5′-PAM Cas12a the staggered cut is ~18–23 nt *distal* to the
    /// PAM-side protospacer end, a positive offset.
    pub cut_offset: i32,
    /// The chemistry of the break.
    pub cut_ends: CutEnds,
    /// A one-line description of the variant and its niche.
    pub notes: String,
}

impl Nuclease {
    /// PAM length in nucleotides.
    pub fn pam_len(&self) -> usize {
        self.pam.len()
    }

    /// `true` when this nuclease edits DNA (Cas9 / Cas12 classes).
    pub fn edits_dna(&self) -> bool {
        self.class.edits_dna()
    }

    /// Builds a [`valenx_genomics`] `PamSpec` from this nuclease so the
    /// guide-design module can hand it directly to that crate's
    /// PAM-scanning code — reusing the scanner, not duplicating it.
    pub fn pam_spec(&self) -> PamSpec {
        PamSpec {
            motif: self.pam.clone(),
            protospacer_len: self.guide_len,
            side: if self.pam_three_prime {
                PamSide::ThreePrime
            } else {
                PamSide::FivePrime
            },
        }
    }
}

/// Looks up a catalogued [`Nuclease`] by [`NucleaseId`].
///
/// The parameters are representative published values; see the
/// per-entry `notes`. They are sufficient to scan a locus, place a
/// cut and describe the break — they are not a substitute for the
/// primary literature when an exotic engineered variant is in play.
pub fn nuclease(id: NucleaseId) -> Nuclease {
    match id {
        NucleaseId::SpCas9 => Nuclease {
            id,
            name: "SpCas9".to_string(),
            class: NucleaseClass::Cas9,
            pam: "NGG".to_string(),
            pam_three_prime: true,
            guide_len: 20,
            cut_offset: -3,
            cut_ends: CutEnds::Blunt,
            notes: "Canonical S. pyogenes Cas9; blunt DSB 3 bp 5' of \
                    the NGG PAM."
                .to_string(),
        },
        NucleaseId::SpCas9Ng => Nuclease {
            id,
            name: "SpCas9-NG".to_string(),
            class: NucleaseClass::Cas9,
            pam: "NG".to_string(),
            pam_three_prime: true,
            guide_len: 20,
            cut_offset: -3,
            cut_ends: CutEnds::Blunt,
            notes: "Engineered SpCas9 with a relaxed NG PAM; widens \
                    targetable sites at some cost to activity."
                .to_string(),
        },
        NucleaseId::SaCas9 => Nuclease {
            id,
            name: "SaCas9".to_string(),
            class: NucleaseClass::Cas9,
            pam: "NNGRRT".to_string(),
            pam_three_prime: true,
            guide_len: 21,
            cut_offset: -3,
            cut_ends: CutEnds::Blunt,
            notes: "Compact S. aureus Cas9 (fits a single AAV with a \
                    promoter); NNGRRT PAM, 21 nt guide."
                .to_string(),
        },
        NucleaseId::Cas12a => Nuclease {
            id,
            name: "Cas12a (Cpf1)".to_string(),
            class: NucleaseClass::Cas12,
            pam: "TTTV".to_string(),
            pam_three_prime: false,
            guide_len: 23,
            cut_offset: 18,
            cut_ends: CutEnds::Staggered { overhang: 5 },
            notes: "Type-V Cas12a; T-rich 5' PAM, staggered DSB with a \
                    ~5 nt 5' overhang distal to the protospacer."
                .to_string(),
        },
        NucleaseId::Cas12f => Nuclease {
            id,
            name: "Cas12f (Cas14)".to_string(),
            class: NucleaseClass::Cas12,
            pam: "TTTR".to_string(),
            pam_three_prime: false,
            guide_len: 20,
            cut_offset: 18,
            cut_ends: CutEnds::Staggered { overhang: 5 },
            notes: "Ultra-compact type-V nuclease (~400-700 aa); T-rich \
                    5' PAM, staggered cut. Small size aids delivery."
                .to_string(),
        },
        NucleaseId::Cas13 => Nuclease {
            id,
            name: "Cas13".to_string(),
            class: NucleaseClass::Cas13,
            pam: "".to_string(),
            pam_three_prime: true,
            guide_len: 28,
            cut_offset: 0,
            cut_ends: CutEnds::RnaCleavage,
            notes: "Type-VI RNA-targeting effector; no DNA PAM (a PFS \
                    flanking-sequence preference instead). Used for \
                    transcript knockdown, not genome editing."
                .to_string(),
        },
        NucleaseId::XCas9 => Nuclease {
            id,
            name: "xCas9 3.7".to_string(),
            class: NucleaseClass::Cas9,
            pam: "NG".to_string(),
            pam_three_prime: true,
            guide_len: 20,
            cut_offset: -3,
            cut_ends: CutEnds::Blunt,
            notes: "Phage-assisted-evolved SpCas9 with broadened PAM \
                    recognition (NG / GAA / GAT); modelled here with \
                    its dominant NG PAM."
                .to_string(),
        },
    }
}

/// The full nuclease database, in [`NucleaseId::all`] order.
pub fn all_nucleases() -> Vec<Nuclease> {
    NucleaseId::all().into_iter().map(nuclease).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn database_has_seven_entries() {
        let db = all_nucleases();
        assert_eq!(db.len(), 7);
        // Each entry's id round-trips through the lookup.
        for n in &db {
            assert_eq!(nuclease(n.id).id, n.id);
        }
    }

    #[test]
    fn spcas9_parameters() {
        let n = nuclease(NucleaseId::SpCas9);
        assert_eq!(n.pam, "NGG");
        assert!(n.pam_three_prime);
        assert_eq!(n.guide_len, 20);
        assert_eq!(n.cut_offset, -3);
        assert_eq!(n.cut_ends, CutEnds::Blunt);
        assert!(n.edits_dna());
    }

    #[test]
    fn cas12a_is_five_prime_pam_and_staggered() {
        let n = nuclease(NucleaseId::Cas12a);
        assert!(!n.pam_three_prime);
        assert_eq!(n.pam, "TTTV");
        assert!(matches!(n.cut_ends, CutEnds::Staggered { overhang: 5 }));
        assert_eq!(n.class, NucleaseClass::Cas12);
    }

    #[test]
    fn cas13_targets_rna_not_dna() {
        let n = nuclease(NucleaseId::Cas13);
        assert_eq!(n.class, NucleaseClass::Cas13);
        assert!(!n.edits_dna());
        assert_eq!(n.cut_ends, CutEnds::RnaCleavage);
        assert!(n.pam.is_empty(), "Cas13 has no DNA PAM");
    }

    #[test]
    fn pam_spec_bridges_to_genomics() {
        let n = nuclease(NucleaseId::SpCas9);
        let spec = n.pam_spec();
        assert_eq!(spec.motif, "NGG");
        assert_eq!(spec.protospacer_len, 20);
        assert_eq!(spec.side, PamSide::ThreePrime);

        let c = nuclease(NucleaseId::Cas12a);
        assert_eq!(c.pam_spec().side, PamSide::FivePrime);
    }
}
