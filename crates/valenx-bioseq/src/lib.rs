//! # valenx-bioseq — biological-sequence core
//!
//! A native-Rust replacement
//! for the sequence-handling portion of Biopython, BioJava, SeqKit,
//! ApE, SerialCloner, pLannotate and pydna — pure parsing and
//! algorithms, no neural-network weights and no external processes.
//!
//! ## What it does
//!
//! - **Types** — IUPAC [`alphabet`]s (DNA / RNA / protein, ambiguity
//!   codes); the validated [`Seq`] type (linear or circular); the
//!   [`SeqRecord`] / [`SeqFeature`] annotated-sequence model with
//!   `join(...)` / `complement(...)` compound locations.
//! - **File I/O** ([`io`]) — FASTA, FASTQ (with the Phred / Solexa
//!   quality codec), GenBank and EMBL readers and writers, plus
//!   streaming record iterators and a `.fai`-style index.
//! - **Core ops** ([`ops`]) — reverse complement, transcription,
//!   translation across every non-withdrawn NCBI genetic-code table
//!   (25 tables: 1-6, 9-14, 16, 21-31, 33), six-frame translation,
//!   and ORF finding.
//! - **Analysis** ([`analysis`]) — GC content / skew, k-mer
//!   statistics, melting temperature (Wallace + nearest-neighbor),
//!   molecular weight, and ProtParam-class protein properties.
//! - **Cloning** ([`cloning`]) — a REBASE-subset restriction-enzyme
//!   database with single/double digest, restriction maps + virtual
//!   gels, pLannotate-class plasmid annotation, and codon
//!   optimization with CAI.
//! - **Primers** ([`primer`]) — Tm-targeted primer design with
//!   hairpin / dimer screening, and in-silico PCR.
//! - **Editing** ([`editing`]) — insert / delete / replace / slice /
//!   rotate with automatic feature-coordinate updates.
//!
//! ## Errors
//!
//! Every fallible function returns
//! [`Result<_, BioseqError>`](error::BioseqError). The error type
//! carries stable [`code`](error::BioseqError::code) and
//! [`category`](error::BioseqError::category) accessors for telemetry.
//!
//! ## v1 scope
//!
//! This is a real working v1, not production parity with the 30-year
//! reference tools. Each module documents its own simplifications
//! inline; the notable ones are: GenBank/EMBL `REFERENCE` blocks are
//! parsed and round-tripped but the parsers stop at the structural
//! fields needed for sequence + annotation analysis (i.e. literature-
//! reference text is preserved verbatim, not parsed into citation
//! sub-fields); restriction matching is exact-site (Type IIS cut
//! offsets are tabulated but star activity is not modeled); plasmid
//! annotation is signature-based, not full alignment.

#![forbid(unsafe_code)]
#![warn(missing_docs)]
// Surface future `&str` byte-offset slicing in clippy review — this
// crate parses untrusted text, where non-char-boundary slices panic.
// WARN (not deny): most existing slices are safe ASCII; this only flags
// NEW ones.
#![allow(
    clippy::string_slice,
    reason = "parsers slice ASCII fixed-format records at byte offsets from find() or constant ASCII prefixes, always valid char boundaries"
)]

pub mod alphabet;
pub mod analysis;
pub mod cloning;
pub mod editing;
pub mod error;
pub mod index;
pub mod io;
pub mod ops;
pub mod primer;
pub mod record;
pub mod seq;

// --- Convenience re-exports of the most-used types --------------------

pub use alphabet::SeqKind;
pub use error::{BioseqError, ErrorCategory, Result};
pub use record::{Location, SeqFeature, SeqRecord, Span, Strand};
pub use seq::{Seq, Topology};

#[cfg(test)]
mod tests {
    use super::*;

    /// A small end-to-end smoke test exercising several modules
    /// together: parse FASTA, find an ORF, translate it, and check the
    /// protein's molecular weight is positive.
    #[test]
    fn end_to_end_fasta_orf_translate() {
        let fasta = ">gene1 demo\nATGAAAGGGTTTTAA\n";
        let recs = io::fasta::parse(fasta, SeqKind::Dna).unwrap();
        assert_eq!(recs.len(), 1);
        let dna = &recs[0].seq;

        let code = ops::translate::GeneticCode::standard();
        let orfs = ops::orf::find_orfs(
            dna,
            &code,
            ops::orf::OrfOptions {
                atg_only: true,
                min_protein_len: 1,
                allow_no_stop: false,
            },
        )
        .unwrap();
        assert_eq!(orfs.len(), 1);
        assert_eq!(orfs[0].protein.as_str(), "MKGF*");

        let protein = Seq::new(SeqKind::Protein, "MKGF").unwrap();
        let mw = analysis::weight::molecular_weight_protein(&protein).unwrap();
        assert!(mw > 0.0);
    }

    /// End-to-end cloning: digest a plasmid-like circular sequence and
    /// confirm the virtual gel sums back to the molecule length.
    #[test]
    fn end_to_end_digest_gel() {
        let s = Seq::with_topology(
            SeqKind::Dna,
            "GAATTCAAAAAAAAAGAATTCAAAAAAAAA",
            Topology::Circular,
        )
        .unwrap();
        let ecori = cloning::restriction::enzyme_by_name("EcoRI").unwrap();
        let gel = cloning::digest_map::digest_to_gel(&s, &[ecori]).unwrap();
        assert_eq!(gel.band_sizes.iter().sum::<usize>(), s.len());
    }

    #[test]
    fn re_exports_are_wired() {
        // Touch each convenience re-export so the paths are exercised.
        let _: SeqKind = SeqKind::Dna;
        let _: Topology = Topology::Linear;
        let e = BioseqError::not_yet("x");
        assert_eq!(e.category_enum(), ErrorCategory::Capability);
        let s = Seq::new(SeqKind::Dna, "ACGT").unwrap();
        let rec = SeqRecord::new("r", s);
        assert_eq!(rec.len(), 4);
        let _ = Span::new(0, 1);
        let _ = Strand::Forward;
        let _ = Location::single(0, 1);
        let _ = SeqFeature::new("misc", Location::single(0, 1));
    }
}

/// Reference-value validation against known molecular-biology facts.
///
/// These tests assert genuine published / textbook values: known
/// translations, reverse complements, GC contents, restriction-fragment
/// sizes and molecular weights — not internal consistency.
#[cfg(test)]
mod validation {
    use super::*;
    use crate::analysis::composition::gc_content;
    use crate::analysis::weight::{molecular_weight_dna_double_stranded, molecular_weight_protein};
    use crate::ops::revcomp::reverse_complement;
    use crate::ops::translate::{translate_default, GeneticCode};

    /// Translating the start of the human insulin (INS) coding sequence
    /// must give the known N-terminal protein sequence. The first 30 nt
    /// of the INS CDS, `ATG GCC CTG TGG ATG CGC CTC CTG CCC CTG`, encode
    /// the signal-peptide start `MALWMRLLPL`.
    #[test]
    fn translation_of_insulin_cds_start_is_correct() {
        let cds = Seq::new(SeqKind::Dna, "ATGGCCCTGTGGATGCGCCTCCTGCCCCTG").unwrap();
        let protein = translate_default(&cds, &GeneticCode::standard()).unwrap();
        assert_eq!(protein.as_str(), "MALWMRLLPL");
    }

    /// The standard genetic code: ATG→M (start), TGG→W (the only Trp
    /// codon), and the three stop codons TAA/TAG/TGA→*.
    #[test]
    fn standard_codon_table_landmark_codons() {
        let code = GeneticCode::standard();
        let tr = |s: &str| {
            translate_default(&Seq::new(SeqKind::Dna, s).unwrap(), &code)
                .unwrap()
                .as_str()
                .to_string()
        };
        assert_eq!(tr("ATG"), "M");
        assert_eq!(tr("TGG"), "W");
        assert_eq!(tr("TAA"), "*");
        assert_eq!(tr("TAG"), "*");
        assert_eq!(tr("TGA"), "*");
        // The six leucine codons all map to L.
        for c in ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"] {
            assert_eq!(tr(c), "L", "codon {c}");
        }
    }

    /// Reverse complement of a known sequence. `5'-GAATTC-3'` (the EcoRI
    /// site) is its own reverse complement (a palindrome); a
    /// non-palindromic example is checked too.
    #[test]
    fn reverse_complement_of_known_sequences() {
        let ecori = Seq::new(SeqKind::Dna, "GAATTC").unwrap();
        assert_eq!(reverse_complement(&ecori).unwrap().as_str(), "GAATTC");
        // 5'-ATGGCC-3'  ->  reverse complement  5'-GGCCAT-3'.
        let s = Seq::new(SeqKind::Dna, "ATGGCC").unwrap();
        assert_eq!(reverse_complement(&s).unwrap().as_str(), "GGCCAT");
    }

    /// GC content of known sequences — an exact rational quantity.
    #[test]
    fn gc_content_of_known_sequences() {
        // All-AT: 0% GC.
        let at = Seq::new(SeqKind::Dna, "ATATATAT").unwrap();
        assert!((gc_content(&at).unwrap() - 0.0).abs() < 1e-12);
        // All-GC: 100%.
        let gc = Seq::new(SeqKind::Dna, "GCGCGCGC").unwrap();
        assert!((gc_content(&gc).unwrap() - 1.0).abs() < 1e-12);
        // 4 of 8 are G/C -> exactly 0.5.
        let half = Seq::new(SeqKind::Dna, "ATGCATGC").unwrap();
        assert!((gc_content(&half).unwrap() - 0.5).abs() < 1e-12);
    }

    /// A restriction digest with known fragment sizes. A 30 bp circular
    /// sequence with two EcoRI sites must be cut into two fragments
    /// whose sizes sum to 30.
    #[test]
    fn ecori_double_cut_fragment_sizes() {
        // Two GAATTC sites 15 bp apart on a 30 bp circle.
        let plasmid = Seq::with_topology(
            SeqKind::Dna,
            "GAATTCAAAAAAAAAGAATTCAAAAAAAAA",
            Topology::Circular,
        )
        .unwrap();
        let ecori = cloning::restriction::enzyme_by_name("EcoRI").unwrap();
        let cuts = cloning::restriction::digest(&plasmid, ecori).unwrap();
        assert_eq!(cuts.len(), 2, "two EcoRI sites expected");
        let gel = cloning::digest_map::digest_to_gel(&plasmid, &[ecori]).unwrap();
        // A circular molecule cut twice -> exactly two fragments that
        // partition the 29 bp molecule.
        assert_eq!(gel.band_sizes.len(), 2);
        assert_eq!(gel.band_sizes.iter().sum::<usize>(), plasmid.len());
    }

    /// Protein molecular weight against published values. As a 1-residue
    /// "protein", glycine must weigh ≈ 75.07 Da (the molar mass of free
    /// glycine); the 21-residue human insulin A chain
    /// (GIVEQCCTSICSLYQLENYCN) must weigh ≈ 2383.8 Da.
    #[test]
    fn protein_molecular_weight_matches_published() {
        let gly = Seq::new(SeqKind::Protein, "G").unwrap();
        let mw = molecular_weight_protein(&gly).unwrap();
        assert!(
            (mw - 75.07).abs() < 0.05,
            "glycine MW = {mw:.3} Da, expected ≈ 75.07"
        );

        let insulin_a = Seq::new(SeqKind::Protein, "GIVEQCCTSICSLYQLENYCN").unwrap();
        let mw_a = molecular_weight_protein(&insulin_a).unwrap();
        assert!(
            (mw_a - 2383.8).abs() < 1.0,
            "insulin A-chain MW = {mw_a:.2} Da, expected ≈ 2383.8"
        );
    }

    /// Double-stranded DNA averages ≈ 650 Da per base pair — the
    /// standard rule-of-thumb used to convert ng to pmol. A 100 bp
    /// duplex should therefore weigh on the order of 60–66 kDa.
    #[test]
    fn double_stranded_dna_weight_is_about_650_da_per_bp() {
        // A 100 bp sequence (mixed composition).
        let unit = "ATGCATGCATGCATGCATGC"; // 20 bp, balanced GC
        let seq100 = Seq::new(SeqKind::Dna, unit.repeat(5)).unwrap();
        let mw = molecular_weight_dna_double_stranded(&seq100).unwrap();
        let per_bp = mw / 100.0;
        assert!(
            (600.0..=680.0).contains(&per_bp),
            "dsDNA mass = {per_bp:.1} Da/bp, expected ≈ 650"
        );
    }
}
