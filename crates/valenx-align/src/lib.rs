//! # valenx-align — sequence alignment and homology search
//!
//! Round 6 Block 2 of the Valenx roadmap. A native-Rust replacement
//! for the alignment and search core of BLAST+, BWA, Bowtie2,
//! minimap2, ClustalΩ, MUSCLE, MAFFT, T-Coffee, HMMER, DIAMOND and
//! MMseqs2 — pure dynamic-programming algorithms, no neural-network
//! weights and no external processes.
//!
//! It builds on [`valenx_bioseq`] (Block 6.1): the
//! [`Seq`](valenx_bioseq::Seq) type and its IUPAC alphabets are the
//! natural input for the [`bio`] convenience layer, while the
//! performance-critical DP works directly on `&[u8]`.
//!
//! ## What it does
//!
//! - **Scoring** ([`matrix`]) — the BLOSUM / PAM substitution-matrix
//!   families, `NUC.4.4`, and an affine [`matrix::GapCost`] /
//!   [`matrix::ScoringScheme`].
//! - **Pairwise** ([`pairwise`]) — Needleman-Wunsch, Gotoh affine,
//!   Smith-Waterman, semi-global / overlap, banded and Hirschberg
//!   linear-space alignment, all yielding an
//!   [`pairwise::Alignment`] with CIGAR and statistics.
//! - **Search** ([`search`]) — a k-mer index, BLAST-class
//!   seed-and-extend with Karlin-Altschul E-values, an FM-index / BWT
//!   for exact matching, minimizer sketches and anchor chaining.
//! - **Multiple** ([`msa`]) — distance matrices, UPGMA / NJ guide
//!   trees, progressive alignment, MUSCLE-class iterative refinement
//!   and column-conservation analysis.
//! - **HMM** ([`hmm`]) — a pair HMM (forward + Viterbi), a Plan7-style
//!   profile HMM built from an MSA, and PSSM motif scanning.
//! - **I/O & misc** ([`io`], [`util`]) — Clustal / Stockholm / PHYLIP
//!   / aligned-FASTA / MSF readers and writers, Levenshtein and Myers
//!   edit distance, dot-plots, SAM CIGAR conversion and a
//!   seed-and-extend read mapper.
//!
//! ## Errors
//!
//! Every fallible function returns
//! [`Result<_, AlignError>`](error::AlignError); the error type
//! carries stable [`code`](error::AlignError::code) and
//! [`category`](error::AlignError::category) accessors for telemetry.
//!
//! ## v1 scope
//!
//! This is a real working v1, not production parity with the 30-year
//! reference tools. Each module documents its own simplifications
//! inline; the notable ones are: seed-and-extend is a real
//! heuristic but not NCBI BLAST+'s production indexing; the FM-index
//! stores an uncompressed Occ table and matches exactly (the
//! heuristic search covers inexact); Karlin-Altschul `λ`/`K` are the
//! published presets rather than solved per scoring system; the
//! profile-HMM transition probabilities use fixed Plan7 defaults
//! rather than counted state paths; and the read mapper handles the
//! forward strand single-end only.

#![forbid(unsafe_code)]
// Surface future `&str` byte-offset slicing in clippy review — this
// crate parses untrusted text alignment formats (PHYLIP/Clustal/etc.),
// where non-char-boundary slices panic. WARN (not deny): most existing
// slices are safe ASCII; this only flags NEW ones.
#![warn(clippy::string_slice)]

pub mod bio;
pub mod error;
pub mod hmm;
pub mod io;
pub mod limits;
pub mod matrix;
pub mod msa;
pub mod pairwise;
pub mod search;
pub mod util;

// --- Convenience re-exports of the most-used types --------------------

pub use error::{AlignError, ErrorCategory, Result};
pub use limits::{check_dp_size, MAX_DP_CELLS};
pub use matrix::{GapCost, ScoringScheme, SubstitutionMatrix};
pub use msa::{align as align_msa, Msa};
pub use pairwise::{Alignment, AlignStats, Cigar, CigarOp};
pub use pairwise::global::{gotoh, needleman_wunsch};
pub use pairwise::local::smith_waterman;
pub use search::{KmerIndex, SeedSearch};

#[cfg(test)]
mod tests {
    use super::*;

    /// End-to-end pairwise: score two DNA sequences globally and check
    /// CIGAR + stats wiring. The two input sequences differ in a single
    /// base, so the optimal global alignment under an affine scheme is
    /// ungapped: 11 matches + 1 mismatch.
    #[test]
    fn end_to_end_pairwise() {
        use crate::matrix::{GapCost, SubstitutionMatrix};
        let scheme = ScoringScheme::dna_default(); // NUC.4.4, 10/1 affine
        let a = b"ACGTACGTACGT";
        let b = b"ACGTTCGTACGT";

        // The affine-gap routine respects the 10/1 open/extend cost, so
        // a gap-pair (cost 22+) never beats the single mismatch (-4):
        // the optimal alignment is the ungapped 11M+1X = 11*5 - 4 = 51.
        let go = gotoh(a, b, &scheme).unwrap();
        assert_eq!(go.score, 51, "affine global score of an ungapped 11M1X");
        assert_eq!(go.len(), 12);

        // The plain linear-gap Needleman-Wunsch ignores the `open` term
        // by design; to make it agree it must be given a linear scheme
        // (open = 0) — then both routines reach the same ungapped score.
        let linear =
            ScoringScheme::new(SubstitutionMatrix::nuc44(), GapCost::new(0, 5));
        let nw = needleman_wunsch(a, b, &linear).unwrap();
        let go_lin = gotoh(a, b, &linear).unwrap();
        assert_eq!(nw.score, go_lin.score, "linear NW == Gotoh under a linear scheme");
        assert_eq!(nw.len(), 12);

        // CIGAR + stats are wired through.
        let cigar = go.cigar();
        assert_eq!(cigar.ref_len(), 12);
        let stats = go.stats(&scheme.matrix);
        assert!(stats.identities >= 10);
    }

    /// End-to-end search: index a database, run seed-and-extend, and
    /// confirm the planted homology is found with a sane E-value.
    #[test]
    fn end_to_end_search() {
        use search::{KarlinAltschul, SeedParams};
        let db: Vec<&[u8]> = vec![b"TTTTTGATTACACATAGCATGGGGG"];
        let idx = KmerIndex::build_many(&db, 7).unwrap();
        let scheme = ScoringScheme::new(
            SubstitutionMatrix::dna_simple(2, -3),
            GapCost::new(5, 2),
        );
        let search = SeedSearch::new(&idx, db.clone(), &scheme, SeedParams::default())
            .with_stats(KarlinAltschul::dna_ungapped());
        let hsps = search.search(b"GATTACACATAGCAT");
        assert!(!hsps.is_empty());
        assert!(hsps[0].e_value.unwrap() >= 0.0);
    }

    /// End-to-end MSA: align three related sequences, refine, and
    /// confirm the consensus and conservation analysis run.
    #[test]
    fn end_to_end_msa() {
        use msa::{analysis, refine, RefineParams};
        let scheme = ScoringScheme::new(
            SubstitutionMatrix::dna_simple(2, -1),
            GapCost::new(4, 1),
        );
        let seqs: &[&[u8]] = &[b"ACGTACGTACGT", b"ACGTACGTACGT", b"ACGTCGTACGT"];
        let msa = align_msa(seqs, &scheme).unwrap();
        assert_eq!(msa.depth(), 3);
        let refined = refine::refine(&msa, &scheme, RefineParams::default()).unwrap();
        assert!(refined.sum_of_pairs(&scheme) >= msa.sum_of_pairs(&scheme));
        let cons = analysis::consensus(&refined, analysis::ConsensusOptions::default());
        assert_eq!(cons.len(), refined.width());
    }

    #[test]
    fn re_exports_are_wired() {
        let _ = SubstitutionMatrix::blosum62();
        let _ = GapCost::default();
        let _ = ScoringScheme::blosum62_default();
        let e = AlignError::not_yet("x");
        assert_eq!(e.category_enum(), ErrorCategory::Capability);
        let _: Cigar = Cigar::new();
        let _ = CigarOp::Match;
        let al = needleman_wunsch(b"ACGT", b"ACGT", &ScoringScheme::dna_default()).unwrap();
        let _: AlignStats = al.stats(&SubstitutionMatrix::nuc44());
    }
}
