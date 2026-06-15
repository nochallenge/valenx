//! Ready-made **illustrative** scoring matrices.
//!
//! These are *teaching* matrices: they encode documented qualitative anchor
//! preferences so the engine is usable out of the box, but they are **not**
//! quantitative experimental matrices (IEDB / SYFPEITHI / NetMHCpan). Treat any
//! result as illustrative, never as a clinical immunogenicity assessment.

use crate::aa::{aa_index, N_AA};
use crate::matrix::Pssm;

fn reward(row: &mut [f64; N_AA], residue: u8, w: f64) {
    if let Some(i) = aa_index(residue) {
        row[i] = w;
    }
}

/// An **illustrative** 9-mer class-I matrix in the style of HLA-A\*02:01.
///
/// HLA-A\*02:01 is the textbook example of anchor-dominated class-I binding:
/// position 2 favours the aliphatic residues Leu/Met and the C-terminal
/// position 9 favours Val/Leu (Falk *et al.*, 1991). This matrix rewards those
/// anchors and leaves the other positions neutral — enough to rank
/// anchor-bearing peptides above non-anchor peptides, and nothing more.
///
/// It is deliberately simple and **not** a quantitative experimental matrix.
pub fn illustrative_hla_a0201() -> Pssm {
    let mut rows: Vec<[f64; N_AA]> = vec![[0.0; N_AA]; 9];
    // P2 anchor (index 1): Leu/Met strongly favoured, Ile weakly tolerated.
    reward(&mut rows[1], b'L', 1.0);
    reward(&mut rows[1], b'M', 1.0);
    reward(&mut rows[1], b'I', 0.5);
    // P-Omega anchor (index 8): Val/Leu favoured.
    reward(&mut rows[8], b'V', 1.0);
    reward(&mut rows[8], b'L', 1.0);
    Pssm::new("HLA-A*02:01 (illustrative anchor motif)", rows)
        .expect("hand-built weights are all finite")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scan::{scan, top_n};

    #[test]
    fn matrix_is_nine_long() {
        assert_eq!(illustrative_hla_a0201().length(), 9);
    }

    #[test]
    fn anchor_bearing_peptide_outscores_non_anchor() {
        let m = illustrative_hla_a0201();
        // P2 = L, P9 = V  -> both anchors satisfied.
        let good = m.score("ALAAAAAAV").unwrap();
        // P2 = P, P9 = P  -> neither anchor satisfied.
        let bad = m.score("APAAAAAAP").unwrap();
        assert!(good > bad);
        assert!((good - 2.0).abs() < 1e-12);
        assert!(bad.abs() < 1e-12);
    }

    #[test]
    fn scan_finds_the_anchor_window_as_top_hit() {
        let m = illustrative_hla_a0201();
        // Embed one strong-anchor 9-mer ("ALAAAAAAV") inside a longer protein.
        let protein = "GGGALAAAAAAVGGG";
        let hits = scan(&m, protein).unwrap();
        let best = top_n(hits, 1);
        assert_eq!(best[0].peptide, "ALAAAAAAV");
        assert!((best[0].score - 2.0).abs() < 1e-12);
    }
}
