//! Pairwise structure alignment (v1) and TM-score.
//!
//! The aligner is a **sequence-anchored iterative superposition** in
//! the CE / TM-align family:
//!
//! 1. Seed a residue correspondence by a fast sequence alignment of
//!    the two chains' Cα sequences (a Needleman-Wunsch with an
//!    identity-style matrix).
//! 2. Superpose on the seed's matched Cα pairs ([`kabsch`]).
//! 3. Re-derive the correspondence: keep matched pairs whose
//!    superposed Cα–Cα distance is within a cutoff, and re-superpose.
//! 4. Iterate until the matched set stabilises.
//!
//! This is a real working structural aligner. It is **not** a
//! full CE (no fragment-pair combinatorial search) or a full
//! TM-align (no dynamic-programming rotation refinement); a
//! sequence-dissimilar but structurally identical pair will align
//! less well than the reference tools. The [`tm_score`] function is
//! the exact published TM-score formula.

use crate::error::{BiostructError, Result};
use crate::structure::Chain;
use crate::superpose::{kabsch, RigidTransform, Superposition};
use nalgebra::Point3;

/// A pairwise structure-alignment result.
#[derive(Clone, Debug, PartialEq)]
pub struct StructureAlignment {
    /// Aligned residue-index pairs `(i_in_a, j_in_b)`.
    pub matched: Vec<(usize, usize)>,
    /// The transform mapping chain A's coordinates onto chain B.
    pub transform: RigidTransform,
    /// RMSD over the aligned Cα pairs, ångström.
    pub rmsd: f64,
    /// TM-score of the alignment, normalised by the *shorter* chain.
    pub tm_score: f64,
    /// Number of aligned pairs.
    pub aligned_length: usize,
}

/// Cα coordinates of a chain's amino-acid residues, paired with each
/// residue's index in the chain.
fn ca_trace(chain: &Chain) -> Vec<(usize, Point3<f64>)> {
    chain
        .residues
        .iter()
        .enumerate()
        .filter_map(|(i, r)| {
            if r.is_amino_acid() {
                r.ca().map(|a| (i, a.coord))
            } else {
                None
            }
        })
        .collect()
}

/// Align two protein chains by sequence-anchored iterative
/// superposition.
///
/// `distance_cutoff` is the Cα–Cα distance (ångström) below which a
/// superposed pair is kept in the correspondence; `4.0` is a sensible
/// default in the spirit of CE.
pub fn align_chains(
    chain_a: &Chain,
    chain_b: &Chain,
    distance_cutoff: f64,
) -> Result<StructureAlignment> {
    if distance_cutoff <= 0.0 || distance_cutoff.is_nan() {
        return Err(BiostructError::invalid(
            "distance_cutoff",
            "must be positive",
        ));
    }
    let trace_a = ca_trace(chain_a);
    let trace_b = ca_trace(chain_b);
    if trace_a.len() < 3 || trace_b.len() < 3 {
        return Err(BiostructError::invalid(
            "chain",
            "each chain needs at least 3 Cα atoms to align",
        ));
    }

    // --- 1. seed correspondence by sequence alignment --------------
    let seq_a = chain_a.observed_sequence();
    let seq_b = chain_b.observed_sequence();
    // The observed_sequence covers amino acids AND nucleotides, but
    // the traces only have amino acids; rebuild the amino-acid-only
    // one-letter strings parallel to the traces.
    let aa_a: Vec<char> = trace_a
        .iter()
        .map(|(i, _)| crate::structure::residue_one_letter(&chain_a.residues[*i].name))
        .collect();
    let aa_b: Vec<char> = trace_b
        .iter()
        .map(|(i, _)| crate::structure::residue_one_letter(&chain_b.residues[*i].name))
        .collect();
    let _ = (&seq_a, &seq_b); // observed_sequence kept for callers

    let mut matched_trace = needleman_wunsch(&aa_a, &aa_b);
    if matched_trace.len() < 3 {
        // Sequences too dissimilar to seed — fall back to aligning
        // the first min(len) residues positionally.
        let m = trace_a.len().min(trace_b.len());
        matched_trace = (0..m).map(|k| (k, k)).collect();
    }

    // --- 2-4. iterative superposition ------------------------------
    // Each iteration superposes on the current matched set, then
    // re-derives the matched set from the full sequence alignment
    // keeping only superposed pairs within `distance_cutoff`. The
    // final transform / metrics are recomputed once below from the
    // converged matched set, so the per-iteration transform is
    // purely loop-local.
    let mut prev_count = 0;
    for _iter in 0..20 {
        if matched_trace.len() < 3 {
            break;
        }
        let mob: Vec<Point3<f64>> = matched_trace.iter().map(|(a, _)| trace_a[*a].1).collect();
        let refp: Vec<Point3<f64>> = matched_trace.iter().map(|(_, b)| trace_b[*b].1).collect();
        let sup: Superposition = kabsch(&mob, &refp)?;
        let iter_transform = sup.transform;

        // Re-derive the correspondence from the full sequence
        // alignment, keeping only close pairs.
        let seq_pairs = needleman_wunsch(&aa_a, &aa_b);
        let mut kept = Vec::new();
        for (a, b) in seq_pairs {
            let moved = iter_transform.apply(&trace_a[a].1);
            if (moved - trace_b[b].1).norm() <= distance_cutoff {
                kept.push((a, b));
            }
        }
        if kept.len() < 3 {
            // keep the previous matched set; cannot improve.
            break;
        }
        if kept.len() == prev_count {
            matched_trace = kept;
            break;
        }
        prev_count = kept.len();
        matched_trace = kept;
    }

    // Final metrics.
    let mob: Vec<Point3<f64>> = matched_trace.iter().map(|(a, _)| trace_a[*a].1).collect();
    let refp: Vec<Point3<f64>> = matched_trace.iter().map(|(_, b)| trace_b[*b].1).collect();
    let final_sup = kabsch(&mob, &refp)?;
    let transform = final_sup.transform.clone();

    // Map trace indices back to residue indices for the public result.
    let matched: Vec<(usize, usize)> = matched_trace
        .iter()
        .map(|(a, b)| (trace_a[*a].0, trace_b[*b].0))
        .collect();

    let moved = transform.apply_all(&mob);
    let norm_len = trace_a.len().min(trace_b.len());
    let tm = tm_score(&moved, &refp, norm_len)?;

    Ok(StructureAlignment {
        aligned_length: matched.len(),
        matched,
        transform,
        rmsd: final_sup.rmsd,
        tm_score: tm,
    })
}

/// TM-score of two already-superposed, equal-length Cα point sets.
///
/// `TM = (1/L_norm) · Σᵢ 1 / (1 + (dᵢ/d₀)²)` with the published
/// length-dependent scale `d₀(L) = 1.24·∛(L−15) − 1.8` (floored at
/// `0.5`). `length_norm` is the chain length to normalise by — the
/// TM-align convention is the target / shorter chain.
pub fn tm_score(
    superposed: &[Point3<f64>],
    reference: &[Point3<f64>],
    length_norm: usize,
) -> Result<f64> {
    if superposed.len() != reference.len() {
        return Err(BiostructError::invalid(
            "points",
            "TM-score needs equal-length point sets",
        ));
    }
    if length_norm == 0 {
        return Err(BiostructError::invalid(
            "length_norm",
            "normalisation length must be positive",
        ));
    }
    let d0 = tm_d0(length_norm);
    let mut sum = 0.0;
    for (p, q) in superposed.iter().zip(reference) {
        let d = (p - q).norm();
        sum += 1.0 / (1.0 + (d / d0).powi(2));
    }
    Ok(sum / length_norm as f64)
}

/// The TM-score length-dependent distance scale `d₀(L)`.
pub fn tm_d0(length: usize) -> f64 {
    if length <= 15 {
        return 0.5;
    }
    let d0 = 1.24 * ((length as f64 - 15.0).cbrt()) - 1.8;
    d0.max(0.5)
}

/// A minimal Needleman-Wunsch global alignment over two one-letter
/// sequences, returning the matched (non-gap) index pairs.
///
/// Scoring: `+2` identical residue, `−1` mismatch, `−2` gap. This is
/// only the seeding aligner — a full substitution matrix lives in the
/// `valenx-align` crate.
fn needleman_wunsch(a: &[char], b: &[char]) -> Vec<(usize, usize)> {
    let (n, m) = (a.len(), b.len());
    if n == 0 || m == 0 {
        return Vec::new();
    }
    const MATCH: i32 = 2;
    const MISMATCH: i32 = -1;
    const GAP: i32 = -2;

    let mut score = vec![vec![0i32; m + 1]; n + 1];
    for (i, row) in score.iter_mut().enumerate() {
        row[0] = i as i32 * GAP;
    }
    for j in 0..=m {
        score[0][j] = j as i32 * GAP;
    }
    for i in 1..=n {
        for j in 1..=m {
            let s = if a[i - 1] == b[j - 1] && a[i - 1] != 'X' {
                MATCH
            } else {
                MISMATCH
            };
            let diag = score[i - 1][j - 1] + s;
            let up = score[i - 1][j] + GAP;
            let left = score[i][j - 1] + GAP;
            score[i][j] = diag.max(up).max(left);
        }
    }

    // Traceback.
    let mut pairs = Vec::new();
    let (mut i, mut j) = (n, m);
    while i > 0 && j > 0 {
        let s = if a[i - 1] == b[j - 1] && a[i - 1] != 'X' {
            MATCH
        } else {
            MISMATCH
        };
        if score[i][j] == score[i - 1][j - 1] + s {
            pairs.push((i - 1, j - 1));
            i -= 1;
            j -= 1;
        } else if score[i][j] == score[i - 1][j] + GAP {
            i -= 1;
        } else {
            j -= 1;
        }
    }
    pairs.reverse();
    pairs
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::{Atom, Residue};
    use nalgebra::Matrix3;

    /// Build a protein chain with Cα atoms at the given coordinates.
    fn chain_from_cas(name: &str, cas: &[Point3<f64>], seq: &[char]) -> Chain {
        let mut c = Chain::new(name);
        for (k, ca) in cas.iter().enumerate() {
            let resn = match seq.get(k) {
                Some('A') => "ALA",
                Some('G') => "GLY",
                Some('S') => "SER",
                Some('V') => "VAL",
                Some('L') => "LEU",
                _ => "ALA",
            };
            let mut r = Residue::new(resn, k as i32 + 1);
            r.atoms.push(Atom::new("CA", "C", *ca));
            c.residues.push(r);
        }
        c
    }

    fn helix_cas(n: usize) -> Vec<Point3<f64>> {
        (0..n)
            .map(|i| {
                let t = i as f64 * 100.0_f64.to_radians();
                Point3::new(2.3 * t.cos(), 2.3 * t.sin(), i as f64 * 1.5)
            })
            .collect()
    }

    #[test]
    fn d0_increases_with_length() {
        assert!((tm_d0(10) - 0.5).abs() < 1e-12);
        assert!(tm_d0(100) > tm_d0(50));
        assert!(tm_d0(50) > 0.5);
    }

    #[test]
    fn tm_score_of_identical_is_one() {
        let cas = helix_cas(40);
        let tm = tm_score(&cas, &cas, 40).unwrap();
        assert!((tm - 1.0).abs() < 1e-9, "tm of identical = {tm}");
    }

    #[test]
    fn needleman_wunsch_aligns_identical() {
        let a: Vec<char> = "AGSVL".chars().collect();
        let pairs = needleman_wunsch(&a, &a);
        assert_eq!(pairs.len(), 5);
        assert_eq!(pairs[0], (0, 0));
        assert_eq!(pairs[4], (4, 4));
    }

    #[test]
    fn needleman_wunsch_handles_gap() {
        let a: Vec<char> = "AGSVL".chars().collect();
        let b: Vec<char> = "AGVL".chars().collect(); // S deleted
        let pairs = needleman_wunsch(&a, &b);
        // 4 residues of b should all be matched.
        assert_eq!(pairs.len(), 4);
    }

    #[test]
    fn aligns_a_chain_to_itself() {
        let cas = helix_cas(30);
        let seq: Vec<char> = std::iter::repeat_n('A', 30).collect();
        let chain = chain_from_cas("A", &cas, &seq);
        let aln = align_chains(&chain, &chain, 4.0).unwrap();
        assert_eq!(aln.aligned_length, 30);
        assert!(aln.rmsd < 1e-6, "self-alignment rmsd {}", aln.rmsd);
        assert!(aln.tm_score > 0.99, "self-alignment tm {}", aln.tm_score);
    }

    #[test]
    fn aligns_a_rotated_copy() {
        // Rotate + translate a helix; the aligner must recover a
        // near-zero RMSD and TM ~ 1.
        let cas = helix_cas(35);
        let rot = Matrix3::new(0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0);
        let shift = nalgebra::Vector3::new(10.0, 5.0, -7.0);
        let moved: Vec<_> = cas
            .iter()
            .map(|p| Point3::from(rot * p.coords + shift))
            .collect();
        let seq: Vec<char> = std::iter::repeat_n('A', 35).collect();
        let a = chain_from_cas("A", &cas, &seq);
        let b = chain_from_cas("B", &moved, &seq);
        let aln = align_chains(&a, &b, 4.0).unwrap();
        assert!(aln.rmsd < 1e-4, "rotated-copy rmsd {}", aln.rmsd);
        assert!(aln.tm_score > 0.99, "rotated-copy tm {}", aln.tm_score);
    }

    #[test]
    fn rejects_short_chains() {
        let cas = helix_cas(2);
        let seq: Vec<char> = vec!['A', 'A'];
        let chain = chain_from_cas("A", &cas, &seq);
        assert!(align_chains(&chain, &chain, 4.0).is_err());
    }

    #[test]
    fn rejects_bad_cutoff() {
        let cas = helix_cas(10);
        let seq: Vec<char> = std::iter::repeat_n('A', 10).collect();
        let chain = chain_from_cas("A", &cas, &seq);
        assert!(align_chains(&chain, &chain, -1.0).is_err());
    }
}
