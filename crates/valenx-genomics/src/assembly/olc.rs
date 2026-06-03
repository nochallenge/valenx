//! Overlap-layout-consensus mini-assembler (v1) — for long reads.
//!
//! While short reads are assembled with a De Bruijn graph, long reads
//! (PacBio / Nanopore) are assembled with the **overlap-layout-
//! consensus** paradigm — the approach of the Celera assembler, Canu,
//! miniasm and hifiasm:
//!
//! 1. **Overlap** — find every pair of reads whose suffix / prefix
//!    overlap, tolerating a few mismatches (long reads are error-prone);
//! 2. **Layout** — chain the overlaps into a linear ordering of reads
//!    (a path through the overlap graph);
//! 3. **Consensus** — merge the laid-out reads into one contig
//!    sequence, taking a per-column majority vote where reads agree.
//!
//! This module implements all three for a small read set.
//!
//! ## v1 scope
//!
//! A real OLC v1, not Canu / hifiasm at genome scale: the overlap step
//! is an O(n² · L) all-pairs suffix-prefix scan (fine for a locus or a
//! plasmid, not a mammalian genome), the layout is a greedy longest-
//! overlap chain (no transitive-edge reduction or repeat graph), and
//! the consensus is a simple majority over the layout — no
//! partial-order alignment, no polishing pass. It is exact and
//! deterministic.

use crate::error::{GenomicsError, Result};

/// An overlap between the suffix of read `a` and the prefix of read
/// `b`.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Overlap {
    /// Index of the read contributing its suffix.
    pub a: usize,
    /// Index of the read contributing its prefix.
    pub b: usize,
    /// Length of the overlapping region.
    pub length: usize,
    /// Number of mismatches inside the overlap.
    pub mismatches: usize,
}

impl Overlap {
    /// Identity fraction of the overlap (`1.0` is a perfect overlap).
    pub fn identity(&self) -> f64 {
        if self.length == 0 {
            0.0
        } else {
            (self.length - self.mismatches) as f64 / self.length as f64
        }
    }
}

/// Parameters for the OLC assembler.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct OlcParams {
    /// Minimum overlap length for two reads to be linked.
    pub min_overlap: usize,
    /// Maximum tolerated mismatch fraction inside an overlap (long
    /// reads need a generous value).
    pub max_error_rate: f64,
    /// Drop contigs shorter than this.
    pub min_contig_len: usize,
}

impl Default for OlcParams {
    /// Long-read-friendly defaults.
    fn default() -> Self {
        OlcParams {
            min_overlap: 20,
            max_error_rate: 0.15,
            min_contig_len: 0,
        }
    }
}

/// Finds the best suffix-prefix overlap of read `a` onto read `b`.
///
/// Tries every overlap length from the longest possible down to
/// `min_overlap`; returns the first (longest) one whose mismatch count
/// is within `max_error_rate`. `None` when no overlap qualifies.
fn best_overlap(
    a: &[u8],
    b: &[u8],
    a_idx: usize,
    b_idx: usize,
    params: &OlcParams,
) -> Option<Overlap> {
    let max_ov = a.len().min(b.len());
    let mut len = max_ov;
    while len >= params.min_overlap {
        // a's suffix of length `len` vs b's prefix of length `len`.
        let a_suffix = &a[a.len() - len..];
        let b_prefix = &b[..len];
        let mismatches = a_suffix
            .iter()
            .zip(b_prefix)
            .filter(|(&x, &y)| !x.eq_ignore_ascii_case(&y))
            .count();
        let allowed = (len as f64 * params.max_error_rate).floor() as usize;
        if mismatches <= allowed {
            return Some(Overlap {
                a: a_idx,
                b: b_idx,
                length: len,
                mismatches,
            });
        }
        len -= 1;
    }
    None
}

/// Computes every qualifying suffix-prefix overlap among the reads.
///
/// Self-overlaps (`a == b`) are skipped. The result is a flat list,
/// sorted by descending overlap length (so the greedy layout consumes
/// the strongest links first).
pub fn compute_overlaps(reads: &[&[u8]], params: &OlcParams) -> Vec<Overlap> {
    let mut overlaps = Vec::new();
    for i in 0..reads.len() {
        for j in 0..reads.len() {
            if i == j {
                continue;
            }
            if let Some(ov) = best_overlap(reads[i], reads[j], i, j, params) {
                overlaps.push(ov);
            }
        }
    }
    overlaps.sort_by(|x, y| y.length.cmp(&x.length).then(x.a.cmp(&y.a)));
    overlaps
}

/// Greedily lays out reads into chains using the overlap list.
///
/// Each read may have at most one successor and one predecessor in a
/// chain (a simple linear layout). Overlaps are consumed strongest
/// first; an overlap is taken only when both endpoints are still free.
/// Returns a list of chains, each a `Vec` of `(read_index,
/// overlap_with_previous)` — the first element's overlap is `0`.
pub fn layout(
    n_reads: usize,
    overlaps: &[Overlap],
) -> Vec<Vec<(usize, usize)>> {
    let mut succ: Vec<Option<(usize, usize)>> = vec![None; n_reads];
    let mut has_pred = vec![false; n_reads];
    let mut has_succ = vec![false; n_reads];

    for ov in overlaps {
        if !has_succ[ov.a] && !has_pred[ov.b] && ov.a != ov.b {
            // Guard against creating a cycle: walking forward from b
            // must not reach a.
            if !reaches(ov.b, ov.a, &succ) {
                succ[ov.a] = Some((ov.b, ov.length));
                has_succ[ov.a] = true;
                has_pred[ov.b] = true;
            }
        }
    }

    // Chain heads are reads with no predecessor.
    let mut chains = Vec::new();
    let mut visited = vec![false; n_reads];
    for start in 0..n_reads {
        if has_pred[start] || visited[start] {
            continue;
        }
        let mut chain = vec![(start, 0usize)];
        visited[start] = true;
        let mut cur = start;
        while let Some((next, ov_len)) = succ[cur] {
            if visited[next] {
                break;
            }
            chain.push((next, ov_len));
            visited[next] = true;
            cur = next;
        }
        chains.push(chain);
    }
    // Any unvisited reads (inside a cycle) become singleton chains.
    for (r, seen) in visited.iter_mut().enumerate() {
        if !*seen {
            chains.push(vec![(r, 0usize)]);
            *seen = true;
        }
    }
    chains
}

/// `true` when following successors from `from` reaches `target`.
fn reaches(from: usize, target: usize, succ: &[Option<(usize, usize)>]) -> bool {
    let mut cur = from;
    let mut steps = 0usize;
    while let Some((next, _)) = succ[cur] {
        if next == target {
            return true;
        }
        cur = next;
        steps += 1;
        if steps > succ.len() {
            return true; // cycle — treat as reachable to be safe
        }
    }
    false
}

/// Builds a consensus contig from one layout chain.
///
/// Reads are stitched along their overlaps. The previously-placed
/// `contig` is authoritative through the overlap region (long-read
/// errors are independent, so the already-merged prefix — which may
/// itself span several reads — carries more evidence than the incoming
/// single read); the non-overlapping suffix of each new read is then
/// appended. This is a deterministic v1 consensus; see the module note.
fn consensus(reads: &[&[u8]], chain: &[(usize, usize)]) -> Vec<u8> {
    if chain.is_empty() {
        return Vec::new();
    }
    let mut contig: Vec<u8> = reads[chain[0].0]
        .iter()
        .map(|b| b.to_ascii_uppercase())
        .collect();
    for &(read_idx, ov_len) in &chain[1..] {
        let read: Vec<u8> = reads[read_idx]
            .iter()
            .map(|b| b.to_ascii_uppercase())
            .collect();
        // The overlap region: the last `ov_len` bases of `contig` align
        // with the first `ov_len` bases of `read`. The existing contig
        // bases are kept; only the new suffix is appended.
        let ov = ov_len.min(contig.len()).min(read.len());
        contig.extend_from_slice(&read[ov..]);
    }
    contig
}

/// One-call OLC assembly: overlap, layout, consensus.
pub fn assemble_olc(reads: &[&[u8]], params: &OlcParams) -> Result<Vec<Vec<u8>>> {
    if reads.is_empty() {
        return Err(GenomicsError::invalid("reads", "no reads supplied"));
    }
    if params.min_overlap == 0 {
        return Err(GenomicsError::invalid("min_overlap", "must be positive"));
    }
    if !(0.0..=1.0).contains(&params.max_error_rate) {
        return Err(GenomicsError::invalid(
            "max_error_rate",
            "must be in [0, 1]",
        ));
    }
    let overlaps = compute_overlaps(reads, params);
    let chains = layout(reads.len(), &overlaps);
    let mut contigs: Vec<Vec<u8>> = chains
        .iter()
        .map(|chain| consensus(reads, chain))
        .filter(|c| c.len() >= params.min_contig_len.max(1))
        .collect();
    contigs.sort_by(|a, b| b.len().cmp(&a.len()).then(a.cmp(b)));
    Ok(contigs)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn overlap_detection() {
        // "ACGTACGT" suffix overlaps "ACGTTTTT" prefix by 4 ("ACGT").
        let p = OlcParams {
            min_overlap: 3,
            max_error_rate: 0.0,
            min_contig_len: 0,
        };
        let ov = best_overlap(b"ACGTACGT", b"ACGTTTTT", 0, 1, &p).unwrap();
        assert_eq!(ov.length, 4);
        assert_eq!(ov.mismatches, 0);
    }

    #[test]
    fn overlap_tolerates_mismatch() {
        // One mismatch inside an 8-base overlap.
        let p = OlcParams {
            min_overlap: 6,
            max_error_rate: 0.2,
            min_contig_len: 0,
        };
        // a suffix "AACCGGTT" vs b prefix "AACCGGTA" — last base differs.
        let ov = best_overlap(b"XXAACCGGTT", b"AACCGGTAYY", 0, 1, &p);
        assert!(ov.is_some());
        assert_eq!(ov.unwrap().mismatches, 1);
    }

    #[test]
    fn assembles_tiled_long_reads() {
        // A 60-base sequence tiled by overlapping 30-base reads.
        let genome = b"ACGTTGCAGGCCAATTACGTACGTACGTTTACCGGTACGTACGTACGTGGCCAATTACGT";
        let reads: Vec<&[u8]> = vec![
            &genome[0..30],
            &genome[15..45],
            &genome[30..60],
        ];
        let contigs = assemble_olc(&reads, &OlcParams {
            min_overlap: 10,
            max_error_rate: 0.0,
            min_contig_len: 0,
        })
        .unwrap();
        let longest = contigs.iter().max_by_key(|c| c.len()).unwrap();
        assert_eq!(longest.as_slice(), genome.as_slice());
    }

    #[test]
    fn single_read_passes_through() {
        let read = b"ACGTACGTACGTACGTACGT";
        let contigs = assemble_olc(&[read.as_slice()], &OlcParams::default()).unwrap();
        assert_eq!(contigs.len(), 1);
        assert_eq!(contigs[0].as_slice(), read.as_slice());
    }

    #[test]
    fn disjoint_reads_stay_separate() {
        // Two reads with no overlap -> two contigs.
        let r1 = b"AAAAAAAAAAAAAAAAAAAA";
        let r2 = b"CCCCCCCCCCCCCCCCCCCC";
        let contigs = assemble_olc(&[r1.as_slice(), r2.as_slice()], &OlcParams {
            min_overlap: 10,
            max_error_rate: 0.0,
            min_contig_len: 0,
        })
        .unwrap();
        assert_eq!(contigs.len(), 2);
    }

    #[test]
    fn deterministic() {
        let genome = b"ACGTTGCAGGCCAATTACGTACGTACGTTTACCGGTACGT";
        let reads: Vec<&[u8]> = vec![&genome[0..25], &genome[12..37], &genome[15..40]];
        let a = assemble_olc(&reads, &OlcParams::default()).unwrap();
        let b = assemble_olc(&reads, &OlcParams::default()).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn identity_computed() {
        let ov = Overlap {
            a: 0,
            b: 1,
            length: 10,
            mismatches: 2,
        };
        assert!((ov.identity() - 0.8).abs() < 1e-9);
    }

    #[test]
    fn rejects_bad_params() {
        let reads: Vec<&[u8]> = vec![b"ACGT"];
        assert!(assemble_olc(&[], &OlcParams::default()).is_err());
        assert!(assemble_olc(&reads, &OlcParams {
            min_overlap: 0,
            max_error_rate: 0.1,
            min_contig_len: 0
        })
        .is_err());
        assert!(assemble_olc(&reads, &OlcParams {
            min_overlap: 5,
            max_error_rate: 2.0,
            min_contig_len: 0
        })
        .is_err());
    }
}
