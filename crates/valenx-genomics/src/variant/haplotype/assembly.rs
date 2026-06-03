//! Local haplotype reassembly inside an
//! [`ActiveRegion`](super::active::ActiveRegion).
//!
//! Within each active region, GATK HaplotypeCaller rebuilds candidate
//! haplotypes from the supporting reads using a small **De Bruijn
//! graph** that is seeded with the reference sequence. Path enumeration
//! from the leftmost `(k−1)`-mer of the reference to the rightmost
//! yields candidate haplotypes — variants are simply the bases on a
//! path that diverge from the reference path.
//!
//! This module implements that step on top of the crate's general
//! [`crate::assembly::debruijn::DeBruijnGraph`] by
//! building a fresh small graph per region (so the global assembler
//! used for whole-genome assembly is reused, not duplicated). Then it
//! enumerates source→sink paths up to a bounded number — this is the
//! discipline that lets the caller stay tractable on repetitive
//! regions.
//!
//! ## v1 scope
//!
//! - The graph is built from the *reads* covering the region; the
//!   reference is added as one extra "read" so the reference path is
//!   always present (the same trick GATK uses).
//! - A cycle-bounded BFS enumerates source-to-sink paths; the total
//!   number of haplotypes returned is capped by
//!   [`LocalAssemblyParams::max_haplotypes`].
//! - If the graph contains no source / sink (e.g. the reference's
//!   `(k−1)`-mer flanks do not occur in the read pool), the local
//!   assembler falls back to returning just the reference haplotype —
//!   this is the documented "noop" fall-through that keeps the caller
//!   safe on degenerate inputs.

use crate::assembly::debruijn::DeBruijnGraph;
use std::collections::{HashMap, HashSet, VecDeque};

/// One bounded-BFS frame: `(current node, reconstructed sequence,
/// edges already used on this path)`. Aliased so the `VecDeque` we
/// push into stays readable.
type EnumFrame = (Vec<u8>, Vec<u8>, HashSet<(Vec<u8>, Vec<u8>)>);

/// One candidate haplotype reassembled from the reads in a region.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Haplotype {
    /// The full haplotype sequence (upper-case ASCII bases).
    pub bases: Vec<u8>,
    /// `true` when this haplotype is the unmodified reference path.
    pub is_reference: bool,
}

impl Haplotype {
    /// `true` when the haplotype sequence equals the reference.
    pub fn matches_reference(&self, reference: &[u8]) -> bool {
        self.bases.as_slice() == reference
    }
}

/// Tunables for [`assemble_local_haplotypes`].
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct LocalAssemblyParams {
    /// k-mer length for the per-region De Bruijn graph. Must be `>= 2`.
    pub k: usize,
    /// Maximum number of haplotypes returned per region. Hard cap on
    /// path enumeration to keep the enumeration bounded.
    pub max_haplotypes: usize,
    /// Maximum *expansion* in node-visit count beyond the reference
    /// path length. Keeps cycles bounded.
    pub max_path_expansion: usize,
}

impl Default for LocalAssemblyParams {
    fn default() -> Self {
        LocalAssemblyParams {
            k: 10,
            max_haplotypes: 64,
            max_path_expansion: 128,
        }
    }
}

/// Reassembles candidate haplotypes for a region.
///
/// `reference` is the reference sub-sequence covering the region.
/// `reads` is the set of supporting reads (already trimmed to the
/// region). The reference is always returned as the first haplotype;
/// any genuine reassembled paths follow.
pub fn assemble_local_haplotypes(
    reference: &[u8],
    reads: &[&[u8]],
    params: &LocalAssemblyParams,
) -> Vec<Haplotype> {
    if reference.is_empty() {
        return Vec::new();
    }
    let k = params.k.max(2);
    // The reference must be long enough to seed a (k-1)-mer.
    if reference.len() < k {
        return vec![reference_only(reference)];
    }

    // Build the graph from reads ∪ reference (the reference path is
    // always present in the graph, so the reference haplotype always
    // exists in the path enumeration).
    let upper_ref: Vec<u8> = reference
        .iter()
        .map(|b| b.to_ascii_uppercase())
        .collect();
    let mut pool: Vec<Vec<u8>> = reads
        .iter()
        .map(|r| r.iter().map(|b| b.to_ascii_uppercase()).collect::<Vec<u8>>())
        .filter(|r| r.len() >= k)
        .collect();
    pool.push(upper_ref.clone());
    let pool_refs: Vec<&[u8]> = pool.iter().map(|v| v.as_slice()).collect();

    let graph = match DeBruijnGraph::build(&pool_refs, k) {
        Ok(g) => g,
        Err(_) => return vec![reference_only(reference)],
    };

    // Source = first (k-1)-mer of the reference; sink = last (k-1)-mer.
    let source: Vec<u8> = upper_ref[..k - 1].to_vec();
    let sink: Vec<u8> = upper_ref[upper_ref.len() - (k - 1)..].to_vec();

    let mut bases_of: Vec<Vec<u8>> =
        enumerate_paths(&graph, &source, &sink, params);
    if bases_of.is_empty() {
        return vec![reference_only(reference)];
    }
    // Dedup and tag the reference haplotype.
    bases_of.sort();
    bases_of.dedup();

    let mut out: Vec<Haplotype> = Vec::new();
    // Reference first, always present.
    out.push(reference_only(reference));
    for b in bases_of {
        if b == upper_ref {
            continue; // already in the output
        }
        out.push(Haplotype {
            bases: b,
            is_reference: false,
        });
        if out.len() >= params.max_haplotypes {
            break;
        }
    }
    out
}

fn reference_only(reference: &[u8]) -> Haplotype {
    Haplotype {
        bases: reference
            .iter()
            .map(|b| b.to_ascii_uppercase())
            .collect::<Vec<u8>>(),
        is_reference: true,
    }
}

/// Bounded BFS enumeration of paths from `source` to `sink`.
///
/// The path returns the *reconstructed sequence* (the `(k−1)`-mer plus
/// one base per edge). Cycles are bounded by an overall expansion cap:
/// a path is dropped when its node count exceeds `reference_len +
/// max_path_expansion` (so a self-loop or a long bubble cannot blow up
/// the enumeration).
fn enumerate_paths(
    graph: &DeBruijnGraph,
    source: &[u8],
    sink: &[u8],
    params: &LocalAssemblyParams,
) -> Vec<Vec<u8>> {
    let cap = params.max_haplotypes.saturating_add(1).max(2);
    let visit_cap = params.max_path_expansion.saturating_add(source.len()).max(64);

    let edges = graph_edges(graph);

    let mut out: Vec<Vec<u8>> = Vec::new();
    // BFS queue of EnumFrame: (current node, reconstructed sequence so
    // far, edges visited on this path).
    let mut q: VecDeque<EnumFrame> = VecDeque::new();
    q.push_back((source.to_vec(), source.to_vec(), HashSet::new()));

    while let Some((node, seq, used)) = q.pop_front() {
        if seq.len() > visit_cap {
            continue;
        }
        if node == sink {
            out.push(seq.clone());
            if out.len() >= cap {
                break;
            }
            // Keep going — there may be more paths.
        }
        if let Some(succs) = edges.get(&node) {
            for succ in succs {
                let edge_key = (node.clone(), succ.clone());
                if used.contains(&edge_key) {
                    // Already used this edge on the current path — a
                    // cycle. The expansion cap will eventually stop us
                    // but skipping is cheaper.
                    continue;
                }
                let mut new_used = used.clone();
                new_used.insert(edge_key);
                let mut new_seq = seq.clone();
                new_seq.push(*succ.last().unwrap());
                q.push_back((succ.clone(), new_seq, new_used));
                if q.len() > 4 * cap.max(64) {
                    // Prevent runaway queue growth on dense graphs.
                    break;
                }
            }
        }
    }
    out
}

/// Materialises the graph's adjacency for path enumeration. We rebuild
/// a `HashMap<Vec<u8>, Vec<Vec<u8>>>` because [`DeBruijnGraph`] keeps
/// its edges private — but the simplified-graph public surface
/// ([`DeBruijnGraph::unitigs`]) loses the per-node adjacency we need.
/// Cheap on small per-region graphs.
fn graph_edges(graph: &DeBruijnGraph) -> HashMap<Vec<u8>, Vec<Vec<u8>>> {
    // We reconstruct adjacency by walking the unitigs and the contained
    // k-mers — but that is fragile. Instead, we depend on a tiny public
    // helper added to the assembly module: see `DeBruijnGraph::adjacency`.
    graph.adjacency()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn upper(s: &[u8]) -> Vec<u8> {
        s.iter().map(|b| b.to_ascii_uppercase()).collect()
    }

    #[test]
    fn reference_haplotype_always_emitted() {
        let reference = b"ACGTACGTACGTACGTACGT";
        let reads: Vec<&[u8]> = vec![reference];
        let haps = assemble_local_haplotypes(
            reference,
            &reads,
            &LocalAssemblyParams::default(),
        );
        assert!(!haps.is_empty());
        assert!(haps[0].is_reference);
        assert_eq!(haps[0].bases, upper(reference));
    }

    #[test]
    fn snv_haplotype_recovered() {
        // A 30 bp non-repetitive reference; reads carry one SNV in
        // the middle.
        let reference = b"GCATCGATCGATGCATCGATCGATCGATGC";
        let mut alt: Vec<u8> = reference.to_vec();
        alt[15] = b'T'; // was C
        // Many reads carry the alt; a few carry the reference.
        let mut reads: Vec<&[u8]> = Vec::new();
        let ref_slice: &[u8] = reference;
        reads.extend(std::iter::repeat_n(alt.as_slice(), 12));
        reads.extend(std::iter::repeat_n(ref_slice, 3));
        let haps = assemble_local_haplotypes(
            reference,
            &reads,
            &LocalAssemblyParams {
                k: 8,
                ..LocalAssemblyParams::default()
            },
        );
        assert!(
            haps.iter().any(|h| h.bases == upper(reference)),
            "missing reference haplotype: {haps:?}"
        );
        assert!(
            haps.iter().any(|h| h.bases == upper(&alt)),
            "missing alt haplotype: {haps:?}"
        );
    }

    #[test]
    fn insertion_haplotype_recovered() {
        // Reference and a 2 bp insertion. Use a non-repetitive context
        // so the (k-1)-mer flanks unambiguously seed the path
        // enumeration.
        let reference = b"AAGTCGATGCCGTAACGGTCAGGTAACCGATTCGAATGCC";
        let mut alt: Vec<u8> = reference[..20].to_vec();
        alt.extend_from_slice(b"TT");
        alt.extend_from_slice(&reference[20..]);
        let mut reads: Vec<&[u8]> = Vec::new();
        for _ in 0..10 {
            reads.push(&alt);
        }
        let haps = assemble_local_haplotypes(
            reference,
            &reads,
            &LocalAssemblyParams {
                k: 8,
                ..LocalAssemblyParams::default()
            },
        );
        assert!(
            haps.iter().any(|h| h.bases == upper(&alt)),
            "insertion haplotype missing in {:?}",
            haps.iter()
                .map(|h| String::from_utf8_lossy(&h.bases).into_owned())
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn deletion_haplotype_recovered() {
        let reference = b"AAGTCGATGCCGTAACGGTCAGGTAACCGATTCGAATGCC";
        // Delete one base at position 20.
        let mut alt: Vec<u8> = reference[..20].to_vec();
        alt.extend_from_slice(&reference[21..]);
        let mut reads: Vec<&[u8]> = Vec::new();
        for _ in 0..10 {
            reads.push(&alt);
        }
        let haps = assemble_local_haplotypes(
            reference,
            &reads,
            &LocalAssemblyParams {
                k: 8,
                ..LocalAssemblyParams::default()
            },
        );
        assert!(
            haps.iter().any(|h| h.bases == upper(&alt)),
            "deletion haplotype missing in {:?}",
            haps.iter().map(|h| String::from_utf8_lossy(&h.bases).into_owned()).collect::<Vec<_>>()
        );
    }

    #[test]
    fn short_reference_falls_back_to_reference_only() {
        let reference = b"AC"; // shorter than k
        let reads: Vec<&[u8]> = vec![];
        let haps = assemble_local_haplotypes(
            reference,
            &reads,
            &LocalAssemblyParams::default(),
        );
        assert_eq!(haps.len(), 1);
        assert!(haps[0].is_reference);
    }

    #[test]
    fn empty_reference_yields_nothing() {
        let reference: &[u8] = &[];
        let reads: Vec<&[u8]> = vec![];
        let haps = assemble_local_haplotypes(
            reference,
            &reads,
            &LocalAssemblyParams::default(),
        );
        assert!(haps.is_empty());
    }

    #[test]
    fn haplotype_count_bounded() {
        let reference = b"GCATCGATCGATGCATCGATCGATCGATGCATCGATCGAT";
        let reads: Vec<&[u8]> = vec![reference];
        let params = LocalAssemblyParams {
            k: 7,
            max_haplotypes: 3,
            ..LocalAssemblyParams::default()
        };
        let haps = assemble_local_haplotypes(reference, &reads, &params);
        assert!(haps.len() <= 3);
    }

    #[test]
    fn duplicate_haplotypes_removed() {
        let reference = b"GCATCGATCGATGCATCGATCGAT";
        // Duplicate reads should not multiply the same haplotype.
        let reads: Vec<&[u8]> = vec![reference, reference, reference, reference];
        let haps = assemble_local_haplotypes(
            reference,
            &reads,
            &LocalAssemblyParams {
                k: 7,
                ..LocalAssemblyParams::default()
            },
        );
        let n = haps
            .iter()
            .filter(|h| h.bases == upper(reference))
            .count();
        assert_eq!(n, 1, "reference must appear exactly once");
    }
}
