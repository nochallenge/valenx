//! Simple De Bruijn graph assembler (v1).
//!
//! De Bruijn graph assembly is the paradigm behind every short-read
//! assembler — Velvet, SPAdes, ABySS, MEGAHIT. The recipe:
//!
//! 1. break every read into k-mers;
//! 2. build a graph whose nodes are the `(k−1)`-mers and whose edges
//!    are the k-mers (each k-mer connects its prefix to its suffix);
//! 3. simplify the graph — collapse unambiguous (non-branching) paths
//!    into unitigs, clip short dead-end tips, pop simple bubbles;
//! 4. emit the unitig sequences as contigs.
//!
//! This module implements all four steps as a real working v1.
//!
//! ## v1 scope
//!
//! A correct graph-algorithm assembler, not SPAdes at scale: it builds
//! the graph in memory (fine for a small genome / a single locus), the
//! bubble-popping handles **simple** two-path bubbles only, and it does
//! not do paired-end scaffolding, error-corrected multi-k assembly, or
//! repeat resolution beyond unitig collapse. It is exact and
//! deterministic.

use crate::error::{GenomicsError, Result};
use std::collections::{BTreeMap, HashMap, HashSet};

/// A De Bruijn graph built from a set of reads.
///
/// Nodes are `(k−1)`-mers; an edge `u → v` exists when some k-mer has
/// `u` as its prefix and `v` as its suffix. Edge multiplicity (the
/// k-mer count) drives tip and bubble decisions.
#[derive(Clone, Debug)]
pub struct DeBruijnGraph {
    /// The k-mer length the graph was built with.
    pub k: usize,
    /// Adjacency: node → list of `(successor, multiplicity)`.
    edges: BTreeMap<Vec<u8>, Vec<(Vec<u8>, u32)>>,
    /// In-degree of every node (counting edge multiplicity as 1 per
    /// distinct edge).
    indegree: HashMap<Vec<u8>, usize>,
}

impl DeBruijnGraph {
    /// Builds the graph from a set of reads with k-mer length `k`.
    ///
    /// k-mers containing an ambiguous base are skipped. Returns
    /// [`GenomicsError::Invalid`] when `k < 2` (a `(k−1)`-mer node
    /// needs `k ≥ 2`).
    pub fn build(reads: &[&[u8]], k: usize) -> Result<Self> {
        if k < 2 {
            return Err(GenomicsError::invalid("k", "k must be >= 2"));
        }
        // Count k-mers.
        let mut kmer_counts: HashMap<Vec<u8>, u32> = HashMap::new();
        for read in reads {
            let upper: Vec<u8> = read.iter().map(|b| b.to_ascii_uppercase()).collect();
            if upper.len() < k {
                continue;
            }
            for w in upper.windows(k) {
                if w.iter().any(|&b| !matches!(b, b'A' | b'C' | b'G' | b'T')) {
                    continue;
                }
                *kmer_counts.entry(w.to_vec()).or_insert(0) += 1;
            }
        }
        // Build the edge map: each k-mer prefix -> suffix.
        let mut edges: BTreeMap<Vec<u8>, Vec<(Vec<u8>, u32)>> = BTreeMap::new();
        let mut indegree: HashMap<Vec<u8>, usize> = HashMap::new();
        for (kmer, count) in &kmer_counts {
            let prefix = kmer[..k - 1].to_vec();
            let suffix = kmer[1..].to_vec();
            edges
                .entry(prefix.clone())
                .or_default()
                .push((suffix.clone(), *count));
            *indegree.entry(suffix).or_insert(0) += 1;
            indegree.entry(prefix).or_insert(0);
        }
        Ok(DeBruijnGraph { k, edges, indegree })
    }

    /// Number of distinct nodes (`(k−1)`-mers).
    pub fn node_count(&self) -> usize {
        let mut nodes: HashSet<&Vec<u8>> = self.edges.keys().collect();
        for succs in self.edges.values() {
            for (s, _) in succs {
                nodes.insert(s);
            }
        }
        nodes.len()
    }

    /// Number of distinct directed edges.
    pub fn edge_count(&self) -> usize {
        self.edges.values().map(|v| v.len()).sum()
    }

    /// Adjacency snapshot — for every node, the list of successor
    /// `(k−1)`-mers. Cheap clone; intended for downstream graph
    /// algorithms that need the raw structure (e.g. the local
    /// haplotype reassembler in [`crate::variant::haplotype`]).
    pub fn adjacency(&self) -> HashMap<Vec<u8>, Vec<Vec<u8>>> {
        let mut out: HashMap<Vec<u8>, Vec<Vec<u8>>> = HashMap::with_capacity(self.edges.len());
        for (k, vs) in &self.edges {
            out.insert(k.clone(), vs.iter().map(|(s, _)| s.clone()).collect());
        }
        out
    }

    /// Out-degree of a node.
    fn outdeg(&self, node: &[u8]) -> usize {
        self.edges.get(node).map(|v| v.len()).unwrap_or(0)
    }

    /// In-degree of a node.
    fn indeg(&self, node: &[u8]) -> usize {
        self.indegree.get(node).copied().unwrap_or(0)
    }

    /// Collapses the graph into **unitigs** — maximal non-branching
    /// paths — and returns their sequences.
    ///
    /// A unitig starts at any node that is *not* the interior of a
    /// simple path (in-degree ≠ 1 or out-degree ≠ 1) and extends while
    /// the next node has in-degree 1 and the current node out-degree 1.
    /// Isolated simple cycles are emitted once each.
    pub fn unitigs(&self) -> Vec<Vec<u8>> {
        let mut contigs = Vec::new();
        let mut used_edges: HashSet<(Vec<u8>, Vec<u8>)> = HashSet::new();

        // All nodes.
        let mut nodes: Vec<Vec<u8>> = self.all_nodes();
        nodes.sort();

        // Start unitigs at non-simple-interior nodes.
        for start in &nodes {
            let is_branch_start = self.indeg(start) != 1 || self.outdeg(start) != 1;
            if !is_branch_start {
                continue;
            }
            if let Some(succs) = self.edges.get(start) {
                for (succ, _) in succs {
                    if used_edges.contains(&(start.clone(), succ.clone())) {
                        continue;
                    }
                    let seq = self.walk_unitig(start, succ, &mut used_edges);
                    if !seq.is_empty() {
                        contigs.push(seq);
                    }
                }
            }
        }

        // Any edges still unused belong to pure cycles — emit them.
        for start in &nodes {
            if let Some(succs) = self.edges.get(start) {
                for (succ, _) in succs {
                    if !used_edges.contains(&(start.clone(), succ.clone())) {
                        let seq = self.walk_unitig(start, succ, &mut used_edges);
                        if !seq.is_empty() {
                            contigs.push(seq);
                        }
                    }
                }
            }
        }

        // A graph with nodes but no edges (e.g. one read shorter than
        // any extension): emit each isolated node.
        if contigs.is_empty() {
            for n in &nodes {
                if self.outdeg(n) == 0 && self.indeg(n) == 0 {
                    contigs.push(n.clone());
                }
            }
        }
        contigs
    }

    /// Walks a single unitig forward from `start` along edge
    /// `start → first`, consuming edges into `used_edges`.
    fn walk_unitig(
        &self,
        start: &[u8],
        first: &[u8],
        used_edges: &mut HashSet<(Vec<u8>, Vec<u8>)>,
    ) -> Vec<u8> {
        let mut seq: Vec<u8> = start.to_vec();
        let mut cur = start.to_vec();
        let mut next = first.to_vec();
        loop {
            if used_edges.contains(&(cur.clone(), next.clone())) {
                break;
            }
            used_edges.insert((cur.clone(), next.clone()));
            // Append the last base of `next` — the new base this edge
            // contributes.
            seq.push(*next.last().unwrap());
            // Continue only along a simple path.
            if self.indeg(&next) == 1 && self.outdeg(&next) == 1 {
                let succ = &self.edges[&next][0].0;
                cur = next.clone();
                next = succ.clone();
                // Cycle guard: a self-returning simple cycle.
                if used_edges.contains(&(cur.clone(), next.clone())) {
                    break;
                }
            } else {
                break;
            }
        }
        seq
    }

    /// Returns every distinct node, sorted.
    fn all_nodes(&self) -> Vec<Vec<u8>> {
        let mut nodes: HashSet<Vec<u8>> = self.edges.keys().cloned().collect();
        for succs in self.edges.values() {
            for (s, _) in succs {
                nodes.insert(s.clone());
            }
        }
        nodes.into_iter().collect()
    }

    /// Removes short dead-end **tips** — a tip is a low-coverage branch
    /// that leaves a node and dead-ends within `max_tip_len` edges.
    /// Returns the number of tip edges removed. A tip is usually an
    /// uncorrected sequencing error.
    pub fn clip_tips(&mut self, max_tip_len: usize) -> usize {
        let mut removed = 0usize;
        let nodes = self.all_nodes();
        let mut to_remove: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();

        for node in &nodes {
            // A node that branches (out-degree > 1) may emit tips.
            if self.outdeg(node) <= 1 {
                continue;
            }
            let succs = self.edges[node].clone();
            // Find the strongest edge to keep.
            let max_mult = succs.iter().map(|(_, m)| *m).max().unwrap_or(0);
            for (succ, mult) in &succs {
                if *mult >= max_mult {
                    continue; // keep the strongest
                }
                // Does this branch dead-end quickly?
                if self.is_short_tip(succ, max_tip_len) {
                    to_remove.push((node.clone(), succ.clone()));
                }
            }
        }
        for (u, v) in to_remove {
            if let Some(succs) = self.edges.get_mut(&u) {
                let before = succs.len();
                succs.retain(|(s, _)| s != &v);
                if succs.len() < before {
                    removed += 1;
                    if let Some(d) = self.indegree.get_mut(&v) {
                        *d = d.saturating_sub(1);
                    }
                }
            }
        }
        removed
    }

    /// `true` when a forward walk from `node` dead-ends within
    /// `max_len` edges (following only out-degree-1 steps).
    fn is_short_tip(&self, node: &[u8], max_len: usize) -> bool {
        let mut cur = node.to_vec();
        for _ in 0..max_len {
            match self.outdeg(&cur) {
                0 => return true, // dead end
                1 => {
                    cur = self.edges[&cur][0].0.clone();
                }
                _ => return false, // re-branches — not a tip
            }
        }
        // Still going after max_len — not a short tip.
        self.outdeg(&cur) == 0
    }

    /// Pops **simple bubbles** — a pair of parallel paths of equal
    /// length between the same two nodes, keeping the higher-coverage
    /// path. Returns the number of bubbles popped.
    ///
    /// Only the simplest bubble topology is handled: a source node with
    /// exactly two out-edges that re-converge at one sink after one
    /// intermediate node each.
    pub fn pop_bubbles(&mut self) -> usize {
        let mut popped = 0usize;
        let nodes = self.all_nodes();
        let mut edges_to_drop: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();

        for src in &nodes {
            if self.outdeg(src) != 2 {
                continue;
            }
            let succs = self.edges[src].clone();
            let (mid_a, mult_a) = (&succs[0].0, succs[0].1);
            let (mid_b, mult_b) = (&succs[1].0, succs[1].1);
            // Each intermediate node must have in-degree 1, out-degree 1.
            if self.indeg(mid_a) != 1
                || self.outdeg(mid_a) != 1
                || self.indeg(mid_b) != 1
                || self.outdeg(mid_b) != 1
            {
                continue;
            }
            let sink_a = &self.edges[mid_a][0].0;
            let sink_b = &self.edges[mid_b][0].0;
            // The two paths must re-converge.
            if sink_a != sink_b {
                continue;
            }
            // Drop the lower-coverage side.
            let (loser_src_edge, loser_mid) = if mult_a < mult_b {
                ((src.clone(), mid_a.clone()), mid_a.clone())
            } else {
                ((src.clone(), mid_b.clone()), mid_b.clone())
            };
            edges_to_drop.push(loser_src_edge);
            edges_to_drop.push((loser_mid.clone(), sink_a.clone()));
            popped += 1;
        }
        for (u, v) in edges_to_drop {
            if let Some(succs) = self.edges.get_mut(&u) {
                let before = succs.len();
                succs.retain(|(s, _)| s != &v);
                if succs.len() < before {
                    if let Some(d) = self.indegree.get_mut(&v) {
                        *d = d.saturating_sub(1);
                    }
                }
            }
        }
        popped
    }
}

/// Parameters for the one-call [`assemble`] entry point.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct AssemblyParams {
    /// k-mer length.
    pub k: usize,
    /// Maximum length (in edges) of a dead-end tip to clip.
    pub max_tip_len: usize,
    /// Whether to pop simple bubbles.
    pub pop_bubbles: bool,
    /// Drop contigs shorter than this from the final output.
    pub min_contig_len: usize,
}

impl Default for AssemblyParams {
    /// Reasonable defaults for a small assembly.
    fn default() -> Self {
        AssemblyParams {
            k: 21,
            max_tip_len: 5,
            pop_bubbles: true,
            min_contig_len: 0,
        }
    }
}

/// One-call De Bruijn assembly: build the graph, simplify it, emit
/// contigs.
pub fn assemble(reads: &[&[u8]], params: &AssemblyParams) -> Result<Vec<Vec<u8>>> {
    let mut graph = DeBruijnGraph::build(reads, params.k)?;
    graph.clip_tips(params.max_tip_len);
    if params.pop_bubbles {
        graph.pop_bubbles();
    }
    let mut contigs: Vec<Vec<u8>> = graph
        .unitigs()
        .into_iter()
        .filter(|c| c.len() >= params.min_contig_len.max(1))
        .collect();
    contigs.sort_by(|a, b| b.len().cmp(&a.len()).then(a.cmp(b)));
    Ok(contigs)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn assembles_a_single_sequence_from_overlapping_reads() {
        // A 40 bp "genome", tiled by overlapping 20-mer reads. The
        // genome is non-repetitive at k = 11 (every 11-mer and every
        // 10-mer node is distinct), so a De Bruijn assembler resolves
        // it into a single unitig — a repeat longer than k would
        // necessarily fragment the assembly.
        let genome = b"GCCGAATAGGGATATAGGCAACGACATGTGCGGCGACCCT";
        let reads: Vec<&[u8]> = vec![
            &genome[0..20],
            &genome[5..25],
            &genome[10..30],
            &genome[15..35],
            &genome[20..40],
        ];
        let contigs = assemble(
            &reads,
            &AssemblyParams {
                k: 11,
                max_tip_len: 3,
                pop_bubbles: true,
                min_contig_len: 0,
            },
        )
        .unwrap();
        assert!(!contigs.is_empty());
        // The longest contig should reconstruct the whole genome.
        let longest = contigs.iter().max_by_key(|c| c.len()).unwrap();
        assert_eq!(longest.as_slice(), genome.as_slice());
    }

    #[test]
    fn reconstructs_original_string_from_its_kmers() {
        // GROUND TRUTH (Eulerian-path reconstruction): a string whose
        // De Bruijn graph is a single non-branching path must reassemble
        // EXACTLY from the (unordered) set of its own k-mers — this is the
        // defining correctness property of De Bruijn assembly.
        //
        // "GATTACAGGCTA" at k=4 has all 9 four-mers distinct and all
        // 3-mer nodes of in/out-degree ≤ 1 (verified: zero branching
        // nodes, one source, one sink), so the unique Eulerian path is
        // the original string and the assembler must emit it verbatim.
        //
        // Unlike the existing read-tiling tests, this drives assembly from
        // the explicit k-mer multiset (each k-mer fed as its own read) and
        // shuffles their order, proving order-independence of the result.
        let original: &[u8] = b"GATTACAGGCTA";
        let k = 4usize;
        // Enumerate the k-mers, then present them in a deliberately
        // non-sequential order to show reconstruction does not depend on
        // input ordering.
        let mut kmers: Vec<&[u8]> = original.windows(k).collect();
        let n = kmers.len();
        assert_eq!(n, original.len() - k + 1, "expected {n} k-mers");
        kmers.reverse();
        if n >= 3 {
            kmers.swap(0, n / 2); // perturb the order further
        }

        let contigs = assemble(
            &kmers,
            &AssemblyParams {
                k,
                max_tip_len: 0, // nothing to clip — the path is clean
                pop_bubbles: false,
                min_contig_len: 0,
            },
        )
        .unwrap();

        // Exactly one unitig, equal to the original string (not merely a
        // substring or rotation — the path is linear and unique).
        assert_eq!(
            contigs.len(),
            1,
            "single clean path must yield one contig, got {contigs:?}"
        );
        assert_eq!(
            contigs[0].as_slice(),
            original,
            "reconstructed {:?} != original {:?}",
            String::from_utf8_lossy(&contigs[0]),
            String::from_utf8_lossy(original),
        );
    }

    #[test]
    fn rejects_small_k() {
        let reads: Vec<&[u8]> = vec![b"ACGT"];
        assert!(DeBruijnGraph::build(&reads, 1).is_err());
    }

    #[test]
    fn graph_has_nodes_and_edges() {
        let reads: Vec<&[u8]> = vec![b"ACGTACGT"];
        let g = DeBruijnGraph::build(&reads, 4).unwrap();
        assert!(g.node_count() > 0);
        assert!(g.edge_count() > 0);
    }

    #[test]
    fn linear_sequence_one_contig() {
        // A genuinely non-repetitive sequence (every 7-mer and 6-mer
        // distinct) with one tiling read — the previous sequence
        // repeated "ACGTTGCA", which a De Bruijn graph cannot collapse
        // into a single contig.
        let seq = b"AAGCCCAATAAACCACTCTGACTG";
        let reads: Vec<&[u8]> = vec![seq];
        let contigs = assemble(
            &reads,
            &AssemblyParams {
                k: 7,
                max_tip_len: 3,
                pop_bubbles: false,
                min_contig_len: 0,
            },
        )
        .unwrap();
        // The whole read reassembles into one contig.
        assert!(contigs.iter().any(|c| c.as_slice() == seq.as_slice()));
    }

    #[test]
    fn tip_clipping_removes_error_branch() {
        // A clean path covered many times, plus a single erroneous read
        // creating a low-coverage tip.
        let clean: &[u8] = b"AAAACCCCGGGGTTTTACGT";
        let mut reads: Vec<&[u8]> = vec![clean; 10];
        // An erroneous read that branches off then dead-ends.
        let err: &[u8] = b"AAAACCCCGGGGTTTTAAAA";
        reads.push(err);
        let mut g = DeBruijnGraph::build(&reads, 7).unwrap();
        let removed = g.clip_tips(5);
        // The error branch should have been clipped.
        assert!(removed > 0, "expected a tip to be clipped");
    }

    #[test]
    fn pop_bubbles_runs() {
        // Two near-identical reads forming a SNP bubble.
        let allele_a: &[u8] = b"ACGTACGTACGTACAT";
        let allele_b: &[u8] = b"ACGTACGTACGTACGT";
        let mut reads: Vec<&[u8]> = vec![allele_a; 5];
        reads.extend(std::iter::repeat_n(allele_b, 3));
        let mut g = DeBruijnGraph::build(&reads, 5).unwrap();
        // Bubble popping must not panic and returns a count.
        let _ = g.pop_bubbles();
        let contigs = g.unitigs();
        assert!(!contigs.is_empty());
    }

    #[test]
    fn deterministic_output() {
        let genome = b"ACGTTGCAACGATGCAGGCCAATTACGTACGT";
        let reads: Vec<&[u8]> = vec![&genome[0..16], &genome[8..24], &genome[16..32]];
        let a = assemble(
            &reads,
            &AssemblyParams {
                k: 9,
                ..AssemblyParams::default()
            },
        )
        .unwrap();
        let b = assemble(
            &reads,
            &AssemblyParams {
                k: 9,
                ..AssemblyParams::default()
            },
        )
        .unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn min_contig_length_filters() {
        let genome = b"ACGTTGCAACGATGCAGGCCAATTACGTACGT";
        let reads: Vec<&[u8]> = vec![genome];
        let big = assemble(
            &reads,
            &AssemblyParams {
                k: 9,
                max_tip_len: 3,
                pop_bubbles: false,
                min_contig_len: 1000,
            },
        )
        .unwrap();
        // Nothing reaches 1000 bp.
        assert!(big.is_empty());
    }
}
