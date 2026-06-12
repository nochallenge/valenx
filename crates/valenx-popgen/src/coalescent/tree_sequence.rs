//! A succinct tree sequence — the `tskit`-class data structure.
//!
//! A *tree sequence* records the genealogy of a sample along a
//! recombining chromosome as four tables:
//!
//! - **Node table** — every genome ever ancestral to the sample, with
//!   a `time` (generations before present) and a `is_sample` flag.
//! - **Edge table** — a parent/child relationship valid over a genomic
//!   interval `[left, right)`. Because nearby positions usually share
//!   the same genealogy, one edge covers a whole non-recombined
//!   stretch; the chromosome's local tree only changes at edge
//!   boundaries. This is the succinctness `tskit` is famous for.
//! - **Site table** — the segregating sites (positions + ancestral
//!   state).
//! - **Mutation table** — which node a mutation arose on, at which
//!   site.
//!
//! [`TreeSequence`] stores all four and can extract the *local tree* at
//! any genomic position as a [`valenx_phylo::Tree`], so the whole of
//! `valenx-phylo`'s rendering and comparison machinery applies to a
//! simulated genealogy for free.

use crate::error::{PopgenError, Result};
use serde::{Deserialize, Serialize};
use valenx_phylo::tree::{Node as PhyloNode, Tree};

/// A node of a tree sequence — one ancestral genome.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TsNode {
    /// Time of the node, in generations before the present.
    pub time: f64,
    /// `true` if the node is part of the sampled set.
    pub is_sample: bool,
}

/// An edge: `child` inherits from `parent` over `[left, right)`.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Edge {
    /// Parent node id.
    pub parent: usize,
    /// Child node id.
    pub child: usize,
    /// Inclusive left genomic coordinate.
    pub left: f64,
    /// Exclusive right genomic coordinate.
    pub right: f64,
}

/// A segregating site of a tree sequence.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TsSite {
    /// Genomic position.
    pub position: f64,
    /// Ancestral allelic state (`0` by convention).
    pub ancestral_state: u8,
}

/// A mutation: a state change on a `node` at a `site`.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TsMutation {
    /// Index into the site table.
    pub site: usize,
    /// Node on which the mutation arose; all descendants over the
    /// site's covering edge interval inherit the derived state.
    pub node: usize,
    /// Derived allelic state introduced by the mutation.
    pub derived_state: u8,
}

/// A succinct tree sequence: node + edge + site + mutation tables.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TreeSequence {
    nodes: Vec<TsNode>,
    edges: Vec<Edge>,
    sites: Vec<TsSite>,
    mutations: Vec<TsMutation>,
    sequence_length: f64,
    finalized: bool,
}

impl TreeSequence {
    /// Creates an empty tree sequence over a chromosome of the given
    /// length.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] if `sequence_length <= 0`.
    pub fn new(sequence_length: f64) -> Result<Self> {
        if sequence_length <= 0.0 {
            return Err(PopgenError::invalid("sequence_length", "must be positive"));
        }
        Ok(TreeSequence {
            nodes: Vec::new(),
            edges: Vec::new(),
            sites: Vec::new(),
            mutations: Vec::new(),
            sequence_length,
            finalized: false,
        })
    }

    /// Appends a node and returns its id.
    pub fn add_node(&mut self, time: f64, is_sample: bool) -> usize {
        self.nodes.push(TsNode { time, is_sample });
        self.nodes.len() - 1
    }

    /// Appends an edge.
    ///
    /// # Errors
    /// [`PopgenError::Model`] if either endpoint is out of range or the
    /// interval is empty / outside the chromosome.
    pub fn add_edge(&mut self, edge: Edge) -> Result<()> {
        if edge.parent >= self.nodes.len() || edge.child >= self.nodes.len() {
            return Err(PopgenError::model("edge references a missing node"));
        }
        if edge.right <= edge.left {
            return Err(PopgenError::model("edge interval is empty"));
        }
        if edge.left < 0.0 || edge.right > self.sequence_length + 1e-9 {
            return Err(PopgenError::model("edge interval outside the chromosome"));
        }
        self.edges.push(edge);
        Ok(())
    }

    /// Appends a site and returns its id.
    pub fn add_site(&mut self, position: f64, ancestral_state: u8) -> usize {
        self.sites.push(TsSite {
            position,
            ancestral_state,
        });
        self.sites.len() - 1
    }

    /// Appends a mutation.
    ///
    /// # Errors
    /// [`PopgenError::Model`] if the node or site is out of range.
    pub fn add_mutation(&mut self, mutation: TsMutation) -> Result<()> {
        if mutation.node >= self.nodes.len() {
            return Err(PopgenError::model("mutation references a missing node"));
        }
        if mutation.site >= self.sites.len() {
            return Err(PopgenError::model("mutation references a missing site"));
        }
        self.mutations.push(mutation);
        Ok(())
    }

    /// Finalises the tree sequence: validates the parent/child time
    /// ordering (a parent must be strictly older than its child) and
    /// sorts the edge table by left coordinate.
    ///
    /// # Errors
    /// [`PopgenError::Model`] if any edge has a parent no older than
    /// its child.
    pub fn finalize(&mut self) -> Result<()> {
        for e in &self.edges {
            if self.nodes[e.parent].time <= self.nodes[e.child].time {
                return Err(PopgenError::model(
                    "edge parent is not older than its child",
                ));
            }
        }
        self.edges.sort_by(|a, b| {
            a.left
                .partial_cmp(&b.left)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        self.finalized = true;
        Ok(())
    }

    /// Number of nodes.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Number of edges.
    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }

    /// Number of sites.
    pub fn site_count(&self) -> usize {
        self.sites.len()
    }

    /// Number of mutations.
    pub fn mutation_count(&self) -> usize {
        self.mutations.len()
    }

    /// Number of sample nodes.
    pub fn sample_count(&self) -> usize {
        self.nodes.iter().filter(|n| n.is_sample).count()
    }

    /// Sample node ids.
    pub fn samples(&self) -> Vec<usize> {
        (0..self.nodes.len())
            .filter(|&i| self.nodes[i].is_sample)
            .collect()
    }

    /// The node table.
    pub fn nodes(&self) -> &[TsNode] {
        &self.nodes
    }

    /// The edge table.
    pub fn edges(&self) -> &[Edge] {
        &self.edges
    }

    /// The site table.
    pub fn sites(&self) -> &[TsSite] {
        &self.sites
    }

    /// The mutation table.
    pub fn mutations(&self) -> &[TsMutation] {
        &self.mutations
    }

    /// Chromosome length.
    pub fn sequence_length(&self) -> f64 {
        self.sequence_length
    }

    /// The number of *distinct local trees*: the chromosome is
    /// partitioned by the set of edge endpoints into intervals, each
    /// with one marginal tree.
    pub fn tree_count(&self) -> usize {
        let mut breaks: Vec<f64> = vec![0.0, self.sequence_length];
        for e in &self.edges {
            breaks.push(e.left);
            breaks.push(e.right);
        }
        breaks.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        breaks.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
        breaks.len().saturating_sub(1)
    }

    /// Extracts the marginal genealogy covering genomic `position` as a
    /// [`valenx_phylo::Tree`].
    ///
    /// Only edges whose interval contains `position` participate. The
    /// resulting parent map is assembled into a phylo arena; node
    /// labels are `n<id>` and branch lengths are time differences.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] if `position` is outside the
    /// chromosome; [`PopgenError::Model`] if the local genealogy is a
    /// forest with more than one root (an incompletely coalesced
    /// sample).
    pub fn local_tree(&self, position: f64) -> Result<Tree> {
        if position < 0.0 || position >= self.sequence_length {
            return Err(PopgenError::invalid("position", "outside the chromosome"));
        }
        // Parent of each ts-node under the local tree.
        let mut parent_of: Vec<Option<usize>> = vec![None; self.nodes.len()];
        for e in &self.edges {
            if position >= e.left && position < e.right {
                parent_of[e.child] = Some(e.parent);
            }
        }
        // Which ts-nodes participate: samples and their ancestors.
        let mut used = vec![false; self.nodes.len()];
        for s in self.samples() {
            let mut cur = Some(s);
            while let Some(c) = cur {
                if used[c] {
                    break;
                }
                used[c] = true;
                cur = parent_of[c];
            }
        }
        let participating: Vec<usize> = (0..self.nodes.len()).filter(|&i| used[i]).collect();
        if participating.is_empty() {
            return Err(PopgenError::model("local tree has no nodes"));
        }
        // Map ts-node -> phylo arena index.
        let mut arena_of = vec![usize::MAX; self.nodes.len()];
        for (idx, &ts) in participating.iter().enumerate() {
            arena_of[ts] = idx;
        }
        // Build phylo nodes.
        let mut phylo_nodes: Vec<PhyloNode> = participating
            .iter()
            .map(|&ts| PhyloNode {
                label: Some(format!("n{ts}")),
                branch_length: None,
                parent: None,
                children: Vec::new(),
            })
            .collect();
        let mut roots = Vec::new();
        for &ts in &participating {
            let me = arena_of[ts];
            match parent_of[ts] {
                Some(p) if used[p] => {
                    let pa = arena_of[p];
                    phylo_nodes[me].parent = Some(pa);
                    phylo_nodes[me].branch_length =
                        Some((self.nodes[p].time - self.nodes[ts].time).abs());
                    phylo_nodes[pa].children.push(me);
                }
                _ => roots.push(me),
            }
        }
        if roots.len() != 1 {
            return Err(PopgenError::model(format!(
                "local tree is a forest with {} roots (sample not fully coalesced)",
                roots.len()
            )));
        }
        Tree::new(phylo_nodes, roots[0], true).map_err(|e| PopgenError::model(e.to_string()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Builds a tiny fixed tree sequence: 3 samples, one MRCA, no
    /// recombination (a single tree).
    fn tiny() -> TreeSequence {
        let mut ts = TreeSequence::new(100.0).unwrap();
        let a = ts.add_node(0.0, true); // 0
        let b = ts.add_node(0.0, true); // 1
        let c = ts.add_node(0.0, true); // 2
        let ab = ts.add_node(1.0, false); // 3
        let root = ts.add_node(2.0, false); // 4
        ts.add_edge(Edge {
            parent: ab,
            child: a,
            left: 0.0,
            right: 100.0,
        })
        .unwrap();
        ts.add_edge(Edge {
            parent: ab,
            child: b,
            left: 0.0,
            right: 100.0,
        })
        .unwrap();
        ts.add_edge(Edge {
            parent: root,
            child: ab,
            left: 0.0,
            right: 100.0,
        })
        .unwrap();
        ts.add_edge(Edge {
            parent: root,
            child: c,
            left: 0.0,
            right: 100.0,
        })
        .unwrap();
        ts.finalize().unwrap();
        ts
    }

    #[test]
    fn tables_have_expected_sizes() {
        let ts = tiny();
        assert_eq!(ts.node_count(), 5);
        assert_eq!(ts.edge_count(), 4);
        assert_eq!(ts.sample_count(), 3);
        assert_eq!(ts.samples(), vec![0, 1, 2]);
        assert_eq!(ts.tree_count(), 1);
    }

    #[test]
    fn edge_validation_catches_bad_endpoints() {
        let mut ts = TreeSequence::new(50.0).unwrap();
        ts.add_node(0.0, true);
        // Edge to a missing node 9.
        assert!(ts
            .add_edge(Edge {
                parent: 9,
                child: 0,
                left: 0.0,
                right: 10.0,
            })
            .is_err());
    }

    #[test]
    fn finalize_rejects_inverted_times() {
        let mut ts = TreeSequence::new(50.0).unwrap();
        let young = ts.add_node(0.0, true);
        let old = ts.add_node(5.0, false);
        // Parent younger than child.
        ts.add_edge(Edge {
            parent: young,
            child: old,
            left: 0.0,
            right: 50.0,
        })
        .unwrap();
        assert!(ts.finalize().is_err());
    }

    #[test]
    fn local_tree_reconstructs_the_genealogy() {
        let ts = tiny();
        let tree = ts.local_tree(50.0).unwrap();
        assert_eq!(tree.leaf_count(), 3);
        assert!(tree.validate().is_ok());
        // Total length: a->ab (1) + b->ab (1) + ab->root (1) +
        // c->root (2) = 5.
        assert!((tree.total_length() - 5.0).abs() < 1e-9);
    }

    #[test]
    fn local_tree_rejects_out_of_range_position() {
        let ts = tiny();
        assert!(ts.local_tree(-1.0).is_err());
        assert!(ts.local_tree(100.0).is_err());
    }

    #[test]
    fn recombination_splits_into_two_trees() {
        // Two samples, two MRCAs over disjoint halves of the genome.
        let mut ts = TreeSequence::new(100.0).unwrap();
        let a = ts.add_node(0.0, true);
        let b = ts.add_node(0.0, true);
        let r1 = ts.add_node(1.0, false);
        let r2 = ts.add_node(3.0, false);
        // Left half coalesces at r1.
        ts.add_edge(Edge {
            parent: r1,
            child: a,
            left: 0.0,
            right: 50.0,
        })
        .unwrap();
        ts.add_edge(Edge {
            parent: r1,
            child: b,
            left: 0.0,
            right: 50.0,
        })
        .unwrap();
        // Right half coalesces at r2.
        ts.add_edge(Edge {
            parent: r2,
            child: a,
            left: 50.0,
            right: 100.0,
        })
        .unwrap();
        ts.add_edge(Edge {
            parent: r2,
            child: b,
            left: 50.0,
            right: 100.0,
        })
        .unwrap();
        ts.finalize().unwrap();
        assert_eq!(ts.tree_count(), 2);
        // The two local trees differ in height.
        let left = ts.local_tree(25.0).unwrap();
        let right = ts.local_tree(75.0).unwrap();
        assert!(right.total_length() > left.total_length());
    }
}
