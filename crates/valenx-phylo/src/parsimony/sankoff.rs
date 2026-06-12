//! Sankoff weighted parsimony (Sankoff 1975).
//!
//! Sankoff's algorithm generalises [Fitch](super::fitch): instead of
//! every change costing 1, a [`CostMatrix`] gives an arbitrary cost for
//! each ordered state-to-state transition. The minimum total cost is
//! found by dynamic programming over the tree.
//!
//! For one character, a single bottom-up (postorder) pass fills, at
//! every node, a vector `g[s]` = the minimum subtree cost given that
//! the node is in state `s`:
//!
//! - **Leaf:** `g[observed] = 0`, every other `g[s] = +∞` (a gap makes
//!   all `g[s] = 0`).
//! - **Internal:** `g[s] = Σ_children min_t ( cost(s→t) + child.g[t] )`.
//!
//! The character's parsimony cost is `min_s root.g[s]`. A top-down pass
//! then back-traces one minimum-cost ancestral state per node.
//!
//! With a unit cost matrix (0 on the diagonal, 1 elsewhere) Sankoff
//! reproduces the Fitch score exactly — a property the tests check.

use crate::error::{PhyloError, Result};
use crate::tree::{NodeId, Tree};

/// A square matrix of state-transition costs.
///
/// `cost(from, to)` must be `0` when `from == to` and non-negative
/// everywhere; it need not be symmetric.
#[derive(Debug, Clone, PartialEq)]
pub struct CostMatrix {
    n: usize,
    data: Vec<f64>,
}

impl CostMatrix {
    /// Builds a cost matrix from a flat row-major `n × n` buffer.
    ///
    /// # Errors
    /// [`PhyloError::Dimension`] if `data.len() != n²`;
    /// [`PhyloError::Invalid`] if a diagonal entry is non-zero or any
    /// entry is negative.
    pub fn new(n: usize, data: Vec<f64>) -> Result<Self> {
        if data.len() != n * n {
            return Err(PhyloError::dimension(n * n, data.len(), "cost matrix"));
        }
        for i in 0..n {
            if data[i * n + i] != 0.0 {
                return Err(PhyloError::invalid(
                    "cost_matrix",
                    "diagonal entries must be zero",
                ));
            }
        }
        if data.iter().any(|&c| c < 0.0) {
            return Err(PhyloError::invalid(
                "cost_matrix",
                "costs must be non-negative",
            ));
        }
        Ok(CostMatrix { n, data })
    }

    /// A unit cost matrix: `0` on the diagonal, `1` off it. With this
    /// matrix, Sankoff parsimony equals Fitch parsimony.
    pub fn unit(n: usize) -> Self {
        let mut data = vec![1.0; n * n];
        for i in 0..n {
            data[i * n + i] = 0.0;
        }
        CostMatrix { n, data }
    }

    /// A transition / transversion matrix for 4 nucleotide states
    /// (A,C,G,T order): transitions A↔G and C↔T cost `ti`, all
    /// transversions cost `tv`.
    pub fn transition_transversion(ti: f64, tv: f64) -> Self {
        // A=0, C=1, G=2, T=3. Purine = {0,2}, pyrimidine = {1,3}.
        let mut data = vec![tv; 16];
        for i in 0..4 {
            data[i * 4 + i] = 0.0;
        }
        // Index (row*4 + col): A=0, C=1, G=2, T=3.
        data[2] = ti; // A->G
        data[8] = ti; // G->A
        data[7] = ti; // C->T
        data[13] = ti; // T->C
        CostMatrix { n: 4, data }
    }

    /// Number of states.
    pub fn n_states(&self) -> usize {
        self.n
    }

    /// Cost of changing `from` → `to`.
    ///
    /// # Panics
    /// If either index is out of range.
    pub fn cost(&self, from: usize, to: usize) -> f64 {
        self.data[from * self.n + to]
    }
}

/// Result of a Sankoff weighted-parsimony analysis.
#[derive(Debug, Clone, PartialEq)]
pub struct SankoffResult {
    /// Total weighted parsimony cost over every column.
    pub cost: f64,
    /// Per-column cost (length = alignment width).
    pub site_costs: Vec<f64>,
    /// Back-traced ancestral states, indexed `[node][column]`.
    pub ancestral: Vec<Vec<u8>>,
}

/// Runs Sankoff weighted parsimony on a tree, an alignment and a cost
/// matrix.
///
/// `alignment` maps a leaf label to its row of `u8` state indices; a
/// state byte of `u8::MAX` is a gap / missing wildcard.
///
/// # Errors
/// - [`PhyloError::Invalid`] on an empty alignment or a leaf with no
///   row.
/// - [`PhyloError::Dimension`] if the rows differ in width.
pub fn sankoff_parsimony(
    tree: &Tree,
    alignment: &[(String, Vec<u8>)],
    costs: &CostMatrix,
) -> Result<SankoffResult> {
    if alignment.is_empty() {
        return Err(PhyloError::invalid("alignment", "no sequences supplied"));
    }
    let width = alignment[0].1.len();
    if width == 0 {
        return Err(PhyloError::invalid("alignment", "zero-width alignment"));
    }
    for (_, row) in alignment {
        if row.len() != width {
            return Err(PhyloError::dimension(width, row.len(), "alignment rows"));
        }
    }
    let k = costs.n_states();
    let n = tree.node_count();

    let row_for = |id: NodeId| -> Result<&Vec<u8>> {
        let label = tree
            .node(id)
            .label
            .as_deref()
            .ok_or_else(|| PhyloError::invalid("tree", "leaf without a label"))?;
        alignment
            .iter()
            .find(|(name, _)| name == label)
            .map(|(_, row)| row)
            .ok_or_else(|| PhyloError::invalid("alignment", format!("no row for leaf `{label}`")))
    };

    let post = tree.postorder();
    let pre = tree.preorder();
    let mut site_costs = vec![0.0; width];
    let mut ancestral = vec![vec![u8::MAX; width]; n];

    // `col` indexes several per-column arrays (site_costs, ancestral
    // rows, the alignment rows) — a range loop is the clearest form.
    #[allow(clippy::needless_range_loop)]
    for col in 0..width {
        // g[node][state] — min subtree cost given the node's state.
        let mut g = vec![vec![f64::INFINITY; k]; n];
        for &id in &post {
            let node = tree.node(id);
            if node.is_leaf() {
                let s = row_for(id)?[col];
                if s == u8::MAX || s as usize >= k {
                    // Gap: any state is free.
                    g[id].iter_mut().for_each(|v| *v = 0.0);
                } else {
                    g[id][s as usize] = 0.0;
                }
            } else {
                for s in 0..k {
                    let mut total = 0.0;
                    for &c in &node.children {
                        let mut best = f64::INFINITY;
                        for t in 0..k {
                            let v = costs.cost(s, t) + g[c][t];
                            if v < best {
                                best = v;
                            }
                        }
                        total += best;
                    }
                    g[id][s] = total;
                }
            }
        }
        let root = tree.root();
        let (root_state, root_cost) =
            g[root]
                .iter()
                .enumerate()
                .fold((0usize, f64::INFINITY), |(bs, bc), (s, &c)| {
                    if c < bc {
                        (s, c)
                    } else {
                        (bs, bc)
                    }
                });
        site_costs[col] = root_cost;
        ancestral[root][col] = root_state as u8;

        // Top-down back-trace: each node picks the state minimising
        // cost(parent_state -> state) + g[node][state].
        for &id in &pre {
            if id == root {
                continue;
            }
            let parent = tree.node(id).parent.expect("non-root has a parent");
            let ps = ancestral[parent][col] as usize;
            let (best_state, _) = (0..k).fold((0usize, f64::INFINITY), |(bs, bc), t| {
                let v = costs.cost(ps, t) + g[id][t];
                if v < bc {
                    (t, v)
                } else {
                    (bs, bc)
                }
            });
            ancestral[id][col] = best_state as u8;
        }
    }

    let cost = site_costs.iter().sum();
    Ok(SankoffResult {
        cost,
        site_costs,
        ancestral,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::newick::read_newick;
    use crate::parsimony::fitch::fitch_parsimony;

    fn row(label: &str, states: &[u8]) -> (String, Vec<u8>) {
        (label.to_string(), states.to_vec())
    }

    #[test]
    fn cost_matrix_validation() {
        assert!(CostMatrix::new(2, vec![0.0, 1.0, 1.0, 0.0]).is_ok());
        // Non-zero diagonal.
        assert!(CostMatrix::new(2, vec![1.0, 1.0, 1.0, 0.0]).is_err());
        // Negative cost.
        assert!(CostMatrix::new(2, vec![0.0, -1.0, 1.0, 0.0]).is_err());
        // Wrong shape.
        assert!(CostMatrix::new(2, vec![0.0, 1.0]).is_err());
    }

    #[test]
    fn unit_cost_matches_fitch() {
        let tree = read_newick("((A,B),(C,D));").unwrap();
        let aln = vec![
            row("A", &[0, 0]),
            row("B", &[0, 1]),
            row("C", &[1, 0]),
            row("D", &[1, 1]),
        ];
        let fitch = fitch_parsimony(&tree, &aln, 4).unwrap();
        let sankoff = sankoff_parsimony(&tree, &aln, &CostMatrix::unit(4)).unwrap();
        assert!((sankoff.cost - fitch.score as f64).abs() < 1e-9);
        for (s, f) in sankoff.site_costs.iter().zip(&fitch.site_scores) {
            assert!((s - *f as f64).abs() < 1e-9);
        }
    }

    #[test]
    fn weighted_costs_change_the_score() {
        let tree = read_newick("((A,B),(C,D));").unwrap();
        // A,B in state 0 (A); C,D in state 2 (G) — one A<->G change.
        let aln = vec![
            row("A", &[0]),
            row("B", &[0]),
            row("C", &[2]),
            row("D", &[2]),
        ];
        // Transitions cheap (1), transversions expensive (5).
        let ts = CostMatrix::transition_transversion(1.0, 5.0);
        let r = sankoff_parsimony(&tree, &aln, &ts).unwrap();
        // One transition => cost 1.
        assert!((r.cost - 1.0).abs() < 1e-9);

        // Now A,B = A(0), C,D = C(1) — a transversion.
        let aln2 = vec![
            row("A", &[0]),
            row("B", &[0]),
            row("C", &[1]),
            row("D", &[1]),
        ];
        let r2 = sankoff_parsimony(&tree, &aln2, &ts).unwrap();
        assert!((r2.cost - 5.0).abs() < 1e-9);
    }

    #[test]
    fn ancestral_back_trace_is_complete() {
        let tree = read_newick("((A,B),(C,D));").unwrap();
        let aln = vec![
            row("A", &[0]),
            row("B", &[0]),
            row("C", &[1]),
            row("D", &[1]),
        ];
        let r = sankoff_parsimony(&tree, &aln, &CostMatrix::unit(4)).unwrap();
        for id in 0..tree.node_count() {
            assert_ne!(r.ancestral[id][0], u8::MAX);
        }
    }

    #[test]
    fn gaps_are_free() {
        let tree = read_newick("((A,B),(C,D));").unwrap();
        let aln = vec![
            row("A", &[0]),
            row("B", &[0]),
            row("C", &[0]),
            row("D", &[u8::MAX]),
        ];
        let r = sankoff_parsimony(&tree, &aln, &CostMatrix::unit(4)).unwrap();
        assert!(r.cost.abs() < 1e-9);
    }
}
