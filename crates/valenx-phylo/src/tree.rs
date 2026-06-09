//! The [`Tree`] data model — the central type of `valenx-phylo`.
//!
//! A phylogenetic tree is an *arena* of [`Node`]s: every node carries
//! its parent index and an ordered list of child indices, so the
//! structure is a forest-free single tree with one designated root.
//! An arena (`Vec<Node>` + `usize` ids) is used rather than `Rc`/`Box`
//! pointers because phylogenetic algorithms (NNI, SPR, pruning) need
//! cheap whole-tree copies and index-keyed scratch arrays.
//!
//! Both **rooted** and **unrooted** trees are represented with the same
//! arena; the [`Tree::rooted`] flag records the caller's intent. An
//! unrooted tree is stored rooted at an arbitrary node (conventionally
//! a trifurcating root), and topology-comparison code that must be
//! root-agnostic works from the bipartition set, not the node parents.
//!
//! Branch lengths live on the *child* node ([`Node::branch_length`]) —
//! the length of the edge connecting the node to its parent. The root
//! has no incoming edge, so its branch length is `None`.

use crate::error::{PhyloError, Result};
use serde::{Deserialize, Serialize};

/// Index of a node within a [`Tree`]'s arena.
pub type NodeId = usize;

/// A single node of a [`Tree`].
///
/// Leaves have an empty [`children`](Node::children) list and (almost
/// always) a non-empty [`label`](Node::label). Internal nodes may or
/// may not be labelled; an internal label is conventionally a clade
/// name or a support value written by [`crate::io::newick`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Node {
    /// Node label — a taxon name for leaves, an optional clade name or
    /// support annotation for internal nodes.
    pub label: Option<String>,
    /// Length of the edge from this node up to its parent, if known.
    /// `None` on the root and on trees with no branch lengths.
    pub branch_length: Option<f64>,
    /// Parent node id, or `None` for the root.
    pub parent: Option<NodeId>,
    /// Child node ids, in left-to-right order.
    pub children: Vec<NodeId>,
}

impl Node {
    /// `true` if this node has no children (a tip / terminal taxon).
    pub fn is_leaf(&self) -> bool {
        self.children.is_empty()
    }

    /// `true` if this node has children (an internal / Hennigian node).
    pub fn is_internal(&self) -> bool {
        !self.children.is_empty()
    }
}

/// A phylogenetic tree: an arena of [`Node`]s with one root.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Tree {
    /// All nodes; index = [`NodeId`].
    nodes: Vec<Node>,
    /// Root node id. Always a valid index once the tree is non-empty.
    root: NodeId,
    /// Caller-declared rooting. `false` means the stored root is an
    /// arbitrary computational anchor with no biological meaning.
    pub rooted: bool,
}

impl Tree {
    /// Builds a tree from a pre-validated arena and root.
    ///
    /// # Errors
    /// [`PhyloError::InvalidTree`] if `root` is out of range, if the
    /// parent/child links are inconsistent, or if the graph contains a
    /// cycle or is disconnected.
    pub fn new(nodes: Vec<Node>, root: NodeId, rooted: bool) -> Result<Self> {
        if nodes.is_empty() {
            return Err(PhyloError::invalid_tree("tree has no nodes"));
        }
        if root >= nodes.len() {
            return Err(PhyloError::invalid_tree(format!(
                "root id {root} out of range (n = {})",
                nodes.len()
            )));
        }
        let tree = Tree {
            nodes,
            root,
            rooted,
        };
        tree.validate()?;
        Ok(tree)
    }

    /// Creates a single-leaf tree (one labelled node, itself the root).
    pub fn leaf(label: impl Into<String>) -> Self {
        Tree {
            nodes: vec![Node {
                label: Some(label.into()),
                branch_length: None,
                parent: None,
                children: Vec::new(),
            }],
            root: 0,
            rooted: true,
        }
    }

    /// Builds a star tree: one internal root with `labels.len()` leaves,
    /// all branch lengths set to `branch_length`.
    ///
    /// # Errors
    /// [`PhyloError::Invalid`] if fewer than two labels are supplied.
    pub fn star(labels: &[String], branch_length: f64) -> Result<Self> {
        if labels.len() < 2 {
            return Err(PhyloError::invalid(
                "labels",
                "a star tree needs at least two leaves",
            ));
        }
        let mut nodes = Vec::with_capacity(labels.len() + 1);
        nodes.push(Node {
            label: None,
            branch_length: None,
            parent: None,
            children: (1..=labels.len()).collect(),
        });
        for name in labels {
            nodes.push(Node {
                label: Some(name.clone()),
                branch_length: Some(branch_length),
                parent: Some(0),
                children: Vec::new(),
            });
        }
        Tree::new(nodes, 0, false)
    }

    /// Checks that parent/child links are mutually consistent, that the
    /// graph is acyclic, and that every node is reachable from the root.
    ///
    /// # Errors
    /// [`PhyloError::InvalidTree`] on the first inconsistency found.
    pub fn validate(&self) -> Result<()> {
        let n = self.nodes.len();
        // Parent of root is None; every other node's declared parent
        // must list it as a child.
        for (id, node) in self.nodes.iter().enumerate() {
            match node.parent {
                None => {
                    if id != self.root {
                        return Err(PhyloError::invalid_tree(format!(
                            "node {id} has no parent but is not the root"
                        )));
                    }
                }
                Some(p) => {
                    if p >= n {
                        return Err(PhyloError::invalid_tree(format!(
                            "node {id} parent {p} out of range"
                        )));
                    }
                    if !self.nodes[p].children.contains(&id) {
                        return Err(PhyloError::invalid_tree(format!(
                            "node {id} parent {p} does not list it as a child"
                        )));
                    }
                }
            }
            for &c in &node.children {
                if c >= n {
                    return Err(PhyloError::invalid_tree(format!(
                        "node {id} child {c} out of range"
                    )));
                }
                if self.nodes[c].parent != Some(id) {
                    return Err(PhyloError::invalid_tree(format!(
                        "node {id} child {c} does not point back to it"
                    )));
                }
            }
        }
        // Reachability + acyclicity: a DFS from the root must visit
        // every node exactly once.
        let mut seen = vec![false; n];
        let mut stack = vec![self.root];
        while let Some(id) = stack.pop() {
            if seen[id] {
                return Err(PhyloError::invalid_tree(format!(
                    "cycle through node {id}"
                )));
            }
            seen[id] = true;
            stack.extend(self.nodes[id].children.iter().copied());
        }
        if let Some(orphan) = seen.iter().position(|&v| !v) {
            return Err(PhyloError::invalid_tree(format!(
                "node {orphan} is not reachable from the root"
            )));
        }
        Ok(())
    }

    /// Number of nodes (leaves + internal).
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// `true` if the tree has a single node only.
    pub fn is_trivial(&self) -> bool {
        self.nodes.len() <= 1
    }

    /// Root node id.
    pub fn root(&self) -> NodeId {
        self.root
    }

    /// Immutable access to a node.
    ///
    /// # Panics
    /// If `id` is out of range — `id`s obtained from this tree's own
    /// accessors are always valid.
    pub fn node(&self, id: NodeId) -> &Node {
        &self.nodes[id]
    }

    /// Mutable access to a node (e.g. to set a branch length).
    ///
    /// Editing `parent` / `children` directly can corrupt the arena —
    /// call [`Tree::validate`] afterwards, or prefer the structural
    /// methods in [`crate::compare`].
    ///
    /// # Panics
    /// If `id` is out of range.
    pub fn node_mut(&mut self, id: NodeId) -> &mut Node {
        &mut self.nodes[id]
    }

    /// All nodes, in arena order.
    pub fn nodes(&self) -> &[Node] {
        &self.nodes
    }

    /// Replaces the root flag (does not re-root — see
    /// [`crate::compare::reroot`]).
    pub fn set_rooted(&mut self, rooted: bool) {
        self.rooted = rooted;
    }

    /// Ids of all leaf nodes, in ascending id order.
    pub fn leaves(&self) -> Vec<NodeId> {
        (0..self.nodes.len())
            .filter(|&i| self.nodes[i].is_leaf())
            .collect()
    }

    /// Number of leaves (terminal taxa).
    pub fn leaf_count(&self) -> usize {
        self.nodes.iter().filter(|n| n.is_leaf()).count()
    }

    /// Number of internal (non-leaf) nodes — the complement of
    /// [`leaf_count`](Self::leaf_count); together they sum to
    /// [`node_count`](Self::node_count).
    pub fn internal_count(&self) -> usize {
        self.nodes.iter().filter(|n| n.is_internal()).count()
    }

    /// Number of **cherries** — internal nodes whose children are exactly two leaves
    /// (a pair of sibling tips). A tree-balance statistic, bounded above by
    /// `leaf_count / 2`.
    pub fn cherry_count(&self) -> usize {
        self.nodes
            .iter()
            .filter(|n| {
                n.children.len() == 2 && n.children.iter().all(|&c| self.nodes[c].is_leaf())
            })
            .count()
    }

    /// Sorted list of leaf labels. Unlabelled leaves are skipped.
    pub fn leaf_labels(&self) -> Vec<String> {
        let mut v: Vec<String> = self
            .nodes
            .iter()
            .filter(|n| n.is_leaf())
            .filter_map(|n| n.label.clone())
            .collect();
        v.sort();
        v
    }

    /// Finds a node by exact label match. Returns the first match in
    /// arena order, or `None`.
    pub fn find(&self, label: &str) -> Option<NodeId> {
        self.nodes
            .iter()
            .position(|n| n.label.as_deref() == Some(label))
    }

    /// Preorder (root-first / NLR) traversal of node ids.
    ///
    /// Children are visited in their stored left-to-right order. This
    /// is the order used to fill top-down scratch arrays.
    pub fn preorder(&self) -> Vec<NodeId> {
        let mut out = Vec::with_capacity(self.nodes.len());
        let mut stack = vec![self.root];
        while let Some(id) = stack.pop() {
            out.push(id);
            // Push children reversed so the leftmost is popped first.
            for &c in self.nodes[id].children.iter().rev() {
                stack.push(c);
            }
        }
        out
    }

    /// Postorder (children-first / LRN) traversal of node ids.
    ///
    /// Every node appears after all of its descendants — the order
    /// required by Fitch parsimony and Felsenstein pruning.
    pub fn postorder(&self) -> Vec<NodeId> {
        // Reverse of a "root, then children" preorder where children
        // are pushed in natural order yields a valid postorder.
        let mut tmp = Vec::with_capacity(self.nodes.len());
        let mut stack = vec![self.root];
        while let Some(id) = stack.pop() {
            tmp.push(id);
            for &c in &self.nodes[id].children {
                stack.push(c);
            }
        }
        tmp.reverse();
        tmp
    }

    /// Path of node ids from `id` up to and including the root.
    pub fn path_to_root(&self, mut id: NodeId) -> Vec<NodeId> {
        let mut path = vec![id];
        while let Some(p) = self.nodes[id].parent {
            path.push(p);
            id = p;
        }
        path
    }

    /// Lowest common ancestor of two nodes.
    ///
    /// Both ids must belong to this tree; the LCA always exists (the
    /// root, at worst).
    pub fn lca(&self, a: NodeId, b: NodeId) -> NodeId {
        let path_a: std::collections::HashSet<NodeId> =
            self.path_to_root(a).into_iter().collect();
        let mut x = b;
        loop {
            if path_a.contains(&x) {
                return x;
            }
            match self.nodes[x].parent {
                Some(p) => x = p,
                None => return self.root,
            }
        }
    }

    /// Patristic distance between two nodes: the sum of branch lengths
    /// on the path connecting them. Missing branch lengths count as 0.
    pub fn patristic_distance(&self, a: NodeId, b: NodeId) -> f64 {
        let lca = self.lca(a, b);
        let up = |mut x: NodeId| -> f64 {
            let mut d = 0.0;
            while x != lca {
                d += self.nodes[x].branch_length.unwrap_or(0.0);
                x = self.nodes[x].parent.expect("path interrupted before lca");
            }
            d
        };
        up(a) + up(b)
    }

    /// Total tree length: the sum of every edge's branch length.
    pub fn total_length(&self) -> f64 {
        self.nodes
            .iter()
            .filter_map(|n| n.branch_length)
            .sum()
    }

    /// Sum of all leaves' depths in edge count (an unweighted height).
    pub fn max_depth(&self) -> usize {
        let mut depth = vec![0usize; self.nodes.len()];
        for &id in &self.preorder() {
            if let Some(p) = self.nodes[id].parent {
                depth[id] = depth[p] + 1;
            }
        }
        depth.into_iter().max().unwrap_or(0)
    }

    /// **Sackin index** — the sum of every leaf's depth (in edges from the root). A
    /// tree-balance statistic: lower is more balanced. Distinct from
    /// [`max_depth`](Self::max_depth) (which is the *maximum* leaf depth, not the *sum*);
    /// a star tree gives `sackin = leaf_count` (all leaves at depth 1).
    pub fn sackin_index(&self) -> usize {
        let mut depth = vec![0usize; self.nodes.len()];
        for &id in &self.preorder() {
            if let Some(p) = self.nodes[id].parent {
                depth[id] = depth[p] + 1;
            }
        }
        self.leaves().iter().map(|&leaf_id| depth[leaf_id]).sum()
    }

    /// `true` if the tree is **binary** (bifurcating) — every internal node has
    /// exactly two children. Leaves and a single-node tree are vacuously binary.
    /// Cannot be inferred from the leaf/internal counts alone (a 4-leaf star and a
    /// caterpillar tree can share counts yet differ in this property).
    pub fn is_binary(&self) -> bool {
        self.nodes
            .iter()
            .all(|n| n.is_leaf() || n.children.len() == 2)
    }

    /// Descendant leaf ids of `id` (just `[id]` if `id` is itself a
    /// leaf).
    pub fn descendant_leaves(&self, id: NodeId) -> Vec<NodeId> {
        let mut out = Vec::new();
        let mut stack = vec![id];
        while let Some(x) = stack.pop() {
            if self.nodes[x].is_leaf() {
                out.push(x);
            } else {
                stack.extend(self.nodes[x].children.iter().copied());
            }
        }
        out
    }

    /// Adds a node to the arena and returns its id. The caller is
    /// responsible for wiring `parent`/`children` consistently; this is
    /// a low-level builder primitive used by the I/O and inference
    /// modules.
    pub(crate) fn push_node(&mut self, node: Node) -> NodeId {
        self.nodes.push(node);
        self.nodes.len() - 1
    }

    /// Constructs an *empty-arena* tree for incremental building. The
    /// returned tree is **not** valid until nodes are added and
    /// [`Tree::finish_building`] is called.
    pub(crate) fn building() -> Tree {
        Tree {
            nodes: Vec::new(),
            root: 0,
            rooted: true,
        }
    }

    /// Finalises an incrementally-built tree: sets the root and runs
    /// [`Tree::validate`].
    pub(crate) fn finish_building(mut self, root: NodeId, rooted: bool) -> Result<Tree> {
        self.root = root;
        self.rooted = rooted;
        self.validate()?;
        Ok(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Builds the rooted tree `((A,B),C);` — two cherries' worth of
    /// nodes used across several tests.
    pub(crate) fn sample_tree() -> Tree {
        let nodes = vec![
            // 0: root
            Node {
                label: None,
                branch_length: None,
                parent: None,
                children: vec![1, 4],
            },
            // 1: internal (A,B)
            Node {
                label: None,
                branch_length: Some(0.3),
                parent: Some(0),
                children: vec![2, 3],
            },
            // 2: A
            Node {
                label: Some("A".into()),
                branch_length: Some(0.1),
                parent: Some(1),
                children: vec![],
            },
            // 3: B
            Node {
                label: Some("B".into()),
                branch_length: Some(0.2),
                parent: Some(1),
                children: vec![],
            },
            // 4: C
            Node {
                label: Some("C".into()),
                branch_length: Some(0.5),
                parent: Some(0),
                children: vec![],
            },
        ];
        Tree::new(nodes, 0, true).unwrap()
    }

    #[test]
    fn validate_accepts_a_good_tree() {
        let t = sample_tree();
        assert!(t.validate().is_ok());
        assert_eq!(t.leaf_count(), 3);
        assert_eq!(t.node_count(), 5);
    }

    #[test]
    fn internal_count_complements_leaf_count() {
        let t = sample_tree();
        // The sample tree has 5 nodes and 3 leaves → 2 internal nodes.
        assert_eq!(t.internal_count(), 2);
        // Non-tautological invariant: every node is exactly one of leaf or internal.
        assert_eq!(t.leaf_count() + t.internal_count(), t.node_count());
    }

    #[test]
    fn cherry_count_on_sample_tree() {
        // sample tree ((A,B),C): node 1 (children A,B) is the only cherry.
        let t = sample_tree();
        assert_eq!(t.cherry_count(), 1);
        // Non-tautological bounds: a cherry consumes 2 leaves and is itself internal.
        assert!(t.cherry_count() <= t.leaf_count() / 2);
        assert!(t.cherry_count() <= t.internal_count());
    }

    #[test]
    fn is_binary_detects_bifurcation() {
        // sample tree ((A,B),C): every internal node has two children → binary.
        let t = sample_tree();
        assert!(t.is_binary());
        // For a rooted binary tree, internal_count == leaf_count − 1 (non-tautological thread).
        assert_eq!(t.internal_count(), t.leaf_count() - 1);
        // A 4-leaf star (root with four children) is NOT binary.
        let labels: Vec<String> = ["w", "x", "y", "z"].iter().map(|s| s.to_string()).collect();
        assert!(!Tree::star(&labels, 1.0).unwrap().is_binary());
        // A single leaf is vacuously binary (no internal nodes).
        assert!(Tree::leaf("X").is_binary());
    }

    #[test]
    fn sackin_index_sums_leaf_depths() {
        // sample tree ((A,B),C): leaves A,B at depth 2 and C at depth 1 → Sackin = 5.
        let t = sample_tree();
        assert_eq!(t.sackin_index(), 5);
        // Distinct from max_depth (the MAX leaf depth, here 2) — Sackin is the SUM.
        assert_eq!(t.max_depth(), 2);
        // A 4-leaf star has every leaf at depth 1 → Sackin = leaf_count.
        let labels: Vec<String> = ["w", "x", "y", "z"].iter().map(|s| s.to_string()).collect();
        let star = Tree::star(&labels, 1.0).unwrap();
        assert_eq!(star.sackin_index(), star.leaf_count());
        // Non-tautological: a non-trivial tree's Sackin ≥ leaf_count; a lone leaf is 0.
        assert!(t.sackin_index() >= t.leaf_count());
        assert_eq!(Tree::leaf("X").sackin_index(), 0);
    }

    #[test]
    fn validate_rejects_a_dangling_child() {
        let mut nodes = sample_tree().nodes;
        nodes[0].children.push(99);
        assert!(Tree::new(nodes, 0, true).is_err());
    }

    #[test]
    fn traversals_have_the_right_shape() {
        let t = sample_tree();
        let pre = t.preorder();
        let post = t.postorder();
        assert_eq!(pre.len(), 5);
        assert_eq!(post.len(), 5);
        // Preorder: root first.
        assert_eq!(pre[0], 0);
        // Postorder: root last, each child before its parent.
        assert_eq!(*post.last().unwrap(), 0);
        let pos = |id: NodeId| post.iter().position(|&x| x == id).unwrap();
        assert!(pos(2) < pos(1) && pos(3) < pos(1) && pos(1) < pos(0));
    }

    #[test]
    fn lca_and_distances() {
        let t = sample_tree();
        let a = t.find("A").unwrap();
        let b = t.find("B").unwrap();
        let c = t.find("C").unwrap();
        assert_eq!(t.lca(a, b), 1);
        assert_eq!(t.lca(a, c), 0);
        // A..B = 0.1 + 0.2
        assert!((t.patristic_distance(a, b) - 0.3).abs() < 1e-12);
        // A..C = 0.1 + 0.3 + 0.5
        assert!((t.patristic_distance(a, c) - 0.9).abs() < 1e-12);
        assert!((t.total_length() - 1.1).abs() < 1e-12);
    }

    #[test]
    fn star_tree_is_valid_and_flat() {
        let labels: Vec<String> = ["w", "x", "y", "z"].iter().map(|s| s.to_string()).collect();
        let t = Tree::star(&labels, 1.0).unwrap();
        assert_eq!(t.leaf_count(), 4);
        assert_eq!(t.max_depth(), 1);
        assert!(!t.rooted);
    }

    #[test]
    fn descendant_leaves_collects_a_clade() {
        let t = sample_tree();
        let mut d = t.descendant_leaves(1);
        d.sort();
        assert_eq!(d, vec![2, 3]);
    }
}
