//! Structural tree editing.
//!
//! - [`reroot`] re-roots a tree on the edge above a chosen node — the
//!   classic "outgroup rooting".
//! - [`midpoint_root`] roots at the midpoint of the longest leaf-to-leaf
//!   path (the standard rooting when no outgroup is known).
//! - [`ladderize`] reorders children so subtrees sort by size — purely
//!   cosmetic, the FigTree "increasing node order" view.
//! - [`prune_taxa`] removes a set of leaves and suppresses the
//!   resulting degree-2 nodes.
//! - [`subtree`] extracts the clade rooted at a node as a standalone
//!   tree.
//!
//! Re-rooting rebuilds the arena: the chosen edge gets a fresh root
//! node, and every parent pointer between the old and new root is
//! reversed.

use crate::error::{PhyloError, Result};
use crate::tree::{Node, NodeId, Tree};
use std::collections::HashSet;

/// Re-roots `tree` on the edge above the node labelled `outgroup`.
///
/// A new root is inserted in the middle of that edge; if the node had a
/// branch length it is split evenly between the two new edges. The
/// result is flagged [`Tree::rooted`]` = true`.
///
/// # Errors
/// [`PhyloError::Invalid`] if no leaf carries `outgroup`, or it is
/// already the (only) child of the root.
pub fn reroot(tree: &Tree, outgroup: &str) -> Result<Tree> {
    let target = tree
        .find(outgroup)
        .ok_or_else(|| PhyloError::invalid("outgroup", format!("no node `{outgroup}`")))?;
    reroot_on_node(tree, target)
}

/// Re-roots on the edge above `target` (by node id).
fn reroot_on_node(tree: &Tree, target: NodeId) -> Result<Tree> {
    reroot_on_edge(tree, target, None)
}

/// Reroots `tree` on the edge above `target`. `dist_from_target`, when
/// `Some`, places the new root that far from `target` along the edge
/// (the remainder of the edge length goes above the root, toward the old
/// parent); `None` places it at the edge midpoint. This is what
/// midpoint-rooting needs — the root can sit at an arbitrary point along
/// an edge, not only at its half.
fn reroot_on_edge(
    tree: &Tree,
    target: NodeId,
    dist_from_target: Option<f64>,
) -> Result<Tree> {
    let old_root = tree.root();
    if target == old_root {
        return Err(PhyloError::invalid(
            "outgroup",
            "node is already the root",
        ));
    }
    let parent = tree
        .node(target)
        .parent
        .ok_or_else(|| PhyloError::invalid("outgroup", "node has no parent"))?;

    // Rebuild: a fresh root with two children — `target` and the rest
    // of the tree (the old graph with edges above `target` reversed).
    let mut nodes: Vec<Node> = tree
        .nodes()
        .iter()
        .map(|n| Node {
            label: n.label.clone(),
            branch_length: n.branch_length,
            parent: None,
            children: Vec::new(),
        })
        .collect();

    // Reproduce the original adjacency as an undirected graph.
    let mut adj: Vec<Vec<NodeId>> = vec![Vec::new(); tree.node_count()];
    for (id, node) in tree.nodes().iter().enumerate() {
        for &c in &node.children {
            adj[id].push(c);
            adj[c].push(id);
        }
    }
    // Branch length is stored on the child; index by the unordered
    // edge (min,max).
    let edge_len = |a: NodeId, b: NodeId| -> Option<f64> {
        if tree.node(a).parent == Some(b) {
            tree.node(a).branch_length
        } else if tree.node(b).parent == Some(a) {
            tree.node(b).branch_length
        } else {
            None
        }
    };

    let new_root = nodes.len();
    // Split the target↔parent edge at `dist_from_target` (clamped into
    // the edge), defaulting to the midpoint. `below` is the new branch
    // length from the root down to `target`; `above` goes up to the old
    // parent — together they preserve the original edge length.
    let edge = edge_len(target, parent).unwrap_or(0.0);
    let below = match dist_from_target {
        Some(d) => d.clamp(0.0, edge),
        None => edge / 2.0,
    };
    let above = (edge - below).max(0.0);
    nodes.push(Node {
        label: None,
        branch_length: None,
        parent: None,
        children: vec![target, parent],
    });

    // Orient the tree away from `new_root` by a DFS over `adj`, but
    // route the `target`-`parent` edge through `new_root`.
    let mut visited: HashSet<NodeId> = HashSet::new();
    visited.insert(new_root);
    nodes[target].parent = Some(new_root);
    nodes[target].branch_length = Some(below);
    nodes[parent].parent = Some(new_root);
    nodes[parent].branch_length = Some(above);

    // BFS from `target` and `parent`, excluding the direct
    // target<->parent link (handled by the new root).
    let mut stack: Vec<NodeId> = vec![target, parent];
    visited.insert(target);
    visited.insert(parent);
    while let Some(cur) = stack.pop() {
        for &nb in &adj[cur] {
            // Skip the original edge that the new root replaced.
            if (cur == target && nb == parent) || (cur == parent && nb == target)
            {
                continue;
            }
            if visited.insert(nb) {
                nodes[nb].parent = Some(cur);
                nodes[nb].branch_length = edge_len(cur, nb);
                nodes[cur].children.push(nb);
                stack.push(nb);
            }
        }
    }

    // The old root may now be a redundant degree-2 node — suppress it.
    let mut tree = Tree::new(nodes, new_root, true)?;
    suppress_unifurcations(&mut tree)?;
    Ok(tree)
}

/// Roots `tree` at the midpoint of its longest leaf-to-leaf path.
///
/// # Errors
/// [`PhyloError::Invalid`] if the tree has fewer than two leaves.
pub fn midpoint_root(tree: &Tree) -> Result<Tree> {
    let leaves = tree.leaves();
    if leaves.len() < 2 {
        return Err(PhyloError::invalid("tree", "need at least two leaves"));
    }
    // Find the most-distant leaf pair (the tree's "diameter").
    let (mut best_a, mut best_b, mut best_d) = (leaves[0], leaves[1], -1.0);
    for (i, &a) in leaves.iter().enumerate() {
        for &b in &leaves[i + 1..] {
            let d = tree.patristic_distance(a, b);
            if d > best_d {
                best_d = d;
                best_a = a;
                best_b = b;
            }
        }
    }
    if best_d <= 0.0 {
        // All branch lengths zero — fall back to outgroup-rooting on
        // the first leaf.
        return reroot_on_node(tree, leaves[0]);
    }
    let half = best_d / 2.0;
    // Walk the path from `best_a` toward `best_b`, accumulating the
    // length of each consecutive edge, until the cumulative distance
    // first reaches the midpoint. The diameter midpoint generally falls
    // PART-WAY along an edge, so the tree is rerooted at that exact
    // position — not merely on the nearest node.
    let path = path_between(tree, best_a, best_b);
    let mut acc = 0.0;
    for w in path.windows(2) {
        let (prev, node) = (w[0], w[1]);
        // The edge between two consecutive path nodes — its length is
        // stored on whichever node is the child of the other.
        let (child, bl) = if tree.node(node).parent == Some(prev) {
            (node, tree.node(node).branch_length.unwrap_or(0.0))
        } else {
            (prev, tree.node(prev).branch_length.unwrap_or(0.0))
        };
        if acc + bl >= half {
            // The midpoint is on this edge. `child` carries the edge;
            // place the root at the midpoint, measured from `child`.
            // Distance from `child` to the midpoint depends on whether
            // the walk crosses the edge from `child` first (the edge is
            // `prev -> node` and stored on `child`).
            let into_edge = half - acc; // distance into the edge from `prev`
            let from_child = if child == node {
                // walking prev -> node: `node` is the far end of the
                // edge, the root sits `bl - into_edge` above `node`.
                (bl - into_edge).clamp(0.0, bl)
            } else {
                // edge stored on `prev` (`prev` is `node`'s child in
                // the original tree): the root sits `into_edge` above
                // `prev`.
                into_edge.clamp(0.0, bl)
            };
            return reroot_on_edge(tree, child, Some(from_child));
        }
        acc += bl;
    }
    // Numerical edge case — root on the far leaf.
    reroot_on_node(tree, best_b)
}

/// Node-id path between two nodes, inclusive of both endpoints.
fn path_between(tree: &Tree, a: NodeId, b: NodeId) -> Vec<NodeId> {
    let lca = tree.lca(a, b);
    let mut up_a = Vec::new();
    let mut x = a;
    while x != lca {
        up_a.push(x);
        x = tree.node(x).parent.expect("path to lca");
    }
    up_a.push(lca);
    let mut up_b = Vec::new();
    let mut y = b;
    while y != lca {
        up_b.push(y);
        y = tree.node(y).parent.expect("path to lca");
    }
    up_b.reverse();
    up_a.extend(up_b);
    up_a
}

/// Reorders every node's children so that subtrees sort by leaf count
/// (smallest first). Cosmetic only — topology is unchanged.
pub fn ladderize(tree: &Tree) -> Tree {
    let mut t = tree.clone();
    // Precompute the leaf count under each node.
    let mut size = vec![1usize; t.node_count()];
    for &id in &t.postorder() {
        if t.node(id).is_internal() {
            size[id] = t.node(id).children.iter().map(|&c| size[c]).sum();
        }
    }
    for id in 0..t.node_count() {
        if t.node(id).is_internal() {
            let mut kids = t.node(id).children.clone();
            kids.sort_by_key(|&c| size[c]);
            t.node_mut(id).children = kids;
        }
    }
    t
}

/// Removes the leaves whose labels appear in `taxa` and suppresses any
/// degree-2 nodes the removal creates.
///
/// # Errors
/// [`PhyloError::Invalid`] if pruning would leave fewer than two
/// leaves; [`PhyloError::InvalidTree`] if the rebuilt arena is
/// inconsistent.
pub fn prune_taxa(tree: &Tree, taxa: &[String]) -> Result<Tree> {
    let drop: HashSet<&str> = taxa.iter().map(|s| s.as_str()).collect();
    let keep_leaves: Vec<NodeId> = tree
        .leaves()
        .into_iter()
        .filter(|&l| {
            tree.node(l)
                .label
                .as_deref()
                .map(|n| !drop.contains(n))
                .unwrap_or(true)
        })
        .collect();
    if keep_leaves.len() < 2 {
        return Err(PhyloError::invalid(
            "taxa",
            "pruning would leave fewer than two leaves",
        ));
    }
    // Mark every node on a kept-leaf-to-root path as retained.
    let mut retained: HashSet<NodeId> = HashSet::new();
    for &leaf in &keep_leaves {
        for node in tree.path_to_root(leaf) {
            retained.insert(node);
        }
    }
    // Rebuild keeping only retained nodes, then suppress unifurcations.
    let rebuilt = rebuild_retained(tree, &retained)?;
    Ok(rebuilt)
}

/// Extracts the clade rooted at the node labelled `clade_root` as an
/// independent [`Tree`].
///
/// # Errors
/// [`PhyloError::Invalid`] if no node carries `clade_root`.
pub fn subtree(tree: &Tree, clade_root: &str) -> Result<Tree> {
    let root = tree
        .find(clade_root)
        .ok_or_else(|| PhyloError::invalid("clade_root", "node not found"))?;
    let keep: HashSet<NodeId> = tree.descendant_set(root);
    let mut retained = keep.clone();
    retained.insert(root);
    // Rebuild with `root` as the new root.
    let mut id_map = std::collections::HashMap::new();
    let mut nodes = Vec::new();
    // Preorder over the subtree.
    let mut stack = vec![root];
    let mut order = Vec::new();
    while let Some(x) = stack.pop() {
        order.push(x);
        for &c in tree.node(x).children.iter().rev() {
            stack.push(c);
        }
    }
    for &old in &order {
        let new = nodes.len();
        id_map.insert(old, new);
        nodes.push(Node {
            label: tree.node(old).label.clone(),
            branch_length: tree.node(old).branch_length,
            parent: None,
            children: Vec::new(),
        });
    }
    for &old in &order {
        let new = id_map[&old];
        if let Some(p) = tree.node(old).parent {
            if let Some(&np) = id_map.get(&p) {
                nodes[new].parent = Some(np);
            }
        }
        for &c in &tree.node(old).children {
            if let Some(&nc) = id_map.get(&c) {
                nodes[new].children.push(nc);
            }
        }
    }
    nodes[0].parent = None;
    nodes[0].branch_length = None;
    Tree::new(nodes, 0, tree.rooted)
}

/// Rebuilds a tree keeping only the nodes in `retained`, then
/// suppresses degree-2 nodes.
fn rebuild_retained(tree: &Tree, retained: &HashSet<NodeId>) -> Result<Tree> {
    let mut id_map = std::collections::HashMap::new();
    let mut nodes = Vec::new();
    for old in 0..tree.node_count() {
        if retained.contains(&old) {
            let new = nodes.len();
            id_map.insert(old, new);
            nodes.push(Node {
                label: tree.node(old).label.clone(),
                branch_length: tree.node(old).branch_length,
                parent: None,
                children: Vec::new(),
            });
        }
    }
    for old in 0..tree.node_count() {
        if let Some(&new) = id_map.get(&old) {
            if let Some(p) = tree.node(old).parent {
                if let Some(&np) = id_map.get(&p) {
                    nodes[new].parent = Some(np);
                }
            }
            for &c in &tree.node(old).children {
                if let Some(&nc) = id_map.get(&c) {
                    nodes[new].children.push(nc);
                }
            }
        }
    }
    let new_root = id_map[&tree.root()];
    nodes[new_root].parent = None;
    nodes[new_root].branch_length = None;
    let mut t = Tree::new(nodes, new_root, tree.rooted)?;
    suppress_unifurcations(&mut t)?;
    Ok(t)
}

/// Collapses every degree-2 (one-child) node, merging its incoming and
/// outgoing branch lengths. Rebuilds the arena because node ids shift.
fn suppress_unifurcations(tree: &mut Tree) -> Result<()> {
    loop {
        // Find a non-root node with exactly one child.
        let victim = (0..tree.node_count()).find(|&id| {
            id != tree.root() && tree.node(id).children.len() == 1
        });
        let Some(v) = victim else { break };
        let parent = tree.node(v).parent.expect("non-root has parent");
        let child = tree.node(v).children[0];
        let merged_len = tree.node(v).branch_length.unwrap_or(0.0)
            + tree.node(child).branch_length.unwrap_or(0.0);
        // Re-point: parent -> child directly.
        let slot = tree
            .node(parent)
            .children
            .iter()
            .position(|&c| c == v)
            .expect("parent lists v");
        tree.node_mut(parent).children[slot] = child;
        tree.node_mut(child).parent = Some(parent);
        tree.node_mut(child).branch_length = Some(merged_len);
        // `v` is now detached; rebuild to compact the arena.
        *tree = compact(tree);
    }
    // A degree-1 root (after suppression) is also collapsed.
    if tree.node(tree.root()).children.len() == 1 {
        let only = tree.node(tree.root()).children[0];
        let mut t = Tree::building();
        // Re-root at `only`.
        let mut id_map = std::collections::HashMap::new();
        let keep: HashSet<NodeId> = {
            let mut s = tree.descendant_set(only);
            s.insert(only);
            s
        };
        let mut order = Vec::new();
        let mut stack = vec![only];
        while let Some(x) = stack.pop() {
            order.push(x);
            for &c in tree.node(x).children.iter().rev() {
                stack.push(c);
            }
        }
        for &old in &order {
            if !keep.contains(&old) {
                continue;
            }
            let new = t.push_node(Node {
                label: tree.node(old).label.clone(),
                branch_length: tree.node(old).branch_length,
                parent: None,
                children: Vec::new(),
            });
            id_map.insert(old, new);
        }
        for &old in &order {
            let Some(&new) = id_map.get(&old) else {
                continue;
            };
            if let Some(p) = tree.node(old).parent {
                if let Some(&np) = id_map.get(&p) {
                    t.node_mut(new).parent = Some(np);
                }
            }
            for &c in &tree.node(old).children {
                if let Some(&nc) = id_map.get(&c) {
                    t.node_mut(new).children.push(nc);
                }
            }
        }
        let nr = id_map[&only];
        t.node_mut(nr).parent = None;
        t.node_mut(nr).branch_length = None;
        *tree = t.finish_building(nr, tree.rooted)?;
    }
    Ok(())
}

/// Rebuilds `tree` into a fresh compact arena (drops detached nodes).
fn compact(tree: &Tree) -> Tree {
    let mut id_map = std::collections::HashMap::new();
    let mut nodes = Vec::new();
    // Reachable set from the root.
    let mut reach = HashSet::new();
    let mut stack = vec![tree.root()];
    while let Some(x) = stack.pop() {
        if reach.insert(x) {
            stack.extend(tree.node(x).children.iter().copied());
        }
    }
    for old in 0..tree.node_count() {
        if reach.contains(&old) {
            let new = nodes.len();
            id_map.insert(old, new);
            nodes.push(Node {
                label: tree.node(old).label.clone(),
                branch_length: tree.node(old).branch_length,
                parent: None,
                children: Vec::new(),
            });
        }
    }
    for old in 0..tree.node_count() {
        if let Some(&new) = id_map.get(&old) {
            if let Some(p) = tree.node(old).parent {
                if let Some(&np) = id_map.get(&p) {
                    nodes[new].parent = Some(np);
                }
            }
            for &c in &tree.node(old).children {
                if let Some(&nc) = id_map.get(&c) {
                    nodes[new].children.push(nc);
                }
            }
        }
    }
    let new_root = id_map[&tree.root()];
    Tree::new(nodes, new_root, tree.rooted).expect("compaction preserves validity")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::newick::read_newick;

    #[test]
    fn reroot_changes_the_root_but_keeps_leaves() {
        let t = read_newick("((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);").unwrap();
        let r = reroot(&t, "A").unwrap();
        assert_eq!(r.leaf_labels(), t.leaf_labels());
        assert!(r.validate().is_ok());
        // A is now adjacent to the root.
        let a = r.find("A").unwrap();
        assert_eq!(r.node(a).parent, Some(r.root()));
    }

    #[test]
    fn reroot_rejects_unknown_outgroup() {
        let t = read_newick("((A,B),(C,D));").unwrap();
        assert!(reroot(&t, "Z").is_err());
    }

    #[test]
    fn midpoint_root_balances_the_diameter() {
        // A long A branch makes the midpoint fall on A's side.
        let t = read_newick("((A:2.0,B:0.1):0.1,(C:0.1,D:0.1):0.1);").unwrap();
        let r = midpoint_root(&t).unwrap();
        assert!(r.validate().is_ok());
        assert_eq!(r.leaf_count(), 4);
        // Root's two subtrees should have nearly equal height.
        let root = r.root();
        let children = &r.node(root).children;
        let h: Vec<f64> = children
            .iter()
            .map(|&c| {
                r.descendant_leaves(c)
                    .iter()
                    .map(|&l| r.patristic_distance(c, l))
                    .fold(0.0, f64::max)
                    + r.node(c).branch_length.unwrap_or(0.0)
            })
            .collect();
        assert!((h[0] - h[1]).abs() < 0.3, "unbalanced: {h:?}");
    }

    #[test]
    fn ladderize_sorts_children_by_size() {
        // ((A,B),C): the (A,B) clade (2 leaves) and C (1 leaf).
        let t = read_newick("((A,B),C);").unwrap();
        let l = ladderize(&t);
        let root = l.root();
        let sizes: Vec<usize> = l
            .node(root)
            .children
            .iter()
            .map(|&c| l.descendant_leaves(c).len())
            .collect();
        // Sorted ascending.
        assert!(sizes.windows(2).all(|w| w[0] <= w[1]), "{sizes:?}");
    }

    #[test]
    fn prune_removes_taxa_and_suppresses_nodes() {
        let t = read_newick("(((A,B),C),(D,E));").unwrap();
        let p = prune_taxa(&t, &["B".to_string()]).unwrap();
        assert_eq!(p.leaf_count(), 4);
        assert!(p.find("B").is_none());
        assert!(p.find("A").is_some());
        assert!(p.validate().is_ok());
        // No degree-2 nodes remain.
        for id in 0..p.node_count() {
            if id != p.root() {
                assert_ne!(p.node(id).children.len(), 1, "unifurcation left");
            }
        }
    }

    #[test]
    fn prune_rejects_emptying_the_tree() {
        let t = read_newick("(A,B);").unwrap();
        assert!(prune_taxa(&t, &["A".to_string()]).is_err());
    }

    #[test]
    fn subtree_extracts_a_clade() {
        let t = read_newick("((A,B)clade1,(C,D)clade2);").unwrap();
        let s = subtree(&t, "clade1").unwrap();
        let mut leaves = s.leaf_labels();
        leaves.sort();
        assert_eq!(leaves, vec!["A", "B"]);
        assert!(s.validate().is_ok());
    }
}
