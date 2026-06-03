//! Study tree — Salome's central data model.
//!
//! A study is an ordered list of [`StudyNode`]s, each carrying:
//! - a numeric `id` (assigned at insertion time),
//! - a kind label,
//! - a vector of upstream dependency ids.
//!
//! Operations support:
//! - **insertion** ([`Study::add_solid`] / `add_mesh` / `add_fem`
//!   / `add_result`),
//! - **dependency query** ([`Study::depends`]),
//! - **rebuild** ([`Study::rebuild`]) — re-run a node and every
//!   downstream node in topological order.

use serde::{Deserialize, Serialize};

use crate::error::SalomeError;

/// Numeric node id (monotonic — never reused even after deletion).
pub type NodeId = u32;

/// Kind tag of a study node. The actual payload lives in the
/// bridged module (geom/mesh/analysis) and is referenced by name —
/// the study keeps only metadata.
#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum NodeKind {
    /// A solid produced by the geom module.
    Solid {
        /// Free-form name of the geom-module object.
        name: String,
    },
    /// A mesh produced by the mesh module.
    Mesh {
        /// Free-form name of the mesh.
        name: String,
    },
    /// An FEM analysis case.
    FemAnalysis {
        /// Free-form analysis name.
        name: String,
    },
    /// A computed result (post-processed field).
    Result {
        /// Free-form result name.
        name: String,
    },
}

/// One study node — id + kind + dependencies + dirty flag.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct StudyNode {
    /// Numeric id.
    pub id: NodeId,
    /// Kind metadata.
    pub kind: NodeKind,
    /// Upstream dependency ids.
    pub dependencies: Vec<NodeId>,
    /// Has been (re)evaluated since last upstream change.
    pub up_to_date: bool,
}

/// The full study.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Study {
    /// All nodes in insertion order.
    pub nodes: Vec<StudyNode>,
    /// Next id to assign.
    next_id: NodeId,
}

impl Study {
    /// Empty study.
    pub fn new() -> Self {
        Self::default()
    }

    fn allocate(&mut self) -> NodeId {
        let id = self.next_id;
        self.next_id += 1;
        id
    }

    /// Append a solid node.
    pub fn add_solid(&mut self, name: impl Into<String>, deps: Vec<NodeId>) -> NodeId {
        let id = self.allocate();
        self.nodes.push(StudyNode {
            id,
            kind: NodeKind::Solid { name: name.into() },
            dependencies: deps,
            up_to_date: false,
        });
        id
    }

    /// Append a mesh node.
    pub fn add_mesh(&mut self, name: impl Into<String>, deps: Vec<NodeId>) -> NodeId {
        let id = self.allocate();
        self.nodes.push(StudyNode {
            id,
            kind: NodeKind::Mesh { name: name.into() },
            dependencies: deps,
            up_to_date: false,
        });
        id
    }

    /// Append an FEM-analysis node.
    pub fn add_fem(&mut self, name: impl Into<String>, deps: Vec<NodeId>) -> NodeId {
        let id = self.allocate();
        self.nodes.push(StudyNode {
            id,
            kind: NodeKind::FemAnalysis { name: name.into() },
            dependencies: deps,
            up_to_date: false,
        });
        id
    }

    /// Append a result node.
    pub fn add_result(&mut self, name: impl Into<String>, deps: Vec<NodeId>) -> NodeId {
        let id = self.allocate();
        self.nodes.push(StudyNode {
            id,
            kind: NodeKind::Result { name: name.into() },
            dependencies: deps,
            up_to_date: false,
        });
        id
    }

    /// Lookup a node by id.
    pub fn node(&self, id: NodeId) -> Result<&StudyNode, SalomeError> {
        self.nodes
            .iter()
            .find(|n| n.id == id)
            .ok_or(SalomeError::UnknownNode(id))
    }

    /// Direct upstream dependencies of `id`.
    pub fn depends(&self, id: NodeId) -> Result<Vec<NodeId>, SalomeError> {
        Ok(self.node(id)?.dependencies.clone())
    }

    /// Rebuild — mark `id` dirty plus every transitive *downstream*
    /// node, then walk them in topological order marking up_to_date.
    ///
    /// Returns the ordered list of node ids that were rebuilt.
    pub fn rebuild(&mut self, id: NodeId) -> Result<Vec<NodeId>, SalomeError> {
        // Validate.
        self.node(id)?;

        // Downstream set — anything that lists `id` (directly or
        // transitively) as a dependency.
        let downstream = self.downstream_of(id);

        // Topologically sort just the dirty set.
        let order = self.topo_sort_subset(&downstream)?;

        for nid in &order {
            // Find the node and flip its flag.
            for n in &mut self.nodes {
                if n.id == *nid {
                    n.up_to_date = true;
                }
            }
        }
        Ok(order)
    }

    fn downstream_of(&self, root: NodeId) -> Vec<NodeId> {
        let mut visited = std::collections::BTreeSet::new();
        visited.insert(root);
        let mut grew = true;
        while grew {
            grew = false;
            for n in &self.nodes {
                if !visited.contains(&n.id)
                    && n.dependencies.iter().any(|d| visited.contains(d))
                {
                    visited.insert(n.id);
                    grew = true;
                }
            }
        }
        visited.into_iter().collect()
    }

    fn topo_sort_subset(&self, ids: &[NodeId]) -> Result<Vec<NodeId>, SalomeError> {
        let id_set: std::collections::BTreeSet<NodeId> = ids.iter().copied().collect();
        let mut in_deg: std::collections::BTreeMap<NodeId, usize> =
            id_set.iter().map(|id| (*id, 0)).collect();
        for n in &self.nodes {
            if id_set.contains(&n.id) {
                for d in &n.dependencies {
                    if id_set.contains(d) {
                        *in_deg.entry(n.id).or_insert(0) += 1;
                    }
                }
            }
        }
        let mut q: std::collections::VecDeque<NodeId> = in_deg
            .iter()
            .filter(|(_, k)| **k == 0)
            .map(|(id, _)| *id)
            .collect();
        let mut out = Vec::new();
        while let Some(id) = q.pop_front() {
            out.push(id);
            for n in &self.nodes {
                if id_set.contains(&n.id) && n.dependencies.contains(&id) {
                    let e = in_deg.entry(n.id).or_insert(0);
                    if *e > 0 {
                        *e -= 1;
                        if *e == 0 {
                            q.push_back(n.id);
                        }
                    }
                }
            }
        }
        if out.len() != id_set.len() {
            // Find first node with a non-zero remaining in-degree.
            for n in &self.nodes {
                if id_set.contains(&n.id) && in_deg.get(&n.id).copied().unwrap_or(0) > 0 {
                    return Err(SalomeError::Cycle(n.id));
                }
            }
        }
        Ok(out)
    }
}
