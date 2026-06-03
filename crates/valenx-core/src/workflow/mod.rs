//! Workflow DAG — the orchestration layer that composes adapter
//! invocations into a reproducible sequence.
//!
//! Spec: [ARCHITECTURE.md § 6](../../ARCHITECTURE.md). A `Workflow`
//! is a set of nodes (adapter calls) connected by typed edges
//! carrying canonical data (`Geometry`, `Mesh`, `Case`, `Results`).
//! The engine topologically sorts and executes; a future concurrent
//! executor can replace the sequential one without changing the
//! data model.

use std::collections::{BTreeMap, BTreeSet, VecDeque};

use serde::{Deserialize, Serialize};
use thiserror::Error;

/// A port name on a workflow node — an edge flows from a producing
/// port on one node into a consuming port on another.
pub type PortName = String;

/// Canonical type carried across an edge. Matches the canonical
/// crates (`valenx-geo`, `valenx-mesh`, `valenx-fields`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum PortType {
    Geometry,
    Mesh,
    Case,
    Results,
    /// Untyped blob — escape hatch for adapter-specific data that
    /// hasn't been canonicalised yet.
    Raw,
}

/// One node = one adapter invocation.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct WorkflowNode {
    pub id: String,
    /// Adapter ID (matches `AdapterInfo::id`).
    pub adapter_id: String,
    /// Ports this node consumes.
    #[serde(default)]
    pub inputs: BTreeMap<PortName, PortType>,
    /// Ports this node produces.
    #[serde(default)]
    pub outputs: BTreeMap<PortName, PortType>,
    /// Free-form config passed to the adapter (serialised TOML
    /// table). Adapters validate against their own schema.
    #[serde(default)]
    pub config: BTreeMap<String, toml::Value>,
}

/// A typed edge from one node's output port to another's input port.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct WorkflowEdge {
    pub from: String,
    pub from_port: PortName,
    pub to: String,
    pub to_port: PortName,
}

/// A workflow — a DAG. Serialised into the project file so runs
/// reproduce.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Workflow {
    #[serde(default)]
    pub nodes: Vec<WorkflowNode>,
    #[serde(default)]
    pub edges: Vec<WorkflowEdge>,
}

/// Errors raised by [`Workflow::validate`].
#[derive(Debug, Error)]
pub enum WorkflowError {
    /// An edge endpoint named a node id that doesn't exist.
    #[error("edge references unknown node {0:?}")]
    UnknownNode(String),

    /// Edge connects two ports with incompatible [`PortType`]s.
    #[error(
        "edge {from}:{from_port} -> {to}:{to_port} type mismatch: \
         {from_type:?} vs {to_type:?}"
    )]
    PortTypeMismatch {
        /// Source node id.
        from: String,
        /// Source port name.
        from_port: String,
        /// Sink node id.
        to: String,
        /// Sink port name.
        to_port: String,
        /// Output port's declared type.
        from_type: PortType,
        /// Input port's declared type.
        to_type: PortType,
    },

    /// Edge referenced a port name not declared on the named node.
    #[error("edge references unknown port {node}:{port}")]
    UnknownPort {
        /// Node id.
        node: String,
        /// Port name.
        port: String,
    },

    /// The directed graph contains at least one cycle; the variant
    /// names one node on that cycle.
    #[error("workflow has a cycle involving node {0:?}")]
    Cycle(String),
}

impl Workflow {
    /// New, empty workflow (no nodes, no edges).
    pub fn new() -> Self {
        Self::default()
    }

    /// Validate node IDs, port names, port types, and acyclicity.
    pub fn validate(&self) -> Result<(), WorkflowError> {
        let ids: BTreeSet<&str> = self.nodes.iter().map(|n| n.id.as_str()).collect();
        let by_id: BTreeMap<&str, &WorkflowNode> =
            self.nodes.iter().map(|n| (n.id.as_str(), n)).collect();

        for edge in &self.edges {
            let from = by_id
                .get(edge.from.as_str())
                .ok_or_else(|| WorkflowError::UnknownNode(edge.from.clone()))?;
            let to = by_id
                .get(edge.to.as_str())
                .ok_or_else(|| WorkflowError::UnknownNode(edge.to.clone()))?;

            let from_type = from.outputs.get(&edge.from_port).copied().ok_or_else(|| {
                WorkflowError::UnknownPort {
                    node: edge.from.clone(),
                    port: edge.from_port.clone(),
                }
            })?;
            let to_type = to.inputs.get(&edge.to_port).copied().ok_or_else(|| {
                WorkflowError::UnknownPort {
                    node: edge.to.clone(),
                    port: edge.to_port.clone(),
                }
            })?;

            if from_type != to_type {
                return Err(WorkflowError::PortTypeMismatch {
                    from: edge.from.clone(),
                    from_port: edge.from_port.clone(),
                    to: edge.to.clone(),
                    to_port: edge.to_port.clone(),
                    from_type,
                    to_type,
                });
            }
        }

        // Topological sort — reject cycles.
        let _ = self.topo_order()?;
        let _ = ids; // keep the binding so future checks can use it
        Ok(())
    }

    /// Topological order of node IDs. Uses Kahn's algorithm.
    pub fn topo_order(&self) -> Result<Vec<String>, WorkflowError> {
        let mut in_degree: BTreeMap<&str, usize> =
            self.nodes.iter().map(|n| (n.id.as_str(), 0)).collect();
        for edge in &self.edges {
            *in_degree.entry(edge.to.as_str()).or_insert(0) += 1;
        }

        let mut queue: VecDeque<&str> = in_degree
            .iter()
            .filter(|(_, d)| **d == 0)
            .map(|(k, _)| *k)
            .collect();

        let mut out: Vec<String> = Vec::new();
        while let Some(id) = queue.pop_front() {
            out.push(id.to_string());
            for edge in self.edges.iter().filter(|e| e.from == id) {
                let d = in_degree.entry(edge.to.as_str()).or_insert(0);
                *d = d.saturating_sub(1);
                if *d == 0 {
                    queue.push_back(edge.to.as_str());
                }
            }
        }

        if out.len() != self.nodes.len() {
            let stuck = self
                .nodes
                .iter()
                .find(|n| !out.contains(&n.id))
                .map(|n| n.id.clone())
                .unwrap_or_else(|| "<unknown>".to_string());
            return Err(WorkflowError::Cycle(stuck));
        }

        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn port_map(pairs: &[(&str, PortType)]) -> BTreeMap<PortName, PortType> {
        pairs.iter().map(|(k, v)| ((*k).to_string(), *v)).collect()
    }

    fn sample_dag() -> Workflow {
        Workflow {
            nodes: vec![
                WorkflowNode {
                    id: "cad".into(),
                    adapter_id: "freecad".into(),
                    inputs: Default::default(),
                    outputs: port_map(&[("geom", PortType::Geometry)]),
                    config: Default::default(),
                },
                WorkflowNode {
                    id: "mesh".into(),
                    adapter_id: "gmsh".into(),
                    inputs: port_map(&[("geom", PortType::Geometry)]),
                    outputs: port_map(&[("mesh", PortType::Mesh)]),
                    config: Default::default(),
                },
                WorkflowNode {
                    id: "solve".into(),
                    adapter_id: "openfoam".into(),
                    inputs: port_map(&[("mesh", PortType::Mesh), ("case", PortType::Case)]),
                    outputs: port_map(&[("results", PortType::Results)]),
                    config: Default::default(),
                },
            ],
            edges: vec![
                WorkflowEdge {
                    from: "cad".into(),
                    from_port: "geom".into(),
                    to: "mesh".into(),
                    to_port: "geom".into(),
                },
                WorkflowEdge {
                    from: "mesh".into(),
                    from_port: "mesh".into(),
                    to: "solve".into(),
                    to_port: "mesh".into(),
                },
            ],
        }
    }

    #[test]
    fn valid_dag_validates() {
        let w = sample_dag();
        w.validate().unwrap();
    }

    #[test]
    fn topo_order_respects_edges() {
        let order = sample_dag().topo_order().unwrap();
        let cad_i = order.iter().position(|n| n == "cad").unwrap();
        let mesh_i = order.iter().position(|n| n == "mesh").unwrap();
        let solve_i = order.iter().position(|n| n == "solve").unwrap();
        assert!(cad_i < mesh_i);
        assert!(mesh_i < solve_i);
    }

    #[test]
    fn unknown_edge_node_fails() {
        let mut w = sample_dag();
        w.edges.push(WorkflowEdge {
            from: "cad".into(),
            from_port: "geom".into(),
            to: "ghost".into(),
            to_port: "x".into(),
        });
        assert!(matches!(
            w.validate(),
            Err(WorkflowError::UnknownNode(id)) if id == "ghost"
        ));
    }

    #[test]
    fn port_type_mismatch_fails() {
        let mut w = sample_dag();
        // Add a bogus mesh → case edge (mesh port outputs Mesh, solve.case expects Case).
        w.edges.push(WorkflowEdge {
            from: "mesh".into(),
            from_port: "mesh".into(),
            to: "solve".into(),
            to_port: "case".into(),
        });
        assert!(matches!(
            w.validate(),
            Err(WorkflowError::PortTypeMismatch { .. })
        ));
    }

    #[test]
    fn cycle_detected() {
        let mut w = sample_dag();
        // solve.results has no consumer; make a cycle by feeding it
        // back into cad.
        w.nodes[0].inputs.insert("back".into(), PortType::Results);
        w.edges.push(WorkflowEdge {
            from: "solve".into(),
            from_port: "results".into(),
            to: "cad".into(),
            to_port: "back".into(),
        });
        assert!(matches!(w.validate(), Err(WorkflowError::Cycle(_))));
    }
}
