//! Genetic-circuit model and Boolean logic simulation — feature 25
//! (Cello-class v1).
//!
//! Cello compiles a Boolean specification into a genetic circuit: it
//! assigns characterised logic gates from a part library to the nodes
//! of a logic netlist and verifies the assignment by simulating the
//! circuit's truth table. This module provides the Cello workflow's
//! tractable v1:
//!
//! - a [`Gate`] library — each gate is a NOT / NOR-class repressor
//!   device with a Hill-style response curve (low / high output
//!   levels and a switching threshold);
//! - a [`Circuit`] — a directed acyclic netlist of [`LogicNode`]s
//!   wired from primary [`InputSignal`]s to a single output;
//! - [`Circuit::simulate_boolean`] — evaluates the **digital** truth
//!   table by thresholding; and
//! - [`Circuit::simulate_analog`] — propagates **continuous** signal
//!   levels through the gate response curves, the more realistic model
//!   that catches a circuit that is logically correct but has no noise
//!   margin.
//!
//! ## v1 caveats
//!
//! This *evaluates and scores* a gate assignment a caller supplies;
//! it does not perform Cello's combinatorial gate-assignment search
//! (simulated annealing over a large gate library) or its full
//! signal-noise "predicted circuit score". The supported gate types
//! are the NOT / NOR repressor family plus the derived AND / OR / NAND
//! (via De Morgan) — the genuine Cello v1 gate set; arbitrary
//! many-input gates are out of scope.

use std::collections::HashMap;

use crate::error::{Result, SysbioError};

/// A characterised genetic logic gate — a repressor device with a
/// sigmoidal transfer function.
///
/// The transfer function is a Hill curve: a high repressor input
/// drives the output low, a low input lets it go high. `ymin` /
/// `ymax` are the saturating output levels, `k` the input level at
/// half-maximal response and `n` the Hill steepness.
#[derive(Debug, Clone, PartialEq)]
pub struct Gate {
    /// Gate identifier (e.g. a repressor name like `"AmtR"`).
    pub id: String,
    /// Minimum (fully repressed) output level.
    pub ymin: f64,
    /// Maximum (de-repressed) output level.
    pub ymax: f64,
    /// Half-maximal input level `K`.
    pub k: f64,
    /// Hill coefficient `n`.
    pub n: f64,
}

impl Gate {
    /// A gate with the given transfer-function parameters.
    pub fn new(id: impl Into<String>, ymin: f64, ymax: f64, k: f64, n: f64) -> Self {
        Gate {
            id: id.into(),
            ymin,
            ymax,
            k,
            n,
        }
    }

    /// The repressor response: output level for a given total input
    /// repressor level. A NOT gate — high input ⇒ low output.
    pub fn response(&self, input: f64) -> f64 {
        let x = input.max(0.0);
        let ratio = (x / self.k).powf(self.n);
        // Inverted Hill: ymax at x=0, ymin as x -> infinity.
        self.ymin + (self.ymax - self.ymin) / (1.0 + ratio)
    }

    /// The digital threshold — the geometric mean of the output rails,
    /// the canonical "halfway in log space" decision level.
    pub fn threshold(&self) -> f64 {
        (self.ymin.max(1e-12) * self.ymax.max(1e-12)).sqrt()
    }
}

/// The Boolean operation a [`LogicNode`] computes over its inputs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LogicOp {
    /// Identity buffer (one input).
    Buffer,
    /// Logical NOT (one input).
    Not,
    /// Logical AND.
    And,
    /// Logical OR.
    Or,
    /// Logical NAND.
    Nand,
    /// Logical NOR.
    Nor,
    /// Exclusive OR.
    Xor,
}

impl LogicOp {
    /// Compute the operation over a slice of Boolean inputs.
    pub fn compute(&self, inputs: &[bool]) -> bool {
        match self {
            LogicOp::Buffer => inputs.first().copied().unwrap_or(false),
            LogicOp::Not => !inputs.first().copied().unwrap_or(false),
            LogicOp::And => inputs.iter().all(|&b| b),
            LogicOp::Or => inputs.iter().any(|&b| b),
            LogicOp::Nand => !inputs.iter().all(|&b| b),
            LogicOp::Nor => !inputs.iter().any(|&b| b),
            LogicOp::Xor => inputs.iter().filter(|&&b| b).count() % 2 == 1,
        }
    }
}

/// A reference to a signal feeding a logic node — either a primary
/// circuit input or the output of an earlier node.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Wire {
    /// A primary input, by index into [`Circuit::inputs`].
    Input(usize),
    /// The output of node `index` in [`Circuit::nodes`].
    Node(usize),
}

/// One logic node of a circuit: an operation, its input wires and the
/// genetic [`Gate`] assigned to realise it.
#[derive(Debug, Clone, PartialEq)]
pub struct LogicNode {
    /// Node identifier.
    pub id: String,
    /// The Boolean operation.
    pub op: LogicOp,
    /// Input wires (order matters for non-commutative ops).
    pub inputs: Vec<Wire>,
    /// Index of the assigned gate in the gate library, if any.
    pub gate: Option<usize>,
}

/// A named primary input signal with its low / high analog levels.
#[derive(Debug, Clone, PartialEq)]
pub struct InputSignal {
    /// Input identifier.
    pub id: String,
    /// Analog level representing logic `0`.
    pub low: f64,
    /// Analog level representing logic `1`.
    pub high: f64,
}

/// A genetic logic circuit — a DAG of [`LogicNode`]s.
#[derive(Debug, Clone, PartialEq)]
pub struct Circuit {
    /// Circuit identifier.
    pub id: String,
    /// Primary input signals.
    pub inputs: Vec<InputSignal>,
    /// Logic nodes in topological-friendly order (a node may only
    /// reference earlier nodes — checked by [`Circuit::validate`]).
    pub nodes: Vec<LogicNode>,
    /// The gate library available for assignment.
    pub gates: Vec<Gate>,
    /// Index of the node whose output is the circuit output.
    pub output_node: usize,
}

/// The result of one Boolean truth-table row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TruthRow {
    /// The Boolean input combination (parallel to [`Circuit::inputs`]).
    pub inputs: Vec<bool>,
    /// The resulting circuit output.
    pub output: bool,
}

/// The result of one analog evaluation.
#[derive(Debug, Clone, PartialEq)]
pub struct AnalogRow {
    /// The Boolean input combination that was applied.
    pub inputs: Vec<bool>,
    /// Continuous output level.
    pub output_level: f64,
    /// The digital value obtained by thresholding `output_level`.
    pub output_digital: bool,
}

impl Circuit {
    /// An empty circuit shell.
    pub fn new(id: impl Into<String>) -> Self {
        Circuit {
            id: id.into(),
            inputs: Vec::new(),
            nodes: Vec::new(),
            gates: Vec::new(),
            output_node: 0,
        }
    }

    /// Structural validation: acyclic wiring, in-range references, a
    /// valid output node.
    pub fn validate(&self) -> Result<()> {
        if self.nodes.is_empty() {
            return Err(SysbioError::invalid_model(
                "circuit",
                "circuit has no nodes",
            ));
        }
        if self.output_node >= self.nodes.len() {
            return Err(SysbioError::invalid_model(
                "circuit",
                "output node index out of range",
            ));
        }
        for (i, node) in self.nodes.iter().enumerate() {
            for w in &node.inputs {
                match w {
                    Wire::Input(idx) => {
                        if *idx >= self.inputs.len() {
                            return Err(SysbioError::invalid_model(
                                "circuit",
                                format!("node `{}` wires to missing input {idx}", node.id),
                            ));
                        }
                    }
                    Wire::Node(idx) => {
                        // DAG constraint: only reference earlier nodes.
                        if *idx >= i {
                            return Err(SysbioError::invalid_model(
                                "circuit",
                                format!(
                                    "node `{}` wires forward / to itself (node {idx})",
                                    node.id
                                ),
                            ));
                        }
                    }
                }
            }
            if let Some(g) = node.gate {
                if g >= self.gates.len() {
                    return Err(SysbioError::invalid_model(
                        "circuit",
                        format!("node `{}` assigned missing gate {g}", node.id),
                    ));
                }
            }
        }
        Ok(())
    }

    /// Evaluate the circuit's Boolean output for one input vector.
    fn run_boolean(&self, input_bits: &[bool]) -> bool {
        let mut node_out: Vec<bool> = Vec::with_capacity(self.nodes.len());
        for node in &self.nodes {
            let ins: Vec<bool> = node
                .inputs
                .iter()
                .map(|w| match w {
                    Wire::Input(i) => input_bits[*i],
                    Wire::Node(i) => node_out[*i],
                })
                .collect();
            node_out.push(node.op.compute(&ins));
        }
        node_out[self.output_node]
    }

    /// Simulate the full digital truth table (feature 25, digital).
    ///
    /// Enumerates all `2^#inputs` input combinations and records the
    /// circuit output of each.
    pub fn simulate_boolean(&self) -> Result<Vec<TruthRow>> {
        self.validate()?;
        let n = self.inputs.len();
        if n > 20 {
            return Err(SysbioError::invalid(
                "inputs",
                "truth-table enumeration capped at 20 inputs",
            ));
        }
        let mut rows = Vec::with_capacity(1 << n);
        for combo in 0..(1u32 << n) {
            let bits: Vec<bool> = (0..n).map(|b| (combo >> b) & 1 == 1).collect();
            rows.push(TruthRow {
                output: self.run_boolean(&bits),
                inputs: bits,
            });
        }
        Ok(rows)
    }

    /// Simulate the circuit with continuous signal levels (feature 25,
    /// analog).
    ///
    /// Each primary input contributes its `low` / `high` analog level
    /// per the Boolean combination; every node with an assigned gate
    /// propagates that level through the gate's transfer function. A
    /// multi-input node combines its inputs by summing (the
    /// repressor-loading model for NOR-class gates) before applying
    /// the gate. The continuous output is thresholded to a digital
    /// value so the caller can compare it against the Boolean result
    /// and detect a lost noise margin.
    ///
    /// Every node on the active path must have a gate assigned —
    /// otherwise this returns [`SysbioError::InvalidModel`].
    pub fn simulate_analog(&self) -> Result<Vec<AnalogRow>> {
        self.validate()?;
        let n = self.inputs.len();
        if n > 16 {
            return Err(SysbioError::invalid(
                "inputs",
                "analog enumeration capped at 16 inputs",
            ));
        }
        let mut rows = Vec::with_capacity(1 << n);
        for combo in 0..(1u32 << n) {
            let bits: Vec<bool> = (0..n).map(|b| (combo >> b) & 1 == 1).collect();
            let mut level: Vec<f64> = Vec::with_capacity(self.nodes.len());
            for node in &self.nodes {
                // Combine the input analog levels.
                let in_levels: Vec<f64> = node
                    .inputs
                    .iter()
                    .map(|w| match w {
                        Wire::Input(i) => {
                            if bits[*i] {
                                self.inputs[*i].high
                            } else {
                                self.inputs[*i].low
                            }
                        }
                        Wire::Node(i) => level[*i],
                    })
                    .collect();
                let combined: f64 = in_levels.iter().sum();
                let gate = node.gate.ok_or_else(|| {
                    SysbioError::invalid_model(
                        "circuit",
                        format!("node `{}` has no gate for analog simulation", node.id),
                    )
                })?;
                level.push(self.gates[gate].response(combined));
            }
            let out = level[self.output_node];
            let thr = self.nodes[self.output_node]
                .gate
                .map(|g| self.gates[g].threshold())
                .unwrap_or(0.5);
            rows.push(AnalogRow {
                inputs: bits,
                output_level: out,
                output_digital: out >= thr,
            });
        }
        Ok(rows)
    }

    /// Compare an analog simulation against the Boolean truth table
    /// and report the smallest output **separation** (ratio of the
    /// lowest logic-1 output to the highest logic-0 output) — a
    /// Cello-style on/off margin. A separation `> 1` means every
    /// logic-1 output beats every logic-0 output; the larger, the more
    /// robust the circuit.
    pub fn noise_margin(&self) -> Result<f64> {
        let digital = self.simulate_boolean()?;
        let analog = self.simulate_analog()?;
        let mut min_on = f64::INFINITY;
        let mut max_off = 0.0_f64;
        for (d, a) in digital.iter().zip(&analog) {
            if d.output {
                min_on = min_on.min(a.output_level);
            } else {
                max_off = max_off.max(a.output_level);
            }
        }
        if !min_on.is_finite() || max_off <= 0.0 {
            // Output is constant — no meaningful margin.
            return Ok(f64::INFINITY);
        }
        Ok(min_on / max_off)
    }

    /// Verify that the analog simulation reproduces the intended
    /// Boolean truth table after thresholding. Returns the indices of
    /// any mismatching rows (empty ⇒ the gate assignment is correct).
    pub fn verify_assignment(&self) -> Result<Vec<usize>> {
        let digital = self.simulate_boolean()?;
        let analog = self.simulate_analog()?;
        Ok(digital
            .iter()
            .zip(&analog)
            .enumerate()
            .filter(|(_, (d, a))| d.output != a.output_digital)
            .map(|(i, _)| i)
            .collect())
    }
}

/// Build a NOT-gate truth-table reference for a single input — a
/// helper for tests and quick checks.
pub fn not_truth_table() -> HashMap<bool, bool> {
    [(false, true), (true, false)].into_iter().collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A repressor gate with a clean 4-decade on/off ratio.
    fn strong_gate(id: &str) -> Gate {
        Gate::new(id, 0.01, 10.0, 1.0, 2.5)
    }

    /// A single-NOT circuit: one input, one NOT node.
    fn not_circuit() -> Circuit {
        let mut c = Circuit::new("not");
        c.inputs.push(InputSignal {
            id: "a".into(),
            low: 0.02,
            high: 8.0,
        });
        c.gates.push(strong_gate("R1"));
        c.nodes.push(LogicNode {
            id: "n0".into(),
            op: LogicOp::Not,
            inputs: vec![Wire::Input(0)],
            gate: Some(0),
        });
        c.output_node = 0;
        c
    }

    #[test]
    fn gate_response_is_inverting() {
        let g = strong_gate("R");
        // Low input -> near ymax; high input -> near ymin.
        assert!(g.response(0.0) > 9.0);
        assert!(g.response(100.0) < 0.1);
        // Monotone decreasing.
        assert!(g.response(0.5) > g.response(2.0));
    }

    #[test]
    fn boolean_not_truth_table() {
        let c = not_circuit();
        let rows = c.simulate_boolean().unwrap();
        assert_eq!(rows.len(), 2);
        let r0 = rows.iter().find(|r| r.inputs == vec![false]).unwrap();
        let r1 = rows.iter().find(|r| r.inputs == vec![true]).unwrap();
        assert!(r0.output); // NOT 0 = 1
        assert!(!r1.output); // NOT 1 = 0
    }

    #[test]
    fn analog_not_reproduces_digital() {
        let c = not_circuit();
        let mismatches = c.verify_assignment().unwrap();
        assert!(mismatches.is_empty(), "mismatched rows: {mismatches:?}");
    }

    #[test]
    fn strong_gate_has_good_noise_margin() {
        let c = not_circuit();
        let margin = c.noise_margin().unwrap();
        // logic-1 output well above logic-0 output.
        assert!(margin > 5.0, "margin {margin}");
    }

    #[test]
    fn weak_gate_loses_noise_margin() {
        // A gate with rails too close together and a shallow curve.
        let mut c = not_circuit();
        c.gates[0] = Gate::new("weak", 1.0, 1.5, 5.0, 1.0);
        let margin = c.noise_margin().unwrap();
        assert!(margin < 2.0, "expected poor margin, got {margin}");
    }

    #[test]
    fn nor_circuit_truth_table() {
        // Two inputs, one NOR node — the fundamental Cello gate.
        let mut c = Circuit::new("nor");
        c.inputs.push(InputSignal {
            id: "a".into(),
            low: 0.02,
            high: 8.0,
        });
        c.inputs.push(InputSignal {
            id: "b".into(),
            low: 0.02,
            high: 8.0,
        });
        c.gates.push(strong_gate("R1"));
        c.nodes.push(LogicNode {
            id: "n0".into(),
            op: LogicOp::Nor,
            inputs: vec![Wire::Input(0), Wire::Input(1)],
            gate: Some(0),
        });
        c.output_node = 0;
        let rows = c.simulate_boolean().unwrap();
        // NOR is high only when both inputs are low.
        for r in &rows {
            let expected = !(r.inputs[0] || r.inputs[1]);
            assert_eq!(r.output, expected, "{:?}", r.inputs);
        }
    }

    #[test]
    fn two_level_circuit_evaluates_in_topological_order() {
        // (a NOR b) feeding a NOT -> equals (a OR b).
        let mut c = Circuit::new("nor_not");
        for id in ["a", "b"] {
            c.inputs.push(InputSignal {
                id: id.into(),
                low: 0.02,
                high: 8.0,
            });
        }
        c.gates.push(strong_gate("R1"));
        c.gates.push(strong_gate("R2"));
        c.nodes.push(LogicNode {
            id: "nor".into(),
            op: LogicOp::Nor,
            inputs: vec![Wire::Input(0), Wire::Input(1)],
            gate: Some(0),
        });
        c.nodes.push(LogicNode {
            id: "inv".into(),
            op: LogicOp::Not,
            inputs: vec![Wire::Node(0)],
            gate: Some(1),
        });
        c.output_node = 1;
        let rows = c.simulate_boolean().unwrap();
        for r in &rows {
            assert_eq!(r.output, r.inputs[0] || r.inputs[1]);
        }
    }

    #[test]
    fn validate_rejects_forward_reference() {
        let mut c = not_circuit();
        c.nodes[0].inputs = vec![Wire::Node(0)]; // self-reference
        assert!(c.validate().is_err());
    }

    #[test]
    fn validate_rejects_missing_gate() {
        let mut c = not_circuit();
        c.nodes[0].gate = Some(99);
        assert!(c.validate().is_err());
    }

    #[test]
    fn logic_op_compute_matches_truth() {
        assert!(LogicOp::Xor.compute(&[true, false]));
        assert!(!LogicOp::Xor.compute(&[true, true]));
        assert!(LogicOp::Nand.compute(&[true, false]));
        assert!(!LogicOp::Nand.compute(&[true, true]));
        assert!(!not_truth_table()[&true]);
    }
}
