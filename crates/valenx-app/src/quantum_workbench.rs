//! The right-side **Quantum circuit** workbench — an in-house state-vector
//! quantum-circuit simulator over the [`valenx_quantum`] crate.
//!
//! The user picks a qubit count (1..=4), appends gates (Hadamard `H` and
//! Pauli-X `X` on a selectable qubit; `CNOT` on selectable control/target),
//! and presses **Run**. The recorded program is replayed on a fresh
//! `|0…0>` register via [`valenx_quantum::Circuit::run`] and the readout
//! lists the Born-rule `StateVector::probabilities` per computational
//! basis state.
//!
//! It mirrors the other real-time workbenches (`brep_workbench`): a
//! [`crate::workbench_chrome::workbench_shell`] panel gated on
//! [`crate::ValenxApp::show_quantum_workbench`], toggled from the View menu
//! and openable by the agent bridge under the workbench id `"quantum"`. The
//! bridge can set the controls (`agent_set` / `agent_control_names`), read
//! a status line (`agent_readout`), and fire the run via the RunCommand id
//! `quantum.run`.

use eframe::egui;

use valenx_quantum::Circuit;

use crate::ValenxApp;

/// Max qubits the workbench selector offers (the panel readout lists all
/// `2^n` basis states, so a small bound keeps the panel readable).
const MAX_PANEL_QUBITS: usize = 4;

// ---------------------------------------------------------------------------
// Workbench state
// ---------------------------------------------------------------------------

/// One recorded gate in the circuit. A `Copy` record (the [`Circuit`]
/// builder consumes `self`, so the panel keeps its own replayable list and
/// rebuilds the circuit on each Run).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Op {
    /// Hadamard on the given qubit.
    H(usize),
    /// Pauli-X (NOT) on the given qubit.
    X(usize),
    /// CNOT with `(control, target)`.
    Cnot(usize, usize),
}

impl Op {
    /// A short human/agent-readable label for the op list.
    fn label(self) -> String {
        match self {
            Op::H(q) => format!("H q{q}"),
            Op::X(q) => format!("X q{q}"),
            Op::Cnot(c, t) => format!("CNOT c{c}\u{2192}t{t}"),
        }
    }
}

/// The result of a Run: the per-basis-state probability distribution.
struct QuantumResult {
    /// Probability of each computational basis state `0..2^n`, in index order.
    probs: Vec<f64>,
    /// Qubit count the result was computed for (so the basis labels match).
    num_qubits: usize,
}

/// Persistent state for the Quantum circuit workbench: the qubit count, the
/// gate-edit selectors, the recorded op list, and the latest run result.
pub struct QuantumWorkbenchState {
    /// Number of qubits (1..=`MAX_PANEL_QUBITS`).
    pub num_qubits: usize,
    /// Currently selected qubit for the single-qubit gate buttons (H / X).
    pub sel_qubit: usize,
    /// Currently selected control qubit for CNOT.
    pub ctrl_qubit: usize,
    /// Currently selected target qubit for CNOT.
    pub tgt_qubit: usize,

    /// The recorded gate program, replayed on Run.
    ops: Vec<Op>,
    /// Latest run result, or `None` before the first run.
    result: Option<QuantumResult>,
    /// Last error string (shown in the panel), cleared on a good run.
    error: Option<String>,
}

impl Default for QuantumWorkbenchState {
    fn default() -> Self {
        // A known-good seed: the 2-qubit Bell state H(0) then CNOT(0->1),
        // whose run yields p(|00>) = p(|11>) = 0.5. Run is NOT executed
        // until the user (or the bridge) presses Run.
        Self {
            num_qubits: 2,
            sel_qubit: 0,
            ctrl_qubit: 0,
            tgt_qubit: 1,
            ops: vec![Op::H(0), Op::Cnot(0, 1)],
            result: None,
            error: None,
        }
    }
}

impl QuantumWorkbenchState {
    /// Clamp the gate-edit qubit selectors into `0..num_qubits` (called
    /// after the qubit count shrinks so stale indices can't escape range).
    fn clamp_selectors(&mut self) {
        let top = self.num_qubits.saturating_sub(1);
        self.sel_qubit = self.sel_qubit.min(top);
        self.ctrl_qubit = self.ctrl_qubit.min(top);
        self.tgt_qubit = self.tgt_qubit.min(top);
    }

    /// Drop any recorded op that references a qubit `>= num_qubits` (called
    /// when the qubit count shrinks so a later Run can't fail on a stale op).
    fn prune_ops(&mut self) {
        let n = self.num_qubits;
        self.ops.retain(|op| match *op {
            Op::H(q) | Op::X(q) => q < n,
            Op::Cnot(c, t) => c < n && t < n,
        });
    }

    /// Build a [`Circuit`] from the recorded ops, run it, and store the
    /// probability distribution (or an error). Factored out so the in-panel
    /// **Run** button and the `quantum.run` bridge id share one path.
    fn run_now(&mut self) {
        match self.try_run() {
            Ok(res) => {
                self.result = Some(res);
                self.error = None;
            }
            Err(e) => {
                self.error = Some(e);
                // Keep any previous result on screen; the error line says why.
            }
        }
    }

    /// The fallible run pipeline. Separated so `run_now` can map a typed
    /// error into the panel's error line.
    fn try_run(&self) -> Result<QuantumResult, String> {
        let mut circuit = Circuit::new(self.num_qubits);
        for op in &self.ops {
            circuit = match *op {
                Op::H(q) => circuit.h(q),
                Op::X(q) => circuit.x(q),
                Op::Cnot(c, t) => circuit.cnot(c, t),
            };
        }
        let sv = circuit.run().map_err(|e| e.to_string())?;
        Ok(QuantumResult {
            probs: sv.probabilities(),
            num_qubits: self.num_qubits,
        })
    }

    /// The user-visible captions of every control the agent bridge can set
    /// via `SetControl` (returned by `ListControls`). Order follows the form.
    pub fn agent_control_names() -> &'static [&'static str] {
        &["Qubits", "Gate qubit", "CNOT control", "CNOT target"]
    }

    /// Set one labelled control by its caption, for the agent `SetControl`
    /// bridge. Fail-loud on an unknown caption / wrong type / out-of-range;
    /// no state is written on error and nothing panics. All four controls
    /// read an integer (via an `f64` value); `Qubits` is clamped to
    /// `1..=MAX_PANEL_QUBITS` and the qubit selectors to `0..num_qubits`.
    pub fn agent_set(
        &mut self,
        name: &str,
        value: &crate::agent_commands::AgentValue,
    ) -> Result<(), String> {
        /// Read a non-negative integer from an `AgentValue::Float`.
        fn as_index(v: f64, what: &str) -> Result<usize, String> {
            if !v.is_finite() || v < 0.0 || v.fract() != 0.0 || v > 1.0e6 {
                return Err(format!("{what} must be a non-negative integer, got {v}"));
            }
            Ok(v as usize)
        }
        match name {
            "Qubits" => {
                let n = as_index(value.as_f64()?, "Qubits")?;
                if !(1..=MAX_PANEL_QUBITS).contains(&n) {
                    return Err(format!("Qubits must be in 1..={MAX_PANEL_QUBITS}, got {n}"));
                }
                self.num_qubits = n;
                self.clamp_selectors();
                self.prune_ops();
            }
            "Gate qubit" => {
                let q = as_index(value.as_f64()?, "Gate qubit")?;
                if q >= self.num_qubits {
                    return Err(format!(
                        "Gate qubit must be < Qubits ({}), got {q}",
                        self.num_qubits
                    ));
                }
                self.sel_qubit = q;
            }
            "CNOT control" => {
                let q = as_index(value.as_f64()?, "CNOT control")?;
                if q >= self.num_qubits {
                    return Err(format!(
                        "CNOT control must be < Qubits ({}), got {q}",
                        self.num_qubits
                    ));
                }
                self.ctrl_qubit = q;
            }
            "CNOT target" => {
                let q = as_index(value.as_f64()?, "CNOT target")?;
                if q >= self.num_qubits {
                    return Err(format!(
                        "CNOT target must be < Qubits ({}), got {q}",
                        self.num_qubits
                    ));
                }
                self.tgt_qubit = q;
            }
            other => return Err(format!("unknown quantum control: {other:?}")),
        }
        Ok(())
    }

    /// The current readout text for the agent `ReadReadout` bridge: the
    /// circuit (qubit count + op count) plus the per-basis-state
    /// probabilities once a run exists. `Some` once run, `None` before the
    /// first run (or after a run error with no prior result).
    pub fn agent_readout(&self) -> Option<String> {
        if let Some(err) = &self.error {
            if self.result.is_none() {
                return Some(format!("Quantum run failed: {err}"));
            }
        }
        let r = self.result.as_ref()?;
        let mut parts = Vec::new();
        for (b, p) in r.probs.iter().enumerate() {
            // Only list the basis states that carry meaningful probability,
            // so the line stays compact for larger registers.
            if *p > 1.0e-9 {
                parts.push(format!("|{:0width$b}>={:.4}", b, p, width = r.num_qubits));
            }
        }
        Some(format!(
            "Quantum \u{00B7} {} qubits \u{00B7} {} gates \u{00B7} P: {}",
            r.num_qubits,
            self.ops.len(),
            parts.join(" "),
        ))
    }
}

// ---------------------------------------------------------------------------
// Bridge run action (run)
// ---------------------------------------------------------------------------

/// Run the circuit (the in-panel **Run** action). Factored out so the
/// button and the `quantum.run` bridge id share one path.
pub(crate) fn run(app: &mut ValenxApp) {
    app.quantum.run_now();
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Draw the Quantum circuit workbench. A no-op unless toggled on via
/// View → Quantum circuit.
pub fn draw_quantum_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_quantum_workbench {
        return;
    }
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_quantum_workbench",
        "Quantum circuit \u{2014} state-vector simulator",
        quantum_workbench_body,
    );
    if close {
        app.show_quantum_workbench = false;
    }
}

// ---------------------------------------------------------------------------
// Workbench body
// ---------------------------------------------------------------------------

fn quantum_workbench_body(app: &mut ValenxApp, ui: &mut egui::Ui) {
    egui::ScrollArea::vertical()
        .auto_shrink([false, false])
        .show(ui, |ui| quantum_workbench_body_inner(app, ui));
}

fn quantum_workbench_body_inner(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.label(
        egui::RichText::new(
            "In-house state-vector quantum-circuit simulator \u{00B7} build a circuit (Hadamard / \
             Pauli-X on a selectable qubit; CNOT on selectable control/target), then Run replays it \
             on a fresh |0\u{2026}0> register and reports the probability of each computational \
             basis state.",
        )
        .weak()
        .small(),
    );
    ui.separator();

    let s = &mut app.quantum;

    // --- Qubit count --------------------------------------------------------
    egui::Grid::new("quantum_counts")
        .num_columns(2)
        .striped(true)
        .show(ui, |ui| {
            let lbl = ui.label("Qubits");
            if ui
                .add(
                    egui::DragValue::new(&mut s.num_qubits)
                        .speed(0.1)
                        .range(1..=MAX_PANEL_QUBITS),
                )
                .labelled_by(lbl.id)
                .on_hover_text("Number of qubits in the register (1..=4).")
                .changed()
            {
                s.clamp_selectors();
                s.prune_ops();
            }
            ui.end_row();

            let lbl = ui.label("Gate qubit");
            ui.add(
                egui::DragValue::new(&mut s.sel_qubit)
                    .speed(0.1)
                    .range(0..=s.num_qubits.saturating_sub(1)),
            )
            .labelled_by(lbl.id)
            .on_hover_text("Target qubit for the single-qubit gate buttons (H / X).");
            ui.end_row();

            let lbl = ui.label("CNOT control");
            ui.add(
                egui::DragValue::new(&mut s.ctrl_qubit)
                    .speed(0.1)
                    .range(0..=s.num_qubits.saturating_sub(1)),
            )
            .labelled_by(lbl.id)
            .on_hover_text("Control qubit for the CNOT gate.");
            ui.end_row();

            let lbl = ui.label("CNOT target");
            ui.add(
                egui::DragValue::new(&mut s.tgt_qubit)
                    .speed(0.1)
                    .range(0..=s.num_qubits.saturating_sub(1)),
            )
            .labelled_by(lbl.id)
            .on_hover_text("Target qubit for the CNOT gate (must differ from the control).");
            ui.end_row();
        });

    ui.add_space(6.0);

    // --- Gate buttons -------------------------------------------------------
    ui.label(egui::RichText::new("Append gate").strong());
    ui.horizontal(|ui| {
        if ui
            .button("H")
            .on_hover_text("Append a Hadamard on the selected gate qubit.")
            .clicked()
        {
            s.ops.push(Op::H(s.sel_qubit));
        }
        if ui
            .button("X")
            .on_hover_text("Append a Pauli-X (NOT) on the selected gate qubit.")
            .clicked()
        {
            s.ops.push(Op::X(s.sel_qubit));
        }
        if ui
            .button("CNOT")
            .on_hover_text("Append a CNOT with the selected control & target qubits.")
            .clicked()
        {
            if s.ctrl_qubit == s.tgt_qubit {
                s.error = Some("CNOT control and target must differ".to_string());
            } else {
                s.ops.push(Op::Cnot(s.ctrl_qubit, s.tgt_qubit));
            }
        }
        if ui
            .button("\u{232B} Clear")
            .on_hover_text("Remove all gates from the circuit.")
            .clicked()
        {
            s.ops.clear();
            s.result = None;
        }
    });

    // --- Circuit list -------------------------------------------------------
    ui.add_space(4.0);
    if s.ops.is_empty() {
        ui.label(egui::RichText::new("(empty circuit)").italics().weak());
    } else {
        let listing: Vec<String> = s.ops.iter().map(|op| op.label()).collect();
        ui.label(format!("Circuit: {}", listing.join("  \u{00B7}  ")));
    }

    // --- Run ----------------------------------------------------------------
    ui.add_space(6.0);
    if ui
        .button("\u{25B6} Run")
        .on_hover_text("Replay the circuit on a fresh |0\u{2026}0> register and report the basis-state probabilities.")
        .clicked()
    {
        s.run_now();
    }

    ui.add_space(6.0);
    ui.separator();
    draw_result(s, ui);
}

// ---------------------------------------------------------------------------
// Result render
// ---------------------------------------------------------------------------

fn draw_result(s: &QuantumWorkbenchState, ui: &mut egui::Ui) {
    if let Some(err) = &s.error {
        ui.label(
            egui::RichText::new(format!("Run error: {err}"))
                .color(egui::Color32::from_rgb(220, 110, 90))
                .strong(),
        );
        ui.add_space(4.0);
    }

    let Some(r) = s.result.as_ref() else {
        ui.label(
            egui::RichText::new("No result yet \u{2014} build a circuit, then press Run.")
                .italics()
                .weak(),
        );
        return;
    };

    ui.label(
        egui::RichText::new(format!(
            "Probabilities ({} basis states of {} qubits)",
            r.probs.len(),
            r.num_qubits
        ))
        .strong()
        .color(egui::Color32::from_rgb(150, 200, 230)),
    );
    egui::Grid::new("quantum_probs")
        .num_columns(2)
        .striped(true)
        .show(ui, |ui| {
            for (b, p) in r.probs.iter().enumerate() {
                ui.monospace(format!("|{:0width$b}>", b, width = r.num_qubits));
                ui.monospace(format!("{p:.5}"));
                ui.end_row();
            }
        });
}

// ---------------------------------------------------------------------------
// Tests (unit)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_seed_is_a_bell_pair() {
        let s = QuantumWorkbenchState::default();
        assert_eq!(s.num_qubits, 2);
        assert_eq!(s.ops.len(), 2);
        assert!(s.result.is_none());
    }

    #[test]
    fn run_produces_bell_distribution_and_readout() {
        let mut s = QuantumWorkbenchState::default();
        s.run_now();
        assert!(s.error.is_none(), "Bell seed should run: {:?}", s.error);
        let r = s.result.as_ref().expect("a result after run");
        assert_eq!(r.probs.len(), 4, "2 qubits => 4 basis states");
        // Bell state: p(|00>) = p(|11>) = 0.5, p(|01>) = p(|10>) = 0.
        assert!((r.probs[0] - 0.5).abs() < 1e-9, "p00: {}", r.probs[0]);
        assert!((r.probs[3] - 0.5).abs() < 1e-9, "p11: {}", r.probs[3]);
        assert!(r.probs[1].abs() < 1e-9 && r.probs[2].abs() < 1e-9);
        let readout = s.agent_readout().expect("readout after run");
        assert!(readout.contains("Quantum"), "readout: {readout}");
        assert!(
            readout.contains("|11>"),
            "readout names a basis state: {readout}"
        );
    }

    #[test]
    fn x_gate_flips_the_qubit() {
        let mut s = QuantumWorkbenchState {
            num_qubits: 1,
            sel_qubit: 0,
            ctrl_qubit: 0,
            tgt_qubit: 0,
            ops: vec![Op::X(0)],
            result: None,
            error: None,
        };
        s.run_now();
        let r = s.result.as_ref().expect("result");
        // X|0> = |1>: all probability on basis state 1.
        assert!(r.probs[1] > 0.999, "X should put p(|1>)=1: {:?}", r.probs);
    }

    #[test]
    fn agent_set_controls_and_clamps() {
        use crate::agent_commands::AgentValue;
        let mut s = QuantumWorkbenchState::default();
        s.agent_set("Qubits", &AgentValue::Float(3.0)).expect("q");
        assert_eq!(s.num_qubits, 3);
        s.agent_set("Gate qubit", &AgentValue::Float(2.0))
            .expect("gq");
        assert_eq!(s.sel_qubit, 2);
        // Shrinking the register clamps selectors and prunes stale ops.
        s.agent_set("Qubits", &AgentValue::Float(1.0)).expect("q1");
        assert_eq!(s.num_qubits, 1);
        assert_eq!(s.sel_qubit, 0, "selector clamped on shrink");
    }

    #[test]
    fn agent_set_rejects_bad_values() {
        use crate::agent_commands::AgentValue;
        let mut s = QuantumWorkbenchState::default();
        assert!(s.agent_set("Qubits", &AgentValue::Float(0.0)).is_err());
        assert!(s.agent_set("Qubits", &AgentValue::Float(9.0)).is_err());
        assert!(s.agent_set("Qubits", &AgentValue::Float(2.5)).is_err());
        assert!(s.agent_set("Gate qubit", &AgentValue::Float(5.0)).is_err());
        assert!(s.agent_set("bogus", &AgentValue::Float(1.0)).is_err());
    }

    #[test]
    fn prune_drops_out_of_range_ops() {
        let mut s = QuantumWorkbenchState {
            num_qubits: 4,
            ops: vec![Op::H(3), Op::Cnot(0, 1)],
            ..Default::default()
        };
        s.num_qubits = 1;
        s.prune_ops();
        // Only ops fully inside 0..1 survive.
        assert_eq!(s.ops.len(), 0, "all seeded ops referenced qubit >= 1");
    }

    #[test]
    fn readout_is_none_before_run() {
        let s = QuantumWorkbenchState::default();
        assert!(s.agent_readout().is_none(), "no readout before a run");
    }

    #[test]
    fn run_bridge_helper_runs_through_app() {
        let mut app = ValenxApp::default();
        run(&mut app);
        assert!(
            app.quantum.result.is_some(),
            "the quantum.run bridge helper should produce a result"
        );
    }

    #[test]
    fn control_names_are_listed() {
        let names = QuantumWorkbenchState::agent_control_names();
        for c in ["Qubits", "Gate qubit", "CNOT control", "CNOT target"] {
            assert!(names.contains(&c), "missing control name {c}");
        }
    }

    #[test]
    fn gate_enum_is_reachable() {
        // The panel uses the convenience h()/x()/cnot() Circuit builders;
        // confirm the underlying re-exported Gate type is still reachable.
        use valenx_quantum::Gate;
        assert_eq!(Gate::H.matrix().len(), 2);
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use egui::accesskit::{Node, NodeId, Role};

    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_quantum_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    fn has_named_node(nodes: &[(NodeId, Node)], name: &str) -> bool {
        nodes.iter().any(|(_, n)| n.name() == Some(name))
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_quantum_workbench);
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_quantum_workbench(&mut app, ctx);
        });
        // No panic = pass.
    }

    #[test]
    fn workbench_draws_when_shown() {
        let mut app = ValenxApp::default();
        app.show_quantum_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);
        assert!(!nodes.is_empty(), "a shown workbench produces a11y nodes");
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // Every numeric DragValue is a SpinButton and must be `labelled_by`
        // its caption so an AI / screen reader can find it by caption text.
        let mut app = ValenxApp::default();
        app.show_quantum_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let by_id: std::collections::HashMap<NodeId, &Node> =
            nodes.iter().map(|(id, n)| (*id, n)).collect();

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        // Qubits + Gate qubit + CNOT control + CNOT target = 4 spin buttons.
        assert!(
            spin_buttons.len() >= 4,
            "expected the 4 numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every DragValue must be labelled_by a caption (AI-drivable name)"
        );
        assert!(
            spin_buttons.iter().all(|n| {
                n.labelled_by()
                    .iter()
                    .any(|id| by_id.get(id).is_some_and(|t| t.name().is_some()))
            }),
            "every DragValue's labelled_by must point at a named caption node"
        );

        for caption in ["Qubits", "Gate qubit", "CNOT control", "CNOT target"] {
            assert!(
                has_named_node(&nodes, caption),
                "caption '{caption}' must be a named node in the a11y tree"
            );
        }
    }
}
