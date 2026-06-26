//! The right-side **Node Graph** workbench — a native, in-house **visual
//! node-graph editor**. The user drops nodes onto a canvas, drags them around,
//! and wires an **output port -> input port** to flow values between them; an
//! evaluation pass walks the graph in topological order and computes every
//! node's value. Wiring `Constant(2) + Constant(2) -> Add -> Output` makes the
//! Output node read `4`.
//!
//! This is the foundation for later **pipeline / systems** editing (bond
//! graphs, CAD -> FEA pipelines, signal flow): the node-type system is an
//! extensible [`NodeKind`] enum with a single `NodeKind::compute` dispatch and
//! typed in/out [`PortType`]s, so a new kind is one enum variant plus one
//! match arm — no canvas / wiring / topo-sort changes.
//!
//! ## Why in-house (crate-first decision)
//!
//! The canonical `egui_node_graph2` crate that fits valenx's egui **0.28** pin
//! (its `0.6.0` release) is a heavyweight, generics-first *framework* that owns
//! its own state, response and node-template trait model and renders its own
//! internal widget tree. valenx's hard requirements pull the other way: an
//! **extensible `NodeKind` enum + a compute fn** that *we* own, and — the #1
//! standing release gate — **AI-drivability** via directly-placed, accessibly-
//! named egui widgets the agent bridge can find by caption. A focused in-house
//! editor (a few hundred lines: node rects, bezier wires, port hit-testing,
//! click-drag wiring, a topological pass) fits those requirements far better
//! than bending the framework around them, so this follows the crate-first
//! policy's "adopt only if it fits" — here it genuinely doesn't.
//!
//! ## AI-drivable surface
//!
//! Mirrors the other workbenches: a [`crate::workbench_chrome::workbench_shell`]
//! panel gated on [`crate::ValenxApp::show_nodegraph_workbench`], toggled from
//! the View menu and openable by the agent bridge under the workbench id
//! `"nodegraph"`. The bridge can:
//! - set the labelled controls (`agent_set` / `agent_control_names`): the
//!   **selected node's constant value**, and a textual **add-node** command
//!   (e.g. `"Add node" = "add"`) so a graph can be built headlessly;
//! - read a status line (`agent_readout`): node count, edge count and the
//!   Output value(s);
//! - fire the evaluation via the `RunCommand` id `nodegraph.eval`, which runs
//!   the SAME topological pass the in-panel **Evaluate** button calls.

use std::collections::HashMap;

use eframe::egui;

use crate::ValenxApp;

// ---------------------------------------------------------------------------
// Graph model — extensible node kinds + typed ports
// ---------------------------------------------------------------------------

/// The data type carried by a port. Today every value is a scalar
/// [`PortType::Number`]; the enum exists so later kinds (bond-graph efforts /
/// flows, vectors, meshes, …) can add a variant and the wiring layer can refuse
/// to connect mismatched types without any other change.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PortType {
    /// A scalar `f64`.
    Number,
}

impl PortType {
    /// Whether a wire may run from an output of `self` into an input of `other`.
    /// Same-type only — the single rule the wiring layer enforces.
    fn connectable_to(self, other: PortType) -> bool {
        self == other
    }
}

/// The kind of a node — the extensible heart of the editor. A new node is one
/// variant here plus an arm in `NodeKind::inputs` / `NodeKind::compute`
/// (and, if it has a stored parameter, the constant-value plumbing). Bond-graph
/// elements and pipeline stages slot in the same way.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NodeKind {
    /// A source: emits its stored constant on its single output.
    Constant,
    /// Sums its two number inputs onto its output.
    Add,
    /// Multiplies its two number inputs onto its output.
    Multiply,
    /// A sink/display: shows whatever reaches its single input (no output).
    Output,
}

impl NodeKind {
    /// Every node kind, in palette order.
    pub const ALL: [NodeKind; 4] = [
        NodeKind::Constant,
        NodeKind::Add,
        NodeKind::Multiply,
        NodeKind::Output,
    ];

    /// Short, stable id used by the agent bridge / `add` command.
    pub fn id(self) -> &'static str {
        match self {
            NodeKind::Constant => "constant",
            NodeKind::Add => "add",
            NodeKind::Multiply => "multiply",
            NodeKind::Output => "output",
        }
    }

    /// Human caption shown on the node body and in menus.
    pub fn label(self) -> &'static str {
        match self {
            NodeKind::Constant => "Constant",
            NodeKind::Add => "Add",
            NodeKind::Multiply => "Multiply",
            NodeKind::Output => "Output",
        }
    }

    /// Parse a kind from its [`id`](Self::id) (case-insensitive). A few natural
    /// aliases are accepted so an agent can say `"+"` / `"sum"` etc.
    pub fn from_id(s: &str) -> Option<NodeKind> {
        match s.trim().to_ascii_lowercase().as_str() {
            "constant" | "const" | "value" | "number" => Some(NodeKind::Constant),
            "add" | "sum" | "+" => Some(NodeKind::Add),
            "multiply" | "mul" | "product" | "*" => Some(NodeKind::Multiply),
            "output" | "display" | "out" | "sink" => Some(NodeKind::Output),
            _ => None,
        }
    }

    /// The input ports of this kind (their captions + types), left-to-right.
    fn inputs(self) -> &'static [(&'static str, PortType)] {
        match self {
            NodeKind::Constant => &[],
            NodeKind::Add | NodeKind::Multiply => {
                &[("a", PortType::Number), ("b", PortType::Number)]
            }
            NodeKind::Output => &[("in", PortType::Number)],
        }
    }

    /// The output ports of this kind. `Output` is a pure sink (no outputs).
    fn outputs(self) -> &'static [(&'static str, PortType)] {
        match self {
            NodeKind::Output => &[],
            _ => &[("out", PortType::Number)],
        }
    }

    /// Whether this kind stores an editable constant (only `Constant` does).
    fn has_constant(self) -> bool {
        matches!(self, NodeKind::Constant)
    }

    /// Compute this node's output value from its already-evaluated input values
    /// (`None` = an input port with nothing wired in, treated as `0.0`). For
    /// `Constant` the inputs are empty and `param` is the stored constant. The
    /// single dispatch point new kinds extend.
    fn compute(self, inputs: &[Option<f64>], param: f64) -> f64 {
        let v = |i: usize| inputs.get(i).copied().flatten().unwrap_or(0.0);
        match self {
            NodeKind::Constant => param,
            NodeKind::Add => v(0) + v(1),
            NodeKind::Multiply => v(0) * v(1),
            // A sink simply passes its input through (so a downstream readout /
            // future chained pipeline can observe it).
            NodeKind::Output => v(0),
        }
    }
}

/// Stable handle for a node within one [`NodeGraph`]. A monotonically-rising
/// id (never reused) so wires stay valid as nodes are added.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct NodeId(pub u64);

/// One node instance: its kind, canvas position, stored constant and last
/// computed value.
#[derive(Clone, Debug)]
pub struct Node {
    /// This node's kind (drives ports + the compute step).
    pub kind: NodeKind,
    /// Top-left position on the canvas, in canvas-local points.
    pub pos: egui::Pos2,
    /// The editable constant (used only by [`NodeKind::Constant`]).
    pub value: f64,
    /// The value computed by the last [`NodeGraph::evaluate`] (for display).
    pub computed: Option<f64>,
}

/// A directed wire: the `from_node`'s output port `from_out` feeds the
/// `to_node`'s input port `to_in` (ports identified by index).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Edge {
    /// Source node.
    pub from_node: NodeId,
    /// Source output-port index.
    pub from_out: usize,
    /// Destination node.
    pub to_node: NodeId,
    /// Destination input-port index.
    pub to_in: usize,
}

/// The graph: a node table + an edge list + the id allocator. Pure data and
/// logic — no egui — so it is fully unit-testable headless.
#[derive(Clone, Debug, Default)]
pub struct NodeGraph {
    /// All nodes, keyed by id.
    pub nodes: HashMap<NodeId, Node>,
    /// All wires.
    pub edges: Vec<Edge>,
    /// Next id to hand out.
    next_id: u64,
}

impl NodeGraph {
    /// Add a node of `kind` at `pos`, returning its fresh id. A `Constant`
    /// starts at `0.0`.
    pub fn add_node(&mut self, kind: NodeKind, pos: egui::Pos2) -> NodeId {
        let id = NodeId(self.next_id);
        self.next_id += 1;
        self.nodes.insert(
            id,
            Node {
                kind,
                pos,
                value: 0.0,
                computed: None,
            },
        );
        id
    }

    /// Remove a node and every wire touching it.
    pub fn remove_node(&mut self, id: NodeId) {
        self.nodes.remove(&id);
        self.edges.retain(|e| e.from_node != id && e.to_node != id);
    }

    /// Try to add a wire `from_node.from_out -> to_node.to_in`. Refuses (returns
    /// `false`, changing nothing) when: either endpoint is missing, a port index
    /// is out of range, the port **types** differ, the endpoints are the same
    /// node, or the destination input is already occupied (each input takes one
    /// wire). A duplicate identical edge is also refused. Replacing the input's
    /// existing source is the caller's job (disconnect first).
    pub fn try_connect(
        &mut self,
        from_node: NodeId,
        from_out: usize,
        to_node: NodeId,
        to_in: usize,
    ) -> bool {
        if from_node == to_node {
            return false;
        }
        let (Some(src), Some(dst)) = (self.nodes.get(&from_node), self.nodes.get(&to_node)) else {
            return false;
        };
        let Some((_, out_ty)) = src.kind.outputs().get(from_out).copied() else {
            return false;
        };
        let Some((_, in_ty)) = dst.kind.inputs().get(to_in).copied() else {
            return false;
        };
        if !out_ty.connectable_to(in_ty) {
            return false;
        }
        // One wire per input port.
        if self
            .edges
            .iter()
            .any(|e| e.to_node == to_node && e.to_in == to_in)
        {
            return false;
        }
        self.edges.push(Edge {
            from_node,
            from_out,
            to_node,
            to_in,
        });
        true
    }

    /// Drop any wire feeding `(to_node, to_in)`.
    pub fn disconnect_input(&mut self, to_node: NodeId, to_in: usize) {
        self.edges
            .retain(|e| !(e.to_node == to_node && e.to_in == to_in));
    }

    /// Evaluate the whole graph in topological order, writing each node's
    /// [`Node::computed`]. On a cycle (no valid topo order) the nodes still on
    /// the frontier are left `None` (fail-soft: a cyclic graph shows blanks
    /// rather than hanging or panicking). Returns nothing; read results via
    /// [`Node::computed`] / [`Self::output_values`].
    pub fn evaluate(&mut self) {
        // Kahn's algorithm over node dependencies (an edge to_node depends on
        // from_node). in_deg counts prerequisite wires per node.
        let ids: Vec<NodeId> = {
            let mut v: Vec<NodeId> = self.nodes.keys().copied().collect();
            v.sort_unstable(); // deterministic order
            v
        };
        let mut in_deg: HashMap<NodeId, usize> = ids.iter().map(|&id| (id, 0usize)).collect();
        for e in &self.edges {
            if self.nodes.contains_key(&e.from_node) && self.nodes.contains_key(&e.to_node) {
                *in_deg.entry(e.to_node).or_insert(0) += 1;
            }
        }
        // Frontier: zero-in-degree nodes.
        let mut frontier: Vec<NodeId> = ids
            .iter()
            .copied()
            .filter(|id| in_deg.get(id).copied().unwrap_or(0) == 0)
            .collect();
        let mut order: Vec<NodeId> = Vec::with_capacity(ids.len());
        while let Some(id) = frontier.pop() {
            order.push(id);
            for e in &self.edges {
                if e.from_node == id {
                    if let Some(d) = in_deg.get_mut(&e.to_node) {
                        *d -= 1;
                        if *d == 0 {
                            frontier.push(e.to_node);
                        }
                    }
                }
            }
        }

        // Reset everything, then fill in topo order. Nodes left out of `order`
        // (in a cycle) stay `None`.
        let mut computed: HashMap<NodeId, f64> = HashMap::new();
        for &id in &ids {
            if let Some(n) = self.nodes.get_mut(&id) {
                n.computed = None;
            }
        }
        for &id in &order {
            let (kind, param) = {
                let n = &self.nodes[&id];
                (n.kind, n.value)
            };
            // Gather this node's input values by port index from incoming wires.
            let n_inputs = kind.inputs().len();
            let mut inputs = vec![None; n_inputs];
            for e in &self.edges {
                if e.to_node == id && e.to_in < n_inputs {
                    inputs[e.to_in] = computed.get(&e.from_node).copied();
                }
            }
            let out = kind.compute(&inputs, param);
            computed.insert(id, out);
            if let Some(n) = self.nodes.get_mut(&id) {
                n.computed = Some(out);
            }
        }
    }

    /// The computed values of every [`NodeKind::Output`] node, in node-id order
    /// — the "results" of the graph. `None` for an Output left in a cycle.
    pub fn output_values(&self) -> Vec<(NodeId, Option<f64>)> {
        let mut outs: Vec<(NodeId, Option<f64>)> = self
            .nodes
            .iter()
            .filter(|(_, n)| n.kind == NodeKind::Output)
            .map(|(id, n)| (*id, n.computed))
            .collect();
        outs.sort_by_key(|(id, _)| *id);
        outs
    }

    /// Number of nodes / edges (for the readout).
    fn counts(&self) -> (usize, usize) {
        (self.nodes.len(), self.edges.len())
    }
}

// ---------------------------------------------------------------------------
// Workbench state
// ---------------------------------------------------------------------------

/// Where, on each node, a port circle sits (resolved at draw time). Cached per
/// frame so click-drag wiring can hit-test ports.
#[derive(Clone, Copy)]
struct PortHit {
    node: NodeId,
    index: usize,
    center: egui::Pos2,
}

/// An in-progress wire drag from an output port (set on press, cleared on
/// release).
#[derive(Clone, Copy)]
struct PendingWire {
    from_node: NodeId,
    from_out: usize,
}

/// Persistent state for the Node Graph workbench: the graph itself plus the
/// transient canvas interaction (selection, the node being dragged, the wire
/// being dragged).
pub struct NodeGraphWorkbenchState {
    /// The edited graph.
    pub graph: NodeGraph,
    /// The currently-selected node (its constant is the agent-settable value).
    pub selected: Option<NodeId>,
    /// The kind chosen in the "add node" combo.
    pub add_kind: NodeKind,
    /// Node currently being dragged (and the grab offset), transient per drag.
    dragging: Option<(NodeId, egui::Vec2)>,
    /// A wire being dragged out of an output port, transient per drag.
    pending: Option<PendingWire>,
}

impl Default for NodeGraphWorkbenchState {
    fn default() -> Self {
        // Seed a tiny, self-explanatory graph: Constant(2) + Constant(2) -> Add
        // -> Output, pre-wired so the very first Evaluate shows 4 (the canonical
        // demo). Users can clear it and build their own.
        let mut graph = NodeGraph::default();
        let c1 = graph.add_node(NodeKind::Constant, egui::pos2(20.0, 40.0));
        let c2 = graph.add_node(NodeKind::Constant, egui::pos2(20.0, 140.0));
        let add = graph.add_node(NodeKind::Add, egui::pos2(200.0, 80.0));
        let out = graph.add_node(NodeKind::Output, egui::pos2(360.0, 90.0));
        if let Some(n) = graph.nodes.get_mut(&c1) {
            n.value = 2.0;
        }
        if let Some(n) = graph.nodes.get_mut(&c2) {
            n.value = 2.0;
        }
        graph.try_connect(c1, 0, add, 0);
        graph.try_connect(c2, 0, add, 1);
        graph.try_connect(add, 0, out, 0);
        graph.evaluate();
        Self {
            graph,
            selected: Some(c1),
            add_kind: NodeKind::Constant,
            dragging: None,
            pending: None,
        }
    }
}

impl NodeGraphWorkbenchState {
    /// Add a node of `kind` at a default spawn spot, select it, and return its
    /// id. Shared by the in-panel **Add** button and the agent `add` command.
    fn add_node(&mut self, kind: NodeKind) -> NodeId {
        // Stagger spawns so successive adds don't stack exactly.
        let n = self.graph.nodes.len() as f32;
        let pos = egui::pos2(40.0 + (n % 5.0) * 24.0, 40.0 + (n % 5.0) * 24.0);
        let id = self.graph.add_node(kind, pos);
        self.selected = Some(id);
        id
    }

    /// Run the topological evaluation now. Factored out so the in-panel
    /// **Evaluate** button and the `nodegraph.eval` bridge id share one path.
    fn run_eval(&mut self) {
        self.graph.evaluate();
    }

    /// Captions of every control the agent bridge can `SetControl`.
    pub fn agent_control_names() -> &'static [&'static str] {
        &["Constant value", "Add node"]
    }

    /// Set one labelled control by caption for the agent `SetControl` bridge.
    /// Fail-loud on an unknown caption / wrong type / unusable target; no state
    /// is written on error and nothing panics.
    ///
    /// - **`Constant value`** (number): sets the *selected* node's stored
    ///   constant. Errors if there is no selection or the selected node is not a
    ///   `Constant`.
    /// - **`Add node`** (string): a node kind id (`constant` / `add` /
    ///   `multiply` / `output`, plus aliases) — adds that node and selects it.
    pub fn agent_set(
        &mut self,
        name: &str,
        value: &crate::agent_commands::AgentValue,
    ) -> Result<(), String> {
        match name {
            "Constant value" => {
                let v = value.as_f64()?;
                if !v.is_finite() {
                    return Err(format!("Constant value must be finite, got {v}"));
                }
                let Some(id) = self.selected else {
                    return Err("Constant value: no node is selected".to_string());
                };
                let Some(node) = self.graph.nodes.get_mut(&id) else {
                    return Err("Constant value: the selected node no longer exists".to_string());
                };
                if !node.kind.has_constant() {
                    return Err(format!(
                        "Constant value: selected node is a {}, not a Constant",
                        node.kind.label()
                    ));
                }
                node.value = v;
                Ok(())
            }
            "Add node" => {
                let s = value.as_str()?;
                match NodeKind::from_id(s) {
                    Some(kind) => {
                        self.add_node(kind);
                        Ok(())
                    }
                    None => Err(format!(
                        "Add node: unknown kind {s:?} (use constant/add/multiply/output)"
                    )),
                }
            }
            other => Err(format!("unknown nodegraph control: {other:?}")),
        }
    }

    /// Readout for the agent `ReadReadout` bridge: node count, edge count and
    /// the Output value(s). Always `Some` (the graph always exists); reports
    /// "no Output node" / "not evaluated" honestly.
    pub fn agent_readout(&self) -> Option<String> {
        let (n, e) = self.graph.counts();
        let outs = self.graph.output_values();
        let outs_txt = if outs.is_empty() {
            "no Output node".to_string()
        } else {
            outs.iter()
                .map(|(id, v)| match v {
                    Some(x) => format!("#{}={:.4}", id.0, x),
                    None => format!("#{}=(not evaluated)", id.0),
                })
                .collect::<Vec<_>>()
                .join(", ")
        };
        Some(format!(
            "Node graph \u{00B7} {n} nodes \u{00B7} {e} edges \u{00B7} outputs: {outs_txt}"
        ))
    }
}

// ---------------------------------------------------------------------------
// Bridge run action (evaluate)
// ---------------------------------------------------------------------------

/// Run the topological evaluation (the in-panel **Evaluate** action). Factored
/// out so the button and the `nodegraph.eval` bridge id share one path.
pub(crate) fn run(app: &mut ValenxApp) {
    app.nodegraph.run_eval();
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Draw the Node Graph workbench. A no-op unless toggled on via View -> Node
/// Graph.
pub fn draw_nodegraph_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_nodegraph_workbench {
        return;
    }
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_nodegraph_workbench",
        "Node Graph (visual node editor)",
        nodegraph_workbench_body,
    );
    if close {
        app.show_nodegraph_workbench = false;
    }
}

// ---------------------------------------------------------------------------
// Workbench body
// ---------------------------------------------------------------------------

fn nodegraph_workbench_body(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.label(
        egui::RichText::new(
            "In-house visual node-graph editor [add nodes, drag them, and wire an output port \
             -> an input port; press Evaluate to flow values through the graph in topological \
             order. The seeded demo Constant(2) + Constant(2) -> Add -> Output reads 4. The \
             node-type system is extensible — the foundation for bond graphs and CAD->FEA \
             pipelines].",
        )
        .weak()
        .small(),
    );
    ui.separator();

    // --- Toolbar: add a node / evaluate / clear -----------------------------
    ui.horizontal(|ui| {
        let lbl = ui.label("Add node");
        egui::ComboBox::from_id_source("nodegraph_add_kind")
            .selected_text(app.nodegraph.add_kind.label())
            .show_ui(ui, |ui| {
                for k in NodeKind::ALL {
                    ui.selectable_value(&mut app.nodegraph.add_kind, k, k.label());
                }
            })
            .response
            .labelled_by(lbl.id);
        if ui
            .button("\u{2795} Add")
            .on_hover_text("Add a node of the chosen kind to the canvas.")
            .clicked()
        {
            let k = app.nodegraph.add_kind;
            app.nodegraph.add_node(k);
        }
        ui.separator();
        if ui
            .button("\u{25B6} Evaluate")
            .on_hover_text("Walk the graph in topological order and compute every node's value.")
            .clicked()
        {
            app.nodegraph.run_eval();
        }
        if ui
            .button("\u{1F5D1} Clear")
            .on_hover_text("Delete all nodes and wires.")
            .clicked()
        {
            app.nodegraph.graph = NodeGraph::default();
            app.nodegraph.selected = None;
        }
    });

    // --- Selected-node inspector (the agent-settable Constant value) --------
    ui.horizontal(|ui| match app.nodegraph.selected {
        Some(id) if app.nodegraph.graph.nodes.contains_key(&id) => {
            let node = app.nodegraph.graph.nodes.get_mut(&id).unwrap();
            ui.label(format!("Selected: {} #{}", node.kind.label(), id.0));
            if node.kind.has_constant() {
                ui.separator();
                let lbl = ui.label("Constant value");
                ui.add(
                    egui::DragValue::new(&mut node.value)
                        .speed(0.1)
                        .max_decimals(4),
                )
                .labelled_by(lbl.id)
                .on_hover_text("This Constant node's emitted value.");
            }
            if ui
                .button("\u{2716} Delete node")
                .on_hover_text("Remove the selected node and its wires.")
                .clicked()
            {
                app.nodegraph.graph.remove_node(id);
                app.nodegraph.selected = None;
            }
        }
        _ => {
            ui.label(egui::RichText::new("No node selected — click a node to select it.").weak());
        }
    });

    ui.separator();

    // --- The canvas ---------------------------------------------------------
    draw_canvas(app, ui);

    // --- Output readout below the canvas ------------------------------------
    ui.add_space(4.0);
    let outs = app.nodegraph.graph.output_values();
    if outs.is_empty() {
        ui.label(egui::RichText::new("Outputs: (add an Output node)").weak());
    } else {
        for (id, v) in outs {
            let txt = match v {
                Some(x) => format!("Output #{}  =  {x:.4}", id.0),
                None => format!("Output #{}  =  (press Evaluate)", id.0),
            };
            ui.label(egui::RichText::new(txt).strong());
        }
    }
}

// ---------------------------------------------------------------------------
// Canvas rendering + interaction
// ---------------------------------------------------------------------------

/// Fixed node-body width / port spacing (canvas-local points).
const NODE_W: f32 = 120.0;
const HEADER_H: f32 = 22.0;
const PORT_DY: f32 = 22.0;
const PORT_R: f32 = 5.0;

/// Total height of a node body given its larger port column.
fn node_height(kind: NodeKind) -> f32 {
    let rows = kind.inputs().len().max(kind.outputs().len()).max(1);
    HEADER_H + rows as f32 * PORT_DY + 8.0
}

/// Centre of input port `i` for a node whose top-left is `origin`.
fn input_pos(origin: egui::Pos2, i: usize) -> egui::Pos2 {
    egui::pos2(origin.x, origin.y + HEADER_H + PORT_DY * (i as f32 + 0.5))
}

/// Centre of output port `i` for a node whose top-left is `origin`.
fn output_pos(origin: egui::Pos2, i: usize) -> egui::Pos2 {
    egui::pos2(
        origin.x + NODE_W,
        origin.y + HEADER_H + PORT_DY * (i as f32 + 0.5),
    )
}

fn draw_canvas(app: &mut ValenxApp, ui: &mut egui::Ui) {
    let desired = egui::vec2(
        ui.available_width(),
        320.0_f32.max(ui.available_height() - 40.0),
    );
    let (resp, painter) = ui.allocate_painter(desired, egui::Sense::click_and_drag());
    let origin = resp.rect.min.to_vec2(); // canvas-local → screen offset

    // Background.
    painter.rect_filled(resp.rect, 4.0, egui::Color32::from_rgb(24, 26, 32));
    painter.rect_stroke(
        resp.rect,
        4.0,
        egui::Stroke::new(1.0, egui::Color32::from_rgb(50, 55, 68)),
    );

    let s = &mut app.nodegraph;

    // Collect port screen positions for wires + hit-testing.
    let mut input_hits: Vec<PortHit> = Vec::new();
    let mut output_hits: Vec<PortHit> = Vec::new();
    let ids: Vec<NodeId> = {
        let mut v: Vec<NodeId> = s.graph.nodes.keys().copied().collect();
        v.sort_unstable();
        v
    };
    for &id in &ids {
        let node = &s.graph.nodes[&id];
        let scr = node.pos + origin; // top-left on screen (Pos2 + Vec2 = Pos2)
        for i in 0..node.kind.inputs().len() {
            input_hits.push(PortHit {
                node: id,
                index: i,
                center: input_pos(scr, i),
            });
        }
        for i in 0..node.kind.outputs().len() {
            output_hits.push(PortHit {
                node: id,
                index: i,
                center: output_pos(scr, i),
            });
        }
    }
    let port_at = |hits: &[PortHit], p: egui::Pos2| -> Option<PortHit> {
        hits.iter()
            .find(|h| h.center.distance(p) <= PORT_R * 2.2)
            .copied()
    };

    // --- Draw existing wires (cubic beziers, left→right) --------------------
    for e in &s.graph.edges {
        let (Some(src), Some(dst)) = (
            s.graph.nodes.get(&e.from_node),
            s.graph.nodes.get(&e.to_node),
        ) else {
            continue;
        };
        let a = output_pos(src.pos + origin, e.from_out);
        let b = input_pos(dst.pos + origin, e.to_in);
        draw_wire(&painter, a, b, egui::Color32::from_rgb(120, 170, 255));
    }

    // --- Pointer interaction ------------------------------------------------
    let pointer = resp.interact_pointer_pos();

    // Begin interactions on press.
    if resp.drag_started() {
        if let Some(p) = pointer {
            // Output port press → start a wire.
            if let Some(h) = port_at(&output_hits, p) {
                s.pending = Some(PendingWire {
                    from_node: h.node,
                    from_out: h.index,
                });
            } else if let Some(h) = port_at(&input_hits, p) {
                // Input port press → detach any existing wire (lets the user
                // re-route).
                s.graph.disconnect_input(h.node, h.index);
            } else if let Some(id) = node_at(s, p - origin) {
                // Node body press → select + begin moving it. `grab` is the
                // pointer-to-node-origin offset (Pos2 - Pos2 = Vec2).
                s.selected = Some(id);
                let grab = (p - origin) - s.graph.nodes[&id].pos;
                s.dragging = Some((id, grab));
            }
        }
    }

    // Continue a node drag. `p - origin` is a Pos2; subtracting the Vec2 grab
    // offset yields the node's new canvas-local Pos2.
    if let (Some((id, grab)), Some(p)) = (s.dragging, pointer) {
        if let Some(node) = s.graph.nodes.get_mut(&id) {
            node.pos = (p - origin) - grab;
        }
    }

    // Draw an in-progress wire to the cursor.
    if let (Some(pw), Some(p)) = (s.pending, pointer) {
        if let Some(src) = s.graph.nodes.get(&pw.from_node) {
            let a = output_pos(src.pos + origin, pw.from_out);
            draw_wire(&painter, a, p, egui::Color32::from_rgb(200, 200, 120));
        }
    }

    // Finish interactions on release.
    if resp.drag_stopped() {
        if let (Some(pw), Some(p)) = (s.pending, pointer) {
            // Drop on an input port → connect.
            if let Some(h) = port_at(&input_hits, p) {
                s.graph
                    .try_connect(pw.from_node, pw.from_out, h.node, h.index);
            }
        }
        s.pending = None;
        s.dragging = None;
    }

    // Plain click (no drag) on a node selects it.
    if resp.clicked() {
        if let Some(p) = pointer {
            if let Some(id) = node_at(s, p - origin) {
                s.selected = Some(id);
            }
        }
    }

    // --- Draw nodes (on top of wires) ---------------------------------------
    for &id in &ids {
        let node = &s.graph.nodes[&id];
        let scr = node.pos + origin; // Pos2 + Vec2 = Pos2
        let h = node_height(node.kind);
        let body = egui::Rect::from_min_size(scr, egui::vec2(NODE_W, h));
        let selected = s.selected == Some(id);

        // Body.
        painter.rect_filled(body, 5.0, egui::Color32::from_rgb(40, 44, 54));
        painter.rect_stroke(
            body,
            5.0,
            egui::Stroke::new(
                if selected { 2.0 } else { 1.0 },
                if selected {
                    egui::Color32::from_rgb(255, 210, 90)
                } else {
                    egui::Color32::from_rgb(70, 78, 95)
                },
            ),
        );
        // Header.
        let header = egui::Rect::from_min_size(scr, egui::vec2(NODE_W, HEADER_H));
        painter.rect_filled(header, 5.0, header_color(node.kind));
        painter.text(
            header.center(),
            egui::Align2::CENTER_CENTER,
            node.kind.label(),
            egui::FontId::proportional(12.0),
            egui::Color32::from_gray(235),
        );

        // Constant value / computed value line.
        let val_txt = if node.kind.has_constant() {
            format!("= {:.3}", node.value)
        } else if let Some(c) = node.computed {
            format!("-> {c:.3}")
        } else {
            String::new()
        };
        if !val_txt.is_empty() {
            painter.text(
                egui::pos2(scr.x + NODE_W * 0.5, scr.y + h - 8.0),
                egui::Align2::CENTER_BOTTOM,
                val_txt,
                egui::FontId::monospace(10.0),
                egui::Color32::from_gray(190),
            );
        }

        // Input ports (left) + captions.
        for (i, (cap, _)) in node.kind.inputs().iter().enumerate() {
            let c = input_pos(scr, i);
            painter.circle_filled(c, PORT_R, egui::Color32::from_rgb(150, 200, 150));
            painter.text(
                egui::pos2(c.x + 8.0, c.y),
                egui::Align2::LEFT_CENTER,
                *cap,
                egui::FontId::proportional(10.0),
                egui::Color32::from_gray(170),
            );
        }
        // Output ports (right) + captions.
        for (i, (cap, _)) in node.kind.outputs().iter().enumerate() {
            let c = output_pos(scr, i);
            painter.circle_filled(c, PORT_R, egui::Color32::from_rgb(200, 180, 120));
            painter.text(
                egui::pos2(c.x - 8.0, c.y),
                egui::Align2::RIGHT_CENTER,
                *cap,
                egui::FontId::proportional(10.0),
                egui::Color32::from_gray(170),
            );
        }
    }
}

/// Header tint per node kind (visual grouping).
fn header_color(kind: NodeKind) -> egui::Color32 {
    match kind {
        NodeKind::Constant => egui::Color32::from_rgb(60, 90, 120),
        NodeKind::Add => egui::Color32::from_rgb(80, 100, 60),
        NodeKind::Multiply => egui::Color32::from_rgb(100, 80, 60),
        NodeKind::Output => egui::Color32::from_rgb(100, 60, 90),
    }
}

/// Draw a left->right cubic bezier wire with horizontal tangents.
fn draw_wire(painter: &egui::Painter, a: egui::Pos2, b: egui::Pos2, color: egui::Color32) {
    let dx = ((b.x - a.x).abs() * 0.5).max(30.0);
    let bezier = egui::epaint::CubicBezierShape::from_points_stroke(
        [a, egui::pos2(a.x + dx, a.y), egui::pos2(b.x - dx, b.y), b],
        false,
        egui::Color32::TRANSPARENT,
        egui::Stroke::new(2.0, color),
    );
    painter.add(bezier);
}

/// Topmost node whose body contains canvas-local point `p` (later-drawn /
/// higher-id nodes win, matching draw order).
fn node_at(s: &NodeGraphWorkbenchState, p: egui::Pos2) -> Option<NodeId> {
    let mut ids: Vec<NodeId> = s.graph.nodes.keys().copied().collect();
    ids.sort_unstable();
    let mut found = None;
    for id in ids {
        let node = &s.graph.nodes[&id];
        let r = egui::Rect::from_min_size(node.pos, egui::vec2(NODE_W, node_height(node.kind)));
        if r.contains(p) {
            found = Some(id); // keep the last (topmost) hit
        }
    }
    found
}

// ---------------------------------------------------------------------------
// Tests (unit)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::agent_commands::AgentValue;

    #[test]
    fn constant_add_output_evaluates_to_four() {
        // The canonical demo: Constant(2) + Constant(2) -> Add -> Output == 4.
        let mut g = NodeGraph::default();
        let c1 = g.add_node(NodeKind::Constant, egui::pos2(0.0, 0.0));
        let c2 = g.add_node(NodeKind::Constant, egui::pos2(0.0, 50.0));
        let add = g.add_node(NodeKind::Add, egui::pos2(100.0, 25.0));
        let out = g.add_node(NodeKind::Output, egui::pos2(200.0, 25.0));
        g.nodes.get_mut(&c1).unwrap().value = 2.0;
        g.nodes.get_mut(&c2).unwrap().value = 2.0;
        assert!(g.try_connect(c1, 0, add, 0));
        assert!(g.try_connect(c2, 0, add, 1));
        assert!(g.try_connect(add, 0, out, 0));
        g.evaluate();
        assert_eq!(g.nodes[&out].computed, Some(4.0));
        assert_eq!(g.output_values(), vec![(out, Some(4.0))]);
    }

    #[test]
    fn multiply_chains_through_topo_order() {
        // (3 + 4) * 5 == 35, regardless of node insertion order.
        let mut g = NodeGraph::default();
        let mul = g.add_node(NodeKind::Multiply, egui::pos2(0.0, 0.0));
        let out = g.add_node(NodeKind::Output, egui::pos2(0.0, 0.0));
        let a = g.add_node(NodeKind::Constant, egui::pos2(0.0, 0.0));
        let b = g.add_node(NodeKind::Constant, egui::pos2(0.0, 0.0));
        let add = g.add_node(NodeKind::Add, egui::pos2(0.0, 0.0));
        let five = g.add_node(NodeKind::Constant, egui::pos2(0.0, 0.0));
        g.nodes.get_mut(&a).unwrap().value = 3.0;
        g.nodes.get_mut(&b).unwrap().value = 4.0;
        g.nodes.get_mut(&five).unwrap().value = 5.0;
        assert!(g.try_connect(a, 0, add, 0));
        assert!(g.try_connect(b, 0, add, 1));
        assert!(g.try_connect(add, 0, mul, 0));
        assert!(g.try_connect(five, 0, mul, 1));
        assert!(g.try_connect(mul, 0, out, 0));
        g.evaluate();
        assert_eq!(g.nodes[&out].computed, Some(35.0));
    }

    #[test]
    fn unwired_input_is_zero() {
        // Add with only one wire -> the other input defaults to 0.
        let mut g = NodeGraph::default();
        let c = g.add_node(NodeKind::Constant, egui::pos2(0.0, 0.0));
        let add = g.add_node(NodeKind::Add, egui::pos2(0.0, 0.0));
        g.nodes.get_mut(&c).unwrap().value = 7.0;
        assert!(g.try_connect(c, 0, add, 0));
        g.evaluate();
        assert_eq!(g.nodes[&add].computed, Some(7.0));
    }

    #[test]
    fn type_and_arity_rules_reject_bad_wires() {
        let mut g = NodeGraph::default();
        let c = g.add_node(NodeKind::Constant, egui::pos2(0.0, 0.0));
        let out = g.add_node(NodeKind::Output, egui::pos2(0.0, 0.0));
        // Output has no output port -> cannot be a wire source.
        assert!(!g.try_connect(out, 0, c, 0));
        // Same-node + Constant-has-no-input both refused.
        assert!(!g.try_connect(c, 0, c, 0));
        // Out-of-range port index refused.
        let add = g.add_node(NodeKind::Add, egui::pos2(0.0, 0.0));
        assert!(!g.try_connect(c, 0, add, 9));
        // Valid wire accepted, second wire into the same input refused.
        assert!(g.try_connect(c, 0, add, 0));
        assert!(!g.try_connect(c, 0, add, 0));
    }

    #[test]
    fn cycle_is_fail_soft_not_a_hang() {
        // Two Adds feeding each other: no topo order -> both stay None, no panic.
        let mut g = NodeGraph::default();
        let a = g.add_node(NodeKind::Add, egui::pos2(0.0, 0.0));
        let b = g.add_node(NodeKind::Add, egui::pos2(0.0, 0.0));
        assert!(g.try_connect(a, 0, b, 0));
        assert!(g.try_connect(b, 0, a, 0));
        g.evaluate();
        assert_eq!(g.nodes[&a].computed, None);
        assert_eq!(g.nodes[&b].computed, None);
    }

    #[test]
    fn disconnect_input_clears_the_wire() {
        let mut g = NodeGraph::default();
        let c = g.add_node(NodeKind::Constant, egui::pos2(0.0, 0.0));
        let out = g.add_node(NodeKind::Output, egui::pos2(0.0, 0.0));
        assert!(g.try_connect(c, 0, out, 0));
        assert_eq!(g.edges.len(), 1);
        g.disconnect_input(out, 0);
        assert!(g.edges.is_empty());
    }

    #[test]
    fn remove_node_drops_incident_edges() {
        let mut g = NodeGraph::default();
        let c = g.add_node(NodeKind::Constant, egui::pos2(0.0, 0.0));
        let out = g.add_node(NodeKind::Output, egui::pos2(0.0, 0.0));
        assert!(g.try_connect(c, 0, out, 0));
        g.remove_node(c);
        assert!(g.edges.is_empty());
        assert!(!g.nodes.contains_key(&c));
    }

    #[test]
    fn default_state_demo_reads_four() {
        // The seeded workbench graph evaluates to 4 on construction.
        let s = NodeGraphWorkbenchState::default();
        let outs = s.graph.output_values();
        assert_eq!(outs.len(), 1);
        assert_eq!(outs[0].1, Some(4.0));
    }

    #[test]
    fn agent_set_constant_value_targets_selection() {
        let mut s = NodeGraphWorkbenchState::default();
        // The default selection is a Constant; set it to 10.
        s.agent_set("Constant value", &AgentValue::Float(10.0))
            .expect("set constant on the selected Constant node");
        s.run_eval();
        // Constant(10) + Constant(2) = 12.
        assert_eq!(s.graph.output_values()[0].1, Some(12.0));
    }

    #[test]
    fn agent_set_constant_value_rejects_non_constant_selection() {
        let mut s = NodeGraphWorkbenchState::default();
        // Add an Add node -> it becomes the selection -> setting a constant fails.
        let id = s.add_node(NodeKind::Add);
        assert_eq!(s.selected, Some(id));
        let err = s
            .agent_set("Constant value", &AgentValue::Float(1.0))
            .unwrap_err();
        assert!(err.contains("not a Constant"), "got: {err}");
    }

    #[test]
    fn agent_set_add_node_builds_the_graph() {
        let mut s = NodeGraphWorkbenchState::default();
        let before = s.graph.nodes.len();
        s.agent_set("Add node", &AgentValue::Str("multiply".to_string()))
            .expect("add a multiply node by id");
        assert_eq!(s.graph.nodes.len(), before + 1);
        assert_eq!(s.graph.nodes[&s.selected.unwrap()].kind, NodeKind::Multiply);
        // Unknown kind is fail-loud.
        assert!(s
            .agent_set("Add node", &AgentValue::Str("frobnicate".to_string()))
            .is_err());
    }

    #[test]
    fn agent_set_unknown_control_is_loud() {
        let mut s = NodeGraphWorkbenchState::default();
        assert!(s.agent_set("Nope", &AgentValue::Float(1.0)).is_err());
    }

    #[test]
    fn control_names_are_listed() {
        let names = NodeGraphWorkbenchState::agent_control_names();
        assert!(names.contains(&"Constant value"));
        assert!(names.contains(&"Add node"));
    }

    #[test]
    fn readout_reports_counts_and_output() {
        let s = NodeGraphWorkbenchState::default();
        let r = s.agent_readout().expect("readout always present");
        assert!(r.contains("4 nodes"), "got: {r}");
        assert!(r.contains("3 edges"), "got: {r}");
        assert!(r.contains("=4.0000"), "output value in readout; got: {r}");
    }

    #[test]
    fn run_bridge_helper_evaluates_through_app() {
        let mut app = ValenxApp::default();
        // Mutate a constant, then evaluate via the bridge `run`.
        app.nodegraph
            .agent_set("Constant value", &AgentValue::Float(5.0))
            .unwrap();
        run(&mut app);
        assert_eq!(app.nodegraph.graph.output_values()[0].1, Some(7.0));
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
            draw_nodegraph_workbench(app, ctx);
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
        assert!(!app.show_nodegraph_workbench);
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_nodegraph_workbench(&mut app, ctx);
        });
        // No panic = pass.
    }

    #[test]
    fn workbench_draws_when_shown() {
        let mut app = ValenxApp::default();
        app.show_nodegraph_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);
        assert!(!nodes.is_empty(), "a shown workbench produces a11y nodes");
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // The Constant-value DragValue is a SpinButton and must be `labelled_by`
        // its caption so an AI / screen reader can find it by caption text. The
        // default selection is a Constant, so the spin button is present.
        let mut app = ValenxApp::default();
        app.show_nodegraph_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let by_id: std::collections::HashMap<NodeId, &Node> =
            nodes.iter().map(|(id, n)| (*id, n)).collect();

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        assert!(
            !spin_buttons.is_empty(),
            "expected the Constant-value numeric control as a spin button"
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
        assert!(
            has_named_node(&nodes, "Constant value"),
            "the 'Constant value' caption is a named node"
        );
        assert!(
            has_named_node(&nodes, "Add node"),
            "the 'Add node' caption is a named node"
        );
    }
}
