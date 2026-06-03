//! UI panel envelope — DAG-visualisation state.

use crate::history::History;

/// Layout hint for a single node in the DAG visualisation.
#[derive(Clone, Debug, PartialEq)]
pub struct NodeLayout {
    /// Entry index.
    pub index: usize,
    /// Layer (depth from root = longest dependency chain).
    pub layer: usize,
    /// Slot within the layer.
    pub slot: usize,
}

/// Workbench panel state.
#[derive(Default)]
pub struct ParamHistPanelState {
    /// History being browsed.
    pub history: History,
    /// Cached layout from the last [`Self::relayout`] call.
    pub layout: Vec<NodeLayout>,
    /// User-selected entry index.
    pub selection: Option<usize>,
    /// Status message.
    pub last_status: Option<String>,
    /// Error message.
    pub last_error: Option<String>,
}

impl ParamHistPanelState {
    /// New empty.
    pub fn new() -> Self {
        Self::default()
    }

    /// Recompute layer + slot positions for every entry.
    pub fn relayout(&mut self) {
        let n = self.history.entries.len();
        let mut layer = vec![0usize; n];
        for (i, e) in self.history.entries.iter().enumerate() {
            let max_dep_layer = e
                .dependencies
                .iter()
                .map(|d| layer[*d] + 1)
                .max()
                .unwrap_or(0);
            layer[i] = max_dep_layer;
        }
        // Slot per layer.
        let mut slot_counter: std::collections::BTreeMap<usize, usize> = Default::default();
        let mut out = Vec::with_capacity(n);
        for (i, l) in layer.iter().enumerate() {
            let slot = *slot_counter.entry(*l).or_insert(0);
            *slot_counter.entry(*l).or_insert(0) += 1;
            out.push(NodeLayout {
                index: i,
                layer: *l,
                slot,
            });
        }
        self.layout = out;
    }

    /// Select a node.
    pub fn select(&mut self, i: Option<usize>) {
        self.selection = i;
        if let Some(i) = i {
            self.last_status = Some(format!("selected entry {i}"));
            self.last_error = None;
        }
    }

    /// Status setter.
    pub fn set_status(&mut self, s: impl Into<String>) {
        self.last_status = Some(s.into());
        self.last_error = None;
    }

    /// Error setter.
    pub fn set_error(&mut self, s: impl Into<String>) {
        self.last_error = Some(s.into());
        self.last_status = None;
    }
}
