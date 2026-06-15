//! Flow-arrangement enum shared by the LMTD and effectiveness-NTU
//! methods.

use serde::{Deserialize, Serialize};

/// Relative direction of the two streams in a two-stream heat
/// exchanger.
///
/// Only the two arrangements with closed-form, single-parameter
/// effectiveness relations are modelled here:
///
/// - [`FlowArrangement::Counterflow`] — the hot and cold streams run in
///   opposite directions. For identical terminal temperatures this
///   yields the largest log-mean temperature difference and the highest
///   attainable effectiveness of any two-stream layout.
/// - [`FlowArrangement::ParallelFlow`] — the streams run in the same
///   direction (also called co-current). Both outlet temperatures
///   approach a common intermediate value, capping the effectiveness
///   below the counterflow value.
///
/// Cross-flow and multi-pass shell-and-tube layouts (which need a
/// correction factor `F` or two-parameter relations) are intentionally
/// out of scope.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum FlowArrangement {
    /// Streams flow in opposite directions (counter-current).
    Counterflow,
    /// Streams flow in the same direction (co-current).
    ParallelFlow,
}

impl FlowArrangement {
    /// Short UI / log label.
    pub fn label(self) -> &'static str {
        match self {
            Self::Counterflow => "Counterflow",
            Self::ParallelFlow => "ParallelFlow",
        }
    }
}
