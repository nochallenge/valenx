//! Per-operation `generate` functions. Each submodule emits a
//! [`crate::Toolpath`] for a single op kind.
//!
//! ## Cross-cutting conventions
//!
//! Every op generator follows the same conventions so emitted
//! toolpaths chain together cleanly:
//!
//! - **Safe-Z** — the first move in every op is a rapid to `stock
//!   .top_z() + safe_z_clearance`. The op never assumes the tool is
//!   already at safe-Z when it starts.
//! - **Plunge-rate** — all Z-decreasing moves use [`MoveKind::Plunge`]
//!   and the params' `plunge_feed` (typically slower than the cut
//!   feed). All XY-only moves at depth use [`MoveKind::Cut`] and the
//!   `feed_mm_per_min`.
//! - **Rapid traversal** — moves above the stock between cutting
//!   regions are [`MoveKind::Rapid`].
//! - **Final rapid** — every op ends with a rapid up to safe-Z so
//!   the next op in a chain can start from a known state.
//!
//! The [`safe_z_for`] helper exposes the safe-Z computation so
//! callers (e.g. the host UI) can preview clearance heights.
//!
//! The [`crate::Toolpath::concatenate`] method chains the per-op
//! toolpaths into a single program for postprocessing.

pub mod drill;
pub mod face;
pub mod pocket;
pub mod profile;

// Phase 17A — Adaptive clearing + new entry primitives.
pub mod adaptive_clearing;
pub mod adaptive_constant_engagement;
pub mod helix_bore;
pub mod peck_drill_full;
pub mod plunge_rough;
pub mod ramp_entry;

// Phase 17B — More 2D + 3D ops.
pub mod contour_2d;
pub mod contour_3d;
pub mod engrave;
pub mod rest_machining;
pub mod scribe;
pub mod slot;
pub mod spiral_pocket;
pub mod thread_mill;
pub mod trochoidal_slot;
pub mod waterline_3d;

// Phase 17E — 5-axis ops.
pub mod fixed_axis_indexing;
pub mod tcp_5ax_contour;

use crate::{stock::Stock, toolpath::MoveKind};

/// Compute the safe-Z plane for an op: `stock.top_z() +
/// safe_z_clearance`. Every op uses this as the start/end Z.
pub fn safe_z_for(stock: &Stock, safe_z_clearance: f64) -> f64 {
    stock.top_z() + safe_z_clearance
}

/// Returns `true` if the given [`MoveKind`] is a Z-decreasing
/// (plunge) move that should use the plunge feed rather than the
/// XY cut feed.
pub fn is_plunge_move(kind: MoveKind) -> bool {
    matches!(kind, MoveKind::Plunge)
}

/// Hard cap on the number of Z-stepping passes a single op may
/// generate. Real toolpaths sit in the 1-100 pass range; 10 000 is
/// far beyond any plausible cut and small enough to defend against
/// the `step_down = f64::MIN_POSITIVE` attack — without the cap,
/// `(depth / MIN_POSITIVE).ceil() as usize` saturates to
/// `usize::MAX`, then `for k in 1..=n_passes` would loop ~2^64
/// times before terminating.
///
/// Op generators that use the `depth / step_down` formula must call
/// `compute_n_passes` to apply this cap consistently.
pub const MAX_N_PASSES: usize = 10_000;

/// Compute the number of Z-stepping passes from `depth` and
/// `step_down`, returning an error when the ratio would exceed
/// [`MAX_N_PASSES`]. Inputs are assumed to have already passed the
/// op's own finite+positive validation; this function only adds the
/// cap-vs-ratio check that the bare arithmetic would miss.
pub(crate) fn compute_n_passes(
    depth: f64,
    step_down: f64,
    op_name: &'static str,
) -> Result<usize, crate::CamError> {
    let ratio = depth / step_down;
    if ratio > MAX_N_PASSES as f64 {
        return Err(crate::CamError::BadOperation {
            name: op_name.into(),
            reason: format!(
                "depth / step_down ratio {ratio} exceeds {MAX_N_PASSES} pass cap — \
                 step_down is implausibly small relative to depth"
            ),
        });
    }
    Ok((ratio.ceil() as usize).max(1))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stock::Stock;

    #[test]
    fn safe_z_adds_clearance() {
        let s = Stock::default();
        // Default stock top is +10.0.
        assert!((safe_z_for(&s, 5.0) - 15.0).abs() < 1e-9);
    }

    #[test]
    fn is_plunge_classifies_correctly() {
        assert!(is_plunge_move(MoveKind::Plunge));
        assert!(!is_plunge_move(MoveKind::Cut));
        assert!(!is_plunge_move(MoveKind::Rapid));
    }
}
