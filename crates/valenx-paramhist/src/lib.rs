//! # valenx-paramhist
//!
//! NaroCAD-style parametric history alternatives — Phase 67.
//!
//! Modules:
//! - [`history::History`] / [`history::HistEntry`] — dependency-
//!   carrying entry list.
//! - [`history::topological_sort`] — Kahn's algorithm.
//! - [`history::dirty_set`] — downstream propagation.
//! - [`history::partial_rebuild`] — re-run only dirty entries.
//! - [`History::insert_before`] / [`History::move_to`] —
//!   DAG-preserving editing.
//! - [`panel::ParamHistPanelState`] — UI envelope with DAG layout.
//! - [`error`] — typed [`ParamHistError`].

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod history;
pub mod panel;

pub use error::{ErrorCategory, ParamHistError};
pub use history::{dirty_set, partial_rebuild, topological_sort, HistEntry, History, RebuildResult};
pub use panel::{NodeLayout, ParamHistPanelState};

#[cfg(test)]
mod tests {
    use super::*;

    fn build_simple() -> History {
        let mut h = History::new();
        h.push(HistEntry::new("a", vec![])).unwrap();
        h.push(HistEntry::new("b", vec![0])).unwrap();
        h.push(HistEntry::new("c", vec![1])).unwrap();
        h.push(HistEntry::new("d", vec![0])).unwrap();
        h
    }

    #[test]
    fn push_records_indices() {
        let h = build_simple();
        assert_eq!(h.entries.len(), 4);
        assert_eq!(h.entries[2].dependencies, vec![1]);
    }

    #[test]
    fn push_rejects_forward_reference() {
        let mut h = History::new();
        let r = h.push(HistEntry::new("a", vec![0]));
        assert!(matches!(r, Err(ParamHistError::InvalidMove(_))));
    }

    #[test]
    fn topological_sort_acyclic() {
        let h = build_simple();
        let order = topological_sort(&h.entries).unwrap();
        // 'a' first; 'c' after 'b'.
        let pos = |name: &str| {
            order
                .iter()
                .position(|i| h.entries[*i].name == name)
                .unwrap()
        };
        assert!(pos("a") < pos("b"));
        assert!(pos("b") < pos("c"));
        assert!(pos("a") < pos("d"));
    }

    #[test]
    fn topological_sort_handles_duplicate_dependency() {
        // A duplicate dependency previously made in_deg double-count, so the
        // node never reached in-degree 0 → a false Cycle for an acyclic graph.
        let mut h = History::new();
        h.push(HistEntry::new("a", vec![])).unwrap();
        h.push(HistEntry::new("b", vec![0, 0])).unwrap(); // dep 0 listed twice
        let order =
            topological_sort(&h.entries).expect("acyclic graph must not report a cycle");
        assert_eq!(order.len(), 2);
    }

    #[test]
    fn relayout_tolerates_dangling_dependency() {
        // A dangling dependency index (from a corrupt project / direct
        // construction) must not panic relayout's `layer[*d]`.
        let mut s = crate::panel::ParamHistPanelState::new();
        s.history.entries = vec![HistEntry::new("x", vec![99])]; // dep 99 out of range
        s.relayout(); // must not panic
    }

    #[test]
    fn dirty_set_propagates_downstream() {
        let h = build_simple();
        let d = dirty_set(&h.entries, 0);
        assert!(d.contains(&0));
        assert!(d.contains(&1));
        assert!(d.contains(&2));
        assert!(d.contains(&3));
        let d2 = dirty_set(&h.entries, 1);
        assert!(d2.contains(&1));
        assert!(d2.contains(&2));
        assert!(!d2.contains(&3));
    }

    #[test]
    fn partial_rebuild_marks_up_to_date() {
        let mut h = build_simple();
        let d = dirty_set(&h.entries, 1);
        let r = partial_rebuild(&mut h.entries, &d).unwrap();
        // Two entries dirty (1 and 2), in topo order.
        assert_eq!(r.len(), 2);
        for entry in &h.entries {
            if d.contains(&h.entries.iter().position(|x| x.name == entry.name).unwrap()) {
                assert!(entry.up_to_date);
            }
        }
    }

    #[test]
    fn insert_before_bumps_dependencies() {
        let mut h = build_simple();
        h.insert_before(1, HistEntry::new("x", vec![0])).unwrap();
        // Now 'x' is at index 1; 'b' moved to 2, its dep should be 0
        // (originally referenced 0, not bumped); 'c' at 3 should
        // reference 2.
        assert_eq!(h.entries[1].name, "x");
        assert_eq!(h.entries[2].name, "b");
        assert_eq!(h.entries[2].dependencies, vec![0]);
        assert_eq!(h.entries[3].dependencies, vec![2]);
    }

    #[test]
    fn insert_before_rejects_forward_dep() {
        let mut h = build_simple();
        let r = h.insert_before(1, HistEntry::new("x", vec![1]));
        assert!(matches!(r, Err(ParamHistError::InvalidMove(_))));
    }

    #[test]
    fn move_to_same_position_no_op() {
        let mut h = build_simple();
        h.move_to(2, 2).unwrap();
        assert_eq!(h.entries[2].name, "c");
    }

    #[test]
    fn move_to_rejects_back_edge() {
        let mut h = build_simple();
        // Move 'b' (idx 1) to 3 would put it after 'c' (idx 2 -> 1).
        // 'c' depends on 'b' so it's fine. But 'a' (idx 0) moving
        // to 2 would put it after 'b' which depends on 'a' --> bad.
        let r = h.move_to(0, 2);
        assert!(matches!(r, Err(ParamHistError::InvalidMove(_))));
    }

    #[test]
    fn move_to_remaps_dependencies() {
        let mut h = build_simple();
        // Swap 'b' and 'd' — both depend on 'a' so they can re-order.
        h.move_to(3, 1).unwrap();
        assert_eq!(h.entries[1].name, "d");
        // 'b' is now at index 2; still depends on 'a' (index 0).
        assert_eq!(h.entries[2].name, "b");
        assert_eq!(h.entries[2].dependencies, vec![0]);
    }

    #[test]
    fn panel_relayout_assigns_layers() {
        let mut p = ParamHistPanelState::new();
        p.history = build_simple();
        p.relayout();
        // a -> layer 0; b -> 1; c -> 2; d -> 1.
        assert_eq!(p.layout[0].layer, 0);
        assert_eq!(p.layout[1].layer, 1);
        assert_eq!(p.layout[2].layer, 2);
        assert_eq!(p.layout[3].layer, 1);
    }

    #[test]
    fn panel_select_updates_status() {
        let mut p = ParamHistPanelState::new();
        p.select(Some(0));
        assert_eq!(p.selection, Some(0));
        assert!(p.last_status.is_some());
    }
}
