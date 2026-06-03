//! Generic undo / redo snapshot history.
//!
//! Panels with stateful text-edit or numeric inputs (Sequence editor,
//! Alignment input, RNA Designer wizard, Sketcher) wrap their edit
//! state in a [`History<T>`]. Each call to [`History::record`]
//! pushes a snapshot of the current value; [`History::undo`] pops
//! and returns the previous snapshot; [`History::redo`] reapplies a
//! previously-undone snapshot.
//!
//! The implementation deliberately uses simple `Vec` stacks rather
//! than the textbook "patch + delta" command pattern. The state
//! objects in question are small — a sequence string is ≤ ~100 kB
//! in normal use, a sketch is tens of kilobytes — so even a 64-entry
//! cap on each stack costs at most a few megabytes per panel. Patch
//! deltas would buy memory at the cost of code surface that nothing
//! else needs.
//!
//! ## Coverage
//!
//! The panels with reversible edits wire the history struct in;
//! the read-only result-display panels (composition results, ORF
//! tables, PCR amplicons) don't. The cheat-sheet overlay's per-panel
//! row column declares which panels have undo/redo so users aren't
//! surprised when Ctrl+Z does nothing in (e.g.) the Genomics VCF
//! summary.

/// Bounded undo / redo stack over a `Clone` state value.
///
/// Each `record(state)` push moves the current value onto the undo
/// stack and clears the redo stack (the conventional editor
/// semantics: making a new edit branches off the timeline). Stacks
/// are capped at [`MAX_DEPTH`] entries — when full, the oldest
/// snapshot is discarded so memory stays bounded.
#[derive(Clone, Debug, Default)]
pub struct History<T: Clone + PartialEq> {
    undo_stack: Vec<T>,
    redo_stack: Vec<T>,
}

/// Per-stack snapshot cap. 64 deep is enough for a long editing
/// session without bounding RAM use even for fat values (a 100 kB
/// sequence × 64 × 2 = ~13 MB worst case, comfortably under any
/// realistic limit).
pub const MAX_DEPTH: usize = 64;

impl<T: Clone + PartialEq> History<T> {
    /// Empty history — no undo / no redo available.
    pub fn new() -> Self {
        Self {
            undo_stack: Vec::new(),
            redo_stack: Vec::new(),
        }
    }

    /// Push `before` onto the undo stack, then clear the redo stack.
    ///
    /// Callers invoke this **before** mutating their state so the
    /// snapshot captures the value the user is about to leave.
    /// Duplicate snapshots (identical to the most-recent entry) are
    /// dropped to keep undo-Z from making the user mash Z eight
    /// times to escape a typing burst that produced no net change.
    pub fn record(&mut self, before: T) {
        if let Some(top) = self.undo_stack.last() {
            if top == &before {
                return;
            }
        }
        self.undo_stack.push(before);
        if self.undo_stack.len() > MAX_DEPTH {
            // Drop the oldest entry (front of the vec). pop_front
            // is O(n) on a Vec — but at MAX_DEPTH = 64 that's a
            // 64-element memmove, negligible compared to the user's
            // typing speed.
            self.undo_stack.remove(0);
        }
        // A new edit invalidates the redo timeline — match the
        // canonical editor behaviour.
        self.redo_stack.clear();
    }

    /// Pop the last undo snapshot and push `current` onto the redo
    /// stack, returning the snapshot so the caller can swap their
    /// state back to it. Returns `None` if there's nothing to undo.
    pub fn undo(&mut self, current: T) -> Option<T> {
        let prev = self.undo_stack.pop()?;
        self.redo_stack.push(current);
        Some(prev)
    }

    /// Pop the last redo snapshot, pushing `current` onto the undo
    /// stack so a subsequent undo can reverse the redo. Returns
    /// `None` if there's nothing to redo.
    pub fn redo(&mut self, current: T) -> Option<T> {
        let next = self.redo_stack.pop()?;
        self.undo_stack.push(current);
        Some(next)
    }

    /// `true` if an undo would do something.
    pub fn can_undo(&self) -> bool {
        !self.undo_stack.is_empty()
    }

    /// `true` if a redo would do something.
    pub fn can_redo(&self) -> bool {
        !self.redo_stack.is_empty()
    }

    /// Drop every undo / redo entry. Called when a panel loads a
    /// fresh document so the user can't "undo back into" the previous
    /// document's history.
    pub fn clear(&mut self) {
        self.undo_stack.clear();
        self.redo_stack.clear();
    }

    /// Number of undo snapshots stored. Used by tests + the
    /// debug-tools panel.
    pub fn undo_depth(&self) -> usize {
        self.undo_stack.len()
    }

    /// Number of redo snapshots stored.
    pub fn redo_depth(&self) -> usize {
        self.redo_stack.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_history_has_no_undo_or_redo() {
        let h: History<String> = History::new();
        assert!(!h.can_undo());
        assert!(!h.can_redo());
        assert_eq!(h.undo_depth(), 0);
        assert_eq!(h.redo_depth(), 0);
    }

    #[test]
    fn record_pushes_snapshots() {
        let mut h: History<String> = History::new();
        h.record("a".into());
        h.record("b".into());
        h.record("c".into());
        assert_eq!(h.undo_depth(), 3);
        assert!(h.can_undo());
    }

    #[test]
    fn record_drops_consecutive_duplicates() {
        // Identical consecutive snapshots collapse so a long typing
        // burst that returns to the same string doesn't pollute the
        // undo stack.
        let mut h: History<String> = History::new();
        h.record("a".into());
        h.record("a".into());
        h.record("a".into());
        assert_eq!(h.undo_depth(), 1);
    }

    #[test]
    fn undo_returns_previous_snapshot_and_populates_redo() {
        let mut h: History<String> = History::new();
        h.record("a".into());
        h.record("b".into());
        // We are currently at "c" — undoing returns "b" and pushes
        // "c" onto the redo stack.
        let prev = h.undo("c".into()).expect("undo returns prev");
        assert_eq!(prev, "b");
        assert_eq!(h.undo_depth(), 1);
        assert_eq!(h.redo_depth(), 1);
    }

    #[test]
    fn redo_round_trips_with_undo() {
        let mut h: History<String> = History::new();
        h.record("a".into());
        // current = "b"; undo → state = "a" and redo stack has "b".
        let after_undo = h.undo("b".into()).unwrap();
        assert_eq!(after_undo, "a");
        // Redo from "a" gives us "b" back; undo stack now has "a".
        let after_redo = h.redo("a".into()).unwrap();
        assert_eq!(after_redo, "b");
        assert_eq!(h.undo_depth(), 1);
        assert_eq!(h.redo_depth(), 0);
    }

    #[test]
    fn recording_after_undo_clears_redo_stack() {
        let mut h: History<String> = History::new();
        h.record("a".into());
        let _ = h.undo("b".into()); // redo stack now has "b"
        assert_eq!(h.redo_depth(), 1);
        // A fresh edit branches off the timeline — redo is invalid.
        h.record("a-prime".into());
        assert_eq!(h.redo_depth(), 0);
    }

    #[test]
    fn undo_on_empty_returns_none() {
        let mut h: History<String> = History::new();
        assert!(h.undo("current".into()).is_none());
    }

    #[test]
    fn stack_is_bounded_to_max_depth() {
        let mut h: History<u32> = History::new();
        for i in 0..(MAX_DEPTH * 2) {
            h.record(i as u32);
        }
        assert_eq!(h.undo_depth(), MAX_DEPTH);
        // The earliest values were dropped — undo should give us
        // one of the later ones.
        let v = h.undo(u32::MAX).unwrap();
        assert!(v >= MAX_DEPTH as u32);
    }

    #[test]
    fn clear_drops_both_stacks() {
        let mut h: History<String> = History::new();
        h.record("a".into());
        let _ = h.undo("b".into());
        assert!(h.can_undo() || h.can_redo());
        h.clear();
        assert!(!h.can_undo());
        assert!(!h.can_redo());
    }
}
