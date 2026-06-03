//! Per-case run + sweep history persistence. The browser tree shows
//! ✓/✗ badges next to each case based on the last recorded outcome,
//! and those badges have to survive an app restart — so the maps
//! land on disk after every run / sweep completes.

use std::io::Read;

use crate::settings_io::MAX_STATE_FILE_BYTES;
use crate::state_paths::{atomic_write, run_history_path, sweep_history_path};
use crate::types::{RunHistoryEntry, SweepHistoryEntry};

/// Load the persisted per-case `run_history` map from disk. Returns
/// `None` if the state file doesn't exist or isn't readable / parseable.
/// Treated as "no history" rather than fatal — silently empty is
/// always safer than crashing on a corrupted state file.
///
/// Round-8 hardening: bounded read via
/// [`crate::settings_io::MAX_STATE_FILE_BYTES`]. A run-history JSON
/// past the cap is rejected as "no history" so a hostile / corrupted
/// state file can't OOM the host on startup.
pub fn load_run_history_from_state_dir(
) -> Option<std::collections::BTreeMap<String, RunHistoryEntry>> {
    let path = run_history_path()?;
    let meta = std::fs::metadata(&path).ok()?;
    if meta.len() > MAX_STATE_FILE_BYTES as u64 {
        return None;
    }
    let mut buf = Vec::with_capacity(meta.len().min(MAX_STATE_FILE_BYTES as u64) as usize);
    let mut file = std::fs::File::open(&path).ok()?;
    file.by_ref()
        .take(MAX_STATE_FILE_BYTES as u64)
        .read_to_end(&mut buf)
        .ok()?;
    serde_json::from_slice(&buf).ok()
}

/// Persist the per-case `run_history` map to disk. Best-effort —
/// failures (no state dir, can't create, can't write) are swallowed
/// silently because losing history-badge persistence isn't worth
/// surfacing as a UI error every time it happens (which on a fresh
/// install would be every single run until the user happens to
/// notice).
pub fn save_run_history_to_state_dir(
    history: &std::collections::BTreeMap<String, RunHistoryEntry>,
) {
    let Some(path) = run_history_path() else {
        return;
    };
    if let Ok(text) = serde_json::to_string_pretty(history) {
        // atomic_write writes to <path>.tmp and renames, so a crash
        // mid-write doesn't truncate the previous history map.
        let _ = atomic_write(&path, &text);
    }
}

/// Load the persisted per-case `sweep_history` map from disk.
/// Returns `None` when the file doesn't exist or fails to parse —
/// mirrors `load_run_history_from_state_dir`'s lenient policy so a
/// fresh-install / corrupted-history doesn't crash the launcher.
///
/// Round-8 hardening: bounded read via
/// [`crate::settings_io::MAX_STATE_FILE_BYTES`] — same OOM-guard as
/// the run-history loader.
pub fn load_sweep_history_from_state_dir(
) -> Option<std::collections::BTreeMap<String, SweepHistoryEntry>> {
    let path = sweep_history_path()?;
    let meta = std::fs::metadata(&path).ok()?;
    if meta.len() > MAX_STATE_FILE_BYTES as u64 {
        return None;
    }
    let mut buf = Vec::with_capacity(meta.len().min(MAX_STATE_FILE_BYTES as u64) as usize);
    let mut file = std::fs::File::open(&path).ok()?;
    file.by_ref()
        .take(MAX_STATE_FILE_BYTES as u64)
        .read_to_end(&mut buf)
        .ok()?;
    serde_json::from_slice(&buf).ok()
}

/// Persist the per-case `sweep_history` map to disk. Best-effort —
/// matches the run-history pattern.
pub fn save_sweep_history_to_state_dir(
    history: &std::collections::BTreeMap<String, SweepHistoryEntry>,
) {
    let Some(path) = sweep_history_path() else {
        return;
    };
    if let Ok(text) = serde_json::to_string_pretty(history) {
        // atomic_write writes to <path>.tmp and renames, so a crash
        // mid-write doesn't truncate the previous sweep history.
        let _ = atomic_write(&path, &text);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn run_history_entry_serializes_round_trip() {
        // Persistence depends on RunHistoryEntry serialising cleanly
        // through serde_json. Lock that down here so a future field
        // addition can't silently break the saved-state file format.
        let mut map = std::collections::BTreeMap::new();
        map.insert(
            "cfd-steady".to_string(),
            RunHistoryEntry {
                succeeded: true,
                wall_time: std::time::Duration::from_millis(1234),
                converged: Some(true),
            },
        );
        let json = serde_json::to_string(&map).expect("serialize");
        let parsed: std::collections::BTreeMap<String, RunHistoryEntry> =
            serde_json::from_str(&json).expect("deserialize");
        assert_eq!(parsed, map);
    }

    #[test]
    fn sweep_history_entry_round_trips_through_serde() {
        let entry = SweepHistoryEntry {
            planned: 32,
            succeeded: 24,
            failed: 8,
            workdir: PathBuf::from("/tmp/valenx-sweep-airfoil-1700000000"),
            completed_at: "2026-04-25T12:00:00Z".into(),
        };
        let json = serde_json::to_string(&entry).unwrap();
        let parsed: SweepHistoryEntry = serde_json::from_str(&json).unwrap();
        assert_eq!(parsed, entry);
    }

    #[test]
    fn save_and_load_sweep_history_round_trips_through_disk() {
        // Skip when no state dir is available (rare on CI hosts
        // without HOME / APPDATA env vars). The path-resolver
        // returns None in that case and the helpers no-op.
        let Some(_) = sweep_history_path() else {
            return;
        };
        let mut map = std::collections::BTreeMap::new();
        let entry = SweepHistoryEntry {
            planned: 5,
            succeeded: 4,
            failed: 1,
            workdir: PathBuf::from("/tmp/sweep-test"),
            completed_at: "2026-04-25T13:00:00Z".into(),
        };
        map.insert(
            format!(
                "test-case-{}",
                std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .unwrap()
                    .as_nanos()
            ),
            entry.clone(),
        );
        save_sweep_history_to_state_dir(&map);
        let loaded = load_sweep_history_from_state_dir().expect("must load");
        // The on-disk file gets shared across test runs; check that
        // the entry we just saved is at least present in the
        // round-trip without asserting exact map equality.
        assert!(loaded.values().any(|v| v == &entry));
    }
}
