//! # valenx-test-utils
//!
//! Tiny test-only helpers shared across Valenx's 120+ adapter unit
//! tests. Lives as its own crate (rather than as a `pub mod` in
//! `valenx-core`) so the helpers can be pulled in via
//! `[dev-dependencies]` without dragging core dependencies into the
//! adapter's runtime graph.
//!
//! Currently the only export is [`tempdir`] — every adapter's
//! `case_input.rs` test module used to hand-roll the same
//! `std::env::temp_dir().join("valenx-<name>-{nanos}")` boilerplate.
//! Centralising it here keeps the per-adapter test files focused on
//! schema parsing rather than scratch-directory plumbing.

#![forbid(unsafe_code)]

use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

/// Atomic counter that disambiguates `tempdir` calls landing on the
/// same nanosecond. Without it, two tests from the same adapter that
/// happen to schedule on a coarse-resolution timer (Windows hits
/// 100 ns granularity in some configurations) can collide and one
/// run will fail with `AlreadyExists` from `create_dir_all`.
static SEQ: AtomicU64 = AtomicU64::new(0);

/// Build (and create) a fresh per-test scratch directory underneath
/// the OS temp dir. Returns the directory path; the caller owns
/// cleanup (most adapter tests `let _ = std::fs::remove_dir_all(&d);`
/// at the end of the test body).
///
/// The `label` is embedded into the directory name so simultaneous
/// runs across adapters don't clobber each other's cases. Convention
/// is to pass the adapter id (`"bwa"`, `"clustalo"`, ...).
///
/// Filename shape:
/// `valenx-<label>-<unix_nanos>-<process_pid>-<atomic_seq>`
///
/// The atomic sequence number guarantees uniqueness even when the
/// system clock has coarse granularity (Windows) or when many tests
/// run in quick succession on the same thread.
///
/// Panics if the directory cannot be created — same contract as the
/// per-adapter helpers it replaces.
///
/// Callers should pass a non-empty label that doesn't contain path
/// separators. In debug builds these are checked via `debug_assert!`
/// to catch foot-guns; in release builds the label is used as-is.
pub fn tempdir(label: &str) -> PathBuf {
    debug_assert!(
        !label.is_empty(),
        "tempdir label must not be empty (use the adapter id, e.g. \"bwa\")"
    );
    debug_assert!(
        !label.contains('/') && !label.contains('\\'),
        "tempdir label `{label}` must not contain path separators"
    );
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system clock before UNIX_EPOCH")
        .as_nanos();
    let seq = SEQ.fetch_add(1, Ordering::Relaxed);
    let pid = std::process::id();
    let d = std::env::temp_dir().join(format!("valenx-{label}-{nanos}-{pid}-{seq}"));
    std::fs::create_dir_all(&d).expect("create scratch tempdir");
    d
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tempdir_creates_unique_directory_per_call() {
        let a = tempdir("unit-test");
        let b = tempdir("unit-test");
        assert_ne!(a, b, "two tempdir() calls must produce distinct paths");
        assert!(a.is_dir(), "first tempdir must exist on disk");
        assert!(b.is_dir(), "second tempdir must exist on disk");
        let _ = std::fs::remove_dir_all(&a);
        let _ = std::fs::remove_dir_all(&b);
    }

    #[test]
    fn tempdir_includes_label_in_path() {
        let d = tempdir("my-adapter");
        let name = d
            .file_name()
            .expect("dir has a name")
            .to_string_lossy()
            .to_string();
        assert!(
            name.starts_with("valenx-my-adapter-"),
            "path `{name}` must include the label between the prefix and the disambiguators"
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "must not be empty")]
    fn tempdir_rejects_empty_label_in_debug_builds() {
        let _ = tempdir("");
    }

    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "must not contain path separators")]
    fn tempdir_rejects_label_with_forward_slash_in_debug_builds() {
        let _ = tempdir("foo/bar");
    }

    #[test]
    fn tempdir_handles_rapid_sequential_calls_without_collision() {
        // Hammer the helper. Without the atomic seq, Windows' coarse
        // clock granularity used to collide here.
        let mut paths = Vec::with_capacity(100);
        for _ in 0..100 {
            paths.push(tempdir("hammer"));
        }
        let mut seen = std::collections::HashSet::new();
        for p in &paths {
            assert!(seen.insert(p.clone()), "duplicate path produced: {p:?}");
        }
        for p in paths {
            let _ = std::fs::remove_dir_all(&p);
        }
    }
}
