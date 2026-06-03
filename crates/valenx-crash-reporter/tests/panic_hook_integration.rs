//! Integration test: install the panic hook in a subprocess that
//! deliberately panics, then assert a [`CrashReport`] lands on
//! disk in the configured directory.
//!
//! Why a subprocess: panic hooks are global state and the parent
//! test process needs to keep running. Spawning a tiny example
//! binary (built at test time via the `examples/` mechanism) lets
//! us isolate the panic to a child whose exit status we observe.

use std::path::PathBuf;
use std::process::Command;

use valenx_crash_reporter::CrashReport;

fn example_binary() -> PathBuf {
    // Cargo populates CARGO_BIN_EXE_<name> for `[[bin]]` targets;
    // for `[[example]]` targets it sets CARGO_BIN_EXE_<example>
    // when --example builds run. Because integration tests don't
    // automatically build examples, we use a `[[bin]]` shim
    // instead — see `src/bin/panic-bin.rs`.
    PathBuf::from(env!("CARGO_BIN_EXE_panic-bin"))
}

#[test]
fn panic_hook_writes_report_to_configured_dir() {
    // Pick a fresh per-test crashes dir so we can assert exactly
    // one report appears.
    let dir = std::env::temp_dir().join(format!(
        "valenx-crash-hook-test-{}",
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));

    let out = Command::new(example_binary())
        .env("VALENX_CRASH_DIR", &dir)
        .output()
        .expect("spawn panic-bin");

    // The child panicked — exit status must be non-zero. (Rust's
    // default panic behaviour is abort or unwind-then-exit-101.)
    assert!(
        !out.status.success(),
        "child exited cleanly; stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    // Exactly one report on disk.
    let reports = CrashReport::load_all(&dir).expect("load reports");
    assert_eq!(reports.len(), 1, "expected one report; got: {reports:#?}");
    let report = &reports[0].1;
    // Confirm the panic message survived.
    assert!(
        report.message.contains("deliberate test panic"),
        "got: {report:?}"
    );
    // Location captured.
    assert!(report.location.is_some());
    // Valenx version surfaced from the bin's CARGO_PKG_VERSION.
    assert!(
        !report.valenx_version.is_empty(),
        "version missing: {report:?}"
    );

    let _ = std::fs::remove_dir_all(&dir);
}
