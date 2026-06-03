//! Test-only binary: install the Valenx panic hook + deliberately
//! panic. The integration test in
//! `tests/panic_hook_integration.rs` spawns this and verifies a
//! crash report lands on disk.
//!
//! Reads the crashes directory from `VALENX_CRASH_DIR` so the
//! test driver picks a fresh path per run.

use std::path::PathBuf;

fn main() {
    let dir = std::env::var_os("VALENX_CRASH_DIR")
        .map(PathBuf::from)
        .expect("VALENX_CRASH_DIR not set — this binary is test-only");
    valenx_crash_reporter::install_panic_hook(dir, env!("CARGO_PKG_VERSION").into());
    // Trigger a panic the test can match against.
    panic!("deliberate test panic from panic-bin");
}
