//! Integration tests for the `valenx-blast` binary.
//!
//! Build the binary via cargo's test fixture (CARGO_BIN_EXE_*) and
//! spawn it as a subprocess. We can NOT actually run BLAST+ in CI:
//! the runner won't have the binary or a database to point `--db`
//! at, and shipping a small bundled database is out of scope for
//! v0.1. Coverage stops at the parse / probe / spawn boundary:
//!
//! 1. `--help` / `-V` exit 0 with the expected output.
//! 2. Missing-binary path: with a stripped `PATH` env var the
//!    BLAST+ probe fails and we exit 3 with an actionable message.
//! 3. Empty-query path: a zero-byte `--db` argument plus an empty
//!    query exits 1 (content error) before any spawn.
//!
//! Happy-path BLAST+ runs land in the BLAST+ adapter's smoke tests
//! once that adapter exists.

use std::path::PathBuf;
use std::process::Command;

/// Path to the compiled `valenx-blast` binary, set by cargo at test
/// build time.
fn bin_path() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_valenx-blast"))
}

#[test]
fn help_prints_usage_and_exits_zero() {
    let out = Command::new(bin_path())
        .arg("--help")
        .output()
        .expect("spawn binary");
    assert!(
        out.status.success(),
        "expected exit 0, got {:?}",
        out.status
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(stdout.contains("valenx-blast"));
    assert!(stdout.contains("USAGE"));
    assert!(stdout.contains("--db"));
    assert!(stdout.contains("--evalue"));
}

#[test]
fn version_prints_and_exits_zero() {
    let out = Command::new(bin_path())
        .arg("-V")
        .output()
        .expect("spawn binary");
    assert!(
        out.status.success(),
        "expected exit 0, got {:?}",
        out.status
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(stdout.contains("valenx-blast"));
    assert!(stdout.contains('v'), "expected version marker: {stdout}");
}

#[test]
fn missing_binary_exits_three_with_actionable_message() {
    // Strip PATH (and PATHEXT on Windows) entirely so `find_on_path`
    // returns None even on dev machines that have BLAST+ installed.
    // Empty PATH is the cleanest "binary not present" signal.
    let tmp = std::env::temp_dir().join(format!(
        "valenx-blast-cli-noblast-{}.fasta",
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::write(&tmp, ">x\nMEEPQSDPSV\n").expect("write fixture");
    let out = Command::new(bin_path())
        .arg(&tmp)
        .arg("--db")
        .arg("nr")
        .env("PATH", "")
        .env("PATHEXT", "")
        .output()
        .expect("spawn binary");
    assert_eq!(
        out.status.code(),
        Some(3),
        "expected exit 3, got {:?} (stderr: {})",
        out.status,
        String::from_utf8_lossy(&out.stderr)
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("BLAST+") || stderr.contains("blast"),
        "stderr should mention BLAST+: {stderr}"
    );
    assert!(
        stderr.contains("PATH"),
        "stderr should hint about PATH: {stderr}"
    );
    let _ = std::fs::remove_file(&tmp);
}

#[test]
fn empty_query_exits_one_before_spawn() {
    // A zero-byte query trips load_query's "query is empty" branch
    // before we even attempt to find_on_path. Exit code 1 (content)
    // — this is the expected shape on bad input.
    let tmp = std::env::temp_dir().join(format!(
        "valenx-blast-cli-empty-{}.fasta",
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::write(&tmp, "").expect("write empty fixture");
    let out = Command::new(bin_path())
        .arg(&tmp)
        .arg("--db")
        .arg("nr")
        .output()
        .expect("spawn binary");
    assert_eq!(
        out.status.code(),
        Some(1),
        "expected exit 1, got {:?} (stderr: {})",
        out.status,
        String::from_utf8_lossy(&out.stderr)
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("empty"),
        "stderr should mention empty query: {stderr}"
    );
    let _ = std::fs::remove_file(&tmp);
}
