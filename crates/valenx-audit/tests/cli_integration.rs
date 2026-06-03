//! Integration tests for the `valenx-audit` binary.
//!
//! Build the binary via cargo's CARGO_BIN_EXE_* env var, spawn it
//! as a subprocess, and verify text + JSON output for both `verify`
//! and `tail` subcommands. End-to-end coverage that the in-binary
//! unit tests can't reach (those exercise parse_args / library
//! calls; these exercise the actual exe + stdout / stderr / exit
//! pipeline).

use std::path::PathBuf;
use std::process::Command;

use valenx_audit::{AuditActor, AuditEntry, AuditWriter};

/// Path to the compiled binary, set by cargo at test build time.
fn bin_path() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_valenx-audit"))
}

/// Write a small valid audit log to a unique tempfile and return
/// (workdir, log_path). Caller cleans up.
fn write_fixture_log(n_entries: usize) -> (PathBuf, PathBuf) {
    let tmp = std::env::temp_dir().join(format!(
        "valenx-audit-cli-{}",
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir_all(&tmp).unwrap();
    let log_path = tmp.join("audit.log.jsonl");
    let writer = AuditWriter::new(log_path.clone());
    for i in 0..n_entries {
        let entry = AuditEntry {
            timestamp: format!("2026-04-28T00:00:{i:02}Z"),
            actor: AuditActor {
                id: "alice".to_string(),
                session_id: None,
            },
            action: format!("test.action.{i}"),
            target: serde_json::json!({"kind": "case", "case": format!("case-{i}")}),
            context: serde_json::json!({}),
            prev_hash: String::new(), // writer fills this in
        };
        writer.append(entry).expect("append");
    }
    (tmp, log_path)
}

#[test]
fn help_prints_usage_and_exits_zero() {
    let out = Command::new(bin_path())
        .arg("help")
        .output()
        .expect("spawn binary");
    assert!(out.status.success());
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(stdout.contains("valenx-audit"));
    assert!(stdout.contains("verify"));
    assert!(stdout.contains("tail"));
}

#[test]
fn verify_succeeds_on_well_formed_log() {
    let (tmp, log) = write_fixture_log(3);
    let out = Command::new(bin_path())
        .arg("verify")
        .arg(&log)
        .output()
        .expect("spawn binary");
    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(stdout.contains("verified"));
    assert!(stdout.contains("3 entries"));
    let _ = std::fs::remove_dir_all(&tmp);
}

#[test]
fn verify_json_emits_parseable_object() {
    let (tmp, log) = write_fixture_log(2);
    let out = Command::new(bin_path())
        .arg("verify")
        .arg("--json")
        .arg(&log)
        .output()
        .expect("spawn binary");
    assert!(out.status.success());
    let stdout = String::from_utf8(out.stdout).expect("utf-8");
    let v: serde_json::Value = serde_json::from_str(&stdout).expect("parseable JSON");
    assert_eq!(v["ok"], true);
    assert_eq!(v["entries_verified"], 2);
    let _ = std::fs::remove_dir_all(&tmp);
}

#[test]
fn tail_text_prints_one_line_per_entry() {
    let (tmp, log) = write_fixture_log(5);
    let out = Command::new(bin_path())
        .arg("tail")
        .arg("-n")
        .arg("3")
        .arg(&log)
        .output()
        .expect("spawn binary");
    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines.len(), 3, "expected 3 lines, got {lines:?}");
    // Last 3 entries are test.action.2, .3, .4 (in chronological order).
    assert!(lines[0].contains("test.action.2"));
    assert!(lines[1].contains("test.action.3"));
    assert!(lines[2].contains("test.action.4"));
    // Each line includes the actor + the kind.
    assert!(lines[0].contains("alice"));
    assert!(lines[0].contains("case"));
    let _ = std::fs::remove_dir_all(&tmp);
}

#[test]
fn tail_default_count_is_fifty() {
    // 5 entries < 50 default -> all printed.
    let (tmp, log) = write_fixture_log(5);
    let out = Command::new(bin_path())
        .arg("tail")
        .arg(&log)
        .output()
        .expect("spawn binary");
    assert!(out.status.success());
    let stdout = String::from_utf8_lossy(&out.stdout);
    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines.len(), 5);
    let _ = std::fs::remove_dir_all(&tmp);
}

#[test]
fn tail_json_emits_parseable_array() {
    let (tmp, log) = write_fixture_log(3);
    let out = Command::new(bin_path())
        .arg("tail")
        .arg("--json")
        .arg(&log)
        .output()
        .expect("spawn binary");
    assert!(out.status.success());
    let stdout = String::from_utf8(out.stdout).expect("utf-8");
    let v: serde_json::Value = serde_json::from_str(&stdout).expect("parseable JSON");
    assert!(v.is_array(), "expected array, got: {stdout}");
    let arr = v.as_array().unwrap();
    assert_eq!(arr.len(), 3);
    assert_eq!(arr[0]["actor"]["id"], "alice");
    let _ = std::fs::remove_dir_all(&tmp);
}

#[test]
fn tail_since_drops_entries_before_cutoff() {
    // write_fixture_log produces timestamps `2026-04-28T00:00:NNZ`
    // for NN = 00, 01, 02, …. Cut at :02 to keep 02, 03, 04.
    let (tmp, log) = write_fixture_log(5);
    let out = Command::new(bin_path())
        .arg("tail")
        .arg("--since")
        .arg("2026-04-28T00:00:02Z")
        .arg(&log)
        .output()
        .expect("spawn binary");
    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    let lines: Vec<&str> = stdout.lines().collect();
    // Three entries on or after the cutoff: test.action.2, .3, .4.
    assert_eq!(lines.len(), 3, "got {lines:?}");
    assert!(lines[0].contains("test.action.2"));
    assert!(lines[2].contains("test.action.4"));
    // Earlier entries dropped.
    assert!(
        !stdout.contains("test.action.0"),
        "pre-cutoff entry leaked: {stdout}"
    );
    let _ = std::fs::remove_dir_all(&tmp);
}

#[test]
fn tail_since_combines_with_n_filters_first_truncates_second() {
    // 5 entries, cutoff at :02 (keeps 3), -n 1 → keep just the last
    // surviving entry (test.action.4).
    let (tmp, log) = write_fixture_log(5);
    let out = Command::new(bin_path())
        .arg("tail")
        .arg("-n")
        .arg("1")
        .arg("--since")
        .arg("2026-04-28T00:00:02Z")
        .arg(&log)
        .output()
        .expect("spawn binary");
    assert!(out.status.success());
    let stdout = String::from_utf8_lossy(&out.stdout);
    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines.len(), 1);
    assert!(lines[0].contains("test.action.4"), "got: {lines:?}");
    let _ = std::fs::remove_dir_all(&tmp);
}

#[test]
fn tail_missing_file_exits_one() {
    let out = Command::new(bin_path())
        .arg("tail")
        .arg("/tmp/this-does-not-exist-valenx-audit-cli.jsonl")
        .output()
        .expect("spawn binary");
    assert_eq!(out.status.code(), Some(1));
}

#[test]
fn tail_invalid_count_exits_two() {
    let (tmp, log) = write_fixture_log(1);
    let out = Command::new(bin_path())
        .arg("tail")
        .arg("-n")
        .arg("not-a-number")
        .arg(&log)
        .output()
        .expect("spawn binary");
    assert_eq!(out.status.code(), Some(2));
    let _ = std::fs::remove_dir_all(&tmp);
}
