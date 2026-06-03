//! Integration tests for the `valenx-validate` binary. Spawn the
//! compiled CLI against the bundled `tests/fixtures/minimal.valenx`
//! fixture and assert exit code + stdout shape.
//!
//! These tests use `CARGO_BIN_EXE_valenx-validate`, which Cargo
//! populates for every `[[bin]]` target in the same package. The
//! sibling `validate_cli.rs` tests inside the binary cover argument
//! parsing; this file covers actually running the binary.

use std::path::PathBuf;
use std::process::Command;

/// Locate the validate binary via the env var Cargo sets at test
/// time. Fails the build if the binary doesn't exist (it should,
/// because Cargo always builds bins before integration tests).
fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_valenx-validate"))
}

/// Absolute path to `tests/fixtures/minimal.valenx`. Same shape as
/// `project_roundtrip::fixture_path`; duplicated here to keep the
/// two files independent.
fn fixture_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .expect("workspace root")
        .join("tests")
        .join("fixtures")
        .join("minimal.valenx")
}

#[test]
fn exit_zero_on_clean_fixture_with_text_output() {
    let out = Command::new(binary())
        .arg(fixture_path())
        .output()
        .expect("spawn validate");
    assert!(
        out.status.success(),
        "validate failed; stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8(out.stdout).expect("utf-8");
    // Header line names the project; one line per case.
    assert!(stdout.contains("project `minimal`"), "got: {stdout}");
    // Every case in the fixture's order list should appear in the
    // text output.
    for name in [
        "box-mesh",
        "cfd-steady",
        "cfd-transient",
        "fea-cantilever",
        "heat-cube",
        "netgen-cylinder",
    ] {
        assert!(
            stdout.contains(name),
            "stdout missing case `{name}`:\n{stdout}"
        );
    }
}

#[test]
fn exit_zero_on_clean_fixture_with_json_output() {
    let out = Command::new(binary())
        .arg(fixture_path())
        .arg("--format")
        .arg("json")
        .output()
        .expect("spawn validate");
    assert!(out.status.success(), "validate failed");
    let stdout = String::from_utf8(out.stdout).expect("utf-8");
    let v: serde_json::Value = serde_json::from_str(&stdout)
        .unwrap_or_else(|e| panic!("stdout not JSON: {e}\n---\n{stdout}"));
    assert_eq!(v["ok"], serde_json::Value::Bool(true));
    assert_eq!(v["project"]["name"], "minimal");
    let cases = v["cases"].as_array().expect("cases array");
    assert_eq!(cases.len(), 6, "expected 6 cases; got {}", cases.len());
}

#[test]
fn exit_three_when_path_does_not_exist() {
    let out = Command::new(binary())
        .arg("/does/not/exist/hopefully-/.valenx-not-a-real-path")
        .output()
        .expect("spawn validate");
    let code = out.status.code().expect("native exit code");
    assert_eq!(code, 3, "expected exit 3 (IO); got {code}");
}

#[test]
fn exit_one_when_project_toml_is_missing() {
    let tmp = std::env::temp_dir().join(format!(
        "valenx-validate-empty-{}",
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir_all(&tmp).unwrap();
    let out = Command::new(binary())
        .arg(&tmp)
        .output()
        .expect("spawn validate");
    let code = out.status.code().expect("native exit code");
    let _ = std::fs::remove_dir_all(&tmp);
    assert_eq!(
        code, 1,
        "expected exit 1 (structural — missing manifest); got {code}"
    );
}

#[test]
fn exit_two_on_usage_error() {
    let out = Command::new(binary())
        .arg("--bogus")
        .output()
        .expect("spawn validate");
    let code = out.status.code().expect("native exit code");
    assert_eq!(code, 2, "expected exit 2 (usage); got {code}");
}

#[test]
fn json_error_envelope_when_load_fails() {
    let out = Command::new(binary())
        .arg("/does/not/exist/hopefully-/.valenx-also-not-real")
        .arg("--format")
        .arg("json")
        .output()
        .expect("spawn validate");
    let stdout = String::from_utf8(out.stdout).expect("utf-8");
    let v: serde_json::Value = serde_json::from_str(&stdout)
        .unwrap_or_else(|e| panic!("stdout not JSON: {e}\n---\n{stdout}"));
    assert_eq!(v["ok"], serde_json::Value::Bool(false));
    // Error message echoes the path so users can see what went wrong.
    assert!(v["error"].as_str().is_some());
}
