//! Integration tests for the `valenx-pdb-info` binary.
//!
//! Build the binary via cargo's test fixture (CARGO_BIN_EXE_*),
//! spawn it as a subprocess, and verify text + JSON output shapes
//! against the bundled `1ubq-tiny.pdb` fixture (the same one the
//! pdb_round_trip integration test consumes).

use std::path::PathBuf;
use std::process::Command;

/// Path to the compiled `valenx-pdb-info` binary, set by cargo at
/// test build time.
fn bin_path() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_valenx-pdb-info"))
}

/// Path to the bundled `tests/fixtures/biology/1ubq-tiny.pdb`
/// fixture (workspace-rooted, shared with `pdb_round_trip`).
fn fixture() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .ancestors()
        .nth(2)
        .expect("workspace root")
        .join("tests")
        .join("fixtures")
        .join("biology")
        .join("1ubq-tiny.pdb")
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
    assert!(stdout.contains("valenx-pdb-info"));
    assert!(stdout.contains("USAGE"));
    assert!(stdout.contains("--format"));
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
    assert!(stdout.contains("valenx-pdb-info"));
    assert!(stdout.contains('v'), "expected version marker: {stdout}");
}

#[test]
fn text_format_prints_chain_summary() {
    let out = Command::new(bin_path())
        .arg(fixture())
        .output()
        .expect("spawn binary");
    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    // Header line: id (1ubq-tiny), chain count (1), residues (2), atoms (9).
    assert!(stdout.contains("1ubq-tiny"), "id missing: {stdout}");
    assert!(stdout.contains("1 chain"), "chain count missing: {stdout}");
    assert!(
        stdout.contains("2 residues"),
        "residue count missing: {stdout}"
    );
    assert!(stdout.contains("9 atoms"), "atom count missing: {stdout}");
    // Per-chain row.
    assert!(stdout.contains("chain A"), "chain A row missing: {stdout}");
    // Element tally — fixture has 4 C, 2 N, 2 O, 1 (the second MET CB
    // is a C). Recount: MET has N, CA, C, O, CB — 1N 3C 1O. GLN has
    // N, CA, C, O — 1N 2C 1O. Totals: 5C 2N 2O.
    assert!(
        stdout.contains("elements:"),
        "elements line missing: {stdout}"
    );
    assert!(
        stdout.contains(" C") || stdout.contains("C,"),
        "carbon entry missing: {stdout}"
    );
    // Residue range row.
    assert!(
        stdout.contains("residue range: 1..=2"),
        "residue range missing: {stdout}"
    );
}

#[test]
fn json_format_emits_parseable_object() {
    let out = Command::new(bin_path())
        .arg(fixture())
        .arg("--format")
        .arg("json")
        .output()
        .expect("spawn binary");
    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8(out.stdout).expect("UTF-8 output");
    let v: serde_json::Value = serde_json::from_str(&stdout).expect("stdout is parseable JSON");
    assert_eq!(v["id"], "1ubq-tiny");
    assert_eq!(v["atom_count"], 9);
    assert_eq!(v["residue_count"], 2);
    let chains = v["chains"].as_array().expect("chains is array");
    assert_eq!(chains.len(), 1);
    assert_eq!(chains[0]["id"], "A");
    assert_eq!(chains[0]["residues"], 2);
    assert_eq!(chains[0]["atoms"], 9);
}

#[test]
fn reads_pdb_from_stdin_when_path_is_dash() {
    // Pipe the fixture through stdin via `-`. Mirrors the
    // valenx-mesh-info stdin-mode test.
    let pdb_text = std::fs::read(fixture()).expect("read fixture");

    let mut child = Command::new(bin_path())
        .arg("-")
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .expect("spawn binary");
    {
        use std::io::Write;
        let stdin = child.stdin.as_mut().expect("stdin");
        stdin.write_all(&pdb_text).expect("write stdin");
    }
    let out = child.wait_with_output().expect("wait");
    assert!(
        out.status.success(),
        "stdin mode failed; stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    // When reading from stdin the id falls back to "stdin".
    assert!(stdout.contains("stdin"), "stdin id missing: {stdout}");
    assert!(stdout.contains("9 atoms"), "atom count missing: {stdout}");
}
