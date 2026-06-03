//! Integration tests for the `valenx-fasta` binary.
//!
//! Build the binary via cargo's test fixture (CARGO_BIN_EXE_*),
//! spawn it as a subprocess, and verify each subcommand's output.
//! End-to-end coverage that the unit tests inside the binary file
//! can't reach (they exercise the parser + helpers, not the actual
//! exe + stdout pipeline).

use std::path::PathBuf;
use std::process::Command;

/// Path to the compiled `valenx-fasta` binary, set by cargo at test
/// build time.
fn bin_path() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_valenx-fasta"))
}

/// Path to the bundled `tests/fixtures/biology/sample.fasta` fixture.
/// Lives at the workspace root so the same file serves the
/// integration tests AND example users running
/// `valenx-fasta inspect tests/fixtures/biology/sample.fasta`.
fn fixture() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .ancestors()
        .nth(2)
        .expect("workspace root")
        .join("tests")
        .join("fixtures")
        .join("biology")
        .join("sample.fasta")
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
    assert!(stdout.contains("valenx-fasta"));
    assert!(stdout.contains("USAGE"));
    assert!(stdout.contains("inspect"));
    assert!(stdout.contains("validate"));
    assert!(stdout.contains("extract"));
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
    assert!(stdout.contains("valenx-fasta"));
    // Format is `<name> v<CARGO_PKG_VERSION>`.
    assert!(stdout.contains('v'), "expected version marker: {stdout}");
}

#[test]
fn inspect_text_output_lists_records() {
    let out = Command::new(bin_path())
        .arg("inspect")
        .arg(fixture())
        .output()
        .expect("spawn binary");
    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    // Header row should mention "2 records" — the fixture has 2.
    assert!(
        stdout.contains("2 records"),
        "expected `2 records` in output: {stdout}"
    );
    // Both record names should appear.
    assert!(stdout.contains("p53"), "p53 missing: {stdout}");
    assert!(stdout.contains("ubq"), "ubq missing: {stdout}");
    // Both alphabets should be detected as protein.
    assert!(
        stdout.contains("protein"),
        "protein label missing: {stdout}"
    );
    // 76 residues (ubiquitin) appears in output.
    assert!(stdout.contains("76"), "expected 76 in output: {stdout}");
}

#[test]
fn inspect_json_output_parses_as_jsonarray() {
    let out = Command::new(bin_path())
        .arg("inspect")
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
    assert!(v["path"].is_string());
    let records = v["records"].as_array().expect("records is array");
    assert_eq!(records.len(), 2);
    assert_eq!(records[0]["alphabet"], "protein");
    assert_eq!(records[1]["length"], 76);
}

#[test]
fn validate_rejects_protein_input_with_dna_alphabet() {
    // The bundled fixture is protein. Forcing alphabet=dna should fail
    // because `M` (or rather the first non-DNA byte; `E` for glutamate)
    // is invalid in DNA. Exit code 1 (content error).
    let out = Command::new(bin_path())
        .arg("validate")
        .arg(fixture())
        .arg("--alphabet")
        .arg("dna")
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
        stderr.contains("invalid byte") || stderr.contains("position"),
        "stderr should mention the offending byte: {stderr}"
    );
    assert!(
        stderr.contains("dna"),
        "stderr should mention dna: {stderr}"
    );
}

#[test]
fn extract_named_record_prints_fasta_shape() {
    let out = Command::new(bin_path())
        .arg("extract")
        .arg(fixture())
        .arg("--name")
        .arg("ubq|ubiquitin")
        .output()
        .expect("spawn binary");
    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    // First line is `>name`, second line is the body.
    let lines: Vec<&str> = stdout.lines().collect();
    assert!(lines.len() >= 2, "expected at least 2 lines: {stdout}");
    assert_eq!(lines[0], ">ubq|ubiquitin");
    // Body must contain the expected residues from the fixture's
    // ubiquitin sequence (joined to one line, 76 residues).
    assert!(lines[1].starts_with("MQIFVKTLTG"), "body: {}", lines[1]);
    assert_eq!(lines[1].len(), 76, "body length should be 76");
}
