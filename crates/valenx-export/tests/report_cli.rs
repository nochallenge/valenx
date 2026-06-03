//! Integration tests for the `valenx-report` binary. Build a
//! synthetic Results in memory, persist it, then spawn the binary
//! and assert HTML / CSV land on disk with a reasonable shape.

use std::path::PathBuf;
use std::process::Command;

use valenx_fields::{
    artifact::{Artifact, ArtifactKind},
    provenance::{Provenance, Sha256Hex},
    scalar::ScalarRecord,
    units::DIMENSIONLESS,
    Results, TimeKey,
};

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_valenx-report"))
}

fn tempdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!(
        "valenx-report-cli-{tag}-{}",
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn synthetic_provenance() -> Provenance {
    Provenance {
        adapter: "test".into(),
        adapter_version: "0.0.0".into(),
        tool: "test-tool".into(),
        tool_version: "1.2.3".into(),
        case_hash: Sha256Hex(String::new()),
        mesh_hash: Sha256Hex(String::new()),
        input_hash: Sha256Hex(String::new()),
        tools_lock_hash: Sha256Hex(String::new()),
        run_id: "00000000-0000-0000-0000-000000000000".into(),
        wall_time_seconds: 2.5,
        completed_at: "2026-04-28T00:00:00Z".into(),
        ancestors: Vec::new(),
    }
}

fn write_results(path: &std::path::Path) {
    let mut r = Results::empty("smoke", synthetic_provenance());
    r.scalars.insert(ScalarRecord {
        name: "T_final".into(),
        value: 298.15,
        units: DIMENSIONLESS,
        time: TimeKey::Steady,
        source: valenx_fields::scalar::ScalarSource::Extracted,
        description: None,
    });
    r.artifacts.push(Artifact {
        path: PathBuf::from("flow.vtu"),
        kind: ArtifactKind::VizData,
        checksum: None,
        label: "synthetic".into(),
    });
    let s = serde_json::to_string_pretty(&r).unwrap();
    std::fs::write(path, s).unwrap();
}

#[test]
fn writes_html_when_asked() {
    let d = tempdir("html");
    let input = d.join("results.json");
    let html = d.join("report.html");
    write_results(&input);

    let out = Command::new(binary())
        .arg(&input)
        .arg("--html")
        .arg(&html)
        .output()
        .expect("spawn report");
    assert!(
        out.status.success(),
        "report failed; stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(html.is_file(), "html file should have been created");
    let body = std::fs::read_to_string(&html).unwrap();
    // Sanity-check the HTML contains the case + a scalar name.
    assert!(body.contains("smoke"), "html missing case id; got: {body}");
    assert!(
        body.contains("T_final"),
        "html missing scalar name; got: {body}"
    );
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn writes_csv_when_asked() {
    let d = tempdir("csv");
    let input = d.join("results.json");
    let csv = d.join("scalars.csv");
    write_results(&input);

    let out = Command::new(binary())
        .arg(&input)
        .arg("--csv")
        .arg(&csv)
        .output()
        .expect("spawn report");
    assert!(out.status.success(), "report failed");
    assert!(csv.is_file(), "csv file should have been created");
    let body = std::fs::read_to_string(&csv).unwrap();
    assert!(
        body.contains("T_final"),
        "csv missing scalar name; got: {body}"
    );
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn writes_both_html_and_csv() {
    let d = tempdir("both");
    let input = d.join("results.json");
    let html = d.join("r.html");
    let csv = d.join("r.csv");
    write_results(&input);

    let out = Command::new(binary())
        .arg(&input)
        .arg("--html")
        .arg(&html)
        .arg("--csv")
        .arg(&csv)
        .output()
        .expect("spawn report");
    assert!(out.status.success(), "report failed");
    assert!(html.is_file());
    assert!(csv.is_file());
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn writes_markdown_when_asked() {
    let d = tempdir("md");
    let input = d.join("results.json");
    let md = d.join("report.md");
    write_results(&input);

    let out = Command::new(binary())
        .arg(&input)
        .arg("--markdown")
        .arg(&md)
        .output()
        .expect("spawn report");
    assert!(
        out.status.success(),
        "report failed; stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(md.is_file(), "markdown file should have been created");
    let body = std::fs::read_to_string(&md).unwrap();
    // Markdown structural sanity-check: header + a section heading +
    // pipe-table markers + the synthetic scalar name.
    assert!(body.starts_with("# Valenx run report"), "got: {body}");
    assert!(body.contains("## Provenance"), "got: {body}");
    assert!(body.contains("## Scalars"), "got: {body}");
    assert!(body.contains("`T_final`"), "missing scalar: {body}");
    assert!(body.contains("|---|"), "missing table separator: {body}");
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn writes_all_three_formats_in_one_run() {
    let d = tempdir("all3");
    let input = d.join("results.json");
    let html = d.join("r.html");
    let md = d.join("r.md");
    let csv = d.join("r.csv");
    write_results(&input);

    let out = Command::new(binary())
        .arg(&input)
        .arg("--html")
        .arg(&html)
        .arg("--markdown")
        .arg(&md)
        .arg("--csv")
        .arg(&csv)
        .output()
        .expect("spawn report");
    assert!(out.status.success(), "report failed");
    assert!(html.is_file());
    assert!(md.is_file());
    assert!(csv.is_file());
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn reads_results_from_stdin_when_input_is_dash() {
    let d = tempdir("stdin");
    let html = d.join("r.html");

    // Build a synthetic Results in memory and pipe it via stdin.
    let mut r = Results::empty("smoke", synthetic_provenance());
    r.scalars.insert(ScalarRecord {
        name: "T_final".into(),
        value: 298.15,
        units: DIMENSIONLESS,
        time: TimeKey::Steady,
        source: valenx_fields::scalar::ScalarSource::Extracted,
        description: None,
    });
    let json = serde_json::to_string_pretty(&r).unwrap();

    let mut child = Command::new(binary())
        .arg("-")
        .arg("--html")
        .arg(&html)
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .expect("spawn report");
    {
        use std::io::Write;
        let stdin = child.stdin.as_mut().expect("stdin");
        stdin.write_all(json.as_bytes()).expect("write stdin");
    }
    let out = child.wait_with_output().expect("wait");
    assert!(
        out.status.success(),
        "stdin path failed; stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(html.is_file(), "html file should have been created");
    let body = std::fs::read_to_string(&html).unwrap();
    assert!(body.contains("smoke"), "html missing case id; got: {body}");
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn exits_two_when_no_destination_given() {
    let d = tempdir("nodst");
    let input = d.join("results.json");
    write_results(&input);

    let out = Command::new(binary())
        .arg(&input)
        .output()
        .expect("spawn report");
    let code = out.status.code().expect("exit code");
    let _ = std::fs::remove_dir_all(&d);
    assert_eq!(
        code, 2,
        "expected exit 2 (usage — needs --html or --csv); got {code}"
    );
}

#[test]
fn exits_one_on_malformed_json() {
    let d = tempdir("bad");
    let input = d.join("results.json");
    let html = d.join("r.html");
    std::fs::write(&input, b"{ this is not JSON").unwrap();

    let out = Command::new(binary())
        .arg(&input)
        .arg("--html")
        .arg(&html)
        .output()
        .expect("spawn report");
    let code = out.status.code().expect("exit code");
    let _ = std::fs::remove_dir_all(&d);
    assert_eq!(code, 1, "expected exit 1 (parse); got {code}");
}

#[test]
fn exits_three_when_input_missing() {
    let d = tempdir("missing-in");
    let html = d.join("r.html");
    let out = Command::new(binary())
        .arg(d.join("does-not-exist.json"))
        .arg("--html")
        .arg(&html)
        .output()
        .expect("spawn report");
    let code = out.status.code().expect("exit code");
    let _ = std::fs::remove_dir_all(&d);
    assert_eq!(code, 3, "expected exit 3 (IO); got {code}");
}
