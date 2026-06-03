//! Integration tests for the `valenx-results` binary. Build a
//! synthetic `Results` in code, persist it as `results.json` in a
//! temp dir, then spawn the compiled CLI against it and assert the
//! exit code + stdout shape.

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
    PathBuf::from(env!("CARGO_BIN_EXE_valenx-results"))
}

fn tempdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!(
        "valenx-results-cli-{tag}-{}",
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
        wall_time_seconds: 1.5,
        completed_at: "2026-04-28T00:00:00Z".into(),
        ancestors: Vec::new(),
    }
}

/// Build a results.json with a known shape and write it to disk.
fn write_synthetic_results(path: &std::path::Path) {
    let mut r = Results::empty("smoke", synthetic_provenance());
    r.scalars.insert(ScalarRecord {
        name: "T_final".into(),
        value: 298.15,
        units: DIMENSIONLESS,
        time: TimeKey::Steady,
        source: valenx_fields::scalar::ScalarSource::Extracted,
        description: None,
    });
    r.scalars.insert(ScalarRecord {
        name: "X_N2".into(),
        value: 0.78,
        units: DIMENSIONLESS,
        time: TimeKey::Steady,
        source: valenx_fields::scalar::ScalarSource::Extracted,
        description: None,
    });
    r.artifacts.push(Artifact {
        path: PathBuf::from("flow.vtu"),
        kind: ArtifactKind::VizData,
        checksum: None,
        label: "synthetic VTU".into(),
    });
    let s = serde_json::to_string_pretty(&r).unwrap();
    std::fs::write(path, s).unwrap();
}

#[test]
fn text_mode_lists_scalars_and_artifacts() {
    let d = tempdir("text");
    let p = d.join("results.json");
    write_synthetic_results(&p);

    let out = Command::new(binary())
        .arg(&p)
        .output()
        .expect("spawn results");
    assert!(
        out.status.success(),
        "results failed; stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8(out.stdout).expect("utf-8");
    assert!(stdout.contains("smoke"), "got: {stdout}");
    assert!(stdout.contains("test-tool"), "got: {stdout}");
    assert!(stdout.contains("T_final"), "got: {stdout}");
    assert!(stdout.contains("X_N2"), "got: {stdout}");
    assert!(stdout.contains("flow.vtu"), "got: {stdout}");
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn json_mode_pretty_prints_full_envelope() {
    let d = tempdir("json");
    let p = d.join("results.json");
    write_synthetic_results(&p);

    let out = Command::new(binary())
        .arg(&p)
        .arg("--format")
        .arg("json")
        .output()
        .expect("spawn results");
    assert!(out.status.success(), "results failed");
    let stdout = String::from_utf8(out.stdout).expect("utf-8");
    let v: serde_json::Value = serde_json::from_str(&stdout).expect("stdout JSON");
    // Round-trip — the JSON envelope must contain the same case_id
    // and provenance tool name.
    assert_eq!(v["meta"]["case_id"], "smoke");
    assert_eq!(v["provenance"]["tool"], "test-tool");
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn reads_results_from_stdin_when_path_is_dash() {
    // Synthesise a Results, pipe it as JSON via stdin, and assert
    // the binary prints the same case + scalar names. Lets users
    // chain `kubectl get …` / `aws s3 cp s3://… -` style sources
    // into the inspector without staging to a temp file.
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
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .expect("spawn results");
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
    let stdout = String::from_utf8(out.stdout).expect("utf-8");
    assert!(stdout.contains("smoke"), "got: {stdout}");
    assert!(stdout.contains("T_final"), "got: {stdout}");
    // Display path swaps to <stdin> rather than the literal `-`
    // since the latter is just a sentinel.
    assert!(
        stdout.contains("<stdin>"),
        "expected <stdin> marker; got: {stdout}"
    );
}

#[test]
fn exit_three_on_missing_file() {
    let out = Command::new(binary())
        .arg("/does/not/exist/hopefully-/.results.json")
        .output()
        .expect("spawn results");
    let code = out.status.code().expect("native exit code");
    assert_eq!(code, 3, "expected exit 3 (IO); got {code}");
}

#[test]
fn exit_one_on_malformed_json() {
    let d = tempdir("bad");
    let p = d.join("results.json");
    std::fs::write(&p, b"{ this is not JSON").unwrap();
    let out = Command::new(binary())
        .arg(&p)
        .output()
        .expect("spawn results");
    let code = out.status.code().expect("native exit code");
    let _ = std::fs::remove_dir_all(&d);
    assert_eq!(code, 1, "expected exit 1 (parse); got {code}");
}

#[test]
fn exit_two_on_usage_error() {
    let out = Command::new(binary())
        .arg("--bogus")
        .output()
        .expect("spawn results");
    let code = out.status.code().expect("native exit code");
    assert_eq!(code, 2, "expected exit 2 (usage); got {code}");
}
