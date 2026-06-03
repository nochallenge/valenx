//! Integration tests for the `valenx-mesh-info` binary.
//!
//! Build the binary via cargo's test fixture (CARGO_BIN_EXE_*),
//! spawn it as a subprocess, and verify both text + JSON output
//! shapes against a known mesh. End-to-end coverage that the unit
//! tests inside the binary file can't reach (they exercise the
//! library functions, not the actual exe + stdout pipeline).

use std::path::PathBuf;
use std::process::Command;

use nalgebra::Vector3;
use valenx_mesh::{ElementBlock, ElementType, Mesh};

/// Path to the compiled `valenx-mesh-info` binary, set by cargo at
/// test build time.
fn bin_path() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_valenx-mesh-info"))
}

/// Write a one-element right-isoceles Tri3 mesh to a tempdir as
/// canonical JSON. Returns (workdir, mesh_path) — caller cleans up.
fn write_fixture_mesh() -> (PathBuf, PathBuf) {
    let tmp = std::env::temp_dir().join(format!(
        "valenx-mesh-info-cli-{}",
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir_all(&tmp).unwrap();

    let mut m = Mesh::new("right-iso");
    m.nodes = vec![
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
    ];
    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = vec![0, 1, 2];
    m.element_blocks.push(block);

    let mesh_path = tmp.join("mesh.json");
    std::fs::write(&mesh_path, serde_json::to_string(&m).unwrap()).unwrap();
    (tmp, mesh_path)
}

#[test]
fn help_prints_usage_and_exits_zero() {
    let out = Command::new(bin_path())
        .arg("help")
        .output()
        .expect("spawn binary");
    assert!(
        out.status.success(),
        "expected exit 0, got {:?}",
        out.status
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(stdout.contains("valenx-mesh-info"));
    assert!(stdout.contains("USAGE"));
    assert!(stdout.contains("--format"));
}

#[test]
fn no_args_prints_usage_to_stderr_and_exits_two() {
    let out = Command::new(bin_path()).output().expect("spawn binary");
    assert_eq!(
        out.status.code(),
        Some(2),
        "expected exit 2, got {:?}",
        out.status
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("missing mesh file"));
    assert!(stderr.contains("USAGE"));
}

#[test]
fn missing_file_exits_one_with_io_error() {
    let out = Command::new(bin_path())
        .arg("/tmp/this-path-does-not-exist-valenx-cli-test.json")
        .output()
        .expect("spawn binary");
    assert_eq!(out.status.code(), Some(1));
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("error"));
}

#[test]
fn text_format_prints_quality_report_lines() {
    let (tmp, mesh_path) = write_fixture_mesh();
    let out = Command::new(bin_path())
        .arg(&mesh_path)
        .output()
        .expect("spawn binary");
    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(stdout.contains("Mesh quality report"));
    assert!(stdout.contains("elements: 1"));
    assert!(stdout.contains("max aspect"));
    assert!(stdout.contains("max skew"));
    // Right-isoceles triangle has skewness 0.25, aspect sqrt(2).
    assert!(
        stdout.contains("0.250"),
        "expected 0.250 in output: {stdout}"
    );
    assert!(stdout.contains("Aspect-ratio histogram"));
    assert!(stdout.contains("Skewness histogram"));
    let _ = std::fs::remove_dir_all(&tmp);
}

#[test]
fn json_format_emits_parseable_object() {
    let (tmp, mesh_path) = write_fixture_mesh();
    let out = Command::new(bin_path())
        .arg(&mesh_path)
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
    assert_eq!(v["quality"]["element_count"], 1);
    assert_eq!(v["quality"]["inverted_count"], 0);
    // Skewness 0.25 lands in bucket 0 (≤ 0.25).
    assert_eq!(v["skewness_histogram"]["counts"][0], 1);
    assert!(v["aspect_histogram"]["buckets"].is_array());
    assert!(v["path"].is_string());
    let _ = std::fs::remove_dir_all(&tmp);
}

#[test]
fn unknown_format_exits_two() {
    let (tmp, mesh_path) = write_fixture_mesh();
    let out = Command::new(bin_path())
        .arg(&mesh_path)
        .arg("--format")
        .arg("xml")
        .output()
        .expect("spawn binary");
    assert_eq!(out.status.code(), Some(2));
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("xml"), "stderr: {stderr}");
    let _ = std::fs::remove_dir_all(&tmp);
}

#[test]
fn check_passes_when_thresholds_met_exits_zero() {
    // Right-isoceles fixture: skew=0.25, aspect=sqrt(2), inverted=0.
    // Set permissive thresholds — all should pass.
    let (tmp, mesh_path) = write_fixture_mesh();
    let out = Command::new(bin_path())
        .arg(&mesh_path)
        .arg("--check")
        .arg("max-skew=0.5")
        .arg("--check")
        .arg("max-aspect=10")
        .arg("--check")
        .arg("inverted=0")
        .output()
        .expect("spawn binary");
    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let _ = std::fs::remove_dir_all(&tmp);
}

#[test]
fn check_fails_when_threshold_exceeded_exits_four() {
    // Skew = 0.25 against threshold 0.1 -> failure.
    let (tmp, mesh_path) = write_fixture_mesh();
    let out = Command::new(bin_path())
        .arg(&mesh_path)
        .arg("--check")
        .arg("max-skew=0.1")
        .output()
        .expect("spawn binary");
    assert_eq!(
        out.status.code(),
        Some(4),
        "expected exit 4, got {:?} (stderr: {})",
        out.status,
        String::from_utf8_lossy(&out.stderr)
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("max-skew"), "stderr: {stderr}");
    assert!(stderr.contains("0.1"), "stderr: {stderr}");
    // Stdout should still contain the report — checks run AFTER
    // printing so users see the numbers either way.
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(stdout.contains("max skew"), "stdout: {stdout}");
    let _ = std::fs::remove_dir_all(&tmp);
}

#[test]
fn check_unknown_metric_exits_two() {
    let (tmp, mesh_path) = write_fixture_mesh();
    let out = Command::new(bin_path())
        .arg(&mesh_path)
        .arg("--check")
        .arg("fancy-metric=42")
        .output()
        .expect("spawn binary");
    assert_eq!(out.status.code(), Some(2));
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("fancy-metric"));
    let _ = std::fs::remove_dir_all(&tmp);
}

#[test]
fn reads_mesh_from_stdin_when_path_is_dash() {
    // Build the same right-isoceles fixture mesh in memory and pipe
    // it to the binary via stdin. End-to-end coverage that the new
    // `-` directive correctly substitutes stdin for the file path.
    let mut m = Mesh::new("right-iso-stdin");
    m.nodes = vec![
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
    ];
    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = vec![0, 1, 2];
    m.element_blocks.push(block);
    let json = serde_json::to_string(&m).unwrap();

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
        stdin.write_all(json.as_bytes()).expect("write stdin");
    }
    let out = child.wait_with_output().expect("wait");
    assert!(
        out.status.success(),
        "stdin path failed; stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    // Display path swaps to <stdin>; mesh stats appear regardless.
    assert!(
        stdout.contains("<stdin>"),
        "expected <stdin> marker; got: {stdout}"
    );
    assert!(
        stdout.contains("elements:"),
        "expected `elements:` line in stats; got: {stdout}"
    );
}
