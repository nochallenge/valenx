use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Stdio};

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_valenx-fastq"))
}

const FQ: &str = "\
@r1
ACGTACGT
+
IIIIIIII
@r2
NNNN
+
!!!!
";

#[test]
fn inspect_text_mode_prints_record_count() {
    let mut child = Command::new(bin())
        .arg("inspect")
        .arg("-")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("spawn");
    child
        .stdin
        .as_mut()
        .unwrap()
        .write_all(FQ.as_bytes())
        .unwrap();
    let out = child.wait_with_output().unwrap();
    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(stdout.contains("records: 2"), "got: {stdout}");
    assert!(
        stdout.contains("total bases: 12") || stdout.contains("total_bases: 12"),
        "got: {stdout}"
    );
}

#[test]
fn inspect_json_mode_emits_record_count() {
    let mut child = Command::new(bin())
        .arg("inspect")
        .arg("--format")
        .arg("json")
        .arg("-")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("spawn");
    child
        .stdin
        .as_mut()
        .unwrap()
        .write_all(FQ.as_bytes())
        .unwrap();
    let out = child.wait_with_output().unwrap();
    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let v: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(v["records"], 2);
    assert_eq!(v["total_bases"], 12);
}
