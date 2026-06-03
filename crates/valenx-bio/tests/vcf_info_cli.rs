use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Stdio};

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_valenx-vcf-info"))
}

const VCF: &str = "\
##fileformat=VCFv4.2
##source=test
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\t100\trs1\tA\tG\t40.0\tPASS\tDP=20\tGT:DP\t0/1:20
chr1\t200\t.\tA\tT\t.\t.\tDP=10\tGT:DP\t0/1:10
chr1\t300\t.\tA\t.\t.\tLowQual\tDP=5\tGT:DP\t0/0:5
";

#[test]
fn text_mode_lists_pass_fail_split() {
    let mut child = Command::new(bin())
        .arg("-")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("spawn");
    child
        .stdin
        .as_mut()
        .unwrap()
        .write_all(VCF.as_bytes())
        .unwrap();
    let out = child.wait_with_output().unwrap();
    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    // 3 records: rs1 (PASS), .200 (empty filter = pass), .300 (LowQual = fail)
    assert!(stdout.contains("records: 3"), "got: {stdout}");
    assert!(stdout.contains("pass: 2"), "got: {stdout}");
    assert!(stdout.contains("fail: 1"), "got: {stdout}");
    // 1 record (.300) has no ALT (.).
    assert!(stdout.contains("no-alt: 1"), "got: {stdout}");
    assert!(stdout.contains("samples: 1"));
    assert!(stdout.contains("header lines: 2"));
}

#[test]
fn json_mode_round_trips() {
    let mut child = Command::new(bin())
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
        .write_all(VCF.as_bytes())
        .unwrap();
    let out = child.wait_with_output().unwrap();
    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let v: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(v["records"], 3);
    assert_eq!(v["pass"], 2);
    assert_eq!(v["fail"], 1);
    assert_eq!(v["no_alt"], 1);
    assert_eq!(v["samples"], 1);
    assert_eq!(v["header_lines"], 2);
}
