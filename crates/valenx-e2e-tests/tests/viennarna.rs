//! E2E integration test for the ViennaRNA (`RNAfold`) adapter.
//!
//! Spawns a real `RNAfold` subprocess via the adapter's
//! prepare/run/collect pipeline. Skipped automatically when the
//! upstream binary isn't on PATH.
//!
//! RNAfold reads a FASTA-style input and writes the minimum-free-energy
//! structure to stdout. The adapter redirects stdout to the
//! user-specified output filename in the workdir; `collect()` should
//! surface that as the canonical structure artifact.

use valenx_adapter_viennarna::ViennaRnaAdapter;
use valenx_core::{Adapter, CancellationToken, Case, LogLevel, LogSink, ProgressSink, RunContext};
use valenx_test_utils::tempdir;

struct NoopProgress;
impl ProgressSink for NoopProgress {
    fn report(&self, _pct: f32, _message: &str) {}
}
struct NoopLog;
impl LogSink for NoopLog {
    fn log_line(&self, _level: LogLevel, _line: &str) {}
}

fn skip_if_missing(adapter: &dyn Adapter) -> bool {
    match adapter.probe() {
        Ok(report) if report.ok => false,
        _ => {
            eprintln!(
                "Skipping E2E test — `{}` upstream binary not installed on PATH.",
                adapter.info().id
            );
            true
        }
    }
}

/// A short RNA sequence that folds into a clean stem-loop (hairpin).
/// 30 nt, easy for RNAfold to process in milliseconds. The sequence
/// has a self-complementary 5'/3' region so the MFE structure is a
/// short hairpin.
const TINY_RNA_FASTA: &str = "\
>hairpin
GGGAAACCCAAAGGGAAACCCAAAGGGAAACCC
";

const MINIMAL_CASE_TOML: &str = r#"[case]
physics = "bio"
solver  = "viennarna.fold"

[bio.viennarna]
input  = "rna.fa"
output = "fold.out"
"#;

#[test]
fn viennarna_runs_end_to_end_against_real_binary() {
    let adapter = ViennaRnaAdapter::new();
    if skip_if_missing(&adapter) {
        return;
    }

    let case_dir = tempdir("e2e-viennarna");
    std::fs::write(case_dir.join("case.toml"), MINIMAL_CASE_TOML).unwrap();
    std::fs::write(case_dir.join("rna.fa"), TINY_RNA_FASTA).unwrap();

    let workdir = case_dir.join(".valenx-workdir");
    std::fs::create_dir_all(&workdir).unwrap();

    let case = Case {
        id: "e2e-viennarna".into(),
        path: case_dir.clone(),
    };

    let prepared = adapter.prepare(&case, &workdir).expect("prepare");

    let cancel = CancellationToken::new();
    let mut ctx = RunContext {
        cancel: &cancel,
        progress: Box::new(NoopProgress),
        log: Box::new(NoopLog),
    };
    let report = adapter.run(&prepared, &mut ctx).expect("run");
    assert_eq!(report.exit_code, 0, "RNAfold exited non-zero: {report:?}");

    let results = adapter.collect(&prepared).expect("collect");
    assert!(
        !results.artifacts.is_empty(),
        "collect found no artifacts in {}",
        workdir.display()
    );
    // The output file the user named in case.toml.
    assert!(
        workdir.join("fold.out").is_file(),
        "expected fold.out in workdir, found artifacts {:?}",
        results.artifacts
    );

    let _ = std::fs::remove_dir_all(&case_dir);
}
