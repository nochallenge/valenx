//! E2E integration test for the Clustal Omega adapter.
//!
//! Spawns a real `clustalo` subprocess via the adapter's
//! prepare/run/collect pipeline. Skipped automatically when the
//! upstream binary isn't on PATH.

use valenx_adapter_clustalo::ClustaloAdapter;
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

/// Three short, similar protein sequences for a fast progressive
/// alignment.
const TINY_FASTA: &str = "\
>seq1
ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY
>seq2
ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY
>seq3
ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY
";

const MINIMAL_CASE_TOML: &str = r#"[case]
physics = "bio"
solver  = "clustalo.align"

[bio.clustalo]
input           = "input.fa"
output_basename = "alignment"
threads         = 1
"#;

#[test]
fn clustalo_runs_end_to_end_against_real_binary() {
    let adapter = ClustaloAdapter::new();
    if skip_if_missing(&adapter) {
        return;
    }

    let case_dir = tempdir("e2e-clustalo");
    std::fs::write(case_dir.join("case.toml"), MINIMAL_CASE_TOML).unwrap();
    std::fs::write(case_dir.join("input.fa"), TINY_FASTA).unwrap();

    let workdir = case_dir.join(".valenx-workdir");
    std::fs::create_dir_all(&workdir).unwrap();

    let case = Case {
        id: "e2e-clustalo".into(),
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
    assert_eq!(report.exit_code, 0, "clustalo exited non-zero: {report:?}");

    let results = adapter.collect(&prepared).expect("collect");
    assert!(
        !results.artifacts.is_empty(),
        "collect found no artifacts in {}",
        workdir.display()
    );

    let _ = std::fs::remove_dir_all(&case_dir);
}
