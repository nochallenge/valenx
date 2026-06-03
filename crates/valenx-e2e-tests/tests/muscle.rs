//! E2E integration test for the MUSCLE adapter.
//!
//! Spawns a real `muscle` subprocess via the adapter's
//! prepare/run/collect pipeline. Skipped automatically when the
//! upstream binary isn't on PATH.

use valenx_adapter_muscle::MuscleAdapter;
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

/// Three short, similar DNA sequences. MUSCLE 5's `align` path picks
/// the small-input branch and finishes in milliseconds.
const TINY_FASTA: &str = "\
>seq1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>seq2
ACGTACGTACGTACGTACGTACGGACGTACGTACGTACGTACGTACGTACGT
>seq3
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAAATACGTACGT
";

const MINIMAL_CASE_TOML: &str = r#"[case]
physics = "bio"
solver  = "muscle.msa"

[bio.muscle]
input = "seqs.fa"
"#;

#[test]
fn muscle_runs_end_to_end_against_real_binary() {
    let adapter = MuscleAdapter::new();
    if skip_if_missing(&adapter) {
        return;
    }

    let case_dir = tempdir("e2e-muscle");
    std::fs::write(case_dir.join("case.toml"), MINIMAL_CASE_TOML).unwrap();
    std::fs::write(case_dir.join("seqs.fa"), TINY_FASTA).unwrap();

    let workdir = case_dir.join(".valenx-workdir");
    std::fs::create_dir_all(&workdir).unwrap();

    let case = Case {
        id: "e2e-muscle".into(),
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
    assert_eq!(report.exit_code, 0, "muscle exited non-zero: {report:?}");

    let results = adapter.collect(&prepared).expect("collect");
    assert!(
        !results.artifacts.is_empty(),
        "collect found no artifacts in {}",
        workdir.display()
    );

    let _ = std::fs::remove_dir_all(&case_dir);
}
