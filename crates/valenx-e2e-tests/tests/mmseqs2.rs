//! E2E integration test for the MMseqs2 adapter.
//!
//! Spawns a real `mmseqs easy-search` subprocess via the adapter's
//! prepare/run/collect pipeline. MMseqs2's `easy-search` workflow
//! handles DB creation internally — feed it FASTA files for both
//! query and target and it builds whatever scratch DBs it needs.
//! Skipped automatically when the upstream binary isn't on PATH.

use valenx_adapter_mmseqs2::Mmseqs2Adapter;
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

/// A small protein FASTA shared as both query and target — the
/// self-search produces guaranteed hits without needing an external
/// database fixture.
const TINY_PROTEIN_FASTA: &str = "\
>seq1
ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY
>seq2
ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY
";

const MINIMAL_CASE_TOML: &str = r#"[case]
physics = "bio"
solver  = "mmseqs2.search"

[bio.mmseqs2]
action = "easy-search"
query  = "query.fa"
target = "target.fa"
output = "hits.m8"
threads = 1
"#;

#[test]
fn mmseqs2_runs_end_to_end_against_real_binary() {
    let adapter = Mmseqs2Adapter::new();
    if skip_if_missing(&adapter) {
        return;
    }

    let case_dir = tempdir("e2e-mmseqs2");
    std::fs::write(case_dir.join("case.toml"), MINIMAL_CASE_TOML).unwrap();
    std::fs::write(case_dir.join("query.fa"), TINY_PROTEIN_FASTA).unwrap();
    std::fs::write(case_dir.join("target.fa"), TINY_PROTEIN_FASTA).unwrap();

    let workdir = case_dir.join(".valenx-workdir");
    std::fs::create_dir_all(&workdir).unwrap();

    let case = Case {
        id: "e2e-mmseqs2".into(),
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
    assert_eq!(report.exit_code, 0, "mmseqs2 exited non-zero: {report:?}");

    let results = adapter.collect(&prepared).expect("collect");
    assert!(
        !results.artifacts.is_empty(),
        "collect found no artifacts in {}",
        workdir.display()
    );

    let _ = std::fs::remove_dir_all(&case_dir);
}
