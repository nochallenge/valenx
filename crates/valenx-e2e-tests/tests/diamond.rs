//! E2E integration test for the DIAMOND adapter.
//!
//! Spawns a real `diamond makedb` subprocess via the adapter's
//! prepare/run/collect pipeline. We pick `makedb` rather than
//! `blastp` because makedb requires only a FASTA — the search path
//! would need a pre-built `.dmnd` database, which the adapter validates
//! as an existing file before run() and would require two adapter
//! invocations to bootstrap from scratch. Skipped automatically when
//! the upstream binary isn't on PATH.

use valenx_adapter_diamond::DiamondAdapter;
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

/// Two short protein sequences. DIAMOND's `makedb` is happy with any
/// valid FASTA — we don't need a realistic-sized proteome to exercise
/// the prepare/run/collect path.
const TINY_PROTEIN_FASTA: &str = "\
>seq1
ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY
>seq2
ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY
";

const MINIMAL_CASE_TOML: &str = r#"[case]
physics = "bio"
solver  = "diamond.makedb"

[bio.diamond]
action   = "makedb"
query    = "proteins.fa"
database = "proteins"
output   = "ignored.m8"
threads  = 1
"#;

#[test]
fn diamond_runs_end_to_end_against_real_binary() {
    let adapter = DiamondAdapter::new();
    if skip_if_missing(&adapter) {
        return;
    }

    let case_dir = tempdir("e2e-diamond");
    std::fs::write(case_dir.join("case.toml"), MINIMAL_CASE_TOML).unwrap();
    std::fs::write(case_dir.join("proteins.fa"), TINY_PROTEIN_FASTA).unwrap();

    let workdir = case_dir.join(".valenx-workdir");
    std::fs::create_dir_all(&workdir).unwrap();

    let case = Case {
        id: "e2e-diamond".into(),
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
    assert_eq!(report.exit_code, 0, "diamond exited non-zero: {report:?}");

    let results = adapter.collect(&prepared).expect("collect");
    assert!(
        !results.artifacts.is_empty(),
        "collect found no artifacts in {}",
        workdir.display()
    );

    let _ = std::fs::remove_dir_all(&case_dir);
}
