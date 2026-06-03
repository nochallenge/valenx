//! E2E integration test for the RAxML-NG adapter.
//!
//! Spawns a real `raxml-ng` subprocess via the adapter's
//! prepare/run/collect pipeline. Skipped automatically when the
//! upstream binary isn't on PATH.

use valenx_adapter_raxml_ng::RaxmlNgAdapter;
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

/// A 4-taxon nucleotide alignment, same shape as the IQ-TREE test.
const TINY_ALIGNMENT_FA: &str = "\
>taxA
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>taxB
ACGTACGTAGGTACGTACGTACGTACGTACGTACGTACGGACGTACGTACGTACGTACGT
>taxC
ACGTACGTAGGTACGTACGTACATACGTACGTACGTACGTACGTACGTACGTACGTACGT
>taxD
ACGTACGTACGTACGTACGCACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
";

const MINIMAL_CASE_TOML: &str = r#"[case]
physics = "bio"
solver  = "raxml-ng.tree"

[bio.raxml-ng]
alignment = "aln.fa"
model     = "JC"
mode      = "search"
threads   = 1
prefix    = "run1"
"#;

#[test]
fn raxml_ng_runs_end_to_end_against_real_binary() {
    let adapter = RaxmlNgAdapter::new();
    if skip_if_missing(&adapter) {
        return;
    }

    let case_dir = tempdir("e2e-raxml-ng");
    std::fs::write(case_dir.join("case.toml"), MINIMAL_CASE_TOML).unwrap();
    std::fs::write(case_dir.join("aln.fa"), TINY_ALIGNMENT_FA).unwrap();

    let workdir = case_dir.join(".valenx-workdir");
    std::fs::create_dir_all(&workdir).unwrap();

    let case = Case {
        id: "e2e-raxml-ng".into(),
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
    assert_eq!(report.exit_code, 0, "raxml-ng exited non-zero: {report:?}");

    let results = adapter.collect(&prepared).expect("collect");
    assert!(
        !results.artifacts.is_empty(),
        "collect found no artifacts in {}",
        workdir.display()
    );

    let _ = std::fs::remove_dir_all(&case_dir);
}
