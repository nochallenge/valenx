//! E2E integration test for the ART (read simulator) adapter.
//!
//! Spawns a real `art_illumina` subprocess via the adapter's
//! prepare/run/collect pipeline at very low coverage so the run
//! finishes in seconds. Skipped automatically when the upstream
//! binary isn't on PATH.

use valenx_adapter_art::ArtAdapter;
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

/// 400 bp reference — long enough for ART to draw a few 50 bp reads
/// at 1x coverage without complaining about an under-length input.
const TINY_REFERENCE_FA: &str = "\
>ref
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
";

const MINIMAL_CASE_TOML: &str = r#"[case]
physics = "bio"
solver  = "art.simulate"

[bio.art]
reference         = "ref.fa"
output_prefix     = "sim"
sequencing_system = "MinS"
read_length       = 50
fold_coverage     = 1.0
"#;

#[test]
fn art_runs_end_to_end_against_real_binary() {
    let adapter = ArtAdapter::new();
    if skip_if_missing(&adapter) {
        return;
    }

    let case_dir = tempdir("e2e-art");
    std::fs::write(case_dir.join("case.toml"), MINIMAL_CASE_TOML).unwrap();
    std::fs::write(case_dir.join("ref.fa"), TINY_REFERENCE_FA).unwrap();

    let workdir = case_dir.join(".valenx-workdir");
    std::fs::create_dir_all(&workdir).unwrap();

    let case = Case {
        id: "e2e-art".into(),
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
    assert_eq!(
        report.exit_code, 0,
        "art_illumina exited non-zero: {report:?}"
    );

    let results = adapter.collect(&prepared).expect("collect");
    assert!(
        !results.artifacts.is_empty(),
        "collect found no artifacts in {}",
        workdir.display()
    );

    let _ = std::fs::remove_dir_all(&case_dir);
}
