//! E2E integration test for the Cas-OFFinder adapter.
//!
//! Spawns a real `cas-offinder` subprocess via the adapter's
//! prepare/run/collect pipeline against a tiny synthetic reference.
//! Skipped automatically when the upstream binary isn't on PATH.

use valenx_adapter_cas_offinder::CasOffinderAdapter;
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

/// A synthetic 100 bp reference. We only need enough sequence for
/// Cas-OFFinder to scan, not enough for a biologically plausible
/// guide hit.
const TINY_REFERENCE_FA: &str = "\
>chr1
ACGTACGTACGTACGTACGTACGTACGTACGGCGGACGTACGTACGTACGTACGTACGTACGT\
ACGTACGTACGTACGTACGTACGTACGTACGT
";

/// Cas-OFFinder input format:
///   line 1 — reference path (file or directory of `.fa` files)
///   line 2 — pattern (N for any, last few bases = PAM)
///   line 3+ — sgRNA sequence (same length as pattern) + mismatch budget
const TINY_CASOFF_INPUT_TEMPLATE: &str = "{REFERENCE}
NNNNNNNNNNNNNNNNNNNNNGG
ACGTACGTACGTACGTACGTACG 3
";

const MINIMAL_CASE_TOML: &str = r#"[case]
physics = "bio"
solver  = "cas-offinder.search"

[bio.cas_offinder]
input   = "query.in"
output  = "hits.tsv"
backend = "C"
"#;

#[test]
fn cas_offinder_runs_end_to_end_against_real_binary() {
    let adapter = CasOffinderAdapter::new();
    if skip_if_missing(&adapter) {
        return;
    }

    let case_dir = tempdir("e2e-cas-offinder");
    std::fs::write(case_dir.join("case.toml"), MINIMAL_CASE_TOML).unwrap();
    std::fs::write(case_dir.join("ref.fa"), TINY_REFERENCE_FA).unwrap();

    // Cas-OFFinder accepts either a single FASTA file or a directory
    // of FASTAs as the reference. Use the absolute path of `ref.fa`
    // so the input file's reference line resolves correctly when
    // Cas-OFFinder runs from its working directory.
    let reference_abs = case_dir.join("ref.fa");
    let input_text =
        TINY_CASOFF_INPUT_TEMPLATE.replace("{REFERENCE}", &reference_abs.display().to_string());
    std::fs::write(case_dir.join("query.in"), input_text).unwrap();

    let workdir = case_dir.join(".valenx-workdir");
    std::fs::create_dir_all(&workdir).unwrap();

    let case = Case {
        id: "e2e-cas-offinder".into(),
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
        "cas-offinder exited non-zero: {report:?}"
    );

    let _results = adapter.collect(&prepared).expect("collect");
    // Cas-OFFinder writes an empty hits.tsv when no off-targets are
    // found — that's still a valid run. Verifying exit code 0 above
    // is the load-bearing assertion.

    let _ = std::fs::remove_dir_all(&case_dir);
}
