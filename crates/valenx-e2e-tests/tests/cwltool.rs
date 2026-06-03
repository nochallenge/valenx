//! E2E integration test for the cwltool adapter.
//!
//! Spawns a real `cwltool` subprocess via the adapter's prepare/run/collect
//! pipeline with a trivial CommandLineTool. Skipped automatically when
//! the upstream binary isn't on PATH.

use valenx_adapter_cwltool::CwltoolAdapter;
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

/// cwltool's probe accepts `python3` as a fallback when the `cwltool`
/// console script isn't on PATH (the adapter probe returns ok with a
/// "python only" warning). But prepare() needs the real `cwltool`
/// binary. So we treat any probe warning containing "cwltool not
/// found" as a skip, matching the spirit of the other adapter tests.
fn skip_if_missing(adapter: &dyn Adapter) -> bool {
    match adapter.probe() {
        Ok(report) if report.ok => {
            let python_only = report.warnings.iter().any(|w| {
                w.contains("cwltool not found") || w.contains("install via `pip install cwltool`")
            });
            if python_only {
                eprintln!(
                    "Skipping E2E test — `{}` console script not installed (Python-only fallback).",
                    adapter.info().id
                );
                true
            } else {
                false
            }
        }
        _ => {
            eprintln!(
                "Skipping E2E test — `{}` upstream binary not installed on PATH.",
                adapter.info().id
            );
            true
        }
    }
}

/// A minimal CWL v1.2 CommandLineTool that runs `echo` and captures
/// stdout to a file. No Docker requirement, no external tools beyond
/// the system `echo`.
const TINY_CWL: &str = "cwlVersion: v1.2
class: CommandLineTool
baseCommand: echo
arguments: [\"hello\"]
inputs: []
outputs:
  out:
    type: stdout
stdout: out.txt
";

const MINIMAL_CASE_TOML: &str = r#"[case]
physics = "bio"
solver  = "cwltool.run"

[bio.cwltool]
workflow   = "tool.cwl"
output_dir = "results"
"#;

#[test]
fn cwltool_runs_end_to_end_against_real_binary() {
    let adapter = CwltoolAdapter::new();
    if skip_if_missing(&adapter) {
        return;
    }

    let case_dir = tempdir("e2e-cwltool");
    std::fs::write(case_dir.join("case.toml"), MINIMAL_CASE_TOML).unwrap();
    std::fs::write(case_dir.join("tool.cwl"), TINY_CWL).unwrap();

    let workdir = case_dir.join(".valenx-workdir");
    std::fs::create_dir_all(&workdir).unwrap();

    let case = Case {
        id: "e2e-cwltool".into(),
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
    assert_eq!(report.exit_code, 0, "cwltool exited non-zero: {report:?}");

    // cwltool's collect() returns the output_dir contents plus any
    // top-level *.log files; on a successful trivial workflow run
    // there'll be at least the workflow's `out` output staged.
    let _results = adapter.collect(&prepared).expect("collect");

    let _ = std::fs::remove_dir_all(&case_dir);
}
