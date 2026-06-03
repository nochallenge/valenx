//! E2E integration test for the xTB adapter.
//!
//! Spawns a real `xtb` subprocess via the adapter's prepare/run/collect
//! pipeline doing a single-point energy of a tiny H2O geometry. Skipped
//! automatically when the upstream binary isn't on PATH.

use valenx_adapter_xtb::XtbAdapter;
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

/// A minimal H2O geometry in XYZ format. Line 1: atom count.
/// Line 2: comment. Lines 3+: element + x y z (Å). xtb finishes the
/// single-point in <1 second.
const TINY_XYZ: &str = "3
water single-point
O  0.000000  0.000000  0.000000
H  0.758602  0.000000  0.504284
H -0.758602  0.000000  0.504284
";

const MINIMAL_CASE_TOML: &str = r#"[case]
physics = "bio"
solver  = "xtb.compute"

[bio.xtb]
input = "molecule.xyz"
"#;

#[test]
fn xtb_runs_end_to_end_against_real_binary() {
    let adapter = XtbAdapter::new();
    if skip_if_missing(&adapter) {
        return;
    }

    let case_dir = tempdir("e2e-xtb");
    std::fs::write(case_dir.join("case.toml"), MINIMAL_CASE_TOML).unwrap();
    std::fs::write(case_dir.join("molecule.xyz"), TINY_XYZ).unwrap();

    let workdir = case_dir.join(".valenx-workdir");
    std::fs::create_dir_all(&workdir).unwrap();

    let case = Case {
        id: "e2e-xtb".into(),
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
    assert_eq!(report.exit_code, 0, "xtb exited non-zero: {report:?}");

    let results = adapter.collect(&prepared).expect("collect");
    // xtb writes xtb.log (stdout capture) plus charges / wbo / ... files
    // depending on mode. Single-point produces at least the log.
    assert!(
        !results.artifacts.is_empty(),
        "collect found no artifacts in {}",
        workdir.display()
    );

    let _ = std::fs::remove_dir_all(&case_dir);
}
