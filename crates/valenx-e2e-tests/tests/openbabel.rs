//! E2E integration test for the Open Babel adapter.
//!
//! Spawns a real `obabel` subprocess via the adapter's prepare/run/collect
//! pipeline doing a simple SMILES → SDF format conversion. Skipped
//! automatically when the upstream binary isn't on PATH.

use valenx_adapter_openbabel::OpenBabelAdapter;
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

/// A trivial SMILES file — water + ethanol on two lines. Open Babel
/// handles SMILES → SDF natively without external dependencies.
const TINY_SMILES: &str = "O\twater
CCO\tethanol
";

const MINIMAL_CASE_TOML: &str = r#"[case]
physics = "bio"
solver  = "openbabel.convert"

[bio.openbabel]
input  = "input.smi"
output = "output.sdf"
"#;

#[test]
fn openbabel_runs_end_to_end_against_real_binary() {
    let adapter = OpenBabelAdapter::new();
    if skip_if_missing(&adapter) {
        return;
    }

    let case_dir = tempdir("e2e-openbabel");
    std::fs::write(case_dir.join("case.toml"), MINIMAL_CASE_TOML).unwrap();
    std::fs::write(case_dir.join("input.smi"), TINY_SMILES).unwrap();

    let workdir = case_dir.join(".valenx-workdir");
    std::fs::create_dir_all(&workdir).unwrap();

    let case = Case {
        id: "e2e-openbabel".into(),
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
    assert_eq!(report.exit_code, 0, "obabel exited non-zero: {report:?}");

    let _results = adapter.collect(&prepared).expect("collect");
    // The exact artefacts depend on Open Babel's output filter; the
    // important assertion is that the convert subprocess succeeded.

    let _ = std::fs::remove_dir_all(&case_dir);
}
