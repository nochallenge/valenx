//! E2E integration test for the Smoldyn adapter.
//!
//! Spawns a real `smoldyn` subprocess via the adapter's prepare/run/collect
//! pipeline with a trivial 1D reaction-diffusion config. Skipped
//! automatically when the upstream binary isn't on PATH.

use valenx_adapter_smoldyn::SmoldynAdapter;
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

/// A tiny 1D Smoldyn config: 10 particles of species `A` diffusing in
/// a 1D box over 10 time units. `graphics none` turns off the OpenGL
/// viewer (Smoldyn would otherwise try to open a window). The `end_file`
/// directive isn't required — Smoldyn ends after `time_stop`.
const TINY_SMOLDYN_TXT: &str = "# Tiny 1D Smoldyn smoke-test config.
graphics none
dim 1
species A
difc A 1
time_start 0
time_stop 10
time_step 0.1
boundaries 0 0 10
max_mol 100
mol 10 A u
end_file
";

const MINIMAL_CASE_TOML: &str = r#"[case]
physics = "bio"
solver  = "smoldyn.simulate"

[bio.smoldyn]
config = "system.txt"
"#;

#[test]
fn smoldyn_runs_end_to_end_against_real_binary() {
    let adapter = SmoldynAdapter::new();
    if skip_if_missing(&adapter) {
        return;
    }

    let case_dir = tempdir("e2e-smoldyn");
    std::fs::write(case_dir.join("case.toml"), MINIMAL_CASE_TOML).unwrap();
    std::fs::write(case_dir.join("system.txt"), TINY_SMOLDYN_TXT).unwrap();

    let workdir = case_dir.join(".valenx-workdir");
    std::fs::create_dir_all(&workdir).unwrap();

    let case = Case {
        id: "e2e-smoldyn".into(),
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
    assert_eq!(report.exit_code, 0, "smoldyn exited non-zero: {report:?}");

    // Smoldyn's collect surfaces whatever files the user's config wrote
    // out into the workdir. A bare diffusion run with no output
    // directives produces nothing collectable — verifying the run
    // completed cleanly is the main assertion here.
    let _results = adapter.collect(&prepared).expect("collect");

    let _ = std::fs::remove_dir_all(&case_dir);
}
