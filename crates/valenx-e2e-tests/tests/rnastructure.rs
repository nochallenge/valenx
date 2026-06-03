//! E2E integration test for the RNAstructure (`Fold`) adapter.
//!
//! Spawns a real `Fold` subprocess via the adapter's prepare/run/collect
//! pipeline. Skipped automatically when the upstream binary isn't on
//! PATH.

use valenx_adapter_rnastructure::RnaStructureAdapter;
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

/// A short RNA in RNAstructure's native `.seq` format: comment line,
/// title line, sequence (1-letter U code), terminator. The sequence
/// has a self-complementary 5'/3' region so Fold produces a clean
/// hairpin MFE in milliseconds.
const TINY_RNA_SEQ: &str = ";
hairpin
GGGAAACCCAAAGGGAAACCCAAAGGGAAACCC1
";

const MINIMAL_CASE_TOML: &str = r#"[case]
physics = "bio"
solver  = "rnastructure.fold"

[bio.rnastructure]
input  = "rna.seq"
output = "fold.ct"
"#;

#[test]
fn rnastructure_runs_end_to_end_against_real_binary() {
    let adapter = RnaStructureAdapter::new();
    if skip_if_missing(&adapter) {
        return;
    }

    let case_dir = tempdir("e2e-rnastructure");
    std::fs::write(case_dir.join("case.toml"), MINIMAL_CASE_TOML).unwrap();
    std::fs::write(case_dir.join("rna.seq"), TINY_RNA_SEQ).unwrap();

    let workdir = case_dir.join(".valenx-workdir");
    std::fs::create_dir_all(&workdir).unwrap();

    let case = Case {
        id: "e2e-rnastructure".into(),
        path: case_dir.clone(),
    };

    let prepared = adapter.prepare(&case, &workdir).expect("prepare");

    let cancel = CancellationToken::new();
    let mut ctx = RunContext {
        cancel: &cancel,
        progress: Box::new(NoopProgress),
        log: Box::new(NoopLog),
    };
    let report = match adapter.run(&prepared, &mut ctx) {
        Ok(r) => r,
        Err(e) => {
            // On case-insensitive filesystems (Windows, macOS HFS+) the
            // probe can match GNU coreutils `fold` instead of
            // RNAstructure's `Fold`. The wrong binary surfaces an
            // "unknown option" error — skip rather than fail the test
            // when that happens.
            let msg = format!("{e}");
            if msg.contains("unknown option") || msg.contains("unrecognized option") {
                eprintln!(
                    "Skipping E2E test — `Fold` on PATH appears to be GNU `fold` rather \
                     than RNAstructure's Fold (case-insensitive FS collision): {msg}"
                );
                let _ = std::fs::remove_dir_all(&case_dir);
                return;
            }
            panic!("run: {e:?}");
        }
    };
    assert_eq!(report.exit_code, 0, "Fold exited non-zero: {report:?}");

    let results = adapter.collect(&prepared).expect("collect");
    assert!(
        !results.artifacts.is_empty(),
        "collect found no artifacts in {}",
        workdir.display()
    );

    let _ = std::fs::remove_dir_all(&case_dir);
}
