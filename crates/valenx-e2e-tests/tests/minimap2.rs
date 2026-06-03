//! E2E integration test for the minimap2 adapter.
//!
//! Spawns a real `minimap2` subprocess via the adapter's
//! prepare/run/collect pipeline. Skipped automatically when the
//! upstream binary isn't on PATH.

use valenx_adapter_minimap2::Minimap2Adapter;
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

/// 200 bp of synthetic reference — minimap2 builds its index on the
/// fly from the FASTA, so the only requirement is that the reference
/// be longer than the k-mer seed length.
const TINY_REFERENCE_FA: &str = "\
>ref
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
";

/// A short Illumina-style read that matches the start of the
/// reference. Uniform high-quality scores (`I` = 40 in Sanger).
const TINY_READS_FQ: &str = "\
@read1
ACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
";

const MINIMAL_CASE_TOML: &str = r#"[case]
physics = "bio"
solver  = "minimap2.align"

[bio.minimap2]
reference = "ref.fa"
reads     = ["reads.fq"]
preset    = "sr"
threads   = 1
"#;

#[test]
fn minimap2_runs_end_to_end_against_real_binary() {
    let adapter = Minimap2Adapter::new();
    if skip_if_missing(&adapter) {
        return;
    }

    let case_dir = tempdir("e2e-minimap2");
    std::fs::write(case_dir.join("case.toml"), MINIMAL_CASE_TOML).unwrap();
    std::fs::write(case_dir.join("ref.fa"), TINY_REFERENCE_FA).unwrap();
    std::fs::write(case_dir.join("reads.fq"), TINY_READS_FQ).unwrap();

    let workdir = case_dir.join(".valenx-workdir");
    std::fs::create_dir_all(&workdir).unwrap();

    let case = Case {
        id: "e2e-minimap2".into(),
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
    assert_eq!(report.exit_code, 0, "minimap2 exited non-zero: {report:?}");

    let results = adapter.collect(&prepared).expect("collect");
    assert!(
        !results.artifacts.is_empty(),
        "collect found no artifacts in {}",
        workdir.display()
    );

    let _ = std::fs::remove_dir_all(&case_dir);
}
