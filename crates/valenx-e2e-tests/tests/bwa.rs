//! E2E integration test for the BWA adapter.
//!
//! Spawns a real `bwa` subprocess via the adapter's prepare/run/collect
//! pipeline. Skipped automatically when the upstream binary isn't on
//! PATH — so local developers without the tool installed don't see
//! false failures.

use valenx_adapter_bwa::BwaAdapter;
use valenx_core::{Adapter, CancellationToken, Case, LogLevel, LogSink, ProgressSink, RunContext};
use valenx_test_utils::tempdir;

// ---------------------------------------------------------------------------
// Helpers (inline rather than shared — each `tests/*.rs` compiles into
// its own binary, so a shared module would need `mod` declarations
// and a path attribute. Inline keeps the test self-contained and
// matches the pattern in valenx-app/tests/pipeline_e2e.rs).
// ---------------------------------------------------------------------------

struct NoopProgress;
impl ProgressSink for NoopProgress {
    fn report(&self, _pct: f32, _message: &str) {}
}
struct NoopLog;
impl LogSink for NoopLog {
    fn log_line(&self, _level: LogLevel, _line: &str) {}
}

/// Returns true if the upstream binary isn't installed on PATH. Tests
/// call this first and `return` on `true` so we don't panic on local
/// machines that lack the tool.
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

// ---------------------------------------------------------------------------
// Test
// ---------------------------------------------------------------------------

/// A tiny synthetic reference: 200 bp of fixed-ish nucleotides. BWA's
/// `bwa index` refuses references shorter than its seed length, so we
/// keep this comfortably above the seed-cutoff floor.
const TINY_REFERENCE_FA: &str = "\
>ref
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
";

/// Two tiny single-end reads that match the start of the reference.
/// Q-scores set to a uniform high quality (`I` = 40 in Sanger).
const TINY_READS_FQ: &str = "\
@read1
ACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read2
ACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
";

const MINIMAL_CASE_TOML: &str = r#"[case]
physics = "bio"
solver  = "bwa.mem"

[bio.bwa]
reference = "ref.fa"
reads     = ["reads.fq"]
threads   = 1
"#;

#[test]
fn bwa_runs_end_to_end_against_real_binary() {
    let adapter = BwaAdapter::new();
    if skip_if_missing(&adapter) {
        return;
    }

    let case_dir = tempdir("e2e-bwa");
    std::fs::write(case_dir.join("case.toml"), MINIMAL_CASE_TOML).unwrap();
    std::fs::write(case_dir.join("ref.fa"), TINY_REFERENCE_FA).unwrap();
    std::fs::write(case_dir.join("reads.fq"), TINY_READS_FQ).unwrap();

    let workdir = case_dir.join(".valenx-workdir");
    std::fs::create_dir_all(&workdir).unwrap();

    let case = Case {
        id: "e2e-bwa".into(),
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
    assert_eq!(report.exit_code, 0, "bwa exited non-zero: {report:?}");

    let results = adapter.collect(&prepared).expect("collect");
    assert!(
        !results.artifacts.is_empty(),
        "collect found no artifacts in {}",
        workdir.display()
    );
    // BWA's canonical aligned-reads file is `out.sam` — check it shows up.
    let has_sam = results.artifacts.iter().any(|a| {
        a.path
            .extension()
            .and_then(|s| s.to_str())
            .map(|s| s.eq_ignore_ascii_case("sam"))
            .unwrap_or(false)
    });
    assert!(
        has_sam,
        "expected an out.sam artifact, got {:?}",
        results.artifacts
    );

    let _ = std::fs::remove_dir_all(&case_dir);
}
