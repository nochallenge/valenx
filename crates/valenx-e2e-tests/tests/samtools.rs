//! E2E integration test for the samtools adapter.
//!
//! Spawns a real `samtools flagstat` subprocess via the adapter's
//! prepare/run/collect pipeline. Skipped automatically when the
//! upstream binary isn't on PATH.
//!
//! We pick the `flagstat` action because it's the simplest end-to-end
//! shape: feed a tiny but valid SAM, the adapter redirects stdout to a
//! pinned `flagstat.txt`, and `collect()` should surface that one
//! artifact.

use valenx_adapter_samtools::SamtoolsAdapter;
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

/// Minimal valid SAM: header + two unmapped reads. samtools flagstat
/// is the canonical "smoke test" command — it reads any well-formed
/// SAM/BAM and prints alignment counts. Even unmapped-only records
/// are perfectly valid input.
///
/// Format (tab-separated): QNAME FLAG RNAME POS MAPQ CIGAR RNEXT
/// PNEXT TLEN SEQ QUAL. FLAG 4 marks the read as unmapped, which is
/// the lowest-friction valid record we can construct.
const TINY_SAM: &str = "@HD\tVN:1.6\tSO:unsorted\n\
@SQ\tSN:chr1\tLN:1000\n\
read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n\
read2\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n";

const MINIMAL_CASE_TOML: &str = r#"[case]
physics = "bio"
solver  = "samtools.flagstat"

[bio.samtools]
action = "flagstat"
input  = "aligned.sam"
"#;

#[test]
fn samtools_runs_end_to_end_against_real_binary() {
    let adapter = SamtoolsAdapter::new();
    if skip_if_missing(&adapter) {
        return;
    }

    let case_dir = tempdir("e2e-samtools");
    std::fs::write(case_dir.join("case.toml"), MINIMAL_CASE_TOML).unwrap();
    std::fs::write(case_dir.join("aligned.sam"), TINY_SAM).unwrap();

    let workdir = case_dir.join(".valenx-workdir");
    std::fs::create_dir_all(&workdir).unwrap();

    let case = Case {
        id: "e2e-samtools".into(),
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
    assert_eq!(report.exit_code, 0, "samtools exited non-zero: {report:?}");

    let results = adapter.collect(&prepared).expect("collect");
    assert!(
        !results.artifacts.is_empty(),
        "collect found no artifacts in {}",
        workdir.display()
    );
    // The adapter pins flagstat's stdout to `flagstat.txt`.
    assert!(
        workdir.join("flagstat.txt").is_file(),
        "expected flagstat.txt in workdir, found artifacts {:?}",
        results.artifacts
    );

    let _ = std::fs::remove_dir_all(&case_dir);
}
