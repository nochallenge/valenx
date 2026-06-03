//! E2E integration test for the MAFFT adapter.
//!
//! Spawns a real `mafft` subprocess via the adapter's prepare/run/collect
//! pipeline. Skipped automatically when the upstream binary isn't on
//! PATH.
//!
//! We feed MAFFT three short, similar sequences so the "auto" strategy
//! picks a fast progressive-alignment path (FFT-NS-2). The adapter
//! redirects stdout to `aligned.fa` in the workdir; `collect()` should
//! surface that as a FASTA artifact.

use valenx_adapter_mafft::MafftAdapter;
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

/// Three similar but distinct DNA sequences, enough to exercise
/// progressive alignment without being trivially identical (MAFFT
/// handles identical inputs but a real alignment is what we want to
/// smoke-test). Sequences are ~60 bp each.
const TINY_FASTA: &str = "\
>seq1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>seq2
ACGTACGTACGTACGTACGTACGGACGTACGTACGTACGTACGTACGTACGTACGTACGT
>seq3
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAAATACGTACGTACGTACGT
";

const MINIMAL_CASE_TOML: &str = r#"[case]
physics = "bio"
solver  = "mafft.msa"

[bio.mafft]
input    = "seqs.fa"
strategy = "auto"
threads  = 1
"#;

#[test]
fn mafft_runs_end_to_end_against_real_binary() {
    let adapter = MafftAdapter::new();
    if skip_if_missing(&adapter) {
        return;
    }

    let case_dir = tempdir("e2e-mafft");
    std::fs::write(case_dir.join("case.toml"), MINIMAL_CASE_TOML).unwrap();
    std::fs::write(case_dir.join("seqs.fa"), TINY_FASTA).unwrap();

    let workdir = case_dir.join(".valenx-workdir");
    std::fs::create_dir_all(&workdir).unwrap();

    let case = Case {
        id: "e2e-mafft".into(),
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
    assert_eq!(report.exit_code, 0, "mafft exited non-zero: {report:?}");

    let results = adapter.collect(&prepared).expect("collect");
    assert!(
        !results.artifacts.is_empty(),
        "collect found no artifacts in {}",
        workdir.display()
    );
    // The canonical MAFFT output is `aligned.fa` in the workdir.
    assert!(
        workdir.join("aligned.fa").is_file(),
        "expected aligned.fa in workdir, found artifacts {:?}",
        results.artifacts
    );

    let _ = std::fs::remove_dir_all(&case_dir);
}
