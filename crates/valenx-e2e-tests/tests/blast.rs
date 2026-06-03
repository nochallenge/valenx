//! E2E integration test for the BLAST+ adapter.
//!
//! Spawns a real `blastn` subprocess via the adapter's
//! prepare/run/collect pipeline. The adapter only wraps the search
//! programs (`blastn`/`blastp`/`blastx`/`tblastn`/`tblastx`) — it
//! does not call `makeblastdb`. So we bootstrap a tiny BLAST database
//! in the test setup with the canonical `makeblastdb` invocation
//! (also part of the BLAST+ package, so the presence check covers it
//! for free). Skipped automatically when the upstream binary isn't on
//! PATH.

use std::process::Command;

use valenx_adapter_blast::BlastAdapter;
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

/// 200 bp synthetic nucleotide reference — long enough for makeblastdb
/// to accept and produce a valid `.nhr/.nin/.nsq` triple.
const TINY_REFERENCE_FA: &str = "\
>ref
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
";

/// Query that should align perfectly to the reference start.
const TINY_QUERY_FA: &str = "\
>query1
ACGTACGTACGTACGTACGTACGTACGTACGTACGT
";

const MINIMAL_CASE_TOML: &str = r#"[case]
physics = "bio"
solver  = "blast.search"

[bio.blast]
program  = "blastn"
query    = "query.fa"
database = "refdb"
evalue   = 10.0
outfmt   = 6
threads  = 1
"#;

#[test]
fn blast_runs_end_to_end_against_real_binary() {
    let adapter = BlastAdapter::new();
    if skip_if_missing(&adapter) {
        return;
    }

    let case_dir = tempdir("e2e-blast");

    let reference_path = case_dir.join("ref.fa");
    std::fs::write(&reference_path, TINY_REFERENCE_FA).unwrap();
    std::fs::write(case_dir.join("query.fa"), TINY_QUERY_FA).unwrap();

    // Bootstrap a BLAST nucleotide database via `makeblastdb`. If this
    // fails we skip rather than panic — matches the contract of every
    // other E2E test on this crate.
    let build = Command::new("makeblastdb")
        .args(["-dbtype", "nucl", "-in"])
        .arg(&reference_path)
        .args(["-out"])
        .arg(case_dir.join("refdb"))
        .output();
    match build {
        Ok(out) if out.status.success() => {}
        Ok(out) => {
            eprintln!(
                "Skipping E2E test — `makeblastdb` failed (exit {:?}): {}",
                out.status.code(),
                String::from_utf8_lossy(&out.stderr)
            );
            let _ = std::fs::remove_dir_all(&case_dir);
            return;
        }
        Err(e) => {
            eprintln!("Skipping E2E test — could not spawn `makeblastdb`: {e}");
            let _ = std::fs::remove_dir_all(&case_dir);
            return;
        }
    }

    std::fs::write(case_dir.join("case.toml"), MINIMAL_CASE_TOML).unwrap();

    let workdir = case_dir.join(".valenx-workdir");
    std::fs::create_dir_all(&workdir).unwrap();

    let case = Case {
        id: "e2e-blast".into(),
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
    assert_eq!(report.exit_code, 0, "blast exited non-zero: {report:?}");

    let results = adapter.collect(&prepared).expect("collect");
    assert!(
        !results.artifacts.is_empty(),
        "collect found no artifacts in {}",
        workdir.display()
    );

    let _ = std::fs::remove_dir_all(&case_dir);
}
