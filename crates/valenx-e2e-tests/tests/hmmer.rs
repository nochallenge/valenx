//! E2E integration test for the HMMER adapter.
//!
//! Spawns a real `hmmsearch` subprocess via the adapter's
//! prepare/run/collect pipeline. Skipped automatically when the
//! upstream binaries (`hmmsearch` + `hmmbuild`) aren't on PATH.
//!
//! HMMER's adapter only wraps `hmmsearch` / `hmmscan` — both require an
//! existing `.hmm` profile. So we bootstrap a tiny profile in the test
//! setup via `hmmbuild` (also part of the HMMER package, so the
//! presence check covers it for free). The point of the E2E test is to
//! confirm `hmmsearch` actually runs through prepare/run/collect, not
//! to exercise `hmmbuild`.

use std::process::Command;

use valenx_adapter_hmmer::HmmerAdapter;
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

/// A tiny 3-sequence protein alignment in Stockholm format — the
/// canonical input shape for `hmmbuild`. All three sequences are short
/// (15 residues) and trivially similar, so `hmmbuild` finishes in
/// milliseconds.
const TINY_STOCKHOLM: &str = "\
# STOCKHOLM 1.0
seq1            ACDEFGHIKLMNPQR
seq2            ACDEFGHIKLMNPQR
seq3            ACDEFGHIKLMNPQR
//
";

/// A tiny protein sequence database (single sequence, matches the
/// profile). `hmmsearch` reads this as the search target.
const TINY_PROTEIN_FASTA: &str = "\
>target1
ACDEFGHIKLMNPQRACDEFGHIKLMNPQRACDEFGHIKLMNPQR
";

const MINIMAL_CASE_TOML: &str = r#"[case]
physics = "bio"
solver  = "hmmer.search"

[bio.hmmer]
tool      = "hmmsearch"
profile   = "profile.hmm"
sequences = "targets.fa"
cpus      = 1
evalue    = 10.0
"#;

#[test]
fn hmmer_runs_end_to_end_against_real_binary() {
    let adapter = HmmerAdapter::new();
    if skip_if_missing(&adapter) {
        return;
    }

    let case_dir = tempdir("e2e-hmmer");

    // Stage the Stockholm-format alignment + target FASTA next to the
    // case.toml. The profile.hmm we build below ends up here too.
    let stockholm_path = case_dir.join("seed.sto");
    let profile_path = case_dir.join("profile.hmm");
    let targets_path = case_dir.join("targets.fa");
    std::fs::write(&stockholm_path, TINY_STOCKHOLM).unwrap();
    std::fs::write(&targets_path, TINY_PROTEIN_FASTA).unwrap();

    // Bootstrap the profile with `hmmbuild`. If this fails we skip
    // (rather than panic) — same contract as the adapter's own probe.
    let build = Command::new("hmmbuild")
        .arg(&profile_path)
        .arg(&stockholm_path)
        .output();
    match build {
        Ok(out) if out.status.success() => {}
        Ok(out) => {
            eprintln!(
                "Skipping E2E test — `hmmbuild` failed (exit {:?}): {}",
                out.status.code(),
                String::from_utf8_lossy(&out.stderr)
            );
            let _ = std::fs::remove_dir_all(&case_dir);
            return;
        }
        Err(e) => {
            eprintln!("Skipping E2E test — could not spawn `hmmbuild`: {e}");
            let _ = std::fs::remove_dir_all(&case_dir);
            return;
        }
    }

    std::fs::write(case_dir.join("case.toml"), MINIMAL_CASE_TOML).unwrap();

    let workdir = case_dir.join(".valenx-workdir");
    std::fs::create_dir_all(&workdir).unwrap();

    let case = Case {
        id: "e2e-hmmer".into(),
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
    assert_eq!(report.exit_code, 0, "hmmsearch exited non-zero: {report:?}");

    let results = adapter.collect(&prepared).expect("collect");
    assert!(
        !results.artifacts.is_empty(),
        "collect found no artifacts in {}",
        workdir.display()
    );
    // The adapter pins HMMER outputs to `tblout.txt` + `hmmer.out`.
    assert!(
        workdir.join("tblout.txt").is_file() || workdir.join("hmmer.out").is_file(),
        "expected tblout.txt or hmmer.out in workdir, found artifacts {:?}",
        results.artifacts
    );

    let _ = std::fs::remove_dir_all(&case_dir);
}
