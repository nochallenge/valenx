//! E2E integration test for the MrBayes adapter.
//!
//! Spawns a real `mb` subprocess via the adapter's prepare/run/collect
//! pipeline with a trivially small (1000-generation) MCMC chain so the
//! test finishes in seconds. Skipped automatically when the upstream
//! binary isn't on PATH.

use valenx_adapter_mrbayes::MrBayesAdapter;
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

/// A minimal NEXUS file with 4 short taxa + a MrBayes block running a
/// 1000-generation MCMC. The embedded `set autoclose=yes nowarn=yes`
/// + `quit` ensures the run exits cleanly without prompting on stdin.
const TINY_NEXUS: &str = "#NEXUS
begin data;
dimensions ntax=4 nchar=20;
format datatype=dna interleave=no gap=-;
matrix
taxonA ACGTACGTACGTACGTACGT
taxonB ACGTACGTACGTACGTACGC
taxonC ACGTAGGTACGTACGTACGT
taxonD ACGTACGTACGTACGTACGG
;
end;

begin mrbayes;
   set autoclose=yes nowarn=yes;
   lset nst=1 rates=equal;
   mcmcp ngen=1000 samplefreq=100 printfreq=100 nchains=1 nruns=1;
   mcmc;
   quit;
end;
";

const MINIMAL_CASE_TOML: &str = r#"[case]
physics = "bio"
solver  = "mrbayes.mcmc"

[bio.mrbayes]
nexus = "data.nex"
"#;

#[test]
fn mrbayes_runs_end_to_end_against_real_binary() {
    let adapter = MrBayesAdapter::new();
    if skip_if_missing(&adapter) {
        return;
    }

    let case_dir = tempdir("e2e-mrbayes");
    std::fs::write(case_dir.join("case.toml"), MINIMAL_CASE_TOML).unwrap();
    std::fs::write(case_dir.join("data.nex"), TINY_NEXUS).unwrap();

    let workdir = case_dir.join(".valenx-workdir");
    std::fs::create_dir_all(&workdir).unwrap();

    let case = Case {
        id: "e2e-mrbayes".into(),
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
    assert_eq!(report.exit_code, 0, "mb exited non-zero: {report:?}");

    let results = adapter.collect(&prepared).expect("collect");
    assert!(
        !results.artifacts.is_empty(),
        "collect found no artifacts in {}",
        workdir.display()
    );

    let _ = std::fs::remove_dir_all(&case_dir);
}
