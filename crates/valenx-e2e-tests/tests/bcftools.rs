//! E2E integration test for the bcftools adapter.
//!
//! Spawns a real `bcftools view` subprocess via the adapter's
//! prepare/run/collect pipeline. Skipped automatically when the
//! upstream binary isn't on PATH.
//!
//! We pick the `view` action on a tiny VCF — the simplest end-to-end
//! shape that exercises real subprocess execution. `view` accepts a
//! VCF and writes a VCF; no reference FASTA required, no multi-input
//! complications.

use valenx_adapter_bcftools::BcftoolsAdapter;
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

/// Minimal valid VCF v4.2 with a single SNV record. bcftools view
/// will round-trip this happily.
///
/// VCF columns (tab-separated): CHROM POS ID REF ALT QUAL FILTER INFO.
const TINY_VCF: &str = "##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=1000>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t100\t.\tA\tG\t40\tPASS\tDP=30\n\
chr1\t200\t.\tC\tT\t40\tPASS\tDP=25\n";

const MINIMAL_CASE_TOML: &str = r#"[case]
physics = "bio"
solver  = "bcftools.view"

[bio.bcftools]
action = "view"
input  = "in.vcf"
output = "out.vcf"
"#;

#[test]
fn bcftools_runs_end_to_end_against_real_binary() {
    let adapter = BcftoolsAdapter::new();
    if skip_if_missing(&adapter) {
        return;
    }

    let case_dir = tempdir("e2e-bcftools");
    std::fs::write(case_dir.join("case.toml"), MINIMAL_CASE_TOML).unwrap();
    std::fs::write(case_dir.join("in.vcf"), TINY_VCF).unwrap();

    let workdir = case_dir.join(".valenx-workdir");
    std::fs::create_dir_all(&workdir).unwrap();

    let case = Case {
        id: "e2e-bcftools".into(),
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
    assert_eq!(report.exit_code, 0, "bcftools exited non-zero: {report:?}");

    let results = adapter.collect(&prepared).expect("collect");
    assert!(
        !results.artifacts.is_empty(),
        "collect found no artifacts in {}",
        workdir.display()
    );

    let _ = std::fs::remove_dir_all(&case_dir);
}
