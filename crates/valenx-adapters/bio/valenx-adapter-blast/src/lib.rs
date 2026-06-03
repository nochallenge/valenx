//! # valenx-adapter-blast
//!
//! Adapter for [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/) — the
//! foundational sequence-search suite. BLAST (Basic Local Alignment
//! Search Tool) is the most widely used bioinformatics tool on the
//! planet: it finds local alignments between a query sequence and a
//! database, scoring matches by E-value (the expected number of hits of
//! equal-or-better score under a random model). The 2.x line of BLAST+
//! ships from NCBI as a small fleet of binaries — `blastn`, `blastp`,
//! `blastx`, `tblastn`, `tblastx` — covering every nucleotide/protein
//! query/database combination.
//!
//! **Phase 18.7 — subprocess wrapper around the BLAST+ binaries.** The
//! user picks the program via `[bio.blast].program` in `case.toml` and
//! supplies a FASTA query plus a database prefix path. `prepare()`
//! resolves both against the case directory, looks up the matching
//! BLAST+ binary on PATH, and composes the invocation with a pinned
//! output filename so `collect()` can find the results
//! deterministically.
//!
//! Because BLAST writes its results to a file (`-out blast_results.txt`)
//! and is otherwise quiet on stdout/stderr, `run()` uses the shared
//! [`valenx_core::subprocess::run`] runner with a minimal line handler
//! that only lifts `[E::` / `ERROR` lines as warnings.
//!
//! NCBI BLAST is US Government Public Domain — produced by employees of
//! the federal government in the course of their official duties — so
//! it ships with no usage restrictions.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;
pub mod native;

use std::ffi::OsString;
use std::fs;
use std::path::Path;
use std::time::Duration;

use semver::Version;

use valenx_core::{
    adapter_helpers::{detect_tool_version_semver, find_on_path, live_provenance},
    error::RunPhase,
    subprocess, Adapter, AdapterError, AdapterInfo, Capabilities, Case, LicenseMode, Physics,
    PreparedJob, ProbeReport, RunContext, RunReport, VersionRange,
};
use valenx_fields::{
    artifact::{Artifact, ArtifactKind},
    Results,
};

use crate::case_input::BlastInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(BlastAdapter::new())
}

pub struct BlastAdapter;

impl BlastAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for BlastAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "blast";

/// BLAST+'s probe binaries. `blastn` and `blastp` ship with every BLAST+
/// install (one or both is always present); finding either on PATH is
/// enough to call the suite installed. The actual binary used at run
/// time is looked up in `prepare()` based on the `program` field.
const PROBE_BINARIES: &[&str] = &["blastn", "blastp"];

/// The results filename written by `-out`. Pinned so the `prepare()`
/// invocation, the `collect()` walk, and the artifact label all agree
/// on what to look for.
const OUT_RESULTS: &str = "blast_results.txt";

impl Adapter for BlastAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "BLAST+",
            // BLAST+ 2.10 (2019) is the long-running stable line every
            // distro ships; major output-format additions (SAM=17,
            // VCF=18) landed in the 2.10 series. Floor at 2.10.0
            // covers every reasonably modern install; upper bound 3.0
            // reserves room for the next major.
            version_range: VersionRange {
                min_inclusive: Version::new(2, 10, 0),
                max_exclusive: Version::new(3, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            // NCBI BLAST is produced by US federal employees and is
            // therefore in the public domain (17 U.S.C. § 105). No
            // copyright applies.
            tool_license: "Public Domain",
            docs_url: "https://www.ncbi.nlm.nih.gov/books/NBK279690/",
            homepage_url: "https://blast.ncbi.nlm.nih.gov/Blast.cgi",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(PROBE_BINARIES) {
            Some(binary_path) => {
                // `blastn -version` (and the matching command on every
                // BLAST+ binary) prints "blastn: 2.13.0+" on stdout —
                // the trailing `+` is a quirk of the BLAST+ versioning
                // scheme that semver's parser handles fine after the
                // generic detector strips it.
                let found_version =
                    detect_tool_version_semver(&binary_path, &["-version", "--version"]);
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    warnings: Vec::new(),
                    required_env: Vec::new(),
                })
            }
            // Native Rust fallback: BLAST-class seed-and-extend via valenx-align.
            // Requires the database to be (or have a companion) FASTA file.
            None => Ok(ProbeReport {
                ok: true,
                found_version: None,
                binary_path: None,
                warnings: vec![
                    "BLAST+ binary not found; using native Rust seed-and-extend search \
                     (valenx-align). Install BLAST+ 2.10+ via apt/brew/conda for the full \
                     BLAST+ suite (formatted databases, translated search, XML output). \
                     Native mode requires the database to be a FASTA file (or have a \
                     companion .fa/.fasta/.fna/.faa file next to the database prefix)."
                        .to_string(),
                ],
                required_env: Vec::new(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = BlastInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve the query FASTA against the case directory if relative.
        // Same convention as every other Phase 17/18 bio adapter —
        // `query = "query.fa"` next to `case.toml`.
        let source_query = if input.query.is_absolute() {
            input.query.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(
            &case.path,
            &input.query,
        )?
        };
        if !source_query.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.blast].query `{}` not found (resolved {})",
                    input.query.display(),
                    source_query.display()
                ),
            });
        }

        // Resolve the database prefix against the case directory if
        // relative. The prefix is NOT a single file — BLAST databases
        // are sets of three (`.nhr/.nin/.nsq` for nucleotide,
        // `.phr/.pin/.psq` for protein) — so we validate that the
        // *directory* containing the prefix exists, but don't require
        // any specific extension to be present. BLAST itself will
        // produce a clean error message at run time if the database
        // files are missing.
        // Round-9 hardening: relative database prefixes flow into the
        // `-db` flag of the subprocess, so a `..`-traversing or
        // absolute relative form would let a hostile case point BLAST
        // at arbitrary files. Wrap with `confined_join`; absolute
        // paths remain a separate explicit branch since the prefix can
        // legitimately point at a shared sysadmin-managed location
        // (`/var/lib/blast/nr`).
        let source_database = if input.database.is_absolute() {
            input.database.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.database)?
        };
        let database_dir = source_database.parent().unwrap_or_else(|| Path::new("."));
        if !database_dir.is_dir() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.blast].database `{}` directory `{}` not found",
                    input.database.display(),
                    database_dir.display()
                ),
            });
        }

        // Write native_params.toml for the native path. Resolve the
        // database FASTA by looking for a companion FASTA file.
        let db_fasta_path = native::find_db_fasta(&source_database)
            .map(|p| p.to_string_lossy().into_owned())
            .unwrap_or_else(|| source_database.to_string_lossy().into_owned());

        let native_params = native::NativeBlastParams {
            query_path: source_query
                .to_str()
                .ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "query path is not valid UTF-8: {}",
                        source_query.display()
                    ))
                })?
                .to_string(),
            db_fasta_path,
            evalue: input.evalue,
            max_hits: 500,
            output_name: OUT_RESULTS.to_string(),
        };
        native::write_params(workdir, &native_params)?;

        // Look up the binary matching the user-selected `program`.
        // Fall back to native sentinel when not found.
        let native_command: Vec<OsString> =
            match find_on_path(&[input.program.as_str()]) {
                Some(binary_path) => {
                    let mut cmd: Vec<OsString> = vec![
                        binary_path.into_os_string(),
                        OsString::from("-query"),
                        source_query.into_os_string(),
                        OsString::from("-db"),
                        source_database.into_os_string(),
                        OsString::from("-out"),
                        OsString::from(OUT_RESULTS),
                        OsString::from("-evalue"),
                        OsString::from(format_evalue(input.evalue)),
                        OsString::from("-outfmt"),
                        OsString::from(input.outfmt.to_string()),
                        OsString::from("-num_threads"),
                        OsString::from(input.threads.to_string()),
                    ];
                    for arg in &input.extra_args {
                        cmd.push(OsString::from(arg));
                    }
                    cmd
                }
                None => vec![OsString::from(native::NATIVE_SENTINEL)],
            };

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Single-protein-vs-Pfam runs in seconds; a query against
            // nt (~150 GB) on a single node is hours. 4 hours covers
            // the long tail without being absurd.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        // Native Rust path: seed-and-extend search, no subprocess.
        if job.native_command.first().map(|s| s.as_os_str())
            == Some(native::NATIVE_SENTINEL.as_ref())
        {
            return native::run_native(&job.workdir, ctx);
        }

        let report = subprocess::run(job, ctx, "starting BLAST", |line| {
            let mut hint = subprocess::Hint::default();
            // BLAST is status-quiet on stderr — it just writes results
            // to the `-out` file. The only lines worth lifting are
            // explicit error markers; we tag those as warnings so the
            // UI surfaces them but the run still completes.
            if line.contains("[E::") || line.contains("ERROR") {
                hint.warning = Some(line.trim().to_string());
            }
            hint
        })?;
        Ok(RunReport {
            exit_code: report.exit_code,
            wall_time: report.wall_time,
            converged: Some(true),
            residual_history: Vec::new(),
            warnings: report.warnings,
            final_phase: Some(RunPhase::Shutdown),
        })
    }

    fn collect(&self, job: &PreparedJob) -> Result<Results, AdapterError> {
        // Provenance: hash the staged results file if present (the
        // canonical run output). Falls back to case.toml when the search
        // hasn't produced results yet — keeps the provenance block
        // well-formed for partial / failed runs.
        let case_hash_input = {
            let results = job.workdir.join(OUT_RESULTS);
            if results.is_file() {
                results
            } else {
                job.workdir.join("case.toml")
            }
        };
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "BLAST+",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk the workdir top-level. BLAST writes the canonical
        // `blast_results.txt` plus optional `.log` files if the user
        // redirected stderr.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-blast", ?e, "workdir read failed");
                return Ok(results);
            }
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if !path.is_file() {
                continue;
            }
            let file_name = path
                .file_name()
                .and_then(|s| s.to_str())
                .map(|s| s.to_string());
            let ext = path
                .extension()
                .and_then(|s| s.to_str())
                .map(|s| s.to_ascii_lowercase());

            let (kind, label) = match (file_name.as_deref(), ext.as_deref()) {
                (Some(name), _) if name == OUT_RESULTS => {
                    (ArtifactKind::Tabular, "BLAST search results".to_string())
                }
                (_, Some("log")) => (ArtifactKind::Log, "BLAST log".to_string()),
                _ => continue,
            };
            artefacts.push(Artifact {
                path,
                kind,
                checksum: None,
                label,
            });
        }
        artefacts.sort_by(|a, b| a.path.cmp(&b.path));
        results.artifacts = artefacts;
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        // The bio-specific Capability variants land in a follow-up
        // task; ribbon contributions are already enough for the
        // registry to surface the adapter.
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.blast.search"],
        }
    }
}

/// Format an E-value for the BLAST CLI. BLAST accepts any C-style
/// floating-point literal, but it's nicer to pass `1e-5` rather than
/// `0.00001` so the user can recognise their own input in the recorded
/// command. `{e}` formatting in Rust produces the canonical scientific
/// form (e.g. `1e-5`); for whole-number defaults like `10.0` the
/// `{e}` form would be `1e1` which is uglier than just `10` — fall
/// back to the default float formatting in that case.
fn format_evalue(evalue: f64) -> String {
    if evalue.fract() == 0.0 && (1.0..=1e6).contains(&evalue) {
        format!("{evalue}")
    } else {
        format!("{evalue:e}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = BlastAdapter::new().info();
        assert_eq!(info.id, "blast");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "Public Domain");
        assert_eq!(info.display_name, "BLAST+");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = BlastAdapter::new().info();
        // BLAST+ 2.10.x is the de facto stable line every distro
        // ships; 3.0 reserves room for an eventual major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(2, 10, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(3, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = BlastAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.blast.search"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = BlastAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    #[test]
    fn probe_always_succeeds() {
        let report = BlastAdapter::new().probe().unwrap();
        assert!(report.ok, "probe() must return ok=true in all modes");
    }

    #[test]
    fn native_sentinel_is_stable() {
        assert_eq!(native::NATIVE_SENTINEL, "valenx:native:blast");
    }

    /// Round-9 RED→GREEN: `[bio.blast].database` used to be joined
    /// with bare `case.path.join`, which let a hostile case bundle
    /// supply `database = "../../etc/passwd"` and have the BLAST
    /// subprocess open whatever file the user could reach. The fix
    /// wraps the relative branch with `confined_join`; absolute paths
    /// remain explicit (admin-managed shared DBs like `/var/lib/blast/nr`).
    #[test]
    fn prepare_rejects_database_traversing_outside_case_dir() {
        use valenx_test_utils::tempdir;
        let d = tempdir("blast-database-trav");
        std::fs::write(d.join("query.fa"), b">x\nACGT\n").unwrap();
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "blast.search"

[bio.blast]
program  = "blastn"
query    = "query.fa"
database = "../../etc/passwd"
"#,
        )
        .unwrap();
        let case = Case {
            id: "blast-database-trav".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = BlastAdapter::new().prepare(&case, &workdir).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("..") || msg.contains("stay within") || msg.contains("escape"),
            "expected confined_join rejection, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
