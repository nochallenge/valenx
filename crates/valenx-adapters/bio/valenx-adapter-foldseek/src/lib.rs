//! # valenx-adapter-foldseek
//!
//! Adapter for [FoldSeek](https://github.com/steineggerlab/foldseek) —
//! the Steinegger lab's structure-vs-structure protein search engine.
//! FoldSeek is the 3D analogue of MMseqs2 (also Steinegger lab):
//! instead of searching amino-acid sequences against a sequence
//! database, it represents each protein structure as a string in the
//! "3Di" alphabet (a 20-letter encoding of local backbone geometry)
//! and runs MMseqs2-style fast prefilter + alignment over the 3Di
//! sequences. The result is structure-search at orders of magnitude
//! the speed of TM-align or DALI, with comparable accuracy.
//!
//! **Phase 17.7 — subprocess wrapper around `foldseek easy-search`.**
//! The user supplies a structure query (PDB or mmCIF) plus a FoldSeek
//! database prefix in `[bio.foldseek]`. `prepare()` resolves both
//! against the case directory, looks up `foldseek` on PATH, and
//! composes the canonical invocation:
//!
//! ```text
//! foldseek easy-search <query> <database> <output_basename>.m8 \
//!     tmp_<output_basename> --threads <N> [extras...]
//! ```
//!
//! The trailing `tmp_<output_basename>` directory is workdir-relative
//! scratch space FoldSeek requires for intermediate database files;
//! `collect()` deliberately does **not** surface it as an artifact —
//! it's pure intermediate state.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

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

use crate::case_input::FoldseekInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(FoldseekAdapter::new())
}

pub struct FoldseekAdapter;

impl FoldseekAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for FoldseekAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "foldseek";
/// FoldSeek's binary candidates. The on-disk binary is `foldseek` —
/// Bioconda, Homebrew, and source builds all install under that name.
const BINARIES: &[&str] = &["foldseek"];

impl Adapter for FoldseekAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "FoldSeek",
            // FoldSeek 8.x is the long-running stable line (commit-hash
            // suffixed releases like 8-ef4e960); the spec pins the
            // supported band at 8.0.0..10.0.0.
            version_range: VersionRange {
                min_inclusive: Version::new(8, 0, 0),
                max_exclusive: Version::new(10, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-3.0",
            docs_url: "https://github.com/steineggerlab/foldseek",
            homepage_url: "https://search.foldseek.com/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `foldseek version` (positional, MMseqs2-style — they
                // share a CLI ancestor) prints the version banner on
                // stdout; the combined scanner picks it up.
                let found_version =
                    detect_tool_version_semver(&binary_path, &["version", "--version"]);
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    warnings: Vec::new(),
                    required_env: Vec::new(),
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "FoldSeek 8+ required; install via `brew install foldseek` \
                       or `conda install -c bioconda foldseek`"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = FoldseekInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.foldseek].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Resolve the query structure against the case directory if
        // relative. Same convention as every other Phase 17 bio
        // adapter — `query = "query.pdb"` next to `case.toml`.
        let source_query = if input.query.is_absolute() {
            input.query.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.query)?
        };
        if !source_query.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.foldseek].query `{}` not found (resolved {})",
                    input.query.display(),
                    source_query.display()
                ),
            });
        }

        // Resolve the database prefix against the case directory if
        // relative. The prefix is NOT a single file — FoldSeek
        // databases are sets named `<prefix>`, `<prefix>.dbtype`,
        // `<prefix>.index`, `<prefix>_ss`, `<prefix>_ca`, etc. — so
        // we validate the *directory* containing the prefix exists
        // but don't require any specific extension to be present.
        // FoldSeek itself produces a clean error message at run time
        // if the database files are missing. (Same approach as the
        // BLAST adapter.)
        // Round-9 hardening: relative prefixes flow into the foldseek
        // command line; wrap with `confined_join`.
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
                    "[bio.foldseek].database `{}` directory `{}` not found",
                    input.database.display(),
                    database_dir.display()
                ),
            });
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "FoldSeek 8+ required; install via `brew install foldseek` \
                       or `conda install -c bioconda foldseek`"
                .into(),
        })?;

        // Pinned output filename and tmp dir, both derived from
        // `output_basename` so `collect()` can find the result file
        // and ignore the scratch directory by name. Both stay
        // workdir-relative — the subprocess runs inside the workdir
        // (see PreparedJob.workdir) so relative paths land in the
        // right place.
        let output_file = format!("{}.m8", input.output_basename);
        let tmp_dir = format!("tmp_{}", input.output_basename);

        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from("easy-search"),
            source_query.into_os_string(),
            source_database.into_os_string(),
            OsString::from(&output_file),
            OsString::from(&tmp_dir),
            OsString::from("--threads"),
            OsString::from(input.threads.to_string()),
        ];
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // A single-protein vs PDB100 search runs in seconds, but
            // searches against AlphaFold-DB (~200M structures) on a
            // single node take hours. 4 hours mirrors the rest of the
            // bio search adapters' generous default.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting FoldSeek", |line| {
            let mut hint = subprocess::Hint::default();
            // FoldSeek inherits MMseqs2's logging conventions — it
            // prints "Step N of M" markers and per-stage timing
            // summaries; the "Time for processing" line marks the
            // end of the run. (Same shape as the MMseqs2 adapter.)
            if line.contains("Time for processing") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("Step ") && line.contains(" of ") {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.contains("Error") || line.contains("ERROR") {
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
        // Re-derive the output filename and tmp directory name from the
        // prepared command so collect() can find the canonical result
        // file (and skip the scratch dir) without re-reading case.toml.
        // Layout in `native_command`:
        //   [bin, "easy-search", query, db, <output_file>, <tmp_dir>, ...]
        let output_file = job
            .native_command
            .get(4)
            .and_then(|s| s.to_str())
            .map(|s| s.to_string());
        let tmp_dir_name = job
            .native_command
            .get(5)
            .and_then(|s| s.to_str())
            .map(|s| s.to_string());

        // Provenance: hash the staged results file if present.
        // Falls back to case.toml when the search hasn't produced
        // results yet — keeps the provenance block well-formed for
        // partial / failed runs.
        let case_hash_input = match &output_file {
            Some(name) => {
                let p = job.workdir.join(name);
                if p.is_file() {
                    p
                } else {
                    job.workdir.join("case.toml")
                }
            }
            None => job.workdir.join("case.toml"),
        };
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "FoldSeek",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk the workdir top-level. Pick up the canonical
        // `<output_basename>.m8` results file plus any `*.log` files;
        // explicitly skip the `tmp_<output_basename>` scratch
        // directory — it's intermediate state that would clutter the
        // artifact list.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-foldseek", ?e, "workdir read failed");
                return Ok(results);
            }
        };
        for entry in entries.flatten() {
            let path = entry.path();
            // Skip the scratch directory entirely (it's the only
            // entry the adapter consciously hides).
            if path.is_dir() {
                continue;
            }
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

            let is_results = match (&file_name, &output_file) {
                (Some(name), Some(want)) => name == want,
                _ => false,
            };

            let (kind, label) = if is_results {
                (ArtifactKind::Tabular, "FoldSeek search results".to_string())
            } else if ext.as_deref() == Some("log") {
                (ArtifactKind::Log, "FoldSeek log".to_string())
            } else {
                continue;
            };
            artefacts.push(Artifact {
                path,
                kind,
                checksum: None,
                label,
            });
        }
        // The `tmp_dir_name` derivation above is intentional — it
        // documents in code that the scratch directory exists and
        // is deliberately excluded from artifacts. Drop it now that
        // collect() has recorded its presence implicitly via the
        // is_dir() skip above.
        let _ = tmp_dir_name;
        artefacts.sort_by(|a, b| a.path.cmp(&b.path));
        results.artifacts = artefacts;
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.foldseek.search"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = FoldseekAdapter::new().info();
        assert_eq!(info.id, "foldseek");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "GPL-3.0");
        assert_eq!(info.display_name, "FoldSeek");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = FoldseekAdapter::new().info();
        assert_eq!(info.version_range.min_inclusive, Version::new(8, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(10, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = FoldseekAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.foldseek.search"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = FoldseekAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    /// Round-9 RED→GREEN: `[bio.foldseek].database` relative form
    /// used bare `case.path.join`. Wrap with `confined_join`.
    #[test]
    fn prepare_rejects_database_traversing_outside_case_dir() {
        use valenx_test_utils::tempdir;
        let d = tempdir("foldseek-db-trav");
        std::fs::write(d.join("query.pdb"), b"HEADER PLACEHOLDER\n").unwrap();
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "foldseek.search"

[bio.foldseek]
query           = "query.pdb"
database        = "../../etc/passwd"
output_basename = "results"
"#,
        )
        .unwrap();
        let case = Case {
            id: "foldseek-db-trav".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = FoldseekAdapter::new().prepare(&case, &workdir).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("..") || msg.contains("stay within") || msg.contains("escape"),
            "expected confined_join rejection, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
