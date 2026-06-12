//! # valenx-adapter-omegafold
//!
//! Adapter for [OmegaFold](https://github.com/HeliXonProtein/OmegaFold) —
//! HelixonAI's single-sequence protein structure prediction network.
//! Apache-2.0 licensed; sister to ESMFold / OpenFold / AlphaFold 2/3 /
//! RoseTTAFold from Phase 17.5–17.7. Unlike AlphaFold or ESMFold,
//! OmegaFold does **not** require an MSA — it works on a single query
//! sequence — which makes it the predictor of choice for orphan or
//! synthetic proteins where MSAs are weak or unavailable.
//!
//! **Phase 17.7 — wrapper around the bundled `omegafold` CLI.** Unlike
//! ESMFold or RoseTTAFold (which we drive through user-supplied predict
//! scripts because they have no canonical CLI), OmegaFold ships as an
//! installed Python package with its own command-line entry point:
//!
//! ```text
//! omegafold <fasta> <output_dir> [--model <model_dir>]
//! ```
//!
//! `prepare()` resolves the FASTA (passing the absolute path through
//! verbatim — large query batches don't need staging) and composes the
//! invocation. If the standalone `omegafold` binary isn't on PATH at
//! prepare time, we fall back to `<python> -m omegafold ...` so the
//! adapter still works when only the Python package is installed (e.g.
//! in a conda env where the entry-point script wasn't shimmed onto the
//! system PATH).
//!
//! On `collect()` we walk one level deep into the workdir's
//! `<output_basename>/` subdir for `*.pdb` (predicted structures) and
//! `*.json` (per-residue metadata). Top-level `*.log` files in the
//! workdir get surfaced as logs.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Duration;

use semver::Version;

use valenx_core::{
    adapter_helpers::{find_on_path, live_provenance},
    error::RunPhase,
    subprocess, Adapter, AdapterError, AdapterInfo, Capabilities, Case, LicenseMode, Physics,
    PreparedJob, ProbeReport, RunContext, RunReport, VersionRange,
};
use valenx_fields::{
    artifact::{Artifact, ArtifactKind},
    Results,
};

use crate::case_input::OmegaFoldInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(OmegaFoldAdapter::new())
}

pub struct OmegaFoldAdapter;

impl OmegaFoldAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for OmegaFoldAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "omegafold";
/// Probe binary candidates. `omegafold` first because that's the
/// canonical install-name; `python3` / `python` cover the fallback
/// path (`python -m omegafold ...`) when only the package is reachable.
const PROBE_BINARIES: &[&str] = &["omegafold", "python3", "python"];
/// Python interpreter candidates for the `python -m omegafold` fallback.
const PYTHON_BINARIES: &[&str] = &["python3", "python"];

impl Adapter for OmegaFoldAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "OmegaFold",
            // OmegaFold is on its 1.x line as of writing; the upstream
            // project hasn't tagged formal releases but the bundled
            // module exposes a 1.x version. Upper bound 2.0 reserves
            // room for an upcoming major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(1, 0, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "Apache-2.0",
            docs_url: "https://github.com/HeliXonProtein/OmegaFold",
            homepage_url: "https://github.com/HeliXonProtein/OmegaFold",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        // Prefer the standalone `omegafold` CLI; fall back to Python so
        // `python -m omegafold` runs still work. The combined PATH walk
        // matches the spec convention; the per-binary lookups below
        // tease apart the warning condition (omegafold missing but
        // Python present).
        match find_on_path(PROBE_BINARIES) {
            Some(binary_path) => {
                let omegafold_present = find_on_path(&["omegafold"]).is_some();
                let mut warnings = Vec::new();
                if !omegafold_present {
                    // Python is on PATH (we got *some* hit) but the
                    // dedicated CLI shim isn't — the `python -m
                    // omegafold` fallback will kick in at run time,
                    // warn the user so they know why their PATH lookup
                    // didn't surface the binary.
                    warnings.push(
                        "OmegaFold CLI not found on PATH; install via \
                         pip install git+https://github.com/HeliXonProtein/OmegaFold.git"
                            .into(),
                    );
                }
                Ok(ProbeReport {
                    ok: true,
                    found_version: None,
                    binary_path: Some(binary_path),
                    warnings,
                    required_env: Vec::new(),
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "OmegaFold CLI or Python 3.8+ with the OmegaFold \
                       package required; install via \
                       `pip install git+https://github.com/HeliXonProtein/OmegaFold.git` \
                       and ensure `omegafold` (or `python3`) is on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = OmegaFoldInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.omegafold].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Resolve the FASTA path against the case directory if relative.
        // We do **not** stage it into the workdir — OmegaFold reads it
        // by absolute path and large multi-sequence FASTAs don't need
        // duplicating. (This is the same convention as the BWA adapter
        // for its reference / read files.)
        let resolved_fasta = if input.fasta.is_absolute() {
            input.fasta.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.fasta)?
        };
        if !resolved_fasta.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.omegafold].fasta `{}` not found (resolved {})",
                    input.fasta.display(),
                    resolved_fasta.display()
                ),
            });
        }

        // Pin the output dir to a workdir-relative path with the
        // configured basename. OmegaFold's CLI takes the output dir as
        // its second positional; we let the runtime resolve it relative
        // to the workdir (which is the subprocess runner's cwd).
        let output_dir_arg = OsString::from(&input.output_basename);

        // Prefer the standalone `omegafold` binary; fall back to the
        // `python -m omegafold` invocation so the adapter still works
        // when only the package is installed (typical conda layout).
        let mut native_command: Vec<OsString> = if let Some(bin) = find_on_path(&["omegafold"]) {
            vec![
                bin.into_os_string(),
                resolved_fasta.into_os_string(),
                output_dir_arg,
            ]
        } else {
            // Resolve the Python binary. Same logic as every other
            // Phase 17 Python-script adapter: bare `python` / `python3`
            // walks PATH; absolute paths or pinned interpreters are
            // honored verbatim.
            // Round-4 security: validate python interpreter spec
            // against the allow-list AND resolve to a real binary
            // in one step. Closes the arbitrary-binary-exec class
            // that round-3 only patched in 8 of the 48 affected
            // adapters.
            let binary_path = valenx_core::adapter_helpers::resolve_python_binary(
                &input.python,
                PYTHON_BINARIES,
            )
            .map_err(|e| AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: format!("neither `omegafold` nor Python interpreter the configured value reachable on PATH: {e}"),
            })?;
            vec![
                binary_path.into_os_string(),
                OsString::from("-m"),
                OsString::from("omegafold"),
                resolved_fasta.into_os_string(),
                output_dir_arg,
            ]
        };

        // Optional `--model <model_dir>` flag — only emitted when the
        // user explicitly configured a checkpoint directory.
        if let Some(model_dir) = input.model_dir.as_ref() {
            native_command.push(OsString::from("--model"));
            native_command.push(model_dir.clone().into_os_string());
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // OmegaFold runtime is dominated by the model forward pass
            // and scales with sequence length; small targets finish in
            // seconds on a GPU, large multi-sequence FASTAs run for
            // hours. 2 hours is a generous default; long runs override
            // through their own progress reporting.
            estimated_runtime: Some(Duration::from_secs(2 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting OmegaFold", |line| {
            let mut hint = subprocess::Hint::default();
            // OmegaFold prints `Predicting <name>` lines per FASTA
            // entry; lift them to a 50% progress tick so the UI shows
            // forward motion. Tracebacks / errors propagate as
            // warnings so the run report surfaces them next to the
            // exit code.
            if line.contains("Predicting") {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.contains("Traceback") || line.contains("Error") {
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
        // Provenance: hash the configured output dir if present (the
        // canonical run output). Falls back to case.toml when the
        // prediction hasn't produced any files yet — keeps the
        // provenance block well-formed for partial / failed runs.
        let output_dir = first_subdir(&job.workdir);
        let case_hash_input = output_dir
            .clone()
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "OmegaFold",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk one level deep into the output basename subdir for
        // PDBs and metadata JSONs. We don't know the user's configured
        // basename without reparsing case.toml — every directory inside
        // the workdir is a candidate. Top-level workdir files (e.g.
        // `omegafold.log`) get classified as Log artifacts directly.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-omegafold", ?e, "workdir read failed");
                return Ok(results);
            }
        };
        let mut subdirs: Vec<PathBuf> = Vec::new();
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                subdirs.push(path);
                continue;
            }
            if !path.is_file() {
                continue;
            }
            let ext = path
                .extension()
                .and_then(|s| s.to_str())
                .map(|s| s.to_ascii_lowercase());
            if let Some("log") = ext.as_deref() {
                artefacts.push(Artifact {
                    path,
                    kind: ArtifactKind::Log,
                    checksum: None,
                    label: "OmegaFold log".to_string(),
                });
            }
        }
        subdirs.sort();
        for subdir in subdirs {
            let inner = match fs::read_dir(&subdir) {
                Ok(e) => e,
                Err(e) => {
                    tracing::warn!(
                        target: "valenx-omegafold",
                        ?e,
                        subdir = %subdir.display(),
                        "subdir read failed"
                    );
                    continue;
                }
            };
            let mut pdb_paths: Vec<PathBuf> = Vec::new();
            let mut json_paths: Vec<PathBuf> = Vec::new();
            for entry in inner.flatten() {
                let path = entry.path();
                if !path.is_file() {
                    continue;
                }
                let ext = path
                    .extension()
                    .and_then(|s| s.to_str())
                    .map(|s| s.to_ascii_lowercase());
                match ext.as_deref() {
                    Some("pdb") => pdb_paths.push(path),
                    Some("json") => json_paths.push(path),
                    _ => continue,
                }
            }
            pdb_paths.sort();
            json_paths.sort();
            for path in pdb_paths {
                artefacts.push(Artifact {
                    path,
                    kind: ArtifactKind::Native,
                    checksum: None,
                    label: "OmegaFold predicted structure".to_string(),
                });
            }
            for path in json_paths {
                artefacts.push(Artifact {
                    path,
                    kind: ArtifactKind::Log,
                    checksum: None,
                    label: "OmegaFold metadata".to_string(),
                });
            }
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
            ribbon_contributions: vec!["bio.omegafold.predict"],
        }
    }
}

/// Lift the lexicographically-first directory inside the workdir for
/// provenance hashing. Returns `None` when the workdir contains no
/// subdirectories (e.g. before the OmegaFold CLI has produced output).
fn first_subdir(workdir: &Path) -> Option<PathBuf> {
    let entries = fs::read_dir(workdir).ok()?;
    let mut hits: Vec<PathBuf> = entries
        .flatten()
        .map(|e| e.path())
        .filter(|p| p.is_dir())
        .collect();
    hits.sort();
    hits.into_iter().next()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = OmegaFoldAdapter::new().info();
        assert_eq!(info.id, "omegafold");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "Apache-2.0");
        assert_eq!(info.display_name, "OmegaFold");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = OmegaFoldAdapter::new().info();
        assert_eq!(info.version_range.min_inclusive, Version::new(1, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(2, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = OmegaFoldAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.omegafold.predict"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = OmegaFoldAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
