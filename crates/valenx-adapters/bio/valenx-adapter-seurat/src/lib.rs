//! # valenx-adapter-seurat
//!
//! Adapter for [Seurat](https://satijalab.org/seurat/) — the
//! dominant R-based toolkit for single-cell genomics: QC,
//! normalisation, dimensionality reduction (PCA / UMAP),
//! clustering, differential expression, multi-modal integration.
//!
//! **Phase 19.6 — Rscript subprocess wrapper for user-provided R
//! scripts.** Seurat is an R package (not a CLI), so the adapter
//! itself doesn't generate R; the user authors an `analysis.R`
//! referenced from `[bio.seurat].script` in `case.toml` that does
//! `library(Seurat)` and the actual data work. `prepare()` stages
//! the script (and an optional input matrix) into the workdir,
//! drops a flat `valenx_params.json` next to it so the script can
//! read parsed knobs without re-parsing case.toml, and `run()`
//! invokes `Rscript <script>` via the shared subprocess runner.
//!
//! ## `valenx_params.json`
//!
//! ```json
//! {
//!   "output_basename": "analysis",
//!   "input_data": "matrix.h5"
//! }
//! ```
//!
//! `input_data` is omitted entirely (not `null`) when the user did
//! not supply one. Scripts read with
//! `jsonlite::fromJSON("valenx_params.json")` and resolve the
//! filename relative to the cwd (the workdir).
//!
//! On `collect()` we walk the workdir for `<basename>*.rds`
//! (Seurat objects), `<basename>*.csv` (tables), and
//! `<basename>*.png` (plots), plus any `*.log` files Rscript
//! emits.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Duration;

use semver::Version;

use valenx_core::{
    adapter_helpers::{confined_join, find_on_path, live_provenance, validate_rscript_binary},
    error::RunPhase,
    subprocess, Adapter, AdapterError, AdapterInfo, Capabilities, Case, LicenseMode, Physics,
    PreparedJob, ProbeReport, RunContext, RunReport, VersionRange,
};
use valenx_fields::{
    artifact::{Artifact, ArtifactKind},
    Results,
};

use crate::case_input::SeuratInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(SeuratAdapter::new())
}

pub struct SeuratAdapter;

impl SeuratAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for SeuratAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "seurat";
/// We probe `Rscript` (the non-interactive R driver). Seurat itself
/// is an R package — confirming it's installed requires running R,
/// which we deliberately don't do at probe time. The script's
/// `library(Seurat)` call surfaces a missing-package error at run
/// time with R's own (clear) diagnostic.
const RSCRIPT_BINARIES: &[&str] = &["Rscript"];

impl Adapter for SeuratAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "Seurat",
            // Seurat 4.0 (Apr 2021) introduced the v3-assay
            // architecture that the modern Seurat workflow leans on;
            // 5.x continues to extend it. Upper bound 6.0 reserves
            // room for the next major.
            version_range: VersionRange {
                min_inclusive: Version::new(4, 0, 0),
                max_exclusive: Version::new(6, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "MIT",
            docs_url: "https://satijalab.org/seurat/",
            homepage_url: "https://satijalab.org/seurat/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(RSCRIPT_BINARIES) {
            Some(binary_path) => Ok(ProbeReport {
                ok: true,
                // We don't run R at probe time — confirming Seurat
                // is installed requires `library(Seurat)`, which is
                // expensive and surfaces a clear error at run time
                // anyway. Leave found_version unknown.
                found_version: None,
                binary_path: Some(binary_path),
                warnings: Vec::new(),
                required_env: Vec::new(),
            }),
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "Seurat 4.0+ required; install R then \
                       `install.packages('Seurat')` from CRAN"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = SeuratInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.seurat].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Stage the user-supplied R script into the workdir.
        // `confined_join` rejects absolute paths and `..` traversal so
        // the staged copy stays confined to the case directory.
        let source_script = confined_join(&case.path, &input.script)?;
        if !source_script.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.seurat].script `{}` not found (resolved {})",
                    input.script.display(),
                    source_script.display()
                ),
            });
        }
        let script_filename =
            input
                .script
                .file_name()
                .ok_or_else(|| AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.seurat].script path `{}` has no filename",
                        input.script.display()
                    ),
                })?;
        let dest_script = workdir.join(script_filename);
        if source_script != dest_script {
            fs::copy(&source_script, &dest_script)?;
        }

        // Optionally stage the input data file so the script can
        // resolve it via a bare filename inside the workdir.
        let staged_input_data: Option<String> = match input.input_data.as_ref() {
            Some(data_path) => {
                let source_data = confined_join(&case.path, data_path)?;
                if !source_data.is_file() {
                    return Err(AdapterError::InvalidCase {
                        case_path: case.path.join("case.toml"),
                        reason: format!(
                            "[bio.seurat].input_data `{}` not found (resolved {})",
                            data_path.display(),
                            source_data.display()
                        ),
                    });
                }
                let data_filename =
                    data_path
                        .file_name()
                        .ok_or_else(|| AdapterError::InvalidCase {
                            case_path: case.path.join("case.toml"),
                            reason: format!(
                                "[bio.seurat].input_data path `{}` has no filename",
                                data_path.display()
                            ),
                        })?;
                let dest_data = workdir.join(data_filename);
                if source_data != dest_data {
                    fs::copy(&source_data, &dest_data)?;
                }
                Some(data_filename.to_string_lossy().to_string())
            }
            None => None,
        };

        // Drop a flat `valenx_params.json` into the workdir so the
        // user's R script can read the parsed `[bio.seurat]` knobs
        // without having to reparse case.toml itself. Built by hand
        // to avoid pulling in serde_json for a 2-key flat object.
        // When `input_data` is absent we omit the key entirely
        // (matching the spec — no `null` literal in the JSON).
        let mut params = String::new();
        params.push_str("{\n");
        params.push_str("  \"output_basename\": ");
        params.push_str(&json_string(&input.output_basename));
        if let Some(name) = staged_input_data.as_deref() {
            params.push_str(",\n  \"input_data\": ");
            params.push_str(&json_string(name));
        }
        params.push_str("\n}\n");
        valenx_core::io_caps::atomic_write_str(&workdir.join("valenx_params.json"), &params)?;

        // Round-5: validate the user-supplied rscript spec against the
        // allow-list BEFORE doing anything else. A hostile case.toml
        // setting `rscript = "/usr/bin/curl"` is rejected here as
        // InvalidCase (programmer-typed-a-bad-binary) rather than
        // executed silently as arbitrary code.
        let validated = validate_rscript_binary(&input.rscript).map_err(|e| {
            AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "Rscript interpreter rejected ('{}'). See \
                     `valenx_core::adapter_helpers::validate_rscript_binary` \
                     allow-list. Error: {e}",
                    input.rscript
                ),
            }
        })?;
        // Resolve the validated Rscript binary. Bare `Rscript` walks
        // PATH; absolute / relative allow-listed paths the user pinned
        // are honored verbatim if they exist, with a final PATH
        // fallback.
        let binary_path = if validated.is_absolute() && validated.is_file() {
            Some(validated)
        } else if input.rscript == "Rscript" {
            find_on_path(RSCRIPT_BINARIES)
        } else {
            find_on_path(&[input.rscript.as_str()]).or_else(|| find_on_path(RSCRIPT_BINARIES))
        }
        .ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: format!(
                "could not locate Rscript binary `{}` — install R and \
                 ensure Rscript is on PATH",
                input.rscript
            ),
        })?;

        let native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from(script_filename),
        ];

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Seurat workflows can run for a long time —
            // integration of millions of cells, label-transfer
            // across atlases, differential-expression sweeps. 60
            // minutes covers most interactive analyses; the user
            // can re-run for longer batches.
            estimated_runtime: Some(Duration::from_secs(60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting Seurat", |line| {
            let mut hint = subprocess::Hint::default();
            // Convention: the user-supplied script can emit a
            // sentinel line `[valenx] seurat done` to signal
            // completion before exit; lift to a 95% progress tick.
            if line.contains("[valenx] seurat done") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("Error in ") || line.contains("Error:") {
                // R surfaces user-visible errors as `Error in <fn>:`
                // and `Error: ...` — both worth a warning hint.
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
        // Provenance: hash the staged R script as the canonical
        // input descriptor. Falls back to case.toml when no script
        // is staged yet (partial / failed runs).
        let script_path = first_script_in_workdir(&job.workdir);
        let case_hash_input = script_path
            .clone()
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "Seurat",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Restrict typed outputs (.rds / .csv / .png) to those whose
        // stem starts with the configured `output_basename`. .log
        // files are accepted regardless — R's `sink()` / system
        // logs aren't typically prefixed.
        let basename = read_output_basename(&job.workdir);

        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-seurat", ?e, "workdir read failed");
                return Ok(results);
            }
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if !path.is_file() {
                continue;
            }
            let ext = path
                .extension()
                .and_then(|s| s.to_str())
                .map(|s| s.to_ascii_lowercase());
            let stem = path
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("")
                .to_string();
            let stem_matches_basename = match basename.as_deref() {
                Some(b) => stem.starts_with(b),
                None => true,
            };
            match ext.as_deref() {
                Some("rds") => {
                    if !stem_matches_basename {
                        continue;
                    }
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Native,
                        checksum: None,
                        label: "Seurat object (RDS)".to_string(),
                    });
                }
                Some("csv") => {
                    if !stem_matches_basename {
                        continue;
                    }
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Tabular,
                        checksum: None,
                        label: "Seurat output table".to_string(),
                    });
                }
                Some("png") => {
                    if !stem_matches_basename {
                        continue;
                    }
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Native,
                        checksum: None,
                        label: "Seurat plot".to_string(),
                    });
                }
                Some("log") => {
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Log,
                        checksum: None,
                        label: "Seurat log".to_string(),
                    });
                }
                _ => continue,
            }
        }
        artefacts.sort_by(|a, b| a.path.cmp(&b.path));
        results.artifacts = artefacts;
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.seurat.analyze"],
        }
    }
}

/// Escape a string for embedding inside a JSON string literal.
fn json_string(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 2);
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{0008}' => out.push_str("\\b"),
            '\u{000C}' => out.push_str("\\f"),
            c if (c as u32) < 0x20 => {
                out.push_str(&format!("\\u{:04x}", c as u32));
            }
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

/// Lift the staged R script out of the workdir for provenance
/// hashing. Returns the lexicographically-first `.R` file (case-
/// insensitive) at the top level, or `None` if none exists yet.
fn first_script_in_workdir(workdir: &Path) -> Option<PathBuf> {
    let entries = fs::read_dir(workdir).ok()?;
    let mut hits: Vec<PathBuf> = entries
        .flatten()
        .map(|e| e.path())
        .filter(|p| {
            p.extension()
                .and_then(|s| s.to_str())
                .map(|s| s.eq_ignore_ascii_case("R"))
                .unwrap_or(false)
        })
        .collect();
    hits.sort();
    hits.into_iter().next()
}

/// Pull `output_basename` out of our own hand-emitted
/// `valenx_params.json` for collect()-time output filtering.
fn read_output_basename(workdir: &Path) -> Option<String> {
    let text = valenx_core::io_caps::read_capped_to_string(
        &workdir.join("valenx_params.json"),
        valenx_core::io_caps::MAX_ADAPTER_PARAMS_BYTES as usize,
    )
    .ok()?;
    extract_json_string(&text, "output_basename")
}

/// Pull a flat string field out of our own hand-emitted
/// `valenx_params.json`. We wrote the file ourselves so we know
/// its shape; a full JSON parser would be overkill.
fn extract_json_string(text: &str, key: &str) -> Option<String> {
    let needle = format!("\"{key}\":");
    let idx = text.find(&needle)?;
    let rest = &text[idx + needle.len()..];
    let start = rest.find('"')? + 1;
    let body = &rest[start..];
    let mut out = String::new();
    let mut chars = body.chars();
    while let Some(c) = chars.next() {
        match c {
            '"' => return Some(out),
            '\\' => match chars.next()? {
                '"' => out.push('"'),
                '\\' => out.push('\\'),
                'n' => out.push('\n'),
                'r' => out.push('\r'),
                't' => out.push('\t'),
                'b' => out.push('\u{0008}'),
                'f' => out.push('\u{000C}'),
                other => out.push(other),
            },
            c => out.push(c),
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = SeuratAdapter::new().info();
        assert_eq!(info.id, "seurat");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "MIT");
        assert_eq!(info.display_name, "Seurat");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = SeuratAdapter::new().info();
        // Seurat 4.0 (Apr 2021) introduced the v3-assay architecture
        // that the modern workflow leans on; 5.x continues to extend
        // it. Upper bound 6.0 reserves room for the next major.
        assert_eq!(info.version_range.min_inclusive, Version::new(4, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(6, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = SeuratAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.seurat.analyze"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = SeuratAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
