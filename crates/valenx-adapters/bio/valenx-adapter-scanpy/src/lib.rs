//! # valenx-adapter-scanpy
//!
//! Adapter for [Scanpy](https://scanpy.readthedocs.io/) — the
//! reference single-cell analysis library in Python (Wolf, Angerer &
//! Theis, 2018). Wraps the canonical `read → preprocess → PCA →
//! neighbours → UMAP → Leiden` recipe that powers most modern scRNA-seq
//! analyses.
//!
//! **Phase 19.5 — subprocess wrapper for user-provided scripts.** The
//! user supplies an `analyse.py` (or whatever filename) referenced
//! from `[bio.scanpy].script` in `case.toml` plus an input AnnData
//! `.h5ad`. `prepare()` stages the script + h5ad into the workdir and
//! `run()` invokes `python <script>` via the shared subprocess
//! runner. The script is responsible for reading
//! `valenx_params.json`, loading the input AnnData, running its
//! recipe, and writing the named output `.h5ad`.
//!
//! ## `valenx_params.json`
//!
//! Scanpy has no canonical CLI — every analysis script takes its own
//! knobs as call arguments. Rather than guess at a flag layout,
//! `prepare()` writes a flat JSON file `valenx_params.json` into the
//! workdir alongside the staged script and AnnData, containing the
//! parsed `[bio.scanpy]` knobs:
//!
//! ```json
//! {
//!   "input_h5ad":   "raw.h5ad",
//!   "output_h5ad":  "annotated.h5ad",
//!   "n_top_genes":  2000,
//!   "n_pcs":        50,
//!   "n_neighbors":  15,
//!   "resolution":   1.0
//! }
//! ```
//!
//! User scripts read it with `json.load(open("valenx_params.json"))`
//! and pass the values through to Scanpy themselves. This keeps the
//! adapter free of upstream API churn and means `case.toml` knobs
//! actually reach the recipe.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Duration;

use semver::Version;

use valenx_core::{
    adapter_helpers::{confined_join, find_on_path, live_provenance},
    error::RunPhase,
    subprocess, Adapter, AdapterError, AdapterInfo, Capabilities, Case, LicenseMode, Physics,
    PreparedJob, ProbeReport, RunContext, RunReport, VersionRange,
};
use valenx_fields::{
    artifact::{Artifact, ArtifactKind},
    Results,
};

use crate::case_input::ScanpyInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(ScanpyAdapter::new())
}

pub struct ScanpyAdapter;

impl ScanpyAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for ScanpyAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "scanpy";
/// Python interpreter candidates. `python3` first because on Linux
/// `python` may still be Python 2 on legacy distros; on Windows
/// `python` typically resolves to the Windows Store / 3.x install.
const PYTHON_BINARIES: &[&str] = &["python3", "python"];

impl Adapter for ScanpyAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "Scanpy",
            // Scanpy 1.10 (Apr 2024) is the first release with the
            // stable `sc.tl.leiden` flavor='igraph' default and the
            // typed AnnData 0.10 surface modern scripts target. Upper
            // bound 2.0 reserves room for an upcoming major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(1, 10, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "BSD-3-Clause",
            docs_url: "https://scanpy.readthedocs.io/",
            homepage_url: "https://scanpy.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(PYTHON_BINARIES) {
            Some(binary_path) => {
                // Try `import scanpy; print(scanpy.__version__)` —
                // confirms the `scanpy` package is importable from
                // the chosen interpreter (vs. just having Python on
                // PATH). Fall back to a "couldn't import" warning so
                // the probe still surfaces a useful state.
                let found_version = detect_scanpy_version(&binary_path);
                let mut warnings = Vec::new();
                if found_version.is_none() {
                    warnings.push(
                        "probe found `python` on PATH but could not import \
                         `scanpy` — install Scanpy with `pip install scanpy` \
                         (or `conda install -c conda-forge scanpy`) for runs \
                         to succeed"
                            .into(),
                    );
                }
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    warnings,
                    required_env: Vec::new(),
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "Python 3.9+ with Scanpy installed; \
                       `pip install scanpy` (or `conda install -c \
                       conda-forge scanpy`) after ensuring python3 \
                       is on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = ScanpyInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Stage the user-supplied Python script. Resolved against the
        // case directory; `confined_join` rejects absolute paths and
        // `..` traversal so the staged copy stays confined to the
        // case directory.
        let source_script = confined_join(&case.path, &input.script)?;
        if !source_script.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.scanpy].script `{}` not found (resolved {})",
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
                        "[bio.scanpy].script path `{}` has no filename",
                        input.script.display()
                    ),
                })?;
        let dest_script = workdir.join(script_filename);
        if source_script != dest_script {
            fs::copy(&source_script, &dest_script)?;
        }

        // Stage the input AnnData if a relative path was given.
        // Absolute paths stay where the user pointed — single-cell
        // matrices can be huge and copying them every run is wasteful.
        let input_h5ad_filename = if input.input_h5ad.is_absolute() {
            // Use the absolute path directly via the JSON params; the
            // user's script reads it verbatim from valenx_params.
            input.input_h5ad.display().to_string()
        } else {
            let source_h5ad = confined_join(&case.path, &input.input_h5ad)?;
            if !source_h5ad.is_file() {
                return Err(AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.scanpy].input_h5ad `{}` not found (resolved {})",
                        input.input_h5ad.display(),
                        source_h5ad.display()
                    ),
                });
            }
            let h5ad_name =
                input
                    .input_h5ad
                    .file_name()
                    .ok_or_else(|| AdapterError::InvalidCase {
                        case_path: case.path.join("case.toml"),
                        reason: format!(
                            "[bio.scanpy].input_h5ad path `{}` has no filename",
                            input.input_h5ad.display()
                        ),
                    })?;
            let dest_h5ad = workdir.join(h5ad_name);
            if source_h5ad != dest_h5ad {
                fs::copy(&source_h5ad, &dest_h5ad)?;
            }
            h5ad_name.to_string_lossy().into_owned()
        };

        // Drop a flat `valenx_params.json` into the workdir so the
        // user's analysis script can read the parsed `[bio.scanpy]`
        // knobs without reparsing case.toml. Built by hand to avoid
        // pulling in a serde_json dep for a flat object.
        let params_json = format!(
            "{{\n  \"input_h5ad\": {},\n  \"output_h5ad\": {},\n  \
             \"n_top_genes\": {},\n  \"n_pcs\": {},\n  \
             \"n_neighbors\": {},\n  \"resolution\": {}\n}}\n",
            json_string(&input_h5ad_filename),
            json_string(&input.output_h5ad),
            input.n_top_genes,
            input.n_pcs,
            input.n_neighbors,
            format_f64(input.resolution),
        );
        valenx_core::io_caps::atomic_write_str(&workdir.join("valenx_params.json"), &params_json)?;

        // Resolve the Python binary. Same logic as every other
        // Phase 17/17.5 Python-script adapter: bare `python` /
        // `python3` walks PATH; absolute paths or pinned interpreters
        // are honored verbatim.
        // Round-4 security: validate python interpreter spec
        // against the allow-list AND resolve to a real binary
        // in one step. Closes the arbitrary-binary-exec class
        // that round-3 only patched in 8 of the 48 affected
        // adapters.
        let binary_path = valenx_core::adapter_helpers::resolve_python_binary(
            &input.python,
            PYTHON_BINARIES,
        )
        // Round-5: do NOT rewrap as ToolNotInstalled — the resolver
        // returns InvalidCase for allow-list rejections (which a hint
        // string would have hidden) and a clear Other for PATH lookup
        // failures. Pass the error through unchanged.
        ?;

        let native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from(script_filename),
        ];

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Single-cell analyses scale with cell count: a few-thousand
            // cell PBMC tutorial finishes in seconds; a million-cell
            // atlas can run for hours. 2 hours is a generous default;
            // long runs override through their own progress reporting.
            estimated_runtime: Some(Duration::from_secs(2 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting Scanpy", |line| {
            let mut hint = subprocess::Hint::default();
            // Convention: the user-supplied script can emit a sentinel
            // line `[valenx] scanpy done` to signal completion before
            // exit; lift to a 95% progress tick.
            if line.contains("[valenx] scanpy done") {
                hint.progress = Some((95.0, line.to_string()));
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
        // Provenance: hash the staged Python script (the canonical
        // "this case is configured this way" input).
        let script_path = first_script_in_workdir(&job.workdir);
        let case_hash_input = script_path
            .clone()
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "Scanpy",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);

        // Walk the workdir top level and classify Scanpy outputs.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-scanpy", ?e, "workdir read failed");
                return Ok(results);
            }
        };
        let mut artefacts: Vec<Artifact> = Vec::new();
        for entry in entries.flatten() {
            let path = entry.path();
            if !path.is_file() {
                continue;
            }
            let ext = path
                .extension()
                .and_then(|s| s.to_str())
                .map(|s| s.to_ascii_lowercase());
            let (kind, label) = match ext.as_deref() {
                Some("h5ad") => (ArtifactKind::Native, "Scanpy AnnData output".to_string()),
                Some("png") | Some("pdf") => (ArtifactKind::Native, "Scanpy plot".to_string()),
                Some("csv") | Some("tsv") => (ArtifactKind::Tabular, "Scanpy table".to_string()),
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
            ribbon_contributions: vec!["bio.scanpy.analyse"],
        }
    }
}

/// Escape a string for embedding inside a JSON string literal.
/// Avoids pulling in a serde_json dependency for the trivially small
/// `valenx_params.json` we emit.
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

/// Format an f64 for inclusion in the params JSON. Always emits a
/// fractional digit so `1.0` round-trips as `1.0` rather than `1`,
/// keeping the user-side type stable across Python's `json.load`.
fn format_f64(v: f64) -> String {
    if v.is_finite() && v == v.trunc() {
        format!("{v:.1}")
    } else {
        format!("{v}")
    }
}

/// Lift the staged Python script out of the workdir for provenance
/// hashing. Returns the lexicographically-first `.py` file at the
/// top level, or `None` if none exists yet.
fn first_script_in_workdir(workdir: &Path) -> Option<PathBuf> {
    let entries = fs::read_dir(workdir).ok()?;
    let mut hits: Vec<PathBuf> = entries
        .flatten()
        .map(|e| e.path())
        .filter(|p| {
            p.extension()
                .and_then(|s| s.to_str())
                .map(|s| s.eq_ignore_ascii_case("py"))
                .unwrap_or(false)
        })
        .collect();
    hits.sort();
    hits.into_iter().next()
}

/// Run `python -c "import scanpy; print(scanpy.__version__)"` and
/// parse a `semver::Version` out of stdout. Returns `None` on any
/// failure (interpreter unusable, scanpy not importable, version
/// string malformed); `probe()` falls back to a "scanpy not
/// importable" warning in that case.
fn detect_scanpy_version(python_binary: &Path) -> Option<Version> {
    let output = std::process::Command::new(python_binary)
        .arg("-c")
        .arg("import scanpy; print(scanpy.__version__)")
        .stdin(std::process::Stdio::null())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .output()
        .ok()?;
    if !output.status.success() {
        return None;
    }
    let stdout = String::from_utf8(output.stdout).ok()?;
    let raw = valenx_core::adapter_helpers::extract_semver(&stdout)?;
    let dots = raw.chars().filter(|c| *c == '.').count();
    let normalised: String = match dots {
        0 => format!("{raw}.0.0"),
        1 => format!("{raw}.0"),
        _ => raw,
    };
    Version::parse(&normalised).ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = ScanpyAdapter::new().info();
        assert_eq!(info.id, "scanpy");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "BSD-3-Clause");
        assert_eq!(info.display_name, "Scanpy");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = ScanpyAdapter::new().info();
        assert_eq!(info.version_range.min_inclusive, Version::new(1, 10, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(2, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = ScanpyAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.scanpy.analyse"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = ScanpyAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
