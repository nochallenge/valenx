//! # valenx-adapter-scvi
//!
//! Adapter for [scvi-tools](https://scvi-tools.org/) — the
//! probabilistic-modelling toolbox for single-cell omics (Lopez et
//! al., 2018; Gayoso et al., 2022). Wraps the canonical SCVI / SCANVI
//! / TOTALVI / LinearSCVI training loop: load AnnData, set up the
//! model, fit, write a denoised latent representation back into
//! `obsm` and persist the trained model checkpoint.
//!
//! **Phase 19.5 — subprocess wrapper for user-provided scripts.** The
//! user supplies a `train.py` (or whatever filename) referenced from
//! `[bio.scvi].script` in `case.toml` plus an input AnnData `.h5ad`.
//! `prepare()` stages the script + h5ad into the workdir and `run()`
//! invokes `python <script>` via the shared subprocess runner. The
//! script is responsible for reading `valenx_params.json`, loading
//! the AnnData, running its training recipe, and writing the named
//! output `.h5ad` plus any model checkpoints.
//!
//! ## `valenx_params.json`
//!
//! scvi-tools' Python API evolves between minor releases — pin the
//! adapter to a flat JSON contract instead of a CLI flag layout. The
//! file we drop into the workdir contains:
//!
//! ```json
//! {
//!   "input_h5ad":  "raw.h5ad",
//!   "output_h5ad": "with_latent.h5ad",
//!   "model":       "scvi",
//!   "n_latent":    10,
//!   "n_hidden":    128,
//!   "n_layers":    2,
//!   "max_epochs":  400
//! }
//! ```
//!
//! User scripts read it with `json.load(open("valenx_params.json"))`
//! and pass the values through to scvi-tools themselves. This keeps
//! the adapter free of upstream API churn and means `case.toml`
//! knobs actually reach the model.

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

use crate::case_input::ScviInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(ScviAdapter::new())
}

pub struct ScviAdapter;

impl ScviAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for ScviAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "scvi";
/// Python interpreter candidates. `python3` first because on Linux
/// `python` may still be Python 2 on legacy distros; on Windows
/// `python` typically resolves to the Windows Store / 3.x install.
const PYTHON_BINARIES: &[&str] = &["python3", "python"];

impl Adapter for ScviAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "scvi-tools",
            // scvi-tools 1.1 (Jan 2024) is the first release with the
            // stable `setup_anndata` typed surface and Lightning 2.x
            // training loop. Upper bound 2.0 reserves room for an
            // upcoming major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(1, 1, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "BSD-3-Clause",
            docs_url: "https://docs.scvi-tools.org/",
            homepage_url: "https://scvi-tools.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(PYTHON_BINARIES) {
            Some(binary_path) => {
                // Try `import scvi; print(scvi.__version__)` —
                // confirms the `scvi-tools` PyPI package (which
                // imports as `scvi`) is importable. Fall back to a
                // "couldn't import" warning so the probe still
                // surfaces a useful state.
                let found_version = detect_scvi_version(&binary_path);
                let mut warnings = Vec::new();
                if found_version.is_none() {
                    warnings.push(
                        "probe found `python` on PATH but could not import \
                         `scvi` — install scvi-tools with \
                         `pip install scvi-tools` (or `conda install -c \
                         conda-forge scvi-tools`) for runs to succeed"
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
                hint: "Python 3.9+ with scvi-tools installed; \
                       `pip install scvi-tools` (or `conda install -c \
                       conda-forge scvi-tools`) after ensuring python3 \
                       is on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = ScviInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Stage the user-supplied Python script.
        // `confined_join` rejects absolute paths and `..` traversal so
        // the staged copy stays confined to the case directory.
        let source_script = confined_join(&case.path, &input.script)?;
        if !source_script.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.scvi].script `{}` not found (resolved {})",
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
                        "[bio.scvi].script path `{}` has no filename",
                        input.script.display()
                    ),
                })?;
        let dest_script = workdir.join(script_filename);
        if source_script != dest_script {
            fs::copy(&source_script, &dest_script)?;
        }

        // Stage the input AnnData if a relative path was given.
        // Absolute paths stay where the user pointed — scvi training
        // matrices can be tens of GB and copying them is wasteful.
        let input_h5ad_filename = if input.input_h5ad.is_absolute() {
            input.input_h5ad.display().to_string()
        } else {
            let source_h5ad = confined_join(&case.path, &input.input_h5ad)?;
            if !source_h5ad.is_file() {
                return Err(AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.scvi].input_h5ad `{}` not found (resolved {})",
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
                            "[bio.scvi].input_h5ad path `{}` has no filename",
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
        // user's training script can read the parsed `[bio.scvi]`
        // knobs. Built by hand to avoid pulling in serde_json.
        let params_json = format!(
            "{{\n  \"input_h5ad\": {},\n  \"output_h5ad\": {},\n  \
             \"model\": {},\n  \"n_latent\": {},\n  \"n_hidden\": {},\n  \
             \"n_layers\": {},\n  \"max_epochs\": {}\n}}\n",
            json_string(&input_h5ad_filename),
            json_string(&input.output_h5ad),
            json_string(&input.model),
            input.n_latent,
            input.n_hidden,
            input.n_layers,
            input.max_epochs,
        );
        valenx_core::io_caps::atomic_write_str(&workdir.join("valenx_params.json"), &params_json)?;

        // Resolve the Python binary.
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
            // scvi-tools training scales with cell count, epoch budget,
            // and GPU. A few-thousand cell tutorial finishes in
            // minutes; an atlas-scale model can train for many hours.
            // 4 hours is a generous default; long runs override
            // through their own progress reporting.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting scvi-tools", |line| {
            let mut hint = subprocess::Hint::default();
            // Convention: the user-supplied script can emit a sentinel
            // line `[valenx] scvi done` to signal completion before
            // exit; lift to a 95% progress tick.
            if line.contains("[valenx] scvi done") {
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
        // Provenance: hash the staged Python script.
        let script_path = first_script_in_workdir(&job.workdir);
        let case_hash_input = script_path
            .clone()
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "scvi-tools",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);

        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-scvi", ?e, "workdir read failed");
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
                Some("h5ad") => (
                    ArtifactKind::Native,
                    "scvi-tools AnnData output".to_string(),
                ),
                Some("pt") | Some("pkl") => (
                    ArtifactKind::Native,
                    "scvi-tools model checkpoint".to_string(),
                ),
                Some("png") => (ArtifactKind::Native, "scvi-tools plot".to_string()),
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
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.scvi.train"],
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

/// Run `python -c "import scvi; print(scvi.__version__)"` and parse a
/// `semver::Version` out of stdout. Returns `None` on any failure
/// (interpreter unusable, scvi not importable, version string
/// malformed); `probe()` falls back to a "scvi not importable"
/// warning in that case.
fn detect_scvi_version(python_binary: &Path) -> Option<Version> {
    let output = std::process::Command::new(python_binary)
        .arg("-c")
        .arg("import scvi; print(scvi.__version__)")
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
        let info = ScviAdapter::new().info();
        assert_eq!(info.id, "scvi");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "BSD-3-Clause");
        assert_eq!(info.display_name, "scvi-tools");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = ScviAdapter::new().info();
        assert_eq!(info.version_range.min_inclusive, Version::new(1, 1, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(2, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = ScviAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.scvi.train"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = ScviAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
