//! # valenx-adapter-esmc
//!
//! Adapter for [ESM Cambrian (ESMC)](https://github.com/evolutionaryscale/esm)
//! — EvolutionaryScale's open-weight protein representation model. Where
//! ESM3 is a generative multi-modal system and ESMFold targets structure
//! prediction, ESMC is the workhorse: it produces high-quality
//! per-residue (or pooled) embeddings that downstream classifiers and
//! regressors consume directly. Two open checkpoints are exposed:
//!
//! - `esmc-300m` — small / fast, fits on a consumer GPU
//! - `esmc-600m` — larger / better representations
//!
//! **Phase 27.6 — subprocess wrapper for user-provided scripts.** The
//! user supplies an `embed_esmc.py` (or whatever filename) referenced
//! from `[bio.esmc].script` in `case.toml` plus an input FASTA.
//! `prepare()` stages the script + FASTA into the workdir and `run()`
//! invokes `python <script>` via the shared subprocess runner. The
//! script is responsible for invoking ESMC with the parsed knobs and
//! writing embeddings under `<output_basename>.{npy|npz|parquet}`.
//!
//! ## `valenx_params.json`
//!
//! ESMC's API surface evolves with each `esm` release. Rather than
//! guess at a fixed call shape, `prepare()` writes a flat JSON file
//! `valenx_params.json` into the workdir containing the parsed
//! `[bio.esmc]` knobs:
//!
//! ```json
//! {
//!   "input_fasta":     "query.fasta",
//!   "model_variant":   "esmc-300m",
//!   "pooling":         "per-residue",
//!   "output_basename": "embeddings"
//! }
//! ```
//!
//! User scripts read it with `json.load(open("valenx_params.json"))`
//! and pass the values through to ESMC themselves.
//!
//! On `collect()` we walk the workdir for `<output_basename>.{npy,
//! npz, parquet}` and surface each as a `Tabular` artifact with the
//! "ESMC embeddings" label.

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

use crate::case_input::EsmcInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(EsmcAdapter::new())
}

pub struct EsmcAdapter;

impl EsmcAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for EsmcAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "esmc";
/// Python interpreter candidates. `python3` first because on Linux
/// `python` may still be Python 2 on legacy distros; on Windows
/// `python` typically resolves to the Windows Store / 3.x install.
const PYTHON_BINARIES: &[&str] = &["python3", "python"];

impl Adapter for EsmcAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "ESM Cambrian",
            // ESMC ships in EvolutionaryScale's `esm` package; its
            // Cambrian-line releases are tagged as 1.x. Upper bound
            // 2.0 reserves room for an upcoming major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(1, 0, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "Cambrian-Open-License",
            docs_url: "https://github.com/evolutionaryscale/esm",
            homepage_url: "https://www.evolutionaryscale.ai/blog/esm-cambrian",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(PYTHON_BINARIES) {
            Some(binary_path) => {
                // Confirm the `esm` package is importable from the
                // chosen interpreter (vs. just having Python on PATH).
                let found_version = detect_esm_version(&binary_path);
                let mut warnings = Vec::new();
                if found_version.is_none() {
                    warnings.push(
                        "probe found `python` on PATH but could not import \
                         `esm` — install EvolutionaryScale's ESM package \
                         from https://github.com/evolutionaryscale/esm and \
                         ensure it's importable from the chosen interpreter \
                         for runs to succeed"
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
                hint: "Python 3.10+ with EvolutionaryScale's `esm` package \
                       installed; install from \
                       https://github.com/evolutionaryscale/esm after \
                       ensuring python3 is on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = EsmcInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.esmc].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Stage the user-supplied Python script. `confined_join`
        // rejects absolute paths and `..` traversal so the staged copy
        // stays confined to the case directory.
        let source_script = confined_join(&case.path, &input.script)?;
        if !source_script.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.esmc].script `{}` not found (resolved {})",
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
                        "[bio.esmc].script path `{}` has no filename",
                        input.script.display()
                    ),
                })?;
        let dest_script = workdir.join(script_filename);
        if source_script != dest_script {
            fs::copy(&source_script, &dest_script)?;
        }

        // Stage the input FASTA alongside the script.
        let source_fasta = confined_join(&case.path, &input.input_fasta)?;
        if !source_fasta.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.esmc].input_fasta `{}` not found (resolved {})",
                    input.input_fasta.display(),
                    source_fasta.display()
                ),
            });
        }
        let fasta_filename =
            input
                .input_fasta
                .file_name()
                .ok_or_else(|| AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.esmc].input_fasta path `{}` has no filename",
                        input.input_fasta.display()
                    ),
                })?;
        let dest_fasta = workdir.join(fasta_filename);
        if source_fasta != dest_fasta {
            fs::copy(&source_fasta, &dest_fasta)?;
        }

        // Drop a flat `valenx_params.json` into the workdir so the
        // user's embed script can read the parsed `[bio.esmc]` knobs
        // without having to reparse case.toml itself. Built by hand
        // to avoid pulling in a serde_json dep.
        let params_json = format!(
            "{{\n  \"input_fasta\": {},\n  \"model_variant\": {},\n  \"pooling\": {},\n  \"output_basename\": {}\n}}\n",
            json_string(&fasta_filename.to_string_lossy()),
            json_string(&input.model_variant),
            json_string(&input.pooling),
            json_string(&input.output_basename),
        );
        valenx_core::io_caps::atomic_write_str(&workdir.join("valenx_params.json"), &params_json)?;

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
            // ESMC inference is much faster than full structure
            // prediction; embedding a few hundred sequences typically
            // finishes in minutes on a modest GPU. 1 hour is a
            // generous default.
            estimated_runtime: Some(Duration::from_secs(60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting ESMC", |line| {
            let mut hint = subprocess::Hint::default();
            // Convention: the user-supplied script can emit a sentinel
            // line `[valenx] esmc done` to signal completion before
            // exit; lift to a 95% progress tick.
            if line.contains("[valenx] esmc done") {
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
        // Provenance: hash the staged FASTA (the canonical "this case
        // is configured this way" input). Falls back to the staged
        // script, then case.toml, when the FASTA isn't present yet.
        let fasta_path = first_fasta_in_workdir(&job.workdir);
        let script_path = first_script_in_workdir(&job.workdir);
        let case_hash_input = fasta_path
            .clone()
            .or_else(|| script_path.clone())
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "ESM Cambrian",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        let basename = read_params(&job.workdir);

        if let Some(p) = fasta_path {
            artefacts.push(Artifact {
                path: p,
                kind: ArtifactKind::Other,
                checksum: None,
                label: "ESMC input FASTA".to_string(),
            });
        }
        if let Some(p) = script_path {
            artefacts.push(Artifact {
                path: p,
                kind: ArtifactKind::Other,
                checksum: None,
                label: "ESMC script".to_string(),
            });
        }

        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-esmc", ?e, "workdir read failed");
                return Ok(results);
            }
        };
        let mut embedding_paths: Vec<PathBuf> = Vec::new();
        for entry in entries.flatten() {
            let path = entry.path();
            if !path.is_file() {
                continue;
            }
            let ext = path
                .extension()
                .and_then(|s| s.to_str())
                .map(|s| s.to_ascii_lowercase());
            match ext.as_deref() {
                Some("npy") | Some("npz") | Some("parquet") => {
                    let stem_ok = match basename.as_ref() {
                        Some(b) => {
                            let stem = path.file_stem().and_then(|s| s.to_str()).unwrap_or("");
                            stem == b.as_str()
                        }
                        None => true,
                    };
                    if stem_ok {
                        embedding_paths.push(path);
                    }
                }
                Some("log") => artefacts.push(Artifact {
                    path,
                    kind: ArtifactKind::Log,
                    checksum: None,
                    label: "ESMC log".to_string(),
                }),
                _ => continue,
            }
        }
        embedding_paths.sort();
        for path in embedding_paths {
            artefacts.push(Artifact {
                path,
                kind: ArtifactKind::Tabular,
                checksum: None,
                label: "ESMC embeddings".to_string(),
            });
        }

        artefacts.sort_by(|a, b| a.path.cmp(&b.path));
        results.artifacts = artefacts;
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.esmc.embed"],
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

/// Lift the staged FASTA out of the workdir for provenance hashing.
fn first_fasta_in_workdir(workdir: &Path) -> Option<PathBuf> {
    let entries = fs::read_dir(workdir).ok()?;
    let mut hits: Vec<PathBuf> = entries
        .flatten()
        .map(|e| e.path())
        .filter(|p| {
            p.extension()
                .and_then(|s| s.to_str())
                .map(|s| {
                    let s = s.to_ascii_lowercase();
                    matches!(s.as_str(), "fasta" | "fa" | "faa" | "fna")
                })
                .unwrap_or(false)
        })
        .collect();
    hits.sort();
    hits.into_iter().next()
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

/// Read just the `output_basename` field from `valenx_params.json` so
/// `collect()` can restrict embedding pickup to the configured stem.
fn read_params(workdir: &Path) -> Option<String> {
    let text = valenx_core::io_caps::read_capped_to_string(
        &workdir.join("valenx_params.json"),
        valenx_core::io_caps::MAX_ADAPTER_PARAMS_BYTES as usize,
    )
    .ok()?;
    extract_json_string(&text, "output_basename")
}

/// Pull a flat string field out of our own hand-emitted
/// `valenx_params.json`.
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

/// Run `python -c "import esm; print(esm.__version__)"` and parse a
/// `semver::Version` out of stdout.
fn detect_esm_version(python_binary: &Path) -> Option<Version> {
    let output = std::process::Command::new(python_binary)
        .arg("-c")
        .arg("import esm; print(esm.__version__)")
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
        let info = EsmcAdapter::new().info();
        assert_eq!(info.id, "esmc");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "Cambrian-Open-License");
        assert_eq!(info.display_name, "ESM Cambrian");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = EsmcAdapter::new().info();
        assert_eq!(info.version_range.min_inclusive, Version::new(1, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(2, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = EsmcAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.esmc.embed"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = EsmcAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
