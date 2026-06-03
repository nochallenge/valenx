//! # valenx-adapter-forecast
//!
//! Adapter for [FORECasT](https://partslab.sanger.ac.uk/FORECasT) —
//! Felicity Allen / Wellcome Sanger's CRISPR/Cas9 indel pattern
//! predictor (Allen et al., Nature Biotechnology). FORECasT is the
//! Sanger sister to inDelphi: given a target window and a guide RNA it
//! predicts the indel-frequency distribution after Cas9 cleavage and
//! NHEJ repair, ranking insertions and deletions by predicted
//! likelihood.
//!
//! **Phase 35.6 — Python subprocess wrapper for user-provided
//! scripts.** FORECasT is published under the SelfTarget name; the
//! Python module is `selftarget` (github.com/felicityallen/SelfTarget).
//! The user authors a `predict.py` referenced from
//! `[bio.forecast].script` in `case.toml` that does
//! `import selftarget` and the actual prediction logic. `prepare()`
//! stages the script (and an optional input `.fa` template) into the
//! workdir, drops a flat `valenx_params.json` next to it, and `run()`
//! invokes `python <script>` via the shared subprocess runner.
//!
//! ## `valenx_params.json`
//!
//! ```json
//! {
//!   "output_basename": "output",
//!   "input_fasta": "target.fa"
//! }
//! ```
//!
//! `input_fasta` is omitted entirely (not `null`) when the user did
//! not supply one.
//!
//! On `collect()` we walk the workdir for `<basename>*.csv` (predicted
//! indel tables), `<basename>*.txt` (summary), plus any `*.log` files
//! Python emits.

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

use crate::case_input::ForecastInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(ForecastAdapter::new())
}

pub struct ForecastAdapter;

impl ForecastAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for ForecastAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "forecast";
const PYTHON_BINARIES: &[&str] = &["python3", "python"];

impl Adapter for ForecastAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "FORECasT",
            // FORECasT / SelfTarget 1.x is the upstream line; 2.0
            // reserves room for the next major.
            version_range: VersionRange {
                min_inclusive: Version::new(1, 0, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "Apache-2.0",
            docs_url: "https://github.com/felicityallen/SelfTarget",
            homepage_url: "https://partslab.sanger.ac.uk/FORECasT",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(PYTHON_BINARIES) {
            Some(binary_path) => {
                let import_ok = selftarget_importable(&binary_path);
                let mut warnings = Vec::new();
                if !import_ok {
                    warnings.push(
                        "probe found `python` on PATH but could not import \
                         `selftarget` — install from upstream \
                         https://github.com/felicityallen/SelfTarget for \
                         runs to succeed (FORECasT is published as the \
                         `selftarget` Python package)"
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
                hint: "Python 3.8+ with selftarget installed; install \
                       from upstream \
                       https://github.com/felicityallen/SelfTarget after \
                       ensuring python3 is on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = ForecastInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.forecast].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        let source_script = confined_join(&case.path, &input.script)?;
        if !source_script.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.forecast].script `{}` not found (resolved {})",
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
                        "[bio.forecast].script path `{}` has no filename",
                        input.script.display()
                    ),
                })?;
        let dest_script = workdir.join(script_filename);
        if source_script != dest_script {
            fs::copy(&source_script, &dest_script)?;
        }

        let staged_input_fasta: Option<String> = match input.input_fasta.as_ref() {
            Some(fa_path) => {
                let source_fa = confined_join(&case.path, fa_path)?;
                if !source_fa.is_file() {
                    return Err(AdapterError::InvalidCase {
                        case_path: case.path.join("case.toml"),
                        reason: format!(
                            "[bio.forecast].input_fasta `{}` not found (resolved {})",
                            fa_path.display(),
                            source_fa.display()
                        ),
                    });
                }
                let fa_name = fa_path
                    .file_name()
                    .ok_or_else(|| AdapterError::InvalidCase {
                        case_path: case.path.join("case.toml"),
                        reason: format!(
                            "[bio.forecast].input_fasta path `{}` has no filename",
                            fa_path.display()
                        ),
                    })?;
                let dest_fa = workdir.join(fa_name);
                if source_fa != dest_fa {
                    fs::copy(&source_fa, &dest_fa)?;
                }
                Some(fa_name.to_string_lossy().to_string())
            }
            None => None,
        };

        let mut params = String::new();
        params.push_str("{\n");
        params.push_str("  \"output_basename\": ");
        params.push_str(&json_string(&input.output_basename));
        if let Some(name) = staged_input_fasta.as_deref() {
            params.push_str(",\n  \"input_fasta\": ");
            params.push_str(&json_string(name));
        }
        params.push_str("\n}\n");
        valenx_core::io_caps::atomic_write_str(&workdir.join("valenx_params.json"), &params)?;

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
            // FORECasT prediction is CPU-light: a logistic-regression
            // family scoring per-window indel distributions; 30
            // minutes covers batch jobs over many target sites.
            estimated_runtime: Some(Duration::from_secs(30 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting FORECasT", |line| {
            let mut hint = subprocess::Hint::default();
            if line.contains("[valenx] FORECasT done") {
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
        let script_path = first_script_in_workdir(&job.workdir);
        let case_hash_input = script_path
            .clone()
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "FORECasT",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        let basename = read_output_basename(&job.workdir);

        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-forecast", ?e, "workdir read failed");
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
                Some("csv") => {
                    if !stem_matches_basename {
                        continue;
                    }
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Tabular,
                        checksum: None,
                        label: "FORECasT indel predictions".to_string(),
                    });
                }
                Some("txt") => {
                    if !stem_matches_basename {
                        continue;
                    }
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Tabular,
                        checksum: None,
                        label: "FORECasT summary".to_string(),
                    });
                }
                Some("log") => {
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Log,
                        checksum: None,
                        label: "FORECasT log".to_string(),
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
            ribbon_contributions: vec!["bio.forecast.predict"],
        }
    }
}

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

fn read_output_basename(workdir: &Path) -> Option<String> {
    let text = valenx_core::io_caps::read_capped_to_string(
        &workdir.join("valenx_params.json"),
        valenx_core::io_caps::MAX_ADAPTER_PARAMS_BYTES as usize,
    )
    .ok()?;
    extract_json_string(&text, "output_basename")
}

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

fn selftarget_importable(python_binary: &Path) -> bool {
    let output = std::process::Command::new(python_binary)
        .arg("-c")
        .arg("import selftarget")
        .stdin(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .output();
    match output {
        Ok(o) => o.status.success(),
        Err(_) => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = ForecastAdapter::new().info();
        assert_eq!(info.id, "forecast");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "Apache-2.0");
        assert_eq!(info.display_name, "FORECasT");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = ForecastAdapter::new().info();
        assert_eq!(info.version_range.min_inclusive, Version::new(1, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(2, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = ForecastAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.forecast.predict"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = ForecastAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
