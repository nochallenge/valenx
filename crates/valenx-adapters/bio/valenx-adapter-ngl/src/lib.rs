//! # valenx-adapter-ngl
//!
//! Adapter for [NGL Viewer](https://nglviewer.org/) — the Rose lab's
//! WebGL framework for molecular visualization (MIT licensed). NGL is
//! a JavaScript browser library that's widely embedded in Jupyter
//! notebooks via the `nglview` Python binding. This adapter wraps
//! `nglview` so a user-authored `viewer.py` can build an NGL view,
//! export a self-contained HTML page, render a still PNG, or dump
//! the underlying state JSON, all from the standard subprocess
//! pattern.
//!
//! **Phase 42 — Python subprocess wrapper for user-provided
//! scripts.** NGL (via the `nglview` Python binding) is a library, so
//! the adapter itself doesn't generate Python; the user authors a
//! `viewer.py` referenced from `[bio.ngl].script` in `case.toml`
//! that does `import nglview` and the actual viewer / render logic.
//! `prepare()` stages the script (and an optional input `.pdb` /
//! `.cif` structure) into the workdir, drops a flat
//! `valenx_params.json` next to it so the script can read parsed
//! knobs without re-parsing case.toml, and `run()` invokes
//! `python <script>` via the shared subprocess runner.
//!
//! ## `valenx_params.json`
//!
//! ```json
//! {
//!   "output_basename": "viewer",
//!   "input_structure": "structure.pdb"
//! }
//! ```
//!
//! `input_structure` is omitted entirely (not `null`) when the user
//! did not supply one. Scripts read with
//! `json.load(open("valenx_params.json"))` and resolve the
//! filename relative to the cwd (the workdir).
//!
//! On `collect()` we walk the workdir for `<basename>*.html`
//! (self-contained viewer pages), `<basename>*.png` (rendered
//! images), and `<basename>*.json` (NGL state JSON), plus any `*.log`
//! files Python emits.

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

use crate::case_input::NglInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(NglAdapter::new())
}

pub struct NglAdapter;

impl NglAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for NglAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "ngl";
/// Python interpreter candidates. `python3` first because on Linux
/// `python` may still be Python 2 on legacy distros; on Windows
/// `python` typically resolves to the Windows Store / 3.x install.
const PYTHON_BINARIES: &[&str] = &["python3", "python"];

impl Adapter for NglAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "NGL Viewer",
            // NGL 3.x is the modern actively-maintained line; 5.0
            // reserves room for the next two majors.
            version_range: VersionRange {
                min_inclusive: Version::new(3, 0, 0),
                max_exclusive: Version::new(5, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "MIT",
            docs_url: "https://nglviewer.org/ngl/api/",
            homepage_url: "https://nglviewer.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(PYTHON_BINARIES) {
            Some(binary_path) => {
                // Try `import nglview` — confirms the `nglview` package
                // is importable from the chosen interpreter (vs. just
                // having Python on PATH). Fall back to a "couldn't
                // import" warning so the probe still surfaces a useful
                // state — a missing-package case can still validate.
                let import_ok = nglview_importable(&binary_path);
                let mut warnings = Vec::new();
                if !import_ok {
                    warnings.push(
                        "probe found `python` on PATH but could not import \
                         `nglview` — install with `pip install nglview` for \
                         runs to succeed"
                            .into(),
                    );
                }
                Ok(ProbeReport {
                    ok: true,
                    // We intentionally don't surface a found_version
                    // — probing the package version requires running
                    // Python, which is enough work that a missing
                    // package is the more useful signal.
                    found_version: None,
                    binary_path: Some(binary_path),
                    warnings,
                    required_env: Vec::new(),
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "Python 3.8+ with the nglview Python binding \
                       installed; `pip install nglview` after ensuring \
                       python3 is on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = NglInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.ngl].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Stage the user-supplied Python script into the workdir.
        // `confined_join` rejects absolute paths and `..` traversal so
        // the staged copy stays confined to the case directory.
        let source_script = confined_join(&case.path, &input.script)?;
        if !source_script.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.ngl].script `{}` not found (resolved {})",
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
                        "[bio.ngl].script path `{}` has no filename",
                        input.script.display()
                    ),
                })?;
        let dest_script = workdir.join(script_filename);
        if source_script != dest_script {
            fs::copy(&source_script, &dest_script)?;
        }

        // Optionally stage the input structure file so the script can
        // resolve it via a bare filename inside the workdir.
        let staged_input_structure: Option<String> = match input.input_structure.as_ref() {
            Some(struct_path) => {
                let source_struct = confined_join(&case.path, struct_path)?;
                if !source_struct.is_file() {
                    return Err(AdapterError::InvalidCase {
                        case_path: case.path.join("case.toml"),
                        reason: format!(
                            "[bio.ngl].input_structure `{}` not found (resolved {})",
                            struct_path.display(),
                            source_struct.display()
                        ),
                    });
                }
                let struct_name =
                    struct_path
                        .file_name()
                        .ok_or_else(|| AdapterError::InvalidCase {
                            case_path: case.path.join("case.toml"),
                            reason: format!(
                                "[bio.ngl].input_structure path `{}` has no filename",
                                struct_path.display()
                            ),
                        })?;
                let dest_struct = workdir.join(struct_name);
                if source_struct != dest_struct {
                    fs::copy(&source_struct, &dest_struct)?;
                }
                Some(struct_name.to_string_lossy().to_string())
            }
            None => None,
        };

        // Drop a flat `valenx_params.json` into the workdir so the
        // user's Python script can read the parsed `[bio.ngl]` knobs
        // without having to reparse case.toml itself. Built by hand
        // to avoid pulling in serde_json for a 2-key flat object.
        // When `input_structure` is absent we omit the key entirely
        // (matching the spec — no `null` literal in the JSON).
        let mut params = String::new();
        params.push_str("{\n");
        params.push_str("  \"output_basename\": ");
        params.push_str(&json_string(&input.output_basename));
        if let Some(name) = staged_input_structure.as_deref() {
            params.push_str(",\n  \"input_structure\": ");
            params.push_str(&json_string(name));
        }
        params.push_str("\n}\n");
        valenx_core::io_caps::atomic_write_str(&workdir.join("valenx_params.json"), &params)?;

        // Resolve the Python binary. Bare `python` / `python3` walks
        // PATH; absolute / relative paths the user pinned are honored
        // verbatim if they exist, with a final PATH fallback.
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
            // NGL viewer / render scripts are typically dominated by
            // structure parsing + a single render pass; PNG export
            // walks the GPU pipeline once. 30 minutes is a generous
            // default; multi-frame trajectory exports might approach
            // it.
            estimated_runtime: Some(Duration::from_secs(30 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting ngl", |line| {
            let mut hint = subprocess::Hint::default();
            // Convention: the user-supplied script can emit a sentinel
            // line `[valenx] ngl done` to signal completion before
            // exit; lift to a 95% progress tick.
            if line.contains("[valenx] ngl done") {
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
        // Provenance: hash the staged Python script as the canonical
        // input descriptor. Falls back to case.toml when no script
        // is staged yet (partial / failed runs).
        let script_path = first_script_in_workdir(&job.workdir);
        let case_hash_input = script_path
            .clone()
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "ngl",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Restrict typed outputs (.html / .png / .json) to those whose
        // stem starts with the configured `output_basename`. .log
        // files are accepted regardless — Python's logging output
        // isn't typically prefixed.
        let basename = read_output_basename(&job.workdir);

        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-ngl", ?e, "workdir read failed");
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
                Some("html") => {
                    if !stem_matches_basename {
                        continue;
                    }
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Native,
                        checksum: None,
                        label: "NGL viewer HTML".to_string(),
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
                        label: "NGL rendered image".to_string(),
                    });
                }
                Some("json") => {
                    if !stem_matches_basename {
                        continue;
                    }
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Tabular,
                        checksum: None,
                        label: "NGL state JSON".to_string(),
                    });
                }
                Some("log") => {
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Log,
                        checksum: None,
                        label: "NGL log".to_string(),
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
            ribbon_contributions: vec!["bio.ngl.view"],
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

/// Lift the staged Python script out of the workdir for provenance
/// hashing. Returns the lexicographically-first `.py` file (case-
/// insensitive) at the top level, or `None` if none exists yet.
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

/// Run `python -c "import nglview"` and return whether it succeeded.
/// Returns `false` on any failure (interpreter unusable, nglview not
/// importable); `probe()` lifts that to a "nglview not importable"
/// warning rather than a hard error so a missing-package case can
/// still validate.
fn nglview_importable(python_binary: &Path) -> bool {
    let output = std::process::Command::new(python_binary)
        .arg("-c")
        .arg("import nglview")
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
        let info = NglAdapter::new().info();
        assert_eq!(info.id, "ngl");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "MIT");
        assert_eq!(info.display_name, "NGL Viewer");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = NglAdapter::new().info();
        // NGL 3.x is the modern actively-maintained line; 5.0
        // reserves room for the next two majors.
        assert_eq!(info.version_range.min_inclusive, Version::new(3, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(5, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = NglAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.ngl.view"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = NglAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
