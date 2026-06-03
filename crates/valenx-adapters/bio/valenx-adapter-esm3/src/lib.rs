//! # valenx-adapter-esm3
//!
//! Adapter for [ESM3](https://github.com/evolutionaryscale/esm) —
//! EvolutionaryScale's flagship generative multi-modal protein model.
//! Where ESMFold and ESM-IF each tackle a single direction (sequence
//! to structure, structure to sequence), ESM3 reasons jointly over
//! sequence, structure, and function tracks and can be conditioned
//! on any subset to fill in the rest. Common modes covered here:
//!
//! - `design`        — unconditional / partially-conditioned generation
//! - `inverse-fold`  — sample sequences for a given backbone PDB
//! - `scaffold`      — fill in masked regions of a structure
//! - `predict`       — sequence-conditioned structure prediction
//!
//! **Phase 27.6 — subprocess wrapper for user-provided scripts.** The
//! user supplies a `run_esm3.py` (or whatever filename) referenced
//! from `[bio.esm3].script` in `case.toml` plus optional input PDB /
//! FASTA. `prepare()` stages everything into the workdir and `run()`
//! invokes `python <script>` via the shared subprocess runner. The
//! script is responsible for invoking ESM3 with the parsed knobs and
//! writing per-sample outputs under `<output_basename>*`.
//!
//! ## `valenx_params.json`
//!
//! ESM3's API surface evolves rapidly with each `esm` release —
//! pinning a specific call shape would break on every upstream version
//! bump. Instead, `prepare()` writes a flat JSON file
//! `valenx_params.json` into the workdir containing the parsed
//! `[bio.esm3]` knobs:
//!
//! ```json
//! {
//!   "model_variant":   "open",
//!   "mode":            "inverse-fold",
//!   "num_samples":     4,
//!   "input_pdb":       "scaffold.pdb",
//!   "input_fasta":     null,
//!   "temperature":     1.0,
//!   "output_basename": "esm3_run"
//! }
//! ```
//!
//! User scripts read it with `json.load(open("valenx_params.json"))`
//! and pass the values through to ESM3 themselves.
//!
//! On `collect()` we walk the workdir for `<output_basename>*.pdb`
//! (generated structures) and `<output_basename>*.fa` (generated
//! sequences) and surface each with the appropriate artifact kind +
//! label.

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

use crate::case_input::Esm3Input;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(Esm3Adapter::new())
}

pub struct Esm3Adapter;

impl Esm3Adapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for Esm3Adapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "esm3";
/// Python interpreter candidates. `python3` first because on Linux
/// `python` may still be Python 2 on legacy distros; on Windows
/// `python` typically resolves to the Windows Store / 3.x install.
const PYTHON_BINARIES: &[&str] = &["python3", "python"];

impl Adapter for Esm3Adapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "ESM3",
            // ESM3 ships in EvolutionaryScale's `esm` package; the 3.x
            // line is the first to expose the unified generative API.
            // Upper bound 4.0 reserves room for an upcoming major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(3, 0, 0),
                max_exclusive: Version::new(4, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "Cambrian-Open-License",
            docs_url: "https://github.com/evolutionaryscale/esm",
            homepage_url: "https://www.evolutionaryscale.ai/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(PYTHON_BINARIES) {
            Some(binary_path) => {
                // Confirm the `esm` package is importable from the
                // chosen interpreter (vs. just having Python on PATH).
                // ESM3 lives in the same EvolutionaryScale `esm`
                // package as ESM-IF / ESM Cambrian.
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
        let input = Esm3Input::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.esm3].output_basename",
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
                    "[bio.esm3].script `{}` not found (resolved {})",
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
                        "[bio.esm3].script path `{}` has no filename",
                        input.script.display()
                    ),
                })?;
        let dest_script = workdir.join(script_filename);
        if source_script != dest_script {
            fs::copy(&source_script, &dest_script)?;
        }

        // Optionally stage input_pdb (relative paths only — absolute
        // paths reference verbatim, so the `.pdb` file we record in
        // params is the original absolute path).
        let input_pdb_param: Option<String> = stage_optional_input(
            input.input_pdb.as_deref(),
            &case.path,
            workdir,
            "[bio.esm3].input_pdb",
        )?;

        let input_fasta_param: Option<String> = stage_optional_input(
            input.input_fasta.as_deref(),
            &case.path,
            workdir,
            "[bio.esm3].input_fasta",
        )?;

        // Drop a flat `valenx_params.json` into the workdir so the
        // user's script can read the parsed `[bio.esm3]` knobs without
        // having to reparse case.toml itself. Built by hand to avoid
        // pulling in a serde_json dep.
        let params_json = format!(
            "{{\n  \"model_variant\": {},\n  \"mode\": {},\n  \"num_samples\": {},\n  \"input_pdb\": {},\n  \"input_fasta\": {},\n  \"temperature\": {},\n  \"output_basename\": {}\n}}\n",
            json_string(&input.model_variant),
            json_string(&input.mode),
            input.num_samples,
            json_string_or_null(input_pdb_param.as_deref()),
            json_string_or_null(input_fasta_param.as_deref()),
            format_finite_float(input.temperature),
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
            // ESM3 inference can be lengthy on the open multi-billion
            // parameter checkpoints, especially when sampling many
            // designs per input. 2 hours mirrors ESMFold's generous
            // default for big GPU jobs.
            estimated_runtime: Some(Duration::from_secs(2 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting ESM3", |line| {
            let mut hint = subprocess::Hint::default();
            // Convention: the user-supplied script can emit a sentinel
            // line `[valenx] esm3 done` to signal completion before
            // exit; lift to a 95% progress tick.
            if line.contains("[valenx] esm3 done") {
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
        // Provenance: hash the staged script (canonical "this case
        // was run this way" input, since ESM3 inputs are mode-
        // dependent). Falls back to case.toml when the script isn't
        // present yet.
        let script_path = first_script_in_workdir(&job.workdir);
        let case_hash_input = script_path
            .clone()
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "ESM3",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        let basename = read_params(&job.workdir);

        if let Some(p) = script_path {
            artefacts.push(Artifact {
                path: p,
                kind: ArtifactKind::Other,
                checksum: None,
                label: "ESM3 script".to_string(),
            });
        }

        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-esm3", ?e, "workdir read failed");
                return Ok(results);
            }
        };
        let mut pdb_paths: Vec<PathBuf> = Vec::new();
        let mut fasta_paths: Vec<PathBuf> = Vec::new();
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
            // Only pick files whose stem starts with the configured
            // basename — keeps stray inputs out of the outputs list.
            // When we couldn't read the params, accept everything.
            let stem_ok = match basename.as_ref() {
                Some(b) => stem.starts_with(b.as_str()),
                None => true,
            };
            match ext.as_deref() {
                Some("pdb") if stem_ok => pdb_paths.push(path),
                Some("fa") | Some("fasta") | Some("faa") if stem_ok => fasta_paths.push(path),
                Some("log") => artefacts.push(Artifact {
                    path,
                    kind: ArtifactKind::Log,
                    checksum: None,
                    label: "ESM3 log".to_string(),
                }),
                _ => continue,
            }
        }
        pdb_paths.sort();
        fasta_paths.sort();

        for path in pdb_paths {
            artefacts.push(Artifact {
                path,
                kind: ArtifactKind::Native,
                checksum: None,
                label: "ESM3 generated structure".to_string(),
            });
        }
        for path in fasta_paths {
            artefacts.push(Artifact {
                path,
                kind: ArtifactKind::Tabular,
                checksum: None,
                label: "ESM3 generated sequence".to_string(),
            });
        }

        artefacts.sort_by(|a, b| a.path.cmp(&b.path));
        results.artifacts = artefacts;
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.esm3.generate"],
        }
    }
}

/// Stage an optional input file (PDB / FASTA) into the workdir.
/// Relative paths resolve against the case dir and copy in; absolute
/// paths reference verbatim and are returned unchanged so
/// `valenx_params.json` records the original location.
fn stage_optional_input(
    input: Option<&Path>,
    case_dir: &Path,
    workdir: &Path,
    field: &str,
) -> Result<Option<String>, AdapterError> {
    let Some(rel) = input else {
        return Ok(None);
    };
    if rel.is_absolute() {
        // Absolute path — reference verbatim. Validate existence so
        // failures surface during prepare instead of mid-run.
        if !rel.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case_dir.join("case.toml"),
                reason: format!("{field} `{}` not found (absolute path)", rel.display()),
            });
        }
        return Ok(Some(rel.display().to_string()));
    }
    let source = confined_join(case_dir, rel)?;
    if !source.is_file() {
        return Err(AdapterError::InvalidCase {
            case_path: case_dir.join("case.toml"),
            reason: format!(
                "{field} `{}` not found (resolved {})",
                rel.display(),
                source.display()
            ),
        });
    }
    let filename = rel.file_name().ok_or_else(|| AdapterError::InvalidCase {
        case_path: case_dir.join("case.toml"),
        reason: format!("{field} path `{}` has no filename", rel.display()),
    })?;
    let dest = workdir.join(filename);
    if source != dest {
        fs::copy(&source, &dest)?;
    }
    Ok(Some(filename.to_string_lossy().into_owned()))
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

/// Emit either the JSON-quoted string or the literal `null` so the
/// optional `input_pdb` / `input_fasta` fields round-trip cleanly to
/// Python.
fn json_string_or_null(s: Option<&str>) -> String {
    match s {
        Some(value) => json_string(value),
        None => "null".to_string(),
    }
}

/// Format a finite float for embedding in JSON. The case-input
/// validator already guarantees finiteness, but we belt-and-brace it
/// here to keep the JSON output well-formed even if a future caller
/// reaches in with a non-finite value.
fn format_finite_float(f: f64) -> String {
    if f.is_finite() {
        format!("{f}")
    } else {
        "0".to_string()
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

/// Read just the `output_basename` field from `valenx_params.json` so
/// `collect()` can restrict pickup to the configured stem.
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
/// `semver::Version` out of stdout. Returns `None` on any failure
/// (interpreter unusable, esm not importable, version string
/// malformed); `probe()` falls back to an "esm not importable"
/// warning in that case.
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
        let info = Esm3Adapter::new().info();
        assert_eq!(info.id, "esm3");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "Cambrian-Open-License");
        assert_eq!(info.display_name, "ESM3");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = Esm3Adapter::new().info();
        assert_eq!(info.version_range.min_inclusive, Version::new(3, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(4, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = Esm3Adapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.esm3.generate"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = Esm3Adapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
