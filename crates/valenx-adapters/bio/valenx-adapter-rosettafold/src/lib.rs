//! # valenx-adapter-rosettafold
//!
//! Adapter for [RoseTTAFold](https://github.com/RosettaCommons/RoseTTAFold) —
//! the Baker lab's original three-track structure prediction network
//! that predicts protein structure from sequence using coupled 1D, 2D,
//! and 3D representations. MIT-licensed; sister to ESMFold / OpenFold /
//! AlphaFold 2/3 from Phase 17.5.
//!
//! **Phase 17.7 — subprocess wrapper for user-provided scripts.** The
//! user supplies a `predict.py` (or whatever filename) referenced from
//! `[bio.rosettafold].script` in `case.toml`. `prepare()` stages the
//! script into the workdir alongside the FASTA query and `run()`
//! invokes `python <script>` via the shared subprocess runner. The
//! script is responsible for reading `valenx_params.json` and writing
//! its outputs under the configured `output_basename` prefix.
//!
//! ## `valenx_params.json`
//!
//! RoseTTAFold has no canonical CLI — every install drives the model
//! through its own per-site predict script. Rather than guess at a flag
//! layout, `prepare()` writes a flat JSON file `valenx_params.json`
//! into the workdir alongside the staged script and FASTA, containing
//! the parsed `[bio.rosettafold]` knobs:
//!
//! ```json
//! {
//!   "output_basename": "predicted",
//!   "fasta":           "query.fasta"
//! }
//! ```
//!
//! User scripts read it with `json.load(open("valenx_params.json"))`
//! and pass the values through to RoseTTAFold themselves. This keeps
//! the adapter free of upstream API churn and means `case.toml` knobs
//! actually reach the model.
//!
//! On `collect()` we walk the workdir for files matching the configured
//! `<output_basename>` prefix:
//!
//! - `<output_basename>*.pdb` → predicted structures (`Native`).
//! - `<output_basename>*.npz` → confidence / distogram arrays (`Native`).
//! - `*.log`                  → run logs (`Log`).
//!
//! No PDB parsing — the model writes a variety of intermediate and
//! final structures (reference, refined, ensemble members) and the user
//! is the source of truth for which to inspect downstream.

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

use crate::case_input::RoseTTAFoldInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(RoseTTAFoldAdapter::new())
}

pub struct RoseTTAFoldAdapter;

impl RoseTTAFoldAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for RoseTTAFoldAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "rosettafold";
/// Python interpreter candidates. `python3` first because on Linux
/// `python` may still be Python 2 on legacy distros; on Windows
/// `python` typically resolves to the Windows Store / 3.x install.
const PYTHON_BINARIES: &[&str] = &["python3", "python"];

impl Adapter for RoseTTAFoldAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "RoseTTAFold",
            // Per spec: 1.0.0..3.0.0. The Baker lab repo isn't versioned
            // as a pip package — this band tracks the conceptual major
            // generations of the network architecture.
            version_range: VersionRange {
                min_inclusive: Version::new(1, 0, 0),
                max_exclusive: Version::new(3, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "MIT",
            docs_url: "https://github.com/RosettaCommons/RoseTTAFold",
            homepage_url: "https://www.bakerlab.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(PYTHON_BINARIES) {
            Some(binary_path) => {
                // RoseTTAFold isn't a pip package — it's a research
                // codebase the user clones, installs into a conda env,
                // and drives via per-site predict scripts. We can't
                // `import rosettafold` to confirm install state, so the
                // probe just verifies Python is reachable and surfaces a
                // bundling caveat as a warning.
                let warnings = vec![
                    "RoseTTAFold model weights + dependencies are not bundled — \
                     clone https://github.com/RosettaCommons/RoseTTAFold and \
                     follow the install README"
                        .into(),
                ];
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
                hint: "Python 3.8+ with RoseTTAFold installed; clone \
                       https://github.com/RosettaCommons/RoseTTAFold and \
                       follow the install README, then ensure python3 is on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = RoseTTAFoldInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.rosettafold].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Enforce `.py` extension on the script. RoseTTAFold has no
        // bundled CLI — the user is expected to supply a Python predict
        // script, and a non-`.py` path here is almost always a typo.
        let script_ext_ok = input
            .script
            .extension()
            .and_then(|s| s.to_str())
            .map(|s| s.eq_ignore_ascii_case("py"))
            .unwrap_or(false);
        if !script_ext_ok {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.rosettafold].script `{}` must have a `.py` extension",
                    input.script.display()
                ),
            });
        }

        // Stage the user-supplied Python script. Resolved against the
        // case directory; absolute paths and `..` traversal are
        // rejected by `confined_join` so the staged copy stays
        // confined to the case directory.
        let source_script = confined_join(&case.path, &input.script)?;
        if !source_script.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.rosettafold].script `{}` not found (resolved {})",
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
                        "[bio.rosettafold].script path `{}` has no filename",
                        input.script.display()
                    ),
                })?;
        let dest_script = workdir.join(script_filename);
        if source_script != dest_script {
            fs::copy(&source_script, &dest_script)?;
        }

        // Stage the FASTA query alongside the script.
        let source_fasta = confined_join(&case.path, &input.fasta)?;
        if !source_fasta.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.rosettafold].fasta `{}` not found (resolved {})",
                    input.fasta.display(),
                    source_fasta.display()
                ),
            });
        }
        let fasta_filename = input
            .fasta
            .file_name()
            .ok_or_else(|| AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.rosettafold].fasta path `{}` has no filename",
                    input.fasta.display()
                ),
            })?;
        let dest_fasta = workdir.join(fasta_filename);
        if source_fasta != dest_fasta {
            fs::copy(&source_fasta, &dest_fasta)?;
        }

        // Drop a flat `valenx_params.json` into the workdir so the
        // user's predict script can read the parsed `[bio.rosettafold]`
        // knobs without having to reparse case.toml itself. Built by
        // hand to avoid pulling in a serde_json dep for a 2-key flat
        // object.
        let params_json = format!(
            "{{\n  \"output_basename\": {},\n  \"fasta\": {}\n}}\n",
            json_string(&input.output_basename),
            json_string(&fasta_filename.to_string_lossy()),
        );
        valenx_core::io_caps::atomic_write_str(&workdir.join("valenx_params.json"), &params_json)?;

        // Resolve the Python binary. Same logic as every other
        // Phase 17 Python-script adapter (Biopython / OpenMM / RDKit /
        // MDAnalysis): bare `python` / `python3` walks PATH; absolute
        // paths or pinned interpreters are honored verbatim.
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

        // No positional arguments — RoseTTAFold predict scripts read
        // their inputs via `valenx_params.json` (the contract this
        // adapter establishes).
        let native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from(script_filename),
        ];

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // RoseTTAFold inference is dominated by MSA generation and
            // model forward passes; on a single GPU a small target
            // finishes in minutes, large complexes can run for an hour
            // plus. 4 hours is a generous default; long runs override
            // through their own progress reporting.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting RoseTTAFold", |line| {
            let mut hint = subprocess::Hint::default();
            // Convention: the user-supplied script can emit a sentinel
            // line `[valenx] rosettafold done` to signal completion
            // before exit; lift to a 95% progress tick.
            if line.contains("[valenx] rosettafold done") {
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
        // Provenance: hash the staged FASTA (the canonical
        // "this case is configured this way" input). Falls back to the
        // staged script, then case.toml, when the FASTA isn't present
        // yet.
        let fasta_path = first_fasta_in_workdir(&job.workdir);
        let script_path = first_script_in_workdir(&job.workdir);
        let case_hash_input = fasta_path
            .clone()
            .or_else(|| script_path.clone())
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "RoseTTAFold",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Lift the configured output_basename out of valenx_params.json
        // so collect() classifies model outputs vs. unrelated PDBs the
        // user might have dropped into the workdir. Falls back to the
        // empty prefix (matches every PDB / NPZ) if the params file is
        // missing.
        let output_basename = read_output_basename(&job.workdir).unwrap_or_default();

        // Surface the staged FASTA so the user can find their query
        // sequence next to the predictions.
        if let Some(p) = fasta_path {
            artefacts.push(Artifact {
                path: p,
                kind: ArtifactKind::Other,
                checksum: None,
                label: "RoseTTAFold input FASTA".to_string(),
            });
        }
        // Surface the staged script as well — it's the canonical
        // record of which model variant was actually run.
        if let Some(p) = script_path {
            artefacts.push(Artifact {
                path: p,
                kind: ArtifactKind::Other,
                checksum: None,
                label: "RoseTTAFold script".to_string(),
            });
        }

        // Walk the workdir and classify outputs by extension. PDBs and
        // NPZs get the `<output_basename>*` filter; logs are always
        // collected regardless of name.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-rosettafold", ?e, "workdir read failed");
                return Ok(results);
            }
        };
        let mut pdb_paths: Vec<PathBuf> = Vec::new();
        let mut npz_paths: Vec<PathBuf> = Vec::new();
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
                .file_name()
                .and_then(|s| s.to_str())
                .unwrap_or_default();
            match ext.as_deref() {
                Some("pdb") if stem.starts_with(&output_basename) => pdb_paths.push(path),
                Some("npz") if stem.starts_with(&output_basename) => npz_paths.push(path),
                Some("log") => artefacts.push(Artifact {
                    path,
                    kind: ArtifactKind::Log,
                    checksum: None,
                    label: "RoseTTAFold log".to_string(),
                }),
                _ => continue,
            }
        }
        pdb_paths.sort();
        npz_paths.sort();
        for path in pdb_paths {
            artefacts.push(Artifact {
                path,
                kind: ArtifactKind::Native,
                checksum: None,
                label: "RoseTTAFold predicted structure".to_string(),
            });
        }
        for path in npz_paths {
            artefacts.push(Artifact {
                path,
                kind: ArtifactKind::Native,
                checksum: None,
                label: "RoseTTAFold confidence arrays".to_string(),
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
            ribbon_contributions: vec!["bio.rosettafold.predict"],
        }
    }
}

/// Escape a string for embedding inside a JSON string literal.
/// Handles the JSON string-escape set (backslash, double-quote, the
/// standard `\n` / `\r` / `\t` / `\b` / `\f` and any other ASCII
/// control as `\u00XX`). Avoids pulling in a serde_json dependency
/// for the trivially small `valenx_params.json` we emit.
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
/// Returns the lexicographically-first `.fasta` / `.fa` / `.faa` /
/// `.fna` file at the top level, or `None` if none exists yet.
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

/// Read the configured `output_basename` back out of the
/// `valenx_params.json` we wrote during prepare(). Returns `None` if
/// the file is missing or unparseable; collect() falls back to an
/// empty prefix in that case.
fn read_output_basename(workdir: &Path) -> Option<String> {
    let text = valenx_core::io_caps::read_capped_to_string(
        &workdir.join("valenx_params.json"),
        valenx_core::io_caps::MAX_ADAPTER_PARAMS_BYTES as usize,
    )
    .ok()?;
    // Trivial single-key extractor — matches the shape we emit.
    // Handles the embedded escapes we produce (backslash + double
    // quote) but isn't a general JSON parser; that's intentional, the
    // file is a contract we control end-to-end.
    let key = "\"output_basename\"";
    let pos = text.find(key)?;
    let after_key = &text[pos + key.len()..];
    let colon = after_key.find(':')?;
    let after_colon = &after_key[colon + 1..];
    let quote = after_colon.find('"')?;
    let after_open = &after_colon[quote + 1..];
    let mut out = String::new();
    let mut chars = after_open.chars();
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
        let info = RoseTTAFoldAdapter::new().info();
        assert_eq!(info.id, "rosettafold");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "MIT");
        assert_eq!(info.display_name, "RoseTTAFold");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = RoseTTAFoldAdapter::new().info();
        assert_eq!(info.version_range.min_inclusive, Version::new(1, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(3, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = RoseTTAFoldAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.rosettafold.predict"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = RoseTTAFoldAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
