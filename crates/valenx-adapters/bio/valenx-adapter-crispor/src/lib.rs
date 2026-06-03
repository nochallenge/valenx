//! # valenx-adapter-crispor
//!
//! Adapter for [CRISPOR](https://crispor.org/) — Maximilian Haeussler's
//! CRISPR guide-RNA design + off-target prediction tool. CRISPOR's
//! distinguishing feature is the rigorous off-target pass: it scores
//! candidate guides against a reference genome assembly with the CFD
//! scoring model and reports an MIT-style specificity score per guide.
//! It powers the public crispor.org service and is also distributed as
//! a standalone Python script for batch / pipeline use.
//!
//! **Phase 35 — subprocess wrapper for the Python-distributed tool.**
//! The user supplies `crispor.py` (the upstream entry point or a thin
//! wrapper) referenced from `[bio.crispor].script` in `case.toml`,
//! plus a target FASTA and the design knobs (`genome`, `pam`,
//! optional `batch_id`, `output_basename`). `prepare()` stages the
//! script + FASTA into the workdir and writes a flat
//! `valenx_params.json` so the script can read the parsed knobs
//! without re-parsing `case.toml`; `run()` invokes `python <script>`
//! via the shared subprocess runner.
//!
//! ## `valenx_params.json`
//!
//! ```json
//! {
//!   "target":          "target.fa",
//!   "genome":          "hg38",
//!   "pam":             "NGG",
//!   "batch_id":        null,
//!   "output_basename": "guides"
//! }
//! ```
//!
//! `batch_id` is emitted as either a JSON string or the literal `null`
//! depending on whether the user supplied one. Scripts can therefore
//! always do `params["batch_id"]` without an `in` check.
//!
//! On `collect()` we walk the workdir for `<output_basename>*.tsv`
//! guide rankings and `<output_basename>*.txt` log files — the two
//! conventional formats CRISPOR's batch script writes.

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

use crate::case_input::CrisporInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(CrisporAdapter::new())
}

pub struct CrisporAdapter;

impl CrisporAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for CrisporAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "crispor";
/// Python interpreter candidates. `python3` first because on Linux
/// `python` may still be Python 2 on legacy distros; on Windows
/// `python` typically resolves to the Windows Store / 3.x install.
const PYTHON_BINARIES: &[&str] = &["python3", "python"];

impl Adapter for CrisporAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "CRISPOR",
            // CRISPOR's tagged release line is 5.x; the modern Python
            // 3 / batch-mode rewrite landed in 5.0. Upper bound 6.0
            // reserves room for an eventual major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(5, 0, 0),
                max_exclusive: Version::new(6, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-3.0",
            docs_url: "https://crispor.org/",
            homepage_url: "https://crispor.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(PYTHON_BINARIES) {
            Some(binary_path) => {
                // Try `import crispor` to confirm the Python
                // distribution is installed in the chosen
                // interpreter (vs. just having Python on PATH).
                // When the import fails we still return `ok: true`
                // with a warning so the probe surfaces a useful
                // state — the user might have CRISPOR cloned to a
                // PYTHONPATH-controlled directory and intend to
                // invoke it via the user-supplied script rather
                // than `import crispor` directly.
                let importable = crispor_importable(&binary_path);
                let mut warnings = Vec::new();
                if !importable {
                    warnings.push(
                        "probe found `python` on PATH but could not import \
                         `crispor` — install CRISPOR from \
                         https://github.com/maximilianh/crisporWebsite (clone \
                         the source repo and add it to PYTHONPATH) for runs \
                         to succeed"
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
                hint: "Python 3.9+ with CRISPOR installed; clone \
                       https://github.com/maximilianh/crisporWebsite and \
                       ensure python3 is on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = CrisporInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.crispor].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Stage the user-supplied Python script. Resolved against the
        // case directory — same convention as every other Phase 17/27
        // bio Python adapter. `confined_join` rejects absolute paths
        // and `..` traversal so a malicious case bundle can't smuggle
        // arbitrary host files into the workdir.
        let source_script = confined_join(&case.path, &input.script)?;
        if !source_script.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.crispor].script `{}` not found (resolved {})",
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
                        "[bio.crispor].script path `{}` has no filename",
                        input.script.display()
                    ),
                })?;
        let dest_script = workdir.join(script_filename);
        if source_script != dest_script {
            fs::copy(&source_script, &dest_script)?;
        }

        // Stage the target FASTA alongside the script.
        let source_target = confined_join(&case.path, &input.target)?;
        if !source_target.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.crispor].target `{}` not found (resolved {})",
                    input.target.display(),
                    source_target.display()
                ),
            });
        }
        let target_filename =
            input
                .target
                .file_name()
                .ok_or_else(|| AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.crispor].target path `{}` has no filename",
                        input.target.display()
                    ),
                })?;
        let dest_target = workdir.join(target_filename);
        if source_target != dest_target {
            fs::copy(&source_target, &dest_target)?;
        }

        // Drop a flat `valenx_params.json` into the workdir so the
        // user's script can read the parsed `[bio.crispor]` knobs
        // without having to reparse case.toml itself. `batch_id` is
        // emitted as JSON null when the user didn't supply one so
        // scripts can always do `params["batch_id"]` without
        // checking for membership first. Built by hand to avoid
        // pulling in serde_json for a 5-key flat object.
        let batch_id_json = match &input.batch_id {
            Some(s) => json_string(s),
            None => "null".to_string(),
        };
        let params_json = format!(
            "{{\n  \"target\": {},\n  \"genome\": {},\n  \"pam\": {},\n  \"batch_id\": {},\n  \"output_basename\": {}\n}}\n",
            json_string(&target_filename.to_string_lossy()),
            json_string(&input.genome),
            json_string(&input.pam),
            batch_id_json,
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
            // CRISPOR's off-target pass walks a BWA / bowtie index
            // of the reference genome — minutes for a small target,
            // up to an hour for whole-genome scans on large
            // assemblies. 2 hours covers the long tail.
            estimated_runtime: Some(Duration::from_secs(2 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting CRISPOR", |line| {
            let mut hint = subprocess::Hint::default();
            // Convention: the user-supplied script can emit a
            // sentinel line `[valenx] crispor done` to signal
            // completion before exit; lift to a 95% progress tick.
            if line.contains("[valenx] crispor done") {
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
        // Provenance: hash the staged target FASTA. Falls back to
        // the script, then case.toml, when the FASTA isn't present
        // yet (partial / failed runs).
        let target_path = first_target_in_workdir(&job.workdir);
        let script_path = first_script_in_workdir(&job.workdir);
        let case_hash_input = target_path
            .clone()
            .or_else(|| script_path.clone())
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "CRISPOR",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Read the staged params back out so we can restrict the
        // collected outputs to those matching `output_basename`.
        // Failure to read the params is non-fatal — collect still
        // surfaces every TSV / TXT.
        let basename = read_output_basename(&job.workdir);

        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-crispor", ?e, "workdir read failed");
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
            // Restrict to outputs whose stem starts with the
            // configured `output_basename`. If params couldn't be
            // read, accept everything (best-effort).
            let stem_ok = match basename.as_deref() {
                Some(b) => stem.starts_with(b),
                None => true,
            };
            if !stem_ok {
                continue;
            }
            let (kind, label) = match ext.as_deref() {
                Some("tsv") => (ArtifactKind::Tabular, "CRISPOR guide rankings".to_string()),
                Some("txt") => (ArtifactKind::Log, "CRISPOR log".to_string()),
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
            ribbon_contributions: vec!["bio.crispor.design"],
        }
    }
}

/// Probe whether `import crispor` succeeds in the given Python
/// interpreter. Returns true only when Python exits 0 — any failure
/// (interpreter missing, module not on sys.path, import-time
/// exception) is treated as "not importable" and surfaced via a
/// probe warning rather than blocking adapter use.
fn crispor_importable(python_binary: &Path) -> bool {
    let output = std::process::Command::new(python_binary)
        .arg("-c")
        .arg("import crispor")
        .stdin(std::process::Stdio::null())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .output();
    matches!(output, Ok(o) if o.status.success())
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

/// Lift the staged target FASTA out of the workdir for provenance
/// hashing. Prefers the file referenced by `valenx_params.json` when
/// present; if the params can't be read, falls back to the
/// lexicographically-first `.fa` / `.fasta` / `.fna` file at the top
/// level.
fn first_target_in_workdir(workdir: &Path) -> Option<PathBuf> {
    if let Some(name) = read_target_filename(workdir) {
        let candidate = workdir.join(&name);
        if candidate.is_file() {
            return Some(candidate);
        }
    }
    let entries = fs::read_dir(workdir).ok()?;
    let mut hits: Vec<PathBuf> = entries
        .flatten()
        .map(|e| e.path())
        .filter(|p| {
            p.extension()
                .and_then(|s| s.to_str())
                .map(|s| {
                    let l = s.to_ascii_lowercase();
                    matches!(l.as_str(), "fa" | "fasta" | "fna" | "ffn")
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

/// Pull `output_basename` out of our own hand-emitted
/// `valenx_params.json` for collect()-time filtering.
fn read_output_basename(workdir: &Path) -> Option<String> {
    let text = valenx_core::io_caps::read_capped_to_string(
        &workdir.join("valenx_params.json"),
        valenx_core::io_caps::MAX_ADAPTER_PARAMS_BYTES as usize,
    )
    .ok()?;
    extract_json_string(&text, "output_basename")
}

/// Pull `target` (the staged FASTA filename) out of our own
/// hand-emitted `valenx_params.json` for provenance.
fn read_target_filename(workdir: &Path) -> Option<String> {
    let text = valenx_core::io_caps::read_capped_to_string(
        &workdir.join("valenx_params.json"),
        valenx_core::io_caps::MAX_ADAPTER_PARAMS_BYTES as usize,
    )
    .ok()?;
    extract_json_string(&text, "target")
}

/// Pull a flat string field out of our own hand-emitted
/// `valenx_params.json`. Trivially small — we wrote the file
/// ourselves so we know its shape; a full JSON parser would be
/// overkill (and pulling in serde_json just for collect()'s
/// side-band metadata would bloat the dep tree).
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
        let info = CrisporAdapter::new().info();
        assert_eq!(info.id, "crispor");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "GPL-3.0");
        assert_eq!(info.display_name, "CRISPOR");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = CrisporAdapter::new().info();
        // CRISPOR 5.x is the modern Python-3 / batch-mode line;
        // upper bound 6.0 reserves room for an eventual major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(5, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(6, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = CrisporAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.crispor.design"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = CrisporAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
