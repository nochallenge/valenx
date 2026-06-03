//! # valenx-adapter-esm-if
//!
//! Adapter for [Meta ESM-IF](https://github.com/facebookresearch/esm/tree/main/examples/inverse_folding)
//! — the inverse-folding head of the ESM family. Where ESMFold takes a
//! sequence and predicts the structure, ESM-IF goes the other way:
//! given a backbone (e.g. one RFdiffusion or Chroma just sampled), it
//! samples amino-acid sequences that are likely to fold to it.
//! ESM-IF ships as part of the `esm` Python package — same install as
//! ESMFold — so the probe checks `import esm` rather than a separate
//! package.
//!
//! **Phase 27.5 — subprocess wrapper for user-provided scripts.** The
//! user supplies a `design_esmif.py` (or whatever filename) referenced
//! from `[bio.esm-if].script` in `case.toml` plus an input PDB.
//! `prepare()` stages the script + PDB into the workdir and `run()`
//! invokes `python <script>` via the shared subprocess runner. The
//! script is responsible for invoking ESM-IF with the parsed knobs and
//! writing a FASTA file under `<output_basename>.fa`.
//!
//! ## `valenx_params.json`
//!
//! ESM-IF's API surface evolves with each `esm` release — pinning a
//! specific call shape would break on every upstream version bump.
//! Instead, `prepare()` writes a flat JSON file `valenx_params.json`
//! into the workdir containing the parsed `[bio.esm-if]` knobs:
//!
//! ```json
//! {
//!   "input_pdb":       "backbone.pdb",
//!   "model":           "esm_if1_gvp4_t16_142M_UR50",
//!   "temperature":     1.0,
//!   "num_samples":     8,
//!   "output_basename": "design"
//! }
//! ```
//!
//! User scripts read it with `json.load(open("valenx_params.json"))`
//! and pass the values through to ESM-IF themselves.
//!
//! On `collect()` we walk the workdir for `<output_basename>.fa` and
//! parse it via [`valenx_bio::format::fasta::read`]. Successful parses
//! get the richer `"ESM-IF · N sequences"` label; otherwise the
//! generic `"ESM-IF designed sequences"` falls through.

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

use crate::case_input::EsmIfInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(EsmIfAdapter::new())
}

pub struct EsmIfAdapter;

impl EsmIfAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for EsmIfAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "esm-if";
/// Python interpreter candidates. `python3` first because on Linux
/// `python` may still be Python 2 on legacy distros; on Windows
/// `python` typically resolves to the Windows Store / 3.x install.
const PYTHON_BINARIES: &[&str] = &["python3", "python"];

impl Adapter for EsmIfAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "ESM-IF",
            // ESM-IF lives inside the `esm` package; the inverse-
            // folding head landed in 0.4 and the whole family has
            // stayed pre-1.0 for a while. Lower bound 0.4 captures
            // the first ESM-IF-bearing release; upper bound 2.0
            // reserves room for an upcoming major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(0, 4, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "MIT",
            docs_url: "https://github.com/facebookresearch/esm/tree/main/examples/inverse_folding",
            homepage_url: "https://github.com/facebookresearch/esm",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(PYTHON_BINARIES) {
            Some(binary_path) => {
                // ESM-IF ships in the same `esm` package as ESMFold —
                // confirms the package is importable from the chosen
                // interpreter (vs. just having Python on PATH).
                let found_version = detect_esm_version(&binary_path);
                let mut warnings = Vec::new();
                if found_version.is_none() {
                    warnings.push(
                        "probe found `python` on PATH but could not import \
                         `esm` — install Meta's ESM package (which ships the \
                         inverse-folding head as a sub-module) from \
                         https://github.com/facebookresearch/esm and ensure \
                         it's importable from the chosen interpreter for \
                         runs to succeed"
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
                hint: "Python 3.9+ with the `esm` package installed; install \
                       from https://github.com/facebookresearch/esm after \
                       ensuring python3 is on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = EsmIfInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.esm-if].output_basename",
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
                    "[bio.esm-if].script `{}` not found (resolved {})",
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
                        "[bio.esm-if].script path `{}` has no filename",
                        input.script.display()
                    ),
                })?;
        let dest_script = workdir.join(script_filename);
        if source_script != dest_script {
            fs::copy(&source_script, &dest_script)?;
        }

        // Stage the input PDB alongside the script.
        let source_pdb = confined_join(&case.path, &input.input_pdb)?;
        if !source_pdb.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.esm-if].input_pdb `{}` not found (resolved {})",
                    input.input_pdb.display(),
                    source_pdb.display()
                ),
            });
        }
        let pdb_filename =
            input
                .input_pdb
                .file_name()
                .ok_or_else(|| AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.esm-if].input_pdb path `{}` has no filename",
                        input.input_pdb.display()
                    ),
                })?;
        let dest_pdb = workdir.join(pdb_filename);
        if source_pdb != dest_pdb {
            fs::copy(&source_pdb, &dest_pdb)?;
        }

        // Drop a flat `valenx_params.json` into the workdir so the
        // user's design script can read the parsed `[bio.esm-if]`
        // knobs without having to reparse case.toml itself. Built by
        // hand to avoid pulling in a serde_json dep for a 5-key flat
        // object.
        let params_json = format!(
            "{{\n  \"input_pdb\": {},\n  \"model\": {},\n  \"temperature\": {},\n  \"num_samples\": {},\n  \"output_basename\": {}\n}}\n",
            json_string(&pdb_filename.to_string_lossy()),
            json_string(&input.model),
            format_finite_float(input.temperature),
            input.num_samples,
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
            // ESM-IF inference is much faster than diffusion-based
            // design — typical runs finish in seconds to a few
            // minutes. 30 minutes is a generous default, mirroring
            // ProteinMPNN.
            estimated_runtime: Some(Duration::from_secs(30 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting ESM-IF", |line| {
            let mut hint = subprocess::Hint::default();
            // Convention: the user-supplied script can emit a sentinel
            // line `[valenx] esm-if done` to signal completion before
            // exit; lift to a 95% progress tick.
            if line.contains("[valenx] esm-if done") {
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
        // Provenance: hash the staged input PDB. Falls back to the
        // script, then case.toml, when the PDB isn't present yet.
        let pdb_path = first_pdb_in_workdir(&job.workdir);
        let script_path = first_script_in_workdir(&job.workdir);
        let case_hash_input = pdb_path
            .clone()
            .or_else(|| script_path.clone())
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "ESM-IF",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Read the staged params back so we can locate the expected
        // `<output_basename>.fa`. Failure is non-fatal — collect
        // still surfaces every FASTA it can find.
        let basename = read_params(&job.workdir);

        if let Some(p) = pdb_path {
            artefacts.push(Artifact {
                path: p,
                kind: ArtifactKind::Other,
                checksum: None,
                label: "ESM-IF input PDB".to_string(),
            });
        }
        if let Some(p) = script_path {
            artefacts.push(Artifact {
                path: p,
                kind: ArtifactKind::Other,
                checksum: None,
                label: "ESM-IF script".to_string(),
            });
        }

        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-esm-if", ?e, "workdir read failed");
                return Ok(results);
            }
        };
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
            match ext.as_deref() {
                Some("fa") | Some("fasta") | Some("faa") => {
                    // If we know the configured basename, restrict to
                    // matching outputs — otherwise accept everything
                    // (best-effort).
                    let stem_ok = match basename.as_ref() {
                        Some(b) => {
                            let stem = path.file_stem().and_then(|s| s.to_str()).unwrap_or("");
                            stem == b.as_str()
                        }
                        None => true,
                    };
                    if stem_ok {
                        fasta_paths.push(path);
                    }
                }
                Some("log") => artefacts.push(Artifact {
                    path,
                    kind: ArtifactKind::Log,
                    checksum: None,
                    label: "ESM-IF log".to_string(),
                }),
                _ => continue,
            }
        }
        fasta_paths.sort();
        for path in fasta_paths {
            // Soft-validate FASTA but never fail — a partial fasta
            // mid-run shouldn't make collect() blow up. Successful
            // parses get the richer "N sequences" label.
            //
            // Round-23 named finding: bound the FASTA read at
            // MAX_BIO_CLI_BYTES (256 MiB) so a poisoned or runaway
            // designed-sequences file can't OOM the renderer before
            // the parser counts sequences. Same magnitude as the
            // bio-CLI inspector cap, generous for chromosome-scale
            // FASTA but refuses the cat /dev/zero denial of service.
            let label = match valenx_core::io_caps::read_capped_to_string(
                &path,
                valenx_core::io_caps::MAX_BIO_CLI_BYTES as usize,
            ) {
                Ok(text) => {
                    match valenx_bio::format::fasta::read(&text, valenx_bio::Alphabet::Protein) {
                        Ok(seqs) => format!(
                            "ESM-IF \u{00b7} {} sequence{}",
                            seqs.len(),
                            if seqs.len() == 1 { "" } else { "s" }
                        ),
                        Err(_) => "ESM-IF designed sequences".to_string(),
                    }
                }
                Err(_) => "ESM-IF designed sequences".to_string(),
            };
            artefacts.push(Artifact {
                path,
                kind: ArtifactKind::Tabular,
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
            ribbon_contributions: vec!["bio.esm-if.design"],
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

/// Lift the staged input PDB out of the workdir for provenance
/// hashing. Prefers the PDB referenced by `valenx_params.json` when
/// present; if the params can't be read, falls back to the
/// lexicographically-first `.pdb` file at the top level.
fn first_pdb_in_workdir(workdir: &Path) -> Option<PathBuf> {
    if let Some(input_pdb_name) = read_params_input_pdb(workdir) {
        let candidate = workdir.join(&input_pdb_name);
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
                .map(|s| s.eq_ignore_ascii_case("pdb"))
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
/// `collect()` can restrict FASTA pickup to the configured stem.
fn read_params(workdir: &Path) -> Option<String> {
    let text = valenx_core::io_caps::read_capped_to_string(
        &workdir.join("valenx_params.json"),
        valenx_core::io_caps::MAX_ADAPTER_PARAMS_BYTES as usize,
    )
    .ok()?;
    extract_json_string(&text, "output_basename")
}

/// Read just the `input_pdb` field from `valenx_params.json`.
fn read_params_input_pdb(workdir: &Path) -> Option<String> {
    let text = valenx_core::io_caps::read_capped_to_string(
        &workdir.join("valenx_params.json"),
        valenx_core::io_caps::MAX_ADAPTER_PARAMS_BYTES as usize,
    )
    .ok()?;
    extract_json_string(&text, "input_pdb")
}

/// Pull a flat string field out of our own hand-emitted
/// `valenx_params.json`. Trivially small — we wrote the file
/// ourselves so we know its shape.
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
/// malformed); `probe()` falls back to a "esm not importable" warning
/// in that case.
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

    /// A minimal multi-record FASTA so the collect-test can exercise
    /// the "N sequences" label path.
    const SAMPLE_FASTA: &str = ">design_0\nACDEF\n>design_1\nGHIKL\n";

    #[test]
    fn info_is_bio_domain() {
        let info = EsmIfAdapter::new().info();
        assert_eq!(info.id, "esm-if");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "MIT");
        assert_eq!(info.display_name, "ESM-IF");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = EsmIfAdapter::new().info();
        assert_eq!(info.version_range.min_inclusive, Version::new(0, 4, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(2, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = EsmIfAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.esm-if.design"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = EsmIfAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    /// `collect()` should walk the workdir for the configured
    /// `<output_basename>.fa`, parse it via the canonical FASTA
    /// reader, and surface it as a `Tabular` artifact with the
    /// "ESM-IF \u{00b7} N sequences" label.
    #[test]
    fn collect_walks_workdir_and_classifies_outputs() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-esm-if-collect-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        // Stage the params so collect() can find the configured
        // basename.
        fs::write(
            tmp.join("valenx_params.json"),
            "{\n  \"input_pdb\": \"backbone.pdb\",\n  \"model\": \"esm_if1_gvp4_t16_142M_UR50\",\n  \"temperature\": 1.0,\n  \"num_samples\": 2,\n  \"output_basename\": \"design\"\n}\n",
        )
        .unwrap();
        fs::write(tmp.join("backbone.pdb"), b"ATOM      1  N   ALA A   1\n").unwrap();
        fs::write(tmp.join("design.py"), b"# placeholder").unwrap();
        fs::write(tmp.join("design.fa"), SAMPLE_FASTA).unwrap();
        // A stray fasta that shouldn't be picked up as the design.
        fs::write(tmp.join("unrelated.fa"), SAMPLE_FASTA).unwrap();
        fs::write(tmp.join("run.log"), b"ESM-IF run log\n").unwrap();

        let job = PreparedJob {
            workdir: tmp.clone(),
            native_command: vec![],
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: true,
        };
        let results = EsmIfAdapter::new().collect(&job).unwrap();

        // Exactly one design FASTA picked up.
        let designs: Vec<_> = results
            .artifacts
            .iter()
            .filter(|a| a.kind == ArtifactKind::Tabular)
            .collect();
        assert_eq!(
            designs.len(),
            1,
            "expected 1 design FASTA, got {}: {:?}",
            designs.len(),
            results.artifacts
        );
        let design = designs[0];
        assert!(
            design.label.contains("ESM-IF"),
            "label was: {}",
            design.label
        );
        assert!(
            design.label.contains("2 sequences"),
            "label was: {}",
            design.label
        );

        // Input PDB tagged as Other with the documented label.
        let input_art = results
            .artifacts
            .iter()
            .find(|a| {
                a.path
                    .file_name()
                    .and_then(|s| s.to_str())
                    .map(|s| s == "backbone.pdb")
                    .unwrap_or(false)
            })
            .expect("input PDB artifact present");
        assert_eq!(input_art.kind, ArtifactKind::Other);
        assert_eq!(input_art.label, "ESM-IF input PDB");

        let py_art = results
            .artifacts
            .iter()
            .find(|a| a.path.extension().is_some_and(|e| e == "py"))
            .expect("script artifact present");
        assert_eq!(py_art.kind, ArtifactKind::Other);
        assert_eq!(py_art.label, "ESM-IF script");

        let log_art = results
            .artifacts
            .iter()
            .find(|a| a.path.extension().is_some_and(|e| e == "log"))
            .expect("log artifact present");
        assert_eq!(log_art.kind, ArtifactKind::Log);

        let _ = fs::remove_dir_all(&tmp);
    }
}
