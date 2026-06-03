//! # valenx-adapter-proteinmpnn
//!
//! Adapter for [ProteinMPNN](https://github.com/dauparas/ProteinMPNN)
//! — Justas Dauparas et al.'s graph-neural-network sequence designer.
//! The pair-counterpart of RFdiffusion: take a protein backbone (e.g.
//! one RFdiffusion just sampled) and design amino-acid sequences that
//! fold to it. ProteinMPNN ships three weight families (vanilla,
//! soluble-protein-biased, and CA-only) and a temperature knob that
//! trades sequence diversity against likelihood.
//!
//! **Phase 27 — subprocess wrapper for user-provided scripts.** The
//! user supplies a `design_proteinmpnn.py` (or whatever filename)
//! referenced from `[bio.proteinmpnn].script` in `case.toml` plus an
//! input PDB. `prepare()` stages the script + PDB into the workdir
//! and `run()` invokes `python <script>` via the shared subprocess
//! runner. The script is responsible for invoking ProteinMPNN with
//! the parsed knobs and writing a FASTA file under
//! `<output_basename>.fa`.
//!
//! ## `valenx_params.json`
//!
//! ProteinMPNN's CLI is a thin wrapper over its own argparse — fast
//! moving and not a stable contract for downstream automation.
//! Instead, `prepare()` writes a flat JSON file `valenx_params.json`
//! into the workdir containing the parsed `[bio.proteinmpnn]` knobs:
//!
//! ```json
//! {
//!   "model_variant":      "vanilla",
//!   "temperature":        0.1,
//!   "num_seq_per_target": 8,
//!   "output_basename":    "design",
//!   "input_pdb":          "backbone.pdb"
//! }
//! ```
//!
//! User scripts read it with `json.load(open("valenx_params.json"))`
//! and pass the values through to ProteinMPNN themselves.
//!
//! On `collect()` we walk the workdir for `<output_basename>.fa` and
//! parse it via [`valenx_bio::format::fasta::read`]. Each is surfaced
//! as a typed [`ArtifactKind::Tabular`] artifact with a
//! "ProteinMPNN designed sequences" label (or richer
//! `"ProteinMPNN · N sequences"` when the FASTA parses cleanly).

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

use crate::case_input::ProteinMpnnInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(ProteinMpnnAdapter::new())
}

pub struct ProteinMpnnAdapter;

impl ProteinMpnnAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for ProteinMpnnAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "proteinmpnn";
/// Python interpreter candidates. `python3` first because on Linux
/// `python` may still be Python 2 on legacy distros; on Windows
/// `python` typically resolves to the Windows Store / 3.x install.
const PYTHON_BINARIES: &[&str] = &["python3", "python"];

impl Adapter for ProteinMpnnAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "ProteinMPNN",
            // ProteinMPNN's first tagged release is the 1.0 weights /
            // inference code drop. Upper bound 2.0 reserves room for
            // an upcoming major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(1, 0, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "MIT",
            docs_url: "https://github.com/dauparas/ProteinMPNN",
            homepage_url: "https://github.com/dauparas/ProteinMPNN",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(PYTHON_BINARIES) {
            Some(binary_path) => {
                // ProteinMPNN exports its package as `ProteinMPNN`
                // (PascalCase). Some installs run as a script in
                // place rather than installing a package — fall back
                // gracefully to a "couldn't import" warning so the
                // probe still surfaces a useful state.
                let found_version = detect_proteinmpnn_version(&binary_path);
                let mut warnings = Vec::new();
                if found_version.is_none() {
                    warnings.push(
                        "probe found `python` on PATH but could not import \
                         `ProteinMPNN` — clone \
                         https://github.com/dauparas/ProteinMPNN and ensure \
                         the package is importable from the chosen \
                         interpreter for runs to succeed"
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
                hint: "Python 3.9+ with ProteinMPNN installed; clone \
                       https://github.com/dauparas/ProteinMPNN and follow \
                       the install steps after ensuring python3 is on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = ProteinMpnnInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.proteinmpnn].output_basename",
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
                    "[bio.proteinmpnn].script `{}` not found (resolved {})",
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
                        "[bio.proteinmpnn].script path `{}` has no filename",
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
                    "[bio.proteinmpnn].input_pdb `{}` not found (resolved {})",
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
                        "[bio.proteinmpnn].input_pdb path `{}` has no filename",
                        input.input_pdb.display()
                    ),
                })?;
        let dest_pdb = workdir.join(pdb_filename);
        if source_pdb != dest_pdb {
            fs::copy(&source_pdb, &dest_pdb)?;
        }

        // Drop a flat `valenx_params.json` into the workdir so the
        // user's design script can read the parsed
        // `[bio.proteinmpnn]` knobs without having to reparse
        // case.toml itself. Built by hand to avoid pulling in a
        // serde_json dep for a 5-key flat object.
        let params_json = format!(
            "{{\n  \"model_variant\": {},\n  \"temperature\": {},\n  \"num_seq_per_target\": {},\n  \"output_basename\": {},\n  \"input_pdb\": {}\n}}\n",
            json_string(&input.model_variant),
            format_finite_float(input.temperature),
            input.num_seq_per_target,
            json_string(&input.output_basename),
            json_string(&pdb_filename.to_string_lossy()),
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
            // ProteinMPNN inference is much faster than diffusion-
            // based design — typical runs finish in seconds to a few
            // minutes. 30 minutes is a generous default.
            estimated_runtime: Some(Duration::from_secs(30 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting ProteinMPNN", |line| {
            let mut hint = subprocess::Hint::default();
            // Convention: the user-supplied script can emit a sentinel
            // line `[valenx] proteinmpnn done` to signal completion
            // before exit; lift to a 95% progress tick.
            if line.contains("[valenx] proteinmpnn done") {
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
            "ProteinMPNN",
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
        let params = read_params(&job.workdir);

        if let Some(p) = pdb_path {
            artefacts.push(Artifact {
                path: p,
                kind: ArtifactKind::Other,
                checksum: None,
                label: "ProteinMPNN input PDB".to_string(),
            });
        }
        if let Some(p) = script_path {
            artefacts.push(Artifact {
                path: p,
                kind: ArtifactKind::Other,
                checksum: None,
                label: "ProteinMPNN script".to_string(),
            });
        }

        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-proteinmpnn", ?e, "workdir read failed");
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
                    let stem_ok = match params.as_ref() {
                        Some(basename) => {
                            let stem = path.file_stem().and_then(|s| s.to_str()).unwrap_or("");
                            stem == basename
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
                    label: "ProteinMPNN log".to_string(),
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
            // Round-22 M2: cap the read at MAX_PDB_FILE_BYTES (256
            // MiB) so a poisoned workdir with a multi-GB `.fa` can't
            // OOM `collect()` before the parser runs.
            let label = match valenx_core::io_caps::read_capped_to_string(
                &path,
                valenx_core::io_caps::MAX_PDB_FILE_BYTES as usize,
            ) {
                Ok(text) => {
                    match valenx_bio::format::fasta::read(&text, valenx_bio::Alphabet::Protein) {
                        Ok(seqs) => format!(
                            "ProteinMPNN \u{00b7} {} sequence{}",
                            seqs.len(),
                            if seqs.len() == 1 { "" } else { "s" }
                        ),
                        Err(_) => "ProteinMPNN designed sequences".to_string(),
                    }
                }
                Err(_) => "ProteinMPNN designed sequences".to_string(),
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
            ribbon_contributions: vec!["bio.proteinmpnn.design"],
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
        // The default Display for f64 emits "0" for integer-valued
        // floats, which JSON parsers happily accept as a number.
        format!("{f}")
    } else {
        // Defensive: emit 0 rather than NaN/Infinity (which JSON
        // doesn't recognise). Validation above should mean this
        // branch is never hit.
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

/// Run `python -c "import ProteinMPNN; print(ProteinMPNN.__version__)"`
/// and parse a `semver::Version` out of stdout. Returns `None` on any
/// failure (interpreter unusable, ProteinMPNN not importable, version
/// string malformed); `probe()` falls back to a "ProteinMPNN not
/// importable" warning in that case.
fn detect_proteinmpnn_version(python_binary: &Path) -> Option<Version> {
    let output = std::process::Command::new(python_binary)
        .arg("-c")
        .arg("import ProteinMPNN; print(ProteinMPNN.__version__)")
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
        let info = ProteinMpnnAdapter::new().info();
        assert_eq!(info.id, "proteinmpnn");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "MIT");
        assert_eq!(info.display_name, "ProteinMPNN");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = ProteinMpnnAdapter::new().info();
        assert_eq!(info.version_range.min_inclusive, Version::new(1, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(2, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = ProteinMpnnAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.proteinmpnn.design"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = ProteinMpnnAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    /// `collect()` should walk the workdir for the configured
    /// `<output_basename>.fa`, parse it via the canonical FASTA
    /// reader, and surface it as a `Tabular` artifact with the
    /// "ProteinMPNN \u{00b7} N sequences" label.
    #[test]
    fn collect_walks_workdir_and_classifies_outputs() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-proteinmpnn-collect-{}",
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
            "{\n  \"model_variant\": \"vanilla\",\n  \"temperature\": 0.1,\n  \"num_seq_per_target\": 2,\n  \"output_basename\": \"design\",\n  \"input_pdb\": \"backbone.pdb\"\n}\n",
        )
        .unwrap();
        fs::write(tmp.join("backbone.pdb"), b"ATOM      1  N   ALA A   1\n").unwrap();
        fs::write(tmp.join("design.py"), b"# placeholder").unwrap();
        fs::write(tmp.join("design.fa"), SAMPLE_FASTA).unwrap();
        // A stray fasta that shouldn't be picked up as the design.
        fs::write(tmp.join("unrelated.fa"), SAMPLE_FASTA).unwrap();
        fs::write(tmp.join("run.log"), b"ProteinMPNN run log\n").unwrap();

        let job = PreparedJob {
            workdir: tmp.clone(),
            native_command: vec![],
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: true,
        };
        let results = ProteinMpnnAdapter::new().collect(&job).unwrap();

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
            design.label.contains("ProteinMPNN"),
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
        assert_eq!(input_art.label, "ProteinMPNN input PDB");

        let py_art = results
            .artifacts
            .iter()
            .find(|a| a.path.extension().is_some_and(|e| e == "py"))
            .expect("script artifact present");
        assert_eq!(py_art.kind, ArtifactKind::Other);
        assert_eq!(py_art.label, "ProteinMPNN script");

        let log_art = results
            .artifacts
            .iter()
            .find(|a| a.path.extension().is_some_and(|e| e == "log"))
            .expect("log artifact present");
        assert_eq!(log_art.kind, ArtifactKind::Log);

        let _ = fs::remove_dir_all(&tmp);
    }

    /// A malformed FASTA shouldn't crash collect — it should degrade
    /// to the generic "ProteinMPNN designed sequences" label so the
    /// UI still surfaces the raw file.
    #[test]
    fn collect_fasta_parse_failure_degrades_gracefully() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-proteinmpnn-bad-fa-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        fs::write(
            tmp.join("valenx_params.json"),
            "{\n  \"model_variant\": \"vanilla\",\n  \"temperature\": 0.1,\n  \"num_seq_per_target\": 1,\n  \"output_basename\": \"design\",\n  \"input_pdb\": \"backbone.pdb\"\n}\n",
        )
        .unwrap();
        // Body contains characters outside the protein alphabet.
        fs::write(tmp.join("design.fa"), b">x\n!!!!\n").unwrap();

        let job = PreparedJob {
            workdir: tmp.clone(),
            native_command: vec![],
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: true,
        };
        let results = ProteinMpnnAdapter::new().collect(&job).unwrap();
        let fa_art = results
            .artifacts
            .iter()
            .find(|a| a.path.extension().is_some_and(|e| e == "fa"))
            .expect("artifact still surfaced");
        assert_eq!(fa_art.kind, ArtifactKind::Tabular);
        assert_eq!(fa_art.label, "ProteinMPNN designed sequences");
        let _ = fs::remove_dir_all(&tmp);
    }
}
