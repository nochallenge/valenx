//! # valenx-adapter-chroma
//!
//! Adapter for [Chroma](https://github.com/generatebio/chroma) — Generate
//! Biomedicines' diffusion-based protein design model. Chroma jointly
//! samples backbone geometry + sequence under a structure-conditioned
//! diffusion process, exposing knobs over sample count, design length,
//! and temperature (akin to a softmax-temperature trade-off between
//! likelihood and diversity).
//!
//! **Phase 27.5 — subprocess wrapper for user-provided scripts.** The
//! user supplies a `design_chroma.py` (or whatever filename) referenced
//! from `[bio.chroma].script` in `case.toml`. `prepare()` stages the
//! script into the workdir and `run()` invokes `python <script>` via
//! the shared subprocess runner. The script is responsible for
//! invoking Chroma with the parsed knobs and writing PDB outputs under
//! `<output_basename>*.pdb` (and optionally `<output_basename>*.fa`
//! for the paired sequence).
//!
//! ## `valenx_params.json`
//!
//! Chroma's Python API surface evolves with each weights drop — pinning
//! a specific call shape would break on every upstream release.
//! Instead, `prepare()` writes a flat JSON file `valenx_params.json`
//! into the workdir containing the parsed `[bio.chroma]` knobs:
//!
//! ```json
//! {
//!   "num_samples":     4,
//!   "length":          100,
//!   "temperature":     1.0,
//!   "output_basename": "design"
//! }
//! ```
//!
//! User scripts read it with `json.load(open("valenx_params.json"))`
//! and pass the values through to Chroma themselves.
//!
//! On `collect()` we walk the workdir for `<output_basename>*.pdb`
//! design backbones (surfaced as `Native` artifacts) and any
//! `<output_basename>*.fa` paired sequences (surfaced as `Tabular`).

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

use crate::case_input::ChromaInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(ChromaAdapter::new())
}

pub struct ChromaAdapter;

impl ChromaAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for ChromaAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "chroma";
/// Python interpreter candidates. `python3` first because on Linux
/// `python` may still be Python 2 on legacy distros; on Windows
/// `python` typically resolves to the Windows Store / 3.x install.
const PYTHON_BINARIES: &[&str] = &["python3", "python"];

impl Adapter for ChromaAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "Chroma",
            // Chroma's first tagged release is the 1.0 weights /
            // inference code drop. Upper bound 2.0 reserves room for
            // an upcoming major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(1, 0, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "Apache-2.0",
            docs_url: "https://github.com/generatebio/chroma",
            homepage_url: "https://generatebiomedicines.com/chroma",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(PYTHON_BINARIES) {
            Some(binary_path) => {
                // Try `import chroma; print(chroma.__version__)` —
                // confirms the package is importable from the chosen
                // interpreter (vs. just having Python on PATH). Some
                // Chroma checkouts don't expose `__version__` — fall
                // back to a "couldn't import" warning so the probe
                // still surfaces a useful state.
                let found_version = detect_chroma_version(&binary_path);
                let mut warnings = Vec::new();
                if found_version.is_none() {
                    warnings.push(
                        "probe found `python` on PATH but could not import \
                         `chroma` — install Chroma from \
                         https://github.com/generatebio/chroma and ensure \
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
                hint: "Python 3.9+ with Chroma installed; clone \
                       https://github.com/generatebio/chroma and follow \
                       the install steps after ensuring python3 is on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = ChromaInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.chroma].output_basename",
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
                    "[bio.chroma].script `{}` not found (resolved {})",
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
                        "[bio.chroma].script path `{}` has no filename",
                        input.script.display()
                    ),
                })?;
        let dest_script = workdir.join(script_filename);
        if source_script != dest_script {
            fs::copy(&source_script, &dest_script)?;
        }

        // Drop a flat `valenx_params.json` into the workdir so the
        // user's design script can read the parsed `[bio.chroma]`
        // knobs without having to reparse case.toml itself. Built by
        // hand to avoid pulling in a serde_json dep for a 4-key flat
        // object.
        let params_json = format!(
            "{{\n  \"num_samples\": {},\n  \"length\": {},\n  \"temperature\": {},\n  \"output_basename\": {}\n}}\n",
            input.num_samples,
            input.length,
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
            // Chroma sampling on a consumer GPU runs minutes-to-hours
            // depending on `num_samples` / `length`. 4 hours mirrors
            // the RFdiffusion default — generous enough that long
            // runs aren't pre-empted.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting Chroma", |line| {
            let mut hint = subprocess::Hint::default();
            // Convention: the user-supplied script can emit a sentinel
            // line `[valenx] chroma done` to signal completion before
            // exit; lift to a 95% progress tick.
            if line.contains("[valenx] chroma done") {
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
        // Provenance: hash the staged Python script (Chroma takes no
        // input PDB). Falls back to case.toml when the script isn't
        // present yet.
        let script_path = first_script_in_workdir(&job.workdir);
        let case_hash_input = script_path
            .clone()
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "Chroma",
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
        // surfaces every PDB / FASTA it can find.
        let basename = read_params(&job.workdir);

        if let Some(p) = script_path {
            artefacts.push(Artifact {
                path: p,
                kind: ArtifactKind::Other,
                checksum: None,
                label: "Chroma script".to_string(),
            });
        }

        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-chroma", ?e, "workdir read failed");
                return Ok(results);
            }
        };
        let mut design_paths: Vec<PathBuf> = Vec::new();
        let mut sequence_paths: Vec<PathBuf> = Vec::new();
        for entry in entries.flatten() {
            let path = entry.path();
            if !path.is_file() {
                continue;
            }
            let ext = path
                .extension()
                .and_then(|s| s.to_str())
                .map(|s| s.to_ascii_lowercase());
            // Stem-match against `output_basename` so foreign PDBs / FAs
            // don't pollute the results when params are present.
            let stem_ok = match basename.as_ref() {
                Some(b) => {
                    let stem = path.file_stem().and_then(|s| s.to_str()).unwrap_or("");
                    stem.starts_with(b.as_str())
                }
                None => true,
            };
            match ext.as_deref() {
                Some("pdb") => {
                    if stem_ok {
                        design_paths.push(path);
                    }
                }
                Some("fa") | Some("fasta") | Some("faa") => {
                    if stem_ok {
                        sequence_paths.push(path);
                    }
                }
                Some("log") => artefacts.push(Artifact {
                    path,
                    kind: ArtifactKind::Log,
                    checksum: None,
                    label: "Chroma log".to_string(),
                }),
                _ => continue,
            }
        }
        design_paths.sort();
        for path in design_paths {
            let stem = path
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("design")
                .to_string();
            // Round-22 M2: cap the per-PDB read at MAX_PDB_FILE_BYTES
            // (256 MiB) so a poisoned workdir with a multi-GB `.pdb`
            // can't OOM `collect()` before the parser runs.
            let label = match valenx_core::io_caps::read_capped_to_string(
                &path,
                valenx_core::io_caps::MAX_PDB_FILE_BYTES as usize,
            ) {
                Ok(text) => match valenx_bio::format::pdb::read(&stem, &text) {
                    Ok(structure) => format!(
                        "Chroma design `{}` ({} atoms, {} residues)",
                        stem,
                        structure.atom_count(),
                        structure.residue_count()
                    ),
                    Err(_) => "Chroma design".to_string(),
                },
                Err(_) => "Chroma design".to_string(),
            };
            artefacts.push(Artifact {
                path,
                kind: ArtifactKind::Native,
                checksum: None,
                label,
            });
        }
        sequence_paths.sort();
        for path in sequence_paths {
            artefacts.push(Artifact {
                path,
                kind: ArtifactKind::Tabular,
                checksum: None,
                label: "Chroma sequence".to_string(),
            });
        }

        artefacts.sort_by(|a, b| a.path.cmp(&b.path));
        results.artifacts = artefacts;
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.chroma.design"],
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

/// Run `python -c "import chroma; print(chroma.__version__)"`
/// and parse a `semver::Version` out of stdout. Returns `None` on any
/// failure (interpreter unusable, chroma not importable, version
/// string malformed); `probe()` falls back to a "chroma not
/// importable" warning in that case.
fn detect_chroma_version(python_binary: &Path) -> Option<Version> {
    let output = std::process::Command::new(python_binary)
        .arg("-c")
        .arg("import chroma; print(chroma.__version__)")
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

    /// Minimal valid PDB record covering one residue.
    const SAMPLE_PDB: &str = "\
ATOM      1  N   ALA A   1      11.104  13.207   2.063  1.00  0.00           N
ATOM      2  CA  ALA A   1      11.804  13.793   3.215  1.00  0.00           C
ATOM      3  C   ALA A   1      11.072  15.058   3.668  1.00  0.00           C
ATOM      4  O   ALA A   1       9.835  15.117   3.586  1.00  0.00           O
ATOM      5  CB  ALA A   1      11.916  12.789   4.357  1.00  0.00           C
END
";

    #[test]
    fn info_is_bio_domain() {
        let info = ChromaAdapter::new().info();
        assert_eq!(info.id, "chroma");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "Apache-2.0");
        assert_eq!(info.display_name, "Chroma");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = ChromaAdapter::new().info();
        assert_eq!(info.version_range.min_inclusive, Version::new(1, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(2, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = ChromaAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.chroma.design"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = ChromaAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    /// `collect()` should walk the workdir for `<output_basename>*.pdb`
    /// design files (Native) and `<output_basename>*.fa` paired
    /// sequences (Tabular).
    #[test]
    fn collect_walks_workdir_and_classifies_outputs() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-chroma-collect-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        fs::write(
            tmp.join("valenx_params.json"),
            "{\n  \"num_samples\": 2,\n  \"length\": 100,\n  \"temperature\": 1.0,\n  \"output_basename\": \"design\"\n}\n",
        )
        .unwrap();
        fs::write(tmp.join("design.py"), b"# placeholder").unwrap();
        fs::write(tmp.join("design_0.pdb"), SAMPLE_PDB).unwrap();
        fs::write(tmp.join("design_1.pdb"), SAMPLE_PDB).unwrap();
        fs::write(tmp.join("design_0.fa"), b">design_0\nACDEF\n").unwrap();
        // A stray pdb that shouldn't be picked up as a design.
        fs::write(tmp.join("unrelated.pdb"), SAMPLE_PDB).unwrap();
        fs::write(tmp.join("run.log"), b"Chroma run log\n").unwrap();

        let job = PreparedJob {
            workdir: tmp.clone(),
            native_command: vec![],
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: true,
        };
        let results = ChromaAdapter::new().collect(&job).unwrap();

        let designs: Vec<_> = results
            .artifacts
            .iter()
            .filter(|a| a.kind == ArtifactKind::Native)
            .collect();
        assert_eq!(
            designs.len(),
            2,
            "expected 2 design PDBs, got {}: {:?}",
            designs.len(),
            results.artifacts
        );
        for d in &designs {
            assert!(d.label.contains("Chroma design"), "label was: {}", d.label);
        }

        let sequences: Vec<_> = results
            .artifacts
            .iter()
            .filter(|a| a.kind == ArtifactKind::Tabular)
            .collect();
        assert_eq!(sequences.len(), 1);
        assert_eq!(sequences[0].label, "Chroma sequence");

        let py_art = results
            .artifacts
            .iter()
            .find(|a| a.path.extension().is_some_and(|e| e == "py"))
            .expect("script artifact present");
        assert_eq!(py_art.kind, ArtifactKind::Other);
        assert_eq!(py_art.label, "Chroma script");

        let log_art = results
            .artifacts
            .iter()
            .find(|a| a.path.extension().is_some_and(|e| e == "log"))
            .expect("log artifact present");
        assert_eq!(log_art.kind, ArtifactKind::Log);

        let _ = fs::remove_dir_all(&tmp);
    }
}
