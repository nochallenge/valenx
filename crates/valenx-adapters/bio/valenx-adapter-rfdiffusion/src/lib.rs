//! # valenx-adapter-rfdiffusion
//!
//! Adapter for [RFdiffusion](https://github.com/RosettaCommons/RFdiffusion)
//! — RosettaCommons' diffusion-model-based de novo protein backbone
//! generator. Where ColabFold / ESMFold / OpenFold / AlphaFold predict
//! the structure of a known sequence, RFdiffusion goes the other way:
//! it samples novel protein backbones conditioned on the user's design
//! intent (motif scaffolding, binder design against a target,
//! unconditional sampling, or partial diffusion of an existing
//! backbone).
//!
//! **Phase 27 — subprocess wrapper for user-provided scripts.** The
//! user supplies a `design_rfdiffusion.py` (or whatever filename)
//! referenced from `[bio.rfdiffusion].script` in `case.toml` plus an
//! input PDB and a `mode` choice. `prepare()` stages the script + PDB
//! into the workdir and `run()` invokes `python <script>` via the
//! shared subprocess runner. The script is responsible for invoking
//! RFdiffusion's `inference.run_inference` (or the underlying Hydra
//! config) with the parsed knobs and writing PDB outputs under
//! `<output_basename>*.pdb`.
//!
//! ## `valenx_params.json`
//!
//! RFdiffusion's CLI is a Hydra config tree — there's no stable flag
//! shape we can pin without breaking on the next refactor. Instead,
//! `prepare()` writes a flat JSON file `valenx_params.json` into the
//! workdir containing the parsed `[bio.rfdiffusion]` knobs:
//!
//! ```json
//! {
//!   "mode":            "motif",
//!   "num_designs":     8,
//!   "diffusion_steps": 50,
//!   "output_basename": "design",
//!   "input_pdb":       "scaffold.pdb"
//! }
//! ```
//!
//! User scripts read it with `json.load(open("valenx_params.json"))`
//! and pass the values through to RFdiffusion themselves. This keeps
//! the adapter free of upstream API churn and means `case.toml`
//! knobs actually reach the model.
//!
//! On `collect()` we walk the workdir for `<output_basename>*.pdb`
//! files and parse each via [`valenx_bio::format::pdb::read`]. Each is
//! surfaced as a typed [`ArtifactKind::Native`] artifact with an
//! "RFdiffusion design" label.

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

use crate::case_input::RfDiffusionInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(RfDiffusionAdapter::new())
}

pub struct RfDiffusionAdapter;

impl RfDiffusionAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for RfDiffusionAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "rfdiffusion";
/// Python interpreter candidates. `python3` first because on Linux
/// `python` may still be Python 2 on legacy distros; on Windows
/// `python` typically resolves to the Windows Store / 3.x install.
const PYTHON_BINARIES: &[&str] = &["python3", "python"];

impl Adapter for RfDiffusionAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "RFdiffusion",
            // RFdiffusion's first tagged release is the 1.1 weights /
            // inference code drop. Upper bound 2.0 reserves room for
            // an upcoming major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(1, 1, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "BSD-3-Clause",
            docs_url: "https://github.com/RosettaCommons/RFdiffusion",
            homepage_url: "https://github.com/RosettaCommons/RFdiffusion",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(PYTHON_BINARIES) {
            Some(binary_path) => {
                // Try `import rfdiffusion; print(rfdiffusion.__version__)`
                // — confirms the package is importable from the chosen
                // interpreter (vs. just having Python on PATH). Some
                // RFdiffusion checkouts don't expose `__version__` —
                // fall back to a "couldn't import" warning so the probe
                // still surfaces a useful state.
                let found_version = detect_rfdiffusion_version(&binary_path);
                let mut warnings = Vec::new();
                if found_version.is_none() {
                    warnings.push(
                        "probe found `python` on PATH but could not import \
                         `rfdiffusion` — install RFdiffusion from \
                         https://github.com/RosettaCommons/RFdiffusion and \
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
                hint: "Python 3.9+ with RFdiffusion installed; clone \
                       https://github.com/RosettaCommons/RFdiffusion and \
                       follow the install steps after ensuring python3 is \
                       on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = RfDiffusionInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.rfdiffusion].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Stage the user-supplied Python script. Resolved against the
        // case directory; absolute paths and `..` traversal are
        // rejected by `confined_join` so the staged copy stays
        // confined to the case directory.
        let source_script = confined_join(&case.path, &input.script)?;
        if !source_script.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.rfdiffusion].script `{}` not found (resolved {})",
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
                        "[bio.rfdiffusion].script path `{}` has no filename",
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
                    "[bio.rfdiffusion].input_pdb `{}` not found (resolved {})",
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
                        "[bio.rfdiffusion].input_pdb path `{}` has no filename",
                        input.input_pdb.display()
                    ),
                })?;
        let dest_pdb = workdir.join(pdb_filename);
        if source_pdb != dest_pdb {
            fs::copy(&source_pdb, &dest_pdb)?;
        }

        // Drop a flat `valenx_params.json` into the workdir so the
        // user's design script can read the parsed `[bio.rfdiffusion]`
        // knobs without having to reparse case.toml itself. Built by
        // hand to avoid pulling in a serde_json dep for a 5-key flat
        // object.
        let params_json = format!(
            "{{\n  \"mode\": {},\n  \"num_designs\": {},\n  \"diffusion_steps\": {},\n  \"output_basename\": {},\n  \"input_pdb\": {}\n}}\n",
            json_string(&input.mode),
            input.num_designs,
            input.diffusion_steps,
            json_string(&input.output_basename),
            json_string(&pdb_filename.to_string_lossy()),
        );
        valenx_core::io_caps::atomic_write_str(&workdir.join("valenx_params.json"), &params_json)?;

        // Resolve the Python binary. Same logic as every other
        // Phase 17 Python-script adapter (Biopython / OpenMM / RDKit /
        // MDAnalysis / OpenFold / ESMFold): bare `python` / `python3`
        // walks PATH; absolute paths or pinned interpreters are
        // honored verbatim.
        // Round-3 security fix (round-12 sweep): validate + resolve
        // via the shared helper.
        let binary_path =
            valenx_core::adapter_helpers::resolve_python_binary(&input.python, PYTHON_BINARIES)
                .map_err(|e| AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!("[bio.rfdiffusion].python: {e}"),
                })?;

        let native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from(script_filename),
        ];

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // RFdiffusion sampling can take hours per design on a
            // consumer GPU; with `num_designs` typically 8 – 64 the
            // tail is long. 4 hours is a generous default; long runs
            // override through their own progress reporting.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting RFdiffusion", |line| {
            let mut hint = subprocess::Hint::default();
            // Convention: the user-supplied script can emit a sentinel
            // line `[valenx] rfdiffusion done` to signal completion
            // before exit; lift to a 95% progress tick.
            if line.contains("[valenx] rfdiffusion done") {
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
            "RFdiffusion",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Read the staged params back out so we can restrict the
        // collected output PDBs to those matching `output_basename`
        // and label the input PDB distinctly. Failure to read the
        // params is non-fatal — collect still surfaces every PDB.
        let params = read_params(&job.workdir);

        if let Some(p) = pdb_path.clone() {
            artefacts.push(Artifact {
                path: p,
                kind: ArtifactKind::Other,
                checksum: None,
                label: "RFdiffusion input PDB".to_string(),
            });
        }
        if let Some(p) = script_path {
            artefacts.push(Artifact {
                path: p,
                kind: ArtifactKind::Other,
                checksum: None,
                label: "RFdiffusion script".to_string(),
            });
        }

        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-rfdiffusion", ?e, "workdir read failed");
                return Ok(results);
            }
        };
        let mut design_paths: Vec<PathBuf> = Vec::new();
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
                Some("pdb") => {
                    // Skip the input PDB — it's already surfaced above.
                    if Some(&path) == pdb_path.as_ref() {
                        continue;
                    }
                    // Restrict to designs whose stem starts with the
                    // configured `output_basename`. If params couldn't
                    // be read, accept everything (best-effort).
                    let stem_ok = match params.as_ref() {
                        Some((basename, _input_pdb_name)) => {
                            let stem = path.file_stem().and_then(|s| s.to_str()).unwrap_or("");
                            stem.starts_with(basename)
                        }
                        None => true,
                    };
                    if stem_ok {
                        design_paths.push(path);
                    }
                }
                Some("log") => artefacts.push(Artifact {
                    path,
                    kind: ArtifactKind::Log,
                    checksum: None,
                    label: "RFdiffusion log".to_string(),
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
                        "RFdiffusion design `{}` ({} atoms, {} residues)",
                        stem,
                        structure.atom_count(),
                        structure.residue_count()
                    ),
                    Err(e) => format!(
                        "RFdiffusion design `{}` (parse warning: {})",
                        stem,
                        e.to_string().lines().next().unwrap_or("invalid")
                    ),
                },
                Err(_) => format!("RFdiffusion design `{stem}`"),
            };
            artefacts.push(Artifact {
                path,
                kind: ArtifactKind::Native,
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
            ribbon_contributions: vec!["bio.rfdiffusion.design"],
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

/// Lift the staged input PDB out of the workdir for provenance hashing.
/// Prefers the PDB referenced by `valenx_params.json` when present; if
/// the params can't be read, falls back to the lexicographically-first
/// `.pdb` file at the top level.
fn first_pdb_in_workdir(workdir: &Path) -> Option<PathBuf> {
    if let Some((_basename, input_pdb_name)) = read_params(workdir) {
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

/// Read `valenx_params.json` enough to recover `(output_basename,
/// input_pdb_filename)` for collect-time filtering. Returns `None` on
/// any failure — callers degrade gracefully.
fn read_params(workdir: &Path) -> Option<(String, String)> {
    let text = valenx_core::io_caps::read_capped_to_string(
        &workdir.join("valenx_params.json"),
        valenx_core::io_caps::MAX_ADAPTER_PARAMS_BYTES as usize,
    )
    .ok()?;
    let basename = extract_json_string(&text, "output_basename")?;
    let input_pdb = extract_json_string(&text, "input_pdb")?;
    Some((basename, input_pdb))
}

/// Pull a flat string field out of our own hand-emitted
/// `valenx_params.json`. Trivially small — we wrote the file
/// ourselves so we know its shape; a full JSON parser would be
/// overkill (and pulling in serde_json just for collect()'s side-band
/// metadata would bloat the dep tree).
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

/// Run `python -c "import rfdiffusion; print(rfdiffusion.__version__)"`
/// and parse a `semver::Version` out of stdout. Returns `None` on any
/// failure (interpreter unusable, rfdiffusion not importable, version
/// string malformed); `probe()` falls back to a "rfdiffusion not
/// importable" warning in that case.
fn detect_rfdiffusion_version(python_binary: &Path) -> Option<Version> {
    let output = std::process::Command::new(python_binary)
        .arg("-c")
        .arg("import rfdiffusion; print(rfdiffusion.__version__)")
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
        let info = RfDiffusionAdapter::new().info();
        assert_eq!(info.id, "rfdiffusion");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "BSD-3-Clause");
        assert_eq!(info.display_name, "RFdiffusion");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = RfDiffusionAdapter::new().info();
        assert_eq!(info.version_range.min_inclusive, Version::new(1, 1, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(2, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = RfDiffusionAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.rfdiffusion.design"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = RfDiffusionAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    /// `collect()` should walk the workdir for `<output_basename>*.pdb`
    /// design files, surface them as Native artifacts with parsed
    /// labels, and tag the input PDB + script as auxiliary `Other`
    /// artifacts.
    #[test]
    fn collect_walks_workdir_and_filters_designs() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-rfdiffusion-collect-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        // Stage the params + input PDB + script + designs.
        fs::write(
            tmp.join("valenx_params.json"),
            "{\n  \"mode\": \"motif\",\n  \"num_designs\": 2,\n  \"diffusion_steps\": 50,\n  \"output_basename\": \"design\",\n  \"input_pdb\": \"scaffold.pdb\"\n}\n",
        )
        .unwrap();
        fs::write(tmp.join("scaffold.pdb"), SAMPLE_PDB).unwrap();
        fs::write(tmp.join("design.py"), b"# placeholder").unwrap();
        fs::write(tmp.join("design_0.pdb"), SAMPLE_PDB).unwrap();
        fs::write(tmp.join("design_1.pdb"), SAMPLE_PDB).unwrap();
        // A stray pdb that shouldn't be picked up as a design.
        fs::write(tmp.join("unrelated.pdb"), SAMPLE_PDB).unwrap();
        fs::write(tmp.join("run.log"), b"RFdiffusion run log\n").unwrap();

        let job = PreparedJob {
            workdir: tmp.clone(),
            native_command: vec![],
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: true,
        };
        let results = RfDiffusionAdapter::new().collect(&job).unwrap();

        // Two design PDBs picked up.
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
            assert!(
                d.label.contains("RFdiffusion design"),
                "label was: {}",
                d.label
            );
            assert!(
                d.label.contains("5 atoms") && d.label.contains("1 residues"),
                "label was: {}",
                d.label
            );
        }

        // Input PDB tagged as Other with the documented label.
        let input_art = results
            .artifacts
            .iter()
            .find(|a| {
                a.path
                    .file_name()
                    .and_then(|s| s.to_str())
                    .map(|s| s == "scaffold.pdb")
                    .unwrap_or(false)
            })
            .expect("input PDB artifact present");
        assert_eq!(input_art.kind, ArtifactKind::Other);
        assert_eq!(input_art.label, "RFdiffusion input PDB");

        let py_art = results
            .artifacts
            .iter()
            .find(|a| a.path.extension().is_some_and(|e| e == "py"))
            .expect("script artifact present");
        assert_eq!(py_art.kind, ArtifactKind::Other);

        let log_art = results
            .artifacts
            .iter()
            .find(|a| a.path.extension().is_some_and(|e| e == "log"))
            .expect("log artifact present");
        assert_eq!(log_art.kind, ArtifactKind::Log);

        let _ = fs::remove_dir_all(&tmp);
    }
}
