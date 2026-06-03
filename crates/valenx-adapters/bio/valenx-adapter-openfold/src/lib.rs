//! # valenx-adapter-openfold
//!
//! Adapter for [OpenFold](https://github.com/aqlaboratory/openfold) —
//! the open-source PyTorch reimplementation of AlphaFold 2 by the
//! AQ Lab. Same protein-structure-prediction job as ColabFold + AF2,
//! but the entire pipeline is in PyTorch (faster on modern GPUs;
//! easier to fine-tune). Supports the full AF2 model preset family
//! (`model_1` … `model_5_ptm`) plus the multimer v3 weights for
//! complex prediction.
//!
//! **Phase 17.5 — subprocess wrapper for user-provided scripts.** The
//! user supplies a `predict_openfold.py` (or whatever filename)
//! referenced from `[bio.openfold].script` in `case.toml` plus a
//! FASTA query and a `model_preset` choice. `prepare()` stages the
//! script + FASTA into the workdir and `run()` invokes
//! `python <script> <fasta>` via the shared subprocess runner.
//! The script is responsible for invoking
//! `openfold.run_pretrained_openfold` with the chosen preset and
//! writing PDB outputs.
//!
//! ## `valenx_params.json`
//!
//! OpenFold's `run_pretrained_openfold.py` evolves rapidly and there
//! is no stable command-line shape we can pin. Instead, `prepare()`
//! writes a flat JSON file `valenx_params.json` into the workdir
//! containing the parsed `[bio.openfold]` knobs:
//!
//! ```json
//! {
//!   "model_preset":  "model_1_ptm",
//!   "use_templates": false,
//!   "num_recycles":  3,
//!   "query_fasta":   "query.fasta"
//! }
//! ```
//!
//! User scripts read it with `json.load(open("valenx_params.json"))`
//! and pass the values through to OpenFold themselves. This keeps
//! the adapter free of upstream API churn and means `case.toml`
//! knobs actually reach the model.
//!
//! On `collect()` we walk the workdir for `*.pdb` files and parse
//! each via [`valenx_bio::format::pdb::read`]. Each is surfaced as a
//! typed [`ArtifactKind::Native`] artifact with an "OpenFold
//! prediction" label.

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

use crate::case_input::OpenFoldInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(OpenFoldAdapter::new())
}

pub struct OpenFoldAdapter;

impl OpenFoldAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for OpenFoldAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "openfold";
/// Python interpreter candidates. `python3` first because on Linux
/// `python` may still be Python 2 on legacy distros; on Windows
/// `python` typically resolves to the Windows Store / 3.x install.
const PYTHON_BINARIES: &[&str] = &["python3", "python"];

impl Adapter for OpenFoldAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "OpenFold",
            // OpenFold 1.x is the first stable release with the
            // `run_pretrained_openfold` surface we lean on. Upper
            // bound 2.0 reserves room for an upcoming major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(1, 0, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "Apache-2.0",
            docs_url: "https://openfold.readthedocs.io/",
            homepage_url: "https://github.com/aqlaboratory/openfold",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(PYTHON_BINARIES) {
            Some(binary_path) => {
                // Try `import openfold; print(openfold.__version__)`
                // first. Some OpenFold builds don't ship `__version__`
                // — fall back gracefully to "couldn't import" so the
                // probe still surfaces a useful state.
                let found_version = detect_openfold_version(&binary_path);
                let mut warnings = Vec::new();
                if found_version.is_none() {
                    warnings.push(
                        "probe found `python` on PATH but could not import \
                         `openfold` — install OpenFold from \
                         https://github.com/aqlaboratory/openfold and ensure \
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
                hint: "Python 3.9+ with OpenFold installed; clone \
                       https://github.com/aqlaboratory/openfold and follow \
                       the install steps after ensuring python3 is on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = OpenFoldInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Stage the user-supplied Python script. Resolved against the
        // case directory; `confined_join` rejects absolute paths and
        // `..` traversal so the staged copy stays confined to the
        // case directory.
        let source_script = confined_join(&case.path, &input.script)?;
        if !source_script.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.openfold].script `{}` not found (resolved {})",
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
                        "[bio.openfold].script path `{}` has no filename",
                        input.script.display()
                    ),
                })?;
        let dest_script = workdir.join(script_filename);
        if source_script != dest_script {
            fs::copy(&source_script, &dest_script)?;
        }

        // Stage the FASTA query alongside the script.
        let source_fasta = confined_join(&case.path, &input.query_fasta)?;
        if !source_fasta.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.openfold].query_fasta `{}` not found (resolved {})",
                    input.query_fasta.display(),
                    source_fasta.display()
                ),
            });
        }
        let fasta_filename =
            input
                .query_fasta
                .file_name()
                .ok_or_else(|| AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.openfold].query_fasta path `{}` has no filename",
                        input.query_fasta.display()
                    ),
                })?;
        let dest_fasta = workdir.join(fasta_filename);
        if source_fasta != dest_fasta {
            fs::copy(&source_fasta, &dest_fasta)?;
        }

        // Drop a flat `valenx_params.json` into the workdir so the
        // user's predict script can read the parsed `[bio.openfold]`
        // knobs without having to reparse case.toml itself. OpenFold's
        // CLI surface evolves rapidly and there's no stable place to
        // pass these via flags. Built by hand to avoid pulling in a
        // serde_json dep for a 4-key flat object.
        let params_json = format!(
            "{{\n  \"model_preset\": {},\n  \"use_templates\": {},\n  \"num_recycles\": {},\n  \"query_fasta\": {}\n}}\n",
            json_string(&input.model_preset),
            input.use_templates,
            input.num_recycles,
            json_string(&fasta_filename.to_string_lossy()),
        );
        valenx_core::io_caps::atomic_write_str(&workdir.join("valenx_params.json"), &params_json)?;

        // Resolve the Python binary.
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

        // Pass the FASTA filename as the script's first positional
        // argument so the script can find its input by name.
        let native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from(script_filename),
            OsString::from(fasta_filename),
        ];

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // OpenFold runtime varies enormously by sequence length
            // and GPU. 4 hours is a generous default; long runs
            // override through their own progress reporting.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting OpenFold", |line| {
            let mut hint = subprocess::Hint::default();
            // Convention: the user-supplied script can emit a sentinel
            // line `[valenx] openfold done` to signal completion before
            // exit; lift to a 95% progress tick.
            if line.contains("[valenx] openfold done") {
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
        // Provenance: hash the staged FASTA. Falls back to the script,
        // then case.toml, when the FASTA isn't present yet.
        let fasta_path = first_fasta_in_workdir(&job.workdir);
        let script_path = first_script_in_workdir(&job.workdir);
        let case_hash_input = fasta_path
            .clone()
            .or_else(|| script_path.clone())
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "OpenFold",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        if let Some(p) = fasta_path {
            artefacts.push(Artifact {
                path: p,
                kind: ArtifactKind::Other,
                checksum: None,
                label: "OpenFold input FASTA".to_string(),
            });
        }
        if let Some(p) = script_path {
            artefacts.push(Artifact {
                path: p,
                kind: ArtifactKind::Other,
                checksum: None,
                label: "OpenFold script".to_string(),
            });
        }

        // Walk the workdir top-level + a `predictions/` subdirectory
        // (OpenFold's customary output location). Soft-validate each
        // PDB via the canonical reader.
        let mut pdb_paths: Vec<PathBuf> = Vec::new();
        for dir in [job.workdir.clone(), job.workdir.join("predictions")] {
            let entries = match fs::read_dir(&dir) {
                Ok(e) => e,
                Err(_) => continue,
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
                match ext.as_deref() {
                    Some("pdb") => pdb_paths.push(path),
                    Some("log") => artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Log,
                        checksum: None,
                        label: "OpenFold log".to_string(),
                    }),
                    _ => continue,
                }
            }
        }
        pdb_paths.sort();
        for path in pdb_paths {
            let stem = path
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("prediction")
                .to_string();
            // Round-23 named finding: bound the per-PDB collect-label
            // read at MAX_PDB_FILE_BYTES (256 MiB) so a poisoned or
            // runaway prediction can't OOM the renderer before the
            // parser sees the first ATOM line. Same cap as the R22
            // bio-collect-label sweep.
            let label = match valenx_core::io_caps::read_capped_to_string(
                &path,
                valenx_core::io_caps::MAX_PDB_FILE_BYTES as usize,
            ) {
                Ok(text) => match valenx_bio::format::pdb::read(&stem, &text) {
                    Ok(structure) => format!(
                        "OpenFold prediction `{}` ({} atoms, {} residues)",
                        stem,
                        structure.atom_count(),
                        structure.residue_count()
                    ),
                    Err(e) => format!(
                        "OpenFold prediction `{}` (parse warning: {})",
                        stem,
                        e.to_string().lines().next().unwrap_or("invalid")
                    ),
                },
                Err(_) => format!("OpenFold prediction `{stem}`"),
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
            ribbon_contributions: vec!["bio.openfold.predict"],
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
/// hashing.
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

/// Run `python -c "import openfold; print(openfold.__version__)"` and
/// parse a `semver::Version` out of stdout. Returns `None` on any
/// failure (interpreter unusable, openfold not importable, version
/// string malformed); `probe()` falls back to a "openfold not
/// importable" warning in that case.
fn detect_openfold_version(python_binary: &Path) -> Option<Version> {
    let output = std::process::Command::new(python_binary)
        .arg("-c")
        .arg("import openfold; print(openfold.__version__)")
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
        let info = OpenFoldAdapter::new().info();
        assert_eq!(info.id, "openfold");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "Apache-2.0");
        assert_eq!(info.display_name, "OpenFold");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = OpenFoldAdapter::new().info();
        assert_eq!(info.version_range.min_inclusive, Version::new(1, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(2, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = OpenFoldAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.openfold.predict"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = OpenFoldAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    /// `collect()` must walk the workdir + the customary `predictions/`
    /// subdir for PDB outputs and surface them with parsed labels.
    #[test]
    fn collect_walks_workdir_and_predictions_subdir() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-openfold-collect-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        fs::write(tmp.join("query.fasta"), b">x\nACDEF\n").unwrap();
        fs::write(tmp.join("predict.py"), b"# placeholder").unwrap();
        let predictions = tmp.join("predictions");
        fs::create_dir_all(&predictions).unwrap();
        fs::write(predictions.join("model_1_ptm.pdb"), SAMPLE_PDB).unwrap();
        fs::write(tmp.join("run.log"), b"OpenFold run log\n").unwrap();

        let job = PreparedJob {
            workdir: tmp.clone(),
            native_command: vec![],
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: true,
        };
        let results = OpenFoldAdapter::new().collect(&job).unwrap();

        let pdb_art = results
            .artifacts
            .iter()
            .find(|a| a.path.extension().is_some_and(|e| e == "pdb"))
            .expect("PDB artifact present from predictions/ subdir");
        assert_eq!(pdb_art.kind, ArtifactKind::Native);
        assert!(
            pdb_art.label.contains("OpenFold prediction"),
            "label was: {}",
            pdb_art.label
        );
        assert!(
            pdb_art.label.contains("5 atoms") && pdb_art.label.contains("1 residues"),
            "label was: {}",
            pdb_art.label
        );

        let fasta_art = results
            .artifacts
            .iter()
            .find(|a| a.path.extension().is_some_and(|e| e == "fasta"))
            .expect("FASTA artifact present");
        assert_eq!(fasta_art.kind, ArtifactKind::Other);
        assert_eq!(fasta_art.label, "OpenFold input FASTA");

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
