//! # valenx-adapter-deepchem
//!
//! Adapter for [DeepChem](https://deepchem.io/) — a PyTorch-backed
//! cheminformatics framework: graph neural networks for molecular
//! property prediction, generative models for molecule design,
//! reinforcement learning over chemical space, and the bread-and-butter
//! tabular pipelines (Morgan / RDKit fingerprints feeding scikit-learn
//! head models).
//!
//! **Phase 24 — subprocess wrapper for user-provided scripts.** The
//! adapter doesn't generate Python; the user supplies a `train.py` (or
//! whatever filename) referenced from `[bio.deepchem].script` in
//! `case.toml`. `prepare()` stages the script (plus optional dataset
//! CSV and checkpoint) into the workdir, drops a `valenx_params.json`
//! with the parsed knobs, and `run()` invokes `python <script>` via
//! the shared subprocess runner.
//!
//! ## `valenx_params.json`
//!
//! DeepChem has no canonical CLI — every site's training script reads
//! its own knobs. Rather than guess at a flag layout, `prepare()`
//! writes a flat JSON file alongside the staged script:
//!
//! ```json
//! {
//!   "smiles":      ["CCO", "c1ccccc1"],
//!   "dataset_csv": "molecules.csv",
//!   "checkpoint":  "best.pt"
//! }
//! ```
//!
//! Scripts read it with `json.load(open("valenx_params.json"))` and
//! pass values through to DeepChem themselves. Same convention as
//! ESMFold / OpenFold.
//!
//! On `collect()` we walk the workdir for `.csv` (analysis output),
//! `.png` (plots), and `.pkl` / `.pt` (model checkpoints).

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

use crate::case_input::DeepChemInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(DeepChemAdapter::new())
}

pub struct DeepChemAdapter;

impl DeepChemAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for DeepChemAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "deepchem";
/// Python interpreter candidates. `python3` first because on Linux
/// `python` may still be Python 2 on legacy distros; on Windows
/// `python` typically resolves to the Windows Store / 3.x install.
const PYTHON_BINARIES: &[&str] = &["python3", "python"];

impl Adapter for DeepChemAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "DeepChem",
            // DeepChem 2.7.x is the first release that targets the
            // PyTorch 2.x line cleanly; upper bound 3.0 reserves room
            // for an eventual API-breaking major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(2, 7, 0),
                max_exclusive: Version::new(3, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "MIT",
            docs_url: "https://deepchem.readthedocs.io/",
            homepage_url: "https://deepchem.io/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(PYTHON_BINARIES) {
            Some(binary_path) => {
                // Try `import deepchem; print(deepchem.__version__)`
                // first — that's the only string that confirms DeepChem
                // itself is installed (vs. just having Python on PATH).
                // If the import fails, return a successful probe with a
                // warning so the user gets a useful install hint.
                let found_version = detect_deepchem_version(&binary_path);
                let mut warnings = Vec::new();
                if found_version.is_none() {
                    warnings.push(
                        "probe found `python` on PATH but could not import \
                         `deepchem` — install DeepChem with `pip install \
                         deepchem` (also requires PyTorch) for runs to \
                         succeed"
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
                hint: "Python 3.9+ with DeepChem installed; `pip install \
                       deepchem` after ensuring python3 is on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = DeepChemInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Stage the user-supplied Python script. `confined_join`
        // rejects absolute paths and `..` traversal so the staged copy
        // stays confined to the case directory.
        let source_script = confined_join(&case.path, &input.script)?;
        if !source_script.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.deepchem].script `{}` not found (resolved {})",
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
                        "[bio.deepchem].script path `{}` has no filename",
                        input.script.display()
                    ),
                })?;
        let dest_script = workdir.join(script_filename);
        if source_script != dest_script {
            fs::copy(&source_script, &dest_script)?;
        }

        // Stage the optional dataset CSV. Resolves against the case
        // directory when relative.
        let dataset_filename: Option<OsString> = if let Some(ref dataset) = input.dataset_csv {
            let source_dataset = confined_join(&case.path, dataset)?;
            if !source_dataset.is_file() {
                return Err(AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.deepchem].dataset_csv `{}` not found (resolved {})",
                        dataset.display(),
                        source_dataset.display()
                    ),
                });
            }
            let fname = dataset
                .file_name()
                .ok_or_else(|| AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.deepchem].dataset_csv path `{}` has no filename",
                        dataset.display()
                    ),
                })?
                .to_os_string();
            let dest = workdir.join(&fname);
            if source_dataset != dest {
                fs::copy(&source_dataset, &dest)?;
            }
            Some(fname)
        } else {
            None
        };

        // Drop `valenx_params.json` so the script can read parsed
        // knobs without reparsing case.toml. Hand-rolled JSON keeps the
        // dep list aligned with the RDKit / ESMFold sibling crates.
        let smiles_json = format_smiles_array(&input.smiles);
        let dataset_field = match &dataset_filename {
            Some(f) => json_string(&f.to_string_lossy()),
            None => "null".to_string(),
        };
        let checkpoint_field = match &input.checkpoint {
            Some(p) => json_string(&p.display().to_string()),
            None => "null".to_string(),
        };
        let params_json = format!(
            "{{\n  \"smiles\": {smiles_json},\n  \"dataset_csv\": {dataset_field},\n  \"checkpoint\": {checkpoint_field}\n}}\n"
        );
        valenx_core::io_caps::atomic_write_str(&workdir.join("valenx_params.json"), &params_json)?;

        // Resolve the Python binary. Same logic as every other Python-
        // script adapter (RDKit / ESMFold / OpenFold / ...): bare
        // `python` / `python3` walks PATH; absolute paths or pinned
        // interpreters are honored verbatim.
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
            // DeepChem workloads vary widely — fingerprint pipelines
            // finish in seconds, GNN training over QM9-class datasets
            // can run for hours on consumer GPUs. 4 hours is a generous
            // default that covers the long tail without being absurd.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting DeepChem", |line| {
            let mut hint = subprocess::Hint::default();
            // Convention: the user-supplied script can emit a sentinel
            // line `[valenx] deepchem done` to signal completion before
            // exit; lift to a 95% progress tick.
            if line.contains("[valenx] deepchem done") {
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
        // Provenance: hash the staged script (the canonical "this case
        // is configured this way" input). We don't know the user's
        // mesh / lock files so leave those empty.
        let script_path = first_script_in_workdir(&job.workdir);
        let case_hash_input = script_path
            .clone()
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "DeepChem",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);

        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-deepchem", ?e, "workdir read failed");
                return Ok(results);
            }
        };
        let mut artefacts: Vec<Artifact> = Vec::new();
        for entry in entries.flatten() {
            let path = entry.path();
            if !path.is_file() {
                continue;
            }
            let ext = path
                .extension()
                .and_then(|s| s.to_str())
                .map(|s| s.to_ascii_lowercase());
            let (kind, label) = match ext.as_deref() {
                Some("csv") => (
                    ArtifactKind::Tabular,
                    "DeepChem analysis output".to_string(),
                ),
                Some("png") => (ArtifactKind::Native, "DeepChem plot".to_string()),
                Some("pkl") | Some("pt") => (
                    ArtifactKind::Native,
                    "DeepChem model checkpoint".to_string(),
                ),
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
            ribbon_contributions: vec!["bio.deepchem.script"],
        }
    }
}

/// Escape a string for embedding inside a JSON string literal. Mirrors
/// the helper used by ESMFold / OpenFold so we don't have to pull in a
/// serde_json dep just to emit a 3-key flat object.
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

/// Render a `Vec<String>` as a JSON array literal of escaped strings.
fn format_smiles_array(smiles: &[String]) -> String {
    if smiles.is_empty() {
        return "[]".to_string();
    }
    let parts: Vec<String> = smiles.iter().map(|s| json_string(s)).collect();
    format!("[{}]", parts.join(", "))
}

/// Lift the staged Python script out of the workdir for provenance
/// hashing. Returns the lexicographically-first `.py` file at the top
/// level, or `None` if none exists yet.
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

/// Run `python -c "import deepchem; print(deepchem.__version__)"` and
/// parse a `semver::Version` out of stdout. Returns `None` on any
/// failure (interpreter unusable, deepchem not importable, version
/// string malformed); `probe()` falls back to a warning in that case.
fn detect_deepchem_version(python_binary: &Path) -> Option<Version> {
    let output = std::process::Command::new(python_binary)
        .arg("-c")
        .arg("import deepchem; print(deepchem.__version__)")
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
        let info = DeepChemAdapter::new().info();
        assert_eq!(info.id, "deepchem");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "MIT");
        assert_eq!(info.display_name, "DeepChem");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = DeepChemAdapter::new().info();
        // DeepChem 2.7.x is the first PyTorch-2-clean release; upper
        // bound 3.0 leaves room for a future major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(2, 7, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(3, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = DeepChemAdapter::new().capabilities();
        // Capability variants land in a future task; ribbon
        // contributions are already enough for the registry to surface
        // the adapter.
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.deepchem.script"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = DeepChemAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
