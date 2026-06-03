//! # valenx-adapter-rdkit
//!
//! Adapter for [RDKit](https://www.rdkit.org/) — the open-source
//! cheminformatics toolkit: SMILES / SDF / MOL parsing, fingerprints,
//! similarity search, 2D depiction, conformer generation.
//!
//! **Phase 17 — subprocess wrapper for user-provided scripts.** The
//! adapter doesn't generate Python; the user supplies a `screen.py`
//! (or whatever filename) referenced from `[bio.rdkit].script` in
//! `case.toml`. `prepare()` stages the script into the workdir and
//! `run()` invokes `python <script>` via the shared subprocess
//! runner. Inline SMILES lists are surfaced for scripts that prefer
//! configuration over input files.
//!
//! Future work: bundled-recipe variants (read SDF -> compute Morgan
//! fingerprints, ligand-docking pipelines, etc.) that generate the
//! Python the way `valenx-adapter-cantera` does, plus structured
//! Molecule parsing once the `valenx-bio` canonical Molecule type
//! lands.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::fs;
use std::path::Path;
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

use crate::case_input::RdkitInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(RdkitAdapter::new())
}

pub struct RdkitAdapter;

impl RdkitAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for RdkitAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "rdkit";
/// Python interpreter candidates. `python3` first because on Linux
/// `python` may still be Python 2 on legacy distros; on Windows
/// `python` typically resolves to the Windows Store / 3.x install.
const PYTHON_BINARIES: &[&str] = &["python3", "python"];

impl Adapter for RdkitAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "RDKit",
            // RDKit ships on a calendar-versioning cadence
            // (`YYYY.MM.patch`). 2023.9 is the first release that
            // ships wheels for Python 3.12 and the modern
            // `rdkit.Chem` API surface scripts target. Upper bound
            // 2026.0 leaves room for the 2025.x line and bumps when
            // 2026.x lands.
            version_range: VersionRange {
                min_inclusive: Version::new(2023, 9, 0),
                max_exclusive: Version::new(2026, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "BSD-3-Clause",
            docs_url: "https://www.rdkit.org/docs/",
            homepage_url: "https://www.rdkit.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(PYTHON_BINARIES) {
            Some(binary_path) => {
                // Try `import rdkit; print(rdkit.__version__)` first
                // — that's the only string that confirms RDKit
                // itself is installed (vs. just having Python on
                // PATH). Fall back to `python --version` if the
                // import fails so we can still surface a useful
                // probe state.
                let found_version = detect_rdkit_version(&binary_path);
                let mut warnings = Vec::new();
                if found_version.is_none() {
                    warnings.push(
                        "probe found `python` on PATH but could not import \
                         `rdkit` — install RDKit with `pip install rdkit` \
                         (or `conda install -c conda-forge rdkit`) for runs \
                         to succeed"
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
                hint: "Python 3.9+ with RDKit installed; `pip install rdkit` \
                       (or `conda install -c conda-forge rdkit`) after \
                       ensuring python3 is on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = RdkitInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve the script path against the case directory. The
        // user authors `script = "screen.py"` and expects it to live
        // alongside `case.toml`. `confined_join` rejects absolute paths
        // and `..` traversal so a malicious case bundle can't smuggle
        // arbitrary host files into the workdir.
        let source_script = confined_join(&case.path, &input.script)?;
        if !source_script.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.rdkit].script `{}` not found (resolved {})",
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
                        "[bio.rdkit].script path `{}` has no filename",
                        input.script.display()
                    ),
                })?;
        let dest_script = workdir.join(script_filename);
        if source_script != dest_script {
            fs::copy(&source_script, &dest_script)?;
        }

        // Resolve the Python binary. If the user pinned a specific
        // interpreter via `python = "..."`, honor it; otherwise
        // walk PATH.
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
            // RDKit cheminformatics workloads vary widely — a SMILES
            // sanity check finishes in seconds, conformer generation
            // for thousands of compounds runs for minutes-to-hours.
            // Pick 30 minutes as a reasonable default; long-running
            // cases override via their own progress reporting.
            estimated_runtime: Some(Duration::from_secs(30 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting RDKit", |line| {
            let mut hint = subprocess::Hint::default();
            // Convention: the user-supplied script can emit a
            // sentinel line `[valenx] rdkit done` to signal
            // completion before exit; we lift that to a 95% progress
            // tick so the UI doesn't sit at "indeterminate" until
            // the process actually exits.
            if line.contains("[valenx] rdkit done") {
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
        // Provenance: hash the staged script (the canonical
        // "this case is configured this way" input). We don't
        // know the user's mesh / lock files so leave those empty.
        let script_path = first_script_in_workdir(&job.workdir);
        let case_hash_input = script_path
            .clone()
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "RDKit",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);

        // Walk the workdir top level and classify any output files
        // the user's script left behind. We deliberately don't
        // recurse — RDKit scripts that need nested output
        // directories will surface their key artefacts at the top
        // level by convention.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-rdkit", ?e, "workdir read failed");
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
                // .sdf — multi-molecule format with embedded
                // 2D/3D coordinates; the bread-and-butter RDKit
                // output. Structured Molecule parsing lands in a
                // follow-up; for now classify as Native.
                Some("sdf") => (ArtifactKind::Native, "SDF molecule set".to_string()),
                // .mol — single-molecule MDL MOL file.
                Some("mol") => (ArtifactKind::Native, "MOL molecule".to_string()),
                Some("csv") | Some("tsv") => (ArtifactKind::Tabular, "Tabular output".to_string()),
                Some("json") => (ArtifactKind::Other, "JSON output".to_string()),
                Some("py") => (ArtifactKind::Other, "RDKit script".to_string()),
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
        // The bio-specific Capability variants land in a follow-up
        // task; for now we publish an empty capability vector and
        // a single ribbon contribution so the registry can wire the
        // adapter in without crashing the UI's capability-index
        // builder.
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.rdkit.script"],
        }
    }
}

/// Lift the staged Python script out of a workdir for provenance
/// hashing. Returns the lexicographically-first `.py` file at the
/// top level, or `None` if none exists yet.
fn first_script_in_workdir(workdir: &Path) -> Option<std::path::PathBuf> {
    let entries = fs::read_dir(workdir).ok()?;
    let mut hits: Vec<std::path::PathBuf> = entries
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

/// Run `python -c "import rdkit; print(rdkit.__version__)"` and
/// parse a `semver::Version` out of stdout. RDKit version strings
/// look like `2024.03.1`; we treat the first three dot-separated
/// components as the semver triple (treating the year as a major
/// version is fine for ordering since it monotonically increases).
/// Returns `None` on any failure (interpreter unusable, RDKit not
/// importable, version string malformed) — `probe()` falls back to
/// a "RDKit not importable" warning in that case.
fn detect_rdkit_version(python_binary: &Path) -> Option<Version> {
    let output = std::process::Command::new(python_binary)
        .arg("-c")
        .arg("import rdkit; print(rdkit.__version__)")
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
        let info = RdkitAdapter::new().info();
        assert_eq!(info.id, "rdkit");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "BSD-3-Clause");
        assert_eq!(info.display_name, "RDKit");
    }

    #[test]
    fn info_version_range_matches_calendar_versioning() {
        let info = RdkitAdapter::new().info();
        // RDKit uses YYYY.MM.patch — we treat the year as the
        // major component for semver ordering. Min 2023.9 ships
        // wheels for Python 3.12; max-exclusive 2026.0 bumps when
        // 2026.x lands.
        assert_eq!(info.version_range.min_inclusive, Version::new(2023, 9, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(2026, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = RdkitAdapter::new().capabilities();
        // Capability variants land in a future task; ribbon
        // contributions are already enough for the registry to
        // surface the adapter.
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.rdkit.script"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = RdkitAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    #[test]
    fn first_script_in_workdir_picks_lexicographic_first() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-rdkit-script-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        fs::write(tmp.join("z_late.py"), b"# placeholder").unwrap();
        fs::write(tmp.join("a_first.py"), b"# placeholder").unwrap();
        fs::write(tmp.join("not_python.txt"), b"placeholder").unwrap();
        let f = first_script_in_workdir(&tmp).expect("found");
        assert!(f.ends_with("a_first.py"));
        let _ = fs::remove_dir_all(&tmp);
    }

    #[test]
    fn first_script_in_workdir_returns_none_when_empty() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-rdkit-empty-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        fs::write(tmp.join("notes.md"), b"placeholder").unwrap();
        assert!(first_script_in_workdir(&tmp).is_none());
        let _ = fs::remove_dir_all(&tmp);
    }

    /// `collect()` runs against a workdir directly — exercise the
    /// SDF / MOL / CSV classification paths so a regression in the
    /// extension-dispatch table doesn't slip past CI.
    #[test]
    fn collect_classifies_rdkit_outputs() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-rdkit-collect-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        fs::write(tmp.join("screen.py"), b"# placeholder").unwrap();
        fs::write(tmp.join("hits.sdf"), b"# placeholder SDF").unwrap();
        fs::write(tmp.join("benzene.mol"), b"# placeholder MOL").unwrap();
        fs::write(tmp.join("scores.csv"), b"smiles,score\nCCO,0.5\n").unwrap();
        fs::write(tmp.join("ignore.bin"), b"...").unwrap();

        let job = PreparedJob {
            workdir: tmp.clone(),
            native_command: vec![],
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: true,
        };
        let results = RdkitAdapter::new().collect(&job).unwrap();
        let labels: Vec<&str> = results.artifacts.iter().map(|a| a.label.as_str()).collect();
        assert!(labels.contains(&"SDF molecule set"));
        assert!(labels.contains(&"MOL molecule"));
        assert!(labels.contains(&"Tabular output"));
        assert!(labels.contains(&"RDKit script"));
        // .bin must not surface — guards the deny-by-default path.
        assert!(!results
            .artifacts
            .iter()
            .any(|a| a.path.extension().is_some_and(|e| e == "bin")));
        let _ = fs::remove_dir_all(&tmp);
    }

    /// SDF and MOL classify as `Native` (binary-ish solver-format
    /// output); CSV classifies as `Tabular`. Pin the contract.
    #[test]
    fn collect_artifact_kinds_pin_native_vs_tabular() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-rdkit-kinds-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        fs::write(tmp.join("a.sdf"), b"").unwrap();
        fs::write(tmp.join("b.mol"), b"").unwrap();
        fs::write(tmp.join("c.csv"), b"").unwrap();

        let job = PreparedJob {
            workdir: tmp.clone(),
            native_command: vec![],
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: true,
        };
        let results = RdkitAdapter::new().collect(&job).unwrap();
        for art in &results.artifacts {
            let ext = art.path.extension().and_then(|e| e.to_str()).unwrap();
            match ext {
                "sdf" | "mol" => assert_eq!(art.kind, ArtifactKind::Native),
                "csv" => assert_eq!(art.kind, ArtifactKind::Tabular),
                _ => panic!("unexpected extension {ext}"),
            }
        }
        let _ = fs::remove_dir_all(&tmp);
    }
}
