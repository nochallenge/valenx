//! # valenx-adapter-openmm
//!
//! Adapter for [OpenMM](https://openmm.org/) — the dev-friendly,
//! Python-native molecular-dynamics engine: high-performance GPU
//! kernels exposed through a clean Python API. Distinct from the
//! GROMACS / LAMMPS adapters already in the workspace; OpenMM
//! shines when scripting custom integrators, free-energy protocols,
//! or quick-turn protein minimisation runs.
//!
//! **Phase 17 — subprocess wrapper for user-provided scripts.** The
//! adapter doesn't generate Python; the user supplies a `run.py`
//! (or whatever filename) referenced from `[bio.openmm].script` in
//! `case.toml`. `prepare()` stages the script into the workdir and
//! `run()` invokes `python <script>` via the shared subprocess
//! runner. The script is responsible for writing the named PDB +
//! DCD outputs.
//!
//! On `collect()` the named PDB is parsed via
//! [`valenx_bio::format::pdb::read`] and surfaced as a typed
//! `Native` artifact (we hold onto the raw file too so the user can
//! re-open it in PyMOL / VMD without round-tripping through Valenx).
//! Full DCD trajectory parsing lands in Task 11 alongside the
//! MDAnalysis adapter.

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

use crate::case_input::OpenMmInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(OpenMmAdapter::new())
}

pub struct OpenMmAdapter;

impl OpenMmAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for OpenMmAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "openmm";
/// Python interpreter candidates. `python3` first because on Linux
/// `python` may still be Python 2 on legacy distros; on Windows
/// `python` typically resolves to the Windows Store / 3.x install.
const PYTHON_BINARIES: &[&str] = &["python3", "python"];

impl Adapter for OpenMmAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "OpenMM",
            // OpenMM 8.0 is the first release that ships the
            // restructured `openmm` top-level package (the legacy
            // `simtk.openmm` namespace is gone). We expect to
            // revisit the upper bound when 9.x lands.
            version_range: VersionRange {
                min_inclusive: Version::new(8, 0, 0),
                max_exclusive: Version::new(9, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "MIT",
            docs_url: "http://docs.openmm.org/",
            homepage_url: "https://openmm.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(PYTHON_BINARIES) {
            Some(binary_path) => {
                // Try `import openmm; print(openmm.__version__)`
                // first — confirms OpenMM is importable from the
                // chosen interpreter (vs. just having Python on
                // PATH). Fall back to a "couldn't import" warning
                // so the probe still surfaces a useful state.
                let found_version = detect_openmm_version(&binary_path);
                let mut warnings = Vec::new();
                if found_version.is_none() {
                    warnings.push(
                        "probe found `python` on PATH but could not import \
                         `openmm` — install OpenMM with `conda install -c \
                         conda-forge openmm` (or `pip install openmm` on \
                         supported platforms) for runs to succeed"
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
                hint: "Python 3.9+ with OpenMM installed; `conda install -c \
                       conda-forge openmm` after ensuring python3 is on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = OpenMmInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve the script path against the case directory.
        // `confined_join` rejects absolute paths and `..` traversal so
        // the staged copy stays confined to the case directory.
        let source_script = confined_join(&case.path, &input.script)?;
        if !source_script.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.openmm].script `{}` not found (resolved {})",
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
                        "[bio.openmm].script path `{}` has no filename",
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
            // OpenMM MD runs span minutes (small minimisation) to
            // hours / days (production simulations on slow
            // hardware). Pick a generous 4-hour default; long
            // runs override via their own progress reporting.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting OpenMM", |line| {
            let mut hint = subprocess::Hint::default();
            // Convention: the user-supplied script can emit a
            // sentinel line `[valenx] openmm done` to signal
            // completion before exit; we lift that to a 95% progress
            // tick so the UI doesn't sit at "indeterminate" until
            // the process actually exits.
            if line.contains("[valenx] openmm done") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.starts_with("Step ") || line.contains("\"Step\"") {
                // OpenMM's StateDataReporter emits `Step,...` /
                // `Step ` lines as the simulation marches. We
                // don't know the total step count from outside the
                // script, so we can't bump a real progress
                // percentage — but the line is already being
                // forwarded through `ctx.log()` by the runner
                // (this handler runs *before* the log push), which
                // keeps the spinner alive. No-op.
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
        // Re-parse the case-input so we know where to look for the
        // PDB / DCD outputs. The case directory isn't on the
        // PreparedJob, but the staged script is the canonical
        // hashable input — same convention as the Biopython adapter.
        let script_path = first_script_in_workdir(&job.workdir);
        let case_hash_input = script_path
            .clone()
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "OpenMM",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk the workdir and pick up any PDB / DCD files. We
        // don't strictly enforce the case-input filenames at this
        // layer — if the script wrote `equilibrated.pdb` instead
        // of `minimised.pdb` we still want the user to see it. Any
        // PDB we find gets parsed; failed parses degrade to a
        // raw `Native` artifact with a parse-warning label.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-openmm", ?e, "workdir read failed");
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
            let (kind, label) = match ext.as_deref() {
                Some("pdb") => {
                    // Soft-validate via the canonical PDB reader
                    // — never fail collect() on a partial / weird
                    // PDB; surface the warning in the artifact
                    // label so the UI can show it. The PDB id is
                    // taken from the file stem so multiple PDBs
                    // (post-min, post-equil, …) get distinct
                    // labels.
                    let stem = path
                        .file_stem()
                        .and_then(|s| s.to_str())
                        .unwrap_or("openmm")
                        .to_string();
                    // Round-22 M2: cap the per-PDB read at
                    // MAX_PDB_FILE_BYTES (256 MiB) so a poisoned
                    // workdir with a multi-GB `.pdb` can't OOM
                    // `collect()` before the parser runs.
                    let label = match valenx_core::io_caps::read_capped_to_string(
                        &path,
                        valenx_core::io_caps::MAX_PDB_FILE_BYTES as usize,
                    ) {
                        Ok(text) => match valenx_bio::format::pdb::read(&stem, &text) {
                            Ok(structure) => format!(
                                "OpenMM PDB `{}` ({} atoms, {} residues)",
                                stem,
                                structure.atom_count(),
                                structure.residue_count()
                            ),
                            Err(e) => format!(
                                "OpenMM PDB `{}` (parse warning: {})",
                                stem,
                                e.to_string().lines().next().unwrap_or("invalid")
                            ),
                        },
                        Err(_) => format!("OpenMM PDB `{stem}`"),
                    };
                    (ArtifactKind::Native, label)
                }
                Some("dcd") => {
                    // Full DCD parsing lands in Task 11 alongside
                    // MDAnalysis; for now surface as a Native
                    // artifact so the user can open it in VMD /
                    // mdtraj / nglview without round-tripping
                    // through Valenx.
                    (ArtifactKind::Native, "OpenMM trajectory (DCD)".to_string())
                }
                Some("py") => (ArtifactKind::Other, "OpenMM script".to_string()),
                Some("csv") | Some("tsv") => (
                    ArtifactKind::Tabular,
                    "OpenMM StateDataReporter log".to_string(),
                ),
                Some("log") => (ArtifactKind::Log, "OpenMM log".to_string()),
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
            ribbon_contributions: vec!["bio.openmm.script"],
        }
    }
}

/// Lift the staged Python script out of a workdir for provenance
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

/// Run `python -c "import openmm; print(openmm.__version__)"` and
/// parse a `semver::Version` out of stdout. OpenMM version strings
/// look like `8.1.1` or `8.2.0.dev` (pre-releases trim cleanly via
/// `extract_semver`). Returns `None` on any failure (interpreter
/// unusable, OpenMM not importable, version string malformed) —
/// `probe()` falls back to an "OpenMM not importable" warning in
/// that case.
fn detect_openmm_version(python_binary: &Path) -> Option<Version> {
    let output = std::process::Command::new(python_binary)
        .arg("-c")
        .arg("import openmm; print(openmm.__version__)")
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

    /// Minimal valid PDB record covering one residue. Fixed-width
    /// columns matter — the reader is column-indexed. Hand-built
    /// from the wwPDB ATOM record spec.
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
        let info = OpenMmAdapter::new().info();
        assert_eq!(info.id, "openmm");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "MIT");
        assert_eq!(info.display_name, "OpenMM");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = OpenMmAdapter::new().info();
        // We support OpenMM >= 8.0 (top-level `openmm` namespace);
        // expect to revisit upper bound when 9.x lands.
        assert_eq!(info.version_range.min_inclusive, Version::new(8, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(9, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = OpenMmAdapter::new().capabilities();
        // Capability variants land in a future task; ribbon
        // contributions are already enough for the registry to
        // surface the adapter.
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.openmm.script"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = OpenMmAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    #[test]
    fn first_script_in_workdir_picks_lexicographic_first() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-openmm-script-{}",
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

    /// `collect()` exercises the full PDB + DCD classification path.
    /// PDB content here parses cleanly via valenx_bio so the label
    /// must include the "atoms / residues" summary text.
    #[test]
    fn collect_parses_pdb_and_lists_dcd() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-openmm-collect-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        fs::write(tmp.join("run.py"), b"# placeholder").unwrap();
        fs::write(tmp.join("minimised.pdb"), SAMPLE_PDB).unwrap();
        fs::write(tmp.join("trajectory.dcd"), b"\x00\x00\x00\x00fake DCD").unwrap();
        fs::write(tmp.join("state.csv"), b"step,energy\n0,-1.0\n").unwrap();

        let job = PreparedJob {
            workdir: tmp.clone(),
            native_command: vec![],
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: true,
        };
        let results = OpenMmAdapter::new().collect(&job).unwrap();

        // Find the parsed PDB artifact and inspect its label —
        // it must include atom + residue counts from the parser.
        let pdb_art = results
            .artifacts
            .iter()
            .find(|a| a.path.extension().is_some_and(|e| e == "pdb"))
            .expect("PDB artifact present");
        assert_eq!(pdb_art.kind, ArtifactKind::Native);
        assert!(
            pdb_art.label.contains("5 atoms"),
            "label was: {}",
            pdb_art.label
        );
        assert!(
            pdb_art.label.contains("1 residues"),
            "label was: {}",
            pdb_art.label
        );

        // DCD is listed as Native with the documented label.
        let dcd_art = results
            .artifacts
            .iter()
            .find(|a| a.path.extension().is_some_and(|e| e == "dcd"))
            .expect("DCD artifact present");
        assert_eq!(dcd_art.kind, ArtifactKind::Native);
        assert_eq!(dcd_art.label, "OpenMM trajectory (DCD)");

        // CSV (StateDataReporter) classifies as Tabular.
        let csv_art = results
            .artifacts
            .iter()
            .find(|a| a.path.extension().is_some_and(|e| e == "csv"))
            .expect("CSV artifact present");
        assert_eq!(csv_art.kind, ArtifactKind::Tabular);

        let _ = fs::remove_dir_all(&tmp);
    }

    /// A malformed PDB shouldn't crash collect — it should degrade
    /// to a "parse warning" label so the UI can still surface the
    /// raw file.
    #[test]
    fn collect_pdb_parse_failure_degrades_gracefully() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-openmm-bad-pdb-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        // ATOM lines must be >= 78 cols; this one is far too short.
        fs::write(tmp.join("broken.pdb"), b"ATOM      1  N   ALA A   1\n").unwrap();

        let job = PreparedJob {
            workdir: tmp.clone(),
            native_command: vec![],
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: true,
        };
        let results = OpenMmAdapter::new().collect(&job).unwrap();
        let pdb_art = results
            .artifacts
            .iter()
            .find(|a| a.path.extension().is_some_and(|e| e == "pdb"))
            .expect("artifact still surfaced");
        assert_eq!(pdb_art.kind, ArtifactKind::Native);
        assert!(
            pdb_art.label.contains("parse warning"),
            "label was: {}",
            pdb_art.label
        );
        let _ = fs::remove_dir_all(&tmp);
    }

    /// Empty-PDB edge case (no ATOM lines): the parser returns an
    /// empty Structure; collect should still surface the file with
    /// a 0-atom label rather than a parse warning.
    #[test]
    fn collect_empty_pdb_succeeds_with_zero_atoms() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-openmm-empty-pdb-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        fs::write(tmp.join("nothing.pdb"), b"REMARK only\nEND\n").unwrap();

        let job = PreparedJob {
            workdir: tmp.clone(),
            native_command: vec![],
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: true,
        };
        let results = OpenMmAdapter::new().collect(&job).unwrap();
        let pdb_art = results
            .artifacts
            .iter()
            .find(|a| a.path.extension().is_some_and(|e| e == "pdb"))
            .expect("artifact present");
        assert!(
            pdb_art.label.contains("0 atoms"),
            "label: {}",
            pdb_art.label
        );
        assert!(
            !pdb_art.label.contains("parse warning"),
            "label: {}",
            pdb_art.label
        );
        let _ = fs::remove_dir_all(&tmp);
    }
}
