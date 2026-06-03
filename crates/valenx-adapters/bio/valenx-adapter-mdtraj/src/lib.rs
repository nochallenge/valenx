//! # valenx-adapter-mdtraj
//!
//! Adapter for [MDTraj](https://mdtraj.org/) — the Pande / VanderSpoel
//! / Beauchamp lab's Python library for molecular-dynamics trajectory
//! analysis. Sister to MDAnalysis with wider format support (XTC /
//! TRR / DCD / NetCDF / HDF5 / LAMMPS / TNG / multi-model PDB / ...)
//! and tighter OpenMM integration. RMSD / RMSF / SASA / radius of
//! gyration / hydrogen bonds / dihedrals / contact maps all live in
//! a unified `Trajectory` API on top of NumPy arrays.
//!
//! **Phase 5.7 — subprocess wrapper for user-provided MDTraj
//! scripts.** MDTraj is a library (not a CLI), so the adapter itself
//! doesn't generate Python; the user authors an `analyze.py`
//! referenced from `[bio.mdtraj].script` in `case.toml` plus the
//! trajectory + topology to operate on and an output basename.
//! `prepare()` stages the script + trajectory + topology into the
//! workdir, drops a flat `valenx_params.json` next to it so the
//! script can read the parsed knobs without re-parsing case.toml,
//! and `run()` invokes `python <script>` via the shared subprocess
//! runner.
//!
//! ## `valenx_params.json`
//!
//! ```json
//! {
//!   "output_basename": "analysis",
//!   "trajectory":      "md.xtc",
//!   "topology":        "md.pdb"
//! }
//! ```
//!
//! User scripts read it with `json.load(open("valenx_params.json"))`
//! and resolve the bare filenames against the cwd (the workdir).
//!
//! On `collect()` we walk the workdir for `<basename>*.csv`
//! (analysis tables), `<basename>*.npz` (NumPy archives MDTraj
//! commonly emits), `<basename>*.h5` (MDTraj's native HDF5
//! trajectory format), and `<basename>*.png` (plots), plus any
//! `*.log` files Python emits.

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

use crate::case_input::MdtrajInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(MdtrajAdapter::new())
}

pub struct MdtrajAdapter;

impl MdtrajAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for MdtrajAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "mdtraj";
/// Python interpreter candidates. `python3` first because on Linux
/// `python` may still be Python 2 on legacy distros; on Windows
/// `python` typically resolves to the Windows Store / 3.x install.
const PYTHON_BINARIES: &[&str] = &["python3", "python"];

impl Adapter for MdtrajAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "MDTraj",
            // MDTraj 1.9+ (2018+) is the modern stable line that
            // shipped the unified `Trajectory` API + the wider
            // format-reader matrix downstream scripts target. 2.0 is
            // reserved for an upcoming major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(1, 9, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "LGPL-2.1",
            docs_url: "https://mdtraj.org/",
            homepage_url: "https://mdtraj.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(PYTHON_BINARIES) {
            Some(binary_path) => {
                // Try `import mdtraj` — confirms the `mdtraj` package
                // is importable from the chosen interpreter (vs. just
                // having Python on PATH). Failure is non-fatal; we
                // surface a warning so the user can still configure a
                // different interpreter via `python` and a missing-
                // package case can still validate.
                let import_ok = mdtraj_importable(&binary_path);
                let mut warnings = Vec::new();
                if !import_ok {
                    warnings.push(
                        "probe found `python` on PATH but could not import \
                         `mdtraj` — install MDTraj with `pip install \
                         mdtraj` (or `conda install -c conda-forge mdtraj`) \
                         for runs to succeed"
                            .into(),
                    );
                }
                Ok(ProbeReport {
                    ok: true,
                    // We intentionally don't surface a found_version —
                    // probing the package version requires running
                    // Python, which is enough work that a missing
                    // package is the more useful signal.
                    found_version: None,
                    binary_path: Some(binary_path),
                    warnings,
                    required_env: Vec::new(),
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "Python 3.9+ with MDTraj installed; \
                       `pip install mdtraj` (or `conda install -c \
                       conda-forge mdtraj`) after ensuring python3 \
                       is on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = MdtrajInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.mdtraj].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Stage the user-supplied Python script into the workdir.
        // `confined_join` rejects absolute paths and `..` traversal so
        // the staged copy stays confined to the case directory.
        let source_script = confined_join(&case.path, &input.script)?;
        if !source_script.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.mdtraj].script `{}` not found (resolved {})",
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
                        "[bio.mdtraj].script path `{}` has no filename",
                        input.script.display()
                    ),
                })?;
        let dest_script = workdir.join(script_filename);
        if source_script != dest_script {
            fs::copy(&source_script, &dest_script)?;
        }

        // Stage the trajectory so the script can resolve it via a bare
        // relative filename inside the workdir. MDTraj sniffs the
        // format from the extension, so we preserve the user-supplied
        // basename verbatim.
        let source_traj = confined_join(&case.path, &input.trajectory)?;
        if !source_traj.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.mdtraj].trajectory `{}` not found (resolved {})",
                    input.trajectory.display(),
                    source_traj.display()
                ),
            });
        }
        let traj_filename =
            input
                .trajectory
                .file_name()
                .ok_or_else(|| AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.mdtraj].trajectory path `{}` has no filename",
                        input.trajectory.display()
                    ),
                })?;
        let dest_traj = workdir.join(traj_filename);
        if source_traj != dest_traj {
            fs::copy(&source_traj, &dest_traj)?;
        }

        // Stage the topology similarly. MDTraj's `Trajectory.load`
        // takes the topology as a separate argument because most
        // trajectory formats (XTC / TRR / DCD / NetCDF) carry only
        // coordinates — bonds / residues / chains live in the
        // topology file (`.pdb` / `.psf` / `.prmtop` / `.gro` / ...).
        let source_top = confined_join(&case.path, &input.topology)?;
        if !source_top.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.mdtraj].topology `{}` not found (resolved {})",
                    input.topology.display(),
                    source_top.display()
                ),
            });
        }
        let top_filename = input
            .topology
            .file_name()
            .ok_or_else(|| AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.mdtraj].topology path `{}` has no filename",
                    input.topology.display()
                ),
            })?;
        let dest_top = workdir.join(top_filename);
        if source_top != dest_top {
            fs::copy(&source_top, &dest_top)?;
        }

        // Drop a flat `valenx_params.json` into the workdir so the
        // user's Python script can read the parsed `[bio.mdtraj]`
        // knobs without having to reparse case.toml itself. Built by
        // hand to avoid pulling in serde_json for a 3-key flat object.
        // `trajectory` and `topology` are bare basenames pointing into
        // the workdir (not the original full paths).
        let mut params = String::new();
        params.push_str("{\n");
        params.push_str("  \"output_basename\": ");
        params.push_str(&json_string(&input.output_basename));
        params.push_str(",\n  \"trajectory\": ");
        params.push_str(&json_string(&traj_filename.to_string_lossy()));
        params.push_str(",\n  \"topology\": ");
        params.push_str(&json_string(&top_filename.to_string_lossy()));
        params.push_str("\n}\n");
        valenx_core::io_caps::atomic_write_str(&workdir.join("valenx_params.json"), &params)?;

        // Resolve the Python binary. Bare `python` / `python3` walks
        // PATH; absolute / relative paths the user pinned are honored
        // verbatim if they exist, with a final PATH fallback.
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
            // MDTraj passes vary widely — RMSD on a small protein
            // finishes in seconds; SASA / hbond surveys against
            // multi-microsecond trajectories run for hours. 2 hours
            // is a generous default that covers the long tail.
            estimated_runtime: Some(Duration::from_secs(2 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting MDTraj", |line| {
            let mut hint = subprocess::Hint::default();
            // Convention: the user-supplied script can emit a sentinel
            // line `[valenx] mdtraj done` to signal completion before
            // exit; lift to a 95% progress tick.
            if line.contains("[valenx] mdtraj done") {
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
        // Provenance: hash the staged Python script as the canonical
        // input descriptor. Falls back to case.toml when no script is
        // staged yet (partial / failed runs).
        let script_path = first_script_in_workdir(&job.workdir);
        let case_hash_input = script_path
            .clone()
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "MDTraj",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Restrict typed outputs (.csv / .npz / .h5 / .png) to those
        // whose stem starts with the configured `output_basename`.
        // .log files are accepted regardless — Python's logging output
        // isn't typically prefixed.
        let basename = read_output_basename(&job.workdir);

        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-mdtraj", ?e, "workdir read failed");
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
            let stem_matches_basename = match basename.as_deref() {
                Some(b) => stem.starts_with(b),
                None => true,
            };
            match ext.as_deref() {
                Some("csv") => {
                    if !stem_matches_basename {
                        continue;
                    }
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Tabular,
                        checksum: None,
                        label: "MDTraj analysis table".to_string(),
                    });
                }
                Some("npz") => {
                    if !stem_matches_basename {
                        continue;
                    }
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Native,
                        checksum: None,
                        label: "MDTraj numpy archive".to_string(),
                    });
                }
                Some("h5") => {
                    if !stem_matches_basename {
                        continue;
                    }
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Native,
                        checksum: None,
                        label: "MDTraj HDF5 output".to_string(),
                    });
                }
                Some("png") => {
                    if !stem_matches_basename {
                        continue;
                    }
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Native,
                        checksum: None,
                        label: "MDTraj plot".to_string(),
                    });
                }
                Some("log") => {
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Log,
                        checksum: None,
                        label: "MDTraj log".to_string(),
                    });
                }
                _ => continue,
            }
        }
        artefacts.sort_by(|a, b| a.path.cmp(&b.path));
        results.artifacts = artefacts;
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.mdtraj.analyze"],
        }
    }
}

/// Escape a string for embedding inside a JSON string literal.
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

/// Lift the staged Python script out of the workdir for provenance
/// hashing. Returns the lexicographically-first `.py` file
/// (case-insensitive) at the top level, or `None` if none exists yet.
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
/// `valenx_params.json` for collect()-time output filtering.
fn read_output_basename(workdir: &Path) -> Option<String> {
    let text = valenx_core::io_caps::read_capped_to_string(
        &workdir.join("valenx_params.json"),
        valenx_core::io_caps::MAX_ADAPTER_PARAMS_BYTES as usize,
    )
    .ok()?;
    extract_json_string(&text, "output_basename")
}

/// Pull a flat string field out of our own hand-emitted
/// `valenx_params.json`. We wrote the file ourselves so we know its
/// shape; a full JSON parser would be overkill.
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

/// Run `python -c "import mdtraj"` and return whether it succeeded.
/// Returns `false` on any failure (interpreter unusable, mdtraj not
/// importable); `probe()` lifts that to a "mdtraj not importable"
/// warning rather than a hard error so a missing-package case can
/// still validate.
fn mdtraj_importable(python_binary: &Path) -> bool {
    let output = std::process::Command::new(python_binary)
        .arg("-c")
        .arg("import mdtraj")
        .stdin(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .output();
    match output {
        Ok(o) => o.status.success(),
        Err(_) => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = MdtrajAdapter::new().info();
        assert_eq!(info.id, "mdtraj");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "LGPL-2.1");
        assert_eq!(info.display_name, "MDTraj");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = MdtrajAdapter::new().info();
        // MDTraj 1.9+ shipped the modern `Trajectory` API + the wider
        // format-reader matrix downstream scripts target. 2.0 reserves
        // room for the upcoming major.
        assert_eq!(info.version_range.min_inclusive, Version::new(1, 9, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(2, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = MdtrajAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.mdtraj.analyze"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = MdtrajAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
