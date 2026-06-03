//! # valenx-adapter-mdanalysis
//!
//! Adapter for [MDAnalysis](https://www.mdanalysis.org/) — the Python
//! library for analysing molecular-dynamics trajectories: RMSD,
//! radius-of-gyration, hydrogen-bond surveys, contact maps, free-energy
//! reweighting helpers, and a unified reader API across every common
//! trajectory format (DCD / XTC / TRR / NetCDF / …).
//!
//! **Phase 17 — subprocess wrapper for user-provided scripts.** The
//! adapter doesn't generate Python; the user supplies an
//! `analyse_traj.py` (or whatever filename) referenced from
//! `[bio.mdanalysis].script` in `case.toml`. `prepare()` stages the
//! script into the workdir and `run()` invokes `python <script>` via
//! the shared subprocess runner.
//!
//! On `collect()` we walk the workdir for the customary outputs:
//! `.dcd` trajectories (parsed via [`valenx_bio::format::dcd::read`]
//! into a typed [`valenx_bio::Trajectory`] for the artifact label) and
//! `.csv` analysis tables. PDB / XYZ structure outputs aren't part
//! of MDAnalysis's typical write surface, so we skip them; CSV /
//! TSV pickup catches the bulk of analysis-script outputs.

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

use crate::case_input::MdAnalysisInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(MdAnalysisAdapter::new())
}

pub struct MdAnalysisAdapter;

impl MdAnalysisAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for MdAnalysisAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "mdanalysis";
/// Python interpreter candidates. `python3` first because on Linux
/// `python` may still be Python 2 on legacy distros; on Windows
/// `python` typically resolves to the Windows Store / 3.x install.
const PYTHON_BINARIES: &[&str] = &["python3", "python"];

impl Adapter for MdAnalysisAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "MDAnalysis",
            // MDAnalysis 2.7 (2024) is the first release with stable
            // `analysis.AtomicDistances` + the modernised reader API
            // we rely on. Upper bound 3.0 reserves room for an
            // upcoming major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(2, 7, 0),
                max_exclusive: Version::new(3, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-2.0-or-later",
            docs_url: "https://docs.mdanalysis.org/",
            homepage_url: "https://www.mdanalysis.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(PYTHON_BINARIES) {
            Some(binary_path) => {
                // The package import is `MDAnalysis` (camel-case) —
                // `mdanalysis` lower-cased import fails. Mirroring
                // the case in adapter id ("mdanalysis") vs probe
                // ("MDAnalysis") deliberately.
                let found_version = detect_mdanalysis_version(&binary_path);
                let mut warnings = Vec::new();
                if found_version.is_none() {
                    warnings.push(
                        "probe found `python` on PATH but could not import \
                         `MDAnalysis` — install with `pip install MDAnalysis` \
                         (or `conda install -c conda-forge mdanalysis`) for \
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
                hint: "Python 3.9+ with MDAnalysis installed; \
                       `pip install MDAnalysis` (or `conda install -c \
                       conda-forge mdanalysis`) after ensuring python3 \
                       is on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = MdAnalysisInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve the script path against the case directory. Same
        // convention as every other Phase 17 bio adapter —
        // `script = "analyse.py"` next to `case.toml`. `confined_join`
        // rejects absolute paths and `..` traversal so a malicious case
        // bundle can't smuggle arbitrary host files into the workdir.
        let source_script = confined_join(&case.path, &input.script)?;
        if !source_script.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.mdanalysis].script `{}` not found (resolved {})",
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
                        "[bio.mdanalysis].script path `{}` has no filename",
                        input.script.display()
                    ),
                })?;
        let dest_script = workdir.join(script_filename);
        if source_script != dest_script {
            fs::copy(&source_script, &dest_script)?;
        }

        // Resolve the Python binary. Same pinning logic as the other
        // Python-script adapters (Biopython / OpenMM / RDKit).
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
            // MDAnalysis trajectory passes vary widely — a quick RMSD
            // sweep finishes in seconds; exhaustive contact-map
            // surveys against multi-microsecond trajectories run for
            // hours. Pick 30 minutes as a reasonable default.
            estimated_runtime: Some(Duration::from_secs(30 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting MDAnalysis", |line| {
            let mut hint = subprocess::Hint::default();
            // Convention: the user-supplied script can emit a
            // sentinel `[valenx] mdanalysis done` to signal
            // completion before exit; lift to a 95% progress tick.
            if line.contains("[valenx] mdanalysis done") {
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
        // Provenance: hash the staged Python script (the canonical
        // "this case is configured this way" input). Falls back to
        // case.toml when the script isn't present yet.
        let script_path = first_script_in_workdir(&job.workdir);
        let case_hash_input = script_path
            .clone()
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "MDAnalysis",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk the workdir top-level. MDAnalysis scripts that need
        // nested output directories will surface their key artefacts
        // at the top level by convention.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-mdanalysis", ?e, "workdir read failed");
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
                // .dcd — the bread-and-butter MDAnalysis output.
                // Try the Valenx DCD reader; on success surface the
                // frame + atom counts in the label so the UI can
                // show the trajectory's shape without re-opening it
                // in Python. Failed parses degrade to a generic
                // label so the user can still find the file.
                Some("dcd") => {
                    // Round-22 M2: cap the per-DCD read at
                    // MAX_DCD_FRAME_FILE_BYTES (4 GiB) so a poisoned
                    // workdir with a truly pathological binary
                    // trajectory can't OOM `collect()` before the
                    // parser runs. 4 GiB accommodates honest long
                    // production runs (large nframes × natoms).
                    let label = match valenx_core::io_caps::read_capped_to_bytes(
                        &path,
                        valenx_core::io_caps::MAX_DCD_FRAME_FILE_BYTES,
                    ) {
                        Ok(bytes) => {
                            let stem = path
                                .file_stem()
                                .and_then(|s| s.to_str())
                                .unwrap_or("trajectory");
                            match valenx_bio::format::dcd::read(&bytes, stem) {
                                Ok(traj) => format!(
                                    "MDAnalysis trajectory · {} frames · {} atoms",
                                    traj.frame_count(),
                                    traj.atom_count().unwrap_or(0)
                                ),
                                Err(_) => "MDAnalysis trajectory (DCD)".to_string(),
                            }
                        }
                        Err(_) => "MDAnalysis trajectory (DCD)".to_string(),
                    };
                    (ArtifactKind::Native, label)
                }
                Some("csv") | Some("tsv") => (
                    ArtifactKind::Tabular,
                    "MDAnalysis analysis table".to_string(),
                ),
                Some("py") => (ArtifactKind::Other, "MDAnalysis script".to_string()),
                Some("log") => (ArtifactKind::Log, "MDAnalysis log".to_string()),
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
        // task; ribbon contributions are already enough for the
        // registry to surface the adapter.
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.mdanalysis.script"],
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

/// Run `python -c "import MDAnalysis as mda; print(mda.__version__)"`
/// and parse a `semver::Version` out of stdout. The package import
/// case matters — `MDAnalysis` capitalised. Returns `None` on any
/// failure (interpreter unusable, MDAnalysis not importable, version
/// string malformed); `probe()` falls back to a "MDAnalysis not
/// importable" warning in that case.
fn detect_mdanalysis_version(python_binary: &Path) -> Option<Version> {
    let output = std::process::Command::new(python_binary)
        .arg("-c")
        .arg("import MDAnalysis as mda; print(mda.__version__)")
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

    /// Synthesise a minimal valid 1-frame, 7-atom DCD so collect()
    /// can exercise the real-reader success path with a label that
    /// surfaces non-zero frame + atom counts. Same wrapper layout
    /// used by the dcd_round_trip integration tests in valenx-bio.
    fn synth_one_frame_dcd() -> Vec<u8> {
        fn wrap(body: &[u8]) -> Vec<u8> {
            let mut out = Vec::with_capacity(body.len() + 8);
            out.extend_from_slice(&(body.len() as u32).to_le_bytes());
            out.extend_from_slice(body);
            out.extend_from_slice(&(body.len() as u32).to_le_bytes());
            out
        }
        let mut header = Vec::with_capacity(84);
        header.extend_from_slice(b"CORD");
        header.extend_from_slice(&1i32.to_le_bytes()); // nframes
        for _ in 0..8 {
            header.extend_from_slice(&0i32.to_le_bytes());
        }
        header.extend_from_slice(&0.001f32.to_le_bytes());
        header.extend_from_slice(&0i32.to_le_bytes()); // cell flag
        for _ in 0..8 {
            header.extend_from_slice(&0i32.to_le_bytes());
        }
        header.extend_from_slice(&24i32.to_le_bytes());
        debug_assert_eq!(header.len(), 84);
        let mut titles = Vec::new();
        titles.extend_from_slice(&1i32.to_le_bytes());
        titles.extend_from_slice(&[b' '; 80]);
        // 7-atom coordinate streams (X, Y, Z) packed as f32.
        let zeros = [0f32; 7];
        let mut xs = Vec::new();
        let mut ys = Vec::new();
        let mut zs = Vec::new();
        for v in zeros.iter() {
            xs.extend_from_slice(&v.to_le_bytes());
            ys.extend_from_slice(&v.to_le_bytes());
            zs.extend_from_slice(&v.to_le_bytes());
        }
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&wrap(&header));
        bytes.extend_from_slice(&wrap(&titles));
        bytes.extend_from_slice(&wrap(&7i32.to_le_bytes()));
        bytes.extend_from_slice(&wrap(&xs));
        bytes.extend_from_slice(&wrap(&ys));
        bytes.extend_from_slice(&wrap(&zs));
        bytes
    }

    #[test]
    fn info_is_bio_domain() {
        let info = MdAnalysisAdapter::new().info();
        assert_eq!(info.id, "mdanalysis");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "GPL-2.0-or-later");
        assert_eq!(info.display_name, "MDAnalysis");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = MdAnalysisAdapter::new().info();
        assert_eq!(info.version_range.min_inclusive, Version::new(2, 7, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(3, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = MdAnalysisAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.mdanalysis.script"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = MdAnalysisAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    #[test]
    fn first_script_in_workdir_picks_lexicographic_first() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-mdanalysis-script-{}",
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
            "valenx-mdanalysis-empty-{}",
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

    /// The real DCD reader landed alongside this adapter, so collect()
    /// must surface frame + atom counts in the artifact label when a
    /// valid DCD is present in the workdir.
    #[test]
    fn collect_parses_dcd_metadata_into_label() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-mdanalysis-dcd-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        fs::write(tmp.join("analyse.py"), b"# placeholder").unwrap();
        fs::write(tmp.join("traj.dcd"), synth_one_frame_dcd()).unwrap();
        fs::write(tmp.join("rmsd.csv"), b"frame,rmsd\n0,0.0\n").unwrap();

        let job = PreparedJob {
            workdir: tmp.clone(),
            native_command: vec![],
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: true,
        };
        let results = MdAnalysisAdapter::new().collect(&job).unwrap();

        let dcd_art = results
            .artifacts
            .iter()
            .find(|a| a.path.extension().is_some_and(|e| e == "dcd"))
            .expect("DCD artifact present");
        assert_eq!(dcd_art.kind, ArtifactKind::Native);
        // Real reader → label includes "frames" + "atoms" rather than
        // the generic fallback.
        assert!(
            dcd_art.label.contains("frames") && dcd_art.label.contains("atoms"),
            "label was: {}",
            dcd_art.label
        );
        // The synthetic DCD declares 1 frame * 7 atoms — verify
        // both numbers reach the label.
        assert!(
            dcd_art.label.contains("1 frames"),
            "label: {}",
            dcd_art.label
        );
        assert!(
            dcd_art.label.contains("7 atoms"),
            "label: {}",
            dcd_art.label
        );

        // CSV classifies as Tabular.
        let csv_art = results
            .artifacts
            .iter()
            .find(|a| a.path.extension().is_some_and(|e| e == "csv"))
            .expect("CSV artifact present");
        assert_eq!(csv_art.kind, ArtifactKind::Tabular);
        assert_eq!(csv_art.label, "MDAnalysis analysis table");

        let _ = fs::remove_dir_all(&tmp);
    }

    /// A malformed DCD shouldn't crash collect — it should degrade
    /// to the generic "MDAnalysis trajectory (DCD)" label so the UI
    /// can still surface the raw file.
    #[test]
    fn collect_dcd_parse_failure_degrades_gracefully() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-mdanalysis-bad-dcd-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        // Tiny garbage payload — no chance of parsing as DCD.
        fs::write(tmp.join("broken.dcd"), b"\x00\x00\x00\x00fake DCD").unwrap();

        let job = PreparedJob {
            workdir: tmp.clone(),
            native_command: vec![],
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: true,
        };
        let results = MdAnalysisAdapter::new().collect(&job).unwrap();
        let dcd_art = results
            .artifacts
            .iter()
            .find(|a| a.path.extension().is_some_and(|e| e == "dcd"))
            .expect("artifact still surfaced");
        assert_eq!(dcd_art.kind, ArtifactKind::Native);
        assert_eq!(dcd_art.label, "MDAnalysis trajectory (DCD)");
        let _ = fs::remove_dir_all(&tmp);
    }

    /// .bin / .dat files outside the recognised extension list must
    /// not surface — guards the deny-by-default classification path.
    #[test]
    fn collect_skips_unknown_extensions() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-mdanalysis-skip-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        fs::write(tmp.join("foo.bin"), b"...").unwrap();
        fs::write(tmp.join("bar.dat"), b"...").unwrap();
        // One recognised file so collect doesn't bail before the
        // walk; we only assert the unrecognised ones are skipped.
        fs::write(tmp.join("ok.csv"), b"a,b\n").unwrap();

        let job = PreparedJob {
            workdir: tmp.clone(),
            native_command: vec![],
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: true,
        };
        let results = MdAnalysisAdapter::new().collect(&job).unwrap();
        assert!(!results
            .artifacts
            .iter()
            .any(|a| a.path.extension().is_some_and(|e| e == "bin")));
        assert!(!results
            .artifacts
            .iter()
            .any(|a| a.path.extension().is_some_and(|e| e == "dat")));
        let _ = fs::remove_dir_all(&tmp);
    }
}
