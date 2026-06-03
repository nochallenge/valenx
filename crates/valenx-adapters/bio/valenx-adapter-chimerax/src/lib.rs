//! # valenx-adapter-chimerax
//!
//! Adapter for [ChimeraX](https://www.cgl.ucsf.edu/chimerax/) — UCSF's
//! next-generation 3D molecular visualisation tool, the successor to
//! Chimera. ChimeraX renders proteins, nucleic acids, ligands, and
//! complexes from PDB / CIF / many other structural formats; in
//! headless mode it executes `.cxc` command scripts to produce
//! publication-quality images, session files, and re-saved
//! structures.
//!
//! **Phase 17 — subprocess wrapper for user-provided scripts.** The
//! adapter doesn't generate ChimeraX commands; the user supplies a
//! `render.cxc` (or whatever filename) referenced from
//! `[bio.chimerax].script` in `case.toml`. `prepare()` stages the
//! script into the workdir and `run()` invokes
//! `chimerax --nogui --script <script>` via the shared subprocess
//! runner. `--nogui` keeps the run headless by default; flip the
//! `nogui` case-input flag to false for the (rare) interactive path.
//!
//! On `collect()` we walk the workdir for ChimeraX's customary output
//! mix: `.png` rendered images, `.cxs` session files, and any
//! `.pdb` / `.cif` structures the script wrote out. Images surface as
//! `Image` artifacts; sessions and structure outputs surface as
//! `Native` artifacts so the user can re-open them in ChimeraX (or
//! PyMOL / VMD for the structure files).

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Duration;

use semver::Version;

use valenx_core::{
    adapter_helpers::{confined_join, detect_tool_version_semver, find_on_path, live_provenance},
    error::RunPhase,
    subprocess, Adapter, AdapterError, AdapterInfo, Capabilities, Case, LicenseMode, Physics,
    PreparedJob, ProbeReport, RunContext, RunReport, VersionRange,
};
use valenx_fields::{
    artifact::{Artifact, ArtifactKind},
    Results,
};

use crate::case_input::ChimeraXInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(ChimeraXAdapter::new())
}

pub struct ChimeraXAdapter;

impl ChimeraXAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for ChimeraXAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "chimerax";
/// ChimeraX binary candidates. The lowercase form is the standard
/// install on Linux / macOS; the CamelCase form matches Windows
/// installers and some macOS app-bundle wrappers.
const BINARIES: &[&str] = &["chimerax", "ChimeraX"];

impl Adapter for ChimeraXAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "ChimeraX",
            // ChimeraX 1.7 (late-2023) is the first release with the
            // stable `--script` headless invocation we lean on; the
            // upper bound bumps when a 2.x line lands.
            version_range: VersionRange {
                min_inclusive: Version::new(1, 7, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "University of California — free for non-commercial / academic use",
            docs_url: "https://www.cgl.ucsf.edu/chimerax/docs/user/index.html",
            homepage_url: "https://www.cgl.ucsf.edu/chimerax/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // ChimeraX prints its version on `--version`.
                let found_version =
                    detect_tool_version_semver(&binary_path, &["--version", "-version"]);
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    warnings: Vec::new(),
                    required_env: Vec::new(),
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "ChimeraX 1.7+ required; download from \
                       https://www.cgl.ucsf.edu/chimerax/download.html"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = ChimeraXInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve the script path against the case directory. The user
        // authors `script = "render.cxc"` and expects it to live
        // alongside `case.toml`. `confined_join` rejects absolute paths
        // and `..` traversal so a malicious case bundle can't smuggle
        // arbitrary host files into the workdir.
        let source_script = confined_join(&case.path, &input.script)?;
        if !source_script.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.chimerax].script `{}` not found (resolved {})",
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
                        "[bio.chimerax].script path `{}` has no filename",
                        input.script.display()
                    ),
                })?;
        let dest_script = workdir.join(script_filename);
        if source_script != dest_script {
            fs::copy(&source_script, &dest_script)?;
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "ChimeraX 1.7+ required; download from \
                       https://www.cgl.ucsf.edu/chimerax/download.html"
                .into(),
        })?;

        // Build the command. `--nogui` keeps the run headless; the
        // rare interactive path (recording a session on a workstation)
        // drops the flag entirely so ChimeraX brings up its main
        // window.
        let mut native_command: Vec<OsString> = vec![binary_path.into_os_string()];
        if input.nogui {
            native_command.push(OsString::from("--nogui"));
        }
        native_command.push(OsString::from("--script"));
        native_command.push(OsString::from(script_filename));

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // ChimeraX rendering jobs typically finish in seconds to
            // a couple of minutes; complex animations can stretch to
            // tens of minutes. 30 minutes is a generous default that
            // covers the long tail without being absurd.
            estimated_runtime: Some(Duration::from_secs(30 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting ChimeraX", |line| {
            let mut hint = subprocess::Hint::default();
            // ChimeraX banners are loose — the script can emit
            // arbitrary text. We pick three weak signals as
            // best-effort progress hints:
            //   * "Saving" — the script reached an output write
            //   * "Rendering" — a render job started
            //   * "exit" — ChimeraX is shutting down
            // These are heuristics; mismatches just leave the
            // spinner alone.
            if line.contains("exit") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("Saving") {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.contains("Rendering") {
                hint.progress = Some((75.0, line.to_string()));
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
        // Provenance: hash the staged .cxc script (the canonical
        // "this case is configured this way" input). We don't know
        // the user's mesh / lock files so leave those empty.
        let script_path = first_script_in_workdir(&job.workdir);
        let case_hash_input = script_path
            .clone()
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "ChimeraX",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);

        // Walk the workdir top level. ChimeraX scripts conventionally
        // write outputs to the working directory; deeply nested
        // outputs are unusual and surface via the script's own
        // explicit `cd` / `save` paths.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-chimerax", ?e, "workdir read failed");
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
                // Rendered images — ChimeraX `save` with a .png
                // extension is the canonical "give me a picture"
                // path. Surfaces as Image so the UI can preview.
                Some("png") => (ArtifactKind::Image, "ChimeraX render".to_string()),
                // .cxs — ChimeraX session file. Re-openable via
                // `open session.cxs` in another ChimeraX run.
                Some("cxs") => (ArtifactKind::Native, "ChimeraX session".to_string()),
                // Structure outputs the script wrote (typically via
                // `save out.pdb` or `save out.cif`). Soft-validate
                // the PDB but never fail collect — surface a parse
                // warning in the label if the file is broken.
                Some("pdb") => {
                    let stem = path
                        .file_stem()
                        .and_then(|s| s.to_str())
                        .unwrap_or("chimerax")
                        .to_string();
                    // Round-22 M2: cap the read at MAX_PDB_FILE_BYTES (256
                    // MiB) so a poisoned workdir with a multi-GB `.pdb`
                    // can't OOM `collect()` before the parser runs.
                    let label = match valenx_core::io_caps::read_capped_to_string(
                        &path,
                        valenx_core::io_caps::MAX_PDB_FILE_BYTES as usize,
                    ) {
                        Ok(text) => match valenx_bio::format::pdb::read(&stem, &text) {
                            Ok(structure) => format!(
                                "ChimeraX PDB `{}` ({} atoms, {} residues)",
                                stem,
                                structure.atom_count(),
                                structure.residue_count()
                            ),
                            Err(e) => format!(
                                "ChimeraX PDB `{}` (parse warning: {})",
                                stem,
                                e.to_string().lines().next().unwrap_or("invalid")
                            ),
                        },
                        Err(_) => format!("ChimeraX PDB `{stem}`"),
                    };
                    (ArtifactKind::Native, label)
                }
                Some("cif") => (ArtifactKind::Native, "ChimeraX CIF structure".to_string()),
                Some("cxc") => (ArtifactKind::Other, "ChimeraX command script".to_string()),
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
            ribbon_contributions: vec!["bio.chimerax.script"],
        }
    }
}

/// Lift the staged ChimeraX script out of a workdir for provenance
/// hashing. Returns the lexicographically-first `.cxc` file at the
/// top level, or `None` if none exists yet.
fn first_script_in_workdir(workdir: &Path) -> Option<PathBuf> {
    let entries = fs::read_dir(workdir).ok()?;
    let mut hits: Vec<PathBuf> = entries
        .flatten()
        .map(|e| e.path())
        .filter(|p| {
            p.extension()
                .and_then(|s| s.to_str())
                .map(|s| s.eq_ignore_ascii_case("cxc"))
                .unwrap_or(false)
        })
        .collect();
    hits.sort();
    hits.into_iter().next()
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
        let info = ChimeraXAdapter::new().info();
        assert_eq!(info.id, "chimerax");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert!(info.tool_license.starts_with("University of California"));
        assert_eq!(info.display_name, "ChimeraX");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = ChimeraXAdapter::new().info();
        // ChimeraX >= 1.7 (stable --script headless mode); upper
        // bound 2.0 leaves room for the 1.x line and bumps when 2.x
        // lands.
        assert_eq!(info.version_range.min_inclusive, Version::new(1, 7, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(2, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = ChimeraXAdapter::new().capabilities();
        // Capability variants land in a future task; ribbon
        // contributions are already enough for the registry to
        // surface the adapter.
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.chimerax.script"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = ChimeraXAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    #[test]
    fn first_script_in_workdir_picks_lexicographic_first() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-chimerax-script-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        fs::write(tmp.join("z_late.cxc"), b"# placeholder").unwrap();
        fs::write(tmp.join("a_first.cxc"), b"# placeholder").unwrap();
        fs::write(tmp.join("not_cxc.txt"), b"placeholder").unwrap();
        let f = first_script_in_workdir(&tmp).expect("found");
        assert!(f.ends_with("a_first.cxc"));
        let _ = fs::remove_dir_all(&tmp);
    }

    /// `collect()` runs against a workdir directly — exercise the
    /// PNG / CXS / PDB / CIF classification paths so a regression in
    /// the extension-dispatch table doesn't slip past CI.
    #[test]
    fn collect_classifies_chimerax_outputs() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-chimerax-collect-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        fs::write(tmp.join("render.cxc"), b"# placeholder").unwrap();
        fs::write(tmp.join("snapshot.png"), b"\x89PNG\r\n\x1a\n").unwrap();
        fs::write(tmp.join("scene.cxs"), b"fake session bytes").unwrap();
        fs::write(tmp.join("model.pdb"), SAMPLE_PDB).unwrap();
        fs::write(tmp.join("model.cif"), b"# placeholder CIF").unwrap();
        fs::write(tmp.join("ignore.bin"), b"...").unwrap();

        let job = PreparedJob {
            workdir: tmp.clone(),
            native_command: vec![],
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: true,
        };
        let results = ChimeraXAdapter::new().collect(&job).unwrap();
        let labels: Vec<&str> = results.artifacts.iter().map(|a| a.label.as_str()).collect();
        assert!(labels.contains(&"ChimeraX render"));
        assert!(labels.contains(&"ChimeraX session"));
        assert!(labels.iter().any(|l| l.starts_with("ChimeraX PDB `model`")));
        assert!(labels.contains(&"ChimeraX CIF structure"));
        assert!(labels.contains(&"ChimeraX command script"));
        // .bin must not surface — guards the deny-by-default path.
        assert!(!results
            .artifacts
            .iter()
            .any(|a| a.path.extension().is_some_and(|e| e == "bin")));
        let _ = fs::remove_dir_all(&tmp);
    }

    /// PNG renders classify as `Image` (not Native); CXS sessions
    /// and structure files (PDB / CIF) classify as `Native`. Pin the
    /// contract.
    #[test]
    fn collect_artifact_kinds_pin_image_vs_native() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-chimerax-kinds-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        fs::write(tmp.join("a.png"), b"").unwrap();
        fs::write(tmp.join("b.cxs"), b"").unwrap();
        fs::write(tmp.join("c.pdb"), SAMPLE_PDB).unwrap();
        fs::write(tmp.join("d.cif"), b"").unwrap();

        let job = PreparedJob {
            workdir: tmp.clone(),
            native_command: vec![],
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: true,
        };
        let results = ChimeraXAdapter::new().collect(&job).unwrap();
        for art in &results.artifacts {
            let ext = art.path.extension().and_then(|e| e.to_str()).unwrap();
            match ext {
                "png" => assert_eq!(art.kind, ArtifactKind::Image),
                "cxs" | "pdb" | "cif" => assert_eq!(art.kind, ArtifactKind::Native),
                _ => panic!("unexpected extension {ext}"),
            }
        }
        let _ = fs::remove_dir_all(&tmp);
    }

    /// PDB output from a ChimeraX `save out.pdb` should parse via
    /// the canonical `valenx_bio` reader — verify the label includes
    /// the atom + residue summary.
    #[test]
    fn collect_pdb_parses_via_valenx_bio() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-chimerax-pdb-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        fs::write(tmp.join("rendered.pdb"), SAMPLE_PDB).unwrap();

        let job = PreparedJob {
            workdir: tmp.clone(),
            native_command: vec![],
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: true,
        };
        let results = ChimeraXAdapter::new().collect(&job).unwrap();
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
        let _ = fs::remove_dir_all(&tmp);
    }
}
