//! # valenx-adapter-openradioss
//!
//! Adapter for OpenRadioss (Altair) — explicit nonlinear structural
//! dynamics for crash and impact.
//!
//! **Phase 3 — live for the engine phase.** `prepare()` parses
//! `[fea.openradioss]` from case.toml, stages the engine deck (and
//! any siblings) into the workdir, and builds an `engine_<arch> -i
//! <deck> -nspmd N -nthread M` invocation. `run()` spawns the
//! engine via the shared subprocess runner. `collect()` walks the
//! workdir for `.h3d` / `.anim` / `.rst` outputs.
//!
//! Scope today: the engine phase only. The starter→engine
//! conversion is left to the user because it's a one-time step
//! that runs on a workstation while the engine phase is what gets
//! queued on a cluster — splitting them matches actual workflows.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::fs;
use std::path::Path;
use std::time::Duration;

use semver::Version;

use valenx_core::{
    adapter_helpers::{find_on_path, first_workdir_match},
    error::RunPhase,
    subprocess, Adapter, AdapterError, AdapterInfo, Capabilities, Capability, Case, LicenseMode,
    Physics, PreparedJob, ProbeReport, RunContext, RunReport, VersionRange,
};
use valenx_fields::{
    artifact::{Artifact, ArtifactKind},
    Results,
};

use crate::case_input::OpenRadiossInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(OpenRadiossAdapter::new())
}

pub struct OpenRadiossAdapter;

impl OpenRadiossAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for OpenRadiossAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "openradioss";
/// Engine binary — the simulation phase, the part the adapter
/// drives. Common names across distributions.
const ENGINE_BINARIES: &[&str] = &[
    "engine_linux64_gf",
    "engine_linux64",
    "engine_win64",
    "engine_macos64",
    "openradioss",
];

impl Adapter for OpenRadiossAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "OpenRadioss",
            version_range: VersionRange {
                min_inclusive: Version::new(2023, 0, 0),
                max_exclusive: Version::new(2030, 0, 0),
            },
            physics: &[Physics::Fea],
            license_mode: LicenseMode::Subprocess,
            tool_license: "AGPL-3.0-only",
            docs_url: "https://openradioss.atlassian.net/wiki/",
            homepage_url: "https://www.openradioss.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(ENGINE_BINARIES) {
            Some(binary_path) => {
                let found_version = valenx_core::adapter_helpers::detect_tool_version_semver(
                    &binary_path,
                    &["--version", "-v"],
                );
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    warnings: vec!["openradioss adapter runs the engine phase only — \
                         user is responsible for the starter -> engine \
                         conversion (typically `starter_<arch> -i \
                         model_0000.rad`)"
                        .into()],
                    required_env: Vec::new(),
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "OpenRadioss engine required; install from openradioss.org \
                       (binary names like `engine_linux64_gf` / `engine_win64`)"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = OpenRadiossInput::from_case_dir(&case.path)?;
        fs::create_dir_all(workdir)?;

        // Stage the engine deck. Resolve relative path against the
        // case directory.
        let source = valenx_core::adapter_helpers::confined_join(
            &case.path,
            &input.engine_input,
        )?;
        if !source.is_file() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "engine_input not found at {} (resolve relative to case dir)",
                source.display()
            )));
        }
        let staged_name = source
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("engine_input.rad")
            .to_string();
        let staged = workdir.join(&staged_name);
        fs::copy(&source, &staged)
            .map_err(|e| AdapterError::Other(anyhow::anyhow!("stage {}: {e}", source.display())))?;
        // Best-effort: stage any sibling restart / include files
        // sitting in the same directory as the engine deck (e.g.
        // `<root>_0001.rst`). Errors here are warnings — the engine
        // tells us which extra files it actually needs.
        if let Some(parent) = source.parent() {
            stage_sibling_inputs(parent, &staged_name, workdir);
        }

        // Find the engine binary. Probe already checked it exists,
        // but we re-find here so the resolved path is the one we
        // actually invoke (PATH could change between probe and run).
        let binary_path =
            find_on_path(ENGINE_BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "no OpenRadioss engine binary on PATH".into(),
            })?;

        let native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from("-i"),
            OsString::from(&staged_name),
            OsString::from("-nspmd"),
            OsString::from(input.nspmd.to_string()),
            OsString::from("-nthread"),
            OsString::from(input.nthread.to_string()),
        ];

        // Crash sims of single drop tests run in seconds; full
        // vehicle sims can take hours. 30 minutes is a generous
        // default UI ceiling for the cancellation timer.
        let estimated_runtime = Some(Duration::from_secs(30 * 60));

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            estimated_runtime,
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting OpenRadioss engine", |line| {
            let mut hint = subprocess::Hint::default();
            if let Some(pct) = openradioss_progress_hint(line) {
                hint.progress = Some((pct, line.to_string()));
            }
            // OpenRadioss emits ERROR / WARNING tokens on its own
            // diagnostics; surface them as warnings so the GUI
            // residual panel picks them up.
            if line.contains("ERROR ID") || line.contains("WARNING ID") {
                hint.warning = Some(line.trim().to_string());
            }
            hint
        })?;
        Ok(RunReport {
            exit_code: report.exit_code,
            wall_time: report.wall_time,
            // OpenRadioss is explicit time integration — no
            // residual-based convergence concept. We mark
            // converged=true on a clean exit; a non-zero exit
            // already routes through subprocess::run's error path.
            converged: Some(true),
            residual_history: Vec::new(),
            warnings: report.warnings,
            final_phase: Some(RunPhase::Shutdown),
        })
    }

    fn collect(&self, job: &PreparedJob) -> Result<Results, AdapterError> {
        // case_path = first .rad keyword deck — typically the
        // engine deck we staged in prepare(). live_provenance
        // hashes its content for the audit trail.
        let case_path = first_workdir_match(&job.workdir, &["rad"])
            .unwrap_or_else(|| job.workdir.join("(no-rad-found)"));
        let prov = valenx_core::adapter_helpers::live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "OpenRadioss",
            "unknown",
            &case_path,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);

        // Walk the workdir top-level for OpenRadioss outputs.
        // Classification: .h3d / .anim are visualisation dumps,
        // .rst is the binary restart, .out / .log is text log.
        let classifications: &[(&str, ArtifactKind, &str)] = &[
            ("h3d", ArtifactKind::Native, "OpenRadioss H3D animation"),
            ("anim", ArtifactKind::Native, "OpenRadioss anim file"),
            ("rst", ArtifactKind::Native, "OpenRadioss restart"),
            ("out", ArtifactKind::Log, "OpenRadioss text output"),
            ("log", ArtifactKind::Log, "OpenRadioss log"),
            ("rad", ArtifactKind::Other, "OpenRadioss keyword deck"),
        ];
        if let Ok(entries) = fs::read_dir(&job.workdir) {
            for entry in entries.flatten() {
                let path = entry.path();
                let ext = path
                    .extension()
                    .and_then(|s| s.to_str())
                    .map(|s| s.to_ascii_lowercase());
                let Some(ext) = ext else { continue };
                if let Some(&(_, kind, label)) =
                    classifications.iter().find(|(e, ..)| *e == ext.as_str())
                {
                    results.artifacts.push(Artifact {
                        path,
                        kind,
                        checksum: None,
                        label: label.to_string(),
                    });
                }
            }
        }
        results.artifacts.sort_by(|a, b| a.path.cmp(&b.path));
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: vec![Capability::FeaTransient, Capability::FeaContact],
            ribbon_contributions: vec!["fea.openradioss.crash", "fea.openradioss.impact"],
        }
    }
}

/// Best-effort stage of files that sit alongside the engine deck.
/// Skips the engine deck itself (already copied), subdirectories,
/// and our own preview / log artifacts. Errors surface as tracing
/// warnings rather than failing prepare — the engine will tell us
/// at run time which files it actually needs.
fn stage_sibling_inputs(source_dir: &Path, staged_name: &str, workdir: &Path) {
    let Ok(entries) = fs::read_dir(source_dir) else {
        return;
    };
    for entry in entries.flatten() {
        let p = entry.path();
        let Some(name) = p.file_name().and_then(|s| s.to_str()) else {
            continue;
        };
        if name == staged_name || name == "case.toml" {
            continue;
        }
        if entry.file_type().map(|t| !t.is_file()).unwrap_or(true) {
            continue;
        }
        let dst = workdir.join(name);
        if dst.exists() {
            continue;
        }
        if let Err(e) = fs::copy(&p, &dst) {
            tracing::warn!(target: "valenx.openradioss", ?p, ?dst, %e, "stage sibling failed");
        }
    }
}

/// Coarse progress hints for OpenRadioss engine stdout. Based on the
/// banners the engine emits in default verbosity. The percentages
/// approximate "where in the sim are we"; non-monotonic messages
/// are deliberately not promoted to a hint so the progress bar
/// only moves forward.
fn openradioss_progress_hint(line: &str) -> Option<f32> {
    if line.contains("STARTER ENDED") || line.contains("STARTING ENGINE") {
        Some(5.0)
    } else if line.contains("ANALYSIS STARTED") {
        Some(10.0)
    } else if line.contains("CYCLE") && line.contains("TIME") {
        // Mid-run ticks. Keep at 50 so the bar doesn't oscillate
        // on every cycle line — the wall-time elapsed is the more
        // useful signal during the run anyway.
        Some(50.0)
    } else if line.contains("NORMAL TERMINATION") {
        Some(98.0)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_core::adapter_helpers::sha256_hex_file;
    use valenx_test_utils::tempdir;

    #[test]
    fn info_is_agpl() {
        let info = OpenRadiossAdapter::new().info();
        assert_eq!(info.tool_license, "AGPL-3.0-only");
    }

    #[test]
    fn collect_uses_live_provenance_with_real_case_hash() {
        let workdir = tempdir("openradioss-collect");
        let case_path = workdir.join("model_0000.rad");
        let case_bytes = b"#RADIOSS STARTER\n/BEGIN\n/END\n";
        std::fs::write(&case_path, case_bytes).expect("write .rad");

        let job = PreparedJob {
            workdir: workdir.clone(),
            native_command: Vec::new(),
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: false,
        };
        let results = OpenRadiossAdapter::new().collect(&job).expect("collect");
        let prov = &results.provenance;

        assert_eq!(prov.adapter, INFO_ID);
        assert!(!prov.adapter_version.is_empty());
        assert_eq!(prov.tool, "OpenRadioss");
        assert!(!prov.run_id.is_empty(), "run_id empty — stub still wired?");
        assert_eq!(prov.case_hash, sha256_hex_file(&case_path));

        cleanup(&workdir);
    }

    fn cleanup(d: &std::path::Path) {
        let _ = std::fs::remove_dir_all(d);
    }

    #[test]
    fn progress_hints_are_monotonic() {
        // Coarse banners we recognise should produce non-decreasing
        // percentages so the GUI bar only moves forward.
        let pts = [
            openradioss_progress_hint("STARTER ENDED OK"),
            openradioss_progress_hint("ANALYSIS STARTED"),
            openradioss_progress_hint("CYCLE 100   TIME 0.123E-03"),
            openradioss_progress_hint("NORMAL TERMINATION"),
        ];
        let mut last = 0.0_f32;
        for (i, p) in pts.iter().enumerate() {
            let v = p.unwrap_or_else(|| panic!("step {i} returned None"));
            assert!(v >= last, "step {i}: {last} -> {v}");
            last = v;
        }
    }

    #[test]
    fn collect_classifies_h3d_anim_rst_outputs() {
        // Drop a pile of fake outputs in the workdir; collect()
        // should classify each by extension and return them sorted
        // by path.
        let workdir = tempdir("openradioss-orad-collect");
        for (name, content) in [
            ("model.h3d", &b"binary"[..]),
            ("model.anim", &b"binary"[..]),
            ("model_0001.rst", &b"binary"[..]),
            ("model_0001.rad", &b"#RADIOSS\n"[..]),
            ("model.out", &b"text log"[..]),
            ("ignored.txt", &b"unrelated"[..]),
        ] {
            std::fs::write(workdir.join(name), content).unwrap();
        }
        let job = PreparedJob {
            workdir: workdir.clone(),
            native_command: Vec::new(),
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: false,
        };
        let results = OpenRadiossAdapter::new().collect(&job).expect("collect");
        // h3d, anim, rst, rad, out — 5 known artifacts; ignored.txt
        // should be skipped.
        assert_eq!(results.artifacts.len(), 5);
        let labels: Vec<&str> = results.artifacts.iter().map(|a| a.label.as_str()).collect();
        assert!(labels.iter().any(|l| l.contains("H3D animation")));
        assert!(labels.iter().any(|l| l.contains("anim file")));
        assert!(labels.iter().any(|l| l.contains("restart")));
        assert!(labels.iter().any(|l| l.contains("text output")));
        assert!(labels.iter().any(|l| l.contains("keyword deck")));
        cleanup(&workdir);
    }

    #[test]
    fn prepare_stages_engine_deck_into_workdir() {
        // Set up a case dir with a real engine_input file + case.toml,
        // then call prepare(). Need an `engine_<arch>` binary on PATH
        // for find_on_path to succeed; if none is present, the test
        // is skipped (find_on_path returns None) so this works on
        // CI without OpenRadioss installed.
        if find_on_path(ENGINE_BINARIES).is_none() {
            // No engine binary present — exercise the error path.
            let case_dir = tempdir("openradioss-orad-prepare-no-bin");
            std::fs::write(
                case_dir.join("case.toml"),
                r#"
[case]
physics = "fea"

[fea.openradioss]
engine_input = "deck_0001.rad"
"#,
            )
            .unwrap();
            std::fs::write(case_dir.join("deck_0001.rad"), b"#RADIOSS\n").unwrap();
            let workdir = tempdir("openradioss-orad-prepare-no-bin-wd");
            let case = Case {
                id: "orad-test".into(),
                path: case_dir.clone(),
            };
            let r = OpenRadiossAdapter::new().prepare(&case, &workdir);
            assert!(matches!(r, Err(AdapterError::ToolNotInstalled { .. })));
            cleanup(&case_dir);
            cleanup(&workdir);
            return;
        }
        // Engine binary IS on PATH — prepare should stage the deck
        // and return a runnable PreparedJob.
        let case_dir = tempdir("openradioss-orad-prepare");
        std::fs::write(
            case_dir.join("case.toml"),
            r#"
[case]
physics = "fea"

[fea.openradioss]
engine_input = "deck_0001.rad"
nspmd = 2
nthread = 4
"#,
        )
        .unwrap();
        std::fs::write(case_dir.join("deck_0001.rad"), b"#RADIOSS\n").unwrap();
        // A sibling restart file should also get staged.
        std::fs::write(case_dir.join("deck_0001.rst"), b"binary").unwrap();

        let workdir = tempdir("openradioss-orad-prepare-wd");
        let case = Case {
            id: "orad-test".into(),
            path: case_dir.clone(),
        };
        let job = OpenRadiossAdapter::new()
            .prepare(&case, &workdir)
            .expect("prepare");
        // Engine deck staged.
        assert!(workdir.join("deck_0001.rad").is_file());
        // Sibling restart staged.
        assert!(workdir.join("deck_0001.rst").is_file());
        // Native command structure: binary -i deck -nspmd 2 -nthread 4.
        let cmd_strs: Vec<String> = job
            .native_command
            .iter()
            .map(|s| s.to_string_lossy().into_owned())
            .collect();
        assert!(cmd_strs.contains(&"-i".into()));
        assert!(cmd_strs.contains(&"deck_0001.rad".into()));
        assert!(cmd_strs.contains(&"-nspmd".into()));
        assert!(cmd_strs.contains(&"2".into()));
        assert!(cmd_strs.contains(&"-nthread".into()));
        assert!(cmd_strs.contains(&"4".into()));

        cleanup(&case_dir);
        cleanup(&workdir);
    }
}
