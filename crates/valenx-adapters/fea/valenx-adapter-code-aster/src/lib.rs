//! # valenx-adapter-code-aster
//!
//! Adapter for EDF's Code_Aster — industrial-strength FEA
//! (thermomechanics, dynamics, fatigue, contact, crack propagation).
//!
//! **Phase 3 — live for `as_run` on a user-provided .export.**
//! `prepare()` parses `[fea.code_aster]` from case.toml, stages
//! the .export plus all companion files (.comm / .mmed / .med /
//! .py) into the workdir, and builds `as_run case.export`.
//! `run()` spawns it. `collect()` walks the workdir for
//! result outputs.

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

use crate::case_input::CodeAsterInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(CodeAsterAdapter::new())
}

pub struct CodeAsterAdapter;

impl CodeAsterAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for CodeAsterAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "code-aster";
const BINARIES: &[&str] = &["as_run", "run_aster", "aster"];

impl Adapter for CodeAsterAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "Code_Aster",
            version_range: VersionRange {
                // Code_Aster's canonical long-support releases track 17.x.
                min_inclusive: Version::new(16, 4, 0),
                max_exclusive: Version::new(19, 0, 0),
            },
            physics: &[Physics::Fea, Physics::MultiPhysics],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-3.0-or-later",
            docs_url: "https://www.code-aster.org/V2/doc/default/en/",
            homepage_url: "https://www.code-aster.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                let found_version = valenx_core::adapter_helpers::detect_tool_version_semver(
                    &binary_path,
                    &["--version", "-v"],
                );
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
                hint: "Code_Aster 16.4+ required; install via Salome-Meca or \
                       your distribution"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = CodeAsterInput::from_case_dir(&case.path)?;
        fs::create_dir_all(workdir)?;

        let export_source = valenx_core::adapter_helpers::confined_join(
            &case.path,
            &input.export,
        )?;
        if !export_source.is_file() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "export not found at {} (resolve relative to case dir)",
                export_source.display()
            )));
        }
        let export_name = export_source
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("case.export")
            .to_string();
        fs::copy(&export_source, workdir.join(&export_name)).map_err(|e| {
            AdapterError::Other(anyhow::anyhow!("stage {}: {e}", export_source.display()))
        })?;

        // Stage every other file in the case dir — .comm files, .mmed
        // / .med meshes, additional Python helpers. The .export file
        // references most of these by basename so they need to be
        // workdir-local.
        if let Some(parent) = export_source.parent() {
            stage_companion_files(parent, &export_name, workdir);
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "no `as_run` / `run_aster` / `aster` on PATH".into(),
        })?;

        let native_command: Vec<OsString> =
            vec![binary_path.into_os_string(), OsString::from(&export_name)];

        // Industrial Code_Aster runs span minutes-to-hours. 4-hour
        // default ceiling for the GUI cancellation timer.
        let estimated_runtime = Some(Duration::from_secs(4 * 60 * 60));

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            estimated_runtime,
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting Code_Aster as_run", |line| {
            let mut hint = subprocess::Hint::default();
            if let Some(pct) = code_aster_progress_hint(line) {
                hint.progress = Some((pct, line.to_string()));
            }
            // Code_Aster's hard failures show up as "<F>" or
            // "<EXCEPTION>" messages; surface them as warnings.
            if line.contains("<F>") || line.contains("<EXCEPTION>") || line.contains("Erreur") {
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
        // case_path: prefer .comm (the script the export references);
        // fall back to the .export itself.
        let case_path = first_workdir_match(&job.workdir, &["comm"])
            .or_else(|| first_workdir_match(&job.workdir, &["export"]))
            .unwrap_or_else(|| job.workdir.join("(no-comm-found)"));
        let mesh_path = first_workdir_match(&job.workdir, &["med", "mmed"]);
        let prov = valenx_core::adapter_helpers::live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "Code_Aster",
            "unknown",
            &case_path,
            mesh_path.as_deref(),
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);

        // Workdir scan. Code_Aster's typical output set:
        //   .resu / .resu.med — result files
        //   .med              — input or output mesh
        //   .mess             — text run summary
        //   .erre             — text error log
        //   .o*               — text output (per-step)
        let classifications: &[(&str, ArtifactKind, &str)] = &[
            ("resu", ArtifactKind::Native, "Code_Aster .resu result"),
            ("med", ArtifactKind::Native, "Code_Aster MED file"),
            ("mmed", ArtifactKind::Native, "Code_Aster MMED file"),
            ("comm", ArtifactKind::Other, "Code_Aster command file"),
            (
                "export",
                ArtifactKind::Other,
                "Code_Aster export descriptor",
            ),
            ("py", ArtifactKind::Other, "Code_Aster Python helper"),
            ("mess", ArtifactKind::Log, "Code_Aster .mess summary"),
            ("erre", ArtifactKind::Log, "Code_Aster .erre error log"),
            ("log", ArtifactKind::Log, "Code_Aster log"),
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
            capabilities: vec![
                Capability::FeaLinearStatic,
                Capability::FeaNonlinearStatic,
                Capability::FeaModal,
                Capability::FeaHarmonic,
                Capability::FeaTransient,
                Capability::FeaThermal,
                Capability::FeaContact,
            ],
            ribbon_contributions: vec![
                "fea.codeaster.static",
                "fea.codeaster.dynamic",
                "fea.codeaster.fatigue",
            ],
        }
    }
}

/// Stage every file alongside the .export (.comm, .mmed, .med,
/// helper .py) into the workdir. The .export references most of
/// these by basename so they need to be workdir-local.
fn stage_companion_files(source_dir: &Path, export_name: &str, workdir: &Path) {
    let Ok(entries) = fs::read_dir(source_dir) else {
        return;
    };
    for entry in entries.flatten() {
        let p = entry.path();
        let Some(name) = p.file_name().and_then(|s| s.to_str()) else {
            continue;
        };
        if name == export_name || name == "case.toml" {
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
            tracing::warn!(target: "valenx.code-aster", ?p, ?dst, %e, "stage companion failed");
        }
    }
}

/// Coarse progress hints for `as_run` stdout. Code_Aster prints
/// recognisable section banners + per-step diagnostics.
fn code_aster_progress_hint(line: &str) -> Option<f32> {
    if line.contains("EXECUTION_COMMENCEE") || line.contains("DEBUT EXECUTION") {
        Some(10.0)
    } else if line.contains("LECTURE DU MAILLAGE") || line.contains("MAIL_") {
        Some(25.0)
    } else if line.contains("ASSEMBLAGE") || line.contains("CALC_MATR") {
        Some(40.0)
    } else if line.contains("RESOLUTION") || line.contains("Newton") {
        Some(60.0)
    } else if line.contains("INST :") || line.contains("Instant") {
        // Per-time-step banner — stay flat at 70 % during the loop.
        Some(70.0)
    } else if line.contains("SAUVEGARDE") || line.contains("ECRITURE") {
        Some(90.0)
    } else if line.contains("EXECUTION_TERMINEE") || line.contains("FIN EXECUTION") {
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
    fn info_mentions_multi_physics() {
        let info = CodeAsterAdapter::new().info();
        assert!(info.physics.contains(&Physics::MultiPhysics));
    }

    #[test]
    fn collect_uses_live_provenance_with_real_case_hash() {
        // Drop a known .comm file in a tempdir, run collect(), and
        // verify the provenance block carries: real adapter id + the
        // crate version, real tool name, a non-empty run id, and a
        // case_hash that matches sha256 of the .comm bytes. This
        // catches anyone reverting collect() back to stub_provenance.
        let workdir = tempdir("code-aster-collect");
        let comm_path = workdir.join("dummy.comm");
        let comm_bytes = b"DEBUT()\nFIN()\n";
        std::fs::write(&comm_path, comm_bytes).expect("write .comm");

        let job = PreparedJob {
            workdir: workdir.clone(),
            native_command: Vec::new(),
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: false,
        };
        let results = CodeAsterAdapter::new().collect(&job).expect("collect");
        let prov = &results.provenance;

        assert_eq!(prov.adapter, INFO_ID);
        assert!(!prov.adapter_version.is_empty(), "adapter_version empty");
        assert_eq!(prov.tool, "Code_Aster");
        assert!(!prov.run_id.is_empty(), "run_id empty — stub still wired?");
        let expected_hash = sha256_hex_file(&comm_path);
        assert_eq!(
            prov.case_hash, expected_hash,
            "case_hash should match SHA-256 of the .comm file"
        );

        cleanup(&workdir);
    }

    fn cleanup(d: &std::path::Path) {
        let _ = std::fs::remove_dir_all(d);
    }

    #[test]
    fn progress_hints_are_monotonic() {
        let pts = [
            code_aster_progress_hint("EXECUTION_COMMENCEE"),
            code_aster_progress_hint("LECTURE DU MAILLAGE"),
            code_aster_progress_hint("ASSEMBLAGE"),
            code_aster_progress_hint("RESOLUTION"),
            code_aster_progress_hint("INST : 0.5"),
            code_aster_progress_hint("ECRITURE"),
            code_aster_progress_hint("EXECUTION_TERMINEE"),
        ];
        let mut last = 0.0_f32;
        for (i, p) in pts.iter().enumerate() {
            let v = p.unwrap_or_else(|| panic!("step {i} returned None"));
            assert!(v >= last, "step {i}: {last} -> {v}");
            last = v;
        }
    }

    #[test]
    fn collect_classifies_resu_med_and_logs() {
        let workdir = tempdir("code-aster-collect-class");
        for (name, content) in [
            ("case.resu", &b"binary"[..]),
            ("case.med", &b"binary"[..]),
            ("case.comm", &b"DEBUT()\n"[..]),
            ("case.mess", &b"text summary"[..]),
            ("case.erre", &b"text errors"[..]),
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
        let results = CodeAsterAdapter::new().collect(&job).expect("collect");
        // 5 known artifacts; ignored.txt skipped.
        assert_eq!(results.artifacts.len(), 5);
        let labels: Vec<&str> = results.artifacts.iter().map(|a| a.label.as_str()).collect();
        assert!(labels.iter().any(|l| l.contains(".resu result")));
        assert!(labels.iter().any(|l| l.contains("MED")));
        assert!(labels.iter().any(|l| l.contains(".mess")));
        assert!(labels.iter().any(|l| l.contains(".erre")));
        cleanup(&workdir);
    }

    #[test]
    fn prepare_stages_export_and_companions() {
        let case_dir = tempdir("code-aster-aster-prepare");
        std::fs::write(
            case_dir.join("case.toml"),
            "[fea.code_aster]\nexport = \"sim.export\"\n",
        )
        .unwrap();
        std::fs::write(case_dir.join("sim.export"), b"P actions make_etude\n").unwrap();
        std::fs::write(case_dir.join("sim.comm"), b"DEBUT()\nFIN()\n").unwrap();
        std::fs::write(case_dir.join("mesh.med"), b"binary mesh").unwrap();
        let workdir = tempdir("code-aster-aster-prepare-wd");
        let case = Case {
            id: "aster-test".into(),
            path: case_dir.clone(),
        };
        let r = CodeAsterAdapter::new().prepare(&case, &workdir);
        if find_on_path(BINARIES).is_none() {
            assert!(matches!(r, Err(AdapterError::ToolNotInstalled { .. })));
            cleanup(&case_dir);
            cleanup(&workdir);
            return;
        }
        let job = r.expect("prepare");
        // Export + .comm + .med all staged.
        assert!(workdir.join("sim.export").is_file());
        assert!(workdir.join("sim.comm").is_file());
        assert!(workdir.join("mesh.med").is_file());
        // Native command: as_run sim.export
        let cmd: Vec<String> = job
            .native_command
            .iter()
            .map(|s| s.to_string_lossy().into_owned())
            .collect();
        assert!(cmd.iter().any(|s| s == "sim.export"));
        cleanup(&case_dir);
        cleanup(&workdir);
    }

    #[test]
    fn prepare_missing_export_is_actionable() {
        let case_dir = tempdir("code-aster-aster-no-export");
        std::fs::write(
            case_dir.join("case.toml"),
            "[fea.code_aster]\nexport = \"missing.export\"\n",
        )
        .unwrap();
        let workdir = tempdir("code-aster-aster-no-export-wd");
        let case = Case {
            id: "aster-test".into(),
            path: case_dir.clone(),
        };
        assert!(CodeAsterAdapter::new().prepare(&case, &workdir).is_err());
        cleanup(&case_dir);
        cleanup(&workdir);
    }
}
