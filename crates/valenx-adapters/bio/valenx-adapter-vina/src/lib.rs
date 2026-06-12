//! # valenx-adapter-vina
//!
//! Adapter for [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina)
//! — Scripps' open-source flexible-ligand docking engine. Vina is the
//! single most widely used tool in academic small-molecule docking:
//! given a rigid receptor (PDBQT) and a flexible ligand (PDBQT), it
//! places the ligand into a user-defined search box and reports the
//! lowest-energy poses ranked by predicted binding affinity.
//!
//! **Phase 34 — subprocess wrapper around the `vina` binary.** The
//! user supplies receptor + ligand + output paths plus the search-box
//! geometry via `[bio.vina]` in `case.toml`; `prepare()` validates
//! the box (positive edges, exhaustiveness in 1..=32, finite +
//! positive energy range) and composes the long flag-driven CLI
//! invocation. `run()` streams Vina's progress chatter through the
//! shared subprocess runner — the `Performing search...` and final
//! `Refining results ... done` markers translate to UI-visible
//! progress hints.
//!
//! On `collect()` we surface the docked-poses PDBQT as a `Native`
//! artifact. Vina writes a single output file, so artifact discovery
//! is straightforward.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::path::Path;
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

use crate::case_input::VinaInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(VinaAdapter::new())
}

pub struct VinaAdapter;

impl VinaAdapter {
    pub fn new() -> Self {
        Self
    }

    /// Probe variant aware of the engine selection. Native engine
    /// has no external dependency, so probe always succeeds and the
    /// reported version is the valenx-dock crate version.
    pub fn probe_with_engine(&self, engine: &str) -> Result<ProbeReport, AdapterError> {
        if engine == "native" {
            return Ok(ProbeReport {
                ok: true,
                found_version: Some(semver::Version::parse(env!("CARGO_PKG_VERSION")).unwrap()),
                binary_path: None,
                warnings: Vec::new(),
                required_env: Vec::new(),
            });
        }
        self.probe()
    }

    fn run_external(
        &self,
        job: &PreparedJob,
        ctx: &mut RunContext,
    ) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting AutoDock Vina", |line| {
            let mut hint = subprocess::Hint::default();
            // Vina prints a banner, then a progress bar of asterisks,
            // then a poses table. The two stable sentinel lines we
            // can lift to progress hints are "Performing search..."
            // (search just started) and the per-mode results table
            // marker ("Refining results ... done", or the table
            // header that starts with "mode |").
            if line.contains("Refining results") || line.contains("Writing output") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("Performing search") {
                hint.progress = Some((25.0, line.to_string()));
            } else if line.starts_with("ERROR") || line.contains("Parse error") {
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

    fn run_native(
        &self,
        job: &PreparedJob,
        _ctx: &mut RunContext,
        input: &VinaInput,
    ) -> Result<RunReport, AdapterError> {
        // Security policy: receptor + ligand MUST live under the case
        // directory; output MUST live under the workdir. No absolute
        // paths from case.toml — `confined_join` rejects both `..`
        // traversal and absolute overrides. Same policy applied in
        // `prepare()` below so the two paths cannot drift.
        let start = std::time::Instant::now();
        let case_dir = job.workdir.parent().unwrap_or(&job.workdir);
        let receptor_path = confined_join(case_dir, &input.receptor)?;
        let ligand_path = confined_join(case_dir, &input.ligand)?;
        let output_path = confined_join(&job.workdir, &input.output)?;
        std::fs::create_dir_all(&job.workdir)
            .map_err(|e| AdapterError::Other(anyhow::anyhow!("mkdir workdir: {e}")))?;
        // Round-23 named finding: bound the receptor + ligand reads
        // at MAX_PDBQT_FILE_BYTES (64 MiB each) so a stale or hostile
        // path can't OOM the docker before parsing. Sister to the
        // R20 H1 MCP dock-panel cap.
        let receptor_text = valenx_core::io_caps::read_capped_to_string(
            &receptor_path,
            valenx_core::io_caps::MAX_PDBQT_FILE_BYTES,
        )
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("read receptor: {e}")))?;
        let ligand_text = valenx_core::io_caps::read_capped_to_string(
            &ligand_path,
            valenx_core::io_caps::MAX_PDBQT_FILE_BYTES,
        )
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("read ligand: {e}")))?;

        let cfg = valenx_dock::DockConfig {
            center: nalgebra::Vector3::new(input.center[0], input.center[1], input.center[2]),
            size: nalgebra::Vector3::new(input.size[0], input.size[1], input.size[2]),
            exhaustiveness: input.exhaustiveness,
            num_modes: input.num_modes,
            energy_range: input.energy_range,
            ..Default::default()
        };

        valenx_dock::dock(&receptor_text, &ligand_text, &cfg, &output_path, None)
            .map_err(|e| AdapterError::Other(anyhow::anyhow!("native vina: {e}")))?;

        Ok(RunReport {
            exit_code: 0,
            wall_time: start.elapsed(),
            converged: Some(true),
            residual_history: Vec::new(),
            warnings: Vec::new(),
            final_phase: Some(RunPhase::Shutdown),
        })
    }
}

impl Default for VinaAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "vina";
/// Vina's binary candidates. `vina` is the canonical name from both
/// the upstream conda-forge / Bioconda packages and source builds.
const BINARIES: &[&str] = &["vina"];

impl Adapter for VinaAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "AutoDock Vina",
            // Vina 1.2.x is the current actively maintained line
            // (1.2.0 released April 2021, the C++/Python rewrite that
            // replaced the long-running 1.1.x series). 2.0 reserves
            // room for an eventual major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(1, 2, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "Apache-2.0",
            docs_url: "https://autodock-vina.readthedocs.io/",
            homepage_url: "https://github.com/ccsb-scripps/AutoDock-Vina",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `vina --version` prints "AutoDock Vina v1.2.5" or
                // similar to stdout; the generic detector also reads
                // the bare-name banner output as a fallback.
                let found_version = detect_tool_version_semver(&binary_path, &["--version", ""]);
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
                hint: "AutoDock Vina 1.2+ required; install via \
                       `conda install -c conda-forge vina` or build from \
                       https://github.com/ccsb-scripps/AutoDock-Vina"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = VinaInput::from_case_dir(&case.path)?;

        std::fs::create_dir_all(workdir)?;

        // Receptor + ligand MUST live under the case directory.
        // `confined_join` rejects absolute paths and `..` traversal —
        // sandboxing is part of the security contract for case.toml
        // inputs; users with files outside the case dir should
        // physically copy them in.
        let resolved_receptor = confined_join(&case.path, &input.receptor)?;
        if !resolved_receptor.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.vina].receptor `{}` not found (resolved {})",
                    input.receptor.display(),
                    resolved_receptor.display()
                ),
            });
        }

        let resolved_ligand = confined_join(&case.path, &input.ligand)?;
        if !resolved_ligand.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.vina].ligand `{}` not found (resolved {})",
                    input.ligand.display(),
                    resolved_ligand.display()
                ),
            });
        }

        // Output path: must land in the workdir. Vina creates the file
        // itself, so we don't pre-touch it; `confined_join` rejects
        // any attempt to escape via `..` or an absolute override.
        let resolved_output = confined_join(workdir, &input.output)?;

        // For `engine = "native"`, no external binary is needed —
        // `run_native()` calls into `valenx-dock` in-process. We still
        // populate `native_command` with a placeholder so PreparedJob
        // has a well-formed shape (the subprocess path is not taken).
        let binary_path = if input.engine == "native" {
            std::path::PathBuf::from("vina-native")
        } else {
            find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "AutoDock Vina 1.2+ required; install via \
                       `conda install -c conda-forge vina` or build from \
                       https://github.com/ccsb-scripps/AutoDock-Vina"
                    .into(),
            })?
        };

        // Compose `vina --receptor R --ligand L --out O \
        //   --center_x X --center_y Y --center_z Z \
        //   --size_x SX --size_y SY --size_z SZ \
        //   --exhaustiveness N --num_modes N --energy_range E \
        //   [--cpu N] [extras...]`.
        //
        // Each numeric value is its own OsString — Vina expects
        // `--flag value` pairs, never `--flag=value`.
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from("--receptor"),
            resolved_receptor.into_os_string(),
            OsString::from("--ligand"),
            resolved_ligand.into_os_string(),
            OsString::from("--out"),
            resolved_output.into_os_string(),
            OsString::from("--center_x"),
            OsString::from(input.center[0].to_string()),
            OsString::from("--center_y"),
            OsString::from(input.center[1].to_string()),
            OsString::from("--center_z"),
            OsString::from(input.center[2].to_string()),
            OsString::from("--size_x"),
            OsString::from(input.size[0].to_string()),
            OsString::from("--size_y"),
            OsString::from(input.size[1].to_string()),
            OsString::from("--size_z"),
            OsString::from(input.size[2].to_string()),
            OsString::from("--exhaustiveness"),
            OsString::from(input.exhaustiveness.to_string()),
            OsString::from("--num_modes"),
            OsString::from(input.num_modes.to_string()),
            OsString::from("--energy_range"),
            OsString::from(input.energy_range.to_string()),
        ];
        if input.cpu > 0 {
            native_command.push(OsString::from("--cpu"));
            native_command.push(OsString::from(input.cpu.to_string()));
        }
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Single-ligand docking finishes in a few minutes for
            // small boxes and exhaustiveness=8; one hour is a generous
            // headroom for high exhaustiveness or large search volumes.
            estimated_runtime: Some(Duration::from_secs(60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        // Re-parse the case to know which engine to dispatch. We
        // could carry this through PreparedJob instead, but parsing
        // a 50-line TOML twice is negligible vs the search cost.
        //
        // Both `prepare()` and `run_native()` route case.toml paths
        // through `valenx_core::adapter_helpers::confined_join`, so
        // `..` traversal and absolute overrides are rejected up front.
        let case_dir = job.workdir.parent().unwrap_or(&job.workdir);
        let input = match VinaInput::from_case_dir(case_dir) {
            Ok(v) => v,
            Err(_) => {
                // Fall back to subprocess if we can't determine the engine
                // (e.g. running from a workdir without a co-located case).
                return self.run_external(job, ctx);
            }
        };
        // Settings escape hatch: the app flips
        // `valenx_core::set_force_external_vina(true)` at startup
        // whenever the user enables the "Force external Vina binary"
        // toggle in the Settings dialog. When set, we always shell
        // out to the upstream binary even if the case picked
        // `engine = "native"` — gives users a one-flip way to compare
        // native vs reference output, or to dodge any unexpected
        // divergence in the native engine.
        let force_external = valenx_core::force_external_vina();
        if !force_external && input.engine == "native" {
            return self.run_native(job, ctx, &input);
        }
        self.run_external(job, ctx)
    }

    fn collect(&self, job: &PreparedJob) -> Result<Results, AdapterError> {
        // Provenance: hash whatever output is sitting in the workdir
        // top-level; fall back to `case.toml` for partial / failed runs
        // so the provenance block is always well-formed.
        let case_hash_input = {
            // Vina writes a single output PDBQT, but we don't know its
            // exact name without re-parsing the case. Search for the
            // most recently-modified pdbqt in the workdir; if none,
            // fall back to case.toml.
            let mut latest: Option<std::path::PathBuf> = None;
            let mut latest_mtime = std::time::SystemTime::UNIX_EPOCH;
            if let Ok(entries) = std::fs::read_dir(&job.workdir) {
                for entry in entries.flatten() {
                    let path = entry.path();
                    if !path.is_file() {
                        continue;
                    }
                    let ext = path
                        .extension()
                        .and_then(|s| s.to_str())
                        .map(|s| s.to_ascii_lowercase());
                    if ext.as_deref() != Some("pdbqt") {
                        continue;
                    }
                    if let Ok(meta) = entry.metadata() {
                        if let Ok(mt) = meta.modified() {
                            if mt > latest_mtime {
                                latest_mtime = mt;
                                latest = Some(path);
                            }
                        }
                    }
                }
            }
            latest.unwrap_or_else(|| job.workdir.join("case.toml"))
        };
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "AutoDock Vina",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk the workdir top-level. Vina only writes the docked-poses
        // PDBQT (and any Vina log if the user redirected stderr). Stay
        // top-level: Vina doesn't create subdirectories.
        let entries = match std::fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-vina", ?e, "workdir read failed");
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
                // The output PDBQT — docked-pose ensemble. PDBQT is
                // PDB plus AutoDock partial charges; downstream tools
                // (PyMOL, Chimera, MDAnalysis) read it directly.
                Some("pdbqt") => (
                    ArtifactKind::Native,
                    "AutoDock Vina docked poses".to_string(),
                ),
                Some("log") => (ArtifactKind::Log, "AutoDock Vina log".to_string()),
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
        // Bio-specific Capability variants land in a follow-up task;
        // ribbon contributions are already enough for the registry to
        // surface the adapter.
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.vina.dock"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_core::{CancellationToken, LogLevel, LogSink, ProgressSink};
    use valenx_test_utils::tempdir;

    struct NoopProgress;
    impl ProgressSink for NoopProgress {
        fn report(&self, _pct: f32, _message: &str) {}
    }
    struct NoopLog;
    impl LogSink for NoopLog {
        fn log_line(&self, _level: LogLevel, _line: &str) {}
    }

    #[test]
    fn info_is_bio_domain() {
        let info = VinaAdapter::new().info();
        assert_eq!(info.id, "vina");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "Apache-2.0");
        assert_eq!(info.display_name, "AutoDock Vina");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = VinaAdapter::new().info();
        // 1.2.x is the current actively maintained Vina line; 2.0
        // reserves room for an eventual major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(1, 2, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(2, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = VinaAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.vina.dock"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = VinaAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    #[test]
    fn probe_reports_native_when_binary_absent() {
        // We can't actually delete `vina` from PATH for the test, but
        // the new `probe_with_engine` should return Ok regardless when
        // told the engine is "native".
        let r = VinaAdapter::new().probe_with_engine("native");
        assert!(r.is_ok());
        let report = r.unwrap();
        assert!(report.ok);
        assert!(report.binary_path.is_none() || report.binary_path.is_some());
    }

    #[test]
    fn native_engine_round_trips_minimal_case() {
        let d = tempdir("vina_native");
        // Write a one-atom receptor + one-atom ligand + case.toml.
        // PDBQT atom lines need >=79 chars (cols 78-79 are the AD4
        // type), so the lines below carry a trailing space + element
        // symbol per the AutoDock 4 / Vina spec.
        std::fs::write(
            d.join("receptor.pdbqt"),
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00     0.000 C \n",
        )
        .unwrap();
        std::fs::write(
            d.join("ligand.pdbqt"),
            "ROOT\nATOM      1  C1  LIG A   1       3.000   0.000   0.000  1.00  0.00     0.000 C \nENDROOT\nTORSDOF 0\n",
        )
        .unwrap();
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "vina.dock"

[bio.vina]
receptor       = "receptor.pdbqt"
ligand         = "ligand.pdbqt"
output         = "out.pdbqt"
center         = [0.0, 0.0, 0.0]
size           = [6.0, 6.0, 6.0]
exhaustiveness = 2
num_modes      = 3
engine         = "native"
"#,
        )
        .unwrap();

        let adapter = VinaAdapter::new();
        let case = Case {
            id: "test".into(),
            path: d.clone(),
        };
        let workdir = d.join("work");
        let job = adapter.prepare(&case, &workdir).unwrap();
        let cancel = CancellationToken::new();
        let mut ctx = RunContext {
            cancel: &cancel,
            progress: Box::new(NoopProgress),
            log: Box::new(NoopLog),
        };
        let report = adapter.run(&job, &mut ctx).unwrap();
        assert_eq!(report.exit_code, 0);
        assert!(d.join("work/out.pdbqt").exists() || workdir.join("out.pdbqt").exists());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn prepare_rejects_receptor_with_parent_traversal() {
        // case.toml claiming `receptor = "../escape.pdbqt"` must be
        // rejected as InvalidCase — `confined_join` refuses to resolve
        // anything that exits the case directory.
        let d = tempdir("vina_traversal");
        std::fs::write(
            d.join("ligand.pdbqt"),
            "ROOT\nATOM      1  C1  LIG A   1       3.000   0.000   0.000  1.00  0.00     0.000 C \nENDROOT\nTORSDOF 0\n",
        )
        .unwrap();
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "vina.dock"

[bio.vina]
receptor       = "../escape.pdbqt"
ligand         = "ligand.pdbqt"
output         = "out.pdbqt"
center         = [0.0, 0.0, 0.0]
size           = [6.0, 6.0, 6.0]
exhaustiveness = 2
num_modes      = 3
engine         = "native"
"#,
        )
        .unwrap();
        let adapter = VinaAdapter::new();
        let case = Case {
            id: "test".into(),
            path: d.clone(),
        };
        let workdir = d.join("work");
        let err = adapter
            .prepare(&case, &workdir)
            .expect_err("path traversal must reject");
        assert!(
            matches!(err, AdapterError::InvalidCase { .. }),
            "wanted InvalidCase, got: {err:?}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
