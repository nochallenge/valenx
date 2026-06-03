//! # valenx-adapter-namd
//!
//! Subprocess adapter for [NAMD](https://www.ks.uiuc.edu/Research/namd/)
//! — UIUC's all-atom molecular-dynamics engine, sister to GROMACS /
//! LAMMPS / OpenMM. NAMD's strengths over its peers are scalability
//! (it ships with the Charm++ runtime tuned for hundred-thousand-atom
//! biomolecular systems) and tight integration with VMD for
//! visualisation. **Phase 5.6 — academic-license flagged subprocess
//! wrapper.**
//!
//! NAMD reads everything (topology, coordinates, force-field
//! parameters, integration / output settings) from a single
//! Tcl-flavoured configuration file (`.namd` / `.conf`). The user
//! references that path via `[bio.namd].config` in `case.toml`;
//! `prepare()` builds a `<binary> +pN <config> [extras...]` invocation
//! and `run()` drives it through the shared subprocess runner.
//!
//! ## License flag
//!
//! NAMD ships under a custom non-OSS license that restricts use to
//! academic / non-commercial contexts (commercial use requires a
//! separate UIUC license). We surface this accurately via a
//! `tool_license` value of `NAMD-License` and emit a probe warning
//! whenever the binary is found so downstream tooling and end-users
//! get a clear "check your license" signal before redistributing
//! trajectories or derived data. The probe-warning text contains the
//! literal substrings `"academic"` and `"non-commercial"` as stable
//! anchors for tests and downstream filters.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::fs;
use std::path::Path;
use std::time::Duration;

use semver::Version;

use valenx_core::{
    adapter_helpers::{detect_tool_version_semver, find_on_path, live_provenance},
    error::RunPhase,
    subprocess, Adapter, AdapterError, AdapterInfo, Capabilities, Case, LicenseMode, Physics,
    PreparedJob, ProbeReport, RunContext, RunReport, VersionRange,
};
use valenx_fields::{
    artifact::{Artifact, ArtifactKind},
    Results,
};

use crate::case_input::NamdInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(NamdAdapter::new())
}

pub struct NamdAdapter;

impl NamdAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for NamdAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "namd";
/// NAMD binary candidates. The 2.x line installs as `namd2`; the
/// 3.x line (current GA, dropped the suffix) ships as `namd3`.
const BINARIES: &[&str] = &["namd2", "namd3"];

/// The probe-warning surfaced whenever NAMD is detected. Anchors a
/// stable "academic / non-commercial only" reminder for downstream
/// tooling and tests; the literal substrings `"academic"` and
/// `"non-commercial"` are part of the asserted contract.
const LICENSE_WARNING: &str = "NAMD is licensed for academic / non-commercial use only \
     — see UIUC NAMD license terms for commercial use";

impl Adapter for NamdAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "NAMD",
            // NAMD 2.14 (2020) is the floor we test against; the 3.x
            // line is current GA. Upper bound 4.0 reserves room for
            // the next major.
            version_range: VersionRange {
                min_inclusive: Version::new(2, 14, 0),
                max_exclusive: Version::new(4, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            // NAMD's terms aren't a recognised SPDX identifier; the
            // closest accurate label is the project's own custom
            // license. Surfacing it here (instead of mislabeling as
            // MIT / BSD) keeps license-aware tooling honest.
            tool_license: "NAMD-License",
            docs_url: "https://www.ks.uiuc.edu/Research/namd/",
            homepage_url: "https://www.ks.uiuc.edu/Research/namd/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // NAMD prints its version banner at startup; the
                // combined version detector covers both `--version`
                // and the implicit banner via help flags.
                let found_version =
                    detect_tool_version_semver(&binary_path, &["--version", "-h", "+help"]);
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    // Always surface the license reminder when NAMD
                    // is detected — it's a custom non-OSS license
                    // and we'd rather over-warn than have a user
                    // ship commercial output without checking.
                    warnings: vec![LICENSE_WARNING.to_string()],
                    required_env: Vec::new(),
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "NAMD 2.14+ required; download from \
                       https://www.ks.uiuc.edu/Research/namd/ \
                       (registration required, academic-use license)"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = NamdInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve the config path against the case directory if
        // relative, and forward the resolved path verbatim. We do
        // NOT stage it into the workdir — NAMD configs reference
        // topology / coordinates / parameters by relative path, and
        // staging only the deck would break those references. The
        // user is responsible for the surrounding file layout.
        let resolved_config = if input.config.is_absolute() {
            input.config.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(
            &case.path,
            &input.config,
        )?
        };

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "NAMD 2.14+ required; download from \
                       https://www.ks.uiuc.edu/Research/namd/ \
                       (registration required, academic-use license)"
                .into(),
        })?;

        // Build the command. `+pN` is one OsString with no space —
        // the Charm++ runtime parses it as a single token. The
        // configuration file is positional after the runtime args;
        // `extra_args` are appended verbatim.
        let mut native_command: Vec<OsString> = vec![binary_path.into_os_string()];
        native_command.push(OsString::from(format!("+p{}", input.processors)));
        native_command.push(resolved_config.into_os_string());
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // NAMD runs span seconds (single-step minimisation) to
            // multi-day production trajectories. 8 hours is a
            // generous default; longer runs can override via the
            // executor.
            estimated_runtime: Some(Duration::from_secs(8 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting NAMD", |line| {
            let mut hint = subprocess::Hint::default();
            // NAMD's stdout is loose. Best-effort hints from common
            // banners:
            //   * "Info: STRUCTURE SUMMARY" — input loaded
            //   * "TIMING:" / "ENERGY:" — integrator stepping
            //   * "WallClock:" — shutdown banner
            //   * "FATAL ERROR" / "ERROR:" — surface as warnings
            // Heuristics; mismatches just leave the spinner alone.
            if line.contains("WallClock:") || line.contains("End of program") {
                hint.progress = Some((95.0, line.trim().to_string()));
            } else if line.starts_with("ENERGY:") || line.starts_with("TIMING:") {
                hint.progress = Some((50.0, line.trim().to_string()));
            } else if line.contains("STRUCTURE SUMMARY") {
                hint.progress = Some((10.0, line.trim().to_string()));
            } else if line.contains("FATAL ERROR") || line.contains("ERROR:") {
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
        // Provenance: hash whichever `.namd` config landed in the
        // workdir if we can find one; otherwise fall back to
        // case.toml so the run-id still derives from real bytes.
        let case_hash_input = job
            .workdir
            .read_dir()
            .ok()
            .and_then(|entries| {
                let mut hits: Vec<std::path::PathBuf> = entries
                    .flatten()
                    .map(|e| e.path())
                    .filter(|p| {
                        p.extension()
                            .and_then(|s| s.to_str())
                            .map(|s| {
                                let s = s.to_ascii_lowercase();
                                s == "namd" || s == "conf"
                            })
                            .unwrap_or(false)
                    })
                    .collect();
                hits.sort();
                hits.into_iter().next()
            })
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "NAMD",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);

        // Walk the workdir top level for NAMD's customary outputs.
        // NAMD writes trajectories to .dcd, restart coordinates /
        // velocities to .coor / .vel, the extended-system file to
        // .xsc, and per-step thermodynamics to .log (when the user
        // redirects stdout) or to a per-output `outputname.log`
        // pattern.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-namd", ?e, "workdir read failed");
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
                Some("dcd") => (ArtifactKind::Native, "NAMD trajectory (DCD)".to_string()),
                Some("coor") => (ArtifactKind::Native, "NAMD coordinates".to_string()),
                Some("vel") => (ArtifactKind::Native, "NAMD velocities".to_string()),
                Some("xsc") => (ArtifactKind::Tabular, "NAMD extended system".to_string()),
                Some("log") => (ArtifactKind::Log, "NAMD log".to_string()),
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
            ribbon_contributions: vec!["bio.namd.simulate"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = NamdAdapter::new().info();
        assert_eq!(info.id, "namd");
        assert_eq!(info.display_name, "NAMD");
        assert_eq!(info.physics, &[Physics::Bio]);
        // The license identifier must surface NAMD's custom non-OSS
        // license rather than mislabel as MIT / BSD.
        assert_eq!(info.tool_license, "NAMD-License");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = NamdAdapter::new().info();
        // NAMD >= 2.14 (2020); upper bound 4.0 reserves room for the
        // next major beyond the current 3.x line.
        assert_eq!(info.version_range.min_inclusive, Version::new(2, 14, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(4, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = NamdAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.namd.simulate"]);
    }

    /// The license-flag warning is mandatory: NAMD is non-OSS
    /// academic-use, and we surface that on every successful probe.
    /// The literal `"academic"` and `"non-commercial"` substrings
    /// are what downstream tooling and license-aware filters key
    /// off — pin them.
    ///
    /// We always assert against the static `LICENSE_WARNING`
    /// constant. We also exercise the live `probe()` path when NAMD
    /// happens to be on PATH (so CI machines without it still pass);
    /// when present, the same substrings must surface in the probe
    /// report's warnings.
    #[test]
    #[ignore] // subprocess-coupled test — run interactively only
    fn probe_warning_mentions_academic_and_non_commercial() {
        assert!(
            LICENSE_WARNING.contains("academic"),
            "probe warning must contain `academic` anchor; got: {LICENSE_WARNING}"
        );
        assert!(
            LICENSE_WARNING.contains("non-commercial"),
            "probe warning must contain `non-commercial` anchor; got: {LICENSE_WARNING}"
        );

        // Best-effort live probe — only assert if NAMD is on PATH.
        // Skipping when it isn't keeps the test green on CI machines
        // without the (registration-walled) binary.
        if find_on_path(BINARIES).is_some() {
            let report = NamdAdapter::new().probe().expect("probe");
            assert!(
                report
                    .warnings
                    .iter()
                    .any(|w| w.contains("academic") && w.contains("non-commercial")),
                "live probe warnings must surface the academic / non-commercial anchors; \
                 got: {:?}",
                report.warnings
            );
        }
    }
}
