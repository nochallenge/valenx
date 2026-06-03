//! # valenx-adapter-cas-offinder
//!
//! Adapter for [Cas-OFFinder](http://www.rgenome.net/cas-offinder/) —
//! the Bae / Park / Kim group's CRISPR off-target searching tool from
//! Hanyang / Seoul National University. Cas-OFFinder is a fast,
//! OpenCL-accelerated scanner: given a list of guide sequences + PAM
//! patterns + mismatch budget in a plain-text input file, it walks a
//! reference genome and reports every position whose sequence matches
//! one of the guides within the configured Hamming distance. It's the
//! workhorse off-target scanner sitting under most CRISPR design web
//! services (CRISPOR, CRISPRdirect, …) and pipelines.
//!
//! **Phase 35 — subprocess wrapper around the `cas-offinder` binary.**
//! The CLI is fixed-shape:
//!
//! ```sh
//! cas-offinder <input> {C|G|A} <output> [extras...]
//! ```
//!
//! `<input>` is a 3+-line text file with the reference genome path,
//! the PAM pattern, and one guide-sequence row per query. The
//! middle positional argument selects the OpenCL device class — `C`
//! (CPU), `G` (GPU), or `A` (auto-pick fastest at runtime).
//!
//! `prepare()` resolves both paths against the case directory (when
//! relative) and composes the invocation. `run()` streams the run via
//! the shared subprocess runner — Cas-OFFinder is mostly silent on
//! stdout / stderr, so the standard line handler picks up `Loading` /
//! `Finished` markers for free.
//!
//! On `collect()` we report the configured `output` file as a single
//! `Tabular` artifact labeled `"Cas-OFFinder off-target hits"`.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::fs;
use std::path::{Path, PathBuf};
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

use crate::case_input::CasOffinderInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(CasOffinderAdapter::new())
}

pub struct CasOffinderAdapter;

impl CasOffinderAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for CasOffinderAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "cas-offinder";
/// Cas-OFFinder's binary candidates. Conda / source / Bioconda /
/// Homebrew all install under the canonical `cas-offinder` name.
const BINARIES: &[&str] = &["cas-offinder"];

impl Adapter for CasOffinderAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "Cas-OFFinder",
            // Cas-OFFinder's tagged release line is 2.x; the modern
            // OpenCL device-selection CLI stabilised at 2.4. Upper
            // bound 3.0 reserves room for an eventual major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(2, 4, 0),
                max_exclusive: Version::new(3, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "BSD-3-Clause",
            docs_url: "http://www.rgenome.net/cas-offinder/",
            homepage_url: "http://www.rgenome.net/cas-offinder/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // Cas-OFFinder prints a banner including the version
                // when invoked without args (or with `--version` on
                // newer builds). The combined scanner reads stderr
                // and stdout so a bare-name invocation works for
                // older builds that ignore unknown flags.
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
                hint: "Cas-OFFinder 2.4+ required; install via \
                       `conda install -c bioconda cas-offinder` or build \
                       from source at https://github.com/snugel/cas-offinder"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = CasOffinderInput::from_case_dir(&case.path)?;

        // Round-10 H3: `output` is `PathBuf` and pre-fix flowed into
        // `workdir.join(&input.output)`. Validate as a basename
        // before the join so `output = "../etc/passwd"` is rejected.
        if let Some(s) = input.output.to_str() {
            valenx_core::adapter_helpers::validate_output_basename(
                s,
                "[bio.cas_offinder].output",
            )
            .map_err(|e| AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!("{e}"),
            })?;
        } else {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: "[bio.cas_offinder].output: non-UTF-8 path rejected"
                    .into(),
            });
        }

        fs::create_dir_all(workdir)?;

        // Stage `case.toml` into the workdir so collect() can recover
        // the configured `output_basename` for prefix-filtering output
        // artifacts. Without this stage, the basename filter silently
        // degrades to "match everything".
        let staged_case_toml = workdir.join("case.toml");
        let source_case_toml = case.path.join("case.toml");
        if source_case_toml.is_file() {
            fs::copy(&source_case_toml, &staged_case_toml)
                .map_err(|e| AdapterError::Other(anyhow::anyhow!("stage case.toml: {e}")))?;
        }

        // Resolve the input file against the case directory if
        // relative — same convention as every other Phase 18 binary
        // adapter.
        let source_input = if input.input.is_absolute() {
            input.input.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(
            &case.path,
            &input.input,
        )?
        };
        if !source_input.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.cas_offinder].input `{}` not found (resolved {})",
                    input.input.display(),
                    source_input.display()
                ),
            });
        }

        // The output path is *generated*, not consumed. Resolve
        // relative paths into the workdir so the artifact lands next
        // to whatever else Cas-OFFinder writes (logs, tmp index
        // files).
        let output_path: PathBuf = if input.output.is_absolute() {
            input.output.clone()
        } else {
            workdir.join(&input.output)
        };

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "Cas-OFFinder 2.4+ required; install via \
                       `conda install -c bioconda cas-offinder` or build \
                       from source at https://github.com/snugel/cas-offinder"
                .into(),
        })?;

        // Compose `cas-offinder <input> <backend> <output> [extras...]`.
        // Cas-OFFinder's CLI is purely positional — no `-i` / `-o`
        // flags, the order is fixed.
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            source_input.into_os_string(),
            OsString::from(&input.backend),
            output_path.into_os_string(),
        ];
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Cas-OFFinder on a small genome with a handful of guides
            // finishes in seconds; whole-genome scans of tens of
            // thousands of guides routinely run for an hour or more.
            // 4 hours covers the long tail.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting Cas-OFFinder", |line| {
            let mut hint = subprocess::Hint::default();
            // Cas-OFFinder is mostly silent; the few markers it does
            // emit ("Loading...", "Finished") map cleanly to
            // start-up / completion ticks.
            if line.contains("Finished") || line.contains("finished") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("Loading") {
                hint.progress = Some((10.0, line.to_string()));
            } else if line.contains("Error") || line.contains("error") {
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
        // Re-derive the output path from the staged case so collect()
        // reports exactly the file Cas-OFFinder wrote (the user picks
        // the filename via `output = "..."`).
        let case_input = CasOffinderInput::from_case_dir(&job.workdir).ok();

        // Provenance: hash the input file when present (the
        // canonical "this case is configured this way" descriptor).
        // Falls back to case.toml on a partial / failed run.
        let case_hash_input = match &case_input {
            Some(ci) => {
                let p = if ci.input.is_absolute() {
                    ci.input.clone()
                } else {
                    job.workdir.join(&ci.input)
                };
                if p.is_file() {
                    p
                } else {
                    job.workdir.join("case.toml")
                }
            }
            None => job.workdir.join("case.toml"),
        };
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "Cas-OFFinder",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        if let Some(ci) = case_input {
            let output_path: PathBuf = if ci.output.is_absolute() {
                ci.output.clone()
            } else {
                job.workdir.join(&ci.output)
            };
            if output_path.is_file() {
                artefacts.push(Artifact {
                    path: output_path,
                    kind: ArtifactKind::Tabular,
                    checksum: None,
                    label: "Cas-OFFinder off-target hits".to_string(),
                });
            }
        }

        artefacts.sort_by(|a, b| a.path.cmp(&b.path));
        results.artifacts = artefacts;
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.cas-offinder.search"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = CasOffinderAdapter::new().info();
        assert_eq!(info.id, "cas-offinder");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "BSD-3-Clause");
        assert_eq!(info.display_name, "Cas-OFFinder");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = CasOffinderAdapter::new().info();
        // Cas-OFFinder 2.4+ is the modern OpenCL-CLI line; upper
        // bound 3.0 reserves room for an eventual major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(2, 4, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(3, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = CasOffinderAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.cas-offinder.search"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = CasOffinderAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    /// Round-10 H3 RED→GREEN: `output` flowed into
    /// `workdir.join(&input.output)` with no validation. Hostile
    /// `output = "../etc/passwd"` is rejected pre-spawn.
    #[test]
    fn prepare_rejects_output_path_traversal() {
        use valenx_test_utils::tempdir;
        let d = tempdir("cas-offinder-output-trav");
        std::fs::write(d.join("input.txt"), b"test").unwrap();
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "cas-offinder.search"

[bio.cas_offinder]
input   = "input.txt"
output  = "../etc/passwd"
backend = "C"
"#,
        )
        .unwrap();
        let case = Case {
            id: "trav".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = CasOffinderAdapter::new()
            .prepare(&case, &workdir)
            .unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("[bio.cas_offinder].output"),
            "expected [bio.cas_offinder].output in error, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
