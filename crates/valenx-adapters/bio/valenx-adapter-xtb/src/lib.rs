//! # valenx-adapter-xtb
//!
//! Adapter for [xTB](https://github.com/grimme-lab/xtb) — Stefan
//! Grimme's extended tight-binding semiempirical quantum chemistry
//! package. xTB exposes the GFN0 / GFN1 / GFN2 family of
//! parameterisations: orders-of-magnitude faster than ab initio QM at
//! a small accuracy cost, fast enough to drive routine geometry
//! optimisation, conformer search, frequencies, and MD on
//! drug-sized molecules.
//!
//! **Phase 25 — subprocess wrapper around the `xtb` binary.** The
//! user supplies a `.xyz` geometry, picks a run mode (`single-point`,
//! `opt`, `ohess`, `hess`, `md`), and optionally specifies charge,
//! unpaired electrons, GFN parameter set, and an ALPB implicit
//! solvent via `[bio.xtb]` in `case.toml`. `prepare()` resolves the
//! geometry, composes the `xtb <input> --gfn <N> --chrg <q> --uhf <s>
//! [--<mode>] [--alpb <solvent>] [extras...]` invocation, and stages
//! everything in the workdir.
//!
//! xTB writes its run output to **stdout** (energies, gradients,
//! convergence chatter — line-oriented) and drops a stack of
//! satellite files in cwd: `xtbopt.xyz` (the optimised geometry),
//! `xtbopt.log` (the per-step optimisation trajectory), `gradient`
//! (final gradient in Turbomole format), `hessian` (Hessian matrix),
//! and similar. Like MAFFT, we use a custom `run()` that hands the
//! output file's FD to the child as stdout so the line-oriented
//! output lands straight on disk.
//!
//! On `collect()` we surface the canonical artefacts: `xtbopt.xyz`
//! as `Native` (the optimised geometry), `xtbopt.log` / `xtb.log` as
//! `Log`, and `gradient` / `hessian` as `Native`.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::fs;
use std::io::BufReader;
use std::path::Path;
use std::process::{Command, Stdio};
use std::sync::mpsc;
use std::thread;
use std::time::{Duration, Instant};

use semver::Version;

use valenx_core::{
    adapter::LogLevel,
    adapter_helpers::{detect_tool_version_semver, find_on_path, live_provenance},
    error::RunPhase,
    io_caps::read_capped_lines_bounded,
    subprocess::{KillOnDropChild, MAX_LINE_BYTES, SUBPROCESS_CHANNEL_CAPACITY},
    Adapter, AdapterError, AdapterInfo, Capabilities, Case, LicenseMode, Physics, PreparedJob,
    ProbeReport, RunContext, RunReport, VersionRange,
};
use valenx_fields::{
    artifact::{Artifact, ArtifactKind},
    Results,
};

use crate::case_input::XtbInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(XtbAdapter::new())
}

pub struct XtbAdapter;

impl XtbAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for XtbAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "xtb";
/// xTB's binary candidate. Conda / source / Bioconda all install
/// under the canonical `xtb` name.
const BINARIES: &[&str] = &["xtb"];

/// The stdout-redirect target. xTB has no `-o` flag — its run report
/// goes to stdout. Pinned so prepare(), run(), and collect() agree.
const OUT_LOG: &str = "xtb.log";

impl Adapter for XtbAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "xTB",
            // xtb 6.6+ is the floor we test against (current: 6.7+).
            // The 6.x line has been the stable series since 2020;
            // upper bound 7.0 reserves room for an eventual major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(6, 6, 0),
                max_exclusive: Version::new(7, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "LGPL-3.0",
            docs_url: "https://xtb-docs.readthedocs.io/",
            homepage_url: "https://github.com/grimme-lab/xtb",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `xtb --version` prints something like
                // "xtb version 6.6.1" on stdout.
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
                hint: "xTB 6.6+ required; install via \
                       `conda install -c conda-forge xtb` or download a \
                       static binary from https://github.com/grimme-lab/xtb/releases"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = XtbInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve the input .xyz against the case directory if
        // relative. Same convention as every other Phase 17/18 bio
        // adapter — `input = "molecule.xyz"` next to `case.toml`.
        let source_input = if input.input.is_absolute() {
            input.input.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.input)?
        };
        if !source_input.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.xtb].input `{}` not found (resolved {})",
                    input.input.display(),
                    source_input.display()
                ),
            });
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "xTB 6.6+ required; install via \
                       `conda install -c conda-forge xtb` or download a \
                       static binary from https://github.com/grimme-lab/xtb/releases"
                .into(),
        })?;

        // Compose `xtb <input> --gfn <gfn> --chrg <charge> --uhf <uhf>
        //              [--<mode> if not single-point]
        //              [--alpb <solvent> if Some]
        //              [extras...]`.
        //
        // `single-point` is xtb's default run type so it gets no
        // flag; every other mode maps to `--<mode>`. Charge,
        // multiplicity, and the GFN parameter set are always
        // emitted so the invocation is unambiguous regardless of
        // whether xtb's own defaults match ours.
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            source_input.into_os_string(),
            OsString::from("--gfn"),
            OsString::from(input.gfn.to_string()),
            OsString::from("--chrg"),
            OsString::from(input.charge.to_string()),
            OsString::from("--uhf"),
            OsString::from(input.uhf.to_string()),
        ];
        if input.mode != "single-point" {
            native_command.push(OsString::from(format!("--{}", input.mode)));
        }
        if let Some(solvent) = &input.solvent {
            native_command.push(OsString::from("--alpb"));
            native_command.push(OsString::from(solvent));
        }
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Single-point on a small molecule finishes in seconds;
            // ohess on a drug-sized system runs for a few minutes;
            // long MD trajectories can run for hours. 4 hours is a
            // generous default that covers the typical long tail.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        if job.native_command.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "PreparedJob.native_command is empty — prepare() should \
                 have populated it"
            )));
        }

        // Open `xtb.log` for write — xtb writes its run report to
        // stdout, so we hand the file's FD to the child as stdout.
        // Truncates any prior content from a previous run.
        let out_path = job.workdir.join(OUT_LOG);
        let out_file = std::fs::File::create(&out_path).map_err(|e| {
            AdapterError::Other(anyhow::anyhow!(
                "open {} for write: {e}",
                out_path.display()
            ))
        })?;

        let program = &job.native_command[0];
        let args: Vec<&OsString> = job.native_command.iter().skip(1).collect();

        ctx.report_progress(0.0, "starting xTB");
        ctx.log(
            LogLevel::Info,
            &format!(
                "spawning {} with {} arg(s) in {}",
                program.to_string_lossy(),
                args.len(),
                job.workdir.display()
            ),
        );

        let mut cmd = Command::new(program);
        for a in &args {
            cmd.arg(a);
        }
        for (k, v) in &job.environment {
            cmd.env(k, v);
        }
        cmd.current_dir(&job.workdir)
            .stdin(Stdio::null())
            .stdout(Stdio::from(out_file))
            .stderr(Stdio::piped());

        let raw_child = cmd.spawn().map_err(|e| AdapterError::Run {
            exit_code: -1,
            stderr: format!("failed to spawn {}: {e}", program.to_string_lossy()),
            phase: RunPhase::Startup,
        })?;
        // Round-24 H2: KillOnDropChild guard.
        let mut kill_guard = KillOnDropChild::new(raw_child, true);

        // stdout is going to xtb.log; only stderr needs a reader.
        let stderr = kill_guard
            .inner_mut()
            .stderr
            .take()
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("child stderr not captured")))?;
        // Round-24 H2: bounded sync_channel + capped lines.
        let (tx, rx) = mpsc::sync_channel::<String>(SUBPROCESS_CHANNEL_CAPACITY);
        let se_thread = thread::spawn(move || {
            let reader = BufReader::new(stderr);
            for line in read_capped_lines_bounded(reader, MAX_LINE_BYTES) {
                let bytes = match line {
                    Ok(b) => b,
                    Err(_) => break,
                };
                let mut s = String::from_utf8_lossy(&bytes).into_owned();
                if s.ends_with('\n') {
                    s.pop();
                    if s.ends_with('\r') {
                        s.pop();
                    }
                }
                if tx.send(s).is_err() {
                    break;
                }
            }
        });

        let start = Instant::now();
        let mut warnings: Vec<String> = Vec::new();
        let mut stderr_tail: Vec<String> = Vec::new();
        const STDERR_TAIL_MAX: usize = 64;

        loop {
            if ctx.check_cancel().is_err() {
                let _ = kill_guard.inner_mut().kill();
                let _ = kill_guard.inner_mut().wait();
                return Err(AdapterError::Cancelled);
            }
            match rx.recv_timeout(Duration::from_millis(200)) {
                Ok(line) => {
                    process_stderr(&line, ctx, &mut warnings);
                    if stderr_tail.len() >= STDERR_TAIL_MAX {
                        stderr_tail.remove(0);
                    }
                    stderr_tail.push(line);
                }
                Err(mpsc::RecvTimeoutError::Timeout) => {
                    if let Some(status) = kill_guard
                        .inner_mut()
                        .try_wait()
                        .map_err(AdapterError::Io)?
                    {
                        for line in rx.try_iter() {
                            process_stderr(&line, ctx, &mut warnings);
                            if stderr_tail.len() >= STDERR_TAIL_MAX {
                                stderr_tail.remove(0);
                            }
                            stderr_tail.push(line);
                        }
                        let _ = se_thread.join();
                        return finalize(status, start.elapsed(), warnings, stderr_tail);
                    }
                }
                Err(mpsc::RecvTimeoutError::Disconnected) => {
                    let status = kill_guard.inner_mut().wait().map_err(AdapterError::Io)?;
                    let _ = se_thread.join();
                    return finalize(status, start.elapsed(), warnings, stderr_tail);
                }
            }
        }
    }

    fn collect(&self, job: &PreparedJob) -> Result<Results, AdapterError> {
        // Provenance: hash xtbopt.xyz if the optimisation produced
        // one, else xtb.log if the run produced one, else case.toml
        // (so the provenance block is well-formed even on partial /
        // failed runs).
        let case_hash_input = {
            let opt_xyz = job.workdir.join("xtbopt.xyz");
            let log = job.workdir.join(OUT_LOG);
            if opt_xyz.is_file() {
                opt_xyz
            } else if log.is_file() {
                log
            } else {
                job.workdir.join("case.toml")
            }
        };
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "xTB",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk the workdir top-level. xTB writes a stack of fixed-name
        // satellite files: `xtbopt.xyz` for optimised geometries,
        // `xtbopt.log` for the per-step trajectory, our redirected
        // `xtb.log` for the run report, and `gradient` / `hessian`
        // (Turbomole-format text files, no extension) for derivatives.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-xtb", ?e, "workdir read failed");
                return Ok(results);
            }
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if !path.is_file() {
                continue;
            }
            let name = path
                .file_name()
                .and_then(|s| s.to_str())
                .map(|s| s.to_ascii_lowercase());
            // Match by full filename first — xtb's signature outputs
            // are pinned filenames, not extensions. Fall through to
            // extension-only matches for the run log.
            let (kind, label) = match name.as_deref() {
                Some("xtbopt.xyz") => (ArtifactKind::Native, "xTB optimised geometry".to_string()),
                Some("xtbopt.log") => {
                    (ArtifactKind::Log, "xTB optimisation trajectory".to_string())
                }
                Some("xtb.log") => (ArtifactKind::Log, "xTB run log".to_string()),
                Some("gradient") => (
                    ArtifactKind::Native,
                    "xTB gradient (Turbomole format)".to_string(),
                ),
                Some("hessian") => (
                    ArtifactKind::Native,
                    "xTB Hessian (Turbomole format)".to_string(),
                ),
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
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.xtb.compute"],
        }
    }
}

/// Mirror of MAFFT's stderr handler — log every line and lift xTB's
/// progress markers to coarse UI ticks.
fn process_stderr(line: &str, ctx: &mut RunContext<'_>, warnings: &mut Vec<String>) {
    ctx.log(LogLevel::Info, line);
    // xTB's stderr is sparse (most output goes to stdout). The
    // visible markers are version banners at startup and any
    // failure messages.
    if line.contains("normal termination") {
        ctx.report_progress(95.0, line);
    } else if line.contains("xtb version") {
        ctx.report_progress(5.0, line);
    }
    if line.contains("ERROR") || line.contains("Error:") || line.contains("FATAL") {
        warnings.push(line.trim().to_string());
    }
}

/// Bridge a child's exit status into the adapter's RunReport / Run
/// error, mirroring MAFFT's `finalize`.
fn finalize(
    status: std::process::ExitStatus,
    wall_time: Duration,
    warnings: Vec<String>,
    stderr_tail: Vec<String>,
) -> Result<RunReport, AdapterError> {
    let exit_code = status.code().unwrap_or(-1);
    if !status.success() {
        let stderr = if stderr_tail.is_empty() {
            format!("xTB exited {exit_code} with no stderr output")
        } else {
            stderr_tail.join("\n")
        };
        return Err(AdapterError::Run {
            exit_code,
            stderr,
            phase: RunPhase::Solve,
        });
    }
    Ok(RunReport {
        exit_code,
        wall_time,
        converged: Some(true),
        residual_history: Vec::new(),
        warnings,
        final_phase: Some(RunPhase::Shutdown),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = XtbAdapter::new().info();
        assert_eq!(info.id, "xtb");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "LGPL-3.0");
        assert_eq!(info.display_name, "xTB");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = XtbAdapter::new().info();
        assert_eq!(info.version_range.min_inclusive, Version::new(6, 6, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(7, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = XtbAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.xtb.compute"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = XtbAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
