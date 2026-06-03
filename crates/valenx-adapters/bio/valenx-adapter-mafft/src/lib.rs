//! # valenx-adapter-mafft
//!
//! Adapter for [MAFFT](https://mafft.cbrc.jp/alignment/software/) —
//! a multiple-sequence-alignment package built around fast Fourier
//! transforms over reduced amino-acid alphabets, plus several
//! accuracy-tuned iterative-refinement strategies (L-INS-i, G-INS-i,
//! E-INS-i). MAFFT is the de-facto MSA workhorse in modern
//! bioinformatics pipelines — comparable accuracy to MUSCLE on
//! medium-sized inputs but scales to thousands of sequences without
//! collapsing.
//!
//! **Phase 18 — subprocess wrapper around `mafft`.** The user
//! supplies a multi-FASTA via `[bio.mafft]` in `case.toml`.
//! `prepare()` resolves it against the case directory, picks the
//! strategy (`auto` by default), and composes the `mafft` invocation.
//! `run()` spawns MAFFT and **redirects its stdout to
//! `aligned.fa` in the workdir** — MAFFT writes the aligned FASTA
//! to stdout and has no `-o` flag, so the standard subprocess runner
//! (which pipes stdout through a line handler) does not fit cleanly.
//!
//! ## Why a custom run() instead of `subprocess::run`
//!
//! Every other Phase 17/18 bio adapter uses [`valenx_core::subprocess::run`],
//! which captures stdout line-by-line for in-flight progress reporting.
//! That works when stdout is chatty progress output and the
//! "real" outputs land as files on disk (BWA's `out.sam`, biopython's
//! `analyse.py` outputs, etc.).
//!
//! MAFFT inverts that contract: stdout *is* the aligned FASTA — the
//! only run-output the user cares about — and stderr carries the
//! progress chatter. Routing stdout through the line handler would
//! force us to either reconstruct the FASTA from a `Vec<String>` after
//! the fact (lossy on whitespace edge cases) or buffer every byte into
//! memory (a large MSA can exceed the input size by 5-10x). Forking
//! the runner to support stdout-to-file would balloon its surface for
//! one caller. Cleanest path: invoke `Command::new("mafft")` directly
//! from `run()` and redirect stdout to a file via `File::create()` +
//! `Stdio::from(file)`.
//!
//! On `collect()` we surface the canonical `aligned.fa` plus any
//! `.log` the user redirected stderr to.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;
pub mod native;

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

use crate::case_input::MafftInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(MafftAdapter::new())
}

pub struct MafftAdapter;

impl MafftAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for MafftAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "mafft";
/// MAFFT's binary candidates. The canonical name is `mafft` —
/// distros, Bioconda, and Homebrew all install under that exact
/// name; the strategy-specific aliases (`linsi`, `ginsi`, etc.) are
/// thin shell wrappers around it.
const BINARIES: &[&str] = &["mafft"];

/// The aligned-FASTA filename written from MAFFT's stdout. Pinned so
/// `prepare()` (which records the command), `run()` (which redirects
/// stdout), and `collect()` (which labels the artifact) all agree.
const OUT_FA: &str = "aligned.fa";

impl Adapter for MafftAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "MAFFT",
            // MAFFT's 7.x line has been the stable series since 2013
            // and ships frequent point releases (current: 7.520+).
            // Floor at 7.500 covers every reasonably modern install;
            // upper bound 8.0 reserves room for the next major.
            version_range: VersionRange {
                min_inclusive: Version::new(7, 500, 0),
                max_exclusive: Version::new(8, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "BSD-3-Clause",
            docs_url: "https://mafft.cbrc.jp/alignment/software/manual/manual.html",
            homepage_url: "https://mafft.cbrc.jp/alignment/software/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `mafft --version` prints something like
                // "v7.520 (2023/Mar/13)" on stderr; our combined
                // stdout+stderr scanner picks it up either way.
                let found_version = detect_tool_version_semver(&binary_path, &["--version", ""]);
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    warnings: Vec::new(),
                    required_env: Vec::new(),
                })
            }
            // Native Rust fallback via valenx-align progressive+iterative MSA.
            None => Ok(ProbeReport {
                ok: true,
                found_version: None,
                binary_path: None,
                warnings: vec![
                    "mafft binary not found; using native Rust progressive+iterative MSA \
                     (valenx-align). Install MAFFT 7.500+ via apt/brew/conda for the full \
                     MAFFT alignment strategies (L-INS-i, G-INS-i, etc.)."
                        .to_string(),
                ],
                required_env: Vec::new(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = MafftInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve the input FASTA against the case directory if
        // relative. Same convention as every other Phase 17/18 bio
        // adapter — `input = "seqs.fa"` next to `case.toml`.
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
                    "[bio.mafft].input `{}` not found (resolved {})",
                    input.input.display(),
                    source_input.display()
                ),
            });
        }

        // Write native_params.toml for the native path.
        let native_params = native::NativeMsaParams {
            input_path: source_input
                .to_str()
                .ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "input path is not valid UTF-8: {}",
                        source_input.display()
                    ))
                })?
                .to_string(),
            output_name: OUT_FA.to_string(),
            refine: true,
            max_iterations: 8,
        };
        native::write_params(workdir, &native_params)?;

        let native_command: Vec<OsString> = match find_on_path(BINARIES) {
            Some(binary_path) => {
                let mut cmd: Vec<OsString> = vec![binary_path.into_os_string()];
                if input.strategy == "auto" {
                    cmd.push(OsString::from("--auto"));
                } else {
                    cmd.push(OsString::from(format!("--{}", input.strategy)));
                }
                cmd.push(OsString::from("--thread"));
                cmd.push(OsString::from(input.threads.to_string()));
                cmd.push(source_input.into_os_string());
                for arg in &input.extra_args {
                    cmd.push(OsString::from(arg));
                }
                cmd
            }
            None => vec![OsString::from(native::NATIVE_SENTINEL)],
        };

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
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

        // Native Rust path: progressive + iterative MSA, no subprocess.
        if job.native_command[0] == native::NATIVE_SENTINEL {
            return native::run_native(&job.workdir, ctx);
        }

        // Build the alignment-output sink. MAFFT writes its FASTA to
        // stdout, so we open `aligned.fa` for write and hand its FD
        // to the child as stdout. Any prior content from a previous
        // run gets truncated.
        let out_path = job.workdir.join(OUT_FA);
        let out_file = std::fs::File::create(&out_path).map_err(|e| {
            AdapterError::Other(anyhow::anyhow!(
                "open {} for write: {e}",
                out_path.display()
            ))
        })?;

        let program = &job.native_command[0];
        let args: Vec<&OsString> = job.native_command.iter().skip(1).collect();

        ctx.report_progress(0.0, "starting MAFFT");
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
        // Round-24 H2: wrap in KillOnDropChild so an early return
        // (cancel, drain-thread panic, IO error mid-loop) always
        // reaps the child. Sister to subprocess::run's guard.
        let mut kill_guard = KillOnDropChild::new(raw_child, true);

        // stdout is going to the file; only stderr needs a reader.
        let stderr = kill_guard
            .inner_mut()
            .stderr
            .take()
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("child stderr not captured")))?;
        // Round-24 H2: bounded sync_channel back-pressures a runaway
        // stderr producer (pre-fix `mpsc::channel` accepted unlimited
        // pending items, an OOM vector) and `read_capped_lines_bounded`
        // caps each line at MAX_LINE_BYTES so a 4 GiB no-newline stderr
        // stream can't OOM the drain thread before the cap fires.
        let (tx, rx) = mpsc::sync_channel::<String>(SUBPROCESS_CHANNEL_CAPACITY);
        let se_thread = thread::spawn(move || {
            let reader = BufReader::new(stderr);
            for line in read_capped_lines_bounded(reader, MAX_LINE_BYTES) {
                let bytes = match line {
                    Ok(b) => b,
                    Err(_) => break, // cap or IO error → stop draining
                };
                // Strip trailing \n if present; lossy UTF-8 covers
                // mojibake that some MAFFT builds emit on stderr.
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
                    if let Some(status) = kill_guard.inner_mut().try_wait().map_err(AdapterError::Io)? {
                        // Drain remaining lines so nothing's lost.
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
        // Provenance: hash the staged aligned.fa if present, falling
        // back to case.toml when the run hasn't produced one yet.
        let case_hash_input = {
            let fa = job.workdir.join(OUT_FA);
            if fa.is_file() {
                fa
            } else {
                job.workdir.join("case.toml")
            }
        };
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "MAFFT",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk the workdir top-level. MAFFT produces `aligned.fa`
        // (our stdout-redirect target); a `.log` file may appear if
        // future cases configure stderr redirection.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-mafft", ?e, "workdir read failed");
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
                // The aligned-FASTA we redirected stdout into. `fa`
                // is the canonical extension; some pipelines use
                // `fasta` / `aln` instead.
                Some("fa") | Some("fasta") | Some("aln") => {
                    (ArtifactKind::Native, "MAFFT alignment (FASTA)".to_string())
                }
                Some("log") => (ArtifactKind::Log, "MAFFT log".to_string()),
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
        // ribbon contributions are already enough for the registry
        // to surface the adapter.
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.mafft.msa"],
        }
    }
}

/// Mirror of the line-by-line stderr handling that
/// `subprocess::run` does for stdout. Logs the line and lifts
/// MAFFT's progress markers ("Strategy: <name>", "Pre-aligning...")
/// to coarse UI ticks.
fn process_stderr(line: &str, ctx: &mut RunContext<'_>, warnings: &mut Vec<String>) {
    ctx.log(LogLevel::Info, line);
    // MAFFT's progress chatter on stderr: typically a "Strategy:"
    // banner at startup, "Making a distance matrix .." mid-run,
    // and a "Done." at the end. Lift them to coarse UI ticks.
    if line.contains("Strategy:") {
        ctx.report_progress(10.0, line);
    } else if line.contains("Making a distance matrix")
        || line.contains("constructing a UPGMA tree")
        || line.contains("Pre-aligning")
    {
        ctx.report_progress(40.0, line);
    } else if line.contains("Iterative refinement") || line.contains("Aligning") {
        ctx.report_progress(70.0, line);
    } else if line.starts_with("Done.") {
        ctx.report_progress(95.0, line);
    }
    if line.contains("ERROR") || line.contains("Error:") {
        warnings.push(line.trim().to_string());
    }
}

/// Bridge a child's exit status into the adapter's RunReport / Run
/// error, mirroring `subprocess::finalize` so MAFFT's failure mode
/// matches every other Subprocess-mode adapter.
fn finalize(
    status: std::process::ExitStatus,
    wall_time: Duration,
    warnings: Vec<String>,
    stderr_tail: Vec<String>,
) -> Result<RunReport, AdapterError> {
    let exit_code = status.code().unwrap_or(-1);
    if !status.success() {
        let stderr = if stderr_tail.is_empty() {
            format!("MAFFT exited {exit_code} with no stderr output")
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
        let info = MafftAdapter::new().info();
        assert_eq!(info.id, "mafft");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "BSD-3-Clause");
        assert_eq!(info.display_name, "MAFFT");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = MafftAdapter::new().info();
        // 7.500 is the floor we test against; 8.0 reserves room
        // for an eventual major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(7, 500, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(8, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = MafftAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.mafft.msa"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = MafftAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    #[test]
    fn probe_always_succeeds() {
        let report = MafftAdapter::new().probe().unwrap();
        assert!(report.ok, "probe() must return ok=true in all modes");
    }

    #[test]
    fn native_sentinel_is_stable() {
        assert_eq!(native::NATIVE_SENTINEL, "valenx:native:msa");
    }
}
