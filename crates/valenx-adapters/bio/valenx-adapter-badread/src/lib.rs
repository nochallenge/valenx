//! # valenx-adapter-badread
//!
//! Adapter for [Badread](https://github.com/rrwick/Badread) — Ryan
//! Wick's long-read simulator with realistic Nanopore (and PacBio
//! CLR) error profiles. Badread's per-platform error models are
//! calibrated against actual sequencer output: random / chimeric /
//! adapter / glitch read types, junk-read injection, identity drift,
//! and length distributions that match what users see from a live
//! flowcell. The de-facto choice for stress-testing long-read
//! pipelines under realistic conditions.
//!
//! **Phase 31 — subprocess wrapper around `badread simulate`.** The
//! user supplies a reference FASTA plus the desired quantity and
//! error-model parameters via `[bio.badread]` in `case.toml`.
//! `prepare()` resolves the reference, validates the quantity
//! literal and identity / length distribution, and composes the
//! `badread simulate` invocation. Badread writes its simulated FASTQ
//! to **stdout** (no `-o` flag), so `run()` borrows MAFFT's
//! stdout-redirect-to-file pattern: spawn the child directly,
//! attach stdout to a `File`, stream stderr through the line
//! handler.
//!
//! On `collect()` we surface the user-chosen output path as the
//! single Tabular artifact ("Badread simulated reads").
//!
//! ## Why a custom run() instead of `subprocess::run`
//!
//! Mirror of MAFFT's reasoning: stdout *is* the only run-output the
//! user cares about (the simulated FASTQ), and stderr carries the
//! progress chatter. Routing stdout through the shared
//! [`valenx_core::subprocess::run`] line handler would force us to
//! reconstruct the FASTQ from a `Vec<String>` after the fact (lossy
//! on whitespace edge cases) or buffer every byte into memory (a
//! 5 Gb simulation produces a 5+ Gb FASTQ). Cleanest path: invoke
//! `Command::new("badread")` directly from `run()` and redirect
//! stdout to the user's chosen file via `File::create()` +
//! `Stdio::from(file)`.

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

use crate::case_input::BadreadInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(BadreadAdapter::new())
}

pub struct BadreadAdapter;

impl BadreadAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for BadreadAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "badread";
/// Badread's binary candidate. `badread` is the canonical entrypoint
/// from PyPI / Bioconda installs; it dispatches `simulate`, `error_model`,
/// `qscore_model`, etc. via subcommands.
const BINARIES: &[&str] = &["badread"];

impl Adapter for BadreadAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "Badread",
            // Badread's 0.4.x line is the long-running stable
            // series; a 1.0 cut hasn't happened yet but we
            // reserve room for it. Floor at 0.4.0 covers every
            // reasonable install.
            version_range: VersionRange {
                min_inclusive: Version::new(0, 4, 0),
                max_exclusive: Version::new(1, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-3.0",
            docs_url: "https://github.com/rrwick/Badread",
            homepage_url: "https://github.com/rrwick/Badread",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `badread --version` prints "Badread v0.4.x" on
                // stdout; the detector picks it up.
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
                hint: "Badread 0.4+ required; install via `pip install badread` \
                       or `conda install -c bioconda badread`"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = BadreadInput::from_case_dir(&case.path)?;

        // Round-10 H3: `output` is `PathBuf` from `case.toml` and
        // flowed into `workdir.join(&input.output)` in run() with no
        // validation. A hostile case `output = "../etc/passwd"`
        // wrote the simulated FASTQ outside the workdir. Validate as
        // a basename — Badread writes a single FASTQ file, not a
        // directory tree.
        if let Some(s) = input.output.to_str() {
            valenx_core::adapter_helpers::validate_output_basename(s, "[bio.badread].output")
                .map_err(|e| AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!("{e}"),
                })?;
        } else {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: "[bio.badread].output: non-UTF-8 path rejected".into(),
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

        // Resolve the reference path against the case directory if
        // relative. Same convention as every other Phase 17/18 bio
        // adapter — `reference = "ref.fa"` next to `case.toml`.
        let source_reference = if input.reference.is_absolute() {
            input.reference.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.reference)?
        };
        if !source_reference.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.badread].reference `{}` not found (resolved {})",
                    input.reference.display(),
                    source_reference.display()
                ),
            });
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "Badread 0.4+ required; install via `pip install badread` \
                       or `conda install -c bioconda badread`"
                .into(),
        })?;

        // Compose `badread simulate --reference <ref> --quantity <q>
        // --error_model <m> --identity <id> --length <mean>,<sd>
        // [extras...]`. Badread writes the FASTQ to stdout; our
        // run() redirects stdout to `<workdir>/<output>`.
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from("simulate"),
            OsString::from("--reference"),
            source_reference.into_os_string(),
            OsString::from("--quantity"),
            OsString::from(&input.quantity),
            OsString::from("--error_model"),
            OsString::from(&input.error_model),
            OsString::from("--identity"),
            OsString::from(format!("{}", input.identity_mean)),
            OsString::from("--length"),
            OsString::from(format!("{},{}", input.length_mean, input.length_sd)),
        ];
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Small simulations (~100 Mb) finish in seconds; deep
            // whole-genome runs (~50 Gb) can take an hour. 4 hours
            // is generous enough for the long tail.
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

        // Re-read the case to discover the user-chosen output path.
        // The workdir always contains a fresh case.toml (the
        // executor stages it), so this is the canonical place to
        // look. Falling back to "reads.fq" if the case is missing
        // would silently mask a staging bug — fail loudly instead.
        let input = BadreadInput::from_case_dir(&job.workdir).map_err(|e| {
            AdapterError::Other(anyhow::anyhow!(
                "badread run() could not re-read case.toml in {}: {e}",
                job.workdir.display()
            ))
        })?;
        let out_path = if input.output.is_absolute() {
            input.output.clone()
        } else {
            job.workdir.join(&input.output)
        };
        let out_file = std::fs::File::create(&out_path).map_err(|e| {
            AdapterError::Other(anyhow::anyhow!(
                "open {} for write: {e}",
                out_path.display()
            ))
        })?;

        let program = &job.native_command[0];
        let args: Vec<&OsString> = job.native_command.iter().skip(1).collect();

        ctx.report_progress(0.0, "starting Badread");
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
        // Round-24 H2: KillOnDropChild guard — see mafft for rationale.
        let mut kill_guard = KillOnDropChild::new(raw_child, true);

        // stdout is going to the file; only stderr needs a reader.
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
        // Provenance: hash the case.toml as the canonical input —
        // the FASTQ filename is user-chosen and the FASTQ itself
        // can be huge, so case.toml is the stable choice.
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "Badread",
            "unknown",
            &job.workdir.join("case.toml"),
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);

        // Re-resolve the user-chosen output path against the
        // workdir. Badread writes there because we redirected
        // stdout in run(); reporting it directly (rather than
        // walking the workdir for `*.fq`) keeps the artifact's
        // label tied to the user's case.toml choice.
        let input = match BadreadInput::from_case_dir(&job.workdir) {
            Ok(i) => i,
            Err(e) => {
                tracing::warn!(target: "valenx-badread", ?e, "collect could not re-read case.toml");
                return Ok(results);
            }
        };
        let out_path = if input.output.is_absolute() {
            input.output.clone()
        } else {
            job.workdir.join(&input.output)
        };
        results.artifacts = vec![Artifact {
            path: out_path,
            kind: ArtifactKind::Tabular,
            checksum: None,
            label: "Badread simulated reads".to_string(),
        }];
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        // Bio-specific Capability variants land in a follow-up task;
        // ribbon contributions are already enough for the registry
        // to surface the adapter.
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.badread.simulate"],
        }
    }
}

/// Mirror of the line-by-line stderr handling that
/// `subprocess::run` does for stdout. Logs the line and lifts
/// Badread's progress markers ("Generating reads", "Loading error
/// model", etc.) to coarse UI ticks.
fn process_stderr(line: &str, ctx: &mut RunContext<'_>, warnings: &mut Vec<String>) {
    ctx.log(LogLevel::Info, line);
    // Badread's stderr chatter: a "Loading error model" /
    // "Loading qscore model" banner at startup, "Generating reads"
    // mid-run, and a per-read counter trailing toward completion.
    // Lift the obvious markers to coarse UI ticks.
    if line.contains("Loading error model") || line.contains("Loading qscore model") {
        ctx.report_progress(10.0, line);
    } else if line.contains("Generating reads") {
        ctx.report_progress(40.0, line);
    } else if line.contains("simulated") || line.contains("read length") {
        ctx.report_progress(70.0, line);
    } else if line.starts_with("Done") || line.contains("simulation complete") {
        ctx.report_progress(95.0, line);
    }
    if line.to_ascii_lowercase().contains("error")
        || line.to_ascii_lowercase().contains("traceback")
    {
        warnings.push(line.trim().to_string());
    }
}

/// Bridge a child's exit status into the adapter's RunReport / Run
/// error, mirroring `subprocess::finalize` so Badread's failure mode
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
            format!("Badread exited {exit_code} with no stderr output")
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
        let info = BadreadAdapter::new().info();
        assert_eq!(info.id, "badread");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "GPL-3.0");
        assert_eq!(info.display_name, "Badread");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = BadreadAdapter::new().info();
        // 0.4 is the long-running stable line; 1.0 reserves room
        // for an eventual major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(0, 4, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(1, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = BadreadAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.badread.simulate"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = BadreadAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    /// Round-10 H3 RED→GREEN: `output` is `PathBuf` from case.toml
    /// and pre-fix flowed into `workdir.join(&input.output)`. A
    /// hostile `output = "../etc/passwd"` wrote the simulated FASTQ
    /// outside the workdir. Validate now rejects in `prepare()`.
    #[test]
    fn prepare_rejects_output_path_traversal() {
        use valenx_test_utils::tempdir;
        let d = tempdir("badread-output-trav");
        std::fs::write(d.join("ref.fa"), b">x\nACGT\n").unwrap();
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "badread.simulate"

[bio.badread]
reference   = "ref.fa"
output      = "../etc/passwd"
quantity    = "100M"
error_model = "nanopore2023"
"#,
        )
        .unwrap();
        let case = Case {
            id: "badread-output-trav".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = BadreadAdapter::new().prepare(&case, &workdir).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("[bio.badread].output"),
            "expected [bio.badread].output in error, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
