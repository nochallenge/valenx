//! # valenx-adapter-linearfold
//!
//! Adapter for [LinearFold](https://github.com/LinearFold/LinearFold)
//! — the Baidu/OSU linear-time RNA secondary-structure folder. The
//! first algorithm to break the cubic-time barrier of classic
//! dynamic-programming folders, making genome-scale single-sequence
//! folding tractable. Sister to LinearDesign in the upstream group's
//! linear-time-RNA portfolio (LinearDesign = co-design; LinearFold =
//! folding-only).
//!
//! **Phase 44.5 — subprocess wrapper around `linearfold` with stdin /
//! stdout redirect.** LinearFold reads the sequence from **stdin** and
//! writes the predicted structure to **stdout**; the CLI has no
//! `-i` / `-o` flags. The adapter mirrors the MAFFT pattern: spawn the
//! child directly, redirect stdin from the user's sequence file via
//! `Stdio::from(File::open(...))`, redirect stdout to
//! `<output_basename>.txt` via `Stdio::from(File::create(...))`, and
//! drain stderr line-by-line for progress chatter.
//!
//! ## Why a custom run() instead of `subprocess::run`
//!
//! The shared `subprocess::run` helper closes the child's stdin
//! (`Stdio::null()`) and pipes stdout through a line handler. LinearFold
//! inverts both contracts: stdin *is* the sequence input, and stdout
//! *is* the structure output the user cares about. Routing either
//! through the line handler would be lossy or nonsensical. Cleanest
//! path: invoke `Command::new("linearfold")` directly from `run()`,
//! redirect stdin from the source sequence file and stdout to the
//! output file, and let stderr carry the chatter.
//!
//! On `collect()` we surface the captured `<output_basename>.txt` as
//! the canonical structure artifact, plus any `*.log` the user
//! redirected stderr to.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::fs::{self, File};
use std::io::BufReader;
use std::path::Path;
use std::process::{Command, Stdio};
use std::sync::mpsc;
use std::thread;
use std::time::{Duration, Instant};

use semver::Version;

use valenx_core::{
    adapter::LogLevel,
    adapter_helpers::{find_on_path, live_provenance},
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

use crate::case_input::LinearFoldInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(LinearFoldAdapter::new())
}

pub struct LinearFoldAdapter;

impl LinearFoldAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for LinearFoldAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "linearfold";
/// LinearFold's binary candidate. The upstream repo ships a single
/// `linearfold` driver; conda-forge / source builds expose the same
/// canonical lowercase name.
const BINARIES: &[&str] = &["linearfold"];
/// Python-interpreter candidates probed only when `linearfold` is
/// missing — surfaces a more useful "you have Python but not the
/// LinearFold repo" hint for the common case where the user has
/// half-installed the tool.
const PYTHON_BINARIES: &[&str] = &["python3", "python"];

/// Sentinel env var used to thread the configured `output_basename`
/// from `prepare()` to `run()` and `collect()`. The custom run()
/// strips this var before spawning so LinearFold never sees it — the
/// env table here is purely an adapter-internal channel. Same scratch-
/// var pattern ViennaRNA uses for its stdout-redirect filename.
const OUTPUT_BASENAME_ENV_VAR: &str = "VALENX_LINEARFOLD_OUTPUT_BASENAME";
/// Sentinel env var carrying the absolute path to the source sequence
/// file. LinearFold reads from stdin; run() opens this path with
/// `File::open(...)` and hands the FD to the child as stdin. Stripped
/// before spawn.
const SOURCE_SEQUENCE_ENV_VAR: &str = "VALENX_LINEARFOLD_SOURCE_SEQUENCE";

impl Adapter for LinearFoldAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "LinearFold",
            // LinearFold 1.0 is the canonical release line tagged in
            // the upstream repo (it tracks the published Bioinformatics
            // paper). Upper bound 2.0 reserves room for the next major
            // bump.
            version_range: VersionRange {
                min_inclusive: Version::new(1, 0, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "Apache-2.0",
            docs_url: "https://github.com/LinearFold/LinearFold",
            homepage_url: "https://github.com/LinearFold/LinearFold",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => Ok(ProbeReport {
                ok: true,
                // Upstream `linearfold` doesn't expose a `--version`
                // banner; skip detection rather than emit a wrong
                // value. The version_range still gates compatibility
                // when the user supplies one explicitly.
                found_version: None,
                binary_path: Some(binary_path),
                warnings: Vec::new(),
                required_env: Vec::new(),
            }),
            None => {
                // If Python is on PATH but the LinearFold driver isn't,
                // surface the half-installed-tool hint as a warning
                // rather than a hard ToolNotInstalled — the user still
                // benefits from validation flowing through even when
                // the binary's missing. Sister pattern to LinearDesign.
                if find_on_path(PYTHON_BINARIES).is_some() {
                    Ok(ProbeReport {
                        ok: false,
                        found_version: None,
                        binary_path: None,
                        warnings: vec!["LinearFold not found on PATH; clone \
                             https://github.com/LinearFold/LinearFold and \
                             add the bin directory to PATH"
                            .into()],
                        required_env: Vec::new(),
                    })
                } else {
                    Err(AdapterError::ToolNotInstalled {
                        name: INFO_ID,
                        hint: "LinearFold 1.0+ required; clone \
                               https://github.com/LinearFold/LinearFold and \
                               add the bin directory to PATH"
                            .into(),
                    })
                }
            }
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = LinearFoldInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.linearfold].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Resolve the sequence file against the case directory if
        // relative. LinearFold reads from stdin; we don't stage the
        // file into the workdir, just validate it exists so the failure
        // is fast and obvious.
        let source_sequence = if input.sequence.is_absolute() {
            input.sequence.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.sequence)?
        };
        if !source_sequence.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.linearfold].sequence `{}` not found (resolved {})",
                    input.sequence.display(),
                    source_sequence.display()
                ),
            });
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "LinearFold 1.0+ required; clone \
                       https://github.com/LinearFold/LinearFold and add \
                       the bin directory to PATH"
                .into(),
        })?;

        // Compose the LinearFold invocation. Model selection is via
        // `-V` (ViennaRNA) or `-C` (CONTRAfold); the beam size follows
        // as a separate token. LinearFold's CLI accepts both
        // `-b <size>` and `--beamsize <size>`; the short form survives
        // quoting better.
        let mut native_command: Vec<OsString> = vec![binary_path.into_os_string()];
        match input.model.as_str() {
            "V" => native_command.push(OsString::from("-V")),
            "C" => native_command.push(OsString::from("-C")),
            // case_input.rs validation rejects anything else; this
            // arm is unreachable but keeps the match exhaustive.
            _ => native_command.push(OsString::from("-C")),
        }
        native_command.push(OsString::from("-b"));
        native_command.push(OsString::from(input.beam_size.to_string()));
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        // Stash the output basename and source sequence path under
        // sentinel env vars so run() can recover them. The custom run()
        // strips both vars before spawning the child so LinearFold
        // never sees them — they're an adapter-internal channel.
        let environment: Vec<(OsString, OsString)> = vec![
            (
                OsString::from(OUTPUT_BASENAME_ENV_VAR),
                OsString::from(&input.output_basename),
            ),
            (
                OsString::from(SOURCE_SEQUENCE_ENV_VAR),
                source_sequence.into_os_string(),
            ),
        ];

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment,
            // A short tRNA finishes in milliseconds; long mRNAs and
            // viral genomes can run for tens of minutes at large beam
            // sizes. 30 minutes is generous for the long tail.
            estimated_runtime: Some(Duration::from_secs(30 * 60)),
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

        // Recover the output basename and source sequence path that
        // prepare() stashed, then strip both vars from the env table
        // so LinearFold doesn't see them.
        let mut output_basename: Option<OsString> = None;
        let mut source_sequence: Option<OsString> = None;
        let mut filtered_env: Vec<(OsString, OsString)> = Vec::with_capacity(job.environment.len());
        for (k, v) in &job.environment {
            if k == OsString::from(OUTPUT_BASENAME_ENV_VAR).as_os_str() {
                output_basename = Some(v.clone());
            } else if k == OsString::from(SOURCE_SEQUENCE_ENV_VAR).as_os_str() {
                source_sequence = Some(v.clone());
            } else {
                filtered_env.push((k.clone(), v.clone()));
            }
        }
        let output_basename = output_basename.ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!(
                "missing {OUTPUT_BASENAME_ENV_VAR} in PreparedJob.environment — \
                 prepare() should have populated it"
            ))
        })?;
        let source_sequence = source_sequence.ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!(
                "missing {SOURCE_SEQUENCE_ENV_VAR} in PreparedJob.environment — \
                 prepare() should have populated it"
            ))
        })?;

        // Build the stdin / stdout sinks. LinearFold reads sequence
        // from stdin and writes structure to stdout; we open the
        // source file for read and the output file for write, then
        // hand the FDs to the child via Stdio::from(file). Mirror of
        // MAFFT's stdout-redirect plus CTFFIND's stdin-redirect, but
        // doing both simultaneously.
        let in_file = File::open(&source_sequence).map_err(|e| {
            AdapterError::Other(anyhow::anyhow!(
                "open {} for read: {e}",
                std::path::Path::new(&source_sequence).display()
            ))
        })?;

        let mut out_filename = output_basename.clone();
        out_filename.push(".txt");
        let out_path = job.workdir.join(std::path::PathBuf::from(&out_filename));
        let out_file = File::create(&out_path).map_err(|e| {
            AdapterError::Other(anyhow::anyhow!(
                "open {} for write: {e}",
                out_path.display()
            ))
        })?;

        let program = &job.native_command[0];
        let args: Vec<&OsString> = job.native_command.iter().skip(1).collect();

        ctx.report_progress(0.0, "starting LinearFold");
        ctx.log(
            LogLevel::Info,
            &format!(
                "spawning {} with {} arg(s) in {}; stdin <- {}; stdout -> {}",
                program.to_string_lossy(),
                args.len(),
                job.workdir.display(),
                std::path::Path::new(&source_sequence).display(),
                out_path.display()
            ),
        );

        let mut cmd = Command::new(program);
        for a in &args {
            cmd.arg(a);
        }
        for (k, v) in &filtered_env {
            cmd.env(k, v);
        }
        cmd.current_dir(&job.workdir)
            .stdin(Stdio::from(in_file))
            .stdout(Stdio::from(out_file))
            .stderr(Stdio::piped());

        let raw_child = cmd.spawn().map_err(|e| AdapterError::Run {
            exit_code: -1,
            stderr: format!("failed to spawn {}: {e}", program.to_string_lossy()),
            phase: RunPhase::Startup,
        })?;
        // Round-24 H2: KillOnDropChild guard.
        let mut kill_guard = KillOnDropChild::new(raw_child, true);

        // stdin is reading from the file; stdout is going to the file;
        // only stderr needs a reader.
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
        // Recover the configured output basename so collect() filters
        // typed outputs to that stem. Same env-channel pattern run()
        // uses; keeps prepare/run/collect in lockstep without re-parsing
        // case.toml.
        let basename = job
            .environment
            .iter()
            .find(|(k, _)| k == OUTPUT_BASENAME_ENV_VAR)
            .and_then(|(_, v)| v.to_str().map(|s| s.to_string()));

        // Provenance: hash the captured structure file when present,
        // falling back to case.toml when the run hasn't produced one
        // yet — keeps the provenance block well-formed for partial /
        // failed runs.
        let case_hash_input = match basename.as_deref() {
            Some(b) => {
                let p = job.workdir.join(format!("{b}.txt"));
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
            "linearfold",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk the workdir top-level. LinearFold's only canonical
        // output is the redirected `<output_basename>.txt`; we also
        // surface any `.log` file (typical when the user redirected
        // stderr) so users get a complete artefact list.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-linearfold", ?e, "workdir read failed");
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
            let stem = path
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("")
                .to_string();
            let stem_matches_basename = match basename.as_deref() {
                Some(b) => stem.starts_with(b),
                None => true,
            };
            match ext.as_deref() {
                Some("txt") => {
                    if !stem_matches_basename {
                        continue;
                    }
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Tabular,
                        checksum: None,
                        label: "LinearFold structure output".to_string(),
                    });
                }
                Some("log") => {
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Log,
                        checksum: None,
                        label: "LinearFold log".to_string(),
                    });
                }
                _ => continue,
            }
        }
        artefacts.sort_by(|a, b| a.path.cmp(&b.path));
        results.artifacts = artefacts;
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.linearfold.fold"],
        }
    }
}

/// Mirror of MAFFT's stderr-line handler. Logs the line and lifts
/// LinearFold's progress markers to coarse UI ticks. LinearFold's
/// stderr is sparse — typically a startup banner and a final
/// "Free Energy" line — so the heuristics are best-effort.
fn process_stderr(line: &str, ctx: &mut RunContext<'_>, warnings: &mut Vec<String>) {
    ctx.log(LogLevel::Info, line);
    if line.contains("LinearFold") || line.contains("linearfold") {
        ctx.report_progress(5.0, line);
    } else if line.contains("Free Energy") || line.contains("free energy") {
        ctx.report_progress(95.0, line);
    } else if line.contains("Folding") || line.contains("folding") {
        ctx.report_progress(50.0, line);
    }
    if line.contains("ERROR") || line.contains("Error:") {
        warnings.push(line.trim().to_string());
    }
}

/// Bridge a child's exit status into RunReport / Run-error, mirroring
/// MAFFT's finalize so LinearFold's failure mode matches every other
/// Subprocess-mode adapter.
fn finalize(
    status: std::process::ExitStatus,
    wall_time: Duration,
    warnings: Vec<String>,
    stderr_tail: Vec<String>,
) -> Result<RunReport, AdapterError> {
    let exit_code = status.code().unwrap_or(-1);
    if !status.success() {
        let stderr = if stderr_tail.is_empty() {
            format!("LinearFold exited {exit_code} with no stderr output")
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
        let info = LinearFoldAdapter::new().info();
        assert_eq!(info.id, "linearfold");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "Apache-2.0");
        assert_eq!(info.display_name, "LinearFold");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = LinearFoldAdapter::new().info();
        // LinearFold 1.0 is the canonical release line tagged in the
        // upstream repo; 2.0 reserves room for the next major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(1, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(2, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = LinearFoldAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.linearfold.fold"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = LinearFoldAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
