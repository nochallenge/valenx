//! # valenx-adapter-samtools
//!
//! Adapter for [samtools](https://www.htslib.org/) — Heng Li's
//! Swiss-army knife for SAM / BAM / CRAM data, the htslib companion
//! tool that powers nearly every short- and long-read sequencing
//! pipeline. This adapter wraps the four most common subcommands:
//! `view` (SAM↔BAM conversion), `sort` (coordinate sort), `index`
//! (BAM index sidecar), and `flagstat` (alignment QC summary).
//!
//! **Phase 18 — subprocess wrapper around `samtools <action>`.** The
//! user picks the subcommand via `action` in `[bio.samtools]` and
//! supplies an input SAM/BAM/CRAM. `prepare()` dispatches on the
//! action: `view` / `sort` / `index` use the shared subprocess runner
//! since their outputs go to files; `flagstat` writes its summary to
//! stdout, so we adopt the MAFFT-style stdout-capture pattern and
//! redirect to `flagstat.txt` in the workdir.
//!
//! ## Why a custom run() for flagstat
//!
//! Three of the four wrapped actions (`view`, `sort`, `index`) write
//! their output to a file — `view`/`sort` via `-o <file>`, `index`
//! via the `<input>.bai` sidecar — so the shared
//! [`valenx_core::subprocess::run`] runner that streams stdout
//! through a line handler works fine. `flagstat` doesn't have an
//! `-o` flag; its summary always goes to stdout. Routing that
//! through the line-handler runner would force us to either rebuild
//! the file from logged lines (lossy on the formatting characters
//! samtools uses to align columns) or buffer the entire output in
//! memory. Cleanest path mirrors MAFFT: a custom `run()` that
//! redirects stdout to a file via `Stdio::from(File::create(...))`
//! and pipes stderr through a per-line handler.
//!
//! Each-action command construction is factored into [`build_command`]
//! so the dispatch is one self-contained match without polluting the
//! Adapter::prepare body.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::fs;
use std::io::BufReader;
use std::path::{Path, PathBuf};
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
    subprocess::{self, KillOnDropChild, MAX_LINE_BYTES, SUBPROCESS_CHANNEL_CAPACITY},
    Adapter, AdapterError, AdapterInfo, Capabilities, Case, LicenseMode, Physics, PreparedJob,
    ProbeReport, RunContext, RunReport, VersionRange,
};
use valenx_fields::{
    artifact::{Artifact, ArtifactKind},
    Results,
};

use crate::case_input::SamtoolsInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(SamtoolsAdapter::new())
}

pub struct SamtoolsAdapter;

impl SamtoolsAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for SamtoolsAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "samtools";
/// samtools' binary candidates. Bioconda, Homebrew, and source
/// builds all install under the canonical name.
const BINARIES: &[&str] = &["samtools"];

/// The flagstat-output filename the adapter pins. samtools writes
/// flagstat to stdout; we redirect via `Stdio::from(File::create(...))`.
const OUT_FLAGSTAT: &str = "flagstat.txt";

/// How to drive a particular samtools action. The two shapes the
/// runner cares about: "stdout is the artifact, capture it to a
/// file" (flagstat) vs "outputs go to files already, just run"
/// (view / sort / index).
#[derive(Clone, Debug, PartialEq)]
pub enum Action {
    /// Standard subprocess execution — outputs land on disk via the
    /// command's own flags.
    Direct,
    /// Capture stdout to the given filename in the workdir.
    CaptureStdout(&'static str),
}

impl Adapter for SamtoolsAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "samtools",
            // samtools 1.x has been the stable line forever; 1.17
            // (2023) is the floor we test against — it carries the
            // modernised CLI we recommend (consistent `-@` for
            // threads, `-o` for output, `--threads` long form). The
            // upper bound 2.0 reserves room for the next major.
            version_range: VersionRange {
                min_inclusive: Version::new(1, 17, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "MIT",
            docs_url: "http://www.htslib.org/doc/samtools.html",
            homepage_url: "https://www.htslib.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `samtools --version` prints "samtools 1.17 ..." on
                // stdout; the combined scanner picks the version
                // up cleanly.
                let found_version = detect_tool_version_semver(&binary_path, &["--version"]);
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
                hint: "samtools 1.17+ required; install via `apt install samtools`, \
                       `brew install samtools`, or `conda install -c bioconda samtools`"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = SamtoolsInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve the input path against the case directory if
        // relative. Same convention as every other Phase 17/18 bio
        // adapter — `input = "aligned.sam"` next to `case.toml`.
        let source_input = if input.input.is_absolute() {
            input.input.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.input)?
        };
        if !source_input.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.samtools].input `{}` not found (resolved {})",
                    input.input.display(),
                    source_input.display()
                ),
            });
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "samtools 1.17+ required; install via `apt install samtools`, \
                       `brew install samtools`, or `conda install -c bioconda samtools`"
                .into(),
        })?;

        let (native_command, _action) = build_command(&binary_path, &source_input, &input);

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // BAM sort / view on a typical aligned-reads file runs
            // in minutes; large pipelines (whole-genome sort, CRAM
            // recompression) can run for an hour or more. 4 hours
            // is a generous default that covers the long tail.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        // Re-derive the action from the prepared command. Inspecting
        // the command preserves the action info across the
        // prepare/run boundary without baking action-specific state
        // into PreparedJob (which is shared across every adapter).
        // The first arg after the binary is always the subcommand
        // name.
        let action_kind = match job.native_command.get(1).and_then(|s| s.to_str()) {
            Some("flagstat") => Action::CaptureStdout(OUT_FLAGSTAT),
            _ => Action::Direct,
        };

        match action_kind {
            Action::Direct => {
                let report = subprocess::run(job, ctx, "starting samtools", |line| {
                    let mut hint = subprocess::Hint::default();
                    // samtools prints sparse progress on stderr (sort
                    // emits "[bam_sort_core] merging from N files..."
                    // when streaming temp files). The shared runner
                    // routes stderr through Warn-level logging on its
                    // own; this stdout handler just lifts a couple of
                    // markers if they appear.
                    if line.contains("Real time:") || line.contains("CPU time:") {
                        hint.progress = Some((95.0, line.to_string()));
                    } else if line.contains("[E::") || line.contains("ERROR") {
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
            Action::CaptureStdout(filename) => run_capture_stdout(job, ctx, filename),
        }
    }

    fn collect(&self, job: &PreparedJob) -> Result<Results, AdapterError> {
        // Re-derive the action from the prepared command (same
        // technique as `run()`).
        let action = job
            .native_command
            .get(1)
            .and_then(|s| s.to_str())
            .unwrap_or("");

        // Provenance: hash whichever output the action produces.
        // Falls back to case.toml when the run hasn't produced
        // anything yet — keeps the provenance block well-formed for
        // partial / failed runs.
        let case_hash_input = primary_output_path(job, action)
            .filter(|p| p.is_file())
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "samtools",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        match action {
            "view" | "sort" => {
                // The output path was recorded as the value after the
                // `-o` flag in the prepared command. It may live
                // outside the workdir (we accept absolute paths as
                // the case-relative resolution rule allows), but we
                // still surface it as an artifact.
                if let Some(out) = output_after_flag(job, "-o") {
                    if out.is_file() {
                        let label = match action {
                            "view" => format!("samtools view output ({})", out.display()),
                            "sort" => format!("samtools sort output ({})", out.display()),
                            _ => unreachable!(),
                        };
                        let kind = if out
                            .extension()
                            .and_then(|s| s.to_str())
                            .map(|s| s.eq_ignore_ascii_case("sam"))
                            .unwrap_or(false)
                        {
                            ArtifactKind::Tabular
                        } else {
                            ArtifactKind::Native
                        };
                        artefacts.push(Artifact {
                            path: out,
                            kind,
                            checksum: None,
                            label,
                        });
                    }
                }
            }
            "index" => {
                // `samtools index <input>` writes `<input>.bai` next
                // to the input. Locate the input (the last positional
                // before extra_args) and look for the sidecar.
                if let Some(input_path) = positional_input(job) {
                    let bai = with_extra_extension(&input_path, "bai");
                    if bai.is_file() {
                        artefacts.push(Artifact {
                            path: bai,
                            kind: ArtifactKind::Native,
                            checksum: None,
                            label: "samtools BAM index (.bai)".to_string(),
                        });
                    }
                }
            }
            "flagstat" => {
                // We pinned the captured stdout to `flagstat.txt` in
                // the workdir. Surface as a Log artifact.
                let flagstat = job.workdir.join(OUT_FLAGSTAT);
                if flagstat.is_file() {
                    artefacts.push(Artifact {
                        path: flagstat,
                        kind: ArtifactKind::Log,
                        checksum: None,
                        label: "samtools flagstat".to_string(),
                    });
                }
            }
            _ => {}
        }

        artefacts.sort_by(|a, b| a.path.cmp(&b.path));
        results.artifacts = artefacts;
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        // The bio-specific Capability variants land in a follow-up
        // task; ribbon contributions are already enough for the
        // registry to surface the adapter.
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.samtools.view"],
        }
    }
}

/// Compose the `samtools` invocation for the given action.
///
/// Returns `(native_command, Action)` — the command vector to run and
/// a discriminator that tells `run()` whether the standard subprocess
/// runner is enough (`Action::Direct`) or whether stdout needs to be
/// redirected to a file (`Action::CaptureStdout`).
///
/// Each action gets its own command shape:
///
/// - `view` / `sort`: `samtools <action> -@ N -o <output> <input> [extras...]`
/// - `index`:         `samtools index <input> [extras...]`
/// - `flagstat`:      `samtools flagstat <input> [extras...]` (stdout
///   captured into `flagstat.txt`)
pub fn build_command(
    binary_path: &Path,
    source_input: &Path,
    case: &SamtoolsInput,
) -> (Vec<OsString>, Action) {
    let mut cmd: Vec<OsString> = vec![
        binary_path.as_os_str().to_owned(),
        OsString::from(&case.action),
    ];

    match case.action.as_str() {
        "view" | "sort" => {
            cmd.push(OsString::from("-@"));
            cmd.push(OsString::from(case.threads.to_string()));
            // `output` is required by case_input parsing for view /
            // sort, so unwrap is safe here. Any future relaxation
            // would also need to update build_command.
            let output = case
                .output
                .as_ref()
                .expect("case_input enforces output for view/sort");
            cmd.push(OsString::from("-o"));
            cmd.push(output.as_os_str().to_owned());
            cmd.push(source_input.as_os_str().to_owned());
            for arg in &case.extra_args {
                cmd.push(OsString::from(arg));
            }
            (cmd, Action::Direct)
        }
        "index" => {
            // `samtools index` doesn't expose a `-@` flag in the same
            // shape — the index step is fast and single-threaded for
            // typical BAMs; pass extra_args through if the user wants
            // multi-threaded indexing (`-@`) explicitly.
            cmd.push(source_input.as_os_str().to_owned());
            for arg in &case.extra_args {
                cmd.push(OsString::from(arg));
            }
            (cmd, Action::Direct)
        }
        "flagstat" => {
            cmd.push(source_input.as_os_str().to_owned());
            for arg in &case.extra_args {
                cmd.push(OsString::from(arg));
            }
            (cmd, Action::CaptureStdout(OUT_FLAGSTAT))
        }
        // Defensive — case_input rejects unknown actions, but a
        // future schema-only change shouldn't reach `unreachable!()`.
        other => {
            panic!("samtools build_command: unsupported action `{other}` slipped past case_input")
        }
    }
}

/// MAFFT-style stdout-capture run: spawn the prepared command,
/// redirect stdout to `<filename>` in the workdir, pipe stderr
/// through a line handler that lifts samtools markers and tracks
/// warnings. Used only for the `flagstat` action.
fn run_capture_stdout(
    job: &PreparedJob,
    ctx: &mut RunContext,
    filename: &str,
) -> Result<RunReport, AdapterError> {
    if job.native_command.is_empty() {
        return Err(AdapterError::Other(anyhow::anyhow!(
            "PreparedJob.native_command is empty — prepare() should \
             have populated it"
        )));
    }

    // Open the output sink. samtools writes flagstat to stdout, so
    // we open `<filename>` for write and hand its FD to the child.
    // Any prior content from a previous run gets truncated.
    let out_path = job.workdir.join(filename);
    let out_file = std::fs::File::create(&out_path).map_err(|e| {
        AdapterError::Other(anyhow::anyhow!(
            "open {} for write: {e}",
            out_path.display()
        ))
    })?;

    let program = &job.native_command[0];
    let args: Vec<&OsString> = job.native_command.iter().skip(1).collect();

    ctx.report_progress(0.0, "starting samtools flagstat");
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

/// Mirror of the line-by-line stderr handling that
/// `subprocess::run` does for stdout. Logs the line and lifts
/// samtools' progress / error markers to coarse UI ticks.
fn process_stderr(line: &str, ctx: &mut RunContext<'_>, warnings: &mut Vec<String>) {
    ctx.log(LogLevel::Info, line);
    if line.contains("[E::") || line.contains("ERROR") {
        warnings.push(line.trim().to_string());
    }
}

/// Bridge a child's exit status into the adapter's RunReport / Run
/// error, mirroring `subprocess::finalize` so flagstat's failure
/// mode matches every other Subprocess-mode adapter.
fn finalize(
    status: std::process::ExitStatus,
    wall_time: Duration,
    warnings: Vec<String>,
    stderr_tail: Vec<String>,
) -> Result<RunReport, AdapterError> {
    let exit_code = status.code().unwrap_or(-1);
    if !status.success() {
        let stderr = if stderr_tail.is_empty() {
            format!("samtools exited {exit_code} with no stderr output")
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

/// Walk the prepared command for the value following `flag`. Used
/// from `collect()` to recover the `-o` output path so we can
/// surface it as an artifact.
fn output_after_flag(job: &PreparedJob, flag: &str) -> Option<PathBuf> {
    let mut iter = job.native_command.iter().peekable();
    while let Some(arg) = iter.next() {
        if arg.to_str() == Some(flag) {
            if let Some(val) = iter.next() {
                return Some(PathBuf::from(val));
            }
        }
    }
    None
}

/// Locate the positional input argument in a samtools command. We
/// know which slot it's in — the first non-flag positional after the
/// subcommand and any flag/value pairs. For the actions this adapter
/// produces, the input always lands as the first non-`-@` /
/// non-`-o` positional in the command.
fn positional_input(job: &PreparedJob) -> Option<PathBuf> {
    // Skip the binary, the subcommand, and any flag/value pairs.
    // For `index` the command is just `samtools index <input>
    // [extras]`, so the first arg after `index` is the input.
    let mut iter = job.native_command.iter().skip(2).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        if s.starts_with('-') {
            // For known flags that take a value (we don't expect any
            // here for `index`'s typical shape, but defensive: skip
            // the value if it isn't another flag), peek ahead.
            if matches!(s, "-@" | "-o") {
                let _ = iter.next();
            }
            continue;
        }
        return Some(PathBuf::from(arg));
    }
    None
}

/// Pick the primary output path to hash for provenance, given the
/// action. `view`/`sort` produce a file at `-o`; `index` produces a
/// `.bai` next to the input; `flagstat` produces `flagstat.txt` in
/// the workdir.
fn primary_output_path(job: &PreparedJob, action: &str) -> Option<PathBuf> {
    match action {
        "view" | "sort" => output_after_flag(job, "-o"),
        "index" => positional_input(job).map(|p| with_extra_extension(&p, "bai")),
        "flagstat" => Some(job.workdir.join(OUT_FLAGSTAT)),
        _ => None,
    }
}

/// Append `.<ext>` to a path. Used to derive `<input>.bai` from
/// `<input>` for the `index` action's sidecar.
fn with_extra_extension(path: &Path, ext: &str) -> PathBuf {
    let mut s = path.as_os_str().to_owned();
    s.push(".");
    s.push(ext);
    PathBuf::from(s)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = SamtoolsAdapter::new().info();
        assert_eq!(info.id, "samtools");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "MIT");
        assert_eq!(info.display_name, "samtools");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = SamtoolsAdapter::new().info();
        // 1.17 is the floor we test against; 2.0 reserves room for
        // an eventual major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(1, 17, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(2, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = SamtoolsAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.samtools.view"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = SamtoolsAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    #[test]
    fn build_command_dispatches_per_action() {
        // Confirm the Action discriminator threads correctly: view /
        // sort / index → Direct; flagstat → CaptureStdout. This is
        // the seam between prepare() and run() — getting it wrong
        // would silently route flagstat through the line-handler
        // runner and lose the artifact.
        let bin = PathBuf::from("/usr/bin/samtools");
        let src = PathBuf::from("/data/aligned.bam");

        let view = SamtoolsInput {
            action: "view".to_string(),
            input: src.clone(),
            output: Some(PathBuf::from("out.bam")),
            threads: 1,
            extra_args: Vec::new(),
        };
        let (_, kind) = build_command(&bin, &src, &view);
        assert_eq!(kind, Action::Direct);

        let sort = SamtoolsInput {
            action: "sort".to_string(),
            output: Some(PathBuf::from("sorted.bam")),
            ..view.clone()
        };
        let (_, kind) = build_command(&bin, &src, &sort);
        assert_eq!(kind, Action::Direct);

        let index = SamtoolsInput {
            action: "index".to_string(),
            output: None,
            ..view.clone()
        };
        let (_, kind) = build_command(&bin, &src, &index);
        assert_eq!(kind, Action::Direct);

        let flagstat = SamtoolsInput {
            action: "flagstat".to_string(),
            output: None,
            ..view
        };
        let (_, kind) = build_command(&bin, &src, &flagstat);
        assert_eq!(kind, Action::CaptureStdout(OUT_FLAGSTAT));
    }
}
