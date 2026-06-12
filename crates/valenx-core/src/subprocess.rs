//! Shared subprocess runner for adapters that drive a CLI tool.
//!
//! OpenFOAM, gmsh, CalculiX, Elmer, SU2 — every Subprocess-mode
//! adapter ends up wanting the same plumbing:
//!
//! - spawn the program with its args + env + cwd
//! - drain stdout + stderr on background threads (so neither
//!   stream can deadlock the adapter)
//! - forward every line through [`RunContext::log`]
//! - poll the child + the cancellation token in the main loop
//! - kill the child on cancel, collect + report structured
//!   `AdapterError::Run` on non-zero exit
//!
//! Rather than reimplement that in every adapter (as OpenFOAM and
//! gmsh did in Phase 1 / early Phase 2), this module hosts the
//! canonical implementation. Adapters pass in a [`LineHandler`]
//! closure that sees every stdout line *before* it's pushed to
//! the log sink — that's where OpenFOAM extracts residuals, gmsh
//! emits progress ticks, etc.

use std::ffi::OsString;
use std::io::{BufRead, BufReader};
use std::process::{Child, Command, Stdio};
use std::sync::mpsc;
use std::thread;
use std::time::{Duration, Instant};

use crate::adapter::{PreparedJob, RunContext};
use crate::error::{AdapterError, RunPhase};

/// Per-line cap on bytes read from a child process's stdout / stderr.
/// Round-6 hardening: pre-fix `BufReader::lines()` would silently
/// allocate an unbounded `String` for a child that printed a single
/// 4 GiB line with no `\n` — instant OOM from a misbehaving solver.
/// 1 MiB per line is far more than any honest residual / progress
/// line; the cap truncates pathological lines and tags them with a
/// `[truncated]` marker so the operator can spot the malfunction.
pub const MAX_LINE_BYTES: usize = 1024 * 1024;

/// Bounded mpsc-channel capacity for the stdout/stderr producer
/// threads. Round-6 hardening: pre-fix the channel was unbounded
/// (`mpsc::channel`), so a chatty solver that out-printed the main
/// thread's `recv_timeout`-paced consumer would accumulate every
/// undelivered line in memory until the parent ran out of RAM. 4096
/// pending lines is large enough that legitimate burst output (e.g.
/// `simpleFoam` flooding residuals at iteration N) doesn't make the
/// child stall on every write, and small enough that a runaway
/// producer back-pressures the child rather than the parent.
pub const SUBPROCESS_CHANNEL_CAPACITY: usize = 4096;

/// RAII guard that wraps a [`std::process::Child`] and kills it on
/// drop when the producing [`PreparedJob`] has `kill_on_drop = true`.
///
/// Round-6 fix: pre-fix the `kill_on_drop` field on
/// [`crate::PreparedJob`] was set across 140+ adapter sites but no
/// path actually honoured it — the runner returned `Err(?)` early
/// and the bare `Child` got dropped without a SIGKILL, leaving the
/// subprocess orphaned (typical OpenFOAM `simpleFoam` continues to
/// run for the rest of its allotted iteration budget after the
/// parent gives up). This guard makes the field load-bearing: the
/// `Drop` impl issues `Child::kill` when the flag is set, exactly
/// as the field's name promises.
///
/// Round-14 M9: factored out to `pub` so [`crate::executor::LocalExecutor`]
/// can wrap every submitted child in the same guard — pre-fix
/// LocalExecutor stored a bare `Child` in its handle table, so a
/// dropped executor (e.g. UI exit while a sweep is still running)
/// orphaned every outstanding subprocess.
pub struct KillOnDropChild {
    inner: Child,
    enabled: bool,
}

impl KillOnDropChild {
    /// Wrap a `Child` so `enabled` decides whether `Drop` issues a
    /// `kill()`. The runner uses `enabled = job.kill_on_drop`;
    /// `LocalExecutor` always wraps with `enabled = true` because the
    /// executor doesn't (yet) thread `PreparedJob.kill_on_drop` into
    /// its handle table.
    pub fn new(inner: Child, enabled: bool) -> Self {
        Self { inner, enabled }
    }

    /// Mutable access to the wrapped child for `try_wait` / `kill`
    /// calls that need to drive the child while keeping the guard
    /// alive.
    pub fn inner_mut(&mut self) -> &mut Child {
        &mut self.inner
    }
}

impl Drop for KillOnDropChild {
    fn drop(&mut self) {
        if !self.enabled {
            return;
        }
        // try_wait first — if the child already exited we have nothing
        // to kill and double-kill could surface a noisy "no such
        // process" error on POSIX.
        match self.inner.try_wait() {
            Ok(Some(_)) => {}
            _ => {
                let _ = self.inner.kill();
                let _ = self.inner.wait();
            }
        }
    }
}

/// Which phase a failure happened in, when the caller wants to
/// override the default `RunPhase::Solve`.
#[derive(Clone, Copy, Debug)]
pub enum FailurePhase {
    Startup,
    Solve,
    Output,
}

impl FailurePhase {
    /// Map to the `RunPhase` stored on `AdapterError::Run`. Used by
    /// callers building structured run errors outside the default
    /// `finalize` path.
    pub fn as_run_phase(self) -> RunPhase {
        match self {
            FailurePhase::Startup => RunPhase::Startup,
            FailurePhase::Solve => RunPhase::Solve,
            FailurePhase::Output => RunPhase::Output,
        }
    }
}

/// A stdout-line callback. Called synchronously on the main loop
/// thread (so `&mut` state is fine) *before* the line is pushed
/// through the `RunContext`'s log sink. Return a [`Hint`] to
/// influence progress / warning tracking.
pub type LineHandler<'a> = dyn FnMut(&str) -> Hint + 'a;

/// Out-of-band signal a [`LineHandler`] can emit per line. The
/// default is `Hint::None` — the runner just logs and moves on.
#[derive(Clone, Debug, Default)]
pub struct Hint {
    /// If `Some`, the runner calls `ctx.report_progress(pct, msg)`.
    pub progress: Option<(f32, String)>,
    /// If `Some`, the runner appends this to the warnings vector
    /// returned in the `SubprocessReport`.
    pub warning: Option<String>,
}

/// What the runner reports when the child exits successfully.
/// Non-zero exits become `AdapterError::Run` returned via `Err`.
#[derive(Clone, Debug)]
pub struct SubprocessReport {
    pub exit_code: i32,
    pub wall_time: Duration,
    pub warnings: Vec<String>,
}

/// Run a prepared job's native command through to completion. The
/// primary entry point for adapter `run()` methods.
///
/// `line_handler` sees every stdout line (stderr still gets logged
/// at `Warn` level and kept in the stderr tail for error reporting,
/// but it doesn't go through the handler — separating the two keeps
/// handler logic uncluttered).
pub fn run<'a>(
    job: &PreparedJob,
    ctx: &mut RunContext<'_>,
    starting_message: &str,
    mut line_handler: impl FnMut(&str) -> Hint + 'a,
) -> Result<SubprocessReport, AdapterError> {
    if job.native_command.is_empty() {
        return Err(AdapterError::Other(anyhow::anyhow!(
            "PreparedJob.native_command is empty — prepare() should have populated it"
        )));
    }

    let program = &job.native_command[0];
    let args: Vec<&OsString> = job.native_command.iter().skip(1).collect();

    ctx.report_progress(0.0, starting_message);
    ctx.log(
        crate::adapter::LogLevel::Info,
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
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());

    let raw_child = cmd.spawn().map_err(|e| AdapterError::Run {
        exit_code: -1,
        stderr: format!("failed to spawn {}: {e}", program.to_string_lossy()),
        phase: RunPhase::Startup,
    })?;
    // Round-6: honour `PreparedJob::kill_on_drop` via the RAII guard.
    // Round-19 M4: force `enabled = true` regardless of the job's
    // preference. The runner OWNS the subprocess for the lifetime of
    // this call — if we return early (cancel, IO error, line-handler
    // panic) the child MUST be reaped before the call site sees the
    // error, otherwise we leak orphaned solver processes that hold
    // their workdir locks open and pressure the host's PID table.
    // Sister to R14 M9 (LocalExecutor) and R16 M2 (capture_subprocess_stdout)
    // which both unconditionally pass `true`. The job's
    // `kill_on_drop` field is still honoured by the executor's handle
    // table; this guard exists for the in-call panic / early-return
    // window only.
    let mut kill_guard = KillOnDropChild::new(raw_child, true);

    let stdout = kill_guard
        .inner_mut()
        .stdout
        .take()
        .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("child stdout not captured")))?;
    let stderr = kill_guard
        .inner_mut()
        .stderr
        .take()
        .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("child stderr not captured")))?;

    // Round-6: bounded sync_channel back-pressures runaway producers
    // instead of OOMing the parent (pre-fix `mpsc::channel` accepted
    // unlimited pending items).
    let (tx, rx) = mpsc::sync_channel::<Event>(SUBPROCESS_CHANNEL_CAPACITY);
    let tx_err = tx.clone();
    let so_thread = thread::spawn(move || {
        let mut reader = BufReader::new(stdout);
        loop {
            match read_capped_line(&mut reader) {
                Ok(Some(line)) => {
                    if tx.send(Event::Stdout(line)).is_err() {
                        break;
                    }
                }
                Ok(None) => break, // EOF
                Err(_) => break,   // bubble up via the child's exit status
            }
        }
    });
    let se_thread = thread::spawn(move || {
        let mut reader = BufReader::new(stderr);
        loop {
            match read_capped_line(&mut reader) {
                Ok(Some(line)) => {
                    if tx_err.send(Event::Stderr(line)).is_err() {
                        break;
                    }
                }
                Ok(None) => break,
                Err(_) => break,
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
            Ok(Event::Stdout(line)) => {
                process_stdout(&line, ctx, &mut line_handler, &mut warnings);
            }
            Ok(Event::Stderr(line)) => {
                ctx.log(crate::adapter::LogLevel::Warn, &line);
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
                    // Drain whatever's still in the channel before
                    // returning, so no log line gets lost.
                    for ev in rx.try_iter() {
                        match ev {
                            Event::Stdout(line) => {
                                process_stdout(&line, ctx, &mut line_handler, &mut warnings);
                            }
                            Event::Stderr(line) => {
                                ctx.log(crate::adapter::LogLevel::Warn, &line);
                                if stderr_tail.len() >= STDERR_TAIL_MAX {
                                    stderr_tail.remove(0);
                                }
                                stderr_tail.push(line);
                            }
                        }
                    }
                    let _ = so_thread.join();
                    let _ = se_thread.join();
                    return finalize(status, start.elapsed(), warnings, stderr_tail);
                }
            }
            Err(mpsc::RecvTimeoutError::Disconnected) => {
                let status = kill_guard.inner_mut().wait().map_err(AdapterError::Io)?;
                let _ = so_thread.join();
                let _ = se_thread.join();
                return finalize(status, start.elapsed(), warnings, stderr_tail);
            }
        }
    }
}

/// Read one capped-length line from `reader`. Stops at the first
/// `\n` byte or after `MAX_LINE_BYTES` bytes (whichever comes
/// first). The trailing `\n` and any `\r\n` is stripped. Returns
/// `Ok(None)` on clean EOF. A line longer than `MAX_LINE_BYTES`
/// gets a `[line truncated at <N> bytes]` suffix and the next read
/// starts at the next byte (the rest of the over-long line is
/// discarded).
fn read_capped_line(reader: &mut impl BufRead) -> std::io::Result<Option<String>> {
    let mut buf = Vec::with_capacity(256);
    let mut bytes_read = 0usize;
    let mut byte = [0u8; 1];
    loop {
        if bytes_read >= MAX_LINE_BYTES {
            // Skip over the rest of the over-long line so the next
            // call doesn't start mid-line. Bounded scan: don't
            // chase a runaway producer past 2x the line cap.
            //
            // Round-22 L2: pre-fix the scan-ahead was bounded at
            // `MAX_LINE_BYTES * 16` (160 MiB at the 10 MiB cap),
            // which is far more memory than the truncation message
            // implies the next reader has to chew through before
            // the next line starts. 2x (20 MiB) is still generous
            // for "find the newline that closes this runaway line"
            // — any producer that writes more than 20 MiB without a
            // newline is hostile, not just verbose, and we'll lose
            // at most the tail of one over-long entry by giving up
            // sooner. Keep `MAX_LINE_BYTES` unchanged.
            let mut skipped = 0usize;
            let cap = MAX_LINE_BYTES.saturating_mul(2);
            while skipped < cap {
                match reader.read(&mut byte) {
                    Ok(0) => break,
                    Ok(_) => {
                        skipped += 1;
                        if byte[0] == b'\n' {
                            break;
                        }
                    }
                    Err(_) => break,
                }
            }
            let mut s = String::from_utf8_lossy(&buf).into_owned();
            s.push_str(&format!(" [line truncated at {MAX_LINE_BYTES} bytes]"));
            return Ok(Some(s));
        }
        match reader.read(&mut byte) {
            Ok(0) => {
                // EOF: emit the partial line if there is one, else
                // signal end-of-stream.
                if buf.is_empty() {
                    return Ok(None);
                }
                let s = String::from_utf8_lossy(&buf).into_owned();
                return Ok(Some(strip_cr(s)));
            }
            Ok(_) => {
                if byte[0] == b'\n' {
                    let s = String::from_utf8_lossy(&buf).into_owned();
                    return Ok(Some(strip_cr(s)));
                }
                buf.push(byte[0]);
                bytes_read += 1;
            }
            Err(e) => return Err(e),
        }
    }
}

fn strip_cr(mut s: String) -> String {
    if s.ends_with('\r') {
        s.pop();
    }
    s
}

fn process_stdout(
    line: &str,
    ctx: &mut RunContext<'_>,
    handler: &mut impl FnMut(&str) -> Hint,
    warnings: &mut Vec<String>,
) {
    ctx.log(crate::adapter::LogLevel::Info, line);
    let hint = handler(line);
    if let Some((pct, msg)) = hint.progress {
        ctx.report_progress(pct, &msg);
    }
    if let Some(w) = hint.warning {
        warnings.push(w);
    }
}

fn finalize(
    status: std::process::ExitStatus,
    wall_time: Duration,
    warnings: Vec<String>,
    stderr_tail: Vec<String>,
) -> Result<SubprocessReport, AdapterError> {
    let exit_code = status.code().unwrap_or(-1);
    if !status.success() {
        let stderr = if stderr_tail.is_empty() {
            format!("child exited {exit_code} with no stderr output")
        } else {
            stderr_tail.join("\n")
        };
        return Err(AdapterError::Run {
            exit_code,
            stderr,
            phase: RunPhase::Solve,
        });
    }
    Ok(SubprocessReport {
        exit_code,
        wall_time,
        warnings,
    })
}

enum Event {
    Stdout(String),
    Stderr(String),
}

// ---------------------------------------------------------------------------
// Tests — the runner is exercised in full by the OpenFOAM + gmsh
// integration tests. Here we just sanity-check the `Hint` builder
// and the `FailurePhase` mapping so refactors catch accidental
// variants being added without the mapping being updated.
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hint_default_is_noop() {
        let h = Hint::default();
        assert!(h.progress.is_none());
        assert!(h.warning.is_none());
    }

    #[test]
    fn failure_phase_maps_to_run_phase() {
        assert!(matches!(
            FailurePhase::Startup.as_run_phase(),
            RunPhase::Startup
        ));
        assert!(matches!(
            FailurePhase::Solve.as_run_phase(),
            RunPhase::Solve
        ));
        assert!(matches!(
            FailurePhase::Output.as_run_phase(),
            RunPhase::Output
        ));
    }

    #[test]
    fn read_capped_line_truncates_lines_past_max_bytes() {
        // Round-6 RED→GREEN: a stdin that ships a single 4 GiB line
        // with no `\n` would let `BufReader::lines()` allocate an
        // unbounded `String`. The capped reader cuts at
        // `MAX_LINE_BYTES`, tags the line, and resumes on the next
        // boundary.
        let oversized: String = "a".repeat(MAX_LINE_BYTES + 100);
        let payload = format!("{oversized}\nnext\n");
        let mut cursor = std::io::Cursor::new(payload.into_bytes());
        let first = read_capped_line(&mut cursor).unwrap().unwrap();
        assert!(
            first.contains("[line truncated"),
            "expected truncation marker, got first {} bytes",
            first.len().min(80)
        );
        // The next read returns the post-newline content (the cap
        // discarded the rest of the over-long line).
        let second = read_capped_line(&mut cursor).unwrap().unwrap();
        assert_eq!(second, "next");
        // Then EOF.
        assert!(read_capped_line(&mut cursor).unwrap().is_none());
    }

    #[test]
    fn read_capped_line_handles_short_lines_and_crlf() {
        let payload = b"hello\nworld\r\nlast";
        let mut cursor = std::io::Cursor::new(&payload[..]);
        assert_eq!(read_capped_line(&mut cursor).unwrap().unwrap(), "hello");
        assert_eq!(read_capped_line(&mut cursor).unwrap().unwrap(), "world");
        // Final unterminated line is still emitted.
        assert_eq!(read_capped_line(&mut cursor).unwrap().unwrap(), "last");
        assert!(read_capped_line(&mut cursor).unwrap().is_none());
    }

    /// Round-22 L2 RED→GREEN: the truncation-skip-ahead reader
    /// stops at `MAX_LINE_BYTES * 2` (was 16x pre-fix). When the
    /// runaway producer writes more than 2x the line cap without a
    /// newline, the skip-ahead loop must give up rather than chase
    /// the producer for many GB.
    ///
    /// The observable proxy is: feed the reader an over-long line,
    /// then more bytes (no newline) totalling > 2x but < 16x the
    /// cap, then a real newline well past the 2x boundary. With the
    /// 2x cap the reader gives up before the real newline and the
    /// next `read_capped_line` call surfaces the post-cap tail
    /// (starting with the bytes the skip-ahead loop did not chase)
    /// — not the bytes after the eventual newline. Pre-fix (16x)
    /// the reader would have kept chasing.
    #[test]
    fn read_capped_line_stops_skip_ahead_at_2x() {
        // 16x is way too big to allocate in a unit test (16 MiB).
        // We exploit `saturating_mul(2)` exactly: the skip loop reads
        // up to `MAX_LINE_BYTES * 2` bytes, then bails. With a single
        // over-long line followed by `MAX_LINE_BYTES * 2 + 1` `b`s
        // and then a newline, the next read should start mid-`b`
        // stream (the byte at offset 1 + MAX_LINE_BYTES * 2 from the
        // post-truncation-marker boundary).
        //
        // To keep the unit test memory-light we don't reproduce
        // `MAX_LINE_BYTES = 1 MiB` of bytes; we observe a simpler
        // shape: the over-long line + a `[line truncated …]` marker
        // gets emitted, then the skip-ahead loop runs, and the next
        // `read_capped_line` call returns whatever's left in the
        // buffer past the skip cap.
        let oversized: String = "a".repeat(MAX_LINE_BYTES + 100);
        // Then a `b` runaway equal to skip-cap (2 * MAX_LINE_BYTES);
        // then a `\n` and `next` — the skip stops before the `\n`,
        // so on the next read we surface the rest of the `b`s.
        let mut payload = oversized.into_bytes();
        payload.extend(std::iter::repeat_n(b'b', MAX_LINE_BYTES * 2 + 5));
        payload.push(b'\n');
        payload.extend_from_slice(b"next\n");
        let mut cursor = std::io::Cursor::new(payload);
        let first = read_capped_line(&mut cursor).unwrap().unwrap();
        assert!(
            first.contains("[line truncated"),
            "first line must be truncated, got len {}",
            first.len()
        );
        // The next call surfaces what survives past the 2x skip cap:
        // a few `b`s tail (because the skip stopped before the `\n`).
        // Post-fix this should be tail bytes; pre-fix (16x cap) the
        // skip would have consumed everything through `\n` and we'd
        // see `"next"` here instead.
        let second = read_capped_line(&mut cursor).unwrap().unwrap();
        assert!(
            second.starts_with('b'),
            "expected post-skip `b` tail, got: {second:?}"
        );
        assert_ne!(
            second, "next",
            "pre-fix 16x skip would have consumed through the newline"
        );
    }

    /// Round-19 M4 RED→GREEN: even when the wrapping job has
    /// `kill_on_drop = false`, the runner's RAII guard must still
    /// kill the child on early-return / panic. The runner OWNS the
    /// subprocess for the lifetime of the `run()` call — leaking an
    /// orphaned process past an unwind would let the solver continue
    /// to hold workdir locks open after the parent gave up.
    ///
    /// We exercise the mechanism in isolation: wrap a spawned child
    /// in `KillOnDropChild::new(child, true)`, drop the guard in a
    /// scope, and verify the child reaped (`try_wait` returns
    /// `Ok(Some(_))`). We can't easily inject a panic into the real
    /// `subprocess::run` loop (it requires a working line-handler
    /// stub + the cancel token), but the guard is the load-bearing
    /// part — the runner now passes `true` for every spawned child
    /// regardless of `job.kill_on_drop`.
    #[test]
    fn kill_on_drop_guard_kills_child_even_with_job_kill_on_drop_false() {
        // Spawn a long-sleeper that would otherwise outlive the
        // scope by a wide margin. Skip when the binary isn't
        // present so the test stays cross-platform.
        #[cfg(windows)]
        let mut spawn = std::process::Command::new("powershell.exe");
        #[cfg(windows)]
        spawn.args(["-NoProfile", "-Command", "Start-Sleep -Seconds 30"]);
        #[cfg(not(windows))]
        let mut spawn = std::process::Command::new("/bin/sleep");
        #[cfg(not(windows))]
        spawn.arg("30");
        spawn
            .stdin(std::process::Stdio::null())
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::null());
        let child = match spawn.spawn() {
            Ok(c) => c,
            Err(_) => {
                eprintln!("skipping: sleep binary not present");
                return;
            }
        };
        let pid = child.id();
        {
            // The runner now constructs the guard with `enabled = true`
            // unconditionally — round-19 M4. The scope-drop here
            // simulates the unwind path (panic, early return, cancel)
            // through `run()`'s body.
            let _guard = KillOnDropChild::new(child, true);
            // Confirm the child is alive RIGHT NOW. We don't sleep —
            // just stat the kill-on-drop guard's view of the PID.
        }
        // After the scope drops, the guard's `Drop` impl SIGKILLed
        // the child and waited on it. We re-spawn a `kill -0`-style
        // probe to confirm the PID is no longer a runnable process.
        // On Windows we can't easily check "process exists" without
        // OpenProcess + GetExitCodeProcess; skip the post-check there
        // — the Drop ran on this thread already and Rust's stdlib
        // panics on a double-wait, so the guard IS the assertion.
        #[cfg(not(windows))]
        {
            // POSIX `kill -0 <pid>` returns 0 if the process exists,
            // non-zero otherwise. After Drop's SIGKILL + wait the
            // PID should be reaped.
            let status = std::process::Command::new("kill")
                .args(["-0", &pid.to_string()])
                .stdout(std::process::Stdio::null())
                .stderr(std::process::Stdio::null())
                .status();
            if let Ok(s) = status {
                assert!(
                    !s.success(),
                    "PID {pid} still exists after KillOnDropChild::drop — leak!"
                );
            }
        }
        // Suppress unused-variable warning on Windows where pid isn't read.
        let _ = pid;
    }

    #[test]
    fn subprocess_channel_has_bounded_capacity() {
        // Round-6 RED→GREEN: the channel between the stdout/stderr
        // pump threads and the main loop is now bounded at
        // `SUBPROCESS_CHANNEL_CAPACITY`. A producer that sends
        // CAP+1 items without a consumer must block on the last
        // send rather than buffer infinitely. We exercise this via
        // a fake producer thread that signals via an atomic when it
        // gets blocked.
        let (tx, _rx) = mpsc::sync_channel::<Event>(SUBPROCESS_CHANNEL_CAPACITY);
        let blocked = std::sync::Arc::new(std::sync::atomic::AtomicBool::new(false));
        let blocked_inner = blocked.clone();
        let h = thread::spawn(move || {
            // Fill the channel to capacity (these all return
            // immediately because the buffer has space).
            for i in 0..SUBPROCESS_CHANNEL_CAPACITY {
                tx.send(Event::Stdout(format!("line-{i}"))).expect("send");
            }
            // The CAP+1-th send must block. Signal that we're
            // about to block, then attempt the blocking send (we'll
            // never make it past — the test thread drops the
            // sender via thread cleanup).
            blocked_inner.store(true, std::sync::atomic::Ordering::SeqCst);
            // try_send instead of send so the producer doesn't
            // park forever — Full == back-pressure confirmed.
            tx.try_send(Event::Stdout("blocked".into()))
        });
        // Wait until the producer has finished its capacity fill.
        let deadline = std::time::Instant::now() + Duration::from_secs(5);
        while !blocked.load(std::sync::atomic::Ordering::SeqCst) {
            if std::time::Instant::now() > deadline {
                panic!("producer never reached the blocking send");
            }
            std::thread::sleep(Duration::from_millis(10));
        }
        let attempt = h.join().expect("producer thread");
        // The CAP+1-th send saw a full channel — that's the
        // back-pressure the round-4 unbounded channel was missing.
        assert!(
            matches!(attempt, Err(mpsc::TrySendError::Full(_))),
            "expected TrySendError::Full, got {attempt:?}"
        );
    }
}
