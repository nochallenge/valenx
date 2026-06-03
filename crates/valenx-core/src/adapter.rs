//! The `Adapter` trait every integrated tool implements.
//!
//! Spec: [RFC 0002](../../rfcs/0002-adapter-contract.md). This file is
//! the authoritative Rust surface; the RFC is the plain-English
//! version.

use std::ffi::OsString;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use std::time::Duration;

use semver::Version;
use serde::{Deserialize, Serialize};

use valenx_fields::Results;

use crate::error::{AdapterError, RunPhase};
use crate::license::LicenseMode;
use crate::physics::{Capability, Physics};

/// Static metadata an adapter publishes.
#[derive(Clone, Debug)]
pub struct AdapterInfo {
    pub id: &'static str,
    pub display_name: &'static str,
    pub version_range: VersionRange,
    pub physics: &'static [Physics],
    pub license_mode: LicenseMode,
    pub tool_license: &'static str,
    pub docs_url: &'static str,
    pub homepage_url: &'static str,
}

/// An inclusive-lower, exclusive-upper SemVer range the adapter has
/// been tested against.
#[derive(Clone, Debug)]
pub struct VersionRange {
    pub min_inclusive: Version,
    pub max_exclusive: Version,
}

impl VersionRange {
    /// `true` if `v` falls inside the half-open interval
    /// `[min_inclusive, max_exclusive)`.
    pub fn contains(&self, v: &Version) -> bool {
        v >= &self.min_inclusive && v < &self.max_exclusive
    }
}

/// What `probe()` reports.
#[derive(Clone, Debug)]
pub struct ProbeReport {
    pub ok: bool,
    pub found_version: Option<Version>,
    pub binary_path: Option<PathBuf>,
    pub warnings: Vec<String>,
    pub required_env: Vec<(&'static str, String)>,
}

impl ProbeReport {
    /// Sentinel "tool isn't installed" result — `ok = false`, all
    /// optional fields cleared.
    pub fn not_found() -> Self {
        Self {
            ok: false,
            found_version: None,
            binary_path: None,
            warnings: Vec::new(),
            required_env: Vec::new(),
        }
    }
}

/// A prepared-but-not-yet-run unit of work.
#[derive(Clone, Debug)]
pub struct PreparedJob {
    pub workdir: PathBuf,
    pub native_command: Vec<OsString>,
    pub environment: Vec<(OsString, OsString)>,
    pub estimated_runtime: Option<Duration>,
    /// When `true`, the [`subprocess::run`](crate::subprocess::run)
    /// path wraps the spawned child in a SIGKILL-on-Drop guard so an
    /// early `return Err(?)` (or a panic that unwinds through the
    /// runner) does not orphan the subprocess. Round-6 fix: pre-fix
    /// this field was set across 140+ adapter sites but no path
    /// actually honoured it; the `std::process::Child` returned by
    /// `Command::spawn` just got dropped on early-return, leaving
    /// the OS process running until it exited on its own. The
    /// runner now consumes the field at spawn time. The
    /// [`Executor`](crate::Executor) path is unaffected — it
    /// already kills children explicitly via the
    /// `submit/poll/cancel` handle lifecycle, so the field there
    /// remains advisory.
    pub kill_on_drop: bool,
}

/// Report from `run()`. Residual samples and warnings travel inside.
#[derive(Clone, Debug, Default)]
pub struct RunReport {
    pub exit_code: i32,
    pub wall_time: Duration,
    pub converged: Option<bool>,
    pub residual_history: Vec<ResidualSample>,
    pub warnings: Vec<String>,
    pub final_phase: Option<RunPhase>,
}

/// One sample of residual history — one iteration, one field,
/// one value.
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct ResidualSample {
    pub iteration: u64,
    pub field: &'static str,
    pub value: f64,
}

/// Cooperative cancellation passed into `run()`.
#[derive(Clone, Debug, Default)]
pub struct CancellationToken(Arc<AtomicBool>);

impl CancellationToken {
    /// New uncancelled token.
    pub fn new() -> Self {
        Self::default()
    }
    /// Flip the token to "cancelled". Subsequent
    /// [`Self::is_cancelled`] calls return `true`. Idempotent.
    pub fn cancel(&self) {
        self.0.store(true, Ordering::SeqCst);
    }
    /// `true` once any clone of the token has had [`Self::cancel`]
    /// called on it. Cheap atomic load; safe to poll from the adapter
    /// loop.
    pub fn is_cancelled(&self) -> bool {
        self.0.load(Ordering::SeqCst)
    }
}

/// Progress sink — the adapter calls `report` periodically.
pub trait ProgressSink: Send + Sync {
    fn report(&self, pct: f32, message: &str);
}

/// Log sink — the adapter forwards the underlying tool's output.
pub trait LogSink: Send + Sync {
    fn log_line(&self, level: LogLevel, line: &str);
}

/// Log severity.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LogLevel {
    Trace,
    Debug,
    Info,
    Warn,
    Error,
}

/// Runtime context passed to `run()`.
pub struct RunContext<'a> {
    pub cancel: &'a CancellationToken,
    pub progress: Box<dyn ProgressSink + 'a>,
    pub log: Box<dyn LogSink + 'a>,
}

impl RunContext<'_> {
    /// Propagate `AdapterError::Cancelled` if the token is set.
    pub fn check_cancel(&self) -> Result<(), AdapterError> {
        if self.cancel.is_cancelled() {
            Err(AdapterError::Cancelled)
        } else {
            Ok(())
        }
    }

    /// Report progress between 0.0 and 100.0 with a short message.
    pub fn report_progress(&self, pct: f32, message: &str) {
        self.progress.report(pct.clamp(0.0, 100.0), message);
    }

    /// Emit a log line at the given level.
    pub fn log(&self, level: LogLevel, line: &str) {
        self.log.log_line(level, line);
    }
}

/// UI-facing capability map — informs ribbon composition and menu
/// enablement without peeking at each adapter's internals.
#[derive(Clone, Debug, Default)]
pub struct Capabilities {
    pub capabilities: Vec<Capability>,
    pub ribbon_contributions: Vec<&'static str>,
}

// ---------------------------------------------------------------------------
// The contract
// ---------------------------------------------------------------------------

/// Opaque handle for a canonical `Case` — resolved to a concrete type
/// in `valenx-core::project` once that lands.
#[derive(Clone, Debug)]
pub struct Case {
    pub id: String,
    pub path: PathBuf,
}

/// The trait every integrated tool implements.
pub trait Adapter: Send + Sync {
    fn info(&self) -> AdapterInfo;
    fn probe(&self) -> Result<ProbeReport, AdapterError>;

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError>;

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError>;

    fn collect(&self, job: &PreparedJob) -> Result<Results, AdapterError>;

    /// Optional quick validation before queueing. Default: OK.
    fn validate(&self, case: &Case) -> Result<(), AdapterError> {
        let _ = case;
        Ok(())
    }

    /// Capability map. Default: empty.
    fn capabilities(&self) -> Capabilities {
        Capabilities::default()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn version_range_contains() {
        let r = VersionRange {
            min_inclusive: Version::parse("1.0.0").unwrap(),
            max_exclusive: Version::parse("2.0.0").unwrap(),
        };
        assert!(r.contains(&Version::parse("1.5.0").unwrap()));
        assert!(r.contains(&Version::parse("1.0.0").unwrap()));
        assert!(!r.contains(&Version::parse("2.0.0").unwrap()));
        assert!(!r.contains(&Version::parse("0.9.9").unwrap()));
    }

    #[test]
    fn cancellation_token() {
        let tok = CancellationToken::new();
        assert!(!tok.is_cancelled());
        tok.cancel();
        assert!(tok.is_cancelled());
    }
}
