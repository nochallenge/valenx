//! License-mode declarations and the host-side enforcement hook.
//!
//! Every adapter declares one of three modes. The workspace-level
//! `cargo-deny` policy enforces the mode at build time; `valenx-core`
//! provides the runtime sanity-check hook adapters call before
//! touching the OS.

use serde::{Deserialize, Serialize};

/// How an adapter runs its external tool.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum LicenseMode {
    /// Tool shipped inside the Valenx binary (statically linked, or
    /// dynamically linked with a permissive license). OK for Apache /
    /// BSD / MIT / ISC / MPL code.
    Bundled,

    /// Tool linked at runtime to an OS-level `.so` / `.dll` / `.dylib`.
    /// OK for LGPL code; user may substitute their own build.
    DynamicLinked,

    /// Adapter is implemented entirely in Rust inside the Valenx
    /// process — no external binary, no FFI library. Probe is
    /// trivially true; run dispatches directly into a native crate.
    Native,

    /// Tool runs as a child process. The only interface is argv,
    /// environment, stdin/stdout/stderr, and files in a work dir. This
    /// is the mode GPL tools live in — we never link against them.
    Subprocess,
}

/// Strongly-typed guard on the process-spawn path: a `Bundled`
/// adapter must not shell out; a `Subprocess` adapter should.
///
/// Call this in the host helper that wraps `std::process::Command`
/// or `tokio::process::Command`; it bubbles up a structured error if
/// the declared mode disagrees with what's happening.
pub fn assert_spawn_allowed(mode: LicenseMode) -> Result<(), LicenseModeViolation> {
    match mode {
        LicenseMode::Bundled => Err(LicenseModeViolation::SpawnInBundled),
        LicenseMode::DynamicLinked | LicenseMode::Native | LicenseMode::Subprocess => Ok(()),
    }
}

/// What went wrong when the declared license mode and the code path
/// disagreed.
#[derive(Debug, thiserror::Error)]
pub enum LicenseModeViolation {
    #[error("adapter declared `Bundled` but attempted to spawn a child process")]
    SpawnInBundled,
}
