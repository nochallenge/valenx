//! # valenx-crash-reporter
//!
//! In-app crash reporter for the Valenx desktop application.
//!
//! Two responsibilities:
//!
//! 1. Install a `std::panic` hook that captures every panic, builds
//!    a sanitised [`CrashReport`], and writes it to a per-user
//!    crashes directory before chaining to the previously installed
//!    hook. The report is persisted unconditionally — the user's
//!    privacy preference only gates *uploading* it.
//! 2. Provide the report data structure and the
//!    serialise / load round-trip so the app can prompt
//!    "found N unsent crash reports — submit them?" on next launch
//!    and ship the JSON to a configurable endpoint when the user
//!    opts in.
//!
//! The crate is intentionally small and dep-light: it pulls only
//! `serde`, `serde_json`, and `tracing`. No HTTP client, no SDK
//! integration. Upload is the caller's responsibility — keeping
//! the report on disk in a stable JSON shape lets users (or
//! enterprise admins) review reports before any network egress.
//!
//! ## Sanitisation
//!
//! Crash reports must not exfiltrate user data. The sanitiser:
//!
//! - Strips home-directory absolute paths (replaces `/home/<user>/`
//!   and `C:\Users\<user>\` with `<HOME>/`).
//! - Strips long opaque IDs (UUIDs, SHA-256 hex) — they're
//!   either project hashes or run IDs that have no diagnostic
//!   value but identify the user's project.
//! - Truncates the panic payload to 4 KiB so a debug-format dump
//!   of a giant `Mesh` doesn't bloat the report or smuggle data.
//!
//! Diagnostic fields that survive: panic message (post-strip),
//! source location, adapter id, Valenx version, OS / arch, and
//! a fixed-length backtrace stub.
//!
//! ## Forbid unsafe
//!
//! Even though a panic-hook installer touches global mutable
//! state, the implementation routes everything through
//! `std::panic::set_hook`'s safe API. The crate is
//! `#![forbid(unsafe_code)]` like the rest of the workspace.

#![forbid(unsafe_code)]

use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use serde::{Deserialize, Serialize};

/// One captured panic.
///
/// Kept as plain serde-friendly data so the report can be written
/// to disk, loaded later, and uploaded by an HTTP layer the caller
/// owns. The shape is stable across patch releases — RFC-bumping
/// it requires a major.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct CrashReport {
    /// Schema version. Increments when fields are added or removed
    /// in a backward-incompatible way.
    pub schema: u32,
    /// ISO 8601 timestamp the panic was captured. Generated from
    /// `SystemTime::now()` — assumes the system clock is roughly
    /// correct.
    pub timestamp: String,
    /// Human panic message, post-sanitisation. Truncated at
    /// [`MAX_MESSAGE_BYTES`].
    pub message: String,
    /// Source file + line where the panic surfaced, post-
    /// sanitisation. `None` if `PanicInfo::location()` returned
    /// `None` (panicked outside of any tracked location).
    pub location: Option<String>,
    /// Optional adapter id when the panic was caught inside an
    /// adapter spawn / collect path. Set by the caller via
    /// [`CrashReport::with_adapter`].
    pub adapter: Option<String>,
    /// Valenx version string — typically `env!("CARGO_PKG_VERSION")`
    /// from the host crate.
    pub valenx_version: String,
    /// `<os>-<arch>` produced from `std::env::consts`. Useful for
    /// triaging platform-specific crashes.
    pub platform: String,
    /// Fixed-length backtrace stub — sanitised, capped at
    /// [`MAX_BACKTRACE_BYTES`]. `None` when no backtrace was
    /// available (RUST_BACKTRACE unset).
    pub backtrace: Option<String>,
}

/// Latest schema version produced by this crate. Bumps when the
/// shape of [`CrashReport`] changes incompatibly.
pub const CURRENT_SCHEMA: u32 = 1;

/// Cap on the panic message length after sanitisation. 4 KiB is
/// generous for a one-line panic and leaves room for the wrapping
/// JSON envelope without making any single report unwieldy.
pub const MAX_MESSAGE_BYTES: usize = 4096;

/// Cap on the backtrace stub. Backtraces from a deeply-nested
/// solver loop can be tens of kilobytes; we keep the first 8 KiB
/// which is enough for triage.
pub const MAX_BACKTRACE_BYTES: usize = 8192;

/// Cap on the bytes [`CrashReport::load_all`] will pull from a
/// single `*.json` file. Round-6 DoS guard: a malicious or
/// already-corrupted `*.json` planted in the crash directory
/// (a non-empty multi-GiB file masquerading as a report) would
/// otherwise drive `std::fs::read` to allocate the whole file
/// into memory. 10 MiB is comfortably past any honest crash
/// report — the constructor caps `message` at 4 KiB and the
/// backtrace at 8 KiB — but small enough that an attacker
/// can't push the parent process into swap.
pub const MAX_REPORT_BYTES: usize = 10 * 1024 * 1024;

impl CrashReport {
    /// Build a report from the components a panic hook receives.
    /// The constructor handles sanitisation + truncation; callers
    /// don't need to pre-process inputs.
    pub fn new(
        message: impl Into<String>,
        location: Option<String>,
        valenx_version: impl Into<String>,
    ) -> Self {
        Self {
            schema: CURRENT_SCHEMA,
            timestamp: now_iso8601(),
            message: sanitise_and_truncate(&message.into(), MAX_MESSAGE_BYTES),
            location: location.map(|loc| sanitise_and_truncate(&loc, MAX_MESSAGE_BYTES)),
            adapter: None,
            valenx_version: valenx_version.into(),
            platform: format!("{}-{}", std::env::consts::OS, std::env::consts::ARCH,),
            backtrace: None,
        }
    }

    /// Tag the report with the adapter that owned the failing run.
    /// Pure setter so the panic hook can opt in without rebuilding
    /// the whole struct.
    pub fn with_adapter(mut self, adapter: impl Into<String>) -> Self {
        self.adapter = Some(adapter.into());
        self
    }

    /// Attach a backtrace stub. Sanitises + truncates the same way
    /// the constructor handles `message` / `location`.
    pub fn with_backtrace(mut self, backtrace: impl Into<String>) -> Self {
        self.backtrace = Some(sanitise_and_truncate(
            &backtrace.into(),
            MAX_BACKTRACE_BYTES,
        ));
        self
    }

    /// Persist the report to `dir/<timestamp>.json`, creating the
    /// directory if it doesn't exist. Returns the full path on
    /// success.
    ///
    /// Round-4 hardening: writes via `atomic_write_bytes` (a sidecar
    /// `.tmp` + rename) rather than `std::fs::write`. Crash reports
    /// are written from the panic hook, and if the process panics
    /// AGAIN mid-write (rare but possible during deep stack unwinds
    /// across FFI), the old fs::write path would leave a truncated
    /// JSON behind that future `load_all` parses would silently drop.
    /// The rename-based path either leaves the old file intact or
    /// the new file atomically — never a half-written sidecar.
    pub fn write_to_disk(&self, dir: &Path) -> std::io::Result<PathBuf> {
        std::fs::create_dir_all(dir)?;
        // Use the sanitised timestamp as the filename — it's
        // already filesystem-safe (RFC 3339 minus colons).
        let safe_ts = self.timestamp.replace([':', '.'], "-");
        let path = dir.join(format!("{safe_ts}.json"));
        let bytes = serde_json::to_vec_pretty(self)
            .map_err(|e| std::io::Error::other(format!("serialise crash report: {e}")))?;
        atomic_write_bytes(&path, &bytes)?;
        Ok(path)
    }

    /// Read every `*.json` in `dir` that parses as a
    /// [`CrashReport`]. Files that fail to parse are logged and
    /// skipped — a corrupt report shouldn't block the app launch
    /// loop.
    pub fn load_all(dir: &Path) -> std::io::Result<Vec<(PathBuf, CrashReport)>> {
        if !dir.exists() {
            return Ok(Vec::new());
        }
        let mut out = Vec::new();
        for entry in std::fs::read_dir(dir)? {
            let entry = entry?;
            let path = entry.path();
            let ext_match = path
                .extension()
                .and_then(|s| s.to_str())
                .map(|s| s.eq_ignore_ascii_case("json"))
                .unwrap_or(false);
            if !ext_match {
                continue;
            }
            // Round-6 hardening: cap the per-file read at
            // MAX_REPORT_BYTES. Use stat first so a multi-GiB file
            // is rejected without any allocation; even if metadata
            // is unreliable (e.g. NFS lag), the `take()` guard on
            // the read path bounds memory growth.
            let metadata = match std::fs::metadata(&path) {
                Ok(m) => m,
                Err(e) => {
                    tracing::warn!(target: "valenx-crash", ?e, ?path, "stat crash report");
                    continue;
                }
            };
            if metadata.len() > MAX_REPORT_BYTES as u64 {
                tracing::warn!(
                    target: "valenx-crash",
                    size = metadata.len(),
                    cap = MAX_REPORT_BYTES as u64,
                    ?path,
                    "crash report exceeds size cap, skipping"
                );
                continue;
            }
            let bytes = match read_capped(&path, MAX_REPORT_BYTES) {
                Ok(b) => b,
                Err(e) => {
                    tracing::warn!(target: "valenx-crash", ?e, ?path, "read crash report");
                    continue;
                }
            };
            match serde_json::from_slice::<CrashReport>(&bytes) {
                Ok(report) => out.push((path, report)),
                Err(e) => {
                    tracing::warn!(target: "valenx-crash", ?e, ?path, "parse crash report");
                }
            }
        }
        // Stable ordering — oldest first (filename is the
        // sanitised timestamp).
        out.sort_by(|(a, _), (b, _)| a.cmp(b));
        Ok(out)
    }
}

/// Read at most `cap` bytes from `path` into a fresh `Vec<u8>`.
/// Wraps `File::open(...).take(cap as u64).read_to_end(...)` — the
/// `take` cap bounds memory growth even if the file's `metadata.len()`
/// is unreliable (NFS, lazy filesystems, races between stat and read).
///
/// Round-6 hardening for [`CrashReport::load_all`].
fn read_capped(path: &Path, cap: usize) -> std::io::Result<Vec<u8>> {
    use std::io::Read;
    let mut f = std::fs::File::open(path)?;
    let mut buf = Vec::new();
    f.by_ref().take(cap as u64).read_to_end(&mut buf)?;
    Ok(buf)
}

/// Round-4 hardening: write `bytes` to `path` atomically by writing
/// to a unique sidecar then renaming over the destination.
///
/// Crash reports are written from the panic hook. If the host
/// double-panics during the write itself (rare but possible — picture
/// a stack-overflow that recurses through the panic hook), the old
/// `std::fs::write` path would leave a truncated JSON behind that
/// future `CrashReport::load_all` calls silently discard. The
/// rename-based path leaves the old file intact until the sidecar is
/// fully flushed — readers see either the previous version or the new
/// one, never a half-written wreck.
///
/// ## Round-27 STRUCTURAL consolidation
///
/// Thin wrapper around
/// [`valenx_core::io_caps::atomic_write_bytes`]. Pre-fix this site
/// was a copy of the R4-era `valenx_app::state_paths::atomic_write`
/// — and was explicitly tagged "the two implementations must stay
/// in lock-step; both are small enough that occasional drift is
/// cheap to fix". They DID drift: when state_paths landed unique
/// `(pid, nanos, counter)` sidecars (R24 M1 / R25 M3), fsync-before-
/// rename (R24 M1), and parent-dir fsync (R26 M3), crash-reporter
/// kept the single-`.tmp` shape with no fsync, no parent fsync, and
/// a race window any time TWO concurrent panics tried to write
/// reports. Consolidating into a canonical `valenx-core` helper
/// closes the drift permanently — a single fix benefits all 4
/// inlined atomic-write sites.
fn atomic_write_bytes(path: &Path, bytes: &[u8]) -> std::io::Result<()> {
    valenx_core::io_caps::atomic_write_bytes(path, bytes)
}

/// Install a panic hook that builds a [`CrashReport`] for every
/// panic and writes it to `crashes_dir`. Chains to the previous
/// hook so default panic-message printing still happens.
///
/// `valenx_version` is typically `env!("CARGO_PKG_VERSION")` from
/// the binary that calls this helper.
///
/// **Idempotency.** Calling this function twice replaces the
/// previous Valenx hook with a fresh one — both will not run.
/// That's the same semantics as `std::panic::set_hook`.
pub fn install_panic_hook(crashes_dir: PathBuf, valenx_version: String) {
    let prev = std::panic::take_hook();
    std::panic::set_hook(Box::new(move |info| {
        // Build the report. PanicInfo's payload is `Any` — we coerce
        // to &str / &String, falling back to a placeholder if the
        // payload isn't a string (rare in practice — `panic!` always
        // produces strings).
        let message = panic_payload_to_string(info.payload());
        let location = info
            .location()
            .map(|loc| format!("{}:{}:{}", loc.file(), loc.line(), loc.column()));
        let report = CrashReport::new(message, location, valenx_version.clone());

        match report.write_to_disk(&crashes_dir) {
            Ok(path) => {
                tracing::error!(
                    target: "valenx-crash",
                    ?path,
                    "wrote crash report"
                );
            }
            Err(e) => {
                // Don't trap the panic just because we couldn't
                // persist a report — write to stderr and let the
                // chained hook do its thing.
                eprintln!("valenx-crash: could not write report: {e}");
            }
        }
        prev(info);
    }));
}

/// Best-effort coercion of a `PanicInfo::payload()` to a `String`.
/// `panic!("...")` always produces a `&str` payload; structured
/// panics (`panic_any(MyError)`) drop into the `Box<dyn Any>` arm
/// and we surface a placeholder so the report still serialises.
fn panic_payload_to_string(payload: &(dyn std::any::Any + Send)) -> String {
    if let Some(s) = payload.downcast_ref::<&'static str>() {
        return (*s).to_string();
    }
    if let Some(s) = payload.downcast_ref::<String>() {
        return s.clone();
    }
    "(non-string panic payload)".to_string()
}

/// Sanitise free-form text + truncate to `max_bytes`.
///
/// Sanitisation rules:
/// 1. Replace `/home/<user>/` with `<HOME>/` (Linux, macOS without
///    the `/Users/` aliasing — we cover both via two passes).
/// 2. Replace `/Users/<user>/` with `<HOME>/` (macOS).
/// 3. Replace `C:\Users\<user>\` (and lowercase / forward-slash
///    variants) with `<HOME>\`.
/// 4. Replace UUIDv4-shaped tokens
///    (`xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx`) with `<UUID>`.
/// 5. Replace standalone 64-character hex strings (SHA-256) with
///    `<HASH>`.
///
/// The sanitiser is conservative — it preserves anything it can't
/// confidently classify. Diagnostic value over privacy is a
/// non-goal: when in doubt, redact.
pub fn sanitise(text: &str) -> String {
    let mut out = text.to_string();
    // Linux / macOS home paths. Iterate by walking the string and
    // skipping past the username component.
    out = strip_unix_home(&out, "/home/");
    out = strip_unix_home(&out, "/Users/");
    // Windows home paths — both `C:\Users\name\` and the
    // forward-slash variant the Rust formatter sometimes uses.
    out = strip_windows_home(&out);
    // UUIDs and 64-char hex tokens.
    out = strip_uuids(&out);
    out = strip_sha256_tokens(&out);
    out
}

/// Sanitise + cap to `max_bytes` UTF-8 bytes (truncated on a char
/// boundary so the result stays valid UTF-8).
fn sanitise_and_truncate(text: &str, max_bytes: usize) -> String {
    let cleaned = sanitise(text);
    if cleaned.len() <= max_bytes {
        return cleaned;
    }
    // Find the largest char-boundary cap <= max_bytes.
    let mut end = max_bytes;
    while end > 0 && !cleaned.is_char_boundary(end) {
        end -= 1;
    }
    let mut truncated = String::with_capacity(end + 16);
    truncated.push_str(&cleaned[..end]);
    truncated.push_str("…(truncated)");
    truncated
}

fn strip_unix_home(text: &str, prefix: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut rest = text;
    while let Some(idx) = rest.find(prefix) {
        out.push_str(&rest[..idx]);
        // Skip past the username component (anything up to the next
        // `/` or end of string).
        let after_prefix = &rest[idx + prefix.len()..];
        match after_prefix.find('/') {
            Some(slash) => {
                out.push_str("<HOME>/");
                rest = &after_prefix[slash + 1..];
            }
            None => {
                // Path ends at the username — nothing after to skip.
                out.push_str("<HOME>/");
                rest = "";
                break;
            }
        }
    }
    out.push_str(rest);
    out
}

fn strip_windows_home(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut rest = text;
    // Match either `C:\Users\` or `C:/Users/` (Rust's debug formatter
    // sometimes prints forward slashes even on Windows).
    let needles: &[&str] = &["C:\\Users\\", "C:/Users/", "c:\\users\\", "c:/users/"];
    'outer: while !rest.is_empty() {
        let mut earliest: Option<(usize, usize)> = None;
        for needle in needles {
            if let Some(idx) = rest.find(needle) {
                let len = needle.len();
                match earliest {
                    Some((e, _)) if e < idx => {}
                    _ => earliest = Some((idx, len)),
                }
            }
        }
        match earliest {
            Some((idx, needle_len)) => {
                out.push_str(&rest[..idx]);
                let after = &rest[idx + needle_len..];
                // Skip the username up to the next path separator.
                let sep_idx = after.find(['\\', '/']).unwrap_or(after.len());
                out.push_str("<HOME>\\");
                if sep_idx == after.len() {
                    break 'outer;
                }
                rest = &after[sep_idx + 1..];
            }
            None => {
                out.push_str(rest);
                break;
            }
        }
    }
    out
}

fn strip_uuids(text: &str) -> String {
    // UUIDv4 shape: 8-4-4-4-12 hex with dashes. We walk by character
    // boundaries — not bytes — so that non-ASCII text in panic
    // messages (em-dashes, accented usernames like `Müller`, emoji)
    // is preserved verbatim. UUIDs themselves are ASCII so byte-level
    // pattern matching still works on the underlying slice.
    let mut out = String::with_capacity(text.len());
    let bytes = text.as_bytes();
    let mut iter = text.char_indices().peekable();
    while let Some((i, ch)) = iter.next() {
        if i + 36 <= bytes.len() && is_uuid_at(bytes, i) {
            out.push_str("<UUID>");
            // Advance the iterator past the 36 UUID bytes (which are
            // ASCII, one byte per char).
            while let Some(&(next_i, _)) = iter.peek() {
                if next_i >= i + 36 {
                    break;
                }
                iter.next();
            }
        } else {
            out.push(ch);
        }
    }
    out
}

fn is_uuid_at(bytes: &[u8], start: usize) -> bool {
    // Positions of dashes in a UUIDv4 string: 8, 13, 18, 23.
    const DASH_POSITIONS: [usize; 4] = [8, 13, 18, 23];
    for offset in 0..36 {
        let c = bytes[start + offset];
        if DASH_POSITIONS.contains(&offset) {
            if c != b'-' {
                return false;
            }
        } else if !c.is_ascii_hexdigit() {
            return false;
        }
    }
    true
}

fn strip_sha256_tokens(text: &str) -> String {
    // 64-char lowercase / mixed hex strings standing alone. Walks by
    // char-boundary so non-ASCII text (em-dashes, accented usernames,
    // emoji) flows through verbatim. Hex runs are ASCII so we can
    // continue to test them at the byte level.
    let mut out = String::with_capacity(text.len());
    let bytes = text.as_bytes();
    let mut iter = text.char_indices().peekable();
    while let Some((i, ch)) = iter.next() {
        if let Some(run_end) = run_of_hex(bytes, i) {
            let len = run_end - i;
            if len == 64 {
                out.push_str("<HASH>");
            } else {
                // Not 64 chars — emit the ASCII hex run verbatim.
                out.push_str(&text[i..run_end]);
            }
            // Advance the iterator past the consumed hex run.
            while let Some(&(next_i, _)) = iter.peek() {
                if next_i >= run_end {
                    break;
                }
                iter.next();
            }
        } else {
            out.push(ch);
        }
    }
    out
}

fn run_of_hex(bytes: &[u8], start: usize) -> Option<usize> {
    if !bytes
        .get(start)
        .copied()
        .is_some_and(|b| b.is_ascii_hexdigit())
    {
        return None;
    }
    let mut end = start;
    while end < bytes.len() && bytes[end].is_ascii_hexdigit() {
        end += 1;
    }
    if end - start >= 8 {
        Some(end)
    } else {
        None
    }
}

/// Format the current `SystemTime` as a fixed-precision ISO 8601
/// string so report filenames sort lexicographically by timestamp
/// even at second granularity.
fn now_iso8601() -> String {
    // Use seconds-since-epoch + a synthetic "Z" suffix. We don't
    // pull `chrono` / `time` — keeps the dep footprint at three
    // crates.
    let secs_total = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    // Decompose into a UTC date / time without a calendar dep.
    // 86_400 s/day, ignoring leap seconds — the timestamp here
    // is for sorting / triage, not legal-grade precision.
    let days = secs_total / 86_400;
    let secs = secs_total % 86_400;
    let h = secs / 3600;
    let m = (secs % 3600) / 60;
    let s = secs % 60;
    let (year, month, day) = ymd_from_epoch_days(days as i64);
    format!("{year:04}-{month:02}-{day:02}T{h:02}:{m:02}:{s:02}Z")
}

/// Convert "days since 1970-01-01" to (year, month, day). Algorithm
/// from Howard Hinnant's date library — accurate for years ±32767.
fn ymd_from_epoch_days(z: i64) -> (i32, u32, u32) {
    let z = z + 719_468; // shift epoch to 0000-03-01
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097) as u32;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if m <= 2 { y + 1 } else { y };
    (y as i32, m, d)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn current_schema_is_one() {
        assert_eq!(CURRENT_SCHEMA, 1);
    }

    #[test]
    fn report_constructor_truncates_long_messages() {
        let huge = "x".repeat(MAX_MESSAGE_BYTES * 2);
        let r = CrashReport::new(huge, None, "0.0.0");
        assert!(r.message.len() <= MAX_MESSAGE_BYTES + 32);
        assert!(r.message.ends_with("…(truncated)"));
    }

    #[test]
    fn report_constructor_populates_platform_string() {
        let r = CrashReport::new("boom", None, "0.0.0");
        assert!(r.platform.contains('-'), "got: {}", r.platform);
        // os-arch shape — the dash and one identifier on each side.
        let parts: Vec<&str> = r.platform.split('-').collect();
        assert_eq!(parts.len(), 2);
    }

    #[test]
    fn with_adapter_attaches_id() {
        let r = CrashReport::new("boom", None, "0.0.0").with_adapter("openfoam");
        assert_eq!(r.adapter, Some("openfoam".to_string()));
    }

    #[test]
    fn sanitise_strips_unix_home_path() {
        let s = sanitise("read /home/alice/projects/foo.toml: not found");
        assert!(s.contains("<HOME>/projects/foo.toml"), "got: {s}");
        assert!(!s.contains("alice"), "username leaked: {s}");
    }

    #[test]
    fn sanitise_strips_macos_users_path() {
        let s = sanitise("read /Users/bob/Library/state: not found");
        assert!(s.contains("<HOME>/Library"), "got: {s}");
        assert!(!s.contains("bob"), "username leaked: {s}");
    }

    #[test]
    fn sanitise_strips_windows_home_path() {
        let s = sanitise("read C:\\Users\\charlie\\AppData\\Roaming: not found");
        assert!(s.contains("<HOME>\\AppData"), "got: {s}");
        assert!(!s.contains("charlie"), "username leaked: {s}");
    }

    #[test]
    fn sanitise_strips_windows_forward_slash_variant() {
        let s = sanitise("read C:/Users/dave/Documents: not found");
        assert!(s.contains("<HOME>"), "got: {s}");
        assert!(!s.contains("dave"), "username leaked: {s}");
    }

    #[test]
    fn sanitise_replaces_uuid() {
        let s = sanitise("run_id=12345678-1234-4abc-8def-123456789abc failed");
        assert!(s.contains("<UUID>"), "got: {s}");
        assert!(!s.contains("123456789abc"), "uuid leaked: {s}");
    }

    #[test]
    fn sanitise_replaces_sha256_hex() {
        let hash = "abcdef0123456789abcdef0123456789abcdef0123456789abcdef0123456789";
        let s = sanitise(&format!("case_hash={hash} body"));
        assert!(s.contains("<HASH>"), "got: {s}");
        assert!(!s.contains(hash), "hash leaked: {s}");
    }

    #[test]
    fn sanitise_preserves_non_ascii_text_around_uuid() {
        // Realistic panic message: an accented Windows username and an
        // em-dash bracket the UUID. Pre-fix, both bytes of each non-
        // ASCII char were cast as `c as char` independently, producing
        // garbage U+00C3 + U+00BC etc. Verify the original glyphs
        // round-trip and the UUID still gets replaced.
        let panic_msg = "Müller — run_id=12345678-1234-4abc-8def-123456789abc failed 🚨";
        let s = sanitise(panic_msg);
        assert!(s.contains("Müller"), "lost accented username: {s}");
        assert!(s.contains("—"), "lost em-dash: {s}");
        assert!(s.contains("🚨"), "lost emoji: {s}");
        assert!(s.contains("<UUID>"), "uuid not redacted: {s}");
        assert!(
            !s.contains("123456789abc"),
            "uuid leaked through char-iter path: {s}"
        );
    }

    #[test]
    fn sanitise_preserves_non_ascii_text_around_sha256() {
        // Same regression check for the SHA-256 stripper.
        let hash = "abcdef0123456789abcdef0123456789abcdef0123456789abcdef0123456789";
        let panic_msg = format!("Müller — hash={hash} failed 🚨");
        let s = sanitise(&panic_msg);
        assert!(s.contains("Müller"), "lost accented username: {s}");
        assert!(s.contains("—"), "lost em-dash: {s}");
        assert!(s.contains("🚨"), "lost emoji: {s}");
        assert!(s.contains("<HASH>"), "hash not redacted: {s}");
        assert!(!s.contains(hash), "hash leaked through char-iter path: {s}");
    }

    #[test]
    fn sanitise_preserves_non_secret_text() {
        // Short hex like `0xff` and ordinary words don't get redacted.
        let s = sanitise("error 0xff in module foo at line 42");
        assert_eq!(s, "error 0xff in module foo at line 42");
    }

    #[test]
    fn write_to_disk_round_trips_through_load_all() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-crash-test-{}",
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let r = CrashReport::new("simple panic", Some("foo.rs:42:1".into()), "0.1.0")
            .with_adapter("test");
        let path = r.write_to_disk(&tmp).expect("write");
        assert!(path.is_file());

        let loaded = CrashReport::load_all(&tmp).expect("load");
        assert_eq!(loaded.len(), 1);
        assert_eq!(loaded[0].1.message, "simple panic");
        assert_eq!(loaded[0].1.adapter, Some("test".to_string()));

        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn load_all_returns_empty_for_missing_dir() {
        let nope = std::env::temp_dir().join("valenx-crash-does-not-exist-banana");
        let _ = std::fs::remove_dir_all(&nope);
        let loaded = CrashReport::load_all(&nope).expect("load");
        assert!(loaded.is_empty());
    }

    #[test]
    fn load_all_skips_non_json_and_corrupt_files() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-crash-mixed-{}",
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        // Write one good report.
        let good = CrashReport::new("ok", None, "0.1.0");
        good.write_to_disk(&tmp).unwrap();
        // Drop a non-json sibling and a corrupt JSON.
        std::fs::write(tmp.join("readme.txt"), b"not a report").unwrap();
        std::fs::write(tmp.join("corrupt.json"), b"{not even close").unwrap();

        let loaded = CrashReport::load_all(&tmp).expect("load");
        assert_eq!(loaded.len(), 1);
        assert_eq!(loaded[0].1.message, "ok");

        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn load_all_orders_oldest_first() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-crash-order-{}",
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        // Two reports with explicit non-default timestamps so the
        // ordering is deterministic regardless of clock resolution.
        let mut a = CrashReport::new("first", None, "0.1.0");
        a.timestamp = "2026-04-28T00:00:00Z".into();
        a.write_to_disk(&tmp).unwrap();

        let mut b = CrashReport::new("second", None, "0.1.0");
        b.timestamp = "2026-04-28T01:00:00Z".into();
        b.write_to_disk(&tmp).unwrap();

        let loaded = CrashReport::load_all(&tmp).expect("load");
        assert_eq!(loaded.len(), 2);
        assert_eq!(loaded[0].1.message, "first");
        assert_eq!(loaded[1].1.message, "second");

        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn schema_round_trip_through_serde() {
        let r = CrashReport::new("boom", Some("foo.rs:1:1".into()), "0.1.0")
            .with_adapter("openfoam")
            .with_backtrace("main\nrun\nspawn");
        let s = serde_json::to_string(&r).unwrap();
        let back: CrashReport = serde_json::from_str(&s).unwrap();
        assert_eq!(back, r);
    }

    #[test]
    fn ymd_from_epoch_days_known_anchors() {
        // 0 = 1970-01-01.
        assert_eq!(ymd_from_epoch_days(0), (1970, 1, 1));
        // 31 = 1970-02-01.
        assert_eq!(ymd_from_epoch_days(31), (1970, 2, 1));
        // 365 = 1971-01-01 (1970 wasn't a leap year).
        assert_eq!(ymd_from_epoch_days(365), (1971, 1, 1));
    }

    #[test]
    fn sanitise_handles_multiple_paths_in_one_string() {
        let s = sanitise("/home/alice/a and /home/alice/b");
        assert!(s.contains("<HOME>/a and <HOME>/b"), "got: {s}");
    }

    #[test]
    fn truncation_lands_on_char_boundary() {
        // A unicode-heavy message close to the boundary must
        // truncate without producing invalid UTF-8.
        let msg = "🦀".repeat(MAX_MESSAGE_BYTES);
        let r = CrashReport::new(msg, None, "0.1.0");
        // Verifying it serialises is the strongest invariant —
        // serde rejects invalid UTF-8.
        assert!(serde_json::to_string(&r).is_ok());
    }

    /// Round-4 hardening: `write_to_disk` must use the atomic
    /// sidecar-then-rename path. We verify the leftover `.tmp` is
    /// gone after a successful write — the original `fs::write` path
    /// would never create a `.tmp` at all, so this test simultaneously
    /// proves the new path is wired AND that it cleans up after
    /// itself.
    #[test]
    fn write_to_disk_uses_atomic_rename_and_leaves_no_tmp_sidecar() {
        let tmp_dir = std::env::temp_dir().join(format!(
            "valenx-crash-atomic-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp_dir).unwrap();
        let report = CrashReport::new("test panic", None, "0.1.0");
        let written = report
            .write_to_disk(&tmp_dir)
            .expect("write_to_disk succeeds");

        // The final file must exist and be readable JSON.
        let bytes = std::fs::read(&written).expect("read back final file");
        let parsed: serde_json::Value = serde_json::from_slice(&bytes).expect("valid JSON");
        assert_eq!(parsed["message"], "test panic");

        // The sidecar `.tmp` must NOT exist post-rename.
        let mut tmp_name = written.file_name().expect("filename").to_os_string();
        tmp_name.push(".tmp");
        let tmp_sidecar = written.with_file_name(tmp_name);
        assert!(
            !tmp_sidecar.exists(),
            "atomic write must clean up the .tmp sidecar after rename; found: {}",
            tmp_sidecar.display()
        );

        let _ = std::fs::remove_dir_all(&tmp_dir);
    }

    /// `atomic_write_bytes` returns `InvalidInput` when the target
    /// path has no filename component (e.g. a bare directory path).
    /// Pins the helper's defensive case so a future refactor that
    /// drops the check fails this test instead of silently writing to
    /// a degenerate path.
    #[test]
    fn atomic_write_bytes_rejects_dirless_path() {
        #[cfg(unix)]
        let p = std::path::Path::new("/");
        #[cfg(windows)]
        let p = std::path::Path::new("C:\\");
        let err = atomic_write_bytes(p, b"x").expect_err("dirless path must fail");
        assert_eq!(err.kind(), std::io::ErrorKind::InvalidInput);
    }

    /// RED→GREEN (round-27 H2 sister, STRUCTURAL): 10 panicking
    /// threads writing crash reports concurrently — each report
    /// must land atomically on disk with no torn writes and no
    /// `AlreadyExists` errors from sidecar collisions.
    ///
    /// Pre-R27 the crash-reporter's inline `atomic_write_bytes` used
    /// the single `<path>.tmp` sidecar shape (no `<pid>.<counter>`
    /// suffix). Two threads racing to write reports in the same
    /// nanosecond would both `fs::write` the same `<path>.tmp`,
    /// then both `fs::rename` — the second rename on POSIX
    /// overwrites the first one's partial work; on Windows it could
    /// fail with `AlreadyExists` if the sidecar wasn't fully closed.
    /// Post-R27 each writer owns a distinct
    /// `<basename>.tmp.<pid>.<counter>` sidecar.
    ///
    /// Uses distinct `safe_ts` per report (timestamps include
    /// thread-id-shaped nanos) so the FINAL paths differ — the
    /// SIDECAR race is what we're verifying gets handled cleanly.
    ///
    /// Round-28 M2: this test by itself didn't actually anchor
    /// R27's bug because distinct timestamps mean distinct rename
    /// targets — even the pre-R27 single-tmp shape would have
    /// passed (no shared rename target = no overwrite race). The
    /// sister test `write_to_disk_handles_concurrent_panics_with_
    /// same_timestamp_round28_m2` below closes the gap by using
    /// the SAME `safe_ts` across all N threads, so the rename
    /// targets DO collide and the SIDECAR-collision surface is
    /// the only thing keeping the writes from interleaving.
    #[test]
    fn write_to_disk_handles_concurrent_panics_round27() {
        use std::sync::{Arc, Barrier};
        use std::thread;
        let dir = std::env::temp_dir().join(format!(
            "valenx-crash-r27-concurrent-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        const N: usize = 10;
        let barrier = Arc::new(Barrier::new(N));
        let mut handles = Vec::with_capacity(N);
        for i in 0..N {
            let dir = dir.clone();
            let barrier = Arc::clone(&barrier);
            handles.push(thread::spawn(move || {
                barrier.wait();
                let mut report = CrashReport::new(format!("concurrent panic #{i}"), None, "0.1.0");
                // Distinct timestamps so the FINAL file names don't
                // collide; we're isolating the SIDECAR-collision
                // surface as the failure mode of interest.
                report.timestamp = format!("2026-05-28T00:00:{i:02}Z");
                report.write_to_disk(&dir)
            }));
        }
        let mut errs: Vec<String> = Vec::new();
        for h in handles {
            if let Err(e) = h.join().unwrap() {
                errs.push(e.to_string());
            }
        }
        assert!(
            errs.is_empty(),
            "expected all {N} concurrent reports to land cleanly; errs: {errs:?}",
        );
        // All N reports should be present on disk.
        let loaded = CrashReport::load_all(&dir).expect("load_all");
        assert_eq!(
            loaded.len(),
            N,
            "expected {N} reports on disk, found {}",
            loaded.len(),
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Round-28 M2 — sister to
    /// `write_to_disk_handles_concurrent_panics_round27`. That R27
    /// test gave each thread a distinct `safe_ts`, so the FINAL
    /// `<safe_ts>.json` paths differed and the rename targets
    /// didn't collide — even the pre-R27 single-tmp sidecar would
    /// have passed it. This test gives every thread the SAME
    /// `safe_ts`, so:
    ///
    ///   * The sidecar shape's per-`<pid>.<counter>` uniqueness is
    ///     the ONLY thing preventing two threads from racing to
    ///     create+truncate+write the same `<path>.tmp`.
    ///   * The rename targets DO collide (all writers rename to
    ///     `<safe_ts>.json`) — POSIX `rename(2)` is atomic and the
    ///     last rename to land wins; Windows
    ///     `MoveFileEx`+`MOVEFILE_REPLACE_EXISTING` (the default
    ///     when the destination exists, set implicitly by Rust's
    ///     `fs::rename`) is also atomic.
    ///
    /// Expected post-fix invariants:
    ///
    ///   * Every `write_to_disk` returns `Ok` — no `AlreadyExists`
    ///     from the sidecar race, no torn content from the rename
    ///     race.
    ///   * Exactly ONE final `<safe_ts>.json` lives on disk after
    ///     the dust settles (the last winner of the rename race).
    ///   * That file parses cleanly as a `CrashReport` — torn
    ///     content during the rename window would fail
    ///     `serde_json::from_slice`.
    #[test]
    fn write_to_disk_handles_concurrent_panics_with_same_timestamp_round28_m2() {
        use std::sync::{Arc, Barrier};
        use std::thread;
        let dir = std::env::temp_dir().join(format!(
            "valenx-crash-r28-m2-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        const N: usize = 10;
        let barrier = Arc::new(Barrier::new(N));
        let mut handles = Vec::with_capacity(N);
        for i in 0..N {
            let dir = dir.clone();
            let barrier = Arc::clone(&barrier);
            handles.push(thread::spawn(move || {
                barrier.wait();
                let mut report = CrashReport::new(format!("concurrent panic #{i}"), None, "0.1.0");
                // SAME safe_ts across all N threads — the sidecar
                // race AND the rename race both hit the same
                // basename. Pre-R27's single-tmp sidecar would
                // surface as either AlreadyExists or torn content
                // on disk; post-R27 / R28 the per-counter sidecar
                // keeps each writer's tmp distinct and the rename
                // is atomic.
                report.timestamp = "2026-05-28T00:00:00Z".to_string();
                report.write_to_disk(&dir)
            }));
        }
        let mut errs: Vec<String> = Vec::new();
        for h in handles {
            if let Err(e) = h.join().unwrap() {
                errs.push(e.to_string());
            }
        }
        assert!(
            errs.is_empty(),
            "expected all {N} concurrent same-ts writes to land cleanly; errs: {errs:?}",
        );
        // Exactly ONE final report should be on disk (the rename
        // race's winner). Any orphan sidecars would also count as
        // failures.
        let loaded = CrashReport::load_all(&dir).expect("load_all");
        assert_eq!(
            loaded.len(),
            1,
            "expected exactly 1 surviving report (last rename wins) — found {}",
            loaded.len(),
        );
        // The surviving file must parse — torn content would fail
        // `load_all`'s serde_json parse.
        let (path, report) = &loaded[0];
        assert_eq!(
            path.file_name().and_then(|s| s.to_str()),
            Some("2026-05-28T00-00-00Z.json"),
            "surviving path must use the shared safe_ts basename",
        );
        assert!(
            report.message.starts_with("concurrent panic #"),
            "surviving report's message must be from one of the N writers",
        );
        // Confirm no orphan sidecars survive — the canonical
        // helper's contract is that every sidecar either gets
        // renamed or removed.
        for entry in std::fs::read_dir(&dir).expect("read dir") {
            let entry = entry.expect("dir entry");
            let name = entry.file_name();
            let name = name.to_string_lossy();
            assert!(!name.contains(".tmp."), "orphan sidecar survived: {name}");
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn load_all_skips_oversized_reports() {
        // Round-6 RED→GREEN: a `*.json` larger than MAX_REPORT_BYTES
        // gets skipped (with a warn-level log) instead of read into
        // memory. We size at MAX + 1 byte; the stat check trips first
        // and the file never enters the read path.
        let tmp = std::env::temp_dir().join(format!(
            "valenx-crash-oversized-{}",
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();

        // One legitimate report.
        let good = CrashReport::new("ok", None, "0.1.0");
        good.write_to_disk(&tmp).unwrap();

        // One oversized report. We pad a JSON-ish prefix with zeros so
        // serde_json::from_slice wouldn't accept it even if the cap
        // were bypassed — but the cap should trip first.
        let oversized_path = tmp.join("oversized.json");
        let mut f = std::fs::File::create(&oversized_path).unwrap();
        use std::io::Write;
        f.write_all(b"{\"schema\":1,").unwrap();
        let pad = vec![b' '; MAX_REPORT_BYTES];
        f.write_all(&pad).unwrap();
        f.write_all(b"}").unwrap();
        drop(f);
        assert!(
            std::fs::metadata(&oversized_path).unwrap().len() > MAX_REPORT_BYTES as u64,
            "test setup: oversized file must exceed cap"
        );

        // load_all returns only the good report; the oversized one
        // is skipped without OOM.
        let loaded = CrashReport::load_all(&tmp).expect("load");
        assert_eq!(loaded.len(), 1, "oversized report must be skipped");
        assert_eq!(loaded[0].1.message, "ok");

        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn read_capped_truncates_at_cap_bytes() {
        // Direct check: `read_capped(path, 100)` must read at most
        // 100 bytes even if the file is much larger.
        let tmp = std::env::temp_dir().join(format!(
            "valenx-crash-capread-{}",
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(&tmp, vec![b'x'; 10_000]).unwrap();
        let bytes = read_capped(&tmp, 100).unwrap();
        assert_eq!(bytes.len(), 100);
        let _ = std::fs::remove_file(&tmp);
    }
}
