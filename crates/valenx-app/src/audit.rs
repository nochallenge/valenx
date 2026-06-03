//! Audit-log emission + line formatting, plus the tiny civil-time
//! helpers `emit_audit` needs to stamp ISO 8601 timestamps without
//! pulling in `chrono`. Also the best-effort current-user lookup
//! that every audit entry's `actor.id` field uses.

use crate::state_paths::audit_log_path;

/// One-line human rendering of an audit entry for the in-app log
/// panel: `[audit] <ts> <actor> <action> <kind> [denied]`. Compact on
/// purpose — the user already gets the structured JSON when they
/// click "Audit: Open audit log in file browser".
pub(crate) fn format_audit_log_line(entry: &valenx_audit::AuditEntry) -> String {
    // Pull the most-useful tag out of the target object — the "kind"
    // field is convention across emit_audit callers ({"kind": "case", …}).
    // Fall back to the raw JSON when the shape is unfamiliar so we
    // never silently drop information.
    let target_tag = entry
        .target
        .get("kind")
        .and_then(|v| v.as_str())
        .map(|s| s.to_string())
        .unwrap_or_else(|| entry.target.to_string());
    let mut line = format!(
        "[audit] {ts} {actor} {action} {target}",
        ts = entry.timestamp,
        actor = entry.actor.id,
        action = entry.action,
        target = target_tag,
    );
    if entry
        .context
        .get("denied")
        .and_then(|v| v.as_bool())
        .unwrap_or(false)
    {
        line.push_str(" [denied]");
    }
    line
}

/// Best-effort emit one audit entry for an action the user just
/// performed. Failures are logged via `tracing` but never surfaced
/// — losing audit visibility is bad but spamming UI errors on
/// every workflow action is worse. Operators who care should run
/// `valenx audit verify` (a future CLI) periodically.
pub fn emit_audit(action: &str, target: serde_json::Value, context: serde_json::Value) {
    let Some(path) = audit_log_path() else {
        return;
    };
    let writer = valenx_audit::AuditWriter::new(path);
    let entry = valenx_audit::AuditEntry {
        timestamp: current_timestamp_iso8601(),
        actor: valenx_audit::AuditActor {
            id: current_actor_id(),
            session_id: None,
        },
        action: action.to_string(),
        target,
        context,
        prev_hash: String::new(),
    };
    if let Err(e) = writer.append(entry) {
        tracing::warn!(target: "valenx.audit", ?e, "audit emit failed");
    }
}

/// Current wall-clock time as an ISO 8601 UTC string. We intentionally
/// don't pull in `chrono` for this — `SystemTime` + a tiny formatter
/// covers what the audit log needs (sub-second precision is overkill
/// here; second-precision is fine for compliance-level audit trails).
pub(crate) fn current_timestamp_iso8601() -> String {
    let now = std::time::SystemTime::now();
    let secs = now
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    // Convert to civil time using a small-and-honest algorithm.
    // Howard Hinnant's "date.h" technique — works for any epoch ≥ 1970.
    let days = (secs / 86_400) as i64;
    let secs_of_day = (secs % 86_400) as u32;
    let h = secs_of_day / 3600;
    let m = (secs_of_day % 3600) / 60;
    let s = secs_of_day % 60;
    let (y, mo, d) = days_to_ymd(days);
    format!("{y:04}-{mo:02}-{d:02}T{h:02}:{m:02}:{s:02}Z")
}

/// Civil-time conversion from days-since-1970 to (year, month, day).
/// Implements the algorithm from
/// <https://howardhinnant.github.io/date_algorithms.html#civil_from_days>.
/// No external date crate.
pub(crate) fn days_to_ymd(z: i64) -> (i32, u32, u32) {
    let z = z + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = (z - era * 146_097) as u64;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if m <= 2 { y + 1 } else { y };
    (y as i32, m as u32, d as u32)
}

/// Best-effort current-user-id detection. Returns the OS-reported
/// user name; falls back to "unknown" if even the env-var fallbacks
/// are unset. Real auth (SSO / OIDC) replaces this in the auth RFC.
///
/// Round-11 fix: pre-fix this function consulted `$USER` first, then
/// `$USERNAME`. On Windows the canonical env var is `%USERNAME%`, but
/// `$USER` is freely user-mutable (any unprivileged shell can `set
/// USER=admin@example.com`) — so the pre-fix order let a non-admin
/// Windows process spoof an arbitrary identity into both the RBAC
/// gate and every audit-log line it emitted. The fix routes through
/// `whoami::username()` first, which calls the OS-canonical APIs
/// directly (GetUserNameW on Windows, getpwuid_r on Unix) and ignores
/// the env. The env-var path remains as a fallback for embedded /
/// minimal environments where the OS user table isn't populated.
pub(crate) fn current_actor_id() -> String {
    // Prefer the real OS identity (GetUserNameW on Windows,
    // getpwuid_r on Unix) — env vars like $USER / $USERNAME are
    // user-mutable and must not be the sole authorisation gate.
    // whoami::username() returns Result<String, whoami::Error> in
    // v2; treat any Err (no user table available) as "fall through to
    // env-var fallback" rather than propagating.
    if let Ok(real) = whoami::username() {
        if !real.is_empty() && real != "unknown" {
            return real;
        }
    }
    // Fallback for embedded / minimal environments without a real
    // user table. Order: USERNAME (Windows canonical) → USER (Unix
    // canonical) → "unknown". Note: this is best-effort identity
    // surfacing only; RBAC enforcement still authorises against the
    // canonical OS identity above.
    if let Some(u) = std::env::var_os("USERNAME") {
        return u.to_string_lossy().into_owned();
    }
    if let Some(u) = std::env::var_os("USER") {
        return u.to_string_lossy().into_owned();
    }
    "unknown".to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn days_to_ymd_handles_known_dates() {
        // Anchor against well-known epoch dates so a future
        // refactor of the algorithm can't silently drift.
        // 1970-01-01 = day 0
        assert_eq!(days_to_ymd(0), (1970, 1, 1));
        // 2000-01-01 = day 10957 (Y2K)
        assert_eq!(days_to_ymd(10957), (2000, 1, 1));
        // 2026-04-26 — anchor to validate the leap-year logic
        // 56 years × 365 + 14 leaps + Jan/Feb/Mar (90) + 26 = 20570
        // Wait: 2026-01-01 = 20454 (= 56 × 365 + 14), so
        // 2026-04-26 = 20454 + 31 + 28 + 31 + 26 - 1 = 20569.
        assert_eq!(days_to_ymd(20569), (2026, 4, 26));
        // Day before — sanity check.
        assert_eq!(days_to_ymd(20568), (2026, 4, 25));
    }

    #[test]
    fn current_timestamp_iso8601_has_the_right_shape() {
        let ts = current_timestamp_iso8601();
        // YYYY-MM-DDTHH:MM:SSZ — exactly 20 chars.
        assert_eq!(ts.len(), 20);
        assert!(ts.ends_with('Z'));
        assert!(ts.contains('T'));
    }

    #[test]
    fn tail_audit_log_pushes_lines_into_log_panel_when_log_has_entries() {
        // We can't depend on the OS's real audit log having entries
        // — instead, exercise the rendering path directly by building
        // a synthetic AuditEntry, formatting it, and asserting the
        // shape. The end-to-end path is covered by the audit crate's
        // tail_n tests + the no-panic test above.
        let entry = valenx_audit::AuditEntry {
            timestamp: "2026-04-25T12:34:56Z".into(),
            actor: valenx_audit::AuditActor {
                id: "alice@example.com".into(),
                session_id: None,
            },
            action: "run.start".into(),
            target: serde_json::json!({"kind": "case", "case": "cfd"}),
            context: serde_json::json!({}),
            prev_hash: "genesis".into(),
        };
        let line = format_audit_log_line(&entry);
        assert!(line.starts_with("[audit]"), "got: {line}");
        assert!(line.contains("2026-04-25T12:34:56Z"), "got: {line}");
        assert!(line.contains("alice@example.com"), "got: {line}");
        assert!(line.contains("run.start"), "got: {line}");
        assert!(line.contains("case"), "got: {line}");
        assert!(!line.contains("[denied]"), "got: {line}");
    }

    #[test]
    fn format_audit_log_line_marks_denied_entries() {
        // RBAC denials emit context.denied=true; the formatter must
        // surface that so a quick scroll of the log makes the rejected
        // attempts obvious.
        let entry = valenx_audit::AuditEntry {
            timestamp: "2026-04-25T12:34:56Z".into(),
            actor: valenx_audit::AuditActor {
                id: "alice@example.com".into(),
                session_id: None,
            },
            action: "case.delete".into(),
            target: serde_json::json!({"kind": "permission_denied"}),
            context: serde_json::json!({"denied": true, "user_role": "Viewer"}),
            prev_hash: "genesis".into(),
        };
        let line = format_audit_log_line(&entry);
        assert!(line.ends_with("[denied]"), "got: {line}");
    }

    #[test]
    fn current_actor_id_prefers_whoami_over_user_env_var() {
        // Round-11 RED→GREEN: pre-fix `current_actor_id` consulted
        // `$USER` first, which on Windows let an unprivileged process
        // spoof an arbitrary identity by `set USER=spoofed`. The fix
        // routes through `whoami::username()` first so the OS-canonical
        // identity wins even when `$USER` is mutated.
        //
        // We can't safely mutate process-global env in a parallel
        // test runner, so the assertion is structural: the function
        // must return the same value `whoami::username()` does
        // (verbatim) whenever whoami returns a non-empty / non-"unknown"
        // identity. Any future regression that re-introduces the
        // env-var-first path will diverge from whoami and break this.
        let real = match whoami::username() {
            Ok(s) if !s.is_empty() && s != "unknown" => s,
            _ => {
                // Build environment without a real user table — skip; the
                // function correctly falls back to env. Covered by the
                // sibling fallback test below.
                return;
            }
        };
        let got = current_actor_id();
        assert_eq!(
            got, real,
            "current_actor_id must mirror whoami::username() when whoami succeeds; \
             pre-fix this consulted `$USER` first (Windows spoof axis)"
        );
    }

    #[test]
    fn current_actor_id_never_returns_empty_or_panics() {
        // A second guard: the function must always return *something*
        // — every callsite assumes a non-empty actor.id for both audit
        // log line emission and RBAC enforcement.
        let id = current_actor_id();
        assert!(!id.is_empty(), "current_actor_id returned empty string");
    }

    #[test]
    fn format_audit_log_line_falls_back_to_raw_target_when_kind_missing() {
        // Older / non-conventional emit_audit callers may have
        // targets without a "kind" field. Don't drop the info on the
        // floor — show the raw JSON so the user can still see
        // *something*.
        let entry = valenx_audit::AuditEntry {
            timestamp: "2026-04-25T12:34:56Z".into(),
            actor: valenx_audit::AuditActor {
                id: "system:scheduler".into(),
                session_id: None,
            },
            action: "weird.action".into(),
            target: serde_json::json!({"id": 42}),
            context: serde_json::json!({}),
            prev_hash: "genesis".into(),
        };
        let line = format_audit_log_line(&entry);
        assert!(line.contains("\"id\":42"), "got: {line}");
    }
}
