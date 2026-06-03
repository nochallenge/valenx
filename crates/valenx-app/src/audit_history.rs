//! Audit-log + history + report-export methods on [`ValenxApp`].
//! Split out of `lib.rs` as part of the structural refactor.

use valenx_core::LogLevel;

use crate::audit::{emit_audit, format_audit_log_line};
use crate::file_browser::{open_path_or_copy, POPUP_DISABLED_PREFIX};
use crate::rbac_io::rbac_check;
use crate::residuals;
use crate::state_paths::{audit_log_path, export_csv_path, run_history_path, sweep_history_path};
use crate::types::BottomTab;
use crate::ValenxApp;

impl ValenxApp {
    /// Walk the audit log file end-to-end, verifying every entry's
    /// `prev_hash` chains correctly to the previous line. Surfaces the
    /// result in `status` (success path) or `last_error` (broken
    /// chain / unreadable file). RFC 0013 §15.2.
    ///
    /// Compliance use-case: an auditor opens the app on a frozen
    /// machine image and clicks "Verify audit log integrity"; the
    /// status bar reports `"verified N entries from <first_ts> to
    /// <last_ts>, head <hash>"` or names the first broken line.
    ///
    /// No RBAC gate — verification is read-only and any user with
    /// access to the log file should be able to check it. The result
    /// is itself audit-emitted so the verification attempt becomes
    /// part of the chain.
    pub fn verify_audit_log(&mut self) {
        let Some(path) = audit_log_path() else {
            self.last_error = Some("No state directory — audit log location unknown.".into());
            return;
        };
        if !path.exists() {
            self.status = Some(format!(
                "Audit log empty (no file at {}) — nothing to verify yet.",
                path.display()
            ));
            self.last_error = None;
            return;
        }
        match valenx_audit::verify_chain_report(&path) {
            Ok(report) => {
                let head = report
                    .head_hash
                    .as_deref()
                    .map(|h| &h[..8.min(h.len())])
                    .unwrap_or("(empty)");
                self.status = Some(format!(
                    "Audit log verified: {} entries from {} to {} (head {}…)",
                    report.entries_verified,
                    report.first_timestamp.as_deref().unwrap_or("(none)"),
                    report.last_timestamp.as_deref().unwrap_or("(none)"),
                    head,
                ));
                self.last_error = None;
                emit_audit(
                    "audit.verify",
                    serde_json::json!({"kind": "audit_log"}),
                    serde_json::json!({
                        "entries_verified": report.entries_verified,
                        "head_hash": report.head_hash,
                    }),
                );
            }
            Err(e) => {
                self.last_error = Some(format!("Audit log verification failed: {e}"));
                emit_audit(
                    "audit.verify",
                    serde_json::json!({"kind": "audit_log"}),
                    serde_json::json!({
                        "denied": true,
                        "error": e.to_string(),
                    }),
                );
            }
        }
    }

    /// Clear the in-memory run-history map AND the persisted
    /// `<state_dir>/run-history.json` file. After this call the
    /// case browser's ✓ / ✗ badges go away and previously-converged
    /// cases look "never run" until the user re-runs them.
    ///
    /// RBAC-gated on PrepareCase (any Runner); the action is
    /// reversible (just re-run the cases) so it doesn't need
    /// Admin-level gating. Audit-emits "history.clear" so the
    /// clear is itself part of the chain.
    pub fn clear_run_history(&mut self) {
        if let Err(e) = rbac_check(
            valenx_rbac::Action::PrepareCase,
            self.project_rbac_override.as_ref(),
        ) {
            self.last_error = Some(format!("{e}"));
            return;
        }
        let n = self.run_history.len();
        self.run_history.clear();
        // Delete the on-disk file so the cleared state survives a
        // restart. Empty-map serialise + write would also work but
        // would leave a `{}` file on disk.
        if let Some(path) = run_history_path() {
            let _ = std::fs::remove_file(&path);
        }
        self.status = Some(format!("Cleared run history ({n} entries)"));
        self.last_error = None;
        emit_audit(
            "history.clear",
            serde_json::json!({"kind": "run_history"}),
            serde_json::json!({"entries_removed": n}),
        );
    }

    /// Clear the in-memory sweep-history map AND the persisted
    /// `<state_dir>/sweep-history.json` file. Same RBAC + audit
    /// shape as [`Self::clear_run_history`].
    pub fn clear_sweep_history(&mut self) {
        if let Err(e) = rbac_check(
            valenx_rbac::Action::PrepareCase,
            self.project_rbac_override.as_ref(),
        ) {
            self.last_error = Some(format!("{e}"));
            return;
        }
        let n = self.sweep_history.len();
        self.sweep_history.clear();
        if let Some(path) = sweep_history_path() {
            let _ = std::fs::remove_file(&path);
        }
        self.status = Some(format!("Cleared sweep history ({n} entries)"));
        self.last_error = None;
        emit_audit(
            "history.clear",
            serde_json::json!({"kind": "sweep_history"}),
            serde_json::json!({"entries_removed": n}),
        );
    }

    /// Render the most recent run's `Results` bundle as a self-
    /// contained HTML report and drop it under
    /// `<state_dir>/exports/report-<unix>.html`. No external CSS /
    /// JS / images — the file works offline + can be emailed
    /// directly.
    ///
    /// No RBAC gate — exporting data the user already sees in the
    /// Results pane isn't elevated. Audit-emits "report.export" so
    /// the export becomes part of the chain.
    pub fn export_html_report(&mut self) {
        let Some(results) = self.last_run_results.as_deref() else {
            self.last_error =
                Some("No run results yet — finish a run before exporting a report.".into());
            return;
        };
        let target = export_csv_path("report").with_extension("html");
        match valenx_export::write_html_report(results, &target) {
            Ok(()) => {
                self.status = Some(format!("HTML report written: {}", target.display()));
                self.last_error = None;
                emit_audit(
                    "report.export",
                    serde_json::json!({"kind": "html_report"}),
                    serde_json::json!({
                        "path": target.display().to_string(),
                        "case_id": results.meta.case_id,
                    }),
                );
            }
            Err(e) => {
                self.last_error = Some(format!(
                    "Couldn't write HTML report to {}: {e}",
                    target.display()
                ));
            }
        }
    }

    /// Export the current run's residual time-series as a long-format
    /// CSV file. Lands in `<state_dir>/exports/residuals-<ts>.csv`
    /// (or the system temp dir if no state dir resolves) so the user
    /// can pull it into pandas / Excel without screen-scraping the
    /// chart.
    ///
    /// No RBAC gate — exporting data the user already sees in the
    /// chart isn't elevated. Audit-emits "residuals.export" so the
    /// export is itself part of the chain.
    pub fn export_residuals_csv(&mut self) {
        if self.residuals.is_empty() {
            self.last_error = Some("No residuals to export — run a case first.".into());
            return;
        }
        let target = export_csv_path("residuals");
        match residuals::write_residuals_csv(&self.residuals, &target) {
            Ok(()) => {
                self.status = Some(format!("Residuals exported: {}", target.display()));
                self.last_error = None;
                emit_audit(
                    "residuals.export",
                    serde_json::json!({"kind": "residuals"}),
                    serde_json::json!({
                        "path": target.display().to_string(),
                        "fields": self
                            .residuals
                            .by_field
                            .keys()
                            .copied()
                            .collect::<Vec<_>>(),
                    }),
                );
            }
            Err(e) => {
                self.last_error = Some(format!(
                    "Couldn't write residual CSV to {}: {e}",
                    target.display()
                ));
            }
        }
    }

    /// Open the audit log file in the host's default JSON / text
    /// viewer (or surface its parent directory in the file browser
    /// if no JSONL viewer is registered). Useful for compliance
    /// auditors who want to scroll the raw chain entries directly
    /// rather than going through the GUI's filtered view.
    ///
    /// No RBAC gate — opening a file the user already has read
    /// access to via the OS isn't elevated. Audit-emits the
    /// "audit.open" action so the open is itself part of the
    /// chain.
    pub fn open_audit_log(&mut self) {
        let Some(path) = audit_log_path() else {
            self.last_error = Some("No state directory — audit log location unknown.".into());
            return;
        };
        // Prefer to surface the file directly. If it doesn't exist
        // yet (fresh install, no actions performed) fall back to
        // the parent directory so the user lands somewhere useful
        // instead of getting an "open path failed" error.
        let target = if path.is_file() {
            path.clone()
        } else if let Some(parent) = path.parent() {
            parent.to_path_buf()
        } else {
            path.clone()
        };
        let disable = self.settings.disable_file_browser_popups;
        match open_path_or_copy(&target, disable) {
            Ok(()) => {
                self.status = Some(format!("Opened audit log location: {}", target.display()));
                self.last_error = None;
                emit_audit(
                    "audit.open",
                    serde_json::json!({"kind": "audit_log"}),
                    serde_json::json!({
                        "path": target.display().to_string(),
                    }),
                );
            }
            Err(reason) if reason.starts_with(POPUP_DISABLED_PREFIX) => {
                // Kill-switch path: surface as a neutral status line
                // and still emit the audit entry — the user *did*
                // ask to know the audit log location, the
                // distinction is just delivery mechanism.
                self.status = Some(reason);
                self.last_error = None;
                emit_audit(
                    "audit.open",
                    serde_json::json!({"kind": "audit_log"}),
                    serde_json::json!({
                        "path": target.display().to_string(),
                    }),
                );
            }
            Err(reason) => {
                self.last_error = Some(format!(
                    "Couldn't open audit log: {reason}. Path: {}",
                    target.display()
                ));
            }
        }
    }

    /// Read the last `n` audit-log entries and surface them in the
    /// bottom dock's Log tab so the user can review recent activity
    /// without leaving the GUI.
    ///
    /// Each entry becomes one log line of shape
    /// `[audit] <timestamp> <actor> <action> <target_kind>` plus an
    /// optional `denied=true` suffix if the entry was an RBAC denial.
    /// Switches the bottom tab to Log so the lines are immediately
    /// visible.
    ///
    /// No RBAC gate — viewing data the user already has read access
    /// to via the on-disk file isn't elevated. The view is itself
    /// audit-emitted so the inspection becomes part of the chain
    /// (compliance use-case: detect `"operator scrolled the audit log
    /// at <ts>"` patterns).
    pub fn tail_audit_log(&mut self, n: usize) {
        let Some(path) = audit_log_path() else {
            self.last_error = Some("No state directory — audit log location unknown.".into());
            return;
        };
        match valenx_audit::tail_n(&path, n) {
            Ok(entries) if entries.is_empty() => {
                self.status = Some(format!(
                    "Audit log empty (no entries at {}).",
                    path.display()
                ));
                self.last_error = None;
            }
            Ok(entries) => {
                let count = entries.len();
                self.log.push(
                    LogLevel::Info,
                    format!(
                        "[audit] --- last {count} {} ---",
                        if count == 1 { "entry" } else { "entries" },
                    ),
                );
                for entry in &entries {
                    self.log.push(LogLevel::Info, format_audit_log_line(entry));
                }
                self.bottom_tab = BottomTab::Log;
                self.status = Some(format!(
                    "Showed last {count} audit {} in Log tab.",
                    if count == 1 { "entry" } else { "entries" },
                ));
                self.last_error = None;
                emit_audit(
                    "audit.tail",
                    serde_json::json!({"kind": "audit_log"}),
                    serde_json::json!({"requested": n, "shown": count}),
                );
            }
            Err(e) => {
                self.last_error = Some(format!("Audit log tail failed: {e}"));
                emit_audit(
                    "audit.tail",
                    serde_json::json!({"kind": "audit_log"}),
                    serde_json::json!({
                        "denied": true,
                        "requested": n,
                        "error": e.to_string(),
                    }),
                );
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::types::RunHistoryEntry;
    use crate::ValenxApp;

    #[test]
    fn export_html_report_without_run_errors_cleanly() {
        let mut app = ValenxApp::default();
        app.export_html_report();
        assert!(
            app.last_error
                .as_ref()
                .is_some_and(|e| e.contains("No run results")),
            "got {:?}",
            app.last_error
        );
    }

    #[test]
    fn clear_run_history_empties_in_memory_map() {
        let mut app = ValenxApp::default();
        app.run_history.insert(
            "case-a".into(),
            RunHistoryEntry {
                succeeded: true,
                wall_time: std::time::Duration::from_secs(1),
                converged: Some(true),
            },
        );
        assert_eq!(app.run_history.len(), 1);
        app.clear_run_history();
        assert!(app.run_history.is_empty());
        assert!(
            app.status.as_ref().is_some_and(|s| s.contains("Cleared")),
            "got status: {:?}",
            app.status
        );
    }

    #[test]
    fn clear_run_history_on_empty_app_is_a_no_op() {
        let mut app = ValenxApp::default();
        app.clear_run_history();
        assert!(app.run_history.is_empty());
        assert!(app.last_error.is_none());
    }

    #[test]
    fn export_residuals_csv_without_data_errors_cleanly() {
        // No run yet -> nothing to export -> friendly error
        // surfaces in last_error (no panic, no half-written file).
        let mut app = ValenxApp::default();
        app.export_residuals_csv();
        assert!(
            app.last_error
                .as_ref()
                .is_some_and(|e| e.contains("residuals")),
            "got {:?}",
            app.last_error
        );
    }

    #[test]
    fn open_audit_log_doesnt_panic_when_no_state_dir_or_file() {
        // Mirrors verify_audit_log_handles_missing_log_gracefully:
        // a fresh-default app should produce either a status update
        // or a state-dir error rather than panicking.
        let mut app = ValenxApp::default();
        app.open_audit_log();
        let either = app.status.is_some()
            || app
                .last_error
                .as_ref()
                .is_some_and(|e| e.contains("state directory") || e.contains("Couldn't open"));
        assert!(
            either,
            "expected status or known error; got status={:?} last_error={:?}",
            app.status, app.last_error
        );
    }

    #[test]
    fn verify_audit_log_handles_missing_log_gracefully() {
        // No state dir / no log file = a friendly status message,
        // not an error. Single-user-machine first-run case.
        //
        // CI hosts and dev machines may have a real audit log that
        // happens to be malformed (corrupted entries from prior
        // crashes, partial writes from interrupted runs). In that
        // case the verifier surfaces a structured AuditError that
        // we treat as "expected behaviour" — the goal is to confirm
        // the verifier doesn't panic on real-world input, not to
        // assert the user's machine is in a clean state.
        let mut app = ValenxApp::default();
        app.verify_audit_log();
        let either_outcome_is_handled = app.status.is_some()
            || app.last_error.as_ref().is_some_and(|e| {
                e.contains("state directory") || e.contains("Audit log verification failed")
            });
        assert!(
            either_outcome_is_handled,
            "expected status or state-dir / verification error, got status={:?} last_error={:?}",
            app.status, app.last_error
        );
    }

    #[test]
    fn tail_audit_log_handles_missing_log_gracefully() {
        // No state dir / empty log = a friendly status message, not
        // an error. The CI host's real audit log may have entries —
        // in that case we just want "no panic, status updated".
        let mut app = ValenxApp::default();
        app.tail_audit_log(10);
        let handled = app.status.is_some()
            || app
                .last_error
                .as_ref()
                .is_some_and(|e| e.contains("state directory") || e.contains("Audit log tail"));
        assert!(
            handled,
            "expected status or known error; got status={:?} last_error={:?}",
            app.status, app.last_error
        );
    }

    #[test]
    fn tail_audit_log_zero_request_is_a_no_op() {
        // n=0 should not push anything into the log panel even when
        // the audit log has entries. Mirrors valenx_audit::tail_n's
        // n=0 short-circuit semantics.
        let mut app = ValenxApp::default();
        let lines_before = app.log.lines.len();
        app.tail_audit_log(0);
        assert_eq!(
            app.log.lines.len(),
            lines_before,
            "n=0 should not add any log lines"
        );
    }
}
