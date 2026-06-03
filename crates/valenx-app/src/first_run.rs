//! First-launch wizard — egui shell that consumes the state
//! machine in `valenx-first-run`.
//!
//! The wizard fires on the first run after install. It probes
//! every registered adapter, renders a per-adapter row with the
//! result, and lets the user dismiss the dialog ("Done"). The
//! resulting [`FirstRunDecision`] persists at
//! `<state_dir>/first-run.json` so subsequent launches skip the
//! wizard.
//!
//! ## Why minimal
//!
//! v0.1.0 targets the "what's installed?" disclosure. Per-adapter
//! enable/disable checkboxes are a v0.2.0 follow-up — the data
//! model in `valenx-first-run::FirstRunDecision::adapter_enabled`
//! already supports them; the UI just doesn't render the controls
//! yet.
//!
//! ## Loading / saving
//!
//! Persistence routes through [`load_first_run_from_state_dir`] /
//! [`save_first_run_to_state_dir`] mirroring how `Settings` works,
//! so users get the same per-platform path resolution and the
//! same "fall back to no-op when state_dir is None" semantics.

use std::path::PathBuf;

use valenx_first_run::{
    install_hint_for, AdapterAvailability, AdapterStatus, EnvironmentReport, FirstRunDecision,
};

use crate::state_paths::atomic_write;

/// Visible-row filter for the wizard's adapter table.
///
/// Persisted across frames via `egui::Memory::data_mut()` so toggling
/// the filter mid-session survives the wizard repaint loop. Defaults
/// to `All` on first open.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum WizardFilter {
    #[default]
    All,
    Installed,
    Missing,
}

impl WizardFilter {
    fn matches(&self, a: &AdapterStatus) -> bool {
        match self {
            Self::All => true,
            Self::Installed => matches!(a.availability, AdapterAvailability::Installed),
            Self::Missing => matches!(a.availability, AdapterAvailability::Missing),
        }
    }
}

/// Build the environment report from the current adapter registry.
///
/// `probe_results` is a vector of `(id, display_name, availability,
/// detected_version)` — the GUI's existing adapter probe loop already
/// produces this shape, so wiring is a one-liner.
pub fn build_report(
    probe_results: &[(String, String, AdapterAvailability, Option<String>)],
) -> EnvironmentReport {
    let adapters: Vec<AdapterStatus> = probe_results
        .iter()
        .map(|(id, display, availability, version)| {
            let install_hint = if availability.is_ready() {
                None
            } else {
                install_hint_for(id)
            };
            AdapterStatus {
                id: id.clone(),
                display_name: display.clone(),
                availability: *availability,
                detected_version: version.clone(),
                install_hint,
            }
        })
        .collect();
    EnvironmentReport { adapters }
}

/// Render the wizard. Returns `true` when the user dismisses (clicks
/// Done or Skip), so the caller can flip the in-memory decision's
/// `completed` flag and trigger persistence.
///
/// Pure rendering — does no I/O. Persistence is the caller's job.
///
/// `catalogue` supplies localised strings. The per-adapter status
/// labels still come from `availability.label()` (English-only at
/// v0.1.0); a follow-up routes those through the catalogue too.
pub fn show(
    ctx: &egui::Context,
    open: &mut bool,
    report: &EnvironmentReport,
    decision: &mut FirstRunDecision,
    catalogue: &valenx_i18n::LocaleCatalogue,
) -> bool {
    let mut dismissed = false;
    let title = catalogue.lookup("dialog.first-run.title-window");
    // Filter state persists across frames via egui's per-id memory
    // store. Keyed off a stable id so re-opening the wizard restores
    // the last filter pick within the session. Reset to `All` at
    // process exit (memory isn't serialised here).
    let filter_id = egui::Id::new("valenx_first_run_filter");
    let mut filter: WizardFilter = ctx.data_mut(|d| d.get_temp(filter_id).unwrap_or_default());
    egui::Window::new(title)
        .open(open)
        .collapsible(false)
        .resizable(true)
        .default_size([880.0, 640.0])
        .min_width(720.0)
        .show(ctx, |ui| {
            // Headline summary — always reflects totals, never the
            // current filter, so the user always sees the "Detected
            // N of M adapters" headline truthfully.
            ui.heading(catalogue.lookup("dialog.first-run.title"));
            ui.add_space(4.0);
            ui.label(catalogue.format_with(
                "dialog.first-run.summary-ready",
                &[
                    ("count", &report.ready_count().to_string()),
                    ("total", &report.adapters.len().to_string()),
                ],
            ));
            ui.add_space(6.0);

            // Filter selector — three radio-style choices in a row.
            ui.horizontal(|ui| {
                ui.label("Show:");
                ui.selectable_value(&mut filter, WizardFilter::All, "All");
                ui.selectable_value(&mut filter, WizardFilter::Installed, "Installed");
                ui.selectable_value(&mut filter, WizardFilter::Missing, "Missing");
            });
            ui.add_space(6.0);

            // Per-adapter table. Filter is applied here; totals above
            // are untouched.
            egui::ScrollArea::vertical()
                .max_height(480.0)
                .show(ui, |ui| {
                    egui::Grid::new("first-run-adapter-table")
                        .num_columns(3)
                        .spacing([16.0, 4.0])
                        .striped(true)
                        .show(ui, |ui| {
                            ui.strong(catalogue.lookup("dialog.first-run.col-adapter"));
                            ui.strong(catalogue.lookup("dialog.first-run.col-status"));
                            ui.strong(catalogue.lookup("dialog.first-run.col-hint"));
                            ui.end_row();
                            for a in report.adapters.iter().filter(|a| filter.matches(a)) {
                                ui.label(&a.display_name);
                                ui.label(a.availability.label());
                                if let Some(hint) = &a.install_hint {
                                    ui.label(hint);
                                } else {
                                    ui.label("—");
                                }
                                ui.end_row();
                            }
                        });
                });
            ui.add_space(8.0);

            ui.horizontal(|ui| {
                if ui
                    .button(catalogue.lookup("dialog.first-run.done"))
                    .clicked()
                {
                    decision.mark_completed();
                    dismissed = true;
                }
                if ui
                    .button(catalogue.lookup("dialog.first-run.skip"))
                    .clicked()
                {
                    // Skip still marks completed — the wizard
                    // doesn't auto-reopen. Users can re-launch it
                    // from the command palette.
                    decision.mark_completed();
                    dismissed = true;
                }
            });
        });
    // Persist the (possibly mutated) filter back into egui memory so
    // the next frame's read picks up the user's selection.
    ctx.data_mut(|d| d.insert_temp(filter_id, filter));
    dismissed
}

/// Path the [`FirstRunDecision`] persists at.
///
/// Mirrors `settings_path()` in lib.rs — same per-platform state
/// dir, same `Option` semantics so a HOME-less host degrades to
/// "decision lives only in memory" instead of crashing.
pub fn first_run_path() -> Option<PathBuf> {
    crate::state_dir().map(|d| d.join("first-run.json"))
}

/// Load the persisted decision, or `None` on first launch / when
/// the file is unreadable.
///
/// Round-8 hardening: bounded read via
/// [`crate::settings_io::MAX_STATE_FILE_BYTES`] — sister-loader to
/// settings/run-history/sweep-history. A first-run JSON past the
/// cap is rejected as "no decision" so a hostile / corrupted state
/// file can't OOM the host on startup.
pub fn load_first_run_from_state_dir() -> Option<FirstRunDecision> {
    use std::io::Read;
    let path = first_run_path()?;
    if !path.is_file() {
        return None;
    }
    let meta = std::fs::metadata(&path).ok()?;
    if meta.len() > crate::settings_io::MAX_STATE_FILE_BYTES as u64 {
        return None;
    }
    let mut buf = Vec::with_capacity(
        meta.len()
            .min(crate::settings_io::MAX_STATE_FILE_BYTES as u64) as usize,
    );
    let mut file = std::fs::File::open(&path).ok()?;
    file.by_ref()
        .take(crate::settings_io::MAX_STATE_FILE_BYTES as u64)
        .read_to_end(&mut buf)
        .ok()?;
    serde_json::from_slice(&buf).ok()
}

/// Persist the decision. Best-effort — write failures log but
/// don't bubble up so the GUI's main loop stays smooth. Uses
/// [`crate::state_paths::atomic_write`] (write to `<path>.tmp` then
/// rename) so a crash mid-write keeps the previous decision file
/// rather than producing a partial / empty `first-run.json`.
pub fn save_first_run_to_state_dir(decision: &FirstRunDecision) {
    let Some(path) = first_run_path() else {
        return;
    };
    match serde_json::to_string_pretty(decision) {
        Ok(text) => {
            if let Err(e) = atomic_write(&path, &text) {
                tracing::warn!(target: "valenx-first-run", ?e, ?path, "write decision");
            }
        }
        Err(e) => {
            tracing::warn!(target: "valenx-first-run", ?e, "serialise decision");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn build_report_attaches_install_hints_only_to_unready_adapters() {
        let probes = vec![
            (
                "openfoam".to_string(),
                "OpenFOAM".to_string(),
                AdapterAvailability::Installed,
                Some("v2406".to_string()),
            ),
            (
                "gmsh".to_string(),
                "gmsh".to_string(),
                AdapterAvailability::Missing,
                None,
            ),
        ];
        let report = build_report(&probes);
        assert_eq!(report.adapters.len(), 2);
        // Installed → no hint surfaced.
        assert!(report.adapters[0].install_hint.is_none());
        // Missing → hint pulled from the per-OS catalogue.
        assert!(report.adapters[1].install_hint.is_some());
    }

    #[test]
    fn build_report_preserves_input_order() {
        let probes = vec![
            (
                "z".to_string(),
                "Z".into(),
                AdapterAvailability::Missing,
                None,
            ),
            (
                "a".to_string(),
                "A".into(),
                AdapterAvailability::Missing,
                None,
            ),
            (
                "m".to_string(),
                "M".into(),
                AdapterAvailability::Missing,
                None,
            ),
        ];
        let report = build_report(&probes);
        let ids: Vec<&str> = report.adapters.iter().map(|a| a.id.as_str()).collect();
        assert_eq!(ids, vec!["z", "a", "m"]);
    }

    #[test]
    fn wizard_filter_default_is_all() {
        assert_eq!(WizardFilter::default(), WizardFilter::All);
    }

    #[test]
    fn wizard_filter_matches_by_availability() {
        let installed = AdapterStatus {
            id: "x".into(),
            display_name: "X".into(),
            availability: AdapterAvailability::Installed,
            detected_version: None,
            install_hint: None,
        };
        let missing = AdapterStatus {
            id: "y".into(),
            display_name: "Y".into(),
            availability: AdapterAvailability::Missing,
            detected_version: None,
            install_hint: None,
        };
        // All accepts everything.
        assert!(WizardFilter::All.matches(&installed));
        assert!(WizardFilter::All.matches(&missing));
        // Installed accepts only Installed.
        assert!(WizardFilter::Installed.matches(&installed));
        assert!(!WizardFilter::Installed.matches(&missing));
        // Missing accepts only Missing.
        assert!(!WizardFilter::Missing.matches(&installed));
        assert!(WizardFilter::Missing.matches(&missing));
    }

    #[test]
    fn first_run_path_is_state_dir_relative() {
        // We can't assert the exact path (state_dir() varies by
        // OS + env), but if HOME / APPDATA / XDG_STATE_HOME is set
        // the path's filename is canonical.
        if let Some(p) = first_run_path() {
            assert_eq!(
                p.file_name().and_then(|s| s.to_str()),
                Some("first-run.json")
            );
        }
    }
}
