//! Shared chrome for the right-side workbench panels.
//!
//! Every workbench draws its own [`egui::SidePanel`] with the same header
//! shape: a title, an optional weak subtitle, and a separator. This module
//! factors that out and adds a right-aligned "✕" close button so any
//! workbench can be dismissed without hunting through the View menu.
//!
//! The helper does not own the visibility state — it just reports a click.
//! The caller clears its own `show_<name>_workbench` flag when
//! [`header`] returns `true`:
//!
//! ```ignore
//! if crate::workbench_ui::header(ui, "Springs", "native helical-spring design") {
//!     app.show_springs_workbench = false;
//! }
//! ```

use eframe::egui;

/// Draw a workbench panel header: the `title` heading on the left and a
/// "✕" close button pinned to the right, an optional `subtitle` line
/// beneath (skipped when empty), then a separator.
///
/// Returns `true` on the frame the close button is clicked, so the caller
/// can hide the workbench by clearing its visibility flag.
pub fn header(ui: &mut egui::Ui, title: &str, subtitle: &str) -> bool {
    let mut close = false;
    ui.horizontal(|ui| {
        ui.heading(title);
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            if ui
                .button("✕")
                .on_hover_text("Close this workbench")
                .clicked()
            {
                close = true;
            }
        });
    });
    if !subtitle.is_empty() {
        ui.label(egui::RichText::new(subtitle).weak().small());
    }
    ui.separator();
    close
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn header_draws_headless_and_stays_open_without_a_click() {
        let ctx = egui::Context::default();
        let mut closed = true;
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            egui::CentralPanel::default().show(ctx, |ui| {
                // No pointer input is synthesised, so the close button is
                // never clicked and the workbench stays open.
                closed = header(ui, "Test Workbench", "a subtitle");
            });
        });
        assert!(!closed, "header reports close only when ✕ is clicked");
    }

    #[test]
    fn header_with_empty_subtitle_does_not_panic() {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            egui::CentralPanel::default().show(ctx, |ui| {
                let _ = header(ui, "No Subtitle", "");
            });
        });
    }
}
