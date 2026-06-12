//! Keyboard-shortcut cheat-sheet overlay.
//!
//! A floating egui window mounted by [`render_cheatsheet`] that
//! enumerates [`ShortcutAction::ALL`](crate::shortcuts::ShortcutAction)
//! as a three-column grid: **binding · label · description**.
//!
//! Triggered by the `?` key (handled in `update.rs` via
//! [`crate::shortcuts::ShortcutAction::ToggleKeyboardHelp`]) and by
//! the Help menu entry.

use eframe::egui;

use crate::shortcuts::ShortcutAction;

/// Render the keyboard cheat-sheet window. Closes when the user
/// clicks the `[x]` or presses Esc.
pub fn render_cheatsheet(ctx: &egui::Context, open: &mut bool) {
    egui::Window::new("Keyboard shortcuts")
        .open(open)
        .collapsible(false)
        .resizable(false)
        .default_width(440.0)
        .show(ctx, |ui| {
            ui.label(
                egui::RichText::new(
                    "Press ? at any time to toggle this overlay. Most shortcuts \
                     are disabled while the command palette is open.",
                )
                .weak()
                .small(),
            );
            ui.add_space(6.0);

            egui::Grid::new("kbd_help_grid")
                .num_columns(3)
                .spacing([16.0, 4.0])
                .striped(true)
                .show(ui, |ui| {
                    // Header row.
                    ui.label(egui::RichText::new("Shortcut").strong());
                    ui.label(egui::RichText::new("Action").strong());
                    ui.label(egui::RichText::new("Description").strong());
                    ui.end_row();

                    for a in ShortcutAction::ALL {
                        ui.label(egui::RichText::new(a.binding()).monospace().strong());
                        ui.label(a.label());
                        ui.label(egui::RichText::new(a.description()).weak().small());
                        ui.end_row();
                    }
                });

            ui.add_space(8.0);
            ui.separator();
            ui.label(
                egui::RichText::new("F1 opens contextual help for whichever panel is in focus.")
                    .weak()
                    .small(),
            );
        });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cheatsheet_renders_without_panic() {
        let ctx = egui::Context::default();
        let mut open = true;
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            render_cheatsheet(ctx, &mut open);
        });
    }

    #[test]
    fn cheatsheet_lists_every_action() {
        // Smoke check: the rendered cheat-sheet must contain a row
        // per action — we can't measure that headlessly, but we can
        // ensure no action's label / binding is empty, which would
        // produce an unreadable row in the grid.
        for a in ShortcutAction::ALL {
            assert!(!a.label().is_empty());
            assert!(!a.binding().is_empty());
        }
    }
}
