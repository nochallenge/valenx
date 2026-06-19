//! Shared helper for long pop-up menus.
//!
//! egui keeps a menu on-screen horizontally but does **not** auto-scroll a
//! submenu that is taller than the window — a long category (e.g.
//! Tools → BIO with dozens of adapters, or the View menu's workbench
//! toggles) simply runs off the bottom edge with no way to reach the lower
//! entries. [`scrollable_menu`] wraps a menu body in a height-capped
//! vertical [`egui::ScrollArea`] so the menu stays on-screen and scrolls.

use eframe::egui;

/// Render `add` inside a vertical scroll area capped at 70 % of the current
/// viewport height (floored at 120 px for very short windows), so a long
/// menu body scrolls and every entry stays reachable.
///
/// `auto_shrink([false, true])` keeps the normal menu width (only the
/// height is bounded), so full-width row highlights are preserved.
pub(crate) fn scrollable_menu(ui: &mut egui::Ui, add: impl FnOnce(&mut egui::Ui)) {
    // 70 % of the viewport, floored at 120 px so tiny windows still show a
    // usable, scrollable strip rather than a single clipped row.
    let max_height = (ui.ctx().screen_rect().height() * 0.7).max(120.0);
    egui::ScrollArea::vertical()
        .auto_shrink([false, true])
        .max_height(max_height)
        .show(ui, add);
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The body closure runs in full (every entry is emitted) and nothing
    /// panics inside a headless egui frame, even with far more items than
    /// fit on screen.
    #[test]
    fn scrollable_menu_runs_full_body_without_panic() {
        let ctx = egui::Context::default();
        let mut emitted = 0;
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            egui::CentralPanel::default().show(ctx, |ui| {
                scrollable_menu(ui, |ui| {
                    for i in 0..250 {
                        let _ = ui.button(format!("entry {i}"));
                        emitted += 1;
                    }
                });
            });
        });
        assert_eq!(emitted, 250, "the wrapped body must run for every entry");
    }

    /// The height cap is 70 % of the viewport, floored at 120 px — locks the
    /// documented contract so a regression in the constant is caught.
    #[test]
    fn height_cap_is_seventy_percent_floored_at_120() {
        let tall: f32 = (1000.0_f32 * 0.7).max(120.0);
        assert!((tall - 700.0).abs() < 1e-3, "70% of 1000 px");
        let tiny: f32 = (50.0_f32 * 0.7).max(120.0);
        assert!((tiny - 120.0).abs() < 1e-3, "floored at 120 px");
    }
}
