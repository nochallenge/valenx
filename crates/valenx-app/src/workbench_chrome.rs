//! Shared chrome (collapse / pop-out / close) for the right-side workbench
//! panels.
//!
//! Every right-side workbench used to draw its own [`egui::SidePanel`] with a
//! bespoke header. This module factors out a single [`workbench_shell`] that
//! wraps a workbench body in one of three containers — a docked
//! [`egui::SidePanel`], a floating in-app [`egui::Window`], or a *separate OS
//! window* (an immediate child viewport) — and stamps every one with the same
//! header: a title plus `collapse / pop-out / ✕` controls.
//!
//! The pop-out menu offers both "Float (in-app window)" and "Pop out (new OS
//! window)", so a workbench can be torn off the right edge into its own egui
//! window or its own native window and docked back at any time. The chrome
//! state for each panel (collapsed? which mode?) lives in a map on the
//! [`ValenxApp`] keyed by the panel's stable id string, so it survives across
//! frames without each workbench having to carry its own flags.
//!
//! The separate-OS-window mode uses [`egui::Context::show_viewport_immediate`],
//! an immediate-mode child viewport that runs its UI callback inline and can
//! therefore borrow `&mut ValenxApp` just like the docked path. On a backend
//! that cannot open extra native windows (or a headless test context),
//! `show_viewport_immediate` transparently falls back to embedding the callback
//! in the current viewport, so the call is always safe.

use eframe::egui;

use crate::ValenxApp;

/// Where a workbench panel is currently shown.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum PanelMode {
    /// Docked to the right edge as an [`egui::SidePanel`] (the classic layout).
    #[default]
    Docked,
    /// Floating as a draggable in-app [`egui::Window`].
    Floating,
    /// Popped out into a separate OS window (an immediate child viewport).
    Window,
}

/// Per-panel chrome state: whether the body is collapsed and where it is shown.
#[derive(Clone, Copy, Debug, Default)]
pub struct PanelChromeState {
    /// When `true` the header is drawn but the body is hidden.
    pub collapsed: bool,
    /// Docked / floating / separate-OS-window.
    pub mode: PanelMode,
}

/// Draw the header row: title + `[collapse][pop-out menu][✕]`. Mutates the
/// passed flags; sets `*close` when ✕ is clicked.
fn panel_header(
    ui: &mut egui::Ui,
    title: &str,
    collapsed: &mut bool,
    mode: &mut PanelMode,
    close: &mut bool,
) {
    ui.horizontal(|ui| {
        ui.strong(title);
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            if ui
                .small_button("✕")
                .on_hover_text("Close (reopen from the View / Tools menu)")
                .clicked()
            {
                *close = true;
            }
            ui.menu_button("⧉", |ui| {
                if ui.button("Dock right").clicked() {
                    *mode = PanelMode::Docked;
                    ui.close_menu();
                }
                if ui.button("Float (in-app window)").clicked() {
                    *mode = PanelMode::Floating;
                    ui.close_menu();
                }
                if ui.button("Pop out (new OS window)").clicked() {
                    *mode = PanelMode::Window;
                    ui.close_menu();
                }
            })
            .response
            .on_hover_text("Pop out");
            let (glyph, tip) = if *collapsed {
                ("▸", "Expand")
            } else {
                ("▾", "Collapse")
            };
            if ui.small_button(glyph).on_hover_text(tip).clicked() {
                *collapsed = !*collapsed;
            }
        });
    });
    ui.separator();
}

/// Wrap a workbench body in the right container (docked [`egui::SidePanel`] /
/// floating [`egui::Window`] / separate OS viewport) with the standard header.
/// Returns `true` if the user clicked ✕ (the caller flips its own show-flag).
///
/// `id` must be a stable unique string — reuse the panel's existing
/// `SidePanel` id. The header's collapse and pop-out controls mutate the
/// chrome state stored on `app` under `id`, so the choice persists across
/// frames.
pub fn workbench_shell(
    app: &mut ValenxApp,
    ctx: &egui::Context,
    id: &str,
    title: &str,
    body: impl FnOnce(&mut ValenxApp, &mut egui::Ui),
) -> bool {
    let st = app.workbench_chrome.get(id).copied().unwrap_or_default();
    let mut collapsed = st.collapsed;
    let mut mode = st.mode;
    let mut close = false;
    match mode {
        PanelMode::Docked => {
            // `SidePanel::right` borrows a `&str` id for the whole `show`
            // closure, which would force `id: &'static str`; hash it into an
            // owned `egui::Id` up front so the non-static `id` doesn't escape.
            egui::SidePanel::right(egui::Id::new(id))
                .resizable(true)
                .show(ctx, |ui| {
                    panel_header(ui, title, &mut collapsed, &mut mode, &mut close);
                    if !collapsed {
                        body(app, ui);
                    }
                });
        }
        PanelMode::Floating => {
            egui::Window::new(title)
                .id(egui::Id::new(id))
                .resizable(true)
                .default_width(340.0)
                .show(ctx, |ui| {
                    panel_header(ui, title, &mut collapsed, &mut mode, &mut close);
                    if !collapsed {
                        body(app, ui);
                    }
                });
        }
        PanelMode::Window => {
            let vid = egui::ViewportId::from_hash_of(id);
            ctx.show_viewport_immediate(
                vid,
                egui::ViewportBuilder::default()
                    .with_title(title)
                    .with_inner_size([380.0, 680.0]),
                |vctx, _class| {
                    egui::CentralPanel::default().show(vctx, |ui| {
                        panel_header(ui, title, &mut collapsed, &mut mode, &mut close);
                        if !collapsed {
                            body(app, ui);
                        }
                    });
                    // The OS-window close button (the title-bar ✕) re-docks the
                    // panel rather than hiding the workbench, so the user never
                    // loses it by closing the torn-off window.
                    if vctx.input(|i| i.viewport().close_requested()) {
                        mode = PanelMode::Docked;
                    }
                },
            );
        }
    }
    app.workbench_chrome
        .insert(id.to_string(), PanelChromeState { collapsed, mode });
    close
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_mode_is_docked() {
        assert_eq!(PanelMode::default(), PanelMode::Docked);
    }

    #[test]
    fn default_state_is_expanded_and_docked() {
        let st = PanelChromeState::default();
        assert!(!st.collapsed, "a fresh panel is not collapsed");
        assert_eq!(st.mode, PanelMode::Docked);
    }
}
