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

/// The three header controls. Each is rendered by [`icon_button`] with
/// [`egui::Painter`] primitives rather than a font glyph, so the chrome looks
/// identical on every platform and never falls back to a "tofu" box when the
/// active font lacks a symbol code-point.
#[derive(Clone, Copy)]
enum Icon {
    /// A single horizontal bar — collapse the body to just the header.
    Minimize,
    /// Three dots in a row — open the dock / float / pop-out menu.
    More,
    /// Two crossing diagonals — hide the workbench.
    Close,
    /// A filled triangle pointing left — collapse a left-edge panel.
    ArrowLeft,
    /// A filled triangle pointing right — expand a collapsed left-edge panel.
    ArrowRight,
    /// A filled triangle pointing up — expand a collapsed bottom panel.
    ArrowUp,
    /// A filled triangle pointing down — collapse a bottom panel.
    ArrowDown,
}

/// Which edge a collapsible panel sits on, picking the arrow pair a
/// [`collapse_arrow_button`] draws: a **left**-edge panel (the Browser) points
/// `◂` when open / `▸` when collapsed; a **bottom**-edge panel (Residuals/Log)
/// points `▾` when open / `▴` when collapsed.
#[derive(Clone, Copy)]
pub(crate) enum PanelEdge {
    /// Docked to the left edge (collapses leftward).
    Left,
    /// Docked to the bottom edge (collapses downward).
    Bottom,
}

/// The clickable footprint of a header icon button (a crisp square).
const ICON_BUTTON_SIZE: f32 = 18.0;

/// Draw a painter-rendered close (`✕`) button — two crossing diagonals — at
/// the standard icon footprint, with the same hover background + tooltip as the
/// workbench header controls. Exposed so other chrome (e.g. the project-tab
/// strip's per-tab close) shares the exact same crisp ✕ instead of falling back
/// to a font glyph that can render as a "tofu" box.
pub(crate) fn close_x_button(ui: &mut egui::Ui, tip: &str) -> egui::Response {
    icon_button(ui, Icon::Close, tip)
}

/// Draw one painter-rendered header icon button with a subtle rounded hover
/// background and a tooltip. No font glyphs are used, so the control renders
/// the same regardless of the loaded font set.
fn icon_button(ui: &mut egui::Ui, icon: Icon, tip: &str) -> egui::Response {
    let (rect, resp) = ui.allocate_exact_size(
        egui::vec2(ICON_BUTTON_SIZE, ICON_BUTTON_SIZE),
        egui::Sense::click(),
    );
    // Copy the small `Copy` visuals out before borrowing the painter, so the
    // immutable `ui` borrow from `style()` does not overlap `ui.painter()`.
    let vis = *ui.style().interact(&resp);
    let painter = ui.painter();
    if resp.hovered() {
        painter.rect_filled(rect, 3.0, vis.bg_fill);
    }
    let c = rect.center();
    let col = vis.fg_stroke.color;
    let stroke = egui::Stroke::new(1.5, col);
    let r = 4.0;
    match icon {
        Icon::Close => {
            painter.line_segment([c + egui::vec2(-r, -r), c + egui::vec2(r, r)], stroke);
            painter.line_segment([c + egui::vec2(-r, r), c + egui::vec2(r, -r)], stroke);
        }
        Icon::Minimize => {
            painter.line_segment([c + egui::vec2(-r, 0.0), c + egui::vec2(r, 0.0)], stroke);
        }
        Icon::More => {
            for dx in [-5.0_f32, 0.0, 5.0] {
                painter.circle_filled(c + egui::vec2(dx, 0.0), 1.5, col);
            }
        }
        // Filled triangles (chevron-style) for the panel collapse toggles. Each
        // is three points around the centre; the apex points in the icon's
        // direction. Drawn as a filled convex polygon so it reads as a crisp
        // arrowhead at the 18px footprint regardless of the loaded font.
        Icon::ArrowLeft => {
            painter.add(egui::Shape::convex_polygon(
                vec![
                    c + egui::vec2(-r, 0.0),
                    c + egui::vec2(r, -r),
                    c + egui::vec2(r, r),
                ],
                col,
                egui::Stroke::NONE,
            ));
        }
        Icon::ArrowRight => {
            painter.add(egui::Shape::convex_polygon(
                vec![
                    c + egui::vec2(r, 0.0),
                    c + egui::vec2(-r, -r),
                    c + egui::vec2(-r, r),
                ],
                col,
                egui::Stroke::NONE,
            ));
        }
        Icon::ArrowUp => {
            painter.add(egui::Shape::convex_polygon(
                vec![
                    c + egui::vec2(0.0, -r),
                    c + egui::vec2(-r, r),
                    c + egui::vec2(r, r),
                ],
                col,
                egui::Stroke::NONE,
            ));
        }
        Icon::ArrowDown => {
            painter.add(egui::Shape::convex_polygon(
                vec![
                    c + egui::vec2(0.0, r),
                    c + egui::vec2(-r, -r),
                    c + egui::vec2(r, -r),
                ],
                col,
                egui::Stroke::NONE,
            ));
        }
    }
    resp.on_hover_text(tip)
}

/// Draw a painter-rendered **collapse / expand arrow toggle** for a docked side
/// panel, matching the workbench header chrome (no font glyph → never a "tofu"
/// box). The arrow direction follows `edge` + `open`:
///
/// - [`PanelEdge::Left`] (Browser): `◂` when open (click collapses) / `▸` when
///   collapsed (click expands);
/// - [`PanelEdge::Bottom`] (Residuals/Log): `▾` when open / `▴` when collapsed.
///
/// `name` is attached as the widget's **accessible label** (via
/// [`egui::WidgetInfo::labeled`]) so an external AI / screen-reader driver finds
/// and Invokes it by a stable Name (e.g. `"Collapse Browser"`) even though the
/// control is painter-drawn with no text. The same `name` is the hover tooltip.
/// Returns the [`egui::Response`] so the caller flips its own collapse flag on
/// `.clicked()`.
pub(crate) fn collapse_arrow_button(
    ui: &mut egui::Ui,
    edge: PanelEdge,
    open: bool,
    name: &str,
) -> egui::Response {
    let icon = match (edge, open) {
        (PanelEdge::Left, true) => Icon::ArrowLeft,
        (PanelEdge::Left, false) => Icon::ArrowRight,
        (PanelEdge::Bottom, true) => Icon::ArrowDown,
        (PanelEdge::Bottom, false) => Icon::ArrowUp,
    };
    let resp = icon_button(ui, icon, name);
    // Stable accessible Name independent of the arrow direction, so AI / screen
    // readers always find this control by the same label.
    resp.widget_info(|| egui::WidgetInfo::labeled(egui::WidgetType::Button, true, name));
    resp
}

/// Draw the header row: the title on the left and, pinned to the far right, a
/// uniform `−  ⋯  ✕` cluster (Close rightmost). The buttons are painter-drawn
/// so they never render as missing-glyph boxes. Mutates the passed flags; sets
/// `*close` when the close button is clicked.
fn panel_header(
    ui: &mut egui::Ui,
    title: &str,
    collapsed: &mut bool,
    mode: &mut PanelMode,
    close: &mut bool,
) {
    ui.horizontal(|ui| {
        // Keep the controls a fixed, equal distance apart on every panel.
        ui.spacing_mut().item_spacing.x = 4.0;
        // Title, vertically centred against the icon cluster.
        ui.with_layout(egui::Layout::left_to_right(egui::Align::Center), |ui| {
            ui.strong(title);
        });
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            // Right-to-left layout lays the first-added widget out rightmost,
            // so add Close first to get the visual order `−  ⋯  ✕`.
            if icon_button(ui, Icon::Close, "Close — reopen from the View menu").clicked() {
                *close = true;
            }
            // "More" opens a small menu via a toggled popup keyed on the
            // button's own id, so the trigger stays a painter-drawn icon.
            let more = icon_button(ui, Icon::More, "Dock, float, or pop out");
            let popup_id = more.id.with("workbench_chrome_more");
            if more.clicked() {
                ui.memory_mut(|m| m.toggle_popup(popup_id));
            }
            egui::popup::popup_below_widget(
                ui,
                popup_id,
                &more,
                egui::popup::PopupCloseBehavior::CloseOnClick,
                |ui| {
                    ui.set_min_width(170.0);
                    if ui.button("Dock right").clicked() {
                        *mode = PanelMode::Docked;
                    }
                    if ui.button("Float (in-app window)").clicked() {
                        *mode = PanelMode::Floating;
                    }
                    if ui.button("Pop out (new window)").clicked() {
                        *mode = PanelMode::Window;
                    }
                },
            );
            let tip = if *collapsed { "Expand" } else { "Minimize" };
            if icon_button(ui, Icon::Minimize, tip).clicked() {
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
    // Historical signature: docked at egui's default width. Delegates to the
    // sized variant with no width override so every existing caller is byte-for
    // -byte unchanged.
    workbench_shell_sized(app, ctx, id, title, None, body)
}

/// [`workbench_shell`] plus an optional **docked balanced width**. When
/// `docked_width` is `Some((default, min))`, the docked-right [`egui::SidePanel`]
/// opens at `default` px (clamped to at least `min`) and can't be dragged
/// narrower than `min`. This is FIX B for calculator-style tabs: the classic
/// Assistant panel gets ~40% of the window (min ~320px) so the workbench form
/// in the centre and the Assistant chat are **both** readable, instead of the
/// Assistant collapsing to a sliver. `None` keeps egui's default width (the
/// behaviour every other workbench keeps). The floating / pop-out modes are
/// unaffected — the balanced width only makes sense for the docked edge.
pub fn workbench_shell_sized(
    app: &mut ValenxApp,
    ctx: &egui::Context,
    id: &str,
    title: &str,
    docked_width: Option<(f32, f32)>,
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
            let mut panel = egui::SidePanel::right(egui::Id::new(id)).resizable(true);
            if let Some((default_w, min_w)) = docked_width {
                let default_w = default_w.max(min_w);
                panel = panel
                    .default_width(default_w)
                    .width_range(min_w..=f32::INFINITY);
            }
            panel.show(ctx, |ui| {
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

    #[test]
    fn icon_buttons_draw_headless_without_panicking() {
        // Each icon is painter-drawn (no font glyph), so it must render on a
        // bare headless context regardless of the loaded font set.
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            egui::CentralPanel::default().show(ctx, |ui| {
                let _ = icon_button(ui, Icon::Minimize, "Minimize");
                let _ = icon_button(ui, Icon::More, "More");
                let _ = icon_button(ui, Icon::Close, "Close");
                // The four collapse arrows are painter-drawn polygons; they must
                // also render headlessly.
                let _ = icon_button(ui, Icon::ArrowLeft, "left");
                let _ = icon_button(ui, Icon::ArrowRight, "right");
                let _ = icon_button(ui, Icon::ArrowUp, "up");
                let _ = icon_button(ui, Icon::ArrowDown, "down");
            });
        });
    }

    #[test]
    fn collapse_arrow_button_draws_both_edges_and_states() {
        // The arrow toggle must render for each (edge, open) combination on a
        // bare headless context — it's a painter-drawn polygon plus an explicit
        // accessible Name, no font glyph.
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            egui::CentralPanel::default().show(ctx, |ui| {
                let _ = collapse_arrow_button(ui, PanelEdge::Left, true, "Collapse Browser");
                let _ = collapse_arrow_button(ui, PanelEdge::Left, false, "Expand Browser");
                let _ = collapse_arrow_button(ui, PanelEdge::Bottom, true, "Collapse log");
                let _ = collapse_arrow_button(ui, PanelEdge::Bottom, false, "Expand log");
            });
        });
    }

    #[test]
    fn panel_header_draws_and_stays_open_without_a_click() {
        // With no synthesised pointer input the close button is never clicked,
        // so the header reports the panel should stay open.
        let ctx = egui::Context::default();
        let mut collapsed = false;
        let mut mode = PanelMode::Docked;
        let mut close = true;
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            egui::CentralPanel::default().show(ctx, |ui| {
                close = false;
                panel_header(ui, "Test", &mut collapsed, &mut mode, &mut close);
            });
        });
        assert!(!close, "header reports close only when ✕ is clicked");
        assert!(!collapsed);
        assert_eq!(mode, PanelMode::Docked);
    }
}
