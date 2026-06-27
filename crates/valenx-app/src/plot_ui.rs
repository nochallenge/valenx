//! Shared **managed-plot** UX for the workbenches.
//!
//! Many workbenches embed an [`egui_plot::Plot`] inside a long scrolling
//! parameter form. Left to default behaviour the plot scrolls away with the
//! form, and its scroll-wheel zoom is swallowed by the surrounding
//! [`egui::ScrollArea`]. [`managed_plot`] is the reusable fix first developed
//! for the Rocket workbench: it pins a fully *maneuverable*, vertically
//! resizable plot with a discoverable, AI-drivable **"Reset view"** control.
//!
//! Adopt it from any workbench by:
//! 1. rendering the plot *above* / outside the form's `ScrollArea` (so its
//!    wheel-zoom isn't swallowed), and
//! 2. routing the draw through [`managed_plot`] with a per-panel `id_salt`,
//!    a caller-owned one-shot `reset` flag, and a `build` closure.

use eframe::egui;
use egui_plot::{Plot, PlotUi};

/// A fully *maneuverable*, vertically-resizable plot with a discoverable
/// "Reset view" control — the managed-plot pattern the workbenches pin above
/// their scrolling forms.
///
/// Draws a small header row holding a named **"Reset view"** button (so the
/// reset is AI-drivable / screen-reader-named, not just a hidden double-click),
/// then the plot inside an [`egui::Resize`] frame so the user can drag its
/// bottom edge to grow / shrink it. The [`Plot`] is configured for full
/// manipulation — drag-to-pan, scroll-to-zoom while hovered, box-zoom (drag
/// with the secondary button), and double-click-to-reset — and, because it is
/// rendered *outside* the form's `ScrollArea`, its scroll-wheel zoom is not
/// swallowed by a surrounding scroll container.
///
/// `reset` is a caller-owned one-shot flag: the "Reset view" button sets it,
/// and this function consumes it (forcing the plot back to auto-bounds on that
/// frame, then clearing the flag). `default_height` is the plot's initial
/// height in points; the user-dragged height persists across frames under
/// `id_salt` via egui's `Resize` state.
///
/// Reusable: other workbenches adopt the same pinned-maneuverable-plot UX by
/// calling this with their own `id_salt` and a `build` closure that draws their
/// lines / points.
pub(crate) fn managed_plot(
    ui: &mut egui::Ui,
    id_salt: &str,
    default_height: f32,
    reset: &mut bool,
    build: impl FnOnce(&mut PlotUi),
) {
    // The plain form: no extra Plot configuration (the identity configurator).
    managed_plot_cfg(ui, id_salt, default_height, reset, |plot| plot, build);
}

/// As [`managed_plot`], but with a `configure` hook that runs on the [`Plot`]
/// builder *before* the manipulation flags and `.show()` — so callers can keep
/// their per-panel `Plot` setup (axis labels, [`egui_plot::Legend`], data
/// aspect, …) while still getting the pinned / pan / zoom / box-zoom / reset /
/// resize behaviour. The manipulation flags are applied **after** `configure`,
/// so they always win (a panel can't accidentally disable pan/zoom).
pub(crate) fn managed_plot_cfg(
    ui: &mut egui::Ui,
    id_salt: &str,
    default_height: f32,
    reset: &mut bool,
    configure: impl FnOnce(Plot) -> Plot,
    build: impl FnOnce(&mut PlotUi),
) {
    // Discoverable, named reset control (double-click on the plot also resets).
    ui.horizontal(|ui| {
        if ui
            .button("⟲ Reset view")
            .on_hover_text(
                "Restore the auto-fit bounds after panning / zooming (or double-click the plot).",
            )
            .clicked()
        {
            *reset = true;
        }
        ui.label(
            egui::RichText::new(
                "drag = pan · scroll = zoom · right-drag = box-zoom · drag bottom edge = resize",
            )
            .weak()
            .small(),
        );
    });

    let do_reset = *reset;
    egui::Resize::default()
        .id_source(id_salt)
        .default_height(default_height)
        .min_height(120.0)
        // Vertical resize only — width follows the panel.
        .resizable([false, true])
        .show(ui, |ui| {
            let mut plot = configure(Plot::new(id_salt))
                .allow_drag(true) // drag-to-pan
                .allow_zoom(true) // scroll / pinch zoom
                .allow_scroll(true) // wheel zoom while hovered (not swallowed: pinned outside the ScrollArea)
                .allow_boxed_zoom(true) // secondary-drag box zoom
                .allow_double_click_reset(true); // double-click restores auto-bounds
            if do_reset {
                plot = plot.reset();
            }
            plot.show(ui, build);
        });
    *reset = false;
}

/// As [`managed_plot`], but the one-shot reset flag is stored in egui's
/// per-widget temporary memory keyed by `id_salt` instead of a caller-owned
/// `&mut bool`. This lets call sites that have **no convenient state struct**
/// (e.g. free `draw_*` helpers that take borrowed data) still get the full
/// pinned / maneuverable / resettable plot UX without threading a `plot_reset`
/// field through. Each plot needs a unique `id_salt`.
pub(crate) fn managed_plot_mem(
    ui: &mut egui::Ui,
    id_salt: &str,
    default_height: f32,
    build: impl FnOnce(&mut PlotUi),
) {
    managed_plot_mem_cfg(ui, id_salt, default_height, |plot| plot, build);
}

/// [`managed_plot_mem`] with the [`managed_plot_cfg`] `configure` hook (axis
/// labels, legend, …). Reset state lives in egui memory keyed by `id_salt`.
pub(crate) fn managed_plot_mem_cfg(
    ui: &mut egui::Ui,
    id_salt: &str,
    default_height: f32,
    configure: impl FnOnce(Plot) -> Plot,
    build: impl FnOnce(&mut PlotUi),
) {
    let key = egui::Id::new(("managed_plot_reset", id_salt));
    let mut reset = ui.data_mut(|d| d.get_temp::<bool>(key).unwrap_or(false));
    managed_plot_cfg(ui, id_salt, default_height, &mut reset, configure, build);
    ui.data_mut(|d| d.insert_temp(key, reset));
}
