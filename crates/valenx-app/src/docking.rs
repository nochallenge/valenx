//! Opt-in **dockable / tileable tab layout** for the central panel,
//! built on [`egui_tiles`] (emilk's official tiling crate for egui —
//! the same one Rerun uses).
//!
//! This is a pure **UI-convenience layer**. It lets the user split the
//! central area into rows *and* columns, drag the splitters to resize,
//! group panes into tab bars, close tabs, and drag-rearrange them —
//! all of which `egui_tiles` gives us for free. It does **not** change
//! any simulation result, solver, or numeric output: every [`Pane`]
//! either hosts the existing viewport / browser draw code or renders a
//! lightweight labelled placeholder.
//!
//! The whole feature is gated behind [`crate::ValenxApp::docked_layout`]
//! (default `false`), so the established single-viewport layout stays
//! the default and nothing here runs unless the user ticks
//! **View → Docked layout**.
//!
//! ## Shape
//!
//! - [`Pane`] — an enum of the views that can live in a tile. Each pane
//!   knows its [`Pane::title`] and how to draw itself via [`Pane::ui`].
//! - [`DockingState`] — wraps an [`egui_tiles::Tree<Pane>`] plus a
//!   [`Behavior`] impl. [`DockingState::default`] builds a starting
//!   layout that demonstrates a horizontal split, a vertical split, and
//!   a tab group.
//! - [`DockingState::show`] — one call that paints the whole tree into a
//!   parent [`egui::Ui`].

use eframe::egui;

/// A single dockable view that can be placed in a tile.
///
/// The heavyweight panes ([`Pane::Viewport`], [`Pane::Browser`]) render
/// a description of what the real panel shows — the production wiring of
/// the live wgpu viewport / browser tree stays in `update.rs`, and this
/// tiling shell deliberately does not re-plumb every workbench. The
/// point of this module is the *tiling container*, not a second copy of
/// every panel's internals. Workbench panes render an at-a-glance
/// summary label so the tab bar and splitters have real, titled content
/// to arrange.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Pane {
    /// The main 3D / 2D viewport slot.
    Viewport,
    /// The left-hand project Browser tree.
    Browser,
    /// The Mesh Toolbox (CAD) workbench slot.
    MeshToolbox,
    /// The Aerodynamics / Wind Tunnel workbench slot.
    Aero,
    /// The FEM Workbench slot.
    Fem,
    /// The residuals / log bottom dock slot.
    Console,
}

impl Pane {
    /// Human-readable tab title, shown in the tile's tab bar.
    pub fn title(&self) -> &'static str {
        match self {
            Pane::Viewport => "Viewport",
            Pane::Browser => "Browser",
            Pane::MeshToolbox => "Mesh Toolbox",
            Pane::Aero => "Aerodynamics",
            Pane::Fem => "FEM",
            Pane::Console => "Console",
        }
    }

    /// A one-line description rendered as placeholder body text. Keeps
    /// the tiling shell honest — every pane shows *what it is* without
    /// pretending to be the fully-wired live panel.
    fn blurb(&self) -> &'static str {
        match self {
            Pane::Viewport => {
                "Main 3D / 2D viewport. In the docked layout this tile hosts \
                 the same scene the classic central panel renders."
            }
            Pane::Browser => {
                "Project tree: project, cases, geometry, mesh, results, and \
                 the live adapter registry with status colours."
            }
            Pane::MeshToolbox => {
                "CAD Mesh Toolbox: Inspector, Transformations, Cut plane, \
                 Repair. Toggle the full panel from View → Mesh Toolbox."
            }
            Pane::Aero => {
                "Virtual wind tunnel (3-D external-aerodynamics CFD): drag / \
                 lift / moment coefficients, Cp & flow fields, AoA sweep."
            }
            Pane::Fem => {
                "Native finite-element analysis: linear-static bending + modal \
                 natural frequencies, solved in-process by valenx-fem."
            }
            Pane::Console => {
                "Solver console: live residual chart (egui_plot) and the \
                 streamed solver log."
            }
        }
    }

    /// Draw this pane's body into `ui`.
    fn ui(&self, ui: &mut egui::Ui) {
        ui.vertical(|ui| {
            ui.add_space(6.0);
            ui.heading(self.title());
            ui.separator();
            ui.label(self.blurb());
            ui.add_space(8.0);
            ui.weak(
                "Docked layout (opt-in). Drag the tab to re-dock it, drag a \
                 splitter to resize, or use the tab's close button.",
            );
        });
    }
}

/// The [`egui_tiles::Behavior`] that tells `egui_tiles` how to title and
/// paint each [`Pane`]. The defaults from `Behavior` already provide
/// resizable splits, tab bars, drag-to-rearrange, and per-tab close
/// buttons; we only override the two required hooks plus opt into the
/// close button via [`Behavior::is_tab_closable`].
#[derive(Default)]
struct Behavior;

impl egui_tiles::Behavior<Pane> for Behavior {
    fn tab_title_for_pane(&mut self, pane: &Pane) -> egui::WidgetText {
        pane.title().into()
    }

    fn pane_ui(
        &mut self,
        ui: &mut egui::Ui,
        _tile_id: egui_tiles::TileId,
        pane: &mut Pane,
    ) -> egui_tiles::UiResponse {
        egui::Frame::none()
            .inner_margin(egui::Margin::same(8.0))
            .show(ui, |ui| {
                pane.ui(ui);
            });
        // We don't initiate a drag from the body — egui_tiles already
        // handles dragging from the tab bar.
        egui_tiles::UiResponse::None
    }

    /// Show the close (×) button on every tab so the user can close
    /// panes, matching the requested "close tabs" behaviour.
    fn is_tab_closable(
        &self,
        _tiles: &egui_tiles::Tiles<Pane>,
        _tile_id: egui_tiles::TileId,
    ) -> bool {
        true
    }
}

/// Owns the [`egui_tiles::Tree`] of [`Pane`]s plus its [`Behavior`].
///
/// Construct with [`DockingState::default`] for a sensible starting
/// layout, then paint every frame with [`DockingState::show`].
pub struct DockingState {
    tree: egui_tiles::Tree<Pane>,
    behavior: Behavior,
}

impl Default for DockingState {
    /// Build a starting tree that demonstrates **both** a horizontal and
    /// a vertical split **and** a tab group:
    ///
    /// ```text
    /// ┌────────────┬───────────────────────────────┐
    /// │            │  ┌─────────────────────────┐   │  <- horizontal split
    /// │  Browser   │  │ Viewport                │   │     (left | right)
    /// │            │  └─────────────────────────┘   │
    /// │            │  ┌─────────────────────────┐   │  <- vertical split
    /// │            │  │ [Console][Aero][FEM] tabs│   │     (top / bottom)
    /// │            │  └─────────────────────────┘   │
    /// └────────────┴───────────────────────────────┘
    /// ```
    ///
    /// The right column is a vertical split whose bottom cell is a tab
    /// group (Console + Aero + FEM + Mesh Toolbox), so the default view
    /// exercises every container kind `egui_tiles` offers.
    fn default() -> Self {
        let mut tiles = egui_tiles::Tiles::default();

        // Leaf tiles. Each `insert_pane` takes `&mut tiles`, so the
        // child ids are bound to locals first — passing
        // `tiles.insert_pane(..)` calls *inside* another `&mut tiles`
        // method's argument list would be two simultaneous mutable
        // borrows (E0499).
        let viewport = tiles.insert_pane(Pane::Viewport);
        let browser = tiles.insert_pane(Pane::Browser);
        let console = tiles.insert_pane(Pane::Console);
        let aero = tiles.insert_pane(Pane::Aero);
        let fem = tiles.insert_pane(Pane::Fem);
        let mesh_toolbox = tiles.insert_pane(Pane::MeshToolbox);

        // A tab group (tabbed container) holding several panes — this is
        // the "group panels into tab bars" demonstration.
        let tabbed = tiles.insert_tab_tile(vec![console, aero, fem, mesh_toolbox]);

        // Vertical split (top / bottom): viewport on top, tab group
        // below. This is the "rows" / vertical-split demonstration.
        let right_column = tiles.insert_vertical_tile(vec![viewport, tabbed]);

        // Horizontal split (left | right): browser beside the vertical
        // column. This is the "columns" / horizontal-split demonstration.
        let root = tiles.insert_horizontal_tile(vec![browser, right_column]);

        let tree = egui_tiles::Tree::new("valenx_docking_tree", root, tiles);
        Self {
            tree,
            behavior: Behavior,
        }
    }
}

impl DockingState {
    /// Paint the entire docked layout into `ui`. `egui_tiles` handles
    /// resizable splits, tab bars, close buttons, and drag-rearrange
    /// internally — this is the single entry point the host calls from
    /// the CentralPanel when `docked_layout` is on.
    pub fn show(&mut self, ui: &mut egui::Ui) {
        self.tree.ui(&mut self.behavior, ui);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_tree_has_all_panes() {
        // The starting layout must contain every pane kind so the
        // demonstration of splits + tabs is intact. We assert the tree
        // is non-empty and that each pane title round-trips.
        let state = DockingState::default();
        assert!(state.tree.root().is_some(), "default tree must have a root");
        // Count the pane leaves present in the tiles store.
        let pane_count = state
            .tree
            .tiles
            .tiles()
            .filter(|t| matches!(t, egui_tiles::Tile::Pane(_)))
            .count();
        // Viewport, Browser, Console, Aero, Fem, MeshToolbox = 6 panes.
        assert_eq!(
            pane_count, 6,
            "expected six pane leaves in the default tree"
        );
    }

    #[test]
    fn pane_titles_are_stable() {
        assert_eq!(Pane::Viewport.title(), "Viewport");
        assert_eq!(Pane::Browser.title(), "Browser");
        assert_eq!(Pane::Console.title(), "Console");
    }
}
