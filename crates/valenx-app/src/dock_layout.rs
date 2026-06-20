//! Opt-in **dockable / tileable layout for the right-side workbench
//! panels**, built on [`egui_tiles`] (emilk's official tiling crate for
//! egui — the same one Rerun and the central-panel [`crate::docking`]
//! layout use).
//!
//! ## What it gives the user
//!
//! With the classic layout, each open workbench is its own stacked
//! [`egui::SidePanel`] on the right edge — they can be collapsed / floated /
//! popped out (via [`crate::workbench_chrome`]) but not *reordered relative
//! to each other*. Ticking **View → "Dockable panel layout (beta)"** flips
//! [`crate::ValenxApp::dock_enabled`] on and replaces that whole run with a
//! single resizable right region hosting every open workbench as a tile in
//! an `egui_tiles` [`Tree`]. `egui_tiles` then provides, for free:
//!
//! - **drag a panel by its tab to reorder** it among the others (they
//!   reflow),
//! - **drop it into a second row** (or column) to split the region,
//! - **group panels into a shared tab bar**, and
//! - **resize the splits** by dragging the boundaries.
//!
//! This is a pure **hosting / layout** layer: every tile renders the *same*
//! `<name>_workbench_body` function the classic [`egui::SidePanel`] renders,
//! so no simulation, solver, or numeric output changes — only *where* the
//! panel is drawn. Turning the toggle back off restores the classic stacked
//! layout exactly.
//!
//! ## Safety / non-regression
//!
//! - Default **off** ([`crate::ValenxApp::dock_enabled`] defaults `false`);
//!   nothing in this module runs unless the user opts in.
//! - The [`Tree`] is **lazily built** and **synced** every frame to the set
//!   of currently-open workbenches: opening a workbench adds a tile, closing
//!   one (here or from the View menu) drops its tile, and when none are open
//!   the tree is dropped and the region paints nothing.
//! - Distinct from [`crate::docking`] (which tiles the *central* viewport):
//!   this tiles the *right-side workbenches*.

use eframe::egui;

use crate::ValenxApp;

/// The right-side workbenches that are wired into the dockable layout, as
/// `(panel_id, human_title)`. The `panel_id` is the **same stable string**
/// each workbench passes to [`crate::workbench_chrome::workbench_shell`], so
/// the dock tab and the classic panel share one identity. The order here is
/// the order tiles are created in when several are already open on the first
/// dock frame.
///
/// To wire another workbench: extract its body into a
/// `pub(crate) fn <name>_workbench_body(app, ui)` (see any entry in
/// [`render_panel_body`]), add its `(id, title)` here *and* its
/// `is_panel_open` arm *and* its `render_panel_body` arm.
pub(crate) const DOCKABLE_PANELS: &[(&str, &str)] = &[
    ("valenx_mesh_toolbox", "Mesh Toolbox"),
    ("valenx_genetics_workbench", "Genetics Workbench"),
    ("valenx_aero_workbench", "Wind Tunnel"),
    ("valenx_fem_workbench", "FEM Workbench"),
    ("valenx_cfd_workbench", "CFD Workbench"),
    ("valenx_astro_workbench", "Astro / Launch"),
    ("valenx_rocket_workbench", "Rocket — design → simulate"),
    ("valenx_engine_workbench", "Engine — design → analyze"),
    ("valenx_car_workbench", "Car — design → simulate"),
    ("valenx_assistant_panel", "Assistant"),
];

/// Tile-id prefix for a **"Workbench + Agent"** unit's empty build-canvas
/// half: `"workspace:<n>"`. Paired with [`AGENT_PREFIX`].
const WORKSPACE_PREFIX: &str = "workspace:";
/// Tile-id prefix for a **"Workbench + Agent"** unit's Claude-chat half:
/// `"agent:<n>"`. Paired with [`WORKSPACE_PREFIX`].
const AGENT_PREFIX: &str = "agent:";

/// Is this a special "Workbench + Agent" tile id (either half)? These are
/// inserted by the launcher rather than coming from [`DOCKABLE_PANELS`], are
/// **not** gated on any `show_*` flag, and must survive [`sync_tree`]'s
/// open-set pruning until the user closes them.
fn is_wb_agent_pane(panel_id: &str) -> bool {
    panel_id.starts_with(WORKSPACE_PREFIX) || panel_id.starts_with(AGENT_PREFIX)
}

/// The human title for a dockable `panel_id`, used for the tile's tab.
///
/// - `"workspace:<n>"` → `"Workspace N"`, `"agent:<n>"` → `"Agent N"` (the
///   paired "Workbench + Agent" tiles).
/// - Otherwise a [`DOCKABLE_PANELS`] title, falling back to the raw id for an
///   unrecognised pane (shouldn't happen — every other pane we insert comes
///   from [`DOCKABLE_PANELS`]).
fn panel_title(panel_id: &str) -> String {
    if let Some(n) = panel_id.strip_prefix(WORKSPACE_PREFIX) {
        return format!("Workspace {n}");
    }
    if let Some(n) = panel_id.strip_prefix(AGENT_PREFIX) {
        return format!("Agent {n}");
    }
    DOCKABLE_PANELS
        .iter()
        .find(|(id, _)| *id == panel_id)
        .map(|(_, title)| (*title).to_string())
        .unwrap_or_else(|| panel_id.to_string())
}

/// Is the workbench identified by `panel_id` currently open (its `show_*`
/// flag set)? Drives the per-frame sync: a tile exists in the dock tree iff
/// this returns `true`. Unknown ids return `false` so a stale pane is
/// dropped rather than rendered as a stub forever.
fn is_panel_open(app: &ValenxApp, panel_id: &str) -> bool {
    match panel_id {
        "valenx_mesh_toolbox" => app.show_mesh_toolbox,
        "valenx_genetics_workbench" => app.show_genetics_workbench,
        "valenx_aero_workbench" => app.show_aero_workbench,
        "valenx_fem_workbench" => app.show_fem_workbench,
        "valenx_cfd_workbench" => app.show_cfd_workbench,
        "valenx_astro_workbench" => app.show_astro_workbench,
        "valenx_rocket_workbench" => app.show_rocket_workbench,
        "valenx_engine_workbench" => app.show_engine_workbench,
        "valenx_car_workbench" => app.show_car_workbench,
        "valenx_assistant_panel" => app.show_assistant_panel,
        _ => false,
    }
}

/// Close the workbench identified by `panel_id` (clear its `show_*` flag),
/// invoked when the user clicks a tab's ✕ in the dock. Mirrors what the
/// classic panel's own ✕ does, so closing behaves identically in both hosts.
fn close_panel(app: &mut ValenxApp, panel_id: &str) {
    match panel_id {
        "valenx_mesh_toolbox" => app.show_mesh_toolbox = false,
        "valenx_genetics_workbench" => app.show_genetics_workbench = false,
        "valenx_aero_workbench" => app.show_aero_workbench = false,
        "valenx_fem_workbench" => app.show_fem_workbench = false,
        "valenx_cfd_workbench" => app.show_cfd_workbench = false,
        "valenx_astro_workbench" => app.show_astro_workbench = false,
        "valenx_rocket_workbench" => app.show_rocket_workbench = false,
        "valenx_engine_workbench" => app.show_engine_workbench = false,
        "valenx_car_workbench" => app.show_car_workbench = false,
        "valenx_assistant_panel" => app.show_assistant_panel = false,
        _ => {}
    }
}

/// Render a workbench's body into a dock tile, dispatching on its
/// `panel_id`. Each arm calls the very same `<name>_workbench_body` the
/// classic [`egui::SidePanel`] path calls, so there is **no duplicated panel
/// logic** — only one source of truth per workbench.
///
/// A `panel_id` that isn't wired yet renders a small graceful notice rather
/// than panicking, so partially-wired states stay usable.
pub(crate) fn render_panel_body(app: &mut ValenxApp, ui: &mut egui::Ui, panel_id: &str) {
    // "Workbench + Agent" tiles: the agent half routes to the single shared
    // Claude chat bridge; the workspace half is an empty build canvas.
    if panel_id.starts_with(AGENT_PREFIX) {
        crate::assistant_workbench::assistant_workbench_body(app, ui);
        return;
    }
    if let Some(n) = panel_id.strip_prefix(WORKSPACE_PREFIX) {
        render_workspace_body(app, ui, n);
        return;
    }
    match panel_id {
        "valenx_mesh_toolbox" => crate::mesh_toolbox::mesh_toolbox_body(app, ui),
        "valenx_genetics_workbench" => crate::genetics_workbench::genetics_workbench_body(app, ui),
        "valenx_aero_workbench" => crate::aero_workbench::aero_workbench_body(app, ui),
        "valenx_fem_workbench" => crate::fem_workbench::fem_workbench_body(app, ui),
        "valenx_cfd_workbench" => crate::cfd_workbench::cfd_workbench_body(app, ui),
        "valenx_astro_workbench" => crate::astro_workbench::astro_workbench_body(app, ui),
        "valenx_rocket_workbench" => crate::rocket_workbench::rocket_workbench_body(app, ui),
        "valenx_engine_workbench" => crate::engine_workbench::engine_workbench_body(app, ui),
        "valenx_car_workbench" => crate::car_workbench::car_workbench_body(app, ui),
        "valenx_assistant_panel" => crate::assistant_workbench::assistant_workbench_body(app, ui),
        _ => {
            ui.label("This panel isn't dockable yet — turn off Dockable layout to use it.");
        }
    }
}

/// Render a `"workspace:<n>"` tile body — the **empty build canvas** half of a
/// "Workbench + Agent" unit, where `n` is the unit number (the suffix after
/// the `"workspace:"` prefix).
///
/// **Honest first cut:** this is a *labelled placeholder canvas*, not a live
/// per-unit 3-D viewport — six independent wgpu scenes is a deep follow-up and
/// out of scope here. It shows a `"Workspace N"` heading, a weak hint that the
/// agent on the right builds here, and — if the assistant feed has reported a
/// build/result/ship line — the latest such status, all inside a bordered,
/// scrollable region.
fn render_workspace_body(app: &mut ValenxApp, ui: &mut egui::Ui, n: &str) {
    ui.heading(format!("Workspace {n}"));
    ui.label(
        egui::RichText::new("Ask the agent on the right to build something here.")
            .weak()
            .small(),
    );
    ui.add_space(6.0);
    // Pull the latest build status *before* the borrow of `ui`'s frame so the
    // bordered canvas can show it (best-effort; `None` → just the placeholder).
    let status = crate::assistant_workbench::latest_build_status(app);
    egui::Frame::group(ui.style())
        .inner_margin(egui::Margin::same(8.0))
        .show(ui, |ui| {
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| match status {
                    Some(line) => {
                        ui.label(
                            egui::RichText::new("Latest from the agent")
                                .strong()
                                .small(),
                        );
                        ui.label(line);
                    }
                    None => {
                        ui.centered_and_justified(|ui| {
                            ui.label(egui::RichText::new("Empty workspace").weak().italics());
                        });
                    }
                });
        });
}

/// After the tree has been drawn, drain the per-workbench deferred requests
/// (3-D mesh loads, field-overlay pushes) the bodies may have set. These run
/// *outside* any panel borrow, exactly as the classic `draw_<x>_workbench`
/// functions drain them after `workbench_shell` returns, so the 3-D /
/// visualization buttons work in the dock host too.
fn drain_workbench_deferred(app: &mut ValenxApp) {
    crate::engine_workbench::drain_deferred(app);
    crate::rocket_workbench::drain_deferred(app);
    crate::fem_workbench::drain_deferred(app);
}

/// The [`egui_tiles::Behavior`] that titles + paints each dock tile and
/// wires the per-tab close button back to the workbench's `show_*` flag. It
/// borrows the whole app so [`render_panel_body`] can mutate workbench state
/// while drawing — see the borrow note on [`ValenxApp::draw_dock_layout`].
struct DockBehavior<'a> {
    app: &'a mut ValenxApp,
}

impl egui_tiles::Behavior<String> for DockBehavior<'_> {
    fn tab_title_for_pane(&mut self, pane: &String) -> egui::WidgetText {
        panel_title(pane).into()
    }

    fn pane_ui(
        &mut self,
        ui: &mut egui::Ui,
        _tile_id: egui_tiles::TileId,
        pane: &mut String,
    ) -> egui_tiles::UiResponse {
        // A little breathing room around the body, matching the central
        // docking layout's pane frame.
        egui::Frame::none()
            .inner_margin(egui::Margin::same(6.0))
            .show(ui, |ui| {
                render_panel_body(self.app, ui, pane);
            });
        // The drag handle is the tab bar (handled by egui_tiles); we never
        // start a drag from the body.
        egui_tiles::UiResponse::None
    }

    /// Every workbench tab gets a ✕ so the user can close a panel straight
    /// from the dock.
    fn is_tab_closable(
        &self,
        _tiles: &egui_tiles::Tiles<String>,
        _tile_id: egui_tiles::TileId,
    ) -> bool {
        true
    }

    /// When a tab's ✕ is clicked, clear the workbench's `show_*` flag so the
    /// close sticks (otherwise the per-frame sync would re-add the tile next
    /// frame). Returning `true` lets `egui_tiles` then remove the tile.
    fn on_tab_close(
        &mut self,
        tiles: &mut egui_tiles::Tiles<String>,
        tile_id: egui_tiles::TileId,
    ) -> bool {
        if let Some(egui_tiles::Tile::Pane(panel_id)) = tiles.get(tile_id) {
            let panel_id = panel_id.clone();
            close_panel(self.app, &panel_id);
        }
        true
    }

    /// Give **every** pane its own draggable tab bar — even a lone pane in a
    /// split. This is what makes the dock actually reorganizable by hand:
    /// without it, panes inside `Linear` (row/column) containers have only
    /// resize handles between them and **no grab handle**, so the user can't
    /// drag to reorder. With it, each pane shows a title tab you grab and
    /// drag — drop it on another pane's edge to split into a new row/column,
    /// onto a tab bar to stack, or along the strip to reorder.
    fn simplification_options(&self) -> egui_tiles::SimplificationOptions {
        egui_tiles::SimplificationOptions {
            all_panes_must_have_tabs: true,
            ..Default::default()
        }
    }
}

impl ValenxApp {
    /// Paint the opt-in dockable workbench layout: sync the tile tree to the
    /// set of currently-open workbenches, then render it inside one resizable
    /// right [`egui::SidePanel`]. Called from `update.rs` in place of the
    /// classic per-workbench `SidePanel` run when
    /// [`ValenxApp::dock_enabled`] is on.
    ///
    /// ### Borrow note
    ///
    /// [`egui_tiles::Tree::ui`] needs `&mut self.dock_tree` *and* a
    /// `Behavior` that holds `&mut self` (so [`render_panel_body`] can mutate
    /// workbench state). Those two `&mut self` borrows can't coexist, so we
    /// [`Option::take`] the tree into a local, build the behavior against
    /// `self`, draw, then put the tree back.
    pub(crate) fn draw_dock_layout(&mut self, ctx: &egui::Context) {
        // 1. Which workbenches are open right now, in the registry's order.
        let open_ids: Vec<String> = DOCKABLE_PANELS
            .iter()
            .filter(|(id, _)| is_panel_open(self, id))
            .map(|(id, _)| (*id).to_string())
            .collect();

        // Are there any "Workbench + Agent" tiles in the current tree? Those
        // are launcher-created and not gated on a `show_*` flag, so the region
        // must stay alive for them even when no `DOCKABLE_PANELS` are open.
        let has_wb_agent_tiles = self
            .dock_tree
            .as_ref()
            .map(|t| {
                t.tiles
                    .tiles()
                    .any(|tile| matches!(tile, egui_tiles::Tile::Pane(id) if is_wb_agent_pane(id)))
            })
            .unwrap_or(false);

        // 2. Nothing open *and* no Workbench+Agent tiles → drop the tree and
        //    paint nothing (the region vanishes, like the classic layout when
        //    every panel is closed).
        if open_ids.is_empty() && !has_wb_agent_tiles {
            self.dock_tree = None;
            return;
        }

        // 3. Lazily build the tree from the open set the first time we need
        //    one (a single horizontal row of panes; the user can then split
        //    it into rows/columns/tabs by dragging). The launcher buttons
        //    build their own tree, so we only reach this branch when a regular
        //    workbench was opened with the dock on but no tree yet.
        if self.dock_tree.is_none() {
            self.dock_tree = Some(egui_tiles::Tree::new_horizontal(
                "valenx_dock_tree",
                open_ids.clone(),
            ));
        }

        // 4. Sync the existing tree to the open set: drop tiles whose panel
        //    was closed, and add a tile for any panel opened since last
        //    frame. Done before drawing so the tree always matches state.
        if let Some(tree) = self.dock_tree.as_mut() {
            sync_tree(tree, &open_ids);
        }

        // 5. Draw. Take the tree out so the behavior can borrow `self`.
        if let Some(mut tree) = self.dock_tree.take() {
            egui::SidePanel::right("valenx_dock_region")
                .resizable(true)
                .default_width(700.0)
                .show(ctx, |ui| {
                    let mut beh = DockBehavior { app: self };
                    tree.ui(&mut beh, ui);
                });
            // Put it back for next frame (preserves the user's layout edits).
            self.dock_tree = Some(tree);
        }

        // 6. Outside the panel borrow: drain any 3-D / overlay requests the
        //    bodies queued this frame.
        drain_workbench_deferred(self);
    }

    /// Ensure the dock tree exists and return it mutably, creating an empty
    /// one (no root) on first use. Used by the launcher helpers below, which
    /// add their tiles into it. Turning the dock *on* is the caller's job.
    fn dock_tree_or_empty(&mut self) -> &mut egui_tiles::Tree<String> {
        self.dock_tree
            .get_or_insert_with(|| egui_tiles::Tree::empty("valenx_dock_tree"))
    }

    /// Launch **one "Workbench + Agent" unit** into the dock: a horizontal
    /// pair `[workspace:n | agent:n]` (empty build canvas on the left, Claude
    /// chat on the right). Turns the dock on, bumps
    /// [`Self::wb_agent_counter`] to the new `n`, and appends the pair as a new
    /// row of the dock's root (wrapping the existing root into a vertical
    /// container so units stack top-to-bottom). Wired to View → "New Workbench
    /// + Agent".
    pub(crate) fn add_workbench_agent_pair(&mut self) {
        self.dock_enabled = true;
        self.wb_agent_counter += 1;
        let n = self.wb_agent_counter;
        let tree = self.dock_tree_or_empty();
        let pair = insert_pair(tree, n);
        attach_unit_to_root(tree, pair);
    }

    /// Launch **six "Workbench + Agent" units in a 3×2 grid** (3 columns ×
    /// 2 rows) — the demo layout. Turns the dock on, hands out six fresh unit
    /// numbers, and *replaces* the dock root with a vertical container of two
    /// rows, each a horizontal container of three units, each unit a
    /// horizontal `[workspace | agent]` pair. Wired to View → "Open 6
    /// Workbench+Agents (3×2 grid)".
    ///
    /// Replacing (rather than appending to) the root keeps the demo's grid
    /// crisp; the user can still drag tiles around afterwards.
    pub(crate) fn open_six_workbench_agents(&mut self) {
        self.dock_enabled = true;
        // Reserve six fresh unit numbers up front (before borrowing the tree),
        // as `[[r0c0, r0c1, r0c2], [r1c0, r1c1, r1c2]]`.
        let nums: [[usize; 3]; 2] = std::array::from_fn(|_row| {
            std::array::from_fn(|_col| {
                self.wb_agent_counter += 1;
                self.wb_agent_counter
            })
        });
        let tree = self.dock_tree_or_empty();
        let rows: Vec<egui_tiles::TileId> = nums
            .iter()
            .map(|row| {
                let units: Vec<egui_tiles::TileId> =
                    row.iter().map(|&n| insert_pair(tree, n)).collect();
                tree.tiles.insert_horizontal_tile(units)
            })
            .collect();
        let grid_root = tree.tiles.insert_vertical_tile(rows);
        tree.root = Some(grid_root);
    }
}

/// Insert one "Workbench + Agent" unit's two panes (`workspace:n`, `agent:n`)
/// and wrap them in a horizontal container `[workspace | agent]`, returning
/// that container's [`egui_tiles::TileId`].
fn insert_pair(tree: &mut egui_tiles::Tree<String>, n: usize) -> egui_tiles::TileId {
    let workspace = tree.tiles.insert_pane(format!("{WORKSPACE_PREFIX}{n}"));
    let agent = tree.tiles.insert_pane(format!("{AGENT_PREFIX}{n}"));
    tree.tiles.insert_horizontal_tile(vec![workspace, agent])
}

/// Attach a freshly-built unit (`new_unit`) to the dock root as a new **row**:
///
/// - empty tree → the unit becomes the root;
/// - root is a vertical container → append the unit as another row;
/// - otherwise → wrap the old root and the unit into a fresh vertical
///   container so successive units stack top-to-bottom.
fn attach_unit_to_root(tree: &mut egui_tiles::Tree<String>, new_unit: egui_tiles::TileId) {
    match tree.root() {
        None => tree.root = Some(new_unit),
        Some(root) => {
            if let Some(egui_tiles::Tile::Container(egui_tiles::Container::Linear(linear))) =
                tree.tiles.get_mut(root)
            {
                if linear.dir == egui_tiles::LinearDir::Vertical {
                    linear.add_child(new_unit);
                    return;
                }
            }
            let new_root = tree.tiles.insert_vertical_tile(vec![root, new_unit]);
            tree.root = Some(new_root);
        }
    }
}

/// Reconcile `tree` with `open_ids`: remove any pane whose panel id is no
/// longer open, and append a new pane for any open id that has no tile yet.
/// New panes are added into the tree's root container so they reflow next to
/// the existing ones; the user's manual splits / order are otherwise left
/// untouched.
fn sync_tree(tree: &mut egui_tiles::Tree<String>, open_ids: &[String]) {
    // 4a. Remove tiles for closed panels. Collect first to avoid mutating
    //     while iterating. "Workbench + Agent" panes (`workspace:` / `agent:`)
    //     are exempt: they're launcher-created, never appear in `open_ids`,
    //     and persist until the user closes them via the tab ✕.
    let to_remove: Vec<egui_tiles::TileId> = tree
        .tiles
        .iter()
        .filter_map(|(tile_id, tile)| match tile {
            egui_tiles::Tile::Pane(panel_id)
                if !is_wb_agent_pane(panel_id) && !open_ids.contains(panel_id) =>
            {
                Some(*tile_id)
            }
            _ => None,
        })
        .collect();
    for tile_id in to_remove {
        tree.tiles.remove(tile_id);
    }

    // 4b. Which panel ids does the tree still host as panes?
    let present: std::collections::HashSet<String> = tree
        .tiles
        .tiles()
        .filter_map(|tile| match tile {
            egui_tiles::Tile::Pane(panel_id) => Some(panel_id.clone()),
            _ => None,
        })
        .collect();

    // 4c. Add a pane for any open id missing from the tree.
    let missing: Vec<String> = open_ids
        .iter()
        .filter(|id| !present.contains(*id))
        .cloned()
        .collect();
    for panel_id in missing {
        let new_pane = tree.tiles.insert_pane(panel_id);
        match tree.root() {
            // Append into the existing root container so it reflows with the
            // others.
            Some(root) => {
                if let Some(egui_tiles::Tile::Container(container)) = tree.tiles.get_mut(root) {
                    container.add_child(new_pane);
                } else {
                    // Root is a lone pane (or empty): wrap both into a fresh
                    // horizontal container so we get a real multi-pane layout.
                    let children = match tree.root() {
                        Some(existing) => vec![existing, new_pane],
                        None => vec![new_pane],
                    };
                    let new_root = tree.tiles.insert_horizontal_tile(children);
                    tree.root = Some(new_root);
                }
            }
            // Empty tree: this pane becomes the root.
            None => {
                tree.root = Some(new_pane);
            }
        }
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod tests {
    use super::*;

    #[test]
    fn registry_title_round_trips() {
        for (id, title) in DOCKABLE_PANELS {
            assert_eq!(panel_title(id), *title);
        }
        // Unknown id falls back to the id itself.
        assert_eq!(panel_title("nope"), "nope");
    }

    #[test]
    fn wb_agent_titles_and_classification() {
        // The "Workbench + Agent" tile ids title as "Workspace N" / "Agent N"
        // and are recognised as launcher panes (exempt from sync pruning).
        assert_eq!(panel_title("workspace:1"), "Workspace 1");
        assert_eq!(panel_title("agent:1"), "Agent 1");
        assert_eq!(panel_title("workspace:42"), "Workspace 42");
        assert_eq!(panel_title("agent:42"), "Agent 42");
        assert!(is_wb_agent_pane("workspace:3"));
        assert!(is_wb_agent_pane("agent:3"));
        assert!(!is_wb_agent_pane("valenx_fem_workbench"));
        assert!(!is_wb_agent_pane("nope"));
    }

    #[test]
    fn every_registry_panel_has_open_close_and_render_wiring() {
        // is_panel_open / close_panel must recognise every registry id (no
        // `_ => false` / `_ => {}` swallowing a real panel), and the default
        // app has them all closed.
        let mut app = ValenxApp::default();
        for (id, _) in DOCKABLE_PANELS {
            assert!(!is_panel_open(&app, id), "{id} should default closed");
            // Flip it on via the matching show flag, confirm is_panel_open
            // observes it, then close_panel turns it back off.
            set_open(&mut app, id, true);
            assert!(is_panel_open(&app, id), "{id} open flag not observed");
            close_panel(&mut app, id);
            assert!(!is_panel_open(&app, id), "{id} close_panel failed");
        }
    }

    /// Test helper: set a registry panel's show flag directly.
    fn set_open(app: &mut ValenxApp, panel_id: &str, on: bool) {
        match panel_id {
            "valenx_mesh_toolbox" => app.show_mesh_toolbox = on,
            "valenx_genetics_workbench" => app.show_genetics_workbench = on,
            "valenx_aero_workbench" => app.show_aero_workbench = on,
            "valenx_fem_workbench" => app.show_fem_workbench = on,
            "valenx_cfd_workbench" => app.show_cfd_workbench = on,
            "valenx_astro_workbench" => app.show_astro_workbench = on,
            "valenx_rocket_workbench" => app.show_rocket_workbench = on,
            "valenx_engine_workbench" => app.show_engine_workbench = on,
            "valenx_car_workbench" => app.show_car_workbench = on,
            "valenx_assistant_panel" => app.show_assistant_panel = on,
            other => panic!("unhandled id {other}"),
        }
    }

    #[test]
    fn sync_tree_adds_and_removes_panes_to_match_open_set() {
        // Start with two open panes.
        let mut tree = egui_tiles::Tree::new_horizontal(
            "t",
            vec![
                "valenx_fem_workbench".to_string(),
                "valenx_cfd_workbench".to_string(),
            ],
        );
        let count_panes = |t: &egui_tiles::Tree<String>| {
            t.tiles
                .tiles()
                .filter(|tile| matches!(tile, egui_tiles::Tile::Pane(_)))
                .count()
        };
        assert_eq!(count_panes(&tree), 2);

        // Close one, open a new one → still two, but the membership changed.
        let open = vec![
            "valenx_cfd_workbench".to_string(),
            "valenx_engine_workbench".to_string(),
        ];
        sync_tree(&mut tree, &open);
        let present: std::collections::HashSet<String> = tree
            .tiles
            .tiles()
            .filter_map(|tile| match tile {
                egui_tiles::Tile::Pane(id) => Some(id.clone()),
                _ => None,
            })
            .collect();
        assert_eq!(present.len(), 2);
        assert!(present.contains("valenx_cfd_workbench"));
        assert!(present.contains("valenx_engine_workbench"));
        assert!(!present.contains("valenx_fem_workbench"));
    }

    /// Test helper: the set of pane (leaf) ids currently in a tree.
    fn pane_ids(tree: &egui_tiles::Tree<String>) -> std::collections::HashSet<String> {
        tree.tiles
            .tiles()
            .filter_map(|tile| match tile {
                egui_tiles::Tile::Pane(id) => Some(id.clone()),
                _ => None,
            })
            .collect()
    }

    #[test]
    fn add_workbench_agent_pair_builds_workspace_and_agent_tiles() {
        // "New Workbench + Agent" turns the dock on, bumps the counter, and
        // adds exactly the paired workspace:1 / agent:1 panes.
        let mut app = ValenxApp::default();
        assert_eq!(app.wb_agent_counter, 0);
        app.add_workbench_agent_pair();
        assert!(app.dock_enabled, "launcher must enable the dock");
        assert_eq!(app.wb_agent_counter, 1);
        let tree = app.dock_tree.as_ref().expect("pair built a tree");
        let panes = pane_ids(tree);
        assert!(panes.contains("workspace:1"));
        assert!(panes.contains("agent:1"));
        assert_eq!(panes.len(), 2);
        assert!(tree.root().is_some());

        // A second pair stacks on without disturbing the first.
        app.add_workbench_agent_pair();
        assert_eq!(app.wb_agent_counter, 2);
        let panes = pane_ids(app.dock_tree.as_ref().unwrap());
        for id in ["workspace:1", "agent:1", "workspace:2", "agent:2"] {
            assert!(panes.contains(id), "{id} missing after second pair");
        }
        assert_eq!(panes.len(), 4);
    }

    #[test]
    fn open_six_workbench_agents_builds_a_3x2_grid_of_pairs() {
        // "Open 6 Workbench+Agents (3x2 grid)" yields six units = twelve panes
        // (workspace/agent for n=1..=6), arranged as a vertical-of-2-rows,
        // each row a horizontal-of-3-units, each unit a horizontal pair.
        let mut app = ValenxApp::default();
        app.open_six_workbench_agents();
        assert!(app.dock_enabled);
        assert_eq!(app.wb_agent_counter, 6);

        let tree = app.dock_tree.as_ref().expect("grid built a tree");
        let panes = pane_ids(tree);
        assert_eq!(panes.len(), 12, "6 units x 2 panes each");
        for n in 1..=6 {
            assert!(panes.contains(&format!("workspace:{n}")));
            assert!(panes.contains(&format!("agent:{n}")));
        }

        // Structure check: root = vertical container with 2 row-children,
        // each row = horizontal container with 3 unit-children, each unit =
        // horizontal pair of 2 panes.
        let root = tree.root().expect("grid has a root");
        let egui_tiles::Tile::Container(egui_tiles::Container::Linear(vroot)) =
            tree.tiles.get(root).unwrap()
        else {
            panic!("root must be a Linear container");
        };
        assert_eq!(vroot.dir, egui_tiles::LinearDir::Vertical);
        let rows: Vec<egui_tiles::TileId> = vroot.children.clone();
        assert_eq!(rows.len(), 2, "two rows");
        for row in rows {
            let egui_tiles::Tile::Container(egui_tiles::Container::Linear(hrow)) =
                tree.tiles.get(row).unwrap()
            else {
                panic!("each row must be a horizontal Linear container");
            };
            assert_eq!(hrow.dir, egui_tiles::LinearDir::Horizontal);
            assert_eq!(hrow.children.len(), 3, "three units per row");
            for unit in &hrow.children {
                let egui_tiles::Tile::Container(egui_tiles::Container::Linear(pair)) =
                    tree.tiles.get(*unit).unwrap()
                else {
                    panic!("each unit must be a horizontal pair");
                };
                assert_eq!(pair.dir, egui_tiles::LinearDir::Horizontal);
                assert_eq!(pair.children.len(), 2, "workspace + agent");
            }
        }
    }

    #[test]
    fn wb_agent_tiles_survive_sync_and_dock_keeps_region_alive() {
        // sync_tree must NOT prune workspace:/agent: panes (they're never in
        // the open set), and draw_dock_layout must keep the region alive while
        // they exist even with no DOCKABLE_PANELS open. Drive it headless.
        let mut app = ValenxApp::default();
        app.open_six_workbench_agents();
        let ctx = egui::Context::default();
        // Several frames: build → sync (with empty open_ids) → persist.
        for _ in 0..3 {
            let _ = ctx.run(egui::RawInput::default(), |ctx| {
                app.draw_dock_layout(ctx);
            });
        }
        let tree = app
            .dock_tree
            .as_ref()
            .expect("region must persist for WB+Agent tiles with nothing else open");
        let panes = pane_ids(tree);
        assert_eq!(panes.len(), 12, "all six pairs survived the sync");
        assert!(panes.contains("workspace:1"));
        assert!(panes.contains("agent:6"));

        // Opening a regular workbench alongside adds its pane without dropping
        // the WB+Agent tiles.
        app.show_assistant_panel = true;
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            app.draw_dock_layout(ctx);
        });
        let panes = pane_ids(app.dock_tree.as_ref().unwrap());
        assert!(panes.contains("valenx_assistant_panel"));
        assert!(panes.contains("workspace:1"));
        assert_eq!(panes.len(), 13);
    }

    #[test]
    fn dock_layout_is_a_noop_when_nothing_open() {
        // With every workbench closed, draw_dock_layout must clear the tree
        // and paint nothing without panicking.
        let mut app = ValenxApp::default();
        app.dock_enabled = true;
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            app.draw_dock_layout(ctx);
        });
        assert!(app.dock_tree.is_none());
    }

    #[test]
    fn dock_layout_builds_tree_and_renders_open_panels_without_panic() {
        // Open several workbenches, then draw the dock layout headless. This
        // exercises the lazy tree build, the per-frame sync, the
        // take()/put-back borrow dance, and pane_ui → render_panel_body for
        // each wired body — none of which may panic, and the tree must
        // persist across frames carrying exactly the open panes.
        let mut app = ValenxApp::default();
        app.dock_enabled = true;
        app.show_engine_workbench = true;
        app.show_fem_workbench = true;
        app.show_assistant_panel = true;

        let ctx = egui::Context::default();
        // Two frames: first builds, second hits the sync + persisted-tree path.
        for _ in 0..2 {
            let _ = ctx.run(egui::RawInput::default(), |ctx| {
                app.draw_dock_layout(ctx);
            });
        }
        let tree = app.dock_tree.as_ref().expect("tree built for open panels");
        let panes: std::collections::HashSet<String> = tree
            .tiles
            .tiles()
            .filter_map(|tile| match tile {
                egui_tiles::Tile::Pane(id) => Some(id.clone()),
                _ => None,
            })
            .collect();
        assert_eq!(panes.len(), 3);
        assert!(panes.contains("valenx_engine_workbench"));
        assert!(panes.contains("valenx_fem_workbench"));
        assert!(panes.contains("valenx_assistant_panel"));

        // Close one workbench; the next dock frame must drop its pane.
        app.show_fem_workbench = false;
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            app.draw_dock_layout(ctx);
        });
        let tree = app.dock_tree.as_ref().expect("tree still has open panels");
        let panes: std::collections::HashSet<String> = tree
            .tiles
            .tiles()
            .filter_map(|tile| match tile {
                egui_tiles::Tile::Pane(id) => Some(id.clone()),
                _ => None,
            })
            .collect();
        assert!(!panes.contains("valenx_fem_workbench"));
        assert!(panes.contains("valenx_engine_workbench"));
    }
}
