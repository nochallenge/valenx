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
use eframe::egui_wgpu;

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

/// Where a **new** "Workbench + Agent" unit should be inserted into the dock
/// grid, chosen from the tab-strip "+ Workbench+Agent" dropdown. Lets the user
/// place the new unit precisely rather than always tacking it onto a new bottom
/// row.
///
/// The grid is a vertical stack of full-width **rows**, each row a horizontal
/// strip of unit pairs. `RowStart`/`RowEnd` carry a **0-based** row index into
/// that stack and add the new unit at the left / right end *within* that row;
/// `NewRowTop`/`NewRowBottom` add the unit as a brand-new first / last row.
///
/// Consumed by [`ValenxApp::add_workbench_agent_pair_at`]. The menu that
/// produces it is built from [`ValenxApp::dock_grid_rows`].
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum UnitAddTarget {
    /// Add the new unit as a brand-new **first** row of the grid.
    NewRowTop,
    /// Add the new unit as a brand-new **last** row of the grid.
    NewRowBottom,
    /// Add the new unit at the **left** end of the existing row at this
    /// 0-based index (falling back to a new bottom row if out of range).
    RowStart(usize),
    /// Add the new unit at the **right** end of the existing row at this
    /// 0-based index (falling back to a new bottom row if out of range).
    RowEnd(usize),
}

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
///
/// `wgpu_renderer` / `render_state` / `pixels_per_point` are only used by the
/// `workspace:<n>` branch (live 3-D render); every other body ignores them.
pub(crate) fn render_panel_body(
    app: &mut ValenxApp,
    ui: &mut egui::Ui,
    panel_id: &str,
    wgpu_renderer: &mut Option<crate::wgpu_renderer::WgpuRenderer>,
    render_state: Option<&egui_wgpu::RenderState>,
    pixels_per_point: f32,
) {
    // "Workbench + Agent" tiles: the agent half is an **independent** Claude
    // chat keyed by the unit number (its own feed file + input buffer — see
    // `assistant_chat_ui`'s `Some(n)` channel); the workspace half is an empty
    // build canvas. Whole-unit relocation is via egui_tiles' per-tile tab drag.
    if let Some(n) = panel_id.strip_prefix(AGENT_PREFIX) {
        if let Ok(n) = n.parse::<usize>() {
            crate::assistant_workbench::assistant_chat_ui(app, ui, Some(n));
        } else {
            // Unparseable unit number (shouldn't happen for ids we insert) →
            // fall back to the classic shared channel rather than panicking.
            crate::assistant_workbench::assistant_chat_ui(app, ui, None);
        }
        return;
    }
    if let Some(n) = panel_id.strip_prefix(WORKSPACE_PREFIX) {
        render_workspace_body(app, ui, n, wgpu_renderer, render_state, pixels_per_point);
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

/// Render a `"workspace:<n>"` tile body — the **empty build-output canvas**
/// half of a "Workbench + Agent" unit, where `n` is the unit number (the suffix
/// after the `"workspace:"` prefix).
///
/// The workspace is the **agent's build-output area**, not a chat mirror: it is
/// deliberately quiet — a small subtle `"Workspace N"` title above a clean
/// bordered area. Once the external agent posts a finished result for this unit
/// (a [`crate::WorkspaceProduct`] via
/// [`crate::agent_commands::AgentCommand::ShowProduct`], keyed by the unit
/// number `n`), that result renders into the bordered area as a **result card**
/// — a bold title heading over one row per line. Until then it shows a faint
/// centered hint. No chat echo and no "move whole unit" header here (that
/// control lives on the agent half).
fn render_workspace_body(
    app: &mut ValenxApp,
    ui: &mut egui::Ui,
    n: &str,
    wgpu_renderer: &mut Option<crate::wgpu_renderer::WgpuRenderer>,
    render_state: Option<&egui_wgpu::RenderState>,
    pixels_per_point: f32,
) {
    ui.label(egui::RichText::new(format!("Workspace {n}")).weak());
    ui.add_space(4.0);
    // This unit's posted result is keyed by its numeric id (the same `n` the
    // agent bridge uses). A non-numeric suffix (shouldn't happen for ids we
    // insert) simply finds nothing → placeholder. We deliberately do NOT bind
    // the product to a long-lived immutable borrow here: the 3-D branch below
    // mutates the product's camera (`get_mut`), so the lookups are kept short.
    let idx = n.parse::<usize>().ok();

    // A live 3-D product (a `show_3d` command set `mesh: Some`) renders as an
    // actual lit viewport — same look as the central viewport. The tile id keys
    // a dedicated per-tile offscreen target so it never aliases the central
    // viewport's or another tile's. If anything required is missing (no GPU
    // ctx, no mesh) it falls through to the text card / placeholder.
    let has_3d = idx
        .and_then(|i| app.workspace_products.get(&i))
        .map(|p| p.mesh.is_some())
        .unwrap_or(false);
    if has_3d {
        // Borrow plan: the offscreen render reads the product's mesh through an
        // immutable `app.workspace_products.get(..)` borrow, but afterwards we
        // mutate that same product's camera via `get_mut`. The two borrows can't
        // overlap, so:
        //   1. clone the camera into an owned `camera` (it's `Clone`, not
        //      `Copy`) — the value we orbit / zoom / frame this frame;
        //   2. run the render inside a scoped immutable borrow of the product
        //      (for its mesh), evaluating to the *owned* `Response` + AABB so no
        //      borrow escapes the block;
        //   3. apply the input deltas to the owned `camera`, then re-borrow
        //      mutably (`get_mut`) and commit it back.
        // `wgpu_renderer` / `render_state` were lifted out of `self` (see the
        // `DockBehavior` doc), so they don't alias `workspace_products` — the
        // mesh borrow and the renderer `&mut` coexist fine.
        let camera = idx
            .and_then(|i| app.workspace_products.get(&i))
            .map(|p| p.camera.clone());
        if let (Some(mut camera), Some(renderer), Some(rs), Some(i)) =
            (camera, wgpu_renderer.as_mut(), render_state, idx)
        {
            let tile_key = format!("workspace:{n}");
            // Scope the mesh borrow to the render call only; it yields owned
            // values (`Response`, AABB) that outlive the borrow.
            let drawn = match app.workspace_products.get(&i).and_then(|p| p.mesh.as_ref()) {
                Some(mesh) => render_tile_mesh_3d(
                    ui,
                    &tile_key,
                    mesh,
                    &camera,
                    renderer,
                    rs,
                    pixels_per_point,
                ),
                None => None,
            };
            if let Some((response, aabb)) = drawn {
                // Mirror the central viewport's input blocks (viewport.rs::show)
                // against this pane's camera clone, then commit it back.
                let mut changed = false;
                if response.dragged_by(egui::PointerButton::Primary)
                    || response.dragged_by(egui::PointerButton::Middle)
                {
                    let delta = response.drag_delta();
                    camera.orbit(delta.x * 0.5, -delta.y * 0.5);
                    changed = true;
                }
                if response.hovered() {
                    let scroll = ui.input(|i| i.raw_scroll_delta.y);
                    if scroll.abs() > f32::EPSILON {
                        camera.zoom(scroll * 0.01);
                        changed = true;
                    }
                }
                if response.double_clicked() {
                    if let Some((min, max)) = aabb {
                        camera.frame_bounds(min, max);
                        changed = true;
                    }
                }
                if changed {
                    if let Some(p) = app.workspace_products.get_mut(&i) {
                        p.camera = camera;
                    }
                    ui.ctx().request_repaint();
                }
                return;
            }
        }
    }

    // Fall-through: text card / placeholder. Re-borrow the product fresh (a
    // short immutable borrow) — it was deliberately not held across the 3-D
    // branch's camera mutation above.
    let product = idx.and_then(|i| app.workspace_products.get(&i));
    egui::Frame::group(ui.style())
        .inner_margin(egui::Margin::same(8.0))
        .show(ui, |ui| {
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| match product {
                    Some(product) => {
                        // Result card: bold title heading, then one row per line.
                        ui.heading(&product.title);
                        if !product.lines.is_empty() {
                            ui.add_space(4.0);
                            for line in &product.lines {
                                ui.label(line);
                            }
                        }
                    }
                    None => {
                        ui.centered_and_justified(|ui| {
                            ui.label(
                                egui::RichText::new("the agent's output will appear here")
                                    .weak()
                                    .italics(),
                            );
                        });
                    }
                });
        });
}

/// An axis-aligned bounding box as `(min, max)` corner coordinates — the shape
/// [`crate::viewport::mesh_aabb`] returns, used here for double-click-to-frame.
type MeshAabb = ([f32; 3], [f32; 3]);

/// Render `mesh` as a lit 3-D view filling the tile's remaining rect, from
/// `camera`, into the per-tile offscreen target keyed by `tile_key`, and blit
/// it. Mirrors [`crate::viewport`]'s `render_wgpu_scene` (the central
/// viewport's path) but: keys its own offscreen target (so two viewports per
/// frame don't alias) and draws **grid-less** (background + shaded mesh only)
/// for a clean, cheap mini-view. Allocates a `click_and_drag` rect (Stage 2 —
/// the caller reads the returned [`egui::Response`] to orbit / zoom / frame
/// this pane's own camera). Returns `Some((response, aabb))` if it drew — the
/// AABB is the mesh's bounds for double-click-to-frame — or `None` if the rect
/// was degenerate or the GPU returned no texture (caller falls back to the
/// text card).
fn render_tile_mesh_3d(
    ui: &mut egui::Ui,
    tile_key: &str,
    mesh: &crate::types::LoadedMesh,
    camera: &valenx_viz::OrbitCamera,
    renderer: &mut crate::wgpu_renderer::WgpuRenderer,
    render_state: &egui_wgpu::RenderState,
    pixels_per_point: f32,
) -> Option<(egui::Response, Option<MeshAabb>)> {
    // Build the shaded surface for the mesh's surface elements (Tri3 / Quad4).
    // A volumetric-only mesh yields nothing → fall back to the card.
    let surface = crate::viewport::mesh_to_triangle_surface(&mesh.mesh)?;
    let vertices = crate::wgpu_renderer::triangles_to_vertices(&surface);
    if vertices.is_empty() {
        return None;
    }

    // Take the whole remaining tile rect; `click_and_drag` so the caller can
    // orbit (drag), zoom (scroll while hovered) and frame (double-click) this
    // pane's camera — same interaction model as the central viewport.
    let avail = ui.available_size();
    let (rect, response) = ui.allocate_exact_size(
        egui::vec2(avail.x.max(16.0), avail.y.max(16.0)),
        egui::Sense::click_and_drag(),
    );
    let w_logical = rect.width();
    let h_logical = rect.height();
    if w_logical < 1.0 || h_logical < 1.0 {
        return None;
    }
    let ppp = pixels_per_point.max(0.5);
    let size_px = [(w_logical * ppp) as u32, (h_logical * ppp) as u32];

    // Camera → MVP / inverse-MVP for this tile's aspect ratio, matching the
    // central viewport's matrix construction exactly.
    let mvp = crate::wgpu_renderer::mvp_from_camera(camera, w_logical, h_logical);
    let inv_mvp = crate::wgpu_renderer::inv_mvp_from_camera(camera, w_logical, h_logical);
    let light_dir = [-0.3f32, -0.8, -0.5];
    let eye = camera.eye();
    let cam_pos = [eye.x, eye.y, eye.z];
    // Grid params are still supplied (the uniform layout is shared) but the
    // grid pass is disabled below via `draw_grid = false`, so they're inert.
    let (minor, blend_t) = valenx_viz::grid_lod_params(camera.distance);
    let grid = [minor, 0.0, 0.0, camera.distance * 14.0];
    let grid2 = [blend_t, 0.0, 0.0, 0.0];

    let mut render_guard = render_state.renderer.write();
    let texture_id = renderer.render_keyed(
        tile_key,
        &mut render_guard,
        size_px,
        mvp,
        inv_mvp,
        light_dir,
        cam_pos,
        grid,
        grid2,
        false, // draw_grid = false: background + shaded mesh only (mini-view)
        &vertices,
    );
    drop(render_guard);

    let texture_id = texture_id?;
    let uv = egui::Rect::from_min_max(egui::pos2(0.0, 0.0), egui::pos2(1.0, 1.0));
    ui.painter_at(rect)
        .image(texture_id, rect, uv, egui::Color32::WHITE);
    // Mesh bounds for the caller's double-click-to-frame.
    let aabb = crate::viewport::mesh_aabb(&mesh.mesh);
    Some((response, aabb))
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
///
/// It also carries the wgpu bits a `workspace:<n>` tile needs to render a
/// live 3-D model: `wgpu_renderer` is the app's renderer *lifted out of
/// `self`* (so it can be held alongside `&mut app` without aliasing — see
/// [`ValenxApp::render_dock_tree_into`]), plus the frame's `render_state`
/// and `pixels_per_point`.
struct DockBehavior<'a> {
    app: &'a mut ValenxApp,
    /// The app's wgpu renderer, taken into a local for the duration of the
    /// draw. `None` when eframe has no wgpu backend → tiles fall back to text.
    wgpu_renderer: &'a mut Option<crate::wgpu_renderer::WgpuRenderer>,
    /// The frame's wgpu render state (egui texture registry + device/queue).
    render_state: Option<&'a egui_wgpu::RenderState>,
    /// Logical-to-physical pixel ratio, for sizing the offscreen target.
    pixels_per_point: f32,
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
        // docking layout's pane frame. Borrow the behavior's fields
        // disjointly so the renderer/render-state can ride along to a
        // `workspace:<n>` tile while `app` is mutably borrowed for the body.
        let renderer = &mut *self.wgpu_renderer;
        let render_state = self.render_state;
        let ppp = self.pixels_per_point;
        egui::Frame::none()
            .inner_margin(egui::Margin::same(6.0))
            .show(ui, |ui| {
                render_panel_body(self.app, ui, pane, renderer, render_state, ppp);
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
    ///
    /// `render_state` + `pixels_per_point` are the wgpu bits a
    /// `workspace:<n>` tile needs to render its live 3-D model; they thread
    /// through to [`Self::render_dock_tree_into`]. `None`/any value is fine
    /// when no GPU backend exists — the tile falls back to its text card.
    pub(crate) fn draw_dock_layout(
        &mut self,
        ctx: &egui::Context,
        render_state: Option<&egui_wgpu::RenderState>,
        pixels_per_point: f32,
    ) {
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

        // 5. Host the dock. With the 3-D viewport visible the dock lives in a
        //    resizable right SidePanel beside it (classic layout). When the
        //    user hides the viewport (its ✕ / View → "Hide 3D viewport"), the
        //    CentralPanel hosts the dock **full-width** instead — so e.g. a
        //    Workbench+Agent grid uses the entire workspace — and we render
        //    nothing here. The CentralPanel in `update.rs` calls
        //    [`Self::render_dock_tree_into`] in that case.
        if !self.viewport_hidden {
            egui::SidePanel::right("valenx_dock_region")
                .resizable(true)
                .default_width(700.0)
                .show(ctx, |ui| {
                    self.render_dock_tree_into(ui, render_state, pixels_per_point);
                });
        }
    }

    /// Render the dock tile-tree into `ui`: the actual draw, the
    /// [`Option::take`]/put-back borrow dance (so [`DockBehavior`] can hold
    /// `&mut self` while [`render_panel_body`] mutates workbench state), and
    /// the post-draw deferred-work drain. Host-agnostic — called both from the
    /// right SidePanel (viewport visible) and from the CentralPanel (viewport
    /// hidden → the dock fills the whole workspace).
    ///
    /// `render_state` + `pixels_per_point` let a `workspace:<n>` tile build a
    /// [`crate::viewport::WgpuCtx`] and render its live 3-D model. The
    /// **same** `Option::take` trick that frees `self.dock_tree` is applied to
    /// `self.wgpu_renderer`: it's lifted into a local so `DockBehavior` can
    /// borrow `&mut self` (for workbench-state mutation) and the renderer
    /// **simultaneously** without aliasing `self`. The renderer is restored
    /// after the draw.
    pub(crate) fn render_dock_tree_into(
        &mut self,
        ui: &mut egui::Ui,
        render_state: Option<&egui_wgpu::RenderState>,
        pixels_per_point: f32,
    ) {
        if let Some(mut tree) = self.dock_tree.take() {
            // Lift the renderer out of `self` too, so the behavior can hold
            // both `&mut self` and `&mut renderer` without two borrows of
            // `self`. Put back below regardless of the draw outcome.
            let mut renderer = self.wgpu_renderer.take();
            let mut beh = DockBehavior {
                app: self,
                wgpu_renderer: &mut renderer,
                render_state,
                pixels_per_point,
            };
            tree.ui(&mut beh, ui);
            // Restore the renderer (it carries the per-tile offscreen targets,
            // so it MUST persist across frames) and the tree (preserves the
            // user's layout edits).
            self.wgpu_renderer = renderer;
            self.dock_tree = Some(tree);
        }
        // Drain any 3-D / overlay requests the bodies queued this frame.
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
        // Same full-workspace behaviour as before; just route through the
        // precise-placement path with the historical default (new bottom row).
        self.add_workbench_agent_pair_at(UnitAddTarget::NewRowBottom);
    }

    /// Read-only snapshot of the dock grid's **row shape**: one entry per row
    /// of the vertical grid, each the number of "Workbench + Agent" **units**
    /// in that row. Used to build the tab-strip "+ Workbench+Agent" placement
    /// dropdown so the menu reflects the live layout.
    ///
    /// The launcher builds each **unit** as a horizontal pair
    /// `[workspace | agent]` (a container) and a multi-unit **row** as a
    /// horizontal `Linear` whose children are those unit-containers; a *single*
    /// unit row is just the pair itself (its children are the two panes). So a
    /// row's unit count is the number of its **container** children, except a
    /// row whose children are bare panes (it *is* one pair) counts as `1` — see
    /// [`row_unit_count`]:
    /// - root is a **vertical** `Linear` → one unit-count per row-child;
    /// - root is **not** a vertical `Linear` (a single row — one unit, or a
    ///   horizontal/tab arrangement) → a one-element vec with that row's count;
    /// - no tree / no root → empty vec.
    ///
    /// Caveat: "unit" assumes the launcher's pair shape. After heavy manual
    /// dragging a row may hold arbitrary tiles, so a count is the best honest
    /// reading of how many side-by-side groups a row has — not a guarantee
    /// every group is a genuine `workspace:`/`agent:` pair.
    pub(crate) fn dock_grid_rows(&self) -> Vec<usize> {
        let Some(tree) = self.dock_tree.as_ref() else {
            return Vec::new();
        };
        let Some(root) = tree.root() else {
            return Vec::new();
        };
        match tree.tiles.get(root) {
            Some(egui_tiles::Tile::Container(egui_tiles::Container::Linear(vroot)))
                if vroot.dir == egui_tiles::LinearDir::Vertical =>
            {
                vroot
                    .children
                    .iter()
                    .map(|&row| row_unit_count(tree, row))
                    .collect()
            }
            // Root isn't a vertical grid: the whole thing is one row.
            _ => vec![row_unit_count(tree, root)],
        }
    }

    /// Launch **one "Workbench + Agent" unit** into the dock at a *chosen*
    /// position (`target`) — the precise-placement counterpart of
    /// [`Self::add_workbench_agent_pair`] (which is now just
    /// `NewRowBottom`). Turns the dock on full-workspace (hides the 3-D
    /// viewport), bumps [`Self::wb_agent_counter`] to the new `n`, builds the
    /// `[workspace:n | agent:n]` pair, then places it per `target`:
    ///
    /// - [`UnitAddTarget::NewRowBottom`] → append as a new last row
    ///   ([`attach_unit_to_root`]).
    /// - [`UnitAddTarget::NewRowTop`] → prepend as a new first row
    ///   (`children.insert(0, …)` when the root is a vertical `Linear`, else
    ///   wrap `[pair, old_root]` into a fresh vertical container and make it the
    ///   root).
    /// - [`UnitAddTarget::RowEnd(i)`] / [`UnitAddTarget::RowStart(i)`] → add the
    ///   pair at the right / left end **within** the existing `i`-th row: if
    ///   that row is a horizontal `Linear`, push / insert into its children; if
    ///   the row is a lone pane / non-horizontal tile, wrap `[row, pair]`
    ///   (RowEnd) or `[pair, row]` (RowStart) into a fresh horizontal container
    ///   and swap it in for that row-child. If the root isn't a vertical
    ///   `Linear`, or `i` is out of range, fall back to `NewRowBottom`.
    ///
    /// Pane ids are positional labels (the agent tiles share one chat bridge),
    /// so insertion never disturbs other units' state.
    pub(crate) fn add_workbench_agent_pair_at(&mut self, target: UnitAddTarget) {
        self.dock_enabled = true;
        // A Workbench+Agent grid is a full-workspace mode: hide the 3-D
        // viewport so the dock fills the whole central area (restore it from
        // the dock's "Show 3D viewport" bar or the View menu).
        self.viewport_hidden = true;
        self.wb_agent_counter += 1;
        let n = self.wb_agent_counter;
        let tree = self.dock_tree_or_empty();
        let pair = insert_pair(tree, n);
        place_unit(tree, pair, target);
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
        // Fill the whole workspace with the grid: hide the 3-D viewport (the
        // dock's "Show 3D viewport" bar / View menu restores it).
        self.viewport_hidden = true;
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

    /// Close the entire dock: drop every tile (workbenches + Workbench+Agent
    /// units) and restore the 3-D viewport. Wired to the dock workspace bar's
    /// **Close all**. Also clears the flag-gated dockable workbenches so
    /// [`sync_tree`] won't re-add them on the next frame.
    pub(crate) fn clear_dock(&mut self) {
        self.dock_tree = None;
        self.viewport_hidden = false;
        for (id, _) in DOCKABLE_PANELS {
            close_panel(self, id);
        }
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

/// Reset a [`egui_tiles::Container::Linear`] container's per-child **shares** to
/// equal, so its children auto-size the same. Used right after a unit is added
/// into a row (or a row into the vertical root), so a freshly-added unit/row
/// gets an equal slice rather than squeezing in at the default `1.0` beside
/// children the user may have resized.
///
/// Mechanism (egui_tiles 0.9.1): [`egui_tiles::Linear::shares`] is a public
/// [`egui_tiles::Shares`] map; its `Index`/`IndexMut` impls **default any
/// missing child to `1.0`**, and `Shares::split` treats an all-equal (or empty)
/// map as an even split. So replacing the map with `Shares::default()` (empty)
/// makes every child weigh an equal `1.0`. A no-op if `id` isn't a `Linear`.
///
/// Deliberately only called **on add** — it doesn't run on manual drags, so a
/// user's hand-tuned split between two adds is preserved until the next add
/// touches that same container.
fn equalize_shares(tree: &mut egui_tiles::Tree<String>, id: egui_tiles::TileId) {
    if let Some(egui_tiles::Tile::Container(egui_tiles::Container::Linear(linear))) =
        tree.tiles.get_mut(id)
    {
        linear.shares = egui_tiles::Shares::default();
    }
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
                    // Re-equalize the row shares so every row (incl. the new
                    // one) is the same height.
                    equalize_shares(tree, root);
                    return;
                }
            }
            let new_root = tree.tiles.insert_vertical_tile(vec![root, new_unit]);
            tree.root = Some(new_root);
        }
    }
}

/// Does the tile `row` host **container** children (i.e. it's a horizontal
/// strip of unit-pairs we can add a sibling unit *into*), as opposed to a lone
/// pair whose children are the two `workspace:`/`agent:` panes? Returns `true`
/// only for a horizontal `Linear` at least one of whose children is itself a
/// container. A bare pair (children are panes) returns `false` so placement
/// wraps it instead of mixing a nested pair in beside loose panes.
fn row_is_unit_strip(tree: &egui_tiles::Tree<String>, row: egui_tiles::TileId) -> bool {
    match tree.tiles.get(row) {
        Some(egui_tiles::Tile::Container(egui_tiles::Container::Linear(hrow)))
            if hrow.dir == egui_tiles::LinearDir::Horizontal =>
        {
            hrow.children
                .iter()
                .any(|&c| matches!(tree.tiles.get(c), Some(egui_tiles::Tile::Container(_))))
        }
        _ => false,
    }
}

/// How many "Workbench + Agent" units does the tile `row` represent? A
/// horizontal strip of unit-containers ([`row_is_unit_strip`]) counts its
/// container children; anything else — a lone pair, a single pane, a tab group,
/// a nested split — is one unit. See [`ValenxApp::dock_grid_rows`].
fn row_unit_count(tree: &egui_tiles::Tree<String>, row: egui_tiles::TileId) -> usize {
    if row_is_unit_strip(tree, row) {
        match tree.tiles.get(row) {
            Some(egui_tiles::Tile::Container(egui_tiles::Container::Linear(hrow))) => hrow
                .children
                .iter()
                .filter(|&&c| matches!(tree.tiles.get(c), Some(egui_tiles::Tile::Container(_))))
                .count(),
            _ => 1,
        }
    } else {
        1
    }
}

/// Place an already-built unit (`new_unit`) into the grid at `target`.
///
/// The placement rules are documented on
/// [`ValenxApp::add_workbench_agent_pair_at`]; this is the tree surgery. Any
/// case that can't be honoured exactly (root not a vertical grid, row index out
/// of range) falls back to attaching the unit as a new bottom row via
/// [`attach_unit_to_root`], so the unit is never dropped.
fn place_unit(
    tree: &mut egui_tiles::Tree<String>,
    new_unit: egui_tiles::TileId,
    target: UnitAddTarget,
) {
    match target {
        // New last row — the historical default.
        UnitAddTarget::NewRowBottom => attach_unit_to_root(tree, new_unit),
        // New first row — prepend into the vertical root (or wrap it).
        UnitAddTarget::NewRowTop => match tree.root() {
            None => tree.root = Some(new_unit),
            Some(root) => {
                let is_vertical_root = matches!(
                    tree.tiles.get(root),
                    Some(egui_tiles::Tile::Container(egui_tiles::Container::Linear(l)))
                        if l.dir == egui_tiles::LinearDir::Vertical
                );
                if is_vertical_root {
                    if let Some(egui_tiles::Tile::Container(egui_tiles::Container::Linear(
                        linear,
                    ))) = tree.tiles.get_mut(root)
                    {
                        linear.children.insert(0, new_unit);
                    }
                    // New top row → re-equalize the vertical root's row heights.
                    equalize_shares(tree, root);
                } else {
                    let new_root = tree.tiles.insert_vertical_tile(vec![new_unit, root]);
                    tree.root = Some(new_root);
                }
            }
        },
        // Into an existing row, at the left (RowStart) or right (RowEnd) end.
        UnitAddTarget::RowStart(i) | UnitAddTarget::RowEnd(i) => {
            let at_start = matches!(target, UnitAddTarget::RowStart(_));
            // The i-th row is only well-defined under a vertical Linear root.
            let row_id = tree.root().and_then(|root| match tree.tiles.get(root) {
                Some(egui_tiles::Tile::Container(egui_tiles::Container::Linear(vroot)))
                    if vroot.dir == egui_tiles::LinearDir::Vertical =>
                {
                    vroot.children.get(i).copied()
                }
                _ => None,
            });
            let Some(row_id) = row_id else {
                // No clean vertical grid, or i out of range → new bottom row.
                attach_unit_to_root(tree, new_unit);
                return;
            };
            // Is that row a horizontal strip of unit-containers we can add a
            // sibling unit *into* directly? (A lone pair — children are panes —
            // is NOT: adding into it would mix a nested pair beside loose panes,
            // so we wrap it instead, keeping the existing pair whole as one
            // unit beside the new one.)
            if row_is_unit_strip(tree, row_id) {
                if let Some(egui_tiles::Tile::Container(egui_tiles::Container::Linear(row))) =
                    tree.tiles.get_mut(row_id)
                {
                    if at_start {
                        row.children.insert(0, new_unit);
                    } else {
                        row.add_child(new_unit);
                    }
                }
                // Re-equalize this row's unit widths so the added unit gets an
                // equal slice rather than squeezing in beside resized siblings.
                equalize_shares(tree, row_id);
            } else {
                // Row is a lone pair / single pane / tab group / vertical split:
                // wrap it and the new unit into a fresh horizontal row and swap
                // it in for the old row-child inside the vertical root.
                let children = if at_start {
                    vec![new_unit, row_id]
                } else {
                    vec![row_id, new_unit]
                };
                let wrapped = tree.tiles.insert_horizontal_tile(children);
                if let Some(root) = tree.root() {
                    if let Some(egui_tiles::Tile::Container(egui_tiles::Container::Linear(vroot))) =
                        tree.tiles.get_mut(root)
                    {
                        if let Some(slot) = vroot.children.get_mut(i) {
                            *slot = wrapped;
                        }
                    }
                }
            }
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
                // Tests run headless (no wgpu backend): pass no render state.
                app.draw_dock_layout(ctx, None, 1.0);
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
            // Tests run headless (no wgpu backend): pass no render state.
            app.draw_dock_layout(ctx, None, 1.0);
        });
        let panes = pane_ids(app.dock_tree.as_ref().unwrap());
        assert!(panes.contains("valenx_assistant_panel"));
        assert!(panes.contains("workspace:1"));
        assert_eq!(panes.len(), 13);
    }

    /// Test helper: collect all pane ids reachable from `root` (a tile id),
    /// walking through any nested containers. Used to assert which row a moved
    /// unit landed in.
    fn panes_under(
        tree: &egui_tiles::Tree<String>,
        root: egui_tiles::TileId,
    ) -> std::collections::HashSet<String> {
        let mut out = std::collections::HashSet::new();
        let mut stack = vec![root];
        while let Some(id) = stack.pop() {
            match tree.tiles.get(id) {
                Some(egui_tiles::Tile::Pane(p)) => {
                    out.insert(p.clone());
                }
                Some(egui_tiles::Tile::Container(c)) => stack.extend(c.children().copied()),
                None => {}
            }
        }
        out
    }

    #[test]
    fn dock_layout_is_a_noop_when_nothing_open() {
        // With every workbench closed, draw_dock_layout must clear the tree
        // and paint nothing without panicking.
        let mut app = ValenxApp::default();
        app.dock_enabled = true;
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            // Tests run headless (no wgpu backend): pass no render state.
            app.draw_dock_layout(ctx, None, 1.0);
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
                // Tests run headless (no wgpu backend): pass no render state.
                app.draw_dock_layout(ctx, None, 1.0);
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
            // Tests run headless (no wgpu backend): pass no render state.
            app.draw_dock_layout(ctx, None, 1.0);
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

    /// Total pane count in a tree (workspace + agent leaves).
    fn pane_count(tree: &egui_tiles::Tree<String>) -> usize {
        tree.tiles
            .tiles()
            .filter(|tile| matches!(tile, egui_tiles::Tile::Pane(_)))
            .count()
    }

    /// The vertical root's row children (panics if the root isn't a vertical
    /// Linear) — used by the placement tests.
    fn grid_rows(tree: &egui_tiles::Tree<String>) -> Vec<egui_tiles::TileId> {
        let root = tree.root().expect("tree has a root");
        let egui_tiles::Tile::Container(egui_tiles::Container::Linear(vroot)) =
            tree.tiles.get(root).unwrap()
        else {
            panic!("root must be a vertical Linear");
        };
        assert_eq!(vroot.dir, egui_tiles::LinearDir::Vertical);
        vroot.children.clone()
    }

    #[test]
    fn dock_grid_rows_reports_the_3x2_grid_shape() {
        // The 3x2 demo grid is two rows of three unit-pairs each.
        let mut app = ValenxApp::default();
        assert!(app.dock_grid_rows().is_empty(), "no tree → empty row shape");
        app.open_six_workbench_agents();
        assert_eq!(app.dock_grid_rows(), vec![3, 3]);
    }

    #[test]
    fn add_pair_at_row_end_grows_that_row_in_place() {
        // Adding a unit at the RIGHT end of row 0 makes it a 4-unit row and
        // takes the grid from 12 to 14 panes (the new pair = 2 panes), without
        // adding a row.
        let mut app = ValenxApp::default();
        app.open_six_workbench_agents();
        assert_eq!(pane_count(app.dock_tree.as_ref().unwrap()), 12);

        app.add_workbench_agent_pair_at(UnitAddTarget::RowEnd(0));
        let tree = app.dock_tree.as_ref().unwrap();
        assert_eq!(pane_count(tree), 14, "the new pair added two panes");
        assert_eq!(app.dock_grid_rows(), vec![4, 3], "row 0 grew to 4 units");

        // The new unit (n = 7) lives in row 0, at its last position.
        let rows = grid_rows(tree);
        let row0 = panes_under(tree, rows[0]);
        assert!(row0.contains("workspace:7"));
        assert!(row0.contains("agent:7"));
    }

    #[test]
    fn adding_a_unit_equalizes_the_row_shares() {
        // After adding a unit at the end of a row, that row's Linear `shares`
        // are reset to equal (empty map → every child defaults to 1.0), so the
        // new unit auto-sizes the same as its siblings instead of squeezing in.
        let mut app = ValenxApp::default();
        app.open_six_workbench_agents();
        let tree = app.dock_tree.as_mut().unwrap();
        // Manually skew row 0's shares to simulate a prior manual drag.
        let row0 = grid_rows(tree)[0];
        if let Some(egui_tiles::Tile::Container(egui_tiles::Container::Linear(row))) =
            tree.tiles.get_mut(row0)
        {
            let first = row.children[0];
            row.shares.set_share(first, 9.0);
            assert!(
                row.shares.iter().any(|(_, &s)| s == 9.0),
                "precondition: skewed share is present"
            );
        }
        // Now add a unit to the right end of row 0 — it should equalize.
        app.add_workbench_agent_pair_at(UnitAddTarget::RowEnd(0));
        let tree = app.dock_tree.as_ref().unwrap();
        let row0 = grid_rows(tree)[0];
        let egui_tiles::Tile::Container(egui_tiles::Container::Linear(row)) =
            tree.tiles.get(row0).unwrap()
        else {
            panic!("row 0 must be a horizontal Linear");
        };
        // Equalized → the share map was cleared (all children index to the
        // default 1.0), so no explicit non-1.0 entries remain.
        assert!(
            row.shares.iter().all(|(_, &s)| s == 1.0),
            "row shares should be equal after add, got {:?}",
            row.shares.iter().collect::<Vec<_>>()
        );
    }

    #[test]
    fn adding_a_new_bottom_row_equalizes_the_root_row_heights() {
        // NewRowBottom into a vertical grid re-equalizes the vertical root's row
        // shares so the rows stay the same height.
        let mut app = ValenxApp::default();
        app.open_six_workbench_agents();
        // Skew the root's row shares.
        let tree = app.dock_tree.as_mut().unwrap();
        let root = tree.root().unwrap();
        if let Some(egui_tiles::Tile::Container(egui_tiles::Container::Linear(vroot))) =
            tree.tiles.get_mut(root)
        {
            let first = vroot.children[0];
            vroot.shares.set_share(first, 5.0);
        }
        app.add_workbench_agent_pair_at(UnitAddTarget::NewRowBottom);
        let tree = app.dock_tree.as_ref().unwrap();
        let root = tree.root().unwrap();
        let egui_tiles::Tile::Container(egui_tiles::Container::Linear(vroot)) =
            tree.tiles.get(root).unwrap()
        else {
            panic!("root must be a vertical Linear");
        };
        assert!(
            vroot.shares.iter().all(|(_, &s)| s == 1.0),
            "root row shares should be equal after adding a new bottom row"
        );
    }

    #[test]
    fn add_pair_at_row_start_inserts_at_the_left_of_that_row() {
        // RowStart(1) puts the new unit at index 0 of row 1, ahead of the
        // existing units there.
        let mut app = ValenxApp::default();
        app.open_six_workbench_agents();
        app.add_workbench_agent_pair_at(UnitAddTarget::RowStart(1));
        let tree = app.dock_tree.as_ref().unwrap();
        assert_eq!(app.dock_grid_rows(), vec![3, 4], "row 1 grew to 4 units");

        // Row 1 is a horizontal Linear; its FIRST child is the new unit (n=7).
        let rows = grid_rows(tree);
        let egui_tiles::Tile::Container(egui_tiles::Container::Linear(row1)) =
            tree.tiles.get(rows[1]).unwrap()
        else {
            panic!("row 1 must be a horizontal Linear");
        };
        let first_unit = *row1.children.first().expect("row 1 has a first unit");
        let first_panes = panes_under(tree, first_unit);
        assert!(
            first_panes.contains("workspace:7"),
            "the new unit must be at index 0 of row 1, got {first_panes:?}"
        );
        assert!(first_panes.contains("agent:7"));
    }

    #[test]
    fn add_pair_new_row_top_becomes_a_single_unit_first_row() {
        // NewRowTop prepends a brand-new row holding only the new unit; the
        // grid grows from [3,3] to [1,3,3].
        let mut app = ValenxApp::default();
        app.open_six_workbench_agents();
        app.add_workbench_agent_pair_at(UnitAddTarget::NewRowTop);
        assert_eq!(app.dock_grid_rows(), vec![1, 3, 3]);

        let tree = app.dock_tree.as_ref().unwrap();
        let rows = grid_rows(tree);
        let top = panes_under(tree, rows[0]);
        assert_eq!(top.len(), 2, "the new top row holds exactly one unit");
        assert!(top.contains("workspace:7"));
        assert!(top.contains("agent:7"));
    }

    #[test]
    fn add_pair_new_row_bottom_matches_the_legacy_add() {
        // NewRowBottom appends a single-unit last row: [3,3] → [3,3,1].
        let mut app = ValenxApp::default();
        app.open_six_workbench_agents();
        app.add_workbench_agent_pair_at(UnitAddTarget::NewRowBottom);
        assert_eq!(app.dock_grid_rows(), vec![3, 3, 1]);

        // And the View-menu helper still delegates to NewRowBottom.
        let mut app2 = ValenxApp::default();
        app2.open_six_workbench_agents();
        app2.add_workbench_agent_pair();
        assert_eq!(app2.dock_grid_rows(), vec![3, 3, 1]);
    }

    #[test]
    fn add_pair_row_index_out_of_range_falls_back_to_new_bottom_row() {
        // RowEnd(99) can't target a real row → it lands as a new bottom row,
        // and no panes are lost.
        let mut app = ValenxApp::default();
        app.open_six_workbench_agents();
        app.add_workbench_agent_pair_at(UnitAddTarget::RowEnd(99));
        assert_eq!(app.dock_grid_rows(), vec![3, 3, 1]);
        assert_eq!(pane_count(app.dock_tree.as_ref().unwrap()), 14);
    }

    #[test]
    fn add_pair_into_a_lone_row_wraps_it_into_a_horizontal_pair() {
        // Start from a single Workbench+Agent unit (root = one horizontal pair,
        // NOT a vertical grid). A lone pair counts as ONE unit (its children are
        // panes, not nested unit-containers). Then RowEnd(0) on the
        // freshly-formed vertical grid keeps every unit and adds the new one
        // beside it.
        let mut app = ValenxApp::default();
        app.add_workbench_agent_pair(); // unit 1, lone horizontal pair = 1 unit
        assert_eq!(app.dock_grid_rows(), vec![1]);

        // A second unit at the bottom forms the vertical grid [1-unit, 1-unit].
        app.add_workbench_agent_pair_at(UnitAddTarget::NewRowBottom);
        assert_eq!(app.dock_grid_rows(), vec![1, 1]);

        // Now grow row 0 (a lone-unit row) to the right — it must wrap into a
        // horizontal row of two units.
        app.add_workbench_agent_pair_at(UnitAddTarget::RowEnd(0));
        assert_eq!(app.dock_grid_rows(), vec![2, 1]);
        assert_eq!(pane_count(app.dock_tree.as_ref().unwrap()), 6);
    }

    #[test]
    fn pane_camera_orbit_zoom_frame_mutate_the_camera() {
        // The workspace pane (render_tile_mesh_3d caller) applies exactly the
        // central viewport's input math to its *own* camera clone: orbit by
        // `drag_delta * 0.5` (y inverted), zoom by `scroll * 0.01`, and frame to
        // a mesh AABB on double-click. A headless egui Response can't be
        // synthesised here, so this guards the underlying mutations the caller
        // performs — that each genuinely moves the camera (no silent no-op).
        let mut camera = valenx_viz::OrbitCamera::default();
        let before = camera.clone();

        // Orbit: a primary/middle drag of (+10, +6) px → orbit(+5.0, -3.0).
        camera.orbit(10.0 * 0.5, -6.0 * 0.5);
        assert!(
            (camera.azimuth_deg - before.azimuth_deg).abs() > f32::EPSILON,
            "orbit should change azimuth"
        );
        assert!(
            (camera.elevation_deg - before.elevation_deg).abs() > f32::EPSILON,
            "orbit should change elevation"
        );

        // Zoom: a scroll of +40 → zoom(0.4), pulling the camera closer.
        let dist_before = camera.distance;
        camera.zoom(40.0 * 0.01);
        assert!(
            camera.distance < dist_before,
            "scroll-in should reduce distance"
        );

        // Frame: double-click frames the mesh AABB → target moves to its center.
        camera.frame_bounds([0.0, 0.0, 0.0], [4.0, 4.0, 4.0]);
        assert!((camera.target.x - 2.0).abs() < 1e-4);
        assert!((camera.target.y - 2.0).abs() < 1e-4);
        assert!((camera.target.z - 2.0).abs() < 1e-4);
    }
}
