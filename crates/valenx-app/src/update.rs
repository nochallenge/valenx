//! The egui [`eframe::App`] impl plus the `draw_browser` helper that
//! renders the left-hand browser tree. Split out of `lib.rs` as part
//! of the structural refactor.

use eframe::egui;

use valenx_viz::ViewDirection;

use crate::adapter_status::{status_color, status_label};
use crate::commands;
use crate::file_browser::{open_path_or_copy, POPUP_DISABLED_PREFIX};
use crate::first_run;
use crate::histograms::{render_aspect_histogram, render_skewness_histogram};
use crate::keyboard_help;
use crate::landing_page;
use crate::log_panel;
use crate::new_project_dialog;
use crate::panel_help;
use crate::residuals;
use crate::settings;
use crate::settings_io::save_settings_to_state_dir;
use crate::setup::crashes_dir;
use crate::shortcuts::{self, ShortcutAction};
use crate::solver_parse::adapter_id_from_solver;
use crate::time_format::format_time_key;
use crate::types::BottomTab;
use crate::viewport::{self, ShadingMode};
use crate::viewport_kind::ViewportKind;
use crate::welcome_tour;
use crate::ValenxApp;

/// The clickable footprint of a viewport-header icon button — matches the
/// workbench-chrome header (`workbench_chrome::ICON_BUTTON_SIZE`) so the
/// central viewport's `−` / `⋯` / `✕` cluster lines up with the right-side
/// panels' chrome.
const VIEWPORT_ICON_BUTTON_SIZE: f32 = 18.0;

/// Which painter-drawn glyph a [`viewport_icon_button`] renders.
///
/// Mirrors the private `workbench_chrome::Icon` painter geometry so the
/// central viewport header looks identical to the workbench panel chrome.
/// Only `−` (minimize) and `⋯` (more) live here; the `✕` reuses the shared
/// [`crate::workbench_chrome::close_x_button`] directly.
#[derive(Clone, Copy)]
enum ViewportHeaderIcon {
    /// A single horizontal bar — collapse the viewport body to just the header.
    Minimize,
    /// Three dots in a row — open the viewport `⋯` menu.
    More,
}

/// Draw one painter-rendered viewport-header icon button (no font glyph, so it
/// never renders as a missing-glyph "tofu" box), with the same subtle rounded
/// hover background + tooltip as the workbench-panel chrome.
fn viewport_icon_button(ui: &mut egui::Ui, icon: ViewportHeaderIcon, tip: &str) -> egui::Response {
    let (rect, resp) = ui.allocate_exact_size(
        egui::vec2(VIEWPORT_ICON_BUTTON_SIZE, VIEWPORT_ICON_BUTTON_SIZE),
        egui::Sense::click(),
    );
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
        ViewportHeaderIcon::Minimize => {
            painter.line_segment([c + egui::vec2(-r, 0.0), c + egui::vec2(r, 0.0)], stroke);
        }
        ViewportHeaderIcon::More => {
            for dx in [-5.0_f32, 0.0, 5.0] {
                painter.circle_filled(c + egui::vec2(dx, 0.0), 1.5, col);
            }
        }
    }
    resp.on_hover_text(tip)
}

/// Draw the central viewport's slim chrome header — a title on the left and a
/// right-aligned `−  ⋯  ✕` cluster (Close rightmost), matching the right-side
/// workbench panels' chrome. Mutates `app`'s viewport-view-pref fields in
/// place:
///
/// * `✕` sets [`ValenxApp::viewport_hidden`] (hide the viewport; reopen from
///   the `⋯` menu or View → "Hide 3D viewport").
/// * `⋯` opens a small menu: **Reset camera** (resets [`ValenxApp::camera`]
///   to [`valenx_viz::OrbitCamera::default`]), a **3D / 2D viewport switch**
///   (mirrors the View-menu radios), and **Hide viewport**.
/// * `−` toggles [`ValenxApp::viewport_collapsed`] (roll the body up to just
///   this header).
fn viewport_chrome_header(app: &mut ValenxApp, ui: &mut egui::Ui) {
    let title = match app.active_viewport {
        ViewportKind::Viewport2dDna => "2D Sketch",
        ViewportKind::Viewport3D => "3D Viewport",
    };
    ui.horizontal(|ui| {
        // Keep the controls a fixed, equal distance apart (matches the
        // workbench header spacing).
        ui.spacing_mut().item_spacing.x = 4.0;
        ui.with_layout(egui::Layout::left_to_right(egui::Align::Center), |ui| {
            ui.strong(title);
        });
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            // Right-to-left: add Close first so the visual order is `−  ⋯  ✕`.
            if crate::workbench_chrome::close_x_button(
                ui,
                "Hide viewport — reopen from the View menu",
            )
            .clicked()
            {
                app.viewport_hidden = true;
            }
            // `⋯` opens a small popup keyed on the button's own id, so the
            // trigger stays a painter-drawn icon (no font-glyph caret).
            let more = viewport_icon_button(
                ui,
                ViewportHeaderIcon::More,
                "Reset camera, switch 2D/3D, hide",
            );
            let popup_id = more.id.with("viewport_chrome_more");
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
                    if ui.button("Reset camera").clicked() {
                        app.camera = valenx_viz::OrbitCamera::default();
                    }
                    ui.separator();
                    // 3D / 2D switch — mirrors the View-menu "Central
                    // viewport" radios.
                    if ui
                        .radio(
                            app.active_viewport == ViewportKind::Viewport3D,
                            "3D Viewport",
                        )
                        .clicked()
                    {
                        app.active_viewport = ViewportKind::Viewport3D;
                    }
                    if ui
                        .radio(
                            app.active_viewport == ViewportKind::Viewport2dDna,
                            "2D Sketch",
                        )
                        .clicked()
                    {
                        app.active_viewport = ViewportKind::for_genetics();
                    }
                    ui.separator();
                    if ui.button("Hide viewport").clicked() {
                        app.viewport_hidden = true;
                    }
                },
            );
            let tip = if app.viewport_collapsed {
                "Expand viewport"
            } else {
                "Collapse viewport"
            };
            if viewport_icon_button(ui, ViewportHeaderIcon::Minimize, tip).clicked() {
                app.viewport_collapsed = !app.viewport_collapsed;
            }
        });
    });
    ui.separator();
}

impl eframe::App for ValenxApp {
    fn update(&mut self, ctx: &egui::Context, frame: &mut eframe::Frame) {
        if !self.theme_applied {
            self.settings.apply_theme(ctx);
            self.theme_applied = true;
        }

        // Agent-bridge heartbeat (UNCONDITIONAL, every frame). egui is
        // reactive, so when valenx is idle/unfocused — the normal case while
        // an external agent drives it from the background — `update()` would
        // otherwise stop being called, and the agent-command poll below
        // (`poll_and_apply_agent_commands`) would never run. Commands appended
        // to a channel's cmd file (including the **global** `new_unit`, which
        // must be processed even when zero Workbench+Agent units exist) would
        // then sit unprocessed. Scheduling a repaint on the same cadence as
        // `agent_commands::POLL_INTERVAL` wakes egui ~once/second with zero
        // input so the poll always fires. This is additive: `request_repaint_after`
        // is a no-op when a sooner repaint is already scheduled (run/sweep,
        // automation, adapter-probe, etc.), so the idle cost is only ~1 fps.
        ctx.request_repaint_after(crate::agent_commands::POLL_INTERVAL);

        self.pump_run_events();
        self.pump_sweep_events();
        // Repaint promptly while a run or sweep is live so residual
        // updates and per-case progress don't lag the channel.
        if self.run_handle.is_some() || self.sweep_handle.is_some() {
            ctx.request_repaint_after(std::time::Duration::from_millis(50));
        }

        // Automation / front-end testing: when the VALENX_AUTOMATION env var
        // is set, keep frames ticking (~10 fps) even when the app is
        // otherwise idle. egui only pushes the accessibility tree and
        // processes queued UI-Automation actions (AccessKit Invoke / Toggle
        // from an external test driver) on a live frame, so without this an
        // idle window leaves a stale tree and dropped actions. Read once and
        // cached; zero effect on the normal interactive build.
        {
            use std::sync::OnceLock;
            static AUTOMATION: OnceLock<bool> = OnceLock::new();
            if *AUTOMATION.get_or_init(|| std::env::var_os("VALENX_AUTOMATION").is_some()) {
                ctx.request_repaint_after(std::time::Duration::from_millis(100));
            }
        }

        // Drag-and-drop
        let dropped = ctx.input(|i| i.raw.dropped_files.clone());
        for f in dropped {
            if let Some(path) = f.path {
                let lower = path.to_string_lossy().to_ascii_lowercase();
                if lower.ends_with(".stl") {
                    self.load_stl(path);
                } else if lower.ends_with(".msh") || lower.ends_with("mesh.canonical.json") {
                    self.load_mesh(path);
                } else if path.is_dir() && lower.ends_with(".valenx") {
                    self.load_project(path);
                } else {
                    self.last_error = Some(format!(
                        "Don't know how to open {:?}",
                        path.file_name().unwrap_or_default()
                    ));
                }
            }
        }

        // Keyboard shortcuts. The polled-action layer in
        // `crate::shortcuts` decodes input into typed
        // `ShortcutAction`s the dispatcher below routes to the
        // right `ValenxApp` method. The view-direction numpad +
        // F (frame) + F5 (run) shortcuts predate the new layer
        // and stay inline — they're viewport-specific and bypass
        // the `Ctrl+digit` workbench-switching bindings since
        // the digit alone (no modifier) means "view direction".
        let mut run_selected_requested = false;
        let palette_open = self.palette.open;
        let shortcut_actions: Vec<ShortcutAction> =
            ctx.input(|i| shortcuts::poll_shortcut(i, palette_open));
        ctx.input(|i| {
            if i.key_pressed(egui::Key::F) && !palette_open {
                // Prefer the mesh's AABB if it's loaded — users who
                // have both usually care about the meshed domain.
                if self.mesh.is_some() {
                    self.frame_current_mesh();
                } else {
                    self.frame_current_stl();
                }
            }
            if i.key_pressed(egui::Key::F5) && !palette_open {
                run_selected_requested = true;
            }
            // View-direction numpad (no modifier). The Ctrl+digit
            // bindings for workbench switching live in shortcuts.rs
            // and don't collide because they require Ctrl.
            if !palette_open && !i.modifiers.ctrl && !i.modifiers.command {
                for (key, dir) in [
                    (egui::Key::Num1, ViewDirection::Front),
                    (egui::Key::Num2, ViewDirection::Back),
                    (egui::Key::Num3, ViewDirection::Right),
                    (egui::Key::Num4, ViewDirection::Left),
                    (egui::Key::Num5, ViewDirection::Top),
                    (egui::Key::Num6, ViewDirection::Bottom),
                    (egui::Key::Num7, ViewDirection::Iso),
                ] {
                    if i.key_pressed(key) {
                        self.camera.set_view(dir);
                    }
                }
            }
        });
        for action in shortcut_actions {
            self.dispatch_shortcut(action);
        }
        if run_selected_requested {
            self.run_selected_case();
        }

        // Apply background adapter-probe results as they stream in (see
        // AdapterRegistry::spawn_probe_all). This keeps startup instant
        // instead of blocking the first frame while 141 external tools are
        // probed. Native adapter paths work regardless of probe state.
        if let Some(rx) = self.adapter_probe_rx.take() {
            let mut still_probing = true;
            loop {
                match rx.try_recv() {
                    Ok((id, status)) => {
                        self.registry.apply_probe_result(id, status);
                    }
                    Err(std::sync::mpsc::TryRecvError::Empty) => break,
                    Err(std::sync::mpsc::TryRecvError::Disconnected) => {
                        still_probing = false;
                        break;
                    }
                }
            }
            if still_probing {
                // Put the receiver back and keep frames coming so results
                // keep applying even if the UI is otherwise idle.
                self.adapter_probe_rx = Some(rx);
                ctx.request_repaint();
            }
        }

        // Persist window geometry on change so Valenx reopens where you
        // left it (e.g. a second monitor). Rounded to whole pixels so a
        // stationary window doesn't re-save on sub-pixel noise; writes
        // settings.json only on an actual move/resize. Gated on the setting.
        if self.settings.remember_window_geometry {
            let (pos, size) = ctx.input(|i| {
                let info = i.viewport();
                (
                    info.outer_rect.map(|r| [r.min.x.round(), r.min.y.round()]),
                    info.inner_rect
                        .map(|r| [r.width().round(), r.height().round()]),
                )
            });
            let mut dirty = false;
            if let Some(p) = pos {
                // Skip off-screen sentinels (minimized windows can report
                // huge negative coordinates on some platforms).
                if p[0] > -30000.0 && p[1] > -30000.0 && self.settings.window_position != Some(p) {
                    self.settings.window_position = Some(p);
                    dirty = true;
                }
            }
            if let Some(s) = size {
                if s[0] >= 1.0 && s[1] >= 1.0 && self.settings.window_size != Some(s) {
                    self.settings.window_size = Some(s);
                    dirty = true;
                }
            }
            if dirty {
                crate::settings_io::save_settings_to_state_dir(&self.settings);
            }
        }

        // Ribbon
        egui::TopBottomPanel::top("valenx_ribbon").show(ctx, |ui| {
            egui::menu::bar(ui, |ui| {
                // Clone once per frame so the immutable borrow doesn't
                // collide with `self.pick_project()` / `self.reprobe()`
                // / etc. inside each menu_button closure (those need
                // `&mut self`). Catalogue is a small BTreeMap so the
                // clone is microseconds.
                let cat = self.catalogue.clone();
                ui.menu_button(cat.lookup("menu.file"), |ui| {
                    if ui.button("New Project…\tCtrl+N").clicked() {
                        ui.close_menu();
                        self.open_new_project_dialog();
                    }
                    if ui.button("Open Project…\tCtrl+O").clicked() {
                        ui.close_menu();
                        self.pick_project();
                    }
                    if ui
                        .add_enabled(
                            self.project.is_some(),
                            egui::Button::new("Close Project (Home)"),
                        )
                        .on_hover_text(
                            "Close the current project and return to the home page \
                             to open or create another.",
                        )
                        .clicked()
                    {
                        ui.close_menu();
                        self.close_project();
                    }
                    ui.separator();
                    if ui.button(cat.lookup("menu.file.import-stl")).clicked() {
                        ui.close_menu();
                        self.pick_stl();
                    }
                    if ui.button(cat.lookup("menu.file.load-mesh")).clicked() {
                        ui.close_menu();
                        self.pick_mesh();
                    }
                    ui.separator();
                    if ui.button(cat.lookup("menu.file.exit")).clicked() {
                        ctx.send_viewport_cmd(egui::ViewportCommand::Close);
                    }
                });
                ui.menu_button(cat.lookup("menu.view"), |ui| {
                    crate::menu_ui::scrollable_menu(ui, |ui| {
                    ui.horizontal(|ui| {
                        ui.label("Text size:");
                        if ui.button("A-").on_hover_text("Smaller text (Ctrl -)").clicked() {
                            let z = (ui.ctx().zoom_factor() - 0.1_f32).max(0.6);
                            ui.ctx().set_zoom_factor(z);
                        }
                        if ui.button("Reset").on_hover_text("Default text size").clicked() {
                            ui.ctx().set_zoom_factor(1.3);
                        }
                        if ui.button("A+").on_hover_text("Bigger text (Ctrl +)").clicked() {
                            let z = (ui.ctx().zoom_factor() + 0.1_f32).min(2.5);
                            ui.ctx().set_zoom_factor(z);
                        }
                    });
                    ui.separator();
                    if ui
                        .radio(
                            self.shading == ShadingMode::Shaded,
                            cat.lookup("menu.view.shaded"),
                        )
                        .clicked()
                    {
                        self.shading = ShadingMode::Shaded;
                        ui.close_menu();
                    }
                    if ui
                        .radio(
                            self.shading == ShadingMode::Wireframe,
                            cat.lookup("menu.view.wireframe"),
                        )
                        .clicked()
                    {
                        self.shading = ShadingMode::Wireframe;
                        ui.close_menu();
                    }
                    ui.separator();
                    if ui.button(cat.lookup("menu.view.frame")).clicked() {
                        self.frame_current_stl();
                        ui.close_menu();
                    }
                    ui.separator();
                    // Toggle the left-side Browser panel — collapse it to
                    // give the viewport the full width.
                    if ui
                        .checkbox(&mut self.show_browser, "Browser")
                        .on_hover_text(
                            "Show / hide the left Browser panel (project tree, \
                             cases, geometry, mesh, results).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.docked_layout, "Docked layout")
                        .on_hover_text(
                            "Opt-in: arrange panels as resizable, tabbable, \
                             closable tiles (rows + columns) in the central area. \
                             Off = classic single-viewport layout.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.dock_enabled, "Dockable panel layout (beta)")
                        .on_hover_text(
                            "Opt-in: host the open right-side workbenches in one \
                             dockable region — grab a panel by its tab and drag to \
                             reorder, drop into a second row, or split. \
                             Off = classic stacked side panels.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // "Workbench + Agent" launchers — each unit is a paired
                    // empty workspace (build canvas) + Claude chat tile in the
                    // dock. Both turn the dockable layout on as a side effect.
                    if ui
                        .button("New Workbench + Agent")
                        .on_hover_text(
                            "Add a paired Workspace + Agent unit to the dockable \
                             region: an empty build canvas on the left and a Claude \
                             chat on the right. Prompt the agent to build into the \
                             workspace. (Turns the dockable layout on.)",
                        )
                        .clicked()
                    {
                        self.add_workbench_agent_pair();
                        ui.close_menu();
                    }
                    if ui
                        .button("Open 6 Workbench+Agents (3x2 grid)")
                        .on_hover_text(
                            "Open six paired Workspace + Agent units arranged as a \
                             3-column x 2-row grid in the dockable region. \
                             (Turns the dockable layout on.)",
                        )
                        .clicked()
                    {
                        self.open_six_workbench_agents();
                        ui.close_menu();
                    }
                    ui.checkbox(&mut self.snap_to_grid, "Snap to grid")
                        .on_hover_text(
                            "Snap the viewport cursor to the nearest ground-grid \
                             node — shows a marker + snapped X/Y/Z readout.",
                        );
                    // Hide / show the central 2D/3D viewport body. Mirrors the
                    // viewport header's ✕ / "Hide viewport" so the user can
                    // always bring it back from the menu.
                    if ui
                        .checkbox(&mut self.viewport_hidden, "Hide 3D viewport")
                        .on_hover_text(
                            "Hide the central 2D/3D viewport body (its render is \
                             skipped); the central area shows a placeholder. \
                             Uncheck — or use the viewport header / its ⋯ menu — \
                             to bring it back.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }

                    ui.separator();
                    ui.label(egui::RichText::new("Workbenches").weak().small())
                        .on_hover_text(
                            "Right-side tool panels: Mesh Toolbox (CAD), Genetics \
                             (15 native computational-biology panels), Aerodynamics \
                             (wind tunnel), and Astro (launch). Toggle whichever you're \
                             working in.",
                        );
                    // Toggle the right-side Mesh Toolbox. Hidden by
                    // default once flipped off; reopen from here or
                    // from the command palette.
                    if ui
                        .checkbox(&mut self.show_mesh_toolbox, "Mesh Toolbox")
                        .on_hover_text(
                            "Show / hide the right-side Mesh Toolbox panel \
                             (Inspector + Transformations + Cut plane + Repair).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Genetics Workbench — the
                    // Round-6 computational-biology toolkits (sequence,
                    // alignment, phylogenetics, RNA folding, MD,
                    // cheminformatics, quantum chemistry, genomics,
                    // docking, gene editing, …). Off by default.
                    if ui
                        .checkbox(&mut self.show_genetics_workbench, "Genetics Workbench")
                        .on_hover_text(
                            "Show / hide the right-side Genetics Workbench — 15 native \
                             computational-biology panels (sequence, alignment, \
                             phylogenetics, population genetics, RNA structure, RNA \
                             design, molecular dynamics, cheminformatics, macromolecular \
                             structure, quantum chemistry, genomics, systems biology, \
                             docking, gene editing, structure prediction).",
                        )
                        .changed()
                    {
                        // When the Genetics Workbench is first opened,
                        // auto-switch the central viewport to the 2D DNA /
                        // plasmid view — the user almost certainly wants to
                        // see their sequence rather than an empty 3D scene.
                        // Only fires when switching *on* and the user hasn't
                        // already manually chosen a non-3D viewport.
                        if self.show_genetics_workbench
                            && self.active_viewport == ViewportKind::Viewport3D
                        {
                            self.active_viewport = ViewportKind::for_genetics();
                        }
                        ui.close_menu();
                    }
                    // Toggle the right-side Aerodynamics / Wind Tunnel
                    // workbench — the native 3-D external-aero CFD
                    // engine (virtual wind tunnel: drag/lift/moment
                    // coefficients, Cp / flow fields, AoA polar sweep).
                    // Off by default.
                    if ui
                        .checkbox(
                            &mut self.show_aero_workbench,
                            "Aerodynamics / Wind Tunnel",
                        )
                        .on_hover_text(
                            "Show / hide the right-side Wind Tunnel workbench — a \
                             virtual wind tunnel (3-D external-aerodynamics CFD): \
                             drag / lift / moment coefficients, Cp and flow-field \
                             visualization, and an angle-of-attack polar sweep.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side FEM Workbench — native
                    // finite-element analysis (linear-static bending +
                    // modal natural frequencies) on valenx-fem's
                    // in-process solvers. Off by default.
                    if ui
                        .checkbox(&mut self.show_fem_workbench, "FEM Workbench")
                        .on_hover_text(
                            "Show / hide the right-side FEM Workbench — native \
                             finite-element analysis (linear-static bending + modal \
                             natural frequencies) on a built-in structured mesh, \
                             solved in-process by valenx-fem (no external solver).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Induction Motor workbench — native
                    // 3-phase induction-motor slip / power (valenx-inductionmotor).
                    if ui
                        .checkbox(&mut self.show_inductionmotor_workbench, "Induction Motor")
                        .on_hover_text(
                            "Show / hide the right-side Induction Motor Workbench — native \
                             3-phase induction machine: synchronous speed, slip, rotor \
                             frequency and the air-gap power split into rotor copper loss and \
                             mechanical power, computed in-process by valenx-inductionmotor.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side CFD Workbench — native 2-D
                    // incompressible laminar CFD (SIMPLE): the lid-driven
                    // cavity + developing channel-flow cases on
                    // valenx-cfd-native. Off by default.
                    if ui
                        .checkbox(&mut self.show_cfd_workbench, "CFD Workbench")
                        .on_hover_text(
                            "Show / hide the right-side CFD Workbench — native 2-D \
                             incompressible laminar CFD (the SIMPLE algorithm): the \
                             lid-driven cavity and developing channel-flow benchmarks, \
                             solved in-process by valenx-cfd-native (no external solver).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Reaction Dynamics workbench —
                    // native ab-initio MD (AIMD): small molecules move and
                    // react in 3-D under first-principles quantum forces
                    // (valenx-reactdyn). Off by default.
                    if ui
                        .checkbox(&mut self.show_reactdyn_workbench, "Reaction Dynamics")
                        .on_hover_text(
                            "Show / hide the right-side Reaction Dynamics workbench — native \
                             ab-initio molecular dynamics (Born-Oppenheimer AIMD): watch small \
                             molecules move and react under first-principles quantum forces, \
                             solved in-process by valenx-reactdyn.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Springs Workbench — native
                    // helical-spring design (valenx-springs). Off by default.
                    if ui
                        .checkbox(&mut self.show_springs_workbench, "Springs")
                        .on_hover_text(
                            "Show / hide the right-side Springs Workbench — native helical-spring \
                             design: spring index, axial stiffness, Wahl factor, the diameters, \
                             pitch, and developed wire length, computed in-process by valenx-springs.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Bearing workbench — native rolling
                    // bearing ISO 281 rating life (valenx-bearing).
                    if ui
                        .checkbox(&mut self.show_bearing_workbench, "Bearing")
                        .on_hover_text(
                            "Show / hide the right-side Bearing Workbench — native rolling-bearing \
                             ISO 281 basic rating life: dynamic equivalent load → L10 in \
                             revolutions and hours, computed in-process by valenx-bearing.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Belt Drive workbench — native flat-belt
                    // drive analysis (valenx-beltdrive).
                    if ui
                        .checkbox(&mut self.show_beltdrive_workbench, "Belt Drive")
                        .on_hover_text(
                            "Show / hide the right-side Belt Drive Workbench — native open \
                             flat-belt drive: speed ratio, belt speed, wrap angles, capstan \
                             tension ratio and transmissible power, computed in-process by \
                             valenx-beltdrive.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Buckling workbench — native Euler /
                    // Johnson column buckling (valenx-buckling).
                    if ui
                        .checkbox(&mut self.show_buckling_workbench, "Buckling")
                        .on_hover_text(
                            "Show / hide the right-side Buckling Workbench — native Euler / Johnson \
                             column buckling: critical load, critical stress, slenderness and the \
                             design load for the end condition, computed in-process by \
                             valenx-buckling.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Brake workbench — native disc-brake
                    // torque and stop dynamics (valenx-brake).
                    if ui
                        .checkbox(&mut self.show_brake_workbench, "Brake")
                        .on_hover_text(
                            "Show / hide the right-side Brake Workbench — native disc-brake \
                             friction torque from clamp force and effective radius, plus the stop \
                             energy / distance / time, computed in-process by valenx-brake.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Fatigue workbench — native high-cycle
                    // stress-life analysis (valenx-fatigue).
                    if ui
                        .checkbox(&mut self.show_fatigue_workbench, "Fatigue")
                        .on_hover_text(
                            "Show / hide the right-side Fatigue Workbench — native high-cycle \
                             stress-life: Goodman / Soderberg / Gerber mean-stress factor of \
                             safety and S-N life, computed in-process by valenx-fatigue.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Gear Tooth workbench — native
                    // spur-gear tooth bending strength (valenx-geartooth).
                    if ui
                        .checkbox(&mut self.show_geartooth_workbench, "Gear Tooth")
                        .on_hover_text(
                            "Show / hide the right-side Gear Tooth Workbench — native spur-gear \
                             tooth bending strength (Lewis equation + AGMA refinement), computed \
                             in-process by valenx-geartooth.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Pharmacokinetics workbench — native
                    // one-compartment IV dosing (valenx-pharmacokinetics).
                    if ui
                        .checkbox(&mut self.show_pharmacokinetics_workbench, "Pharmacokinetics")
                        .on_hover_text(
                            "Show / hide the right-side Pharmacokinetics Workbench — native \
                             one-compartment IV-bolus dosing: elimination rate, half-life, Cmax, \
                             AUC and concentration-time, computed in-process by \
                             valenx-pharmacokinetics.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Pipe Network workbench — native
                    // Hardy-Cross flow balancing (valenx-pipenetwork).
                    if ui
                        .checkbox(&mut self.show_pipenetwork_workbench, "Pipe Network")
                        .on_hover_text(
                            "Show / hide the right-side Pipe Network Workbench — native \
                             pipe-network flow balancing (Hardy-Cross loop solve) for the loop \
                             flows and head losses, computed in-process by valenx-pipenetwork.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side RC Beam workbench — native reinforced
                    // concrete beam flexure (valenx-rcbeam).
                    if ui
                        .checkbox(&mut self.show_rcbeam_workbench, "RC Beam")
                        .on_hover_text(
                            "Show / hide the right-side RC Beam Workbench — native singly- \
                             reinforced concrete beam flexure (Whitney stress block): stress-block \
                             depth, lever arm, nominal and design moment, computed in-process by \
                             valenx-rcbeam.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Marine / Hull Workbench — native
                    // box-form hull hydrostatics (valenx-marine). Off by default.
                    if ui
                        .checkbox(&mut self.show_marine_workbench, "Marine / Hull")
                        .on_hover_text(
                            "Show / hide the right-side Marine / Hull Workbench — box-form hull \
                             hydrostatics: displacement, the centre of buoyancy KB, the metacentric \
                             radius BM and the metacentric height GM with a stability verdict, plus a \
                             3-D hull solid, computed in-process by valenx-marine.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Capacitor workbench — native
                    // parallel-plate electrostatics (valenx-capacitor).
                    if ui
                        .checkbox(&mut self.show_capacitor_workbench, "Capacitor")
                        .on_hover_text(
                            "Show / hide the right-side Capacitor Workbench — native \
                             parallel-plate capacitance (εr·ε0·A/d), stored energy, charge and \
                             capacitive reactance, computed in-process by valenx-capacitor.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Fan Laws workbench — native fan /
                    // blower affinity-law scaling (valenx-fanlaws).
                    if ui
                        .checkbox(&mut self.show_fanlaws_workbench, "Fan Laws")
                        .on_hover_text(
                            "Show / hide the right-side Fan Laws Workbench — native fan / blower \
                             affinity-law scaling (flow ~ N, pressure ~ N^2, power ~ N^3) to a new \
                             speed, computed in-process by valenx-fanlaws.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Creep workbench — native
                    // high-temperature creep / rupture (valenx-creep).
                    if ui
                        .checkbox(&mut self.show_creep_workbench, "Creep")
                        .on_hover_text(
                            "Show / hide the right-side Creep Workbench — native high-temperature \
                             creep: Larson-Miller rupture life and the Norton-Bailey secondary \
                             creep rate, computed in-process by valenx-creep.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Electrochemistry workbench — native
                    // Nernst / Faraday cell analysis (valenx-electrochem).
                    if ui
                        .checkbox(&mut self.show_electrochem_workbench, "Electrochemistry")
                        .on_hover_text(
                            "Show / hide the right-side Electrochemistry Workbench — native Nernst \
                             cell potential and Faraday electrolysis (mass deposited), computed \
                             in-process by valenx-electrochem.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Enzyme Kinetics workbench — native
                    // Michaelis-Menten rate laws (valenx-enzymekinetics).
                    if ui
                        .checkbox(&mut self.show_enzymekinetics_workbench, "Enzyme Kinetics")
                        .on_hover_text(
                            "Show / hide the right-side Enzyme Kinetics Workbench — native \
                             Michaelis-Menten and reversible-inhibition initial-rate laws, \
                             computed in-process by valenx-enzymekinetics.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Gears Workbench — native
                    // involute-gear design (valenx-gears). Off by default.
                    if ui
                        .checkbox(&mut self.show_gears_workbench, "Gears")
                        .on_hover_text(
                            "Show / hide the right-side Gears Workbench — native involute-gear \
                             design: circular pitch, the pitch / base / addendum / dedendum \
                             diameters, and the meshing gear ratio, computed in-process by valenx-gears.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Pneumatics workbench (valenx-pneumatics).
                    if ui
                        .checkbox(&mut self.show_pneumatics_workbench, "Pneumatics")
                        .on_hover_text(
                            "Show / hide the Pneumatics Workbench — native pneumatic cylinder \
                             force, compression ratio, choked-flow check and free-air consumption, \
                             computed in-process by valenx-pneumatics.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Psychrometrics workbench (valenx-psychrometrics).
                    if ui
                        .checkbox(&mut self.show_psychrometrics_workbench, "Psychrometrics")
                        .on_hover_text(
                            "Show / hide the Psychrometrics Workbench — native moist-air state: \
                             humidity ratio, dew point, degree of saturation and enthalpy, \
                             computed in-process by valenx-psychrometrics.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Thermistor workbench (valenx-thermistor).
                    if ui
                        .checkbox(&mut self.show_thermistor_workbench, "Thermistor")
                        .on_hover_text(
                            "Show / hide the Thermistor Workbench — native NTC resistance \
                             ↔ temperature (Beta and Steinhart-Hart models) and temperature \
                             coefficient, computed in-process by valenx-thermistor.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Strain Gauge workbench (valenx-straingauge).
                    if ui
                        .checkbox(&mut self.show_straingauge_workbench, "Strain Gauge")
                        .on_hover_text(
                            "Show / hide the Strain Gauge Workbench — native Wheatstone-bridge \
                             gauge factor → bridge output and stress for quarter / half / full \
                             bridges, computed in-process by valenx-straingauge.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Drone Workbench — native
                    // multirotor hover performance (valenx-drone). Off by default.
                    if ui
                        .checkbox(&mut self.show_drone_workbench, "Drone / Multirotor")
                        .on_hover_text(
                            "Show / hide the right-side Drone Workbench — native multirotor \
                             hover performance: disk loading, induced velocity, ideal / actual \
                             hover power, thrust-to-weight, and battery endurance, computed \
                             in-process by valenx-drone.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Acoustics workbench — native room
                    // reverberation + SPL (valenx-acoustics).
                    if ui
                        .checkbox(&mut self.show_acoustics_workbench, "Acoustics")
                        .on_hover_text(
                            "Show / hide the right-side Acoustics Workbench — native \
                             rectangular-room reverberation (Sabine / Eyring RT60), lowest room \
                             mode and free-field SPL distance drop, computed in-process by \
                             valenx-acoustics.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Acid-Base workbench — native aqueous
                    // pH / buffer equilibria (valenx-acidbase).
                    if ui
                        .checkbox(&mut self.show_acidbase_workbench, "Acid-Base")
                        .on_hover_text(
                            "Show / hide the right-side Acid-Base Workbench — native aqueous pH \
                             for strong / weak acids and bases plus Henderson-Hasselbalch \
                             buffers, computed in-process by valenx-acidbase.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side BJT workbench — native bipolar
                    // transistor DC bias (valenx-bjt).
                    if ui
                        .checkbox(&mut self.show_bjt_workbench, "BJT")
                        .on_hover_text(
                            "Show / hide the right-side BJT Workbench — native bipolar-junction \
                             transistor DC bias Q-point (currents, Vce, region, stability factor) \
                             for divider / fixed-base networks, computed in-process by valenx-bjt.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side BMR / TDEE workbench — native resting
                    // + daily energy expenditure (valenx-bmr).
                    if ui
                        .checkbox(&mut self.show_bmr_workbench, "BMR / TDEE")
                        .on_hover_text(
                            "Show / hide the right-side BMR / TDEE Workbench — native basal \
                             metabolic rate (Mifflin-St Jeor / Harris-Benedict) and total daily \
                             energy expenditure, computed in-process by valenx-bmr.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Bolted Joint workbench — native
                    // preloaded joint mechanics (valenx-bolt).
                    if ui
                        .checkbox(&mut self.show_bolt_workbench, "Bolted Joint")
                        .on_hover_text(
                            "Show / hide the right-side Bolted Joint Workbench — native preloaded \
                             bolted-joint mechanics: preload from torque, bolt / member load \
                             sharing, separation and overload safety, computed in-process by \
                             valenx-bolt.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Geomatics Workbench — native
                    // geodesic calculations (valenx-geomatics). Off by default.
                    if ui
                        .checkbox(&mut self.show_geomatics_workbench, "Geomatics")
                        .on_hover_text(
                            "Show / hide the right-side Geomatics Workbench — native geodesic \
                             calculations: great-circle (haversine) + rhumb-line distance, initial \
                             and final bearings, and cross-track / along-track offsets, computed \
                             in-process by valenx-geomatics.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // EE / DSP workbenches (electronics batch). Each off by default.
                    if ui
                        .checkbox(&mut self.show_opamp_workbench, "Op-Amp")
                        .on_hover_text(
                            "Show / hide the right-side Op-Amp Workbench — native ideal \
                             closed-loop gain and bandwidth (valenx-opamp).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.show_led_workbench, "LED")
                        .on_hover_text(
                            "Show / hide the right-side LED Workbench — native series \
                             current-limiting resistor sizing (valenx-led).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.show_thermocouple_workbench, "Thermocouple")
                        .on_hover_text(
                            "Show / hide the right-side Thermocouple Workbench — native \
                             Seebeck EMF + cold-junction compensation (valenx-thermocouple).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(
                            &mut self.show_transmissionline_workbench,
                            "Transmission Line",
                        )
                        .on_hover_text(
                            "Show / hide the right-side Transmission Line Workbench — native \
                             lossless reflection / VSWR (valenx-transmissionline).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.show_powerfactor_workbench, "Power Factor")
                        .on_hover_text(
                            "Show / hide the right-side Power Factor Workbench — native AC \
                             power triangle + shunt-capacitor correction (valenx-powerfactor).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(
                            &mut self.show_resistornetwork_workbench,
                            "Resistor Network",
                        )
                        .on_hover_text(
                            "Show / hide the right-side Resistor Network Workbench — native \
                             series / parallel / divider analysis (valenx-resistor-network).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.show_rectifier_workbench, "Rectifier")
                        .on_hover_text(
                            "Show / hide the right-side Rectifier Workbench — native ideal \
                             rectifier figures + capacitor-filter ripple (valenx-rectifier).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.show_filter_workbench, "RC / RLC Filter")
                        .on_hover_text(
                            "Show / hide the right-side Filter Workbench — native RC / series-RLC \
                             analog filter response (valenx-filter).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Heat Transfer workbench — native
                    // composite-wall 1-D heat loss (valenx-heat-transfer).
                    if ui
                        .checkbox(&mut self.show_heattransfer_workbench, "Heat Transfer")
                        .on_hover_text(
                            "Show / hide the right-side Heat Transfer Workbench — native \
                             composite-wall heat loss: series conduction + convective-film \
                             resistances, total resistance, heat-loss rate and the overall \
                             U-value, computed in-process by valenx-heat-transfer.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Four-Bar Linkage Workbench — native
                    // planar mechanism kinematics (valenx-kinematics). Off by default.
                    if ui
                        .checkbox(&mut self.show_fourbar_workbench, "Four-Bar Linkage")
                        .on_hover_text(
                            "Show / hide the right-side Four-Bar Linkage Workbench — native planar \
                             four-bar kinematics: Grashof classification and, at a crank angle, the \
                             solved pin positions, coupler / rocker angles and transmission angle, \
                             computed in-process by valenx-kinematics.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Mechanical + civil batch workbenches. Each off by default.
                    if ui
                        .checkbox(&mut self.show_shaftdesign_workbench, "Shaft Design")
                        .on_hover_text(
                            "Show / hide the right-side Shaft Design Workbench — native \
                             combined bending + torsion shaft sizing (valenx-shaftdesign).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.show_screwthread_workbench, "Power Screw")
                        .on_hover_text(
                            "Show / hide the right-side Power Screw Workbench — native \
                             square-thread lead-screw torque + efficiency (valenx-screwthread).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.show_pulley_workbench, "Pulley System")
                        .on_hover_text(
                            "Show / hide the right-side Pulley System Workbench — native \
                             block-and-tackle mechanical advantage (valenx-pulley).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.show_springdesign_workbench, "Spring Design")
                        .on_hover_text(
                            "Show / hide the right-side Spring Design Workbench — native \
                             helical compression spring (Wahl, rate, stress) (valenx-spring-design).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(
                            &mut self.show_springcombination_workbench,
                            "Spring Combination",
                        )
                        .on_hover_text(
                            "Show / hide the right-side Spring Combination Workbench — native \
                             series / parallel spring networks (valenx-springcombination).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.show_vibration_workbench, "Vibration (SDOF)")
                        .on_hover_text(
                            "Show / hide the right-side Vibration Workbench — native single-DOF \
                             forced-vibration response (valenx-vibration).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.show_rivet_workbench, "Riveted Joint")
                        .on_hover_text(
                            "Show / hide the right-side Riveted Joint Workbench — native \
                             rivet-joint strength + failure mode (valenx-rivet).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.show_soilbearing_workbench, "Soil Bearing")
                        .on_hover_text(
                            "Show / hide the right-side Soil Bearing Workbench — native \
                             Terzaghi strip-footing bearing capacity (valenx-soilbearing).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Piping Workbench — native
                    // pipe-section sizing (valenx-piping). Off by default.
                    if ui
                        .checkbox(&mut self.show_piping_workbench, "Piping")
                        .on_hover_text(
                            "Show / hide the right-side Piping Workbench — native pipe-section \
                             sizing: NPS outer / inner diameters, flow + metal cross-section areas, \
                             wetted perimeter, and external surface area, computed in-process by \
                             valenx-piping.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Science / bio / civil batch workbenches. Each off by default.
                    if ui
                        .checkbox(&mut self.show_retainingwall_workbench, "Retaining Wall")
                        .on_hover_text(
                            "Show / hide the right-side Retaining Wall Workbench — native \
                             Rankine lateral earth pressure + thrust (valenx-retainingwall).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.show_openchannel_workbench, "Open Channel")
                        .on_hover_text(
                            "Show / hide the right-side Open Channel Workbench — native \
                             Manning uniform flow + Froude regime (valenx-openchannel).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.show_weir_workbench, "Weir Flow")
                        .on_hover_text(
                            "Show / hide the right-side Weir Flow Workbench — native \
                             sharp-crested rectangular / V-notch discharge (valenx-weir).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(
                            &mut self.show_thermocycle_workbench,
                            "Thermodynamic Cycle",
                        )
                        .on_hover_text(
                            "Show / hide the right-side Thermodynamic Cycle Workbench — native \
                             ideal Carnot / Otto / Diesel / Brayton / Rankine efficiency (valenx-thermocycle).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.show_queueing_workbench, "Queueing (M/M/1)")
                        .on_hover_text(
                            "Show / hide the right-side Queueing Workbench — native \
                             M/M/1 steady-state queue metrics (valenx-queueing).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(
                            &mut self.show_radioactivity_workbench,
                            "Radioactive Decay",
                        )
                        .on_hover_text(
                            "Show / hide the right-side Radioactive Decay Workbench — native \
                             single-nuclide exponential decay + activity (valenx-radioactivity).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.show_osmosis_workbench, "Osmosis / Starling")
                        .on_hover_text(
                            "Show / hide the right-side Osmosis Workbench — native van't Hoff \
                             osmotic pressure + Starling filtration (valenx-osmosis).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.show_thermoreg_workbench, "Thermoregulation")
                        .on_hover_text(
                            "Show / hide the right-side Thermoregulation Workbench — native \
                             single-node human heat balance (valenx-thermoreg).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.show_hemodynamics_workbench, "Hemodynamics")
                        .on_hover_text(
                            "Show / hide the right-side Hemodynamics Workbench — native \
                             cardiac output, Poiseuille flow, Windkessel (valenx-hemodynamics).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(
                            &mut self.show_popdynamics_workbench,
                            "Population Dynamics",
                        )
                        .on_hover_text(
                            "Show / hide the right-side Population Dynamics Workbench — native \
                             SIR / logistic / Lotka-Volterra closed forms (valenx-popdynamics).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Rail / Train Workbench — native
                    // train resistance + tractive effort (valenx-rail). Off by default.
                    if ui
                        .checkbox(&mut self.show_rail_workbench, "Rail / Train")
                        .on_hover_text(
                            "Show / hide the right-side Rail / Train Workbench — native train \
                             running resistance (Davis A + B·v + C·v²), grade resistance, net \
                             tractive force, acceleration, drawbar power, and the constant-effort \
                             balancing speed, computed in-process by valenx-rail.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Bone Mechanics workbench — native
                    // long-bone stress / strength (valenx-bonemech).
                    if ui
                        .checkbox(&mut self.show_bonemech_workbench, "Bone Mechanics")
                        .on_hover_text(
                            "Show / hide the right-side Bone Mechanics Workbench — native \
                             long-bone hollow-shaft bending / axial stress and density-scaled \
                             strength, computed in-process by valenx-bonemech.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Chain Drive workbench — native
                    // roller-chain kinematics (valenx-chaindrive).
                    if ui
                        .checkbox(&mut self.show_chaindrive_workbench, "Chain Drive")
                        .on_hover_text(
                            "Show / hide the right-side Chain Drive Workbench — native \
                             single-stage roller-chain kinematics (speed ratio, chain velocity, \
                             output torque, chain length), computed in-process by valenx-chaindrive.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Clutch workbench — native friction
                    // clutch torque capacity (valenx-clutch).
                    if ui
                        .checkbox(&mut self.show_clutch_workbench, "Clutch")
                        .on_hover_text(
                            "Show / hide the right-side Clutch Workbench — native dry \
                             friction-clutch torque capacity (uniform-wear vs uniform-pressure) \
                             and transmissible power, computed in-process by valenx-clutch.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Solenoid Coil workbench — native
                    // long-solenoid field / inductance (valenx-coil).
                    if ui
                        .checkbox(&mut self.show_coil_workbench, "Solenoid Coil")
                        .on_hover_text(
                            "Show / hide the right-side Solenoid Coil Workbench — native \
                             long-solenoid axial field, inductance and reactance, computed \
                             in-process by valenx-coil.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Steel Column workbench — native
                    // Euler-Johnson axial buckling (valenx-columnsteel).
                    if ui
                        .checkbox(&mut self.show_columnsteel_workbench, "Steel Column")
                        .on_hover_text(
                            "Show / hide the right-side Steel Column Workbench — native \
                             Euler-Johnson (AISC-ASD) axial buckling: slenderness, critical \
                             stress and allowable load, computed in-process by valenx-columnsteel.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Collision Workbench — native AABB
                    // geometry + overlap tests (valenx-collision). Off by default.
                    if ui
                        .checkbox(&mut self.show_collision_workbench, "Collision")
                        .on_hover_text(
                            "Show / hide the right-side Collision Workbench — native AABB \
                             geometry: per-box volume / surface / diagonal / inradius, plus the \
                             pairwise overlap test and L2 separation, computed in-process by \
                             valenx-collision.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Statics workbench (valenx-statics).
                    if ui
                        .checkbox(&mut self.show_statics_workbench, "Statics")
                        .on_hover_text(
                            "Show / hide the Statics Workbench — native simply-supported beam \
                             reactions (ΣF=0, ΣM=0), computed in-process by valenx-statics.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Projectile workbench (valenx-projectile).
                    if ui
                        .checkbox(&mut self.show_projectile_workbench, "Projectile")
                        .on_hover_text(
                            "Show / hide the Projectile Workbench — native ballistic trajectory \
                             (vacuum range / apex / time-of-flight + drag), computed in-process \
                             by valenx-projectile.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Conveyor workbench (valenx-conveyor).
                    if ui
                        .checkbox(&mut self.show_conveyor_workbench, "Conveyor")
                        .on_hover_text(
                            "Show / hide the Conveyor Workbench — native belt-conveyor mass flow, \
                             capacity and lift / drive power, computed in-process by valenx-conveyor.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Fluid Statics workbench (valenx-fluid-statics).
                    if ui
                        .checkbox(&mut self.show_fluidstatics_workbench, "Fluid Statics")
                        .on_hover_text(
                            "Show / hide the Fluid Statics Workbench — native hydrostatic pressure \
                             and submerged-plate resultant force / centre of pressure, computed \
                             in-process by valenx-fluid-statics.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Plate Bending workbench (valenx-plate).
                    if ui
                        .checkbox(&mut self.show_plate_workbench, "Plate Bending")
                        .on_hover_text(
                            "Show / hide the Plate Bending Workbench — native thin circular-plate \
                             flexural rigidity, centre deflection and max bending stress, computed \
                             in-process by valenx-plate.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Strain Rosette workbench (valenx-strainrosette).
                    if ui
                        .checkbox(&mut self.show_strainrosette_workbench, "Strain Rosette")
                        .on_hover_text(
                            "Show / hide the Strain Rosette Workbench — native 0/45/90 rosette \
                             reduction to Cartesian + principal strains and plane-stress, computed \
                             in-process by valenx-strainrosette.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Transformer workbench (valenx-transformer).
                    if ui
                        .checkbox(&mut self.show_transformer_workbench, "Transformer")
                        .on_hover_text(
                            "Show / hide the Transformer Workbench — native ideal two-winding \
                             turns-ratio, EMF equation and efficiency, computed in-process by \
                             valenx-transformer.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Three-Phase workbench (valenx-threephase).
                    if ui
                        .checkbox(&mut self.show_threephase_workbench, "Three-Phase")
                        .on_hover_text(
                            "Show / hide the Three-Phase Workbench — native balanced wye / delta \
                             line-phase voltages / currents and the power triangle, computed \
                             in-process by valenx-threephase.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Truss Workbench — native planar
                    // pin-jointed Warren truss analysis (valenx-truss). Off by default.
                    if ui
                        .checkbox(&mut self.show_truss_workbench, "Truss")
                        .on_hover_text(
                            "Show / hide the right-side Truss Workbench — native planar \
                             pin-jointed Warren truss: determinacy, support reactions, and the \
                             peak tension / compression members by the method of joints, computed \
                             in-process by valenx-truss.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Solar PV Workbench — native
                    // single-diode photovoltaic cell performance (valenx-solarpv).
                    if ui
                        .checkbox(&mut self.show_solarpv_workbench, "Solar PV")
                        .on_hover_text(
                            "Show / hide the right-side Solar PV Workbench — native single-diode \
                             cell performance: short-circuit current, open-circuit voltage, the \
                             maximum-power point, fill factor and conversion efficiency, computed \
                             in-process by valenx-solarpv.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Sheet Metal Workbench — native bend
                    // allowance / deduction (valenx-sheet-metal). Off by default.
                    if ui
                        .checkbox(&mut self.show_sheetmetal_workbench, "Sheet Metal")
                        .on_hover_text(
                            "Show / hide the right-side Sheet Metal Workbench — native bend \
                             allowance and bend deduction (neutral-axis arc, flat-blank \
                             correction) from thickness, inside radius, angle, and k-factor, \
                             computed in-process by valenx-sheet-metal.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Field Statistics Workbench — descriptive
                    // statistics over a pasted number list (valenx-fields). Off by default.
                    if ui
                        .checkbox(&mut self.show_fields_workbench, "Field Statistics")
                        .on_hover_text(
                            "Show / hide the right-side Field Statistics Workbench — descriptive \
                             statistics (mean, median, variance, std dev, rms, skewness, excess \
                             kurtosis, coefficient of variation, min/max) over a pasted list of \
                             numbers, computed in-process by valenx-fields.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Gearbox workbench — native two-stage
                    // compound gear-train analysis (valenx-gearbox).
                    if ui
                        .checkbox(&mut self.show_gearbox_workbench, "Gearbox")
                        .on_hover_text(
                            "Show / hide the right-side Gearbox Workbench — native two-stage \
                             compound gear train: overall ratio, efficiency, output speed and \
                             torque and the input / output power, computed in-process by \
                             valenx-gearbox.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Science batch 5 workbenches. Each off by default.
                    if ui
                        .checkbox(&mut self.show_camdynamics_workbench, "Cam Dynamics")
                        .on_hover_text(
                            "Show / hide the right-side Cam Dynamics Workbench — native \
                             cam-follower rise kinematics (valenx-camdynamics).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.show_batteryecm_workbench, "Battery ECM")
                        .on_hover_text(
                            "Show / hide the right-side Battery ECM Workbench — native \
                             first-order Thevenin cell terminal voltage (valenx-battery-ecm).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.show_diffusion_workbench, "Diffusion")
                        .on_hover_text(
                            "Show / hide the right-side Diffusion Workbench — native Fickian \
                             flux + Gaussian spread (valenx-diffusion).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(
                            &mut self.show_dimensional_workbench,
                            "Dimensionless Numbers",
                        )
                        .on_hover_text(
                            "Show / hide the right-side Dimensionless Numbers Workbench — native \
                             similitude groups + regime classifiers (valenx-dimensional).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.show_fft_workbench, "FFT / Spectrum")
                        .on_hover_text(
                            "Show / hide the right-side FFT / Spectrum Workbench — native DFT \
                             of a synthesized test tone (valenx-fft).",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Fasteners Workbench — ISO 4017 hex-bolt
                    // dimensions (valenx-fasteners). Off by default.
                    if ui
                        .checkbox(&mut self.show_fasteners_workbench, "Fasteners")
                        .on_hover_text(
                            "Show / hide the right-side Fasteners Workbench — ISO 4017 hex-bolt \
                             dimensions: width across flats, head height, pitch diameter, and \
                             tensile stress area, for a standard metric size, from \
                             valenx-fasteners.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Fixed-Wing / Aircraft Workbench —
                    // native preliminary aircraft performance (valenx-fixedwing).
                    if ui
                        .checkbox(&mut self.show_fixedwing_workbench, "Fixed-Wing / Aircraft")
                        .on_hover_text(
                            "Show / hide the right-side Fixed-Wing / Aircraft Workbench — native \
                             preliminary point-performance: wing loading, stall speed, cruise \
                             lift / drag, the best lift-to-drag ratio and the still-air glide \
                             range, computed in-process by valenx-fixedwing.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Combustion workbench (valenx-combustion).
                    if ui
                        .checkbox(&mut self.show_combustion_workbench, "Combustion")
                        .on_hover_text(
                            "Show / hide the Combustion Workbench — native hydrocarbon air-fuel \
                             stoichiometry, equivalence ratio, product moles and adiabatic flame \
                             temperature, computed in-process by valenx-combustion.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Flywheel workbench (valenx-flywheel).
                    if ui
                        .checkbox(&mut self.show_flywheel_workbench, "Flywheel")
                        .on_hover_text(
                            "Show / hide the Flywheel Workbench — native rotor inertia, stored \
                             kinetic energy, usable energy, rim speed and hoop stress, computed \
                             in-process by valenx-flywheel.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Fracture Mechanics workbench (valenx-fracture).
                    if ui
                        .checkbox(&mut self.show_fracture_workbench, "Fracture Mechanics")
                        .on_hover_text(
                            "Show / hide the Fracture Mechanics Workbench — native Mode-I LEFM: \
                             stress intensity K, critical crack length, fracture stress and \
                             plastic-zone radius, computed in-process by valenx-fracture.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Hydraulics workbench (valenx-hydraulics).
                    if ui
                        .checkbox(&mut self.show_hydraulics_workbench, "Hydraulics")
                        .on_hover_text(
                            "Show / hide the Hydraulics Workbench — native cylinder / circuit \
                             force (p·A), flow (A·v), hydraulic power and valve flow, computed \
                             in-process by valenx-hydraulics.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Inclined Plane workbench (valenx-inclinedplane).
                    if ui
                        .checkbox(&mut self.show_inclinedplane_workbench, "Inclined Plane")
                        .on_hover_text(
                            "Show / hide the Inclined Plane Workbench — native ramp statics: \
                             mechanical advantage, normal / slope / friction forces, effort and \
                             self-locking, computed in-process by valenx-inclinedplane.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Insulation workbench (valenx-insulation).
                    if ui
                        .checkbox(&mut self.show_insulation_workbench, "Insulation")
                        .on_hover_text(
                            "Show / hide the Insulation Workbench — native composite-wall R-value, \
                             U-value and steady heat loss, computed in-process by valenx-insulation.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Lead Screw workbench (valenx-leadscrew).
                    if ui
                        .checkbox(&mut self.show_leadscrew_workbench, "Lead Screw")
                        .on_hover_text(
                            "Show / hide the Lead Screw Workbench — native power-screw lead angle, \
                             linear speed, thrust / torque, resolution and self-locking, computed \
                             in-process by valenx-leadscrew.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Lever workbench (valenx-leverage).
                    if ui
                        .checkbox(&mut self.show_leverage_workbench, "Lever")
                        .on_hover_text(
                            "Show / hide the Lever Workbench — native lever mechanical advantage, \
                             class, balancing force and displacement, computed in-process by \
                             valenx-leverage.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Mohr's Circle workbench (valenx-mohr).
                    if ui
                        .checkbox(&mut self.show_mohr_workbench, "Mohr's Circle")
                        .on_hover_text(
                            "Show / hide the Mohr's Circle Workbench — native 2-D stress \
                             transformation: principal stresses, max shear, principal angle and \
                             stress on a plane, computed in-process by valenx-mohr.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side MOSFET workbench (valenx-mosfet).
                    if ui
                        .checkbox(&mut self.show_mosfet_workbench, "MOSFET")
                        .on_hover_text(
                            "Show / hide the MOSFET Workbench — native square-law NMOS bias: \
                             overdrive, operating region, drain current and transconductance, \
                             computed in-process by valenx-mosfet.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Optics workbench (valenx-optics).
                    if ui
                        .checkbox(&mut self.show_optics_workbench, "Optics")
                        .on_hover_text(
                            "Show / hide the Optics Workbench — native geometric optics: lensmaker \
                             focal length, Snell refraction / TIR and thin-lens imaging, computed \
                             in-process by valenx-optics.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Orifice Meter workbench (valenx-orifice).
                    if ui
                        .checkbox(&mut self.show_orifice_workbench, "Orifice Meter")
                        .on_hover_text(
                            "Show / hide the Orifice Meter Workbench — native differential-pressure \
                             flow metering: beta ratio, discharge coefficient, volumetric / mass \
                             flow and permanent loss, computed in-process by valenx-orifice.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Pressure Vessel workbench (valenx-pressure-vessel).
                    if ui
                        .checkbox(&mut self.show_pressurevessel_workbench, "Pressure Vessel")
                        .on_hover_text(
                            "Show / hide the Pressure Vessel Workbench — native thin / thick-wall \
                             hoop, longitudinal, shear and von Mises stresses for cylinders and \
                             spheres, computed in-process by valenx-pressure-vessel.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Torsion workbench (valenx-torsion).
                    if ui
                        .checkbox(&mut self.show_torsion_workbench, "Torsion")
                        .on_hover_text(
                            "Show / hide the Torsion Workbench — native circular-shaft torsion: \
                             polar moment, shear stress, allowable torque, angle of twist and \
                             power, computed in-process by valenx-torsion.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Refrigeration workbench (valenx-refrigeration).
                    if ui
                        .checkbox(&mut self.show_refrigeration_workbench, "Refrigeration")
                        .on_hover_text(
                            "Show / hide the Refrigeration Workbench — native vapor-compression \
                             cycle COP, Carnot bound and second-law efficiency, computed in-process \
                             by valenx-refrigeration.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side DC Motor Workbench — native
                    // brushed-DC-motor performance (valenx-dcmotor). Off by default.
                    if ui
                        .checkbox(&mut self.show_dcmotor_workbench, "DC Motor")
                        .on_hover_text(
                            "Show / hide the right-side DC Motor Workbench — native brushed-DC-\
                             motor performance: stall torque / current, no-load speed, peak \
                             output power, and the operating-point current / torque / power / \
                             efficiency, computed in-process by valenx-dcmotor.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Frames Workbench — structural
                    // cross-section properties (valenx-frames). Off by default.
                    if ui
                        .checkbox(&mut self.show_frames_workbench, "Frames")
                        .on_hover_text(
                            "Show / hide the right-side Frames Workbench — structural \
                             cross-section properties: the area and perimeter of an I-beam, \
                             channel, angle, rectangular / round HSS, or T-beam profile, \
                             computed in-process by valenx-frames.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Thermal Expansion workbench — linear
                    // expansion + constrained stress (valenx-thermalexpansion). Off by default.
                    if ui
                        .checkbox(&mut self.show_thermalexpansion_workbench, "Thermal Expansion")
                        .on_hover_text(
                            "Show / hide the right-side Thermal Expansion workbench — linear / \
                             areal / volumetric expansion, the length change and free thermal \
                             strain, and the fully-constrained thermal stress E·α·ΔT, computed \
                             in-process by valenx-thermalexpansion.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Gas Dynamics workbench — 1-D
                    // compressible-flow relations (valenx-gasdynamics). Off by default.
                    if ui
                        .checkbox(&mut self.show_gasdynamics_workbench, "Gas Dynamics")
                        .on_hover_text(
                            "Show / hide the right-side Gas Dynamics workbench — 1-D \
                             compressible-flow relations for a perfect gas: isentropic \
                             stagnation + area-Mach ratios, normal-shock jumps, \
                             Prandtl-Meyer expansion, and Fanno / Rayleigh flow, computed \
                             in-process by valenx-gasdynamics.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.show_neuro_workbench, "Neural Interface")
                        .on_hover_text(
                            "Show / hide the right-side Neural-Interface workbench — native \
                             BCI stimulation: an electrode's extracellular field drives \
                             Hodgkin-Huxley axons (recruitment, bioheat, electrode impedance), \
                             solved in-process by valenx-neuro.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.show_windturbine_workbench, "Wind Turbine")
                        .on_hover_text(
                            "Show / hide the right-side Wind Turbine workbench — native \
                             actuator-disc power: swept area, wind / Betz / Cp-extracted power, \
                             the cut-in / rated / cut-out power curve, capacity factor and \
                             tip-speed ratio, computed in-process by valenx-windturbine.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.show_cad_workbench, "Parametric CAD")
                        .on_hover_text(
                            "Show / hide the right-side Parametric-CAD workbench — named \
                             parameters with expressions (Fusion-style) driving sketch \
                             geometry through the valenx-solvespace-3d constraint solver.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Antenna workbench — native parabolic
                    // dish gain / beamwidth (valenx-antenna).
                    if ui
                        .checkbox(&mut self.show_antenna_workbench, "Antenna")
                        .on_hover_text(
                            "Show / hide the right-side Antenna Workbench — native parabolic \
                             dish: wavelength, effective aperture, gain (linear and dBi) and a \
                             half-power beamwidth estimate, computed in-process by \
                             valenx-antenna.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.show_draft2d_workbench, "2D Drafting")
                        .on_hover_text(
                            "Show / hide the right-side 2D Drafting workbench — a \
                             LibreCAD-style 2D drawing canvas (lines, circles, arcs, \
                             polylines, text) with DXF export, on valenx-librecad-2d.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.show_reinforcement_workbench, "Reinforcement")
                        .on_hover_text(
                            "Show / hide the right-side Concrete-Reinforcement workbench — \
                             parametric beam/column rebar cages rendered in the 3D \
                             viewport, on valenx-reinforcement.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.show_render_workbench, "Path-Traced Render")
                        .on_hover_text(
                            "Show / hide the right-side Path-Traced Render workbench — a \
                             global-illumination Cornell-box render displayed as an image, \
                             on the native valenx-pathtrace CPU path tracer.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.show_beam_workbench, "Beam")
                        .on_hover_text(
                            "Show / hide the right-side Beam Workbench — native Euler-Bernoulli \
                             beam bending: peak deflection, bending moment and bending stress for \
                             a cantilever / simply-supported / fixed-fixed beam under a point or \
                             distributed load, computed in-process by valenx-beam.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.show_hvac_workbench, "HVAC")
                        .on_hover_text(
                            "Show / hide the right-side HVAC workbench — duct sizing and \
                             Darcy-Weisbach pressure-drop calculation, on valenx-hvac.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.show_reverse_workbench, "Reverse Engineering")
                        .on_hover_text(
                            "Show / hide the right-side Reverse-Engineering workbench — \
                             point-cloud → mesh reconstruction shown in the 3D viewport, \
                             on valenx-reverse.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.show_pump_workbench, "Pump")
                        .on_hover_text(
                            "Show / hide the right-side Pump Workbench — native centrifugal \
                             pump duty point (pump curve ∩ system curve), hydraulic and shaft \
                             power, and the NPSH cavitation margin, computed in-process by \
                             valenx-pump.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.show_interior_workbench, "Interior Design")
                        .on_hover_text(
                            "Show / hide the right-side Interior-Design workbench — a 2D \
                             floor plan with a furniture palette, on valenx-interior.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Heat Pump workbench — native
                    // Carnot COP + second-law derating (valenx-heatpump).
                    if ui
                        .checkbox(&mut self.show_heatpump_workbench, "Heat Pump")
                        .on_hover_text(
                            "Show / hide the right-side Heat Pump Workbench — native Carnot \
                             COP: temperature lift, the ideal heating / cooling COP and the \
                             real-machine COPs after a Carnot-fraction derating, computed \
                             in-process by valenx-heatpump.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.show_animate_workbench, "Animation")
                        .on_hover_text(
                            "Show / hide the right-side Animation workbench — a joint \
                             keyframe timeline with easing curves, on valenx-animate.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Battery Pack workbench — native
                    // series / parallel pack sizing (valenx-batterypack).
                    if ui
                        .checkbox(&mut self.show_batterypack_workbench, "Battery Pack")
                        .on_hover_text(
                            "Show / hide the right-side Battery Pack Workbench — native \
                             series / parallel pack sizing: pack voltage, capacity, energy, \
                             cell count, discharge current and runtime at a C-rate, computed \
                             in-process by valenx-batterypack.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    if ui
                        .checkbox(&mut self.show_variant_effect_workbench, "Variant Effect")
                        .on_hover_text(
                            "Show / hide the right-side Variant-Effect workbench — an HGVS \
                             variant parser (protein / coding substitutions), on \
                             valenx-variant-effect.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Pipe Flow workbench — native
                    // Darcy-Weisbach pipe-flow analysis (valenx-pipeflow).
                    if ui
                        .checkbox(&mut self.show_pipeflow_workbench, "Pipe Flow")
                        .on_hover_text(
                            "Show / hide the right-side Pipe Flow Workbench — native \
                             Darcy-Weisbach analysis: Reynolds number, flow regime, friction \
                             factor, head loss and pressure drop for a straight pipe run, \
                             computed in-process by valenx-pipeflow.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Astro / Launch workbench —
                    // the native launch-vehicle ascent + trajectory
                    // simulator (fly a rocket to orbit) plus the
                    // closed-form mission planners (Hohmann, hoverslam,
                    // rendezvous, launch azimuth). Off by default.
                    if ui
                        .checkbox(&mut self.show_astro_workbench, "Astro / Launch")
                        .on_hover_text(
                            "Show / hide the right-side Astro / Launch workbench — a \
                             launch-vehicle ascent simulator (orbit reached, Δv \
                             budget, max-Q, staging timeline + flight-profile chart) \
                             plus Hohmann / hoverslam / rendezvous / launch-azimuth \
                             planners.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Heat Exchanger workbench — native
                    // effectiveness-NTU exchanger analysis (valenx-heatexchanger).
                    if ui
                        .checkbox(&mut self.show_heatexchanger_workbench, "Heat Exchanger")
                        .on_hover_text(
                            "Show / hide the right-side Heat Exchanger Workbench — native \
                             effectiveness-NTU analysis: capacity ratio, NTU, effectiveness, \
                             duty and both outlet temperatures for counter- or parallel-flow, \
                             computed in-process by valenx-heatexchanger.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Car workbench — vehicle dynamics
                    // (top speed, 0-100 / 0-200 km/h, braking, skidpad) over
                    // valenx-vehicle. Off by default.
                    if ui
                        .checkbox(&mut self.show_car_workbench, "Car")
                        .on_hover_text(
                            "Show / hide the right-side Car workbench — a design → \
                             simulate panel: top speed, 0-100 / 0-200 km/h, braking \
                             distance, and skidpad grip over valenx-vehicle.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Rocket workbench — the coupled
                    // design→simulate pipeline (valenx-rocket-demo): fly the
                    // ascent to orbit, then size the interstage struts against
                    // the trajectory's peak g-load with a live SAFE /
                    // OVER-STRESSED verdict. Off by default.
                    if ui
                        .checkbox(
                            &mut self.show_rocket_workbench,
                            "Rocket — design → simulate",
                        )
                        .on_hover_text(
                            "Show / hide the right-side Rocket workbench — a reactive \
                             design→simulate loop (valenx-rocket-demo): fly the \
                             medium-lift two-stage preset to orbit, then size the \
                             interstage struts against the trajectory's peak g-load \
                             with a live safety-factor verdict.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // Toggle the right-side Engine workbench — the reactive
                    // engine design → analyze → optimize → export loop
                    // (ideal-nozzle performance + Bartz regen-cooling).
                    if ui
                        .checkbox(
                            &mut self.show_engine_workbench,
                            "Engine — design → analyze",
                        )
                        .on_hover_text(
                            "Show / hide the right-side Engine workbench — reactive \
                             ideal-nozzle performance + Bartz regen-cooling, a \
                             chamber-pressure × expansion-ratio optimizer, and a 3-D \
                             nozzle you can export to STL.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // The Assistant activity sidebar — a live in-app feed of
                    // what the AI assistant is designing / simulating /
                    // building, read from a .jsonl activity file. On by
                    // default so the app narrates its own work.
                    if ui
                        .checkbox(&mut self.show_assistant_panel, "Assistant (AI activity)")
                        .on_hover_text(
                            "Show / hide the right-side Assistant sidebar — a live \
                             feed of what the AI assistant is designing / simulating / \
                             building, read from a .jsonl activity file.",
                        )
                        .changed()
                    {
                        ui.close_menu();
                    }
                    // ── Central viewport selector ─────────────────────────
                    // Lets the user swap the 3D viewport for a
                    // domain-appropriate alternative (e.g. the 2D DNA /
                    // plasmid view for biology work) without losing the
                    // 3D camera state.
                    ui.separator();
                    ui.label(egui::RichText::new("Central viewport:").weak());
                    if ui
                        .radio(
                            self.active_viewport == ViewportKind::Viewport3D,
                            ViewportKind::Viewport3D.label(),
                        )
                        .on_hover_text("Blender-style 3D orbit / pan / zoom viewport (CAD, mesh, FEM, aero, astro).")
                        .clicked()
                    {
                        self.active_viewport = ViewportKind::Viewport3D;
                        ui.close_menu();
                    }
                    if ui
                        .radio(
                            self.active_viewport == ViewportKind::Viewport2dDna,
                            ViewportKind::Viewport2dDna.label(),
                        )
                        .on_hover_text("2D linear + circular annotated-sequence / plasmid map (Genetics Workbench).")
                        .clicked()
                    {
                        self.active_viewport = ViewportKind::Viewport2dDna;
                        ui.close_menu();
                    }
                    });
                });
                // ── Part Design ─────────────────────────────────────────────
                // Surfaces valenx's in-house parametric CAD modeler (the
                // valenx-solvespace-3d feature tree) as a first-class top-bar
                // menu, so the Extrude / Revolve / primitive / boolean
                // operations are discoverable without first hunting for the
                // View → Parametric CAD toggle. Every item is a clearly named,
                // AI-drivable widget; each reveals the workbench, mutates the
                // shared feature tree through the same `CadWorkbenchState`
                // methods the panel's own buttons call, then requests a rebuild
                // so the central 3-D viewport updates immediately.
                ui.menu_button("Part Design", |ui| {
                    crate::menu_ui::scrollable_menu(ui, |ui| {
                        if ui
                            .button("Open Part Design workbench")
                            .on_hover_text(
                                "Reveal the right-side parametric CAD modeler — named-parameter \
                                 feature tree (sketch, extrude, revolve, primitives, booleans), \
                                 built in-house on valenx-solvespace-3d + valenx-cad.",
                            )
                            .clicked()
                        {
                            self.show_cad_workbench = true;
                            ui.close_menu();
                        }
                        if ui
                            .button("New part (reset)")
                            .on_hover_text(
                                "Reset the feature tree to a single base solid (a parametric \
                                 box) — a clean start for a new part.",
                            )
                            .clicked()
                        {
                            self.show_cad_workbench = true;
                            self.cad.reset_to_base_solid();
                            self.cad.request_rebuild();
                            ui.close_menu();
                        }
                        ui.separator();
                        ui.label(egui::RichText::new("Add feature").weak().small());
                        if ui
                            .button("Sketch (extrude profile)")
                            .on_hover_text(
                                "Add an Extrude feature whose editable (x, y) profile IS the \
                                 sketch — the in-house sketch → extrude → 3-D-part workflow.",
                            )
                            .clicked()
                        {
                            self.show_cad_workbench = true;
                            self.cad.add_extrude();
                            self.cad.request_rebuild();
                            ui.close_menu();
                        }
                        if ui
                            .button("Extrude")
                            .on_hover_text("Sweep a 2-D profile along +Z into a solid prism.")
                            .clicked()
                        {
                            self.show_cad_workbench = true;
                            self.cad.add_extrude();
                            self.cad.request_rebuild();
                            ui.close_menu();
                        }
                        if ui
                            .button("Revolve")
                            .on_hover_text(
                                "Revolve a half-section profile about the Z axis (lathe / turn) \
                                 into a solid of revolution.",
                            )
                            .clicked()
                        {
                            self.show_cad_workbench = true;
                            self.cad.add_revolve();
                            self.cad.request_rebuild();
                            ui.close_menu();
                        }
                        if ui
                            .button("Box")
                            .on_hover_text("Add an axis-aligned box primitive.")
                            .clicked()
                        {
                            self.show_cad_workbench = true;
                            self.cad.add_box();
                            self.cad.request_rebuild();
                            ui.close_menu();
                        }
                        if ui
                            .button("Cylinder")
                            .on_hover_text("Add a cylinder primitive (axis along +Z).")
                            .clicked()
                        {
                            self.show_cad_workbench = true;
                            self.cad.add_cylinder();
                            self.cad.request_rebuild();
                            ui.close_menu();
                        }
                        if ui
                            .button("Sphere")
                            .on_hover_text("Add a sphere primitive.")
                            .clicked()
                        {
                            self.show_cad_workbench = true;
                            self.cad.add_sphere();
                            self.cad.request_rebuild();
                            ui.close_menu();
                        }
                        if ui
                            .button("Cone")
                            .on_hover_text("Add a (truncated) cone primitive.")
                            .clicked()
                        {
                            self.show_cad_workbench = true;
                            self.cad.add_cone();
                            self.cad.request_rebuild();
                            ui.close_menu();
                        }
                        if ui
                            .button("Torus")
                            .on_hover_text("Add a torus primitive (major axis along Z).")
                            .clicked()
                        {
                            self.show_cad_workbench = true;
                            self.cad.add_torus();
                            self.cad.request_rebuild();
                            ui.close_menu();
                        }
                        ui.separator();
                        ui.menu_button("Boolean", |ui| {
                            ui.label(
                                egui::RichText::new("combine the last feature with the body")
                                    .weak()
                                    .small(),
                            );
                            let has = self.cad.has_steps();
                            if ui
                                .add_enabled(has, egui::Button::new("Union (join)"))
                                .on_hover_text("Weld the last feature onto the running body.")
                                .clicked()
                            {
                                self.show_cad_workbench = true;
                                self.cad.set_last_op_union();
                                self.cad.request_rebuild();
                                ui.close_menu();
                            }
                            if ui
                                .add_enabled(has, egui::Button::new("Cut (subtract)"))
                                .on_hover_text("Subtract the last feature from the running body.")
                                .clicked()
                            {
                                self.show_cad_workbench = true;
                                self.cad.set_last_op_cut();
                                self.cad.request_rebuild();
                                ui.close_menu();
                            }
                            if ui
                                .add_enabled(has, egui::Button::new("Intersect"))
                                .on_hover_text("Keep only the overlap of the last feature and body.")
                                .clicked()
                            {
                                self.show_cad_workbench = true;
                                self.cad.set_last_op_intersect();
                                self.cad.request_rebuild();
                                ui.close_menu();
                            }
                        });
                        ui.separator();
                        if ui
                            .button("Save part…")
                            .on_hover_text(
                                "Open the Part Design workbench, where the feature tree can be \
                                 saved to a .ron document (and exported to STL).",
                            )
                            .clicked()
                        {
                            self.show_cad_workbench = true;
                            ui.close_menu();
                        }
                    });
                });
                ui.menu_button(cat.lookup("menu.run"), |ui| {
                    let running = self.run_handle.is_some();
                    let has_selection = self.selected_case.is_some();
                    if ui
                        .add_enabled(
                            has_selection && !running,
                            egui::Button::new(cat.lookup("menu.run.selected")),
                        )
                        .clicked()
                    {
                        ui.close_menu();
                        self.run_selected_case();
                    }
                    if ui
                        .add_enabled(
                            has_selection && !running,
                            egui::Button::new(cat.lookup("menu.run.prepare")),
                        )
                        .on_hover_text(cat.lookup("menu.run.prepare-tooltip"))
                        .clicked()
                    {
                        ui.close_menu();
                        self.prepare_selected_case();
                    }
                    let has_prepared = self.last_prepared_job.is_some();
                    if ui
                        .add_enabled(
                            has_prepared && !running,
                            egui::Button::new(cat.lookup("menu.run.from-prepared")),
                        )
                        .on_hover_text(cat.lookup("menu.run.from-prepared-tooltip"))
                        .clicked()
                    {
                        ui.close_menu();
                        self.run_from_prepared_workdir();
                    }
                    if ui
                        .add_enabled(
                            has_selection && !running,
                            egui::Button::new(cat.lookup("menu.run.sweep")),
                        )
                        .on_hover_text(cat.lookup("menu.run.sweep-tooltip"))
                        .clicked()
                    {
                        ui.close_menu();
                        self.sweep_selected_case();
                    }
                    if ui
                        .add_enabled(!running, egui::Button::new(cat.lookup("menu.run.first")))
                        .clicked()
                    {
                        ui.close_menu();
                        self.run_first_case();
                    }
                    if ui
                        .add_enabled(running, egui::Button::new(cat.lookup("menu.run.cancel")))
                        .clicked()
                    {
                        ui.close_menu();
                        self.cancel_run();
                    }
                });
                ui.menu_button(cat.lookup("menu.settings"), |ui| {
                    if ui.button(cat.lookup("menu.settings.preferences")).clicked() {
                        self.settings_open = true;
                        ui.close_menu();
                    }
                    if ui.button(cat.lookup("menu.settings.reprobe")).clicked() {
                        ui.close_menu();
                        self.reprobe();
                    }
                });
                // Tools menu (Phase 44b) — every registered adapter
                // gets a "New case from this adapter" entry, grouped
                // by the first segment of its `ribbon_contributions`
                // (`bio.dnachisel.*` → `BIO`). Clicking an entry
                // triggers `new_case_for_adapter`. We borrow `self`
                // immutably to iterate the registry, so writes are
                // deferred via `pending_new_case` until after the
                // closure ends.
                let mut pending_new_case: Option<String> = None;
                let show_non_oss = self.settings.show_non_oss_adapters;
                ui.menu_button("Tools", |ui| {
                    use std::collections::BTreeMap;
                    let mut by_category: BTreeMap<&str, Vec<&valenx_core::registry::AdapterEntry>> =
                        BTreeMap::new();
                    for (_, entry) in self.registry.iter() {
                        // Default-hide adapters wrapping
                        // academic-only / non-commercial tools.
                        if !show_non_oss
                            && !valenx_core::adapter_helpers::is_oss_license(
                                entry.adapter.info().tool_license,
                            )
                        {
                            continue;
                        }
                        let caps = entry.adapter.capabilities();
                        let category = caps
                            .ribbon_contributions
                            .iter()
                            .filter_map(|c| c.split('.').next())
                            .next()
                            .unwrap_or("uncategorised");
                        by_category.entry(category).or_default().push(entry);
                    }
                    crate::menu_ui::scrollable_menu(ui, |ui| {
                        for (cat_name, mut entries) in by_category {
                            entries.sort_by_key(|e| e.adapter.info().display_name);
                            ui.menu_button(cat_name.to_uppercase(), |ui| {
                                crate::menu_ui::scrollable_menu(ui, |ui| {
                                    for entry in entries {
                                        let info = entry.adapter.info();
                                        if ui.button(info.display_name).clicked() {
                                            pending_new_case = Some(info.id.to_string());
                                            ui.close_menu();
                                        }
                                    }
                                });
                            });
                        }
                    });
                });
                if let Some(id) = pending_new_case.take() {
                    self.new_case_for_adapter(&id);
                }
                // Phase 22 — Add-on Manager menu entry. Opens / closes
                // the dedicated panel.
                ui.menu_button("Add-ons", |ui| {
                    if ui
                        .button(format!(
                            "Manager ({} installed)",
                            self.addons.installed().len()
                        ))
                        .clicked()
                    {
                        self.show_addon_manager = !self.show_addon_manager;
                        if self.show_addon_manager {
                            // Refresh on open so the list reflects any
                            // out-of-band install (e.g. user copied a
                            // folder into ~/.valenx/addons manually).
                            let _ = self.addons.refresh();
                        }
                        ui.close_menu();
                    }
                    if ui.button("Install from directory…").clicked() {
                        if let Some(dir) = rfd::FileDialog::new()
                            .set_title("Pick an add-on directory (must contain valenx-addon.toml)")
                            .pick_folder()
                        {
                            match self.addons.install_from_dir(&dir) {
                                Ok(local) => {
                                    self.status = Some(format!(
                                        "Installed add-on {} v{}",
                                        local.manifest.name, local.manifest.version,
                                    ));
                                }
                                Err(e) => {
                                    self.last_error = Some(format!("Add-on install: {e}"));
                                }
                            }
                        }
                        ui.close_menu();
                    }
                });

                // Phase 21 — Macro menu. Always rendered (regardless
                // of registry contents) because macros work without
                // any external adapters.
                ui.menu_button("Macro", |ui| {
                    let recording = self.macro_recorder.is_recording();
                    let action_count = self.macro_recorder.current_action_count();
                    if !recording {
                        if ui.button("Start Recording").clicked() {
                            self.macro_recorder.start_recording("Untitled");
                            ui.close_menu();
                        }
                    } else if ui
                        .button(format!("Stop Recording ({action_count} actions)"))
                        .clicked()
                    {
                        if let Some(m) = self.macro_recorder.stop_recording() {
                            if let Some(dir) = valenx_macro::persist::default_macro_dir() {
                                if let Ok(p) = valenx_macro::persist::save(&m, &dir) {
                                    self.status = Some(format!("Saved macro to {}", p.display()));
                                }
                            }
                        }
                        ui.close_menu();
                    }
                    if ui
                        .button(format!(
                            "Library ({} macros)",
                            self.macro_recorder.library.len()
                        ))
                        .clicked()
                    {
                        // Open the macro library panel via a status
                        // message; the dedicated panel lives in
                        // mesh_toolbox::macro_library and is opened by
                        // the user from the dock.
                        self.status =
                            Some("Macro Library panel: pick a macro to run/edit/delete.".into());
                        ui.close_menu();
                    }
                });

                ui.menu_button(cat.lookup("menu.help"), |ui| {
                    if ui.button(cat.lookup("menu.help.palette")).clicked() {
                        self.palette.request_open();
                        ui.close_menu();
                    }
                    if ui
                        .button("Keyboard shortcuts…")
                        .on_hover_text("Open the keyboard cheat-sheet overlay (? toggles).")
                        .clicked()
                    {
                        self.keyboard_help_open = !self.keyboard_help_open;
                        self.settings.keyboard_shortcuts_overlay_open = self.keyboard_help_open;
                        save_settings_to_state_dir(&self.settings);
                        ui.close_menu();
                    }
                    if ui
                        .button("Panel help (F1)")
                        .on_hover_text(
                            "Pop up contextual help for the focused workbench panel.",
                        )
                        .clicked()
                    {
                        self.dispatch_shortcut(ShortcutAction::OpenContextHelp);
                        ui.close_menu();
                    }
                    if ui
                        .button("Welcome tour…")
                        .on_hover_text("Restart the 3-step welcome walkthrough.")
                        .clicked()
                    {
                        self.welcome_tour_open = true;
                        self.welcome_tour_state = crate::welcome_tour::TourState::new();
                        ui.close_menu();
                    }
                    if ui.button(cat.lookup("menu.help.about")).clicked() {
                        self.about_open = true;
                        ui.close_menu();
                    }
                });
                ui.separator();
                ui.label(format!("Valenx {} — pre-alpha", env!("CARGO_PKG_VERSION")));
                // One-click panel toggles — "make space" controls pinned to
                // the far right of the ribbon (Blender / Fusion style).
                // A toggle is highlighted while its panel is visible.
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    let any_open = self.show_browser || self.show_mesh_toolbox;
                    if ui
                        .selectable_label(!any_open, "Max view")
                        .on_hover_text(
                            "Maximize viewport: hide both side panels \
                             (click again to restore both).",
                        )
                        .clicked()
                    {
                        let new_state = !any_open;
                        self.show_browser = new_state;
                        self.show_mesh_toolbox = new_state;
                    }
                    ui.separator();
                    if ui
                        .selectable_label(self.show_mesh_toolbox, "Toolbox")
                        .on_hover_text("Show / hide the right Mesh Toolbox panel.")
                        .clicked()
                    {
                        self.show_mesh_toolbox = !self.show_mesh_toolbox;
                    }
                    if ui
                        .selectable_label(self.show_browser, "Browser")
                        .on_hover_text("Show / hide the left Browser panel.")
                        .clicked()
                    {
                        self.show_browser = !self.show_browser;
                    }
                });
            });
        });

        // Agent-drives-valenx bridge: apply any new commands an external agent
        // appended to a channel's command file (open/focus/rename/close tabs,
        // switch a workbench, post a chat note) before the tab strip is drawn,
        // so those effects show up this frame. Throttled to ~1/sec internally;
        // clean &mut self here (no dock-tree borrow held).
        crate::agent_commands::poll_and_apply_agent_commands(self);

        // Project tabs (Chrome-style) — a slim strip just below the
        // ribbon. Each tab is a project workbench the user switches
        // between; additive, so it's empty until the user opens a tab.
        crate::project_tabs::draw_tab_strip(self, ctx);

        // Status bar with progress bar if running
        egui::TopBottomPanel::bottom("valenx_status").show(ctx, |ui| {
            ui.horizontal(|ui| {
                if self.run_handle.is_some() {
                    ui.add(
                        egui::ProgressBar::new(self.run_progress.clamp(0.0, 100.0) / 100.0)
                            .desired_width(160.0)
                            .show_percentage(),
                    );
                    ui.label(&self.run_message);
                } else if let Some(err) = &self.last_error {
                    ui.colored_label(egui::Color32::LIGHT_RED, err);
                } else if let Some(err) = &self.last_run_error {
                    ui.colored_label(egui::Color32::LIGHT_RED, format!("Last run failed: {err}"));
                } else if let Some(report) = &self.last_run_report {
                    ui.colored_label(
                        egui::Color32::LIGHT_GREEN,
                        format!(
                            "Run finished — exit {} · {} residual samples",
                            report.exit_code,
                            report.residual_history.len()
                        ),
                    );
                } else if let Some(msg) = &self.status {
                    ui.label(msg);
                } else {
                    ui.label(
                        "Ctrl+P for commands · Ctrl+N for a new project · Ctrl+O to open one.",
                    );
                }

                // Always show which case is selected, so the user
                // doesn't have to guess what Run / F5 will do.
                if let Some(case) = &self.selected_case {
                    ui.separator();
                    ui.label(
                        egui::RichText::new(format!("▶ {case}"))
                            .color(egui::Color32::from_rgb(180, 210, 240)),
                    );
                }
            });
        });

        // Tabbed bottom panel — Residuals or Log. Collapsible: when
        // `self.bottom_panel_collapsed` is set the panel shrinks to
        // just its header strip (tab selectors + toggle) and the
        // content body is skipped, freeing the vertical space. The
        // resize handle / tall default height only apply in the
        // expanded state — collapsed it auto-sizes to the bar.
        let collapsed = self.bottom_panel_collapsed;
        let bottom_dock = egui::TopBottomPanel::bottom("valenx_bottom_dock");
        let bottom_dock = if collapsed {
            // Thin, fixed bar — no resize handle, no minimum that would
            // keep reserving the tall content area.
            bottom_dock.resizable(false)
        } else {
            bottom_dock
                .resizable(true)
                .default_height(240.0)
                .height_range(140.0..=520.0)
        };
        bottom_dock.show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.selectable_value(
                    &mut self.bottom_tab,
                    BottomTab::Residuals,
                    format!("Residuals ({})", self.residuals.by_field.len()),
                );
                ui.selectable_value(
                    &mut self.bottom_tab,
                    BottomTab::Log,
                    format!("Log ({})", self.log.lines.len()),
                );
                ui.separator();
                if !collapsed && self.bottom_tab == BottomTab::Log {
                    log_panel::header(ui, &mut self.log);
                }
                // Collapse / expand toggle on the right edge. The
                // button's TEXT is the accessibility / UI-Automation
                // `Name`, so it is a uniquely- and stably-named,
                // AI-invokable widget: "Collapse panel" when expanded,
                // "Expand panel" when collapsed. egui derives a
                // button's accesskit name from its label text, so the
                // plain-ASCII phrase is what an external driver finds
                // and Invokes by Name.
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    let name = if collapsed {
                        "Expand panel"
                    } else {
                        "Collapse panel"
                    };
                    if ui.button(name).clicked() {
                        self.toggle_bottom_panel();
                    }
                });
            });

            // Collapsed: stop here — render only the header strip and
            // skip the residual plot / log / placeholder body so the
            // panel stays a thin bar.
            if collapsed {
                return;
            }

            ui.separator();
            match self.bottom_tab {
                BottomTab::Residuals => {
                    residuals::show(
                        ui,
                        &self.residuals,
                        self.settings.residual_scale,
                        self.settings.convergence_target,
                    );
                }
                BottomTab::Log => {
                    log_panel::show(ui, &self.log);
                }
            }
        });

        // Browser (left) — shown/hidden by `self.show_browser` (the
        // ribbon / View toggle). When shown it is ALSO collapsible to a
        // thin bar via `self.browser_collapsed`, mirroring the bottom
        // dock: collapsed renders only a narrow vertical strip holding
        // the AI-drivable "Expand panel" button and skips the heavy
        // browser body; expanded renders the full Browser plus a named
        // "Collapse panel" button in the header row. The resize handle /
        // wide default width only apply when expanded.
        if self.show_browser {
            let browser_collapsed = self.browser_collapsed;
            let browser = egui::SidePanel::left("valenx_browser");
            let browser = if browser_collapsed {
                // Thin, fixed strip — no resize handle, no wide default
                // that would keep reserving the browser's full width.
                browser.resizable(false).exact_width(30.0)
            } else {
                browser
                    .resizable(true)
                    .default_width(300.0)
                    .width_range(220.0..=520.0)
            };
            browser.show(ctx, |ui| {
                if browser_collapsed {
                    // Collapsed: only the named "Expand panel" button.
                    // Its label TEXT is the accessibility / UI-Automation
                    // `Name` (egui derives a button's accesskit name from
                    // its label, exactly like the bottom dock's toggle),
                    // so the full plain-ASCII phrase is what an external
                    // driver finds and Invokes by Name. egui clips the
                    // rendered glyphs to the narrow 30px strip while the
                    // accesskit Name stays the complete "Expand panel"
                    // string.
                    if ui.button("Expand panel").clicked() {
                        self.toggle_browser_panel();
                    }
                    return;
                }
                // Expanded: full Browser. The header row carries the
                // right-aligned, named "Collapse panel" toggle alongside
                // the "Browser" title.
                ui.horizontal(|ui| {
                    ui.heading("Browser");
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if ui.button("Collapse panel").clicked() {
                            self.toggle_browser_panel();
                        }
                    });
                });
                ui.separator();
                self.draw_browser(ui);
            });
        }

        // wgpu render state + DPI, grabbed once up-front so BOTH the dockable
        // workbench layout (below — its `workspace:<n>` tiles render live 3-D)
        // and the central viewport (further down) can build a `WgpuCtx`. The
        // clone is an Arc-backed handle (cheap) and owned, so it doesn't pin a
        // borrow on `frame`.
        let wgpu_render_state = frame.wgpu_render_state().cloned();
        let pixels_per_point = ctx.pixels_per_point();

        // Right-side workbench dispatch. Each open workbench is its own
        // `SidePanel::right`; the dispatch ALWAYS runs (in both classic and
        // dock modes) so that a workbench panel can sit on the right alongside
        // the dock in the centre. In a normal (non-product) tab at most the
        // active tab's one workbench flag is set, so exactly one panel renders;
        // in a `new_unit` product tab the tab's linked `workbench_kind` flag is
        // the only one on (see `project_tabs::sync_active`), so its single
        // workbench reserves the right edge BEFORE the CentralPanel/dock below.
        //
        // When the user opts into View → "Dockable panel layout (beta)", the
        // `egui_tiles` dock region is ALSO drawn (further down, before the
        // CentralPanel) — see [`crate::dock_layout::draw_dock_layout`]. The
        // SidePanel(s) are added before it, so egui reserves their right-hand
        // width first and the dock fills the remaining centre: workbench right,
        // dock centre, no overlap.
        {
            // The 10 [`crate::dock_layout::DOCKABLE_PANELS`] (mesh toolbox,
            // genetics, aero, fem, cfd, astro, rocket, engine, car, assistant)
            // are hosted INSIDE the dock tree when the dock is on AND this is
            // not a clean agent product tab (`!dock_agent_only`) — see
            // `draw_dock_layout`'s `open_ids`. In that one case they must NOT
            // also paint as classic right SidePanels (that would duplicate
            // them), so their draws below are gated on this. Every other
            // workbench — and these same panels in a `dock_agent_only` product
            // tab (whose dock hosts only `workspace:n | agent:n`) or in classic
            // mode — always paints as a SidePanel.
            let dock_owns_dockable_panels = self.dock_enabled && !self.dock_agent_only;

            // Mesh Toolbox (right) — only paints when a mesh / STL is
            // loaded AND the user hasn't hidden it. Mounted before the
            // CentralPanel so egui docks it to the right of the viewport.
            if !dock_owns_dockable_panels {
                crate::mesh_toolbox::draw_mesh_toolbox(self, ctx);
            }

            // Genetics Workbench (right) — the Round-6 computational-
            // biology panels. A no-op unless the user has toggled it on
            // via View → Genetics Workbench. Mounted before the
            // CentralPanel so egui docks it to the right (alongside the
            // Mesh Toolbox when both are open).
            if !dock_owns_dockable_panels {
                crate::genetics_workbench::draw_genetics_workbench(self, ctx);
            }

            // Aerodynamics / Wind Tunnel workbench (right) — the
            // valenx-aero virtual-wind-tunnel CFD panels. A no-op unless
            // toggled on via View → Aerodynamics / Wind Tunnel. Mounted
            // before the CentralPanel so egui docks it to the right.
            if !dock_owns_dockable_panels {
                crate::aero_workbench::draw_aero_workbench(self, ctx);
            }

            // FEM Workbench (right) — native linear-static + modal FEA on
            // valenx-fem's in-process solvers. A no-op unless toggled on via
            // View → FEM Workbench. Mounted before the CentralPanel so egui
            // docks it to the right (alongside the other open workbenches).
            if !dock_owns_dockable_panels {
                crate::fem_workbench::draw_fem_workbench(self, ctx);
            }

            // Induction Motor workbench (right) — 3-phase induction-motor slip /
            // power on valenx-inductionmotor. Off unless toggled via View.
            crate::inductionmotor_workbench::draw_inductionmotor_workbench(self, ctx);

            // CFD Workbench (right) — native 2-D incompressible laminar CFD
            // (SIMPLE) on valenx-cfd-native. A no-op unless toggled on via
            // View → CFD Workbench. Mounted before the CentralPanel so egui
            // docks it to the right alongside the other workbenches.
            if !dock_owns_dockable_panels {
                crate::cfd_workbench::draw_cfd_workbench(self, ctx);
            }

            // Reaction Dynamics workbench (right) — native ab-initio MD on
            // valenx-reactdyn. A no-op unless toggled on via View → Reaction
            // Dynamics. Mounted before the CentralPanel so egui docks it to
            // the right alongside the other workbenches.
            crate::reactdyn_workbench::draw_reactdyn_workbench(self, ctx);

            // Springs Workbench (right) — native helical-spring design on
            // valenx-springs. A no-op unless toggled on via View → Springs.
            crate::springs_workbench::draw_springs_workbench(self, ctx);

            // Bearing workbench (right) — rolling-bearing ISO 281 rating life on
            // valenx-bearing. Off unless toggled via View.
            crate::bearing_workbench::draw_bearing_workbench(self, ctx);

            // Belt Drive workbench (right) — flat-belt drive analysis on
            // valenx-beltdrive. Off unless toggled via View.
            crate::beltdrive_workbench::draw_beltdrive_workbench(self, ctx);

            // Buckling workbench (right) — Euler / Johnson column buckling on
            // valenx-buckling. Off unless toggled via View.
            crate::buckling_workbench::draw_buckling_workbench(self, ctx);

            // Brake workbench (right) — disc-brake torque + stop dynamics on
            // valenx-brake. Off unless toggled via View.
            crate::brake_workbench::draw_brake_workbench(self, ctx);

            // Fatigue workbench (right) — high-cycle stress-life on valenx-fatigue.
            // Off unless toggled via View.
            crate::fatigue_workbench::draw_fatigue_workbench(self, ctx);

            // Gear Tooth workbench (right) — spur-gear tooth bending strength on
            // valenx-geartooth. Off unless toggled via View.
            crate::geartooth_workbench::draw_geartooth_workbench(self, ctx);

            // Pharmacokinetics workbench (right) — one-compartment IV dosing on
            // valenx-pharmacokinetics. Off unless toggled via View.
            crate::pharmacokinetics_workbench::draw_pharmacokinetics_workbench(self, ctx);

            // Pipe Network workbench (right) — Hardy-Cross flow balancing on
            // valenx-pipenetwork. Off unless toggled via View.
            crate::pipenetwork_workbench::draw_pipenetwork_workbench(self, ctx);

            // RC Beam workbench (right) — reinforced-concrete beam flexure on
            // valenx-rcbeam. Off unless toggled via View.
            crate::rcbeam_workbench::draw_rcbeam_workbench(self, ctx);

            // Marine / Hull Workbench (right) — native box-form hull
            // hydrostatics on valenx-marine. A no-op unless toggled on.
            crate::marine_workbench::draw_marine_workbench(self, ctx);

            // Capacitor workbench (right) — parallel-plate electrostatics on
            // valenx-capacitor. Off unless toggled via View.
            crate::capacitor_workbench::draw_capacitor_workbench(self, ctx);

            // Fan Laws workbench (right) — fan / blower affinity-law scaling on
            // valenx-fanlaws. Off unless toggled via View.
            crate::fanlaws_workbench::draw_fanlaws_workbench(self, ctx);

            // Creep workbench (right) — Larson-Miller / Norton creep on
            // valenx-creep. Off unless toggled via View.
            crate::creep_workbench::draw_creep_workbench(self, ctx);

            // Electrochemistry workbench (right) — Nernst / Faraday cell analysis
            // on valenx-electrochem. Off unless toggled via View.
            crate::electrochem_workbench::draw_electrochem_workbench(self, ctx);

            // Enzyme Kinetics workbench (right) — Michaelis-Menten rate laws on
            // valenx-enzymekinetics. Off unless toggled via View.
            crate::enzymekinetics_workbench::draw_enzymekinetics_workbench(self, ctx);

            // Gears Workbench (right) — native involute-gear design on
            // valenx-gears. A no-op unless toggled on via View → Gears.
            crate::gears_workbench::draw_gears_workbench(self, ctx);

            // Pneumatics workbench (right) — cylinder force / consumption on
            // valenx-pneumatics. Off unless toggled via View.
            crate::pneumatics_workbench::draw_pneumatics_workbench(self, ctx);

            // Psychrometrics workbench (right) — moist-air state on
            // valenx-psychrometrics. Off unless toggled via View.
            crate::psychrometrics_workbench::draw_psychrometrics_workbench(self, ctx);

            // Thermistor workbench (right) — NTC R<->T on valenx-thermistor.
            // Off unless toggled via View.
            crate::thermistor_workbench::draw_thermistor_workbench(self, ctx);

            // Strain Gauge workbench (right) — Wheatstone-bridge output on
            // valenx-straingauge. Off unless toggled via View.
            crate::straingauge_workbench::draw_straingauge_workbench(self, ctx);

            // Drone Workbench (right) — native multirotor hover performance on
            // valenx-drone. A no-op unless toggled on via View → Drone / Multirotor.
            crate::drone_workbench::draw_drone_workbench(self, ctx);

            // Acoustics workbench (right) — room reverberation + SPL on
            // valenx-acoustics. Off unless toggled via View.
            crate::acoustics_workbench::draw_acoustics_workbench(self, ctx);

            // Acid-Base workbench (right) — aqueous pH / buffer equilibria on
            // valenx-acidbase. Off unless toggled via View.
            crate::acidbase_workbench::draw_acidbase_workbench(self, ctx);

            // BJT workbench (right) — bipolar-transistor DC bias on valenx-bjt.
            // Off unless toggled via View.
            crate::bjt_workbench::draw_bjt_workbench(self, ctx);

            // BMR / TDEE workbench (right) — energy expenditure on valenx-bmr.
            // Off unless toggled via View.
            crate::bmr_workbench::draw_bmr_workbench(self, ctx);

            // Bolted Joint workbench (right) — preloaded joint mechanics on
            // valenx-bolt. Off unless toggled via View.
            crate::bolt_workbench::draw_bolt_workbench(self, ctx);

            // Geomatics Workbench (right) — native geodesic calculations on
            // valenx-geomatics. A no-op unless toggled on via View → Geomatics.
            crate::geomatics_workbench::draw_geomatics_workbench(self, ctx);

            // EE / DSP workbenches (electronics batch). Each a no-op unless toggled
            // on via the View menu.
            crate::opamp_workbench::draw_opamp_workbench(self, ctx);
            crate::led_workbench::draw_led_workbench(self, ctx);
            crate::thermocouple_workbench::draw_thermocouple_workbench(self, ctx);
            crate::transmissionline_workbench::draw_transmissionline_workbench(self, ctx);
            crate::powerfactor_workbench::draw_powerfactor_workbench(self, ctx);
            crate::resistornetwork_workbench::draw_resistornetwork_workbench(self, ctx);
            crate::rectifier_workbench::draw_rectifier_workbench(self, ctx);
            crate::filter_workbench::draw_filter_workbench(self, ctx);

            // Heat Transfer workbench (right) — composite-wall 1-D heat loss on
            // valenx-heat-transfer. Off unless toggled via View.
            crate::heattransfer_workbench::draw_heattransfer_workbench(self, ctx);

            // Four-Bar Linkage Workbench (right) — native planar mechanism
            // kinematics on valenx-kinematics. A no-op unless toggled via View → Four-Bar Linkage.
            crate::fourbar_workbench::draw_fourbar_workbench(self, ctx);

            // Mechanical + civil batch workbenches. Each a no-op unless toggled via View.
            crate::shaftdesign_workbench::draw_shaftdesign_workbench(self, ctx);
            crate::screwthread_workbench::draw_screwthread_workbench(self, ctx);
            crate::pulley_workbench::draw_pulley_workbench(self, ctx);
            crate::springdesign_workbench::draw_springdesign_workbench(self, ctx);
            crate::springcombination_workbench::draw_springcombination_workbench(self, ctx);
            crate::vibration_workbench::draw_vibration_workbench(self, ctx);
            crate::rivet_workbench::draw_rivet_workbench(self, ctx);
            crate::soilbearing_workbench::draw_soilbearing_workbench(self, ctx);

            // Piping Workbench (right) — native pipe-section sizing on
            // valenx-piping. A no-op unless toggled on via View → Piping.
            crate::piping_workbench::draw_piping_workbench(self, ctx);

            // Science / bio / civil batch workbenches. Each a no-op unless toggled via View.
            crate::retainingwall_workbench::draw_retainingwall_workbench(self, ctx);
            crate::openchannel_workbench::draw_openchannel_workbench(self, ctx);
            crate::weir_workbench::draw_weir_workbench(self, ctx);
            crate::thermocycle_workbench::draw_thermocycle_workbench(self, ctx);
            crate::queueing_workbench::draw_queueing_workbench(self, ctx);
            crate::radioactivity_workbench::draw_radioactivity_workbench(self, ctx);
            crate::osmosis_workbench::draw_osmosis_workbench(self, ctx);
            crate::thermoreg_workbench::draw_thermoreg_workbench(self, ctx);
            crate::hemodynamics_workbench::draw_hemodynamics_workbench(self, ctx);
            crate::popdynamics_workbench::draw_popdynamics_workbench(self, ctx);

            // Rail / Train Workbench (right) — native train resistance + tractive
            // effort on valenx-rail. A no-op unless toggled on via View → Rail / Train.
            crate::rail_workbench::draw_rail_workbench(self, ctx);

            // Bone Mechanics workbench (right) — long-bone stress / strength on
            // valenx-bonemech. Off unless toggled via View.
            crate::bonemech_workbench::draw_bonemech_workbench(self, ctx);

            // Chain Drive workbench (right) — roller-chain kinematics on
            // valenx-chaindrive. Off unless toggled via View.
            crate::chaindrive_workbench::draw_chaindrive_workbench(self, ctx);

            // Clutch workbench (right) — friction-clutch torque capacity on
            // valenx-clutch. Off unless toggled via View.
            crate::clutch_workbench::draw_clutch_workbench(self, ctx);

            // Solenoid Coil workbench (right) — long-solenoid field / inductance
            // on valenx-coil. Off unless toggled via View.
            crate::coil_workbench::draw_coil_workbench(self, ctx);

            // Steel Column workbench (right) — Euler-Johnson axial buckling on
            // valenx-columnsteel. Off unless toggled via View.
            crate::columnsteel_workbench::draw_columnsteel_workbench(self, ctx);

            // Collision Workbench (right) — native AABB geometry + overlap tests
            // on valenx-collision. A no-op unless toggled on via View → Collision.
            crate::collision_workbench::draw_collision_workbench(self, ctx);

            // Science / engineering batch 4 (right) — each a no-op unless toggled via View.
            crate::statics_workbench::draw_statics_workbench(self, ctx);
            crate::projectile_workbench::draw_projectile_workbench(self, ctx);
            crate::conveyor_workbench::draw_conveyor_workbench(self, ctx);
            crate::fluidstatics_workbench::draw_fluidstatics_workbench(self, ctx);
            crate::plate_workbench::draw_plate_workbench(self, ctx);
            crate::strainrosette_workbench::draw_strainrosette_workbench(self, ctx);
            crate::transformer_workbench::draw_transformer_workbench(self, ctx);
            crate::threephase_workbench::draw_threephase_workbench(self, ctx);

            // Solar PV Workbench (right) — native single-diode photovoltaic cell
            // performance on valenx-solarpv. A no-op unless toggled via View → Solar PV.
            crate::solarpv_workbench::draw_solarpv_workbench(self, ctx);

            // Sheet Metal Workbench (right) — native bend allowance / deduction
            // on valenx-sheet-metal. A no-op unless toggled on via View → Sheet Metal.
            crate::sheetmetal_workbench::draw_sheetmetal_workbench(self, ctx);

            // Truss Workbench (right) — native planar pin-jointed truss analysis
            // on valenx-truss. A no-op unless toggled on via View → Truss.
            crate::truss_workbench::draw_truss_workbench(self, ctx);

            // Field Statistics Workbench (right) — descriptive statistics over a
            // pasted number list on valenx-fields. A no-op unless toggled on via
            // View → Field Statistics.
            crate::fields_workbench::draw_fields_workbench(self, ctx);

            // Gearbox workbench (right) — two-stage compound gear-train analysis
            // on valenx-gearbox. Off unless toggled via View.
            crate::gearbox_workbench::draw_gearbox_workbench(self, ctx);

            // Science batch 5 workbenches. Each a no-op unless toggled via View.
            crate::camdynamics_workbench::draw_camdynamics_workbench(self, ctx);
            crate::batteryecm_workbench::draw_batteryecm_workbench(self, ctx);
            crate::diffusion_workbench::draw_diffusion_workbench(self, ctx);
            crate::dimensional_workbench::draw_dimensional_workbench(self, ctx);
            crate::fft_workbench::draw_fft_workbench(self, ctx);

            // Fasteners Workbench (right) — ISO 4017 hex-bolt dimensions on
            // valenx-fasteners. A no-op unless toggled on via View → Fasteners.
            crate::fasteners_workbench::draw_fasteners_workbench(self, ctx);

            // Fixed-Wing / Aircraft Workbench (right) — native preliminary aircraft
            // point-performance on valenx-fixedwing. A no-op unless toggled on via
            // View → Fixed-Wing / Aircraft.
            crate::fixedwing_workbench::draw_fixedwing_workbench(self, ctx);

            // Science / engineering batch (right) — each a no-op unless toggled via View.
            crate::combustion_workbench::draw_combustion_workbench(self, ctx);
            crate::flywheel_workbench::draw_flywheel_workbench(self, ctx);
            crate::fracture_workbench::draw_fracture_workbench(self, ctx);
            crate::hydraulics_workbench::draw_hydraulics_workbench(self, ctx);
            crate::inclinedplane_workbench::draw_inclinedplane_workbench(self, ctx);
            crate::insulation_workbench::draw_insulation_workbench(self, ctx);
            crate::leadscrew_workbench::draw_leadscrew_workbench(self, ctx);
            crate::leverage_workbench::draw_leverage_workbench(self, ctx);
            crate::mohr_workbench::draw_mohr_workbench(self, ctx);
            crate::mosfet_workbench::draw_mosfet_workbench(self, ctx);
            crate::optics_workbench::draw_optics_workbench(self, ctx);
            crate::orifice_workbench::draw_orifice_workbench(self, ctx);
            crate::pressurevessel_workbench::draw_pressurevessel_workbench(self, ctx);
            crate::torsion_workbench::draw_torsion_workbench(self, ctx);
            crate::refrigeration_workbench::draw_refrigeration_workbench(self, ctx);

            // Frames Workbench (right) — structural cross-section properties on
            // valenx-frames. A no-op unless toggled on via View → Frames.
            crate::frames_workbench::draw_frames_workbench(self, ctx);

            // DC Motor Workbench (right) — native brushed-DC-motor performance on
            // valenx-dcmotor. A no-op unless toggled on via View → DC Motor.
            crate::dcmotor_workbench::draw_dcmotor_workbench(self, ctx);

            // Gas Dynamics workbench (right) — 1-D compressible-flow relations on
            // valenx-gasdynamics. A no-op unless toggled on via View → Gas Dynamics.
            crate::gasdynamics_workbench::draw_gasdynamics_workbench(self, ctx);

            // Thermal Expansion workbench (right) — linear expansion + constrained
            // stress on valenx-thermalexpansion. A no-op unless toggled via View → Thermal Expansion.
            crate::thermalexpansion_workbench::draw_thermalexpansion_workbench(self, ctx);

            // Neural-Interface workbench (right) — native BCI stimulation on
            // valenx-neuro. A no-op unless toggled on via View → Neural Interface.
            crate::neuro_workbench::draw_neuro_workbench(self, ctx);

            // Wind Turbine workbench (right) — native actuator-disc wind-turbine
            // power on valenx-windturbine. A no-op unless toggled via View → Wind Turbine.
            crate::windturbine_workbench::draw_windturbine_workbench(self, ctx);

            // Parametric-CAD workbench (right) — named parameters driving sketch
            // geometry on valenx-solvespace-3d. Off unless toggled via View.
            crate::cad_workbench::draw_cad_workbench(self, ctx);

            // Antenna workbench (right) — parabolic-dish gain / beamwidth on
            // valenx-antenna. Off unless toggled via View.
            crate::antenna_workbench::draw_antenna_workbench(self, ctx);

            // 2D Drafting workbench (right) — LibreCAD-style 2D canvas on
            // valenx-librecad-2d. Off unless toggled via View.
            crate::draft2d_workbench::draw_draft2d_workbench(self, ctx);

            // Reinforcement workbench (right) — parametric rebar cages pushed to
            // the 3D viewport, on valenx-reinforcement. Off unless toggled.
            crate::reinforcement_workbench::draw_reinforcement_workbench(self, ctx);

            // Path-Traced Render workbench (right) — global-illumination image on
            // valenx-pathtrace. Off unless toggled via View.
            crate::render_workbench::draw_render_workbench(self, ctx);

            // HVAC workbench (right) — duct sizing + pressure drop on valenx-hvac.
            // Off unless toggled via View.
            crate::hvac_workbench::draw_hvac_workbench(self, ctx);

            // Beam Workbench (right) — native Euler-Bernoulli beam bending on
            // valenx-beam. A no-op unless toggled on via View → Beam.
            crate::beam_workbench::draw_beam_workbench(self, ctx);

            // Reverse-Engineering workbench (right) — point-cloud reconstruction
            // pushed to the 3D viewport, on valenx-reverse. Off unless toggled.
            crate::reverse_workbench::draw_reverse_workbench(self, ctx);

            // Pump workbench (right) — centrifugal-pump duty point + NPSH on
            // valenx-pump. Off unless toggled via View.
            crate::pump_workbench::draw_pump_workbench(self, ctx);

            // Interior-Design workbench (right) — 2D floor plan + furniture on
            // valenx-interior. Off unless toggled via View.
            crate::interior_workbench::draw_interior_workbench(self, ctx);

            // Animation workbench (right) — joint keyframe timeline on
            // valenx-animate. Off unless toggled via View.
            crate::animate_workbench::draw_animate_workbench(self, ctx);

            // Variant-Effect workbench (right) — HGVS variant parser on
            // valenx-variant-effect. Off unless toggled via View.
            crate::variant_effect_workbench::draw_variant_effect_workbench(self, ctx);

            // Heat Pump workbench (right) — Carnot COP + second-law derating on
            // valenx-heatpump. Off unless toggled via View.
            crate::heatpump_workbench::draw_heatpump_workbench(self, ctx);

            // Astro / Launch workbench (right) — the valenx-astro launch
            // ascent simulator + mission planners. A no-op unless toggled
            // on via View → Astro / Launch. Mounted before the CentralPanel
            // so egui docks it to the right (alongside the other open
            // workbenches).
            if !dock_owns_dockable_panels {
                crate::astro_workbench::draw_astro_workbench(self, ctx);
            }

            // Pipe Flow workbench (right) — Darcy-Weisbach pipe-flow analysis
            // on valenx-pipeflow. Off unless toggled via View.
            crate::pipeflow_workbench::draw_pipeflow_workbench(self, ctx);

            // Rocket workbench (right) — the valenx-rocket-demo coupled
            // design→simulate pipeline. A no-op unless toggled on via
            // View → Rocket. Mounted before the CentralPanel so egui docks it
            // to the right alongside the other workbenches.
            if !dock_owns_dockable_panels {
                crate::rocket_workbench::draw_rocket_workbench(self, ctx);
            }

            // Battery Pack workbench (right) — series / parallel pack sizing on
            // valenx-batterypack. Off unless toggled via View.
            crate::batterypack_workbench::draw_batterypack_workbench(self, ctx);

            // Engine workbench (right) — reactive engine design → analyze →
            // optimize → export. A no-op unless toggled on via View → Engine.
            if !dock_owns_dockable_panels {
                crate::engine_workbench::draw_engine_workbench(self, ctx);
            }

            // Heat Exchanger workbench (right) — effectiveness-NTU analysis on
            // valenx-heatexchanger. Off unless toggled via View.
            crate::heatexchanger_workbench::draw_heatexchanger_workbench(self, ctx);

            // Car workbench (right) — vehicle dynamics design → simulate over
            // valenx-vehicle. A no-op unless toggled on via View → Car.
            if !dock_owns_dockable_panels {
                crate::car_workbench::draw_car_workbench(self, ctx);
            }

            // Assistant activity sidebar — a live in-app feed of what the AI
            // assistant is building. On by default; toggle via View. Renders
            // ONLY in classic (non-dock) mode: a normal dock tab hosts the
            // Assistant inside its dock tree, and an agent product tab
            // (`dock_agent_only`) has its OWN per-unit Agent chat — in both
            // dock cases the global Assistant panel would duplicate/overlap the
            // dock, so it is suppressed whenever the dock is enabled. (Unlike
            // the dockable workbenches above, `show_assistant_panel` defaults
            // on, so gating it on `!dock_owns_dockable_panels` leaked it into
            // every product tab — the layout glitch this fixes.)
            if !self.dock_enabled {
                crate::assistant_workbench::draw_assistant_workbench(self, ctx);
            }
        }

        // Opt-in dockable workbench layout: one resizable right region hosting
        // every open dock pane (workspace render + agent chat) as a draggable /
        // reorderable / splittable tile. Mounted AFTER the workbench
        // SidePanel(s) above so they reserve the right edge first and this dock
        // fills the centre — the workbench tool panel and the dock coexist
        // side-by-side. The render state + DPI are threaded in so a
        // `workspace:<n>` tile can render its live 3-D model.
        if self.dock_enabled {
            self.draw_dock_layout(ctx, wgpu_render_state.as_ref(), pixels_per_point);
        }

        // Welcome / landing-page action picked this frame. Resolved
        // inside the CentralPanel closure (the landing page renders
        // there) and dispatched against `&mut self` AFTER the closure
        // ends so we don't fight the borrow checker.
        let mut pending_landing_action: Option<landing_page::LandingAction> = None;
        // Show the welcome landing page whenever nothing is loaded —
        // no project, no STL drop, no canonical mesh. The moment any
        // of those load (e.g. user drops a .stl) we fall through to
        // the regular viewport render path. Driven purely by app
        // state — no menu / toggle — so the page exists strictly as
        // an empty-state for the central panel.
        let show_landing = self.project.is_none() && self.stl.is_none() && self.mesh.is_none();
        egui::CentralPanel::default().show(ctx, |ui| {
            // Opt-in docked / tiling layout: when on, the central area is
            // an egui_tiles tree (resizable splits, tabs, close, drag)
            // instead of the classic body. Return early so the rest of
            // the closure (landing page, 2D/3D viewport) is bypassed.
            if self.docked_layout {
                self.docking.show(ui);
                return;
            }

            // Dock-fills-workspace: when the dock is on and the 3-D viewport is
            // hidden (the Workbench+Agent launcher hides it so the units get
            // the whole screen), render the dock HERE in the CentralPanel —
            // ahead of the empty-project landing page and the viewport — with a
            // slim bar to restore the viewport or close the whole dock. Without
            // a dock tree yet, fall through to the normal central content.
            if self.dock_enabled && self.viewport_hidden && self.dock_tree.is_some() {
                let mut do_close_all = false;
                let mut do_show_viewport = false;
                ui.horizontal(|ui| {
                    ui.strong("Workbench workspace");
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if ui
                            .button("Close all")
                            .on_hover_text("Close every dock panel and Workbench+Agent unit")
                            .clicked()
                        {
                            do_close_all = true;
                        }
                        if ui
                            .button("Show 3D viewport")
                            .on_hover_text(
                                "Restore the 3-D viewport (the dock returns to the right).",
                            )
                            .clicked()
                        {
                            do_show_viewport = true;
                        }
                    });
                });
                ui.separator();
                if do_close_all {
                    self.clear_dock();
                } else {
                    if do_show_viewport {
                        self.viewport_hidden = false;
                    }
                    self.render_dock_tree_into(ui, wgpu_render_state.as_ref(), pixels_per_point);
                }
                return;
            }

            if show_landing {
                // CARGO_PKG_REPOSITORY inherits from
                // [workspace.package].repository at build time, so the
                // URL stays in lock-step with Cargo.toml without us
                // hardcoding it here. Falls back to the empty string
                // if (somehow) the manifest is missing the field; the
                // link rows still render but resolve to "/blob/…"
                // suffixes, which are visibly broken (good — that's
                // the kind of regression the build should make loud).
                pending_landing_action = landing_page::render(
                    ui,
                    env!("CARGO_PKG_VERSION"),
                    &self.settings.recent_projects,
                    env!("CARGO_PKG_REPOSITORY"),
                    self.landing_inline_message.as_deref(),
                );
                return;
            }

            // ── Central-viewport chrome header ────────────────────────────
            // A slim `−  ⋯  ✕` header above the 2D/3D body, matching the
            // right-side workbench panels' chrome. Drawn only for the live
            // viewport (not the docked / landing branches, which returned
            // above). The header mutates the global view-prefs below.
            viewport_chrome_header(self, ui);

            // Hidden → render a centered placeholder instead of the viewport
            // and skip the body entirely (no wgpu ctx / no viewport::show).
            if self.viewport_hidden {
                ui.centered_and_justified(|ui| {
                    ui.vertical_centered(|ui| {
                        ui.label(egui::RichText::new("3-D viewport hidden").weak());
                        ui.add_space(8.0);
                        if ui.button("Show viewport").clicked() {
                            self.viewport_hidden = false;
                        }
                    });
                });
                return;
            }
            // Collapsed (but not hidden) → header only; skip the body.
            if self.viewport_collapsed {
                return;
            }

            // ── Viewport dispatch ─────────────────────────────────────────
            // The 2D DNA viewport has no GPU / wgpu dependency and no
            // borrows of camera / mesh / field overlays, so it gets its
            // own clean branch. The 3D path is unchanged — all its
            // state-building code is in the `else` branch below.
            if self.active_viewport == ViewportKind::Viewport2dDna {
                // Show the user's live sequence from the Genetics
                // Sequence panel; falls back to the demo plasmid when
                // the panel is empty / unparseable.
                let live = self.genetics.sequence.to_seq_record();
                crate::viewport_2d::show(ui, &mut self.viewport_2d, live.as_ref());
                return;
            }

            // ── 3D viewport (unchanged) ───────────────────────────────────
            let wgpu_ctx = match (self.wgpu_renderer.as_mut(), wgpu_render_state.as_ref()) {
                (Some(renderer), Some(rs)) => Some(viewport::WgpuCtx {
                    renderer,
                    render_state: rs,
                    pixels_per_point,
                }),
                _ => None,
            };
            // Pick a field overlay for the viewport. The Wind Tunnel
            // workbench's aero field (Cp / velocity / pressure / Q),
            // when present, takes priority — it is paired with the
            // coloured surface / cut-plane mesh the workbench pushed
            // into `self.mesh`. Otherwise fall back to the post-run
            // results overlay below.
            //
            // Prefer the user-selected results field (clicked in the
            // Results pane); fall back to "first scalar OnNode field
            // whose data length matches the current mesh's node count"
            // so freshly-finished runs land coloured immediately. Most
            // CFD runs include `p` (pressure), most thermal runs `T`,
            // both of which satisfy the filter for a clean default.
            //
            // Within the chosen field's time series, `selected_time_index`
            // picks which snapshot to render. The index is clamped
            // here against the actual time-series length so the
            // slider can't outrun reality.
            // Resolved overlay is owned (synthesized magnitude for
            // vector fields, clone for scalars) so its data outlives
            // any borrow into self.last_run_results. Viewport receives
            // an Option<&Field> referencing this local.
            let results_field_overlay: Option<valenx_fields::Field> = self
                .last_run_results
                .as_ref()
                .zip(self.mesh.as_ref())
                .and_then(|(results, loaded)| {
                    let n = loaded.mesh.nodes.len();
                    let renderable = |f: &&valenx_fields::Field| -> bool {
                        if !matches!(f.location, valenx_fields::Location::OnNode) {
                            return false;
                        }
                        let comps = f.kind.components();
                        match f.kind {
                            valenx_fields::FieldKind::Scalar => f.data.len() == n,
                            // 1D / 2D / 3D vectors render as |v|; Voigt
                            // stress (Vector{dim:6}) does too — magnitude_field
                            // handles arbitrary dim.
                            valenx_fields::FieldKind::Vector { .. } if comps > 0 => {
                                f.data.len() == n * comps
                            }
                            // Tensors (e.g. 3x3 stress) render as
                            // Frobenius norm via magnitude_field.
                            valenx_fields::FieldKind::Tensor { .. } if comps > 0 => {
                                f.data.len() == n * comps
                            }
                            _ => false,
                        }
                    };
                    let selected = self.selected_field_name.as_deref();
                    let mut chosen_name: Option<&str> = None;
                    for name in results.fields.names() {
                        let is_selected = selected == Some(name);
                        let has_renderable = results.fields.by_name(name).any(|f| renderable(&f));
                        if has_renderable && (is_selected || chosen_name.is_none()) {
                            chosen_name = Some(name);
                            if is_selected {
                                break;
                            }
                        }
                    }
                    let name = chosen_name?;
                    let series: Vec<&valenx_fields::Field> =
                        results.fields.by_name(name).filter(renderable).collect();
                    if series.is_empty() {
                        return None;
                    }
                    let idx = self.selected_time_index.min(series.len() - 1);
                    // Project to a per-node scalar — the rendering
                    // pipeline below only knows how to colour scalars.
                    // magnitude_field is identity for already-scalar
                    // fields and ||v|| per node for vectors.
                    series[idx].magnitude_field()
                });
            // The aero flow-visualization overlay wins when active:
            // the Wind Tunnel workbench paired it with the coloured
            // body / cut-plane mesh now in `self.mesh`.
            let field_overlay_owned: Option<valenx_fields::Field> =
                self.aero_field_overlay.clone().or(results_field_overlay);
            // Sketch overlay (Phase 1H — Task 51). Forwards the active
            // sketch + selection + pending-click state only when the
            // user has the panel checkbox ticked. The rubber-band line
            // preview reads pending_click_x/y as a stand-in for the
            // real viewport cursor in world space — the numeric click
            // pipeline is the Phase 1 placeholder for a real
            // viewport-click pipeline that lands in Phase 2.
            let sketcher = &self.mesh_toolbox.sketcher;
            let pending_line_start: Option<(f64, f64)> =
                if matches!(sketcher.tool, crate::mesh_toolbox::SketcherTool::Line) {
                    sketcher
                        .pending_first_click
                        .and_then(|id| sketcher.sketch.point_at(id).ok())
                        .map(|p| p.read(&sketcher.sketch.vars))
                } else {
                    None
                };
            let cursor_screen = ctx.pointer_latest_pos();
            let sketch_overlay = if sketcher.show_overlay {
                Some(crate::sketch_overlay::SketchOverlayState {
                    sketch: &sketcher.sketch,
                    selected: &sketcher.selected,
                    pending_line_start,
                    cursor_screen,
                })
            } else {
                None
            };
            // Draft overlay (Phase 4E). Forward the document + cursor
            // + grid spacing only when the user has opted in via the
            // panel checkbox — zero cost when off.
            let draft = &self.mesh_toolbox.draft;
            let draft_overlay = if draft.show_overlay {
                Some(crate::draft_overlay::DraftOverlayState {
                    document: &draft.document,
                    selected: draft.selected_entity,
                    cursor_screen,
                    grid_spacing: draft.grid_spacing,
                })
            } else {
                None
            };
            let state = viewport::ViewportState {
                camera: &mut self.camera,
                stl: self.stl.as_ref(),
                mesh: self.mesh.as_ref(),
                shading: self.shading,
                wgpu: wgpu_ctx,
                field_overlay: field_overlay_owned.as_ref(),
                // Cut-plane preview overlay: forward the Mesh
                // Toolbox's plane only while its "Show cut overlay"
                // checkbox is ticked, so the viewport pays nothing
                // for the feature when it's off.
                cut_overlay: if self.mesh_toolbox.cut_show_overlay {
                    Some((self.mesh_toolbox.cut_point, self.mesh_toolbox.cut_normal))
                } else {
                    None
                },
                sketch_overlay,
                draft_overlay,
                // CAM simulation overlay (Phase 10E). Forward the
                // last-generated toolpath only when the user has
                // pressed "Toggle Simulate Overlay" in the CAM panel.
                cam_overlay: if self.mesh_toolbox.cam.show_overlay {
                    self.mesh_toolbox.cam.last_toolpath.as_ref()
                } else {
                    None
                },
                snap_to_grid: self.snap_to_grid,
            };
            viewport::show(ui, state);
        });

        // Dispatch the landing-page click captured above (if any).
        // Done after the closure exits so the action handlers have
        // exclusive `&mut self` — the click handlers need to call
        // `open_new_project_dialog` / `pick_project` / `load_project`,
        // all of which mutate the same fields the closure borrowed.
        if let Some(action) = pending_landing_action {
            match action {
                landing_page::LandingAction::NewProject => {
                    self.landing_inline_message = None;
                    self.open_new_project_dialog();
                }
                landing_page::LandingAction::OpenProject => {
                    self.landing_inline_message = None;
                    self.pick_project();
                }
                landing_page::LandingAction::OpenRecent(path) => {
                    self.landing_inline_message = None;
                    self.load_project(path);
                }
                landing_page::LandingAction::DropMissingRecent(path) => {
                    // The user clicked a recent-project row whose path
                    // no longer exists on disk. Prune the entry +
                    // persist + render the confirmation inline on the
                    // landing page (cleaner than routing the error to
                    // the top status bar where it doesn't make sense
                    // for an action the user took on this view).
                    let before = self.settings.recent_projects.len();
                    self.settings.recent_projects.retain(|p| p != &path);
                    if self.settings.recent_projects.len() != before {
                        crate::settings_io::save_settings_to_state_dir(&self.settings);
                    }
                    self.landing_inline_message = Some(format!(
                        "Removed missing project from recents: {}",
                        path.display()
                    ));
                }
                landing_page::LandingAction::OpenLink(target) => {
                    if let Err(e) = landing_page::open_link(&target) {
                        self.last_error = Some(format!("Open link: {e}"));
                    }
                }
            }
        }

        // Command palette overlay — split into "decide" + "invoke"
        // so the borrow checker can see that we only touch `palette`
        // while the palette is drawing, then the invoke step has
        // exclusive access to the whole app.
        //
        // `build_visible_commands` allocates ~360 `String`s per call
        // (one per static command + per-adapter dynamic entries + one
        // per workbench template + one per saved project), so running
        // it on every frame is wasteful. We memoise on
        // `(registry.len(), library.content_rev(), show_non_oss_adapters)`:
        // the registry is monotonic (adapters only ever land via
        // `register`, so a change in `len()` covers re-probes and
        // post-startup registration); `content_rev()` is a fingerprint of
        // every project's `(id, name)`, so add/remove/**rename** all flip
        // it (the old `projects.len()` key missed a rename, which keeps the
        // count identical, so the launcher showed the stale name); and
        // the OSS-filter is a bool. The cache invalidates automatically
        // when any of the three changes; a Settings toggle, a new saved
        // project, or a rename reflects on the next frame.
        //
        // Round-8: skip the per-frame Vec clone entirely when the
        // palette is closed (the common case). The clone runs only
        // when the user has the overlay up.
        if self.palette.open {
            let show_non_oss = self.settings.show_non_oss_adapters;
            let registry_len = self.registry.len();
            let projects_rev = self.library.content_rev();
            let cache_valid = matches!(
                &self.palette_cache,
                Some((len, prev, oss, _))
                    if *len == registry_len && *prev == projects_rev && *oss == show_non_oss
            );
            if !cache_valid {
                let built =
                    commands::build_visible_commands(&self.registry, show_non_oss, &self.library);
                self.palette_cache = Some((registry_len, projects_rev, show_non_oss, built));
            }
            // Take the cached Vec by clone — the dispatch step needs
            // exclusive `&mut self`, so we can't hold a borrow into
            // `self.palette_cache` across it. Cloning a `Vec<CommandKind>`
            // is a few hundred enum copies plus per-Dynamic `String`
            // clones — still vastly cheaper than rebuilding the entries
            // from the static catalogue + registry every frame.
            let visible: Vec<commands::CommandKind> = self
                .palette_cache
                .as_ref()
                .map(|(_, _, _, v)| v.clone())
                .unwrap_or_default();
            if let Some(idx) = commands::show(ctx, &mut self.palette, &visible) {
                if let Some(kind) = visible.get(idx) {
                    commands::dispatch(self, kind);
                }
            }
        }

        // Settings window — returns true on change so we know to
        // reapply the theme and optionally re-probe.
        if self.settings_open {
            let was_open = self.settings_open;
            let response = settings::show_with_response(
                ctx,
                &mut self.settings_open,
                &mut self.settings,
                &self.catalogue,
            );
            if response.changed {
                self.settings.apply_theme(ctx);
                // Persist immediately so a force-quit between change
                // and exit doesn't lose the preference.
                save_settings_to_state_dir(&self.settings);
                // Drop the cached first-run wizard report so a flip
                // of `show_non_oss_adapters` reflects the next time
                // the wizard opens. Cheap — the report is rebuilt
                // from the live registry on the next `ensure_*` call.
                self.invalidate_first_run_report();
            }
            if response.open_crashes_folder {
                let dir = crashes_dir();
                let disable = self.settings.disable_file_browser_popups;
                // Pre-create only when the kill-switch is OFF — i.e.
                // when we're about to actually open the folder in the
                // host's file browser. With the kill-switch ON we
                // surface the path as a status string only; there's
                // no value in pre-creating a directory the user
                // explicitly asked us not to "touch by default".
                if !disable {
                    let _ = std::fs::create_dir_all(&dir);
                }
                match open_path_or_copy(&dir, disable) {
                    Ok(()) => {}
                    Err(e) if e.starts_with(POPUP_DISABLED_PREFIX) => {
                        // Kill-switch path: not an error — surface as a
                        // neutral status line so the user gets the
                        // crashes-dir location without an Explorer popup.
                        self.status = Some(e);
                    }
                    Err(e) => {
                        self.last_error = Some(format!("Open crashes folder: {e}"));
                    }
                }
            }
            // Detect "just closed" → re-probe if asked.
            if was_open && !self.settings_open && self.settings.reprobe_on_close {
                self.reprobe();
                self.invalidate_first_run_report();
            }
        }

        // "File → New Project…" modal. Opened from the File menu,
        // the command palette, or Ctrl+N. Rendering mutates the
        // dialog's text-field / template-pick state in place; the
        // outcome enum tells us whether the user clicked Create
        // (load the new project), Cancel (close the dialog), or
        // neither (keep rendering next frame).
        if let Some(mut dialog) = self.new_project_dialog.take() {
            let outcome = new_project_dialog::render(ctx, &mut dialog);
            match outcome {
                new_project_dialog::DialogOutcome::Stay => {
                    self.new_project_dialog = Some(dialog);
                }
                new_project_dialog::DialogOutcome::Cancel => {
                    // Drop the dialog state — `take` already pulled
                    // it out of `self.new_project_dialog`.
                }
                new_project_dialog::DialogOutcome::Created(path) => {
                    // Project was scaffolded successfully. Mirror the
                    // existing Open-project flow so the new project
                    // shows up in the Browser tree.
                    self.load_project(path);
                    // Optional Blender-style starter cube — opt-in via
                    // Settings → New projects, so analysis projects stay
                    // empty unless the user asked for geometry to model.
                    if self.settings.starter_cube_in_new_projects {
                        self.apply_create_primitive(crate::mesh_toolbox::CadPrimitiveKind::Box);
                    }
                }
            }
        }

        // First-launch wizard — auto-opens on a fresh install,
        // re-openable from the command palette afterwards.
        // Borrow split: build the report under a separate scope so
        // the &mut self for `first_run::show` doesn't overlap with
        // the &self in `ensure_first_run_report`.
        if self.first_run_open {
            // Build / refresh the report, then clone it for the
            // render call — the report lives on `self` but
            // `show` takes a borrow + `&mut decision` which would
            // conflict with that lifetime.
            let report = self.ensure_first_run_report().clone();
            let dismissed = first_run::show(
                ctx,
                &mut self.first_run_open,
                &report,
                &mut self.first_run_decision,
                &self.catalogue,
            );
            if dismissed {
                self.first_run_open = false;
                first_run::save_first_run_to_state_dir(&self.first_run_decision);
            }
        }

        // Keyboard cheat-sheet overlay — toggled by `?` and by the
        // Help menu. Persists open state via `settings` so a user
        // who pinned it has it remembered across launches.
        if self.keyboard_help_open {
            let was_open = self.keyboard_help_open;
            keyboard_help::render_cheatsheet(ctx, &mut self.keyboard_help_open);
            if was_open != self.keyboard_help_open {
                self.settings.keyboard_shortcuts_overlay_open = self.keyboard_help_open;
                save_settings_to_state_dir(&self.settings);
            }
        }

        // Per-panel contextual help. F1 (handled by the shortcut
        // dispatcher) populates `panel_help_target` from the active
        // workbench before flipping `panel_help_open`.
        if self.panel_help_open {
            let target = self.panel_help_target.clone();
            panel_help::render_help_window(ctx, &mut self.panel_help_open, &target);
        }

        // First-launch welcome tour — independent from the
        // first-run wizard. Auto-opens on a fresh install (gated by
        // `settings.welcome_tour_completed`); re-openable from the
        // Help menu / command palette.
        if self.welcome_tour_open {
            let dismissed = welcome_tour::render(
                ctx,
                &mut self.welcome_tour_open,
                &mut self.welcome_tour_state,
            );
            if dismissed {
                self.welcome_tour_open = false;
                self.settings.welcome_tour_completed = true;
                save_settings_to_state_dir(&self.settings);
            }
        }

        // About dialog. Strings route through the locale catalogue
        // so a future translation drop changes the rendered text
        // without touching this code.
        if self.about_open {
            let cat = &self.catalogue;
            egui::Window::new(cat.lookup("dialog.about.title"))
                .open(&mut self.about_open)
                .collapsible(false)
                .resizable(false)
                .show(ctx, |ui| {
                    ui.heading(cat.lookup("dialog.about.heading"));
                    ui.label(cat.format_with(
                        "dialog.about.version",
                        &[("version", env!("CARGO_PKG_VERSION"))],
                    ));
                    ui.label(cat.lookup("dialog.about.tagline"));
                    ui.label(cat.lookup("dialog.about.licence"));
                    ui.add_space(8.0);
                    ui.label(cat.lookup("dialog.about.unifies-1"));
                    ui.label(cat.lookup("dialog.about.unifies-2"));
                    ui.label(cat.lookup("dialog.about.unifies-3"));
                    ui.add_space(8.0);
                    ui.label(cat.lookup("dialog.about.tip"));
                });
        }
    }

    /// Window-close handler. Without this, closing the window while a
    /// run (or sweep) is in-flight just drops the `JoinHandle`, leaving
    /// the subprocess (e.g. `simpleFoam.exe`) orphaned — the OS reaps
    /// it eventually but in the meantime it's still writing to the
    /// workdir, holding file handles, and consuming CPU. Cooperatively
    /// cancel both handles here and give the worker a brief grace
    /// window to honour the cancellation token before eframe tears
    /// down the GL context and exits.
    ///
    /// We don't `join()` the threads — egui's `on_exit` runs on the UI
    /// thread and a blocking join could hang the close in pathological
    /// cases. The 200ms sleep is a best-effort pause that's long
    /// enough for `LocalExecutor` to send SIGTERM / `TerminateProcess`
    /// and short enough that an impatient user clicking the X again
    /// doesn't notice the lag.
    fn on_exit(&mut self) {
        if let Some(h) = &self.run_handle {
            h.cancel();
        }
        if let Some(h) = &self.sweep_handle {
            h.cancel();
        }
        // Round-4 fix: in-flight aero solves were previously orphaned
        // on window close (the solver thread held the body mesh + the
        // ray-marcher state, so a long sweep would keep running until
        // process exit). Wire the cancellation token so on_exit asks
        // the aero worker to stop at the next per-angle checkpoint.
        if let Some(h) = &self.aero.run {
            h.cancel();
        }
        // Round-6: the RNA Designer panel can have a `BackgroundRun`
        // in flight (inverse design / LinearDesign / λ-sweep). Cancel
        // it the same way so closing the window doesn't leak a
        // worker thread that's deep in valenx-rnadesign.
        self.genetics.rna_designer.cancel_run();
        if self.run_handle.is_some() || self.sweep_handle.is_some() || self.aero.run.is_some() {
            std::thread::sleep(std::time::Duration::from_millis(200));
        }
    }
}

impl ValenxApp {
    fn draw_browser(&mut self, ui: &mut egui::Ui) {
        // VS-Code-style "Open Editors" list at the very TOP of the Browser:
        // a collapsible "Open Tabs (N)" section that mirrors EVERY open tab
        // (100+) so the user can see + switch + close them all from the
        // organized left panel rather than the cramped horizontal strip.
        // Deliberately separate from (and above) the Projects navigator:
        // Open Tabs = the open docs, Projects = the saved library. See
        // [`crate::project_tabs::draw_open_tabs_list`].
        crate::project_tabs::draw_open_tabs_list(self, ui);
        ui.separator();

        // IDE-style project library navigator below it in the Browser:
        // search + ★ Pinned / 🕘 Recent / 📁 All (folders), with
        // open-as-tab + per-project context actions. See
        // [`crate::project_navigator`].
        crate::project_navigator::draw_navigator(self, ui);
        ui.separator();

        ui.collapsing("Project", |ui| match &self.project {
            None => {
                // The welcome landing page in the central viewport
                // owns the "create / open a project" affordance now;
                // this row just acts as a quiet reminder for users
                // who scan the Browser tree before they look at the
                // centre of the window.
                ui.weak("(no project loaded — see the welcome page)");
            }
            Some(p) => {
                ui.label(&p.project.project.name);
                ui.label(format!("format: {}", p.project.project.format));
                if let Some(path) = &self.project_path {
                    ui.label(format!("path: {}", path.display()));
                }
            }
        });

        let mut pending_select: Option<String> = None;
        let mut pending_run = false;
        ui.collapsing("Cases", |ui| match &self.project {
            None => {
                ui.weak("(empty — load a project to see its cases)");
            }
            Some(p) => {
                if p.case_names().is_empty() {
                    ui.weak("(project has no cases)");
                    return;
                }
                for name in p.case_names() {
                    let is_selected = self.selected_case.as_ref() == Some(name);
                    // Show each case as a selectable row with its
                    // solver id, so users can see which adapter will
                    // handle it at a glance.
                    let case_def = p.cases.get(name);
                    let solver_hint = case_def
                        .map(|cd| cd.case.solver.as_str())
                        .unwrap_or("(no solver)");
                    let row_label = format!("{name}  ·  {solver_hint}");

                    // Resolve the adapter status so we can paint a
                    // coloured ● badge before the row. The convention
                    // is consistent with the Adapters panel under the
                    // "Advanced — external integrations" section below:
                    // green = Ready, gray = Missing, yellow = Outdated,
                    // red = Broken, purple = Disabled, dark gray =
                    // unregistered (no adapter for that solver string).
                    let badge = case_def
                        .map(|cd| adapter_id_from_solver(&cd.case.solver))
                        .and_then(|id| {
                            self.registry
                                .get(id)
                                .map(|e| (id.to_string(), Some(e.status.clone())))
                        });
                    let (badge_id, badge_status) = match badge {
                        Some((id, status)) => (id, status),
                        None => (
                            case_def
                                .map(|cd| adapter_id_from_solver(&cd.case.solver).to_string())
                                .unwrap_or_default(),
                            None,
                        ),
                    };
                    let (dot_color, status_text) = match &badge_status {
                        Some(s) => (status_color(s).0, status_label(s)),
                        None => (egui::Color32::from_rgb(80, 80, 80), "Unregistered"),
                    };
                    let tooltip = if badge_id.is_empty() {
                        "Case has no solver — adapter unknown".to_string()
                    } else {
                        format!("`{badge_id}` adapter: {status_text}")
                    };

                    // Run-history badge: ✓ if last run succeeded
                    // here, ✗ if it failed, blank otherwise. Sits
                    // between the adapter-status dot and the row
                    // label so users get an at-a-glance "have I
                    // run this case yet?" signal.
                    let history = self.run_history.get(name).cloned();
                    let (history_glyph, history_color, history_tooltip) = match &history {
                        Some(entry) if entry.succeeded => (
                            "✓",
                            egui::Color32::from_rgb(120, 220, 140),
                            format!("last run succeeded in {:?}", entry.wall_time),
                        ),
                        Some(_) => (
                            "✗",
                            egui::Color32::from_rgb(230, 120, 120),
                            "last run failed".to_string(),
                        ),
                        None => ("·", egui::Color32::from_gray(80), "not yet run".to_string()),
                    };
                    let response = ui
                        .horizontal(|ui| {
                            ui.colored_label(dot_color, "●");
                            ui.colored_label(history_color, history_glyph)
                                .on_hover_text(&history_tooltip);
                            ui.selectable_label(is_selected, row_label)
                        })
                        .inner
                        .on_hover_text(&tooltip);
                    if response.clicked() {
                        pending_select = Some(name.clone());
                    }
                    if response.double_clicked() {
                        pending_select = Some(name.clone());
                        pending_run = true;
                    }
                }
                ui.separator();
                let has_selection = self.selected_case.is_some();
                let running = self.run_handle.is_some();
                if ui
                    .add_enabled(has_selection && !running, egui::Button::new("Run selected"))
                    .clicked()
                {
                    pending_run = true;
                }
            }
        });
        if let Some(name) = pending_select {
            self.selected_case = Some(name);
        }
        if pending_run {
            self.run_selected_case();
        }

        ui.collapsing("Geometry", |ui| match &self.stl {
            None => {
                ui.weak("(empty — File → Import STL… or drop an .stl here)");
            }
            Some(stl) => {
                ui.label(
                    stl.path
                        .file_name()
                        .unwrap_or_default()
                        .to_string_lossy()
                        .as_ref(),
                );
                ui.label(format!("{} triangles", stl.mesh.triangle_count()));
            }
        });

        let mut load_mesh_requested = false;
        let mut frame_mesh_requested = false;
        ui.collapsing("Mesh", |ui| match &self.mesh {
            None => {
                ui.weak("(empty — File → Load canonical mesh…)");
                ui.weak("Run the gmsh adapter, or drop a .msh / mesh.canonical.json file.");
                if ui.small_button("Load mesh…").clicked() {
                    load_mesh_requested = true;
                }
            }
            Some(loaded) => {
                ui.label(
                    loaded
                        .path
                        .file_name()
                        .unwrap_or_default()
                        .to_string_lossy()
                        .as_ref(),
                );
                ui.label(format!(
                    "{} nodes · {} elements",
                    loaded.mesh.stats.node_count, loaded.mesh.stats.element_count
                ));
                ui.label(format!(
                    "regions: {} · boundary groups: {}",
                    loaded.mesh.stats.region_count, loaded.mesh.stats.boundary_group_count
                ));
                ui.separator();
                ui.label("Quality");
                if let Some(min) = loaded.quality.min_size {
                    ui.label(format!("min size: {min:.4e}"));
                }
                if let Some(max) = loaded.quality.max_size {
                    ui.label(format!("max size: {max:.4e}"));
                }
                if let Some(mean) = loaded.quality.mean_size {
                    ui.label(format!("mean size: {mean:.4e}"));
                }
                if let Some(ar) = loaded.quality.max_aspect_ratio {
                    ui.label(format!("max aspect: {ar:.2}"));
                }
                if let Some(sk) = loaded.quality.max_skewness {
                    ui.label(format!("max skew: {sk:.3}"));
                }
                if let Some(orth) = loaded.quality.min_orthogonality {
                    ui.label(format!("min orthogonality: {orth:.3}"));
                }
                if loaded.quality.inverted_count > 0 {
                    ui.colored_label(
                        egui::Color32::LIGHT_RED,
                        format!("{} inverted elements", loaded.quality.inverted_count),
                    );
                } else if loaded.quality.element_count > 0 {
                    ui.colored_label(egui::Color32::LIGHT_GREEN, "no inverted elements");
                }
                egui::CollapsingHeader::new("Quality histograms")
                    .id_source("mesh_quality_histograms")
                    .default_open(false)
                    .show(ui, |ui| {
                        render_aspect_histogram(ui, &loaded.aspect_hist);
                        ui.separator();
                        render_skewness_histogram(ui, &loaded.skew_hist);
                    });
                ui.horizontal(|ui| {
                    if ui.small_button("Load…").clicked() {
                        load_mesh_requested = true;
                    }
                    if ui.small_button("Frame (F)").clicked() {
                        frame_mesh_requested = true;
                    }
                });
            }
        });
        if load_mesh_requested {
            self.pick_mesh();
        }
        if frame_mesh_requested {
            self.frame_current_mesh();
        }

        let mut open_prepare_requested = false;
        let mut open_run_requested = false;
        ui.collapsing("Results", |ui| {
            // Run report from the last `Run …` invocation. Prepare-only
            // runs don't populate this — they land in the section below.
            match &self.last_run_report {
                None => {
                    ui.weak("(empty — runs populate this section)");
                }
                Some(report) => {
                    ui.label(format!("exit code: {}", report.exit_code));
                    ui.label(format!("wall time: {:?}", report.wall_time));
                    match report.converged {
                        Some(true) => ui.colored_label(egui::Color32::LIGHT_GREEN, "converged"),
                        Some(false) => {
                            ui.colored_label(egui::Color32::LIGHT_RED, "did not converge")
                        }
                        None => ui.label("convergence unknown"),
                    };
                    ui.label(format!(
                        "{} residual samples",
                        report.residual_history.len()
                    ));
                    if !report.warnings.is_empty() {
                        ui.label(format!("{} warnings", report.warnings.len()));
                    }
                }
            }

            // Field catalog from the most recent successful collect().
            // Only shown when there's something to show — empty
            // catalogs (e.g. preCICE meta-runs that don't write VTKs)
            // stay hidden so the pane doesn't grow noise.
            let mut pending_field_select: Option<String> = None;
            if let Some(results) = &self.last_run_results {
                if !results.fields.is_empty() {
                    ui.separator();
                    ui.label(format!(
                        "{} field samples in {} unique fields",
                        results.fields.len(),
                        results.fields.names().count()
                    ));
                    // Each field name is a clickable row — the
                    // currently-overlaid one is highlighted. Click
                    // toggles which field the viewport's colour
                    // overlay shows; the wireframe re-paints next
                    // frame. Non-renderable fields (vector / tensor
                    // / cell-data, etc.) still show but stay weak
                    // so users see they exist without thinking the
                    // viewport will switch to them.
                    let n_mesh_nodes = self.mesh.as_ref().map(|m| m.mesh.nodes.len()).unwrap_or(0);
                    for name in results.fields.names() {
                        let times = results.fields.time_series(name).len();
                        // Same predicate the field_overlay resolution
                        // uses upstream — accept Scalar OnNode plus
                        // any Vector / Tensor OnNode that
                        // magnitude_field can collapse into a scalar.
                        let mut renderable = false;
                        let mut kind_label = "scalar";
                        for f in results.fields.by_name(name) {
                            if !matches!(f.location, valenx_fields::Location::OnNode) {
                                continue;
                            }
                            let comps = f.kind.components();
                            if comps == 0 {
                                continue;
                            }
                            match f.kind {
                                valenx_fields::FieldKind::Scalar
                                    if f.data.len() == n_mesh_nodes =>
                                {
                                    renderable = true;
                                    kind_label = "scalar";
                                    break;
                                }
                                valenx_fields::FieldKind::Vector { .. }
                                    if f.data.len() == n_mesh_nodes * comps =>
                                {
                                    renderable = true;
                                    kind_label = "vector |v|";
                                    break;
                                }
                                valenx_fields::FieldKind::Tensor { .. }
                                    if f.data.len() == n_mesh_nodes * comps =>
                                {
                                    renderable = true;
                                    kind_label = "tensor |T|_F";
                                    break;
                                }
                                _ => {}
                            }
                        }
                        let row_label = if renderable {
                            format!("  · {name}  ({times} timesteps · {kind_label})")
                        } else {
                            format!("  · {name}  ({times} timesteps)")
                        };
                        if renderable {
                            let is_selected = self
                                .selected_field_name
                                .as_deref()
                                .map(|s| s == name)
                                .unwrap_or(false);
                            if ui.selectable_label(is_selected, row_label).clicked() {
                                pending_field_select = Some(name.to_string());
                            }
                        } else {
                            ui.weak(row_label).on_hover_text(
                                "tensor / cell-data fields and shape-mismatched \
                                 vectors aren't renderable on the wireframe yet",
                            );
                        }
                    }
                }
            }
            if let Some(name) = pending_field_select {
                // Toggle: clicking the already-selected field clears
                // the explicit pick (returns to auto-default).
                if self.selected_field_name.as_deref() == Some(name.as_str()) {
                    self.selected_field_name = None;
                } else {
                    self.selected_field_name = Some(name);
                }
                // Different field has a different time series — reset
                // the time-step slider to step 0 so we don't accidentally
                // jump to an invalid index. The clamp in the resolver
                // would catch it anyway, but resetting visibly is the
                // right UX.
                self.selected_time_index = 0;
            }

            // Time-series slider — only meaningful when the chosen
            // field has multiple snapshots (transient runs). Steady
            // runs have one snapshot and the slider would be a no-op
            // taking up space, so it's hidden in that case.
            if let Some(results) = &self.last_run_results {
                let chosen_name = self
                    .selected_field_name
                    .clone()
                    .or_else(|| results.fields.names().next().map(|s| s.to_string()));
                if let Some(name) = chosen_name {
                    let timesteps = results.fields.time_series(&name);
                    if timesteps.len() > 1 {
                        ui.separator();
                        // Clamp before reading so the slider widget
                        // never displays an out-of-range value.
                        let max_idx = timesteps.len() - 1;
                        if self.selected_time_index > max_idx {
                            self.selected_time_index = max_idx;
                        }
                        ui.label(format!(
                            "Time step {} of {} — {}",
                            self.selected_time_index + 1,
                            timesteps.len(),
                            format_time_key(timesteps[self.selected_time_index]),
                        ));
                        ui.add(
                            egui::Slider::new(&mut self.selected_time_index, 0..=max_idx)
                                .text("step"),
                        );
                    }
                }
            }

            // Run workdir — the solver's output landing zone (.vtu /
            // .frd / log.simpleFoam etc.). Shown only after at least
            // one run completes (success or failure), so users can
            // dig into post-mortem artifacts from a failed solve too.
            if let Some(workdir) = &self.last_run_workdir {
                ui.separator();
                ui.label("Run workdir:");
                let s = workdir.display().to_string();
                ui.add(
                    egui::TextEdit::singleline(&mut s.clone())
                        .desired_width(f32::INFINITY)
                        .interactive(false),
                );
                if ui
                    .button("Open in file browser")
                    .on_hover_text("Opens the run workdir in Explorer / Finder / xdg-open.")
                    .clicked()
                {
                    open_run_requested = true;
                }
            }

            // Prepare-only workdir from `Run → Prepare selected case`.
            // Surfaced separately from run results so the user can
            // distinguish "the run produced these residuals" from
            // "I asked for a dry run and these are the dicts".
            if let Some(workdir) = &self.last_prepare_workdir {
                ui.separator();
                ui.label("Prepared workdir:");
                let s = workdir.display().to_string();
                // Read-only text edit so the user can copy the path
                // out without us reaching for clipboard plumbing.
                ui.add(
                    egui::TextEdit::singleline(&mut s.clone())
                        .desired_width(f32::INFINITY)
                        .interactive(false),
                );
                if ui
                    .button("Open in file browser")
                    .on_hover_text("Opens the prepared workdir in Explorer / Finder / xdg-open.")
                    .clicked()
                {
                    open_prepare_requested = true;
                }
            }
        });
        // Defer the actual open calls so the button-click closures
        // don't fight the borrow checker (they capture &mut self
        // indirectly through the `ui` closure).
        if open_prepare_requested {
            self.open_prepare_workdir();
        }
        if open_run_requested {
            self.open_run_workdir();
        }

        // ---- Advanced (external integrations) ---------------------------
        // The Adapters tree lives here, at the bottom of the Browser,
        // default-collapsed. Valenx ships native Rust engines for every
        // major simulation domain (CFD, FEA, MD, RNA folding, quantum
        // chemistry, …), so the adapter surface — which wraps EXTERNAL
        // tools (OpenFOAM, GROMACS, Python scripts, …) — is a
        // power-user feature, not the first impression. The collapsing
        // header makes that visible: "Advanced" reads as "you almost
        // certainly don't need this", and the default-closed state
        // keeps the ~141-adapter list off-screen until someone goes
        // looking.
        let mut reprobe_requested = false;
        let show_non_oss = self.settings.show_non_oss_adapters;
        egui::CollapsingHeader::new("Advanced — external integrations")
            .id_source("browser_advanced_section")
            .default_open(false)
            .show(ui, |ui| {
                ui.weak(
                    "Optional bridges to external simulation tools. Most users \
                     don't need these — Valenx's built-in engines handle CFD, \
                     FEA, MD, and more natively.",
                );
                ui.add_space(4.0);
                let adapter_header = ui.collapsing("Adapters", |ui| {
                    ui.horizontal(|ui| {
                        ui.label(format!("{} registered", self.registry.len()));
                        if ui.small_button("Re-probe").clicked() {
                            reprobe_requested = true;
                        }
                    });
                    ui.separator();
                    if self.registry.is_empty() {
                        ui.label("(none registered)");
                    } else {
                        // Default-hide adapters wrapping academic-only /
                        // non-commercial tools. Count what we suppressed so
                        // we can surface a "(N hidden — flip Settings →
                        // Show non-OSS adapters to see them)" footer line and
                        // the user doesn't think the missing rows are a bug.
                        let mut hidden = 0usize;
                        for (_id, entry) in self.registry.iter() {
                            let info = entry.adapter.info();
                            if !show_non_oss
                                && !valenx_core::adapter_helpers::is_oss_license(info.tool_license)
                            {
                                hidden += 1;
                                continue;
                            }
                            let (dot, label_color) = status_color(&entry.status);
                            ui.horizontal(|ui| {
                                ui.colored_label(dot, "●");
                                ui.colored_label(label_color, info.display_name);
                                ui.weak(status_label(&entry.status));
                            });
                        }
                        if hidden > 0 {
                            ui.weak(format!(
                                "({hidden} hidden — flip Settings → Show non-OSS adapters \
                                 to see them)"
                            ));
                        }
                    }
                });
                adapter_header.header_response.on_hover_text(
                    "Optional bridges to external simulation tools (OpenFOAM, \
                     GROMACS, etc.). Most users don't need these — Valenx's \
                     built-in engines handle CFD, FEA, MD, and more natively.",
                );
            });
        if reprobe_requested {
            self.reprobe();
        }
    }
}

impl ValenxApp {
    /// Route a [`ShortcutAction`] decoded by [`crate::shortcuts::poll_shortcut`]
    /// to the appropriate `ValenxApp` method.
    ///
    /// The dispatcher centralises the shortcut → action mapping so
    /// every binding's behaviour is reviewable in one place, and so
    /// the command palette can call the same dispatcher to trigger
    /// keyboard-equivalent actions from the search-result list (a
    /// follow-up — palette items currently dispatch separately).
    pub(crate) fn dispatch_shortcut(&mut self, action: ShortcutAction) {
        match action {
            ShortcutAction::OpenCommandPalette => {
                self.palette.request_open();
            }
            ShortcutAction::OpenContextHelp => {
                // Resolve the active workbench's panel name so the F1
                // popup shows context-appropriate text. The priority
                // order mirrors what egui paints — aero / genetics
                // before the mesh toolbox since they sit at the
                // right and are usually what the user is staring at.
                let target = if self.show_astro_workbench {
                    "Astro / Launch".to_string()
                } else if self.show_aero_workbench {
                    "Run".to_string()
                } else if self.show_genetics_workbench {
                    self.genetics.active.label().to_string()
                } else if self.show_mesh_toolbox {
                    "Part".to_string()
                } else {
                    "Sequence".to_string()
                };
                self.panel_help_target = target;
                self.panel_help_open = true;
            }
            ShortcutAction::ToggleKeyboardHelp => {
                self.keyboard_help_open = !self.keyboard_help_open;
                self.settings.keyboard_shortcuts_overlay_open = self.keyboard_help_open;
                save_settings_to_state_dir(&self.settings);
            }
            ShortcutAction::NewProject => {
                self.open_new_project_dialog();
            }
            ShortcutAction::OpenProject => {
                self.pick_project();
            }
            ShortcutAction::ShowMeshToolbox => {
                self.show_mesh_toolbox = true;
            }
            ShortcutAction::ShowGeneticsWorkbench => {
                self.show_genetics_workbench = true;
            }
            ShortcutAction::ShowAeroWorkbench => {
                self.show_aero_workbench = true;
            }
            ShortcutAction::ShowAstroWorkbench => {
                self.show_astro_workbench = true;
            }
            ShortcutAction::RunPrimary => {
                // The active workbench owns "primary action". Priority
                // ladder: astro ascent → aero solve → genetics panel
                // run → project case. The astro / aero / genetics
                // workbenches each have a small helper that knows their
                // own primary action; if none is up we fall back to the
                // CFD/FEA case-driven Run-selected. (The astro ascent
                // runs synchronously, so there's no in-flight guard.)
                if self.show_astro_workbench {
                    crate::astro::run::run_ascent_from_shortcut(self);
                } else if self.show_aero_workbench && !self.aero.is_running() {
                    crate::aero::panels::start_run_from_shortcut(self);
                } else if self.show_genetics_workbench {
                    // Genetics panels' Run buttons are panel-internal
                    // — the simplest hook is to route Ctrl+R to the
                    // existing panel-level "run the current sub-tool"
                    // helper, which each panel exposes through a
                    // `run_active` method. For panels without an
                    // explicit Run (read-only result displays), this
                    // is a no-op.
                    crate::genetics::run_active_panel(self);
                } else if self.selected_case.is_some() {
                    self.run_selected_case();
                }
            }
            ShortcutAction::SavePrimary => {
                // The CAD-side and aero workbenches don't yet have a
                // Save concept distinct from "the panel just persists
                // its state on every edit". For the case-driven CFD/
                // FEA pipeline we re-persist settings — which is the
                // user-visible "I saved my preferences" feedback most
                // users expect from Ctrl+S.
                save_settings_to_state_dir(&self.settings);
                self.status = Some("Settings saved.".into());
            }
            ShortcutAction::Undo => {
                // The CAD-side mesh-toolbox sketcher and the
                // RNA-Designer wizard own per-panel undo stacks. They
                // expose `try_undo` methods the host calls — the
                // dispatcher tries each in priority order (active
                // workbench first) and stops at the first stack that
                // had something to undo.
                let _ = self.try_undo_in_active_panel();
            }
            ShortcutAction::Redo => {
                let _ = self.try_redo_in_active_panel();
            }
            ShortcutAction::Cancel => {
                // Esc cascade: close overlays first, then cancel any
                // in-flight long-running operation. Matches what every
                // commercial app does — a stray Esc on a focused
                // viewport shouldn't kill a 40-minute solve.
                if self.welcome_tour_open {
                    self.welcome_tour_open = false;
                    return;
                }
                if self.keyboard_help_open {
                    self.keyboard_help_open = false;
                    self.settings.keyboard_shortcuts_overlay_open = false;
                    save_settings_to_state_dir(&self.settings);
                    return;
                }
                if self.panel_help_open {
                    self.panel_help_open = false;
                    return;
                }
                if self.first_run_open {
                    self.first_run_open = false;
                    return;
                }
                if self.new_project_dialog.is_some() {
                    self.new_project_dialog = None;
                    return;
                }
                if self.run_handle.is_some() {
                    self.cancel_run();
                }
            }
        }
    }

    /// Try to undo in whatever panel currently owns the active
    /// editing surface. Returns `true` if any panel popped an undo
    /// entry, `false` if every undo stack was empty.
    ///
    /// Priority order: the panel-level edits the user is most likely
    /// to want reversed first — sequence text, alignment input, RNA
    /// designer wizard, sketcher. Read-only result panels (PCR
    /// amplicons, residual plot) don't own histories and are skipped.
    ///
    /// Coverage broadened in the second polish pass: 11 more genetics
    /// editor panels + the wind-tunnel workbench. Read-only panels
    /// (molecule viewer) stay skipped.
    pub(crate) fn try_undo_in_active_panel(&mut self) -> bool {
        // Genetics workbench panels — only meaningful when the
        // workbench is visible AND on one of the editor-style panels.
        if self.show_genetics_workbench {
            use crate::genetics_workbench::GeneticsPanel;
            match self.genetics.active {
                GeneticsPanel::Sequence => return self.genetics.sequence.undo_edit(),
                GeneticsPanel::Alignment => return self.genetics.alignment.undo_edit(),
                GeneticsPanel::RnaStructure => return self.genetics.rnastruct.undo_edit(),
                GeneticsPanel::Phylogenetics => return self.genetics.phylogenetics.undo_edit(),
                GeneticsPanel::PopulationGenetics => return self.genetics.popgen.undo_edit(),
                GeneticsPanel::MolecularDynamics => return self.genetics.md.undo_edit(),
                GeneticsPanel::Cheminformatics => return self.genetics.cheminf.undo_edit(),
                GeneticsPanel::MacromolecularStructure => {
                    return self.genetics.biostruct.undo_edit();
                }
                GeneticsPanel::QuantumChemistry => return self.genetics.qchem.undo_edit(),
                GeneticsPanel::Genomics => return self.genetics.genomics.undo_edit(),
                GeneticsPanel::SystemsBiology => return self.genetics.sysbio.undo_edit(),
                GeneticsPanel::Docking => return self.genetics.docking.undo_edit(),
                GeneticsPanel::GeneEditing => return self.genetics.genediting.undo_edit(),
                GeneticsPanel::RnaDesigner => return self.genetics.rna_designer.undo_edit(),
                GeneticsPanel::StructurePrediction => {
                    return self.genetics.structpredict.undo_edit()
                }
            }
        }
        if self.show_astro_workbench {
            return self.astro.undo_edit();
        }
        if self.show_aero_workbench {
            return self.aero.undo_edit();
        }
        // Mesh Toolbox — Sketcher + Part Design panels (polish pass 3).
        // The toolbox renders both panels inside collapsing headers in
        // one Window; try each in turn so Ctrl+Z hits whichever has the
        // most recent edit recorded.
        if self.show_mesh_toolbox {
            if self.mesh_toolbox.sketcher.undo_edit() {
                return true;
            }
            if self.mesh_toolbox.part_design.undo_edit() {
                return true;
            }
        }
        false
    }

    /// Try to redo in the active panel. Mirror of [`Self::try_undo_in_active_panel`].
    pub(crate) fn try_redo_in_active_panel(&mut self) -> bool {
        if self.show_genetics_workbench {
            use crate::genetics_workbench::GeneticsPanel;
            match self.genetics.active {
                GeneticsPanel::Sequence => return self.genetics.sequence.redo_edit(),
                GeneticsPanel::Alignment => return self.genetics.alignment.redo_edit(),
                GeneticsPanel::RnaStructure => return self.genetics.rnastruct.redo_edit(),
                GeneticsPanel::Phylogenetics => return self.genetics.phylogenetics.redo_edit(),
                GeneticsPanel::PopulationGenetics => return self.genetics.popgen.redo_edit(),
                GeneticsPanel::MolecularDynamics => return self.genetics.md.redo_edit(),
                GeneticsPanel::Cheminformatics => return self.genetics.cheminf.redo_edit(),
                GeneticsPanel::MacromolecularStructure => {
                    return self.genetics.biostruct.redo_edit();
                }
                GeneticsPanel::QuantumChemistry => return self.genetics.qchem.redo_edit(),
                GeneticsPanel::Genomics => return self.genetics.genomics.redo_edit(),
                GeneticsPanel::SystemsBiology => return self.genetics.sysbio.redo_edit(),
                GeneticsPanel::Docking => return self.genetics.docking.redo_edit(),
                GeneticsPanel::GeneEditing => return self.genetics.genediting.redo_edit(),
                GeneticsPanel::RnaDesigner => return self.genetics.rna_designer.redo_edit(),
                GeneticsPanel::StructurePrediction => {
                    return self.genetics.structpredict.redo_edit()
                }
            }
        }
        if self.show_astro_workbench {
            return self.astro.redo_edit();
        }
        if self.show_aero_workbench {
            return self.aero.redo_edit();
        }
        if self.show_mesh_toolbox {
            if self.mesh_toolbox.sketcher.redo_edit() {
                return true;
            }
            if self.mesh_toolbox.part_design.redo_edit() {
                return true;
            }
        }
        false
    }
}
