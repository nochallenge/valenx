//! Command palette — Ctrl+P overlay that lists every action the
//! app knows how to perform and runs the one the user picks.
//!
//! A command is identified by a stable `id` (for analytics + for the
//! future user-facing keybinding UI) and carries a human label plus
//! optional shortcut hint. Matching is a simple subsequence search
//! so typing "iso" matches "View: Iso" before "View: Isolate".
//!
//! The palette itself holds no app state; callers pass their
//! `ValenxApp` into `show` and the palette mutates it by invoking
//! the command closure.
//!
//! Two command tiers:
//!
//! * Static — the hand-curated [`Command`] entries returned by
//!   [`static_commands`], one per fixed verb the app supports
//!   (open project, run case, frame mesh, …).
//! * Dynamic — per-adapter run / new-case actions composed from the
//!   live `AdapterRegistry` by [`build_visible_commands`]. Lets the
//!   palette surface every registered adapter without hard-coding
//!   each one in the static slice.
//!
//! Both tiers flow through [`CommandKind`] and [`dispatch`] so the
//! palette UI doesn't need to know which tier produced an entry.

use eframe::egui;

use valenx_core::registry::AdapterRegistry;

use crate::ValenxApp;

/// Stable identifier for a command. Keep this string immutable per
/// command — changing it invalidates any recorded keybindings.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CommandId(pub &'static str);

/// One entry in the command palette — a static action with a stable
/// id, display label, optional keyboard shortcut hint, and a function
/// pointer that mutates the app when invoked.
pub struct Command {
    /// Stable identifier (see [`CommandId`]).
    pub id: CommandId,
    /// Human-readable label shown in the palette.
    pub label: &'static str,
    /// Display-only shortcut hint (the actual binding lives in the
    /// keybinding store).
    pub shortcut: Option<&'static str>,
    /// Mutator that runs when the user picks this command.
    pub invoke: fn(&mut ValenxApp),
}

/// Every static user-facing action. Order is the order they show
/// with an empty query.
pub fn static_commands() -> &'static [Command] {
    &COMMANDS
}

/// Back-compat alias for the previous public API. New callers should
/// use [`static_commands`] (static-only) or [`build_visible_commands`]
/// (static + dynamic per-adapter entries).
pub fn all_commands() -> &'static [Command] {
    static_commands()
}

/// One palette entry. Either points at a static [`Command`] or
/// describes a per-adapter dynamic action composed from the
/// `AdapterRegistry`.
///
/// `Clone` is derived so the palette-rebuild cache (see
/// `ValenxApp::palette_cache`) can hand out an owned `Vec` without
/// re-running the per-frame composition loop.
#[derive(Clone)]
pub enum CommandKind {
    Static(&'static Command),
    DynamicRunAdapter {
        id: String,
        label: String,
        adapter_id: String,
    },
    DynamicNewCase {
        id: String,
        label: String,
        adapter_id: String,
    },
}

impl CommandKind {
    /// Stable id used by analytics + the (future) keybinding UI.
    /// Static entries reuse `Command::id.0`; dynamic entries follow
    /// the conventions documented on each variant.
    pub fn id(&self) -> &str {
        match self {
            CommandKind::Static(cmd) => cmd.id.0,
            CommandKind::DynamicRunAdapter { id, .. } => id,
            CommandKind::DynamicNewCase { id, .. } => id,
        }
    }

    /// User-visible label rendered in the palette row + the future
    /// per-row shortcut indicator.
    pub fn label(&self) -> &str {
        match self {
            CommandKind::Static(cmd) => cmd.label,
            CommandKind::DynamicRunAdapter { label, .. } => label,
            CommandKind::DynamicNewCase { label, .. } => label,
        }
    }

    /// Optional shortcut hint shown in the palette gutter. Only the
    /// hand-curated static commands set this today.
    pub fn shortcut(&self) -> Option<&'static str> {
        match self {
            CommandKind::Static(cmd) => cmd.shortcut,
            CommandKind::DynamicRunAdapter { .. } => None,
            CommandKind::DynamicNewCase { .. } => None,
        }
    }
}

/// Compose the full set of palette-visible commands from the static
/// catalogue + the live adapter registry. Each registered adapter
/// contributes two dynamic entries:
///
/// * `Run: <display_name> on selected case` — runs the selected
///   case forced through this adapter (after a solver-id sanity check).
/// * `New case: <display_name>…` — opens the folder picker and
///   scaffolds a fresh project rooted at this adapter's starter
///   template.
///
/// When `show_non_oss` is `false` (the default), adapters whose
/// `tool_license` doesn't satisfy
/// [`valenx_core::adapter_helpers::is_oss_license`] are omitted from
/// the dynamic entries. The static commands always surface — the
/// filter is per-adapter, not per-command. The flag flows through
/// `Settings::show_non_oss_adapters`; callers in the GUI pass it
/// every frame so a flip in the Settings dialog takes effect on
/// the next palette open.
pub fn build_visible_commands(registry: &AdapterRegistry, show_non_oss: bool) -> Vec<CommandKind> {
    let mut cmds: Vec<CommandKind> = static_commands().iter().map(CommandKind::Static).collect();
    // Sort adapter entries by display name so the palette ordering
    // is deterministic across runs (the registry is a HashMap).
    let mut entries: Vec<_> = registry.iter().collect();
    entries.sort_by_key(|(_, e)| e.adapter.info().display_name);
    for (_, entry) in entries {
        let info = entry.adapter.info();
        if !show_non_oss && !valenx_core::adapter_helpers::is_oss_license(info.tool_license) {
            continue;
        }
        cmds.push(CommandKind::DynamicRunAdapter {
            id: format!("run.adapter.{}", info.id),
            label: format!("Run: {} on selected case", info.display_name),
            adapter_id: info.id.to_string(),
        });
        cmds.push(CommandKind::DynamicNewCase {
            id: format!("new.case.{}", info.id),
            label: format!("New case: {}…", info.display_name),
            adapter_id: info.id.to_string(),
        });
    }
    cmds
}

/// Run the action a [`CommandKind`] describes. Static entries call
/// their stored `fn` pointer; dynamic entries delegate to the
/// matching `ValenxApp` method.
pub fn dispatch(app: &mut ValenxApp, kind: &CommandKind) {
    match kind {
        CommandKind::Static(cmd) => (cmd.invoke)(app),
        CommandKind::DynamicRunAdapter { adapter_id, .. } => {
            app.run_selected_case_with_adapter(adapter_id);
        }
        CommandKind::DynamicNewCase { adapter_id, .. } => {
            app.new_case_for_adapter(adapter_id);
        }
    }
}

/// Default size for the "Audit: Show recent activity" tail action.
/// 50 entries comfortably fits the Log tab's visible window without
/// flooding the user with thousands of lines on a busy machine.
pub const AUDIT_TAIL_DEFAULT_N: usize = 50;

const COMMANDS: [Command; 65] = [
    Command {
        id: CommandId("file.new-project"),
        label: "File: New project…",
        shortcut: Some("Ctrl+N"),
        invoke: |app| app.open_new_project_dialog(),
    },
    Command {
        id: CommandId("file.open-project"),
        label: "File: Open project…",
        shortcut: Some("Ctrl+O"),
        invoke: |app| app.pick_project(),
    },
    Command {
        id: CommandId("file.import-stl"),
        label: "File: Import STL…",
        shortcut: None,
        invoke: |app| app.pick_stl(),
    },
    Command {
        id: CommandId("file.load-mesh"),
        label: "File: Load canonical mesh…",
        shortcut: None,
        invoke: |app| app.pick_mesh(),
    },
    Command {
        id: CommandId("view.frame-mesh"),
        label: "View: Frame mesh (AABB)",
        shortcut: None,
        invoke: |app| app.frame_current_mesh(),
    },
    Command {
        id: CommandId("view.front"),
        label: "View: Front",
        shortcut: Some("1"),
        invoke: |app| app.camera_mut().set_view(valenx_viz::ViewDirection::Front),
    },
    Command {
        id: CommandId("view.back"),
        label: "View: Back",
        shortcut: Some("2"),
        invoke: |app| app.camera_mut().set_view(valenx_viz::ViewDirection::Back),
    },
    Command {
        id: CommandId("view.right"),
        label: "View: Right",
        shortcut: Some("3"),
        invoke: |app| app.camera_mut().set_view(valenx_viz::ViewDirection::Right),
    },
    Command {
        id: CommandId("view.left"),
        label: "View: Left",
        shortcut: Some("4"),
        invoke: |app| app.camera_mut().set_view(valenx_viz::ViewDirection::Left),
    },
    Command {
        id: CommandId("view.top"),
        label: "View: Top",
        shortcut: Some("5"),
        invoke: |app| app.camera_mut().set_view(valenx_viz::ViewDirection::Top),
    },
    Command {
        id: CommandId("view.bottom"),
        label: "View: Bottom",
        shortcut: Some("6"),
        invoke: |app| app.camera_mut().set_view(valenx_viz::ViewDirection::Bottom),
    },
    Command {
        id: CommandId("view.iso"),
        label: "View: Iso",
        shortcut: Some("7"),
        invoke: |app| app.camera_mut().set_view(valenx_viz::ViewDirection::Iso),
    },
    Command {
        id: CommandId("view.frame"),
        label: "View: Frame geometry",
        shortcut: Some("F"),
        invoke: |app| app.frame_current_stl(),
    },
    Command {
        id: CommandId("view.toggle-shading"),
        label: "View: Toggle shaded/wireframe",
        shortcut: None,
        invoke: |app| app.toggle_shading(),
    },
    Command {
        id: CommandId("run.selected-case"),
        label: "Run: Selected case",
        shortcut: Some("F5"),
        invoke: |app| app.run_selected_case(),
    },
    Command {
        id: CommandId("run.prepare-selected-case"),
        label: "Run: Prepare selected case (no execute)",
        shortcut: None,
        invoke: |app| app.prepare_selected_case(),
    },
    Command {
        id: CommandId("run.from-prepared-workdir"),
        label: "Run: From prepared workdir (skip prepare)",
        shortcut: None,
        invoke: |app| app.run_from_prepared_workdir(),
    },
    Command {
        id: CommandId("run.sweep-selected-case"),
        label: "Run: Sweep selected case (materialise only)",
        shortcut: None,
        invoke: |app| app.sweep_selected_case(),
    },
    Command {
        id: CommandId("run.materialised-sweep-via-local-executor"),
        label: "Run: Run materialised sweep (local executor, blocking)",
        shortcut: None,
        invoke: |app| app.run_materialised_sweep_via_local_executor(),
    },
    Command {
        id: CommandId("run.materialised-sweep-async"),
        label: "Run: Run materialised sweep (async, threaded)",
        shortcut: None,
        invoke: |app| app.run_materialised_sweep_async(),
    },
    Command {
        id: CommandId("run.cancel-sweep"),
        label: "Run: Cancel running sweep",
        shortcut: None,
        invoke: |app| app.cancel_sweep(),
    },
    Command {
        id: CommandId("run.assemble-sweep-dataset"),
        label: "Run: Assemble sweep dataset (export to npy + manifest)",
        shortcut: None,
        invoke: |app| app.assemble_sweep_dataset(),
    },
    Command {
        id: CommandId("audit.verify"),
        label: "Audit: Verify audit-log chain integrity",
        shortcut: None,
        invoke: |app| app.verify_audit_log(),
    },
    Command {
        id: CommandId("audit.open"),
        label: "Audit: Open audit log in file browser",
        shortcut: None,
        invoke: |app| app.open_audit_log(),
    },
    Command {
        id: CommandId("audit.tail"),
        label: "Audit: Show recent activity (last 50) in Log tab",
        shortcut: None,
        invoke: |app| app.tail_audit_log(AUDIT_TAIL_DEFAULT_N),
    },
    Command {
        id: CommandId("residuals.export-csv"),
        label: "Residuals: Export current chart as CSV",
        shortcut: None,
        invoke: |app| app.export_residuals_csv(),
    },
    Command {
        id: CommandId("history.clear-runs"),
        label: "History: Clear all run-history badges",
        shortcut: None,
        invoke: |app| app.clear_run_history(),
    },
    Command {
        id: CommandId("history.clear-sweeps"),
        label: "History: Clear all sweep-history entries",
        shortcut: None,
        invoke: |app| app.clear_sweep_history(),
    },
    Command {
        id: CommandId("report.export-html"),
        label: "Report: Export current run as standalone HTML",
        shortcut: None,
        invoke: |app| app.export_html_report(),
    },
    Command {
        id: CommandId("view.open-prepare-workdir"),
        label: "View: Open prepared workdir in file browser",
        shortcut: None,
        invoke: |app| app.open_prepare_workdir(),
    },
    Command {
        id: CommandId("view.open-run-workdir"),
        label: "View: Open run workdir in file browser",
        shortcut: None,
        invoke: |app| app.open_run_workdir(),
    },
    Command {
        id: CommandId("run.case"),
        label: "Run: First case",
        shortcut: None,
        invoke: |app| app.run_first_case(),
    },
    Command {
        id: CommandId("run.cancel"),
        label: "Run: Cancel",
        shortcut: None,
        invoke: |app| app.cancel_run(),
    },
    Command {
        id: CommandId("settings.open"),
        label: "Settings: Preferences…",
        shortcut: None,
        invoke: |app| app.open_settings(),
    },
    Command {
        id: CommandId("settings.reprobe"),
        label: "Settings: Re-probe adapters",
        shortcut: None,
        invoke: |app| app.reprobe(),
    },
    Command {
        id: CommandId("help.about"),
        label: "Help: About Valenx",
        shortcut: None,
        invoke: |app| app.open_about(),
    },
    Command {
        id: CommandId("view.toggle-mesh-toolbox"),
        label: "View: Toggle Mesh Toolbox panel",
        shortcut: None,
        invoke: |app| app.toggle_mesh_toolbox(),
    },
    Command {
        id: CommandId("view.apply-translate-focus"),
        label: "View: Show Mesh Toolbox (focus translate inputs)",
        shortcut: None,
        invoke: |app| {
            // No-op style focus shift: re-show the toolbox so users
            // hitting Ctrl+P → "translate" land where they can type
            // deltas, even when they previously hid the panel.
            app.show_mesh_toolbox = true;
        },
    },
    // Right-side workbench visibility toggles. One entry per panel in
    // the View menu, so the palette honours its "fuzzy-search every
    // action" promise — typing a workbench name (Ctrl+P → "gears")
    // surfaces its toggle. Each flips the same `show_*` flag the menu
    // checkbox drives; the primary four also carry their Ctrl+N hint.
    Command {
        id: CommandId("view.toggle-genetics"),
        label: "View: Toggle Genetics Workbench",
        shortcut: Some("Ctrl+2"),
        invoke: |app| app.show_genetics_workbench = !app.show_genetics_workbench,
    },
    Command {
        id: CommandId("view.toggle-aero"),
        label: "View: Toggle Aerodynamics / Wind Tunnel",
        shortcut: Some("Ctrl+3"),
        invoke: |app| app.show_aero_workbench = !app.show_aero_workbench,
    },
    Command {
        id: CommandId("view.toggle-astro"),
        label: "View: Toggle Astro / Launch",
        shortcut: Some("Ctrl+4"),
        invoke: |app| app.show_astro_workbench = !app.show_astro_workbench,
    },
    Command {
        id: CommandId("view.toggle-fem"),
        label: "View: Toggle FEM Workbench",
        shortcut: None,
        invoke: |app| app.show_fem_workbench = !app.show_fem_workbench,
    },
    Command {
        id: CommandId("view.toggle-cfd"),
        label: "View: Toggle CFD Workbench",
        shortcut: None,
        invoke: |app| app.show_cfd_workbench = !app.show_cfd_workbench,
    },
    Command {
        id: CommandId("view.toggle-reactdyn"),
        label: "View: Toggle Reaction Dynamics",
        shortcut: None,
        invoke: |app| app.show_reactdyn_workbench = !app.show_reactdyn_workbench,
    },
    Command {
        id: CommandId("view.toggle-springs"),
        label: "View: Toggle Springs",
        shortcut: None,
        invoke: |app| app.show_springs_workbench = !app.show_springs_workbench,
    },
    Command {
        id: CommandId("view.toggle-gears"),
        label: "View: Toggle Gears",
        shortcut: None,
        invoke: |app| app.show_gears_workbench = !app.show_gears_workbench,
    },
    Command {
        id: CommandId("view.toggle-geomatics"),
        label: "View: Toggle Geomatics",
        shortcut: None,
        invoke: |app| app.show_geomatics_workbench = !app.show_geomatics_workbench,
    },
    Command {
        id: CommandId("view.toggle-piping"),
        label: "View: Toggle Piping",
        shortcut: None,
        invoke: |app| app.show_piping_workbench = !app.show_piping_workbench,
    },
    Command {
        id: CommandId("view.toggle-collision"),
        label: "View: Toggle Collision",
        shortcut: None,
        invoke: |app| app.show_collision_workbench = !app.show_collision_workbench,
    },
    Command {
        id: CommandId("view.toggle-sheetmetal"),
        label: "View: Toggle Sheet Metal",
        shortcut: None,
        invoke: |app| app.show_sheetmetal_workbench = !app.show_sheetmetal_workbench,
    },
    Command {
        id: CommandId("view.toggle-fields"),
        label: "View: Toggle Field Statistics",
        shortcut: None,
        invoke: |app| app.show_fields_workbench = !app.show_fields_workbench,
    },
    Command {
        id: CommandId("view.toggle-fasteners"),
        label: "View: Toggle Fasteners",
        shortcut: None,
        invoke: |app| app.show_fasteners_workbench = !app.show_fasteners_workbench,
    },
    Command {
        id: CommandId("view.toggle-frames"),
        label: "View: Toggle Frames",
        shortcut: None,
        invoke: |app| app.show_frames_workbench = !app.show_frames_workbench,
    },
    Command {
        id: CommandId("view.toggle-neuro"),
        label: "View: Toggle Neural Interface",
        shortcut: None,
        invoke: |app| app.show_neuro_workbench = !app.show_neuro_workbench,
    },
    Command {
        id: CommandId("view.toggle-cad"),
        label: "View: Toggle Parametric CAD",
        shortcut: None,
        invoke: |app| app.show_cad_workbench = !app.show_cad_workbench,
    },
    Command {
        id: CommandId("view.toggle-draft2d"),
        label: "View: Toggle 2D Drafting",
        shortcut: None,
        invoke: |app| app.show_draft2d_workbench = !app.show_draft2d_workbench,
    },
    Command {
        id: CommandId("view.toggle-reinforcement"),
        label: "View: Toggle Reinforcement",
        shortcut: None,
        invoke: |app| app.show_reinforcement_workbench = !app.show_reinforcement_workbench,
    },
    Command {
        id: CommandId("view.toggle-render"),
        label: "View: Toggle Path-Traced Render",
        shortcut: None,
        invoke: |app| app.show_render_workbench = !app.show_render_workbench,
    },
    Command {
        id: CommandId("view.toggle-hvac"),
        label: "View: Toggle HVAC",
        shortcut: None,
        invoke: |app| app.show_hvac_workbench = !app.show_hvac_workbench,
    },
    Command {
        id: CommandId("view.toggle-reverse"),
        label: "View: Toggle Reverse Engineering",
        shortcut: None,
        invoke: |app| app.show_reverse_workbench = !app.show_reverse_workbench,
    },
    Command {
        id: CommandId("view.toggle-interior"),
        label: "View: Toggle Interior Design",
        shortcut: None,
        invoke: |app| app.show_interior_workbench = !app.show_interior_workbench,
    },
    Command {
        id: CommandId("view.toggle-animate"),
        label: "View: Toggle Animation",
        shortcut: None,
        invoke: |app| app.show_animate_workbench = !app.show_animate_workbench,
    },
    Command {
        id: CommandId("view.toggle-variant-effect"),
        label: "View: Toggle Variant Effect",
        shortcut: None,
        invoke: |app| app.show_variant_effect_workbench = !app.show_variant_effect_workbench,
    },
    Command {
        id: CommandId("file.save-mesh-stl"),
        label: "File: Save mesh as STL…",
        shortcut: None,
        invoke: |app| app.save_mesh_as_stl(),
    },
    Command {
        id: CommandId("file.open-in-freecad"),
        label: "File: Open current mesh in FreeCAD",
        shortcut: None,
        invoke: |app| app.open_in_freecad(),
    },
];

/// Transient state for one open palette instance.
#[derive(Default)]
pub struct CommandPalette {
    pub open: bool,
    pub query: String,
    pub selected: usize,
}

impl CommandPalette {
    /// Open the palette, clear the query, and reset the selection
    /// cursor.
    pub fn request_open(&mut self) {
        self.open = true;
        self.query.clear();
        self.selected = 0;
    }

    /// Close the palette without clearing its state.
    pub fn close(&mut self) {
        self.open = false;
    }
}

/// Draw the palette over the supplied `commands` slice. Returns
/// `Some(idx)` when the user selected a command this frame; the
/// caller dispatches via [`dispatch(app, &commands[idx])`](dispatch).
/// Splitting the "decide what to run" step from the "run it" step
/// is what keeps the borrow checker happy: the palette mutably
/// borrows a field of `ValenxApp` while the dispatch step needs
/// `&mut ValenxApp`.
pub fn show(
    ctx: &egui::Context,
    palette: &mut CommandPalette,
    commands: &[CommandKind],
) -> Option<usize> {
    if !palette.open {
        return None;
    }

    let mut invoked: Option<usize> = None;
    let mut should_close = false;

    // Score every command by relevance to the query, drop misses,
    // sort highest-score-first so the closest match bubbles to the
    // top instead of staying in declaration order. We track the
    // ORIGINAL index in `commands` so the caller can dispatch
    // exactly the entry the user picked.
    let mut scored: Vec<(i32, usize, &CommandKind)> = commands
        .iter()
        .enumerate()
        .filter_map(|(i, c)| fuzzy_score(&palette.query, c.label()).map(|s| (s, i, c)))
        .collect();
    scored.sort_by(|a, b| b.0.cmp(&a.0).then_with(|| a.2.label().cmp(b.2.label())));
    let filtered: Vec<(usize, &CommandKind)> = scored.iter().map(|(_, i, c)| (*i, *c)).collect();

    if palette.selected >= filtered.len().max(1) {
        palette.selected = 0;
    }

    let screen = ctx.screen_rect();
    let width = (screen.width() * 0.5).clamp(320.0, 640.0);
    let palette_rect = egui::Rect::from_center_size(
        egui::pos2(screen.center().x, screen.min.y + 180.0),
        egui::vec2(width, 340.0),
    );

    egui::Area::new(egui::Id::new("valenx_command_palette"))
        .order(egui::Order::Foreground)
        .fixed_pos(palette_rect.min)
        .show(ctx, |ui| {
            egui::Frame::window(&ctx.style()).show(ui, |ui| {
                ui.set_width(width);

                // Keyboard handling
                let arrow_down = ctx.input(|i| i.key_pressed(egui::Key::ArrowDown));
                let arrow_up = ctx.input(|i| i.key_pressed(egui::Key::ArrowUp));
                let page_down = ctx.input(|i| i.key_pressed(egui::Key::PageDown));
                let page_up = ctx.input(|i| i.key_pressed(egui::Key::PageUp));
                let home = ctx.input(|i| i.key_pressed(egui::Key::Home));
                let end = ctx.input(|i| i.key_pressed(egui::Key::End));
                let enter = ctx.input(|i| i.key_pressed(egui::Key::Enter));
                let escape = ctx.input(|i| i.key_pressed(egui::Key::Escape));

                const PAGE_STEP: usize = 8;

                if arrow_down && !filtered.is_empty() {
                    palette.selected = (palette.selected + 1).min(filtered.len() - 1);
                }
                if arrow_up && palette.selected > 0 {
                    palette.selected -= 1;
                }
                if page_down && !filtered.is_empty() {
                    palette.selected = (palette.selected + PAGE_STEP).min(filtered.len() - 1);
                }
                if page_up {
                    palette.selected = palette.selected.saturating_sub(PAGE_STEP);
                }
                if home {
                    palette.selected = 0;
                }
                if end && !filtered.is_empty() {
                    palette.selected = filtered.len() - 1;
                }

                let edit = egui::TextEdit::singleline(&mut palette.query)
                    .hint_text("Search commands…")
                    .desired_width(width - 24.0)
                    .show(ui);
                edit.response.request_focus();

                ui.separator();
                egui::ScrollArea::vertical()
                    .max_height(240.0)
                    .show(ui, |ui| {
                        for (i, (orig_idx, cmd)) in filtered.iter().enumerate() {
                            let selected = i == palette.selected;
                            let response = ui.selectable_label(selected, row_text(cmd));
                            if response.clicked() {
                                invoked = Some(*orig_idx);
                                should_close = true;
                            }
                        }
                    });

                if enter {
                    if let Some((orig_idx, _)) = filtered.get(palette.selected) {
                        invoked = Some(*orig_idx);
                        should_close = true;
                    }
                }
                if escape {
                    should_close = true;
                }
            });
        });

    if should_close {
        palette.close();
    }
    invoked
}

/// Invoke a previously-selected static command on the app by id.
/// Silent no-op if the `id` doesn't match any registered command.
/// Kept for back-compat callers that still hold a [`CommandId`];
/// new code should use [`dispatch`] over a [`CommandKind`].
pub fn invoke(id: CommandId, app: &mut ValenxApp) {
    if let Some(cmd) = COMMANDS.iter().find(|c| c.id == id) {
        (cmd.invoke)(app);
    }
}

fn row_text(cmd: &CommandKind) -> String {
    match cmd.shortcut() {
        Some(s) => format!("{}        [{s}]", cmd.label()),
        None => cmd.label().to_string(),
    }
}

/// Case-insensitive subsequence match: every char of `query` appears
/// in `label` in order. Empty query matches everything. Kept as a
/// thin convenience over [`fuzzy_score`]; the runtime palette uses
/// the score directly so the ranking signal isn't discarded.
#[cfg(test)]
fn fuzzy_match(query: &str, label: &str) -> bool {
    fuzzy_score(query, label).is_some()
}

/// Score a label against a query. Higher = better. `None` means no
/// match. Used by the palette to rank results so the closest match
/// bubbles to the top instead of staying in declaration order.
///
/// Scoring tiers (additive):
/// - +1000  exact match (case-insensitive)
/// - +500   label starts with the query
/// - +200   label contains the query as a contiguous substring
/// - +100   first matched char is the start of a word (after a
///   space / `:` / `-` / `_`)
/// - +N     N consecutive matched chars (rewards "iso" matching
///   "isolate" over "i...s...o...late")
/// - -i     small penalty per char before the first match (rewards
///   early matches over deep ones)
///
/// `None` when even subsequence matching fails. Empty query
/// returns `Some(0)` (matches everything, no ranking signal).
///
/// `pub(crate)` so the project navigator ([`crate::project_navigator`])
/// can reuse the same ranking for its project-search box.
pub(crate) fn fuzzy_score(query: &str, label: &str) -> Option<i32> {
    if query.is_empty() {
        return Some(0);
    }
    let q: Vec<char> = query.chars().flat_map(|c| c.to_lowercase()).collect();
    let l: Vec<char> = label.chars().flat_map(|c| c.to_lowercase()).collect();
    if q.len() > l.len() {
        return None;
    }
    // Walk both arrays, recording matched positions so we can score
    // consecutive runs + start-of-word bonuses.
    let mut matched: Vec<usize> = Vec::with_capacity(q.len());
    let mut qi = 0;
    for (li, &lc) in l.iter().enumerate() {
        if qi < q.len() && lc == q[qi] {
            matched.push(li);
            qi += 1;
        }
        if qi == q.len() {
            break;
        }
    }
    if qi < q.len() {
        return None;
    }
    // Base subsequence-match credit.
    let mut score: i32 = 0;
    // Exact / prefix / contiguous substring tiers.
    if l == q {
        score += 1000;
    } else if l.starts_with(q.as_slice()) {
        score += 500;
    } else if l.windows(q.len()).any(|w| w == q.as_slice()) {
        score += 200;
    }
    // Start-of-word bonus on the first matched char.
    let first = matched[0];
    let prev = if first == 0 { None } else { Some(l[first - 1]) };
    let starts_word = matches!(prev, None | Some(' ' | ':' | '-' | '_' | '.'));
    if starts_word {
        score += 100;
    }
    // Reward consecutive runs of matched chars.
    let mut run: i32 = 1;
    for w in matched.windows(2) {
        if w[1] == w[0] + 1 {
            run += 1;
            score += run;
        } else {
            run = 1;
        }
    }
    // Penalise deep first-match positions.
    score -= first as i32;
    Some(score)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fuzzy_matches_subsequence() {
        assert!(fuzzy_match("", "anything"));
        assert!(fuzzy_match("iso", "View: Iso"));
        assert!(fuzzy_match("impstl", "File: Import STL…"));
        assert!(!fuzzy_match("xyz", "View: Front"));
    }

    #[test]
    fn fuzzy_is_case_insensitive() {
        assert!(fuzzy_match("ISO", "View: Iso"));
        assert!(fuzzy_match("iso", "VIEW: ISO"));
    }

    #[test]
    fn fuzzy_score_exact_match_beats_prefix_beats_substring() {
        // "iso" in three labels — exact match wins, prefix next,
        // substring third, deep subsequence last.
        let exact = fuzzy_score("iso", "iso").unwrap();
        let prefix = fuzzy_score("iso", "iso settings").unwrap();
        let substring = fuzzy_score("iso", "view: iso projection").unwrap();
        let subsequence = fuzzy_score("iso", "i s o").unwrap();
        assert!(exact > prefix, "exact {exact} should beat prefix {prefix}");
        assert!(
            prefix > substring,
            "prefix {prefix} should beat substring {substring}"
        );
        assert!(
            substring > subsequence,
            "substring {substring} should beat subseq {subsequence}"
        );
    }

    #[test]
    fn fuzzy_score_consecutive_run_beats_scattered() {
        // "ru" matches "Run: ..." (consecutive r-u) better than
        // "report: thumbnail" (r in "report", u in "thumbnail").
        let consecutive = fuzzy_score("ru", "Run: case").unwrap();
        let scattered = fuzzy_score("ru", "report: thumbnail").unwrap();
        assert!(
            consecutive > scattered,
            "consec {consecutive} vs scattered {scattered}"
        );
    }

    #[test]
    fn fuzzy_score_start_of_word_bonus() {
        // "verify" matching "Audit: Verify..." (after the colon
        // word boundary) beats matching "supervisor" (mid-word).
        let after_colon = fuzzy_score("verify", "Audit: Verify chain").unwrap();
        let mid_word = fuzzy_score("verify", "supervisor: meaning").unwrap_or(i32::MIN);
        assert!(after_colon > mid_word, "{after_colon} vs {mid_word}");
    }

    #[test]
    fn fuzzy_score_returns_none_for_no_match() {
        assert!(fuzzy_score("xyz", "View: Front").is_none());
        assert!(fuzzy_score("longer than label", "x").is_none());
    }

    #[test]
    fn fuzzy_score_empty_query_returns_zero() {
        assert_eq!(fuzzy_score("", "anything"), Some(0));
    }

    #[test]
    fn every_command_has_unique_id() {
        let mut seen = std::collections::HashSet::new();
        for cmd in all_commands() {
            assert!(seen.insert(cmd.id.0), "duplicate command id: {}", cmd.id.0);
        }
    }

    #[test]
    fn every_command_has_non_empty_label() {
        for cmd in all_commands() {
            assert!(!cmd.label.is_empty(), "empty label on {:?}", cmd.id);
        }
    }

    #[test]
    fn workbench_toggle_commands_each_flip_a_distinct_flag() {
        // The palette advertises "fuzzy-search every action", so it carries
        // one toggle per right-side workbench (mirroring the View menu).
        // Each `view.toggle-*` command must flip exactly ONE workbench
        // flag, and no two may share a flag — this catches a copy-paste
        // error pointing a command at the wrong panel (which the unique-id
        // test cannot see, since the ids stay unique) and locks the
        // coverage count against future drift.
        fn flags(a: &ValenxApp) -> Vec<bool> {
            vec![
                a.show_mesh_toolbox,
                a.show_genetics_workbench,
                a.show_aero_workbench,
                a.show_astro_workbench,
                a.show_fem_workbench,
                a.show_cfd_workbench,
                a.show_reactdyn_workbench,
                a.show_springs_workbench,
                a.show_gears_workbench,
                a.show_geomatics_workbench,
                a.show_piping_workbench,
                a.show_collision_workbench,
                a.show_sheetmetal_workbench,
                a.show_fields_workbench,
                a.show_fasteners_workbench,
                a.show_frames_workbench,
                a.show_neuro_workbench,
                a.show_cad_workbench,
                a.show_draft2d_workbench,
                a.show_reinforcement_workbench,
                a.show_render_workbench,
                a.show_hvac_workbench,
                a.show_reverse_workbench,
                a.show_interior_workbench,
                a.show_animate_workbench,
                a.show_variant_effect_workbench,
            ]
        }

        // The 25 new entries. Exclude two pre-existing `view.toggle-*`
        // commands that aren't workbench flag-flips: `mesh-toolbox` (runs
        // through a method) and `shading` (toggles render mode, not a panel).
        let toggles: Vec<&Command> = static_commands()
            .iter()
            .filter(|c| {
                c.id.0.starts_with("view.toggle-")
                    && c.id.0 != "view.toggle-mesh-toolbox"
                    && c.id.0 != "view.toggle-shading"
            })
            .collect();
        assert_eq!(toggles.len(), 25, "expected 25 workbench-toggle commands");

        // One app, all flags start false; each command flips its own flag
        // false→true exactly once. A wrong-flag closure surfaces as a
        // duplicate `touched` index.
        let mut app = ValenxApp::default();
        let mut touched = std::collections::HashSet::new();
        for cmd in &toggles {
            let before = flags(&app);
            (cmd.invoke)(&mut app);
            let after = flags(&app);
            let changed: Vec<usize> = (0..before.len())
                .filter(|&i| before[i] != after[i])
                .collect();
            assert_eq!(
                changed.len(),
                1,
                "command `{}` should flip exactly one workbench flag, flipped {changed:?}",
                cmd.id.0
            );
            assert!(
                touched.insert(changed[0]),
                "command `{}` flips a flag already owned by another command",
                cmd.id.0
            );
        }
    }

    // ---- Phase 44: dynamic per-adapter commands ----

    #[test]
    fn build_visible_commands_empty_registry_returns_static_commands() {
        let registry = AdapterRegistry::new();
        let cmds = build_visible_commands(&registry, false);
        // Empty registry → exactly the static commands surface.
        assert_eq!(cmds.len(), static_commands().len());
        // And every entry is the Static variant.
        for cmd in &cmds {
            assert!(
                matches!(cmd, CommandKind::Static(_)),
                "expected only Static; got id `{}`",
                cmd.id()
            );
        }
    }

    #[test]
    fn build_visible_commands_populated_registry_includes_per_adapter_entries() {
        // We can't build a real adapter inside a unit test without
        // pulling in a heavy dep; instead use one of the workspace
        // adapters via the existing registration path. The OpenFOAM
        // adapter is plain (no probe side effects) so it's a good
        // smoke choice.
        use std::sync::Arc;
        let mut registry = AdapterRegistry::new();
        registry.register(Arc::new(valenx_adapter_openfoam::OpenFoamAdapter::new()));
        let cmds = build_visible_commands(&registry, false);
        // 2 dynamic entries per registered adapter (run + new-case).
        // OpenFOAM is GPL-3.0-only (OSS) so the OSS-only filter
        // keeps it visible.
        assert_eq!(cmds.len(), static_commands().len() + 2);
        // The labels follow the documented pattern.
        let labels: Vec<&str> = cmds.iter().map(|c| c.label()).collect();
        assert!(
            labels
                .iter()
                .any(|l| l.starts_with("Run: OpenFOAM on selected case")),
            "missing Run entry; labels: {labels:?}"
        );
        assert!(
            labels.iter().any(|l| l.starts_with("New case: OpenFOAM")),
            "missing New case entry; labels: {labels:?}"
        );
        // Stable ids follow the documented `run.adapter.<id>` /
        // `new.case.<id>` shape.
        let ids: Vec<&str> = cmds.iter().map(|c| c.id()).collect();
        assert!(ids.contains(&"run.adapter.openfoam"), "got: {ids:?}");
        assert!(ids.contains(&"new.case.openfoam"), "got: {ids:?}");
    }

    #[test]
    fn build_visible_commands_hides_non_oss_adapter_by_default() {
        // Rosetta's tool_license is "Rosetta-License" (academic /
        // non-commercial). The OSS-only filter must suppress its
        // two dynamic palette entries.
        use std::sync::Arc;
        let mut registry = AdapterRegistry::new();
        registry.register(Arc::new(valenx_adapter_rosetta::RosettaAdapter::new()));
        let cmds = build_visible_commands(&registry, false);
        assert_eq!(
            cmds.len(),
            static_commands().len(),
            "non-OSS adapter must contribute zero dynamic entries when show_non_oss=false"
        );
        // Flipping the flag surfaces both dynamic entries again.
        let cmds_all = build_visible_commands(&registry, true);
        assert_eq!(cmds_all.len(), static_commands().len() + 2);
    }
}
