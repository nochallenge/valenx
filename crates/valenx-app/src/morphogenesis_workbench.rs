//! The right-side **Morphogenesis** workbench — a native, in-house **Turing
//! reaction–diffusion** sandbox that grows organic biological patterns
//! (spots / coral / mitosis / mazes) **live** on a 3-D surface.
//!
//! It is a front-end over [`crate::morphogenesis::Morphogenesis`], the pure
//! Gray–Scott reaction–diffusion core (Turing 1952). The `V` morphogen field is
//! rendered as a **shaded 3-D height-field**: each grid cell becomes a quad
//! lifted by its `V` value, lit with a simple Lambert term and tinted by a
//! viridis-style colormap, drawn in back-to-front isometric order so the user
//! sees a living surface. While playing, the workbench steps the simulation
//! every frame and requests a repaint, so the pattern visibly emerges in real
//! time.
//!
//! Mirrors the other real-time workbenches (`mission_planner_workbench`): a
//! [`crate::workbench_chrome::workbench_shell`] panel gated on
//! [`crate::ValenxApp::show_morphogenesis_workbench`], toggled from the View
//! menu and openable by the agent bridge under the workbench id
//! `"morphogenesis"`. The agent bridge can set the controls
//! (`agent_set` / `agent_control_names`), read a status line (`agent_readout`),
//! and drive playback via the RunCommand ids `morphogenesis.play` /
//! `morphogenesis.pause`.

use eframe::egui;

use crate::morphogenesis::Morphogenesis;
use crate::ValenxApp;

/// The grid edge (cells) the workbench seeds a new field with by default.
const DEFAULT_GRID: u32 = 96;
/// Default explicit-Euler steps advanced per animation frame.
const DEFAULT_STEPS_PER_FRAME: u32 = 8;

// ---------------------------------------------------------------------------
// Workbench state
// ---------------------------------------------------------------------------

/// Persistent state for the Morphogenesis workbench: the live reaction–diffusion
/// field plus the real-time playback controls and the selected preset.
pub struct MorphogenesisWorkbenchState {
    /// The live Gray–Scott field (`U`/`V` morphogens on a toroidal grid).
    pub field: Morphogenesis,
    /// Edge length (cells) of the (square) grid; changing it resizes + reseeds.
    pub grid: u32,
    /// Explicit-Euler steps advanced per animation frame while playing. Higher
    /// = the pattern grows faster (more compute per frame).
    pub steps_per_frame: u32,
    /// Name of the active Gray–Scott preset (`spots` / `coral` / `mitosis` /
    /// `maze`); drives `feed`/`kill`.
    pub preset: String,
    /// Whether real-time playback is running. While `true` the workbench steps
    /// the field every frame and requests a repaint so the pattern animates.
    pub playing: bool,
}

impl Default for MorphogenesisWorkbenchState {
    fn default() -> Self {
        Self {
            field: Morphogenesis::new(DEFAULT_GRID as usize, DEFAULT_GRID as usize),
            grid: DEFAULT_GRID,
            steps_per_frame: DEFAULT_STEPS_PER_FRAME,
            preset: "mitosis".to_string(),
            playing: false,
        }
    }
}

impl MorphogenesisWorkbenchState {
    /// Re-seed the field (restart the pattern from the centred germ), keeping
    /// the current grid size, preset, and reaction parameters.
    fn reseed(&mut self) {
        self.field.reseed();
    }

    /// Resize the field to a `grid × grid` square and reseed.
    fn resize(&mut self) {
        let g = self.grid.clamp(16, 512) as usize;
        self.field.resize(g, g);
    }

    /// Apply the named preset (sets `feed`/`kill` and reseeds).
    fn apply_preset(&mut self, name: &str) {
        self.preset = name.to_string();
        self.field.apply_preset(name);
    }

    /// The user-visible captions of every control the agent bridge can set via
    /// `SetControl` (returned by `ListControls`). Order follows the form.
    pub fn agent_control_names() -> &'static [&'static str] {
        &[
            "Feed rate F",
            "Kill rate k",
            "Grid size",
            "Sim speed (steps/frame)",
            "Preset",
        ]
    }

    /// Set one labelled control by its caption, for the agent `SetControl`
    /// bridge. Fail-loud on an unknown caption / wrong type / out-of-range; no
    /// state is written on error and nothing panics. Numeric controls read
    /// [`crate::agent_commands::AgentValue`] as `f64`/`i64`; `Preset` reads a
    /// string and must name a known régime.
    pub fn agent_set(
        &mut self,
        name: &str,
        value: &crate::agent_commands::AgentValue,
    ) -> Result<(), String> {
        match name {
            "Feed rate F" => {
                let v = value.as_f64()? as f32;
                if !(0.0..=0.12).contains(&v) {
                    return Err(format!("Feed rate F must be in 0..=0.12, got {v}"));
                }
                self.field.feed = v;
            }
            "Kill rate k" => {
                let v = value.as_f64()? as f32;
                if !(0.0..=0.10).contains(&v) {
                    return Err(format!("Kill rate k must be in 0..=0.10, got {v}"));
                }
                self.field.kill = v;
            }
            "Grid size" => {
                let n = value.as_i64()?;
                if !(16..=512).contains(&n) {
                    return Err(format!("Grid size must be in 16..=512, got {n}"));
                }
                self.grid = n as u32;
                self.resize();
            }
            "Sim speed (steps/frame)" => {
                let n = value.as_i64()?;
                if !(0..=200).contains(&n) {
                    return Err(format!(
                        "Sim speed (steps/frame) must be in 0..=200, got {n}"
                    ));
                }
                self.steps_per_frame = n as u32;
            }
            "Preset" => {
                let s = value.as_str()?;
                if !matches!(s, "spots" | "coral" | "mitosis" | "maze") {
                    return Err(format!(
                        "Preset must be one of spots/coral/mitosis/maze, got {s:?}"
                    ));
                }
                self.apply_preset(s);
            }
            other => return Err(format!("unknown morphogenesis control: {other:?}")),
        }
        Ok(())
    }

    /// The current readout text for the agent `ReadReadout` bridge: the preset,
    /// step count, grid size, and mean `V`. Always `Some` (the field always
    /// exists).
    pub fn agent_readout(&self) -> Option<String> {
        Some(format!(
            "Morphogenesis \u{00B7} preset {} (F={:.4}, k={:.4}) \u{00B7} {}x{} grid \u{00B7} \
             {} steps \u{00B7} mean V {:.4} \u{00B7} {}",
            self.preset,
            self.field.feed,
            self.field.kill,
            self.field.w,
            self.field.h,
            self.field.steps,
            self.field.mean_v(),
            if self.playing { "playing" } else { "paused" },
        ))
    }
}

// ---------------------------------------------------------------------------
// Bridge run actions (play / pause)
// ---------------------------------------------------------------------------

/// Start real-time playback (the in-panel Play action). Factored out so the
/// button and the `morphogenesis.play` bridge id share one path.
pub(crate) fn play(app: &mut ValenxApp) {
    app.morphogenesis.playing = true;
}

/// Pause real-time playback (the in-panel Pause action). Factored out so the
/// button and the `morphogenesis.pause` bridge id share one path.
pub(crate) fn pause(app: &mut ValenxApp) {
    app.morphogenesis.playing = false;
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Draw the Morphogenesis workbench. A no-op unless toggled on via
/// View → Morphogenesis.
pub fn draw_morphogenesis_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_morphogenesis_workbench {
        return;
    }
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_morphogenesis_workbench",
        "Morphogenesis (Turing reaction–diffusion, live 3-D)",
        morphogenesis_workbench_body,
    );
    if close {
        app.show_morphogenesis_workbench = false;
    }
}

// ---------------------------------------------------------------------------
// Workbench body
// ---------------------------------------------------------------------------

fn morphogenesis_workbench_body(app: &mut ValenxApp, ui: &mut egui::Ui) {
    egui::ScrollArea::vertical()
        .auto_shrink([false, false])
        .show(ui, |ui| morphogenesis_workbench_body_inner(app, ui));
}

fn morphogenesis_workbench_body_inner(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.label(
        egui::RichText::new(
            "In-house Turing morphogenesis \u{00B7} Gray\u{2013}Scott reaction\u{2013}diffusion \
             [two morphogens U,V on a toroidal grid grow organic patterns \u{2014} spots / coral \
             / mitosis / mazes \u{2014} rendered as a live 3-D surface]",
        )
        .weak()
        .small(),
    );
    ui.separator();

    let s = &mut app.morphogenesis;

    // --- Preset selector ----------------------------------------------------
    ui.horizontal(|ui| {
        let lbl = ui.label("Preset");
        let mut chosen: Option<&'static str> = None;
        egui::ComboBox::from_id_source("morphogenesis_preset")
            .selected_text(s.preset.clone())
            .show_ui(ui, |ui| {
                for name in ["spots", "coral", "mitosis", "maze"] {
                    if ui
                        .selectable_label(s.preset == name, name)
                        .on_hover_text(preset_hint(name))
                        .clicked()
                    {
                        chosen = Some(name);
                    }
                }
            })
            .response
            .labelled_by(lbl.id);
        if let Some(name) = chosen {
            s.apply_preset(name);
        }
    });

    // --- Numeric controls ---------------------------------------------------
    egui::Grid::new("morphogenesis_controls")
        .num_columns(2)
        .striped(true)
        .show(ui, |ui| {
            // Feed rate F.
            let lbl = ui.label("Feed rate F");
            ui.add(
                egui::DragValue::new(&mut s.field.feed)
                    .speed(0.001)
                    .range(0.0..=0.12)
                    .max_decimals(4),
            )
            .labelled_by(lbl.id)
            .on_hover_text("Gray\u{2013}Scott feed rate (replenishes U). Sets the pattern régime.");
            ui.end_row();

            // Kill rate k.
            let lbl = ui.label("Kill rate k");
            ui.add(
                egui::DragValue::new(&mut s.field.kill)
                    .speed(0.001)
                    .range(0.0..=0.10)
                    .max_decimals(4),
            )
            .labelled_by(lbl.id)
            .on_hover_text("Gray\u{2013}Scott kill rate (removes V). Sets the pattern régime.");
            ui.end_row();

            // Grid size (square edge; resizes + reseeds).
            let lbl = ui.label("Grid size");
            let resp = ui
                .add(egui::DragValue::new(&mut s.grid).speed(1).range(16..=512))
                .labelled_by(lbl.id)
                .on_hover_text("Grid edge in cells; changing it resizes the field and reseeds.");
            if resp.changed() {
                s.resize();
            }
            ui.end_row();

            // Sim speed (steps per frame).
            let lbl = ui.label("Sim speed (steps/frame)");
            ui.add(
                egui::DragValue::new(&mut s.steps_per_frame)
                    .speed(0.5)
                    .range(0..=200),
            )
            .labelled_by(lbl.id)
            .on_hover_text(
                "Explicit-Euler steps advanced each frame while playing. Higher = faster growth.",
            );
            ui.end_row();
        });

    // --- Play / Pause / Reseed ----------------------------------------------
    ui.add_space(4.0);
    ui.horizontal(|ui| {
        let toggle_label = if s.playing {
            "\u{23F8} Pause"
        } else {
            "\u{25B6} Play"
        };
        if ui
            .button(toggle_label)
            .on_hover_text("Start / pause real-time growth of the reaction\u{2013}diffusion field.")
            .clicked()
        {
            s.playing = !s.playing;
        }
        if ui
            .button("\u{21BA} Reseed")
            .on_hover_text("Restart the pattern from the centred seed (keeps grid + parameters).")
            .clicked()
        {
            s.reseed();
        }
    });

    // --- Real-time step -----------------------------------------------------
    // While playing, advance the field by `steps_per_frame` and request a
    // repaint so the next frame keeps animating the pattern's growth.
    if s.playing && s.steps_per_frame > 0 {
        s.field.step(s.steps_per_frame);
        ui.ctx().request_repaint();
    }

    // --- Status readout -----------------------------------------------------
    ui.add_space(6.0);
    ui.label(
        egui::RichText::new(format!(
            "Steps: {}   \u{00B7}   {}x{} grid   \u{00B7}   mean V: {:.4}",
            s.field.steps,
            s.field.w,
            s.field.h,
            s.field.mean_v(),
        ))
        .strong()
        .color(egui::Color32::from_rgb(150, 210, 160)),
    );

    // --- Live 3-D surface ---------------------------------------------------
    ui.add_space(6.0);
    ui.separator();
    draw_surface_3d(&s.field, ui);
}

/// One-line description of what each preset grows into (combo tooltips).
fn preset_hint(name: &str) -> &'static str {
    match name {
        "spots" => "Isolated dots (leopard-like spots).",
        "coral" => "Branching coral / fingerprint ridges.",
        "mitosis" => "Self-replicating cell-like blobs.",
        "maze" => "Labyrinthine stripes.",
        _ => "",
    }
}

// ---------------------------------------------------------------------------
// 3-D shaded height-field rendering (egui painter, isometric)
// ---------------------------------------------------------------------------

/// Render the `V` morphogen field as a **shaded 3-D height-field**.
///
/// Each grid cell is projected isometrically with its height proportional to
/// `V`, drawn as a filled quad whose colour comes from a viridis-style colormap
/// (by normalised height) modulated by a simple Lambert shading term from a
/// fixed light. Quads are painted **back-to-front** (far rows first) so nearer
/// ridges correctly occlude farther ones, giving a solid living-surface look.
///
/// To keep the per-frame cost bounded on large grids, the field is sampled on a
/// capped lattice (`MAX_CELLS` per side) rather than one quad per cell.
fn draw_surface_3d(field: &Morphogenesis, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Live 3-D surface (V morphogen height)").strong());

    let available = ui.available_size();
    let (rect, _) = ui.allocate_exact_size(
        egui::vec2(available.x.clamp(240.0, 640.0), 360.0),
        egui::Sense::hover(),
    );
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 2.0, egui::Color32::from_rgb(8, 10, 16));

    if field.w == 0 || field.h == 0 {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "empty field",
            egui::FontId::monospace(12.0),
            egui::Color32::from_gray(120),
        );
        return;
    }

    // Sampling lattice: at most MAX_CELLS samples per axis.
    const MAX_CELLS: usize = 110;
    let nx = field.w.min(MAX_CELLS);
    let ny = field.h.min(MAX_CELLS);

    // Normalisation range for height + colour.
    let (lo, hi) = field.field_minmax();
    let span = (hi - lo).max(1e-4);

    // Sample the field onto an `ny × nx` height grid (normalised to [0,1]).
    let mut heights = vec![0.0f32; nx * ny];
    for gy in 0..ny {
        let sy = gy * (field.h - 1) / ny.max(2).saturating_sub(1).max(1);
        for gx in 0..nx {
            let sx = gx * (field.w - 1) / nx.max(2).saturating_sub(1).max(1);
            let v = field.v_at(sx, sy);
            heights[gy * nx + gx] = ((v - lo) / span).clamp(0.0, 1.0);
        }
    }

    // --- Isometric projection setup -----------------------------------------
    // Grid coordinates (gx in 0..nx, gy in 0..ny) map to screen via an oblique
    // projection; height lifts points up (−y on screen).
    let inner = rect.shrink(20.0);
    let cell = (inner.width() / (nx as f32 + ny as f32 * 0.5)).max(1.0);
    // Vertical exaggeration so the surface reads as 3-D.
    let height_px = cell * (nx.min(ny) as f32) * 0.22;
    // Origin: near the top so the whole skewed grid fits.
    let ox = inner.left();
    let oy = inner.top() + ny as f32 * cell * 0.5 + height_px * 0.5;

    let project = |gx: f32, gy: f32, z: f32| -> egui::Pos2 {
        // Oblique: x shifts right with gx and (half) with gy; y goes down with
        // gy and up with height.
        let px = ox + gx * cell + gy * cell * 0.5;
        let py = oy + gy * cell * 0.5 - z * height_px;
        egui::pos2(px, py)
    };

    // --- Lambert shading from a fixed light ---------------------------------
    // Surface-normal approximation per quad via its height gradient.
    let light = normalize3([-0.4, -0.5, 0.76]);

    // Paint quads back-to-front (far rows = high gy first) so near ridges
    // occlude far ones.
    for gy in (0..ny.saturating_sub(1)).rev() {
        for gx in 0..nx.saturating_sub(1) {
            let h00 = heights[gy * nx + gx];
            let h10 = heights[gy * nx + gx + 1];
            let h01 = heights[(gy + 1) * nx + gx];
            let h11 = heights[(gy + 1) * nx + gx + 1];
            let havg = 0.25 * (h00 + h10 + h01 + h11);

            // Quad corners in screen space.
            let p00 = project(gx as f32, gy as f32, h00);
            let p10 = project(gx as f32 + 1.0, gy as f32, h10);
            let p11 = project(gx as f32 + 1.0, gy as f32 + 1.0, h11);
            let p01 = project(gx as f32, gy as f32 + 1.0, h01);

            // Approximate normal from finite differences of the height field.
            let dzdx = (h10 - h00) + (h11 - h01);
            let dzdy = (h01 - h00) + (h11 - h10);
            let normal = normalize3([-dzdx, -dzdy, 1.0]);
            let lambert = dot3(normal, light).clamp(0.0, 1.0);
            // Ambient floor so valleys aren't pure black.
            let shade = 0.35 + 0.65 * lambert;

            let base = viridis(havg);
            let col = egui::Color32::from_rgb(
                (base.0 as f32 * shade) as u8,
                (base.1 as f32 * shade) as u8,
                (base.2 as f32 * shade) as u8,
            );
            painter.add(egui::Shape::convex_polygon(
                vec![p00, p10, p11, p01],
                col,
                egui::Stroke::NONE,
            ));
        }
    }

    // Frame + caption.
    painter.rect_stroke(
        rect,
        2.0,
        egui::Stroke::new(1.0, egui::Color32::from_rgb(50, 60, 80)),
    );
    painter.text(
        egui::pos2(rect.right() - 5.0, rect.bottom() - 8.0),
        egui::Align2::RIGHT_BOTTOM,
        format!("{nx}x{ny} samples \u{00B7} viridis(V) + Lambert"),
        egui::FontId::monospace(9.0),
        egui::Color32::from_gray(140),
    );
}

/// A compact viridis-style colormap: `t` in `[0,1]` → `(r,g,b)` going
/// dark-purple → blue → teal → green → yellow. Piecewise-linear over a handful
/// of anchors (no external colormap dep).
fn viridis(t: f32) -> (u8, u8, u8) {
    const STOPS: [(f32, f32, f32); 6] = [
        (68.0, 1.0, 84.0),    // dark purple
        (59.0, 82.0, 139.0),  // blue
        (33.0, 145.0, 140.0), // teal
        (94.0, 201.0, 98.0),  // green
        (253.0, 231.0, 37.0), // yellow
        (253.0, 231.0, 37.0), // clamp tail
    ];
    let t = t.clamp(0.0, 1.0) * (STOPS.len() as f32 - 2.0);
    let i = t.floor() as usize;
    let f = t - i as f32;
    let a = STOPS[i];
    let b = STOPS[i + 1];
    (
        (a.0 + (b.0 - a.0) * f) as u8,
        (a.1 + (b.1 - a.1) * f) as u8,
        (a.2 + (b.2 - a.2) * f) as u8,
    )
}

fn normalize3(v: [f32; 3]) -> [f32; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len <= 1e-9 {
        [0.0, 0.0, 1.0]
    } else {
        [v[0] / len, v[1] / len, v[2] / len]
    }
}

fn dot3(a: [f32; 3], b: [f32; 3]) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

// ---------------------------------------------------------------------------
// Tests (unit)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_seeds_a_field() {
        let s = MorphogenesisWorkbenchState::default();
        assert_eq!(s.grid, DEFAULT_GRID);
        assert_eq!(s.field.w, DEFAULT_GRID as usize);
        assert_eq!(s.preset, "mitosis");
        assert!(!s.playing);
    }

    #[test]
    fn agent_set_feed_and_kill() {
        use crate::agent_commands::AgentValue;
        let mut s = MorphogenesisWorkbenchState::default();
        s.agent_set("Feed rate F", &AgentValue::Float(0.04))
            .expect("valid feed");
        assert!((s.field.feed - 0.04).abs() < 1e-6);
        s.agent_set("Kill rate k", &AgentValue::Float(0.06))
            .expect("valid kill");
        assert!((s.field.kill - 0.06).abs() < 1e-6);
    }

    #[test]
    fn agent_set_grid_resizes() {
        use crate::agent_commands::AgentValue;
        let mut s = MorphogenesisWorkbenchState::default();
        s.agent_set("Grid size", &AgentValue::Int(48))
            .expect("valid grid");
        assert_eq!(s.grid, 48);
        assert_eq!(s.field.w, 48);
        assert_eq!(s.field.h, 48);
    }

    #[test]
    fn agent_set_preset_changes_params() {
        use crate::agent_commands::AgentValue;
        let mut s = MorphogenesisWorkbenchState::default();
        s.agent_set("Preset", &AgentValue::Str("maze".to_string()))
            .expect("valid preset");
        assert_eq!(s.preset, "maze");
        assert_eq!((s.field.feed, s.field.kill), (0.029, 0.057));
    }

    #[test]
    fn agent_set_rejects_out_of_range_and_unknown() {
        use crate::agent_commands::AgentValue;
        let mut s = MorphogenesisWorkbenchState::default();
        assert!(s.agent_set("Feed rate F", &AgentValue::Float(9.9)).is_err());
        assert!(s.agent_set("Grid size", &AgentValue::Int(2)).is_err());
        assert!(s.agent_set("Grid size", &AgentValue::Int(9999)).is_err());
        assert!(s
            .agent_set("Preset", &AgentValue::Str("nope".to_string()))
            .is_err());
        assert!(s.agent_set("bogus", &AgentValue::Float(1.0)).is_err());
    }

    #[test]
    fn agent_readout_reports_preset_and_steps() {
        let mut s = MorphogenesisWorkbenchState::default();
        s.field.step(10);
        let r = s.agent_readout().expect("readout always present");
        assert!(r.contains("Morphogenesis"), "readout: {r}");
        assert!(r.contains("preset"), "readout should name the preset: {r}");
        assert!(r.contains("mean V"), "readout should report mean V: {r}");
    }

    #[test]
    fn play_pause_helpers_toggle_flag() {
        let mut app = ValenxApp::default();
        assert!(!app.morphogenesis.playing);
        play(&mut app);
        assert!(app.morphogenesis.playing);
        pause(&mut app);
        assert!(!app.morphogenesis.playing);
    }

    #[test]
    fn control_names_are_listed() {
        let names = MorphogenesisWorkbenchState::agent_control_names();
        assert!(names.contains(&"Feed rate F"));
        assert!(names.contains(&"Kill rate k"));
        assert!(names.contains(&"Grid size"));
        assert!(names.contains(&"Sim speed (steps/frame)"));
        assert!(names.contains(&"Preset"));
    }

    #[test]
    fn viridis_is_monotone_dark_to_bright() {
        // The colormap brightens from t=0 (dark purple) to t=1 (yellow).
        let dark = viridis(0.0);
        let bright = viridis(1.0);
        let lum = |c: (u8, u8, u8)| c.0 as u32 + c.1 as u32 + c.2 as u32;
        assert!(lum(bright) > lum(dark), "viridis should brighten with t");
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use egui::accesskit::{Node, NodeId, Role};

    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_morphogenesis_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    fn has_named_node(nodes: &[(NodeId, Node)], name: &str) -> bool {
        nodes.iter().any(|(_, n)| n.name() == Some(name))
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_morphogenesis_workbench);
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_morphogenesis_workbench(&mut app, ctx);
        });
        // No panic = pass.
    }

    #[test]
    fn workbench_draws_when_shown() {
        let mut app = ValenxApp::default();
        app.show_morphogenesis_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);
        assert!(!nodes.is_empty(), "a shown workbench produces a11y nodes");
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // Every numeric DragValue is a SpinButton and must be `labelled_by` its
        // caption so an AI / screen reader can find it by caption text.
        let mut app = ValenxApp::default();
        app.show_morphogenesis_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let by_id: std::collections::HashMap<NodeId, &Node> =
            nodes.iter().map(|(id, n)| (*id, n)).collect();

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        assert!(
            spin_buttons.len() >= 4,
            "expected the four numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every DragValue must be labelled_by a caption (AI-drivable name)"
        );
        assert!(
            spin_buttons.iter().all(|n| {
                n.labelled_by()
                    .iter()
                    .any(|id| by_id.get(id).is_some_and(|t| t.name().is_some()))
            }),
            "every DragValue's labelled_by must point at a named caption node"
        );

        for caption in [
            "Feed rate F",
            "Kill rate k",
            "Grid size",
            "Sim speed (steps/frame)",
            "Preset",
        ] {
            assert!(
                has_named_node(&nodes, caption),
                "caption '{caption}' must be a named node in the a11y tree"
            );
        }
    }
}
