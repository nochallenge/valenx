//! The right-side **Mission Planner** workbench — a native, in-house **visual
//! mission planner** modelled after ArduPilot Mission Planner / NASA GMAT.
//!
//! Entities sit on a geographic **map** (latitude / longitude) and follow ordered
//! **waypoint routes**; the workbench plays the plan back in **real time** so the
//! user watches the entities move along their legs. It is a front-end over
//! [`valenx_mission_sim::planner`], which owns the geographic frame, the entities
//! and routes, and the pure per-tick movement step.
//!
//! ## Stage 1 scope
//!
//! Movement + routes only: a lat/lon map with a graticule, entity markers,
//! polyline routes, and real-time playback (Play / Pause + a playback-speed
//! multiplier). There is **no** engagement, sensor, or orbit modelling here —
//! those are later stages. Purely defensive / planning posture: entities move
//! along routes, nothing more.
//!
//! Mirrors the other workbenches (`missionsim_workbench`): a
//! [`crate::workbench_chrome::workbench_shell`] panel gated on
//! [`crate::ValenxApp::show_mission_planner_workbench`], toggled from the View
//! menu and openable by the agent bridge under the workbench id
//! `"missionplanner"`. The agent bridge can set the controls
//! (`agent_set` / `agent_control_names`), read the sim-time readout
//! (`agent_readout`), and drive playback via the RunCommand ids
//! `missionplanner.play` / `missionplanner.pause`.

use eframe::egui;
use valenx_mission_sim::planner::{project, PlannerScenario};

use crate::ValenxApp;

/// The number of entities the demo scenario is seeded with by default.
const DEFAULT_ENTITY_COUNT: u32 = 4;

// ---------------------------------------------------------------------------
// Workbench state
// ---------------------------------------------------------------------------

/// Persistent state for the Mission Planner workbench: the live geographic
/// scenario plus the playback controls.
pub struct MissionPlannerWorkbenchState {
    /// The live planner scenario (entities + routes + simulated clock).
    pub scenario: PlannerScenario,
    /// Number of demo entities (re-seeds [`PlannerScenario::demo`] when changed).
    pub entity_count: u32,
    /// Real-time playback multiplier: each frame advances the scenario by
    /// `frame_dt · playback_speed` seconds. `1.0` = wall-clock real time.
    pub playback_speed: f64,
    /// Whether playback is running. While `true` the workbench steps the
    /// scenario every frame and requests a repaint so entities animate live.
    pub playing: bool,
}

impl Default for MissionPlannerWorkbenchState {
    fn default() -> Self {
        Self {
            scenario: PlannerScenario::demo(DEFAULT_ENTITY_COUNT as usize),
            entity_count: DEFAULT_ENTITY_COUNT,
            playback_speed: 30.0,
            playing: false,
        }
    }
}

impl MissionPlannerWorkbenchState {
    /// Re-seed the scenario with `entity_count` demo entities and reset the clock
    /// (also stops playback so the reset is visible before it resumes).
    fn reseed(&mut self) {
        self.scenario = PlannerScenario::demo(self.entity_count.max(1) as usize);
    }

    /// The user-visible captions of every control the agent bridge can set via
    /// `SetControl` (returned by `ListControls`). Order follows the form.
    pub fn agent_control_names() -> &'static [&'static str] {
        &["Entity count", "Playback speed x"]
    }

    /// Set one labelled control by its caption, for the agent `SetControl`
    /// bridge. Fail-loud on an unknown caption / wrong type; nothing is written
    /// on error and nothing panics. `Entity count` reads
    /// [`crate::agent_commands::AgentValue::as_i64`] (and re-seeds the scenario);
    /// `Playback speed x` is an `f64`.
    pub fn agent_set(
        &mut self,
        name: &str,
        value: &crate::agent_commands::AgentValue,
    ) -> Result<(), String> {
        match name {
            "Entity count" => {
                let n = value.as_i64()?;
                if !(1..=64).contains(&n) {
                    return Err(format!("Entity count must be in 1..=64, got {n}"));
                }
                self.entity_count = n as u32;
                self.reseed();
            }
            "Playback speed x" => {
                let v = value.as_f64()?;
                if !(0.0..=10_000.0).contains(&v) {
                    return Err(format!("Playback speed x must be in 0..=10000, got {v}"));
                }
                self.playback_speed = v;
            }
            other => return Err(format!("unknown mission-planner control: {other:?}")),
        }
        Ok(())
    }

    /// The current readout text for the agent `ReadReadout` bridge: the simulated
    /// time, entity count, and how many entities have finished their routes.
    /// Always `Some` (the scenario always exists).
    pub fn agent_readout(&self) -> Option<String> {
        let done = self
            .scenario
            .entities
            .iter()
            .filter(|e| e.is_done())
            .count();
        Some(format!(
            "Sim time {:.1}s \u{00B7} {} entities ({} arrived) \u{00B7} playback {}x \u{00B7} {}",
            self.scenario.sim_time_s,
            self.scenario.entities.len(),
            done,
            self.playback_speed,
            if self.playing { "playing" } else { "paused" },
        ))
    }
}

// ---------------------------------------------------------------------------
// Bridge run actions (play / pause)
// ---------------------------------------------------------------------------

/// Start real-time playback (the in-panel Play action). Factored out so the
/// button and the `missionplanner.play` bridge id share one path.
pub(crate) fn play(app: &mut ValenxApp) {
    app.mission_planner.playing = true;
}

/// Pause real-time playback (the in-panel Pause action). Factored out so the
/// button and the `missionplanner.pause` bridge id share one path.
pub(crate) fn pause(app: &mut ValenxApp) {
    app.mission_planner.playing = false;
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Draw the Mission Planner workbench. A no-op unless toggled on via
/// View → Mission Planner.
///
/// Mirrors [`crate::missionsim_workbench::draw_missionsim_workbench`].
pub fn draw_mission_planner_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_mission_planner_workbench {
        return;
    }
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_mission_planner_workbench",
        "Mission Planner (geographic / waypoint routes)",
        mission_planner_workbench_body,
    );
    if close {
        app.show_mission_planner_workbench = false;
    }
}

// ---------------------------------------------------------------------------
// Workbench body
// ---------------------------------------------------------------------------

fn mission_planner_workbench_body(app: &mut ValenxApp, ui: &mut egui::Ui) {
    egui::ScrollArea::vertical()
        .auto_shrink([false, false])
        .show(ui, |ui| mission_planner_workbench_body_inner(app, ui));
}

fn mission_planner_workbench_body_inner(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.label(
        egui::RichText::new(
            "In-house visual mission planner \u{00B7} valenx-mission-sim::planner \
             [geographic map + waypoint routes + real-time playback; movement only \u{2014} \
             no engagement / sensors / orbits]",
        )
        .weak()
        .small(),
    );
    ui.separator();

    // --- Controls -----------------------------------------------------------
    let s = &mut app.mission_planner;

    egui::Grid::new("mission_planner_controls")
        .num_columns(2)
        .striped(true)
        .show(ui, |ui| {
            // Entity count (re-seeds demo(N)).
            let lbl = ui.label("Entity count");
            let resp = ui
                .add(
                    egui::DragValue::new(&mut s.entity_count)
                        .speed(1)
                        .range(1..=64),
                )
                .labelled_by(lbl.id)
                .on_hover_text("Number of demo entities; changing this re-seeds the scenario.");
            if resp.changed() {
                s.reseed();
            }
            ui.end_row();

            // Playback speed multiplier.
            let lbl = ui.label("Playback speed x");
            ui.add(
                egui::DragValue::new(&mut s.playback_speed)
                    .speed(0.5)
                    .range(0.0..=10_000.0),
            )
            .labelled_by(lbl.id)
            .on_hover_text(
                "Real-time multiplier: each frame advances the sim by frame_dt x this value. \
                 1 = wall-clock; higher = fast-forward.",
            );
            ui.end_row();
        });

    ui.add_space(4.0);
    ui.horizontal(|ui| {
        // Play / Pause toggle (a single button whose label reflects state).
        let toggle_label = if s.playing {
            "\u{23F8} Pause"
        } else {
            "\u{25B6} Play"
        };
        if ui
            .button(toggle_label)
            .on_hover_text("Start / pause real-time playback of the routes.")
            .clicked()
        {
            s.playing = !s.playing;
        }
        if ui
            .button("\u{21BA} Reset")
            .on_hover_text("Re-seed the scenario and reset the simulated clock to 0.")
            .clicked()
        {
            s.reseed();
        }
    });

    // --- Real-time playback step --------------------------------------------
    // While playing, advance the scenario by (frame dt x speed) and request a
    // repaint so the next frame keeps animating. `stable_dt` is egui's smoothed
    // wall-clock frame time.
    if s.playing {
        let frame_dt = ui.input(|i| i.stable_dt) as f64;
        s.scenario.step(frame_dt * s.playback_speed);
        ui.ctx().request_repaint();
        if s.scenario.all_done() {
            // Everything has arrived — stop so we don't spin idle frames.
            s.playing = false;
        }
    }

    // --- Sim-time readout ----------------------------------------------------
    ui.add_space(6.0);
    let done = s.scenario.entities.iter().filter(|e| e.is_done()).count();
    ui.label(
        egui::RichText::new(format!(
            "Sim time (s): {:.1}   \u{00B7}   {} entities ({} arrived)",
            s.scenario.sim_time_s,
            s.scenario.entities.len(),
            done,
        ))
        .strong()
        .color(egui::Color32::from_rgb(120, 200, 140)),
    );

    // --- Map -----------------------------------------------------------------
    ui.add_space(6.0);
    ui.separator();
    draw_map(&s.scenario, ui);
}

// ---------------------------------------------------------------------------
// Map rendering (graticule + routes + entity markers)
// ---------------------------------------------------------------------------

/// Draw the geographic map: a lat/lon **graticule** (grid lines with labels), a
/// **scale hint**, each entity's **route** as a polyline through its waypoints,
/// and each **entity** as a filled circle + name label at its projected
/// position. Equirectangular projection (`x = lon`, `y = -lat`), auto-framed to
/// the data with a margin.
fn draw_map(sc: &PlannerScenario, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Map (equirectangular lat/lon)").strong());

    let available = ui.available_size();
    let (rect, _) = ui.allocate_exact_size(
        egui::vec2(available.x.clamp(200.0, 560.0), 320.0),
        egui::Sense::hover(),
    );
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 2.0, egui::Color32::from_rgb(12, 20, 30));

    // World extent over all waypoints + current entity positions.
    let mut min_lat = f64::INFINITY;
    let mut max_lat = f64::NEG_INFINITY;
    let mut min_lon = f64::INFINITY;
    let mut max_lon = f64::NEG_INFINITY;
    let mut expand = |lat: f64, lon: f64| {
        if lat.is_finite() && lon.is_finite() {
            min_lat = min_lat.min(lat);
            max_lat = max_lat.max(lat);
            min_lon = min_lon.min(lon);
            max_lon = max_lon.max(lon);
        }
    };
    for e in &sc.entities {
        expand(e.lat, e.lon);
        for wp in &e.route {
            expand(wp.lat, wp.lon);
        }
    }
    if !(min_lat.is_finite() && min_lon.is_finite() && max_lat.is_finite() && max_lon.is_finite()) {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "no entities",
            egui::FontId::monospace(12.0),
            egui::Color32::from_gray(120),
        );
        return;
    }

    // Pad the extent a little so markers near the edge stay inside.
    let span_lat = (max_lat - min_lat).max(0.01);
    let span_lon = (max_lon - min_lon).max(0.01);
    let pad_lat = span_lat * 0.12;
    let pad_lon = span_lon * 0.12;
    let (lo_lat, hi_lat) = (min_lat - pad_lat, max_lat + pad_lat);
    let (lo_lon, hi_lon) = (min_lon - pad_lon, max_lon + pad_lon);

    // Fit the world box into the drawing rect with a margin, preserving the
    // equirectangular shape (a single uniform deg->px scale).
    let margin = 28.0_f32;
    let inner = rect.shrink(margin);
    let world_w = (hi_lon - lo_lon) as f32;
    let world_h = (hi_lat - lo_lat) as f32;
    let scale = (inner.width() / world_w).min(inner.height() / world_h);
    // Centre the projected world in `inner`.
    let proj_w = world_w * scale;
    let proj_h = world_h * scale;
    let ox = inner.left() + (inner.width() - proj_w) / 2.0;
    let oy = inner.top() + (inner.height() - proj_h) / 2.0;
    // `project` maps (lat, lon) about an origin pixel for lat=lon=0; build that
    // origin so (lo_lon, hi_lat) lands at the top-left of the centred box.
    let origin = (ox - lo_lon as f32 * scale, oy + hi_lat as f32 * scale);
    let to_px = |lat: f64, lon: f64| -> egui::Pos2 {
        let (x, y) = project(lat, lon, origin, scale);
        egui::pos2(x, y)
    };

    // --- Graticule (grid lines every `step` degrees, with labels) -----------
    let grid_col = egui::Color32::from_rgb(40, 60, 80);
    let label_col = egui::Color32::from_gray(150);
    let step = graticule_step(span_lat.max(span_lon));
    // Longitude lines (vertical).
    let first_lon = (lo_lon / step).ceil() * step;
    let mut lon = first_lon;
    while lon <= hi_lon {
        let top = to_px(hi_lat, lon);
        let bot = to_px(lo_lat, lon);
        painter.line_segment([top, bot], egui::Stroke::new(0.7, grid_col));
        painter.text(
            egui::pos2(bot.x, rect.bottom() - 11.0),
            egui::Align2::CENTER_CENTER,
            format!("{lon:.2}\u{00B0}"),
            egui::FontId::monospace(9.0),
            label_col,
        );
        lon += step;
    }
    // Latitude lines (horizontal).
    let first_lat = (lo_lat / step).ceil() * step;
    let mut lat = first_lat;
    while lat <= hi_lat {
        let left = to_px(lat, lo_lon);
        let right = to_px(lat, hi_lon);
        painter.line_segment([left, right], egui::Stroke::new(0.7, grid_col));
        painter.text(
            egui::pos2(rect.left() + 3.0, left.y),
            egui::Align2::LEFT_CENTER,
            format!("{lat:.2}\u{00B0}"),
            egui::FontId::monospace(9.0),
            label_col,
        );
        lat += step;
    }
    // Frame the region.
    painter.rect_stroke(
        rect,
        2.0,
        egui::Stroke::new(1.0, egui::Color32::from_rgb(60, 90, 120)),
    );

    // --- Scale hint (one grid cell in km) -----------------------------------
    // One degree of latitude ~= 111.32 km; a graticule cell is `step` degrees.
    let cell_km = step * 111.32;
    painter.text(
        egui::pos2(rect.right() - 4.0, rect.top() + 10.0),
        egui::Align2::RIGHT_CENTER,
        format!("grid {step:.2}\u{00B0} \u{2248} {cell_km:.0} km"),
        egui::FontId::monospace(9.0),
        label_col,
    );

    // --- Routes (polylines) + waypoint dots ---------------------------------
    let route_col = egui::Color32::from_rgb(90, 130, 170);
    let wp_col = egui::Color32::from_rgb(150, 180, 210);
    for e in &sc.entities {
        let pts: Vec<egui::Pos2> = e.route.iter().map(|wp| to_px(wp.lat, wp.lon)).collect();
        for w in pts.windows(2) {
            painter.line_segment([w[0], w[1]], egui::Stroke::new(1.2, route_col));
        }
        for p in &pts {
            painter.circle_stroke(*p, 2.0, egui::Stroke::new(1.0, wp_col));
        }
    }

    // --- Entity markers (filled circle + name) ------------------------------
    let entity_col = egui::Color32::from_rgb(240, 180, 70);
    let done_col = egui::Color32::from_rgb(120, 200, 140);
    for e in &sc.entities {
        let c = to_px(e.lat, e.lon);
        let col = if e.is_done() { done_col } else { entity_col };
        painter.circle_filled(c, 4.0, col);
        painter.circle_stroke(c, 4.0, egui::Stroke::new(1.0, egui::Color32::BLACK));
        painter.text(
            c + egui::vec2(6.0, -6.0),
            egui::Align2::LEFT_BOTTOM,
            &e.name,
            egui::FontId::monospace(10.0),
            egui::Color32::from_gray(220),
        );
    }
}

/// Choose a "nice" graticule step (degrees) for a world span, so the map shows a
/// handful of grid lines rather than too many / too few.
fn graticule_step(span_deg: f64) -> f64 {
    // Aim for ~5 divisions across the larger span.
    let raw = (span_deg / 5.0).max(1e-3);
    // Snap up to 1, 2, 5 x 10^k.
    let pow = 10f64.powf(raw.log10().floor());
    let m = raw / pow;
    let nice = if m <= 1.0 {
        1.0
    } else if m <= 2.0 {
        2.0
    } else if m <= 5.0 {
        5.0
    } else {
        10.0
    };
    nice * pow
}

// ---------------------------------------------------------------------------
// Tests (unit)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_seeds_a_scenario() {
        let s = MissionPlannerWorkbenchState::default();
        assert_eq!(s.entity_count, DEFAULT_ENTITY_COUNT);
        assert_eq!(s.scenario.entities.len(), DEFAULT_ENTITY_COUNT as usize);
        assert!(!s.playing);
    }

    #[test]
    fn agent_set_entity_count_reseeds() {
        use crate::agent_commands::AgentValue;
        let mut s = MissionPlannerWorkbenchState::default();
        s.agent_set("Entity count", &AgentValue::Int(7))
            .expect("valid entity count");
        assert_eq!(s.entity_count, 7);
        assert_eq!(s.scenario.entities.len(), 7);
    }

    #[test]
    fn agent_set_rejects_out_of_range_and_unknown() {
        use crate::agent_commands::AgentValue;
        let mut s = MissionPlannerWorkbenchState::default();
        assert!(s.agent_set("Entity count", &AgentValue::Int(0)).is_err());
        assert!(s.agent_set("Entity count", &AgentValue::Int(999)).is_err());
        assert!(s.agent_set("nope", &AgentValue::Float(1.0)).is_err());
    }

    #[test]
    fn agent_set_playback_speed() {
        use crate::agent_commands::AgentValue;
        let mut s = MissionPlannerWorkbenchState::default();
        s.agent_set("Playback speed x", &AgentValue::Float(12.5))
            .expect("valid speed");
        assert!((s.playback_speed - 12.5).abs() < 1e-9);
    }

    #[test]
    fn agent_readout_reports_sim_time() {
        let mut s = MissionPlannerWorkbenchState::default();
        s.scenario.step(10.0);
        let r = s.agent_readout().expect("readout always present");
        assert!(
            r.contains("Sim time"),
            "readout should mention sim time: {r}"
        );
    }

    #[test]
    fn play_pause_helpers_toggle_flag() {
        let mut app = ValenxApp::default();
        assert!(!app.mission_planner.playing);
        play(&mut app);
        assert!(app.mission_planner.playing);
        pause(&mut app);
        assert!(!app.mission_planner.playing);
    }

    #[test]
    fn control_names_are_listed() {
        let names = MissionPlannerWorkbenchState::agent_control_names();
        assert!(names.contains(&"Entity count"));
        assert!(names.contains(&"Playback speed x"));
    }

    #[test]
    fn graticule_step_is_nice() {
        // A ~1.4 deg span -> ~5 divisions -> step around 0.2-0.5.
        let st = graticule_step(1.4);
        assert!(
            st > 0.0 && st <= 0.5,
            "step {st} should be a small nice value"
        );
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
            draw_mission_planner_workbench(app, ctx);
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
        assert!(!app.show_mission_planner_workbench);
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_mission_planner_workbench(&mut app, ctx);
        });
        // No panic = pass.
    }

    #[test]
    fn workbench_draws_when_shown() {
        let mut app = ValenxApp::default();
        app.show_mission_planner_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);
        assert!(!nodes.is_empty(), "a shown workbench produces a11y nodes");
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // Every numeric DragValue is a SpinButton and must be `labelled_by` its
        // caption so an AI / screen reader can find it by caption text.
        let mut app = ValenxApp::default();
        app.show_mission_planner_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let by_id: std::collections::HashMap<NodeId, &Node> =
            nodes.iter().map(|(id, n)| (*id, n)).collect();

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        assert!(
            spin_buttons.len() >= 2,
            "expected the two numeric controls as spin buttons, got {}",
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

        for caption in ["Entity count", "Playback speed x"] {
            assert!(
                has_named_node(&nodes, caption),
                "caption '{caption}' must be a named node in the a11y tree"
            );
        }
    }
}
