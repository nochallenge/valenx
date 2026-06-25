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
//! Movement + routes only: a **real OpenStreetMap tile basemap** (via the
//! in-process `walkers` slippy-map widget) with entity markers, polyline routes,
//! and real-time playback (Play / Pause + a playback-speed multiplier) overlaid
//! at true lat/lon. There is **no** engagement, sensor, or orbit modelling here —
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
use valenx_mission_sim::planner::PlannerScenario;
use walkers::sources::{OpenStreetMap, TileSource};
use walkers::{HttpTiles, Map, MapMemory, Position, Projector};

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
    /// OpenStreetMap raster-tile basemap source (real slippy-map tiles via
    /// [`walkers::HttpTiles`]). Lazily created on the first interactive frame
    /// because [`HttpTiles::new`] needs the live [`egui::Context`]; stays `None`
    /// in headless unit tests so no tile-download thread is ever spawned.
    pub tiles: Option<HttpTiles>,
    /// Pan / zoom state for the [`walkers::Map`] widget. When the user has not
    /// dragged the map, the camera follows the scenario centroid each frame; a
    /// drag detaches it (standard slippy-map behaviour).
    pub map_memory: MapMemory,
}

impl Default for MissionPlannerWorkbenchState {
    fn default() -> Self {
        Self {
            scenario: PlannerScenario::demo(DEFAULT_ENTITY_COUNT as usize),
            entity_count: DEFAULT_ENTITY_COUNT,
            playback_speed: 30.0,
            playing: false,
            tiles: None,
            map_memory: MapMemory::default(),
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
    draw_map(s, ui);
}

// ---------------------------------------------------------------------------
// Map rendering — real OpenStreetMap tile basemap (walkers) + entity overlay
// ---------------------------------------------------------------------------

/// Geographic centroid (mean lat / mean lon) of every entity and waypoint, used
/// as the map camera target. `None` when there is nothing to show.
fn scenario_centroid(sc: &PlannerScenario) -> Option<(f64, f64)> {
    let mut sum_lat = 0.0;
    let mut sum_lon = 0.0;
    let mut n = 0u32;
    let mut add = |lat: f64, lon: f64| {
        if lat.is_finite() && lon.is_finite() {
            sum_lat += lat;
            sum_lon += lon;
            n += 1;
        }
    };
    for e in &sc.entities {
        add(e.lat, e.lon);
        for wp in &e.route {
            add(wp.lat, wp.lon);
        }
    }
    (n > 0).then(|| (sum_lat / n as f64, sum_lon / n as f64))
}

/// The geographic extent (min/max lat & lon) over all entities and waypoints,
/// used to pick a zoom level that frames the whole scenario. `None` when empty.
fn scenario_extent(sc: &PlannerScenario) -> Option<(f64, f64, f64, f64)> {
    let mut min_lat = f64::INFINITY;
    let mut max_lat = f64::NEG_INFINITY;
    let mut min_lon = f64::INFINITY;
    let mut max_lon = f64::NEG_INFINITY;
    let mut seen = false;
    let mut expand = |lat: f64, lon: f64| {
        if lat.is_finite() && lon.is_finite() {
            min_lat = min_lat.min(lat);
            max_lat = max_lat.max(lat);
            min_lon = min_lon.min(lon);
            max_lon = max_lon.max(lon);
            seen = true;
        }
    };
    for e in &sc.entities {
        expand(e.lat, e.lon);
        for wp in &e.route {
            expand(wp.lat, wp.lon);
        }
    }
    seen.then_some((min_lat, max_lat, min_lon, max_lon))
}

/// Pick a slippy-map zoom (web-mercator `z`) whose tiles span the scenario
/// extent across roughly `view_px` pixels. Clamped to `[2, 18]`. Uses the
/// standard relation: at zoom `z` the whole 360° of longitude maps to
/// `256 · 2^z` pixels, so `z ≈ log2(view_px · 360 / (256 · span_lon))`. Returns
/// `f32` to feed [`walkers::MapMemory::set_zoom`] directly.
fn fit_zoom(sc: &PlannerScenario, view_px: f32) -> f32 {
    let Some((min_lat, max_lat, min_lon, max_lon)) = scenario_extent(sc) else {
        return 6.0;
    };
    // Pad the span so markers near the edge are not clipped, and guard a
    // single-point scenario (zero span) with a small floor.
    let span_lon = ((max_lon - min_lon).abs() * 1.3).max(0.02);
    let span_lat = ((max_lat - min_lat).abs() * 1.3).max(0.02);
    // Convert the lat span to an equivalent longitude span at this latitude so
    // the tighter of the two axes drives the zoom (keeps everything in view).
    let mid_lat_rad = (0.5 * (min_lat + max_lat)).to_radians();
    let span_lat_as_lon = span_lat / mid_lat_rad.cos().abs().max(1e-3);
    let span = span_lon.max(span_lat_as_lon);
    let z = ((view_px as f64) * 360.0 / (256.0 * span)).log2();
    z.floor().clamp(2.0, 18.0) as f32
}

/// Draw the geographic map: a **real OpenStreetMap raster-tile basemap** via the
/// [`walkers::Map`] widget, with the scenario's **routes** (polylines through
/// each entity's waypoints), **waypoint dots**, and **entity markers + name
/// labels** overlaid at their true lat/lon through walkers' [`Projector`]. The
/// camera follows the scenario centroid until the user drags the map; OSM
/// attribution is drawn in the corner per the tile-usage policy.
///
/// In headless unit tests `tiles` stays `None` (set lazily only on a live,
/// interactive frame), so the widget renders a plain basemap with the overlay
/// and **no tile-download thread is ever spawned**.
fn draw_map(s: &mut MissionPlannerWorkbenchState, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Map (OpenStreetMap basemap)").strong());

    // Camera target: the scenario centroid (or a sensible default if empty).
    let (centre_lat, centre_lon) = scenario_centroid(&s.scenario).unwrap_or((0.0, 0.0));
    let my_position = Position::from_lon_lat(centre_lon, centre_lat);

    // Allocate a fixed-height map viewport so the surrounding ScrollArea works.
    let available = ui.available_size();
    let map_size = egui::vec2(available.x.clamp(200.0, 720.0), 360.0);

    // Lazily create the OSM tile source on the first interactive frame (needs
    // the live egui Context). Skipped under #[cfg(test)] so headless tests stay
    // hermetic — no network, no Tokio thread. On first creation, frame the map
    // to the scenario by setting an initial fit-zoom on the MapMemory.
    #[cfg(not(test))]
    if s.tiles.is_none() {
        s.tiles = Some(HttpTiles::new(OpenStreetMap, ui.ctx().clone()));
        let z = fit_zoom(&s.scenario, map_size.x);
        // `set_zoom` only fails on an out-of-range argument; `fit_zoom` already
        // clamps to the valid slippy-map range, so this cannot error here.
        let _ = s.map_memory.set_zoom(z);
    }

    // Build the overlay plugin from a cheap snapshot of what to draw, so the
    // closure-free `Plugin` impl borrows nothing from `s` during `ui.add`.
    let overlay = MissionOverlay::from_scenario(&s.scenario);

    ui.allocate_ui(map_size, |ui| {
        let map = Map::new(
            s.tiles.as_mut().map(|t| t as &mut dyn walkers::Tiles),
            &mut s.map_memory,
            my_position,
        )
        .with_plugin(overlay);
        let response = ui.add(map);

        // OSM attribution (required by the tile-usage policy), bottom-left.
        let attr = OpenStreetMap.attribution();
        let painter = ui.painter_at(response.rect);
        let pos = response.rect.left_bottom() + egui::vec2(4.0, -4.0);
        let galley = painter.layout_no_wrap(
            format!("\u{00A9} {}", attr.text),
            egui::FontId::proportional(10.0),
            egui::Color32::from_gray(230),
        );
        let bg = egui::Rect::from_min_size(
            egui::pos2(pos.x - 2.0, pos.y - galley.size().y - 2.0),
            galley.size() + egui::vec2(4.0, 4.0),
        );
        painter.rect_filled(bg, 2.0, egui::Color32::from_black_alpha(150));
        painter.galley(
            egui::pos2(pos.x, pos.y - galley.size().y),
            galley,
            egui::Color32::from_gray(230),
        );
    });

    // Drag the map to detach the camera from the centroid; double-click resets.
    if ui
        .button("\u{1F3AF} Recenter on entities")
        .on_hover_text("Snap the map camera back to the scenario centroid and refit the zoom.")
        .clicked()
    {
        s.map_memory.follow_my_position();
        let z = fit_zoom(&s.scenario, map_size.x);
        let _ = s.map_memory.set_zoom(z);
    }
}

// ---------------------------------------------------------------------------
// Overlay plugin: routes + waypoints + entity markers drawn over the tiles
// ---------------------------------------------------------------------------

/// One entity's draw data, snapshotted from the scenario so the [`walkers::Plugin`]
/// owns its data and borrows nothing live during rendering.
struct OverlayEntity {
    name: String,
    lat: f64,
    lon: f64,
    done: bool,
    route: Vec<(f64, f64)>,
}

/// A [`walkers::Plugin`] that paints the mission scenario — every entity's route
/// polyline, its waypoint dots, and the entity marker + name label — on top of
/// the OSM tile basemap, projecting each lat/lon to screen pixels with walkers'
/// [`Projector`] so the overlay stays pinned to the map as it pans / zooms.
struct MissionOverlay {
    entities: Vec<OverlayEntity>,
}

impl MissionOverlay {
    fn from_scenario(sc: &PlannerScenario) -> Self {
        let entities = sc
            .entities
            .iter()
            .map(|e| OverlayEntity {
                name: e.name.clone(),
                lat: e.lat,
                lon: e.lon,
                done: e.is_done(),
                route: e.route.iter().map(|wp| (wp.lat, wp.lon)).collect(),
            })
            .collect();
        Self { entities }
    }
}

impl walkers::Plugin for MissionOverlay {
    fn run(&mut self, _response: &egui::Response, painter: egui::Painter, projector: &Projector) {
        // Project a (lat, lon) to an absolute screen position on the viewport.
        let to_px = |lat: f64, lon: f64| -> egui::Pos2 {
            projector
                .project(Position::from_lon_lat(lon, lat))
                .to_pos2()
        };

        let route_col = egui::Color32::from_rgb(70, 150, 240);
        let wp_col = egui::Color32::from_rgb(170, 200, 240);
        let entity_col = egui::Color32::from_rgb(240, 180, 70);
        let done_col = egui::Color32::from_rgb(120, 200, 140);

        // Routes (polylines) + waypoint dots.
        for e in &self.entities {
            let pts: Vec<egui::Pos2> = e.route.iter().map(|&(lat, lon)| to_px(lat, lon)).collect();
            for w in pts.windows(2) {
                painter.line_segment([w[0], w[1]], egui::Stroke::new(2.0, route_col));
            }
            for p in &pts {
                painter.circle_stroke(*p, 2.5, egui::Stroke::new(1.2, wp_col));
            }
        }

        // Entity markers (filled circle + black halo + name label).
        for e in &self.entities {
            let c = to_px(e.lat, e.lon);
            let col = if e.done { done_col } else { entity_col };
            painter.circle_filled(c, 5.0, col);
            painter.circle_stroke(c, 5.0, egui::Stroke::new(1.5, egui::Color32::BLACK));
            // Halo behind the label for legibility over busy tiles.
            painter.text(
                c + egui::vec2(8.0, -8.0) + egui::vec2(1.0, 1.0),
                egui::Align2::LEFT_BOTTOM,
                &e.name,
                egui::FontId::monospace(11.0),
                egui::Color32::from_black_alpha(200),
            );
            painter.text(
                c + egui::vec2(8.0, -8.0),
                egui::Align2::LEFT_BOTTOM,
                &e.name,
                egui::FontId::monospace(11.0),
                egui::Color32::WHITE,
            );
        }
    }
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
    fn centroid_and_extent_track_the_demo_scenario() {
        let sc = PlannerScenario::demo(4);
        let (clat, clon) = scenario_centroid(&sc).expect("demo has entities");
        let (min_lat, max_lat, min_lon, max_lon) = scenario_extent(&sc).expect("demo has entities");
        // The centroid must lie inside the extent box.
        assert!(
            (min_lat..=max_lat).contains(&clat),
            "centroid lat in extent"
        );
        assert!(
            (min_lon..=max_lon).contains(&clon),
            "centroid lon in extent"
        );
        assert!(
            max_lat >= min_lat && max_lon >= min_lon,
            "extent is well-ordered"
        );
    }

    #[test]
    fn centroid_is_none_for_empty_scenario() {
        // `demo(n)` clamps to >=1 entity, so build a genuinely empty scenario
        // directly to exercise the empty-extent fallback path.
        let sc = PlannerScenario {
            entities: Vec::new(),
            sim_time_s: 0.0,
        };
        assert!(scenario_centroid(&sc).is_none());
        assert!(scenario_extent(&sc).is_none());
    }

    #[test]
    fn fit_zoom_is_in_slippy_range() {
        let sc = PlannerScenario::demo(4);
        let z = fit_zoom(&sc, 480.0);
        assert!(
            (2.0..=18.0).contains(&z),
            "fit zoom {z} must be a valid slippy zoom"
        );
        // An empty scenario falls back to a sane default zoom, still in range.
        let empty = PlannerScenario {
            entities: Vec::new(),
            sim_time_s: 0.0,
        };
        let z_empty = fit_zoom(&empty, 480.0);
        assert!((2.0..=18.0).contains(&z_empty));
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
