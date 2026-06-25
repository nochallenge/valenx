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
use valenx_mission_sim::los::{line_of_sight, line_of_sight_terrain};
use valenx_mission_sim::planner::{Affiliation, PlannerScenario, M_PER_DEG_LAT};
use valenx_mission_sim::routing::{astar, demo_field, CostGrid};
use valenx_mission_sim::terrain::{demo_terrain, HeightGrid};
use walkers::sources::{OpenStreetMap, TileSource};
use walkers::{HttpTiles, Map, MapMemory, Position, Projector};

use crate::ValenxApp;

/// The number of entities the demo scenario is seeded with by default.
const DEFAULT_ENTITY_COUNT: u32 = 4;

/// Width (columns) of the tactical-routing demo cost grid laid over the map.
const ROUTE_GRID_W: usize = 32;
/// Height (rows) of the tactical-routing demo cost grid laid over the map.
const ROUTE_GRID_H: usize = 20;

/// Base traversal cost of a flat cell when the cost field is derived from terrain
/// slope (`cost = TERRAIN_BASE + TERRAIN_SLOPE_K · slope`).
const TERRAIN_BASE: f32 = 1.0;
/// Slope-to-cost gain: how much each unit of terrain slope (rise-over-run) adds to
/// a cell's traversal cost, so A\* prefers gentle ground.
const TERRAIN_SLOPE_K: f32 = 6.0;
/// Slope above which terrain is **impassable** (a cliff / sheer ridge flank): such
/// cells become `f32::INFINITY` so the route goes around them. Tuned to the demo
/// landscape's steep ridge flanks (~5-6) while leaving valleys (~0.9) traversable.
const TERRAIN_IMPASSABLE_SLOPE: f32 = 4.0;

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
    /// The affiliation currently selected in the "Set all to" combo. Applying it
    /// (the combo `changed`, or the `Affiliation (all)` agent control) overrides
    /// every entity's [`Affiliation`]; purely a display/iconography choice.
    pub set_all_affiliation: Affiliation,

    // --- Tactical routing (in-house A* over a cost field) -------------------
    /// The demo traversal-cost grid (with obstacle walls) the A\* route is
    /// planned over. `ROUTE_GRID_W × ROUTE_GRID_H` cells, mapped onto the current
    /// map extent at draw time so the route overlays the basemap geographically.
    pub route_grid: CostGrid,
    /// Routing start cell `(x, y)` in grid coordinates (`0..W`, `0..H`).
    pub route_start: (usize, usize),
    /// Routing goal cell `(x, y)` in grid coordinates.
    pub route_goal: (usize, usize),
    /// The most recently computed A\* path as grid cells (`start..=goal`), or
    /// `None` if not yet computed or the goal was unreachable.
    pub route: Option<Vec<(usize, usize)>>,
    /// Human-readable status of the last `Compute route` (length in cells / km,
    /// or "no route"). Surfaced in the readout and the agent readout.
    pub route_status: String,

    // --- Line of sight (in-house DDA ray-march over the same cost grid) ------
    /// The observer / sensor cell `(x, y)` line-of-sight is computed FROM, over
    /// the same [`Self::route_grid`] occupancy. Defaults to the route start.
    pub los_observer: (usize, usize),
    /// The last computed visibility to each target cell, paired as
    /// `(target_cell, visible)`: `true` = clear sight line (GREEN), `false` =
    /// masked by an obstacle (RED). Empty until `Compute LoS` runs. The targets
    /// are the route goal plus each entity's current cell.
    pub los_results: Vec<((usize, usize), bool)>,
    /// Human-readable status of the last `Compute LoS` (visible / blocked
    /// counts). Surfaced in the readout and the agent readout.
    pub los_status: String,

    // --- Terrain elevation (in-house procedural heightfield) ----------------
    /// The procedural elevation heightfield (metres) laid under the routes / LoS,
    /// same `ROUTE_GRID_W × ROUTE_GRID_H` extent as [`Self::route_grid`] so the two
    /// overlay cell-for-cell. When [`Self::terrain_on`] the routing cost field is
    /// derived from this terrain's slope and line-of-sight is terrain-masked (2.5-D).
    pub terrain: HeightGrid,
    /// Whether terrain-awareness is active: when `true`, the routing cost grid is
    /// derived from terrain slope (gentle = cheap, steep ridge = impassable) and
    /// LoS uses the 2.5-D elevation ray-march; when `false`, routing/LoS fall back
    /// to the flat obstacle-wall [`demo_field`]. Drawn as a colour shade either way.
    pub terrain_on: bool,
    /// Observer height above the ground in metres, used by terrain-masked LoS (a
    /// taller observer sees over low hills). Display + LoS only.
    pub obs_height_m: f32,
    /// Target height above the ground in metres, used by terrain-masked LoS.
    pub tgt_height_m: f32,
}

impl Default for MissionPlannerWorkbenchState {
    fn default() -> Self {
        // Build the procedural terrain and derive the routing cost field from its
        // slope (terrain-aware by default — this is the capstone tying routing +
        // LoS to elevation). The grid extent matches the route grid cell-for-cell.
        let terrain = demo_terrain(ROUTE_GRID_W, ROUTE_GRID_H);
        let route_grid = CostGrid::from_terrain(
            &terrain,
            TERRAIN_BASE,
            TERRAIN_SLOPE_K,
            TERRAIN_IMPASSABLE_SLOPE,
        );
        Self {
            scenario: PlannerScenario::demo(DEFAULT_ENTITY_COUNT as usize),
            entity_count: DEFAULT_ENTITY_COUNT,
            playback_speed: 30.0,
            playing: false,
            tiles: None,
            map_memory: MapMemory::default(),
            set_all_affiliation: Affiliation::Friendly,
            route_grid,
            // Default cross-ridge planning problem: SW area (NW side of the ridge)
            // to NE area (SE side). The diagonal ridge lies between them, so the
            // slope-aware route must thread the mountain PASS rather than climb the
            // crest — the capstone demo of terrain-aware routing.
            route_start: (1, ROUTE_GRID_H - 2),
            route_goal: (ROUTE_GRID_W - 1, 1),
            route: None,
            route_status: "no route computed".to_string(),
            // Observe from the route start by default; LoS to the goal + entities.
            los_observer: (1, ROUTE_GRID_H - 2),
            los_results: Vec::new(),
            los_status: "no LoS computed".to_string(),
            terrain,
            terrain_on: true,
            obs_height_m: 2.0,
            tgt_height_m: 0.0,
        }
    }
}

impl MissionPlannerWorkbenchState {
    /// Re-seed the scenario with `entity_count` demo entities and reset the clock
    /// (also stops playback so the reset is visible before it resumes).
    fn reseed(&mut self) {
        self.scenario = PlannerScenario::demo(self.entity_count.max(1) as usize);
    }

    /// Rebuild the routing [`Self::route_grid`] to match the current terrain mode:
    /// when [`Self::terrain_on`], derive the cost field from the terrain's slope
    /// (gentle = cheap, steep ridge = impassable); otherwise fall back to the flat
    /// obstacle-wall [`demo_field`]. Invalidates any stale computed route so the
    /// map does not show a path planned over the previous field. Called whenever
    /// the terrain toggle flips.
    fn rebuild_route_grid(&mut self) {
        self.route_grid = if self.terrain_on {
            CostGrid::from_terrain(
                &self.terrain,
                TERRAIN_BASE,
                TERRAIN_SLOPE_K,
                TERRAIN_IMPASSABLE_SLOPE,
            )
        } else {
            demo_field(ROUTE_GRID_W, ROUTE_GRID_H)
        };
        // The old route was planned over a different cost field — clear it.
        self.route = None;
        self.route_status = "no route computed".to_string();
    }

    /// The user-visible captions of every control the agent bridge can set via
    /// `SetControl` (returned by `ListControls`). Order follows the form.
    pub fn agent_control_names() -> &'static [&'static str] {
        &[
            "Entity count",
            "Playback speed x",
            "Affiliation (all)",
            "Route start X",
            "Route start Y",
            "Route goal X",
            "Route goal Y",
            "Observer X",
            "Observer Y",
            "Terrain",
            "Observer height (m)",
            "Target height (m)",
        ]
    }

    /// Compute the A\* tactical route over the demo cost field from
    /// [`Self::route_start`] to [`Self::route_goal`], storing the cell path in
    /// [`Self::route`] and a human-readable summary in [`Self::route_status`].
    /// Shared by the in-panel `Compute route` button and the
    /// `missionplanner.route` bridge id. Pure aside from updating `self`.
    pub(crate) fn compute_route(&mut self) {
        self.route = astar(&self.route_grid, self.route_start, self.route_goal);
        self.route_status = match &self.route {
            Some(path) => {
                // Approx ground length: each step is one cell; map the grid to the
                // current scenario extent to size a cell in km (diagonal = √2·cell).
                let km = self.route_length_km(path);
                format!(
                    "route: {} cells \u{00B7} ~{:.1} km \u{00B7} ({},{}) \u{2192} ({},{})",
                    path.len(),
                    km,
                    self.route_start.0,
                    self.route_start.1,
                    self.route_goal.0,
                    self.route_goal.1,
                )
            }
            None => format!(
                "no route: ({},{}) \u{2192} ({},{}) unreachable / invalid",
                self.route_start.0, self.route_start.1, self.route_goal.0, self.route_goal.1
            ),
        };
    }

    /// Approximate the ground length of a grid path in kilometres by mapping the
    /// cell grid onto the current scenario geographic extent and summing the
    /// great-circle-ish (equirectangular) length of each leg. Returns `0.0` for a
    /// path of fewer than two cells.
    fn route_length_km(&self, path: &[(usize, usize)]) -> f64 {
        if path.len() < 2 {
            return 0.0;
        }
        let bbox = route_bbox(&self.scenario);
        let mut total_m = 0.0;
        for w in path.windows(2) {
            let (a_lat, a_lon) = cell_to_latlon(w[0], &self.route_grid, bbox);
            let (b_lat, b_lon) = cell_to_latlon(w[1], &self.route_grid, bbox);
            let cos_lat = ((a_lat + b_lat) * 0.5).to_radians().cos();
            let north_m = (b_lat - a_lat) * M_PER_DEG_LAT;
            let east_m = (b_lon - a_lon) * M_PER_DEG_LAT * cos_lat;
            total_m += (north_m * north_m + east_m * east_m).sqrt();
        }
        total_m / 1000.0
    }

    /// Compute **line of sight** from [`Self::los_observer`] to a set of target
    /// cells over the same [`Self::route_grid`] occupancy, storing
    /// `(target, visible)` pairs in [`Self::los_results`] and a visible/blocked
    /// summary in [`Self::los_status`]. Shared by the in-panel `Compute LoS`
    /// button and the `missionplanner.los` bridge id.
    ///
    /// Targets are the routing **goal** cell plus each **entity's** current cell
    /// (entity lat/lon snapped to the nearest grid cell over the routing bbox).
    /// Targets that coincide with the observer are skipped (a cell trivially sees
    /// itself). Pure aside from updating `self`.
    pub(crate) fn compute_los(&mut self) {
        let bbox = route_bbox(&self.scenario);
        let mut targets: Vec<(usize, usize)> = Vec::new();
        // The routing goal is always a target of interest.
        targets.push(self.route_goal);
        // Each entity's current position, snapped to a grid cell.
        for e in &self.scenario.entities {
            targets.push(latlon_to_cell(e.lat, e.lon, &self.route_grid, bbox));
        }
        // De-duplicate and drop the observer's own cell (trivially visible).
        targets.sort_unstable();
        targets.dedup();
        targets.retain(|&c| c != self.los_observer);

        // Terrain-aware (2.5-D dead-ground masking via the elevation ray-march)
        // when terrain is on; otherwise flat occupancy over the obstacle field.
        self.los_results = targets
            .iter()
            .map(|&t| {
                let vis = if self.terrain_on {
                    line_of_sight_terrain(
                        &self.terrain,
                        self.los_observer,
                        t,
                        self.obs_height_m,
                        self.tgt_height_m,
                    )
                } else {
                    line_of_sight(&self.route_grid, self.los_observer, t)
                };
                (t, vis)
            })
            .collect();

        let visible = self.los_results.iter().filter(|(_, v)| *v).count();
        let blocked = self.los_results.len() - visible;
        self.los_status = format!(
            "LoS ({}) from ({},{}): {} visible \u{00B7} {} blocked \u{00B7} {} targets",
            if self.terrain_on {
                "terrain-masked 2.5-D"
            } else {
                "flat"
            },
            self.los_observer.0,
            self.los_observer.1,
            visible,
            blocked,
            self.los_results.len(),
        );
    }

    /// Set every entity's APP-6 affiliation at once (the per-side override). Used
    /// by both the GUI combo and the `Affiliation (all)` agent control.
    fn set_all_affiliations(&mut self, a: Affiliation) {
        for e in &mut self.scenario.entities {
            e.affiliation = a;
        }
    }

    /// Count of entities per APP-6 affiliation, in canonical
    /// [`Affiliation::ALL`] order (friendly, hostile, neutral, unknown).
    fn affiliation_counts(&self) -> [usize; 4] {
        let mut counts = [0usize; 4];
        for e in &self.scenario.entities {
            for (i, a) in Affiliation::ALL.iter().enumerate() {
                if e.affiliation == *a {
                    counts[i] += 1;
                }
            }
        }
        counts
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
            "Affiliation (all)" => {
                // Enum-by-name: set every entity to one APP-6 affiliation.
                let s = value.as_str()?;
                let a = Affiliation::from_name(s).ok_or_else(|| {
                    format!(
                        "Affiliation (all) must be one of friendly/hostile/neutral/unknown, got {s:?}"
                    )
                })?;
                self.set_all_affiliations(a);
            }
            "Route start X" => {
                self.route_start.0 = route_coord(value, ROUTE_GRID_W, "Route start X")?
            }
            "Route start Y" => {
                self.route_start.1 = route_coord(value, ROUTE_GRID_H, "Route start Y")?
            }
            "Route goal X" => self.route_goal.0 = route_coord(value, ROUTE_GRID_W, "Route goal X")?,
            "Route goal Y" => self.route_goal.1 = route_coord(value, ROUTE_GRID_H, "Route goal Y")?,
            "Observer X" => self.los_observer.0 = route_coord(value, ROUTE_GRID_W, "Observer X")?,
            "Observer Y" => self.los_observer.1 = route_coord(value, ROUTE_GRID_H, "Observer Y")?,
            "Terrain" => {
                // Bool toggle: rebuild the routing cost field for the new mode so
                // routing + LoS become (or stop being) terrain-aware.
                let on = value.as_bool()?;
                if on != self.terrain_on {
                    self.terrain_on = on;
                    self.rebuild_route_grid();
                }
            }
            "Observer height (m)" => {
                let v = value.as_f64()? as f32;
                if !(0.0..=10_000.0).contains(&v) {
                    return Err(format!("Observer height (m) must be in 0..=10000, got {v}"));
                }
                self.obs_height_m = v;
            }
            "Target height (m)" => {
                let v = value.as_f64()? as f32;
                if !(0.0..=10_000.0).contains(&v) {
                    return Err(format!("Target height (m) must be in 0..=10000, got {v}"));
                }
                self.tgt_height_m = v;
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
        let [fr, ho, ne, un] = self.affiliation_counts();
        let (elev_lo, elev_hi) = self.terrain.minmax();
        Some(format!(
            "Sim time {:.1}s \u{00B7} {} entities ({} arrived) \u{00B7} \
             affiliation F{} H{} N{} U{} \u{00B7} playback {}x \u{00B7} {} \u{00B7} \
             terrain {} (elev {:.0}\u{2013}{:.0} m; routing + LoS are terrain-aware) \u{00B7} \
             {} \u{00B7} {}",
            self.scenario.sim_time_s,
            self.scenario.entities.len(),
            done,
            fr,
            ho,
            ne,
            un,
            self.playback_speed,
            if self.playing { "playing" } else { "paused" },
            if self.terrain_on { "ON" } else { "OFF" },
            elev_lo,
            elev_hi,
            self.route_status,
            self.los_status,
        ))
    }
}

/// Validate an `AgentValue` as a routing grid coordinate in `0..extent` and
/// return it as a `usize`. Fail-loud on a non-integer or out-of-range value so a
/// bad `SetControl` is rejected without mutating state.
fn route_coord(
    value: &crate::agent_commands::AgentValue,
    extent: usize,
    caption: &str,
) -> Result<usize, String> {
    let n = value.as_i64()?;
    if n < 0 || n as usize >= extent {
        return Err(format!("{caption} must be in 0..{extent}, got {n}"));
    }
    Ok(n as usize)
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

/// Compute the A\* tactical route (the in-panel `Compute route` action).
/// Factored out so the button and the `missionplanner.route` bridge id share one
/// path; delegates to [`MissionPlannerWorkbenchState::compute_route`].
pub(crate) fn route(app: &mut ValenxApp) {
    app.mission_planner.compute_route();
}

/// Compute line of sight from the observer to the goal + entities (the in-panel
/// `Compute LoS` action). Factored out so the button and the
/// `missionplanner.los` bridge id share one path; delegates to
/// [`MissionPlannerWorkbenchState::compute_los`].
pub(crate) fn los(app: &mut ValenxApp) {
    app.mission_planner.compute_los();
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

            // APP-6 affiliation: set ALL entities' symbol frame at once.
            let lbl = ui.label("Affiliation (all)");
            let mut chosen = s.set_all_affiliation;
            egui::ComboBox::from_id_source("mission_planner_affiliation")
                .selected_text(affiliation_label(chosen))
                .show_ui(ui, |ui| {
                    for a in Affiliation::ALL {
                        ui.selectable_value(&mut chosen, a, affiliation_label(a));
                    }
                })
                .response
                .labelled_by(lbl.id)
                .on_hover_text(
                    "MIL-STD-2525 / APP-6 standard identity for every entity's map symbol: \
                     Friendly = blue rounded rectangle, Hostile = red diamond, \
                     Neutral = green square, Unknown = yellow quatrefoil. Map iconography only.",
                );
            if chosen != s.set_all_affiliation {
                s.set_all_affiliation = chosen;
                s.set_all_affiliations(chosen);
            }
            ui.end_row();
        });

    // --- Tactical routing (A* over the cost field) --------------------------
    ui.add_space(6.0);
    ui.separator();
    ui.label(
        egui::RichText::new(
            "Tactical routing \u{00B7} in-house A* over a grid cost field with obstacles \
             [defensive route planning \u{2014} finds a least-cost path AROUND obstacles]",
        )
        .weak()
        .small(),
    );
    egui::Grid::new("mission_planner_routing")
        .num_columns(4)
        .striped(true)
        .show(ui, |ui| {
            // Start cell (X, Y).
            let lbl = ui.label("Route start X");
            ui.add(
                egui::DragValue::new(&mut s.route_start.0)
                    .speed(1)
                    .range(0..=ROUTE_GRID_W - 1),
            )
            .labelled_by(lbl.id)
            .on_hover_text("Start cell column (0..grid width) over the cost field.");
            let lbl = ui.label("Route start Y");
            ui.add(
                egui::DragValue::new(&mut s.route_start.1)
                    .speed(1)
                    .range(0..=ROUTE_GRID_H - 1),
            )
            .labelled_by(lbl.id)
            .on_hover_text("Start cell row (0 = north edge).");
            ui.end_row();

            // Goal cell (X, Y).
            let lbl = ui.label("Route goal X");
            ui.add(
                egui::DragValue::new(&mut s.route_goal.0)
                    .speed(1)
                    .range(0..=ROUTE_GRID_W - 1),
            )
            .labelled_by(lbl.id)
            .on_hover_text("Goal cell column (0..grid width).");
            let lbl = ui.label("Route goal Y");
            ui.add(
                egui::DragValue::new(&mut s.route_goal.1)
                    .speed(1)
                    .range(0..=ROUTE_GRID_H - 1),
            )
            .labelled_by(lbl.id)
            .on_hover_text("Goal cell row (grid height - 1 = south edge).");
            ui.end_row();
        });

    if ui
        .button("\u{1F9ED} Compute route")
        .on_hover_text(
            "Run in-house A* over the demo cost field (with obstacle walls) from start to goal; \
             the least-cost path is drawn as an amber polyline on the map.",
        )
        .clicked()
    {
        s.compute_route();
    }
    // Route status readout (cells / approx km, or 'no route').
    ui.label(
        egui::RichText::new(&s.route_status)
            .strong()
            .color(egui::Color32::from_rgb(255, 180, 90)),
    );

    // --- Terrain elevation (in-house procedural heightfield) ----------------
    ui.add_space(6.0);
    ui.separator();
    ui.label(
        egui::RichText::new(
            "Terrain elevation \u{00B7} in-house procedural heightfield \
             [makes routing slope-aware + line-of-sight terrain-masked in 2.5-D; \
             shaded green\u{2192}tan\u{2192}brown by height under the map]",
        )
        .weak()
        .small(),
    );
    let (elev_lo, elev_hi) = s.terrain.minmax();
    egui::Grid::new("mission_planner_terrain")
        .num_columns(4)
        .striped(true)
        .show(ui, |ui| {
            // Terrain on/off toggle: flips routing/LoS between terrain-aware and
            // the flat obstacle field, and rebuilds the cost grid accordingly.
            let lbl = ui.label("Terrain");
            let resp = ui
                .checkbox(&mut s.terrain_on, "slope-aware routing + masked LoS")
                .labelled_by(lbl.id)
                .on_hover_text(
                    "When on, the routing cost field is derived from terrain slope \
                     (gentle = cheap, steep ridge = impassable) and line-of-sight is \
                     terrain-masked in 2.5-D (dead ground behind ridges). When off, \
                     routing/LoS use the flat obstacle field.",
                );
            if resp.changed() {
                s.rebuild_route_grid();
            }
            ui.end_row();

            // Observer height above ground (metres) for terrain-masked LoS.
            let lbl = ui.label("Observer height (m)");
            ui.add(
                egui::DragValue::new(&mut s.obs_height_m)
                    .speed(0.5)
                    .range(0.0..=10_000.0),
            )
            .labelled_by(lbl.id)
            .on_hover_text(
                "Observer eye height above the ground, in metres. A taller observer \
                 (e.g. a mast or high ground) sees OVER low hills in terrain-masked LoS.",
            );
            // Target height above ground (metres).
            let lbl = ui.label("Target height (m)");
            ui.add(
                egui::DragValue::new(&mut s.tgt_height_m)
                    .speed(0.5)
                    .range(0.0..=10_000.0),
            )
            .labelled_by(lbl.id)
            .on_hover_text("Target height above the ground, in metres, for terrain-masked LoS.");
            ui.end_row();
        });
    // Terrain elevation-range readout.
    ui.label(
        egui::RichText::new(format!(
            "Elevation range: {:.0}\u{2013}{:.0} m \u{00B7} {} \u{00B7} {}\u{00D7}{} grid",
            elev_lo,
            elev_hi,
            if s.terrain_on {
                "routing + LoS terrain-aware"
            } else {
                "terrain display only (flat routing/LoS)"
            },
            s.terrain.w,
            s.terrain.h,
        ))
        .strong()
        .color(egui::Color32::from_rgb(180, 200, 150)),
    );

    // --- Line of sight (in-house DDA ray-march over the cost field) ---------
    ui.add_space(6.0);
    ui.separator();
    ui.label(
        egui::RichText::new(
            "Line of sight \u{00B7} in-house 2-D ray-march over the same grid \
             [sensor visibility \u{2014} which targets are SEEN vs MASKED by terrain / obstacles]",
        )
        .weak()
        .small(),
    );
    egui::Grid::new("mission_planner_los")
        .num_columns(4)
        .striped(true)
        .show(ui, |ui| {
            let lbl = ui.label("Observer X");
            ui.add(
                egui::DragValue::new(&mut s.los_observer.0)
                    .speed(1)
                    .range(0..=ROUTE_GRID_W - 1),
            )
            .labelled_by(lbl.id)
            .on_hover_text("Observer / sensor cell column (0..grid width) over the cost field.");
            let lbl = ui.label("Observer Y");
            ui.add(
                egui::DragValue::new(&mut s.los_observer.1)
                    .speed(1)
                    .range(0..=ROUTE_GRID_H - 1),
            )
            .labelled_by(lbl.id)
            .on_hover_text("Observer / sensor cell row (0 = north edge).");
            ui.end_row();
        });
    ui.horizontal(|ui| {
        if ui
            .button("\u{1F441} Compute LoS")
            .on_hover_text(
                "Ray-march line of sight from the observer to the route goal and every entity \
                 over the demo cost field; clear lines draw GREEN, terrain-masked lines RED.",
            )
            .clicked()
        {
            s.compute_los();
        }
        if ui
            .button("Observer = route start")
            .on_hover_text("Move the observer cell to the current routing start cell.")
            .clicked()
        {
            s.los_observer = s.route_start;
        }
    });
    // LoS status readout (visible / blocked counts).
    ui.label(
        egui::RichText::new(&s.los_status)
            .strong()
            .color(egui::Color32::from_rgb(120, 220, 170)),
    );

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

/// The geographic bounding box `(min_lat, max_lat, min_lon, max_lon)` the routing
/// cost grid is laid over: the scenario extent, padded ~20% so the grid frames
/// the action, with a fixed fallback (a ~1°×1° Bay-Area box) when the scenario
/// is empty. The route grid's cell `(0,0)` maps to the box's NW corner.
fn route_bbox(sc: &PlannerScenario) -> (f64, f64, f64, f64) {
    match scenario_extent(sc) {
        Some((min_lat, max_lat, min_lon, max_lon)) => {
            let pad_lat = ((max_lat - min_lat).abs() * 0.2).max(0.05);
            let pad_lon = ((max_lon - min_lon).abs() * 0.2).max(0.05);
            (
                min_lat - pad_lat,
                max_lat + pad_lat,
                min_lon - pad_lon,
                max_lon + pad_lon,
            )
        }
        None => (37.0, 38.0, -122.5, -121.5),
    }
}

/// Map a routing grid cell `(x, y)` to a geographic `(lat, lon)` inside `bbox`
/// `(min_lat, max_lat, min_lon, max_lon)`. Cell centres are placed across the box
/// so the path overlays the basemap; `y = 0` is the **north** (max-lat) edge so
/// the grid reads top-down like the screen. Robust to a 1×1 grid (no divide by
/// zero).
fn cell_to_latlon(cell: (usize, usize), grid: &CostGrid, bbox: (f64, f64, f64, f64)) -> (f64, f64) {
    let (min_lat, max_lat, min_lon, max_lon) = bbox;
    let (x, y) = cell;
    let fx = if grid.w > 1 {
        x as f64 / (grid.w - 1) as f64
    } else {
        0.5
    };
    let fy = if grid.h > 1 {
        y as f64 / (grid.h - 1) as f64
    } else {
        0.5
    };
    let lon = min_lon + fx * (max_lon - min_lon);
    // y grows southward: y=0 -> max_lat (north), y=h-1 -> min_lat (south).
    let lat = max_lat - fy * (max_lat - min_lat);
    (lat, lon)
}

/// Inverse of [`cell_to_latlon`]: snap a geographic `(lat, lon)` to the nearest
/// routing grid cell `(x, y)` inside `bbox`, clamped to the grid bounds. Used to
/// place entities (which live at lat/lon) onto the cost grid for line-of-sight.
/// A degenerate (zero-span) box maps to the grid centre.
fn latlon_to_cell(
    lat: f64,
    lon: f64,
    grid: &CostGrid,
    bbox: (f64, f64, f64, f64),
) -> (usize, usize) {
    let (min_lat, max_lat, min_lon, max_lon) = bbox;
    let span_lon = max_lon - min_lon;
    let span_lat = max_lat - min_lat;
    let fx = if span_lon.abs() > f64::EPSILON {
        (lon - min_lon) / span_lon
    } else {
        0.5
    };
    // y=0 is the north (max_lat) edge, mirroring cell_to_latlon.
    let fy = if span_lat.abs() > f64::EPSILON {
        (max_lat - lat) / span_lat
    } else {
        0.5
    };
    let x = (fx * (grid.w.saturating_sub(1)) as f64).round();
    let y = (fy * (grid.h.saturating_sub(1)) as f64).round();
    // Clamp into bounds (entities outside the padded bbox snap to the edge).
    let x = (x.max(0.0) as usize).min(grid.w.saturating_sub(1));
    let y = (y.max(0.0) as usize).min(grid.h.saturating_sub(1));
    (x, y)
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
    let overlay = MissionOverlay::from_scenario_with_route(s);

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
    affiliation: Affiliation,
    route: Vec<(f64, f64)>,
}

/// Human-readable label for an APP-6 affiliation, used in the combo box and the
/// selected-text (the frame shape is named so the choice is unambiguous).
fn affiliation_label(a: Affiliation) -> &'static str {
    match a {
        Affiliation::Friendly => "Friendly (blue rounded rect)",
        Affiliation::Hostile => "Hostile (red diamond)",
        Affiliation::Neutral => "Neutral (green square)",
        Affiliation::Unknown => "Unknown (yellow quatrefoil)",
    }
}

/// The canonical APP-6 / MIL-STD-2525 **frame colour** for an affiliation
/// (friendly blue/cyan, hostile red, neutral green, unknown yellow).
fn affiliation_color(a: Affiliation) -> egui::Color32 {
    match a {
        Affiliation::Friendly => egui::Color32::from_rgb(80, 200, 255), // cyan/blue
        Affiliation::Hostile => egui::Color32::from_rgb(240, 80, 80),   // red
        Affiliation::Neutral => egui::Color32::from_rgb(110, 210, 130), // green
        Affiliation::Unknown => egui::Color32::from_rgb(245, 215, 70),  // yellow
    }
}

/// Paint the standard **APP-6 affiliation frame** for one entity at projected
/// screen position `c`, sized by `r` (half-extent in pixels): Friendly = blue
/// **rounded rectangle**, Hostile = red **diamond**, Neutral = green **square**,
/// Unknown = yellow **quatrefoil** (a rounded square at icon size). The frame is
/// stroked in the affiliation colour over a translucent dark fill (so it reads
/// over busy tiles), with a thin black halo for contrast — the canonical 2525
/// frame shape + colour. Map iconography only; no engagement semantics.
fn draw_app6_frame(painter: &egui::Painter, c: egui::Pos2, r: f32, a: Affiliation) {
    let col = affiliation_color(a);
    let fill = egui::Color32::from_black_alpha(150);
    let halo = egui::Stroke::new(2.6, egui::Color32::from_black_alpha(180));
    let edge = egui::Stroke::new(1.8, col);
    match a {
        // Friendly: rounded rectangle (wider than tall, per 2525 framing).
        Affiliation::Friendly => {
            let rect = egui::Rect::from_center_size(c, egui::vec2(r * 2.4, r * 1.8));
            let rounding = egui::Rounding::same(r * 0.6);
            painter.rect_filled(rect, rounding, fill);
            painter.rect_stroke(rect, rounding, halo);
            painter.rect_stroke(rect, rounding, edge);
        }
        // Hostile: diamond (square rotated 45°).
        Affiliation::Hostile => {
            let d = r * 1.45;
            let pts = vec![
                c + egui::vec2(0.0, -d),
                c + egui::vec2(d, 0.0),
                c + egui::vec2(0.0, d),
                c + egui::vec2(-d, 0.0),
            ];
            painter.add(egui::Shape::convex_polygon(pts.clone(), fill, halo));
            painter.add(egui::Shape::closed_line(pts, edge));
        }
        // Neutral: upright square.
        Affiliation::Neutral => {
            let rect = egui::Rect::from_center_size(c, egui::vec2(r * 2.0, r * 2.0));
            let sq = egui::Rounding::ZERO;
            painter.rect_filled(rect, sq, fill);
            painter.rect_stroke(rect, sq, halo);
            painter.rect_stroke(rect, sq, edge);
        }
        // Unknown: quatrefoil — approximated at icon size by a heavily rounded
        // square (the 2525 cloverleaf reads as a "blob"); good enough as a frame
        // and cheap to paint with the same stroke style as the others.
        Affiliation::Unknown => {
            let rect = egui::Rect::from_center_size(c, egui::vec2(r * 2.0, r * 2.0));
            let rounding = egui::Rounding::same(r * 1.0);
            painter.rect_filled(rect, rounding, fill);
            painter.rect_stroke(rect, rounding, halo);
            painter.rect_stroke(rect, rounding, edge);
        }
    }
}

/// Paint a **dashed** line segment from `a` to `b` with the given `stroke`,
/// alternating `dash` px of line with `gap` px of space. Used for **masked**
/// line-of-sight rays so a blocked sight line reads distinctly from a clear
/// (solid) one even for colour-blind users.
fn dashed_line(
    painter: &egui::Painter,
    a: egui::Pos2,
    b: egui::Pos2,
    stroke: egui::Stroke,
    dash: f32,
    gap: f32,
) {
    let delta = b - a;
    let len = delta.length();
    if len < f32::EPSILON {
        return;
    }
    let dir = delta / len;
    let period = (dash + gap).max(1.0);
    let mut t = 0.0;
    while t < len {
        let seg_end = (t + dash).min(len);
        painter.line_segment([a + dir * t, a + dir * seg_end], stroke);
        t += period;
    }
}

/// A [`walkers::Plugin`] that paints the mission scenario — every entity's route
/// polyline, its waypoint dots, and the entity marker + name label — on top of
/// the OSM tile basemap, projecting each lat/lon to screen pixels with walkers'
/// [`Projector`] so the overlay stays pinned to the map as it pans / zooms.
/// One terrain shading cell: its NW corner `(lat, lon)`, its SE corner
/// `(lat, lon)`, and its elevation normalized to `[0, 1]` across the field relief.
/// Drawn as a translucent hypsometric-tinted rect under the map overlay.
type TerrainShadeCell = ((f64, f64), (f64, f64), f32);

struct MissionOverlay {
    entities: Vec<OverlayEntity>,
    /// Terrain elevation shading: one [`TerrainShadeCell`] per grid cell. Drawn
    /// first (under everything) as translucent filled rects ramped green (low) →
    /// tan (mid) → brown/white (high), so the user sees the landscape the routing
    /// + LoS are reasoning over. Empty when terrain display is off.
    terrain_cells: Vec<TerrainShadeCell>,
    /// The computed A\* tactical route as geographic `(lat, lon)` points (empty
    /// when no route is computed) — drawn as a DISTINCT-coloured polyline,
    /// separate from the per-entity blue routes.
    tactical_route: Vec<(f64, f64)>,
    /// The geographic centres of the routing grid's **obstacle** cells, drawn as
    /// small hatched markers so the user sees what the route is avoiding.
    obstacles: Vec<(f64, f64)>,
    /// The line-of-sight observer position as `(lat, lon)`, present only when a
    /// LoS result has been computed (the origin of every sight line).
    los_observer: Option<(f64, f64)>,
    /// Each computed sight line as `(target_latlon, visible)`: a GREEN solid line
    /// when `visible`, a RED dashed/dim line when masked by terrain.
    los_lines: Vec<((f64, f64), bool)>,
}

impl MissionOverlay {
    /// Build the overlay from the scenario plus the tactical-routing state: the
    /// A\* path and the cost-grid obstacles are mapped from grid cells to lat/lon
    /// over the routing bbox so they pin to the basemap with the entity routes.
    fn from_scenario_with_route(s: &MissionPlannerWorkbenchState) -> Self {
        let sc = &s.scenario;
        let entities = sc
            .entities
            .iter()
            .map(|e| OverlayEntity {
                name: e.name.clone(),
                lat: e.lat,
                lon: e.lon,
                done: e.is_done(),
                affiliation: e.affiliation,
                route: e.route.iter().map(|wp| (wp.lat, wp.lon)).collect(),
            })
            .collect();

        let bbox = route_bbox(sc);

        // Terrain elevation shading: snapshot each cell as a lat/lon quad plus its
        // height normalized to [0,1] over the field relief. Shown whenever terrain
        // is on (the landscape that routing/LoS reason over). Half-cell offsets give
        // each cell's NW/SE corners so the rects tile the bbox without gaps.
        let mut terrain_cells = Vec::new();
        if s.terrain_on {
            let (elev_lo, elev_hi) = s.terrain.minmax();
            let span = (elev_hi - elev_lo).max(1e-3);
            let (min_lat, max_lat, min_lon, max_lon) = bbox;
            let (gw, gh) = (s.terrain.w, s.terrain.h);
            // Per-cell half-spans in degrees (so a cell rect spans one grid step).
            let half_lon = if gw > 1 {
                0.5 * (max_lon - min_lon) / (gw - 1) as f64
            } else {
                0.5 * (max_lon - min_lon)
            };
            let half_lat = if gh > 1 {
                0.5 * (max_lat - min_lat) / (gh - 1) as f64
            } else {
                0.5 * (max_lat - min_lat)
            };
            terrain_cells.reserve(gw * gh);
            for y in 0..gh {
                for x in 0..gw {
                    // route_grid and terrain share the same extent (cell-for-cell),
                    // so cell_to_latlon (keyed on grid dims) is valid for both.
                    let (clat, clon) = cell_to_latlon((x, y), &s.route_grid, bbox);
                    let nw = (clat + half_lat, clon - half_lon); // north / west
                    let se = (clat - half_lat, clon + half_lon); // south / east
                    let hnorm = ((s.terrain.elevation_at(x, y) - elev_lo) / span).clamp(0.0, 1.0);
                    terrain_cells.push((nw, se, hnorm));
                }
            }
        }

        let tactical_route = s
            .route
            .as_ref()
            .map(|path| {
                path.iter()
                    .map(|&cell| cell_to_latlon(cell, &s.route_grid, bbox))
                    .collect()
            })
            .unwrap_or_default();

        // Snapshot obstacle cell centres (only when a route exists, to keep the
        // idle map uncluttered) for context.
        let mut obstacles = Vec::new();
        if s.route.is_some() {
            for y in 0..s.route_grid.h {
                for x in 0..s.route_grid.w {
                    if s.route_grid.cost_at(x, y).is_infinite() {
                        obstacles.push(cell_to_latlon((x, y), &s.route_grid, bbox));
                    }
                }
            }
        }

        // Line-of-sight: project the observer + each target sight line to lat/lon
        // over the same routing bbox so they pin to the basemap. Present only
        // once `compute_los` has run (non-empty results).
        let (los_observer, los_lines) = if s.los_results.is_empty() {
            (None, Vec::new())
        } else {
            let obs = cell_to_latlon(s.los_observer, &s.route_grid, bbox);
            let lines = s
                .los_results
                .iter()
                .map(|&(cell, vis)| (cell_to_latlon(cell, &s.route_grid, bbox), vis))
                .collect();
            (Some(obs), lines)
        };

        Self {
            entities,
            terrain_cells,
            tactical_route,
            obstacles,
            los_observer,
            los_lines,
        }
    }
}

/// The elevation colour ramp for the terrain shade: `h` in `[0, 1]` maps low → a
/// translucent **green** lowland, mid → **tan**, high → **brown**, top → near
/// **white** (snow-capped), a conventional hypsometric tint. Alpha is kept low so
/// the basemap, routes, and markers stay legible over it.
fn elevation_color(h: f32) -> egui::Color32 {
    let h = h.clamp(0.0, 1.0);
    // Piecewise-linear over four stops: green, tan, brown, white.
    let stops = [
        (0.0_f32, (70.0, 130.0, 70.0)), // low: green
        (0.45, (170.0, 160.0, 95.0)),   // mid: tan
        (0.80, (130.0, 95.0, 70.0)),    // high: brown
        (1.0, (235.0, 235.0, 235.0)),   // peak: near-white
    ];
    let mut col = stops[stops.len() - 1].1;
    for w in stops.windows(2) {
        let (h0, c0) = w[0];
        let (h1, c1) = w[1];
        if h <= h1 {
            let t = if (h1 - h0).abs() > f32::EPSILON {
                (h - h0) / (h1 - h0)
            } else {
                0.0
            };
            col = (
                c0.0 + (c1.0 - c0.0) * t,
                c0.1 + (c1.1 - c0.1) * t,
                c0.2 + (c1.2 - c0.2) * t,
            );
            break;
        }
    }
    // Translucent so the OSM basemap and overlays read through the shade.
    egui::Color32::from_rgba_unmultiplied(col.0 as u8, col.1 as u8, col.2 as u8, 90)
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
        let done_col = egui::Color32::from_rgb(120, 200, 140);
        // The tactical A* route uses a DISTINCT amber/orange, deliberately unlike
        // the blue entity routes and every APP-6 affiliation colour.
        let tac_col = egui::Color32::from_rgb(255, 160, 30);
        let obstacle_col = egui::Color32::from_rgb(200, 60, 60);

        // TERRAIN elevation shade (drawn first, UNDER everything): one translucent
        // filled rect per grid cell, hypsometric-tinted by normalized height
        // (green lowland → tan → brown ridge → white peak). Pinned to the map via
        // the projector so it pans / zooms with the basemap.
        for &(nw, se, hnorm) in &self.terrain_cells {
            let p_nw = to_px(nw.0, nw.1);
            let p_se = to_px(se.0, se.1);
            let rect = egui::Rect::from_two_pos(p_nw, p_se);
            painter.rect_filled(rect, 0.0, elevation_color(hnorm));
        }

        // Cost-field OBSTACLES (drawn first, under everything) as small red X's.
        for &(lat, lon) in &self.obstacles {
            let p = to_px(lat, lon);
            let r = 2.2;
            let st = egui::Stroke::new(1.2, obstacle_col.gamma_multiply(0.85));
            painter.line_segment([p + egui::vec2(-r, -r), p + egui::vec2(r, r)], st);
            painter.line_segment([p + egui::vec2(-r, r), p + egui::vec2(r, -r)], st);
        }

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

        // The tactical A* route: a thick amber polyline with a dark halo so it
        // reads over the basemap, plus a green start dot and a red-ringed goal.
        if self.tactical_route.len() >= 2 {
            let pts: Vec<egui::Pos2> = self
                .tactical_route
                .iter()
                .map(|&(lat, lon)| to_px(lat, lon))
                .collect();
            for w in pts.windows(2) {
                // Halo first, then the bright route line on top.
                painter.line_segment(
                    [w[0], w[1]],
                    egui::Stroke::new(4.5, egui::Color32::from_black_alpha(170)),
                );
                painter.line_segment([w[0], w[1]], egui::Stroke::new(2.6, tac_col));
            }
        }
        // Start / goal markers (shown whenever a route has been computed).
        if let (Some(&start), Some(&goal)) =
            (self.tactical_route.first(), self.tactical_route.last())
        {
            let sp = to_px(start.0, start.1);
            let gp = to_px(goal.0, goal.1);
            painter.circle_filled(sp, 4.0, done_col);
            painter.circle_stroke(sp, 5.0, egui::Stroke::new(1.5, egui::Color32::BLACK));
            painter.circle_stroke(gp, 5.0, egui::Stroke::new(2.0, tac_col));
            painter.circle_filled(gp, 2.5, egui::Color32::from_rgb(240, 80, 80));
        }

        // Line-of-sight rays from the observer to each target: GREEN solid when
        // the sight line is clear, RED dashed/dim when terrain or an obstacle
        // MASKS it. Drawn under the entity markers; the observer is an eye-ringed
        // node. Defensive sensor-visibility analysis only.
        if let Some((olat, olon)) = self.los_observer {
            let op = to_px(olat, olon);
            let vis_col = egui::Color32::from_rgb(90, 230, 130); // clear → green
            let masked_col = egui::Color32::from_rgb(235, 90, 80); // masked → red
            for &((tlat, tlon), visible) in &self.los_lines {
                let tp = to_px(tlat, tlon);
                if visible {
                    // Solid bright green line for a clear sight line.
                    painter.line_segment([op, tp], egui::Stroke::new(2.0, vis_col));
                } else {
                    // Dashed dim red line for a masked sight line.
                    dashed_line(
                        &painter,
                        op,
                        tp,
                        egui::Stroke::new(1.6, masked_col),
                        6.0,
                        5.0,
                    );
                }
                // Small target pip in the line's colour.
                painter.circle_filled(tp, 2.6, if visible { vis_col } else { masked_col });
            }
            // Observer node: a white-ringed dot so the LoS origin is obvious.
            painter.circle_filled(op, 3.4, egui::Color32::from_rgb(250, 250, 250));
            painter.circle_stroke(op, 5.2, egui::Stroke::new(1.6, egui::Color32::BLACK));
            painter.circle_stroke(op, 6.6, egui::Stroke::new(1.4, vis_col));
        }

        // Entity markers: the standard APP-6 / MIL-STD-2525 affiliation FRAME
        // (shape + colour from the entity's affiliation) instead of a plain dot,
        // plus a small centre dot that turns green when the route is complete,
        // plus the entity name label with a halo for legibility over tiles.
        for e in &self.entities {
            let c = to_px(e.lat, e.lon);
            draw_app6_frame(&painter, c, 6.0, e.affiliation);
            // Centre status pip: green once arrived, else the frame colour.
            let pip = if e.done {
                done_col
            } else {
                affiliation_color(e.affiliation)
            };
            painter.circle_filled(c, 1.8, pip);
            // Halo behind the label for legibility over busy tiles.
            painter.text(
                c + egui::vec2(10.0, -10.0) + egui::vec2(1.0, 1.0),
                egui::Align2::LEFT_BOTTOM,
                &e.name,
                egui::FontId::monospace(11.0),
                egui::Color32::from_black_alpha(200),
            );
            painter.text(
                c + egui::vec2(10.0, -10.0),
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
        // The readout also surfaces the APP-6 affiliation counts.
        assert!(
            r.contains("affiliation F"),
            "readout should mention affiliation counts: {r}"
        );
    }

    #[test]
    fn agent_set_affiliation_all_sets_every_entity() {
        use crate::agent_commands::AgentValue;
        let mut s = MissionPlannerWorkbenchState::default();
        // The demo seeds a mix; force all to hostile via the bridge.
        s.agent_set("Affiliation (all)", &AgentValue::Str("hostile".into()))
            .expect("valid affiliation");
        assert!(
            s.scenario
                .entities
                .iter()
                .all(|e| e.affiliation == Affiliation::Hostile),
            "every entity should now be hostile"
        );
        let [fr, ho, ne, un] = s.affiliation_counts();
        assert_eq!((fr, ne, un), (0, 0, 0));
        assert_eq!(ho, s.scenario.entities.len());
    }

    #[test]
    fn agent_set_affiliation_rejects_bad_name_and_type() {
        use crate::agent_commands::AgentValue;
        let mut s = MissionPlannerWorkbenchState::default();
        // Unknown affiliation name is rejected (fail-loud, nothing written).
        let before: Vec<_> = s.scenario.entities.iter().map(|e| e.affiliation).collect();
        assert!(s
            .agent_set("Affiliation (all)", &AgentValue::Str("banana".into()))
            .is_err());
        // A non-string value is a type error.
        assert!(s
            .agent_set("Affiliation (all)", &AgentValue::Int(1))
            .is_err());
        let after: Vec<_> = s.scenario.entities.iter().map(|e| e.affiliation).collect();
        assert_eq!(before, after, "a rejected set must not mutate state");
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
        assert!(names.contains(&"Affiliation (all)"));
        // The tactical-routing controls are also discoverable.
        assert!(names.contains(&"Route start X"));
        assert!(names.contains(&"Route start Y"));
        assert!(names.contains(&"Route goal X"));
        assert!(names.contains(&"Route goal Y"));
    }

    #[test]
    fn compute_route_solves_the_default_demo_field() {
        // The default corner-to-corner problem over the demo field is solvable;
        // computing it stores a contiguous path and a km-bearing status.
        let mut s = MissionPlannerWorkbenchState::default();
        assert!(s.route.is_none());
        s.compute_route();
        let path = s.route.as_ref().expect("demo field is solvable");
        assert_eq!(path.first(), Some(&s.route_start));
        assert_eq!(path.last(), Some(&s.route_goal));
        assert!(s.route_status.contains("route:") && s.route_status.contains("km"));
        // The path must avoid every obstacle in the cost grid.
        for &(x, y) in path {
            assert!(s.route_grid.cost_at(x, y).is_finite());
        }
    }

    #[test]
    fn compute_route_reports_no_route_when_endpoint_is_blocked() {
        let mut s = MissionPlannerWorkbenchState::default();
        // Make the start cell an obstacle → unreachable, fail-loud status.
        s.route_grid.cost[0] = f32::INFINITY;
        s.route_start = (0, 0);
        s.compute_route();
        assert!(s.route.is_none());
        assert!(s.route_status.starts_with("no route"));
    }

    #[test]
    fn agent_set_route_coords_update_endpoints() {
        use crate::agent_commands::AgentValue;
        let mut s = MissionPlannerWorkbenchState::default();
        s.agent_set("Route start X", &AgentValue::Int(3)).unwrap();
        s.agent_set("Route start Y", &AgentValue::Int(4)).unwrap();
        s.agent_set("Route goal X", &AgentValue::Int(10)).unwrap();
        s.agent_set("Route goal Y", &AgentValue::Int(9)).unwrap();
        assert_eq!(s.route_start, (3, 4));
        assert_eq!(s.route_goal, (10, 9));
    }

    #[test]
    fn agent_set_route_coords_reject_out_of_range_and_type() {
        use crate::agent_commands::AgentValue;
        let mut s = MissionPlannerWorkbenchState::default();
        let before = s.route_goal;
        assert!(s.agent_set("Route goal X", &AgentValue::Int(9999)).is_err());
        assert!(s.agent_set("Route start Y", &AgentValue::Int(-1)).is_err());
        // A fractional float is a type error (an integral float like 2.0 is, per
        // the AgentValue::as_i64 convention, accepted as a coordinate).
        assert!(s
            .agent_set("Route start X", &AgentValue::Float(2.5))
            .is_err());
        assert_eq!(s.route_goal, before, "a rejected set must not mutate state");
    }

    #[test]
    fn route_bridge_helper_computes_route() {
        let mut app = ValenxApp::default();
        assert!(app.mission_planner.route.is_none());
        route(&mut app);
        assert!(
            app.mission_planner.route.is_some(),
            "the route() bridge helper computes a path"
        );
    }

    #[test]
    fn readout_includes_route_status() {
        let mut s = MissionPlannerWorkbenchState::default();
        s.compute_route();
        let r = s.agent_readout().expect("readout present");
        assert!(
            r.contains("route:"),
            "readout should surface the route status: {r}"
        );
    }

    #[test]
    fn los_control_names_are_listed() {
        let names = MissionPlannerWorkbenchState::agent_control_names();
        assert!(names.contains(&"Observer X"));
        assert!(names.contains(&"Observer Y"));
    }

    #[test]
    fn compute_los_classifies_targets_over_the_flat_field() {
        // With terrain OFF the LoS runs over the flat obstacle field (demo_field).
        // From the NW corner observer, computing LoS produces one result per
        // distinct target (goal + entities, minus the observer's own cell) and a
        // visible/blocked summary; each stored flag matches a direct flat
        // line_of_sight call, and the SE corner goal is masked by the demo walls.
        let mut s = MissionPlannerWorkbenchState::default();
        s.terrain_on = false;
        s.rebuild_route_grid();
        s.route_goal = (ROUTE_GRID_W - 1, ROUTE_GRID_H - 1); // SE corner (behind walls)
        assert!(s.los_results.is_empty());
        s.los_observer = (0, 0);
        s.compute_los();
        assert!(
            !s.los_results.is_empty(),
            "LoS computes at least the goal target"
        );
        assert!(s.los_status.contains("flat"));
        assert!(s.los_status.contains("visible") && s.los_status.contains("blocked"));
        // Every result's `visible` flag matches a direct (flat) line_of_sight call.
        for &(target, vis) in &s.los_results {
            assert_eq!(
                vis,
                line_of_sight(&s.route_grid, s.los_observer, target),
                "stored visibility must match the pure LoS for {target:?}"
            );
            // The observer's own cell is never a target.
            assert_ne!(target, s.los_observer);
        }
        // The corner-to-corner demo field walls SE-bound sight lines, so the goal
        // (SE corner) must be masked from the NW observer.
        let goal_vis = s
            .los_results
            .iter()
            .find(|(c, _)| *c == s.route_goal)
            .map(|(_, v)| *v);
        assert_eq!(
            goal_vis,
            Some(false),
            "the SE goal is behind the demo walls → masked from the NW observer"
        );
    }

    #[test]
    fn agent_set_observer_coords_update_and_reject() {
        use crate::agent_commands::AgentValue;
        let mut s = MissionPlannerWorkbenchState::default();
        s.agent_set("Observer X", &AgentValue::Int(5)).unwrap();
        s.agent_set("Observer Y", &AgentValue::Int(6)).unwrap();
        assert_eq!(s.los_observer, (5, 6));
        // Out-of-range and wrong-type are fail-loud, leaving state intact.
        let before = s.los_observer;
        assert!(s.agent_set("Observer X", &AgentValue::Int(9999)).is_err());
        assert!(s.agent_set("Observer Y", &AgentValue::Int(-1)).is_err());
        assert!(s.agent_set("Observer X", &AgentValue::Float(1.5)).is_err());
        assert_eq!(
            s.los_observer, before,
            "a rejected set must not mutate state"
        );
    }

    #[test]
    fn los_bridge_helper_computes_los() {
        let mut app = ValenxApp::default();
        assert!(app.mission_planner.los_results.is_empty());
        los(&mut app);
        assert!(
            !app.mission_planner.los_results.is_empty(),
            "the los() bridge helper computes visibility"
        );
    }

    #[test]
    fn readout_includes_los_status() {
        let mut s = MissionPlannerWorkbenchState::default();
        s.compute_los();
        let r = s.agent_readout().expect("readout present");
        assert!(
            r.contains("LoS") && r.contains("from"),
            "readout should surface the LoS status: {r}"
        );
    }

    #[test]
    fn terrain_controls_are_listed() {
        let names = MissionPlannerWorkbenchState::agent_control_names();
        assert!(names.contains(&"Terrain"));
        assert!(names.contains(&"Observer height (m)"));
        assert!(names.contains(&"Target height (m)"));
    }

    #[test]
    fn default_state_is_terrain_aware() {
        // By default terrain is on, the heightfield has real relief, and the cost
        // grid was derived from it (it contains impassable steep cells).
        let s = MissionPlannerWorkbenchState::default();
        assert!(s.terrain_on);
        let (lo, hi) = s.terrain.minmax();
        assert!(hi - lo > 100.0, "default terrain should have relief");
        assert_eq!((s.terrain.w, s.terrain.h), (ROUTE_GRID_W, ROUTE_GRID_H));
        assert!(
            s.route_grid.cost.iter().any(|c| c.is_infinite()),
            "the terrain-derived cost grid should mark steep ridge cells impassable"
        );
    }

    #[test]
    fn toggling_terrain_off_rebuilds_to_the_flat_field() {
        use crate::agent_commands::AgentValue;
        let mut s = MissionPlannerWorkbenchState::default();
        // Turn terrain OFF via the bridge → cost grid becomes the flat demo_field.
        s.agent_set("Terrain", &AgentValue::Bool(false)).unwrap();
        assert!(!s.terrain_on);
        assert_eq!(
            s.route_grid,
            demo_field(ROUTE_GRID_W, ROUTE_GRID_H),
            "terrain off should restore the flat obstacle field"
        );
        // Toggling back on re-derives from terrain (impassable steep cells return).
        s.agent_set("Terrain", &AgentValue::Bool(true)).unwrap();
        assert!(s.terrain_on);
        assert!(s.route_grid.cost.iter().any(|c| c.is_infinite()));
    }

    #[test]
    fn agent_set_observer_target_height_validate() {
        use crate::agent_commands::AgentValue;
        let mut s = MissionPlannerWorkbenchState::default();
        s.agent_set("Observer height (m)", &AgentValue::Float(35.0))
            .unwrap();
        s.agent_set("Target height (m)", &AgentValue::Float(2.0))
            .unwrap();
        assert!((s.obs_height_m - 35.0).abs() < 1e-6);
        assert!((s.tgt_height_m - 2.0).abs() < 1e-6);
        // Out-of-range is fail-loud, leaving state intact.
        let before = (s.obs_height_m, s.tgt_height_m);
        assert!(s
            .agent_set("Observer height (m)", &AgentValue::Float(-1.0))
            .is_err());
        assert!(s
            .agent_set("Target height (m)", &AgentValue::Float(99_999.0))
            .is_err());
        assert_eq!((s.obs_height_m, s.tgt_height_m), before);
    }

    #[test]
    fn terrain_masked_los_masks_dead_ground_and_a_tall_observer_sees_more() {
        // With terrain on, an observer on one side of the diagonal ridge cannot
        // see a target in the dead ground on the far side at ground level; raising
        // the observer high enough recovers at least as many sight lines (it sees
        // over the terrain). Uses an explicit cross-ridge observer/target so the
        // masking is deterministic regardless of the demo entity positions.
        let mut s = MissionPlannerWorkbenchState::default();
        s.los_observer = (1, ROUTE_GRID_H - 2); // NW side of the ridge
        let across = (ROUTE_GRID_W - 1, 1); // SE side, beyond the ridge
        s.obs_height_m = 2.0;
        s.tgt_height_m = 0.0;
        // Directly assert the cross-ridge sight line is dead ground at ground level.
        assert!(
            !line_of_sight_terrain(&s.terrain, s.los_observer, across, 2.0, 0.0),
            "a ground observer must be masked across the ridge (dead ground)"
        );
        // A tall mast lifts the sight line over the ridge crest → now visible.
        assert!(
            line_of_sight_terrain(&s.terrain, s.los_observer, across, 5000.0, 0.0),
            "a very tall observer should see over the ridge"
        );

        // Through the workbench API: the terrain-masked status is reported, and a
        // taller observer never sees fewer of the goal/entity targets.
        s.compute_los();
        assert!(s.los_status.contains("terrain-masked"));
        let visible_low = s.los_results.iter().filter(|(_, v)| *v).count();
        s.obs_height_m = 5000.0;
        s.compute_los();
        let visible_high = s.los_results.iter().filter(|(_, v)| *v).count();
        assert!(
            visible_high >= visible_low,
            "a taller observer must not see fewer targets ({visible_high} < {visible_low})"
        );
    }

    #[test]
    fn readout_mentions_terrain_elevation_range_and_awareness() {
        let s = MissionPlannerWorkbenchState::default();
        let r = s.agent_readout().expect("readout present");
        assert!(
            r.contains("terrain ON"),
            "readout should report terrain mode: {r}"
        );
        assert!(
            r.contains("elev") && r.contains("terrain-aware"),
            "readout should surface the elevation range + terrain-awareness: {r}"
        );
    }

    #[test]
    fn elevation_color_ramps_low_to_high() {
        // The hypsometric ramp must be translucent and vary across the range.
        let low = elevation_color(0.0);
        let mid = elevation_color(0.5);
        let high = elevation_color(1.0);
        assert!(low.a() < 255, "terrain shade must be translucent");
        assert_ne!(low, mid);
        assert_ne!(mid, high);
        // Out-of-range inputs clamp without panicking.
        let _ = elevation_color(-5.0);
        let _ = elevation_color(5.0);
    }

    #[test]
    fn latlon_to_cell_round_trips_cell_centres() {
        // latlon_to_cell is the inverse of cell_to_latlon at cell centres: every
        // grid cell maps to lat/lon and back to itself.
        let grid = demo_field(ROUTE_GRID_W, ROUTE_GRID_H);
        let bbox = (37.0, 38.0, -122.5, -121.5);
        for &cell in &[
            (0, 0),
            (ROUTE_GRID_W - 1, ROUTE_GRID_H - 1),
            (10, 7),
            (5, 0),
        ] {
            let (lat, lon) = cell_to_latlon(cell, &grid, bbox);
            assert_eq!(
                latlon_to_cell(lat, lon, &grid, bbox),
                cell,
                "round trip must recover cell {cell:?}"
            );
        }
    }

    #[test]
    fn cell_to_latlon_maps_corners_into_the_bbox() {
        let grid = demo_field(ROUTE_GRID_W, ROUTE_GRID_H);
        let bbox = (37.0, 38.0, -122.5, -121.5);
        // NW cell (0,0) -> (max_lat, min_lon); SE cell -> (min_lat, max_lon).
        let nw = cell_to_latlon((0, 0), &grid, bbox);
        let se = cell_to_latlon((ROUTE_GRID_W - 1, ROUTE_GRID_H - 1), &grid, bbox);
        assert!((nw.0 - 38.0).abs() < 1e-9 && (nw.1 - (-122.5)).abs() < 1e-9);
        assert!((se.0 - 37.0).abs() < 1e-9 && (se.1 - (-121.5)).abs() < 1e-9);
    }

    #[test]
    fn affiliation_color_is_distinct_per_side() {
        // Each APP-6 affiliation must map to a distinct frame colour.
        let cols: Vec<_> = Affiliation::ALL
            .iter()
            .map(|a| affiliation_color(*a))
            .collect();
        for i in 0..cols.len() {
            for j in (i + 1)..cols.len() {
                assert_ne!(cols[i], cols[j], "affiliation colours must be distinct");
            }
        }
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
    fn overlay_draws_without_panic_for_each_affiliation() {
        // Drawing the workbench paints the MissionOverlay (the APP-6 frames) for
        // whatever affiliation the entities carry. Force every entity to each of
        // the four affiliations in turn and confirm a full frame renders without
        // panicking (the headless run drives the map widget + overlay plugin).
        for a in Affiliation::ALL {
            let mut app = ValenxApp::default();
            app.show_mission_planner_workbench = true;
            for e in &mut app.mission_planner.scenario.entities {
                e.affiliation = a;
            }
            let nodes = draw_and_collect_nodes(&mut app);
            assert!(
                !nodes.is_empty(),
                "workbench with all-{a:?} entities must still render"
            );
        }
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
