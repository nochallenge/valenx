//! # `planner` — a geographic, waypoint-route **mission planner** (Stage 1).
//!
//! A constructive, real-time **visual mission planner** modelled after ArduPilot
//! Mission Planner / NASA GMAT: entities sit on a geographic map (latitude /
//! longitude) and follow ordered **waypoint routes**, advanced in real time so a
//! user can watch them move along their legs.
//!
//! ## Scope (Stage 1 only)
//!
//! This module is deliberately limited to **movement + routes**: a geographic
//! frame, entities, waypoint routes, and a pure [`PlannerScenario::step`] that
//! moves each entity toward its current waypoint. There is **no** engagement, no
//! sensor model, and no orbital mechanics here — those are later stages. The
//! posture is purely defensive / planning: entities move along routes, nothing
//! more.
//!
//! ## Geographic frame & distance model
//!
//! Positions are stored as latitude / longitude in **degrees**. For movement we
//! use a simple **equirectangular (planar) degrees↔metres** model rather than a
//! full great-circle integration, which is accurate over the short legs a
//! mission plan uses and keeps [`PlannerScenario::step`] trivially pure:
//!
//! - one degree of **latitude**  ≈ [`M_PER_DEG_LAT`] metres (`111_320 m`);
//! - one degree of **longitude** ≈ `M_PER_DEG_LAT · cos(latitude)` metres (the
//!   meridian convergence factor), evaluated at the entity's current latitude.
//!
//! Each [`PlannerScenario::step`] converts the entity→waypoint offset to metres
//! with that model, advances the entity along the straight line to the waypoint
//! by `speed_mps · dt_s`, snaps to the waypoint when it would overshoot, and then
//! advances the `leg` index. At the final waypoint the entity **holds**.
//!
//! For *display* a separate equirectangular **projection to screen** is provided
//! by [`project`] (`x = lon`, `y = -lat`, scaled), independent of the movement
//! model above.

/// Metres per degree of latitude (mean Earth radius, `111_320 m/deg`). Also the
/// per-degree scale of longitude at the equator; longitude shrinks by
/// `cos(latitude)` toward the poles.
pub const M_PER_DEG_LAT: f64 = 111_320.0;

/// A single geographic waypoint in latitude / longitude **degrees**.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Waypoint {
    /// Latitude in degrees (`+` north).
    pub lat: f64,
    /// Longitude in degrees (`+` east).
    pub lon: f64,
}

impl Waypoint {
    /// Construct a waypoint from latitude / longitude degrees.
    pub fn new(lat: f64, lon: f64) -> Self {
        Self { lat, lon }
    }
}

/// A moving entity that follows an ordered [`Waypoint`] route across the map.
///
/// `leg` is the index of the route waypoint the entity is currently moving
/// **toward**. When the entity reaches that waypoint, `leg` advances; once
/// `leg` reaches the end of `route` the entity **holds** at the final waypoint.
#[derive(Debug, Clone, PartialEq)]
pub struct Entity {
    /// Stable identifier (unique within a scenario).
    pub id: u32,
    /// Human-readable label shown on the map.
    pub name: String,
    /// Current latitude in degrees.
    pub lat: f64,
    /// Current longitude in degrees.
    pub lon: f64,
    /// Ground speed in metres / second. `0` (or negative) ⇒ the entity does not
    /// move (treated as stationary; never panics).
    pub speed_mps: f64,
    /// Ordered waypoint route the entity follows.
    pub route: Vec<Waypoint>,
    /// Index of the waypoint currently being moved toward (the current leg). At
    /// or beyond `route.len()` the entity is holding at the final waypoint.
    pub leg: usize,
}

impl Entity {
    /// Whether the entity has finished its route (holding at the last waypoint,
    /// i.e. there is no further waypoint to move toward).
    pub fn is_done(&self) -> bool {
        self.leg >= self.route.len()
    }

    /// The waypoint the entity is currently heading for, if any (`None` once the
    /// route is complete).
    pub fn current_target(&self) -> Option<Waypoint> {
        self.route.get(self.leg).copied()
    }

    /// Advance this entity toward its current waypoint by `dt_s` seconds.
    ///
    /// Pure kinematics under the equirectangular metres model documented at the
    /// module level. A non-positive `dt_s` or `speed_mps`, an empty route, or a
    /// completed route all leave the entity unchanged (no movement, no panic).
    /// When the remaining distance to the current waypoint is within one step,
    /// the entity snaps exactly onto the waypoint and `leg` advances by one.
    fn step(&mut self, dt_s: f64) {
        if dt_s <= 0.0 || self.speed_mps <= 0.0 {
            return;
        }
        let Some(target) = self.current_target() else {
            return; // route complete — hold.
        };

        // Offset to the target in metres (equirectangular about current lat).
        let cos_lat = (self.lat.to_radians()).cos();
        let east_m = (target.lon - self.lon) * M_PER_DEG_LAT * cos_lat;
        let north_m = (target.lat - self.lat) * M_PER_DEG_LAT;
        let dist_m = (east_m * east_m + north_m * north_m).sqrt();

        let travel_m = self.speed_mps * dt_s;
        if dist_m <= travel_m || dist_m == 0.0 {
            // Reached (or already at) the waypoint: snap on and advance the leg.
            self.lat = target.lat;
            self.lon = target.lon;
            self.leg += 1;
            return;
        }

        // Move a fraction of the way along the straight leg.
        let frac = travel_m / dist_m;
        let new_east_m = east_m * frac;
        let new_north_m = north_m * frac;
        self.lat += new_north_m / M_PER_DEG_LAT;
        // Guard the pole singularity where cos(lat) -> 0.
        if cos_lat.abs() > 1e-9 {
            self.lon += new_east_m / (M_PER_DEG_LAT * cos_lat);
        }
    }
}

/// A complete planner scenario: a set of route-following entities plus the
/// accumulated simulated time.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct PlannerScenario {
    /// The entities on the map.
    pub entities: Vec<Entity>,
    /// Accumulated simulated time in seconds.
    pub sim_time_s: f64,
}

impl PlannerScenario {
    /// Advance the whole scenario by `dt_s` seconds: bump `sim_time_s` and move
    /// every entity toward its current waypoint (see [`Entity::step`]).
    ///
    /// Pure — no I/O. A `dt_s` of `0` (or negative) is a no-op for movement; the
    /// simulated clock only advances for a positive `dt_s`, so a paused tick does
    /// not drift the time readout.
    pub fn step(&mut self, dt_s: f64) {
        if dt_s <= 0.0 {
            return;
        }
        self.sim_time_s += dt_s;
        for e in &mut self.entities {
            e.step(dt_s);
        }
    }

    /// Whether every entity has finished its route (all holding at their final
    /// waypoints). `true` for an empty scenario.
    pub fn all_done(&self) -> bool {
        self.entities.iter().all(Entity::is_done)
    }

    /// Seed a demonstration scenario with `n` entities, each on a simple
    /// multi-waypoint route over a small geographic region (around the Bay Area,
    /// ~37–38° N, −122° E). Entity `i` is offset north so the routes are visually
    /// separated; every entity starts at the first waypoint with `leg = 0`.
    ///
    /// `n` is clamped to at least 1 so the demo is never empty.
    pub fn demo(n: usize) -> Self {
        let n = n.max(1);
        let mut entities = Vec::with_capacity(n);
        for i in 0..n {
            // Stagger each entity's lane northwards by 0.10° (~11 km).
            let base_lat = 37.40 + i as f64 * 0.10;
            let base_lon = -122.20;
            // A simple 3-leg dog-leg route heading roughly east then north-east.
            let route = vec![
                Waypoint::new(base_lat, base_lon),
                Waypoint::new(base_lat + 0.05, base_lon + 0.25),
                Waypoint::new(base_lat + 0.20, base_lon + 0.40),
                Waypoint::new(base_lat + 0.30, base_lon + 0.70),
            ];
            // Vary speed a little per entity so they don't move in lockstep
            // (140–200 m/s — small-aircraft scale).
            let speed_mps = 140.0 + (i % 4) as f64 * 20.0;
            entities.push(Entity {
                id: i as u32,
                name: format!("E{}", i + 1),
                lat: route[0].lat,
                lon: route[0].lon,
                speed_mps,
                route,
                leg: 0, // heading toward route[0]; the first step snaps + advances.
            });
        }
        Self {
            entities,
            sim_time_s: 0.0,
        }
    }
}

/// Project a geographic latitude / longitude (degrees) to **screen** pixels with
/// a plain equirectangular mapping: `x = lon`, `y = -lat` (screen `+y` is down,
/// so northern latitudes map upward), each scaled by `scale` pixels-per-degree
/// and offset by `origin` (the screen pixel of `lat = lon = 0`).
///
/// Returns `(x, y)` in pixels. This is the **display** projection only and is
/// independent of the metres-based movement model used by [`PlannerScenario::step`].
pub fn project(lat: f64, lon: f64, origin: (f32, f32), scale: f32) -> (f32, f32) {
    let x = origin.0 + lon as f32 * scale;
    let y = origin.1 - lat as f32 * scale;
    (x, y)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Planar distance in metres between two lat/lon points under the same
    /// equirectangular model `step` uses (about the first point's latitude).
    fn dist_m(a: Waypoint, b: Waypoint) -> f64 {
        let cos_lat = a.lat.to_radians().cos();
        let east = (b.lon - a.lon) * M_PER_DEG_LAT * cos_lat;
        let north = (b.lat - a.lat) * M_PER_DEG_LAT;
        (east * east + north * north).sqrt()
    }

    fn two_wp_entity(speed_mps: f64) -> Entity {
        let wp1 = Waypoint::new(40.0, -100.0);
        let wp2 = Waypoint::new(40.0, -99.0); // due east, ~85 km at 40°N.
        Entity {
            id: 0,
            name: "T".to_string(),
            lat: wp1.lat,
            lon: wp1.lon,
            speed_mps,
            route: vec![wp1, wp2],
            leg: 1, // already sitting on wp1 -> head straight to wp2.
        }
    }

    #[test]
    fn reaches_second_waypoint_after_distance_over_speed() {
        let speed = 200.0;
        let mut sc = PlannerScenario {
            entities: vec![two_wp_entity(speed)],
            sim_time_s: 0.0,
        };
        let wp1 = sc.entities[0].route[0];
        let wp2 = sc.entities[0].route[1];
        let expected_t = dist_m(wp1, wp2) / speed;

        // Integrate in small steps and record the time at which the entity
        // first completes its route.
        let dt = 1.0;
        let mut t = 0.0;
        let mut arrival_t = None;
        for _ in 0..1000 {
            sc.step(dt);
            t += dt;
            if arrival_t.is_none() && sc.entities[0].is_done() {
                arrival_t = Some(t);
                break;
            }
        }
        let arrival_t = arrival_t.expect("entity should complete its route");
        let e = &sc.entities[0];
        assert!(e.is_done(), "entity should have completed its route");
        assert!(
            (e.lat - wp2.lat).abs() < 1e-6 && (e.lon - wp2.lon).abs() < 1e-6,
            "entity should be at wp2, got ({}, {})",
            e.lat,
            e.lon
        );
        assert!(
            (sc.sim_time_s - t).abs() < 1e-9,
            "sim_time should track the integrated time"
        );
        // Arrival time is ~ distance/speed (within one step of dt: the final
        // partial step snaps onto the waypoint).
        assert!(
            arrival_t >= expected_t && arrival_t <= expected_t + dt,
            "arrival ~ distance/speed: expected ~{expected_t:.1}s, arrived at {arrival_t:.1}s"
        );
    }

    #[test]
    fn position_is_monotone_along_the_leg() {
        let mut sc = PlannerScenario {
            entities: vec![two_wp_entity(200.0)],
            sim_time_s: 0.0,
        };
        // Moving due east (increasing lon); lon must be non-decreasing each step
        // until the waypoint is reached.
        let mut prev_lon = sc.entities[0].lon;
        for _ in 0..50 {
            sc.step(1.0);
            let lon = sc.entities[0].lon;
            assert!(
                lon + 1e-12 >= prev_lon,
                "longitude must advance monotonically toward the eastern waypoint: {lon} < {prev_lon}"
            );
            prev_lon = lon;
            if sc.entities[0].is_done() {
                break;
            }
        }
    }

    #[test]
    fn holds_at_the_end() {
        let mut sc = PlannerScenario {
            entities: vec![two_wp_entity(1_000.0)],
            sim_time_s: 0.0,
        };
        // Overshoot the route entirely.
        for _ in 0..200 {
            sc.step(10.0);
        }
        let wp2 = sc.entities[0].route[1];
        let (lat_after, lon_after) = (sc.entities[0].lat, sc.entities[0].lon);
        assert!(sc.entities[0].is_done());
        // Stepping further must NOT move the held entity.
        sc.step(10.0);
        assert_eq!(
            (sc.entities[0].lat, sc.entities[0].lon),
            (lat_after, lon_after),
            "a completed entity must hold its position"
        );
        assert!((sc.entities[0].lat - wp2.lat).abs() < 1e-6);
        assert!((sc.entities[0].lon - wp2.lon).abs() < 1e-6);
    }

    #[test]
    fn step_with_zero_dt_is_a_noop() {
        let mut sc = PlannerScenario::demo(3);
        let before = sc.clone();
        sc.step(0.0);
        assert_eq!(
            sc, before,
            "dt = 0 must change nothing (positions or clock)"
        );
        sc.step(-5.0);
        assert_eq!(sc, before, "negative dt must change nothing");
    }

    #[test]
    fn zero_or_negative_speed_does_not_move_or_panic() {
        for speed in [0.0, -50.0] {
            let mut sc = PlannerScenario {
                entities: vec![two_wp_entity(speed)],
                sim_time_s: 0.0,
            };
            let start = (sc.entities[0].lat, sc.entities[0].lon);
            sc.step(10.0);
            assert_eq!(
                (sc.entities[0].lat, sc.entities[0].lon),
                start,
                "a non-positive speed ({speed}) entity must not move"
            );
            // The clock still advances on a positive dt even if nothing moves.
            assert!((sc.sim_time_s - 10.0).abs() < 1e-9);
        }
    }

    #[test]
    fn demo_seeds_n_entities_with_routes() {
        let sc = PlannerScenario::demo(5);
        assert_eq!(sc.entities.len(), 5);
        assert!(sc.entities.iter().all(|e| e.route.len() >= 2));
        // demo(0) is clamped to one entity, never empty.
        assert_eq!(PlannerScenario::demo(0).entities.len(), 1);
    }

    #[test]
    fn projection_maps_north_up_and_east_right() {
        let origin = (100.0, 100.0);
        let scale = 4.0;
        let (x0, y0) = project(0.0, 0.0, origin, scale);
        assert_eq!((x0, y0), origin);
        // Higher latitude -> smaller screen y (up).
        let (_, y_north) = project(10.0, 0.0, origin, scale);
        assert!(y_north < y0, "north should map upward (smaller y)");
        // Higher longitude -> larger screen x (right).
        let (x_east, _) = project(0.0, 10.0, origin, scale);
        assert!(x_east > x0, "east should map rightward (larger x)");
    }
}
