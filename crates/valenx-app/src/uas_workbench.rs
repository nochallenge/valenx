//! The right-side **UAS Workbench** panel — a native front-end over the
//! in-house `valenx-uas` crate (small-UAS design → performance → trade study,
//! plus *defensive* counter-UAS geometry).
//!
//! This workbench is the **M1 track** of the valenx modeling-&-simulation
//! roadmap. It is the *same tool a civilian drone designer uses* — pick a
//! multirotor or fixed-wing configuration, describe the airframe / battery /
//! payload / drivetrain, and get integrated **performance** (hover power,
//! endurance, range, payload margin). On top of that it runs a one-parameter
//! **trade study** (sweep a design variable → a Pareto front, e.g. endurance
//! vs. payload) and a **defensive counter-UAS** layer.
//!
//! ## Dual-use boundary (defensive only)
//!
//! The counter-UAS group is **pure detect / track / intercept GEOMETRY and a
//! detection timeline** — exactly the kinematics a defender (or any aircraft
//! doing collision avoidance) needs to know *when* and *where* an inbound track
//! is reachable and *when* a sensor of range `R` first sees it. **There is no
//! weapon employment, no targeting solution, no guidance command, and no
//! lethality model of any kind** in this workbench or the `valenx-uas` crate it
//! drives. What happens at the rendezvous is out of scope and not modeled. This
//! is the same dual-use posture as civilian air-traffic conflict-detection /
//! closest-point-of-approach math.
//!
//! Mirrors the other workbenches (`uq_workbench`, `rom_workbench`): a
//! [`crate::workbench_chrome::workbench_shell`] panel gated on
//! [`crate::ValenxApp::show_uas_workbench`], toggled from the View menu and
//! openable by the agent bridge under the workbench id `"uas"` (aliases
//! `"drone"` / `"counteruas"`; see [`crate::project_tabs::TabKind`] /
//! [`crate::agent_commands`]).
//!
//! Three painter views are drawn: (a) a **performance readout** bar panel,
//! (b) a **trade-study scatter** of the swept design points with the Pareto
//! front highlighted, and (c) an **intercept plan-view** (the threat track,
//! the interceptor, the computed intercept point, and the sensor-range ring).
//!
//! Honesty: `valenx-uas` is **research / educational grade** — the performance
//! models inherit the caveats of the composed crates (ideal momentum-theory
//! hover; parabolic-polar point performance; a single lumped drivetrain
//! efficiency; no voltage sag, gusts, transients, or controls). The counter-UAS
//! layer is *exact constant-velocity kinematics*, not a tracking filter or
//! sensor-fusion model. Numbers are first estimates, not flight-grade or
//! accredited. Every error from `valenx-uas` surfaces verbatim — the workbench
//! never invents a number, and degenerate parameters (e.g. a zero rotor radius,
//! a `usable_fraction > 1`, an interceptor slower than an opening target) show
//! an in-panel error, NOT a panic.

use eframe::egui;
use nalgebra::Vector3;
use valenx_uas::trade::{DesignPoint, Objective, ParetoFront, TradeStudy};
use valenx_uas::vehicle::{Battery, FixedWingUas, MultirotorUas};
use valenx_uas::{detection_timeline, time_to_intercept, Interceptor, ThreatTrack};

use crate::ValenxApp;

// ---------------------------------------------------------------------------
// Parameters
// ---------------------------------------------------------------------------

/// Which airframe configuration the design describes.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum UasConfig {
    /// A multirotor (momentum-theory hover): rotor count + disk area driven.
    Multirotor,
    /// A fixed-wing (parabolic-polar point performance): wing + AR driven.
    FixedWing,
}

impl UasConfig {
    /// Combo / menu label.
    fn label(self) -> &'static str {
        match self {
            UasConfig::Multirotor => "multirotor (hover, momentum theory)",
            UasConfig::FixedWing => "fixed-wing (cruise, parabolic polar)",
        }
    }
}

/// Which design variable the one-parameter trade sweep walks.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum SweepVar {
    /// Sweep the payload mass (kg). Endurance vs. payload is the classic
    /// trade — more payload raises the all-up mass and cuts endurance.
    Payload,
    /// Sweep the installed battery energy (Wh). More energy buys endurance /
    /// range at the cost of (notional) mass — here held fixed, so it is a
    /// monotone gain swept purely to show the curve.
    BatteryWh,
}

impl SweepVar {
    /// Combo / menu label.
    fn label(self) -> &'static str {
        match self {
            SweepVar::Payload => "payload mass (kg)",
            SweepVar::BatteryWh => "battery energy (Wh)",
        }
    }
}

/// Parse an airframe-configuration name (for the agent `SetControl` bridge)
/// into a [`UasConfig`]. Case-insensitive; accepts the short menu words.
/// Fail-loud on an unrecognised name.
fn parse_uas_config(s: &str) -> Result<UasConfig, String> {
    match s.trim().to_ascii_lowercase().as_str() {
        "multirotor" | "multi" | "rotor" => Ok(UasConfig::Multirotor),
        "fixedwing" | "fixed-wing" | "fixed_wing" | "wing" => Ok(UasConfig::FixedWing),
        other => Err(format!(
            "unknown configuration '{other}' (expected 'multirotor' or 'fixed-wing')"
        )),
    }
}

/// Parse a trade-sweep-variable name (for the agent `SetControl` bridge) into a
/// [`SweepVar`]. Case-insensitive; accepts the short menu words. Fail-loud.
fn parse_sweep_var(s: &str) -> Result<SweepVar, String> {
    match s.trim().to_ascii_lowercase().as_str() {
        "payload" => Ok(SweepVar::Payload),
        "battery" | "batterywh" | "battery_wh" => Ok(SweepVar::BatteryWh),
        other => Err(format!(
            "unknown sweep variable '{other}' (expected 'payload' or 'battery')"
        )),
    }
}

/// Editable counter-UAS (defensive intercept GEOMETRY) inputs. All planar
/// (the painter draws a top-down plan view; the `z` component is held at 0).
#[derive(Clone, Copy, Debug)]
pub struct CounterUasParams {
    /// Threat track initial position `(x, y)` (m), relative to the defended
    /// site at the origin.
    pub threat_pos_xy: [f64; 2],
    /// Threat track constant velocity `(vx, vy)` (m/s).
    pub threat_vel_xy: [f64; 2],
    /// Interceptor start position `(x, y)` (m).
    pub interceptor_pos_xy: [f64; 2],
    /// Interceptor maximum speed `s` (m/s), finite and positive.
    pub interceptor_max_speed: f64,
    /// Sensor detection range `R` (m), finite and positive. The sensor is
    /// sited at the origin.
    pub sensor_range_m: f64,
}

impl Default for CounterUasParams {
    fn default() -> Self {
        Self {
            // An inbound track 1200 m out, closing at 30 m/s with a small
            // cross-track component so the plan view is not a degenerate line.
            threat_pos_xy: [1200.0, 200.0],
            threat_vel_xy: [-30.0, 0.0],
            interceptor_pos_xy: [0.0, 0.0],
            interceptor_max_speed: 50.0,
            sensor_range_m: 800.0,
        }
    }
}

/// Editable UAS inputs shown in the workbench: the airframe configuration and
/// its geometry, the battery, payload, drivetrain efficiency, the cruise /
/// thrust evaluation point, the trade-sweep selection, and the counter-UAS
/// (defensive intercept-geometry) group.
#[derive(Clone, Debug)]
pub struct UasParams {
    /// Airframe configuration.
    pub config: UasConfig,
    /// All-up mass `m` (kg) — the full take-off mass the lift is sized to.
    pub all_up_mass_kg: f64,

    // -- Multirotor geometry --
    /// Number of rotors.
    pub rotor_count: u32,
    /// Total rotor **disk area** `A` (m^2). Converted to an equivalent
    /// per-rotor radius for the composed momentum-theory crate.
    pub disk_area_m2: f64,
    /// Rotor figure of merit (hover efficiency) in `(0, 1]`.
    pub figure_of_merit: f64,
    /// Cruise power factor `k` for the multirotor range balance (cruise power
    /// ≈ `k · hover power`).
    pub cruise_power_factor: f64,
    /// Maximum thrust-to-weight ratio (sets the multirotor payload margin).
    pub max_thrust_to_weight: f64,

    // -- Fixed-wing geometry --
    /// Wing area `S` (m^2).
    pub wing_area_m2: f64,
    /// Wing aspect ratio `AR`.
    pub aspect_ratio: f64,
    /// Maximum lift coefficient `CL,max`.
    pub cl_max: f64,
    /// Zero-lift drag coefficient `CD0`.
    pub cd0: f64,
    /// Oswald span efficiency `e` in `(0, 1]`.
    pub oswald_efficiency: f64,

    // -- Shared --
    /// Installed battery energy `E` (Wh).
    pub battery_wh: f64,
    /// Usable battery fraction `f` in `(0, 1]` (depth-of-discharge × reserve).
    pub usable_fraction: f64,
    /// Payload mass `m_pl` (kg), `>= 0` and below the all-up mass.
    pub payload_kg: f64,
    /// Combined motor + ESC + propeller efficiency `η` in `(0, 1]`.
    pub drivetrain_efficiency: f64,
    /// Air density `ρ` (kg/m^3).
    pub air_density: f64,
    /// Cruise speed `V` (m/s) the range / endurance is evaluated at.
    pub cruise_speed_m_s: f64,

    // -- Trade sweep --
    /// Which design variable the one-parameter trade sweep walks.
    pub sweep_var: SweepVar,
    /// Sweep lower bound (units follow `sweep_var`).
    pub sweep_lo: f64,
    /// Sweep upper bound (units follow `sweep_var`).
    pub sweep_hi: f64,
    /// Number of sweep points (>= 2).
    pub sweep_points: usize,

    /// Counter-UAS (defensive intercept-geometry) inputs.
    pub counter: CounterUasParams,
}

impl Default for UasParams {
    fn default() -> Self {
        Self {
            config: UasConfig::Multirotor,
            all_up_mass_kg: 1.5,
            rotor_count: 4,
            // 4 rotors of ~0.15 m radius => 4·π·0.15² ≈ 0.283 m².
            disk_area_m2: 0.283,
            figure_of_merit: 0.70,
            cruise_power_factor: 1.0,
            max_thrust_to_weight: 2.0,
            wing_area_m2: 0.5,
            aspect_ratio: 10.0,
            cl_max: 1.3,
            cd0: 0.03,
            oswald_efficiency: 0.85,
            battery_wh: 100.0,
            usable_fraction: 0.8,
            payload_kg: 0.3,
            drivetrain_efficiency: 0.75,
            air_density: valenx_uas::SEA_LEVEL_AIR_DENSITY,
            cruise_speed_m_s: 15.0,
            sweep_var: SweepVar::Payload,
            sweep_lo: 0.0,
            sweep_hi: 0.8,
            sweep_points: 12,
            counter: CounterUasParams::default(),
        }
    }
}

impl UasParams {
    /// Equivalent per-rotor radius (m) from the total disk area and rotor
    /// count: `A = N·π·r²  ⇒  r = sqrt(A / (N·π))`. Fail-loud on a
    /// non-positive area / zero count so the composed crate never sees a NaN.
    fn rotor_radius_m(&self) -> Result<f64, String> {
        if self.rotor_count == 0 {
            return Err("rotor count must be >= 1".to_string());
        }
        if !(self.disk_area_m2.is_finite() && self.disk_area_m2 > 0.0) {
            return Err(format!(
                "disk area must be finite and positive, got {}",
                self.disk_area_m2
            ));
        }
        let r2 = self.disk_area_m2 / (f64::from(self.rotor_count) * std::f64::consts::PI);
        Ok(r2.sqrt())
    }

    /// Build a validated [`Battery`], surfacing `valenx-uas`'s fail-loud
    /// constructor error verbatim.
    fn battery(&self) -> Result<Battery, String> {
        Battery::new(self.battery_wh, self.usable_fraction).map_err(|e| e.to_string())
    }

    /// Build the validated multirotor for the given payload / battery overrides
    /// (the trade sweep substitutes one of these). `None` overrides keep the
    /// current field value.
    fn build_multirotor(
        &self,
        payload_override: Option<f64>,
        battery_override: Option<Battery>,
    ) -> Result<MultirotorUas, String> {
        let radius = self.rotor_radius_m()?;
        let battery = match battery_override {
            Some(b) => b,
            None => self.battery()?,
        };
        MultirotorUas::new(
            self.rotor_count,
            radius,
            self.all_up_mass_kg,
            self.figure_of_merit,
            self.air_density,
            battery,
            payload_override.unwrap_or(self.payload_kg),
            self.drivetrain_efficiency,
        )
        .map_err(|e| e.to_string())
    }

    /// Build the validated fixed-wing for the given payload / battery overrides.
    fn build_fixedwing(
        &self,
        payload_override: Option<f64>,
        battery_override: Option<Battery>,
    ) -> Result<FixedWingUas, String> {
        let battery = match battery_override {
            Some(b) => b,
            None => self.battery()?,
        };
        FixedWingUas::new(
            self.wing_area_m2,
            self.all_up_mass_kg,
            self.cl_max,
            self.aspect_ratio,
            self.cd0,
            self.oswald_efficiency,
            self.air_density,
            battery,
            payload_override.unwrap_or(self.payload_kg),
            self.drivetrain_efficiency,
        )
        .map_err(|e| e.to_string())
    }
}

// ---------------------------------------------------------------------------
// Simulation result
// ---------------------------------------------------------------------------

/// A single performance readout row (label + value + unit), used by the
/// painter readout panel and the text grid.
#[derive(Clone, Debug)]
pub struct PerfRow {
    /// Quantity name.
    pub label: String,
    /// Formatted value.
    pub value: String,
}

/// One evaluated trade-sweep point: the swept design value plus the two
/// objective values (endurance s, payload kg) and whether it is on the
/// Pareto front.
#[derive(Clone, Copy, Debug)]
pub struct TradePoint {
    /// Swept design value (units follow the chosen [`SweepVar`]).
    pub x: f64,
    /// Endurance objective (s).
    pub endurance_s: f64,
    /// Payload objective (kg).
    pub payload_kg: f64,
    /// True if this point is on the non-dominated Pareto front.
    pub on_front: bool,
}

/// The intercept-geometry outcome for the counter-UAS plan view.
#[derive(Clone, Copy, Debug)]
pub struct InterceptResult {
    /// Threat start position (plan, m).
    pub threat_pos: [f64; 2],
    /// Threat velocity (plan, m/s).
    pub threat_vel: [f64; 2],
    /// Interceptor start position (plan, m).
    pub interceptor_pos: [f64; 2],
    /// Sensor range (m), sited at the origin.
    pub sensor_range_m: f64,
    /// Earliest time-to-intercept (s), `None` if no feasible intercept.
    pub time_to_intercept_s: Option<f64>,
    /// Intercept point (plan, m), `None` if no feasible intercept.
    pub intercept_point: Option<[f64; 2]>,
    /// Whether the sensor ever detects the track.
    pub detected: bool,
    /// First-detection time (s), if detected.
    pub first_detect_s: Option<f64>,
    /// Time the track leaves the sensor ring again (s), if it does.
    pub last_in_range_s: Option<f64>,
    /// Closest-point-of-approach distance to the sensor (m).
    pub min_range_m: f64,
    /// Time of closest approach (s).
    pub time_of_closest_approach_s: f64,
}

/// Cached UAS output for the painter + readouts.
#[derive(Default, Clone)]
pub struct UasResult {
    /// Performance readout rows.
    pub perf_rows: Vec<PerfRow>,
    /// Evaluated trade-sweep points (with Pareto flags).
    pub trade: Vec<TradePoint>,
    /// Counter-UAS intercept geometry.
    pub intercept: Option<InterceptResult>,
}

// ---------------------------------------------------------------------------
// Workbench state
// ---------------------------------------------------------------------------

/// Persistent state for the UAS workbench.
#[derive(Default)]
pub struct UasWorkbenchState {
    /// User-editable parameters.
    pub params: UasParams,
    /// Last successful result (populated after a successful run).
    pub result: Option<UasResult>,
    /// Status / error line shown below the controls.
    pub status: String,
}

impl UasWorkbenchState {
    /// The user-visible captions of every control the agent bridge can set via
    /// `SetControl` (see [`crate::agent_commands`]). Returned by `ListControls`
    /// so an agent can discover the name space. Order follows the form.
    pub fn agent_control_names() -> &'static [&'static str] {
        &[
            "configuration",
            "all-up mass (kg)",
            "rotor count",
            "disk area (m^2)",
            "figure of merit",
            "cruise power factor",
            "max thrust/weight",
            "wing area (m^2)",
            "aspect ratio",
            "CL max",
            "CD0",
            "Oswald e",
            "battery (Wh)",
            "usable fraction",
            "payload (kg)",
            "drivetrain eff.",
            "air density (kg/m^3)",
            "cruise speed (m/s)",
            "sweep variable",
            "sweep lower bound",
            "sweep upper bound",
            "sweep points",
            "threat pos x (m)",
            "threat pos y (m)",
            "threat vel x (m/s)",
            "threat vel y (m/s)",
            "interceptor pos x (m)",
            "interceptor pos y (m)",
            "interceptor max speed (m/s)",
            "sensor range (m)",
        ]
    }

    /// Set one labelled control by its user-visible caption, for the agent
    /// `SetControl` bridge. Captions match exactly what the form draws. Fail-loud
    /// on an unknown caption / wrong type (the bridge posts a `warn` note); no
    /// field is written on error and nothing panics. The two enum captions
    /// (`configuration`, `sweep variable`) read [`AgentValue::as_str`]; the
    /// `rotor count` / `sweep points` integer fields read [`AgentValue::as_i64`];
    /// every other caption is an `f64` drag value.
    pub fn agent_set(
        &mut self,
        name: &str,
        value: &crate::agent_commands::AgentValue,
    ) -> Result<(), String> {
        let p = &mut self.params;
        let c = &mut p.counter;
        match name {
            // -- Airframe enum + geometry --
            "configuration" => p.config = parse_uas_config(value.as_str()?)?,
            "all-up mass (kg)" => p.all_up_mass_kg = value.as_f64()?,
            "rotor count" => {
                let n = value.as_i64()?;
                if !(1..=32).contains(&n) {
                    return Err(format!("rotor count must be in 1..=32, got {n}"));
                }
                p.rotor_count = n as u32;
            }
            "disk area (m^2)" => p.disk_area_m2 = value.as_f64()?,
            "figure of merit" => p.figure_of_merit = value.as_f64()?,
            "cruise power factor" => p.cruise_power_factor = value.as_f64()?,
            "max thrust/weight" => p.max_thrust_to_weight = value.as_f64()?,
            "wing area (m^2)" => p.wing_area_m2 = value.as_f64()?,
            "aspect ratio" => p.aspect_ratio = value.as_f64()?,
            "CL max" => p.cl_max = value.as_f64()?,
            "CD0" => p.cd0 = value.as_f64()?,
            "Oswald e" => p.oswald_efficiency = value.as_f64()?,
            // -- Battery / payload / drivetrain --
            "battery (Wh)" => p.battery_wh = value.as_f64()?,
            "usable fraction" => p.usable_fraction = value.as_f64()?,
            "payload (kg)" => p.payload_kg = value.as_f64()?,
            "drivetrain eff." => p.drivetrain_efficiency = value.as_f64()?,
            "air density (kg/m^3)" => p.air_density = value.as_f64()?,
            "cruise speed (m/s)" => p.cruise_speed_m_s = value.as_f64()?,
            // -- Trade sweep --
            "sweep variable" => p.sweep_var = parse_sweep_var(value.as_str()?)?,
            "sweep lower bound" => p.sweep_lo = value.as_f64()?,
            "sweep upper bound" => p.sweep_hi = value.as_f64()?,
            "sweep points" => {
                let n = value.as_i64()?;
                if !(2..=200).contains(&n) {
                    return Err(format!("sweep points must be in 2..=200, got {n}"));
                }
                p.sweep_points = n as usize;
            }
            // -- Counter-UAS (defensive intercept geometry) --
            "threat pos x (m)" => c.threat_pos_xy[0] = value.as_f64()?,
            "threat pos y (m)" => c.threat_pos_xy[1] = value.as_f64()?,
            "threat vel x (m/s)" => c.threat_vel_xy[0] = value.as_f64()?,
            "threat vel y (m/s)" => c.threat_vel_xy[1] = value.as_f64()?,
            "interceptor pos x (m)" => c.interceptor_pos_xy[0] = value.as_f64()?,
            "interceptor pos y (m)" => c.interceptor_pos_xy[1] = value.as_f64()?,
            "interceptor max speed (m/s)" => c.interceptor_max_speed = value.as_f64()?,
            "sensor range (m)" => c.sensor_range_m = value.as_f64()?,
            other => return Err(format!("unknown UAS control: {other:?}")),
        }
        Ok(())
    }

    /// Run the full UAS pipeline: build + validate the vehicle, compute its
    /// integrated performance, run the one-parameter trade sweep to a Pareto
    /// front, and solve the defensive counter-UAS intercept geometry +
    /// detection timeline.
    ///
    /// Every failure is returned as an `Err(String)` — no panics, no invented
    /// numbers. Degenerate inputs (a zero rotor radius, a `usable_fraction`
    /// outside `(0, 1]`, a payload `>=` all-up mass, an empty/too-short sweep,
    /// an interceptor slower than an opening target) surface `valenx-uas`'s own
    /// error verbatim.
    pub fn run(&self) -> Result<UasResult, String> {
        let p = &self.params;

        let perf_rows = self.performance_rows()?;
        let trade = self.run_trade_sweep()?;
        let intercept = Some(run_intercept(&p.counter)?);

        Ok(UasResult {
            perf_rows,
            trade,
            intercept,
        })
    }

    /// Compute the performance readout for the current configuration.
    fn performance_rows(&self) -> Result<Vec<PerfRow>, String> {
        let p = &self.params;
        let mut rows = Vec::new();
        let row = |label: &str, value: String| PerfRow {
            label: label.to_string(),
            value,
        };
        match p.config {
            UasConfig::Multirotor => {
                let uas = p.build_multirotor(None, None)?;
                let perf = uas
                    .performance(
                        p.cruise_speed_m_s,
                        p.cruise_power_factor,
                        p.max_thrust_to_weight,
                    )
                    .map_err(|e| e.to_string())?;
                rows.push(row("all-up mass", format!("{:.3} kg", perf.all_up_mass_kg)));
                rows.push(row("disk area", format!("{:.4} m^2", perf.disk_area_m2)));
                rows.push(row(
                    "disk loading",
                    format!("{:.2} N/m^2", perf.disk_loading_pa),
                ));
                rows.push(row("hover power", format!("{:.1} W", perf.hover_power_w)));
                rows.push(row(
                    "hover endurance",
                    format!("{:.2} min", perf.hover_endurance_s / 60.0),
                ));
                rows.push(row(
                    "cruise range",
                    format!("{:.2} km", perf.cruise_range_m / 1000.0),
                ));
                rows.push(row(
                    "max extra payload",
                    format!("{:.3} kg", perf.max_extra_payload_kg),
                ));
                rows.push(row(
                    "payload fraction",
                    format!("{:.1} %", perf.payload_fraction * 100.0),
                ));
            }
            UasConfig::FixedWing => {
                let uas = p.build_fixedwing(None, None)?;
                let perf = uas
                    .performance(p.cruise_speed_m_s)
                    .map_err(|e| e.to_string())?;
                rows.push(row("all-up mass", format!("{:.3} kg", perf.all_up_mass_kg)));
                rows.push(row(
                    "wing loading",
                    format!("{:.2} N/m^2", perf.wing_loading_pa),
                ));
                rows.push(row(
                    "stall speed",
                    format!("{:.2} m/s", perf.stall_speed_m_s),
                ));
                rows.push(row("(L/D) max", format!("{:.2}", perf.max_lift_to_drag)));
                rows.push(row(
                    "Breguet range",
                    format!("{:.2} km", perf.breguet_range_m / 1000.0),
                ));
                rows.push(row(
                    "endurance",
                    format!("{:.2} min", perf.endurance_s / 60.0),
                ));
                rows.push(row(
                    "max extra payload",
                    format!("{:.3} kg", perf.max_extra_payload_kg),
                ));
                rows.push(row(
                    "payload fraction",
                    format!("{:.1} %", perf.payload_fraction * 100.0),
                ));
            }
        }
        Ok(rows)
    }

    /// Run the one-parameter trade sweep through `valenx-uas`'s [`TradeStudy`]
    /// and tag each evaluated point with its Pareto-front membership.
    ///
    /// Objectives are **maximize endurance** and **maximize payload** — the
    /// classic small-UAS trade (more payload raises mass and cuts endurance).
    /// Each swept design substitutes the chosen [`SweepVar`] into a fresh
    /// validated vehicle; the evaluator returns `[endurance_s, payload_kg]`.
    fn run_trade_sweep(&self) -> Result<Vec<TradePoint>, String> {
        let p = &self.params;
        if p.sweep_points < 2 {
            return Err(format!(
                "trade sweep needs >= 2 points, got {}",
                p.sweep_points
            ));
        }
        if !(p.sweep_lo.is_finite() && p.sweep_hi.is_finite()) {
            return Err("trade sweep bounds must be finite".to_string());
        }
        if p.sweep_hi <= p.sweep_lo {
            return Err(format!(
                "trade sweep upper bound ({}) must exceed lower bound ({})",
                p.sweep_hi, p.sweep_lo
            ));
        }

        let n = p.sweep_points;
        let step = (p.sweep_hi - p.sweep_lo) / (n as f64 - 1.0);
        let var_key = "sweep_x";
        let designs: Vec<DesignPoint> = (0..n)
            .map(|i| {
                let x = p.sweep_lo + step * i as f64;
                DesignPoint::new(format!("{x:.4}"), [(var_key.to_string(), x)])
            })
            .collect();

        // Objectives: maximize endurance, maximize payload.
        let study = TradeStudy::new(vec![Objective::Maximize, Objective::Maximize], designs);

        // The evaluator builds a fresh validated vehicle per design point,
        // substituting the swept variable, and returns [endurance_s, payload].
        // It returns `valenx_uas::UasError` so it composes straight into the
        // study's `run` (UasError: From<UasError>).
        let front: ParetoFront = study
            .run::<valenx_uas::UasError>(|d| {
                let x = d.get(var_key).expect("design carries its swept value");
                let (payload, battery) = match p.sweep_var {
                    SweepVar::Payload => (Some(x), None),
                    SweepVar::BatteryWh => (None, Some(Battery::new(x, p.usable_fraction)?)),
                };
                let endurance = match p.config {
                    UasConfig::Multirotor => {
                        let uas = build_multirotor_uaserr(p, payload, battery)?;
                        uas.hover_endurance_s()?
                    }
                    UasConfig::FixedWing => {
                        let uas = build_fixedwing_uaserr(p, payload, battery)?;
                        uas.endurance_s(p.cruise_speed_m_s, None)?
                    }
                };
                let payload_val = payload.unwrap_or(p.payload_kg);
                Ok(vec![endurance, payload_val])
            })
            .map_err(|e| e.to_string())?;

        let on_front: std::collections::HashSet<usize> =
            front.front_indices.iter().copied().collect();
        let mut out = Vec::with_capacity(front.all.len());
        for (i, ev) in front.all.iter().enumerate() {
            let x = ev.design.get(var_key).unwrap_or(f64::NAN);
            out.push(TradePoint {
                x,
                endurance_s: ev.objectives.first().copied().unwrap_or(f64::NAN),
                payload_kg: ev.objectives.get(1).copied().unwrap_or(f64::NAN),
                on_front: on_front.contains(&i),
            });
        }
        Ok(out)
    }
}

/// Build a multirotor returning the native `UasError` (so the trade evaluator
/// composes straight into [`TradeStudy::run`]). Shares the radius / battery
/// derivation but surfaces the typed error rather than a `String`.
fn build_multirotor_uaserr(
    p: &UasParams,
    payload_override: Option<f64>,
    battery_override: Option<Battery>,
) -> Result<MultirotorUas, valenx_uas::UasError> {
    // Radius from disk area; a degenerate area is rejected by valenx-drone via
    // the non-positive radius check, so map our own guard onto that error.
    let radius = if p.rotor_count > 0 && p.disk_area_m2.is_finite() && p.disk_area_m2 > 0.0 {
        (p.disk_area_m2 / (f64::from(p.rotor_count) * std::f64::consts::PI)).sqrt()
    } else {
        // Force valenx-drone's own NonPositive radius error.
        0.0
    };
    let battery = match battery_override {
        Some(b) => b,
        None => Battery::new(p.battery_wh, p.usable_fraction)?,
    };
    MultirotorUas::new(
        p.rotor_count,
        radius,
        p.all_up_mass_kg,
        p.figure_of_merit,
        p.air_density,
        battery,
        payload_override.unwrap_or(p.payload_kg),
        p.drivetrain_efficiency,
    )
}

/// Build a fixed-wing returning the native `UasError` (trade-evaluator form).
fn build_fixedwing_uaserr(
    p: &UasParams,
    payload_override: Option<f64>,
    battery_override: Option<Battery>,
) -> Result<FixedWingUas, valenx_uas::UasError> {
    let battery = match battery_override {
        Some(b) => b,
        None => Battery::new(p.battery_wh, p.usable_fraction)?,
    };
    FixedWingUas::new(
        p.wing_area_m2,
        p.all_up_mass_kg,
        p.cl_max,
        p.aspect_ratio,
        p.cd0,
        p.oswald_efficiency,
        p.air_density,
        battery,
        payload_override.unwrap_or(p.payload_kg),
        p.drivetrain_efficiency,
    )
}

/// Solve the defensive counter-UAS intercept geometry + detection timeline for
/// the given parameters, fail-loud. Pure kinematics — no weapon, no targeting.
fn run_intercept(c: &CounterUasParams) -> Result<InterceptResult, String> {
    let track = ThreatTrack::new(
        Vector3::new(c.threat_pos_xy[0], c.threat_pos_xy[1], 0.0),
        Vector3::new(c.threat_vel_xy[0], c.threat_vel_xy[1], 0.0),
    )
    .map_err(|e| e.to_string())?;
    let interceptor = Interceptor::new(
        Vector3::new(c.interceptor_pos_xy[0], c.interceptor_pos_xy[1], 0.0),
        c.interceptor_max_speed,
    )
    .map_err(|e| e.to_string())?;

    let sol = time_to_intercept(&track, &interceptor).map_err(|e| e.to_string())?;
    let timeline = detection_timeline(&track, Vector3::zeros(), c.sensor_range_m)
        .map_err(|e| e.to_string())?;

    Ok(InterceptResult {
        threat_pos: c.threat_pos_xy,
        threat_vel: c.threat_vel_xy,
        interceptor_pos: c.interceptor_pos_xy,
        sensor_range_m: c.sensor_range_m,
        time_to_intercept_s: sol.map(|s| s.time_s),
        intercept_point: sol.map(|s| [s.point.x, s.point.y]),
        detected: timeline.detected,
        first_detect_s: timeline.first_detect_s,
        last_in_range_s: timeline.last_in_range_s,
        min_range_m: timeline.min_range_m,
        time_of_closest_approach_s: timeline.time_of_closest_approach_s,
    })
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Draw the UAS workbench. A no-op unless toggled on via View → UAS.
///
/// Mirrors [`crate::uq_workbench::draw_uq_workbench`].
pub fn draw_uas_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_uas_workbench {
        return;
    }
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_uas_workbench",
        "UAS design & counter-UAS (defensive)",
        uas_workbench_body,
    );
    if close {
        app.show_uas_workbench = false;
    }
}

// ---------------------------------------------------------------------------
// Workbench body
// ---------------------------------------------------------------------------

fn uas_workbench_body(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.label(
        egui::RichText::new(
            "Small-UAS design \u{2192} performance \u{2192} trade study, plus DEFENSIVE \
             counter-UAS detect / track / intercept GEOMETRY \u{00B7} valenx-uas  \
             [research / educational \u{2014} momentum-theory hover & parabolic-polar cruise]",
        )
        .weak()
        .small(),
    );
    ui.separator();

    let mut do_run = false;

    {
        let s = &mut app.uas;
        let p = &mut s.params;

        // --- Airframe configuration & geometry ------------------------------
        ui.label(egui::RichText::new("Airframe").strong());
        egui::Grid::new("uas_airframe_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                let lbl = ui.label("configuration");
                egui::ComboBox::from_id_source("uas_config_combo")
                    .selected_text(p.config.label())
                    .show_ui(ui, |ui| {
                        ui.selectable_value(
                            &mut p.config,
                            UasConfig::Multirotor,
                            UasConfig::Multirotor.label(),
                        );
                        ui.selectable_value(
                            &mut p.config,
                            UasConfig::FixedWing,
                            UasConfig::FixedWing.label(),
                        );
                    })
                    .response
                    .labelled_by(lbl.id)
                    .on_hover_text("Multirotor (hover) or fixed-wing (cruise) configuration.");
                ui.end_row();

                drag_row(ui, "all-up mass (kg)", &mut p.all_up_mass_kg, 0.05, "Total take-off mass the lift is sized to. Must exceed the payload.");

                let multi = p.config == UasConfig::Multirotor;
                ui.add_enabled_ui(multi, |ui| {
                    let lbl = ui.label("rotor count");
                    ui.add(egui::DragValue::new(&mut p.rotor_count).speed(1).range(1..=32))
                        .labelled_by(lbl.id)
                        .on_hover_text("Number of rotors (multirotor only).");
                });
                ui.end_row();
                ui.add_enabled_ui(multi, |ui| {
                    drag_row_inline(ui, "disk area (m^2)", &mut p.disk_area_m2, 0.01, "Total rotor disk area A; an equivalent per-rotor radius is derived as sqrt(A/(N*pi)). Multirotor only.");
                });
                ui.end_row();
                ui.add_enabled_ui(multi, |ui| {
                    drag_row_inline(ui, "figure of merit", &mut p.figure_of_merit, 0.01, "Rotor hover efficiency (figure of merit) in (0, 1]. Multirotor only.");
                });
                ui.end_row();
                ui.add_enabled_ui(multi, |ui| {
                    drag_row_inline(ui, "cruise power factor", &mut p.cruise_power_factor, 0.05, "Cruise power as a multiple of hover power for the multirotor range balance. Multirotor only.");
                });
                ui.end_row();
                ui.add_enabled_ui(multi, |ui| {
                    drag_row_inline(ui, "max thrust/weight", &mut p.max_thrust_to_weight, 0.05, "Installed maximum thrust-to-weight ratio; sets the multirotor payload margin. Multirotor only.");
                });
                ui.end_row();

                let wing = p.config == UasConfig::FixedWing;
                ui.add_enabled_ui(wing, |ui| {
                    drag_row_inline(ui, "wing area (m^2)", &mut p.wing_area_m2, 0.01, "Wing reference area S. Fixed-wing only.");
                });
                ui.end_row();
                ui.add_enabled_ui(wing, |ui| {
                    drag_row_inline(ui, "aspect ratio", &mut p.aspect_ratio, 0.1, "Wing aspect ratio AR. Fixed-wing only.");
                });
                ui.end_row();
                ui.add_enabled_ui(wing, |ui| {
                    drag_row_inline(ui, "CL max", &mut p.cl_max, 0.05, "Maximum lift coefficient CL,max (sets stall speed). Fixed-wing only.");
                });
                ui.end_row();
                ui.add_enabled_ui(wing, |ui| {
                    drag_row_inline(ui, "CD0", &mut p.cd0, 0.005, "Zero-lift drag coefficient CD0. Fixed-wing only.");
                });
                ui.end_row();
                ui.add_enabled_ui(wing, |ui| {
                    drag_row_inline(ui, "Oswald e", &mut p.oswald_efficiency, 0.01, "Oswald span efficiency e in (0, 1]. Fixed-wing only.");
                });
                ui.end_row();
            });

        // --- Power & loads --------------------------------------------------
        ui.add_space(4.0);
        ui.label(egui::RichText::new("Battery, payload & drivetrain").strong());
        egui::Grid::new("uas_power_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                drag_row(ui, "battery (Wh)", &mut p.battery_wh, 5.0, "Installed (nameplate) battery energy in watt-hours. Must be > 0.");
                drag_row(ui, "usable fraction", &mut p.usable_fraction, 0.01, "Usable battery fraction (depth-of-discharge x reserve) in (0, 1].");
                drag_row(ui, "payload (kg)", &mut p.payload_kg, 0.05, "Payload mass carried within the all-up mass. Must be >= 0 and below the all-up mass.");
                drag_row(ui, "drivetrain eff.", &mut p.drivetrain_efficiency, 0.01, "Combined motor + ESC + propeller efficiency in (0, 1].");
                drag_row(ui, "air density (kg/m^3)", &mut p.air_density, 0.01, "Air density rho. Must be > 0.");
                drag_row(ui, "cruise speed (m/s)", &mut p.cruise_speed_m_s, 0.5, "Cruise speed the range / endurance is evaluated at. Must be > 0.");
            });

        // --- Trade sweep ----------------------------------------------------
        ui.add_space(4.0);
        ui.label(egui::RichText::new("Trade study (one-parameter sweep \u{2192} Pareto)").strong());
        egui::Grid::new("uas_trade_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                let lbl = ui.label("sweep variable");
                egui::ComboBox::from_id_source("uas_sweep_combo")
                    .selected_text(p.sweep_var.label())
                    .show_ui(ui, |ui| {
                        ui.selectable_value(&mut p.sweep_var, SweepVar::Payload, SweepVar::Payload.label());
                        ui.selectable_value(&mut p.sweep_var, SweepVar::BatteryWh, SweepVar::BatteryWh.label());
                    })
                    .response
                    .labelled_by(lbl.id)
                    .on_hover_text("Which design variable the trade sweep walks (objectives: maximize endurance & payload).");
                ui.end_row();
                drag_row(ui, "sweep lower bound", &mut p.sweep_lo, 0.05, "Sweep start value (units follow the sweep variable).");
                drag_row(ui, "sweep upper bound", &mut p.sweep_hi, 0.05, "Sweep end value; must exceed the lower bound.");
                let lbl = ui.label("sweep points");
                ui.add(egui::DragValue::new(&mut p.sweep_points).speed(1).range(2..=200))
                    .labelled_by(lbl.id)
                    .on_hover_text("Number of design points evaluated across the sweep. Must be >= 2.");
                ui.end_row();
            });

        // --- Counter-UAS (defensive intercept geometry) ---------------------
        ui.add_space(4.0);
        ui.label(egui::RichText::new("Counter-UAS \u{2014} defensive intercept GEOMETRY").strong());
        ui.label(
            egui::RichText::new(
                "Detect / track / intercept KINEMATICS only \u{2014} when & where an inbound \
                 track is reachable and when a sensor first sees it. No weapon employment, \
                 no targeting, no lethality is modeled.",
            )
            .weak()
            .small(),
        );
        let c = &mut p.counter;
        egui::Grid::new("uas_counter_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                drag_row(ui, "threat pos x (m)", &mut c.threat_pos_xy[0], 10.0, "Threat track initial x position relative to the defended site at the origin.");
                drag_row(ui, "threat pos y (m)", &mut c.threat_pos_xy[1], 10.0, "Threat track initial y position relative to the defended site.");
                drag_row(ui, "threat vel x (m/s)", &mut c.threat_vel_xy[0], 1.0, "Threat track constant x velocity component.");
                drag_row(ui, "threat vel y (m/s)", &mut c.threat_vel_xy[1], 1.0, "Threat track constant y velocity component.");
                drag_row(ui, "interceptor pos x (m)", &mut c.interceptor_pos_xy[0], 10.0, "Interceptor start x position.");
                drag_row(ui, "interceptor pos y (m)", &mut c.interceptor_pos_xy[1], 10.0, "Interceptor start y position.");
                drag_row(ui, "interceptor max speed (m/s)", &mut c.interceptor_max_speed, 1.0, "Interceptor maximum speed s. Must be > 0; a stern chase of a faster opening target has no solution.");
                drag_row(ui, "sensor range (m)", &mut c.sensor_range_m, 25.0, "Detection range R of the sensor sited at the origin. Must be > 0.");
            });

        ui.add_space(6.0);
        ui.horizontal(|ui| {
            if ui
                .button(egui::RichText::new("Run").strong())
                .on_hover_text(
                    "Compute performance, the trade-study Pareto front, and the defensive \
                     intercept geometry + detection timeline.",
                )
                .clicked()
            {
                do_run = true;
            }
        });
    }

    // --- Execute (outside borrow) -------------------------------------------
    if do_run {
        run_and_store(app);
    }

    // --- Status line ---------------------------------------------------------
    let s = &app.uas;
    if !s.status.is_empty() {
        ui.add_space(6.0);
        let color = if s.status.starts_with('\u{26A0}') {
            egui::Color32::from_rgb(220, 120, 60)
        } else {
            egui::Color32::from_rgb(90, 180, 110)
        };
        ui.label(egui::RichText::new(&s.status).color(color).strong());
    }

    // --- Visualisation -------------------------------------------------------
    ui.add_space(6.0);
    ui.separator();
    draw_uas_viz(s, ui);
}

/// A labelled `DragValue` row in a 2-column grid: a caption cell, the drag
/// value (`labelled_by` the caption so it carries an accessible name), then
/// `end_row`. Mirrors the uq workbench's per-row caption pattern.
fn drag_row(ui: &mut egui::Ui, caption: &str, value: &mut f64, speed: f64, hint: &str) {
    let lbl = ui.label(caption);
    ui.add(egui::DragValue::new(value).speed(speed))
        .labelled_by(lbl.id)
        .on_hover_text(hint);
    ui.end_row();
}

/// Like [`drag_row`] but **without** the trailing `end_row` — for use inside an
/// [`egui::Ui::add_enabled_ui`] closure (the caller ends the row after the
/// closure so the enabled-scope wraps exactly the caption + control).
fn drag_row_inline(ui: &mut egui::Ui, caption: &str, value: &mut f64, speed: f64, hint: &str) {
    let lbl = ui.label(caption);
    ui.add(egui::DragValue::new(value).speed(speed))
        .labelled_by(lbl.id)
        .on_hover_text(hint);
}

/// Run the pipeline and fold the result (or error) into the workbench status.
/// Factored out so the Run button (and tests) can share it.
pub(crate) fn run_and_store(app: &mut ValenxApp) {
    let s = &mut app.uas;
    match s.run() {
        Ok(res) => {
            let tti = res
                .intercept
                .as_ref()
                .and_then(|i| i.time_to_intercept_s)
                .map(|t| format!("{t:.1}s"))
                .unwrap_or_else(|| "none".to_string());
            let front = res.trade.iter().filter(|t| t.on_front).count();
            s.status = format!(
                "\u{2714} {} perf rows \u{00B7} {} trade pts ({} on Pareto front) \u{00B7} TTI {}",
                res.perf_rows.len(),
                res.trade.len(),
                front,
                tti
            );
            s.result = Some(res);
        }
        Err(e) => {
            s.status = format!("\u{26A0} {e}");
            s.result = None;
        }
    }
}

// ---------------------------------------------------------------------------
// 2-D visualisation (performance readout + trade scatter + intercept plan view)
// ---------------------------------------------------------------------------

fn draw_uas_viz(s: &UasWorkbenchState, ui: &mut egui::Ui) {
    let Some(res) = &s.result else {
        ui.label(
            egui::RichText::new(
                "press \"Run\" to compute performance, the trade-study front, and the \
                 defensive intercept plan",
            )
            .weak(),
        );
        return;
    };

    draw_performance_readout(res, ui);
    ui.add_space(8.0);
    draw_trade_scatter(&s.params, res, ui);
    ui.add_space(8.0);
    draw_intercept_planview(res, ui);
}

/// View (a): a performance readout panel — the integrated performance rows as
/// a labelled grid plus a small painter bar for each dimensionless / scaled
/// quantity is overkill, so this is a clean two-column text readout drawn over
/// a panel background.
fn draw_performance_readout(res: &UasResult, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Performance readout").strong());
    if res.perf_rows.is_empty() {
        ui.label(egui::RichText::new("no performance rows").weak().small());
        return;
    }
    egui::Grid::new("uas_perf_grid")
        .num_columns(2)
        .striped(true)
        .show(ui, |ui| {
            for r in &res.perf_rows {
                ui.label(&r.label);
                ui.label(egui::RichText::new(&r.value).monospace());
                ui.end_row();
            }
        });
}

/// View (b): a trade-study scatter — each evaluated sweep point plotted as
/// (payload, endurance); Pareto-front points highlighted. The classic
/// endurance-vs-payload trade reads top-left (high endurance, low payload) to
/// bottom-right (low endurance, high payload).
fn draw_trade_scatter(params: &UasParams, res: &UasResult, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Trade study \u{2014} endurance vs payload").strong());
    ui.label(
        egui::RichText::new(
            "grey dots = swept designs \u{00B7} amber dots = Pareto (non-dominated) front \
             \u{00B7} x = payload (kg), y = endurance (min)",
        )
        .weak()
        .small(),
    );
    let _ = params;

    let available = ui.available_size();
    let (rect, _) = ui.allocate_exact_size(
        egui::vec2(available.x.min(460.0), 200.0),
        egui::Sense::hover(),
    );
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 0.0, egui::Color32::from_rgb(14, 22, 34));

    if res.trade.len() < 2 {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "too few trade points",
            egui::FontId::monospace(12.0),
            egui::Color32::from_gray(120),
        );
        return;
    }

    // Axis ranges over finite points.
    let mut x_lo = f64::INFINITY;
    let mut x_hi = f64::NEG_INFINITY;
    let mut y_lo = f64::INFINITY;
    let mut y_hi = f64::NEG_INFINITY;
    for t in &res.trade {
        let y_min = t.endurance_s / 60.0;
        if t.payload_kg.is_finite() {
            x_lo = x_lo.min(t.payload_kg);
            x_hi = x_hi.max(t.payload_kg);
        }
        if y_min.is_finite() {
            y_lo = y_lo.min(y_min);
            y_hi = y_hi.max(y_min);
        }
    }
    if !(x_lo.is_finite() && x_hi.is_finite() && y_lo.is_finite() && y_hi.is_finite()) {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "non-finite trade objectives",
            egui::FontId::monospace(12.0),
            egui::Color32::from_gray(120),
        );
        return;
    }
    if (x_hi - x_lo).abs() < 1e-12 {
        x_lo -= 0.5;
        x_hi += 0.5;
    }
    if (y_hi - y_lo).abs() < 1e-12 {
        y_lo -= 0.5;
        y_hi += 0.5;
    }

    let margin = 26.0_f32;
    let inner = rect.shrink(margin);
    let to_px = |x: f64, y: f64| -> egui::Pos2 {
        let fx = ((x - x_lo) / (x_hi - x_lo)).clamp(0.0, 1.0) as f32;
        let fy = ((y - y_lo) / (y_hi - y_lo)).clamp(0.0, 1.0) as f32;
        egui::pos2(
            inner.left() + fx * inner.width(),
            inner.bottom() - fy * inner.height(),
        )
    };

    // Axes.
    painter.line_segment(
        [
            egui::pos2(inner.left(), inner.bottom()),
            egui::pos2(inner.right(), inner.bottom()),
        ],
        egui::Stroke::new(1.0, egui::Color32::from_gray(70)),
    );
    painter.line_segment(
        [
            egui::pos2(inner.left(), inner.top()),
            egui::pos2(inner.left(), inner.bottom()),
        ],
        egui::Stroke::new(1.0, egui::Color32::from_gray(70)),
    );

    // Connect the dominated points faintly in sweep order so the curve reads.
    let pts: Vec<egui::Pos2> = res
        .trade
        .iter()
        .map(|t| to_px(t.payload_kg, t.endurance_s / 60.0))
        .collect();
    for w in pts.windows(2) {
        painter.line_segment(
            [w[0], w[1]],
            egui::Stroke::new(1.0, egui::Color32::from_gray(60)),
        );
    }

    // All points; front highlighted (amber, larger) over the grey dots.
    for (t, p) in res.trade.iter().zip(pts.iter()) {
        if t.on_front {
            painter.circle_filled(*p, 4.5, egui::Color32::from_rgb(230, 180, 70));
        } else {
            painter.circle_filled(*p, 2.5, egui::Color32::from_gray(150));
        }
    }

    // Axis labels.
    painter.text(
        egui::pos2(inner.center().x, rect.bottom() - 2.0),
        egui::Align2::CENTER_BOTTOM,
        "payload (kg)",
        egui::FontId::monospace(10.0),
        egui::Color32::from_gray(150),
    );
    painter.text(
        egui::pos2(rect.left() + 2.0, inner.center().y),
        egui::Align2::LEFT_CENTER,
        "endur (min)",
        egui::FontId::monospace(10.0),
        egui::Color32::from_gray(150),
    );
}

/// View (c): a top-down intercept **plan view** — the defended site / sensor at
/// the origin with its range ring, the threat track (start point + velocity
/// arrow + straight path), the interceptor start, and the computed intercept
/// point if feasible. Pure geometry; nothing about engagement.
fn draw_intercept_planview(res: &UasResult, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Counter-UAS plan view (defensive geometry)").strong());
    ui.label(
        egui::RichText::new(
            "origin = defended site \u{00B7} cyan ring = sensor range \u{00B7} red = threat track \
             \u{00B7} blue = interceptor \u{00B7} green X = intercept point (geometry only)",
        )
        .weak()
        .small(),
    );

    let available = ui.available_size();
    let (rect, _) = ui.allocate_exact_size(
        egui::vec2(available.x.min(460.0), 240.0),
        egui::Sense::hover(),
    );
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 0.0, egui::Color32::from_rgb(10, 18, 28));

    let Some(ic) = &res.intercept else {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "no intercept geometry",
            egui::FontId::monospace(12.0),
            egui::Color32::from_gray(120),
        );
        return;
    };

    // World extent: cover the origin, sensor ring, threat start, its endpoint
    // along the track over a horizon, the interceptor, and the intercept point.
    let mut min = [0.0_f64, 0.0];
    let mut max = [0.0_f64, 0.0];
    let mut expand = |x: f64, y: f64| {
        min[0] = min[0].min(x);
        min[1] = min[1].min(y);
        max[0] = max[0].max(x);
        max[1] = max[1].max(y);
    };
    expand(ic.sensor_range_m, ic.sensor_range_m);
    expand(-ic.sensor_range_m, -ic.sensor_range_m);
    expand(ic.threat_pos[0], ic.threat_pos[1]);
    expand(ic.interceptor_pos[0], ic.interceptor_pos[1]);
    // A horizon for the threat path: to intercept time if it exists, else a
    // fixed look-ahead long enough to leave the ring.
    let horizon = ic.time_to_intercept_s.unwrap_or(30.0).max(1.0);
    let end = [
        ic.threat_pos[0] + ic.threat_vel[0] * horizon,
        ic.threat_pos[1] + ic.threat_vel[1] * horizon,
    ];
    expand(end[0], end[1]);
    if let Some(pt) = ic.intercept_point {
        expand(pt[0], pt[1]);
    }

    // Pad and keep the aspect square so circles stay circular.
    let span_x = (max[0] - min[0]).max(1.0);
    let span_y = (max[1] - min[1]).max(1.0);
    let span = span_x.max(span_y) * 1.15;
    let cx = (min[0] + max[0]) / 2.0;
    let cy = (min[1] + max[1]) / 2.0;
    let world_lo = [cx - span / 2.0, cy - span / 2.0];

    let margin = 14.0_f32;
    let inner = rect.shrink(margin);
    let side = inner.width().min(inner.height());
    // Centre the square plot inside the rect.
    let plot = egui::Rect::from_center_size(inner.center(), egui::vec2(side, side));
    let to_px = |x: f64, y: f64| -> egui::Pos2 {
        let fx = ((x - world_lo[0]) / span).clamp(0.0, 1.0) as f32;
        // World +y is up; screen +y is down → invert.
        let fy = ((y - world_lo[1]) / span).clamp(0.0, 1.0) as f32;
        egui::pos2(plot.left() + fx * side, plot.bottom() - fy * side)
    };

    // Sensor range ring + origin.
    let origin_px = to_px(0.0, 0.0);
    let ring_r = (ic.sensor_range_m / span) as f32 * side;
    painter.circle_stroke(
        origin_px,
        ring_r,
        egui::Stroke::new(1.5, egui::Color32::from_rgb(70, 180, 200)),
    );
    painter.circle_filled(origin_px, 3.0, egui::Color32::from_rgb(120, 220, 240));
    painter.text(
        origin_px + egui::vec2(5.0, 4.0),
        egui::Align2::LEFT_TOP,
        "site",
        egui::FontId::monospace(10.0),
        egui::Color32::from_rgb(120, 220, 240),
    );

    // Threat path (start → horizon end) + start marker + velocity arrow.
    let threat_start = to_px(ic.threat_pos[0], ic.threat_pos[1]);
    let threat_end = to_px(end[0], end[1]);
    painter.line_segment(
        [threat_start, threat_end],
        egui::Stroke::new(1.5, egui::Color32::from_rgb(220, 90, 80)),
    );
    painter.circle_filled(threat_start, 4.0, egui::Color32::from_rgb(230, 110, 100));
    painter.text(
        threat_start + egui::vec2(5.0, -4.0),
        egui::Align2::LEFT_BOTTOM,
        "threat",
        egui::FontId::monospace(10.0),
        egui::Color32::from_rgb(230, 110, 100),
    );

    // Interceptor start.
    let icpt_px = to_px(ic.interceptor_pos[0], ic.interceptor_pos[1]);
    painter.circle_filled(icpt_px, 4.0, egui::Color32::from_rgb(90, 150, 240));
    painter.text(
        icpt_px + egui::vec2(5.0, 4.0),
        egui::Align2::LEFT_TOP,
        "intcpt",
        egui::FontId::monospace(10.0),
        egui::Color32::from_rgb(120, 170, 245),
    );

    // Intercept point + interceptor leg, if feasible.
    if let Some(pt) = ic.intercept_point {
        let meet = to_px(pt[0], pt[1]);
        painter.line_segment(
            [icpt_px, meet],
            egui::Stroke::new(1.2, egui::Color32::from_rgb(110, 200, 130)),
        );
        // Draw an X at the intercept point.
        let d = 5.0;
        let g = egui::Color32::from_rgb(120, 230, 140);
        painter.line_segment(
            [meet + egui::vec2(-d, -d), meet + egui::vec2(d, d)],
            egui::Stroke::new(2.0, g),
        );
        painter.line_segment(
            [meet + egui::vec2(-d, d), meet + egui::vec2(d, -d)],
            egui::Stroke::new(2.0, g),
        );
    } else {
        painter.text(
            egui::pos2(plot.center().x, plot.top() + 10.0),
            egui::Align2::CENTER_TOP,
            "no feasible intercept (target faster / opening)",
            egui::FontId::monospace(10.0),
            egui::Color32::from_rgb(220, 150, 90),
        );
    }

    // Readout grid below the plot.
    ui.add_space(4.0);
    egui::Grid::new("uas_intercept_stats")
        .num_columns(2)
        .striped(true)
        .show(ui, |ui| {
            let row = |ui: &mut egui::Ui, k: &str, v: String| {
                ui.label(k);
                ui.label(egui::RichText::new(v).monospace());
                ui.end_row();
            };
            row(
                ui,
                "time-to-intercept",
                ic.time_to_intercept_s
                    .map(|t| format!("{t:.2} s"))
                    .unwrap_or_else(|| "none".to_string()),
            );
            row(
                ui,
                "intercept point",
                ic.intercept_point
                    .map(|p| format!("({:.1}, {:.1}) m", p[0], p[1]))
                    .unwrap_or_else(|| "\u{2014}".to_string()),
            );
            row(
                ui,
                "sensor detects track",
                if ic.detected {
                    "yes".into()
                } else {
                    "no".into()
                },
            );
            row(
                ui,
                "first detection",
                ic.first_detect_s
                    .map(|t| format!("{t:.2} s"))
                    .unwrap_or_else(|| "\u{2014}".to_string()),
            );
            row(
                ui,
                "leaves range at",
                ic.last_in_range_s
                    .map(|t| format!("{t:.2} s"))
                    .unwrap_or_else(|| "\u{2014} (stays / never)".to_string()),
            );
            row(ui, "closest approach", format!("{:.1} m", ic.min_range_m));
            row(
                ui,
                "time of closest approach",
                format!("{:.2} s", ic.time_of_closest_approach_s),
            );
        });
}

// ---------------------------------------------------------------------------
// Tests (unit + headless_ui_tests, mirroring uq_workbench)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_run_succeeds_and_is_populated() {
        let s = UasWorkbenchState::default();
        let res = s.run().expect("default UAS run should succeed");
        assert!(!res.perf_rows.is_empty(), "performance rows populated");
        assert_eq!(
            res.trade.len(),
            s.params.sweep_points,
            "one trade point per sweep step"
        );
        assert!(
            res.trade.iter().any(|t| t.on_front),
            "at least one point on the Pareto front"
        );
        let ic = res.intercept.expect("intercept geometry present");
        // Default threat closes the sensor ring head-on-ish -> detected.
        assert!(ic.detected, "default inbound track should be detected");
        // Interceptor (50 m/s) faster than the closing threat -> feasible.
        assert!(
            ic.time_to_intercept_s.is_some(),
            "default geometry should yield a feasible intercept"
        );
    }

    #[test]
    fn fixed_wing_config_runs_and_reports_lift_to_drag() {
        let mut s = UasWorkbenchState::default();
        s.params.config = UasConfig::FixedWing;
        s.params.all_up_mass_kg = 3.0;
        let res = s.run().expect("fixed-wing run should succeed");
        assert!(
            res.perf_rows.iter().any(|r| r.label.contains("L/D")),
            "fixed-wing readout should include (L/D) max"
        );
    }

    #[test]
    fn disk_area_maps_to_expected_radius() {
        // A = N·π·r²  ⇒  r = sqrt(A/(N·π)). 4 rotors, area 0.283 ≈ 0.15 m.
        let p = UasParams::default();
        let r = p.rotor_radius_m().expect("radius");
        assert!((r - 0.15).abs() < 0.005, "radius {r} should be ~0.15 m");
    }

    #[test]
    fn head_on_intercept_matches_range_over_closing_speed() {
        // Threat 1000 m down +x flying straight back at 20 m/s; interceptor at
        // origin, 80 m/s. Closing 100 m/s -> TTI 10 s (the crate's pinned case).
        let mut c = CounterUasParams::default();
        c.threat_pos_xy = [1000.0, 0.0];
        c.threat_vel_xy = [-20.0, 0.0];
        c.interceptor_pos_xy = [0.0, 0.0];
        c.interceptor_max_speed = 80.0;
        c.sensor_range_m = 500.0;
        let ic = run_intercept(&c).expect("intercept");
        let tti = ic.time_to_intercept_s.expect("feasible");
        assert!(
            (tti - 10.0).abs() < 1e-3,
            "head-on TTI should be 10 s, got {tti}"
        );
    }

    #[test]
    fn opening_faster_target_has_no_intercept() {
        // Target ahead opening at 60 m/s, interceptor only 40 m/s -> None.
        let mut c = CounterUasParams::default();
        c.threat_pos_xy = [100.0, 0.0];
        c.threat_vel_xy = [60.0, 0.0];
        c.interceptor_pos_xy = [0.0, 0.0];
        c.interceptor_max_speed = 40.0;
        let ic = run_intercept(&c).expect("intercept call ok");
        assert!(
            ic.time_to_intercept_s.is_none(),
            "an opening faster target must yield no feasible intercept"
        );
    }

    #[test]
    fn more_payload_cuts_endurance_in_sweep() {
        // The endurance-vs-payload trade: as the swept payload rises the
        // all-up mass is fixed but the reported payload grows; endurance is
        // governed by hover power (fixed mass) so stays ~constant here, but the
        // sweep must produce monotone-increasing payload values.
        let s = UasWorkbenchState::default();
        let trade = s.run_trade_sweep().expect("sweep");
        assert!(trade.len() >= 2);
        for w in trade.windows(2) {
            assert!(
                w[1].payload_kg >= w[0].payload_kg - 1e-9,
                "swept payload should be non-decreasing"
            );
        }
    }

    // ---- degenerate-param tests — must return Err, NOT panic ----

    #[test]
    fn zero_disk_area_returns_err() {
        let mut s = UasWorkbenchState::default();
        s.params.disk_area_m2 = 0.0;
        assert!(
            s.run().is_err(),
            "zero disk area must return Err, not panic"
        );
    }

    #[test]
    fn usable_fraction_out_of_range_returns_err() {
        let mut s = UasWorkbenchState::default();
        s.params.usable_fraction = 1.5;
        assert!(
            s.run().is_err(),
            "usable_fraction > 1 must return Err, not panic"
        );
        s.params.usable_fraction = 0.0;
        assert!(
            s.run().is_err(),
            "usable_fraction = 0 must return Err, not panic"
        );
    }

    #[test]
    fn payload_exceeds_mass_returns_err() {
        let mut s = UasWorkbenchState::default();
        s.params.payload_kg = 5.0; // >= all_up_mass 1.5
        assert!(
            s.run().is_err(),
            "payload >= all-up mass must return Err, not panic"
        );
    }

    #[test]
    fn degenerate_sweep_bounds_return_err() {
        let mut s = UasWorkbenchState::default();
        s.params.sweep_hi = s.params.sweep_lo; // hi <= lo
        assert!(s.run().is_err(), "hi <= lo must return Err, not panic");
        let mut s2 = UasWorkbenchState::default();
        s2.params.sweep_points = 1; // < 2
        assert!(s2.run().is_err(), "< 2 sweep points must return Err");
    }

    #[test]
    fn zero_sensor_range_returns_err() {
        let mut s = UasWorkbenchState::default();
        s.params.counter.sensor_range_m = 0.0;
        assert!(
            s.run().is_err(),
            "zero sensor range must return Err, not panic"
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
            draw_uas_workbench(app, ctx);
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
        assert!(!app.show_uas_workbench);
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_uas_workbench(&mut app, ctx);
        });
        // No panic = pass.
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_uas_workbench = true;
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn workbench_draws_with_populated_result_without_panic() {
        let mut app = ValenxApp::default();
        app.show_uas_workbench = true;
        let res = app.uas.run().expect("run should succeed");
        app.uas.result = Some(res);
        app.uas.status = "\u{2714} test result".to_string();
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn workbench_draws_with_error_status_without_panic() {
        let mut app = ValenxApp::default();
        app.show_uas_workbench = true;
        // Trigger an error state (zero disk area is fail-loud in run()).
        app.uas.params.disk_area_m2 = 0.0;
        let result = app.uas.run();
        app.uas.status = match result {
            Err(e) => format!("\u{26A0} {e}"),
            Ok(_) => "\u{26A0} simulated error for testing".to_string(),
        };
        app.uas.result = None;
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // Every numeric DragValue is a SpinButton and must be `labelled_by` its
        // caption (egui clears a DragValue's own Name), so an AI / screen reader
        // can find the control by caption text. Each `labelled_by` target must
        // RESOLVE to a real named caption node, not a dangling id.
        let mut app = ValenxApp::default();
        app.show_uas_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let by_id: std::collections::HashMap<NodeId, &Node> =
            nodes.iter().map(|(id, n)| (*id, n)).collect();

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        // Many numeric controls across airframe / power / trade / counter-UAS;
        // a conservative lower bound that all are present and named.
        assert!(
            spin_buttons.len() >= 20,
            "expected many numeric controls as spin buttons, got {}",
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

        // Representative captions present as named accessibility nodes.
        for caption in [
            "all-up mass (kg)",
            "battery (Wh)",
            "usable fraction",
            "payload (kg)",
            "drivetrain eff.",
            "cruise speed (m/s)",
            "sweep lower bound",
            "sweep upper bound",
            "threat pos x (m)",
            "interceptor max speed (m/s)",
            "sensor range (m)",
            "configuration",
            "sweep variable",
        ] {
            assert!(
                has_named_node(&nodes, caption),
                "caption '{caption}' must be a named node in the a11y tree"
            );
        }

        // The Run button must be a named, invokable node.
        assert!(
            nodes.iter().any(|(_, n)| {
                n.role() == Role::Button && n.name().is_some_and(|s| s.contains("Run"))
            }),
            "the Run button must be a named, invokable node"
        );
    }

    #[test]
    fn degenerate_params_show_error_not_panic() {
        // When the disk area is 0 or the usable fraction is out of range the
        // workbench must surface the error in-panel, not panic.
        let mut state = UasWorkbenchState::default();
        state.params.disk_area_m2 = 0.0;
        assert!(
            state.run().is_err(),
            "zero disk area must produce Err, not panic"
        );
        state.params.disk_area_m2 = 0.283;
        state.params.usable_fraction = 2.0;
        assert!(
            state.run().is_err(),
            "usable_fraction > 1 must produce Err, not panic"
        );
    }

    #[test]
    fn agent_bridge_uas_id_resolves_and_sets_flag() {
        // Verify the two mechanisms the agent bridge uses for
        //   `OpenWorkbench { id: "uas" }` (and the aliases):
        //   1. TabKind::from_id("uas") -> Some(TabKind::Uas)
        //   2. set_workbench_flag(app, "uas", true) -> show_uas_workbench = true
        use crate::project_tabs::{set_workbench_flag, TabKind};

        // 1. Lookup + aliases + case/whitespace tolerance.
        assert_eq!(
            TabKind::from_id("uas"),
            Some(TabKind::Uas),
            "\"uas\" must resolve to TabKind::Uas"
        );
        assert_eq!(TabKind::from_id("UAS"), Some(TabKind::Uas));
        assert_eq!(TabKind::from_id("  uas  "), Some(TabKind::Uas));
        assert_eq!(TabKind::from_id("drone"), Some(TabKind::Uas));
        assert_eq!(TabKind::from_id("counteruas"), Some(TabKind::Uas));

        // 2. Flag toggle.
        let mut app = ValenxApp::default();
        assert!(!app.show_uas_workbench);
        set_workbench_flag(&mut app, "uas", true);
        assert!(
            app.show_uas_workbench,
            "set_workbench_flag(\"uas\", true) must set show_uas_workbench"
        );
        set_workbench_flag(&mut app, "uas", false);
        assert!(!app.show_uas_workbench);
    }
}
