//! State, conversions and formatting for the Wind Tunnel workbench.
//!
//! Everything in this module is **pure, non-UI logic** — the form
//! state structs the panels mutate, the `valenx-aero` request builder
//! that turns those forms into an [`valenx_aero::AeroRequest`], the
//! unit conversions (km/h ↔ m/s, degrees ↔ radians), and the result
//! formatters. Splitting it out keeps the `panels` module purely about
//! egui layout and makes the workbench's logic `#[test]`-coverable
//! without standing up an egui context.

use valenx_aero::{AeroRequest, Air, BoundaryConditions, TunnelSizing, TurbulenceModel};

// ---------------------------------------------------------------------------
// Enumerations the form exposes
// ---------------------------------------------------------------------------

/// Which body the wind tunnel tests.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub enum BodySource {
    /// The CAD solid currently loaded in the app (operand A).
    #[default]
    CurrentCadModel,
    /// A triangle mesh imported from an STL file.
    ImportedStl,
    /// A built-in parametric box — the canonical bluff-body / sanity
    /// case, always available even with nothing else loaded.
    DemoBox,
}

impl BodySource {
    /// Human-readable label.
    pub fn label(self) -> &'static str {
        match self {
            BodySource::CurrentCadModel => "Current CAD model",
            BodySource::ImportedStl => "Imported STL file",
            BodySource::DemoBox => "Built-in demo box",
        }
    }
}

/// Grid-resolution preset — sets the cells-across-body count and the
/// total-cell cap that `valenx-aero`'s [`TunnelSizing`] enforces.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub enum GridResolution {
    /// A fast, coarse pass — good for a first look / iterating geometry.
    Coarse,
    /// The balanced default.
    #[default]
    Medium,
    /// A fine grid — slower, sharper coefficients.
    Fine,
}

impl GridResolution {
    /// All presets in display order.
    pub const ALL: [GridResolution; 3] = [
        GridResolution::Coarse,
        GridResolution::Medium,
        GridResolution::Fine,
    ];

    /// Human-readable label.
    pub fn label(self) -> &'static str {
        match self {
            GridResolution::Coarse => "Coarse (fast)",
            GridResolution::Medium => "Medium (balanced)",
            GridResolution::Fine => "Fine (accurate)",
        }
    }

    /// Target number of cells across the body's smallest dimension.
    pub fn cells_across_body(self) -> usize {
        match self {
            GridResolution::Coarse => 10,
            GridResolution::Medium => 16,
            GridResolution::Fine => 26,
        }
    }

    /// Hard cap on the total cell count for this preset.
    pub fn max_cells(self) -> usize {
        match self {
            GridResolution::Coarse => 400_000,
            GridResolution::Medium => 2_000_000,
            GridResolution::Fine => 6_000_000,
        }
    }
}

/// Which turbulence closure the solver runs.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub enum TurbulenceChoice {
    /// Standard k-ε — robust, the generic default.
    KEpsilon,
    /// Menter k-ω SST — the external-aero standard.
    #[default]
    KOmegaSST,
}

impl TurbulenceChoice {
    /// All choices in display order.
    pub const ALL: [TurbulenceChoice; 2] =
        [TurbulenceChoice::KOmegaSST, TurbulenceChoice::KEpsilon];

    /// Human-readable label.
    pub fn label(self) -> &'static str {
        match self {
            TurbulenceChoice::KEpsilon => "k-epsilon",
            TurbulenceChoice::KOmegaSST => "k-omega SST",
        }
    }

    /// The `valenx-aero` turbulence model this choice maps to.
    pub fn model(self) -> TurbulenceModel {
        match self {
            TurbulenceChoice::KEpsilon => TurbulenceModel::KEpsilon,
            TurbulenceChoice::KOmegaSST => TurbulenceModel::KOmegaSST,
        }
    }
}

/// Solver mode — a single steady operating point or an angle-of-attack
/// polar sweep.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub enum RunMode {
    /// One steady RANS solve at the configured wind direction.
    #[default]
    SteadyPoint,
    /// A sequence of steady solves across an angle-of-attack range,
    /// assembled into a lift / drag polar.
    AngleSweep,
}

impl RunMode {
    /// Human-readable label.
    pub fn label(self) -> &'static str {
        match self {
            RunMode::SteadyPoint => "Steady point",
            RunMode::AngleSweep => "Angle-of-attack sweep",
        }
    }
}

/// Which scalar field a flow-visualization push colours the geometry
/// by. Drives both the surface-mesh build and the colour-legend label.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub enum FlowField {
    /// Surface pressure coefficient `Cp` on the body.
    #[default]
    SurfaceCp,
    /// Surface skin-friction coefficient `Cf` on the body.
    SkinFriction,
    /// Velocity magnitude on a domain cut plane.
    VelocityMagnitude,
    /// A static-pressure cut-plane slice.
    PressureSlice,
    /// The wake / vortex field — Q-criterion on a cut plane.
    VortexQ,
}

impl FlowField {
    /// All fields in display order.
    pub const ALL: [FlowField; 5] = [
        FlowField::SurfaceCp,
        FlowField::SkinFriction,
        FlowField::VelocityMagnitude,
        FlowField::PressureSlice,
        FlowField::VortexQ,
    ];

    /// Short button label.
    pub fn label(self) -> &'static str {
        match self {
            FlowField::SurfaceCp => "Surface Cp",
            FlowField::SkinFriction => "Skin friction Cf",
            FlowField::VelocityMagnitude => "Velocity magnitude",
            FlowField::PressureSlice => "Pressure cut plane",
            FlowField::VortexQ => "Wake / Q-criterion",
        }
    }

    /// The field name shown in the viewport colour-bar legend, with
    /// units where it has them.
    pub fn legend_name(self) -> &'static str {
        match self {
            FlowField::SurfaceCp => "Cp (surface)",
            FlowField::SkinFriction => "Cf (surface)",
            FlowField::VelocityMagnitude => "|U| (m/s)",
            FlowField::PressureSlice => "p (Pa)",
            FlowField::VortexQ => "Q-criterion (1/s2)",
        }
    }

    /// `true` if the field lives on the body surface (so the body
    /// triangle shell is coloured), `false` if it is a volumetric cut
    /// plane (so a flat plane patch is coloured instead).
    pub fn is_surface(self) -> bool {
        matches!(self, FlowField::SurfaceCp | FlowField::SkinFriction)
    }
}

/// Which axis a cut-plane visualization slices through.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub enum CutAxis {
    /// A plane of constant x (a y-z slice).
    X,
    /// A plane of constant y (an x-z slice) — the centreline plane.
    #[default]
    Y,
    /// A plane of constant z (an x-y slice).
    Z,
}

impl CutAxis {
    /// All axes in display order.
    pub const ALL: [CutAxis; 3] = [CutAxis::X, CutAxis::Y, CutAxis::Z];

    /// Short label.
    pub fn label(self) -> &'static str {
        match self {
            CutAxis::X => "X (frontal)",
            CutAxis::Y => "Y (centreline)",
            CutAxis::Z => "Z (plan)",
        }
    }

    /// The `valenx-aero` slice axis.
    pub fn slice_axis(self) -> valenx_aero::SliceAxis {
        match self {
            CutAxis::X => valenx_aero::SliceAxis::X,
            CutAxis::Y => valenx_aero::SliceAxis::Y,
            CutAxis::Z => valenx_aero::SliceAxis::Z,
        }
    }
}

// ---------------------------------------------------------------------------
// Unit conversions
// ---------------------------------------------------------------------------

/// Convert kilometres-per-hour to metres-per-second.
pub fn kmh_to_ms(kmh: f64) -> f64 {
    kmh / 3.6
}

/// Convert metres-per-second to kilometres-per-hour.
pub fn ms_to_kmh(ms: f64) -> f64 {
    ms * 3.6
}

/// Convert degrees to radians.
pub fn deg_to_rad(deg: f64) -> f64 {
    deg * std::f64::consts::PI / 180.0
}

/// Convert radians to degrees.
pub fn rad_to_deg(rad: f64) -> f64 {
    rad * 180.0 / std::f64::consts::PI
}

// ---------------------------------------------------------------------------
// The wind-tunnel form state
// ---------------------------------------------------------------------------

/// The full set of wind-tunnel form inputs the panels mutate.
///
/// One instance lives on [`crate::aero_workbench::AeroWorkbenchState`].
/// Defaults describe a sensible road-car case at 30 m/s — the workbench
/// runs out of the box with no input beyond pressing **Run**.
///
/// Implements `PartialEq` so the `History<WindTunnelForm>` undo stack
/// can deduplicate identical consecutive snapshots — see
/// [`crate::undo::History::record`].
#[derive(Clone, Debug, PartialEq)]
pub struct WindTunnelForm {
    // --- 1. Body -----------------------------------------------------
    /// Which body to test.
    pub body_source: BodySource,
    /// The demo box's dimensions (m) — used when `body_source` is
    /// [`BodySource::DemoBox`].
    pub demo_box_size: [f64; 3],

    // --- 2. Wind conditions -----------------------------------------
    /// Free-stream speed (m/s).
    pub speed_ms: f64,
    /// Wind yaw angle (degrees) — sideslip.
    pub yaw_deg: f64,
    /// Wind pitch angle / angle of attack (degrees).
    pub pitch_deg: f64,
    /// Air density (kg/m3).
    pub air_density: f64,
    /// Air dynamic viscosity (Pa·s).
    pub air_viscosity: f64,
    /// Upstream turbulence intensity (fraction, 0..1).
    pub turbulence_intensity: f64,
    /// Absolute air temperature (K) — feeds the Mach number.
    pub temperature_k: f64,
    /// Apply the Prandtl-Glauert compressibility correction.
    pub apply_compressibility: bool,
    /// Wing aspect ratio AR = span²/area (dimensionless) — for induced drag.
    pub aspect_ratio: f64,
    /// Oswald span-efficiency factor `e` (0 < e ≤ 1) — for induced drag.
    pub span_efficiency: f64,

    // --- 3. Ground & wheels -----------------------------------------
    /// `true` for a moving ground (road) carried at the free-stream
    /// speed; `false` for free-air (no ground).
    pub moving_ground: bool,
    /// `true` to spin the wheels (automotive) — rolling without slip.
    pub rotating_wheels: bool,
    /// Wheel radius (m) — sets the rolling angular speed.
    pub wheel_radius: f64,

    // --- 4. Tunnel & mesh -------------------------------------------
    /// Grid-resolution preset.
    pub resolution: GridResolution,
    /// Override the auto-sized domain clearances.
    pub override_domain: bool,
    /// Upstream clearance (body-lengths) — used when `override_domain`.
    pub upstream_lengths: f64,
    /// Downstream clearance (body-lengths) — the wake needs room.
    pub downstream_lengths: f64,
    /// Lateral / vertical clearance (multiples of body size).
    pub lateral_lengths: f64,

    // --- 5. Solver --------------------------------------------------
    /// Turbulence closure.
    pub turbulence: TurbulenceChoice,
    /// Maximum SIMPLE outer iterations.
    pub max_iterations: usize,
    /// Convergence tolerance on the mass-imbalance residual.
    pub tolerance: f64,
    /// Steady point vs. angle-of-attack sweep.
    pub run_mode: RunMode,
    /// Sweep start angle (degrees) — `RunMode::AngleSweep` only.
    pub sweep_start_deg: f64,
    /// Sweep end angle (degrees).
    pub sweep_end_deg: f64,
    /// Number of sweep points.
    pub sweep_points: usize,

    // --- 8. Flow visualization --------------------------------------
    /// Which field a flow-viz push colours the geometry by.
    pub flow_field: FlowField,
    /// Cut-plane axis for volumetric fields.
    pub cut_axis: CutAxis,
    /// Cut-plane position as a fraction (0..1) along its axis.
    pub cut_fraction: f64,
}

impl Default for WindTunnelForm {
    fn default() -> Self {
        WindTunnelForm {
            body_source: BodySource::default(),
            demo_box_size: [4.2, 1.8, 1.4],
            speed_ms: 30.0,
            yaw_deg: 0.0,
            pitch_deg: 0.0,
            air_density: 1.225,
            air_viscosity: 1.81e-5,
            turbulence_intensity: 0.02,
            temperature_k: 288.0,
            apply_compressibility: false,
            aspect_ratio: 6.0,
            span_efficiency: 0.95,
            moving_ground: true,
            rotating_wheels: false,
            wheel_radius: 0.33,
            resolution: GridResolution::default(),
            override_domain: false,
            upstream_lengths: 3.0,
            downstream_lengths: 8.0,
            lateral_lengths: 3.0,
            turbulence: TurbulenceChoice::default(),
            max_iterations: 200,
            tolerance: 1e-4,
            run_mode: RunMode::default(),
            sweep_start_deg: -4.0,
            sweep_end_deg: 12.0,
            sweep_points: 7,
            flow_field: FlowField::default(),
            cut_axis: CutAxis::default(),
            cut_fraction: 0.5,
        }
    }
}

impl WindTunnelForm {
    /// Build the validated [`Air`] this form describes.
    pub fn air(&self) -> Result<Air, String> {
        Air::new(self.air_density, self.air_viscosity).map_err(|e| e.to_string())
    }

    /// The dynamic pressure `q∞ = ½·ρ·U∞²` (Pa) implied by the form.
    pub fn dynamic_pressure(&self) -> f64 {
        0.5 * self.air_density * self.speed_ms * self.speed_ms
    }

    /// The free-stream Reynolds number for a body of characteristic
    /// length `length_m` — `Re = ρ·U·L/μ`. Returns `0` for a
    /// non-positive viscosity.
    pub fn reynolds_number(&self, length_m: f64) -> f64 {
        if self.air_viscosity <= 0.0 {
            return 0.0;
        }
        self.air_density * self.speed_ms * length_m / self.air_viscosity
    }

    /// The rolling-wheel angular speed (rad/s) for the configured wheel
    /// radius and free-stream speed — `ω = V/R`.
    pub fn wheel_omega(&self) -> f64 {
        valenx_aero::presets::rolling_wheel_omega(self.speed_ms, self.wheel_radius)
    }

    /// The [`BoundaryConditions`] the form's ground choice implies.
    pub fn boundary_conditions(&self) -> BoundaryConditions {
        if self.moving_ground {
            // A moving ground carries the road at the free-stream
            // speed (the car and the road move together).
            BoundaryConditions::automotive(self.speed_ms)
        } else {
            BoundaryConditions::external_aero()
        }
    }

    /// The [`TunnelSizing`] the resolution preset + any override imply.
    pub fn tunnel_sizing(&self) -> TunnelSizing {
        let mut s = TunnelSizing {
            cells_across_body: self.resolution.cells_across_body(),
            max_cells: self.resolution.max_cells(),
            ..TunnelSizing::default()
        };
        if self.override_domain {
            s.upstream = self.upstream_lengths.max(0.5);
            s.downstream = self.downstream_lengths.max(0.5);
            s.lateral = self.lateral_lengths.max(0.5);
        }
        s
    }

    /// Build the [`AeroRequest`] this form describes.
    ///
    /// Returns an error string only for an invalid air model — every
    /// other field is range-bounded by the UI widgets.
    pub fn build_request(&self) -> Result<AeroRequest, String> {
        let air = self.air()?;
        let mut req = AeroRequest::new(self.speed_ms)
            .with_yaw(deg_to_rad(self.yaw_deg))
            .with_angle_of_attack(deg_to_rad(self.pitch_deg))
            .with_air(air)
            .with_turbulence(self.turbulence.model())
            .with_boundary(self.boundary_conditions())
            .with_sizing(self.tunnel_sizing())
            .with_max_iterations(self.max_iterations.max(1))
            .with_compressibility(self.apply_compressibility);
        req.turbulence_intensity = self.turbulence_intensity.clamp(0.0, 1.0);
        req.tolerance = self.tolerance.max(1e-9);
        req.temperature_k = self.temperature_k.max(1.0);
        Ok(req)
    }

    /// The list of angle-of-attack values (radians) a sweep would
    /// solve at, evenly spaced from `sweep_start_deg` to
    /// `sweep_end_deg` with `sweep_points` points.
    pub fn sweep_angles(&self) -> Vec<f64> {
        valenx_aero::sweep::linspace_degrees(
            self.sweep_start_deg,
            self.sweep_end_deg,
            self.sweep_points.max(1),
        )
    }
}

// ---------------------------------------------------------------------------
// Result formatting
// ---------------------------------------------------------------------------

/// Format a Reynolds number in compact scientific notation, e.g.
/// `"1.31e7"`.
pub fn format_reynolds(re: f64) -> String {
    if re <= 0.0 || !re.is_finite() {
        return "—".to_string();
    }
    format!("{re:.2e}")
}

/// Format an aerodynamic coefficient with a sign and four decimals.
pub fn format_coefficient(c: f64) -> String {
    if !c.is_finite() {
        return "—".to_string();
    }
    format!("{c:+.4}")
}

/// Lift-to-drag ratio `L/D = Cl / Cd` — the headline aerodynamic-efficiency
/// number. `None` when drag is ~zero or non-finite (an empty tunnel or a
/// degenerate run), where the ratio is meaningless.
pub fn lift_to_drag(cl: f64, cd: f64) -> Option<f64> {
    if !cl.is_finite() || !cd.is_finite() || cd.abs() < 1e-9 {
        None
    } else {
        Some(cl / cd)
    }
}

/// Format a force in newtons, switching to kN above 1000 N.
pub fn format_force_n(force: f64) -> String {
    if !force.is_finite() {
        return "—".to_string();
    }
    if force.abs() >= 1000.0 {
        format!("{:.2} kN", force / 1000.0)
    } else {
        format!("{force:.1} N")
    }
}

/// Format a pressure in pascals, switching to kPa above 1000 Pa.
pub fn format_pressure_pa(pa: f64) -> String {
    if !pa.is_finite() {
        return "—".to_string();
    }
    if pa.abs() >= 1000.0 {
        format!("{:.3} kPa", pa / 1000.0)
    } else {
        format!("{pa:.1} Pa")
    }
}

/// Format the pressure-drag fraction of the total drag as a percentage
/// string, e.g. `"82% pressure / 18% friction"`.
pub fn format_drag_split(pressure_fraction: f64) -> String {
    let p = (pressure_fraction.clamp(0.0, 1.0) * 100.0).round() as i64;
    let f = 100 - p;
    format!("{p}% pressure / {f}% friction")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unit_conversions_round_trip() {
        // km/h ↔ m/s: 108 km/h is exactly 30 m/s.
        assert!((kmh_to_ms(108.0) - 30.0).abs() < 1e-12);
        assert!((ms_to_kmh(30.0) - 108.0).abs() < 1e-12);
        // A round trip preserves the value.
        let v = 42.7;
        assert!((ms_to_kmh(kmh_to_ms(v)) - v).abs() < 1e-9);
        // deg ↔ rad: 180° is π.
        assert!((deg_to_rad(180.0) - std::f64::consts::PI).abs() < 1e-12);
        assert!((rad_to_deg(std::f64::consts::PI) - 180.0).abs() < 1e-9);
        assert!((rad_to_deg(deg_to_rad(33.0)) - 33.0).abs() < 1e-9);
    }

    #[test]
    fn default_form_is_a_sensible_road_car() {
        let f = WindTunnelForm::default();
        // 30 m/s sea-level air, a moving ground, k-ω SST.
        assert!((f.speed_ms - 30.0).abs() < 1e-12);
        assert!((f.air_density - 1.225).abs() < 1e-12);
        assert!(f.moving_ground);
        assert_eq!(f.turbulence, TurbulenceChoice::KOmegaSST);
        assert_eq!(f.run_mode, RunMode::SteadyPoint);
        // The air model the default builds is valid.
        assert!(f.air().is_ok());
    }

    #[test]
    fn dynamic_pressure_matches_the_textbook_value() {
        let f = WindTunnelForm {
            speed_ms: 40.0,
            air_density: 1.225,
            ..WindTunnelForm::default()
        };
        // q = ½·1.225·40² = 980 Pa.
        assert!((f.dynamic_pressure() - 980.0).abs() < 1e-6);
    }

    #[test]
    fn reynolds_number_is_in_the_road_car_regime() {
        let f = WindTunnelForm::default();
        // A 4 m car at 30 m/s in sea-level air → Re ~ 8e6.
        let re = f.reynolds_number(4.2);
        assert!(re > 1.0e6 && re < 2.0e7, "Re {re} out of range");
        // A non-positive viscosity is guarded.
        let bad = WindTunnelForm {
            air_viscosity: 0.0,
            ..WindTunnelForm::default()
        };
        assert_eq!(bad.reynolds_number(4.0), 0.0);
    }

    #[test]
    fn boundary_conditions_follow_the_ground_choice() {
        use valenx_aero::FaceBc;
        // Moving ground → a moving-wall floor at the free-stream speed.
        let car = WindTunnelForm {
            moving_ground: true,
            speed_ms: 30.0,
            ..WindTunnelForm::default()
        };
        assert_eq!(car.boundary_conditions().z_min, FaceBc::MovingWall(30.0));
        // Free-air → a slip floor.
        let air = WindTunnelForm {
            moving_ground: false,
            ..WindTunnelForm::default()
        };
        assert_eq!(air.boundary_conditions().z_min, FaceBc::Slip);
    }

    #[test]
    fn tunnel_sizing_reflects_the_resolution_preset() {
        // The fine preset asks for more cells across the body than the
        // coarse one and a bigger cap.
        let coarse = WindTunnelForm {
            resolution: GridResolution::Coarse,
            ..WindTunnelForm::default()
        };
        let fine = WindTunnelForm {
            resolution: GridResolution::Fine,
            ..WindTunnelForm::default()
        };
        assert!(fine.tunnel_sizing().cells_across_body > coarse.tunnel_sizing().cells_across_body);
        assert!(fine.tunnel_sizing().max_cells > coarse.tunnel_sizing().max_cells);
    }

    #[test]
    fn tunnel_sizing_override_replaces_the_clearances() {
        let f = WindTunnelForm {
            override_domain: true,
            upstream_lengths: 5.0,
            downstream_lengths: 12.0,
            lateral_lengths: 4.0,
            ..WindTunnelForm::default()
        };
        let s = f.tunnel_sizing();
        assert!((s.upstream - 5.0).abs() < 1e-12);
        assert!((s.downstream - 12.0).abs() < 1e-12);
        assert!((s.lateral - 4.0).abs() < 1e-12);
        // Without the override the defaults stand.
        let d = WindTunnelForm {
            override_domain: false,
            ..WindTunnelForm::default()
        };
        let ds = d.tunnel_sizing();
        let def = TunnelSizing::default();
        assert!((ds.upstream - def.upstream).abs() < 1e-12);
    }

    #[test]
    fn build_request_carries_the_form_through() {
        let f = WindTunnelForm {
            speed_ms: 50.0,
            yaw_deg: 6.0,
            pitch_deg: 3.0,
            turbulence: TurbulenceChoice::KEpsilon,
            max_iterations: 120,
            apply_compressibility: true,
            ..WindTunnelForm::default()
        };
        let req = f.build_request().expect("valid request");
        assert!((req.speed - 50.0).abs() < 1e-12);
        // 6° yaw and 3° pitch made it through as radians.
        assert!((req.yaw - deg_to_rad(6.0)).abs() < 1e-12);
        assert!((req.pitch - deg_to_rad(3.0)).abs() < 1e-12);
        assert_eq!(req.turbulence, TurbulenceModel::KEpsilon);
        assert_eq!(req.max_iterations, 120);
        assert!(req.apply_compressibility);
        // The request builds a valid wind.
        assert!(req.wind().is_ok());
    }

    #[test]
    fn build_request_rejects_a_non_physical_air() {
        let f = WindTunnelForm {
            air_density: -1.0,
            ..WindTunnelForm::default()
        };
        assert!(f.build_request().is_err());
    }

    #[test]
    fn wheel_omega_is_speed_over_radius() {
        let f = WindTunnelForm {
            speed_ms: 33.0,
            wheel_radius: 0.33,
            ..WindTunnelForm::default()
        };
        // ω = 33 / 0.33 = 100 rad/s.
        assert!((f.wheel_omega() - 100.0).abs() < 1e-9);
    }

    #[test]
    fn sweep_angles_span_the_configured_range() {
        let f = WindTunnelForm {
            sweep_start_deg: -4.0,
            sweep_end_deg: 8.0,
            sweep_points: 4,
            ..WindTunnelForm::default()
        };
        let angles = f.sweep_angles();
        assert_eq!(angles.len(), 4);
        assert!((angles[0] - deg_to_rad(-4.0)).abs() < 1e-12);
        assert!((angles[3] - deg_to_rad(8.0)).abs() < 1e-12);
        // Strictly ascending.
        for w in angles.windows(2) {
            assert!(w[1] > w[0]);
        }
    }

    #[test]
    fn turbulence_intensity_is_clamped_into_the_request() {
        // An out-of-range intensity is clamped to [0, 1] so the
        // request builder never produces an ill-posed wind.
        let f = WindTunnelForm {
            turbulence_intensity: 5.0,
            ..WindTunnelForm::default()
        };
        let req = f.build_request().expect("valid");
        assert!((0.0..=1.0).contains(&req.turbulence_intensity));
    }

    #[test]
    fn format_reynolds_uses_scientific_notation() {
        assert_eq!(format_reynolds(13_100_000.0), "1.31e7");
        // Non-positive / non-finite degrade to a dash.
        assert_eq!(format_reynolds(0.0), "—");
        assert_eq!(format_reynolds(f64::NAN), "—");
    }

    #[test]
    fn format_coefficient_shows_sign_and_four_decimals() {
        assert_eq!(format_coefficient(0.31), "+0.3100");
        assert_eq!(format_coefficient(-0.5), "-0.5000");
        assert_eq!(format_coefficient(f64::INFINITY), "—");
    }

    #[test]
    fn lift_to_drag_is_cl_over_cd_and_guards_zero_drag() {
        assert_eq!(lift_to_drag(0.6, 0.3), Some(2.0));
        assert_eq!(lift_to_drag(-0.4, 0.2), Some(-2.0));
        // Zero / non-finite drag → no meaningful ratio.
        assert_eq!(lift_to_drag(0.5, 0.0), None);
        assert_eq!(lift_to_drag(0.5, f64::NAN), None);
        assert_eq!(lift_to_drag(f64::INFINITY, 0.3), None);
    }

    #[test]
    fn format_force_switches_to_kilonewtons() {
        assert_eq!(format_force_n(250.0), "250.0 N");
        assert_eq!(format_force_n(2500.0), "2.50 kN");
        assert_eq!(format_force_n(f64::NAN), "—");
    }

    #[test]
    fn format_pressure_switches_to_kilopascals() {
        assert_eq!(format_pressure_pa(500.0), "500.0 Pa");
        assert_eq!(format_pressure_pa(101_325.0), "101.325 kPa");
    }

    #[test]
    fn format_drag_split_sums_to_one_hundred() {
        assert_eq!(format_drag_split(0.82), "82% pressure / 18% friction");
        assert_eq!(format_drag_split(0.0), "0% pressure / 100% friction");
        // Out-of-range fractions are clamped.
        assert_eq!(format_drag_split(1.5), "100% pressure / 0% friction");
    }

    #[test]
    fn flow_field_surface_classification_is_correct() {
        assert!(FlowField::SurfaceCp.is_surface());
        assert!(FlowField::SkinFriction.is_surface());
        assert!(!FlowField::VelocityMagnitude.is_surface());
        assert!(!FlowField::PressureSlice.is_surface());
        assert!(!FlowField::VortexQ.is_surface());
    }

    #[test]
    fn enum_label_tables_are_populated() {
        // Every enum the UI iterates must hand back a non-empty label
        // for every variant — an empty one would be an invisible
        // control.
        for r in GridResolution::ALL {
            assert!(!r.label().is_empty());
            assert!(r.cells_across_body() > 0);
        }
        for t in TurbulenceChoice::ALL {
            assert!(!t.label().is_empty());
        }
        for f in FlowField::ALL {
            assert!(!f.label().is_empty());
            assert!(!f.legend_name().is_empty());
        }
        for a in CutAxis::ALL {
            assert!(!a.label().is_empty());
        }
        for b in [
            BodySource::CurrentCadModel,
            BodySource::ImportedStl,
            BodySource::DemoBox,
        ] {
            assert!(!b.label().is_empty());
        }
        for m in [RunMode::SteadyPoint, RunMode::AngleSweep] {
            assert!(!m.label().is_empty());
        }
    }
}
