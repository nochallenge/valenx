//! The right-side **CFD Workbench** panel — native 2-D incompressible
//! laminar computational fluid dynamics over `valenx-cfd-native`'s
//! SIMPLE solver (no external solver, no case directory).
//!
//! Distinct from the Aerodynamics / Wind Tunnel workbench (external
//! aerodynamics over `valenx-aero`): this is the *internal*-flow solver
//! — the two canonical CFD textbook cases, the **lid-driven cavity** and
//! **developing channel flow**, solved on a staggered finite-volume grid
//! by the SIMPLE pressure-velocity coupling.
//!
//! Mirrors the FEM / aero / astro workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_cfd_workbench`,
//! toggled from the View menu.

use eframe::egui;
use egui_plot::{Line, Plot, PlotPoints, VLine};

use valenx_cfd_native::{solve_simple, Boundaries, Fluid, Grid, SimpleControls};

use crate::ValenxApp;

/// Which canonical flow case the workbench solves.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
enum CfdCase {
    /// A square cavity whose top lid slides — the classic recirculation
    /// benchmark.
    #[default]
    LidDrivenCavity,
    /// Flow entering a channel at a uniform inlet speed and developing a
    /// parabolic profile.
    ChannelFlow,
}

impl CfdCase {
    fn label(self) -> &'static str {
        match self {
            CfdCase::LidDrivenCavity => "lid-driven cavity",
            CfdCase::ChannelFlow => "channel flow",
        }
    }
}

/// Persistent state for the CFD Workbench.
pub struct CfdWorkbenchState {
    // Staggered grid (cells + domain in metres).
    nx: usize,
    ny: usize,
    lx: f64,
    ly: f64,
    // Fluid: density + kinematic viscosity.
    density: f64,
    viscosity: f64,
    // Drive speed: lid speed (cavity) or inlet speed (channel), m/s.
    speed: f64,
    // SIMPLE outer-iteration cap (defaults low for a responsive UI;
    // the engine default is 4000).
    max_iterations: usize,
    case: CfdCase,
    result: String,
    error: Option<String>,
    /// Vertical centreline velocity profile: [speed (m/s), height y (m)].
    profile: Option<Vec<[f64; 2]>>,
    /// Analytic fully-developed Poiseuille profile for channel flow, overlaid on
    /// the plot for comparison; `None` for the lid-driven cavity.
    analytic_profile: Option<Vec<[f64; 2]>>,
    /// Bulk (mean-throughflow) velocity drawn as a reference line on the channel
    /// profile plot (⅔ of the Poiseuille peak); `None` for the lid-driven cavity.
    bulk_velocity: Option<f64>,
}

impl Default for CfdWorkbenchState {
    fn default() -> Self {
        // A unit cavity at Re = U·L/ν = 100 (U=1, L=1, ν=0.01).
        Self {
            nx: 32,
            ny: 32,
            lx: 1.0,
            ly: 1.0,
            density: 1.0,
            viscosity: 0.01,
            speed: 1.0,
            max_iterations: 500,
            case: CfdCase::LidDrivenCavity,
            result: String::new(),
            error: None,
            profile: None,
            analytic_profile: None,
            bulk_velocity: None,
        }
    }
}

/// The characteristic length used for the Reynolds number: the cavity
/// width for the lid-driven case, the channel height for channel flow.
fn characteristic_length(s: &CfdWorkbenchState) -> f64 {
    match s.case {
        CfdCase::LidDrivenCavity => s.lx,
        CfdCase::ChannelFlow => s.ly,
    }
}

/// Flow regime inferred from the Reynolds number.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum FlowRegime {
    Laminar,
    Transitional,
    Turbulent,
}

impl FlowRegime {
    fn label(self) -> &'static str {
        match self {
            FlowRegime::Laminar => "laminar",
            FlowRegime::Transitional => "transitional",
            FlowRegime::Turbulent => "turbulent",
        }
    }
}

/// Classify the flow regime from the Reynolds number. The channel case uses
/// the nominal internal-flow thresholds (laminar < 2300, transitional < 4000,
/// turbulent above); the lid-driven cavity stays steady/laminar to much higher
/// Re — its first instability is near Re ≈ 8000 — so its thresholds are raised.
/// This is the validity gate for the laminar SIMPLE solver: a turbulent Re
/// means the results are not physical.
fn flow_regime(re: f64, case: CfdCase) -> FlowRegime {
    let (laminar_max, turbulent_min) = match case {
        CfdCase::ChannelFlow => (2300.0, 4000.0),
        CfdCase::LidDrivenCavity => (8000.0, 10000.0),
    };
    if re < laminar_max {
        FlowRegime::Laminar
    } else if re < turbulent_min {
        FlowRegime::Transitional
    } else {
        FlowRegime::Turbulent
    }
}

/// Draw the CFD Workbench right-side panel. A no-op when the
/// `show_cfd_workbench` toggle is off.
pub fn draw_cfd_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_cfd_workbench {
        return;
    }
    egui::SidePanel::right("valenx_cfd_workbench")
        .resizable(true)
        .default_width(360.0)
        .width_range(300.0..=560.0)
        .show(ctx, |ui| {
            ui.heading("CFD Workbench");
            ui.label(
                egui::RichText::new("native 2-D incompressible CFD · valenx-cfd-native")
                    .weak()
                    .small(),
            );
            ui.separator();
            let s = &mut app.cfd;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Case").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.case, CfdCase::LidDrivenCavity, "Lid-driven cavity")
                            .on_hover_text("Square cavity, sliding top lid — the recirculation benchmark.");
                        ui.radio_value(&mut s.case, CfdCase::ChannelFlow, "Channel flow")
                            .on_hover_text("Uniform inlet developing a parabolic profile.");
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Grid — staggered finite volume").strong());
                    ui.horizontal(|ui| {
                        ui.label("nx");
                        ui.add(egui::DragValue::new(&mut s.nx).speed(0.5));
                        ui.label("ny");
                        ui.add(egui::DragValue::new(&mut s.ny).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("Lx");
                        ui.add(egui::DragValue::new(&mut s.lx).speed(0.05));
                        ui.label("Ly");
                        ui.add(egui::DragValue::new(&mut s.ly).speed(0.05));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Fluid").strong());
                    ui.horizontal(|ui| {
                        ui.label("density ρ (kg/m³)");
                        ui.add(egui::DragValue::new(&mut s.density).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("kinematic ν (m²/s)");
                        ui.add(egui::DragValue::new(&mut s.viscosity).speed(0.001));
                    });
                    let drive = match s.case {
                        CfdCase::LidDrivenCavity => "lid speed U (m/s)",
                        CfdCase::ChannelFlow => "inlet speed U (m/s)",
                    };
                    ui.horizontal(|ui| {
                        ui.label(drive);
                        ui.add(egui::DragValue::new(&mut s.speed).speed(0.05));
                    });
                    if s.viscosity > 0.0 {
                        let re = s.speed.abs() * characteristic_length(s) / s.viscosity;
                        let regime = flow_regime(re, s.case);
                        ui.label(
                            egui::RichText::new(format!(
                                "Reynolds number ≈ {re:.1}  ({})",
                                regime.label()
                            ))
                            .weak()
                            .small(),
                        );
                        match regime {
                            FlowRegime::Turbulent => {
                                ui.colored_label(
                                    egui::Color32::from_rgb(220, 90, 90),
                                    "⚠ turbulent — this laminar solver will be unphysical here",
                                );
                            }
                            FlowRegime::Transitional => {
                                ui.colored_label(
                                    egui::Color32::from_rgb(220, 170, 80),
                                    "⚠ transitional — laminar results are approximate",
                                );
                            }
                            FlowRegime::Laminar => {}
                        }
                    }

                    ui.add_space(4.0);
                    ui.horizontal(|ui| {
                        ui.label("max SIMPLE iterations");
                        ui.add(egui::DragValue::new(&mut s.max_iterations).speed(5.0));
                    });
                    ui.label(
                        egui::RichText::new("Solve runs synchronously — a finer grid or more iterations takes longer.")
                            .weak()
                            .small(),
                    );

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Solve flow").strong())
                        .clicked()
                    {
                        run_cfd(s);
                    }

                    if let Some(e) = &s.error {
                        ui.add_space(4.0);
                        ui.colored_label(egui::Color32::from_rgb(220, 90, 90), e);
                    }
                    if !s.result.is_empty() {
                        ui.separator();
                        ui.label(egui::RichText::new("Result").strong());
                        ui.label(egui::RichText::new(&s.result).monospace());
                    }
                    if let Some(profile) = &s.profile {
                        ui.add_space(4.0);
                        ui.label(egui::RichText::new("Centreline speed vs height").strong());
                        Plot::new("cfd_profile_plot")
                            .height(160.0)
                            .x_axis_label("speed (m/s)")
                            .y_axis_label("y (m)")
                            .show(ui, |pui| {
                                pui.line(Line::new(PlotPoints::from(profile.clone())).name("|u|"));
                                if let Some(analytic) = &s.analytic_profile {
                                    pui.line(
                                        Line::new(PlotPoints::from(analytic.clone()))
                                            .name("Poiseuille (analytic)"),
                                    );
                                }
                                if let Some(ub) = s.bulk_velocity {
                                    pui.vline(VLine::new(ub).name("U_bulk (Q/H)"));
                                }
                            });
                    }
                });
        });
}

/// Free-stream dynamic pressure `q = ½ ρ U²` (Pa) — the pressure scale that
/// sizes hydrodynamic loads, from the fluid density and the drive speed.
fn dynamic_pressure(density: f64, speed: f64) -> f64 {
    0.5 * density * speed * speed
}

/// Cell (grid) Reynolds number `Re_cell = U·Δx/ν` for a streamwise cell size
/// `cell_size` — a numerical-resolution diagnostic, not a property of the
/// flow. It equals the global Reynolds number divided by the streamwise cell
/// count, so a coarse grid (large `Re_cell`) under-resolves the convective
/// term and the first-order upwind scheme smears sharp gradients.
fn cell_reynolds(speed: f64, cell_size: f64, viscosity: f64) -> f64 {
    speed.abs() * cell_size / viscosity
}

/// Build the grid + BCs and run the SIMPLE solver. Extracted from the
/// draw closure so it is unit-testable, and it validates every input
/// before calling [`Grid::new`] (which *panics* on a bad grid) — so a
/// bad input surfaces as an error, never a crash.
/// The analytic fully-developed plane-Poiseuille velocity profile for a channel
/// of height `height` driven at bulk speed `inlet_speed`, sampled at the same
/// `ny` cell centres as the CFD profile: `u(y) = 1.5·U·(1 − ((y − h/2)/(h/2))²)`
/// — the parabola the channel solve should converge to far from the inlet.
fn poiseuille_profile(inlet_speed: f64, height: f64, ny: usize) -> Vec<[f64; 2]> {
    let u_max = 1.5 * inlet_speed;
    let half = 0.5 * height;
    (0..ny)
        .map(|j| {
            let y = (j as f64 + 0.5) * height / ny as f64;
            let eta = (y - half) / half;
            [u_max * (1.0 - eta * eta), y]
        })
        .collect()
}

fn run_cfd(s: &mut CfdWorkbenchState) {
    s.error = None;
    s.profile = None;
    s.analytic_profile = None;
    s.bulk_velocity = None;
    if s.nx == 0 || s.ny == 0 {
        s.error = Some("grid must have at least one cell per axis".into());
        return;
    }
    if !(s.lx > 0.0 && s.ly > 0.0 && s.lx.is_finite() && s.ly.is_finite()) {
        s.error = Some("domain dimensions must be positive and finite".into());
        return;
    }
    if !(s.viscosity > 0.0 && s.viscosity.is_finite()) {
        s.error = Some("kinematic viscosity must be positive".into());
        return;
    }
    if !(s.density > 0.0 && s.density.is_finite()) {
        s.error = Some("density must be positive".into());
        return;
    }

    let grid = Grid::new(s.nx, s.ny, s.lx, s.ly);
    let fluid = Fluid::new(s.density, s.viscosity);
    let bcs = match s.case {
        CfdCase::LidDrivenCavity => Boundaries::lid_driven_cavity(s.speed),
        CfdCase::ChannelFlow => Boundaries::channel_flow(s.speed),
    };
    let controls = SimpleControls {
        max_iterations: s.max_iterations,
        ..Default::default()
    };

    let sol = solve_simple(&grid, &fluid, &bcs, &controls);

    let max_speed = sol.max_speed();
    let mean_speed = sol.mean_speed();
    let re = s.speed.abs() * characteristic_length(s) / s.viscosity;
    let regime = flow_regime(re, s.case);
    let q = dynamic_pressure(s.density, s.speed);
    // Streamwise grid spacing Δx = lx / nx (nx > 0 guaranteed by the guard above).
    let re_cell = cell_reynolds(s.speed, s.lx / s.nx as f64, s.viscosity);
    // Static-pressure swing Δp = p_max − p_min (gauge-independent).
    let dp = sol.pressure_range();
    // Total (stagnation) pressure swing Δp₀ = max(p+½ρ|u|²) − min(…) — the Bernoulli
    // total-pressure loss, the irreversible degradation the static Δp cannot see.
    let dp0 = sol.total_pressure_range(s.density);
    // Peak vorticity |∂v/∂x − ∂u/∂y| — the strongest local rotation.
    let vorticity = sol.max_vorticity();
    // Location of that peak — the vortex core (when an interior difference exists).
    let vort_loc = sol
        .peak_vorticity_location()
        .map(|(x, y)| format!("  (core @ {x:.2}, {y:.2} m)"))
        .unwrap_or_else(|| "  (max rotation)".to_string());
    // Suction peak: where the static pressure bottoms out (cavitation / vortex core).
    let suction_loc = sol
        .min_pressure_location()
        .map(|(x, y)| format!("  (suction peak @ {x:.2}, {y:.2} m)"))
        .unwrap_or_default();
    // Stagnation point: where the static pressure peaks (the flow brought to rest).
    let stagnation_loc = sol
        .max_pressure_location()
        .map(|(x, y)| format!("  (stagnation @ {x:.2}, {y:.2} m)"))
        .unwrap_or_default();
    // Where the flow is fastest — the convective hot-spot (CFL-limiting cell).
    let peak_speed_loc = sol
        .max_speed_location()
        .map(|(x, y)| format!("  (peak @ {x:.2}, {y:.2} m)"))
        .unwrap_or_default();

    // Vertical centreline velocity profile: speed vs height.
    let i_mid = s.nx / 2;
    let profile: Vec<[f64; 2]> = (0..s.ny)
        .map(|j| {
            let y = (j as f64 + 0.5) * s.ly / s.ny as f64;
            [sol.speed_at_cell(i_mid, j), y]
        })
        .collect();
    s.profile = Some(profile);
    // For channel flow, the analytic Poiseuille parabola the solve converges to
    // (the lid-driven cavity has no such 1-D profile).
    s.analytic_profile = match s.case {
        CfdCase::ChannelFlow => Some(poiseuille_profile(s.speed, s.ly, s.ny)),
        CfdCase::LidDrivenCavity => None,
    };
    // Bulk velocity reference (channel only): U_bulk = Q_in / H, ⅔ of the peak.
    s.bulk_velocity = match s.case {
        CfdCase::ChannelFlow => Some(sol.bulk_velocity()),
        CfdCase::LidDrivenCavity => None,
    };

    // Through-flow throughput + a global mass-continuity check (channel only;
    // the enclosed cavity has no net inflow, so the relative error is undefined).
    let flow_str = match s.case {
        CfdCase::ChannelFlow => format!(
            "\nflow rate  : {:.5} m²/s  (continuity err {:.2}%)",
            sol.inlet_flow_rate(),
            100.0 * sol.continuity_error()
        ),
        CfdCase::LidDrivenCavity => String::new(),
    };

    // Wall shear stress τ_w = μ·(∂u/∂y)|_wall = ρν·(wall shear rate), bottom + top
    // walls (the top wall is the moving lid in the lid-driven cavity).
    let tau_w = s.density * s.viscosity * sol.bottom_wall_shear_rate();
    let tau_w_top = s.density * s.viscosity * sol.top_wall_shear_rate();
    // Side (vertical) walls: τ_w = ρν·(∂v/∂x)|_{x=0, lx} — the recirculation's drag on
    // the left and right walls, the wall-normal axis the horizontal-wall entries miss.
    let tau_w_left = s.density * s.viscosity * sol.left_wall_shear_rate();
    let tau_w_right = s.density * s.viscosity * sol.right_wall_shear_rate();
    s.result = format!(
        "case       : {}\n\
         grid       : {}×{}  ({:.3} × {:.3} m)\n\
         Reynolds   : {:.1}  ({})\n\
         iterations : {} {}\n\
         residual   : {:.3e}\n\
         max |u|    : {:.5} m/s  (mean {:.5}){peak_speed_loc}\n\
         dynamic q  : {:.4} Pa  (½ρU²; mean KE {:.4}; \u{03A6}_visc {:.3e} W/m)\n\
         cell Re    : {:.2}  (U·Δx/ν; ≳2 ⇒ convection under-resolved)\n\
         pressure Δp: {:.4e} Pa  (p_max−p_min){suction_loc}{stagnation_loc}  \u{00B7}  total \u{0394}p\u{2080} {dp0:.3e} Pa  (Bernoulli loss)\n\
         peak vort  : {:.4} 1/s{vort_loc}  \u{00B7}  Q_max {:.4} 1/s\u{00B2}  \u{00B7}  S_max {:.4} 1/s{flow_str}\n\
         circulation: {:.4} m\u{00B2}/s  (\u{222B}\u{03C9}\u{00B7}dA, signed)  \u{00B7}  enstrophy {:.4} m\u{00B2}/s\u{00B2}  (\u{00BD}\u{222B}\u{03C9}\u{00B2}\u{00B7}dA)  \u{00B7}  palin {:.3e} 1/s\u{00B2}  (\u{00BD}\u{222B}|\u{2207}\u{03C9}|\u{00B2}\u{00B7}dA)\n\
         wall shear : {:.4e} Pa  (\u{03C4}_w, bottom)  \u{00B7}  {tau_w_top:.4e} Pa (top)  \u{00B7}  {tau_w_left:.4e} Pa (left)  \u{00B7}  {tau_w_right:.4e} Pa (right)  \u{00B7}  max|\u{2207}\u{00B7}u| {:.2e} 1/s  (peak local continuity residual)  \u{00B7}  reverse-flow {:.0}%  (recirculating area)  \u{00B7}  \u{03C8}-span {:.4} m\u{00B2}/s  (streamline range)",
        s.case.label(),
        s.nx,
        s.ny,
        s.lx,
        s.ly,
        re,
        regime.label(),
        sol.iterations,
        if sol.converged {
            "(converged)"
        } else {
            "(hit iteration cap)"
        },
        sol.residual,
        max_speed,
        mean_speed,
        q,
        sol.mean_kinetic_energy_density(s.density),
        sol.viscous_dissipation(s.density * s.viscosity),
        re_cell,
        dp,
        vorticity,
        sol.max_q_criterion(),
        sol.max_strain_rate(),
        sol.circulation(),
        sol.enstrophy(),
        sol.palinstrophy(),
        tau_w,
        sol.max_divergence(),
        100.0 * sol.reverse_flow_fraction(),
        sol.stream_function_range(),
    );

    // ── Wall, boundary-layer & flow-rate diagnostics (FlowSolution methods) ──
    // s.viscosity is the kinematic viscosity ν; the dynamic viscosity is μ = ρ·ν.
    let mu = s.density * s.viscosity;
    let re_bulk = sol.bulk_reynolds_number(s.viscosity);
    let tau_w_method = sol.bottom_wall_shear_stress(mu);
    let c_f = sol.skin_friction_coefficient(mu, s.density);
    let eu = sol.euler_number(s.density);
    let u_tau = sol.friction_velocity(s.viscosity);
    let m_in = sol.inlet_mass_flow_rate(s.density);
    s.result.push_str(&format!(
        "\nwall & flow: \u{03C4}_w {tau_w_method:.4e} Pa (method)  \u{00B7}  C_f {c_f:.6}  \u{00B7}  \
         u_\u{03C4} {u_tau:.5} m/s  \u{00B7}  Eu {eu:.4}  \u{00B7}  Re_bulk {re_bulk:.1}  \u{00B7}  \
         \u{1E41}_in {m_in:.4} kg/s"
    ));

    // Field statistics (any case): speed min/CoV, pressure mean/rms, vorticity mean/rms.
    let min_speed = sol.min_speed();
    let cov_speed = sol.speed_coefficient_of_variation();
    let mean_press = sol.mean_pressure();
    let rms_press = sol.rms_pressure();
    let mean_vort = sol.mean_vorticity();
    let rms_vort = sol.rms_vorticity();
    s.result.push_str(&format!(
        "\nfield stats: |u|_min {min_speed:.5} m/s  \u{00B7}  speed CoV {cov_speed:.4}  \u{00B7}  \
         p_mean {mean_press:.4e} Pa  \u{00B7}  p_rms {rms_press:.4e} Pa  \u{00B7}  \
         \u{03C9}_mean {mean_vort:.4} 1/s  \u{00B7}  \u{03C9}_rms {rms_vort:.4} 1/s"
    ));

    // Flow regime & throughput (any case): domain/cell Reynolds, residence time, outlet mass.
    let re_domain = sol.domain_reynolds_number(s.viscosity);
    let re_cell = sol.cell_reynolds_number(s.viscosity);
    let tau_residence = sol.flow_through_time();
    let m_out = sol.outlet_mass_flow_rate(s.density);
    s.result.push_str(&format!(
        "\nflow regime: Re_L {re_domain:.1}  \u{00B7}  Re_h {re_cell:.2}  \u{00B7}  \
         \u{03C4}_res {tau_residence:.4} s  \u{00B7}  \u{1E41}_out {m_out:.4} kg/s"
    ));
    // Shear-layer & energy-dissipation diagnostics (any case).
    let delta_omega = sol.vorticity_thickness();
    let eps = sol.mean_dissipation_rate(s.viscosity);
    s.result.push_str(&format!(
        "\nshear & dissipation: \u{03B4}_\u{03C9} {delta_omega:.5} m  \u{00B7}  \
         \u{03B5} {eps:.3e} m\u{00B2}/s\u{00B3}"
    ));
    if matches!(s.case, CfdCase::ChannelFlow) {
        let u_ref = sol.bulk_velocity();
        if u_ref > 0.0 {
            let ds = sol.displacement_thickness(u_ref);
            let th = sol.momentum_thickness(u_ref);
            let de = sol.energy_thickness(u_ref);
            let h = sol.shape_factor(u_ref);
            let hs = sol.energy_shape_factor(u_ref);
            let re_theta = sol.momentum_thickness_reynolds_number(u_ref, s.viscosity);
            let re_dstar = sol.displacement_thickness_reynolds_number(u_ref, s.viscosity);
            s.result.push_str(&format!(
                "\nboundary layer: \u{03B4}* {ds:.5} m  \u{00B7}  \u{03B8} {th:.5} m  \u{00B7}  \
                 \u{03B4}_E {de:.5} m  \u{00B7}  H {h:.3}  \u{00B7}  H* {hs:.3}  \u{00B7}  \
                 Re_\u{03B8} {re_theta:.0}  \u{00B7}  Re_\u{03B4}* {re_dstar:.0}"
            ));
            // Turbulence cascade scales (kinematic ν): Taylor microscale + Re_λ, Kolmogorov
            // length/time, integral length/time. Diagnostic for the resolved field.
            let lambda = sol.taylor_microscale(s.viscosity);
            let re_lambda = sol.taylor_reynolds_number(s.viscosity);
            let eta = sol.kolmogorov_length_scale(s.viscosity);
            let tau_eta = sol.kolmogorov_time_scale(s.viscosity);
            let integral_l = sol.integral_length_scale(s.viscosity);
            let integral_t = sol.integral_time_scale(s.viscosity);
            s.result.push_str(&format!(
                "\nturbulence scales: \u{03BB} {lambda:.4e} m  \u{00B7}  Re_\u{03BB} {re_lambda:.1}  \
                 \u{00B7}  \u{03B7} {eta:.4e} m  \u{00B7}  \u{03C4}_\u{03B7} {tau_eta:.4e} s  \
                 \u{00B7}  L {integral_l:.4e} m  \u{00B7}  T {integral_t:.4e} s"
            ));
            // Scale ratios (dimensionless) + Kolmogorov velocity.
            let lambda_eta = sol.taylor_to_kolmogorov_ratio(s.viscosity);
            let l_lambda = sol.integral_to_taylor_ratio(s.viscosity);
            let l_eta = sol.integral_to_kolmogorov_ratio(s.viscosity);
            let t_tau = sol.integral_to_kolmogorov_time_ratio(s.viscosity);
            let u_eta = sol.kolmogorov_velocity_scale(s.viscosity);
            s.result.push_str(&format!(
                "\nscale ratios: \u{03BB}/\u{03B7} {lambda_eta:.2}  \u{00B7}  L/\u{03BB} {l_lambda:.2}  \
                 \u{00B7}  L/\u{03B7} {l_eta:.2}  \u{00B7}  T/\u{03C4}_\u{03B7} {t_tau:.2}  \
                 \u{00B7}  u_\u{03B7} {u_eta:.4e} m/s"
            ));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lid_driven_cavity_solves() {
        // Small grid + low iteration cap keep the headless test fast.
        let mut s = CfdWorkbenchState {
            nx: 16,
            ny: 16,
            max_iterations: 100,
            ..Default::default()
        };
        run_cfd(&mut s);
        assert!(s.error.is_none(), "unexpected error: {:?}", s.error);
        assert!(s.result.contains("max |u|"));
        // The recirculating cavity has rotation everywhere → the enstrophy
        // readout is surfaced alongside the circulation.
        assert!(s.result.contains("enstrophy"), "enstrophy in result: {}", s.result);
        // The cavity's primary vortex is rotation-dominated → the Q-criterion
        // readout is surfaced on the peak-vorticity line.
        assert!(s.result.contains("Q_max"), "Q-criterion in result: {}", s.result);
        // The sheared cavity flow also deforms → the peak strain-rate readout is
        // surfaced alongside (the strain companion to peak vorticity).
        assert!(s.result.contains("S_max"), "strain rate in result: {}", s.result);
        // The sheared cavity flow dissipates energy → the viscous-dissipation
        // readout is surfaced on the energy line.
        assert!(s.result.contains("\u{03A6}_visc"), "dissipation in result: {}", s.result);
        // The vorticity field is non-uniform → the palinstrophy readout (the
        // gradient-of-vorticity diagnostic) is surfaced on the cascade line.
        assert!(s.result.contains("palin"), "palinstrophy in result: {}", s.result);
        assert!(s.profile.as_ref().is_some_and(|p| p.len() == 16), "centreline profile sampled");
    }

    #[test]
    fn channel_flow_solves() {
        let mut s = CfdWorkbenchState {
            nx: 20,
            ny: 12,
            max_iterations: 100,
            case: CfdCase::ChannelFlow,
            ..Default::default()
        };
        run_cfd(&mut s);
        assert!(s.error.is_none(), "unexpected error: {:?}", s.error);
        assert!(s.result.contains("channel flow"));
    }

    #[test]
    fn poiseuille_profile_is_a_parabola_peaking_at_centreline() {
        let prof = poiseuille_profile(1.0, 1.0, 16);
        assert_eq!(prof.len(), 16);
        let max_u = prof.iter().map(|p| p[0]).fold(0.0_f64, f64::max);
        // Centreline speed ≈ 1.5× the inlet (bulk) speed.
        assert!((max_u - 1.5).abs() < 0.02, "centreline ~1.5× inlet: {max_u}");
        // Symmetric about the centre, non-negative, slower at the walls.
        assert!((prof[0][0] - prof[15][0]).abs() < 1e-9, "symmetric");
        assert!(prof[0][0] < max_u);
        assert!(prof.iter().all(|p| p[0] >= 0.0));
    }

    #[test]
    fn channel_flow_gets_an_analytic_overlay_cavity_does_not() {
        let mut chan = CfdWorkbenchState {
            nx: 20,
            ny: 12,
            max_iterations: 100,
            case: CfdCase::ChannelFlow,
            ..Default::default()
        };
        run_cfd(&mut chan);
        assert!(
            chan.analytic_profile.as_ref().is_some_and(|p| p.len() == chan.ny),
            "channel flow should carry the analytic Poiseuille overlay"
        );
        // The lid-driven cavity has no 1-D analytic profile.
        let mut cav = CfdWorkbenchState {
            case: CfdCase::LidDrivenCavity,
            ..Default::default()
        };
        run_cfd(&mut cav);
        assert!(cav.analytic_profile.is_none(), "cavity has no analytic overlay");
    }

    #[test]
    fn degenerate_grid_fails_loud() {
        // nx=0 would panic Grid::new — run_cfd must guard it.
        let mut s = CfdWorkbenchState {
            nx: 0,
            ..Default::default()
        };
        run_cfd(&mut s);
        assert!(s.error.is_some(), "nx=0 must surface an error, not panic");
    }

    #[test]
    fn zero_viscosity_fails_loud() {
        let mut s = CfdWorkbenchState {
            viscosity: 0.0,
            ..Default::default()
        };
        run_cfd(&mut s);
        assert!(s.error.is_some(), "zero viscosity must surface an error");
    }

    #[test]
    fn flow_regime_classifies_by_reynolds() {
        // Channel/pipe thresholds: laminar < 2300, transitional < 4000, then turbulent.
        assert_eq!(flow_regime(100.0, CfdCase::ChannelFlow), FlowRegime::Laminar);
        assert_eq!(flow_regime(3000.0, CfdCase::ChannelFlow), FlowRegime::Transitional);
        assert_eq!(flow_regime(5000.0, CfdCase::ChannelFlow), FlowRegime::Turbulent);
        // The lid-driven cavity stays laminar to much higher Re.
        assert_eq!(flow_regime(3000.0, CfdCase::LidDrivenCavity), FlowRegime::Laminar);
        assert_eq!(flow_regime(9000.0, CfdCase::LidDrivenCavity), FlowRegime::Transitional);
        // The default cavity (Re = 100) reports "laminar" in the solve result.
        let mut s = CfdWorkbenchState {
            nx: 12,
            ny: 12,
            max_iterations: 60,
            ..Default::default()
        };
        run_cfd(&mut s);
        assert!(s.result.contains("laminar"), "regime in result: {}", s.result);
    }

    #[test]
    fn dynamic_pressure_is_half_rho_u_squared() {
        assert!((dynamic_pressure(1.0, 10.0) - 50.0).abs() < 1e-12);
        // Quadratic in speed: doubling U quadruples q.
        assert!(
            (dynamic_pressure(1.2, 20.0) - 4.0 * dynamic_pressure(1.2, 10.0)).abs() < 1e-9
        );
        // Linear in density.
        assert!(
            (dynamic_pressure(2.0, 10.0) - 2.0 * dynamic_pressure(1.0, 10.0)).abs() < 1e-9
        );
    }

    #[test]
    fn cell_reynolds_is_speed_times_cell_over_nu() {
        // Re_cell = U·Δx/ν; 10 · 0.1 / 0.01 = 100.
        assert!((cell_reynolds(10.0, 0.1, 0.01) - 100.0).abs() < 1e-9);
        // Equals the global Reynolds number divided by the streamwise cell
        // count: with L = 1, ν = 0.01, U = 10 → Re = 1000; over nx = 10 cells
        // (Δx = 0.1) that is Re/nx = 100.
        let re_global = 10.0_f64 * 1.0 / 0.01;
        assert!((cell_reynolds(10.0, 1.0 / 10.0, 0.01) - re_global / 10.0).abs() < 1e-9);
        // Linear in speed and cell size, inverse in viscosity.
        assert!(
            (cell_reynolds(20.0, 0.1, 0.01) - 2.0 * cell_reynolds(10.0, 0.1, 0.01)).abs() < 1e-9
        );
        assert!(
            (cell_reynolds(10.0, 0.1, 0.02) - 0.5 * cell_reynolds(10.0, 0.1, 0.01)).abs() < 1e-9
        );
        // Sign-independent in the drive direction (uses |U|).
        assert!((cell_reynolds(-10.0, 0.1, 0.01) - cell_reynolds(10.0, 0.1, 0.01)).abs() < 1e-9);
    }
}
