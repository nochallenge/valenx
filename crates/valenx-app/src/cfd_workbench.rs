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
use egui_plot::{Line, PlotPoints, VLine};

use valenx_cfd_native::{solve_simple, Boundaries, FlowSolution, Fluid, Grid, SimpleControls};

use crate::background::{BackgroundJob, JobState};
use crate::plot_ui::managed_plot_mem_cfg;
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
    /// A running background solve, polled each frame. While `Some`, the form is
    /// frozen and a spinner shows until the worker delivers its solution.
    job: Option<BackgroundJob<FlowSolution>>,
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
            job: None,
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

impl CfdWorkbenchState {
    /// The user-visible captions of every control the agent bridge can set via
    /// `SetControl` (see [`crate::agent_commands`]). Returned by `ListControls`;
    /// each string matches exactly the caption the form draws. The drive-speed
    /// control's caption is **case-dependent** (`lid speed U (m/s)` for the
    /// cavity, `inlet speed U (m/s)` for the channel) — both spellings address
    /// the same `speed` field in [`agent_set`](Self::agent_set); the neutral
    /// `drive speed U (m/s)` is listed here.
    pub fn agent_control_names() -> &'static [&'static str] {
        &[
            "nx",
            "ny",
            "Lx",
            "Ly",
            "density ρ (kg/m³)",
            "kinematic ν (m²/s)",
            "drive speed U (m/s)",
            "max SIMPLE iterations",
        ]
    }

    /// Set one labelled control by its user-visible caption, for the agent
    /// `SetControl` bridge. Fail-loud: an unknown caption or a value of the
    /// wrong type / out of range returns `Err(String)` — never a panic. Ranges
    /// mirror `validate_inputs`: grid cells `nx`/`ny >= 1`, domain `Lx`/`Ly`
    /// and fluid `density`/`viscosity` finite `> 0`, the SIMPLE iteration cap
    /// `>= 1`. The drive speed accepts any finite value (a negative drive is a
    /// valid reversed lid/inlet). Both case-dependent speed captions
    /// (`lid speed U (m/s)` / `inlet speed U (m/s)`) plus the neutral
    /// `drive speed U (m/s)` map to the same field.
    pub fn agent_set(
        &mut self,
        name: &str,
        value: &crate::agent_commands::AgentValue,
    ) -> Result<(), String> {
        // A finite, strictly-positive real (domain / fluid props).
        let positive = |v: f64, what: &str| -> Result<f64, String> {
            if v.is_finite() && v > 0.0 {
                Ok(v)
            } else {
                Err(format!("{what} must be > 0, got {v}"))
            }
        };
        match name {
            "nx" => {
                let n = value.as_i64()?;
                if n < 1 {
                    return Err(format!("nx must be >= 1, got {n}"));
                }
                self.nx = n as usize;
            }
            "ny" => {
                let n = value.as_i64()?;
                if n < 1 {
                    return Err(format!("ny must be >= 1, got {n}"));
                }
                self.ny = n as usize;
            }
            "Lx" => self.lx = positive(value.as_f64()?, "Lx")?,
            "Ly" => self.ly = positive(value.as_f64()?, "Ly")?,
            "density ρ (kg/m³)" => self.density = positive(value.as_f64()?, "density")?,
            "kinematic ν (m²/s)" => self.viscosity = positive(value.as_f64()?, "viscosity")?,
            "drive speed U (m/s)" | "lid speed U (m/s)" | "inlet speed U (m/s)" => {
                let v = value.as_f64()?;
                if !v.is_finite() {
                    return Err(format!("drive speed U must be finite, got {v}"));
                }
                self.speed = v;
            }
            "max SIMPLE iterations" => {
                let n = value.as_i64()?;
                if n < 1 {
                    return Err(format!("max SIMPLE iterations must be >= 1, got {n}"));
                }
                self.max_iterations = n as usize;
            }
            other => return Err(format!("unknown CFD control: {other:?}")),
        }
        Ok(())
    }

    /// The current computed-result text for the agent `ReadReadout` bridge (see
    /// [`crate::agent_commands`]): the same `Result` string the panel renders
    /// when a solve has produced one, else the last `error`, else `None` when the
    /// case has not been run yet. Read-only — closes the live-driving loop by
    /// letting an agent read the answer back after a `RunCommand`/solve.
    pub fn agent_readout(&self) -> Option<String> {
        if !self.result.is_empty() {
            Some(self.result.clone())
        } else {
            self.error.clone()
        }
    }
}

/// Draw the CFD Workbench right-side panel. A no-op when the
/// `show_cfd_workbench` toggle is off.
pub fn draw_cfd_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_cfd_workbench {
        return;
    }
    poll_cfd(&mut app.cfd);
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_cfd_workbench",
        "CFD Workbench",
        cfd_workbench_body,
    );
    if close {
        app.show_cfd_workbench = false;
    }
}

/// The CFD workbench body — case picker, grid, fluid props, run controls,
/// residual/convergence plots, and the centreline-profile plot. Extracted
/// from [`draw_cfd_workbench`] so it can be hosted by the classic
/// [`crate::workbench_chrome::workbench_shell`] *or* the opt-in dockable
/// tile layout ([`crate::dock_layout`]). Re-runs the non-blocking
/// background-job poll up front so the dock path stays live too.
pub(crate) fn cfd_workbench_body(app: &mut ValenxApp, ui: &mut egui::Ui) {
    poll_cfd(&mut app.cfd);
    ui.label(
        egui::RichText::new("native 2-D incompressible CFD · valenx-cfd-native")
            .weak()
            .small(),
    );
    ui.separator();
    let s = &mut app.cfd;
    let running = s.job.is_some();
    if running {
        ui.ctx().request_repaint();
    }
    egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    if running {
                        ui.horizontal(|ui| {
                            ui.spinner();
                            ui.label("solving…");
                        });
                        ui.disable();
                    }
                    ui.label(egui::RichText::new("Case").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.case, CfdCase::LidDrivenCavity, "Lid-driven cavity")
                            .on_hover_text("Square cavity, sliding top lid — the recirculation benchmark.");
                        ui.radio_value(&mut s.case, CfdCase::ChannelFlow, "Channel flow")
                            .on_hover_text("Uniform inlet developing a parabolic profile.");
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Grid — staggered finite volume").strong());
                    // Associate each numeric `DragValue` with its caption via
                    // `labelled_by` so the spin button carries the caption as its
                    // accessibility / UI-Automation Name (egui clears a
                    // DragValue's own Name, so without this it is anonymous to a
                    // screen reader / AI driver). The hover text mirrors the
                    // caption for a mouse user.
                    ui.horizontal(|ui| {
                        let nx = ui.label("nx");
                        ui.add(egui::DragValue::new(&mut s.nx).speed(0.5))
                            .labelled_by(nx.id)
                            .on_hover_text("Grid cells along x");
                        let ny = ui.label("ny");
                        ui.add(egui::DragValue::new(&mut s.ny).speed(0.5))
                            .labelled_by(ny.id)
                            .on_hover_text("Grid cells along y");
                    });
                    ui.horizontal(|ui| {
                        let lx = ui.label("Lx");
                        ui.add(egui::DragValue::new(&mut s.lx).speed(0.05))
                            .labelled_by(lx.id)
                            .on_hover_text("Domain length Lx (m)");
                        let ly = ui.label("Ly");
                        ui.add(egui::DragValue::new(&mut s.ly).speed(0.05))
                            .labelled_by(ly.id)
                            .on_hover_text("Domain height Ly (m)");
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Fluid").strong());
                    ui.horizontal(|ui| {
                        let rho = ui.label("density ρ (kg/m³)");
                        ui.add(egui::DragValue::new(&mut s.density).speed(0.1))
                            .labelled_by(rho.id)
                            .on_hover_text("Fluid density ρ (kg/m³)");
                    });
                    ui.horizontal(|ui| {
                        let nu = ui.label("kinematic ν (m²/s)");
                        ui.add(egui::DragValue::new(&mut s.viscosity).speed(0.001))
                            .labelled_by(nu.id)
                            .on_hover_text("Kinematic viscosity ν (m²/s)");
                    });
                    let drive = match s.case {
                        CfdCase::LidDrivenCavity => "lid speed U (m/s)",
                        CfdCase::ChannelFlow => "inlet speed U (m/s)",
                    };
                    ui.horizontal(|ui| {
                        let u = ui.label(drive);
                        ui.add(egui::DragValue::new(&mut s.speed).speed(0.05))
                            .labelled_by(u.id)
                            .on_hover_text("Drive speed U (m/s)");
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
                        let it = ui.label("max SIMPLE iterations");
                        ui.add(egui::DragValue::new(&mut s.max_iterations).speed(5.0))
                            .labelled_by(it.id)
                            .on_hover_text("Outer SIMPLE iteration cap");
                    });
                    ui.label(
                        egui::RichText::new("Runs on a background thread — the UI stays responsive; a finer grid or more iterations takes longer.")
                            .weak()
                            .small(),
                    );

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Solve flow").strong())
                        .clicked()
                    {
                        s.error = None;
                        s.profile = None;
                        s.analytic_profile = None;
                        s.bulk_velocity = None;
                        s.result.clear();
                        match validate_inputs(s) {
                            Ok(inp) => {
                                s.job = Some(BackgroundJob::spawn(move || solve(inp)));
                            }
                            Err(e) => s.error = Some(e),
                        }
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
                        managed_plot_mem_cfg(
                            ui,
                            "cfd_profile_plot",
                            160.0,
                            |plot| plot.x_axis_label("speed (m/s)").y_axis_label("y (m)"),
                            |pui| {
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
                            },
                        );
                    }
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

/// An owned snapshot of the solver inputs, moved into the worker thread.
#[derive(Clone, Copy)]
struct CfdInputs {
    nx: usize,
    ny: usize,
    lx: f64,
    ly: f64,
    density: f64,
    viscosity: f64,
    speed: f64,
    max_iterations: usize,
    case: CfdCase,
}

/// Validate the form and snapshot the inputs. Every guard protects a
/// `valenx-cfd-native` precondition — e.g. [`Grid::new`] *panics* on a zero
/// axis — so a bad input surfaces as an error string, never a crash.
fn validate_inputs(s: &CfdWorkbenchState) -> Result<CfdInputs, String> {
    if s.nx == 0 || s.ny == 0 {
        return Err("grid must have at least one cell per axis".into());
    }
    if !(s.lx > 0.0 && s.ly > 0.0 && s.lx.is_finite() && s.ly.is_finite()) {
        return Err("domain dimensions must be positive and finite".into());
    }
    if !(s.viscosity > 0.0 && s.viscosity.is_finite()) {
        return Err("kinematic viscosity must be positive".into());
    }
    if !(s.density > 0.0 && s.density.is_finite()) {
        return Err("density must be positive".into());
    }
    Ok(CfdInputs {
        nx: s.nx,
        ny: s.ny,
        lx: s.lx,
        ly: s.ly,
        density: s.density,
        viscosity: s.viscosity,
        speed: s.speed,
        max_iterations: s.max_iterations,
        case: s.case,
    })
}

/// Build the grid + BCs from validated inputs and run the SIMPLE solver. This
/// is the heavy step — the UI runs it on a worker thread (see
/// [`draw_cfd_workbench`]) so the window stays responsive while it solves.
fn solve(inp: CfdInputs) -> FlowSolution {
    let grid = Grid::new(inp.nx, inp.ny, inp.lx, inp.ly);
    let fluid = Fluid::new(inp.density, inp.viscosity);
    let bcs = match inp.case {
        CfdCase::LidDrivenCavity => Boundaries::lid_driven_cavity(inp.speed),
        CfdCase::ChannelFlow => Boundaries::channel_flow(inp.speed),
    };
    let controls = SimpleControls {
        max_iterations: inp.max_iterations,
        ..Default::default()
    };
    solve_simple(&grid, &fluid, &bcs, &controls)
}

/// Move a finished background solve into the panel state (or surface a worker
/// crash). Called once per frame at the top of [`draw_cfd_workbench`].
fn poll_cfd(s: &mut CfdWorkbenchState) {
    match s.job.as_mut().map(BackgroundJob::poll) {
        Some(JobState::Done(sol)) => {
            s.job = None;
            apply_solution(s, &sol);
        }
        Some(JobState::Failed) => {
            s.job = None;
            s.error = Some("the CFD solver thread stopped unexpectedly".into());
        }
        Some(JobState::Pending) | None => {}
    }
}

/// Synchronous validate → solve → apply. Kept for the headless tests; the UI
/// instead spawns [`solve`] on a worker thread and calls [`apply_solution`]
/// when the result arrives. Because the form is frozen while a solve runs,
/// `apply_solution` reads the inputs straight off `s` and they match `sol`.
#[cfg(test)]
fn run_cfd(s: &mut CfdWorkbenchState) {
    s.error = None;
    s.profile = None;
    s.analytic_profile = None;
    s.bulk_velocity = None;
    match validate_inputs(s) {
        Ok(inp) => apply_solution(s, &solve(inp)),
        Err(e) => s.error = Some(e),
    }
}

/// Populate the panel readout (result string, centreline profile, analytic
/// overlay, bulk-velocity reference) from a finished [`FlowSolution`]. Reads
/// the inputs off `s`, which the frozen form keeps consistent with `sol`.
fn apply_solution(s: &mut CfdWorkbenchState, sol: &FlowSolution) {
    let max_speed = sol.max_speed();
    let mean_speed = sol.mean_speed();
    let rms_speed = sol.rms_speed();
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
         max |u|    : {:.5} m/s  (mean {:.5}, rms {:.5}){peak_speed_loc}\n\
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
        rms_speed,
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

/// Build the **CFD** result card for the Workbench+Agent bridge — a DATA-ONLY
/// [`crate::WorkspaceProduct`] (`mesh: None`) whose `lines` are the genuine
/// lid-driven-cavity flow diagnostics for the canonical default case (a 32×32
/// unit cavity at Re = 100): iteration count / residual, peak and mean speed,
/// dynamic pressure, cell Reynolds number, vorticity, circulation, wall shear,
/// etc. Registered as the `"cfd"` producer in
/// [`crate::products_registry::lookup`]; the tile renders it as a text card, not
/// a 3-D view.
///
/// Unlike the live panel — which spawns the SIMPLE solve on a worker thread so
/// the window stays responsive — a registry builder must be synchronous, so this
/// runs the same validate → [`solve`] → [`apply_solution`] path inline on the
/// default state (the test-only `run_cfd` does the same). The default 32×32 / 500-
/// iteration cavity solve is bounded and modest (the same order of work as the
/// FEM cantilever builder), so it is acceptable on the bridge call. The default
/// inputs always validate; on the (canonically-unreachable) validation error the
/// card carries that message instead.
pub(crate) fn cfd_product() -> crate::WorkspaceProduct {
    let mut s = CfdWorkbenchState::default();
    let lines = match validate_inputs(&s) {
        Ok(inp) => {
            apply_solution(&mut s, &solve(inp));
            crate::products_registry::lines_from_readout(&s.result)
        }
        Err(e) => vec![format!("CFD setup invalid: {e}")],
    };
    crate::WorkspaceProduct {
        title: "CFD (lid-driven cavity)".into(),
        lines,
        mesh: None,
        vertex_colors: None,
        camera: Default::default(),
        kind2d: None,
        last_export: None,
        image: None,
        image_texture: None,
        animation: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::agent_commands::AgentValue;

    #[test]
    fn agent_set_sets_param_unknown_and_type_mismatch_err() {
        let mut s = CfdWorkbenchState::default();
        // A representative integer-grid set lands in state.
        s.agent_set("nx", &AgentValue::Int(48)).unwrap();
        assert_eq!(s.nx, 48);
        // A float field accepts a real value.
        s.agent_set("density ρ (kg/m³)", &AgentValue::Float(1.225))
            .unwrap();
        assert_eq!(s.density, 1.225);
        // Unknown caption -> Err (no panic).
        assert!(s.agent_set("no such control", &AgentValue::Int(1)).is_err());
        // Type mismatch (string into the integer grid count) -> Err.
        assert!(s.agent_set("nx", &AgentValue::Str("many".into())).is_err());
        // Out-of-range (zero cells) -> Err, field untouched.
        assert!(s.agent_set("nx", &AgentValue::Int(0)).is_err());
        assert_eq!(s.nx, 48, "rejected set leaves field untouched");
    }

    #[test]
    fn background_solve_populates_the_readout() {
        // Exercise the reactive path the UI actually uses: spawn `solve` on a
        // worker thread, then poll until `poll_cfd` applies the result. The
        // readout must come out the same as a synchronous run.
        let mut s = CfdWorkbenchState {
            nx: 16,
            ny: 16,
            max_iterations: 100,
            ..Default::default()
        };
        let inp = validate_inputs(&s).expect("default inputs are valid");
        s.job = Some(BackgroundJob::spawn(move || solve(inp)));
        // Bounded so the test can never hang on a stuck worker.
        for _ in 0..2000 {
            poll_cfd(&mut s);
            if s.job.is_none() {
                break;
            }
            std::thread::sleep(std::time::Duration::from_millis(1));
        }
        assert!(
            s.job.is_none(),
            "the background solve should finish within the poll budget"
        );
        assert!(s.error.is_none(), "unexpected error: {:?}", s.error);
        assert!(
            s.result.contains("max |u|"),
            "readout populated: {}",
            s.result
        );
        assert!(
            s.profile.as_ref().is_some_and(|p| p.len() == 16),
            "centreline profile sampled"
        );
    }

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
        assert!(
            s.result.contains("enstrophy"),
            "enstrophy in result: {}",
            s.result
        );
        // The cavity's primary vortex is rotation-dominated → the Q-criterion
        // readout is surfaced on the peak-vorticity line.
        assert!(
            s.result.contains("Q_max"),
            "Q-criterion in result: {}",
            s.result
        );
        // The sheared cavity flow also deforms → the peak strain-rate readout is
        // surfaced alongside (the strain companion to peak vorticity).
        assert!(
            s.result.contains("S_max"),
            "strain rate in result: {}",
            s.result
        );
        // The sheared cavity flow dissipates energy → the viscous-dissipation
        // readout is surfaced on the energy line.
        assert!(
            s.result.contains("\u{03A6}_visc"),
            "dissipation in result: {}",
            s.result
        );
        // The vorticity field is non-uniform → the palinstrophy readout (the
        // gradient-of-vorticity diagnostic) is surfaced on the cascade line.
        assert!(
            s.result.contains("palin"),
            "palinstrophy in result: {}",
            s.result
        );
        assert!(
            s.profile.as_ref().is_some_and(|p| p.len() == 16),
            "centreline profile sampled"
        );
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
        assert!(
            (max_u - 1.5).abs() < 0.02,
            "centreline ~1.5× inlet: {max_u}"
        );
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
            chan.analytic_profile
                .as_ref()
                .is_some_and(|p| p.len() == chan.ny),
            "channel flow should carry the analytic Poiseuille overlay"
        );
        // The lid-driven cavity has no 1-D analytic profile.
        let mut cav = CfdWorkbenchState {
            case: CfdCase::LidDrivenCavity,
            ..Default::default()
        };
        run_cfd(&mut cav);
        assert!(
            cav.analytic_profile.is_none(),
            "cavity has no analytic overlay"
        );
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
        assert_eq!(
            flow_regime(100.0, CfdCase::ChannelFlow),
            FlowRegime::Laminar
        );
        assert_eq!(
            flow_regime(3000.0, CfdCase::ChannelFlow),
            FlowRegime::Transitional
        );
        assert_eq!(
            flow_regime(5000.0, CfdCase::ChannelFlow),
            FlowRegime::Turbulent
        );
        // The lid-driven cavity stays laminar to much higher Re.
        assert_eq!(
            flow_regime(3000.0, CfdCase::LidDrivenCavity),
            FlowRegime::Laminar
        );
        assert_eq!(
            flow_regime(9000.0, CfdCase::LidDrivenCavity),
            FlowRegime::Transitional
        );
        // The default cavity (Re = 100) reports "laminar" in the solve result.
        let mut s = CfdWorkbenchState {
            nx: 12,
            ny: 12,
            max_iterations: 60,
            ..Default::default()
        };
        run_cfd(&mut s);
        assert!(
            s.result.contains("laminar"),
            "regime in result: {}",
            s.result
        );
    }

    #[test]
    fn dynamic_pressure_is_half_rho_u_squared() {
        assert!((dynamic_pressure(1.0, 10.0) - 50.0).abs() < 1e-12);
        // Quadratic in speed: doubling U quadruples q.
        assert!((dynamic_pressure(1.2, 20.0) - 4.0 * dynamic_pressure(1.2, 10.0)).abs() < 1e-9);
        // Linear in density.
        assert!((dynamic_pressure(2.0, 10.0) - 2.0 * dynamic_pressure(1.0, 10.0)).abs() < 1e-9);
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

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use egui::accesskit::{Node, NodeId, Role};

    /// Draw the panel once in a headless egui context **with accesskit enabled**
    /// and return the emitted accessibility tree nodes — the same tree a screen
    /// reader / AI UI-Automation driver consumes. `accesskit` is re-exported by
    /// egui, so this needs no extra dependency.
    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_cfd_workbench(app, ctx);
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
        assert!(!app.show_cfd_workbench);
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_cfd_workbench(&mut app, ctx);
        });
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_cfd_workbench = true;
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // Every DragValue (Role::SpinButton) must be associated with a caption
        // via `labelled_by` so it is findable by name; egui clears a DragValue's
        // own Name. The case radio buttons + Solve button carry their text Name.
        let mut app = ValenxApp::default();
        app.show_cfd_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        // nx, ny, Lx, Ly, density, viscosity, drive speed, max-iterations.
        assert!(
            spin_buttons.len() >= 8,
            "expected the CFD numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every CFD DragValue must be labelled_by its caption (AI-drivable name)"
        );

        for caption in ["nx", "Lx", "density ρ (kg/m³)", "max SIMPLE iterations"] {
            assert!(
                has_named_node(&nodes, caption),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
        // The case radio buttons are named, selectable nodes.
        assert!(
            has_named_node(&nodes, "Lid-driven cavity"),
            "the case radio buttons keep their text Name"
        );
    }
}
