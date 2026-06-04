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
use egui_plot::{Line, Plot, PlotPoints};

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
                        ui.label(
                            egui::RichText::new(format!("Reynolds number ≈ {re:.1}"))
                                .weak()
                                .small(),
                        );
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
                            });
                    }
                });
        });
}

/// Build the grid + BCs and run the SIMPLE solver. Extracted from the
/// draw closure so it is unit-testable, and it validates every input
/// before calling [`Grid::new`] (which *panics* on a bad grid) — so a
/// bad input surfaces as an error, never a crash.
fn run_cfd(s: &mut CfdWorkbenchState) {
    s.error = None;
    s.profile = None;
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

    let mut max_speed = 0.0_f64;
    for j in 0..s.ny {
        for i in 0..s.nx {
            max_speed = max_speed.max(sol.speed_at_cell(i, j));
        }
    }
    let re = s.speed.abs() * characteristic_length(s) / s.viscosity;

    // Vertical centreline velocity profile: speed vs height.
    let i_mid = s.nx / 2;
    let profile: Vec<[f64; 2]> = (0..s.ny)
        .map(|j| {
            let y = (j as f64 + 0.5) * s.ly / s.ny as f64;
            [sol.speed_at_cell(i_mid, j), y]
        })
        .collect();
    s.profile = Some(profile);

    s.result = format!(
        "case       : {}\n\
         grid       : {}×{}  ({:.3} × {:.3} m)\n\
         Reynolds   : {:.1}\n\
         iterations : {} {}\n\
         residual   : {:.3e}\n\
         max |u|    : {:.5} m/s",
        s.case.label(),
        s.nx,
        s.ny,
        s.lx,
        s.ly,
        re,
        sol.iterations,
        if sol.converged {
            "(converged)"
        } else {
            "(hit iteration cap)"
        },
        sol.residual,
        max_speed,
    );
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
}
