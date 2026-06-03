//! The right-side **FEM Workbench** panel — native finite-element
//! analysis over `valenx-fem`'s in-process solvers (no external solver,
//! no input deck).
//!
//! Mirrors the other workbenches (`aero_workbench`, `astro_workbench`):
//! a resizable [`egui::SidePanel`] gated on
//! `crate::ValenxApp::show_fem_workbench`, toggled from the View menu.
//!
//! v1 surfaces the **linear-static** and **modal** solvers on a built-in
//! structured tetrahedral box mesh: set the box dimensions + material,
//! the `x = 0` face is fixed, and either a tip load is applied to the
//! `x = Lx` face (cantilever bending) or the natural frequencies are
//! extracted. This exercises the real native FEM end-to-end — the same
//! `valenx-fem` solvers that previously had no UI at all.

use eframe::egui;

use valenx_fem::material::FemMaterial;
use valenx_fem::modal_solver::solve_modal;
use valenx_fem::native_solver::{
    solve_linear_static, structured_box_mesh, NodalConstraint, NodalForce,
};

use crate::ValenxApp;

/// Which native solver the workbench runs.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
enum FemSolver {
    #[default]
    LinearStatic,
    Modal,
}

/// Persistent state for the FEM Workbench.
pub struct FemWorkbenchState {
    // Structured box mesh (metres + subdivisions).
    lx: f64,
    ly: f64,
    lz: f64,
    nx: usize,
    ny: usize,
    nz: usize,
    // Isotropic material.
    youngs_gpa: f64,
    poisson: f64,
    density: f64,
    // Linear-static tip load (newtons, applied downward in -Y).
    force_n: f64,
    // Modal: number of modes to extract.
    n_modes: usize,
    solver: FemSolver,
    result: String,
    error: Option<String>,
}

impl Default for FemWorkbenchState {
    fn default() -> Self {
        // A 1 m × 0.1 m × 0.1 m steel bar — a classic cantilever.
        Self {
            lx: 1.0,
            ly: 0.1,
            lz: 0.1,
            nx: 12,
            ny: 3,
            nz: 3,
            youngs_gpa: 205.0,
            poisson: 0.29,
            density: 7850.0,
            force_n: 1000.0,
            n_modes: 6,
            solver: FemSolver::LinearStatic,
            result: String::new(),
            error: None,
        }
    }
}

/// Draw the FEM Workbench right-side panel. A no-op when the
/// `show_fem_workbench` toggle is off.
pub fn draw_fem_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_fem_workbench {
        return;
    }
    egui::SidePanel::right("valenx_fem_workbench")
        .resizable(true)
        .default_width(360.0)
        .width_range(300.0..=560.0)
        .show(ctx, |ui| {
            ui.heading("FEM Workbench");
            ui.label(
                egui::RichText::new("native finite-element analysis · valenx-fem")
                    .weak()
                    .small(),
            );
            ui.separator();
            let s = &mut app.fem;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Geometry — structured box mesh (m)").strong());
                    ui.horizontal(|ui| {
                        ui.label("Lx");
                        ui.add(egui::DragValue::new(&mut s.lx).speed(0.05));
                        ui.label("Ly");
                        ui.add(egui::DragValue::new(&mut s.ly).speed(0.01));
                        ui.label("Lz");
                        ui.add(egui::DragValue::new(&mut s.lz).speed(0.01));
                    });
                    ui.horizontal(|ui| {
                        ui.label("nx");
                        ui.add(egui::DragValue::new(&mut s.nx).speed(0.2));
                        ui.label("ny");
                        ui.add(egui::DragValue::new(&mut s.ny).speed(0.2));
                        ui.label("nz");
                        ui.add(egui::DragValue::new(&mut s.nz).speed(0.2));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Material").strong());
                    ui.horizontal(|ui| {
                        ui.label("E (GPa)");
                        ui.add(egui::DragValue::new(&mut s.youngs_gpa).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("Poisson ν");
                        ui.add(egui::DragValue::new(&mut s.poisson).speed(0.005));
                    });
                    ui.horizontal(|ui| {
                        ui.label("Density (kg/m³)");
                        ui.add(egui::DragValue::new(&mut s.density).speed(10.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Analysis").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.solver, FemSolver::LinearStatic, "Linear static")
                            .on_hover_text("Cantilever bending under a tip load (Cholesky).");
                        ui.radio_value(&mut s.solver, FemSolver::Modal, "Modal")
                            .on_hover_text("Natural frequencies + mode shapes (K φ = λ M φ).");
                    });
                    match s.solver {
                        FemSolver::LinearStatic => {
                            ui.horizontal(|ui| {
                                ui.label("Tip load Fy (N, downward)");
                                ui.add(egui::DragValue::new(&mut s.force_n).speed(50.0));
                            });
                        }
                        FemSolver::Modal => {
                            ui.horizontal(|ui| {
                                ui.label("# modes");
                                ui.add(egui::DragValue::new(&mut s.n_modes).speed(0.2));
                            });
                        }
                    }
                    ui.label(
                        egui::RichText::new(
                            "The x=0 face is fixed; the load / modes are evaluated on the bar.",
                        )
                        .weak()
                        .small(),
                    );

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Run analysis").strong())
                        .clicked()
                    {
                        run_fem(s);
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
                });
        });
}

/// Build the box mesh + boundary conditions and run the selected native
/// solver. Extracted from the draw closure so it is unit-testable.
fn run_fem(s: &mut FemWorkbenchState) {
    s.error = None;
    let mesh = match structured_box_mesh(s.lx, s.ly, s.lz, s.nx, s.ny, s.nz) {
        Ok(m) => m,
        Err(e) => {
            s.error = Some(format!("mesh: {e}"));
            return;
        }
    };
    let material = FemMaterial {
        youngs_modulus: s.youngs_gpa * 1e9,
        poisson_ratio: s.poisson,
        density: s.density,
        ..Default::default()
    };

    // Fix the x = 0 face; collect the x = Lx face as the loaded tip.
    let tol = (s.lx / s.nx.max(1) as f64) * 1e-3 + 1e-9;
    let mut constraints = Vec::new();
    let mut tip_nodes = Vec::new();
    for (i, p) in mesh.nodes.iter().enumerate() {
        if p.x <= tol {
            constraints.push(NodalConstraint::fixed(i));
        } else if p.x >= s.lx - tol {
            tip_nodes.push(i);
        }
    }
    if constraints.is_empty() {
        s.error = Some("no nodes found on the fixed (x=0) face".into());
        return;
    }

    match s.solver {
        FemSolver::LinearStatic => {
            let per = if tip_nodes.is_empty() {
                0.0
            } else {
                -s.force_n / tip_nodes.len() as f64
            };
            let forces: Vec<NodalForce> = tip_nodes
                .iter()
                .map(|&n| NodalForce {
                    node: n,
                    force: [0.0, per, 0.0],
                })
                .collect();
            match solve_linear_static(&mesh, &material, &constraints, &forces) {
                Ok(sol) => {
                    let vm = sol.max_von_mises();
                    s.result = format!(
                        "Linear static  ({} nodes, {} fixed)\n\
                         tip load        : {:.1} N downward\n\
                         max displacement: {:.6e} m\n\
                         max von Mises   : {:.4e} Pa  ({:.3} MPa)",
                        mesh.nodes.len(),
                        constraints.len(),
                        s.force_n,
                        sol.max_displacement(),
                        vm,
                        vm / 1e6,
                    );
                }
                Err(e) => s.error = Some(format!("solve: {e}")),
            }
        }
        FemSolver::Modal => {
            let n_modes = s.n_modes.max(1);
            match solve_modal(&mesh, &material, &constraints, n_modes) {
                Ok(sol) => {
                    let mut out = format!(
                        "Modal  ({} nodes, {} fixed)\nnatural frequencies:\n",
                        mesh.nodes.len(),
                        constraints.len()
                    );
                    for (i, m) in sol.modes.iter().enumerate() {
                        out.push_str(&format!("  mode {:>2}: {:>12.4} Hz\n", i + 1, m.frequency_hz));
                    }
                    s.result = out;
                }
                Err(e) => s.error = Some(format!("solve: {e}")),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn linear_static_runs_on_default_box() {
        let mut s = FemWorkbenchState {
            solver: FemSolver::LinearStatic,
            ..Default::default()
        };
        run_fem(&mut s);
        assert!(s.error.is_none(), "unexpected error: {:?}", s.error);
        assert!(s.result.contains("max displacement"));
    }

    #[test]
    fn modal_runs_on_default_box() {
        let mut s = FemWorkbenchState {
            solver: FemSolver::Modal,
            ..Default::default()
        };
        run_fem(&mut s);
        assert!(s.error.is_none(), "unexpected error: {:?}", s.error);
        assert!(s.result.contains("Hz"));
    }

    #[test]
    fn degenerate_mesh_fails_loud() {
        let mut s = FemWorkbenchState {
            nx: 0,
            ..Default::default()
        };
        run_fem(&mut s);
        assert!(s.error.is_some(), "nx=0 must surface an error, not panic");
    }
}
