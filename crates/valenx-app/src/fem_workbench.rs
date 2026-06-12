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
//! extracted. Results are shown as text, a plot (load–displacement line
//! or frequency spectrum), **and** — for the static case — the
//! **deformed shape coloured by von Mises stress** in the central 3-D
//! viewport via the `(Mesh, Field)` colour-ramp overlay.

use eframe::egui;
use egui_plot::{Line, Plot, PlotPoints, Points};
use nalgebra::Vector3;

use valenx_fem::material::FemMaterial;
use valenx_fem::modal_solver::solve_modal;
use valenx_fem::native_solver::{
    solve_linear_static, structured_box_mesh, NodalConstraint, NodalForce,
};
use valenx_fields::{Field, FieldKind, Location, RegionRef, TimeKey};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Which native solver the workbench runs.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
enum FemSolver {
    #[default]
    LinearStatic,
    Modal,
}

/// Result data to plot from the most recent run.
enum FemPlot {
    /// Natural frequencies (Hz); index + 1 is the mode number.
    Modal(Vec<f64>),
    /// Tip load (N) vs maximum displacement (m).
    LoadDisp(Vec<[f64; 2]>),
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
    /// Yield strength (MPa) for the factor-of-safety readout.
    yield_mpa: f64,
    // Linear-static tip load (newtons, applied downward in -Y).
    force_n: f64,
    // Modal: number of modes to extract.
    n_modes: usize,
    solver: FemSolver,
    result: String,
    /// Factor of safety (σy / peak von-Mises) from the last static run.
    fos: Option<f64>,
    /// Structural mass (kg) = density × box volume, from the last run.
    mass_kg: Option<f64>,
    /// Tip stiffness (N/m) = load / max displacement, from the last static run.
    stiffness_n_per_m: Option<f64>,
    /// Elastic strain energy `U = ½·Σ F·d` (J) from the last static run.
    strain_energy_j: Option<f64>,
    /// Serviceability deflection ratio `L/δ` (span over tip deflection).
    deflection_ratio: Option<f64>,
    error: Option<String>,
    plot: Option<FemPlot>,
    /// Deformed mesh + von-Mises field, pending a push to the 3-D viewport.
    viz: Option<(valenx_mesh::Mesh, Field)>,
    push_viz: bool,
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
            yield_mpa: 250.0,
            force_n: 1000.0,
            n_modes: 6,
            solver: FemSolver::LinearStatic,
            result: String::new(),
            fos: None,
            mass_kg: None,
            stiffness_n_per_m: None,
            strain_energy_j: None,
            deflection_ratio: None,
            error: None,
            plot: None,
            viz: None,
            push_viz: false,
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
                        // Clamp the mesh resolution: node count grows as
                        // (nx+1)(ny+1)(nz+1) and the solve as its cube, so an
                        // unbounded drag could hang the app for minutes / OOM.
                        ui.label("nx");
                        ui.add(egui::DragValue::new(&mut s.nx).speed(0.2).range(1..=40));
                        ui.label("ny");
                        ui.add(egui::DragValue::new(&mut s.ny).speed(0.2).range(1..=20));
                        ui.label("nz");
                        ui.add(egui::DragValue::new(&mut s.nz).speed(0.2).range(1..=20));
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
                    ui.horizontal(|ui| {
                        ui.label("Yield σy (MPa)");
                        ui.add(egui::DragValue::new(&mut s.yield_mpa).speed(5.0));
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
                            "x=0 fixed. Static runs colour the deformed shape by von Mises in the 3D viewport.",
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
                        if let Some(m) = s.mass_kg {
                            ui.label(
                                egui::RichText::new(format!(
                                    "structural mass: {m:.2} kg  (weight {:.1} N)",
                                    m * 9.80665
                                ))
                                .monospace()
                                .small(),
                            );
                        }
                        if let Some(f) = s.fos {
                            let (txt, col) = if f >= 1.0 {
                                (
                                    format!("✔ factor of safety {f:.2} — within yield"),
                                    egui::Color32::from_rgb(80, 220, 120),
                                )
                            } else {
                                (
                                    format!("✖ factor of safety {f:.2} — exceeds yield"),
                                    egui::Color32::from_rgb(220, 90, 90),
                                )
                            };
                            ui.colored_label(col, txt);
                        }
                    }

                    ui.add_space(6.0);
                    egui::CollapsingHeader::new("Closed-form beam reference")
                        .default_open(false)
                        .show(ui, |ui| {
                            // Reactive: recomputed every frame from the tip-load and length
                            // inputs (P = tip load N, L = box length Lx m). Independent of the
                            // FEA solve — a textbook cross-check for the current load/span.
                            let p = s.force_n;
                            let l = s.lx;
                            let w = if l > 0.0 { p / l } else { 0.0 };
                            ui.label(
                                egui::RichText::new(
                                    "P = tip load (N), L = box length Lx (m); w = P/L (N/m)",
                                )
                                .weak()
                                .small(),
                            );
                            let row = |ui: &mut egui::Ui, label: &str, val: f64, unit: &str| {
                                ui.label(
                                    egui::RichText::new(format!("  {label}: {val:.4} {unit}"))
                                        .monospace()
                                        .small(),
                                );
                            };
                            ui.label(egui::RichText::new("point load P").small().strong());
                            row(
                                ui,
                                "cantilever root moment P·L",
                                valenx_fem::cantilever_point_load_root_moment(p, l),
                                "N·m",
                            );
                            row(
                                ui,
                                "propped-cantilever prop reaction 5P/16",
                                valenx_fem::propped_cantilever_central_load_prop_reaction(p, l),
                                "N",
                            );
                            row(
                                ui,
                                "propped-cantilever fixed reaction 11P/16",
                                valenx_fem::propped_cantilever_central_load_fixed_end_reaction(p, l),
                                "N",
                            );
                            row(
                                ui,
                                "propped-cantilever clamp moment 3PL/16",
                                valenx_fem::propped_cantilever_central_load_fixed_end_moment(p, l),
                                "N·m",
                            );
                            row(
                                ui,
                                "two-span centre moment 3PL/32",
                                valenx_fem::two_span_continuous_beam_central_point_load_middle_moment(
                                    p, l,
                                ),
                                "N·m",
                            );
                            row(
                                ui,
                                "two-span centre reaction 11P/16",
                                valenx_fem::two_span_continuous_beam_central_point_load_middle_reaction(
                                    p, l,
                                ),
                                "N",
                            );
                            row(
                                ui,
                                "two-span loaded-span reaction 13P/32",
                                valenx_fem::two_span_continuous_beam_central_point_load_loaded_span_outer_reaction(
                                    p, l,
                                ),
                                "N",
                            );
                            row(
                                ui,
                                "two-span unloaded-span reaction −3P/32",
                                valenx_fem::two_span_continuous_beam_central_point_load_unloaded_span_outer_reaction(
                                    p, l,
                                ),
                                "N",
                            );
                            ui.add_space(3.0);
                            ui.label(egui::RichText::new("equivalent UDL w = P/L").small().strong());
                            row(
                                ui,
                                "propped-cantilever prop reaction 3wL/8",
                                valenx_fem::propped_cantilever_udl_prop_reaction(w, l),
                                "N",
                            );
                            row(
                                ui,
                                "propped-cantilever clamp moment wL²/8",
                                valenx_fem::propped_cantilever_udl_fixed_end_moment(w, l),
                                "N·m",
                            );
                            row(
                                ui,
                                "two-span centre moment wL²/8",
                                valenx_fem::two_span_continuous_beam_udl_middle_moment(w, l),
                                "N·m",
                            );
                            row(
                                ui,
                                "two-span centre reaction 5wL/4",
                                valenx_fem::two_span_continuous_beam_udl_middle_reaction(w, l),
                                "N",
                            );
                            ui.add_space(3.0);
                            let e = s.youngs_gpa * 1e9;
                            let i_sec = s.ly * s.lz.powi(3) / 12.0;
                            ui.label(
                                egui::RichText::new(format!(
                                    "deflection (solid rect. section I = ly·lz³/12 = {i_sec:.4e} m⁴)"
                                ))
                                .small()
                                .strong(),
                            );
                            row(
                                ui,
                                "cantilever tip δ = PL³/3EI",
                                valenx_fem::cantilever_tip_deflection(p, l, e, i_sec) * 1000.0,
                                "mm",
                            );
                            row(
                                ui,
                                "cantilever UDL tip δ = wL⁴/8EI",
                                valenx_fem::cantilever_udl_tip_deflection(w, l, e, i_sec) * 1000.0,
                                "mm",
                            );
                            row(
                                ui,
                                "simply-supported centre δ = PL³/48EI",
                                valenx_fem::simply_supported_center_deflection(p, l, e, i_sec)
                                    * 1000.0,
                                "mm",
                            );
                            row(
                                ui,
                                "simply-supported UDL centre δ = 5wL⁴/384EI",
                                valenx_fem::simply_supported_udl_center_deflection(w, l, e, i_sec)
                                    * 1000.0,
                                "mm",
                            );
                        });

                    ui.add_space(6.0);
                    egui::CollapsingHeader::new("Isotropic elastic constants")
                        .default_open(false)
                        .show(ui, |ui| {
                            // Reactive: recomputed every frame from Young's modulus E (GPa) and
                            // Poisson's ratio ν. K (bulk), G (shear), λ (Lamé first), M (P-wave).
                            let e = s.youngs_gpa * 1e9;
                            let nu = s.poisson;
                            ui.label(
                                egui::RichText::new(
                                    "E = Young's modulus (GPa), ν = Poisson's ratio; K, G, λ, M in GPa",
                                )
                                .weak()
                                .small(),
                            );
                            let row = |ui: &mut egui::Ui, label: &str, val: f64| {
                                ui.label(
                                    egui::RichText::new(format!("  {label}: {val:.4} GPa"))
                                        .monospace()
                                        .small(),
                                );
                            };
                            row(ui, "bulk modulus K", valenx_fem::bulk_modulus(e, nu) / 1e9);
                            row(
                                ui,
                                "shear modulus G",
                                valenx_fem::shear_modulus_from_youngs(e, nu) / 1e9,
                            );
                            row(
                                ui,
                                "Lamé first parameter λ",
                                valenx_fem::lames_first_parameter(e, nu) / 1e9,
                            );
                            row(ui, "P-wave modulus M", valenx_fem::p_wave_modulus(e, nu) / 1e9);
                            ui.label(egui::RichText::new("M = K + 4G/3 = λ + 2G").weak().small());
                        });

                    ui.add_space(6.0);
                    egui::CollapsingHeader::new("Section & buckling")
                        .default_open(false)
                        .show(ui, |ui| {
                            // Reactive: from the solid box section (ly × lz), length Lx, and E.
                            let e = s.youngs_gpa * 1e9;
                            let area = s.ly * s.lz;
                            let i_sec = valenx_fem::rectangular_second_moment_of_area(s.ly, s.lz);
                            let l = s.lx;
                            let c = s.lz / 2.0;
                            let r_gyr = valenx_fem::section_radius_of_gyration(i_sec, area);
                            ui.label(
                                egui::RichText::new(
                                    "solid rect. section ly×lz, length Lx, pinned ends (K=1)",
                                )
                                .weak()
                                .small(),
                            );
                            let row = |ui: &mut egui::Ui, label: &str, val: f64, unit: &str| {
                                ui.label(
                                    egui::RichText::new(format!("  {label}: {val:.4} {unit}"))
                                        .monospace()
                                        .small(),
                                );
                            };
                            row(ui, "area A", area, "m²");
                            row(ui, "second moment I = ly·lz³/12", i_sec, "m⁴");
                            row(
                                ui,
                                "polar second moment J",
                                valenx_fem::rectangular_polar_second_moment_of_area(s.ly, s.lz),
                                "m⁴",
                            );
                            row(
                                ui,
                                "elastic section modulus Z = I/c",
                                valenx_fem::elastic_section_modulus(i_sec, c),
                                "m³",
                            );
                            row(ui, "radius of gyration r = √(I/A)", r_gyr, "m");
                            row(
                                ui,
                                "Euler critical load P_cr",
                                valenx_fem::euler_critical_load(e, i_sec, l, 1.0) / 1e3,
                                "kN",
                            );
                            row(
                                ui,
                                "slenderness λ = KL/r",
                                valenx_fem::slenderness_ratio(l, 1.0, r_gyr),
                                "—",
                            );
                            row(
                                ui,
                                "critical buckling stress σ_cr",
                                valenx_fem::critical_buckling_stress(e, i_sec, l, 1.0, area) / 1e6,
                                "MPa",
                            );
                            let g = valenx_fem::shear_modulus_from_youngs(e, s.poisson);
                            let j_polar =
                                valenx_fem::rectangular_polar_second_moment_of_area(s.ly, s.lz);
                            row(
                                ui,
                                "flexural rigidity E·I",
                                valenx_fem::flexural_rigidity(e, i_sec),
                                "N·m²",
                            );
                            row(
                                ui,
                                "torsional rigidity G·J",
                                valenx_fem::torsional_rigidity(g, j_polar),
                                "N·m²",
                            );
                            let tip_moment =
                                valenx_fem::cantilever_point_load_root_moment(s.force_n, l);
                            row(
                                ui,
                                "cantilever root bending stress σ = M·c/I",
                                valenx_fem::bending_stress(tip_moment, c, i_sec) / 1e6,
                                "MPa",
                            );
                            let q_na = s.ly * s.lz * s.lz / 8.0;
                            row(
                                ui,
                                "cantilever max shear stress τ = VQ/Ib",
                                valenx_fem::beam_transverse_shear_stress(s.force_n, q_na, i_sec, s.ly)
                                    / 1e6,
                                "MPa",
                            );
                        });

                    if let Some(plot) = &s.plot {
                        ui.add_space(4.0);
                        match plot {
                            FemPlot::Modal(freqs) => {
                                ui.label(egui::RichText::new("Natural frequencies").strong());
                                Plot::new("fem_modal_plot").height(150.0).show(ui, |pui| {
                                    let pts: Vec<[f64; 2]> = freqs
                                        .iter()
                                        .enumerate()
                                        .map(|(i, &f)| [(i + 1) as f64, f])
                                        .collect();
                                    pui.line(Line::new(PlotPoints::from(pts.clone())).name("Hz"));
                                    pui.points(Points::new(PlotPoints::from(pts)).radius(3.0));
                                });
                            }
                            FemPlot::LoadDisp(pts) => {
                                ui.label(egui::RichText::new("Load – displacement").strong());
                                Plot::new("fem_loaddisp_plot").height(150.0).show(ui, |pui| {
                                    pui.line(Line::new(PlotPoints::from(pts.clone())).name("tip"));
                                });
                            }
                        }
                    }
                });
        });

    // Deferred (outside the panel borrow): push the deformed-shape field
    // overlay into the central 3-D viewport.
    if app.fem.push_viz {
        app.fem.push_viz = false;
        if let Some((mesh, field)) = app.fem.viz.take() {
            let quality = valenx_mesh::quality_report(&mesh);
            let aspect_hist =
                valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
            let skew_hist =
                valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
            app.stl = None;
            app.mesh = Some(LoadedMesh {
                path: std::path::PathBuf::from("<fem>/deformed"),
                mesh,
                quality,
                aspect_hist,
                skew_hist,
            });
            app.aero_field_overlay = Some(field);
            app.frame_current_mesh();
        }
    }
}

/// Build the box mesh + boundary conditions and run the selected native
/// solver. Extracted from the draw closure so it is unit-testable.
/// Elastic strain energy `U = ½·Σ F·d` (joules). By Clapeyron's theorem the
/// work the applied loads do on a linear elastic body equals the energy it
/// stores, so summing `force · displacement` over the loaded DOFs (with the
/// solved nodal displacements) is the exact strain energy — not the `½·F·δ_max`
/// single-point approximation. Force node indices past the end of the
/// displacement field are skipped.
fn strain_energy_j(forces: &[NodalForce], displacement: &[[f64; 3]]) -> f64 {
    0.5 * forces
        .iter()
        .filter_map(|f| {
            displacement
                .get(f.node)
                .map(|d| f.force[0] * d[0] + f.force[1] * d[1] + f.force[2] * d[2])
        })
        .sum::<f64>()
}

/// Serviceability deflection ratio — the span divided by the deflection, i.e.
/// the `n` in the familiar `L/n` limit (a cantilever is typically held to
/// `L/180`, a floor beam to `L/360`). `None` for a non-positive deflection.
fn span_deflection_ratio(span_m: f64, deflection_m: f64) -> Option<f64> {
    if deflection_m > 0.0 {
        Some(span_m / deflection_m)
    } else {
        None
    }
}

fn run_fem(s: &mut FemWorkbenchState) {
    s.error = None;
    s.plot = None;
    s.viz = None;
    s.push_viz = false;
    s.fos = None;
    s.mass_kg = None;
    s.stiffness_n_per_m = None;
    s.strain_energy_j = None;
    s.deflection_ratio = None;
    let mesh = match structured_box_mesh(s.lx, s.ly, s.lz, s.nx, s.ny, s.nz) {
        Ok(m) => m,
        Err(e) => {
            s.error = Some(format!("mesh: {e}"));
            return;
        }
    };
    // Structural mass of the solid box = density × volume.
    s.mass_kg = Some(s.density * s.lx * s.ly * s.lz);
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
                    let max_principal = sol.max_principal_stress();
                    let min_principal = sol.min_principal_stress();
                    let max_shear = sol.max_shear_stress();
                    // Hydrostatic (mean / volumetric) stress extremes — the
                    // pressure-like part complementary to the deviatoric measures.
                    let max_hydro = sol.max_hydrostatic_stress();
                    let min_hydro = sol.min_hydrostatic_stress();
                    // Stress triaxiality σ_m/σ_vm at the peak-von-Mises node.
                    let triax = sol.peak_triaxiality();
                    // Mean von Mises + stress-concentration factor Kt = peak/mean.
                    let mean_vm = sol.mean_von_mises();
                    let kt = if mean_vm > 0.0 { vm / mean_vm } else { 0.0 };
                    // Coordinates of the peak-von-Mises node (where failure initiates).
                    let peak_str = sol
                        .peak_von_mises_index()
                        .and_then(|i| mesh.nodes.get(i))
                        .map(|p| format!("({:.3}, {:.3}, {:.3}) m", p[0], p[1], p[2]))
                        .unwrap_or_else(|| "n/a".to_string());
                    let max_disp = sol.max_displacement();
                    // Factor of safety = yield strength / peak von-Mises stress.
                    let fos = if vm > 0.0 {
                        Some(s.yield_mpa * 1e6 / vm)
                    } else {
                        None
                    };
                    s.fos = fos;
                    let fos_str = match fos {
                        Some(f) => format!("{f:.2}"),
                        None => "n/a".to_string(),
                    };
                    // Tip stiffness k = F / δ (N/m).
                    let stiffness = if max_disp > 0.0 {
                        Some(s.force_n / max_disp)
                    } else {
                        None
                    };
                    s.stiffness_n_per_m = stiffness;
                    let stiffness_str = match stiffness {
                        Some(k) => format!("{k:.4e}"),
                        None => "n/a".to_string(),
                    };
                    // Elastic strain energy U = ½·Σ F·d (Clapeyron), exact over
                    // the loaded tip DOFs.
                    let energy = strain_energy_j(&forces, &sol.displacement);
                    s.strain_energy_j = Some(energy);
                    // Serviceability deflection ratio L/δ (span over tip deflection).
                    let defl_ratio = span_deflection_ratio(s.lx, max_disp);
                    s.deflection_ratio = defl_ratio;
                    let defl_str = match defl_ratio {
                        Some(n) => format!("L/{n:.0}"),
                        None => "—".to_string(),
                    };
                    s.result = format!(
                        "Linear static  ({} nodes, {} fixed)\n\
                         tip load        : {:.1} N downward\n\
                         max displacement: {:.6e} m\n\
                         deflection ratio: {}  (span/δ)\n\
                         max von Mises   : {:.4e} Pa  ({:.3} MPa, triax {:.2}, Lode {:.2})\n\
                         mean von Mises  : {:.4e} Pa  (Kt {:.1} = peak/mean)\n\
                         max principal   : {:.4e} Pa  (min {:.4e} Pa)\n\
                         max shear       : {:.4e} Pa  (Tresca)\n\
                         hydrostatic     : {:.4e} Pa  (min {:.4e} Pa)\n\
                         peak stress @   : {}\n\
                         tip stiffness   : {} N/m\n\
                         strain energy   : {:.4e} J\n\
                         factor of safety: {} (σy = {:.0} MPa)",
                        mesh.nodes.len(),
                        constraints.len(),
                        s.force_n,
                        max_disp,
                        defl_str,
                        vm,
                        vm / 1e6,
                        triax,
                        sol.peak_lode_parameter(),
                        mean_vm,
                        kt,
                        max_principal,
                        min_principal,
                        max_shear,
                        max_hydro,
                        min_hydro,
                        peak_str,
                        stiffness_str,
                        energy,
                        fos_str,
                        s.yield_mpa,
                    );
                    // Linear FEM → displacement scales linearly with load.
                    let pts = (0..=5)
                        .map(|i| {
                            let f = i as f64 / 5.0;
                            [s.force_n * f, max_disp * f]
                        })
                        .collect();
                    s.plot = Some(FemPlot::LoadDisp(pts));

                    // Deformed shape (scaled for visibility), coloured by von Mises.
                    let scale = if max_disp > 1e-12 {
                        0.1 * s.lx / max_disp
                    } else {
                        0.0
                    };
                    let mut deformed = mesh.clone();
                    for (node, d) in deformed.nodes.iter_mut().zip(&sol.displacement) {
                        *node += Vector3::new(d[0], d[1], d[2]) * scale;
                    }
                    deformed.recompute_stats();
                    let mut field = Field {
                        name: "von Mises".to_string(),
                        kind: FieldKind::Scalar,
                        location: Location::OnNode,
                        region: RegionRef("fem".to_string()),
                        units: valenx_fields::units::PASCAL,
                        time: TimeKey::Steady,
                        data: sol.von_mises.clone(),
                        range: None,
                    };
                    field.recompute_range();
                    s.viz = Some((deformed, field));
                    s.push_viz = true;
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
                        let axis = ["X", "Y", "Z"][m.dominant_translation_axis()];
                        out.push_str(&format!(
                            "  mode {:>2}: {:>12.4} Hz   (T = {:.3} ms, dom {axis})\n",
                            i + 1,
                            m.frequency_hz,
                            m.period_s() * 1000.0,
                        ));
                    }
                    s.result = out;
                    s.plot = Some(FemPlot::Modal(
                        sol.modes.iter().map(|m| m.frequency_hz).collect(),
                    ));
                    // Visualise the fundamental mode shape: deform the mesh by
                    // modes[0].shape (mass-normalised, so scaled for visibility),
                    // coloured by per-node modal amplitude.
                    if let Some(mode) = sol.modes.first() {
                        let amp = mode.max_amplitude();
                        let scale = if amp > 1e-12 { 0.1 * s.lx / amp } else { 0.0 };
                        let mut deformed = mesh.clone();
                        for (node, d) in deformed.nodes.iter_mut().zip(&mode.shape) {
                            *node += Vector3::new(d[0], d[1], d[2]) * scale;
                        }
                        deformed.recompute_stats();
                        let data: Vec<f64> = mode
                            .shape
                            .iter()
                            .map(|d| (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt())
                            .collect();
                        let mut field = Field {
                            name: "mode amplitude".to_string(),
                            kind: FieldKind::Scalar,
                            location: Location::OnNode,
                            region: RegionRef("fem".to_string()),
                            units: valenx_fields::units::DIMENSIONLESS,
                            time: TimeKey::Steady,
                            data,
                            range: None,
                        };
                        field.recompute_range();
                        s.viz = Some((deformed, field));
                        s.push_viz = true;
                    }
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
        // The stress-state pair is surfaced: triaxiality (hydrostatic axis) and
        // its deviatoric companion, the Lode parameter, on the von-Mises line.
        assert!(
            s.result.contains("triax"),
            "triaxiality in result: {}",
            s.result
        );
        assert!(
            s.result.contains("Lode"),
            "Lode parameter in result: {}",
            s.result
        );
        assert!(
            matches!(s.plot, Some(FemPlot::LoadDisp(_))),
            "static run plots a curve"
        );
        // The deformed-shape overlay is built and queued for the viewport.
        assert!(s.push_viz, "static run queues the 3D viz");
        let (mesh, field) = s
            .viz
            .as_ref()
            .expect("static run builds the deformed mesh + field");
        assert_eq!(
            field.data.len(),
            mesh.nodes.len(),
            "one von-Mises value per node"
        );
    }

    #[test]
    fn linear_static_reports_factor_of_safety() {
        let mut s = FemWorkbenchState {
            solver: FemSolver::LinearStatic,
            ..Default::default()
        };
        run_fem(&mut s);
        assert!(s.error.is_none(), "unexpected error: {:?}", s.error);
        let f1 = s.fos.expect("static run computes a factor of safety");
        assert!(f1 > 0.0 && f1.is_finite());
        assert!(s.result.contains("factor of safety"));
        // FoS = σy / peak von-Mises; same load + geometry ⇒ same peak stress,
        // so halving the yield strength halves the factor of safety.
        s.yield_mpa *= 0.5;
        run_fem(&mut s);
        let f2 = s.fos.expect("FoS recomputed");
        assert!((f2 - 0.5 * f1).abs() / f1 < 1e-6, "FoS ∝ σy: {f1} → {f2}");
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
        match &s.plot {
            Some(FemPlot::Modal(freqs)) => assert_eq!(freqs.len(), 6, "six modes plotted"),
            other => panic!(
                "modal run should plot frequencies, got {:?}",
                other.is_some()
            ),
        }
        // The modal run now visualises the fundamental mode shape: a deformed
        // mesh coloured by per-node modal amplitude (one value per node).
        let (mesh, field) = s
            .viz
            .as_ref()
            .expect("modal run builds the fundamental mode-shape overlay");
        assert_eq!(field.name, "mode amplitude");
        assert_eq!(field.data.len(), mesh.nodes.len(), "one amplitude per node");
    }

    #[test]
    fn degenerate_mesh_fails_loud() {
        let mut s = FemWorkbenchState {
            nx: 0,
            ..Default::default()
        };
        run_fem(&mut s);
        assert!(s.error.is_some(), "nx=0 must surface an error, not panic");
        assert!(s.plot.is_none(), "a failed run clears the plot");
        assert!(!s.push_viz, "a failed run does not queue a viz");
    }

    #[test]
    fn run_reports_structural_mass() {
        let mut s = FemWorkbenchState::default();
        run_fem(&mut s);
        // mass = ρ·Lx·Ly·Lz = 7850 · 1 · 0.1 · 0.1 = 78.5 kg.
        let m = s.mass_kg.expect("mass computed on a successful run");
        assert!((m - 78.5).abs() < 1e-9, "mass = {m}");
        // A failed run (degenerate mesh) leaves mass cleared, not stale.
        let mut bad = FemWorkbenchState {
            nx: 0,
            ..Default::default()
        };
        run_fem(&mut bad);
        assert!(bad.mass_kg.is_none(), "a failed run clears the mass");
    }

    #[test]
    fn linear_static_reports_tip_stiffness() {
        let mut s = FemWorkbenchState {
            solver: FemSolver::LinearStatic,
            ..Default::default()
        };
        run_fem(&mut s);
        let k1 = s.stiffness_n_per_m.expect("static run computes stiffness");
        assert!(k1 > 0.0 && k1.is_finite());
        assert!(s.result.contains("tip stiffness"));
        // k = F/δ is a structural property: linear FEM ⇒ doubling the load
        // doubles δ, so k is unchanged.
        s.force_n *= 2.0;
        run_fem(&mut s);
        let k2 = s.stiffness_n_per_m.expect("stiffness recomputed");
        assert!(
            (k2 - k1).abs() / k1 < 1e-6,
            "stiffness is load-independent: {k1} → {k2}"
        );
    }

    #[test]
    fn strain_energy_sums_half_force_dot_displacement() {
        // Two loaded nodes; the unloaded node 1 must not contribute.
        let forces = vec![
            NodalForce {
                node: 0,
                force: [0.0, -10.0, 0.0],
            },
            NodalForce {
                node: 2,
                force: [4.0, 0.0, 0.0],
            },
        ];
        let disp = [[0.0, -2.0, 0.0], [9.0, 9.0, 9.0], [3.0, 0.0, 0.0]];
        // U = ½[(-10)(-2) + 4·3] = ½[20 + 12] = 16 J.
        assert!((strain_energy_j(&forces, &disp) - 16.0).abs() < 1e-12);
        // An out-of-range node index is skipped, not a panic.
        let bad = [NodalForce {
            node: 99,
            force: [1.0, 1.0, 1.0],
        }];
        assert_eq!(strain_energy_j(&bad, &disp), 0.0);
    }

    #[test]
    fn linear_static_reports_strain_energy() {
        let mut s = FemWorkbenchState {
            solver: FemSolver::LinearStatic,
            ..Default::default()
        };
        run_fem(&mut s);
        let u1 = s
            .strain_energy_j
            .expect("static run computes strain energy");
        assert!(
            u1 > 0.0 && u1.is_finite(),
            "a loaded body stores positive energy"
        );
        assert!(s.result.contains("strain energy"));
        // U = ½·Σ F·d with d ∝ F (linear FEM) ⇒ U ∝ F²: doubling the load
        // quadruples the stored energy.
        s.force_n *= 2.0;
        run_fem(&mut s);
        let u2 = s.strain_energy_j.expect("strain energy recomputed");
        assert!(
            (u2 / u1 - 4.0).abs() < 1e-3,
            "U scales with F²: {u1} → {u2}"
        );
        // A failed run clears it rather than leaving a stale value.
        let mut bad = FemWorkbenchState {
            nx: 0,
            ..Default::default()
        };
        run_fem(&mut bad);
        assert!(
            bad.strain_energy_j.is_none(),
            "a failed run clears the energy"
        );
    }

    #[test]
    fn span_deflection_ratio_and_serviceability_readout() {
        // The pure ratio is span / deflection (the n in L/n).
        assert!((span_deflection_ratio(2.0, 0.01).unwrap() - 200.0).abs() < 1e-9);
        // Non-positive deflection → no ratio.
        assert!(span_deflection_ratio(2.0, 0.0).is_none());

        // End to end: the cantilever reports a finite, positive ratio, and
        // because δ ∝ load (linear FEM) the ratio halves when the load doubles.
        let mut s = FemWorkbenchState {
            solver: FemSolver::LinearStatic,
            ..Default::default()
        };
        run_fem(&mut s);
        let r1 = s
            .deflection_ratio
            .expect("static run computes a deflection ratio");
        assert!(r1 > 0.0 && r1.is_finite());
        assert!(s.result.contains("deflection ratio"));
        s.force_n *= 2.0;
        run_fem(&mut s);
        let r2 = s.deflection_ratio.expect("ratio recomputed");
        assert!(
            (r2 - 0.5 * r1).abs() / r1 < 1e-6,
            "doubling the load halves L/δ: {r1} → {r2}"
        );
    }

    #[test]
    fn linear_static_reports_max_principal_stress() {
        let mut s = FemWorkbenchState {
            solver: FemSolver::LinearStatic,
            ..Default::default()
        };
        run_fem(&mut s);
        // The result surfaces the maximum principal (Rankine) stress alongside
        // the von Mises measure. (The eigenvalue maths is validated in
        // valenx-fem; here we just confirm the readout is wired through.)
        assert!(s.result.contains("max principal"), "result: {}", s.result);
        // The Tresca maximum shear stress is reported alongside it.
        assert!(s.result.contains("max shear"), "result: {}", s.result);
    }
}
