//! The right-side **Thermodynamic Cycle Workbench** panel — native
//! closed-form ideal-cycle efficiencies over `valenx-thermocycle`.
//!
//! Mirrors the Refrigeration / Heat Transfer workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_thermocycle_workbench`,
//! toggled from the View menu. A radio selects which cycle to analyze —
//! the reversible **Carnot** bound, the air-standard **Otto**, **Diesel**
//! and **Brayton / Joule** cycles, or a basic ideal **Rankine** steam
//! cycle from supplied enthalpies — and only that cycle's inputs are
//! shown. "Analyze" builds the selected cycle value type, reports its
//! ideal thermal efficiency (and, for Carnot, the reversed-cycle
//! refrigerator and heat-pump coefficients of performance; for Diesel, the
//! cutoff penalty factor), and "Show 3-D" loads a representative labeled
//! block into the central viewport.
//!
//! These are textbook closed-form air-standard / ideal-cycle models: real
//! `c_p`, `c_v` (hence `γ`) vary with temperature and every real cycle has
//! irreversibilities, so these are reversible-limit efficiencies, not
//! measured ones. The numbers reproduce the standard worked examples
//! exactly.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_thermocycle::{Brayton, Carnot, Diesel, HeatCapacityRatio, Otto, Rankine};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Which thermodynamic cycle the workbench analyzes. Selecting one shows
/// only that cycle's inputs.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum CycleType {
    /// The reversible Carnot cycle (the upper bound, plus reversed-cycle
    /// coefficients of performance).
    Carnot,
    /// The air-standard Otto (spark-ignition) cycle.
    Otto,
    /// The air-standard Diesel (compression-ignition) cycle.
    Diesel,
    /// The air-standard Brayton / Joule (gas-turbine) cycle.
    Brayton,
    /// The basic ideal Rankine (steam) cycle from supplied enthalpies.
    Rankine,
}

/// Persistent form + result state for the Thermodynamic Cycle Workbench.
pub struct ThermoCycleWorkbenchState {
    /// Which cycle is being analyzed.
    cycle: CycleType,
    /// Cold-reservoir absolute temperature `Tc` (K) — Carnot.
    t_cold_k: f64,
    /// Hot-reservoir absolute temperature `Th` (K) — Carnot.
    t_hot_k: f64,
    /// Compression ratio `r` — Otto / Diesel.
    compression_ratio: f64,
    /// Cutoff ratio `r_c = V3 / V2` — Diesel.
    cutoff_ratio: f64,
    /// Compressor pressure ratio `r_p` — Brayton.
    pressure_ratio: f64,
    /// Working-fluid heat-capacity ratio `γ = c_p / c_v` — gas cycles.
    gamma: f64,
    /// State-1 enthalpy `h1` (kJ/kg), condenser exit / pump inlet — Rankine.
    h1: f64,
    /// State-2 enthalpy `h2` (kJ/kg), pump exit / boiler inlet — Rankine.
    h2: f64,
    /// State-3 enthalpy `h3` (kJ/kg), boiler exit / turbine inlet — Rankine.
    h3: f64,
    /// State-4 enthalpy `h4` (kJ/kg), turbine exit / condenser inlet —
    /// Rankine.
    h4: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D block (serviced after the panel
    /// draws).
    show_3d_request: bool,
}

impl Default for ThermoCycleWorkbenchState {
    fn default() -> Self {
        // Realistic textbook defaults for each cycle:
        //   Carnot  300 K / 600 K          -> eta = 0.5 exactly.
        //   Otto    r = 8, gamma = 1.4     -> eta ~ 0.5647 (classic petrol).
        //   Diesel  r = 18, r_c = 2        -> eta ~ 0.6316 (Cengel example).
        //   Brayton r_p = 10, gamma = 1.4  -> eta ~ 0.4821 (classic turbine).
        //   Rankine Moran/Cengel enthalpies -> eta ~ 0.4043.
        Self {
            cycle: CycleType::Carnot,
            t_cold_k: 300.0,
            t_hot_k: 600.0,
            compression_ratio: 8.0,
            cutoff_ratio: 2.0,
            pressure_ratio: 10.0,
            gamma: HeatCapacityRatio::AIR,
            h1: 191.8,
            h2: 199.7,
            h3: 3247.6,
            h4: 2007.5,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Thermodynamic Cycle Workbench right-side panel. A no-op when
/// the `show_thermocycle_workbench` toggle is off.
pub fn draw_thermocycle_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_thermocycle_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_thermocycle_workbench",
        "Thermodynamic Cycle",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native closed-form ideal-cycle efficiency · valenx-thermocycle",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.thermocycle;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Cycle").strong());
                    ui.horizontal_wrapped(|ui| {
                        ui.radio_value(&mut s.cycle, CycleType::Carnot, "Carnot");
                        ui.radio_value(&mut s.cycle, CycleType::Otto, "Otto");
                        ui.radio_value(&mut s.cycle, CycleType::Diesel, "Diesel");
                        ui.radio_value(&mut s.cycle, CycleType::Brayton, "Brayton");
                        ui.radio_value(&mut s.cycle, CycleType::Rankine, "Rankine");
                    });

                    ui.add_space(4.0);
                    match s.cycle {
                        CycleType::Carnot => {
                            ui.label(egui::RichText::new("Reservoirs (K)").strong());
                            ui.horizontal(|ui| {
                                ui.label("cold Tc (K)");
                                ui.add(egui::DragValue::new(&mut s.t_cold_k).speed(1.0));
                            });
                            ui.horizontal(|ui| {
                                ui.label("hot Th (K)");
                                ui.add(egui::DragValue::new(&mut s.t_hot_k).speed(1.0));
                            });
                        }
                        CycleType::Otto => {
                            ui.label(egui::RichText::new("Air-standard Otto").strong());
                            ui.horizontal(|ui| {
                                ui.label("compression ratio r");
                                ui.add(egui::DragValue::new(&mut s.compression_ratio).speed(0.1));
                            });
                            ui.horizontal(|ui| {
                                ui.label("gamma γ");
                                ui.add(egui::DragValue::new(&mut s.gamma).speed(0.01));
                            });
                        }
                        CycleType::Diesel => {
                            ui.label(egui::RichText::new("Air-standard Diesel").strong());
                            ui.horizontal(|ui| {
                                ui.label("compression ratio r");
                                ui.add(egui::DragValue::new(&mut s.compression_ratio).speed(0.1));
                            });
                            ui.horizontal(|ui| {
                                ui.label("cutoff ratio rc");
                                ui.add(egui::DragValue::new(&mut s.cutoff_ratio).speed(0.05));
                            });
                            ui.horizontal(|ui| {
                                ui.label("gamma γ");
                                ui.add(egui::DragValue::new(&mut s.gamma).speed(0.01));
                            });
                        }
                        CycleType::Brayton => {
                            ui.label(egui::RichText::new("Air-standard Brayton").strong());
                            ui.horizontal(|ui| {
                                ui.label("pressure ratio rp");
                                ui.add(egui::DragValue::new(&mut s.pressure_ratio).speed(0.1));
                            });
                            ui.horizontal(|ui| {
                                ui.label("gamma γ");
                                ui.add(egui::DragValue::new(&mut s.gamma).speed(0.01));
                            });
                        }
                        CycleType::Rankine => {
                            ui.label(egui::RichText::new("Rankine enthalpies (kJ/kg)").strong());
                            ui.horizontal(|ui| {
                                ui.label("h1 cond. exit");
                                ui.add(egui::DragValue::new(&mut s.h1).speed(0.5));
                            });
                            ui.horizontal(|ui| {
                                ui.label("h2 pump exit");
                                ui.add(egui::DragValue::new(&mut s.h2).speed(0.5));
                            });
                            ui.horizontal(|ui| {
                                ui.label("h3 turbine inlet");
                                ui.add(egui::DragValue::new(&mut s.h3).speed(0.5));
                            });
                            ui.horizontal(|ui| {
                                ui.label("h4 turbine exit");
                                ui.add(egui::DragValue::new(&mut s.h4).speed(0.5));
                            });
                        }
                    }

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_thermocycle(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative labeled block tagged with the selected cycle as a 3-D solid and load it into the central viewport to orbit",
                        )
                        .clicked()
                    {
                        s.show_3d_request = true;
                    }

                    if let Some(e) = &s.error {
                        ui.add_space(4.0);
                        ui.colored_label(egui::Color32::from_rgb(220, 90, 90), e);
                    }

                    if !s.result.is_empty() {
                        ui.separator();
                        ui.label(egui::RichText::new("Efficiency").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_thermocycle_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.thermocycle` borrow is
    // released here): build the cycle's 3-D block and load it.
    if app.thermocycle.show_3d_request {
        app.thermocycle.show_3d_request = false;
        load_cycle_3d(app);
    }
}

/// Validate the form, evaluate the selected cycle and format the readout.
fn run_thermocycle(s: &mut ThermoCycleWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// A short label for the currently selected cycle, used in the readout.
fn cycle_name(s: &ThermoCycleWorkbenchState) -> &'static str {
    match s.cycle {
        CycleType::Carnot => "Carnot",
        CycleType::Otto => "Otto",
        CycleType::Diesel => "Diesel",
        CycleType::Brayton => "Brayton",
        CycleType::Rankine => "Rankine",
    }
}

/// Evaluate the selected cycle and format the full readout, mapping any
/// validation / domain error to a display string. Extracted so it is
/// unit-testable.
fn compute(s: &ThermoCycleWorkbenchState) -> Result<String, String> {
    let name = cycle_name(s);
    match s.cycle {
        CycleType::Carnot => {
            let c = Carnot::new(s.t_cold_k, s.t_hot_k).map_err(|e| e.to_string())?;
            let eta = c.efficiency();
            let cop_ref = c.cop_refrigerator();
            let cop_hp = c.cop_heat_pump();
            Ok(format!(
                "cycle           : {name}\n\
                 reservoirs Tc/Th: {:.2} / {:.2} K\n\n\
                 efficiency η    : {eta:.4}  ({:.2} %)\n\
                 COP refrigerator: {cop_ref:.3}\n\
                 COP heat pump   : {cop_hp:.3}",
                s.t_cold_k,
                s.t_hot_k,
                eta * 100.0,
            ))
        }
        CycleType::Otto => {
            let gamma = HeatCapacityRatio::new(s.gamma).map_err(|e| e.to_string())?;
            let o = Otto::new(s.compression_ratio, gamma).map_err(|e| e.to_string())?;
            let eta = o.efficiency();
            Ok(format!(
                "cycle           : {name}\n\
                 compression r   : {:.3}\n\
                 gamma γ         : {:.4}\n\n\
                 efficiency η    : {eta:.4}  ({:.2} %)",
                s.compression_ratio,
                s.gamma,
                eta * 100.0,
            ))
        }
        CycleType::Diesel => {
            let gamma = HeatCapacityRatio::new(s.gamma).map_err(|e| e.to_string())?;
            let d = Diesel::new(s.compression_ratio, s.cutoff_ratio, gamma)
                .map_err(|e| e.to_string())?;
            let eta = d.efficiency();
            let penalty = d.cutoff_penalty();
            Ok(format!(
                "cycle           : {name}\n\
                 compression r   : {:.3}\n\
                 cutoff ratio rc : {:.3}\n\
                 gamma γ         : {:.4}\n\n\
                 cutoff penalty  : {penalty:.4}\n\
                 efficiency η    : {eta:.4}  ({:.2} %)",
                s.compression_ratio,
                s.cutoff_ratio,
                s.gamma,
                eta * 100.0,
            ))
        }
        CycleType::Brayton => {
            let gamma = HeatCapacityRatio::new(s.gamma).map_err(|e| e.to_string())?;
            let b = Brayton::new(s.pressure_ratio, gamma).map_err(|e| e.to_string())?;
            let eta = b.efficiency();
            Ok(format!(
                "cycle           : {name}\n\
                 pressure r_p    : {:.3}\n\
                 gamma γ         : {:.4}\n\n\
                 efficiency η    : {eta:.4}  ({:.2} %)",
                s.pressure_ratio,
                s.gamma,
                eta * 100.0,
            ))
        }
        CycleType::Rankine => {
            let r = Rankine::new(s.h1, s.h2, s.h3, s.h4).map_err(|e| e.to_string())?;
            let eta = r.efficiency();
            let w_turbine = r.turbine_work();
            let w_pump = r.pump_work();
            let q_in = r.heat_in();
            let w_net = r.net_work();
            let bwr = r.back_work_ratio();
            Ok(format!(
                "cycle           : {name}\n\
                 enthalpies 1..4 : {:.1} / {:.1} / {:.1} / {:.1} kJ/kg\n\n\
                 turbine work    : {w_turbine:.2} kJ/kg\n\
                 pump work       : {w_pump:.2} kJ/kg\n\
                 net work        : {w_net:.2} kJ/kg\n\
                 heat in q_in    : {q_in:.2} kJ/kg\n\
                 back-work ratio : {bwr:.4}\n\
                 efficiency η    : {eta:.4}  ({:.2} %)",
                s.h1,
                s.h2,
                s.h3,
                s.h4,
                eta * 100.0,
            ))
        }
    }
}

/// Append an outward-facing box (centre `c`, half-extents `h`) to the
/// buffers.
fn push_box(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    c: Vector3<f64>,
    h: Vector3<f64>,
) {
    let base = nodes.len();
    let signs = [
        (-1.0, -1.0, -1.0),
        (1.0, -1.0, -1.0),
        (1.0, 1.0, -1.0),
        (-1.0, 1.0, -1.0),
        (-1.0, -1.0, 1.0),
        (1.0, -1.0, 1.0),
        (1.0, 1.0, 1.0),
        (-1.0, 1.0, 1.0),
    ];
    for (sx, sy, sz) in signs {
        nodes.push(c + Vector3::new(sx * h.x, sy * h.y, sz * h.z));
    }
    let faces = [
        [1, 2, 6, 5],
        [0, 4, 7, 3],
        [3, 7, 6, 2],
        [0, 1, 5, 4],
        [4, 5, 6, 7],
        [0, 3, 2, 1],
    ];
    for f in faces {
        tris.extend_from_slice(&[
            base + f[0],
            base + f[1],
            base + f[2],
            base + f[0],
            base + f[2],
            base + f[3],
        ]);
    }
}

/// Build a representative labeled block for the selected cycle as a
/// triangle [`Mesh`] — a main slab on a thin base. Representative geometry
/// (not to scale; the efficiency numbers are the `valenx-thermocycle`
/// result). `None` when the current inputs do not form a valid cycle.
fn cycle_block_mesh(s: &ThermoCycleWorkbenchState) -> Option<Mesh> {
    // Gate the solid on the same validation the readout uses.
    compute(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();
    // Main block.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.6),
        Vector3::new(0.5, 0.5, 0.5),
    );
    // Thin base.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.05),
        Vector3::new(0.7, 0.7, 0.05),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-thermocycle");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D cycle block and load it into the central viewport.
fn load_cycle_3d(app: &mut ValenxApp) {
    let Some(mesh) = cycle_block_mesh(&app.thermocycle) else {
        app.thermocycle.error =
            Some("cycle parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<cycle>/valenx-thermocycle"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_is_idle() {
        let s = ThermoCycleWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_carnot_reports_efficiency_and_cops() {
        let mut s = ThermoCycleWorkbenchState::default();
        run_thermocycle(&mut s);
        assert!(
            s.error.is_none(),
            "default Carnot should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("efficiency η"));
        assert!(s.result.contains("COP refrigerator"));
        assert!(s.result.contains("COP heat pump"));
        // 300 K / 600 K -> η = 0.5 exactly; COP_ref = 1.000, COP_hp = 2.000.
        assert!(s.result.contains("0.5000"));
        assert!(s.result.contains("1.000"));
        assert!(s.result.contains("2.000"));
    }

    #[test]
    fn carnot_ground_truth_one_minus_tc_over_th() {
        // Ground truth, hand-computed: the Carnot efficiency between
        // Tc = 300 K and Th = 600 K is 1 - 300/600 = 0.5 exactly.
        let c = Carnot::new(300.0, 600.0).unwrap();
        let hand = 1.0 - 300.0 / 600.0;
        assert!((c.efficiency() - hand).abs() < 1e-12);
        let pinned = 0.5_f64;
        assert!((c.efficiency() - pinned).abs() < 1e-12);
    }

    #[test]
    fn analyze_each_other_cycle_matches_its_classic_value() {
        // Otto r = 8 -> η ~ 0.5647; Brayton r_p = 10 -> η ~ 0.4821;
        // Diesel r = 18, r_c = 2 -> η ~ 0.6316 with cutoff penalty
        // (r_c^γ - 1) / (γ (r_c - 1)) ~ 1.1707; Rankine default ~ 0.4043.
        for (cycle, needle) in [
            (CycleType::Otto, "0.5647"),
            (CycleType::Brayton, "0.4821"),
            (CycleType::Rankine, "0.4043"),
        ] {
            let mut s = ThermoCycleWorkbenchState {
                cycle,
                ..Default::default()
            };
            run_thermocycle(&mut s);
            assert!(s.error.is_none(), "{cycle:?}: {:?}", s.error);
            assert!(
                s.result.contains(needle),
                "{cycle:?} missing {needle}: {}",
                s.result
            );
        }

        // Diesel at the classic r = 18, r_c = 2 -> η ~ 0.6316 (the shared
        // default compression ratio is 8, Otto's value, so the textbook
        // diesel case sets r = 18 explicitly), plus its cutoff penalty
        // factor (r_c^γ - 1) / (γ (r_c - 1)) ~ 1.1707.
        let mut diesel = ThermoCycleWorkbenchState {
            cycle: CycleType::Diesel,
            compression_ratio: 18.0,
            ..Default::default()
        };
        run_thermocycle(&mut diesel);
        assert!(diesel.error.is_none(), "diesel: {:?}", diesel.error);
        assert!(diesel.result.contains("0.6316"));
        assert!(diesel.result.contains("cutoff penalty"));
        assert!(diesel.result.contains("1.1707"));
    }

    #[test]
    fn analyze_rejects_unordered_carnot_temperatures() {
        // Tc >= Th leaves no temperature drop; the constructor rejects it
        // and the workbench surfaces the error.
        let mut s = ThermoCycleWorkbenchState {
            t_cold_k: 600.0,
            t_hot_k: 300.0,
            ..Default::default()
        };
        run_thermocycle(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn block_mesh_is_valid_for_default_and_none_for_invalid() {
        // A valid default cycle yields a non-empty, in-range triangle mesh.
        let s = ThermoCycleWorkbenchState::default();
        let mesh = cycle_block_mesh(&s).expect("default cycle yields a solid");
        assert!(mesh.nodes.len() > 8, "expected block + base");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }

        // An invalid configuration builds no solid.
        let bad = ThermoCycleWorkbenchState {
            t_cold_k: 600.0,
            t_hot_k: 300.0,
            ..Default::default()
        };
        assert!(cycle_block_mesh(&bad).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_thermocycle_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_thermocycle_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_thermocycle_workbench = true;
        run_thermocycle(&mut app.thermocycle);
        draw_workbench(&mut app);
    }
}
