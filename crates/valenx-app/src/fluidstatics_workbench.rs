//! The right-side **Fluid Statics Workbench** panel — native closed-form
//! hydrostatics over `valenx-fluid-statics`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_fluidstatics_workbench`,
//! toggled from the View menu. The form sets a fluid (water / sea water /
//! mercury) filling an open tank, a depth at which to read the gauge
//! pressure, and a vertical rectangular plate (a tank-wall gate) submerged
//! with its top edge at a chosen depth. "Analyze" reports the gauge
//! pressure `P = rho * g * h` at the probe depth and the resultant
//! hydrostatic force `F = rho * g * h_c * A` and centre-of-pressure depth
//! on the gate, and "Show 3-D" loads a representative water tank with the
//! submerged plate into the central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_fluid_statics::{gauge_pressure_standard, Fluid, RectangularPlate};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// The tank fluid choice exposed in the form.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum TankFluid {
    /// Fresh water (1000 kg/m³).
    Water,
    /// Sea water (1025 kg/m³).
    Seawater,
    /// Mercury (13 534 kg/m³).
    Mercury,
}

impl TankFluid {
    /// The corresponding [`Fluid`] density model.
    fn fluid(self) -> Fluid {
        match self {
            TankFluid::Water => Fluid::water(),
            TankFluid::Seawater => Fluid::seawater(),
            TankFluid::Mercury => Fluid::mercury(),
        }
    }

    /// A short display label for the readout.
    fn label(self) -> &'static str {
        match self {
            TankFluid::Water => "water",
            TankFluid::Seawater => "sea water",
            TankFluid::Mercury => "mercury",
        }
    }
}

/// Persistent form + result state for the Fluid Statics Workbench.
pub struct FluidStaticsWorkbenchState {
    /// The fluid filling the open tank.
    fluid: TankFluid,
    /// Total fluid depth in the tank (m).
    tank_depth_m: f64,
    /// Depth below the free surface at which to read the gauge pressure (m).
    probe_depth_m: f64,
    /// Width of the rectangular gate plate, along the surface line (m).
    plate_width_m: f64,
    /// Height of the rectangular gate plate, down the wall (m).
    plate_height_m: f64,
    /// Depth of the gate's top edge below the free surface (m).
    plate_top_depth_m: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D tank solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for FluidStaticsWorkbenchState {
    fn default() -> Self {
        // A 5 m deep open water tank, gauge read at 3 m, holding back a
        // vertical 2 m x 3 m gate whose top edge sits on the free surface
        // (centroid 1.5 m): P(3 m) ~ 29.4 kPa, gate force ~ 88.3 kN with
        // the centre of pressure at exactly 2/3 of the gate height (2.0 m).
        Self {
            fluid: TankFluid::Water,
            tank_depth_m: 5.0,
            probe_depth_m: 3.0,
            plate_width_m: 2.0,
            plate_height_m: 3.0,
            plate_top_depth_m: 0.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Fluid Statics Workbench right-side panel. A no-op when the
/// `show_fluidstatics_workbench` toggle is off.
pub fn draw_fluidstatics_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_fluidstatics_workbench {
        return;
    }

    egui::SidePanel::right("valenx_fluidstatics_workbench")
        .resizable(true)
        .default_width(360.0)
        .width_range(300.0..=560.0)
        .show(ctx, |ui| {
            if crate::workbench_ui::header(
                ui,
                "Fluid Statics",
                "native closed-form hydrostatics · valenx-fluid-statics",
            ) {
                app.show_fluidstatics_workbench = false;
            }

            let s = &mut app.fluidstatics;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Tank fluid").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.fluid, TankFluid::Water, "water");
                        ui.radio_value(&mut s.fluid, TankFluid::Seawater, "sea water");
                        ui.radio_value(&mut s.fluid, TankFluid::Mercury, "mercury");
                    });
                    ui.horizontal(|ui| {
                        ui.label("tank depth (m)");
                        ui.add(egui::DragValue::new(&mut s.tank_depth_m).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("gauge probe depth (m)");
                        ui.add(egui::DragValue::new(&mut s.probe_depth_m).speed(0.1));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Submerged gate plate").strong());
                    ui.horizontal(|ui| {
                        ui.label("width b (m)");
                        ui.add(egui::DragValue::new(&mut s.plate_width_m).speed(0.05));
                    });
                    ui.horizontal(|ui| {
                        ui.label("height H (m)");
                        ui.add(egui::DragValue::new(&mut s.plate_height_m).speed(0.05));
                    });
                    ui.horizontal(|ui| {
                        ui.label("top-edge depth (m)");
                        ui.add(egui::DragValue::new(&mut s.plate_top_depth_m).speed(0.1));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_fluidstatics(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative open fluid tank with the submerged vertical gate plate as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Hydrostatics").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        });

    // Serviced after the panel draws (the `&mut app.fluidstatics` borrow is
    // released here): build the tank's 3-D solid and load it.
    if app.fluidstatics.show_3d_request {
        app.fluidstatics.show_3d_request = false;
        load_tank_3d(app);
    }
}

/// Validate the form, evaluate the hydrostatics and format the readout.
fn run_fluidstatics(s: &mut FluidStaticsWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Evaluate the gauge pressure at the probe depth and the gate's resultant
/// force and centre of pressure, formatting the full readout and mapping
/// any domain error to a display string. Extracted so it is unit-testable.
fn compute(s: &FluidStaticsWorkbenchState) -> Result<String, String> {
    let fluid = s.fluid.fluid();
    let p_probe = gauge_pressure_standard(&fluid, s.probe_depth_m).map_err(|e| e.to_string())?;

    let plate =
        RectangularPlate::vertical(s.plate_width_m, s.plate_height_m).map_err(|e| e.to_string())?;
    let load = plate
        .load_from_top_edge(
            &fluid,
            valenx_fluid_statics::STANDARD_GRAVITY,
            s.plate_top_depth_m,
        )
        .map_err(|e| e.to_string())?;

    Ok(format!(
        "fluid           : {} ({:.0} kg/m³)\n\
         tank depth      : {:.2} m\n\
         probe depth     : {:.2} m\n\
         gauge pressure  : {:.1} Pa ({:.3} kPa)\n\n\
         gate b x H      : {:.3} m x {:.3} m\n\
         gate area A     : {:.4} m²\n\
         top-edge depth  : {:.2} m\n\
         centroid depth  : {:.4} m\n\
         resultant F     : {:.2} N ({:.3} kN)\n\
         centre of press : {:.4} m deep\n\
         CP below cent.  : {:.4} m",
        s.fluid.label(),
        fluid.density(),
        s.tank_depth_m,
        s.probe_depth_m,
        p_probe,
        p_probe / 1000.0,
        s.plate_width_m,
        s.plate_height_m,
        plate.area_m2(),
        s.plate_top_depth_m,
        load.centroid_depth_m,
        load.force_n,
        load.force_n / 1000.0,
        load.center_of_pressure_depth_m,
        load.cp_below_centroid_m(),
    ))
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

/// Build the fluid tank as a triangle [`Mesh`] — an open box of fluid with
/// a thin vertical rectangular gate plate submerged with its top edge at
/// the chosen depth. Representative geometry (not to scale; the
/// hydrostatic numbers are the `valenx-fluid-statics` result). `None` for
/// an invalid configuration.
fn tank_solid_mesh(s: &FluidStaticsWorkbenchState) -> Option<Mesh> {
    // Gate validity gates the whole solid, mirroring the analyze path.
    let plate = RectangularPlate::vertical(s.plate_width_m, s.plate_height_m).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Fluid tank body (z is depth; top of fluid at z = 0, going down to
    // z = -depth). Centred at half-depth.
    let depth = s.tank_depth_m.max(0.1);
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, -0.5 * depth),
        Vector3::new(1.2, 1.2, 0.5 * depth),
    );

    // Submerged vertical gate plate: thin in x, width b in y, height H in
    // z, with its top edge at the chosen depth (z = -top_depth) and the
    // plate hanging down to z = -(top_depth + H).
    let top = -s.plate_top_depth_m;
    let half_h = 0.5 * plate.height_m();
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, top - half_h),
        Vector3::new(0.04, 0.5 * plate.width_m(), half_h),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-fluid-statics");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D tank solid and load it into the central viewport.
fn load_tank_3d(app: &mut ValenxApp) {
    let Some(mesh) = tank_solid_mesh(&app.fluidstatics) else {
        app.fluidstatics.error =
            Some("tank / gate parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<tank>/valenx-fluid-statics"),
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
        let s = FluidStaticsWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_pressure_force_and_cp() {
        let mut s = FluidStaticsWorkbenchState::default();
        run_fluidstatics(&mut s);
        assert!(
            s.error.is_none(),
            "default tank should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("gauge pressure"));
        assert!(s.result.contains("resultant F"));
        assert!(s.result.contains("centre of press"));
        // Probe at 3 m in water: P = 1000*9.80665*3 = 29419.95 Pa = 29.420 kPa.
        assert!(s.result.contains("29.420 kPa"));
        // Surface-piercing 2x3 vertical gate: F = 1000*9.80665*1.5*6 N,
        // and the centre of pressure sits at exactly 2/3 of the 3 m height.
        assert!(s.result.contains("88259.85 N"));
        assert!(s.result.contains("2.0000 m deep"));
    }

    #[test]
    fn analyze_rejects_zero_plate_width() {
        let mut s = FluidStaticsWorkbenchState {
            plate_width_m: 0.0,
            ..Default::default()
        };
        run_fluidstatics(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn ground_truth_pressure_and_force_match_hand_calc() {
        // Ground truth: gauge pressure P = rho*g*h and resultant force on a
        // surface-piercing vertical gate F = rho*g*h_c*A, hand-computed.
        let rho = 1000.0_f64;
        let g = valenx_fluid_statics::STANDARD_GRAVITY;
        let s = FluidStaticsWorkbenchState::default();

        let p_expected = rho * g * s.probe_depth_m;
        let p_actual = gauge_pressure_standard(&s.fluid.fluid(), s.probe_depth_m).unwrap();
        assert!((p_actual - p_expected).abs() < 1e-6, "p {p_actual}");

        let plate = RectangularPlate::vertical(s.plate_width_m, s.plate_height_m).unwrap();
        let area = s.plate_width_m * s.plate_height_m;
        let h_c = 0.5 * s.plate_height_m; // top edge on the surface
        let f_expected = rho * g * h_c * area;
        let load = plate.load_from_top_edge(&s.fluid.fluid(), g, 0.0).unwrap();
        assert!(
            (load.force_n - f_expected).abs() < 1e-6,
            "F {}",
            load.force_n
        );
        // Surface-piercing vertical rectangle -> CP at 2/3 of the height.
        assert!(
            (load.center_of_pressure_depth_m - 2.0 * s.plate_height_m / 3.0).abs() < 1e-9,
            "cp {}",
            load.center_of_pressure_depth_m
        );
    }

    #[test]
    fn mercury_is_denser_so_pressure_is_higher() {
        let mut water = FluidStaticsWorkbenchState::default();
        let mut mercury = FluidStaticsWorkbenchState {
            fluid: TankFluid::Mercury,
            ..Default::default()
        };
        run_fluidstatics(&mut water);
        run_fluidstatics(&mut mercury);
        let p_water = gauge_pressure_standard(&Fluid::water(), 3.0).unwrap();
        let p_mercury = gauge_pressure_standard(&Fluid::mercury(), 3.0).unwrap();
        assert!(p_mercury > p_water);
        assert!(water.error.is_none() && mercury.error.is_none());
    }

    #[test]
    fn tank_mesh_for_default_is_nonempty_and_in_range() {
        let s = FluidStaticsWorkbenchState::default();
        let mesh = tank_solid_mesh(&s).expect("default tank yields a solid");
        assert!(mesh.nodes.len() > 8, "expected tank + gate plate");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn tank_mesh_none_for_invalid() {
        let s = FluidStaticsWorkbenchState {
            plate_height_m: 0.0,
            ..Default::default()
        };
        assert!(tank_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_fluidstatics_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_fluidstatics_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_fluidstatics_workbench = true;
        run_fluidstatics(&mut app.fluidstatics);
        draw_workbench(&mut app);
    }
}
