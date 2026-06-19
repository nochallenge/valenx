//! The right-side **Hemodynamics Workbench** panel — native closed-form
//! cardiovascular analysis over `valenx-hemodynamics`.
//!
//! Mirrors the Heat Transfer / Heat Exchanger workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_hemodynamics_workbench`,
//! toggled from the View menu. The form sets a heart (rate, stroke
//! volume), a single straight vessel segment (radius, length, blood
//! viscosity, pressure drop) and a 2-element Windkessel diastolic decay
//! (start pressure, peripheral resistance, arterial compliance, elapsed
//! time). "Analyze" reports the cardiac output, the Hagen-Poiseuille
//! segment flow / resistance / wall shear stress, and the Windkessel time
//! constant and decayed pressure, and "Show 3-D vessel" loads a
//! representative blood-vessel tube into the central viewport.
//!
//! These are research/educational closed-form relations — not a clinical
//! tool — exactly as `valenx-hemodynamics` documents.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_hemodynamics::circulation::cardiac_output;
use valenx_hemodynamics::flow::{flow_from_resistance, vascular_resistance, wall_shear_stress};
use valenx_hemodynamics::windkessel::{time_constant, windkessel_pressure};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Representative vessel calibre, selecting only the radius of the 3-D
/// tube preview (the numeric flow result is driven by the form's own
/// `radius_mm`, never by this).
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum VesselScale {
    /// A large conduit artery (e.g. aorta-scale).
    Artery,
    /// A small resistance arteriole.
    Arteriole,
}

impl VesselScale {
    /// Half-radius (in scene units) of the representative 3-D tube.
    fn tube_radius(self) -> f64 {
        match self {
            VesselScale::Artery => 0.18,
            VesselScale::Arteriole => 0.07,
        }
    }

    /// Short human label for the readout.
    fn label(self) -> &'static str {
        match self {
            VesselScale::Artery => "artery",
            VesselScale::Arteriole => "arteriole",
        }
    }
}

/// Persistent form + result state for the Hemodynamics Workbench.
pub struct HemodynamicsWorkbenchState {
    /// Heart rate `HR` (beats/min).
    heart_rate_bpm: f64,
    /// Stroke volume `SV` (mL/beat).
    stroke_volume_ml: f64,
    /// Vessel-segment radius `r` (mm).
    radius_mm: f64,
    /// Vessel-segment length `L` (mm).
    length_mm: f64,
    /// Blood dynamic viscosity `mu` (mPa·s).
    viscosity_mpas: f64,
    /// Pressure drop `dP` along the segment (Pa).
    pressure_drop_pa: f64,
    /// Windkessel start (end-systolic) pressure `P0` (mmHg).
    windkessel_p0_mmhg: f64,
    /// Windkessel peripheral resistance `R` (Pa·s/m^3).
    windkessel_r: f64,
    /// Windkessel arterial compliance `C` (m^3/Pa).
    windkessel_c: f64,
    /// Elapsed diastolic time `t` for the decay (s).
    decay_time_s: f64,
    /// Representative tube calibre for the 3-D preview.
    scale: VesselScale,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D vessel solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for HemodynamicsWorkbenchState {
    fn default() -> Self {
        // A resting adult: HR = 70 bpm, SV = 70 mL -> CO = 4900 mL/min
        // (4.9 L/min). A 2 mm-radius, 50 mm arteriolar segment of blood
        // (mu = 3.5 mPa·s) under a 1000 Pa (~7.5 mmHg) drop. A 2-element
        // Windkessel with tau = R*C = 1.2 s relaxing 80 mmHg over 1.2 s.
        Self {
            heart_rate_bpm: 70.0,
            stroke_volume_ml: 70.0,
            radius_mm: 2.0,
            length_mm: 50.0,
            viscosity_mpas: 3.5,
            pressure_drop_pa: 1000.0,
            windkessel_p0_mmhg: 80.0,
            windkessel_r: 1.2e8,
            windkessel_c: 1.0e-8,
            decay_time_s: 1.2,
            scale: VesselScale::Arteriole,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Hemodynamics Workbench right-side panel. A no-op when the
/// `show_hemodynamics_workbench` toggle is off.
pub fn draw_hemodynamics_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_hemodynamics_workbench {
        return;
    }

    egui::SidePanel::right("valenx_hemodynamics_workbench")
        .resizable(true)
        .default_width(360.0)
        .width_range(300.0..=560.0)
        .show(ctx, |ui| {
            if crate::workbench_ui::header(
                ui,
                "Hemodynamics",
                "native closed-form cardiovascular analysis · valenx-hemodynamics",
            ) {
                app.show_hemodynamics_workbench = false;
            }

            let s = &mut app.hemodynamics;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Heart").strong());
                    ui.horizontal(|ui| {
                        ui.label("heart rate (bpm)");
                        ui.add(egui::DragValue::new(&mut s.heart_rate_bpm).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("stroke volume (mL)");
                        ui.add(egui::DragValue::new(&mut s.stroke_volume_ml).speed(1.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Vessel segment").strong());
                    ui.horizontal(|ui| {
                        ui.label("radius (mm)");
                        ui.add(egui::DragValue::new(&mut s.radius_mm).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("length (mm)");
                        ui.add(egui::DragValue::new(&mut s.length_mm).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("viscosity (mPa·s)");
                        ui.add(egui::DragValue::new(&mut s.viscosity_mpas).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("pressure drop (Pa)");
                        ui.add(egui::DragValue::new(&mut s.pressure_drop_pa).speed(10.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("3-D calibre");
                        ui.radio_value(&mut s.scale, VesselScale::Artery, "artery");
                        ui.radio_value(&mut s.scale, VesselScale::Arteriole, "arteriole");
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Windkessel diastolic decay").strong());
                    ui.horizontal(|ui| {
                        ui.label("start P0 (mmHg)");
                        ui.add(egui::DragValue::new(&mut s.windkessel_p0_mmhg).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("resistance R (Pa·s/m³)");
                        ui.add(egui::DragValue::new(&mut s.windkessel_r).speed(1.0e6));
                    });
                    ui.horizontal(|ui| {
                        ui.label("compliance C (m³/Pa)");
                        ui.add(egui::DragValue::new(&mut s.windkessel_c).speed(1.0e-10));
                    });
                    ui.horizontal(|ui| {
                        ui.label("elapsed t (s)");
                        ui.add(egui::DragValue::new(&mut s.decay_time_s).speed(0.05));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_hemodynamics(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D vessel").strong())
                        .on_hover_text(
                            "Build a representative blood-vessel tube (a straight rigid cylinder, the calibre of the selected vessel scale) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Hemodynamics").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        });

    // Serviced after the panel draws (the `&mut app.hemodynamics` borrow is
    // released here): build the vessel's 3-D solid and load it.
    if app.hemodynamics.show_3d_request {
        app.hemodynamics.show_3d_request = false;
        load_vessel_3d(app);
    }
}

/// Validate the form, evaluate the relations and format the readout.
fn run_hemodynamics(s: &mut HemodynamicsWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Evaluate all the cardiovascular relations and format the full readout,
/// mapping any domain error to a display string. Extracted so it is
/// unit-testable.
///
/// The heart relation uses clinical units (bpm · mL → mL/min, reported as
/// L/min); the single-vessel and Windkessel relations use SI throughout
/// (the segment radius / length are entered in mm and converted to metres
/// here), exactly as `valenx-hemodynamics` documents.
fn compute(s: &HemodynamicsWorkbenchState) -> Result<String, String> {
    // Cardiac output, clinical units: bpm * mL = mL/min.
    let co_ml_min =
        cardiac_output(s.heart_rate_bpm, s.stroke_volume_ml).map_err(|e| e.to_string())?;
    let co_l_min = co_ml_min / 1000.0;

    // Single straight vessel segment, SI units.
    let radius_m = s.radius_mm / 1000.0;
    let length_m = s.length_mm / 1000.0;
    let viscosity = s.viscosity_mpas / 1000.0; // mPa·s -> Pa·s
    let resistance =
        vascular_resistance(viscosity, length_m, radius_m).map_err(|e| e.to_string())?;
    let q = flow_from_resistance(s.pressure_drop_pa, resistance).map_err(|e| e.to_string())?;
    // Wall shear stress uses the |flow| magnitude (a reversed flow shears
    // the wall just as hard); the model itself requires a non-negative flow.
    let tau = wall_shear_stress(viscosity, q.abs(), radius_m).map_err(|e| e.to_string())?;
    let q_ml_min = q * 60.0e6; // m^3/s -> mL/min

    // 2-element Windkessel diastolic decay, SI time.
    let tau_wk = time_constant(s.windkessel_r, s.windkessel_c).map_err(|e| e.to_string())?;
    let p_decay = windkessel_pressure(
        s.windkessel_p0_mmhg,
        s.decay_time_s,
        s.windkessel_r,
        s.windkessel_c,
    )
    .map_err(|e| e.to_string())?;

    Ok(format!(
        "heart rate / SV : {:.0} bpm / {:.0} mL\n\
         cardiac output  : {:.2} L/min\n\n\
         vessel ({})\n\
         radius / length : {:.2} mm / {:.1} mm\n\
         viscosity       : {:.3} mPa·s\n\
         pressure drop   : {:.1} Pa\n\
         resistance R    : {:.3e} Pa·s/m³\n\
         flow Q          : {:.3} mL/min\n\
         wall shear τ_w  : {:.3} Pa\n\n\
         Windkessel\n\
         time constant   : {:.3} s\n\
         P({:.2} s)       : {:.2} mmHg",
        s.heart_rate_bpm,
        s.stroke_volume_ml,
        co_l_min,
        s.scale.label(),
        s.radius_mm,
        s.length_mm,
        s.viscosity_mpas,
        s.pressure_drop_pa,
        resistance,
        q_ml_min,
        tau,
        tau_wk,
        s.decay_time_s,
        p_decay,
    ))
}

/// Append a (double-sided) cylinder whose axis runs along `+z`, spanning
/// `base.z ..= base.z + height` with circle centre `(base.x, base.y)`.
fn push_cyl_z(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    base: Vector3<f64>,
    height: f64,
    r: f64,
    seg: usize,
) {
    use std::f64::consts::TAU;
    let (z0, z1) = (base.z, base.z + height);
    let bot = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(base.x + r * a.cos(), base.y + r * a.sin(), z0));
    }
    let top = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(base.x + r * a.cos(), base.y + r * a.sin(), z1));
    }
    for j in 0..seg {
        let jn = (j + 1) % seg;
        tris.extend_from_slice(&[
            bot + j,
            top + j,
            top + jn,
            bot + j,
            top + jn,
            bot + jn,
            bot + j,
            top + jn,
            top + j,
            bot + j,
            bot + jn,
            top + jn,
        ]);
    }
}

/// Build the vessel as a triangle [`Mesh`] — a representative straight
/// rigid blood-vessel tube (a single `+z` cylinder, the calibre of the
/// selected vessel scale). Representative geometry (not to scale; the flow
/// numbers are the `valenx-hemodynamics` result). `None` for an invalid
/// configuration.
fn vessel_solid_mesh(s: &HemodynamicsWorkbenchState) -> Option<Mesh> {
    // Gate the 3-D build on a consistent vessel configuration — the same
    // domain checks the numeric readout runs.
    let radius_m = s.radius_mm / 1000.0;
    let length_m = s.length_mm / 1000.0;
    let viscosity = s.viscosity_mpas / 1000.0;
    vascular_resistance(viscosity, length_m, radius_m).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // One straight vertical tube, the selected calibre.
    push_cyl_z(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.1),
        1.0,
        s.scale.tube_radius(),
        32,
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-hemodynamics");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D vessel solid and load it into the central viewport.
fn load_vessel_3d(app: &mut ValenxApp) {
    let Some(mesh) = vessel_solid_mesh(&app.hemodynamics) else {
        app.hemodynamics.error =
            Some("vessel parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<vessel>/valenx-hemodynamics"),
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
        let s = HemodynamicsWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_all_sections() {
        let mut s = HemodynamicsWorkbenchState::default();
        run_hemodynamics(&mut s);
        assert!(
            s.error.is_none(),
            "default config should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("cardiac output"));
        assert!(s.result.contains("flow Q"));
        assert!(s.result.contains("wall shear"));
        assert!(s.result.contains("time constant"));
        // CO = 70 bpm * 70 mL = 4900 mL/min = 4.90 L/min ({:.2} L/min).
        assert!(s.result.contains("4.90 L/min"));
        // tau = R*C = 1.2e8 * 1.0e-8 = 1.2 s ({:.3} s).
        assert!(s.result.contains("1.200 s"));
    }

    #[test]
    fn analyze_rejects_zero_radius() {
        let mut s = HemodynamicsWorkbenchState {
            radius_mm: 0.0,
            ..Default::default()
        };
        run_hemodynamics(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn cardiac_output_is_hr_times_sv_ground_truth() {
        // GROUND TRUTH (hand-computed): CO = HR * SV = 70 /min * 70 mL
        // = 4900 mL/min = 4.9 L/min, the textbook resting cardiac output.
        let co_ml_min = cardiac_output(70.0, 70.0).expect("valid");
        assert!((co_ml_min - 4900.0).abs() < 1e-9);
        assert!((co_ml_min / 1000.0 - 4.9).abs() < 1e-12);
    }

    #[test]
    fn segment_flow_obeys_ohm_analogue() {
        // VALIDATE: with R the vessel's hydraulic resistance, Q * R = dP.
        let s = HemodynamicsWorkbenchState::default();
        let radius_m = s.radius_mm / 1000.0;
        let length_m = s.length_mm / 1000.0;
        let viscosity = s.viscosity_mpas / 1000.0;
        let r = vascular_resistance(viscosity, length_m, radius_m).expect("valid");
        let q = flow_from_resistance(s.pressure_drop_pa, r).expect("valid");
        assert!((q * r - s.pressure_drop_pa).abs() < 1e-6 * s.pressure_drop_pa);
        assert!(q > 0.0);
    }

    #[test]
    fn vessel_mesh_for_default_is_nonempty_and_in_range() {
        let s = HemodynamicsWorkbenchState::default();
        let mesh = vessel_solid_mesh(&s).expect("default vessel yields a solid");
        assert!(mesh.nodes.len() > 8, "expected a tube of ring nodes");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn vessel_mesh_none_for_invalid() {
        let s = HemodynamicsWorkbenchState {
            radius_mm: 0.0,
            ..Default::default()
        };
        assert!(vessel_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_hemodynamics_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_hemodynamics_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_hemodynamics_workbench = true;
        run_hemodynamics(&mut app.hemodynamics);
        draw_workbench(&mut app);
    }
}
