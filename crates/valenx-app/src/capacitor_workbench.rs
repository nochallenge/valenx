//! The right-side **Capacitor Workbench** panel — native ideal
//! parallel-plate electrostatics over `valenx-capacitor`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_capacitor_workbench`,
//! toggled from the View menu. The form sets a parallel-plate geometry
//! (overlapping plate area, gap, dielectric relative permittivity) plus a
//! terminal voltage, a drive frequency and a series charging resistance;
//! "Analyze" computes the ideal capacitance `C = eps_r eps0 A / d`, the
//! stored energy `E = 1/2 C V^2`, the stored charge `Q = C V`, the
//! capacitive reactance `X_C = 1 / (2 pi f C)` and the RC charging time
//! constant `tau = R C`, and "Show 3-D capacitor" loads a representative
//! parallel-plate solid (two plates with a dielectric slab between) into the
//! central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_capacitor::parallel_plate;
use valenx_capacitor::reactance;
use valenx_capacitor::transient;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Capacitor Workbench.
pub struct CapacitorWorkbenchState {
    /// Overlapping plate area `A` (m^2).
    area_m2: f64,
    /// Plate separation / dielectric thickness `d` (m).
    gap_m: f64,
    /// Relative permittivity of the dielectric `eps_r` (dimensionless, >= 1).
    eps_r: f64,
    /// Terminal voltage `V` (volts).
    voltage_v: f64,
    /// Drive frequency `f` for the reactance (Hz).
    frequency_hz: f64,
    /// Series charging resistance `R` for the RC time constant (ohms).
    series_resistance_ohm: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D capacitor solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for CapacitorWorkbenchState {
    fn default() -> Self {
        // A 100 cm^2 = 0.01 m^2 parallel-plate capacitor with a 0.1 mm
        // dielectric film of eps_r = 3.5, charged to 50 V and driven at
        // 1 kHz: C ~ 3.10 nF, E ~ 3.87 uJ, Q ~ 155 nC, X_C ~ 51.4 kohm.
        // Charged through a 10 kohm series resistor: tau = R C ~ 31.0 us.
        Self {
            area_m2: 0.01,
            gap_m: 0.0001,
            eps_r: 3.5,
            voltage_v: 50.0,
            frequency_hz: 1000.0,
            series_resistance_ohm: 10_000.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Capacitor Workbench right-side panel. A no-op when the
/// `show_capacitor_workbench` toggle is off.
pub fn draw_capacitor_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_capacitor_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_capacitor_workbench",
        "Capacitor",
        |app, ui| {
            ui.label(
                egui::RichText::new("native ideal parallel-plate electrostatics · valenx-capacitor")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.capacitor;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Geometry").strong());
                    ui.horizontal(|ui| {
                        ui.label("plate area A (m²)");
                        ui.add(egui::DragValue::new(&mut s.area_m2).speed(0.001));
                    });
                    ui.horizontal(|ui| {
                        ui.label("gap d (m)");
                        ui.add(egui::DragValue::new(&mut s.gap_m).speed(0.00001));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Dielectric").strong());
                    ui.horizontal(|ui| {
                        ui.label("relative permittivity εr");
                        ui.add(egui::DragValue::new(&mut s.eps_r).speed(0.1));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Excitation").strong());
                    ui.horizontal(|ui| {
                        ui.label("voltage V (V)");
                        ui.add(egui::DragValue::new(&mut s.voltage_v).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("frequency f (Hz)");
                        ui.add(egui::DragValue::new(&mut s.frequency_hz).speed(10.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("series R (Ω)");
                        ui.add(egui::DragValue::new(&mut s.series_resistance_ohm).speed(100.0));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_capacitor(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D capacitor").strong())
                        .on_hover_text(
                            "Build a representative parallel-plate capacitor (two conductive plates with a dielectric slab between) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Electrostatics").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_capacitor_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.capacitor` borrow is
    // released here): build the capacitor's 3-D solid and load it.
    if app.capacitor.show_3d_request {
        app.capacitor.show_3d_request = false;
        load_capacitor_3d(app);
    }
}

/// Validate the form, evaluate the capacitor and format the readout.
fn run_capacitor(s: &mut CapacitorWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Evaluate the parallel-plate capacitor and format the full readout,
/// mapping any domain error to a display string. Extracted so it is
/// unit-testable.
fn compute(s: &CapacitorWorkbenchState) -> Result<String, String> {
    let c = parallel_plate::capacitance(s.eps_r, s.area_m2, s.gap_m).map_err(|e| e.to_string())?;
    let energy = parallel_plate::stored_energy(c, s.voltage_v).map_err(|e| e.to_string())?;
    let q = parallel_plate::charge(c, s.voltage_v).map_err(|e| e.to_string())?;
    let xc = reactance::reactance(s.frequency_hz, c).map_err(|e| e.to_string())?;
    let tau = transient::time_constant(s.series_resistance_ohm, c).map_err(|e| e.to_string())?;

    Ok(format!(
        "plate area A    : {:.4} m²\n\
         gap d           : {:.6} m\n\
         dielectric εr   : {:.3}\n\
         voltage / freq  : {:.1} V / {:.1} Hz\n\
         series R        : {:.1} Ω\n\n\
         capacitance C   : {:.4} nF\n\
         stored energy E : {:.4} µJ\n\
         stored charge Q : {:.4} nC\n\
         reactance X_C   : {:.2} Ω\n\
         RC constant τ   : {:.4} µs",
        s.area_m2,
        s.gap_m,
        s.eps_r,
        s.voltage_v,
        s.frequency_hz,
        s.series_resistance_ohm,
        c * 1.0e9,
        energy * 1.0e6,
        q * 1.0e9,
        xc,
        tau * 1.0e6,
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

/// Build the parallel-plate capacitor as a triangle [`Mesh`] — two flat
/// conductive plates (thin boxes) separated along the field (`x`) direction
/// with a dielectric slab box filling the gap between them. Representative
/// geometry (not to scale; the electrostatic numbers are the
/// `valenx-capacitor` result). `None` for an invalid configuration.
fn cap_solid_mesh(s: &CapacitorWorkbenchState) -> Option<Mesh> {
    // Gate the geometry on a valid, buildable capacitor object.
    valenx_capacitor::ParallelPlate::new(s.eps_r, s.area_m2, s.gap_m).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Left conductive plate (-x face).
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(-0.11, 0.0, 0.6),
        Vector3::new(0.02, 0.7, 0.5),
    );
    // Right conductive plate (+x face).
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.11, 0.0, 0.6),
        Vector3::new(0.02, 0.7, 0.5),
    );
    // Dielectric slab filling the gap between the plates.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.6),
        Vector3::new(0.085, 0.66, 0.46),
    );
    // Base.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.05),
        Vector3::new(0.2, 0.7, 0.05),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-capacitor");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D capacitor solid and load it into the central viewport.
fn load_capacitor_3d(app: &mut ValenxApp) {
    let Some(mesh) = cap_solid_mesh(&app.capacitor) else {
        app.capacitor.error =
            Some("capacitor parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<capacitor>/valenx-capacitor"),
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
        let s = CapacitorWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_capacitance_energy_and_reactance() {
        let mut s = CapacitorWorkbenchState::default();
        run_capacitor(&mut s);
        assert!(
            s.error.is_none(),
            "default capacitor should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("capacitance C"));
        assert!(s.result.contains("stored energy E"));
        assert!(s.result.contains("reactance X_C"));
        assert!(s.result.contains("RC constant τ"));
        // 0.01 m^2, 0.1 mm, eps_r 3.5 -> C ~ 3.099 nF; 50 V -> E ~ 3.874 uJ;
        // Q ~ 154.95 nC; at 1 kHz X_C ~ 51357.44 ohm; with R = 10 kohm
        // tau = R C ~ 30.9897 us.
        assert!(s.result.contains("3.0990"));
        assert!(s.result.contains("3.8737"));
        assert!(s.result.contains("154.9483"));
        assert!(s.result.contains("51357.44"));
        assert!(s.result.contains("30.9897"));
    }

    #[test]
    fn analyze_rejects_zero_gap() {
        let mut s = CapacitorWorkbenchState {
            gap_m: 0.0,
            ..Default::default()
        };
        run_capacitor(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn compute_matches_closed_form_ground_truth() {
        // Ground truth: C = eps_r * eps0 * A / d, E = 1/2 C V^2,
        // X_C = 1 / (2 pi f C) — recomputed here independently.
        let s = CapacitorWorkbenchState::default();
        let c = parallel_plate::VACUUM_PERMITTIVITY * s.eps_r * s.area_m2 / s.gap_m;
        let energy = 0.5 * c * s.voltage_v * s.voltage_v;
        let xc = 1.0 / (2.0 * std::f64::consts::PI * s.frequency_hz * c);

        let c_api = parallel_plate::capacitance(s.eps_r, s.area_m2, s.gap_m).unwrap();
        let e_api = parallel_plate::stored_energy(c_api, s.voltage_v).unwrap();
        let xc_api = reactance::reactance(s.frequency_hz, c_api).unwrap();

        assert!((c_api - c).abs() < 1e-18);
        assert!((e_api - energy).abs() < 1e-18);
        assert!((xc_api - xc).abs() < 1e-6 * xc);
    }

    #[test]
    fn compute_reports_rc_time_constant_ground_truth() {
        // Ground truth: tau = R * C with C = eps_r * eps0 * A / d.
        // Default: R = 10000, C = 8.8541878128e-12 * 3.5 * 0.01 / 0.0001
        //        = 3.098965734e-9 F  ->  tau = 3.098965734e-5 s = 30.9897 us.
        let s = CapacitorWorkbenchState::default();
        let c = parallel_plate::VACUUM_PERMITTIVITY * s.eps_r * s.area_m2 / s.gap_m;
        let tau = s.series_resistance_ohm * c;

        let c_api = parallel_plate::capacitance(s.eps_r, s.area_m2, s.gap_m).unwrap();
        let tau_api = transient::time_constant(s.series_resistance_ohm, c_api).unwrap();

        // Hand value: 10_000 * 3.0989657344800006e-9 = 3.09896573448e-5 s.
        let tau_expected = 3.098_965_734_48e-5;
        assert!(
            (tau_api - tau).abs() < 1e-18,
            "tau_api {tau_api} vs recomputed {tau}"
        );
        assert!(
            (tau_api - tau_expected).abs() < 1e-12,
            "tau_api {tau_api} vs hand value {tau_expected}"
        );

        // The formatted readout shows it in microseconds to 4 dp.
        let r = compute(&s).unwrap();
        assert!(r.contains("RC constant τ   : 30.9897 µs"), "readout: {r}");
    }

    #[test]
    fn cap_mesh_for_default_is_nonempty_and_in_range() {
        let s = CapacitorWorkbenchState::default();
        let mesh = cap_solid_mesh(&s).expect("default capacitor yields a solid");
        assert!(
            mesh.nodes.len() > 8,
            "expected two plates + dielectric + base"
        );
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn cap_mesh_none_for_invalid() {
        let s = CapacitorWorkbenchState {
            eps_r: 0.5,
            ..Default::default()
        };
        assert!(cap_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_capacitor_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_capacitor_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_capacitor_workbench = true;
        run_capacitor(&mut app.capacitor);
        draw_workbench(&mut app);
    }
}
