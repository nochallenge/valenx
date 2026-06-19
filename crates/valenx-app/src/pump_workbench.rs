//! The right-side **Pump Workbench** panel — native centrifugal-pump
//! duty-point analysis over `valenx-pump`.
//!
//! Mirrors the Solar PV / Wind Turbine / Rail workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_pump_workbench`,
//! toggled from the View menu. The form drives a quadratic pump curve
//! `H = H0 − a·Q²` against a system curve `H = Hs + K·Q²`; "Analyze"
//! reports the operating point (where the two heads balance), the
//! hydraulic and shaft power there, and the NPSH cavitation margin, and
//! "Show 3-D pump" loads a volute-casing-on-a-baseplate solid into the
//! central viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_pump::npsh::{available_npsh_m, is_cavitation_free, npsh_margin_m, SuctionConditions};
use valenx_pump::operating::{operating_point, PumpCurve};
use valenx_pump::power::{hydraulic_power_w, shaft_power_w};
use valenx_pump::system::SystemCurve;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Pump Workbench.
pub struct PumpWorkbenchState {
    /// Pump shut-off head `H0` (m of fluid, head at zero flow).
    shutoff_head_m: f64,
    /// Pump droop coefficient `a > 0` (m per (m³/s)²).
    droop_a: f64,
    /// System static head `Hs` (m, the lift the pump must overcome at
    /// zero flow).
    static_head_m: f64,
    /// System resistance `K >= 0` (m per (m³/s)², the friction term).
    resistance_k: f64,
    /// Fluid density `rho` (kg/m³).
    density_kg_m3: f64,
    /// Pump efficiency `eta` in (0, 1].
    efficiency: f64,
    /// Surface (atmospheric) pressure over the suction source (Pa).
    atmospheric_pa: f64,
    /// Liquid vapour pressure at the pumping temperature (Pa).
    vapor_pressure_pa: f64,
    /// Static suction head `+` flooded / `-` lift (m).
    static_suction_head_m: f64,
    /// Suction-line friction loss (m).
    suction_loss_m: f64,
    /// Pump-required NPSH `NPSHr` (m).
    required_npsh_m: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D pump solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for PumpWorkbenchState {
    fn default() -> Self {
        // A representative end-suction centrifugal pump on cold water at
        // sea level: 50 m shut-off, system 10 m static + K = 4000, which
        // balance near 0.089 m³/s at 42 m.
        Self {
            shutoff_head_m: 50.0,
            droop_a: 1000.0,
            static_head_m: 10.0,
            resistance_k: 4000.0,
            density_kg_m3: 1000.0,
            efficiency: 0.75,
            atmospheric_pa: 101_325.0,
            vapor_pressure_pa: 2_340.0,
            static_suction_head_m: 2.0,
            suction_loss_m: 0.5,
            required_npsh_m: 3.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Pump Workbench right-side panel. A no-op when the
/// `show_pump_workbench` toggle is off.
pub fn draw_pump_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_pump_workbench {
        return;
    }

    egui::SidePanel::right("valenx_pump_workbench")
        .resizable(true)
        .default_width(360.0)
        .width_range(300.0..=560.0)
        .show(ctx, |ui| {
            if crate::workbench_ui::header(
                ui,
                "Pump",
                "native centrifugal-pump duty point + NPSH · valenx-pump",
            ) {
                app.show_pump_workbench = false;
            }

            let s = &mut app.pump;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Pump curve  H = H0 − a·Q²").strong());
                    ui.horizontal(|ui| {
                        ui.label("shut-off head H0 (m)");
                        ui.add(egui::DragValue::new(&mut s.shutoff_head_m).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("droop a (m·s²/m⁶)");
                        ui.add(egui::DragValue::new(&mut s.droop_a).speed(10.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("System curve  H = Hs + K·Q²").strong());
                    ui.horizontal(|ui| {
                        ui.label("static head Hs (m)");
                        ui.add(egui::DragValue::new(&mut s.static_head_m).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("resistance K (m·s²/m⁶)");
                        ui.add(egui::DragValue::new(&mut s.resistance_k).speed(10.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Fluid + drive").strong());
                    ui.horizontal(|ui| {
                        ui.label("density ρ (kg/m³)");
                        ui.add(egui::DragValue::new(&mut s.density_kg_m3).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("efficiency η");
                        ui.add(egui::DragValue::new(&mut s.efficiency).speed(0.01));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Suction / NPSH").strong());
                    ui.horizontal(|ui| {
                        ui.label("atmospheric (Pa)");
                        ui.add(egui::DragValue::new(&mut s.atmospheric_pa).speed(100.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("vapour pressure (Pa)");
                        ui.add(egui::DragValue::new(&mut s.vapor_pressure_pa).speed(50.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("static suction head (m)");
                        ui.add(egui::DragValue::new(&mut s.static_suction_head_m).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("suction loss (m)");
                        ui.add(egui::DragValue::new(&mut s.suction_loss_m).speed(0.05));
                    });
                    ui.horizontal(|ui| {
                        ui.label("required NPSHr (m)");
                        ui.add(egui::DragValue::new(&mut s.required_npsh_m).speed(0.1));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_pump(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D pump").strong())
                        .on_hover_text(
                            "Build a volute casing with suction eye, vertical discharge and motor on a baseplate as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Duty point").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        });

    // Serviced after the panel draws (the `&mut app.pump` borrow is
    // released here): build the pump's 3-D solid and load it.
    if app.pump.show_3d_request {
        app.pump.show_3d_request = false;
        load_pump_3d(app);
    }
}

/// Validate the form, compute the duty point and format the readout.
fn run_pump(s: &mut PumpWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Compute the full duty-point readout, mapping any domain error to a
/// display string. Extracted so it is unit-testable.
fn compute(s: &PumpWorkbenchState) -> Result<String, String> {
    let pump = PumpCurve::new(s.shutoff_head_m, s.droop_a).map_err(|e| e.to_string())?;
    let system = SystemCurve::new(s.static_head_m, s.resistance_k).map_err(|e| e.to_string())?;
    let op = operating_point(&pump, &system).map_err(|e| e.to_string())?;
    let p_hyd =
        hydraulic_power_w(s.density_kg_m3, op.flow_m3s, op.head_m).map_err(|e| e.to_string())?;
    let p_shaft = shaft_power_w(s.density_kg_m3, op.flow_m3s, op.head_m, s.efficiency)
        .map_err(|e| e.to_string())?;

    let suction = SuctionConditions::new(
        s.atmospheric_pa,
        s.vapor_pressure_pa,
        s.density_kg_m3,
        s.static_suction_head_m,
        s.suction_loss_m,
    )
    .map_err(|e| e.to_string())?;
    let npsha = available_npsh_m(&suction);
    let margin = npsh_margin_m(&suction, s.required_npsh_m).map_err(|e| e.to_string())?;
    let cav_free = is_cavitation_free(&suction, s.required_npsh_m).map_err(|e| e.to_string())?;

    Ok(format!(
        "pump  H0 / a    : {:.1} m / {:.0}\n\
         system Hs / K   : {:.1} m / {:.0}\n\n\
         operating Q*    : {:.4} m³/s  ({:.1} L/s)\n\
         operating H*    : {:.2} m\n\
         hydraulic power : {:.2} kW\n\
         shaft power     : {:.2} kW  (η {:.0} %)\n\n\
         NPSH available  : {:.2} m\n\
         NPSH required   : {:.2} m\n\
         NPSH margin     : {:.2} m\n\
         cavitation      : {}",
        s.shutoff_head_m,
        s.droop_a,
        s.static_head_m,
        s.resistance_k,
        op.flow_m3s,
        op.flow_m3s * 1000.0,
        op.head_m,
        p_hyd / 1000.0,
        p_shaft / 1000.0,
        s.efficiency * 100.0,
        npsha,
        s.required_npsh_m,
        margin,
        if cav_free {
            "free (margin ≥ 0)"
        } else {
            "RISK (margin < 0)"
        },
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

/// Append a (double-sided) capped-less cylinder whose axis runs along
/// `+x`, spanning `base.x ..= base.x + length` with circle centre
/// `(base.y, base.z)`.
fn push_cyl_x(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    base: Vector3<f64>,
    length: f64,
    r: f64,
    seg: usize,
) {
    let (x0, x1) = (base.x, base.x + length);
    let lo = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(x0, base.y + r * a.cos(), base.z + r * a.sin()));
    }
    let hi = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(x1, base.y + r * a.cos(), base.z + r * a.sin()));
    }
    for j in 0..seg {
        let jn = (j + 1) % seg;
        tris.extend_from_slice(&[
            lo + j,
            hi + j,
            hi + jn,
            lo + j,
            hi + jn,
            lo + jn,
            lo + j,
            hi + jn,
            hi + j,
            lo + j,
            lo + jn,
            hi + jn,
        ]);
    }
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

/// Build the pump as a triangle [`Mesh`] — a volute casing (fat cylinder
/// on the pump axis), a suction eye, a vertical discharge nozzle, a motor
/// stub and a baseplate. Representative geometry (the duty point is the
/// `valenx-pump` curve intersection). `None` for a configuration that
/// delivers no flow.
fn pump_solid_mesh(s: &PumpWorkbenchState) -> Option<Mesh> {
    // Gate the 3-D build on a pump that actually delivers flow against
    // the system (shut-off head must exceed the static lift).
    let pump = PumpCurve::new(s.shutoff_head_m, s.droop_a).ok()?;
    let system = SystemCurve::new(s.static_head_m, s.resistance_k).ok()?;
    operating_point(&pump, &system).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Baseplate the pump is mounted on.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.05),
        Vector3::new(0.7, 0.5, 0.05),
    );
    // Pedestal lifting the casing off the plate.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.28),
        Vector3::new(0.25, 0.18, 0.22),
    );
    // Volute casing — fat cylinder on the pump axis (+x).
    push_cyl_x(
        &mut nodes,
        &mut tris,
        Vector3::new(-0.2, 0.0, 0.55),
        0.4,
        0.45,
        28,
    );
    // Suction eye — narrower inlet entering the casing front (−x).
    push_cyl_x(
        &mut nodes,
        &mut tris,
        Vector3::new(-0.75, 0.0, 0.55),
        0.55,
        0.18,
        20,
    );
    // Motor / driver stub behind the casing (+x).
    push_cyl_x(
        &mut nodes,
        &mut tris,
        Vector3::new(0.2, 0.0, 0.55),
        0.7,
        0.3,
        24,
    );
    // Vertical discharge nozzle off the top of the casing (+z).
    push_cyl_z(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.9),
        0.55,
        0.16,
        18,
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-pump");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D pump solid and load it into the central viewport.
fn load_pump_3d(app: &mut ValenxApp) {
    let Some(mesh) = pump_solid_mesh(&app.pump) else {
        app.pump.error =
            Some("pump cannot deliver flow against the system — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<pump>/valenx-pump"),
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
        let s = PumpWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_operating_point_and_npsh() {
        let mut s = PumpWorkbenchState::default();
        run_pump(&mut s);
        assert!(
            s.error.is_none(),
            "default pump should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("operating Q*"));
        assert!(s.result.contains("shaft power"));
        assert!(s.result.contains("NPSH margin"));
        // Default duty: Q* = sqrt((50-10)/(1000+4000)) ≈ 0.0894 m³/s.
        assert!(s.result.contains("0.0894"));
    }

    #[test]
    fn analyze_rejects_pump_that_cannot_deliver() {
        // Static lift above the shut-off head → no operating point.
        let mut s = PumpWorkbenchState {
            static_head_m: 60.0,
            ..Default::default()
        };
        run_pump(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn analyze_rejects_efficiency_above_one() {
        let mut s = PumpWorkbenchState {
            efficiency: 1.5,
            ..Default::default()
        };
        run_pump(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn pump_mesh_for_default_is_nonempty_and_in_range() {
        let s = PumpWorkbenchState::default();
        let mesh = pump_solid_mesh(&s).expect("default pump yields a solid");
        assert!(mesh.nodes.len() > 8, "expected casing + pipes + base");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn pump_mesh_none_when_no_flow() {
        let s = PumpWorkbenchState {
            static_head_m: 60.0,
            ..Default::default()
        };
        assert!(pump_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_pump_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_pump_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_pump_workbench = true;
        run_pump(&mut app.pump);
        draw_workbench(&mut app);
    }
}
