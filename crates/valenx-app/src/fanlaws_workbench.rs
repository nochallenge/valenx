//! The right-side **Fan Laws Workbench** panel — native fan / blower
//! affinity-law scaling over `valenx-fanlaws`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_fanlaws_workbench`,
//! toggled from the View menu. The form sets a baseline fan operating
//! point (flow, pressure rise, total efficiency, speed) and a new
//! impeller speed; "Analyze" applies the affinity laws — flow linear in
//! speed, pressure with its square, shaft power with its cube — via
//! [`scale_operating_point`] and reports the scaled operating point and
//! the [`air_power`] / [`shaft_power`] split, and "Show 3-D fan" loads a
//! representative hub-and-blades impeller solid into the central viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_fanlaws::{air_power, scale_operating_point, shaft_power, Efficiency, OperatingPoint};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Number of radial blades drawn on the representative 3-D impeller.
const BLADE_COUNT: usize = 6;

/// Persistent form + result state for the Fan Laws Workbench.
pub struct FanLawsWorkbenchState {
    /// Baseline volumetric flow `Q1` (m^3/s).
    flow_m3s: f64,
    /// Baseline pressure rise `dP1` (Pa).
    pressure_pa: f64,
    /// Fan total efficiency (percent, in (0, 100]); held constant across
    /// the speed step.
    efficiency_percent: f64,
    /// Baseline impeller speed `N1` (rev/min).
    speed1_rpm: f64,
    /// New impeller speed `N2` (rev/min) to scale to.
    speed2_rpm: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D fan solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for FanLawsWorkbenchState {
    fn default() -> Self {
        // A fan delivering 2 m^3/s at a 500 Pa rise (1000 W of air power)
        // at ~66.7% total efficiency draws ~1.5 kW shaft at 1000 rev/min;
        // stepping to 1500 rev/min (r = 1.5) gives Q=3 m^3/s, dP=1125 Pa,
        // and shaft power x 1.5^3 = x3.375 (~5.06 kW).
        Self {
            flow_m3s: 2.0,
            pressure_pa: 500.0,
            efficiency_percent: 66.7,
            speed1_rpm: 1000.0,
            speed2_rpm: 1500.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Fan Laws Workbench right-side panel. A no-op when the
/// `show_fanlaws_workbench` toggle is off.
pub fn draw_fanlaws_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_fanlaws_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_fanlaws_workbench",
        "Fan Laws",
        |app, ui| {
            ui.label(
                egui::RichText::new("native fan / blower affinity-law scaling · valenx-fanlaws")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.fanlaws;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Baseline operating point").strong());
                    ui.horizontal(|ui| {
                        ui.label("flow Q₁ (m³/s)");
                        ui.add(egui::DragValue::new(&mut s.flow_m3s).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("pressure ΔP₁ (Pa)");
                        ui.add(egui::DragValue::new(&mut s.pressure_pa).speed(5.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("efficiency η (%)");
                        ui.add(egui::DragValue::new(&mut s.efficiency_percent).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("speed N₁ (rev/min)");
                        ui.add(egui::DragValue::new(&mut s.speed1_rpm).speed(10.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Scale to").strong());
                    ui.horizontal(|ui| {
                        ui.label("new speed N₂ (rev/min)");
                        ui.add(egui::DragValue::new(&mut s.speed2_rpm).speed(10.0));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_fanlaws(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D fan").strong())
                        .on_hover_text(
                            "Build a representative impeller (central hub + radial blades) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Affinity-law scaling").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_fanlaws_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.fanlaws` borrow is
    // released here): build the fan's 3-D solid and load it.
    if app.fanlaws.show_3d_request {
        app.fanlaws.show_3d_request = false;
        load_fan_3d(app);
    }
}

/// Validate the form, apply the affinity laws and format the readout.
fn run_fanlaws(s: &mut FanLawsWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// The baseline [`OperatingPoint`] and its [`Efficiency`] — the validated
/// inputs both the readout and the 3-D gate need. Extracted so it is
/// unit-testable and shared.
///
/// Density is dimension-agnostic in `valenx-fanlaws` and is held equal
/// across the speed step here (only the speed changes), so the pressure
/// and power scale purely by the speed ratio. A nominal 1.2 (kg/m^3 air)
/// is used for both points.
fn baseline(s: &FanLawsWorkbenchState) -> Result<(OperatingPoint, Efficiency), String> {
    let eta = Efficiency::from_percent(s.efficiency_percent).map_err(|e| e.to_string())?;
    let w1 = shaft_power(s.flow_m3s, s.pressure_pa, eta).map_err(|e| e.to_string())?;
    let point = OperatingPoint::new(s.flow_m3s, s.pressure_pa, w1, s.speed1_rpm, 1.2)
        .map_err(|e| e.to_string())?;
    Ok((point, eta))
}

/// Apply the affinity laws and format the full readout, mapping any
/// domain error to a display string. Extracted so it is unit-testable.
fn compute(s: &FanLawsWorkbenchState) -> Result<String, String> {
    let (p0, eta) = baseline(s)?;
    // Same density across the step: only the speed changes.
    let p1 = scale_operating_point(&p0, s.speed2_rpm, p0.density).map_err(|e| e.to_string())?;
    let r = s.speed2_rpm / s.speed1_rpm;
    let air1 = air_power(p0.flow, p0.pressure).map_err(|e| e.to_string())?;
    let air2 = air_power(p1.flow, p1.pressure).map_err(|e| e.to_string())?;

    Ok(format!(
        "speed N₁ → N₂   : {n1:.0} → {n2:.0} rev/min\n\
         speed ratio r   : {r:.4}\n\
         efficiency η    : {eta:.1} %\n\n\
         flow    Q₁ → Q₂ : {q1:.3} → {q2:.3} m³/s   (× r)\n\
         press   P₁ → P₂ : {p1v:.1} → {p2v:.1} Pa     (× r²)\n\
         air pwr         : {a1:.1} → {a2:.1} W\n\
         shaft   W₁ → W₂ : {w1:.1} → {w2:.1} W      (× r³)",
        n1 = s.speed1_rpm,
        n2 = s.speed2_rpm,
        r = r,
        eta = eta.percent(),
        q1 = p0.flow,
        q2 = p1.flow,
        p1v = p0.pressure,
        p2v = p1.pressure,
        a1 = air1,
        a2 = air2,
        w1 = p0.power,
        w2 = p1.power,
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

/// Append a faceted cylinder of `radius` and half-length `half_len`, with
/// its axis along `x`, centred at the origin. Drawn as a closed
/// `segments`-sided prism (side quads + two end caps).
fn push_x_cylinder(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    radius: f64,
    half_len: f64,
    segments: usize,
) {
    let ring0 = nodes.len();
    for k in 0..segments {
        let theta = TAU * (k as f64) / (segments as f64);
        let (y, z) = (radius * theta.cos(), radius * theta.sin());
        nodes.push(Vector3::new(-half_len, y, z));
        nodes.push(Vector3::new(half_len, y, z));
    }
    // Side quads between consecutive rings (2 verts per segment).
    for k in 0..segments {
        let a = ring0 + 2 * k;
        let b = ring0 + 2 * ((k + 1) % segments);
        // (-x of a, +x of a, +x of b, -x of b)
        tris.extend_from_slice(&[a, a + 1, b + 1, a, b + 1, b]);
    }
    // End-cap centres + triangle fans.
    let cap_neg = nodes.len();
    nodes.push(Vector3::new(-half_len, 0.0, 0.0));
    let cap_pos = nodes.len();
    nodes.push(Vector3::new(half_len, 0.0, 0.0));
    for k in 0..segments {
        let a = ring0 + 2 * k;
        let b = ring0 + 2 * ((k + 1) % segments);
        tris.extend_from_slice(&[cap_neg, b, a]);
        tris.extend_from_slice(&[cap_pos, a + 1, b + 1]);
    }
}

/// Build the fan as a triangle [`Mesh`] — a central hub (a cylinder on
/// the `x` spin axis) with [`BLADE_COUNT`] flat blade boxes radiating
/// around it. Representative geometry (not to scale; the affinity-law
/// numbers are the `valenx-fanlaws` result). `None` for an invalid
/// configuration.
fn fan_solid_mesh(s: &FanLawsWorkbenchState) -> Option<Mesh> {
    // Gate on the real crate object: an invalid baseline yields no solid.
    baseline(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Central hub on the spin (x) axis.
    push_x_cylinder(&mut nodes, &mut tris, 0.18, 0.16, 16);

    // Radial blades: thin boxes built around the +y axis, then rotated
    // about x to their angular station.
    let blade_half = Vector3::new(0.02, 0.5, 0.16);
    let blade_centre_r = 0.5; // mid-radius of the blade box along y
    for i in 0..BLADE_COUNT {
        let phi = TAU * (i as f64) / (BLADE_COUNT as f64);
        let (cphi, sphi) = (phi.cos(), phi.sin());
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
            // Box centred at (0, blade_centre_r, 0) before rotation.
            let local = Vector3::new(
                sx * blade_half.x,
                blade_centre_r + sy * blade_half.y,
                sz * blade_half.z,
            );
            // Rotate about the x (spin) axis by phi.
            let rotated = Vector3::new(
                local.x,
                local.y * cphi - local.z * sphi,
                local.y * sphi + local.z * cphi,
            );
            nodes.push(rotated);
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

    // A short back plate behind the hub for orientation.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(-0.2, 0.0, 0.0),
        Vector3::new(0.03, 0.22, 0.22),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-fanlaws");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D fan solid and load it into the central viewport.
fn load_fan_3d(app: &mut ValenxApp) {
    let Some(mesh) = fan_solid_mesh(&app.fanlaws) else {
        app.fanlaws.error = Some("fan parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<fan>/valenx-fanlaws"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// Agent-bridge product: the canonical fan-laws workbench as a 3-D solid plus
/// its `compute()` readout rows (see [`crate::products_registry`]).
pub(crate) fn fanlaws_product() -> crate::WorkspaceProduct {
    let s = FanLawsWorkbenchState::default();
    let mesh = fan_solid_mesh(&s).expect("canonical fan laws ⇒ fan solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<fanlaws>/valenx-fan");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical fan laws ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Fan laws (affinity scaling)".into(),
        lines,
        mesh: Some(loaded),
        vertex_colors: None,
        camera,
        kind2d: None,
        last_export: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_is_idle() {
        let s = FanLawsWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_scaled_flow_pressure_and_power() {
        let mut s = FanLawsWorkbenchState::default();
        run_fanlaws(&mut s);
        assert!(
            s.error.is_none(),
            "default fan should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("speed ratio r"));
        assert!(s.result.contains("flow    Q₁ → Q₂"));
        assert!(s.result.contains("shaft   W₁ → W₂"));
        // r = 1500/1000 = 1.5; flow 2 -> 3 m^3/s, pressure 500 -> 1125 Pa.
        assert!(s.result.contains("1.5000"));
        assert!(s.result.contains("3.000"));
        assert!(s.result.contains("1125.0"));
    }

    #[test]
    fn analyze_rejects_efficiency_over_unity() {
        let mut s = FanLawsWorkbenchState {
            efficiency_percent: 150.0,
            ..Default::default()
        };
        run_fanlaws(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn affinity_exponents_are_one_two_three() {
        // Ground truth: at speed ratio r the affinity laws scale flow by
        // r, pressure by r^2, and shaft power by r^3.
        let s = FanLawsWorkbenchState::default();
        let (p0, _eta) = baseline(&s).expect("default baseline is valid");
        let p1 = scale_operating_point(&p0, s.speed2_rpm, p0.density).unwrap();
        let r = s.speed2_rpm / s.speed1_rpm;
        assert!((p1.flow / p0.flow - r).abs() < 1e-9);
        assert!((p1.pressure / p0.pressure - r * r).abs() < 1e-9);
        assert!((p1.power / p0.power - r * r * r).abs() < 1e-9);
    }

    #[test]
    fn fan_mesh_for_default_is_nonempty_and_in_range() {
        let s = FanLawsWorkbenchState::default();
        let mesh = fan_solid_mesh(&s).expect("default fan yields a solid");
        assert!(mesh.nodes.len() > 8, "expected hub + blades + plate");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn fan_mesh_none_for_invalid() {
        let s = FanLawsWorkbenchState {
            efficiency_percent: 0.0,
            ..Default::default()
        };
        assert!(fan_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_fanlaws_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_fanlaws_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_fanlaws_workbench = true;
        run_fanlaws(&mut app.fanlaws);
        draw_workbench(&mut app);
    }
}
