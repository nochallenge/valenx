//! The right-side **Hydraulics Workbench** panel — native hydraulic
//! cylinder / circuit analysis over `valenx-hydraulics`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_hydraulics_workbench`,
//! toggled from the View menu. The form sets a double-acting cylinder (bore
//! and rod diameter), a supply pressure and piston speed, and a control
//! valve (flow coefficient, pressure drop, fluid specific gravity).
//! "Analyze" reports the bore / rod areas, the developed force `F = p A`
//! for the selected stroke, the required flow `Q = A v`, the transmitted
//! hydraulic power `Power = p Q`, and the valve flow `Q = Cv sqrt(dP / SG)`;
//! "Show 3-D" loads a representative actuator solid (barrel with a
//! protruding rod) into the central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_hydraulics::{hydraulic_power, valve_flow, Cylinder, Stroke};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Hydraulics Workbench.
pub struct HydraulicsWorkbenchState {
    /// Cylinder bore (piston) diameter `D` (m).
    bore_diameter_m: f64,
    /// Cylinder rod diameter `d` (m).
    rod_diameter_m: f64,
    /// Supply gauge pressure `p` (bar; converted to Pa for the solver).
    pressure_bar: f64,
    /// Commanded piston speed `v` (m/s).
    velocity_m_s: f64,
    /// Which stroke (extend uses the bore area, retract the rod annulus).
    stroke: Stroke,
    /// Valve flow coefficient `Cv` (US units; flow in gpm at 1 psi water).
    valve_cv: f64,
    /// Pressure drop across the valve `dP` (psi).
    valve_dp_psi: f64,
    /// Fluid specific gravity `SG` (dimensionless, water = 1).
    specific_gravity: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D actuator solid (serviced after
    /// the panel draws).
    show_3d_request: bool,
}

impl Default for HydraulicsWorkbenchState {
    fn default() -> Self {
        // A 100 mm bore / 50 mm rod double-acting cylinder at 210 bar
        // driving out at 0.2 m/s, with a Cv = 12 control valve dropping
        // 4 psi on hydraulic oil (SG ~ 0.87). Extend force ~ 164.9 kN.
        Self {
            bore_diameter_m: 0.100,
            rod_diameter_m: 0.050,
            pressure_bar: 210.0,
            velocity_m_s: 0.20,
            stroke: Stroke::Extend,
            valve_cv: 12.0,
            valve_dp_psi: 4.0,
            specific_gravity: 0.87,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Hydraulics Workbench right-side panel. A no-op when the
/// `show_hydraulics_workbench` toggle is off.
pub fn draw_hydraulics_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_hydraulics_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_hydraulics_workbench",
        "Hydraulics",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native cylinder / circuit force-flow-power · valenx-hydraulics",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.hydraulics;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Cylinder").strong());
                    ui.horizontal(|ui| {
                        ui.label("bore D (m)");
                        ui.add(egui::DragValue::new(&mut s.bore_diameter_m).speed(0.005));
                    });
                    ui.horizontal(|ui| {
                        ui.label("rod d (m)");
                        ui.add(egui::DragValue::new(&mut s.rod_diameter_m).speed(0.005));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Operating point").strong());
                    ui.horizontal(|ui| {
                        ui.label("pressure (bar)");
                        ui.add(egui::DragValue::new(&mut s.pressure_bar).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("speed v (m/s)");
                        ui.add(egui::DragValue::new(&mut s.velocity_m_s).speed(0.01));
                    });
                    ui.horizontal(|ui| {
                        ui.label("stroke");
                        ui.radio_value(&mut s.stroke, Stroke::Extend, "extend");
                        ui.radio_value(&mut s.stroke, Stroke::Retract, "retract");
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Control valve").strong());
                    ui.horizontal(|ui| {
                        ui.label("Cv");
                        ui.add(egui::DragValue::new(&mut s.valve_cv).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("ΔP (psi)");
                        ui.add(egui::DragValue::new(&mut s.valve_dp_psi).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("specific gravity");
                        ui.add(egui::DragValue::new(&mut s.specific_gravity).speed(0.01));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_hydraulics(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative double-acting actuator (barrel with a protruding rod) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Performance").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_hydraulics_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.hydraulics` borrow is
    // released here): build the actuator's 3-D solid and load it.
    if app.hydraulics.show_3d_request {
        app.hydraulics.show_3d_request = false;
        load_cylinder_3d(app);
    }
}

/// Validate the form, evaluate the circuit and format the readout.
fn run_hydraulics(s: &mut HydraulicsWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build the [`Cylinder`] from the form, mapping any domain error to a
/// display string. Extracted so it is shared by the readout and the 3-D
/// gate.
fn cylinder(s: &HydraulicsWorkbenchState) -> Result<Cylinder, String> {
    Cylinder::new(s.bore_diameter_m, s.rod_diameter_m).map_err(|e| e.to_string())
}

/// Evaluate the cylinder and valve and format the full readout, mapping
/// any domain error to a display string. Extracted so it is unit-testable.
fn compute(s: &HydraulicsWorkbenchState) -> Result<String, String> {
    let cyl = cylinder(s)?;
    let pressure_pa = s.pressure_bar * 1.0e5;

    let bore_area = cyl.bore_area_m2();
    let rod_area = cyl.rod_area_m2();
    let area = cyl.effective_area_m2(s.stroke);

    let force_n = cyl
        .force_n(pressure_pa, s.stroke)
        .map_err(|e| e.to_string())?;
    let flow_m3_s = cyl
        .flow_m3_s(s.velocity_m_s, s.stroke)
        .map_err(|e| e.to_string())?;
    let power_w = hydraulic_power(pressure_pa, flow_m3_s).map_err(|e| e.to_string())?;
    let valve_q =
        valve_flow(s.valve_cv, s.valve_dp_psi, s.specific_gravity).map_err(|e| e.to_string())?;

    // Convenience conversions for the readout (display only).
    let force_kn = force_n / 1.0e3;
    let flow_l_min = flow_m3_s * 60.0e3;
    let power_kw = power_w / 1.0e3;

    Ok(format!(
        "bore D / rod d  : {:.3} / {:.3} m\n\
         bore area A_b   : {:.6} m²\n\
         rod area A_r    : {:.6} m²\n\
         stroke          : {}\n\
         pressure        : {:.1} bar\n\n\
         effective area  : {:.6} m²\n\
         force F = p·A   : {:.3} kN\n\
         flow Q = A·v    : {:.4} m³/s  ({:.2} L/min)\n\
         power P = p·Q   : {:.3} kW\n\n\
         valve Q = Cv·√(ΔP/SG): {:.3}",
        s.bore_diameter_m,
        s.rod_diameter_m,
        bore_area,
        rod_area,
        s.stroke.label(),
        s.pressure_bar,
        area,
        force_kn,
        flow_m3_s,
        flow_l_min,
        power_kw,
        valve_q,
    ))
}

/// Append a capped cylinder of radius `radius` running along `x` from
/// `x0` to `x1` (so `segments` quad facets around the circumference,
/// triangulated, plus the two fan-triangulated end caps) to the buffers.
fn push_cylinder(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    x0: f64,
    x1: f64,
    radius: f64,
    segments: usize,
) {
    let base = nodes.len();
    // Two rings of `segments` vertices, then the two cap centres.
    for &x in &[x0, x1] {
        for k in 0..segments {
            let theta = std::f64::consts::TAU * (k as f64) / (segments as f64);
            nodes.push(Vector3::new(x, radius * theta.cos(), radius * theta.sin()));
        }
    }
    let centre0 = nodes.len();
    nodes.push(Vector3::new(x0, 0.0, 0.0));
    let centre1 = nodes.len();
    nodes.push(Vector3::new(x1, 0.0, 0.0));

    for k in 0..segments {
        let kn = (k + 1) % segments;
        let a = base + k;
        let b = base + kn;
        let c = base + segments + k;
        let d = base + segments + kn;
        // Side wall (two triangles per facet).
        tris.extend_from_slice(&[a, b, d, a, d, c]);
        // End caps (fan to each centre).
        tris.extend_from_slice(&[centre0, b, a]);
        tris.extend_from_slice(&[centre1, c, d]);
    }
}

/// Build the actuator as a triangle [`Mesh`] — a barrel sized to the bore
/// with a thinner rod sized to the rod diameter protruding from the +x
/// end. Representative geometry (the bore / rod proportions are real; the
/// stroke length is illustrative). `None` for an invalid configuration.
fn cylinder_solid_mesh(s: &HydraulicsWorkbenchState) -> Option<Mesh> {
    let cyl = cylinder(s).ok()?;
    let bore_r = cyl.bore_diameter_m() / 2.0;
    let rod_r = cyl.rod_diameter_m() / 2.0;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Barrel: length ~ 4 bore diameters along x.
    let barrel_len = 4.0 * cyl.bore_diameter_m();
    push_cylinder(&mut nodes, &mut tris, 0.0, barrel_len, bore_r, 24);
    // Rod: protrudes a further ~2 bore diameters from the +x end.
    let rod_len = 2.0 * cyl.bore_diameter_m();
    push_cylinder(
        &mut nodes,
        &mut tris,
        barrel_len,
        barrel_len + rod_len,
        rod_r,
        24,
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-hydraulics");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D actuator solid and load it into the central viewport.
fn load_cylinder_3d(app: &mut ValenxApp) {
    let Some(mesh) = cylinder_solid_mesh(&app.hydraulics) else {
        app.hydraulics.error =
            Some("cylinder parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<actuator>/valenx-hydraulics"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"hydraulics"}`** product: the canonical
/// hydraulic cylinder built as a 3-D solid, paired with the workbench's own
/// `compute()` readout rows, at a fixed 3/4 camera. Registered in
/// [`crate::products_registry`]; the per-tool builder the registry dispatches
/// to. Pure — driven off [`HydraulicsWorkbenchState::default`].
pub(crate) fn hydraulics_product() -> crate::WorkspaceProduct {
    let s = HydraulicsWorkbenchState::default();
    let mesh = cylinder_solid_mesh(&s).expect("canonical hydraulic cylinder ⇒ solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<hydraulics>/valenx-hydraulics");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical hydraulics ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Hydraulic cylinder (circuit)".into(),
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
        let s = HydraulicsWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_force_flow_and_power() {
        let mut s = HydraulicsWorkbenchState::default();
        run_hydraulics(&mut s);
        assert!(
            s.error.is_none(),
            "default cylinder should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("force F = p·A"));
        assert!(s.result.contains("flow Q = A·v"));
        assert!(s.result.contains("power P = p·Q"));
        // 100 mm bore at 210 bar extending: F = 21e6 * pi*0.05^2 ~ 164.9 kN.
        assert!(s.result.contains("164.934"));
    }

    #[test]
    fn analyze_rejects_rod_not_thinner_than_bore() {
        let mut s = HydraulicsWorkbenchState {
            rod_diameter_m: 0.100,
            ..Default::default()
        };
        run_hydraulics(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn force_is_pressure_times_bore_area_on_extend() {
        // Ground truth: on the extend stroke F = p * (pi D^2 / 4),
        // hand-computed for 100 mm bore at 210 bar (21 MPa).
        let s = HydraulicsWorkbenchState::default();
        let cyl = cylinder(&s).unwrap();
        let pressure_pa: f64 = s.pressure_bar * 1.0e5;
        let force = cyl.force_n(pressure_pa, Stroke::Extend).unwrap();
        let expected = pressure_pa * std::f64::consts::PI * 0.100 * 0.100 / 4.0;
        assert!((force - expected).abs() < 1e-6);
        // ~164.9 kN.
        assert!((force - 164_933.614).abs() < 1.0);
    }

    #[test]
    fn extend_force_exceeds_retract_force() {
        let mut extend = HydraulicsWorkbenchState::default();
        run_hydraulics(&mut extend);
        let cyl = cylinder(&extend).unwrap();
        let p: f64 = extend.pressure_bar * 1.0e5;
        let fe = cyl.force_n(p, Stroke::Extend).unwrap();
        let fr = cyl.force_n(p, Stroke::Retract).unwrap();
        assert!(fe > fr, "extend force {fe} should exceed retract {fr}");
    }

    #[test]
    fn cylinder_mesh_for_default_is_nonempty_and_in_range() {
        let s = HydraulicsWorkbenchState::default();
        let mesh = cylinder_solid_mesh(&s).expect("default cylinder yields a solid");
        assert!(mesh.nodes.len() > 8, "expected barrel + rod rings");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn cylinder_mesh_none_for_invalid() {
        let s = HydraulicsWorkbenchState {
            rod_diameter_m: 0.100,
            ..Default::default()
        };
        assert!(cylinder_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_hydraulics_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_hydraulics_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_hydraulics_workbench = true;
        run_hydraulics(&mut app.hydraulics);
        draw_workbench(&mut app);
    }
}
