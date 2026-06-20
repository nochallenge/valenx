//! The right-side **Retaining Wall Workbench** panel — native Rankine
//! lateral earth-pressure analysis over `valenx-retainingwall`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on
//! `crate::ValenxApp::show_retainingwall_workbench`, toggled from the View
//! menu. The form sets a dry, cohesionless soil (friction angle `phi`, unit
//! weight `gamma`) and a wall height `H`, and selects the active or passive
//! limit state; "Analyze" reports the Rankine earth-pressure coefficient,
//! the base pressure, the resultant thrust per metre of wall and its line
//! of action above the base, and "Show 3-D" loads a representative
//! cantilever retaining wall (vertical stem on a base footing) into the
//! central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use valenx_retainingwall::SoilProfile;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Which Rankine limit state to evaluate.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum LimitState {
    /// Active case (wall yields away from the soil): `Ka < 1`.
    Active,
    /// Passive case (wall pushed into the soil): `Kp = 1/Ka > 1`.
    Passive,
}

/// Persistent form + result state for the Retaining Wall Workbench.
pub struct RetainingWallWorkbenchState {
    /// Soil internal friction angle `phi` (degrees), `0 <= phi < 90`.
    phi_deg: f64,
    /// Soil unit weight `gamma` (kN/m^3), strictly positive.
    gamma_kn_per_m3: f64,
    /// Wall height `H` (m), non-negative.
    height_m: f64,
    /// Active vs passive limit state.
    state: LimitState,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D wall solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for RetainingWallWorkbenchState {
    fn default() -> Self {
        // Medium-dense sand (phi = 30 deg, gamma = 18 kN/m^3) retained by a
        // 4 m wall: Ka = 1/3, so the active thrust is
        // 0.5 * (1/3) * 18 * 4^2 = 48 kN/m, acting at H/3 = 1.333 m.
        Self {
            phi_deg: 30.0,
            gamma_kn_per_m3: 18.0,
            height_m: 4.0,
            state: LimitState::Active,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Retaining Wall Workbench right-side panel. A no-op when the
/// `show_retainingwall_workbench` toggle is off.
pub fn draw_retainingwall_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_retainingwall_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_retainingwall_workbench",
        "Retaining Wall",
        |app, ui| {
            ui.label(
                egui::RichText::new("native Rankine lateral earth pressure · valenx-retainingwall")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.retainingwall;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Soil").strong());
                    ui.horizontal(|ui| {
                        ui.label("friction angle φ (°)");
                        ui.add(egui::DragValue::new(&mut s.phi_deg).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("unit weight γ (kN/m³)");
                        ui.add(egui::DragValue::new(&mut s.gamma_kn_per_m3).speed(0.5));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Wall").strong());
                    ui.horizontal(|ui| {
                        ui.label("height H (m)");
                        ui.add(egui::DragValue::new(&mut s.height_m).speed(0.1));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Limit state").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.state, LimitState::Active, "active (Ka)");
                        ui.radio_value(&mut s.state, LimitState::Passive, "passive (Kp)");
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_retainingwall(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative cantilever retaining wall (vertical stem on a base footing) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Earth pressure").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_retainingwall_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.retainingwall` borrow is
    // released here): build the wall's 3-D solid and load it.
    if app.retainingwall.show_3d_request {
        app.retainingwall.show_3d_request = false;
        load_wall_3d(app);
    }
}

/// Validate the form, evaluate the wall and format the readout.
fn run_retainingwall(s: &mut RetainingWallWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Evaluate the Rankine earth pressure for the form and format the full
/// readout, mapping any domain error to a display string. Extracted so it
/// is unit-testable.
fn compute(s: &RetainingWallWorkbenchState) -> Result<String, String> {
    let soil = SoilProfile::new(s.phi_deg, s.gamma_kn_per_m3).map_err(|e| e.to_string())?;

    let (label, coeff, base_pressure, thrust) = match s.state {
        LimitState::Active => (
            "active",
            soil.ka(),
            soil.active_pressure_at(s.height_m)
                .map_err(|e| e.to_string())?,
            soil.active_thrust(s.height_m).map_err(|e| e.to_string())?,
        ),
        LimitState::Passive => (
            "passive",
            soil.kp(),
            soil.passive_pressure_at(s.height_m)
                .map_err(|e| e.to_string())?,
            soil.passive_thrust(s.height_m).map_err(|e| e.to_string())?,
        ),
    };

    Ok(format!(
        "limit state     : {label}\n\
         friction angle φ: {:.1} °\n\
         unit weight γ   : {:.1} kN/m³\n\
         wall height H   : {:.2} m\n\n\
         coefficient K   : {:.4}\n\
         base pressure σ : {:.3} kPa\n\
         thrust P        : {:.3} kN/m\n\
         line of action  : {:.4} m above base",
        s.phi_deg,
        s.gamma_kn_per_m3,
        s.height_m,
        coeff,
        base_pressure,
        thrust.resultant,
        thrust.line_of_action,
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

/// Build a representative cantilever retaining wall as a triangle [`Mesh`]
/// — a vertical stem (height scaled to `H`) sitting on a wider base
/// footing. Representative geometry (not a stability-checked section; the
/// earth-pressure numbers are the `valenx-retainingwall` result). `None`
/// for an invalid configuration.
fn wall_solid_mesh(s: &RetainingWallWorkbenchState) -> Option<Mesh> {
    // Reuse the same validation the analyze path uses.
    SoilProfile::new(s.phi_deg, s.gamma_kn_per_m3).ok()?;

    // Scale the drawn stem height to H, clamped to a sane visual range so a
    // tiny or huge H still renders a recognisable wall.
    let stem_h = (s.height_m * 0.5).clamp(0.4, 4.0);
    let base_th = 0.15;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Base footing: wide, thin slab on the ground.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, base_th),
        Vector3::new(0.6, 0.5, base_th),
    );
    // Vertical stem: thin in the through-thickness (x) direction, rising to
    // stem_h above the footing.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(-0.4, 0.0, 2.0 * base_th + stem_h),
        Vector3::new(0.08, 0.5, stem_h),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-retainingwall");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D wall solid and load it into the central viewport.
fn load_wall_3d(app: &mut ValenxApp) {
    let Some(mesh) = wall_solid_mesh(&app.retainingwall) else {
        app.retainingwall.error =
            Some("wall parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<wall>/valenx-retainingwall"),
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
    use std::f64::consts::PI;

    #[test]
    fn default_state_is_idle() {
        let s = RetainingWallWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_coefficient_thrust_and_line_of_action() {
        let mut s = RetainingWallWorkbenchState::default();
        run_retainingwall(&mut s);
        assert!(
            s.error.is_none(),
            "default wall should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("coefficient K"));
        assert!(s.result.contains("thrust P"));
        assert!(s.result.contains("line of action"));
        // phi = 30 deg => Ka = 1/3 = 0.3333 at the readout's {:.4} precision.
        assert!(s.result.contains("0.3333"));
        // Pa = 0.5 * (1/3) * 18 * 4^2 = 48.000 kN/m at the {:.3} precision.
        assert!(s.result.contains("48.000"));
        // Line of action H/3 = 1.3333 m at the {:.4} precision.
        assert!(s.result.contains("1.3333"));
    }

    #[test]
    fn analyze_passive_reports_passive_coefficient() {
        let mut s = RetainingWallWorkbenchState {
            state: LimitState::Passive,
            ..Default::default()
        };
        run_retainingwall(&mut s);
        assert!(
            s.error.is_none(),
            "passive case should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("passive"));
        // phi = 30 deg => Kp = 3 = 3.0000 at the readout's {:.4} precision.
        assert!(s.result.contains("3.0000"));
    }

    #[test]
    fn analyze_rejects_friction_angle_at_ninety() {
        let mut s = RetainingWallWorkbenchState {
            phi_deg: 90.0,
            ..Default::default()
        };
        run_retainingwall(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn ground_truth_ka_and_active_thrust_at_thirty_degrees() {
        // Ground truth: Ka = tan^2(45 - phi/2) = (1 - sin phi)/(1 + sin phi)
        // = 1/3 at phi = 30 deg, and the active thrust on a wall of height H
        // is the area of the triangular pressure diagram,
        // Pa = 0.5 * Ka * gamma * H^2 = 0.5 * (1/3) * 18 * 16 = 48 kN/m.
        let soil = SoilProfile::new(30.0, 18.0).expect("valid soil");

        let phi_rad = 30.0 * (PI / 180.0);
        let sin_phi = phi_rad.sin();
        let ka_identity = (1.0 - sin_phi) / (1.0 + sin_phi);
        assert!((soil.ka() - ka_identity).abs() < 1e-12);
        assert!((soil.ka() - 1.0 / 3.0).abs() < 1e-12);

        let thrust = soil.active_thrust(4.0).expect("valid height");
        let expected = 0.5 * (1.0 / 3.0) * 18.0 * 16.0;
        assert!((thrust.resultant - expected).abs() < 1e-9);
        assert!((thrust.resultant - 48.0).abs() < 1e-9);
        assert!((thrust.line_of_action - 4.0 / 3.0).abs() < 1e-9);
    }

    #[test]
    fn wall_mesh_for_default_is_nonempty_and_in_range() {
        let s = RetainingWallWorkbenchState::default();
        let mesh = wall_solid_mesh(&s).expect("default wall yields a solid");
        assert!(mesh.nodes.len() > 8, "expected stem + base footing");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn wall_mesh_none_for_invalid() {
        let s = RetainingWallWorkbenchState {
            gamma_kn_per_m3: 0.0,
            ..Default::default()
        };
        assert!(wall_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_retainingwall_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_retainingwall_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_retainingwall_workbench = true;
        run_retainingwall(&mut app.retainingwall);
        draw_workbench(&mut app);
    }
}
