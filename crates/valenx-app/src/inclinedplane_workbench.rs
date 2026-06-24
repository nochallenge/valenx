//! The right-side **Inclined Plane Workbench** panel — native ramp
//! statics over `valenx-inclinedplane`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_inclinedplane_workbench`,
//! toggled from the View menu. The form sets a ramp (slope angle, friction
//! coefficient, block weight) and chooses whether the load is being raised
//! or lowered; "Analyze" runs the closed-form rigid-body equilibrium of a
//! block on the plane and reports the normal / slope / friction forces, the
//! effort to raise and to lower the load, the friction angle, whether the
//! ramp is self-locking, and the ideal vs actual mechanical advantage with
//! the resulting efficiency. "Show 3-D" loads a representative wedge
//! (triangular-prism) ramp solid into the central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_inclinedplane::Ramp;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Whether the reported headline effort is for moving the load up or down
/// the ramp. Both directions are always shown in the readout; this only
/// selects which one the summary line leads with.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum LoadDirection {
    /// Push the block steadily up the slope (effort adds friction).
    Raise,
    /// Lower / restrain the block down the slope (friction subtracts).
    Lower,
}

/// Persistent form + result state for the Inclined Plane Workbench.
pub struct InclinedPlaneWorkbenchState {
    /// Slope angle from the horizontal (degrees), in `(0, 90)`.
    angle_deg: f64,
    /// Coulomb friction coefficient `mu` between block and ramp.
    mu: f64,
    /// Block weight `W` (newtons).
    weight_n: f64,
    /// Target vertical lift `h` (metres) used to report the incline
    /// (hypotenuse) distance the effort acts over, `h / sin(θ)`.
    lift_height_m: f64,
    /// Whether the headline effort is for raising or lowering the load.
    direction: LoadDirection,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D ramp solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for InclinedPlaneWorkbenchState {
    fn default() -> Self {
        // A 30-degree loading ramp, mu = 0.2, 100 N block: not self-locking,
        // effort to raise ~67.3 N, efficiency ~0.74. Lifting 1 m up a 30°
        // ramp travels 1 / sin30 = 2.000 m along the incline.
        Self {
            angle_deg: 30.0,
            mu: 0.2,
            weight_n: 100.0,
            lift_height_m: 1.0,
            direction: LoadDirection::Raise,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Inclined Plane Workbench right-side panel. A no-op when the
/// `show_inclinedplane_workbench` toggle is off.
pub fn draw_inclinedplane_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_inclinedplane_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_inclinedplane_workbench",
        "Inclined Plane",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native ramp statics & mechanical advantage · valenx-inclinedplane",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.inclinedplane;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Ramp").strong());
                    ui.horizontal(|ui| {
                        let cap = ui.label("slope angle (°)");
                        ui.add(egui::DragValue::new(&mut s.angle_deg).speed(0.5))
                            .labelled_by(cap.id);
                    });
                    ui.horizontal(|ui| {
                        let cap = ui.label("friction μ");
                        ui.add(egui::DragValue::new(&mut s.mu).speed(0.01))
                            .labelled_by(cap.id);
                    });
                    ui.horizontal(|ui| {
                        let cap = ui.label("block weight W (N)");
                        ui.add(egui::DragValue::new(&mut s.weight_n).speed(1.0))
                            .labelled_by(cap.id);
                    });
                    ui.horizontal(|ui| {
                        let cap = ui.label("lift height h (m)");
                        ui.add(egui::DragValue::new(&mut s.lift_height_m).speed(0.1))
                            .labelled_by(cap.id);
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Direction").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.direction, LoadDirection::Raise, "raise");
                        ui.radio_value(&mut s.direction, LoadDirection::Lower, "lower");
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_inclinedplane(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative wedge (triangular-prism) ramp solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Statics").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_inclinedplane_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.inclinedplane` borrow is
    // released here): build the ramp's 3-D solid and load it.
    if app.inclinedplane.show_3d_request {
        app.inclinedplane.show_3d_request = false;
        load_ramp_3d(app);
    }
}

/// Validate the form, evaluate the ramp and format the readout.
fn run_inclinedplane(s: &mut InclinedPlaneWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build a validated [`Ramp`] from the form, mapping any domain error to a
/// display string. Extracted so the readout and the 3-D gate share it.
fn build_ramp(s: &InclinedPlaneWorkbenchState) -> Result<Ramp, String> {
    Ramp::new(s.angle_deg.to_radians(), s.mu, s.weight_n).map_err(|e| e.to_string())
}

/// Evaluate the ramp and format the full readout, mapping any domain error
/// to a display string. Extracted so it is unit-testable.
fn compute(s: &InclinedPlaneWorkbenchState) -> Result<String, String> {
    let ramp = build_ramp(s)?;
    let f = ramp.forces();
    let ima = ramp.geometry().ideal_mechanical_advantage();
    let ama = ramp.actual_mechanical_advantage();
    let eta = ramp.efficiency();
    let phi_deg = ramp.friction_angle_rad().to_degrees();
    let theta_deg = ramp.geometry().angle_deg();
    // Incline (hypotenuse) distance the effort travels to gain `h`:
    // slope length = h / sin(θ) = h · ideal MA.
    let lift_h = s.lift_height_m;
    let slope_len = ramp
        .geometry()
        .slope_length_for_height(lift_h)
        .map_err(|e| e.to_string())?;

    let lock = if f.is_self_locking {
        "yes (holds under gravity)"
    } else {
        "no (slides without a hold)"
    };
    let headline = match s.direction {
        LoadDirection::Raise => format!("effort to RAISE   : {:.2} N", f.effort_to_raise),
        LoadDirection::Lower => format!("effort to LOWER   : {:.2} N", f.effort_to_lower),
    };

    Ok(format!(
        "slope angle θ     : {theta_deg:.2}°\n\
         friction μ / φ    : {:.3} / {phi_deg:.2}°\n\
         block weight W    : {:.2} N\n\n\
         normal force N    : {:.2} N\n\
         slope force       : {:.2} N\n\
         max friction μN   : {:.2} N\n\
         {headline}\n\
         effort to raise   : {:.2} N\n\
         effort to lower   : {:.2} N\n\
         self-locking      : {lock}\n\n\
         ideal MA (1/sinθ) : {ima:.3}\n\
         actual MA (W/F↑)  : {ama:.3}\n\
         efficiency η      : {eta:.3}\n\
         incline for h={lift_h:.2} m : {slope_len:.3} m",
        s.mu, s.weight_n, f.normal, f.slope_force, f.friction, f.effort_to_raise, f.effort_to_lower,
    ))
}

/// Append a wedge (right-triangular prism) to the buffers: a ramp whose
/// vertical right-angle face rises at the back, the hypotenuse forms the
/// sloped top, and the prism is extruded by `width` along `y`. `base` is
/// the horizontal run, `rise` is the vertical gain at the back.
fn push_wedge(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    base: f64,
    rise: f64,
    width: f64,
) {
    let i0 = nodes.len();
    let hw = width / 2.0;
    // Two triangular end caps in the x-z plane, offset along y by +/- hw.
    // Cap vertices (per side): A=(0,0) front-bottom, B=(base,0) back-bottom,
    // C=(base,rise) back-top.
    nodes.push(Vector3::new(0.0, -hw, 0.0)); // 0 near A
    nodes.push(Vector3::new(base, -hw, 0.0)); // 1 near B
    nodes.push(Vector3::new(base, -hw, rise)); // 2 near C
    nodes.push(Vector3::new(0.0, hw, 0.0)); // 3 far A
    nodes.push(Vector3::new(base, hw, 0.0)); // 4 far B
    nodes.push(Vector3::new(base, hw, rise)); // 5 far C

    // Faces wound outward. The two triangular caps stay single Tri3; the
    // three quad faces (bottom, vertical back, sloped top) each split into
    // two Tri3.
    // Near triangular cap (-y).
    tris.extend_from_slice(&[i0, i0 + 2, i0 + 1]);
    // Far triangular cap (+y).
    tris.extend_from_slice(&[i0 + 3, i0 + 4, i0 + 5]);
    // Bottom rectangle (z = 0): A B B' A'.
    tris.extend_from_slice(&[i0, i0 + 1, i0 + 4, i0, i0 + 4, i0 + 3]);
    // Back vertical rectangle (x = base): B C C' B'.
    tris.extend_from_slice(&[i0 + 1, i0 + 2, i0 + 5, i0 + 1, i0 + 5, i0 + 4]);
    // Sloped hypotenuse rectangle (top): A C C' A'.
    tris.extend_from_slice(&[i0, i0 + 5, i0 + 2, i0, i0 + 3, i0 + 5]);
}

/// Build the ramp as a triangle [`Mesh`] — a wedge (right-triangular prism)
/// whose slope matches the configured angle. Representative geometry (unit
/// run; the statics numbers are the `valenx-inclinedplane` result). `None`
/// for an invalid configuration.
fn ramp_solid_mesh(s: &InclinedPlaneWorkbenchState) -> Option<Mesh> {
    let ramp = build_ramp(s).ok()?;
    let theta = ramp.geometry().angle_rad();

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Unit horizontal run; the rise follows the slope so the wedge angle is
    // visually the real theta. Pinned locals for the float ops.
    let base: f64 = 1.0;
    let rise: f64 = base * theta.tan();
    let width: f64 = 0.6;
    push_wedge(&mut nodes, &mut tris, base, rise, width);

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-inclinedplane");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D ramp solid and load it into the central viewport.
fn load_ramp_3d(app: &mut ValenxApp) {
    let Some(mesh) = ramp_solid_mesh(&app.inclinedplane) else {
        app.inclinedplane.error =
            Some("ramp parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<ramp>/valenx-inclinedplane"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// Agent-bridge product: the canonical inclinedplane workbench as a 3-D solid plus its
/// `compute()` readout rows (see [`crate::products_registry`]).
pub(crate) fn inclinedplane_product() -> crate::WorkspaceProduct {
    let s = InclinedPlaneWorkbenchState::default();
    let mesh = ramp_solid_mesh(&s).expect("canonical inclinedplane ⇒ ramp solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<inclinedplane>/valenx-ramp");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical inclinedplane ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Inclined plane (friction)".into(),
        lines,
        mesh: Some(loaded),
        vertex_colors: None,
        camera,
        kind2d: None,
        last_export: None,
        image: None,
        image_texture: None,
        animation: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_is_idle() {
        let s = InclinedPlaneWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_forces_and_advantage() {
        let mut s = InclinedPlaneWorkbenchState::default();
        run_inclinedplane(&mut s);
        assert!(
            s.error.is_none(),
            "default ramp should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("normal force N"));
        assert!(s.result.contains("effort to raise"));
        assert!(s.result.contains("efficiency"));
        // 30 deg, mu=0.2, 100 N: raise ~67.32 N, efficiency ~0.743.
        assert!(s.result.contains("67.32"));
        assert!(s.result.contains("0.743"));
        // Not self-locking at 30 deg with mu = 0.2.
        assert!(s.result.contains("no (slides without a hold)"));
    }

    #[test]
    fn lower_direction_changes_headline() {
        let mut s = InclinedPlaneWorkbenchState {
            direction: LoadDirection::Lower,
            ..Default::default()
        };
        run_inclinedplane(&mut s);
        assert!(s.error.is_none());
        assert!(s.result.contains("effort to LOWER"));
    }

    #[test]
    fn analyze_rejects_vertical_angle() {
        let mut s = InclinedPlaneWorkbenchState {
            angle_deg: 90.0,
            ..Default::default()
        };
        run_inclinedplane(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn ideal_ma_is_one_over_sin_theta_ground_truth() {
        // Ground truth: at theta = 30 deg, sin(theta) = 0.5, so the ideal
        // mechanical advantage 1/sin(theta) is exactly 2.
        let s = InclinedPlaneWorkbenchState::default();
        let ramp = build_ramp(&s).expect("default ramp builds");
        let ima = ramp.geometry().ideal_mechanical_advantage();
        let hand = 1.0 / (30.0_f64.to_radians().sin());
        assert!((ima - hand).abs() < 1e-12);
        assert!((ima - 2.0).abs() < 1e-12);
    }

    #[test]
    fn incline_distance_for_lift_height_ground_truth() {
        // Ground truth: lifting h = 1 m up a 30° ramp travels
        // h / sin(30°) = 1 / 0.5 = 2.000 m along the incline, which is
        // exactly h times the ideal mechanical advantage (1/sin θ = 2).
        let mut s = InclinedPlaneWorkbenchState::default();
        run_inclinedplane(&mut s);
        assert!(
            s.error.is_none(),
            "default ramp should analyze: {:?}",
            s.error
        );
        assert!(
            s.result.contains("incline for h=1.00 m : 2.000 m"),
            "expected the 2.000 m incline line, got:\n{result}",
            result = s.result
        );

        // Hand check against the crate, independent of formatting.
        let ramp = build_ramp(&s).expect("default ramp builds");
        let len = ramp
            .geometry()
            .slope_length_for_height(s.lift_height_m)
            .expect("non-negative height");
        let hand = 1.0_f64 / 30.0_f64.to_radians().sin();
        assert!((len - hand).abs() < 1e-12, "len={len}, hand={hand}");
        assert!((len - 2.0).abs() < 1e-12, "len={len}");
    }

    #[test]
    fn ramp_mesh_for_default_is_nonempty_and_in_range() {
        let s = InclinedPlaneWorkbenchState::default();
        let mesh = ramp_solid_mesh(&s).expect("default ramp yields a solid");
        assert_eq!(mesh.nodes.len(), 6, "a wedge prism has six vertices");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn ramp_mesh_none_for_invalid() {
        let s = InclinedPlaneWorkbenchState {
            angle_deg: 90.0,
            ..Default::default()
        };
        assert!(ramp_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use egui::accesskit::{Node, NodeId, Role};

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_inclinedplane_workbench(app, ctx);
        });
    }

    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_inclinedplane_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_inclinedplane_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_inclinedplane_workbench = true;
        run_inclinedplane(&mut app.inclinedplane);
        draw_workbench(&mut app);
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        let mut app = ValenxApp::default();
        app.show_inclinedplane_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);
        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        assert!(
            spin_buttons.len() >= 4,
            "expected the numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every DragValue must be labelled_by its caption (AI-drivable name)"
        );
        for caption in ["slope angle (°)", "block weight W (N)"] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
    }
}
