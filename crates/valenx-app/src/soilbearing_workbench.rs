//! The right-side **Soil Bearing Workbench** panel — native shallow-
//! foundation bearing-capacity analysis over `valenx-soilbearing`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_soilbearing_workbench`,
//! toggled from the View menu. The form sets the soil strength parameters
//! (friction angle, cohesion, unit weight) and a continuous strip-footing
//! geometry (width, founding depth) with a factor of safety; "Analyze"
//! evaluates the classical Terzaghi general-shear bearing capacity
//! (`qult = c*Nc + q*Nq + 0.5*gamma*B*Ngamma`) and reports the
//! bearing-capacity factors, the three term contributions, the ultimate
//! and allowable pressures, and the allowable line load; "Show 3-D" loads
//! a representative footing-on-soil solid into the central viewport.
//!
//! Honest scope follows the crate: a single homogeneous drained stratum
//! under a vertically, concentrically loaded continuous **strip** footing
//! on level ground (no shape/depth/inclination factors, no groundwater,
//! no settlement). The crate exposes only this strip-footing model, so
//! the panel does not offer a footing-shape selector. Research/educational
//! grade — not a substitute for a licensed geotechnical engineer.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_soilbearing::{bearing_capacity, Footing, SoilProperties};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Soil Bearing Workbench.
pub struct SoilBearingWorkbenchState {
    /// Drained angle of internal friction `phi` (degrees).
    friction_angle_deg: f64,
    /// Effective cohesion `c` (kPa).
    cohesion_kpa: f64,
    /// Effective soil unit weight `gamma` (kN/m³).
    unit_weight_kn_m3: f64,
    /// Footing width `B` (m).
    width_m: f64,
    /// Founding depth `Df` below grade (m).
    depth_m: f64,
    /// Global factor of safety applied to obtain the allowable pressure.
    factor_of_safety: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D footing solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for SoilBearingWorkbenchState {
    fn default() -> Self {
        // A 2 m strip footing founded 1 m deep on medium-dense sand with a
        // little cohesion (phi = 30 deg, c = 5 kPa, gamma = 18 kN/m^3) at a
        // factor of safety of 3: qult ~ 885 kPa, qall ~ 295 kPa.
        Self {
            friction_angle_deg: 30.0,
            cohesion_kpa: 5.0,
            unit_weight_kn_m3: 18.0,
            width_m: 2.0,
            depth_m: 1.0,
            factor_of_safety: 3.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Soil Bearing Workbench right-side panel. A no-op when the
/// `show_soilbearing_workbench` toggle is off.
pub fn draw_soilbearing_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_soilbearing_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_soilbearing_workbench",
        "Soil Bearing",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native Terzaghi strip-footing bearing capacity · valenx-soilbearing",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.soilbearing;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    // Associate each numeric `DragValue` with its caption via `labelled_by`, so
                    // the spin button carries the caption as its accessibility / UI-Automation
                    // Name (egui clears a DragValue's own Name otherwise).
                    ui.label(egui::RichText::new("Soil").strong());
                    ui.horizontal(|ui| {
                        let l = ui.label("friction angle φ (deg)");
                        ui.add(egui::DragValue::new(&mut s.friction_angle_deg).speed(0.5))
                            .labelled_by(l.id);
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("cohesion c (kPa)");
                        ui.add(egui::DragValue::new(&mut s.cohesion_kpa).speed(1.0))
                            .labelled_by(l.id);
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("unit weight γ (kN/m³)");
                        ui.add(egui::DragValue::new(&mut s.unit_weight_kn_m3).speed(0.5))
                            .labelled_by(l.id);
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Footing (continuous strip)").strong());
                    ui.horizontal(|ui| {
                        let l = ui.label("width B (m)");
                        ui.add(egui::DragValue::new(&mut s.width_m).speed(0.1))
                            .labelled_by(l.id);
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("founding depth Df (m)");
                        ui.add(egui::DragValue::new(&mut s.depth_m).speed(0.1))
                            .labelled_by(l.id);
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Safety").strong());
                    ui.horizontal(|ui| {
                        let l = ui.label("factor of safety");
                        ui.add(egui::DragValue::new(&mut s.factor_of_safety).speed(0.1))
                            .labelled_by(l.id);
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_soilbearing(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative continuous strip footing (a concrete slab atop a wider soil block) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Bearing capacity").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_soilbearing_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.soilbearing` borrow is
    // released here): build the footing's 3-D solid and load it.
    if app.soilbearing.show_3d_request {
        app.soilbearing.show_3d_request = false;
        load_footing_3d(app);
    }
}

/// Validate the form, evaluate the bearing capacity and format the readout.
fn run_soilbearing(s: &mut SoilBearingWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build the validated soil + footing value types from the form, mapping
/// any domain error to a display string. Shared by the readout and the
/// 3-D gate so an invalid configuration is rejected identically.
fn model(s: &SoilBearingWorkbenchState) -> Result<(SoilProperties, Footing), String> {
    let soil = SoilProperties::new(s.friction_angle_deg, s.cohesion_kpa, s.unit_weight_kn_m3)
        .map_err(|e| e.to_string())?;
    let footing = Footing::new(s.width_m, s.depth_m).map_err(|e| e.to_string())?;
    Ok((soil, footing))
}

/// Evaluate the bearing capacity and format the full readout, mapping any
/// domain error to a display string. Extracted so it is unit-testable.
fn compute(s: &SoilBearingWorkbenchState) -> Result<String, String> {
    let (soil, footing) = model(s)?;
    let r = bearing_capacity(&soil, &footing, s.factor_of_safety).map_err(|e| e.to_string())?;
    // Allowable line load Q_all = qall * B (force per unit run of footing).
    let q_all_line = r.q_allowable * footing.width();

    Ok(format!(
        "friction φ      : {:.1} deg\n\
         cohesion c      : {:.1} kPa\n\
         unit weight γ   : {:.1} kN/m³\n\
         footing B / Df  : {:.2} / {:.2} m\n\
         factor of safety: {:.1}\n\n\
         Nc / Nq / Nγ    : {:.2} / {:.2} / {:.2}\n\
         c·Nc term       : {:.2} kPa\n\
         q·Nq term       : {:.2} kPa\n\
         0.5γB·Nγ term   : {:.2} kPa\n\
         q ultimate      : {:.2} kPa\n\
         q allowable     : {:.2} kPa\n\
         Q allow (line)  : {:.2} kN/m",
        s.friction_angle_deg,
        s.cohesion_kpa,
        s.unit_weight_kn_m3,
        s.width_m,
        s.depth_m,
        s.factor_of_safety,
        r.factors.nc,
        r.factors.nq,
        r.factors.ngamma,
        r.cohesion_term,
        r.surcharge_term,
        r.self_weight_term,
        r.q_ultimate,
        r.q_allowable,
        q_all_line,
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

/// Build the footing-on-soil as a triangle [`Mesh`] — a representative
/// continuous strip footing modelled as a concrete slab (width `B` in `x`,
/// long in `y`) sitting at the founding depth atop a wider soil block.
/// Representative geometry (the bearing numbers are the
/// `valenx-soilbearing` result, not derived from this mesh). `None` for an
/// invalid configuration.
fn footing_solid_mesh(s: &SoilBearingWorkbenchState) -> Option<valenx_mesh::Mesh> {
    let (_soil, footing) = model(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Footing half-width in x scales with B (clamped to keep the figure
    // readable); the soil block is wider and deeper than the footing.
    let half_b = (footing.width() * 0.5).clamp(0.3, 2.5);
    let half_y = (half_b * 1.6).clamp(0.6, 3.0);
    let soil_half_x = half_b * 2.2;
    let soil_half_y = half_y * 1.3;
    let soil_depth = 1.4;
    let slab_thick = 0.25;

    // Soil block (top at z = 0, extending downward).
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, -soil_depth * 0.5),
        Vector3::new(soil_half_x, soil_half_y, soil_depth * 0.5),
    );
    // Concrete strip footing slab resting on the soil top.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, slab_thick * 0.5),
        Vector3::new(half_b, half_y, slab_thick * 0.5),
    );

    let mut block =
        valenx_mesh::element::ElementBlock::new(valenx_mesh::element::ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = valenx_mesh::Mesh::new("valenx-soilbearing");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D footing solid and load it into the central viewport.
fn load_footing_3d(app: &mut ValenxApp) {
    let Some(mesh) = footing_solid_mesh(&app.soilbearing) else {
        app.soilbearing.error =
            Some("footing parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<footing>/valenx-soilbearing"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// Agent-bridge product: the canonical soilbearing workbench as a 3-D solid plus its
/// `compute()` readout rows (see [`crate::products_registry`]).
pub(crate) fn soilbearing_product() -> crate::WorkspaceProduct {
    let s = SoilBearingWorkbenchState::default();
    let mesh = footing_solid_mesh(&s).expect("canonical soilbearing ⇒ footing solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<soilbearing>/valenx-footing");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical soilbearing ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Footing (Terzaghi bearing capacity)".into(),
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
        let s = SoilBearingWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_capacity_and_factors() {
        let mut s = SoilBearingWorkbenchState::default();
        run_soilbearing(&mut s);
        assert!(
            s.error.is_none(),
            "default footing should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("q ultimate"));
        assert!(s.result.contains("q allowable"));
        assert!(s.result.contains("Nc / Nq / Nγ"));
        // Textbook anchor (Das, Principles of Foundation Engineering) for
        // phi = 30 deg: Nc = 30.140 -> printed at {:.2} as "30.14".
        assert!(s.result.contains("30.14"));
    }

    #[test]
    fn analyze_rejects_zero_width() {
        let mut s = SoilBearingWorkbenchState {
            width_m: 0.0,
            ..Default::default()
        };
        run_soilbearing(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn analyze_rejects_factor_of_safety_at_or_below_one() {
        let mut s = SoilBearingWorkbenchState {
            factor_of_safety: 1.0,
            ..Default::default()
        };
        run_soilbearing(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn allowable_is_ultimate_over_fs_ground_truth() {
        // Ground truth: the allowable pressure is exactly the ultimate
        // pressure divided by the factor of safety. Recompute the result
        // independently of the formatted readout and check the identity.
        let s = SoilBearingWorkbenchState::default();
        let (soil, footing) = model(&s).expect("default model is valid");
        let r = bearing_capacity(&soil, &footing, s.factor_of_safety).unwrap();
        let expected_qall = r.q_ultimate / s.factor_of_safety;
        assert!(
            (r.q_allowable - expected_qall).abs() < 1e-9_f64,
            "qall {} != qult/FS {}",
            r.q_allowable,
            expected_qall
        );
    }

    #[test]
    fn factors_at_phi_30_match_textbook() {
        // Independent ground-truth anchor for the bearing-capacity factors
        // at phi = 30 deg (Das): Nq = 18.401, Nc = 30.140, Ngamma = 22.402.
        let s = SoilBearingWorkbenchState::default();
        let (soil, footing) = model(&s).expect("default model is valid");
        let r = bearing_capacity(&soil, &footing, s.factor_of_safety).unwrap();
        assert!((r.factors.nq - 18.401).abs() < 1e-2);
        assert!((r.factors.nc - 30.140).abs() < 1e-2);
        assert!((r.factors.ngamma - 22.402).abs() < 1e-2);
    }

    #[test]
    fn footing_mesh_for_default_is_nonempty_and_in_range() {
        let s = SoilBearingWorkbenchState::default();
        let mesh = footing_solid_mesh(&s).expect("default footing yields a solid");
        assert!(mesh.nodes.len() > 8, "expected soil block + footing slab");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn footing_mesh_none_for_invalid() {
        let s = SoilBearingWorkbenchState {
            width_m: 0.0,
            ..Default::default()
        };
        assert!(footing_solid_mesh(&s).is_none());
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
            draw_soilbearing_workbench(app, ctx);
        });
    }

    /// As [`draw_workbench`], but with accesskit enabled, returning the emitted
    /// accessibility tree nodes — the same tree a screen reader / AI driver
    /// consumes. `accesskit` is re-exported by egui, so no extra dependency.
    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_soilbearing_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_soilbearing_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_soilbearing_workbench = true;
        run_soilbearing(&mut app.soilbearing);
        draw_workbench(&mut app);
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // Every numeric `DragValue` is a SpinButton that must be `labelled_by`
        // its caption (egui clears a DragValue's own Name), so an AI / screen
        // reader can find the control by the caption text.
        let mut app = ValenxApp::default();
        app.show_soilbearing_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        assert!(
            spin_buttons.len() >= 6,
            "expected the numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every DragValue must be labelled_by its caption (AI-drivable name)"
        );

        for caption in ["friction angle φ (deg)", "width B (m)", "factor of safety"] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
        assert!(
            nodes.iter().any(|(_, n)| n.role() == Role::Button
                && n.name().is_some_and(|s| s.contains("Analyze"))),
            "the primary action button is a named, invokable node"
        );
    }
}
