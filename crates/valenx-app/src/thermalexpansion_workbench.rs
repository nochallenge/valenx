//! The right-side **Thermal Expansion Workbench** panel — native linear
//! thermal expansion + constrained thermal stress over
//! `valenx-thermalexpansion`.
//!
//! Mirrors the Beam / DC Motor workbenches: a resizable [`egui::SidePanel`]
//! gated on `crate::ValenxApp::show_thermalexpansion_workbench`, toggled from
//! the View menu. The form takes a linear expansion coefficient, a length, a
//! temperature change and a Young's modulus; "Analyze" reports the linear /
//! areal / volumetric coefficients, the length change and final length, the
//! free thermal strain and the fully-constrained thermal stress `E*alpha*dT`,
//! and "Show 3-D bar" loads a rod-between-fixed-walls solid into the central
//! viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_thermalexpansion::{
    constrained_thermal_stress, free_thermal_strain, linear_expansion, linear_final_length,
    LinearCoefficient, YoungsModulus,
};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Thermal Expansion Workbench.
pub struct ThermalExpansionWorkbenchState {
    /// Linear expansion coefficient `alpha` (per kelvin).
    alpha_per_k: f64,
    /// Original length `L0` (m).
    length_m: f64,
    /// Temperature change `dT` (K).
    delta_t_k: f64,
    /// Young's modulus `E` (GPa) for the constrained stress.
    youngs_gpa: f64,
    /// Formatted readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D bar solid (serviced after the panel
    /// draws).
    show_3d_request: bool,
}

impl Default for ThermalExpansionWorkbenchState {
    fn default() -> Self {
        // A 1 m aluminium bar heated 100 K (alpha 23e-6 /K, E 69 GPa).
        Self {
            alpha_per_k: 23.0e-6,
            length_m: 1.0,
            delta_t_k: 100.0,
            youngs_gpa: 69.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Thermal Expansion Workbench right-side panel. A no-op when the
/// `show_thermalexpansion_workbench` toggle is off.
pub fn draw_thermalexpansion_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_thermalexpansion_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_thermalexpansion_workbench",
        "Thermal Expansion",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native linear expansion + constrained stress · valenx-thermalexpansion",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.thermalexpansion;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Material + part").strong());
                    ui.horizontal(|ui| {
                        let cap = ui.label("linear α (/K)");
                        ui.add(egui::DragValue::new(&mut s.alpha_per_k).speed(1.0e-6))
                            .labelled_by(cap.id);
                    });
                    ui.horizontal(|ui| {
                        let cap = ui.label("length L0 (m)");
                        ui.add(egui::DragValue::new(&mut s.length_m).speed(0.05))
                            .labelled_by(cap.id);
                    });
                    ui.horizontal(|ui| {
                        let cap = ui.label("Young's E (GPa)");
                        ui.add(egui::DragValue::new(&mut s.youngs_gpa).speed(1.0))
                            .labelled_by(cap.id);
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Temperature").strong());
                    ui.horizontal(|ui| {
                        let cap = ui.label("ΔT (K)");
                        ui.add(egui::DragValue::new(&mut s.delta_t_k).speed(1.0))
                            .labelled_by(cap.id);
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_expansion(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D bar").strong())
                        .on_hover_text(
                            "Build a rod between two fixed walls as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Expansion + stress").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_thermalexpansion_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.thermalexpansion` borrow
    // is released here): build the bar's 3-D solid and load it.
    if app.thermalexpansion.show_3d_request {
        app.thermalexpansion.show_3d_request = false;
        load_bar_3d(app);
    }
}

/// Validate the form, compute the expansion + stress and format the readout.
fn run_expansion(s: &mut ThermalExpansionWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Compute the full readout, mapping any domain error to a display string.
/// Extracted so it is unit-testable.
fn compute(s: &ThermalExpansionWorkbenchState) -> Result<String, String> {
    let alpha = LinearCoefficient::new(s.alpha_per_k).map_err(|e| e.to_string())?;
    let dl = linear_expansion(alpha, s.length_m, s.delta_t_k).map_err(|e| e.to_string())?;
    let lf = linear_final_length(alpha, s.length_m, s.delta_t_k).map_err(|e| e.to_string())?;
    let strain = free_thermal_strain(alpha, s.delta_t_k).map_err(|e| e.to_string())?;
    let ym = YoungsModulus::from_gpa(s.youngs_gpa).map_err(|e| e.to_string())?;
    let stress = constrained_thermal_stress(ym, alpha, s.delta_t_k).map_err(|e| e.to_string())?;
    Ok(format!(
        "linear α     : {:.3e} /K\n\
         areal 2α     : {:.3e} /K\n\
         volumetric 3α: {:.3e} /K\n\
         length L0    : {:.4} m\n\
         ΔT           : {:.1} K\n\
         E            : {:.0} GPa\n\n\
         ΔL (linear)  : {:.4} mm\n\
         final length : {:.5} m\n\
         free strain  : {:.4} %\n\
         constrained σ: {:.1} MPa  (E·α·ΔT, both ends fixed)",
        alpha.per_kelvin(),
        alpha.areal(),
        alpha.volumetric(),
        s.length_m,
        s.delta_t_k,
        s.youngs_gpa,
        dl * 1000.0,
        lf,
        strain * 100.0,
        stress / 1.0e6,
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

/// Build the part as a triangle [`Mesh`] — a prismatic rod (length along x)
/// between two fixed end walls (the constrained case the stress assumes).
/// `None` for an invalid coefficient / length.
fn bar_solid_mesh(s: &ThermalExpansionWorkbenchState) -> Option<Mesh> {
    LinearCoefficient::new(s.alpha_per_k).ok()?;
    let l = s.length_m;
    if !(l.is_finite() && l > 0.0) {
        return None;
    }
    let hh = (l * 0.03).max(1e-3);
    let ww = (l * 0.04).max(1e-3);
    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();
    // Rod.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::zeros(),
        Vector3::new(l / 2.0, hh, hh),
    );
    // Fixed end walls.
    for sx in [-(l / 2.0 + ww / 2.0), l / 2.0 + ww / 2.0] {
        push_box(
            &mut nodes,
            &mut tris,
            Vector3::new(sx, 0.0, 0.0),
            Vector3::new(ww / 2.0, hh * 2.2, hh * 2.2),
        );
    }
    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-thermalexpansion");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D bar solid and load it into the central viewport.
fn load_bar_3d(app: &mut ValenxApp) {
    let Some(mesh) = bar_solid_mesh(&app.thermalexpansion) else {
        app.thermalexpansion.error =
            Some("bar parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<bar>/valenx-thermalexpansion"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// Agent-bridge product: the canonical thermal-expansion workbench as a 3-D
/// solid plus its `compute()` readout rows (see [`crate::products_registry`]).
pub(crate) fn thermalexpansion_product() -> crate::WorkspaceProduct {
    let s = ThermalExpansionWorkbenchState::default();
    let mesh = bar_solid_mesh(&s).expect("canonical thermal expansion ⇒ bar solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<thermalexpansion>/valenx-bar");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical thermal expansion ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Thermal expansion (ΔL/stress)".into(),
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
        let s = ThermalExpansionWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_expansion_and_stress() {
        let mut s = ThermalExpansionWorkbenchState::default();
        run_expansion(&mut s);
        assert!(s.error.is_none(), "default bar analyzes: {:?}", s.error);
        assert!(s.result.contains("ΔL"));
        assert!(s.result.contains("free strain"));
        assert!(s.result.contains("constrained"));
    }

    #[test]
    fn analyze_rejects_nonpositive_length() {
        let mut s = ThermalExpansionWorkbenchState {
            length_m: 0.0,
            ..Default::default()
        };
        run_expansion(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn bar_mesh_for_default_is_nonempty_and_in_range() {
        let s = ThermalExpansionWorkbenchState::default();
        let mesh = bar_solid_mesh(&s).expect("default bar yields a solid");
        assert!(mesh.nodes.len() > 8, "expected rod + two walls");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn bar_mesh_none_for_invalid() {
        let s = ThermalExpansionWorkbenchState {
            length_m: 0.0,
            ..Default::default()
        };
        assert!(bar_solid_mesh(&s).is_none());
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
            draw_thermalexpansion_workbench(app, ctx);
        });
    }

    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_thermalexpansion_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_thermalexpansion_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_thermalexpansion_workbench = true;
        run_expansion(&mut app.thermalexpansion);
        draw_workbench(&mut app);
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        let mut app = ValenxApp::default();
        app.show_thermalexpansion_workbench = true;
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
        for caption in ["linear α (/K)", "ΔT (K)"] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
    }
}
