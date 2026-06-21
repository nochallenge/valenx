//! The right-side **Fracture Mechanics Workbench** panel — native Mode-I
//! linear-elastic fracture analysis over `valenx-fracture`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_fracture_workbench`,
//! toggled from the View menu. The form sets a material (`K_Ic`, `σ_y`), a
//! crack geometry (central through-crack or edge crack), an applied remote
//! tension and a crack length; "Analyze" evaluates the four closed-form
//! LEFM relations — stress-intensity factor `K`, critical crack length
//! `a_c`, fracture stress `σ_f` and the Irwin plastic-zone radius — and
//! reports the `K` vs `K_Ic` fast-fracture verdict, and "Show 3-D" loads a
//! representative cracked-plate solid into the central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_fracture::{
    critical_crack_length, fracture_stress, plastic_zone_radius, stress_intensity, Material,
};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// The crack configuration, which fixes the dimensionless geometry factor
/// `Y` used everywhere in the Mode-I formulae.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CrackGeometry {
    /// Central through-crack in a wide plate, `Y = 1.0` (`a` is the
    /// half-length).
    CentralThrough,
    /// Single edge crack in a wide plate, `Y ≈ 1.12` (`a` is the full
    /// depth).
    EdgeCrack,
}

impl CrackGeometry {
    /// The dimensionless geometry factor `Y` for this configuration.
    fn geometry_factor(self) -> f64 {
        match self {
            CrackGeometry::CentralThrough => 1.0,
            CrackGeometry::EdgeCrack => 1.12,
        }
    }
}

/// Persistent form + result state for the Fracture Mechanics Workbench.
pub struct FractureWorkbenchState {
    /// Mode-I plane-strain fracture toughness `K_Ic` (MPa·√m).
    fracture_toughness: f64,
    /// Tensile yield strength `σ_y` (MPa).
    yield_strength: f64,
    /// Crack configuration, fixing the geometry factor `Y`.
    geometry: CrackGeometry,
    /// Remote applied tensile stress `σ` (MPa).
    applied_stress_mpa: f64,
    /// Crack length `a` entered in millimetres (converted to metres for the
    /// `MPa·√m` consistent unit set).
    crack_length_mm: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D cracked-plate solid (serviced after
    /// the panel draws).
    show_3d_request: bool,
}

impl Default for FractureWorkbenchState {
    fn default() -> Self {
        // 7075-T6 aluminium (K_Ic = 24 MPa*sqrt(m), sigma_y = 470 MPa) with
        // a 2 mm edge crack under 150 MPa of remote tension: K ~ 13.3, well
        // below the 24 toughness, so the flaw is stable (a_c ~ 6.5 mm).
        Self {
            fracture_toughness: 24.0,
            yield_strength: 470.0,
            geometry: CrackGeometry::EdgeCrack,
            applied_stress_mpa: 150.0,
            crack_length_mm: 2.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Fracture Mechanics Workbench right-side panel. A no-op when the
/// `show_fracture_workbench` toggle is off.
pub fn draw_fracture_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_fracture_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_fracture_workbench",
        "Fracture Mechanics",
        |app, ui| {
            ui.label(
                egui::RichText::new("native Mode-I linear-elastic fracture · valenx-fracture")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.fracture;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Material").strong());
                    ui.horizontal(|ui| {
                        ui.label("toughness K_Ic (MPa·√m)");
                        ui.add(egui::DragValue::new(&mut s.fracture_toughness).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("yield σ_y (MPa)");
                        ui.add(egui::DragValue::new(&mut s.yield_strength).speed(5.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Crack geometry").strong());
                    ui.radio_value(
                        &mut s.geometry,
                        CrackGeometry::CentralThrough,
                        "central through-crack (Y = 1.00)",
                    );
                    ui.radio_value(
                        &mut s.geometry,
                        CrackGeometry::EdgeCrack,
                        "edge crack (Y = 1.12)",
                    );

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Loading").strong());
                    ui.horizontal(|ui| {
                        ui.label("applied σ (MPa)");
                        ui.add(egui::DragValue::new(&mut s.applied_stress_mpa).speed(2.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("crack length a (mm)");
                        ui.add(egui::DragValue::new(&mut s.crack_length_mm).speed(0.1));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_fracture(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative cracked plate (a slab with an edge slot for the crack) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Fracture assessment").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_fracture_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.fracture` borrow is
    // released here): build the cracked-plate 3-D solid and load it.
    if app.fracture.show_3d_request {
        app.fracture.show_3d_request = false;
        load_specimen_3d(app);
    }
}

/// Validate the form, evaluate the LEFM relations and format the readout.
fn run_fracture(s: &mut FractureWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Evaluate the four closed-form Mode-I relations and format the full
/// readout, mapping any domain error to a display string. Extracted so it is
/// unit-testable.
fn compute(s: &FractureWorkbenchState) -> Result<String, String> {
    let material =
        Material::new(s.fracture_toughness, s.yield_strength).map_err(|e| e.to_string())?;
    let y = s.geometry.geometry_factor();
    let a_m = s.crack_length_mm / 1000.0;
    let sigma = s.applied_stress_mpa;

    let k = stress_intensity(y, sigma, a_m).map_err(|e| e.to_string())?;
    let a_c = critical_crack_length(&material, y, sigma).map_err(|e| e.to_string())?;
    let sigma_f = fracture_stress(&material, y, a_m).map_err(|e| e.to_string())?;
    let r_p = plastic_zone_radius(k, &material).map_err(|e| e.to_string())?;

    let k_ic = material.fracture_toughness;
    let verdict = if k >= k_ic {
        "FAST FRACTURE (K >= K_Ic)"
    } else {
        "stable (K < K_Ic)"
    };

    // Bind the derived / unit-converted quantities so every format
    // placeholder is an inlined identifier (clippy `uninlined_format_args`).
    let sy = s.yield_strength;
    let a_mm = s.crack_length_mm;
    let margin = k_ic / k;
    let a_c_mm = a_c * 1000.0;
    let r_p_mm = r_p * 1000.0;

    Ok(format!(
        "toughness K_Ic : {k_ic:.2} MPa·√m\n\
         yield σ_y      : {sy:.1} MPa\n\
         geometry Y     : {y:.2}\n\
         applied σ      : {sigma:.1} MPa\n\
         crack length a : {a_mm:.3} mm\n\n\
         stress int. K  : {k:.3} MPa·√m\n\
         margin K_Ic/K  : {margin:.2}\n\
         critical a_c   : {a_c_mm:.3} mm\n\
         fracture σ_f   : {sigma_f:.2} MPa\n\
         plastic zone r : {r_p_mm:.4} mm\n\
         verdict        : {verdict}"
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

/// Build the cracked specimen as a triangle [`Mesh`] — a plate slab with a
/// thin edge slot standing in for the crack, plus a base. Representative
/// geometry (not to scale; the LEFM numbers are the `valenx-fracture`
/// result). `None` for an invalid material configuration.
fn specimen_solid_mesh(s: &FractureWorkbenchState) -> Option<Mesh> {
    Material::new(s.fracture_toughness, s.yield_strength).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Plate split into two halves around the edge crack slot so the crack
    // reads as a gap entering from the -y edge.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, -0.42, 0.7),
        Vector3::new(0.5, 0.05, 0.6),
    );
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.27, 0.7),
        Vector3::new(0.5, 0.23, 0.6),
    );
    // The two ligaments either side of the through-thickness crack slot.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(-0.28, -0.05, 0.7),
        Vector3::new(0.22, 0.32, 0.6),
    );
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.28, -0.05, 0.7),
        Vector3::new(0.22, 0.32, 0.6),
    );
    // Base.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.05),
        Vector3::new(0.6, 0.6, 0.05),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-fracture");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D cracked specimen solid and load it into the central
/// viewport.
fn load_specimen_3d(app: &mut ValenxApp) {
    let Some(mesh) = specimen_solid_mesh(&app.fracture) else {
        app.fracture.error =
            Some("material parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<specimen>/valenx-fracture"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// Agent-bridge product: the canonical fracture workbench as a 3-D solid plus its
/// `compute()` readout rows (see [`crate::products_registry`]).
pub(crate) fn fracture_product() -> crate::WorkspaceProduct {
    let s = FractureWorkbenchState::default();
    let mesh = specimen_solid_mesh(&s).expect("canonical fracture ⇒ specimen solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<fracture>/valenx-specimen");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical fracture ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Fracture (stress-intensity K)".into(),
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
    use core::f64::consts::PI;

    #[test]
    fn default_state_is_idle() {
        let s = FractureWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_k_and_verdict() {
        let mut s = FractureWorkbenchState::default();
        run_fracture(&mut s);
        assert!(
            s.error.is_none(),
            "default specimen should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("stress int. K"));
        assert!(s.result.contains("critical a_c"));
        assert!(s.result.contains("plastic zone r"));
        // K = 1.12 * 150 * sqrt(pi * 0.002) ~ 13.317 MPa*sqrt(m), below the
        // 24 toughness, so the default flaw is stable.
        assert!(s.result.contains("13.317"));
        assert!(s.result.contains("stable (K < K_Ic)"));
    }

    #[test]
    fn analyze_rejects_zero_toughness() {
        let mut s = FractureWorkbenchState {
            fracture_toughness: 0.0,
            ..Default::default()
        };
        run_fracture(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn k_equals_ground_truth_closed_form() {
        // Ground truth: K = Y * sigma * sqrt(pi * a), hand-computed for the
        // default edge crack (Y = 1.12, sigma = 150 MPa, a = 0.002 m).
        let y = 1.12_f64;
        let sigma = 150.0_f64;
        let a = 0.002_f64;
        let expected = y * sigma * (PI * a).sqrt();
        let k = stress_intensity(y, sigma, a).unwrap();
        assert!((k - expected).abs() < 1e-12, "k = {k}, expected {expected}");
        // And the hand value itself.
        assert!((k - 13.31677972_f64).abs() < 1e-6, "k = {k}");
    }

    #[test]
    fn specimen_mesh_for_default_is_nonempty_and_in_range() {
        let s = FractureWorkbenchState::default();
        let mesh = specimen_solid_mesh(&s).expect("default specimen yields a solid");
        assert!(mesh.nodes.len() > 8, "expected plate + ligaments + base");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn specimen_mesh_none_for_invalid() {
        let s = FractureWorkbenchState {
            yield_strength: 0.0,
            ..Default::default()
        };
        assert!(specimen_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_fracture_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_fracture_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_fracture_workbench = true;
        run_fracture(&mut app.fracture);
        draw_workbench(&mut app);
    }
}
