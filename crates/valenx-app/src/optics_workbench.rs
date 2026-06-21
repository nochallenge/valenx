//! The right-side **Optics Workbench** panel — native geometric-optics
//! analysis over `valenx-optics`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_optics_workbench`,
//! toggled from the View menu. One of three modes is analyzed at a time:
//!
//! - **Lensmaker** — focal length and power of a thin lens from its
//!   material index and the two surface radii,
//!   `1/f = (n - 1)(1/R1 - 1/R2)`.
//! - **Refraction** — Snell's law at a planar interface, with the critical
//!   angle and total-internal-reflection check.
//! - **Thin-lens imaging** — the Gaussian thin-lens equation: image
//!   distance, magnification, and the real/virtual upright/inverted
//!   classification.
//!
//! "Analyze" formats the readout; "Show 3-D" loads a representative
//! biconvex-lens solid (a solid of revolution) into the central viewport.

use std::f64::consts::PI;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_optics::lensmaker::Lens;
use valenx_optics::refraction::{classify_ray, Interface, RayOutcome};
use valenx_optics::thin_lens::{ImageKind, Orientation, ThinLens};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Which geometric-optics calculation the workbench is set to run.
#[derive(Clone, Copy, Debug, PartialEq)]
enum OpticsMode {
    /// Focal length / power from surface geometry (the lensmaker's equation).
    Lensmaker,
    /// Snell refraction at a planar interface, with critical angle / TIR.
    Refraction,
    /// Thin-lens imaging: image distance, magnification, classification.
    ThinLensImaging,
}

/// Persistent form + result state for the Optics Workbench.
pub struct OpticsWorkbenchState {
    /// Which calculation to run on "Analyze".
    mode: OpticsMode,

    // -- Lensmaker inputs -------------------------------------------------
    /// Refractive index `n` of the lens material.
    lens_n: f64,
    /// Front-surface radius of curvature `R1` (m).
    lens_r1_m: f64,
    /// Back-surface radius of curvature `R2` (m).
    lens_r2_m: f64,

    // -- Refraction inputs ------------------------------------------------
    /// Incident-side refractive index `n1`.
    iface_n1: f64,
    /// Transmitted-side refractive index `n2`.
    iface_n2: f64,
    /// Incidence angle from the normal (degrees).
    incidence_deg: f64,

    // -- Thin-lens imaging inputs ----------------------------------------
    /// Lens focal length `f` (m). Positive converging, negative diverging.
    focal_length_m: f64,
    /// Object distance in front of the lens (m).
    object_distance_m: f64,

    /// Formatted readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D lens solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for OpticsWorkbenchState {
    fn default() -> Self {
        // A symmetric biconvex crown-glass lens (n = 1.5, R1 = +0.10 m,
        // R2 = -0.10 m) has f = 0.10 m. The refraction default is a
        // glass -> air ray below its ~41.8 deg critical angle, and the
        // imaging default places an object at 2f for a 1:1 inverted image.
        Self {
            mode: OpticsMode::Lensmaker,
            lens_n: 1.5,
            lens_r1_m: 0.10,
            lens_r2_m: -0.10,
            iface_n1: 1.5,
            iface_n2: 1.0,
            incidence_deg: 30.0,
            focal_length_m: 0.10,
            object_distance_m: 0.20,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Optics Workbench right-side panel. A no-op when the
/// `show_optics_workbench` toggle is off.
pub fn draw_optics_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_optics_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_optics_workbench",
        "Optics",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native geometric-optics: lensmaker, Snell, thin-lens · valenx-optics",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.optics;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Calculation").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.mode, OpticsMode::Lensmaker, "Lensmaker");
                        ui.radio_value(&mut s.mode, OpticsMode::Refraction, "Refraction");
                    });
                    ui.radio_value(
                        &mut s.mode,
                        OpticsMode::ThinLensImaging,
                        "Thin-lens imaging",
                    );

                    ui.add_space(6.0);
                    match s.mode {
                        OpticsMode::Lensmaker => {
                            ui.label(egui::RichText::new("Lens").strong());
                            ui.horizontal(|ui| {
                                ui.label("index n");
                                ui.add(egui::DragValue::new(&mut s.lens_n).speed(0.01));
                            });
                            ui.horizontal(|ui| {
                                ui.label("R1 (m)");
                                ui.add(egui::DragValue::new(&mut s.lens_r1_m).speed(0.005));
                            });
                            ui.horizontal(|ui| {
                                ui.label("R2 (m)");
                                ui.add(egui::DragValue::new(&mut s.lens_r2_m).speed(0.005));
                            });
                            ui.label(
                                egui::RichText::new(
                                    "biconvex: R1 > 0, R2 < 0; flat surface: very large radius",
                                )
                                .weak()
                                .small(),
                            );
                        }
                        OpticsMode::Refraction => {
                            ui.label(egui::RichText::new("Interface").strong());
                            ui.horizontal(|ui| {
                                ui.label("n1 (incident)");
                                ui.add(egui::DragValue::new(&mut s.iface_n1).speed(0.01));
                            });
                            ui.horizontal(|ui| {
                                ui.label("n2 (transmitted)");
                                ui.add(egui::DragValue::new(&mut s.iface_n2).speed(0.01));
                            });
                            ui.horizontal(|ui| {
                                ui.label("incidence (°)");
                                ui.add(
                                    egui::DragValue::new(&mut s.incidence_deg)
                                        .speed(0.5)
                                        .range(0.0..=90.0),
                                );
                            });
                        }
                        OpticsMode::ThinLensImaging => {
                            ui.label(egui::RichText::new("Thin lens").strong());
                            ui.horizontal(|ui| {
                                ui.label("focal length f (m)");
                                ui.add(egui::DragValue::new(&mut s.focal_length_m).speed(0.005));
                            });
                            ui.horizontal(|ui| {
                                ui.label("object distance (m)");
                                ui.add(
                                    egui::DragValue::new(&mut s.object_distance_m).speed(0.01),
                                );
                            });
                            ui.label(
                                egui::RichText::new(
                                    "f > 0 converging, f < 0 diverging; object distance > 0",
                                )
                                .weak()
                                .small(),
                            );
                        }
                    }

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_optics(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative biconvex lens (a solid of revolution) and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Result").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_optics_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.optics` borrow is
    // released here): build the lens's 3-D solid and load it.
    if app.optics.show_3d_request {
        app.optics.show_3d_request = false;
        load_lens_3d(app);
    }
}

/// Validate the form, run the selected calculation and format the readout.
fn run_optics(s: &mut OpticsWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Run the selected calculation and format the full readout, mapping any
/// domain error to a display string. Extracted so it is unit-testable.
fn compute(s: &OpticsWorkbenchState) -> Result<String, String> {
    match s.mode {
        OpticsMode::Lensmaker => compute_lensmaker(s),
        OpticsMode::Refraction => compute_refraction(s),
        OpticsMode::ThinLensImaging => compute_imaging(s),
    }
}

/// Lensmaker's equation: focal length and power from the surface geometry.
fn compute_lensmaker(s: &OpticsWorkbenchState) -> Result<String, String> {
    let lens = Lens::new(s.lens_n, s.lens_r1_m, s.lens_r2_m).map_err(|e| e.to_string())?;
    let f = lens.focal_length().map_err(|e| e.to_string())?;
    let power = lens.power().map_err(|e| e.to_string())?;
    let kind = if f > 0.0 { "converging" } else { "diverging" };

    Ok(format!(
        "index n         : {:.3}\n\
         R1 / R2         : {:.4} / {:.4} m\n\n\
         focal length f  : {:.4} m\n\
         power P = 1/f   : {:.3} dioptres\n\
         lens type       : {kind}",
        s.lens_n, s.lens_r1_m, s.lens_r2_m, f, power,
    ))
}

/// Snell refraction at a planar interface, with the critical angle / TIR.
fn compute_refraction(s: &OpticsWorkbenchState) -> Result<String, String> {
    let iface = Interface::new(s.iface_n1, s.iface_n2).map_err(|e| e.to_string())?;
    let critical = match iface.critical_angle_deg() {
        Some(c) => format!("{c:.3} °"),
        None => "n/a (n1 <= n2)".to_string(),
    };
    let outcome =
        classify_ray(s.iface_n1, s.iface_n2, s.incidence_deg).map_err(|e| e.to_string())?;
    let ray = match outcome {
        RayOutcome::Refracted { angle_deg } => {
            format!("refracted at  : {angle_deg:.3} ° from normal")
        }
        RayOutcome::TotallyReflected { critical_deg } => {
            format!("total internal reflection (>= {critical_deg:.3} °)")
        }
    };

    Ok(format!(
        "n1 / n2         : {:.3} / {:.3}\n\
         incidence       : {:.3} °\n\
         critical angle  : {critical}\n\n\
         {ray}",
        s.iface_n1, s.iface_n2, s.incidence_deg,
    ))
}

/// Thin-lens imaging: image distance, magnification, classification.
fn compute_imaging(s: &OpticsWorkbenchState) -> Result<String, String> {
    let lens = ThinLens::new(s.focal_length_m).map_err(|e| e.to_string())?;
    let img = lens.image(s.object_distance_m).map_err(|e| e.to_string())?;
    let kind = match img.kind {
        ImageKind::Real => "real",
        ImageKind::Virtual => "virtual",
    };
    let orientation = match img.orientation {
        Orientation::Upright => "upright",
        Orientation::Inverted => "inverted",
    };
    let size = if img.is_magnified() {
        "magnified"
    } else if img.is_reduced() {
        "reduced"
    } else {
        "same size"
    };

    Ok(format!(
        "focal length f  : {:.4} m\n\
         object distance : {:.4} m\n\n\
         image distance  : {:.4} m\n\
         magnification   : {:.4}\n\
         image           : {kind}, {orientation}, {size}",
        s.focal_length_m, s.object_distance_m, img.distance, img.magnification,
    ))
}

/// Append an outward-facing quad (as two triangles, `a-b-c`, `a-c-d`) to
/// the index buffer for nodes already pushed at the given indices.
fn push_quad(tris: &mut Vec<usize>, a: usize, b: usize, c: usize, d: usize) {
    tris.extend_from_slice(&[a, b, c, a, c, d]);
}

/// Build a representative **biconvex lens** as a triangle [`Mesh`]: a solid
/// of revolution formed by two spherical caps (front and back) sharing a
/// circular rim, faceted into rings around the optical axis. Representative
/// geometry (not to scale; the optical numbers are the `valenx-optics`
/// result). `None` for an invalid lens configuration.
fn lens_solid_mesh() -> Option<Mesh> {
    // Fixed representative proportions (metres): rim radius and the bulge
    // (sagitta) of each cap. Both positive, so the apex z-offsets below are
    // well defined.
    let rim_radius = 0.5_f64;
    let sagitta = 0.18_f64;

    let radial_segments = 48usize;
    let rings = 8usize; // rings from rim (j = 0) to apex (j = rings).

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Cap apices on the optical (z) axis: front at +sagitta, back at
    // -sagitta. Pushed first so their indices are known.
    let front_apex = nodes.len();
    nodes.push(Vector3::new(0.0, 0.0, sagitta));
    let back_apex = nodes.len();
    nodes.push(Vector3::new(0.0, 0.0, -sagitta));

    // Rim-to-apex rings for the FRONT cap. ring[j][i] is node index.
    // Radius shrinks from rim_radius (j=0) to ~0 (j=rings) as a quarter
    // cosine; z grows from 0 to +sagitta as a quarter sine — a smooth
    // dome. Symmetric for the back cap with z negated.
    let mut front_rings: Vec<Vec<usize>> = Vec::new();
    let mut back_rings: Vec<Vec<usize>> = Vec::new();
    for j in 0..rings {
        // Fraction from rim (0.0) toward apex (1.0), excluding the apex
        // itself (handled by the apex node).
        let t = j as f64 / rings as f64;
        let radius = rim_radius * (t * PI / 2.0).cos();
        let z = sagitta * (t * PI / 2.0).sin();
        let mut front_row = Vec::with_capacity(radial_segments);
        let mut back_row = Vec::with_capacity(radial_segments);
        for i in 0..radial_segments {
            let phi = (i as f64 / radial_segments as f64) * 2.0 * PI;
            let (sin_phi, cos_phi) = phi.sin_cos();
            let x = radius * cos_phi;
            let y = radius * sin_phi;
            front_row.push(nodes.len());
            nodes.push(Vector3::new(x, y, z));
            back_row.push(nodes.len());
            nodes.push(Vector3::new(x, y, -z));
        }
        front_rings.push(front_row);
        back_rings.push(back_row);
    }

    // Quad bands between successive front rings, and between back rings.
    for j in 0..rings - 1 {
        for i in 0..radial_segments {
            let i_next = (i + 1) % radial_segments;
            push_quad(
                &mut tris,
                front_rings[j][i],
                front_rings[j][i_next],
                front_rings[j + 1][i_next],
                front_rings[j + 1][i],
            );
            push_quad(
                &mut tris,
                back_rings[j][i],
                back_rings[j + 1][i],
                back_rings[j + 1][i_next],
                back_rings[j][i_next],
            );
        }
    }

    // Cap the innermost ring to each apex with a triangle fan.
    let inner = rings - 1;
    for i in 0..radial_segments {
        let i_next = (i + 1) % radial_segments;
        tris.extend_from_slice(&[
            front_rings[inner][i],
            front_rings[inner][i_next],
            front_apex,
        ]);
        tris.extend_from_slice(&[back_rings[inner][i_next], back_rings[inner][i], back_apex]);
    }

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-optics");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D lens solid and load it into the central viewport.
fn load_lens_3d(app: &mut ValenxApp) {
    let Some(mesh) = lens_solid_mesh() else {
        app.optics.error = Some("cannot build the 3-D lens solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<lens>/valenx-optics"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"optics"}`** product: the canonical
/// biconvex lens built as a 3-D solid, paired with the workbench's own
/// `compute()` readout rows, at a fixed 3/4 camera. Registered in
/// [`crate::products_registry`]; the per-tool builder the registry dispatches
/// to. Pure — driven off [`OpticsWorkbenchState::default`]. (The lens solid
/// is the canonical representative geometry — it takes no parameters.)
pub(crate) fn optics_product() -> crate::WorkspaceProduct {
    let s = OpticsWorkbenchState::default();
    let mesh = lens_solid_mesh().expect("canonical lens ⇒ solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<lens>/valenx-optics");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical optics ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Optics (lens)".into(),
        lines,
        mesh: Some(loaded),
        vertex_colors: None,
        camera,
        kind2d: None,
        last_export: None,
        image: None,
        image_texture: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_is_idle() {
        let s = OpticsWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
        assert_eq!(s.mode, OpticsMode::Lensmaker);
    }

    #[test]
    fn lensmaker_default_reports_focal_length_and_power() {
        let mut s = OpticsWorkbenchState::default();
        run_optics(&mut s);
        assert!(
            s.error.is_none(),
            "default lens should analyze: {:?}",
            s.error
        );
        // Symmetric biconvex n=1.5, R1=+0.10, R2=-0.10 -> f = 0.10 m,
        // P = 10 dioptres, converging.
        assert!(s.result.contains("focal length f  : 0.1000 m"));
        assert!(s.result.contains("10.000 dioptres"));
        assert!(s.result.contains("converging"));
    }

    #[test]
    fn lensmaker_ground_truth_lensmaker_equation() {
        // Hand-computed ground truth: 1/f = (n - 1)(1/R1 - 1/R2).
        // n=1.5, R1=+0.10, R2=-0.10:
        // 1/f = 0.5 * (10 - (-10)) = 10  ->  f = 0.10 m exactly.
        let n = 1.5_f64;
        let r1 = 0.10_f64;
        let r2 = -0.10_f64;
        let inv_f = (n - 1.0) * (1.0 / r1 - 1.0 / r2);
        let f_expected = 1.0 / inv_f;
        assert!((f_expected - 0.10).abs() < 1e-12, "f = {f_expected}");
        let lens = Lens::new(n, r1, r2).unwrap();
        assert!((lens.focal_length().unwrap() - f_expected).abs() < 1e-12);
    }

    #[test]
    fn refraction_default_refracts_below_critical() {
        let mut s = OpticsWorkbenchState {
            mode: OpticsMode::Refraction,
            ..Default::default()
        };
        run_optics(&mut s);
        assert!(s.error.is_none(), "glass->air at 30 deg: {:?}", s.error);
        // Glass (1.5) -> air (1.0): critical ~41.810 deg; 30 deg refracts.
        assert!(s.result.contains("critical angle  : 41.810 °"));
        assert!(s.result.contains("refracted at"));
    }

    #[test]
    fn refraction_reports_total_internal_reflection_past_critical() {
        let mut s = OpticsWorkbenchState {
            mode: OpticsMode::Refraction,
            incidence_deg: 60.0,
            ..Default::default()
        };
        run_optics(&mut s);
        // 60 deg is past the ~41.8 deg critical angle for glass -> air.
        assert!(
            s.error.is_none(),
            "TIR is an outcome, not an error: {:?}",
            s.error
        );
        assert!(s.result.contains("total internal reflection"));
    }

    #[test]
    fn imaging_default_object_at_two_f_is_real_inverted_same_size() {
        let mut s = OpticsWorkbenchState {
            mode: OpticsMode::ThinLensImaging,
            ..Default::default()
        };
        run_optics(&mut s);
        assert!(s.error.is_none(), "f=0.10, do=0.20: {:?}", s.error);
        // Object at 2f: image at 2f (0.20 m), m = -1, real inverted.
        assert!(s.result.contains("image distance  : 0.2000 m"));
        assert!(s.result.contains("magnification   : -1.0000"));
        assert!(s.result.contains("real, inverted, same size"));
    }

    #[test]
    fn imaging_rejects_nonpositive_object_distance() {
        let mut s = OpticsWorkbenchState {
            mode: OpticsMode::ThinLensImaging,
            object_distance_m: 0.0,
            ..Default::default()
        };
        run_optics(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn lens_mesh_is_nonempty_and_in_range() {
        let mesh = lens_solid_mesh().expect("biconvex lens yields a solid");
        assert!(mesh.nodes.len() > 8, "expected two faceted caps");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_optics_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_optics_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_optics_workbench = true;
        run_optics(&mut app.optics);
        draw_workbench(&mut app);
    }
}
