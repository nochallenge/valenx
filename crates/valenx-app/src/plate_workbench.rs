//! The right-side **Plate Bending Workbench** panel — native thin
//! circular-plate bending over `valenx-plate`.
//!
//! Mirrors the Heat Transfer / Bearing workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_plate_workbench`,
//! toggled from the View menu. The form sets a uniformly-loaded circular
//! plate (radius, thickness, Young's modulus, Poisson's ratio, pressure)
//! with a clamped or simply-supported rim; "Analyze" evaluates the
//! closed-form Kirchhoff-Love small-deflection results — flexural rigidity
//! `D = E t^3 / (12 (1 - nu^2))`, centre deflection `w = k p a^4 / D`, and
//! the maximum extreme-fibre bending stress — and "Show 3-D" loads a
//! representative circular-plate disc (radius `a`, thickness `t`) into the
//! central viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_plate::{CircularPlate, EdgeSupport, PlateMaterial};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Plate Bending Workbench.
pub struct PlateWorkbenchState {
    /// Plate radius `a` (m).
    radius_m: f64,
    /// Plate thickness `t` (m).
    thickness_m: f64,
    /// Young's modulus `E` (Pa).
    youngs_modulus_pa: f64,
    /// Poisson's ratio `nu` (dimensionless, in the open interval `(-1, 0.5)`).
    poisson_ratio: f64,
    /// Uniform transverse pressure `p` over the face (Pa).
    pressure_pa: f64,
    /// Rim support condition (clamped / simply-supported).
    support: EdgeSupport,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D plate disc (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for PlateWorkbenchState {
    fn default() -> Self {
        // The crate's own steel-disc example: E = 200 GPa, nu = 0.3,
        // t = 5 mm, radius 250 mm (a/t = 50), uniform 20 kPa, clamped rim.
        // D ~ 2289 N.m, centre deflection ~ 0.533 mm, peak stress 37.5 MPa.
        Self {
            radius_m: 0.25,
            thickness_m: 0.005,
            youngs_modulus_pa: 200.0e9,
            poisson_ratio: 0.3,
            pressure_pa: 20.0e3,
            support: EdgeSupport::Clamped,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Plate Bending Workbench right-side panel. A no-op when the
/// `show_plate_workbench` toggle is off.
pub fn draw_plate_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_plate_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_plate_workbench",
        "Plate Bending",
        |app, ui| {
            ui.label(
                egui::RichText::new("native thin circular-plate bending · valenx-plate")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.plate;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Geometry").strong());
                    ui.horizontal(|ui| {
                        ui.label("radius a (m)");
                        ui.add(egui::DragValue::new(&mut s.radius_m).speed(0.01));
                    });
                    ui.horizontal(|ui| {
                        ui.label("thickness t (m)");
                        ui.add(egui::DragValue::new(&mut s.thickness_m).speed(0.001));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Material").strong());
                    ui.horizontal(|ui| {
                        ui.label("Young's E (Pa)");
                        ui.add(
                            egui::DragValue::new(&mut s.youngs_modulus_pa)
                                .speed(1.0e9)
                                .max_decimals(0),
                        );
                    });
                    ui.horizontal(|ui| {
                        ui.label("Poisson's ν");
                        ui.add(egui::DragValue::new(&mut s.poisson_ratio).speed(0.01));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Load").strong());
                    ui.horizontal(|ui| {
                        ui.label("pressure p (Pa)");
                        ui.add(egui::DragValue::new(&mut s.pressure_pa).speed(100.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Edge support").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.support, EdgeSupport::Clamped, "clamped");
                        ui.radio_value(
                            &mut s.support,
                            EdgeSupport::SimplySupported,
                            "simply-supported",
                        );
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_plate(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build the circular plate (radius a, thickness t) as a 3-D disc and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Bending").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_plate_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.plate` borrow is
    // released here): build the plate's 3-D disc and load it.
    if app.plate.show_3d_request {
        app.plate.show_3d_request = false;
        load_plate_3d(app);
    }
}

/// Validate the form, evaluate the plate and format the readout.
fn run_plate(s: &mut PlateWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build the validated [`CircularPlate`] from the form, the input both the
/// readout and the 3-D gate need. Extracted so it is unit-testable and
/// shared, mapping any domain error to a display string.
fn build_plate(s: &PlateWorkbenchState) -> Result<CircularPlate, String> {
    let material = PlateMaterial::new(s.youngs_modulus_pa, s.poisson_ratio, s.thickness_m)
        .map_err(|e| e.to_string())?;
    CircularPlate::new(material, s.radius_m, s.pressure_pa, s.support).map_err(|e| e.to_string())
}

/// Evaluate the plate and format the full readout, mapping any domain error
/// to a display string. Extracted so it is unit-testable.
fn compute(s: &PlateWorkbenchState) -> Result<String, String> {
    let plate = build_plate(s)?;
    let support = match s.support {
        EdgeSupport::Clamped => "clamped",
        EdgeSupport::SimplySupported => "simply-supported",
    };
    let d = plate.flexural_rigidity();
    let k = plate.deflection_coefficient();
    let a_over_t = s.radius_m / s.thickness_m;
    let w_mm = plate.center_deflection() * 1.0e3;
    let sigma_mpa = plate.max_bending_stress() * 1.0e-6;

    Ok(format!(
        "radius / thickness : {:.4} m / {:.4} m\n\
         a / t ratio        : {a_over_t:.1}\n\
         E / ν              : {:.3e} Pa / {:.3}\n\
         pressure p         : {:.1} Pa\n\
         edge support       : {support}\n\n\
         flexural rigidity D: {d:.1} N·m\n\
         deflection coeff k : {k:.6}\n\
         centre deflection w: {w_mm:.4} mm\n\
         max bending stress : {sigma_mpa:.1} MPa",
        s.radius_m, s.thickness_m, s.youngs_modulus_pa, s.poisson_ratio, s.pressure_pa,
    ))
}

/// Append a (double-sided) z-axis disc — a short cylindrical rim spanning
/// `c.z ± h_t` at radius `r`, capped top and bottom by a triangle fan to
/// the centre. The plate lies in the `x`-`y` plane with thickness along
/// `z`.
fn push_z_disc(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    c: Vector3<f64>,
    r: f64,
    h_t: f64,
    seg: usize,
) {
    let (z0, z1) = (c.z - h_t, c.z + h_t);

    // Bottom and top rim rings.
    let lo = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(c.x + r * a.cos(), c.y + r * a.sin(), z0));
    }
    let hi = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(c.x + r * a.cos(), c.y + r * a.sin(), z1));
    }
    // Centre points for the two end-cap fans.
    let c_lo = nodes.len();
    nodes.push(Vector3::new(c.x, c.y, z0));
    let c_hi = nodes.len();
    nodes.push(Vector3::new(c.x, c.y, z1));

    for j in 0..seg {
        let jn = (j + 1) % seg;
        // Side wall (two triangles per segment, double-sided).
        tris.extend_from_slice(&[lo + j, hi + j, hi + jn, lo + j, hi + jn, lo + jn]);
        // Bottom cap fan.
        tris.extend_from_slice(&[c_lo, lo + jn, lo + j]);
        // Top cap fan.
        tris.extend_from_slice(&[c_hi, hi + j, hi + jn]);
    }
}

/// Build the circular plate as a triangle [`Mesh`] — a single thin disc of
/// radius `a` and thickness `t`. Representative geometry (true radius and
/// thickness; the bending numbers are the `valenx-plate` result). `None`
/// for an invalid configuration.
fn plate_disc_mesh(s: &PlateWorkbenchState) -> Option<Mesh> {
    build_plate(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    push_z_disc(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.0),
        s.radius_m,
        s.thickness_m * 0.5,
        48,
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-plate");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D plate disc and load it into the central viewport.
fn load_plate_3d(app: &mut ValenxApp) {
    let Some(mesh) = plate_disc_mesh(&app.plate) else {
        app.plate.error = Some("plate parameters are invalid — cannot build the 3-D disc".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<plate>/valenx-plate"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"plate"}`** product: the representative
/// circular-plate disc (true radius `a`, thickness `t`) built from the
/// canonical steel disc (250 mm radius, 5 mm thick, 200 GPa, 20 kPa, clamped
/// rim), paired with the Kirchhoff-Love bending readout rows (flexural
/// rigidity / centre deflection / max stress), at a fixed 3/4 camera.
/// Registered in [`crate::products_registry`]; the per-tool builder the
/// registry dispatches to. Pure — driven off [`PlateWorkbenchState::default`].
pub(crate) fn plate_product() -> crate::WorkspaceProduct {
    let s = PlateWorkbenchState::default();
    let mesh = plate_disc_mesh(&s).expect("canonical plate ⇒ circular-disc solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<plate>/valenx-plate");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical plate ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Circular plate (Kirchhoff-Love)".into(),
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
        let s = PlateWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_rigidity_deflection_and_stress() {
        let mut s = PlateWorkbenchState::default();
        run_plate(&mut s);
        assert!(
            s.error.is_none(),
            "default plate should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("flexural rigidity D"));
        assert!(s.result.contains("centre deflection w"));
        assert!(s.result.contains("max bending stress"));
        // Steel disc: D ~ 2289.4 N.m, k = 1/64, w ~ 0.533 mm, sigma 37.5 MPa.
        assert!(s.result.contains("2289.4"));
        assert!(s.result.contains("0.015625"));
        assert!(s.result.contains("0.5332"));
        assert!(s.result.contains("37.5"));
    }

    #[test]
    fn analyze_rejects_thick_plate() {
        // radius / thickness = 0.02 / 0.01 = 2 < 5 -> NotThin.
        let mut s = PlateWorkbenchState {
            radius_m: 0.02,
            thickness_m: 0.01,
            ..Default::default()
        };
        run_plate(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn analyze_rejects_bad_poisson_ratio() {
        // nu = 0.5 is the excluded incompressible upper bound.
        let mut s = PlateWorkbenchState {
            poisson_ratio: 0.5,
            ..Default::default()
        };
        run_plate(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn flexural_rigidity_matches_hand_computed_ground_truth() {
        // Ground truth: D = E t^3 / (12 (1 - nu^2)). With E = 200 GPa,
        // t = 5 mm, nu = 0.3: D = 2e11 * 1.25e-7 / (12 * 0.91) = 2289.377...
        let s = PlateWorkbenchState::default();
        let plate = build_plate(&s).expect("default plate builds");
        let e = 200.0e9;
        let t: f64 = 0.005;
        let nu = 0.3;
        let expected = e * t.powi(3) / (12.0 * (1.0 - nu * nu));
        assert!(
            (expected - 2289.377289377289).abs() < 1.0e-6,
            "hand value {expected}"
        );
        assert!((plate.flexural_rigidity() - expected).abs() < 1.0e-9);
    }

    #[test]
    fn plate_mesh_for_default_is_nonempty_and_in_range() {
        let s = PlateWorkbenchState::default();
        let mesh = plate_disc_mesh(&s).expect("default plate yields a disc");
        assert!(mesh.nodes.len() > 8, "expected a tessellated disc");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn plate_mesh_none_for_invalid() {
        let s = PlateWorkbenchState {
            radius_m: 0.02,
            thickness_m: 0.01,
            ..Default::default()
        };
        assert!(plate_disc_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_plate_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_plate_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_plate_workbench = true;
        run_plate(&mut app.plate);
        draw_workbench(&mut app);
    }
}
