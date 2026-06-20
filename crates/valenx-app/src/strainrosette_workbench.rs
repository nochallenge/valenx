//! The right-side **Strain Rosette Workbench** panel — native rectangular
//! 0/45/90 strain-gauge reduction over `valenx-strainrosette`.
//!
//! Mirrors the Heat Transfer / Buckling workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_strainrosette_workbench`,
//! toggled from the View menu. The form sets the three normal-strain gauge
//! readings of a rectangular (45-degree) rosette and an isotropic
//! linear-elastic material (Young's modulus, Poisson's ratio); "Analyze"
//! inverts the strain-transformation law to the Cartesian strain
//! `(eps_x, eps_y, gamma_xy)`, solves the principal strains `eps_1, eps_2`
//! and the major-axis angle, and closes with plane-stress Hooke's law for
//! `sigma_x, sigma_y, tau_xy`; "Show 3-D" loads a representative plate with
//! three angled gauge strips into the central viewport.

use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_strainrosette::{analyze, ElasticMaterial, RosetteReadings};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// How the strain magnitudes are rendered in the readout.
///
/// This is a presentation toggle only — the analysis always runs on the
/// dimensionless engineering strains; this just scales the displayed
/// numbers. (The crate models the rectangular rosette only, so there is
/// no physics mode to switch.)
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum StrainUnit {
    /// Dimensionless engineering strain (for example `0.0008`).
    Strain,
    /// Microstrain — strain times one million (for example `800`).
    Microstrain,
}

/// Persistent form + result state for the Strain Rosette Workbench.
pub struct StrainRosetteWorkbenchState {
    /// Gauge bonded along the x-axis (0 degrees), engineering strain.
    eps_0: f64,
    /// Gauge bonded at 45 degrees from the x-axis, engineering strain.
    eps_45: f64,
    /// Gauge bonded along the y-axis (90 degrees), engineering strain.
    eps_90: f64,
    /// Young's modulus `E` (MPa).
    youngs_modulus_mpa: f64,
    /// Poisson's ratio `nu` (dimensionless).
    poisson_ratio: f64,
    /// How strain magnitudes are displayed in the readout.
    unit: StrainUnit,
    /// Formatted analysis readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D plate solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for StrainRosetteWorkbenchState {
    fn default() -> Self {
        // A rectangular rosette on a steel surface (E = 200 GPa = 200000
        // MPa, nu = 0.3) reading 800 / 400 / 200 microstrain. Then
        // eps_x = 0.0008, eps_y = 0.0002, gamma_xy = 2*0.0004 - 0.0008 -
        // 0.0002 = -0.0002, giving principal strains ~ 8.16e-4 and
        // ~1.84e-4.
        Self {
            eps_0: 0.0008,
            eps_45: 0.0004,
            eps_90: 0.0002,
            youngs_modulus_mpa: 200_000.0,
            poisson_ratio: 0.3,
            unit: StrainUnit::Microstrain,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Strain Rosette Workbench right-side panel. A no-op when the
/// `show_strainrosette_workbench` toggle is off.
pub fn draw_strainrosette_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_strainrosette_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_strainrosette_workbench",
        "Strain Rosette",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native rectangular 0/45/90 rosette reduction · valenx-strainrosette",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.strainrosette;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Gauge readings (engineering strain)").strong());
                    ui.horizontal(|ui| {
                        ui.label("ε₀ (0°)");
                        ui.add(egui::DragValue::new(&mut s.eps_0).speed(1.0e-5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("ε₄₅ (45°)");
                        ui.add(egui::DragValue::new(&mut s.eps_45).speed(1.0e-5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("ε₉₀ (90°)");
                        ui.add(egui::DragValue::new(&mut s.eps_90).speed(1.0e-5));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Material").strong());
                    ui.horizontal(|ui| {
                        ui.label("Young's modulus E (MPa)");
                        ui.add(
                            egui::DragValue::new(&mut s.youngs_modulus_mpa)
                                .speed(1_000.0)
                                .range(0.0..=f64::MAX),
                        );
                    });
                    ui.horizontal(|ui| {
                        ui.label("Poisson's ratio ν");
                        ui.add(egui::DragValue::new(&mut s.poisson_ratio).speed(0.01));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Display units").strong());
                    ui.radio_value(&mut s.unit, StrainUnit::Strain, "strain (dimensionless)");
                    ui.radio_value(&mut s.unit, StrainUnit::Microstrain, "microstrain (×10⁻⁶)");

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_strainrosette(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative instrumented plate with three angled gauge strips (0 / 45 / 90 degrees) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Reduction").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_strainrosette_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.strainrosette` borrow is
    // released here): build the plate's 3-D solid and load it.
    if app.strainrosette.show_3d_request {
        app.strainrosette.show_3d_request = false;
        load_plate_3d(app);
    }
}

/// Validate the form, run the rosette reduction and format the readout.
fn run_strainrosette(s: &mut StrainRosetteWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Reduce the gauges, solve the principal strains and map to plane stress,
/// formatting the full readout. The only fallible step is building the
/// [`ElasticMaterial`]; its domain error is mapped to a display string.
/// Extracted so it is unit-testable.
fn compute(s: &StrainRosetteWorkbenchState) -> Result<String, String> {
    let material =
        ElasticMaterial::new(s.youngs_modulus_mpa, s.poisson_ratio).map_err(|e| e.to_string())?;
    let a = analyze(RosetteReadings::new(s.eps_0, s.eps_45, s.eps_90), &material);

    // Strain display scale + suffix.
    let (scale, strain_unit) = match s.unit {
        StrainUnit::Strain => (1.0, ""),
        StrainUnit::Microstrain => (1.0e6, " µε"),
    };
    let ex = a.strain.eps_x * scale;
    let ey = a.strain.eps_y * scale;
    let gxy = a.strain.gamma_xy * scale;
    let e1 = a.principal.eps_1 * scale;
    let e2 = a.principal.eps_2 * scale;
    let gmax = a.principal.max_shear() * scale;
    let theta_deg = a.principal.theta_p.to_degrees();

    Ok(format!(
        "gauges 0/45/90  : {e0:.6} / {e45:.6} / {e90:.6}\n\
         material E / ν   : {emod:.1} MPa / {nu:.3}\n\n\
         eps_x           : {ex:.3}{strain_unit}\n\
         eps_y           : {ey:.3}{strain_unit}\n\
         gamma_xy        : {gxy:.3}{strain_unit}\n\n\
         principal eps_1 : {e1:.3}{strain_unit}\n\
         principal eps_2 : {e2:.3}{strain_unit}\n\
         max shear γmax  : {gmax:.3}{strain_unit}\n\
         major axis θp   : {theta_deg:.2} °\n\n\
         sigma_x         : {sx:.2} MPa\n\
         sigma_y         : {sy:.2} MPa\n\
         tau_xy          : {txy:.2} MPa",
        e0 = s.eps_0,
        e45 = s.eps_45,
        e90 = s.eps_90,
        emod = s.youngs_modulus_mpa,
        nu = s.poisson_ratio,
        sx = a.stress.sigma_x,
        sy = a.stress.sigma_y,
        txy = a.stress.tau_xy,
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

/// Append a thin gauge strip (a flat box) centred on the plate face,
/// rotated by `theta` radians about the plate normal (z). Models a single
/// bonded foil gauge sitting on top of the plate.
fn push_gauge_strip(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    theta: f64,
    length: f64,
    width: f64,
) {
    // Half-extents of the strip in its own (un-rotated) frame: long in
    // local x, narrow in local y, very thin in z. Sit it just above the
    // plate's top face.
    let (hx, hy, hz) = (0.5 * length, 0.5 * width, 0.01);
    let z_top = 0.06;
    let base = nodes.len();
    let (sin_t, cos_t) = theta.sin_cos();
    let local = [
        (-hx, -hy, -hz),
        (hx, -hy, -hz),
        (hx, hy, -hz),
        (-hx, hy, -hz),
        (-hx, -hy, hz),
        (hx, -hy, hz),
        (hx, hy, hz),
        (-hx, hy, hz),
    ];
    for (lx, ly, lz) in local {
        // Rotate (lx, ly) about z, then lift to sit on the plate top.
        let rx = lx * cos_t - ly * sin_t;
        let ry = lx * sin_t + ly * cos_t;
        nodes.push(Vector3::new(rx, ry, z_top + lz));
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

/// Build the instrumented plate as a triangle [`Mesh`] — a flat plate with
/// three thin gauge strips bonded at 0, 45 and 90 degrees, the physical
/// layout of a rectangular rosette. Representative geometry (not to scale;
/// the strain / stress numbers are the `valenx-strainrosette` result).
/// `None` for an invalid material configuration.
fn plate_solid_mesh(s: &StrainRosetteWorkbenchState) -> Option<Mesh> {
    // Reject the same invalid material the analysis would reject.
    ElasticMaterial::new(s.youngs_modulus_mpa, s.poisson_ratio).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // The plate (thin in z).
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(0.6, 0.6, 0.05),
    );
    // Three bonded gauge strips at 0 / 45 / 90 degrees.
    let strip_len = 0.5;
    let strip_w = 0.08;
    push_gauge_strip(&mut nodes, &mut tris, 0.0, strip_len, strip_w);
    push_gauge_strip(&mut nodes, &mut tris, FRAC_PI_4, strip_len, strip_w);
    push_gauge_strip(&mut nodes, &mut tris, FRAC_PI_2, strip_len, strip_w);

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-strainrosette");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D instrumented-plate solid and load it into the central
/// viewport.
fn load_plate_3d(app: &mut ValenxApp) {
    let Some(mesh) = plate_solid_mesh(&app.strainrosette) else {
        app.strainrosette.error =
            Some("material parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<plate>/valenx-strainrosette"),
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
        let s = StrainRosetteWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_strain_principal_and_stress() {
        let mut s = StrainRosetteWorkbenchState::default();
        run_strainrosette(&mut s);
        assert!(
            s.error.is_none(),
            "default rosette should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("eps_x"));
        assert!(s.result.contains("principal eps_1"));
        assert!(s.result.contains("major axis θp"));
        assert!(s.result.contains("sigma_x"));
        // Default reads in microstrain: eps_x = 0.0008 -> 800 µε.
        assert!(s.result.contains("800.000 µε"), "result was:\n{}", s.result);
    }

    #[test]
    fn analyze_strain_units_switch_changes_readout() {
        let mut s = StrainRosetteWorkbenchState {
            unit: StrainUnit::Strain,
            ..Default::default()
        };
        run_strainrosette(&mut s);
        assert!(s.error.is_none());
        // In dimensionless strain, eps_x = 0.0008 prints as 0.001 at the
        // chosen precision, with no microstrain suffix.
        assert!(!s.result.contains("µε"), "result was:\n{}", s.result);
        assert!(s.result.contains("0.001"), "result was:\n{}", s.result);
    }

    #[test]
    fn analyze_rejects_out_of_range_poisson() {
        let mut s = StrainRosetteWorkbenchState {
            poisson_ratio: 0.5,
            ..Default::default()
        };
        run_strainrosette(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn reduction_matches_hand_computed_ground_truth() {
        // Ground truth for the rectangular 0/45/90 rosette:
        //   eps_x    = eps_0
        //   eps_y    = eps_90
        //   gamma_xy = 2*eps_45 - eps_0 - eps_90
        // With 0.0008 / 0.0004 / 0.0002:
        //   eps_x = 0.0008, eps_y = 0.0002,
        //   gamma_xy = 2*0.0004 - 0.0008 - 0.0002 = -0.0002.
        let e0 = 0.0008_f64;
        let e45 = 0.0004_f64;
        let e90 = 0.0002_f64;
        let a = analyze(
            RosetteReadings::new(e0, e45, e90),
            &ElasticMaterial::new(200_000.0, 0.3).unwrap(),
        );
        assert!((a.strain.eps_x - e0).abs() < 1e-12);
        assert!((a.strain.eps_y - e90).abs() < 1e-12);
        let gamma = 2.0 * e45 - e0 - e90;
        assert!((gamma - (-0.0002)).abs() < 1e-12);
        assert!((a.strain.gamma_xy - gamma).abs() < 1e-12);
        // Principal strains from the hand-built Mohr circle:
        //   mean   = (0.0008 + 0.0002)/2 = 0.0005
        //   radius = sqrt(0.0003^2 + 0.0001^2) = sqrt(1e-7)
        let mean = 0.0005_f64;
        let radius = (0.0003_f64 * 0.0003 + 0.0001 * 0.0001).sqrt();
        assert!((a.principal.eps_1 - (mean + radius)).abs() < 1e-12);
        assert!((a.principal.eps_2 - (mean - radius)).abs() < 1e-12);
    }

    #[test]
    fn plate_mesh_for_default_is_nonempty_and_in_range() {
        let s = StrainRosetteWorkbenchState::default();
        let mesh = plate_solid_mesh(&s).expect("default plate yields a solid");
        // Plate box (8) plus three gauge strips (8 each).
        assert!(mesh.nodes.len() > 8, "expected plate + three strips");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn plate_mesh_none_for_invalid_material() {
        let s = StrainRosetteWorkbenchState {
            poisson_ratio: 0.5,
            ..Default::default()
        };
        assert!(plate_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_strainrosette_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_strainrosette_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_strainrosette_workbench = true;
        run_strainrosette(&mut app.strainrosette);
        draw_workbench(&mut app);
    }
}
