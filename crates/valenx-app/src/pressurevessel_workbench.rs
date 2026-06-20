//! The right-side **Pressure Vessel Workbench** panel — native closed-form
//! pressure-vessel stress analysis over `valenx-pressure-vessel`.
//!
//! Mirrors the Heat Transfer / Flywheel workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_pressurevessel_workbench`,
//! toggled from the View menu. The form sets a vessel (shape — cylinder or
//! sphere; wall model — thin membrane or the exact thick-wall Lamé
//! solution) under a gauge internal pressure; "Analyze" reports the
//! membrane / bore stresses (cylinder hoop `p r / t` and longitudinal
//! `p r / (2 t)`, sphere `p r / (2 t)`, or the thick-wall bore hoop and
//! radial), the maximum shear and the von Mises equivalent stress, plus
//! the thin-wall validity check; "Show 3-D" loads a representative
//! cylindrical vessel body (a wall tube with end caps) into the central
//! viewport to orbit.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;
use std::f64::consts::TAU;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_pressure_vessel::{ThickCylinder, ThinCylinder, ThinSphere};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// The vessel geometry chosen in the form.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum VesselShape {
    /// A closed cylindrical vessel — hoop and longitudinal stresses.
    Cylinder,
    /// A spherical vessel — one isotropic membrane stress.
    Sphere,
}

/// The wall model chosen in the form.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum WallModel {
    /// Thin-wall membrane theory (`r / t > 10`).
    ThinMembrane,
    /// The exact thick-wall Lamé elasticity solution (cylinder only).
    ThickLame,
}

/// Persistent form + result state for the Pressure Vessel Workbench.
pub struct PressureVesselWorkbenchState {
    /// Vessel shape (cylinder or sphere).
    shape: VesselShape,
    /// Wall model (thin membrane or thick-wall Lamé).
    wall: WallModel,
    /// Internal gauge pressure `p` (MPa).
    pressure_mpa: f64,
    /// Internal (bore) radius `r` (m).
    radius_m: f64,
    /// Wall thickness `t` (m) — used by the thin-wall membrane model and
    /// to derive the outer radius for the thick-wall model.
    thickness_m: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D vessel solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for PressureVesselWorkbenchState {
    fn default() -> Self {
        // A 1 m bore, 10 mm wall cylinder at 2 MPa: r/t = 100 (firmly
        // thin-wall), hoop = p*r/t = 200 MPa, longitudinal p*r/(2t) =
        // 100 MPa, von Mises = 200*sqrt(3)/2 ~ 173 MPa.
        Self {
            shape: VesselShape::Cylinder,
            wall: WallModel::ThinMembrane,
            pressure_mpa: 2.0,
            radius_m: 1.0,
            thickness_m: 0.01,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Pressure Vessel Workbench right-side panel. A no-op when the
/// `show_pressurevessel_workbench` toggle is off.
pub fn draw_pressurevessel_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_pressurevessel_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_pressurevessel_workbench",
        "Pressure Vessel",
        |app, ui| {
            ui.label(
                egui::RichText::new("native closed-form vessel stress · valenx-pressure-vessel")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.pressurevessel;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Shape").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.shape, VesselShape::Cylinder, "cylinder");
                        ui.radio_value(&mut s.shape, VesselShape::Sphere, "sphere");
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Wall model").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.wall, WallModel::ThinMembrane, "thin membrane");
                        ui.radio_value(&mut s.wall, WallModel::ThickLame, "thick (Lamé)");
                    });
                    if s.shape == VesselShape::Sphere && s.wall == WallModel::ThickLame {
                        ui.label(
                            egui::RichText::new("thick-wall model is cylinder-only — falls back to the membrane sphere")
                                .weak()
                                .small(),
                        );
                    }

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Loading & geometry").strong());
                    ui.horizontal(|ui| {
                        ui.label("pressure p (MPa)");
                        ui.add(egui::DragValue::new(&mut s.pressure_mpa).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("bore radius r (m)");
                        ui.add(egui::DragValue::new(&mut s.radius_m).speed(0.01));
                    });
                    ui.horizontal(|ui| {
                        ui.label("wall thickness t (m)");
                        ui.add(egui::DragValue::new(&mut s.thickness_m).speed(0.001));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_pressurevessel(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative cylindrical vessel body (a wall tube with end caps) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Stresses").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_pressurevessel_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.pressurevessel` borrow
    // is released here): build the vessel's 3-D solid and load it.
    if app.pressurevessel.show_3d_request {
        app.pressurevessel.show_3d_request = false;
        load_vessel_3d(app);
    }
}

/// Validate the form, evaluate the vessel and format the readout.
fn run_pressurevessel(s: &mut PressureVesselWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Whether the configuration evaluates to a usable thick-wall cylinder:
/// thick model selected, cylindrical shape. (Spheres and the thin model
/// route through the membrane formulae instead.)
fn uses_thick_cylinder(s: &PressureVesselWorkbenchState) -> bool {
    s.shape == VesselShape::Cylinder && s.wall == WallModel::ThickLame
}

/// Evaluate the vessel and format the full readout, mapping any domain
/// error to a display string. Extracted so it is unit-testable.
fn compute(s: &PressureVesselWorkbenchState) -> Result<String, String> {
    let p = s.pressure_mpa;
    let r = s.radius_m;
    let t = s.thickness_m;

    let header = format!(
        "shape / wall    : {} / {}\n\
         pressure p      : {:.3} MPa\n\
         bore radius r   : {:.4} m\n\
         wall thickness t: {:.4} m\n\
         r / t ratio     : {:.2}\n",
        match s.shape {
            VesselShape::Cylinder => "cylinder",
            VesselShape::Sphere => "sphere",
        },
        match s.wall {
            WallModel::ThinMembrane => "thin membrane",
            WallModel::ThickLame => "thick (Lamé)",
        },
        p,
        r,
        t,
        r / t,
    );

    match s.shape {
        VesselShape::Cylinder if s.wall == WallModel::ThickLame => {
            // Exact thick-wall Lamé solution; outer radius from the wall
            // thickness. Internal pressure only (gauge), closed ends.
            let c = ThickCylinder::new(r, r + t, p).map_err(|e| e.to_string())?;
            let bore = c.stress_at(r).map_err(|e| e.to_string())?;
            let hoop_out = c.hoop_stress(r + t).map_err(|e| e.to_string())?;
            let max_shear = c.max_shear_at(r).map_err(|e| e.to_string())?;
            let axial = c.axial_stress_closed();
            Ok(format!(
                "{header}\n\
                 -- thick-wall Lamé (bore = critical) --\n\
                 bore hoop  σθ(a): {:.2} MPa\n\
                 bore radial σr(a): {:.2} MPa\n\
                 outer hoop σθ(b): {:.2} MPa\n\
                 axial (closed)  : {:.2} MPa\n\
                 max shear @ bore: {:.2} MPa",
                bore.hoop, bore.radial, hoop_out, axial, max_shear,
            ))
        }
        VesselShape::Cylinder => {
            // Thin-wall membrane cylinder.
            let c = ThinCylinder::new(p, r, t).map_err(|e| e.to_string())?;
            let validity = if c.is_thin_wall() {
                "OK (r/t ≥ 10)"
            } else {
                "INVALID — r/t < 10, use thick-wall"
            };
            Ok(format!(
                "{header}\n\
                 -- thin-wall membrane cylinder --\n\
                 hoop      σh   : {:.2} MPa\n\
                 longitudinal σl: {:.2} MPa\n\
                 max shear (3-D): {:.2} MPa\n\
                 von Mises      : {:.2} MPa\n\
                 thin-wall check: {validity}",
                c.hoop_stress(),
                c.longitudinal_stress(),
                c.max_shear_stress(),
                c.von_mises_stress(),
            ))
        }
        VesselShape::Sphere => {
            // Thin-wall membrane sphere (the thick model has no sphere, so
            // a thick selection falls back to this with a note above).
            let sph = ThinSphere::new(p, r, t).map_err(|e| e.to_string())?;
            let validity = if sph.is_thin_wall() {
                "OK (r/t ≥ 10)"
            } else {
                "INVALID — r/t < 10"
            };
            Ok(format!(
                "{header}\n\
                 -- thin-wall membrane sphere --\n\
                 membrane  σ    : {:.2} MPa\n\
                 max shear      : {:.2} MPa\n\
                 von Mises      : {:.2} MPa\n\
                 thin-wall check: {validity}",
                sph.membrane_stress(),
                sph.max_shear_stress(),
                sph.von_mises_stress(),
            ))
        }
    }
}

/// Append a swept cylindrical wall tube — an annulus of outer radius
/// `r_out` and inner (bore) radius `r_in`, centred on the `z` axis with
/// axial half-length `hz` — to the buffers as a closed triangulated solid
/// (the two end-cap annuli plus the outer and inner bore walls). `seg` is
/// the number of angular segments.
fn push_cylinder(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    r_out: f64,
    r_in: f64,
    hz: f64,
    seg: usize,
) {
    let base = nodes.len();
    // Ring of nodes at each radius on the top (+z) and bottom (-z) cap.
    // Order per angular step i: outer-top, outer-bottom, inner-top,
    // inner-bottom.
    for i in 0..seg {
        let a = TAU * (i as f64) / (seg as f64);
        let (sin, cos) = a.sin_cos();
        nodes.push(Vector3::new(r_out * cos, r_out * sin, hz));
        nodes.push(Vector3::new(r_out * cos, r_out * sin, -hz));
        nodes.push(Vector3::new(r_in * cos, r_in * sin, hz));
        nodes.push(Vector3::new(r_in * cos, r_in * sin, -hz));
    }
    let mut quad = |a: usize, b: usize, c: usize, d: usize| {
        tris.extend_from_slice(&[base + a, base + b, base + c, base + a, base + c, base + d]);
    };
    for i in 0..seg {
        let j = (i + 1) % seg;
        let (ot, ob, it, ib) = (4 * i, 4 * i + 1, 4 * i + 2, 4 * i + 3);
        let (not_, nob, nit, nib) = (4 * j, 4 * j + 1, 4 * j + 2, 4 * j + 3);
        // Top end-cap annulus, bottom end-cap annulus, outer wall, bore.
        quad(it, ot, not_, nit);
        quad(ob, ib, nib, nob);
        quad(ot, ob, nob, not_);
        quad(ib, it, nit, nib);
    }
}

/// Build the vessel body as a triangle [`Mesh`] — a representative
/// cylindrical wall tube (bore radius `r`, wall thickness `t`) with end
/// caps, swept about the `z` axis. Representative geometry (the bore is
/// exaggerated to a slim shell when very thin; the stress numbers are the
/// `valenx-pressure-vessel` result). `None` for an invalid configuration.
fn vessel_solid_mesh(s: &PressureVesselWorkbenchState) -> Option<Mesh> {
    // Reuse the analysis validation: a configuration that cannot evaluate
    // has no solid. For the thick path the constructor also checks the
    // geometry; the thin/sphere paths gate on positive r, t.
    if uses_thick_cylinder(s) {
        ThickCylinder::new(s.radius_m, s.radius_m + s.thickness_m, s.pressure_mpa).ok()?;
    } else {
        ThinCylinder::new(s.pressure_mpa, s.radius_m, s.thickness_m).ok()?;
    }

    let r_in = s.radius_m;
    let r_out = s.radius_m + s.thickness_m;
    // Representative axial half-length: a few diameters tall, keyed to the
    // bore so the proportions read as a vessel rather than a disc.
    let hz = (1.5 * r_in).clamp(0.05, 5.0);

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();
    push_cylinder(&mut nodes, &mut tris, r_out, r_in, hz, 64);

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-pressure-vessel");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D vessel solid and load it into the central viewport.
fn load_vessel_3d(app: &mut ValenxApp) {
    let Some(mesh) = vessel_solid_mesh(&app.pressurevessel) else {
        app.pressurevessel.error =
            Some("vessel parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<vessel>/valenx-pressure-vessel"),
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
        let s = PressureVesselWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_cylinder_reports_hoop_and_von_mises() {
        let mut s = PressureVesselWorkbenchState::default();
        run_pressurevessel(&mut s);
        assert!(
            s.error.is_none(),
            "default cylinder should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("hoop"));
        assert!(s.result.contains("von Mises"));
        // p=2, r=1, t=0.01 -> hoop = 200 MPa, longitudinal = 100 MPa.
        assert!(s.result.contains("200.00"));
        assert!(s.result.contains("100.00"));
        // Firmly thin-wall (r/t = 100).
        assert!(s.result.contains("OK"));
    }

    #[test]
    fn analyze_sphere_reports_membrane_stress() {
        let mut s = PressureVesselWorkbenchState {
            shape: VesselShape::Sphere,
            ..Default::default()
        };
        run_pressurevessel(&mut s);
        assert!(s.error.is_none(), "sphere should analyze: {:?}", s.error);
        assert!(s.result.contains("membrane"));
        // Sphere membrane = p r / (2 t) = 2*1/(2*0.01) = 100 MPa.
        assert!(s.result.contains("100.00"));
    }

    #[test]
    fn analyze_thick_cylinder_reports_bore_hoop() {
        let mut s = PressureVesselWorkbenchState {
            wall: WallModel::ThickLame,
            ..Default::default()
        };
        run_pressurevessel(&mut s);
        assert!(
            s.error.is_none(),
            "thick cylinder should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("thick-wall Lamé"));
        assert!(s.result.contains("bore hoop"));
    }

    #[test]
    fn analyze_rejects_zero_thickness() {
        let mut s = PressureVesselWorkbenchState {
            thickness_m: 0.0,
            ..Default::default()
        };
        run_pressurevessel(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn membrane_ground_truth_hoop_and_longitudinal() {
        // Ground truth: for a thin cylinder the hoop stress is p r / t and
        // the longitudinal stress is p r / (2 t) — exactly half. Hand-
        // compute for p=2, r=1, t=0.01: hoop=200, longitudinal=100.
        let p = 2.0_f64;
        let r = 1.0_f64;
        let t = 0.01_f64;
        let c = ThinCylinder::new(p, r, t).unwrap();
        let hoop_hand = p * r / t;
        let long_hand = p * r / (2.0 * t);
        assert!((c.hoop_stress() - hoop_hand).abs() < 1e-9);
        assert!((c.longitudinal_stress() - long_hand).abs() < 1e-9);
        assert!((hoop_hand - 200.0).abs() < 1e-9);
        assert!((long_hand - 100.0).abs() < 1e-9);
        // And the defining relation: hoop is exactly twice longitudinal.
        assert!((c.hoop_stress() - 2.0 * c.longitudinal_stress()).abs() < 1e-9);
    }

    #[test]
    fn vessel_mesh_for_default_is_nonempty_and_in_range() {
        let s = PressureVesselWorkbenchState::default();
        let mesh = vessel_solid_mesh(&s).expect("default vessel yields a solid");
        assert!(mesh.nodes.len() > 8, "expected a swept wall tube");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn vessel_mesh_none_for_invalid() {
        let s = PressureVesselWorkbenchState {
            thickness_m: 0.0,
            ..Default::default()
        };
        assert!(vessel_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_pressurevessel_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_pressurevessel_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_pressurevessel_workbench = true;
        run_pressurevessel(&mut app.pressurevessel);
        draw_workbench(&mut app);
    }
}
