//! The right-side **Shaft Design Workbench** panel — native combined
//! bending-and-torsion stress analysis of a solid circular shaft over
//! `valenx-shaftdesign`.
//!
//! Mirrors the Heat Transfer / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_shaftdesign_workbench`,
//! toggled from the View menu. The form sets a solid round shaft
//! (diameter) under a steady bending moment `M` and torque `T` with a
//! material allowable shear stress; "Analyze" reports the bending and
//! torsional fibre stresses, the equivalent torque `T_e = sqrt(M^2 + T^2)`,
//! the combined Tresca maximum shear and principal normal stresses, the
//! section properties (section modulus `Z`, polar section modulus `Z_p`,
//! polar second moment `J`), and the minimum diameter that keeps the
//! combined shear within the allowable; "Show 3-D" loads a representative
//! round shaft solid into the central viewport.

use std::f64::consts::PI;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_shaftdesign::bending::bending_stress_section;
use valenx_shaftdesign::torsion::shear_stress_section;
use valenx_shaftdesign::{
    diameter_for_shear_stress, equivalent_torque, max_normal_stress, max_shear_stress, ShaftSection,
};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// The combined-load criterion used to size the shaft diameter.
#[derive(Clone, Copy, PartialEq, Eq)]
enum SizingCriterion {
    /// Maximum-shear-stress (Tresca) via the equivalent torque.
    Shear,
}

/// Persistent form + result state for the Shaft Design Workbench.
pub struct ShaftDesignWorkbenchState {
    /// Shaft outer diameter `d` (mm; converted to metres for the solver).
    diameter_mm: f64,
    /// Applied steady bending moment `M` at the section (N·m).
    bending_moment_nm: f64,
    /// Applied steady torque `T` at the section (N·m).
    torque_nm: f64,
    /// Material allowable shear stress (MPa; converted to pascals).
    allowable_shear_mpa: f64,
    /// Which combined-load criterion drives the diameter sizing.
    criterion: SizingCriterion,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D shaft solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for ShaftDesignWorkbenchState {
    fn default() -> Self {
        // A 50 mm steel line-shaft carrying M = 600 N·m and T = 800 N·m
        // (equivalent torque T_e = 1000 N·m), checked against a 40 MPa
        // allowable shear: combined Tresca shear ~40.7 MPa, so the
        // section needs ~50.3 mm — just over the modelled 50 mm.
        Self {
            diameter_mm: 50.0,
            bending_moment_nm: 600.0,
            torque_nm: 800.0,
            allowable_shear_mpa: 40.0,
            criterion: SizingCriterion::Shear,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Shaft Design Workbench right-side panel. A no-op when the
/// `show_shaftdesign_workbench` toggle is off.
pub fn draw_shaftdesign_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_shaftdesign_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_shaftdesign_workbench",
        "Shaft Design",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native combined bending + torsion shaft stress · valenx-shaftdesign",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.shaftdesign;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Shaft").strong());
                    ui.horizontal(|ui| {
                        ui.label("diameter d (mm)");
                        ui.add(egui::DragValue::new(&mut s.diameter_mm).speed(0.5));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Loads").strong());
                    ui.horizontal(|ui| {
                        ui.label("bending moment M (N·m)");
                        ui.add(egui::DragValue::new(&mut s.bending_moment_nm).speed(5.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("torque T (N·m)");
                        ui.add(egui::DragValue::new(&mut s.torque_nm).speed(5.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Material").strong());
                    ui.horizontal(|ui| {
                        ui.label("allowable shear (MPa)");
                        ui.add(egui::DragValue::new(&mut s.allowable_shear_mpa).speed(1.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Sizing criterion").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(
                            &mut s.criterion,
                            SizingCriterion::Shear,
                            "max shear (Tresca)",
                        );
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_shaftdesign(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative solid round shaft at the entered diameter as a 3-D solid and load it into the central viewport to orbit",
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
        app.show_shaftdesign_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.shaftdesign` borrow is
    // released here): build the shaft's 3-D solid and load it.
    if app.shaftdesign.show_3d_request {
        app.shaftdesign.show_3d_request = false;
        load_shaft_3d(app);
    }
}

/// Validate the form, evaluate the shaft and format the readout.
fn run_shaftdesign(s: &mut ShaftDesignWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build the validated [`ShaftSection`] from the form, converting the
/// diameter from millimetres to metres. Extracted so it is unit-testable
/// and shared by the readout and the 3-D gate.
fn section(s: &ShaftDesignWorkbenchState) -> Result<ShaftSection, String> {
    ShaftSection::new(s.diameter_mm / 1000.0, s.bending_moment_nm, s.torque_nm)
        .map_err(|e| e.to_string())
}

/// Evaluate the shaft and format the full readout, mapping any domain
/// error to a display string. Extracted so it is unit-testable.
fn compute(s: &ShaftDesignWorkbenchState) -> Result<String, String> {
    let sec = section(s)?;
    let sigma = bending_stress_section(&sec);
    let tau = shear_stress_section(&sec);
    let te = equivalent_torque(s.bending_moment_nm, s.torque_nm).map_err(|e| e.to_string())?;
    let tau_max = max_shear_stress(&sec);
    let sigma_1 = max_normal_stress(&sec);

    // Section properties (in metres) — the geometry behind the stresses:
    // the bending stress is M / Z, the torsional shear is T / Z_p, and the
    // polar second moment J governs the torsional stiffness / angle of
    // twist. Reported in engineering units: moduli in cm³, J in cm⁴.
    let z_modulus_cm3 = sec.section_modulus() * 1.0e6;
    let zp_modulus_cm3 = sec.polar_section_modulus() * 1.0e6;
    let polar_j_cm4 = sec.polar_second_moment() * 1.0e8;

    // Equivalent-torque (Tresca) sizing inverse: the minimum diameter that
    // keeps the combined shear within the allowable.
    let allowable_pa = s.allowable_shear_mpa * 1.0e6;
    let d_req_m = match s.criterion {
        SizingCriterion::Shear => {
            diameter_for_shear_stress(s.bending_moment_nm, s.torque_nm, allowable_pa)
                .map_err(|e| e.to_string())?
        }
    };
    let d_req_mm = d_req_m * 1000.0;
    let verdict = if s.diameter_mm + 1e-9 >= d_req_mm {
        "OK (d >= required)"
    } else {
        "UNDERSIZED (d < required)"
    };

    Ok(format!(
        "diameter d      : {:.1} mm\n\
         moment / torque : {:.1} / {:.1} N·m\n\
         allowable shear : {:.1} MPa\n\n\
         bending stress  : {:.3} MPa\n\
         torsional shear : {:.3} MPa\n\
         equivalent torq : {:.2} N·m\n\
         max shear τmax  : {:.3} MPa\n\
         max normal σ1   : {:.3} MPa\n\
         required d      : {:.2} mm\n\n\
         section mod Z   : {z_modulus_cm3:.4} cm³\n\
         polar mod Z_p   : {zp_modulus_cm3:.4} cm³\n\
         polar 2nd mom J : {polar_j_cm4:.4} cm⁴\n\n\
         verdict         : {verdict}",
        s.diameter_mm,
        s.bending_moment_nm,
        s.torque_nm,
        s.allowable_shear_mpa,
        sigma / 1.0e6,
        tau / 1.0e6,
        te,
        tau_max / 1.0e6,
        sigma_1 / 1.0e6,
        d_req_mm,
    ))
}

/// Append a closed coaxial round cylinder (axis along `z`) of `radius`
/// spanning `z0..z1`, faceted into `segs` sides with flat end caps, to
/// the buffers.
fn push_cylinder(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    radius: f64,
    z0: f64,
    z1: f64,
    segs: usize,
) {
    let bottom_centre = nodes.len();
    nodes.push(Vector3::new(0.0, 0.0, z0));
    let top_centre = nodes.len();
    nodes.push(Vector3::new(0.0, 0.0, z1));

    let ring_bottom = nodes.len();
    for k in 0..segs {
        let a = 2.0 * PI * (k as f64) / (segs as f64);
        nodes.push(Vector3::new(radius * a.cos(), radius * a.sin(), z0));
    }
    let ring_top = nodes.len();
    for k in 0..segs {
        let a = 2.0 * PI * (k as f64) / (segs as f64);
        nodes.push(Vector3::new(radius * a.cos(), radius * a.sin(), z1));
    }

    for k in 0..segs {
        let k1 = (k + 1) % segs;
        let b0 = ring_bottom + k;
        let b1 = ring_bottom + k1;
        let t0 = ring_top + k;
        let t1 = ring_top + k1;
        // Side wall (two triangles per facet).
        tris.extend_from_slice(&[b0, b1, t1, b0, t1, t0]);
        // Bottom cap (inward fan) and top cap (outward fan).
        tris.extend_from_slice(&[bottom_centre, b1, b0]);
        tris.extend_from_slice(&[top_centre, t0, t1]);
    }
}

/// Build the shaft as a triangle [`Mesh`] — a single solid round
/// cylinder at the entered diameter, faceted for display. Representative
/// geometry (the stress numbers are the `valenx-shaftdesign` result).
/// `None` for an invalid configuration.
fn shaft_solid_mesh(s: &ShaftDesignWorkbenchState) -> Option<Mesh> {
    let sec = section(s).ok()?;
    let radius = sec.diameter / 2.0;
    // A representative shaft length: eight diameters long, centred on z.
    let length = 8.0 * sec.diameter;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();
    push_cylinder(
        &mut nodes,
        &mut tris,
        radius,
        -0.5 * length,
        0.5 * length,
        48,
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-shaftdesign");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D shaft solid and load it into the central viewport.
fn load_shaft_3d(app: &mut ValenxApp) {
    let Some(mesh) = shaft_solid_mesh(&app.shaftdesign) else {
        app.shaftdesign.error =
            Some("shaft parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<shaft>/valenx-shaftdesign"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"shaftdesign"}`** product: the
/// representative solid round shaft built from the canonical 50 mm steel
/// line-shaft (M = 600 N·m, T = 800 N·m, 40 MPa allowable shear), paired with
/// the combined bending + torsion stress / section-property readout rows, at a
/// fixed 3/4 camera. Registered in [`crate::products_registry`]; the per-tool
/// builder the registry dispatches to. Pure — driven off
/// [`ShaftDesignWorkbenchState::default`].
pub(crate) fn shaftdesign_product() -> crate::WorkspaceProduct {
    let s = ShaftDesignWorkbenchState::default();
    let mesh = shaft_solid_mesh(&s).expect("canonical shaft ⇒ round-bar solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<shaft>/valenx-shaftdesign");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical shaft ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Shaft (bending + torsion)".into(),
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
        let s = ShaftDesignWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_combined_stresses() {
        let mut s = ShaftDesignWorkbenchState::default();
        run_shaftdesign(&mut s);
        assert!(
            s.error.is_none(),
            "default shaft should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("bending stress"));
        assert!(s.result.contains("torsional shear"));
        assert!(s.result.contains("max shear"));
        // d=50 mm, M=600, T=800: bending σ = 48.892 MPa, torsional
        // τ = 32.595 MPa, equivalent torque 1000.00 N·m, combined
        // τmax = 40.744 MPa, σ1 = 65.190 MPa, required d = 50.31 mm.
        assert!(s.result.contains("48.892"));
        assert!(s.result.contains("32.595"));
        assert!(s.result.contains("1000.00"));
        assert!(s.result.contains("40.744"));
        assert!(s.result.contains("65.190"));
        assert!(s.result.contains("50.31"));
    }

    #[test]
    fn analyze_rejects_zero_diameter() {
        let mut s = ShaftDesignWorkbenchState {
            diameter_mm: 0.0,
            ..Default::default()
        };
        run_shaftdesign(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn bending_stress_matches_32m_over_pi_d_cubed_ground_truth() {
        // Ground truth: outer-fibre bending stress σ = 32 M / (π d³).
        // d = 0.05 m, M = 600 N·m -> σ = 32*600/(π*0.05³) Pa.
        let d_m = 0.05_f64;
        let m = 600.0_f64;
        let sec = ShaftSection::new(d_m, m, 800.0).unwrap();
        let want = 32.0 * m / (PI * d_m.powi(3));
        assert!(
            (bending_stress_section(&sec) - want).abs() < 1e-6 * want,
            "σ {}, want {want}",
            bending_stress_section(&sec)
        );
        // ~48.89 MPa.
        assert!((want / 1.0e6 - 48.892_f64).abs() < 0.01);
    }

    #[test]
    fn section_properties_match_pi_d_powers_ground_truth() {
        // Ground truth for the solid round section (d = 0.05 m):
        //   section modulus     Z   = π d³ / 32  (m³)  → ×1e6 cm³
        //   polar sect. modulus Z_p = π d³ / 16  (m³)  → ×1e6 cm³
        //   polar second moment J   = π d⁴ / 32  (m⁴)  → ×1e8 cm⁴
        let d_m = 0.05_f64;
        let z_cm3 = PI * d_m.powi(3) / 32.0 * 1.0e6;
        let zp_cm3 = PI * d_m.powi(3) / 16.0 * 1.0e6;
        let j_cm4 = PI * d_m.powi(4) / 32.0 * 1.0e8;
        // Hand values: Z = 12.2718 cm³, Z_p = 24.5437 cm³, J = 61.3592 cm⁴.
        assert!((z_cm3 - 12.2718_f64).abs() < 1.0e-3, "Z {z_cm3} cm³");
        assert!((zp_cm3 - 24.5437_f64).abs() < 1.0e-3, "Z_p {zp_cm3} cm³");
        assert!((j_cm4 - 61.3592_f64).abs() < 1.0e-3, "J {j_cm4} cm⁴");
        // Z_p is exactly twice Z (π d³/16 = 2 · π d³/32).
        assert!((zp_cm3 - 2.0 * z_cm3).abs() < 1.0e-9);

        // The compute() readout surfaces all three with .4f precision.
        let s = ShaftDesignWorkbenchState::default();
        let out = compute(&s).expect("default shaft computes");
        assert!(out.contains("section mod Z   : 12.2718 cm³"), "{out}");
        assert!(out.contains("polar mod Z_p   : 24.5437 cm³"), "{out}");
        assert!(out.contains("polar 2nd mom J : 61.3592 cm⁴"), "{out}");
        // Consistency: the bending stress equals M / Z and the torsional
        // shear equals T / Z_p, both already shown in the readout — so the
        // new section moduli are the genuine geometry behind those stresses.
        let z_m3 = z_cm3 / 1.0e6;
        let zp_m3 = zp_cm3 / 1.0e6;
        let sigma_mpa = s.bending_moment_nm / z_m3 / 1.0e6;
        let tau_mpa = s.torque_nm / zp_m3 / 1.0e6;
        assert!((sigma_mpa - 48.892_f64).abs() < 0.01, "σ=M/Z {sigma_mpa}");
        assert!((tau_mpa - 32.595_f64).abs() < 0.01, "τ=T/Z_p {tau_mpa}");
    }

    #[test]
    fn shaft_mesh_for_default_is_nonempty_and_in_range() {
        let s = ShaftDesignWorkbenchState::default();
        let mesh = shaft_solid_mesh(&s).expect("default shaft yields a solid");
        assert!(mesh.nodes.len() > 8, "expected a faceted cylinder");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn shaft_mesh_none_for_invalid() {
        let s = ShaftDesignWorkbenchState {
            diameter_mm: 0.0,
            ..Default::default()
        };
        assert!(shaft_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_shaftdesign_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_shaftdesign_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_shaftdesign_workbench = true;
        run_shaftdesign(&mut app.shaftdesign);
        draw_workbench(&mut app);
    }
}
