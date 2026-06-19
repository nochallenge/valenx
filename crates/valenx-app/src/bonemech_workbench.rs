//! The right-side **Bone Mechanics Workbench** panel — native long-bone
//! shaft stress analysis over `valenx-bonemech`.
//!
//! Mirrors the Heat Transfer / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_bonemech_workbench`,
//! toggled from the View menu. The form idealises a long-bone diaphysis as
//! a **hollow circular tube** (a cortical wall around the medullary canal)
//! loaded by a bending moment and an axial force, with an apparent-density
//! strength scaling. "Analyze" reports the annular section's second moment
//! of area `I`, its section modulus `S`, the Euler-Bernoulli outer-fibre
//! bending stress, the axial stress, the density-scaled ultimate stress,
//! the bending fracture moment, and the load-vs-fracture utilisation /
//! safety factor. "Show 3-D bone" loads a representative hollow-shaft solid
//! into the central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;
use std::f64::consts::{PI, TAU};

use valenx_bonemech::{
    bending_moment_for_stress, bending_stress_mpa, second_moment_hollow_circle_mm4,
    section_modulus_mm3, Bone, PowerLaw, CORTICAL_MODULUS_GPA, CORTICAL_ULTIMATE_STRESS_MPA,
};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Bone Mechanics Workbench.
pub struct BonemechWorkbenchState {
    /// Shaft outer diameter `D_o` (mm).
    outer_d_mm: f64,
    /// Medullary-canal (bore) inner diameter `D_i` (mm).
    inner_d_mm: f64,
    /// Applied bending moment `M` (N·m). Converted to N·mm internally.
    moment_nm: f64,
    /// Applied axial force `F` (N): positive tensile, negative compressive.
    axial_force_n: f64,
    /// Apparent density `rho` of the bone tissue (g/cm^3) for the
    /// Carter-Hayes strength scaling.
    apparent_density: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D bone solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for BonemechWorkbenchState {
    fn default() -> Self {
        // A femoral-shaft-like cortical tube: outer 27 mm, inner 17 mm,
        // under a 60 N·m bending moment and a 1500 N axial load, at the
        // reference apparent density (1.9 g/cm^3 -> 150 MPa ultimate).
        // Bending stress ~36.84 MPa (well below the 150 MPa ultimate),
        // safety factor ~4.07 in bending.
        Self {
            outer_d_mm: 27.0,
            inner_d_mm: 17.0,
            moment_nm: 60.0,
            axial_force_n: 1500.0,
            apparent_density: 1.9,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Bone Mechanics Workbench right-side panel. A no-op when the
/// `show_bonemech_workbench` toggle is off.
pub fn draw_bonemech_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_bonemech_workbench {
        return;
    }

    egui::SidePanel::right("valenx_bonemech_workbench")
        .resizable(true)
        .default_width(360.0)
        .width_range(300.0..=560.0)
        .show(ctx, |ui| {
            if crate::workbench_ui::header(
                ui,
                "Bone Mechanics",
                "native long-bone hollow-shaft bending + axial stress · valenx-bonemech",
            ) {
                app.show_bonemech_workbench = false;
            }

            let s = &mut app.bonemech;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Shaft cross-section").strong());
                    ui.horizontal(|ui| {
                        ui.label("outer Ø (mm)");
                        ui.add(egui::DragValue::new(&mut s.outer_d_mm).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("inner Ø (mm)");
                        ui.add(egui::DragValue::new(&mut s.inner_d_mm).speed(0.5));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Loads").strong());
                    ui.horizontal(|ui| {
                        ui.label("bending moment (N·m)");
                        ui.add(egui::DragValue::new(&mut s.moment_nm).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("axial force (N)");
                        ui.add(egui::DragValue::new(&mut s.axial_force_n).speed(10.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Tissue").strong());
                    ui.horizontal(|ui| {
                        ui.label("apparent density (g/cm³)");
                        ui.add(egui::DragValue::new(&mut s.apparent_density).speed(0.02));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_bonemech(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D bone").strong())
                        .on_hover_text(
                            "Build a representative hollow long-bone shaft (cortical wall around the medullary canal) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Section & strength").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        });

    // Serviced after the panel draws (the `&mut app.bonemech` borrow is
    // released here): build the bone's 3-D solid and load it.
    if app.bonemech.show_3d_request {
        app.bonemech.show_3d_request = false;
        load_bone_3d(app);
    }
}

/// Validate the form, evaluate the section and format the readout.
fn run_bonemech(s: &mut BonemechWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// The annular cross-sectional ring area `A = (pi / 4)(D_o^2 - D_i^2)`, in
/// mm^2 — the area carrying the axial load. Returns the area only for a
/// valid hollow-tube geometry (`0 <= D_i < D_o`, both finite). Extracted so
/// it is unit-testable and shared with the axial-stress path.
fn ring_area_mm2(s: &BonemechWorkbenchState) -> Result<f64, String> {
    // Reuse the crate's geometry validation (rejects D_o<=0, D_i<0,
    // D_i>=D_o, non-finite) so the area and the second moment of area stay
    // consistent about what counts as a valid section.
    second_moment_hollow_circle_mm4(s.outer_d_mm, s.inner_d_mm).map_err(|e| e.to_string())?;
    Ok(PI / 4.0 * (s.outer_d_mm.powi(2) - s.inner_d_mm.powi(2)))
}

/// Evaluate the section and format the full readout, mapping any domain
/// error to a display string. Extracted so it is unit-testable.
fn compute(s: &BonemechWorkbenchState) -> Result<String, String> {
    // Section geometry of the annulus.
    let i_mm4 =
        second_moment_hollow_circle_mm4(s.outer_d_mm, s.inner_d_mm).map_err(|e| e.to_string())?;
    let c_mm = s.outer_d_mm / 2.0; // peak fibre at the outer surface
    let s_mm3 = section_modulus_mm3(i_mm4, c_mm).map_err(|e| e.to_string())?;
    let area_mm2 = ring_area_mm2(s)?;

    // Density-scaled ultimate stress via the Carter-Hayes squared power
    // law, anchored at the cortical reference (1.9 g/cm^3 -> 150 MPa).
    let law =
        PowerLaw::carter_hayes(1.9, CORTICAL_ULTIMATE_STRESS_MPA).map_err(|e| e.to_string())?;
    let ult_mpa = law
        .strength(s.apparent_density)
        .map_err(|e| e.to_string())?;
    let bone = Bone::new(area_mm2, CORTICAL_MODULUS_GPA, ult_mpa).map_err(|e| e.to_string())?;

    // Bending: Euler-Bernoulli outer-fibre stress and the fracture moment
    // that brings that fibre to the ultimate stress.
    let moment_nmm = s.moment_nm * 1000.0; // N·m -> N·mm
    let bend_mpa = bending_stress_mpa(moment_nmm, c_mm, i_mm4).map_err(|e| e.to_string())?;
    let strain = bone.strain(bend_mpa).map_err(|e| e.to_string())?;
    let m_fracture_nmm =
        bending_moment_for_stress(ult_mpa, c_mm, i_mm4).map_err(|e| e.to_string())?;
    let bend_util = bend_mpa.abs() / ult_mpa;
    let bend_sf = if bend_mpa == 0.0 {
        f64::INFINITY
    } else {
        ult_mpa / bend_mpa.abs()
    };

    // Axial: stress, ultimate force and utilisation / safety factor.
    let axial_mpa = bone
        .axial_stress_mpa(s.axial_force_n)
        .map_err(|e| e.to_string())?;
    let f_ult_n = bone.fracture_load_n();
    let axial_util = bone
        .utilisation(s.axial_force_n)
        .map_err(|e| e.to_string())?;
    let axial_sf = bone
        .safety_factor(s.axial_force_n)
        .map_err(|e| e.to_string())?;

    Ok(format!(
        "shaft Ø out/in  : {:.1} / {:.1} mm\n\
         ring area A     : {:.2} mm²\n\
         2nd moment I    : {:.1} mm⁴\n\
         extreme fibre c : {:.2} mm\n\
         section mod S   : {:.1} mm³\n\
         apparent ρ      : {:.2} g/cm³\n\
         ultimate stress : {:.1} MPa\n\n\
         bending moment  : {:.2} N·m\n\
         bending stress  : {:.2} MPa\n\
         bending strain  : {:.6}\n\
         fracture moment : {:.2} N·m\n\
         bend util / SF  : {:.3} / {:.2}\n\n\
         axial force     : {:.1} N\n\
         axial stress    : {:.2} MPa\n\
         fracture load   : {:.1} N\n\
         axial util / SF : {:.3} / {:.2}",
        s.outer_d_mm,
        s.inner_d_mm,
        area_mm2,
        i_mm4,
        c_mm,
        s_mm3,
        s.apparent_density,
        ult_mpa,
        s.moment_nm,
        bend_mpa,
        strain,
        m_fracture_nmm / 1000.0,
        bend_util,
        bend_sf,
        s.axial_force_n,
        axial_mpa,
        f_ult_n,
        axial_util,
        axial_sf,
    ))
}

/// Build the hollow long-bone shaft as a triangle [`Mesh`]: an annular
/// cylinder of outer radius `D_o / 2` and bore `D_i / 2`, extruded along
/// the bone's long (z) axis, with both the outer wall and the inner canal
/// wall double-sided so the hollow is visible when sectioned. Geometry is
/// to scale in the cross-section (the diameters are the real inputs); the
/// stress numbers are the `valenx-bonemech` result. `None` for an invalid
/// configuration (the crate's geometry validation rejects it).
fn bone_solid_mesh(s: &BonemechWorkbenchState) -> Option<Mesh> {
    // Gate on the same geometry validation the analysis uses.
    second_moment_hollow_circle_mm4(s.outer_d_mm, s.inner_d_mm).ok()?;

    let r_out = (s.outer_d_mm / 2.0) as f32;
    let r_in = (s.inner_d_mm / 2.0) as f32;
    // Make the shaft a few outer-diameters long so it reads as a long bone.
    let half_len = (s.outer_d_mm * 2.0) as f32;
    let seg = 48usize;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Two rings (z = -half_len, +half_len) for each of the outer and inner
    // walls, laid out as [outer ring -z, outer ring +z, inner ring -z,
    // inner ring +z].
    let ring = |radius: f32, z: f32, nodes: &mut Vec<Vector3<f64>>| -> usize {
        let base = nodes.len();
        for k in 0..seg {
            let a = TAU as f32 * (k as f32) / (seg as f32);
            nodes.push(Vector3::new(
                (radius * a.cos()) as f64,
                (radius * a.sin()) as f64,
                z as f64,
            ));
        }
        base
    };
    let out_lo = ring(r_out, -half_len, &mut nodes);
    let out_hi = ring(r_out, half_len, &mut nodes);
    let in_lo = ring(r_in, -half_len, &mut nodes);
    let in_hi = ring(r_in, half_len, &mut nodes);

    // Double-sided quad strip between two rings of `seg` nodes each.
    let wall = |lo: usize, hi: usize, tris: &mut Vec<usize>| {
        for k in 0..seg {
            let k1 = (k + 1) % seg;
            let a = lo + k;
            let b = lo + k1;
            let c = hi + k1;
            let d = hi + k;
            // Front and back faces (double-sided).
            tris.extend_from_slice(&[a, b, c, a, c, d]);
            tris.extend_from_slice(&[a, c, b, a, d, c]);
        }
    };
    wall(out_lo, out_hi, &mut tris); // outer cortical wall
    wall(in_lo, in_hi, &mut tris); // medullary canal wall

    // End caps: an annular ring (cortical wall thickness) at each end,
    // double-sided.
    let end_cap = |out_ring: usize, in_ring: usize, tris: &mut Vec<usize>| {
        for k in 0..seg {
            let k1 = (k + 1) % seg;
            let oa = out_ring + k;
            let ob = out_ring + k1;
            let ia = in_ring + k;
            let ib = in_ring + k1;
            tris.extend_from_slice(&[oa, ob, ib, oa, ib, ia]);
            tris.extend_from_slice(&[oa, ib, ob, oa, ia, ib]);
        }
    };
    end_cap(out_lo, in_lo, &mut tris);
    end_cap(out_hi, in_hi, &mut tris);

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-bonemech");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D bone solid and load it into the central viewport.
fn load_bone_3d(app: &mut ValenxApp) {
    let Some(mesh) = bone_solid_mesh(&app.bonemech) else {
        app.bonemech.error = Some("shaft geometry is invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<bone>/valenx-bonemech"),
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
        let s = BonemechWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_section_and_strength() {
        let mut s = BonemechWorkbenchState::default();
        run_bonemech(&mut s);
        assert!(
            s.error.is_none(),
            "default shaft should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("bending stress"));
        assert!(s.result.contains("section mod S"));
        assert!(s.result.contains("axial stress"));
        // Outer 27 / inner 17 mm under 60 N·m: outer-fibre bending stress
        // ~36.84 MPa, at the reference density (150 MPa ultimate) the
        // bending safety factor is ~4.07.
        assert!(s.result.contains("36.84"), "got:\n{}", s.result);
        assert!(s.result.contains("150.0"), "got:\n{}", s.result);
        assert!(s.result.contains("4.07"), "got:\n{}", s.result);
    }

    #[test]
    fn analyze_rejects_bore_exceeding_outer() {
        // Inner diameter >= outer is a geometry inconsistency the crate
        // rejects.
        let mut s = BonemechWorkbenchState {
            outer_d_mm: 20.0,
            inner_d_mm: 22.0,
            ..Default::default()
        };
        run_bonemech(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn bending_stress_matches_flexure_ground_truth() {
        // Ground truth: the Euler-Bernoulli flexure formula sigma = M c / I
        // and the section modulus S = I / c (so sigma = M / S agree). With
        // I = 1000 mm^4, c = 5 mm and M = 20000 N·mm the outer-fibre stress
        // is exactly 100 MPa and S = 200 mm^3.
        let i_mm4: f64 = 1000.0;
        let c_mm: f64 = 5.0;
        let moment_nmm: f64 = 20_000.0;
        let sigma = bending_stress_mpa(moment_nmm, c_mm, i_mm4).unwrap();
        assert!((sigma - 100.0).abs() < 1e-9, "got {sigma}");
        let s_mod = section_modulus_mm3(i_mm4, c_mm).unwrap();
        assert!((s_mod - 200.0).abs() < 1e-9, "got {s_mod}");
        assert!((moment_nmm / s_mod - sigma).abs() < 1e-9);
    }

    #[test]
    fn ring_area_is_annulus_formula() {
        // A = (pi/4)(D_o^2 - D_i^2): outer 27, inner 17 -> pi/4 * (729-289)
        // = pi/4 * 440 ~= 345.5752 mm^2.
        let s = BonemechWorkbenchState::default();
        let area: f64 = ring_area_mm2(&s).unwrap();
        let expected: f64 = PI / 4.0 * (27f64.powi(2) - 17f64.powi(2));
        assert!((area - expected).abs() < 1e-9, "got {area}");
        assert!((area - 345.575_191_894_877_3).abs() < 1e-6, "got {area}");
    }

    #[test]
    fn bone_mesh_for_default_is_nonempty_and_in_range() {
        let s = BonemechWorkbenchState::default();
        let mesh = bone_solid_mesh(&s).expect("default shaft yields a solid");
        assert!(mesh.nodes.len() > 8, "expected an extruded annular tube");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn bone_mesh_none_for_invalid() {
        let s = BonemechWorkbenchState {
            outer_d_mm: 20.0,
            inner_d_mm: 25.0,
            ..Default::default()
        };
        assert!(bone_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_bonemech_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_bonemech_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_bonemech_workbench = true;
        run_bonemech(&mut app.bonemech);
        draw_workbench(&mut app);
    }
}
