//! The right-side **Bolted Joint Workbench** panel — native preloaded
//! bolted-joint mechanics over `valenx-bolt`.
//!
//! Mirrors the Heat Transfer / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_bolt_workbench`,
//! toggled from the View menu. The form sets a bolt (nominal diameter,
//! thread pitch, property grade), a tightening torque and nut factor, the
//! joint stiffness constant `C` and a service load. "Analyze" inverts
//! `T = K F d` to the achieved preload, sizes the ISO tensile-stress area
//! and reports the preload, axial / proof stress, the load-sharing split
//! and the separation and bolt-overload safety factors; "Show 3-D bolt"
//! loads a representative head + shank + nut solid into the central
//! viewport.
//!
//! This is the bolted-*joint mechanics* workbench (preload / torque /
//! stress), distinct from the ISO 4017 hex-bolt *dimensions* workbench.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;
use std::f64::consts::TAU;

use valenx_bolt::stress;
use valenx_bolt::{BoltGrade, BoltedJoint, NutFactor, StiffnessRatio};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Bolted Joint Workbench.
pub struct BoltWorkbenchState {
    /// Nominal bolt diameter `d` (mm).
    diameter_mm: f64,
    /// Thread pitch `P` (mm).
    pitch_mm: f64,
    /// Bolt property grade (ISO 898-1).
    grade: BoltGrade,
    /// Tightening torque `T` (N·m).
    torque_nm: f64,
    /// Nut factor `K` (torque coefficient).
    nut_factor: f64,
    /// Joint stiffness constant `C = kb / (kb + km)`.
    stiffness_c: f64,
    /// External service tensile load `P` (kN).
    service_load_kn: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D bolt solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for BoltWorkbenchState {
    fn default() -> Self {
        // An M10 class-8.8 bolt (d = 10 mm, P = 1.5 mm) tightened to
        // 50 N·m with the textbook K = 0.2: preload F = T / (K d) =
        // 50 / (0.2 * 0.010) = 25 000 N. Stiffness split C = 0.25, a
        // 6 kN service load.
        Self {
            diameter_mm: 10.0,
            pitch_mm: 1.5,
            grade: BoltGrade::Class8_8,
            torque_nm: 50.0,
            nut_factor: NutFactor::STEEL_AS_RECEIVED,
            stiffness_c: 0.25,
            service_load_kn: 6.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Bolted Joint Workbench right-side panel. A no-op when the
/// `show_bolt_workbench` toggle is off.
pub fn draw_bolt_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_bolt_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_bolt_workbench",
        "Bolted Joint",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native preloaded-joint mechanics (preload / torque / stress) · valenx-bolt",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.bolt;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Bolt").strong());
                    ui.horizontal(|ui| {
                        ui.label("nominal d (mm)");
                        ui.add(egui::DragValue::new(&mut s.diameter_mm).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("thread pitch P (mm)");
                        ui.add(egui::DragValue::new(&mut s.pitch_mm).speed(0.05));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Property grade (ISO 898-1)").strong());
                    ui.horizontal_wrapped(|ui| {
                        ui.radio_value(&mut s.grade, BoltGrade::Class4_6, "4.6");
                        ui.radio_value(&mut s.grade, BoltGrade::Class5_8, "5.8");
                        ui.radio_value(&mut s.grade, BoltGrade::Class8_8, "8.8");
                        ui.radio_value(&mut s.grade, BoltGrade::Class10_9, "10.9");
                        ui.radio_value(&mut s.grade, BoltGrade::Class12_9, "12.9");
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Tightening").strong());
                    ui.horizontal(|ui| {
                        ui.label("torque T (N·m)");
                        ui.add(egui::DragValue::new(&mut s.torque_nm).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("nut factor K");
                        ui.add(egui::DragValue::new(&mut s.nut_factor).speed(0.005));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Joint & load").strong());
                    ui.horizontal(|ui| {
                        ui.label("stiffness C = kb/(kb+km)");
                        ui.add(egui::DragValue::new(&mut s.stiffness_c).speed(0.01));
                    });
                    ui.horizontal(|ui| {
                        ui.label("service load P (kN)");
                        ui.add(egui::DragValue::new(&mut s.service_load_kn).speed(0.5));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_bolt(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D bolt").strong())
                        .on_hover_text(
                            "Build a representative bolt (hex head + threaded shank + nut) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Joint analysis").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_bolt_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.bolt` borrow is
    // released here): build the bolt's 3-D solid and load it.
    if app.bolt.show_3d_request {
        app.bolt.show_3d_request = false;
        load_bolt_3d(app);
    }
}

/// Validate the form, evaluate the joint and format the readout.
fn run_bolt(s: &mut BoltWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Evaluate the bolted joint and format the full readout, mapping any
/// domain error to a display string. Extracted so it is unit-testable.
fn compute(s: &BoltWorkbenchState) -> Result<String, String> {
    // SI: mm -> m, kN -> N.
    let d_m = s.diameter_mm / 1000.0;
    let pitch_m = s.pitch_mm / 1000.0;
    let load_n = s.service_load_kn * 1000.0;

    let k = NutFactor::new(s.nut_factor).map_err(|e| e.to_string())?;
    let c = StiffnessRatio::new(s.stiffness_c).map_err(|e| e.to_string())?;
    let joint = BoltedJoint::from_torque(s.torque_nm, k, d_m, c).map_err(|e| e.to_string())?;

    let material = s.grade.material().map_err(|e| e.to_string())?;
    let area_m2 = stress::tensile_stress_area(d_m, pitch_m).map_err(|e| e.to_string())?;
    let proof_n = stress::proof_load(&material, area_m2).map_err(|e| e.to_string())?;
    // Shigley's reuse target preload Fi = 0.75 Sp At — the recommended
    // tightening preload for a reused connection. Compare against the
    // achieved preload F to see how close the torque puts you.
    let rec_preload = stress::recommended_preload(&material, area_m2).map_err(|e| e.to_string())?;

    let preload = joint.preload_n();
    let preload_stress = stress::axial_stress(preload, area_m2).map_err(|e| e.to_string())?;
    let bolt_load = joint.bolt_load_n(load_n).map_err(|e| e.to_string())?;
    let bolt_stress = stress::axial_stress(bolt_load, area_m2).map_err(|e| e.to_string())?;
    let clamp = joint.clamping_force_n(load_n).map_err(|e| e.to_string())?;
    let p_sep = joint.separation_load_n();
    let sep_sf = joint
        .separation_safety_factor(load_n)
        .map_err(|e| e.to_string())?;
    let load_factor = joint
        .bolt_load_factor(proof_n, load_n)
        .map_err(|e| e.to_string())?;

    // Areas in mm^2, stresses in MPa for the readout.
    let area_mm2 = area_m2 * 1.0e6;
    let preload_mpa = preload_stress / 1.0e6;
    let proof_mpa = material.proof_strength_pa / 1.0e6;
    let bolt_mpa = bolt_stress / 1.0e6;
    let preload_frac = if proof_n > 0.0 {
        preload / proof_n * 100.0
    } else {
        0.0
    };
    // Achieved preload as a fraction of the recommended reuse target.
    let preload_util = if rec_preload > 0.0 {
        preload / rec_preload * 100.0
    } else {
        0.0
    };

    Ok(format!(
        "bolt            : M{:.0} grade {} (P = {:.2} mm)\n\
         tensile area At : {:.1} mm²\n\
         torque / K      : {:.1} N·m / {:.3}\n\
         stiffness C     : {:.3}  (bolt picks up {:.1}% of P)\n\n\
         preload F       : {:.0} N  ({:.0} MPa, {:.0}% of proof)\n\
         rec preload Fi  : {:.0} N  (0.75 Sp At, F = {:.0}% of Fi)\n\
         proof load Fp   : {:.0} N  ({:.0} MPa)\n\n\
         service load P  : {:.0} N\n\
         bolt load F+CP  : {:.0} N  ({:.0} MPa)\n\
         clamp F-(1-C)P  : {:.0} N\n\
         separation Psep : {:.0} N\n\
         sep safety n    : {:.2}\n\
         bolt load factor: {:.2}",
        s.diameter_mm,
        s.grade.label(),
        s.pitch_mm,
        area_mm2,
        s.torque_nm,
        s.nut_factor,
        s.stiffness_c,
        s.stiffness_c * 100.0,
        preload,
        preload_mpa,
        preload_frac,
        rec_preload,
        preload_util,
        proof_n,
        proof_mpa,
        load_n,
        bolt_load,
        bolt_mpa,
        clamp,
        p_sep,
        sep_sf,
        load_factor,
    ))
}

/// Append a capped cylinder along the x-axis, double-sided, centred at
/// `center` with half-length `half_len` and `radius`, `seg` segments.
fn push_cyl_x(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    center: Vector3<f64>,
    half_len: f64,
    radius: f64,
    seg: usize,
) {
    let (x0, x1) = (center.x - half_len, center.x + half_len);
    let lo = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(
            x0,
            center.y + radius * a.cos(),
            center.z + radius * a.sin(),
        ));
    }
    let hi = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(
            x1,
            center.y + radius * a.cos(),
            center.z + radius * a.sin(),
        ));
    }
    let cap0 = nodes.len();
    nodes.push(Vector3::new(x0, center.y, center.z));
    let cap1 = nodes.len();
    nodes.push(Vector3::new(x1, center.y, center.z));
    for j in 0..seg {
        let jn = (j + 1) % seg;
        // Side wall (double-sided).
        tris.extend_from_slice(&[
            lo + j,
            hi + j,
            hi + jn,
            lo + j,
            hi + jn,
            lo + jn,
            lo + j,
            hi + jn,
            hi + j,
            lo + j,
            lo + jn,
            hi + jn,
        ]);
        // Caps (double-sided fans).
        tris.extend_from_slice(&[cap0, lo + jn, lo + j, cap0, lo + j, lo + jn]);
        tris.extend_from_slice(&[cap1, hi + j, hi + jn, cap1, hi + jn, hi + j]);
    }
}

/// Build the bolt as a triangle [`Mesh`] — a hex head (a short, fat
/// 6-segment x-cylinder), a long slender threaded shank and a nut (a
/// short, fat 6-segment x-cylinder) run on near the far end.
/// Representative geometry (not to scale; the mechanics are the
/// `valenx-bolt` result). `None` for an invalid configuration.
fn bolt_solid_mesh(s: &BoltWorkbenchState) -> Option<Mesh> {
    // Gate on the same construction the readout uses: only build a solid
    // for a valid joint / thread.
    let d_m = s.diameter_mm / 1000.0;
    let pitch_m = s.pitch_mm / 1000.0;
    let k = NutFactor::new(s.nut_factor).ok()?;
    let c = StiffnessRatio::new(s.stiffness_c).ok()?;
    BoltedJoint::from_torque(s.torque_nm, k, d_m, c).ok()?;
    stress::tensile_stress_area(d_m, pitch_m).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Hex head at the -x end (6 segments => hexagonal prism).
    push_cyl_x(
        &mut nodes,
        &mut tris,
        Vector3::new(-1.85, 0.0, 0.0),
        0.3,
        0.9,
        6,
    );
    // Threaded shank along x.
    push_cyl_x(&mut nodes, &mut tris, Vector3::zeros(), 1.6, 0.5, 24);
    // Nut run on near the +x end (6 segments => hexagonal).
    push_cyl_x(
        &mut nodes,
        &mut tris,
        Vector3::new(1.35, 0.0, 0.0),
        0.25,
        0.85,
        6,
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-bolt");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D bolt solid and load it into the central viewport.
fn load_bolt_3d(app: &mut ValenxApp) {
    let Some(mesh) = bolt_solid_mesh(&app.bolt) else {
        app.bolt.error = Some("bolt parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<bolt>/valenx-bolt"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"bolt"}`** product: the representative bolt
/// (hex head + threaded shank + nut) built from the canonical M10 class-8.8
/// joint (tightened to 50 N·m, K = 0.2, C = 0.25, 6 kN service load), paired
/// with the preloaded-joint mechanics readout rows (preload / stresses /
/// load-sharing / safety factors), at a fixed 3/4 camera. Registered in
/// [`crate::products_registry`]; the per-tool builder the registry dispatches
/// to. Pure — driven off [`BoltWorkbenchState::default`].
pub(crate) fn bolt_product() -> crate::WorkspaceProduct {
    let s = BoltWorkbenchState::default();
    let mesh = bolt_solid_mesh(&s).expect("canonical bolt ⇒ head-shank-nut solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<bolt>/valenx-bolt");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical bolt ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Bolted joint (M10 8.8)".into(),
        lines,
        mesh: Some(loaded),
        vertex_colors: None,
        camera,
        kind2d: None,
        last_export: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_is_idle() {
        let s = BoltWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_preload_and_factors() {
        let mut s = BoltWorkbenchState::default();
        run_bolt(&mut s);
        assert!(
            s.error.is_none(),
            "default joint should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("preload F"));
        assert!(s.result.contains("separation Psep"));
        assert!(s.result.contains("bolt load factor"));
        // F = T / (K d) = 50 / (0.2 * 0.010) = 25 000 N -> "25000 N".
        assert!(s.result.contains("25000 N"), "readout: {}", s.result);
    }

    #[test]
    fn preload_matches_torque_relation_ground_truth() {
        // Ground truth: tightening to torque T at nut factor K on a bolt
        // of diameter d gives preload F = T / (K d). Hand-computed for
        // the default M10 joint: 50 / (0.2 * 0.010) = 25 000 N exactly.
        let s = BoltWorkbenchState::default();
        let d_m = s.diameter_mm / 1000.0;
        let k = NutFactor::new(s.nut_factor).unwrap();
        let c = StiffnessRatio::new(s.stiffness_c).unwrap();
        let joint = BoltedJoint::from_torque(s.torque_nm, k, d_m, c).unwrap();
        let expected = s.torque_nm / (s.nut_factor * d_m);
        assert!((joint.preload_n() - expected).abs() < 1e-6 * expected);
        assert!((joint.preload_n() - 25_000.0).abs() < 1e-6);
    }

    #[test]
    fn recommended_preload_matches_shigley_reuse_rule_ground_truth() {
        // Ground truth: the recommended reuse-joint preload is
        // Fi = 0.75 Sp At (Shigley). For the default M10 class-8.8 joint,
        // Sp = 600 MPa and At = (pi/4)(d - 0.938194 P)^2 with d = 0.010 m,
        // P = 0.0015 m:
        //   eff = 0.010 - 0.938194*0.0015 = 0.008592709 m
        //   At  = (pi/4)*eff^2           = 5.79895969e-5 m^2
        //   Fi  = 0.75 * 600e6 * At      = 26095.3186... N  -> "26095 N".
        // The achieved torque preload is F = 25000 N, so the readout
        // reports F = 96% of Fi (25000 / 26095.3186 * 100 = 95.80% -> 96).
        let mut s = BoltWorkbenchState::default();
        run_bolt(&mut s);
        assert!(
            s.error.is_none(),
            "default joint should analyze: {:?}",
            s.error
        );
        assert!(
            s.result.contains("rec preload Fi  : 26095 N"),
            "readout: {}",
            s.result
        );
        assert!(s.result.contains("F = 96% of Fi"), "readout: {}", s.result);

        // Hand-computed crate cross-check of the same number.
        let area = stress::tensile_stress_area(0.010, 0.0015).unwrap();
        let material = BoltGrade::Class8_8.material().unwrap();
        let fi = stress::recommended_preload(&material, area).unwrap();
        let expected = 0.75 * 600.0e6 * area;
        assert!(
            (fi - expected).abs() < 1e-6,
            "Fi {fi} vs expected {expected}"
        );
        assert!((fi - 26_095.318_605_830_347).abs() < 1e-3, "Fi was {fi}");
    }

    #[test]
    fn analyze_rejects_out_of_range_nut_factor() {
        // K must lie in (0, 1); zero is rejected by valenx-bolt.
        let mut s = BoltWorkbenchState {
            nut_factor: 0.0,
            ..Default::default()
        };
        run_bolt(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn analyze_rejects_degenerate_thread() {
        // Pitch coarser than the diameter -> non-positive effective
        // diameter -> tensile_stress_area errors.
        let mut s = BoltWorkbenchState {
            diameter_mm: 1.0,
            pitch_mm: 2.0,
            ..Default::default()
        };
        run_bolt(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn higher_grade_raises_proof_load() {
        // Same geometry / torque, a stronger class must report a larger
        // proof load (its readout therefore differs).
        let weak = BoltWorkbenchState {
            grade: BoltGrade::Class4_6,
            ..Default::default()
        };
        let strong = BoltWorkbenchState {
            grade: BoltGrade::Class12_9,
            ..Default::default()
        };
        let area = stress::tensile_stress_area(0.010, 0.0015).unwrap();
        let fp_weak = stress::proof_load(&weak.grade.material().unwrap(), area).unwrap();
        let fp_strong = stress::proof_load(&strong.grade.material().unwrap(), area).unwrap();
        assert!(fp_strong > fp_weak);
    }

    #[test]
    fn bolt_mesh_for_default_is_nonempty_and_in_range() {
        let s = BoltWorkbenchState::default();
        let mesh = bolt_solid_mesh(&s).expect("default bolt yields a solid");
        assert!(mesh.nodes.len() > 8, "expected head + shank + nut");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn bolt_mesh_none_for_invalid() {
        let s = BoltWorkbenchState {
            nut_factor: 0.0,
            ..Default::default()
        };
        assert!(bolt_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_bolt_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_bolt_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_bolt_workbench = true;
        run_bolt(&mut app.bolt);
        draw_workbench(&mut app);
    }
}
