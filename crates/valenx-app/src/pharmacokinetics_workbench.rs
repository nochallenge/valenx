//! The right-side **Pharmacokinetics Workbench** panel — native
//! one-compartment IV-bolus dosing analysis over `valenx-pharmacokinetics`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on
//! `crate::ValenxApp::show_pharmacokinetics_workbench`, toggled from the
//! View menu. The form sets a single-compartment IV bolus (dose, apparent
//! volume of distribution `V`, clearance `CL`) plus a sample time and a
//! threshold concentration; "Analyze" reports the elimination rate, peak
//! concentration, half-life, total exposure (AUC) and the concentration /
//! cumulative-AUC at the sample time, and "Show 3-D vial" loads a
//! representative drug-vial solid into the central viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_pharmacokinetics::onecompartment::OneCompartmentIv;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Pharmacokinetics Workbench.
pub struct PharmacokineticsWorkbenchState {
    /// Administered IV-bolus dose (mg).
    dose_mg: f64,
    /// Apparent volume of distribution `V` (L).
    volume_l: f64,
    /// Clearance `CL` (L/h).
    clearance_l_per_h: f64,
    /// Sample time at which to report `C(t)` and the cumulative AUC (h).
    sample_time_h: f64,
    /// Threshold concentration for the time-to-threshold solve (mg/L).
    threshold_mg_per_l: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D vial solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for PharmacokineticsWorkbenchState {
    fn default() -> Self {
        // A 500 mg IV bolus into V = 35 L with CL = 7 L/h:
        // k = CL/V = 0.2 /h, Cmax = dose/V ≈ 14.29 mg/L, t½ = ln2/k ≈
        // 3.47 h, AUC = dose/CL ≈ 71.43 mg·h/L. Sample at t = 4 h, with a
        // 1 mg/L therapeutic-floor threshold.
        Self {
            dose_mg: 500.0,
            volume_l: 35.0,
            clearance_l_per_h: 7.0,
            sample_time_h: 4.0,
            threshold_mg_per_l: 1.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Pharmacokinetics Workbench right-side panel. A no-op when the
/// `show_pharmacokinetics_workbench` toggle is off.
pub fn draw_pharmacokinetics_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_pharmacokinetics_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_pharmacokinetics_workbench",
        "Pharmacokinetics",
        |app, ui| {
            ui.label(
                egui::RichText::new("native one-compartment PK dosing · valenx-pharmacokinetics")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.pharmacokinetics;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Dose").strong());
                    ui.horizontal(|ui| {
                        ui.label("dose (mg)");
                        ui.add(egui::DragValue::new(&mut s.dose_mg).speed(5.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Disposition").strong());
                    ui.horizontal(|ui| {
                        ui.label("volume Vd (L)");
                        ui.add(egui::DragValue::new(&mut s.volume_l).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("clearance CL (L/h)");
                        ui.add(egui::DragValue::new(&mut s.clearance_l_per_h).speed(0.25));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Sampling").strong());
                    ui.horizontal(|ui| {
                        ui.label("sample time (h)");
                        ui.add(egui::DragValue::new(&mut s.sample_time_h).speed(0.25));
                    });
                    ui.horizontal(|ui| {
                        ui.label("threshold (mg/L)");
                        ui.add(egui::DragValue::new(&mut s.threshold_mg_per_l).speed(0.1));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_pk(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D vial").strong())
                        .on_hover_text(
                            "Build a representative drug vial (capped cylinder body, neck cap and base) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Dosing profile").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_pharmacokinetics_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.pharmacokinetics` borrow
    // is released here): build the vial's 3-D solid and load it.
    if app.pharmacokinetics.show_3d_request {
        app.pharmacokinetics.show_3d_request = false;
        load_pk_3d(app);
    }
}

/// Validate the form, evaluate the model and format the readout.
fn run_pk(s: &mut PharmacokineticsWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Evaluate the one-compartment IV model and format the full readout,
/// mapping any domain error to a display string. Extracted so it is
/// unit-testable.
fn compute(s: &PharmacokineticsWorkbenchState) -> Result<String, String> {
    let m = OneCompartmentIv::new(s.dose_mg, s.volume_l, s.clearance_l_per_h)
        .map_err(|e| e.to_string())?;
    let k = m.elimination_rate();
    let cmax = m.cmax();
    let t_half = m.half_life();
    let auc = m.auc();
    let c_t = m
        .concentration(s.sample_time_h)
        .map_err(|e| e.to_string())?;
    let auc_t = m.auc_to(s.sample_time_h).map_err(|e| e.to_string())?;
    let t_thresh = m
        .time_to_threshold(s.threshold_mg_per_l)
        .map_err(|e| e.to_string())?;

    Ok(format!(
        "dose            : {:.1} mg\n\
         volume Vd       : {:.2} L\n\
         clearance CL    : {:.3} L/h\n\n\
         elim. rate k    : {:.4} /h\n\
         half-life t½    : {:.3} h\n\
         Cmax = C(0)     : {:.3} mg/L\n\
         AUC(0→∞)        : {:.3} mg·h/L\n\n\
         at t = {:.2} h\n\
         C(t)            : {:.4} mg/L\n\
         AUC(0→t)        : {:.3} mg·h/L\n\n\
         time to {:.2} mg/L: {:.3} h",
        s.dose_mg,
        s.volume_l,
        s.clearance_l_per_h,
        k,
        t_half,
        cmax,
        auc,
        s.sample_time_h,
        c_t,
        auc_t,
        s.threshold_mg_per_l,
        t_thresh,
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

/// Append a double-sided cylinder of `height` along `+z` rising from `base`
/// (radius `r`, `seg` segments) to the buffers.
fn push_cyl_z(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    base: Vector3<f64>,
    height: f64,
    r: f64,
    seg: usize,
) {
    let (z0, z1) = (base.z, base.z + height);
    let bot = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(base.x + r * a.cos(), base.y + r * a.sin(), z0));
    }
    let top = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(base.x + r * a.cos(), base.y + r * a.sin(), z1));
    }
    for j in 0..seg {
        let jn = (j + 1) % seg;
        tris.extend_from_slice(&[
            bot + j,
            top + j,
            top + jn,
            bot + j,
            top + jn,
            bot + jn,
            bot + j,
            top + jn,
            top + j,
            bot + j,
            bot + jn,
            top + jn,
        ]);
    }
}

/// Build a representative drug vial as a triangle [`Mesh`] — a capped
/// cylinder body (the liquid container) on the `+z` axis, a neck-cap box
/// (the stopper / crimp seal) and a base. Representative geometry (not to
/// scale; the dosing numbers are the `valenx-pharmacokinetics` result).
/// `None` for an invalid configuration.
fn vial_solid_mesh(s: &PharmacokineticsWorkbenchState) -> Option<Mesh> {
    // Gate the 3-D build on a model that constructs (positive V and CL,
    // non-negative dose) — the same validation the readout runs.
    OneCompartmentIv::new(s.dose_mg, s.volume_l, s.clearance_l_per_h).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Vial body (capped cylinder on the +z axis).
    push_cyl_z(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.12),
        0.7,
        0.22,
        32,
    );
    // Neck / crimp-seal cap.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.88),
        Vector3::new(0.1, 0.1, 0.06),
    );
    // Base.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.06),
        Vector3::new(0.26, 0.26, 0.06),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-pharmacokinetics");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D vial solid and load it into the central viewport.
fn load_pk_3d(app: &mut ValenxApp) {
    let Some(mesh) = vial_solid_mesh(&app.pharmacokinetics) else {
        app.pharmacokinetics.error =
            Some("dosing parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<pk>/valenx-pharmacokinetics"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"pharmacokinetics"}`** product: the
/// canonical drug-vial solid (the panel's "Show 3-D vial" geometry) paired
/// with the workbench's own one-compartment dosing-profile headline numbers,
/// at a fixed 3/4 camera. Registered in [`crate::products_registry`]; the
/// per-tool builder the registry dispatches to. Pure — driven off
/// [`PharmacokineticsWorkbenchState::default`].
///
/// The readout rows mirror the panel's `compute()` dosing-profile readout.
pub(crate) fn pharmacokinetics_product() -> crate::WorkspaceProduct {
    let s = PharmacokineticsWorkbenchState::default();
    let mesh = vial_solid_mesh(&s).expect("default IV-bolus dosing ⇒ a 3-D vial");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<pk>/valenx-pharmacokinetics");
    let readout = compute(&s).expect("default IV-bolus dosing ⇒ a valid readout");
    let lines = crate::products_registry::lines_from_readout(&readout);
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Pharmacokinetics".into(),
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
        let s = PharmacokineticsWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_half_life_cmax_and_auc() {
        let mut s = PharmacokineticsWorkbenchState::default();
        run_pk(&mut s);
        assert!(
            s.error.is_none(),
            "default dosing should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("half-life"));
        assert!(s.result.contains("Cmax"));
        assert!(s.result.contains("AUC(0→∞)"));
        assert!(s.result.contains("elim. rate k"));
        // dose 500 / V 35 = Cmax ≈ 14.286 mg/L (printed to 3 decimals).
        assert!(s.result.contains("14.286"));
    }

    #[test]
    fn analyze_rejects_zero_volume() {
        let mut s = PharmacokineticsWorkbenchState {
            volume_l: 0.0,
            ..Default::default()
        };
        run_pk(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn elimination_rate_half_life_and_clearance_identities() {
        // GROUND TRUTH: k = CL / V, t½ = ln(2)/k, AUC = dose / CL, and the
        // clearance identity CL = k · V all hold to machine precision.
        let s = PharmacokineticsWorkbenchState::default();
        let m = OneCompartmentIv::new(s.dose_mg, s.volume_l, s.clearance_l_per_h).unwrap();
        let k = m.elimination_rate();
        assert!((k - s.clearance_l_per_h / s.volume_l).abs() < 1e-9);
        assert!((m.half_life() - std::f64::consts::LN_2 / k).abs() < 1e-9);
        assert!((m.auc() - s.dose_mg / s.clearance_l_per_h).abs() < 1e-9);
        // Clearance recovered from rate and volume.
        assert!((k * s.volume_l - s.clearance_l_per_h).abs() < 1e-9);
    }

    #[test]
    fn vial_mesh_for_default_is_nonempty_and_in_range() {
        let s = PharmacokineticsWorkbenchState::default();
        let mesh = vial_solid_mesh(&s).expect("default vial yields a solid");
        assert!(mesh.nodes.len() > 8, "expected cylinder body + cap + base");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn vial_mesh_none_for_invalid() {
        let s = PharmacokineticsWorkbenchState {
            clearance_l_per_h: 0.0,
            ..Default::default()
        };
        assert!(vial_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_pharmacokinetics_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_pharmacokinetics_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_pharmacokinetics_workbench = true;
        run_pk(&mut app.pharmacokinetics);
        draw_workbench(&mut app);
    }
}
