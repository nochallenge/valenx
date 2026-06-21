//! The right-side **Fatigue Workbench** panel — native high-cycle
//! stress-life analysis over `valenx-fatigue`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_fatigue_workbench`,
//! toggled from the View menu. The form sets a material (endurance limit
//! `Se`, yield `Sy`, ultimate `Su`), an operating point (alternating
//! stress `sigma_a`, mean stress `sigma_m`) and a constant-life criterion
//! (Goodman / Soderberg / Gerber). "Analyze" reports the factor of safety
//! against the chosen constant-life line, the equivalent fully-reversed
//! stress, and the Basquin S-N life it predicts; "Show 3-D specimen"
//! loads a representative round dog-bone fatigue test specimen into the
//! central viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_fatigue::{Life, Material, MeanStressCriterion, SnCurve};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Fatigue Workbench.
pub struct FatigueWorkbenchState {
    /// Corrected endurance limit `Se` (consistent stress unit, e.g. MPa).
    endurance_limit: f64,
    /// Yield strength `Sy` (Soderberg intercept).
    yield_strength: f64,
    /// Ultimate tensile strength `Su` (Goodman / Gerber intercept).
    ultimate_strength: f64,
    /// Applied alternating stress amplitude `sigma_a`.
    sigma_a: f64,
    /// Applied (tensile) mean stress `sigma_m`.
    sigma_m: f64,
    /// Which constant-life criterion to apply.
    criterion: MeanStressCriterion,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D specimen solid (serviced after
    /// the panel draws).
    show_3d_request: bool,
}

impl Default for FatigueWorkbenchState {
    fn default() -> Self {
        // A medium-strength steel (Se = 200, Sy = 350, Su = 500 MPa) under
        // a tensile duty cycle sigma_a = 120, sigma_m = 80 MPa. Goodman:
        // n ~ 1.32 (inside the line), equivalent reversed stress ~ 142.9
        // MPa, which is below Se -> infinite life on the capped S-N curve.
        Self {
            endurance_limit: 200.0,
            yield_strength: 350.0,
            ultimate_strength: 500.0,
            sigma_a: 120.0,
            sigma_m: 80.0,
            criterion: MeanStressCriterion::Goodman,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Fatigue Workbench right-side panel. A no-op when the
/// `show_fatigue_workbench` toggle is off.
pub fn draw_fatigue_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_fatigue_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_fatigue_workbench",
        "Fatigue",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native high-cycle stress-life (S-N) analysis · valenx-fatigue",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.fatigue;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Material strengths").strong());
                    ui.horizontal(|ui| {
                        ui.label("endurance Se (MPa)");
                        ui.add(egui::DragValue::new(&mut s.endurance_limit).speed(2.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("yield Sy (MPa)");
                        ui.add(egui::DragValue::new(&mut s.yield_strength).speed(2.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("ultimate Su (MPa)");
                        ui.add(egui::DragValue::new(&mut s.ultimate_strength).speed(2.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Operating stress").strong());
                    ui.horizontal(|ui| {
                        ui.label("alternating σa (MPa)");
                        ui.add(egui::DragValue::new(&mut s.sigma_a).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("mean σm (MPa)");
                        ui.add(egui::DragValue::new(&mut s.sigma_m).speed(1.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Constant-life criterion").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(
                            &mut s.criterion,
                            MeanStressCriterion::Goodman,
                            "Goodman",
                        );
                        ui.radio_value(
                            &mut s.criterion,
                            MeanStressCriterion::Soderberg,
                            "Soderberg",
                        );
                        ui.radio_value(&mut s.criterion, MeanStressCriterion::Gerber, "Gerber");
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_fatigue(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D specimen").strong())
                        .on_hover_text(
                            "Build a representative round dog-bone fatigue test specimen (wide grip ends, narrow gauge section) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Fatigue").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_fatigue_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.fatigue` borrow is
    // released here): build the specimen's 3-D solid and load it.
    if app.fatigue.show_3d_request {
        app.fatigue.show_3d_request = false;
        load_specimen_3d(app);
    }
}

/// Validate the form, evaluate the fatigue point and format the readout.
fn run_fatigue(s: &mut FatigueWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build the validated [`Material`], the quantity both the readout and
/// the 3-D gate need. Extracted so it is unit-testable and shared.
fn material(s: &FatigueWorkbenchState) -> Result<Material, String> {
    Material::new(s.endurance_limit, s.yield_strength, s.ultimate_strength)
        .map_err(|e| e.to_string())
}

/// Evaluate the operating point and format the full readout, mapping any
/// domain error to a display string. Extracted so it is unit-testable.
///
/// Reports the factor of safety against the chosen constant-life line,
/// the equivalent completely-reversed stress, the allowable alternating
/// stress at the operating mean stress (the Haigh-diagram `n = 1` ceiling,
/// with the operating-point margin below it), and the Basquin S-N life:
/// the curve is fit through `(1e3, 0.9 Su)` and the fatigue-limit point
/// `(1e6, Se)`, then capped at `Se`, so an equivalent reversed stress at
/// or below the endurance limit reads as infinite life.
fn compute(s: &FatigueWorkbenchState) -> Result<String, String> {
    let mat = material(s)?;
    let crit = s.criterion;

    let n = mat
        .factor_of_safety(crit, s.sigma_a, s.sigma_m)
        .map_err(|e| e.to_string())?;
    let sigma_ar = mat
        .equivalent_reversed_stress(crit, s.sigma_a, s.sigma_m)
        .map_err(|e| e.to_string())?;
    // The Haigh-diagram allowable alternating stress at this mean stress on
    // the n = 1 constant-life line: the largest σa the part tolerates at σm
    // before it reaches the failure line. The operating σa headroom below it
    // is the alternating-stress design margin.
    let sigma_a_allow = mat
        .allowable_alternating(crit, s.sigma_m, 1.0)
        .map_err(|e| e.to_string())?;
    let sigma_a_margin = sigma_a_allow - s.sigma_a;

    // Basquin S-N curve through (1e3, 0.9*Su) and (1e6, Se), capped at the
    // endurance limit Se.
    let curve =
        SnCurve::from_two_points(1.0e3, 0.9 * s.ultimate_strength, 1.0e6, s.endurance_limit)
            .map_err(|e| e.to_string())?
            .with_endurance_limit(s.endurance_limit)
            .map_err(|e| e.to_string())?;

    let life = curve
        .cycles_to_failure(sigma_ar)
        .map_err(|e| e.to_string())?;
    let life_str = match life {
        Life::Infinite => "infinite (σ_ar ≤ Se)".to_string(),
        Life::Finite(cycles) => format!("{cycles:.3e} cycles"),
    };
    // The bare power-law life, ignoring the endurance cap, for context.
    let life_unbounded = curve
        .cycles_to_failure_unbounded(sigma_ar)
        .map_err(|e| e.to_string())?;

    let verdict = if n >= 1.0 { "SAFE" } else { "FAIL" };
    let crit_name = match crit {
        MeanStressCriterion::Goodman => "Goodman",
        MeanStressCriterion::Soderberg => "Soderberg",
        MeanStressCriterion::Gerber => "Gerber",
    };

    Ok(format!(
        "criterion       : {crit_name}\n\
         Se / Sy / Su    : {se:.0} / {sy:.0} / {su:.0} MPa\n\
         σa / σm         : {sa:.1} / {sm:.1} MPa\n\n\
         factor of safety: {n:.2}  ({verdict})\n\
         equiv reversed  : {sigma_ar:.1} MPa\n\
         allow σa @ σm   : {sigma_a_allow:.1} MPa  (margin {sigma_a_margin:+.1})\n\
         S-N life        : {life_str}\n\
         (power-law life): {life_unbounded:.3e} cycles",
        se = s.endurance_limit,
        sy = s.yield_strength,
        su = s.ultimate_strength,
        sa = s.sigma_a,
        sm = s.sigma_m,
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

/// Append an axis-aligned cylinder along the `x` axis (centre `c`, radius
/// `r`, half-length `hx`, `seg` facets) to the buffers as a closed
/// triangle tube with end caps.
fn push_x_cylinder(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    c: Vector3<f64>,
    r: f64,
    hx: f64,
    seg: usize,
) {
    let base = nodes.len();
    // Two rings of `seg` vertices at -hx and +hx.
    for ring in 0..2 {
        let x = c.x + if ring == 0 { -hx } else { hx };
        for k in 0..seg {
            let theta = TAU * (k as f64) / (seg as f64);
            nodes.push(Vector3::new(
                x,
                c.y + r * theta.cos(),
                c.z + r * theta.sin(),
            ));
        }
    }
    // Side wall: quad between ring0[k],ring0[k+1],ring1[k+1],ring1[k].
    for k in 0..seg {
        let kn = (k + 1) % seg;
        let a = base + k;
        let b = base + kn;
        let cc = base + seg + kn;
        let d = base + seg + k;
        tris.extend_from_slice(&[a, b, cc, a, cc, d]);
    }
    // End-cap centres.
    let c0 = nodes.len();
    nodes.push(Vector3::new(c.x - hx, c.y, c.z));
    let c1 = nodes.len();
    nodes.push(Vector3::new(c.x + hx, c.y, c.z));
    for k in 0..seg {
        let kn = (k + 1) % seg;
        // -x cap (faces outward in -x).
        tris.extend_from_slice(&[c0, base + kn, base + k]);
        // +x cap (faces outward in +x).
        tris.extend_from_slice(&[c1, base + seg + k, base + seg + kn]);
    }
}

/// Build the fatigue test specimen as a triangle [`Mesh`] — a round
/// dog-bone coupon laid along the `x` axis: a wide grip cylinder at each
/// end and a narrow gauge cylinder in the middle (where fatigue cracks
/// initiate), on a base. Representative geometry (not to scale; the
/// fatigue numbers are the `valenx-fatigue` result). `None` for an
/// invalid material.
fn specimen_solid_mesh(s: &FatigueWorkbenchState) -> Option<Mesh> {
    material(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    let seg = 24;
    // Narrow gauge section in the middle.
    push_x_cylinder(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.4),
        0.12,
        0.45,
        seg,
    );
    // Wide grip at -x end.
    push_x_cylinder(
        &mut nodes,
        &mut tris,
        Vector3::new(-0.65, 0.0, 0.4),
        0.25,
        0.25,
        seg,
    );
    // Wide grip at +x end.
    push_x_cylinder(
        &mut nodes,
        &mut tris,
        Vector3::new(0.65, 0.0, 0.4),
        0.25,
        0.25,
        seg,
    );
    // Base.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.04),
        Vector3::new(0.95, 0.3, 0.04),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-fatigue");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D specimen solid and load it into the central viewport.
fn load_specimen_3d(app: &mut ValenxApp) {
    let Some(mesh) = specimen_solid_mesh(&app.fatigue) else {
        app.fatigue.error =
            Some("material parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<specimen>/valenx-fatigue"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// Agent-bridge product: the canonical fatigue workbench as a 3-D solid plus its
/// `compute()` readout rows (see [`crate::products_registry`]).
pub(crate) fn fatigue_product() -> crate::WorkspaceProduct {
    let s = FatigueWorkbenchState::default();
    let mesh = specimen_solid_mesh(&s).expect("canonical fatigue ⇒ specimen solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<fatigue>/valenx-specimen");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical fatigue ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Fatigue life (S-N curve)".into(),
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
        let s = FatigueWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_safety_and_life() {
        let mut s = FatigueWorkbenchState::default();
        run_fatigue(&mut s);
        assert!(
            s.error.is_none(),
            "default point should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("factor of safety"));
        assert!(s.result.contains("equiv reversed"));
        assert!(s.result.contains("allow σa @ σm"));
        assert!(s.result.contains("S-N life"));
        // Goodman n = 1/(120/200 + 80/500) = 1/0.76 ~ 1.32, inside -> SAFE.
        assert!(s.result.contains("1.32"));
        assert!(s.result.contains("SAFE"));
        // Equivalent reversed sigma_ar = 120/(1 - 80/500) = 142.857 ~ 142.9.
        assert!(s.result.contains("142.9"));
        // Allowable alternating at sm=80 on the Goodman n=1 line:
        // sa = Se*(1 - sm/Su) = 200*(1 - 80/500) = 200*0.84 = 168.0, with
        // the operating margin 168.0 - 120.0 = +48.0 MPa of headroom.
        assert!(s.result.contains("168.0 MPa  (margin +48.0)"));
        // 142.9 MPa is below Se = 200 -> the capped S-N curve is infinite.
        assert!(s.result.contains("infinite"));
    }

    #[test]
    fn analyze_rejects_yield_above_ultimate() {
        let mut s = FatigueWorkbenchState {
            yield_strength: 600.0,
            ultimate_strength: 500.0,
            ..Default::default()
        };
        run_fatigue(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn goodman_factor_of_safety_matches_hand_value() {
        // Ground truth (Shigley): the Goodman factor of safety satisfies
        // 1/n = sigma_a/Se + sigma_m/Su exactly. For Se=200, Su=500,
        // sigma_a=120, sigma_m=80: 1/n = 0.6 + 0.16 = 0.76, n = 1.315789...
        let mat = Material::new(200.0, 350.0, 500.0).unwrap();
        let n = mat
            .factor_of_safety(MeanStressCriterion::Goodman, 120.0, 80.0)
            .unwrap();
        let expected = 1.0 / (120.0 / 200.0 + 80.0 / 500.0);
        assert!((n - expected).abs() < 1e-9);
        assert!((n - 1.3157894736842106).abs() < 1e-9);
    }

    #[test]
    fn allowable_alternating_matches_hand_value_and_readout() {
        // Ground truth (Haigh diagram, Goodman n = 1): the allowable
        // alternating stress at a mean stress sm is sa = Se*(1 - sm/Su).
        // For Se=200, Su=500, sm=80: sa = 200*(1 - 80/500) = 200*0.84 = 168.
        let mat = Material::new(200.0, 350.0, 500.0).unwrap();
        let allow = mat
            .allowable_alternating(MeanStressCriterion::Goodman, 80.0, 1.0)
            .unwrap();
        let expected = 200.0 * (1.0 - 80.0 / 500.0);
        assert!((allow - expected).abs() < 1e-9);
        assert!((allow - 168.0).abs() < 1e-9);
        // And the operating margin below it: 168 - 120 = 48 MPa of headroom.
        let margin = allow - 120.0;
        assert!((margin - 48.0).abs() < 1e-9);
        // The default workbench point shares these numbers; the readout
        // must format them at .1 precision with a signed margin.
        let s = FatigueWorkbenchState::default();
        let r = compute(&s).expect("default point analyzes");
        assert!(r.contains("allow σa @ σm   : 168.0 MPa  (margin +48.0)"));
    }

    #[test]
    fn specimen_mesh_for_default_is_nonempty_and_in_range() {
        let s = FatigueWorkbenchState::default();
        let mesh = specimen_solid_mesh(&s).expect("default material yields a solid");
        assert!(mesh.nodes.len() > 8, "expected gauge + two grips + base");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn specimen_mesh_none_for_invalid() {
        let s = FatigueWorkbenchState {
            yield_strength: 600.0,
            ultimate_strength: 500.0,
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
            draw_fatigue_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_fatigue_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_fatigue_workbench = true;
        run_fatigue(&mut app.fatigue);
        draw_workbench(&mut app);
    }
}
