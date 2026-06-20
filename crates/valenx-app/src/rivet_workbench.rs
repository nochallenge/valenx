//! The right-side **Riveted Joint Workbench** panel — native closed-form
//! rivet-joint strength analysis over `valenx-rivet`.
//!
//! Mirrors the Heat Transfer / Torsion workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_rivet_workbench`,
//! toggled from the View menu. The form sets a rivet group (diameter, how
//! many rivets per row, rows, shear planes), the plate (gross width and
//! thickness) and the permissible shear / bearing / tensile stresses;
//! "Analyze" evaluates the three classic failure modes — rivet shear,
//! plate bearing and net-section tension — and reports each failure load,
//! the governing strength and mode, the solid (un-drilled) plate strength
//! that anchors the efficiency, the joint efficiency, and the factor of
//! safety / utilization against an applied service load. "Show 3-D"
//! loads a representative riveted-lap-joint solid (two overlapping plates
//! pierced by rivets) into the central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_rivet::{Allowables, FailureMode, Joint, Plate, RivetGroup};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Riveted Joint Workbench.
pub struct RivetWorkbenchState {
    /// Rivet (shank / hole) diameter `d` (m).
    diameter_m: f64,
    /// Number of rivets in the critical (tension-carrying) row, `n_row`.
    rivets_per_row: u32,
    /// Number of rivet rows, `n_rows`.
    rows: u32,
    /// Shear planes per rivet, `s` (1 = lap / single-cover butt, 2 =
    /// double-cover butt).
    shear_planes: u32,
    /// Gross plate width at the critical row, `w` (m).
    width_m: f64,
    /// Plate thickness, `t` (m).
    thickness_m: f64,
    /// Permissible shear stress in the rivet shank, `τ` (MPa).
    shear_mpa: f64,
    /// Permissible bearing (crushing) stress, `σ_b` (MPa).
    bearing_mpa: f64,
    /// Permissible tensile stress in the plate, `σ_t` (MPa).
    tension_mpa: f64,
    /// Applied service load to rate the joint against, `P_applied` (kN).
    applied_kn: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D joint solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for RivetWorkbenchState {
    fn default() -> Self {
        // The textbook lap joint: one row of three 20 mm rivets in single
        // shear, joining 10 mm plates 150 mm wide, with mild-steel
        // working stresses (τ = 80, σ_b = 160, σ_t = 100 MPa). Rivet shear
        // governs here at ~75.4 kN (shear 75.4 < tension 90 < bearing 96),
        // an efficiency of ~50 %. Rated against a 40 kN service load.
        Self {
            diameter_m: 0.020,
            rivets_per_row: 3,
            rows: 1,
            shear_planes: 1,
            width_m: 0.150,
            thickness_m: 0.010,
            shear_mpa: 80.0,
            bearing_mpa: 160.0,
            tension_mpa: 100.0,
            applied_kn: 40.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Riveted Joint Workbench right-side panel. A no-op when the
/// `show_rivet_workbench` toggle is off.
pub fn draw_rivet_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_rivet_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(app, ctx, "valenx_rivet_workbench", "Riveted Joint", |app, ui| {
            ui.label(egui::RichText::new("native closed-form rivet-joint strength · valenx-rivet").weak().small());
            ui.separator();

            let s = &mut app.rivet;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Rivet group").strong());
                    ui.horizontal(|ui| {
                        ui.label("diameter d (m)");
                        ui.add(egui::DragValue::new(&mut s.diameter_m).speed(0.001));
                    });
                    ui.horizontal(|ui| {
                        ui.label("rivets per row");
                        ui.add(egui::DragValue::new(&mut s.rivets_per_row).speed(1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("rows");
                        ui.add(egui::DragValue::new(&mut s.rows).speed(1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("shear planes s");
                        ui.radio_value(&mut s.shear_planes, 1, "1 (lap)");
                        ui.radio_value(&mut s.shear_planes, 2, "2 (butt)");
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Plate").strong());
                    ui.horizontal(|ui| {
                        ui.label("width w (m)");
                        ui.add(egui::DragValue::new(&mut s.width_m).speed(0.005));
                    });
                    ui.horizontal(|ui| {
                        ui.label("thickness t (m)");
                        ui.add(egui::DragValue::new(&mut s.thickness_m).speed(0.001));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Allowable stresses").strong());
                    ui.horizontal(|ui| {
                        ui.label("shear τ (MPa)");
                        ui.add(egui::DragValue::new(&mut s.shear_mpa).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("bearing σ_b (MPa)");
                        ui.add(egui::DragValue::new(&mut s.bearing_mpa).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("tension σ_t (MPa)");
                        ui.add(egui::DragValue::new(&mut s.tension_mpa).speed(1.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Applied load").strong());
                    ui.horizontal(|ui| {
                        ui.label("service load P (kN)");
                        ui.add(egui::DragValue::new(&mut s.applied_kn).speed(1.0));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_rivet(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative riveted lap joint (two overlapping plates pierced by a row of rivets) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Joint strength").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        }, );
    if close { app.show_rivet_workbench = false; }

    // Serviced after the panel draws (the `&mut app.rivet` borrow is
    // released here): build the joint's 3-D solid and load it.
    if app.rivet.show_3d_request {
        app.rivet.show_3d_request = false;
        load_joint_3d(app);
    }
}

/// Validate the form, evaluate the joint and format the readout.
fn run_rivet(s: &mut RivetWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build the validated [`Joint`] from the form, mapping any domain error
/// to a display string. Extracted so it is shared by the readout and the
/// 3-D gate. Stresses enter in MPa and are converted to pascals.
fn build_joint(s: &RivetWorkbenchState) -> Result<Joint, String> {
    let group = RivetGroup::new(s.diameter_m, s.rivets_per_row, s.rows, s.shear_planes)
        .map_err(|e| e.to_string())?;
    let plate = Plate::new(s.width_m, s.thickness_m).map_err(|e| e.to_string())?;
    let allow = Allowables::new(
        s.shear_mpa * 1.0e6,
        s.bearing_mpa * 1.0e6,
        s.tension_mpa * 1.0e6,
    )
    .map_err(|e| e.to_string())?;
    Ok(Joint::new(group, plate, allow))
}

/// Evaluate the joint and format the full readout, mapping any domain
/// error to a display string. Extracted so it is unit-testable.
fn compute(s: &RivetWorkbenchState) -> Result<String, String> {
    let joint = build_joint(s)?;
    let r = joint.analyze().map_err(|e| e.to_string())?;
    let applied = s.applied_kn * 1.0e3;
    let fos = joint.factor_of_safety(applied).map_err(|e| e.to_string())?;
    let util = joint.utilization(applied).map_err(|e| e.to_string())?;

    let mode = match r.mode {
        FailureMode::Shear => "rivet shear",
        FailureMode::Bearing => "plate bearing",
        FailureMode::Tension => "net-section tension",
    };

    Ok(format!(
        "rivets          : {n} ({n_row}/row x {rows} row(s))\n\
         diameter d      : {d:.3} m, {s_pl} shear plane(s)\n\
         plate w x t     : {w:.3} x {t:.3} m\n\
         allowables      : tau {tau:.0} / sb {sb:.0} / st {st:.0} MPa\n\n\
         P shear         : {shear:.2} kN\n\
         P bearing       : {bearing:.2} kN\n\
         P tension       : {tension:.2} kN\n\
         governing       : {mode} at {strength:.2} kN\n\
         solid plate     : {solid:.2} kN\n\
         efficiency      : {eff:.1} %\n\n\
         applied load    : {applied:.2} kN\n\
         factor of safety: {fos:.2}\n\
         utilization     : {util:.1} %",
        n = joint.group.total_rivets(),
        n_row = s.rivets_per_row,
        rows = s.rows,
        d = s.diameter_m,
        s_pl = s.shear_planes,
        w = s.width_m,
        t = s.thickness_m,
        tau = s.shear_mpa,
        sb = s.bearing_mpa,
        st = s.tension_mpa,
        shear = r.shear / 1.0e3,
        bearing = r.bearing / 1.0e3,
        tension = r.tension / 1.0e3,
        strength = r.strength / 1.0e3,
        solid = r.solid / 1.0e3,
        eff = r.efficiency * 100.0,
        applied = s.applied_kn,
        util = util * 100.0,
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

/// Append a rivet shank: a closed cylinder of radius `r` swept along the
/// through-thickness `z` axis from `z0` to `z1`, centred at `(cx, cy)`,
/// with `seg` angular segments (both end caps and the wall). Represents a
/// single rivet piercing the overlapped plates.
#[allow(clippy::too_many_arguments)]
fn push_cylinder(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    cx: f64,
    cy: f64,
    z0: f64,
    z1: f64,
    r: f64,
    seg: usize,
) {
    let base = nodes.len();
    // Per angular step i: a top (z1) and a bottom (z0) node on the wall.
    for i in 0..seg {
        let a = std::f64::consts::TAU * (i as f64) / (seg as f64);
        let (sin, cos) = a.sin_cos();
        nodes.push(Vector3::new(cx + r * cos, cy + r * sin, z1));
        nodes.push(Vector3::new(cx + r * cos, cy + r * sin, z0));
    }
    // Two centre nodes for the end-cap fans.
    let top_c = nodes.len();
    nodes.push(Vector3::new(cx, cy, z1));
    let bot_c = nodes.len();
    nodes.push(Vector3::new(cx, cy, z0));

    for i in 0..seg {
        let j = (i + 1) % seg;
        let (t, b) = (2 * i, 2 * i + 1);
        let (nt, nb) = (2 * j, 2 * j + 1);
        // Outer wall quad.
        tris.extend_from_slice(&[
            base + t,
            base + b,
            base + nb,
            base + t,
            base + nb,
            base + nt,
        ]);
        // Top cap fan and bottom cap fan.
        tris.extend_from_slice(&[top_c, base + t, base + nt]);
        tris.extend_from_slice(&[bot_c, base + nb, base + b]);
    }
}

/// Build the riveted lap joint as a triangle [`Mesh`] — two overlapping
/// plates (offset in the through-thickness `z` direction) pierced by a row
/// of rivet cylinders. Representative geometry (not to scale; the strength
/// numbers are the `valenx-rivet` result). `None` for an invalid
/// configuration.
fn joint_solid_mesh(s: &RivetWorkbenchState) -> Option<Mesh> {
    build_joint(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Plate footprint: load axis along x, width along y, thin in z.
    let half_x = 0.6;
    let half_y = 0.35;
    let half_t = 0.04;

    // Lower plate (z around -half_t) and upper plate (z around +half_t),
    // overlapping in the central x region — a lap joint. The lower plate
    // reaches in from -x, the upper from +x; they overlap near x = 0.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(-0.3, 0.0, -half_t),
        Vector3::new(half_x, half_y, half_t),
    );
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.3, 0.0, half_t),
        Vector3::new(half_x, half_y, half_t),
    );

    // A row of rivets through the overlap (the critical tension row),
    // spanning both plate thicknesses. Lay `n_row` rivets across y.
    let n = s.rivets_per_row.max(1) as usize;
    let r_rivet = 0.05;
    let z_lo = -2.0 * half_t;
    let z_hi = 2.0 * half_t;
    for i in 0..n {
        // Evenly spread across the width, inset from the edges.
        let frac = if n == 1 {
            0.5
        } else {
            (i as f64) / ((n - 1) as f64)
        };
        let cy = -0.22 + frac * 0.44;
        push_cylinder(&mut nodes, &mut tris, 0.0, cy, z_lo, z_hi, r_rivet, 32);
    }

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-rivet");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D joint solid and load it into the central viewport.
fn load_joint_3d(app: &mut ValenxApp) {
    let Some(mesh) = joint_solid_mesh(&app.rivet) else {
        app.rivet.error = Some("joint parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<joint>/valenx-rivet"),
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
        let s = RivetWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_modes_and_governing() {
        let mut s = RivetWorkbenchState::default();
        run_rivet(&mut s);
        assert!(
            s.error.is_none(),
            "default joint should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("P shear"));
        assert!(s.result.contains("P bearing"));
        assert!(s.result.contains("P tension"));
        assert!(s.result.contains("governing"));
        assert!(s.result.contains("factor of safety"));
        // The textbook lap joint (three 20 mm single-shear rivets, 150x10
        // mm plate, mild-steel allowables) is rivet-shear governed at
        // 75.40 kN, an efficiency of 50.3 %.
        assert!(s.result.contains("rivet shear"));
        assert!(s.result.contains("75.40 kN"));
        assert!(s.result.contains("50.3 %"));
    }

    #[test]
    fn solid_plate_strength_matches_hand_computed_value() {
        // Ground truth: the solid (un-drilled) plate strength is
        // `P_solid = w · t · σ_t`, the denominator of the joint
        // efficiency. For the default textbook joint that is
        // 0.150 m · 0.010 m · 100e6 Pa = 150 000 N = 150.00 kN, and the
        // reported efficiency must equal governing / solid exactly.
        let s = RivetWorkbenchState::default();
        let joint = build_joint(&s).expect("valid joint");
        let r = joint.analyze().expect("analyzable");
        let expected = 0.150_f64 * 0.010_f64 * 100.0e6;
        assert!(
            (r.solid - expected).abs() < 1.0e-6,
            "solid-plate strength {} should equal w*t*st {expected} N",
            r.solid
        );
        // Efficiency is governing / solid by construction.
        assert!(
            (r.efficiency - r.strength / r.solid).abs() < 1.0e-12,
            "efficiency {} should equal strength/solid {}",
            r.efficiency,
            r.strength / r.solid
        );
        // And it surfaces in the formatted readout as 150.00 kN.
        let out = compute(&s).expect("default joint formats");
        assert!(
            out.contains("solid plate     : 150.00 kN"),
            "readout should report the solid-plate strength: {out}"
        );
    }

    #[test]
    fn analyze_rejects_zero_diameter() {
        let mut s = RivetWorkbenchState {
            diameter_m: 0.0,
            ..Default::default()
        };
        run_rivet(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn analyze_rejects_consumed_net_section() {
        // 5 rivets of 30 mm across a 150 mm plate remove the whole width,
        // so the net-section tension strength is undefined.
        let mut s = RivetWorkbenchState {
            diameter_m: 0.030,
            rivets_per_row: 5,
            width_m: 0.150,
            ..Default::default()
        };
        run_rivet(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn single_rivet_shear_matches_hand_computed_strength() {
        // Ground truth: one rivet in single shear carries τ·π·d²/4. With
        // d = 20 mm and τ = 80 MPa that is 80e6 · π · 0.020² / 4 ≈ 25.13 kN,
        // and the governing strength of a one-rivet shear-only joint (a
        // very wide, thick, strong plate so bearing and tension cannot
        // govern) must equal exactly that single-plane shear load.
        let s = RivetWorkbenchState {
            diameter_m: 0.020,
            rivets_per_row: 1,
            rows: 1,
            shear_planes: 1,
            width_m: 1.0,
            thickness_m: 0.05,
            shear_mpa: 80.0,
            bearing_mpa: 1000.0,
            tension_mpa: 1000.0,
            ..Default::default()
        };
        let joint = build_joint(&s).expect("valid joint");
        let r = joint.analyze().expect("analyzable");
        let expected = 80.0e6 * std::f64::consts::PI * 0.020_f64 * 0.020_f64 / 4.0;
        assert_eq!(r.mode, FailureMode::Shear);
        assert!(
            (r.strength - expected).abs() < 1.0e-6,
            "governing {} should equal single-rivet shear {expected}",
            r.strength
        );
    }

    #[test]
    fn governing_strength_is_minimum_of_the_three_modes() {
        // Ground truth: the joint strength is the minimum of the three
        // failure loads. Pin the literal 1.0 before .min to keep the
        // comparison in f64.
        let s = RivetWorkbenchState::default();
        let joint = build_joint(&s).expect("valid joint");
        let r = joint.analyze().expect("analyzable");
        let lo = r.shear.min(r.bearing).min(r.tension);
        assert!(
            (r.strength - lo).abs() < 1.0e-6,
            "strength {} should equal min {lo}",
            r.strength
        );
        assert!(r.efficiency > 0.0 && r.efficiency < 1.0_f64);
    }

    #[test]
    fn joint_mesh_for_default_is_nonempty_and_in_range() {
        let s = RivetWorkbenchState::default();
        let mesh = joint_solid_mesh(&s).expect("default joint yields a solid");
        assert!(mesh.nodes.len() > 16, "expected two plates + rivets");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn joint_mesh_none_for_invalid() {
        let s = RivetWorkbenchState {
            diameter_m: 0.0,
            ..Default::default()
        };
        assert!(joint_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_rivet_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_rivet_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_rivet_workbench = true;
        run_rivet(&mut app.rivet);
        draw_workbench(&mut app);
    }
}
