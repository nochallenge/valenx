//! The right-side **Bearing Workbench** panel — native rolling-element
//! bearing basic rating-life analysis over `valenx-bearing`.
//!
//! Mirrors the Heat Transfer / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_bearing_workbench`,
//! toggled from the View menu. The form sets a bearing's basic dynamic load
//! rating `C`, a combined radial/axial load with its `X` / `Y` factors, the
//! rolling-element type and the shaft speed; "Analyze" forms the dynamic
//! equivalent load `P = X·Fr + Y·Fa`, evaluates the ISO 281 basic rating
//! life `L10 = (C / P)^p` and converts it to operating hours
//! `L10h = L10 · 1e6 / (60 · n)`. Alongside the dynamic fatigue life it
//! also reports the ISO 76 *static* check — the static equivalent load
//! `P0 = max(X0·Fr + Y0·Fa, Fr)` and the static safety factor
//! `s0 = C0 / P0`, which guards a slow or stationary bearing against
//! brinelling. "Show 3-D bearing" loads a representative bearing solid
//! (outer ring, inner ring and rolling elements) into the central
//! viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_bearing::{BearingType, EquivalentLoad, RatingLife, StaticEquivalentLoad};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Bearing Workbench.
pub struct BearingWorkbenchState {
    /// Basic dynamic load rating `C` (newtons).
    dynamic_load_rating_n: f64,
    /// Radial load component `Fr` (newtons).
    radial_n: f64,
    /// Axial (thrust) load component `Fa` (newtons).
    axial_n: f64,
    /// Dimensionless radial load factor `X`.
    x_factor: f64,
    /// Dimensionless axial load factor `Y`.
    y_factor: f64,
    /// Basic static load rating `C0` (newtons), for the ISO 76 static
    /// safety factor `s0 = C0 / P0`.
    static_load_rating_n: f64,
    /// Dimensionless static radial load factor `X0`.
    x0_factor: f64,
    /// Dimensionless static axial load factor `Y0`.
    y0_factor: f64,
    /// Rolling-element type, which fixes the load-life exponent `p`.
    bearing_type: BearingType,
    /// Shaft speed `n` (revolutions per minute).
    rpm: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D bearing solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for BearingWorkbenchState {
    fn default() -> Self {
        // A deep-groove ball bearing: C = 25 kN dynamic rating carrying
        // Fr = 5 kN radial with no thrust (radial-only: X = 1, Y = 0), so
        // P = 5 kN and L10 = (25/5)^3 = 125 Mrev; at 1500 rpm that is
        // ~1388.9 h.
        //
        // Static side (ISO 76): C0 = 15 kN basic static rating with the
        // usual ball factors X0 = 0.6, Y0 = 0.5. With Fa = 0 the formula
        // gives X0*Fr = 3 kN, below Fr, so the ISO 76 floor takes P0 = Fr
        // = 5 kN and the static safety factor is s0 = 15/5 = 3.0.
        Self {
            dynamic_load_rating_n: 25_000.0,
            radial_n: 5_000.0,
            axial_n: 0.0,
            x_factor: 1.0,
            y_factor: 0.0,
            static_load_rating_n: 15_000.0,
            x0_factor: 0.6,
            y0_factor: 0.5,
            bearing_type: BearingType::Ball,
            rpm: 1500.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Bearing Workbench right-side panel. A no-op when the
/// `show_bearing_workbench` toggle is off.
pub fn draw_bearing_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_bearing_workbench {
        return;
    }

    egui::SidePanel::right("valenx_bearing_workbench")
        .resizable(true)
        .default_width(360.0)
        .width_range(300.0..=560.0)
        .show(ctx, |ui| {
            if crate::workbench_ui::header(
                ui,
                "Bearing",
                "native ISO 281 basic rating-life L10 · valenx-bearing",
            ) {
                app.show_bearing_workbench = false;
            }

            let s = &mut app.bearing;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Bearing").strong());
                    ui.horizontal(|ui| {
                        ui.label("dynamic rating C (N)");
                        ui.add(
                            egui::DragValue::new(&mut s.dynamic_load_rating_n).speed(100.0),
                        );
                    });
                    ui.horizontal(|ui| {
                        ui.label("element type");
                        ui.selectable_value(&mut s.bearing_type, BearingType::Ball, "ball (p=3)");
                        ui.selectable_value(
                            &mut s.bearing_type,
                            BearingType::Roller,
                            "roller (p=10/3)",
                        );
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Load").strong());
                    ui.horizontal(|ui| {
                        ui.label("radial Fr (N)");
                        ui.add(egui::DragValue::new(&mut s.radial_n).speed(50.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("axial Fa (N)");
                        ui.add(egui::DragValue::new(&mut s.axial_n).speed(50.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("factor X");
                        ui.add(egui::DragValue::new(&mut s.x_factor).speed(0.01));
                    });
                    ui.horizontal(|ui| {
                        ui.label("factor Y");
                        ui.add(egui::DragValue::new(&mut s.y_factor).speed(0.01));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Speed").strong());
                    ui.horizontal(|ui| {
                        ui.label("shaft speed (rpm)");
                        ui.add(egui::DragValue::new(&mut s.rpm).speed(10.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Static (ISO 76)").strong());
                    ui.horizontal(|ui| {
                        ui.label("static rating C0 (N)");
                        ui.add(
                            egui::DragValue::new(&mut s.static_load_rating_n).speed(100.0),
                        );
                    });
                    ui.horizontal(|ui| {
                        ui.label("factor X0");
                        ui.add(egui::DragValue::new(&mut s.x0_factor).speed(0.01));
                    });
                    ui.horizontal(|ui| {
                        ui.label("factor Y0");
                        ui.add(egui::DragValue::new(&mut s.y0_factor).speed(0.01));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_bearing(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D bearing").strong())
                        .on_hover_text(
                            "Build a representative rolling-element bearing (outer ring, inner ring and balls) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Rating life").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        });

    // Serviced after the panel draws (the `&mut app.bearing` borrow is
    // released here): build the bearing's 3-D solid and load it.
    if app.bearing.show_3d_request {
        app.bearing.show_3d_request = false;
        load_bearing_3d(app);
    }
}

/// Validate the form, evaluate the bearing and format the readout.
fn run_bearing(s: &mut BearingWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build the validated [`RatingLife`] from the form, the object both the
/// readout and the 3-D gate need. Extracted so it is unit-testable and
/// shared.
fn rating_life(s: &BearingWorkbenchState) -> Result<RatingLife, String> {
    let load = EquivalentLoad::new(s.radial_n, s.axial_n, s.x_factor, s.y_factor)
        .map_err(|e| e.to_string())?;
    RatingLife::from_equivalent_load(s.dynamic_load_rating_n, &load, s.bearing_type)
        .map_err(|e| e.to_string())
}

/// Evaluate the bearing and format the full readout, mapping any domain
/// error to a display string. Extracted so it is unit-testable.
fn compute(s: &BearingWorkbenchState) -> Result<String, String> {
    let load = EquivalentLoad::new(s.radial_n, s.axial_n, s.x_factor, s.y_factor)
        .map_err(|e| e.to_string())?;
    let p = load.value();
    let life = RatingLife::from_equivalent_load(s.dynamic_load_rating_n, &load, s.bearing_type)
        .map_err(|e| e.to_string())?;
    let l10 = life.l10_million_revs();
    let hours = life.life_hours(s.rpm).map_err(|e| e.to_string())?;
    let exponent = s.bearing_type.life_exponent();
    let type_name = match s.bearing_type {
        BearingType::Ball => "ball",
        BearingType::Roller => "roller",
    };

    // ISO 76 static check: the static equivalent load
    // P0 = max(X0*Fr + Y0*Fa, Fr) and the static safety factor
    // s0 = C0 / P0, which guards a slow / stationary bearing against
    // brinelling rather than rolling-contact fatigue.
    let static_load = StaticEquivalentLoad::new(s.radial_n, s.axial_n, s.x0_factor, s.y0_factor)
        .map_err(|e| e.to_string())?;
    let p0 = static_load.value();
    let s0 = static_load
        .safety_factor(s.static_load_rating_n)
        .map_err(|e| e.to_string())?;

    Ok(format!(
        "dynamic rating C: {:.0} N\n\
         element type    : {} (p = {:.3})\n\
         radial Fr       : {:.0} N\n\
         axial Fa        : {:.0} N\n\
         factors X / Y   : {:.2} / {:.2}\n\
         shaft speed n   : {:.0} rpm\n\n\
         equiv load P    : {:.0} N\n\
         C / P ratio     : {:.3}\n\
         L10 (basic life): {:.1} Mrev\n\
         L10h (hours)    : {:.0} h\n\n\
         static rating C0: {:.0} N\n\
         factors X0 / Y0 : {:.2} / {:.2}\n\
         static load P0  : {:.0} N\n\
         static safety s0: {:.2}",
        s.dynamic_load_rating_n,
        type_name,
        exponent,
        s.radial_n,
        s.axial_n,
        s.x_factor,
        s.y_factor,
        s.rpm,
        p,
        s.dynamic_load_rating_n / p,
        l10,
        hours,
        s.static_load_rating_n,
        s.x0_factor,
        s.y0_factor,
        p0,
        s0,
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

/// Append a (double-sided) cylinder whose axis runs along `+x`, spanning
/// `base.x ..= base.x + length` with circle centre `(base.y, base.z)`.
fn push_cyl_x(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    base: Vector3<f64>,
    length: f64,
    r: f64,
    seg: usize,
) {
    let (x0, x1) = (base.x, base.x + length);
    let lo = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(x0, base.y + r * a.cos(), base.z + r * a.sin()));
    }
    let hi = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(x1, base.y + r * a.cos(), base.z + r * a.sin()));
    }
    for j in 0..seg {
        let jn = (j + 1) % seg;
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
    }
}

/// Build the bearing as a triangle [`Mesh`] — an outer ring and an inner
/// ring (two concentric short cylinders along `x`) with a ring of rolling
/// elements between them and a base. Representative geometry (not to scale;
/// the rating-life numbers are the `valenx-bearing` result). `None` for an
/// invalid configuration.
fn bearing_solid_mesh(s: &BearingWorkbenchState) -> Option<Mesh> {
    rating_life(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    let z = 0.6;
    let width = 0.18;
    let outer_r = 0.5;
    let inner_r = 0.28;
    let ball_r = 0.09;
    let pitch_r = (outer_r + inner_r) * 0.5;

    // Outer ring (large concentric x-cylinder).
    push_cyl_x(
        &mut nodes,
        &mut tris,
        Vector3::new(-width * 0.5, 0.0, z),
        width,
        outer_r,
        36,
    );
    // Inner ring (small concentric x-cylinder, sharing the axis).
    push_cyl_x(
        &mut nodes,
        &mut tris,
        Vector3::new(-width * 0.5, 0.0, z),
        width,
        inner_r,
        28,
    );
    // A ring of rolling elements (balls) on the pitch circle, drawn as
    // small boxes spaced around the axis.
    let balls = 8;
    for j in 0..balls {
        let a = j as f64 / balls as f64 * TAU;
        let cy = pitch_r * a.cos();
        let cz = z + pitch_r * a.sin();
        push_box(
            &mut nodes,
            &mut tris,
            Vector3::new(0.0, cy, cz),
            Vector3::new(width * 0.4, ball_r, ball_r),
        );
    }
    // Base.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.05),
        Vector3::new(0.25, outer_r, 0.05),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-bearing");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D bearing solid and load it into the central viewport.
fn load_bearing_3d(app: &mut ValenxApp) {
    let Some(mesh) = bearing_solid_mesh(&app.bearing) else {
        app.bearing.error =
            Some("bearing parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<bearing>/valenx-bearing"),
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
        let s = BearingWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_l10_and_hours() {
        let mut s = BearingWorkbenchState::default();
        run_bearing(&mut s);
        assert!(
            s.error.is_none(),
            "default bearing should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("L10 (basic life)"));
        assert!(s.result.contains("L10h (hours)"));
        // C = 25 kN, P = 5 kN, ball: L10 = 5^3 = 125 Mrev.
        assert!(s.result.contains("125.0 Mrev"));
        // 125e6 / (60 * 1500) = 1388.9 h -> formats to "1389 h".
        assert!(s.result.contains("1389 h"));
    }

    #[test]
    fn analyze_rejects_zero_load() {
        // Zero load components with the default X = 1, Y = 0 give P = 0,
        // which the life formula cannot evaluate.
        let mut s = BearingWorkbenchState {
            radial_n: 0.0,
            axial_n: 0.0,
            ..Default::default()
        };
        run_bearing(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn analyze_reports_iso76_static_safety_factor() {
        // Ground truth (defaults): static side is Fr = 5000 N, Fa = 0,
        // X0 = 0.6, Y0 = 0.5, C0 = 15000 N. The ISO 76 static equivalent
        // load is P0 = max(0.6*5000 + 0.5*0, 5000) = max(3000, 5000)
        // = 5000 N (the formula falls below Fr, so the floor applies), and
        // the static safety factor is s0 = C0 / P0 = 15000 / 5000 = 3.00.
        let mut s = BearingWorkbenchState::default();
        run_bearing(&mut s);
        assert!(
            s.error.is_none(),
            "default bearing should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("static load P0  : 5000 N"));
        assert!(s.result.contains("static safety s0: 3.00"));
    }

    #[test]
    fn static_equivalent_load_uses_formula_above_floor() {
        // Ground truth: with a real thrust the formula clears Fr. Fr =
        // 2000 N, Fa = 10000 N, X0 = 0.6, Y0 = 0.5 give P0 = 0.6*2000 +
        // 0.5*10000 = 1200 + 5000 = 6200 N (> Fr = 2000, so used directly),
        // and with C0 = 31000 N the static safety factor is
        // s0 = 31000 / 6200 = 5.00.
        let s = BearingWorkbenchState {
            radial_n: 2_000.0,
            axial_n: 10_000.0,
            x0_factor: 0.6,
            y0_factor: 0.5,
            static_load_rating_n: 31_000.0,
            ..Default::default()
        };
        let out = compute(&s).expect("valid bearing");
        assert!(out.contains("static load P0  : 6200 N"));
        assert!(out.contains("static safety s0: 5.00"));
    }

    #[test]
    fn l10_is_c_over_p_cubed_for_ball() {
        // Ground truth: for a ball bearing (p = 3) the basic rating life
        // is exactly (C / P)^3. C = 30 kN, radial-only P = 6 kN -> ratio 5,
        // L10 = 5^3 = 125 Mrev.
        let s = BearingWorkbenchState {
            dynamic_load_rating_n: 30_000.0,
            radial_n: 6_000.0,
            axial_n: 0.0,
            x_factor: 1.0,
            y_factor: 0.0,
            bearing_type: BearingType::Ball,
            ..Default::default()
        };
        let life = rating_life(&s).expect("valid bearing");
        let ratio = 30_000.0_f64 / 6_000.0;
        assert!((life.l10_million_revs() - ratio.powi(3)).abs() < 1e-9);
        assert!((life.l10_million_revs() - 125.0).abs() < 1e-9);
    }

    #[test]
    fn bearing_mesh_for_default_is_nonempty_and_in_range() {
        let s = BearingWorkbenchState::default();
        let mesh = bearing_solid_mesh(&s).expect("default bearing yields a solid");
        assert!(mesh.nodes.len() > 8, "expected rings + balls + base");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn bearing_mesh_none_for_invalid() {
        let s = BearingWorkbenchState {
            radial_n: 0.0,
            axial_n: 0.0,
            ..Default::default()
        };
        assert!(bearing_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_bearing_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_bearing_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_bearing_workbench = true;
        run_bearing(&mut app.bearing);
        draw_workbench(&mut app);
    }
}
