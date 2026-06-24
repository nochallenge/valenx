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

use valenx_bearing::{BearingType, EquivalentLoad, RatingLife, StaticEquivalentLoad};
use valenx_mesh::Mesh;

use crate::mesh_prims::MeshBuilder;
use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Brushed-steel raceways.
const STEEL: [f32; 3] = [0.60, 0.62, 0.66];
/// Bright chrome rolling balls.
const CHROME: [f32; 3] = [0.82, 0.84, 0.88];
/// Brass cage.
const BRASS: [f32; 3] = [0.72, 0.60, 0.30];
/// Presentation spin rate of the rolling assembly (inner race + balls + cage),
/// rad/s (~0.6 rev/s) — a readable inspect speed for the orbiting elements.
const SPIN_RAD_PER_S: f32 = 4.0;

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

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_bearing_workbench",
        "Bearing",
        |app, ui| {
            ui.label(
                egui::RichText::new("native ISO 281 basic rating-life L10 · valenx-bearing")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.bearing;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Bearing").strong());
                    // Associate each numeric `DragValue` with its caption via
                    // `labelled_by`, so the spin button carries the caption as
                    // its accessibility / UI-Automation Name (egui clears a
                    // DragValue's own Name, leaving it anonymous to a screen
                    // reader / AI driver otherwise).
                    ui.horizontal(|ui| {
                        let l = ui.label("dynamic rating C (N)");
                        ui.add(egui::DragValue::new(&mut s.dynamic_load_rating_n).speed(100.0))
                            .labelled_by(l.id);
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
                        let l = ui.label("radial Fr (N)");
                        ui.add(egui::DragValue::new(&mut s.radial_n).speed(50.0))
                            .labelled_by(l.id);
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("axial Fa (N)");
                        ui.add(egui::DragValue::new(&mut s.axial_n).speed(50.0))
                            .labelled_by(l.id);
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("factor X");
                        ui.add(egui::DragValue::new(&mut s.x_factor).speed(0.01))
                            .labelled_by(l.id);
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("factor Y");
                        ui.add(egui::DragValue::new(&mut s.y_factor).speed(0.01))
                            .labelled_by(l.id);
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Speed").strong());
                    ui.horizontal(|ui| {
                        let l = ui.label("shaft speed (rpm)");
                        ui.add(egui::DragValue::new(&mut s.rpm).speed(10.0))
                            .labelled_by(l.id);
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Static (ISO 76)").strong());
                    ui.horizontal(|ui| {
                        let l = ui.label("static rating C0 (N)");
                        ui.add(egui::DragValue::new(&mut s.static_load_rating_n).speed(100.0))
                            .labelled_by(l.id);
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("factor X0");
                        ui.add(egui::DragValue::new(&mut s.x0_factor).speed(0.01))
                            .labelled_by(l.id);
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("factor Y0");
                        ui.add(egui::DragValue::new(&mut s.y0_factor).speed(0.01))
                            .labelled_by(l.id);
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
        },
    );
    if close {
        app.show_bearing_workbench = false;
    }

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

/// Build the bearing as a triangle [`Mesh`] **with per-vertex colours** plus a
/// [`crate::RigidPart`] for the rolling assembly. A real rolling-element
/// bearing cross-section on the +z axle:
///
/// - **outer race** and **inner race** are true *annular* [`MeshBuilder::tube`]s
///   (hollow rings with a bore, not solid cylinders) in brushed steel;
/// - the **rolling elements** are bright-chrome [`MeshBuilder::sphere`] balls on
///   the pitch circle (the primitive the old box-balls were faking);
/// - a thin brass **cage** ([`MeshBuilder::torus`]) rides the pitch circle
///   linking the balls.
///
/// The inner race + balls + cage are built consecutively as one contiguous node
/// range that spins about the +z axle (a readable inspect of the orbiting
/// elements); the outer race stays fixed. Representative dimensions (the
/// rating-life numbers are the `valenx-bearing` result). `None` for an invalid
/// configuration.
///
/// Returns `(mesh, colors, parts)` with `colors.len() == 3 × triangle_count`.
fn bearing_solid_mesh_parts(
    s: &BearingWorkbenchState,
) -> Option<(Mesh, Vec<[f32; 3]>, Vec<crate::RigidPart>)> {
    rating_life(s).ok()?;

    // Representative bearing cross-section (axle = +z).
    let width = 0.18;
    let outer_od = 0.5; // outer race outer radius
    let outer_id = 0.40; // outer race bore (race wall = outer_od − outer_id)
    let inner_od = 0.30; // inner race outer radius
    let inner_id = 0.20; // inner race bore (the shaft hole)
    let pitch_r = (outer_id + inner_od) * 0.5; // ball-centre circle
    let ball_r = (outer_id - inner_od) * 0.5 * 0.92; // ball nearly fills the gap
    let cage_minor = ball_r * 0.32; // thin cage ring tube

    let mut b = MeshBuilder::new();

    // Outer race — a steel annulus (hollow ring), fixed.
    b.tube(
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        outer_id,
        outer_od,
        width,
        40,
        STEEL,
    );

    // Rolling assembly (inner race + balls + cage) — one contiguous range that
    // spins about +z. Build the inner race first so the range starts here.
    let roll_start = b.node_count();
    b.tube(
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        inner_id,
        inner_od,
        width,
        32,
        STEEL,
    );
    // Balls on the pitch circle in the mid-plane (z = 0).
    let balls = 9;
    for j in 0..balls {
        let a = j as f64 / balls as f64 * TAU;
        let (sa, ca) = a.sin_cos();
        b.sphere([pitch_r * ca, pitch_r * sa, 0.0], ball_r, 10, 14, CHROME);
    }
    // Brass cage: a thin torus on the pitch circle in the mid-plane.
    b.torus(
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        pitch_r,
        cage_minor,
        40,
        10,
        BRASS,
    );
    let roll_end = b.node_count();

    let (mut mesh, colors) = b.into_mesh_and_colors();
    mesh.id = "valenx-bearing".to_string();

    let parts = vec![crate::RigidPart {
        node_range: roll_start..roll_end,
        axis: [0.0, 0.0, 1.0],
        pivot: [0.0, 0.0, 0.0],
        rad_per_s: SPIN_RAD_PER_S,
    }];
    Some((mesh, colors, parts))
}

/// Build the bearing [`Mesh`] (without the colour / part metadata) for the
/// central viewport. See [`bearing_solid_mesh_parts`].
fn bearing_solid_mesh(s: &BearingWorkbenchState) -> Option<Mesh> {
    bearing_solid_mesh_parts(s).map(|(mesh, _colors, _parts)| mesh)
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

/// The agent-bridge **`show_3d{kind:"bearing"}`** product: the representative
/// rolling-element bearing (outer ring, inner ring, a ring of balls, base)
/// built from the canonical deep-groove ball bearing (C = 25 kN, Fr = 5 kN
/// radial, 1500 rpm), paired with the ISO 281 rating-life + ISO 76 static
/// readout rows, at a fixed 3/4 camera. Registered in
/// [`crate::products_registry`]; the per-tool builder the registry dispatches
/// to. Pure — driven off [`BearingWorkbenchState::default`].
pub(crate) fn bearing_product() -> crate::WorkspaceProduct {
    let s = BearingWorkbenchState::default();
    let (mesh, colors, parts) =
        bearing_solid_mesh_parts(&s).expect("canonical bearing ⇒ ring-and-ball solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<bearing>/valenx-bearing");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical bearing ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Ball bearing (L10 rating life)".into(),
        lines,
        mesh: Some(loaded),
        vertex_colors: Some(colors),
        camera,
        kind2d: None,
        last_export: None,
        image: None,
        image_texture: None,
        // Animated: the inner race + balls + cage spin about the axle while the
        // outer race stays fixed. Paused at t = 0.
        animation: Some(crate::ProductAnimation {
            playing: false,
            speed: 1.0,
            t: 0.0,
            motion: crate::ProductMotion::RigidParts(parts),
        }),
    }
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

    #[test]
    fn bearing_carries_colours_and_spins_the_rolling_assembly() {
        // Per-vertex colours align to the renderer's coloured path (3/triangle),
        // with the three part colours present (steel races / chrome balls /
        // brass cage). The product spins the inner-race+balls+cage assembly about
        // the +z axle while the outer race (built first) stays fixed.
        let s = BearingWorkbenchState::default();
        let (mesh, colors, parts) = bearing_solid_mesh_parts(&s).expect("default bearing builds");
        assert_eq!(
            colors.len(),
            mesh.total_elements() * 3,
            "vertex_colors must equal 3 × triangle count"
        );
        assert!(colors.contains(&STEEL), "race colour present");
        assert!(colors.contains(&CHROME), "ball colour present");
        assert!(colors.contains(&BRASS), "cage colour present");

        assert_eq!(parts.len(), 1, "one rolling assembly");
        let p = &parts[0];
        assert_eq!(p.axis, [0.0, 0.0, 1.0], "spins about the axle");
        assert!(p.rad_per_s.abs() > 0.0, "non-zero spin rate");
        // The outer race is built first, so the rolling assembly starts past 0
        // and ends at the final node (the cage is the last thing built).
        assert!(
            p.node_range.start > 0,
            "outer race precedes the rolling assembly (fixed)"
        );
        assert_eq!(
            p.node_range.end,
            mesh.nodes.len(),
            "the rolling assembly reaches the final node"
        );
    }

    #[test]
    fn bearing_product_carries_colours_and_animation() {
        let product = bearing_product();
        let loaded = product.mesh.as_ref().expect("bearing product has a mesh");
        let colors = product
            .vertex_colors
            .as_ref()
            .expect("bearing product carries vertex_colors");
        assert_eq!(colors.len(), loaded.mesh.total_elements() * 3);
        let anim = product.animation.expect("bearing product is animated");
        assert!(!anim.playing, "starts paused");
        match anim.motion {
            crate::ProductMotion::RigidParts(parts) => assert_eq!(parts.len(), 1),
            crate::ProductMotion::Turntable { .. } => panic!("bearing uses rigid-part spin"),
        }
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use egui::accesskit::{Node, NodeId, Role};

    /// Render the whole workbench panel once in a headless egui context.
    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_bearing_workbench(app, ctx);
        });
    }

    /// As [`draw_workbench`], but with accesskit enabled, returning the emitted
    /// accessibility tree nodes — the same tree a screen reader / AI driver
    /// consumes. `accesskit` is re-exported by egui, so no extra dependency.
    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_bearing_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
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
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_a_result_and_error_without_panic() {
        let mut app = ValenxApp::default();
        app.show_bearing_workbench = true;
        run_bearing(&mut app.bearing);
        app.bearing.error = Some("invalid bearing parameters".to_string());
        draw_workbench(&mut app);
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // The bearing / load / speed / static DragValues are SpinButtons; each
        // must be `labelled_by` its caption (egui clears a DragValue's own Name),
        // so an AI / screen reader can find the control by the caption text.
        let mut app = ValenxApp::default();
        app.show_bearing_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        // C, Fr, Fa, X, Y, rpm, C0, X0, Y0.
        assert!(
            spin_buttons.len() >= 9,
            "expected the bearing numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every bearing DragValue must be labelled_by its caption (AI-drivable name)"
        );

        for caption in [
            "dynamic rating C (N)",
            "radial Fr (N)",
            "static rating C0 (N)",
        ] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
        // The Analyze button stays a named, invokable node.
        assert!(
            nodes.iter().any(|(_, n)| n.role() == Role::Button
                && n.name().is_some_and(|s| s.contains("Analyze"))),
            "the Analyze button is a named, invokable node"
        );
    }
}
