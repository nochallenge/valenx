//! The right-side **Pulley System Workbench** panel — native
//! rope-and-pulley mechanical-advantage analysis over `valenx-pulley`.
//!
//! Mirrors the Belt Drive / Heat Transfer workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_pulley_workbench`,
//! toggled from the View menu. The form picks the arrangement (a single
//! fixed pulley, a single movable pulley, or a block-and-tackle of `n`
//! supporting rope segments), the load weight, a lumped efficiency `eta`
//! and the lift distance `s`; "Analyze" reports the ideal mechanical
//! advantage `MA = n`, the equal velocity ratio, the ideal and
//! friction-aware real effort, the actual mechanical advantage `MA * eta`
//! and the load that real effort raises, plus the energy story for one
//! lift through `s` — the rope length pulled `VR * s`, the useful output
//! work `W * s`, the operator's input work `W * s / eta` and the work lost
//! to friction — and "Show 3-D" loads a representative pulley assembly (a fixed
//! block over a movable block, drawn as sheave cylinders) into the central
//! viewport.

use std::path::PathBuf;

use eframe::egui;

use valenx_pulley::{
    actual_mechanical_advantage, effort_distance, ideal_effort, input_work, load_from_effort,
    output_work, real_effort, work_lost, PulleySystem,
};

use valenx_mesh::Mesh;

use crate::mesh_prims::MeshBuilder;
use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Cast-iron sheave grey.
const CAST_IRON: [f32; 3] = [0.46, 0.47, 0.50];
/// Dark axle/hub colour.
const HUB: [f32; 3] = [0.26, 0.27, 0.30];

/// The pulley arrangement the form is configured for. The numeric
/// mechanical advantage always comes from the resulting
/// [`PulleySystem`]'s supporting-rope count, never from this tag alone.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
enum Arrangement {
    /// A single fixed pulley (`MA = 1`, direction change only).
    Fixed,
    /// A single movable pulley (`MA = 2`).
    Movable,
    /// A block-and-tackle with `supporting_ropes` segments (`MA = n`).
    BlockAndTackle,
}

/// Persistent form + result state for the Pulley System Workbench.
pub struct PulleyWorkbenchState {
    /// The selected pulley arrangement.
    arrangement: Arrangement,
    /// Number of rope segments supporting the load — used only when the
    /// arrangement is a block-and-tackle.
    supporting_ropes: u32,
    /// Load (weight) to raise (N).
    load_n: f64,
    /// Lumped mechanical efficiency `eta` in `(0, 1]`.
    efficiency: f64,
    /// Distance to raise the load (m) — drives the rope-pulled length and
    /// the work / energy readout.
    lift_distance_m: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D pulley assembly (serviced after
    /// the panel draws).
    show_3d_request: bool,
}

impl Default for PulleyWorkbenchState {
    fn default() -> Self {
        // A 4-rope block-and-tackle lifting a 1200 N load at 90%
        // efficiency: MA = VR = 4, ideal effort 300 N, real effort
        // 1200 / 3.6 = 333.33 N, actual MA 3.60. Raising it 2 m pulls
        // 8 m of rope, doing 2400 J on the load for 2666.67 J of input
        // (266.67 J lost to friction).
        Self {
            arrangement: Arrangement::BlockAndTackle,
            supporting_ropes: 4,
            load_n: 1200.0,
            efficiency: 0.9,
            lift_distance_m: 2.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Pulley System Workbench right-side panel. A no-op when the
/// `show_pulley_workbench` toggle is off.
pub fn draw_pulley_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_pulley_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_pulley_workbench",
        "Pulley System",
        |app, ui| {
            ui.label(
                egui::RichText::new("native rope-and-pulley mechanical advantage · valenx-pulley")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.pulley;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Arrangement").strong());
                    ui.radio_value(&mut s.arrangement, Arrangement::Fixed, "fixed pulley (MA = 1)");
                    ui.radio_value(
                        &mut s.arrangement,
                        Arrangement::Movable,
                        "movable pulley (MA = 2)",
                    );
                    ui.radio_value(
                        &mut s.arrangement,
                        Arrangement::BlockAndTackle,
                        "block and tackle (MA = n)",
                    );
                    if s.arrangement == Arrangement::BlockAndTackle {
                        ui.horizontal(|ui| {
                            ui.label("supporting ropes n");
                            ui.add(
                                egui::DragValue::new(&mut s.supporting_ropes)
                                    .speed(1.0)
                                    .range(1..=24),
                            );
                        });
                    }

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Load & efficiency").strong());
                    ui.horizontal(|ui| {
                        ui.label("load W (N)");
                        ui.add(egui::DragValue::new(&mut s.load_n).speed(10.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("efficiency η (0–1]");
                        ui.add(egui::DragValue::new(&mut s.efficiency).speed(0.01));
                    });
                    ui.horizontal(|ui| {
                        ui.label("lift distance s (m)");
                        ui.add(
                            egui::DragValue::new(&mut s.lift_distance_m)
                                .speed(0.1)
                                .range(0.0..=f64::INFINITY),
                        );
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_pulley(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative pulley assembly (a fixed block over a movable block, drawn as sheave cylinders) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Mechanical advantage").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_pulley_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.pulley` borrow is
    // released here): build the assembly's 3-D solid and load it.
    if app.pulley.show_3d_request {
        app.pulley.show_3d_request = false;
        load_pulley_3d(app);
    }
}

/// Validate the form, analyse the system and format the readout.
fn run_pulley(s: &mut PulleyWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build the [`PulleySystem`] the form describes, mapping any domain error
/// to a display string. Extracted so it is shared by the readout and the
/// 3-D gate.
fn system_of(s: &PulleyWorkbenchState) -> Result<PulleySystem, String> {
    match s.arrangement {
        Arrangement::Fixed => Ok(PulleySystem::fixed()),
        Arrangement::Movable => Ok(PulleySystem::movable()),
        Arrangement::BlockAndTackle => {
            PulleySystem::block_and_tackle(s.supporting_ropes).map_err(|e| e.to_string())
        }
    }
}

/// Analyse the configured system and format the full readout, mapping any
/// domain error to a display string. Extracted so it is unit-testable.
fn compute(s: &PulleyWorkbenchState) -> Result<String, String> {
    let sys = system_of(s)?;
    let ma = sys.ideal_mechanical_advantage();
    let vr = sys.velocity_ratio();
    let f_ideal = ideal_effort(sys, s.load_n).map_err(|e| e.to_string())?;
    let f_real = real_effort(sys, s.load_n, s.efficiency).map_err(|e| e.to_string())?;
    let ama = actual_mechanical_advantage(sys, s.efficiency).map_err(|e| e.to_string())?;
    // The load that real effort raises is the round-trip of `real_effort`.
    let raised = load_from_effort(sys, f_real, s.efficiency).map_err(|e| e.to_string())?;

    // Energy story for raising the load through `lift_distance_m`: the
    // operator pulls `VR * s` of rope, doing `W_in` of work for `W_out`
    // useful work on the load, with the difference lost to friction.
    let rope_pulled = effort_distance(sys, s.lift_distance_m).map_err(|e| e.to_string())?;
    let w_out = output_work(s.load_n, s.lift_distance_m).map_err(|e| e.to_string())?;
    let w_in = input_work(s.load_n, s.lift_distance_m, s.efficiency).map_err(|e| e.to_string())?;
    let w_loss = work_lost(s.load_n, s.lift_distance_m, s.efficiency).map_err(|e| e.to_string())?;

    Ok(format!(
        "arrangement     : {arrangement}\n\
         supporting ropes: {n}\n\
         load W          : {load:.1} N\n\
         efficiency η    : {eta:.3}\n\n\
         ideal MA (= n)  : {ma:.2}\n\
         velocity ratio  : {vr:.2}\n\
         ideal effort    : {f_ideal:.2} N\n\
         real effort     : {f_real:.2} N\n\
         actual MA (η·MA): {ama:.2}\n\
         raises load     : {raised:.1} N\n\n\
         lift distance s : {lift:.2} m\n\
         rope pulled     : {rope_pulled:.2} m\n\
         output work W_out: {w_out:.1} J\n\
         input work W_in : {w_in:.1} J\n\
         work lost (frict): {w_loss:.1} J",
        arrangement = arrangement_label(s.arrangement),
        n = sys.supporting_ropes(),
        load = s.load_n,
        eta = s.efficiency,
        lift = s.lift_distance_m,
    ))
}

/// A short human label for an [`Arrangement`].
fn arrangement_label(a: Arrangement) -> &'static str {
    match a {
        Arrangement::Fixed => "fixed pulley",
        Arrangement::Movable => "movable pulley",
        Arrangement::BlockAndTackle => "block and tackle",
    }
}

/// Presentation spin rate of each sheave wheel, rad/s (~0.6 rev/s) — a readable
/// inspect speed for the rope wheels.
const SHEAVE_RAD_PER_S: f32 = 4.0;

/// Append a **grooved V-sheave** (a pulley/belt wheel) to `b`, centred at
/// `centre` with its axle along `axis`, of rim radius `rim_r`, axial half-width
/// `half_w` and a central bore of radius `bore_r`. The rim is a true revolved
/// **V-groove** — the profile dips from the rim radius down to `rim_r −
/// groove_depth` at the mid-plane and back — and a darker hub tube fills the
/// bore. Built by lathing a closed `(r, z)` cross-section
/// ([`MeshBuilder::revolve`]) the way a real sheave is turned, replacing the
/// plain-cylinder fake. Shared with the Belt Drive Workbench.
///
/// Returns the [`std::ops::Range<usize>`] of node indices the **whole sheave
/// (rim + hub)** added, so the caller can record a [`crate::RigidPart`] that
/// spins it about its axle.
#[allow(clippy::too_many_arguments)]
pub(crate) fn push_grooved_sheave(
    b: &mut MeshBuilder,
    centre: [f64; 3],
    axis: [f64; 3],
    rim_r: f64,
    half_w: f64,
    bore_r: f64,
    rim_color: [f32; 3],
    hub_color: [f32; 3],
) -> std::ops::Range<usize> {
    let start = b.node_count();
    let groove_depth = (rim_r - bore_r).max(0.0) * 0.35;
    let groove_r = rim_r - groove_depth;
    // Closed (r, z) cross-section, revolved 360° about the axle. The closing
    // point repeats the first so the lathe produces a watertight solid: the
    // bore wall, the two flat faces, and the V-groove rim.
    let profile = [
        [bore_r, -half_w],   // bore edge, −z face
        [rim_r, -half_w],    // rim, −z face
        [groove_r, 0.0],     // bottom of the V at the mid-plane
        [rim_r, half_w],     // rim, +z face
        [bore_r, half_w],    // bore edge, +z face
        [bore_r, -half_w],   // close the inner wall back to the start
    ];
    b.revolve(&profile, centre, axis, 360.0, 36, rim_color);
    // Darker hub tube filling the bore (a short collar standing slightly proud).
    let hub_inner = bore_r * 0.45;
    b.tube(
        centre,
        axis,
        hub_inner,
        bore_r,
        half_w * 2.2,
        24,
        hub_color,
    );
    start..b.node_count()
}

/// Build the pulley assembly as a triangle [`Mesh`] **with per-vertex colours**
/// plus a [`crate::RigidPart`] per sheave wheel, so every grooved sheave spins
/// about its own pin while the rope/blocks frame stays put. A fixed (upper)
/// block over a movable (lower) block split the `n` supporting rope segments.
/// Each wheel is now a true revolved **V-groove sheave**
/// ([`push_grooved_sheave`]) — cast-iron grey rim, dark hub — not a plain
/// cylinder. `None` for an invalid configuration.
///
/// Returns `(mesh, colors, parts)` with `colors.len() == 3 × triangle_count`.
fn pulley_solid_mesh_parts(
    s: &PulleyWorkbenchState,
) -> Option<(Mesh, Vec<[f32; 3]>, Vec<crate::RigidPart>)> {
    let sys = system_of(s).ok()?;
    let n = sys.supporting_ropes();

    let mut b = MeshBuilder::new();
    let mut parts: Vec<crate::RigidPart> = Vec::new();

    // Split the supporting segments between a fixed (upper) and movable
    // (lower) block: ceil(n/2) sheaves up, floor(n/2) (at least one) down.
    let fixed_sheaves = n.div_ceil(2).max(1);
    let movable_sheaves = (n / 2).max(1);
    let rim_r = 0.18;
    let half_w = 0.06;
    let bore_r = 0.05;
    let pitch = 0.5;

    // Each sheave is its own rotating body about its pin (+x axle through the
    // centre): record its node range + pin centre.
    let add_sheave = |b: &mut MeshBuilder, parts: &mut Vec<crate::RigidPart>, centre: [f64; 3]| {
        let range = push_grooved_sheave(
            b,
            centre,
            [1.0, 0.0, 0.0],
            rim_r,
            half_w,
            bore_r,
            CAST_IRON,
            HUB,
        );
        parts.push(crate::RigidPart {
            node_range: range,
            axis: [1.0, 0.0, 0.0],
            pivot: [centre[0] as f32, centre[1] as f32, centre[2] as f32],
            rad_per_s: SHEAVE_RAD_PER_S,
        });
    };

    // Fixed block: a row of sheaves along y at the top (+z).
    let fixed_span = (fixed_sheaves.saturating_sub(1)) as f64 * pitch;
    for i in 0..fixed_sheaves {
        let y = i as f64 * pitch - 0.5 * fixed_span;
        add_sheave(&mut b, &mut parts, [0.0, y, 0.9]);
    }

    // Movable block: a row of sheaves lower down (the load hangs here).
    let movable_span = (movable_sheaves.saturating_sub(1)) as f64 * pitch;
    for i in 0..movable_sheaves {
        let y = i as f64 * pitch - 0.5 * movable_span;
        add_sheave(&mut b, &mut parts, [0.0, y, 0.1]);
    }

    let (mut mesh, colors) = b.into_mesh_and_colors();
    mesh.id = "valenx-pulley".to_string();
    Some((mesh, colors, parts))
}

/// Build the pulley assembly as a triangle [`Mesh`] (without the colour / part
/// metadata) for the central viewport. See [`pulley_solid_mesh_parts`].
fn pulley_solid_mesh(s: &PulleyWorkbenchState) -> Option<Mesh> {
    pulley_solid_mesh_parts(s).map(|(mesh, _colors, _parts)| mesh)
}

/// Build the 3-D pulley assembly and load it into the central viewport.
fn load_pulley_3d(app: &mut ValenxApp) {
    let Some(mesh) = pulley_solid_mesh(&app.pulley) else {
        app.pulley.error =
            Some("pulley parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<pulley>/valenx-pulley"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"pulley"}`** product: the representative
/// pulley assembly (a fixed block of sheave cylinders over a movable block)
/// built from the canonical 4-rope block-and-tackle (1200 N load, 90 %
/// efficiency, 2 m lift), paired with the mechanical-advantage + energy
/// readout rows, at a fixed 3/4 camera. Registered in
/// [`crate::products_registry`]; the per-tool builder the registry dispatches
/// to. Pure — driven off [`PulleyWorkbenchState::default`].
pub(crate) fn pulley_product() -> crate::WorkspaceProduct {
    let s = PulleyWorkbenchState::default();
    let (mesh, colors, parts) =
        pulley_solid_mesh_parts(&s).expect("canonical tackle ⇒ sheave-block solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<pulley>/valenx-pulley");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical tackle ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Block & tackle (MA = 4)".into(),
        lines,
        mesh: Some(loaded),
        vertex_colors: Some(colors),
        camera,
        kind2d: None,
        last_export: None,
        image: None,
        image_texture: None,
        // Animated: each sheave wheel spins about its own pin. Paused at t = 0.
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
        let s = PulleyWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_advantage_and_effort() {
        let mut s = PulleyWorkbenchState::default();
        run_pulley(&mut s);
        assert!(
            s.error.is_none(),
            "default tackle should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("ideal MA"));
        assert!(s.result.contains("actual MA"));
        // n = 4 -> MA 4.00, ideal effort 1200 / 4 = 300.00 N, real effort
        // 1200 / (4 * 0.9) = 333.33 N, actual MA 4 * 0.9 = 3.60.
        assert!(s.result.contains("4.00"));
        assert!(s.result.contains("300.00"));
        assert!(s.result.contains("333.33"));
        assert!(s.result.contains("3.60"));
    }

    #[test]
    fn fixed_pulley_has_unit_advantage() {
        let mut s = PulleyWorkbenchState {
            arrangement: Arrangement::Fixed,
            load_n: 500.0,
            efficiency: 1.0,
            ..Default::default()
        };
        run_pulley(&mut s);
        assert!(s.error.is_none());
        // MA = 1, ideal effort == load == 500.0 N.
        assert!(s.result.contains("ideal MA (= n)  : 1.00"));
        assert!(s.result.contains("ideal effort    : 500.00 N"));
    }

    #[test]
    fn analyze_rejects_efficiency_above_one() {
        let mut s = PulleyWorkbenchState {
            efficiency: 1.5,
            ..Default::default()
        };
        run_pulley(&mut s);
        assert!(s.error.is_some());
    }

    /// Ground truth: the ideal mechanical advantage of a block-and-tackle
    /// equals its number of supporting rope segments, and the ideal effort
    /// is the load divided by that count — checked against hand values.
    #[test]
    fn ideal_advantage_equals_supporting_ropes() {
        let load = 600.0;
        for n in 1..=8u32 {
            let sys = PulleySystem::block_and_tackle(n).unwrap();
            assert!((sys.ideal_mechanical_advantage() - f64::from(n)).abs() < 1e-12);
            let f = ideal_effort(sys, load).unwrap();
            let expected = load / f64::from(n);
            assert!((f - expected).abs() < 1e-9, "n = {n}");
        }
    }

    /// Ground truth for the energy block surfaced from `valenx-pulley`'s
    /// `effort_distance` / `output_work` / `input_work` / `work_lost`.
    /// Defaults: n = 4 (VR = 4), W = 1200 N, η = 0.9, lifting s = 2 m.
    /// Hand values: rope pulled = VR·s = 4·2 = 8.00 m; output work
    /// W_out = W·s = 1200·2 = 2400 J; input work W_in = W_out/η =
    /// 2400/0.9 = 2666.67 J; work lost = W_in − W_out = 266.67 J.
    #[test]
    fn analyze_default_reports_rope_pulled_and_work() {
        let mut s = PulleyWorkbenchState::default();
        run_pulley(&mut s);
        assert!(
            s.error.is_none(),
            "default tackle should analyze: {:?}",
            s.error
        );
        assert!(
            s.result.contains("rope pulled     : 8.00 m"),
            "{}",
            s.result
        );
        assert!(
            s.result.contains("output work W_out: 2400.0 J"),
            "{}",
            s.result
        );
        assert!(
            s.result.contains("input work W_in : 2666.7 J"),
            "{}",
            s.result
        );
        assert!(
            s.result.contains("work lost (frict): 266.7 J"),
            "{}",
            s.result
        );
    }

    /// Cross-check the energy identity directly against `valenx-pulley`:
    /// the surfaced rope length is VR·s and W_in − W_out equals the work
    /// lost, with no loss at η = 1.
    #[test]
    fn energy_block_matches_crate_ground_truth() {
        let load = 1200.0_f64;
        let s_load = 2.0_f64;
        let eta = 0.9_f64;
        let sys = PulleySystem::block_and_tackle(4).unwrap();

        let rope = effort_distance(sys, s_load).unwrap();
        assert!((rope - 8.0).abs() < 1e-12, "rope = {rope}");

        let w_out = output_work(load, s_load).unwrap();
        let w_in = input_work(load, s_load, eta).unwrap();
        let w_loss = work_lost(load, s_load, eta).unwrap();
        assert!((w_out - 2400.0).abs() < 1e-9, "w_out = {w_out}");
        assert!((w_in - 2400.0 / 0.9).abs() < 1e-9, "w_in = {w_in}");
        assert!((w_loss - (w_in - w_out)).abs() < 1e-9, "w_loss = {w_loss}");

        // No friction loss for the perfect machine.
        let lossless = work_lost(load, s_load, 1.0).unwrap();
        assert!(lossless.abs() < 1e-9, "lossless = {lossless}");
    }

    #[test]
    fn pulley_mesh_for_default_is_nonempty_and_in_range() {
        let s = PulleyWorkbenchState::default();
        let mesh = pulley_solid_mesh(&s).expect("default tackle yields a solid");
        assert!(mesh.nodes.len() > 8, "expected several sheave cylinders");
        let count = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < count));
        }
    }

    #[test]
    fn pulley_mesh_none_for_invalid() {
        let s = PulleyWorkbenchState {
            arrangement: Arrangement::BlockAndTackle,
            supporting_ropes: 0,
            ..Default::default()
        };
        assert!(pulley_solid_mesh(&s).is_none());
    }

    #[test]
    fn pulley_carries_vertex_aligned_colours() {
        // The grooved sheaves ship per-vertex colours aligned to the renderer's
        // coloured path (3/triangle), with both the cast-iron rim and the dark
        // hub colours present.
        let s = PulleyWorkbenchState::default();
        let (mesh, colors, _parts) =
            pulley_solid_mesh_parts(&s).expect("default tackle builds");
        assert_eq!(
            colors.len(),
            mesh.total_elements() * 3,
            "vertex_colors must equal 3 × triangle count"
        );
        assert!(colors.contains(&CAST_IRON), "sheave rim colour present");
        assert!(colors.contains(&HUB), "hub colour present");
        for c in &colors {
            for ch in c {
                assert!(ch.is_finite() && (0.0..=1.0).contains(ch));
            }
        }
    }

    #[test]
    fn pulley_product_spins_each_sheave() {
        // The product carries a RigidParts animation: one rotating part per
        // sheave, each a non-empty node range within the mesh that spins about
        // its pin (+x) at a non-zero rate.
        let product = pulley_product();
        let loaded = product.mesh.as_ref().expect("pulley product has a mesh");
        let node_count = loaded.mesh.nodes.len();
        let anim = product.animation.expect("pulley product is animated");
        assert!(!anim.playing, "starts paused");
        match anim.motion {
            crate::ProductMotion::RigidParts(parts) => {
                assert!(
                    parts.len() >= 2,
                    "default 4-rope tackle splits into ≥2 sheaves"
                );
                for p in &parts {
                    assert!(
                        p.node_range.start < p.node_range.end,
                        "non-empty sheave range"
                    );
                    assert!(
                        p.node_range.end <= node_count,
                        "sheave range within the mesh"
                    );
                    assert_eq!(p.axis, [1.0, 0.0, 0.0], "each sheave spins about its pin");
                    assert!(p.rad_per_s.abs() > 0.0, "non-zero spin rate");
                }
                // Sheave ranges tile contiguously (no gap/overlap).
                for w in parts.windows(2) {
                    assert_eq!(
                        w[0].node_range.end, w[1].node_range.start,
                        "sheave node ranges are contiguous"
                    );
                }
            }
            crate::ProductMotion::Turntable { .. } => {
                panic!("pulley must use per-part rigid motion")
            }
        }
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_pulley_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_pulley_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_pulley_workbench = true;
        run_pulley(&mut app.pulley);
        draw_workbench(&mut app);
    }
}
