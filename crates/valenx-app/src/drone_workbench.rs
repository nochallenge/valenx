//! The right-side **Drone / Multirotor Workbench** panel — native
//! multirotor hover performance over `valenx-drone`.
//!
//! Mirrors the Springs / Marine workbenches: a resizable [`egui::SidePanel`]
//! gated on `crate::ValenxApp::show_drone_workbench`, toggled from the View
//! menu. The form drives a [`valenx_drone::Multirotor`]; "Analyze" reports
//! the disk loading, the hover induced velocity and the ideal / actual
//! hover power plus a battery-limited hover endurance, and "Show 3-D drone"
//! loads a hub + rotor-disk solid into the central viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;

use valenx_drone::{Multirotor, SEA_LEVEL_AIR_DENSITY};
use valenx_mesh::Mesh;

use crate::mesh_prims::MeshBuilder;
use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Presentation spin rate for the rotors (rad/s of animation time). Alternating
/// rotors counter-rotate (the sign flips per rotor) like a real quad.
const ROTOR_RAD_PER_S: f32 = 6.0;

/// Carbon-dark frame / arms.
const FRAME: [f32; 3] = [0.16, 0.17, 0.19];
/// Electronics stack (flight controller) — a cooler grey-blue.
const STACK: [f32; 3] = [0.30, 0.34, 0.42];
/// Anodised motor body.
const MOTOR: [f32; 3] = [0.55, 0.40, 0.20];
/// Propeller blades — light grey so they read against the dark frame.
const PROP: [f32; 3] = [0.78, 0.80, 0.84];
/// Rotor hub (prop nut / bell).
const HUB: [f32; 3] = [0.10, 0.10, 0.12];

/// Persistent form + result state for the Drone / Multirotor Workbench.
pub struct DroneWorkbenchState {
    /// Number of lifting rotors `n`.
    rotor_count: u32,
    /// Rotor radius `R` (m).
    rotor_radius_m: f64,
    /// All-up mass `m` (kg).
    mass_kg: f64,
    /// Rotor figure of merit `FM` in `(0, 1]`.
    figure_of_merit: f64,
    /// Air density `rho` (kg/m^3).
    air_density: f64,
    /// Usable battery energy (Wh) for the hover-endurance estimate.
    battery_wh: f64,
    /// Formatted hover-performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D drone solid and load it into the
    /// central viewport (serviced after the panel draws).
    show_3d_request: bool,
}

impl Default for DroneWorkbenchState {
    fn default() -> Self {
        // A 1.5 kg quadcopter with 15 cm rotors at sea level.
        Self {
            rotor_count: 4,
            rotor_radius_m: 0.15,
            mass_kg: 1.5,
            figure_of_merit: 0.70,
            air_density: SEA_LEVEL_AIR_DENSITY,
            battery_wh: 50.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Drone / Multirotor Workbench right-side panel. A no-op when the
/// `show_drone_workbench` toggle is off.
pub fn draw_drone_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_drone_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_drone_workbench",
        "Drone / Multirotor",
        |app, ui| {
            ui.label(
                egui::RichText::new("native multirotor hover performance · valenx-drone")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.drone;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Rotors + airframe").strong());
                    ui.horizontal(|ui| {
                        ui.label("rotor count n");
                        ui.add(egui::DragValue::new(&mut s.rotor_count).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("rotor radius R (m)");
                        ui.add(egui::DragValue::new(&mut s.rotor_radius_m).speed(0.005));
                    });
                    ui.horizontal(|ui| {
                        ui.label("all-up mass (kg)");
                        ui.add(egui::DragValue::new(&mut s.mass_kg).speed(0.05));
                    });
                    ui.horizontal(|ui| {
                        ui.label("figure of merit");
                        ui.add(egui::DragValue::new(&mut s.figure_of_merit).speed(0.01));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Air + battery").strong());
                    ui.horizontal(|ui| {
                        ui.label("air ρ (kg/m³)");
                        ui.add(egui::DragValue::new(&mut s.air_density).speed(0.01));
                    });
                    ui.horizontal(|ui| {
                        ui.label("battery (Wh)");
                        ui.add(egui::DragValue::new(&mut s.battery_wh).speed(1.0));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_drone(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D drone").strong())
                        .on_hover_text(
                            "Build the hub + rotor disks as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Hover performance").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_drone_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.drone` borrow is
    // released here): build the drone's 3-D solid and load it.
    if app.drone.show_3d_request {
        app.drone.show_3d_request = false;
        load_drone_3d(app);
    }
}

/// Build a validated [`Multirotor`] from the form, mapping the domain error
/// to a display string.
fn build_drone(s: &DroneWorkbenchState) -> Result<Multirotor, String> {
    Multirotor::new(
        s.rotor_count,
        s.rotor_radius_m,
        s.mass_kg,
        s.figure_of_merit,
        s.air_density,
    )
    .map_err(|e| e.to_string())
}

/// Validate the form, compute the hover performance and format the readout.
fn run_drone(s: &mut DroneWorkbenchState) {
    s.error = None;
    match build_drone(s) {
        Ok(m) => {
            let p = m.hover_performance();
            let endurance = m.hover_endurance_minutes(s.battery_wh.max(0.0));
            s.result = format!(
                "rotors n       : {}\n\
                 rotor radius R : {:.3} m\n\
                 all-up mass m  : {:.2} kg\n\
                 figure of merit: {:.2}\n\n\
                 weight         : {:.1} N\n\
                 disk area      : {:.4} m²\n\
                 disk loading   : {:.1} N/m²\n\
                 induced vel vi : {:.2} m/s\n\
                 ideal power    : {:.0} W\n\
                 actual power   : {:.0} W  (P_ideal / FM)\n\
                 hover endurance: {:.1} min  ({:.0} Wh battery)",
                s.rotor_count,
                s.rotor_radius_m,
                s.mass_kg,
                s.figure_of_merit,
                p.weight_n,
                p.disk_area_m2,
                p.disk_loading_pa,
                p.induced_velocity_m_s,
                p.ideal_power_w,
                p.actual_power_w,
                endurance,
                s.battery_wh,
            );
        }
        Err(e) => s.error = Some(e),
    }
}

/// Number of blades on each rotor (a real prop reads as 2 blades, not a disc).
const BLADES_PER_ROTOR: usize = 2;
/// Points around one closed blade airfoil section.
const BLADE_SECTION_POINTS: usize = 10;

/// A closed thin **cambered blade airfoil** loop of chord `chord` and thickness
/// `thick`, in local `[chordwise, thickness]` coordinates centred on the chord.
/// A flattened ellipse-ish lens (upper arc + lower arc) — cheap but reads as a
/// blade cross-section once swept. Counter-clockwise.
fn blade_section(chord: f64, thick: f64) -> Vec<[f64; 2]> {
    let n = BLADE_SECTION_POINTS;
    let half_c = chord * 0.5;
    let mut pts = Vec::with_capacity(2 * n - 2);
    // Upper surface, leading → trailing edge.
    for i in 0..n {
        let t = i as f64 / (n - 1) as f64;
        let x = -half_c + t * chord;
        // Parabolic thickness, zero at both edges, max mid-chord.
        let y = thick * 0.5 * (1.0 - (2.0 * t - 1.0).powi(2));
        pts.push([x, y]);
    }
    // Lower surface, trailing → leading, skipping shared endpoints.
    for i in 1..n - 1 {
        let t = (n - 1 - i) as f64 / (n - 1) as f64;
        let x = -half_c + t * chord;
        let y = -thick * 0.5 * (1.0 - (2.0 * t - 1.0).powi(2));
        pts.push([x, y]);
    }
    pts
}

/// Place a blade airfoil section into 3-D at radial station `radius` from the
/// rotor centre `c`, along blade azimuth `blade_a` (radians), set to pitch
/// angle `pitch` (radians, twisting the chord about the radial axis). The
/// chord lies tangent to the rotor circle; thickness tilts by the pitch so the
/// lofted blade is a real twisted aerofoil. Returns one ring for
/// [`MeshBuilder::loft`].
fn blade_ring(
    profile: &[[f64; 2]],
    c: [f64; 3],
    radius: f64,
    blade_a: f64,
    pitch: f64,
    z: f64,
) -> Vec<[f64; 3]> {
    let (sa, ca) = blade_a.sin_cos();
    // Radial unit (out along the blade) and tangent unit (chordwise) in plane.
    let radial = [ca, sa];
    let tangent = [-sa, ca];
    let (sp, cp) = pitch.sin_cos();
    profile
        .iter()
        .map(|&[chordwise, th]| {
            // Pitch rotates (chordwise, thickness) in the tangent/vertical plane.
            let tang = chordwise * cp - th * sp;
            let vert = chordwise * sp + th * cp;
            [
                c[0] + radial[0] * radius + tangent[0] * tang,
                c[1] + radial[1] * radius + tangent[1] * tang,
                c[2] + z + vert,
            ]
        })
        .collect()
}

/// Append one rotor — a small **hub cylinder** plus [`BLADES_PER_ROTOR`] real
/// **twisted, tapered propeller blades** lofted from a root section to a tip
/// section — centred at `c`, and return the combined node [`Range<usize>`]
/// covering the whole rotor so the caller can spin it rigidly. Hub coloured
/// [`HUB`], blades [`PROP`].
fn push_rotor(b: &mut MeshBuilder, c: [f64; 3], r: f64) -> std::ops::Range<usize> {
    let start = b.node_count();
    // Hub: a short cylinder on the rotor axis, sitting just above the motor.
    let hub_r = r * 0.12;
    let hub_h = r * 0.18;
    b.cylinder([c[0], c[1], c[2] + hub_h * 0.5], [0.0, 0.0, 1.0], hub_r, hub_h, 16, HUB);
    // Blades radiate from the hub. Root near the hub, tip near radius r.
    let root_radius = hub_r * 1.2;
    let tip_radius = r * 0.98;
    let root_chord = r * 0.30;
    let tip_chord = r * 0.12;
    let blade_thick = r * 0.05;
    let z_blade = c[2] + hub_h; // blades sit on top of the hub
    for k in 0..BLADES_PER_ROTOR {
        let blade_a = k as f64 / BLADES_PER_ROTOR as f64 * TAU;
        // Root pitched ~18°, tip washed out to ~6° (typical blade twist).
        let root = blade_ring(
            &blade_section(root_chord, blade_thick),
            c,
            root_radius,
            blade_a,
            18.0_f64.to_radians(),
            z_blade,
        );
        let tip = blade_ring(
            &blade_section(tip_chord, blade_thick * 0.6),
            c,
            tip_radius,
            blade_a,
            6.0_f64.to_radians(),
            z_blade,
        );
        b.loft(&[root, tip], true, PROP);
    }
    start..b.node_count()
}

/// Build the multirotor as a triangle [`Mesh`] **with per-vertex colours** plus
/// one [`crate::RigidPart`] per rotor, so each rotor (hub + blades) spins about
/// its own shaft while the frame stays put. The airframe is a real **frame
/// plate** + electronics stack, the `n` arms are round **cylinder** booms (with
/// a motor can at each end), and each rotor is a hub + [`BLADES_PER_ROTOR`]
/// twisted tapered **propeller blades** (not a flat disc). Alternating rotors
/// counter-rotate. `None` for an invalid configuration.
///
/// Returns `(mesh, colors, parts)` with `colors.len() == 3 × triangle_count`.
fn drone_solid_mesh_parts(
    s: &DroneWorkbenchState,
) -> Option<(Mesh, Vec<[f32; 3]>, Vec<crate::RigidPart>)> {
    let m = build_drone(s).ok()?;
    let n = m.rotor_count as usize;
    let r = m.rotor_radius_m;
    let arm = r * 2.3;

    let mut b = MeshBuilder::new();

    // Central frame plate (thin square) + flight-controller stack on top.
    let body_half = r * 0.55;
    b.cuboid([0.0, 0.0, 0.0], [body_half * 2.0, body_half * 2.0, r * 0.10], FRAME);
    b.cuboid(
        [0.0, 0.0, r * 0.16],
        [body_half * 1.1, body_half * 1.1, r * 0.16],
        STACK,
    );

    // Arms + motors + spinning rotors, evenly spaced around the frame.
    let mut parts = Vec::with_capacity(n);
    for i in 0..n {
        let a = i as f64 / n as f64 * TAU;
        let (sa, ca) = a.sin_cos();
        let motor_c = [arm * ca, arm * sa, 0.0];
        // Round arm boom from the frame edge out to the motor (a real tube).
        let inner = [body_half * 0.8 * ca, body_half * 0.8 * sa, 0.0];
        let mid = [
            (inner[0] + motor_c[0]) * 0.5,
            (inner[1] + motor_c[1]) * 0.5,
            0.0,
        ];
        let arm_len = ((motor_c[0] - inner[0]).powi(2) + (motor_c[1] - inner[1]).powi(2)).sqrt();
        b.cylinder(mid, [ca, sa, 0.0], r * 0.06, arm_len, 12, FRAME);
        // Motor can under each rotor.
        let motor_h = r * 0.22;
        b.cylinder(
            [motor_c[0], motor_c[1], motor_h * 0.5],
            [0.0, 0.0, 1.0],
            r * 0.16,
            motor_h,
            16,
            MOTOR,
        );
        // The rotor (hub + blades) — record its range and counter-rotate odds.
        let rotor_range = push_rotor(&mut b, [motor_c[0], motor_c[1], motor_h], r);
        let dir = if i % 2 == 0 { 1.0 } else { -1.0 };
        parts.push(crate::RigidPart {
            node_range: rotor_range,
            axis: [0.0, 0.0, 1.0],
            pivot: [motor_c[0] as f32, motor_c[1] as f32, (motor_h) as f32],
            rad_per_s: ROTOR_RAD_PER_S * dir,
        });
    }

    let (mut mesh, colors) = b.into_mesh_and_colors();
    mesh.id = "valenx-drone".to_string();
    Some((mesh, colors, parts))
}

/// Build the multirotor [`Mesh`] (without colour / part metadata) for the
/// central viewport. See [`drone_solid_mesh_parts`].
fn drone_solid_mesh(s: &DroneWorkbenchState) -> Option<Mesh> {
    drone_solid_mesh_parts(s).map(|(mesh, _colors, _parts)| mesh)
}

/// Build the 3-D drone solid and load it into the central viewport.
fn load_drone_3d(app: &mut ValenxApp) {
    let Some(mesh) = drone_solid_mesh(&app.drone) else {
        app.drone.error = Some("drone parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<drone>/valenx-drone"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"drone"}`** product: the canonical hub +
/// rotor-disk multirotor solid (the panel's "Show 3-D drone" geometry) paired
/// with the workbench's own hover-performance headline numbers, at a fixed
/// 3/4 camera. Registered in [`crate::products_registry`]; the per-tool
/// builder the registry dispatches to. Pure — driven off
/// [`DroneWorkbenchState::default`].
///
/// The readout rows mirror the panel's `run_drone` "Hover performance"
/// readout (this workbench formats its readout into `result` rather than via a
/// shared `compute()` string fn, so the rows are taken from that here).
pub(crate) fn drone_product() -> crate::WorkspaceProduct {
    let mut s = DroneWorkbenchState::default();
    run_drone(&mut s);
    let (mesh, colors, parts) =
        drone_solid_mesh_parts(&s).expect("default quadcopter ⇒ a 3-D solid");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<drone>/valenx-drone");
    let lines = crate::products_registry::lines_from_readout(&s.result);
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Drone / Multirotor".into(),
        lines,
        mesh: Some(loaded),
        vertex_colors: Some(colors),
        camera,
        kind2d: None,
        last_export: None,
        image: None,
        image_texture: None,
        // Animated: each rotor (hub + blades) spins about its own shaft;
        // alternating rotors counter-rotate like a real quad. Paused at t = 0.
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
        let s = DroneWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_drone_reports_power_and_endurance() {
        let mut s = DroneWorkbenchState::default();
        run_drone(&mut s);
        assert!(s.error.is_none());
        assert!(s.result.contains("actual power"));
        assert!(s.result.contains("hover endurance"));
        assert!(s.result.contains("disk loading"));
    }

    #[test]
    fn analyze_rejects_zero_rotors_and_bad_fm() {
        let mut z = DroneWorkbenchState {
            rotor_count: 0,
            ..Default::default()
        };
        run_drone(&mut z);
        assert!(z.error.is_some());
        let mut fm = DroneWorkbenchState {
            figure_of_merit: 1.5,
            ..Default::default()
        };
        run_drone(&mut fm);
        assert!(fm.error.is_some());
    }

    #[test]
    fn drone_mesh_for_default_is_nonempty_and_in_range() {
        let s = DroneWorkbenchState::default();
        let mesh = drone_solid_mesh(&s).expect("default drone yields a solid");
        assert!(mesh.nodes.len() > 8, "expected hub + rotor disks");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn drone_mesh_none_for_invalid() {
        let s = DroneWorkbenchState {
            rotor_radius_m: 0.0,
            ..Default::default()
        };
        assert!(drone_solid_mesh(&s).is_none());
    }

    #[test]
    fn blade_section_is_a_closed_thin_loop() {
        // A blade lens: closed loop of 2·n−2 points spanning the chord, with a
        // real (non-zero) thickness and symmetric upper/lower surfaces.
        let chord = 0.1;
        let thick = 0.01;
        let pts = blade_section(chord, thick);
        assert_eq!(pts.len(), 2 * BLADE_SECTION_POINTS - 2);
        let max_x = pts.iter().map(|p| p[0]).fold(f64::MIN, f64::max);
        let min_x = pts.iter().map(|p| p[0]).fold(f64::MAX, f64::min);
        assert!((max_x - chord * 0.5).abs() < 1e-9 && (min_x + chord * 0.5).abs() < 1e-9);
        let max_y = pts.iter().map(|p| p[1]).fold(f64::MIN, f64::max);
        // Parabolic thickness peaks at mid-chord (≤ thick/2 at the samples).
        assert!(max_y > 0.0 && max_y <= thick * 0.5 + 1e-9, "real thickness");
        assert!(max_y > thick * 0.4, "near the max-thickness fraction");
        let min_y = pts.iter().map(|p| p[1]).fold(f64::MAX, f64::min);
        assert!((max_y + min_y).abs() < 1e-9, "symmetric about the chord");
    }

    #[test]
    fn drone_carries_vertex_aligned_colours() {
        // Frame + arms + motors + twisted-blade rotors ship per-vertex colours
        // aligned to the renderer's coloured path (3 / triangle), with the
        // frame, prop and motor colours present.
        let s = DroneWorkbenchState::default();
        let (mesh, colors, _parts) =
            drone_solid_mesh_parts(&s).expect("default drone builds coloured");
        assert!(!mesh.nodes.is_empty(), "non-empty mesh");
        assert!(mesh.total_elements() > 0, "mesh has triangles");
        assert_eq!(
            colors.len(),
            mesh.total_elements() * 3,
            "vertex_colors must equal 3 × triangle count"
        );
        assert!(colors.contains(&FRAME), "frame colour present");
        assert!(colors.contains(&PROP), "propeller colour present");
        assert!(colors.contains(&MOTOR), "motor colour present");
        for c in &colors {
            for ch in c {
                assert!(ch.is_finite() && (0.0..=1.0).contains(ch));
            }
        }
    }

    #[test]
    fn drone_product_spins_each_rotor_counter_rotating() {
        // The product is animated: one RigidPart per rotor, each a non-empty
        // node range within the mesh spinning about +z; alternating rotors
        // counter-rotate (signs alternate). The frame is left static.
        let product = drone_product();
        let loaded = product.mesh.as_ref().expect("drone product has a mesh");
        let node_count = loaded.mesh.nodes.len();
        let anim = product.animation.expect("drone product is animated");
        assert!(!anim.playing, "starts paused");
        match anim.motion {
            crate::ProductMotion::RigidParts(parts) => {
                assert_eq!(parts.len(), 4, "four rotors on the default quad");
                for p in &parts {
                    assert!(p.node_range.start < p.node_range.end, "non-empty rotor range");
                    assert!(p.node_range.end <= node_count, "rotor range within the mesh");
                    assert_eq!(p.axis, [0.0, 0.0, 1.0], "spins about the rotor axis");
                    assert!(p.rad_per_s.abs() > 0.0, "rotor spins");
                }
                // Adjacent rotors counter-rotate.
                assert!(
                    parts[0].rad_per_s.signum() != parts[1].rad_per_s.signum(),
                    "adjacent rotors counter-rotate"
                );
                // Rotor ranges don't overlap (they tile in order).
                for w in 0..parts.len() - 1 {
                    assert!(
                        parts[w].node_range.end <= parts[w + 1].node_range.start,
                        "rotor ranges are ordered + disjoint"
                    );
                }
            }
            crate::ProductMotion::Turntable { .. } => panic!("drone must use per-part motion"),
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
            draw_drone_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_drone_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_drone_workbench = true;
        run_drone(&mut app.drone);
        draw_workbench(&mut app);
    }
}
