//! The right-side **Fixed-Wing / Aircraft Workbench** panel — native
//! preliminary aircraft point-performance over `valenx-fixedwing`.
//!
//! Mirrors the Drone / Rail / Marine workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_fixedwing_workbench`,
//! toggled from the View menu. The form drives a [`valenx_fixedwing::Aircraft`];
//! "Analyze" reports the wing loading, the stall speed, the cruise lift /
//! drag, the best lift-to-drag ratio and the still-air glide range, and
//! "Show 3-D aircraft" loads a fuselage + wings + tail solid into the
//! central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_fixedwing::{Aircraft, SEA_LEVEL_AIR_DENSITY};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Fixed-Wing / Aircraft Workbench.
pub struct FixedWingWorkbenchState {
    /// Reference wing area `S` (m^2).
    wing_area_m2: f64,
    /// All-up mass `m` (kg).
    mass_kg: f64,
    /// Maximum (clean) lift coefficient `CLmax`.
    cl_max: f64,
    /// Wing aspect ratio `AR`.
    aspect_ratio: f64,
    /// Zero-lift drag coefficient `CD0`.
    cd0: f64,
    /// Oswald span efficiency `e` in `(0, 1]`.
    oswald_efficiency: f64,
    /// Air density `rho` (kg/m^3).
    air_density: f64,
    /// Cruise airspeed (m/s) the level-flight readout is evaluated at.
    cruise_speed_m_s: f64,
    /// Altitude (m) the still-air glide range is evaluated from.
    glide_altitude_m: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D aircraft solid and load it into the
    /// central viewport (serviced after the panel draws).
    show_3d_request: bool,
}

impl Default for FixedWingWorkbenchState {
    fn default() -> Self {
        // A light single-engine GA aircraft at sea level.
        Self {
            wing_area_m2: 16.0,
            mass_kg: 1100.0,
            cl_max: 1.6,
            aspect_ratio: 7.5,
            cd0: 0.027,
            oswald_efficiency: 0.8,
            air_density: SEA_LEVEL_AIR_DENSITY,
            cruise_speed_m_s: 50.0,
            glide_altitude_m: 3000.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Fixed-Wing / Aircraft Workbench right-side panel. A no-op when
/// the `show_fixedwing_workbench` toggle is off.
pub fn draw_fixedwing_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_fixedwing_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_fixedwing_workbench",
        "Fixed-Wing / Aircraft",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native preliminary aircraft point-performance · valenx-fixedwing",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.fixedwing;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Wing + mass").strong());
                    ui.horizontal(|ui| {
                        ui.label("wing area S (m²)");
                        ui.add(egui::DragValue::new(&mut s.wing_area_m2).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("all-up mass (kg)");
                        ui.add(egui::DragValue::new(&mut s.mass_kg).speed(10.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("aspect ratio AR");
                        ui.add(egui::DragValue::new(&mut s.aspect_ratio).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("CLmax");
                        ui.add(egui::DragValue::new(&mut s.cl_max).speed(0.05));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Drag polar + air").strong());
                    ui.horizontal(|ui| {
                        ui.label("CD0");
                        ui.add(egui::DragValue::new(&mut s.cd0).speed(0.001));
                    });
                    ui.horizontal(|ui| {
                        ui.label("Oswald e");
                        ui.add(egui::DragValue::new(&mut s.oswald_efficiency).speed(0.01));
                    });
                    ui.horizontal(|ui| {
                        ui.label("air ρ (kg/m³)");
                        ui.add(egui::DragValue::new(&mut s.air_density).speed(0.01));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Operating point").strong());
                    ui.horizontal(|ui| {
                        ui.label("cruise v (m/s)");
                        ui.add(egui::DragValue::new(&mut s.cruise_speed_m_s).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("glide altitude (m)");
                        ui.add(egui::DragValue::new(&mut s.glide_altitude_m).speed(50.0));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_aircraft(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D aircraft").strong())
                        .on_hover_text(
                            "Build the fuselage + wings + tail as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Point performance").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_fixedwing_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.fixedwing` borrow is
    // released here): build the aircraft's 3-D solid and load it.
    if app.fixedwing.show_3d_request {
        app.fixedwing.show_3d_request = false;
        load_aircraft_3d(app);
    }
}

/// Build a validated [`Aircraft`] from the form, mapping the domain error to
/// a display string.
fn build_aircraft(s: &FixedWingWorkbenchState) -> Result<Aircraft, String> {
    Aircraft::new(
        s.wing_area_m2,
        s.mass_kg,
        s.cl_max,
        s.aspect_ratio,
        s.cd0,
        s.oswald_efficiency,
        s.air_density,
    )
    .map_err(|e| e.to_string())
}

/// Validate the form, compute the point-performance and format the readout.
fn run_aircraft(s: &mut FixedWingWorkbenchState) {
    s.error = None;
    match build_aircraft(s) {
        Ok(a) => {
            let p = a.performance();
            let vc = s.cruise_speed_m_s.max(0.0);
            let cl_cruise = a.level_flight_cl(vc);
            let cd_cruise = a.drag_coefficient(cl_cruise);
            let ld_cruise = if cd_cruise > 0.0 {
                cl_cruise / cd_cruise
            } else {
                0.0
            };
            let range_km = a.glide_range(s.glide_altitude_m) / 1000.0;
            s.result = format!(
                "wing area S    : {:.1} m²\n\
                 all-up mass m  : {:.0} kg\n\
                 aspect ratio AR: {:.1}\n\
                 CD0 / e        : {:.3} / {:.2}\n\n\
                 weight W       : {:.0} N\n\
                 wing loading   : {:.0} N/m²\n\
                 stall speed Vs : {:.1} m/s ({:.0} km/h)\n\n\
                 at cruise v = {:.1} m/s:\n\
                 level CL       : {:.3}\n\
                 drag coeff CD  : {:.4}\n\
                 lift/drag L/D  : {:.1}\n\n\
                 CL at max L/D  : {:.3}\n\
                 max L/D        : {:.2}\n\
                 best glide     : {:.1} : 1\n\
                 glide range    : {:.1} km  (from {:.0} m)",
                s.wing_area_m2,
                s.mass_kg,
                s.aspect_ratio,
                s.cd0,
                s.oswald_efficiency,
                p.weight_n,
                p.wing_loading_pa,
                p.stall_speed_m_s,
                p.stall_speed_m_s * 3.6,
                vc,
                cl_cruise,
                cd_cruise,
                ld_cruise,
                p.cl_at_max_ld,
                p.max_lift_to_drag,
                a.best_glide_ratio(),
                range_km,
                s.glide_altitude_m,
            );
        }
        Err(e) => s.error = Some(e),
    }
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

/// Append a flat (double-sided) quad `a-b-c-d` to the buffers. Used for the
/// thin lifting surfaces, whose winding need not be tracked.
fn push_quad(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    a: Vector3<f64>,
    b: Vector3<f64>,
    c: Vector3<f64>,
    d: Vector3<f64>,
) {
    let base = nodes.len();
    nodes.push(a);
    nodes.push(b);
    nodes.push(c);
    nodes.push(d);
    tris.extend_from_slice(&[
        base,
        base + 1,
        base + 2,
        base,
        base + 2,
        base + 3,
        base,
        base + 2,
        base + 1,
        base,
        base + 3,
        base + 2,
    ]);
}

/// Build the aircraft as a triangle [`Mesh`] — a box fuselage plus tapered
/// main wings, a horizontal tail and a vertical fin. The planform follows
/// the aspect ratio and wing area (`span = sqrt(AR * S)`), so the geometry
/// tracks the inputs. `None` for an invalid (out-of-domain) aircraft.
/// `+x` is forward (nose), `+y` is the right wing, `+z` is up.
fn aircraft_solid_mesh(s: &FixedWingWorkbenchState) -> Option<Mesh> {
    let a = build_aircraft(s).ok()?;
    let span = (a.aspect_ratio * a.wing_area_m2).sqrt();
    let chord = a.wing_area_m2 / span;
    let semi = span * 0.5;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Fuselage: an elongated box with real height + width (the 3-D anchor).
    let fus_len = span * 0.85;
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::zeros(),
        Vector3::new(fus_len * 0.5, span * 0.05, span * 0.07),
    );

    // Main wing: tapered, swept slightly aft, at mid-fuselage. Root chord
    // 1.5x mean chord, tip chord 0.55x, both centred near the CG.
    let root_c = chord * 1.5;
    let tip_c = chord * 0.55;
    let root_xc = 0.0;
    let tip_xc = -chord * 0.4; // gentle aft sweep
    let zw = 0.0;
    for side in [1.0_f64, -1.0] {
        let yr = side * span * 0.05; // root at fuselage side
        let yt = side * semi;
        push_quad(
            &mut nodes,
            &mut tris,
            Vector3::new(root_xc + root_c * 0.5, yr, zw),
            Vector3::new(root_xc - root_c * 0.5, yr, zw),
            Vector3::new(tip_xc - tip_c * 0.5, yt, zw),
            Vector3::new(tip_xc + tip_c * 0.5, yt, zw),
        );
    }

    // Horizontal tailplane near the rear.
    let tail_x = -fus_len * 0.46;
    let th_root = chord * 0.7;
    let th_tip = chord * 0.35;
    let tail_semi = semi * 0.4;
    for side in [1.0_f64, -1.0] {
        let yt = side * tail_semi;
        push_quad(
            &mut nodes,
            &mut tris,
            Vector3::new(tail_x + th_root * 0.5, 0.0, span * 0.04),
            Vector3::new(tail_x - th_root * 0.5, 0.0, span * 0.04),
            Vector3::new(tail_x - th_tip * 0.5, yt, span * 0.04),
            Vector3::new(tail_x + th_tip * 0.5, yt, span * 0.04),
        );
    }

    // Vertical fin (in the x-z plane at y = 0), swept.
    let fin_h = span * 0.18;
    push_quad(
        &mut nodes,
        &mut tris,
        Vector3::new(tail_x + th_root * 0.5, 0.0, span * 0.07),
        Vector3::new(tail_x - th_root * 0.5, 0.0, span * 0.07),
        Vector3::new(tail_x - th_root * 0.4, 0.0, span * 0.07 + fin_h),
        Vector3::new(tail_x + th_root * 0.1, 0.0, span * 0.07 + fin_h),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-fixedwing");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D aircraft solid and load it into the central viewport.
fn load_aircraft_3d(app: &mut ValenxApp) {
    let Some(mesh) = aircraft_solid_mesh(&app.fixedwing) else {
        app.fixedwing.error =
            Some("aircraft parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<aircraft>/valenx-fixedwing"),
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
        let s = FixedWingWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_aircraft_reports_glide_and_stall() {
        let mut s = FixedWingWorkbenchState::default();
        run_aircraft(&mut s);
        assert!(s.error.is_none());
        assert!(s.result.contains("stall speed"));
        assert!(s.result.contains("max L/D"));
        assert!(s.result.contains("glide range"));
    }

    #[test]
    fn analyze_rejects_out_of_domain_inputs() {
        let mut zero_area = FixedWingWorkbenchState {
            wing_area_m2: 0.0,
            ..Default::default()
        };
        run_aircraft(&mut zero_area);
        assert!(zero_area.error.is_some());
        let mut bad_e = FixedWingWorkbenchState {
            oswald_efficiency: 1.5,
            ..Default::default()
        };
        run_aircraft(&mut bad_e);
        assert!(bad_e.error.is_some());
    }

    #[test]
    fn aircraft_mesh_for_default_is_nonempty_and_in_range() {
        let s = FixedWingWorkbenchState::default();
        let mesh = aircraft_solid_mesh(&s).expect("default aircraft yields a solid");
        assert!(
            mesh.nodes.len() > 8,
            "expected fuselage + wings + tail surfaces"
        );
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn aircraft_mesh_none_for_invalid() {
        let s = FixedWingWorkbenchState {
            wing_area_m2: 0.0,
            ..Default::default()
        };
        assert!(aircraft_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_fixedwing_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_fixedwing_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_fixedwing_workbench = true;
        run_aircraft(&mut app.fixedwing);
        draw_workbench(&mut app);
    }
}
