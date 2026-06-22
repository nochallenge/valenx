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

use valenx_fixedwing::{Aircraft, SEA_LEVEL_AIR_DENSITY};
use valenx_mesh::Mesh;

use crate::mesh_prims::MeshBuilder;
use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Aluminium airframe grey for the fuselage body.
const FUSELAGE: [f32; 3] = [0.74, 0.75, 0.78];
/// Slightly cooler grey for the main wings.
const WING: [f32; 3] = [0.62, 0.65, 0.70];
/// Empennage (tail surfaces) — a touch darker so they read apart from the wings.
const TAIL: [f32; 3] = [0.52, 0.55, 0.60];

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

/// Half the points around one closed airfoil section. The section is sampled
/// as upper-then-lower surface so a [`MeshBuilder::loft`] of two of these skins
/// a real wing with finite thickness.
const AIRFOIL_HALF_POINTS: usize = 16;

/// A closed **symmetric airfoil** loop for a NACA-`00tt` section of chord
/// `chord` and max-thickness fraction `t_frac` (e.g. `0.12` ⇒ 12 % thick),
/// returned as `[chordwise, thickness]` 2-D points. Built from the classic
/// 4-digit half-thickness distribution
/// `yt = 5t(0.2969√x − 0.1260x − 0.3516x² + 0.2843x³ − 0.1015x⁴)` (chord
/// normalised to 1, scaled by `chord`), sampled with a cosine chordwise
/// spacing so the curved leading edge is well resolved. The loop runs along the
/// upper surface from the trailing edge to the leading edge, then back along
/// the lower surface — a single closed counter-clockwise polygon. `cx` is the
/// chordwise origin (the section's local x), so a swept section just shifts it.
fn airfoil_section(chord: f64, t_frac: f64, cx: f64) -> Vec<[f64; 2]> {
    let yt = |x: f64| {
        5.0 * t_frac
            * (0.2969 * x.sqrt() - 0.1260 * x - 0.3516 * x * x + 0.2843 * x * x * x
                - 0.1015 * x * x * x * x)
    };
    let n = AIRFOIL_HALF_POINTS;
    let mut pts = Vec::with_capacity(2 * n);
    // Upper surface: trailing edge (x=1) → leading edge (x=0), cosine spacing.
    for i in 0..n {
        let beta = i as f64 / (n - 1) as f64 * std::f64::consts::PI;
        let x = 0.5 * (1.0 + beta.cos()); // 1 → 0
        pts.push([cx + x * chord, yt(x) * chord]);
    }
    // Lower surface: leading edge (x=0) → trailing edge (x=1), skipping the
    // shared LE/TE endpoints so the polygon has no duplicate vertices.
    for i in 1..n - 1 {
        let beta = i as f64 / (n - 1) as f64 * std::f64::consts::PI;
        let x = 0.5 * (1.0 - beta.cos()); // 0 → 1
        pts.push([cx + x * chord, -yt(x) * chord]);
    }
    pts
}

/// Place a planar `[chordwise, thickness]` airfoil section into 3-D for a
/// **horizontal** lifting surface (wing / tailplane): chord runs along x,
/// thickness along z, and the whole section sits at spanwise station `y` and
/// height `z0`. The result is one ring ready for [`MeshBuilder::loft`].
fn horiz_section(profile: &[[f64; 2]], y: f64, z0: f64) -> Vec<[f64; 3]> {
    profile.iter().map(|&[cx, th]| [cx, y, z0 + th]).collect()
}

/// Place a planar `[chordwise, thickness]` airfoil section into 3-D for the
/// **vertical** fin: chord runs along x, thickness along y, and the section
/// sits at height station `z` (so lofting root→tip sweeps it upward).
fn vert_section(profile: &[[f64; 2]], z: f64) -> Vec<[f64; 3]> {
    profile.iter().map(|&[cx, th]| [cx, th, z]).collect()
}

/// Build the aircraft as a triangle [`Mesh`] **with per-vertex colours** — a
/// **revolved tapered fuselage** (nose + tail fairing) plus **extruded airfoil
/// surfaces** for the main wings, the horizontal tailplane and the vertical
/// fin, each a real thickness/taper/sweep loft of a NACA-`00tt` section (not a
/// flat zero-thickness quad). The planform follows the aspect ratio and wing
/// area (`span = √(AR·S)`), so the geometry tracks the inputs. Colours:
/// fuselage aluminium grey, wings cool grey, empennage darker. `+x` is forward
/// (nose), `+y` is the right wing, `+z` is up. `None` for an invalid aircraft.
///
/// Returns `(mesh, colors)` with `colors.len() == 3 × triangle_count`, ready
/// for [`crate::WorkspaceProduct::vertex_colors`].
fn aircraft_solid_mesh_colored(s: &FixedWingWorkbenchState) -> Option<(Mesh, Vec<[f32; 3]>)> {
    let a = build_aircraft(s).ok()?;
    let span = (a.aspect_ratio * a.wing_area_m2).sqrt();
    let chord = a.wing_area_m2 / span;
    let semi = span * 0.5;
    let fus_len = span * 0.85;
    let max_r = span * 0.06; // peak fuselage radius

    let mut b = MeshBuilder::new();

    // Fuselage: a revolved body of revolution along +x. The half-profile is
    // `[radius, x]` from a pointed nose, swelling to the cabin, then fairing to
    // a small tailcone — a real tapered fuselage instead of a box.
    let nose = fus_len * 0.5;
    let tail = -fus_len * 0.5;
    let profile: Vec<[f64; 2]> = vec![
        [0.0, nose], // nose apex
        [max_r * 0.45, nose - fus_len * 0.10],
        [max_r * 0.85, nose - fus_len * 0.22],
        [max_r, fus_len * 0.10], // cabin (fullest, slightly fwd of centre)
        [max_r * 0.95, 0.0],
        [max_r * 0.70, tail + fus_len * 0.28],
        [max_r * 0.34, tail + fus_len * 0.06],
        [0.0, tail], // tailcone tip
    ];
    b.revolve(
        &profile,
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        360.0,
        24,
        FUSELAGE,
    );

    // Main wing: tapered + swept, finite-thickness airfoil. Root chord 1.5×
    // mean, tip 0.55×, gentle aft sweep, 12 % thick. Lofted root→tip per side.
    let root_c = chord * 1.5;
    let tip_c = chord * 0.55;
    let tip_sweep = -chord * 0.4;
    let root_y = span * 0.05;
    for side in [1.0_f64, -1.0] {
        let root = horiz_section(
            &airfoil_section(root_c, 0.12, -root_c * 0.5),
            side * root_y,
            0.0,
        );
        let tip = horiz_section(
            &airfoil_section(tip_c, 0.10, tip_sweep - tip_c * 0.5),
            side * semi,
            0.0,
        );
        // Order sections outboard-going so the cap winding faces outward.
        let secs = if side > 0.0 {
            vec![root, tip]
        } else {
            vec![tip, root]
        };
        b.loft(&secs, true, WING);
    }

    // Horizontal tailplane near the rear, smaller + thinner.
    let tail_x = -fus_len * 0.46;
    let th_root = chord * 0.7;
    let th_tip = chord * 0.35;
    let tail_semi = semi * 0.4;
    let tail_z = max_r * 0.4;
    for side in [1.0_f64, -1.0] {
        let root = horiz_section(
            &airfoil_section(th_root, 0.10, tail_x - th_root * 0.5),
            0.0,
            tail_z,
        );
        let tip = horiz_section(
            &airfoil_section(th_tip, 0.09, tail_x - chord * 0.15 - th_tip * 0.5),
            side * tail_semi,
            tail_z,
        );
        let secs = if side > 0.0 {
            vec![root, tip]
        } else {
            vec![tip, root]
        };
        b.loft(&secs, true, TAIL);
    }

    // Vertical fin in the x-z plane at y = 0, swept, lofted bottom→top.
    let fin_h = span * 0.18;
    let fin_root_c = th_root;
    let fin_tip_c = th_root * 0.55;
    let fin_root = vert_section(
        &airfoil_section(fin_root_c, 0.10, tail_x - fin_root_c * 0.5),
        max_r * 0.6,
    );
    let fin_tip = vert_section(
        &airfoil_section(fin_tip_c, 0.09, tail_x - chord * 0.2 - fin_tip_c * 0.5),
        max_r * 0.6 + fin_h,
    );
    b.loft(&[fin_root, fin_tip], true, TAIL);

    let (mut mesh, colors) = b.into_mesh_and_colors();
    mesh.id = "valenx-fixedwing".to_string();
    Some((mesh, colors))
}

/// Build the aircraft [`Mesh`] (without the colour metadata) for the central
/// viewport. See [`aircraft_solid_mesh_colored`].
fn aircraft_solid_mesh(s: &FixedWingWorkbenchState) -> Option<Mesh> {
    aircraft_solid_mesh_colored(s).map(|(mesh, _colors)| mesh)
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

/// The agent-bridge **`show_3d{kind:"fixedwing"}`** product: the canonical
/// fuselage + wings + tail aircraft solid (the panel's "Show 3-D aircraft"
/// geometry) paired with the workbench's own point-performance headline
/// numbers, at a fixed 3/4 camera. Registered in
/// [`crate::products_registry`]; the per-tool builder the registry dispatches
/// to. Pure — driven off [`FixedWingWorkbenchState::default`].
///
/// The readout rows mirror the panel's `run_aircraft` "Point performance"
/// readout (this workbench formats its readout into `result` rather than via a
/// shared `compute()` string fn, so the rows are taken from that here).
pub(crate) fn fixedwing_product() -> crate::WorkspaceProduct {
    let mut s = FixedWingWorkbenchState::default();
    run_aircraft(&mut s);
    let (mesh, colors) =
        aircraft_solid_mesh_colored(&s).expect("default GA aircraft ⇒ a 3-D solid");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<aircraft>/valenx-fixedwing");
    let lines = crate::products_registry::lines_from_readout(&s.result);
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Fixed-Wing / Aircraft".into(),
        lines,
        mesh: Some(loaded),
        vertex_colors: Some(colors),
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

    #[test]
    fn airfoil_section_is_a_closed_thick_loop() {
        // A 12 % NACA-00xx section: a closed loop of 2·n−2 vertices spanning the
        // full chord, with a real (non-zero) max thickness ≈ 12 % of chord, and
        // upper/lower surfaces symmetric about the chord line.
        let chord = 2.0;
        let pts = airfoil_section(chord, 0.12, 0.0);
        assert_eq!(pts.len(), 2 * AIRFOIL_HALF_POINTS - 2);
        let max_x = pts.iter().map(|p| p[0]).fold(f64::MIN, f64::max);
        let min_x = pts.iter().map(|p| p[0]).fold(f64::MAX, f64::min);
        assert!(
            (max_x - chord).abs() < 1e-9 && min_x.abs() < 1e-9,
            "spans chord"
        );
        let max_t = pts.iter().map(|p| p[1]).fold(f64::MIN, f64::max);
        // NACA 4-digit peak half-thickness is ≈ 0.6·t·chord above the chord.
        assert!(
            (max_t - 0.12 * chord * 0.6).abs() < 0.02 * chord,
            "≈12% thickness, got {max_t}"
        );
        let min_t = pts.iter().map(|p| p[1]).fold(f64::MAX, f64::min);
        assert!((max_t + min_t).abs() < 1e-9, "symmetric about the chord");
    }

    #[test]
    fn aircraft_carries_vertex_aligned_colours() {
        // The revolved fuselage + lofted airfoil surfaces ship per-vertex
        // colours aligned to the renderer's coloured path (3 / triangle), with
        // the fuselage, wing and tail colours all present and a real 3-D solid
        // (the wings now have thickness, so node count far exceeds the old
        // flat-quad build).
        let s = FixedWingWorkbenchState::default();
        let (mesh, colors) =
            aircraft_solid_mesh_colored(&s).expect("default aircraft builds coloured");
        assert!(!mesh.nodes.is_empty(), "non-empty mesh");
        assert!(mesh.total_elements() > 0, "mesh has triangles");
        assert_eq!(
            colors.len(),
            mesh.total_elements() * 3,
            "vertex_colors must equal 3 × triangle count"
        );
        assert!(colors.contains(&FUSELAGE), "fuselage colour present");
        assert!(colors.contains(&WING), "wing colour present");
        assert!(colors.contains(&TAIL), "tail colour present");
        for c in &colors {
            for ch in c {
                assert!(ch.is_finite() && (0.0..=1.0).contains(ch));
            }
        }
    }

    #[test]
    fn aircraft_product_is_coloured_and_aligned() {
        let product = fixedwing_product();
        let loaded = product.mesh.as_ref().expect("aircraft product has a mesh");
        let colors = product
            .vertex_colors
            .as_ref()
            .expect("aircraft product carries vertex_colors");
        assert_eq!(
            colors.len(),
            loaded.mesh.total_elements() * 3,
            "product colours aligned to the coloured path"
        );
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
