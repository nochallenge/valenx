//! The right-side **Belt Drive Workbench** panel — native open flat-belt
//! drive analysis over `valenx-beltdrive`.
//!
//! Mirrors the Gearbox / Heat Transfer workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_beltdrive_workbench`,
//! toggled from the View menu. The form sets the driver and driven pulley
//! diameters, the driver speed, the belt/pulley friction coefficient, the
//! shaft centre distance and the belt linear density / maximum tension;
//! "Analyze" reports the speed ratio, belt linear speed, driven speed, the
//! open-belt wrap angles, the capstan tension ratio `T1/T2 = exp(mu*theta)`,
//! the centrifugal tension and the slipping-limited power, and "Show 3-D
//! belt drive" loads a representative two-pulley drive (two discs on
//! parallel shafts joined by a belt loop) into the central viewport.

use std::f64::consts::{PI, TAU};
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_beltdrive::{belt_length_open, DriveAnalysis, DriveSpec};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::mesh_prims::MeshBuilder;
use crate::pulley_workbench::push_grooved_sheave;
use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Cast-iron sheave grey (matches the Pulley Workbench).
const PULLEY_GREY: [f32; 3] = [0.46, 0.47, 0.50];
/// Dark hub/axle.
const HUB: [f32; 3] = [0.26, 0.27, 0.30];
/// Dark rubber belt.
const RUBBER: [f32; 3] = [0.13, 0.13, 0.15];

/// Persistent form + result state for the Belt Drive Workbench.
pub struct BeltDriveWorkbenchState {
    /// Driver (input) pulley pitch diameter (mm).
    driver_diameter_mm: f64,
    /// Driven (output) pulley pitch diameter (mm).
    driven_diameter_mm: f64,
    /// Driver rotational speed (rpm).
    driver_speed_rpm: f64,
    /// Belt/pulley coefficient of friction (dimensionless).
    mu: f64,
    /// Shaft-to-shaft centre distance (mm).
    center_distance_mm: f64,
    /// Belt linear mass density (kg/m).
    linear_density: f64,
    /// Maximum allowable tight-side tension (N).
    t1_max: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D belt-drive solid (serviced after
    /// the panel draws).
    show_3d_request: bool,
}

impl Default for BeltDriveWorkbenchState {
    fn default() -> Self {
        // A 100 mm driver at 1450 rpm driving a 250 mm pulley 500 mm away
        // with mu = 0.3: a 2.5:1 reduction, ~7.6 m/s belt, capstan ratio
        // ~2.35 on the small pulley.
        Self {
            driver_diameter_mm: 100.0,
            driven_diameter_mm: 250.0,
            driver_speed_rpm: 1450.0,
            mu: 0.3,
            center_distance_mm: 500.0,
            linear_density: 0.4,
            t1_max: 1200.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Belt Drive Workbench right-side panel. A no-op when the
/// `show_beltdrive_workbench` toggle is off.
pub fn draw_beltdrive_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_beltdrive_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_beltdrive_workbench",
        "Belt Drive",
        |app, ui| {
            ui.label(
                egui::RichText::new("native open flat-belt drive analysis · valenx-beltdrive")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.beltdrive;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Pulleys").strong());
                    ui.horizontal(|ui| {
                        let l = ui.label("driver Ø (mm)");
                        ui.add(egui::DragValue::new(&mut s.driver_diameter_mm).speed(1.0))
                            .labelled_by(l.id);
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("driven Ø (mm)");
                        ui.add(egui::DragValue::new(&mut s.driven_diameter_mm).speed(1.0))
                            .labelled_by(l.id);
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("centre distance (mm)");
                        ui.add(egui::DragValue::new(&mut s.center_distance_mm).speed(5.0))
                            .labelled_by(l.id);
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Drive").strong());
                    ui.horizontal(|ui| {
                        let l = ui.label("driver speed (rpm)");
                        ui.add(egui::DragValue::new(&mut s.driver_speed_rpm).speed(5.0))
                            .labelled_by(l.id);
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("friction μ");
                        ui.add(egui::DragValue::new(&mut s.mu).speed(0.01))
                            .labelled_by(l.id);
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Belt").strong());
                    ui.horizontal(|ui| {
                        let l = ui.label("linear density (kg/m)");
                        ui.add(egui::DragValue::new(&mut s.linear_density).speed(0.01))
                            .labelled_by(l.id);
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("max tension T1 (N)");
                        ui.add(egui::DragValue::new(&mut s.t1_max).speed(10.0))
                            .labelled_by(l.id);
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_beltdrive(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D belt drive").strong())
                        .on_hover_text(
                            "Build a representative open belt drive (driver and driven discs on two parallel shafts joined by a belt loop) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Drive").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_beltdrive_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.beltdrive` borrow is
    // released here): build the belt drive's 3-D solid and load it.
    if app.beltdrive.show_3d_request {
        app.beltdrive.show_3d_request = false;
        load_beltdrive_3d(app);
    }
}

/// Validate the form, evaluate the drive and format the readout.
fn run_beltdrive(s: &mut BeltDriveWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build a validated [`DriveSpec`] from the form (converting the mm pulley
/// dimensions to metres and the rpm to rev/s), mapping any domain error to
/// a display string. Extracted so it is shared by the readout and the 3-D
/// gate.
fn build_spec(s: &BeltDriveWorkbenchState) -> Result<DriveSpec, String> {
    Ok(DriveSpec {
        driver_diameter: s.driver_diameter_mm / 1000.0,
        driven_diameter: s.driven_diameter_mm / 1000.0,
        driver_rev_per_sec: s.driver_speed_rpm / 60.0,
        center_distance: s.center_distance_mm / 1000.0,
        mu: s.mu,
        linear_density: s.linear_density,
        t1_max: s.t1_max,
    })
}

/// Evaluate the drive and format the full readout, mapping any domain error
/// to a display string. Extracted so it is unit-testable.
fn compute(s: &BeltDriveWorkbenchState) -> Result<String, String> {
    let spec = build_spec(s)?;
    let a: DriveAnalysis = spec.analyze().map_err(|e| e.to_string())?;

    let driven_rpm = a.driven_rev_per_sec * 60.0;
    // Exact open-belt length (two straight tangents + the two wrapped
    // arcs): the practical "what length of belt to order" figure. The
    // free function sorts the small/large radii internally, so the
    // driver/driven order does not matter.
    let belt_length = belt_length_open(
        spec.driver_diameter / 2.0,
        spec.driven_diameter / 2.0,
        spec.center_distance,
    )
    .map_err(|e| e.to_string())?;

    Ok(format!(
        "driver / driven : {:.0} / {:.0} mm\n\
         centre distance : {:.0} mm\n\
         friction μ      : {:.2}\n\n\
         speed ratio     : {:.2} : 1\n\
         belt speed      : {:.2} m/s\n\
         driver / driven : {:.0} / {:.0} rpm\n\n\
         wrap (small)    : {:.1}°\n\
         wrap (large)    : {:.1}°\n\
         belt length     : {belt_length:.3} m ({:.0} mm)\n\n\
         capstan T1/T2   : {:.3}\n\
         centrifugal Tc  : {:.1} N\n\
         max power       : {:.3} kW",
        s.driver_diameter_mm,
        s.driven_diameter_mm,
        s.center_distance_mm,
        s.mu,
        a.speed_ratio,
        a.belt_speed,
        s.driver_speed_rpm,
        driven_rpm,
        a.wrap_small.to_degrees(),
        a.wrap_large.to_degrees(),
        belt_length * 1000.0,
        a.tension_ratio,
        a.centrifugal_tension,
        a.max_power / 1000.0,
    ))
}

/// Presentation spin rate of the driver pulley, rad/s (~1.0 rev/s) — a readable
/// inspect speed; the driven pulley follows at the diameter ratio.
const DRIVER_RAD_PER_S: f32 = 6.0;

/// Build the **open (external-tangent) belt** as its own thin swept-band
/// [`Mesh`]. The belt wraps two pulleys whose centres sit in the wheel plane
/// (the constant-`x` y–z plane) at `c1 = (y, z1)` / `c2 = (y, z2)` with belt
/// pitch radii `rho1` / `rho2` (the pulley radius + the belt's radial
/// mid-thickness). The closed belt centreline is the standard open-belt path —
/// a wrap arc on pulley 1, a straight external tangent, a wrap arc on pulley 2,
/// the return tangent — built from the closed-form tangent angle
/// `gamma = asin((rho2 − rho1) / C)`. The centreline is then swept into a solid
/// of `width` (along the axle `x`) and `2·half_t` radial thickness.
///
/// Returns a [`valenx_mesh::Mesh`] (`Tri3`) ready to drop into a [`MeshBuilder`]
/// via [`MeshBuilder::append_tri_mesh`].
#[allow(clippy::too_many_arguments)]
fn belt_band_mesh(
    cy: f64,
    z1: f64,
    rho1: f64,
    z2: f64,
    rho2: f64,
    belt_x: f64,
    width: f64,
    half_t: f64,
) -> Mesh {
    // Centre-to-centre geometry in the wheel plane (axes: z horizontal, y up).
    let dz = z2 - z1;
    let dy = 0.0; // both centres share the same height
    let c_dist = (dz * dz + dy * dy).sqrt().max(1e-9);
    // Direction from pulley 1 → pulley 2 and the in-plane perpendicular.
    let dir = (dz / c_dist, dy / c_dist); // (z, y)
    let perp = (-dir.1, dir.0); // rotate +90° in the (z, y) plane
                                // Open-belt tangent offset angle.
    let gamma = ((rho2 - rho1) / c_dist).clamp(-1.0, 1.0).asin();

    // A point on pulley `i` (centre `(zc, yc)`, radius `rho`) at angle `phi`
    // measured from `dir`, toward `perp`. Returns (z, y).
    let on_circle = |zc: f64, yc: f64, rho: f64, phi: f64| -> (f64, f64) {
        let (s, c) = phi.sin_cos();
        // Local frame: dir = phi 0, perp = phi +90°.
        let z = zc + rho * (c * dir.0 + s * perp.0);
        let y = yc + rho * (c * dir.1 + s * perp.1);
        (z, y)
    };

    // Tangent contact angles (from `dir`): the upper tangent leaves pulley 1 and
    // meets pulley 2 at +(π/2 + gamma); the lower at −(π/2 + gamma). The belt
    // wraps pulley 1 over the MAJOR arc on the far side (from +(π/2+γ) round the
    // back to −(π/2+γ)), and pulley 2 over its arc on the near side.
    let a = PI / 2.0 + gamma; // upper contact, measured from dir
                              // Sample: pulley-1 wrap (the long way round, away from pulley 2), then
                              // tangent to pulley-2 lower contact, pulley-2 wrap (the long way round away
                              // from pulley 1), then tangent back. Build the closed centreline polyline.
    let arc_steps = 24usize;
    let mut center: Vec<(f64, f64)> = Vec::new();

    // Pulley 1 wrap: from +a, increasing through π to (2π − a) i.e. the arc on
    // the side AWAY from pulley 2 (pulley 2 is at phi = 0 direction).
    let p1_start = a;
    let p1_end = TAU - a;
    for k in 0..=arc_steps {
        let t = k as f64 / arc_steps as f64;
        let phi = p1_start + (p1_end - p1_start) * t;
        center.push(on_circle(z1, cy, rho1, phi));
    }
    // Pulley 2 wrap: pulley 1 is in the −dir direction from pulley 2, i.e. at
    // phi = π relative to pulley 2's own `dir`. The belt wraps pulley 2 on the
    // far side: from −a round through 0 to +a (centred on phi = 0, the side away
    // from pulley 1). Use the SAME `dir`/`perp`; pulley 2's far side is phi ∈
    // [−a, +a] going through 0.
    let p2_start = -a;
    let p2_end = a;
    for k in 0..=arc_steps {
        let t = k as f64 / arc_steps as f64;
        let phi = p2_start + (p2_end - p2_start) * t;
        center.push(on_circle(z2, cy, rho2, phi));
    }

    // Sweep the closed centreline into a solid band. Each centreline point gets
    // a 4-corner cross-section: ±half_t radially (in the wheel plane, along the
    // outward normal from the local pulley centre ≈ the centreline tangent's
    // perpendicular) and ±0.5·width along the axle (x). We approximate the
    // radial outward direction by the segment normal in the (z, y) plane.
    let n = center.len();
    let mut mesh = Mesh::new("valenx-belt");
    let mut block = ElementBlock::new(ElementType::Tri3);
    let hw = 0.5 * width;
    // Build 4 corners per station: order [inner-(-x), outer-(-x), outer+(+x), inner+(+x)].
    for i in 0..n {
        let (z, y) = center[i];
        let (zp, yp) = center[(i + 1) % n];
        let (zm, ym) = center[(i + n - 1) % n];
        // Tangent ≈ next − prev; radial normal = perp of tangent in (z, y).
        let (tz, ty) = (zp - zm, yp - ym);
        let tl = (tz * tz + ty * ty).sqrt().max(1e-9);
        let (nz, ny) = (-ty / tl, tz / tl); // outward-ish normal in plane
        let base = mesh.nodes.len() as u32;
        // inner (−normal), −x
        mesh.nodes
            .push(Vector3::new(belt_x - hw, y - ny * half_t, z - nz * half_t));
        // outer (+normal), −x
        mesh.nodes
            .push(Vector3::new(belt_x - hw, y + ny * half_t, z + nz * half_t));
        // outer (+normal), +x
        mesh.nodes
            .push(Vector3::new(belt_x + hw, y + ny * half_t, z + nz * half_t));
        // inner (−normal), +x
        mesh.nodes
            .push(Vector3::new(belt_x + hw, y - ny * half_t, z - nz * half_t));
        let _ = base;
    }
    // Connect consecutive cross-sections (i → i+1, wrapping) into a tube of 4
    // quad faces (inner, outer, −x side, +x side) = 8 triangles per segment.
    let quad = |block: &mut ElementBlock, a: u32, b: u32, c: u32, d: u32| {
        block.connectivity.extend_from_slice(&[a, b, c, a, c, d]);
    };
    for i in 0..n {
        let s0 = (i * 4) as u32;
        let s1 = (((i + 1) % n) * 4) as u32;
        // corners: 0 inner−x, 1 outer−x, 2 outer+x, 3 inner+x
        // −x face (corners 0,1)
        quad(&mut block, s0, s0 + 1, s1 + 1, s1);
        // outer face (corners 1,2)
        quad(&mut block, s0 + 1, s0 + 2, s1 + 2, s1 + 1);
        // +x face (corners 2,3)
        quad(&mut block, s0 + 2, s0 + 3, s1 + 3, s1 + 2);
        // inner face (corners 3,0)
        quad(&mut block, s0 + 3, s0, s1, s1 + 3);
    }
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    mesh
}

/// Build the belt drive as a triangle [`Mesh`] **with per-vertex colours** plus
/// a [`crate::RigidPart`] for each pulley (driver + driven), so both wheels spin
/// about their shafts while the belt + base stay put. Both pulleys are now true
/// **grooved V-sheaves** ([`push_grooved_sheave`], shared with the Pulley
/// Workbench) and the belt is the real **external-tangent loop** wrapping them
/// ([`belt_band_mesh`]) — not two parked boxes. The driver turns the same
/// direction as the driven but the driven is slower by the diameter ratio (open
/// belt). Colours: pulleys cast-iron grey, belt dark rubber. `None` for an
/// invalid config.
///
/// Returns `(mesh, colors, parts)` with `colors.len() == 3 × triangle_count`.
fn beltdrive_solid_mesh_parts(
    s: &BeltDriveWorkbenchState,
) -> Option<(Mesh, Vec<[f32; 3]>, Vec<crate::RigidPart>)> {
    // Gate the geometry on a buildable, analysable drive.
    build_spec(s).ok()?.analyze().ok()?;

    // Scale the real (mm) dimensions to a tidy on-screen size.
    let scale = 1.0 / s.center_distance_mm.max(1.0);
    let r_drv = 0.5 * s.driver_diameter_mm * scale;
    let r_drn = 0.5 * s.driven_diameter_mm * scale;
    let half_c = 0.5; // centre distance maps to 1.0 along z
    let face_x = -0.06;
    let width = 0.12;
    let half_w = 0.5 * width;
    let bore_r = 0.04;
    let lift = r_drv.max(r_drn) + 0.15;

    let mut b = MeshBuilder::new();

    // Driver grooved sheave (small) at z = −half_c — record its node range.
    let drv_centre = [face_x, lift, -half_c];
    let drv_range = push_grooved_sheave(
        &mut b,
        drv_centre,
        [1.0, 0.0, 0.0],
        r_drv,
        half_w,
        bore_r,
        PULLEY_GREY,
        HUB,
    );
    // Driven grooved sheave (large) at z = +half_c — record its node range.
    let drn_centre = [face_x, lift, half_c];
    let drn_range = push_grooved_sheave(
        &mut b,
        drn_centre,
        [1.0, 0.0, 0.0],
        r_drn,
        half_w,
        bore_r,
        PULLEY_GREY,
        HUB,
    );

    // Real open-belt loop wrapping both sheave rims (centreline radius = rim +
    // half the belt thickness), swept into a thin band. Static.
    let belt_x = face_x; // belt sits centred on the sheave plane
    let belt_half_t = 0.02;
    let belt_band = belt_band_mesh(
        lift,
        -half_c,
        r_drv + belt_half_t,
        half_c,
        r_drn + belt_half_t,
        belt_x,
        width * 0.9,
        belt_half_t,
    );
    b.append_tri_mesh(&belt_band, RUBBER);

    // Base (static).
    b.cuboid([0.0, 0.06, 0.0], [1.0, 0.12, 1.6], HUB);

    let (mut mesh, colors) = b.into_mesh_and_colors();
    mesh.id = "valenx-beltdrive".to_string();

    // Open belt ⇒ same direction; the driven pulley is slower by the diameter
    // ratio (ω_driven = ω_driver · d_driver / d_driven). Presentation-scaled.
    let ratio = if r_drn > 0.0 {
        (r_drv / r_drn) as f32
    } else {
        1.0
    };
    let parts = vec![
        crate::RigidPart {
            node_range: drv_range,
            axis: [1.0, 0.0, 0.0],
            pivot: [
                drv_centre[0] as f32,
                drv_centre[1] as f32,
                drv_centre[2] as f32,
            ],
            rad_per_s: DRIVER_RAD_PER_S,
        },
        crate::RigidPart {
            node_range: drn_range,
            axis: [1.0, 0.0, 0.0],
            pivot: [
                drn_centre[0] as f32,
                drn_centre[1] as f32,
                drn_centre[2] as f32,
            ],
            rad_per_s: DRIVER_RAD_PER_S * ratio,
        },
    ];
    Some((mesh, colors, parts))
}

/// Build the belt-drive [`Mesh`] (without the colour / part metadata) for the
/// central viewport. See [`beltdrive_solid_mesh_parts`].
fn beltdrive_solid_mesh(s: &BeltDriveWorkbenchState) -> Option<Mesh> {
    beltdrive_solid_mesh_parts(s).map(|(mesh, _colors, _parts)| mesh)
}

/// Build the 3-D belt-drive solid and load it into the central viewport.
fn load_beltdrive_3d(app: &mut ValenxApp) {
    let Some(mesh) = beltdrive_solid_mesh(&app.beltdrive) else {
        app.beltdrive.error =
            Some("belt-drive parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<belt-drive>/valenx-beltdrive"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// Agent-bridge product: the canonical belt-drive workbench as a 3-D solid plus
/// its `compute()` readout rows (see [`crate::products_registry`]).
pub(crate) fn beltdrive_product() -> crate::WorkspaceProduct {
    let s = BeltDriveWorkbenchState::default();
    let (mesh, colors, parts) =
        beltdrive_solid_mesh_parts(&s).expect("canonical belt drive ⇒ drive solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<beltdrive>/valenx-beltdrive");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical belt drive ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Belt drive (tension/power)".into(),
        lines,
        mesh: Some(loaded),
        vertex_colors: Some(colors),
        camera,
        kind2d: None,
        last_export: None,
        image: None,
        image_texture: None,
        // Animated: both pulleys spin about their shafts (driven slower by the
        // diameter ratio) while the belt runs and base stay put. Paused at t = 0.
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
    use std::f64::consts::PI;

    #[test]
    fn default_state_is_idle() {
        let s = BeltDriveWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_ratio_speed_and_tension() {
        let mut s = BeltDriveWorkbenchState::default();
        run_beltdrive(&mut s);
        assert!(
            s.error.is_none(),
            "default belt drive should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("speed ratio"));
        assert!(s.result.contains("belt speed"));
        assert!(s.result.contains("capstan T1/T2"));
        // 250 mm driven on a 100 mm driver => 2.5:1 reduction.
        assert!(s.result.contains("2.50 : 1"));
        // Belt rim speed pi * 0.1 m * (1450/60) rev/s ~ 7.59 m/s.
        assert!(s.result.contains("7.59 m/s"));
    }

    #[test]
    fn analyze_rejects_zero_diameter() {
        let mut s = BeltDriveWorkbenchState {
            driver_diameter_mm: 0.0,
            ..Default::default()
        };
        run_beltdrive(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn speed_ratio_is_driven_over_driver_ground_truth() {
        // Ground truth: with no slip the surface speeds match, so the
        // transmission ratio is exactly D_driven / D_driver, and the belt
        // rim speed is pi * D_driver * N_driver.
        let s = BeltDriveWorkbenchState::default();
        let a = build_spec(&s).unwrap().analyze().unwrap();
        let expected_ratio = s.driven_diameter_mm / s.driver_diameter_mm;
        assert!((a.speed_ratio - expected_ratio).abs() < 1e-12);
        let expected_v = PI * (s.driver_diameter_mm / 1000.0) * (s.driver_speed_rpm / 60.0);
        assert!((a.belt_speed - expected_v).abs() < 1e-9);
        // Capstan ratio on the small pulley: T1/T2 = exp(mu * theta_small).
        assert!((a.tension_ratio - (s.mu * a.wrap_small).exp()).abs() < 1e-9);
    }

    #[test]
    fn belt_length_is_exact_open_belt_ground_truth() {
        // Ground truth for the default 100/250 mm pulleys 500 mm apart:
        // the exact open-belt length is the two straight tangents plus the
        // two wrapped arcs,
        //   alpha = asin((R_large - R_small) / C),
        //   L = 2*C*cos(alpha) + R_small*(pi - 2*alpha) + R_large*(pi + 2*alpha).
        // With R_small = 0.05 m, R_large = 0.125 m, C = 0.5 m this is
        // 0.98868 + 0.14202 + 0.43034 = 1.56105 m, i.e. "1.561 m".
        let s = BeltDriveWorkbenchState::default();
        let r_small = s.driver_diameter_mm / 1000.0 / 2.0;
        let r_large = s.driven_diameter_mm / 1000.0 / 2.0;
        let c = s.center_distance_mm / 1000.0;
        let alpha = ((r_large - r_small) / c).asin();
        let expected =
            2.0 * c * alpha.cos() + r_small * (PI - 2.0 * alpha) + r_large * (PI + 2.0 * alpha);
        assert!(
            (expected - 1.561_049_951_958_976).abs() < 1e-12,
            "hand check drifted: {expected}"
        );

        // The crate's free function must reproduce that closed form.
        let l = belt_length_open(r_small, r_large, c).unwrap();
        assert!((l - expected).abs() < 1e-12, "belt_length_open: {l}");

        // …and the readout must surface it (3-decimal metres).
        let out = compute(&s).expect("default belt drive computes");
        assert!(
            out.contains("belt length     : 1.561 m"),
            "missing belt length line: {out}"
        );
        assert!(out.contains("1561 mm"), "missing mm form: {out}");
    }

    #[test]
    fn beltdrive_mesh_for_default_is_nonempty_and_in_range() {
        let s = BeltDriveWorkbenchState::default();
        let mesh = beltdrive_solid_mesh(&s).expect("default belt drive yields a solid");
        assert!(mesh.nodes.len() > 8, "expected two pulleys + belt + base");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn beltdrive_mesh_none_for_invalid() {
        // Centre distance too small for the radius difference => degenerate
        // open-belt geometry, so no 3-D solid.
        let s = BeltDriveWorkbenchState {
            center_distance_mm: 10.0,
            ..Default::default()
        };
        assert!(beltdrive_solid_mesh(&s).is_none());
    }

    #[test]
    fn beltdrive_carries_vertex_aligned_colours() {
        // The two grooved sheaves + the swept belt band ship per-vertex colours
        // aligned to the renderer's coloured path (3/triangle), with the pulley
        // grey and the dark rubber belt colours both present.
        let s = BeltDriveWorkbenchState::default();
        let (mesh, colors, _parts) =
            beltdrive_solid_mesh_parts(&s).expect("default belt drive builds");
        assert_eq!(
            colors.len(),
            mesh.total_elements() * 3,
            "vertex_colors must equal 3 × triangle count"
        );
        assert!(colors.contains(&PULLEY_GREY), "pulley rim colour present");
        assert!(colors.contains(&RUBBER), "belt colour present");
        for c in &colors {
            for ch in c {
                assert!(ch.is_finite() && (0.0..=1.0).contains(ch));
            }
        }
    }

    #[test]
    fn belt_band_is_a_closed_nonempty_loop() {
        // The open-belt band mesh is a closed swept loop: non-empty nodes +
        // triangles, all node indices in range. For unequal radii (100 vs 250
        // mm) the band still closes (the tangent angle is finite).
        let band = belt_band_mesh(0.7, -0.5, 0.12, 0.5, 0.27, -0.06, 0.1, 0.02);
        assert!(band.nodes.len() > 8, "belt band has vertices");
        assert!(band.total_elements() > 0, "belt band has triangles");
        let n = band.nodes.len() as u32;
        for blk in &band.element_blocks {
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn beltdrive_product_spins_both_pulleys_same_dir() {
        // The product carries a RigidParts animation: the driver and driven
        // pulleys are two non-empty node ranges within the mesh, both spinning
        // about +x in the SAME direction (open belt), with the larger driven
        // pulley slower (smaller |rate|). The belt + base are left static.
        let product = beltdrive_product();
        let loaded = product
            .mesh
            .as_ref()
            .expect("belt-drive product has a mesh");
        let node_count = loaded.mesh.nodes.len();
        let anim = product.animation.expect("belt-drive product is animated");
        assert!(!anim.playing, "starts paused");
        match anim.motion {
            crate::ProductMotion::RigidParts(parts) => {
                assert_eq!(parts.len(), 2, "two rotating parts: driver + driven");
                for p in &parts {
                    assert!(
                        p.node_range.start < p.node_range.end,
                        "non-empty pulley range"
                    );
                    assert!(
                        p.node_range.end <= node_count,
                        "pulley range within the mesh"
                    );
                    assert_eq!(p.axis, [1.0, 0.0, 0.0], "spins about the shaft axis");
                }
                let (drv, drn) = (&parts[0], &parts[1]);
                assert!(
                    drv.rad_per_s.signum() == drn.rad_per_s.signum(),
                    "open belt ⇒ same direction"
                );
                assert!(drv.rad_per_s.abs() > 0.0, "driver spins");
                assert!(
                    drn.rad_per_s.abs() < drv.rad_per_s.abs(),
                    "the larger driven pulley (250 mm vs 100 mm) turns slower"
                );
                // Belt + base nodes lie beyond the two pulley ranges (static).
                assert!(
                    drn.node_range.end < node_count,
                    "belt + base are not animated"
                );
            }
            crate::ProductMotion::Turntable { .. } => {
                panic!("belt drive must use per-part rigid motion")
            }
        }
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use egui::accesskit::{Node, NodeId, Role};

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_beltdrive_workbench(app, ctx);
        });
    }

    /// As [`draw_workbench`], but with accesskit enabled, returning the emitted
    /// accessibility tree nodes — the same tree a screen reader / AI driver
    /// consumes. `accesskit` is re-exported by egui, so no extra dependency.
    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_beltdrive_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_beltdrive_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_beltdrive_workbench = true;
        run_beltdrive(&mut app.beltdrive);
        draw_workbench(&mut app);
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // Every DragValue is a SpinButton that must be `labelled_by` its caption
        // (egui clears a DragValue's own Name), so an AI / screen reader can find
        // the control by the caption text.
        let mut app = ValenxApp::default();
        app.show_beltdrive_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        assert!(
            spin_buttons.len() >= 7,
            "expected the numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every DragValue must be labelled_by its caption (AI-drivable name)"
        );

        for caption in ["driver Ø (mm)", "friction μ", "max tension T1 (N)"] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
    }
}
