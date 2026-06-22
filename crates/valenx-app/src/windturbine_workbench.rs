//! The right-side **Wind Turbine Workbench** panel — native horizontal-axis
//! wind-turbine power performance over `valenx-windturbine`.
//!
//! Mirrors the Rail / Drone / Marine workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_windturbine_workbench`,
//! toggled from the View menu. The form drives the closed-form
//! `valenx-windturbine` actuator-disc model; "Analyze" reports the swept
//! area, the wind / Betz / Cp-extracted power, the idealised power-curve
//! output and region, the capacity factor and the tip-speed ratio, and
//! "Show 3-D turbine" loads a tower + nacelle + three-blade rotor solid into
//! the central viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_windturbine::{
    available_power, betz_power, extracted_power, rpm_to_rad_per_s, swept_area, tip_speed,
    tip_speed_ratio, validate_cp, PowerCurve, Region, AIR_DENSITY_SEA_LEVEL, BETZ_LIMIT,
};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Wind Turbine Workbench.
pub struct WindTurbineWorkbenchState {
    /// Rotor radius `R` (m).
    rotor_radius_m: f64,
    /// Power coefficient `Cp` in `[0, 16/27]`.
    power_coefficient: f64,
    /// Air density `rho` (kg/m^3).
    air_density: f64,
    /// Operating wind speed `v` (m/s) the readout is evaluated at.
    wind_speed_m_s: f64,
    /// Rotor speed (rev/min) for the tip-speed ratio.
    rotor_rpm: f64,
    /// Power-curve cut-in wind speed (m/s).
    cut_in_m_s: f64,
    /// Power-curve rated wind speed (m/s).
    rated_m_s: f64,
    /// Power-curve cut-out wind speed (m/s).
    cut_out_m_s: f64,
    /// Rated electrical power (MW).
    rated_power_mw: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D turbine solid and load it into the
    /// central viewport (serviced after the panel draws).
    show_3d_request: bool,
}

impl Default for WindTurbineWorkbenchState {
    fn default() -> Self {
        // A 3 MW-class onshore turbine: 90 m rotor, Cp 0.45, sea-level air.
        Self {
            rotor_radius_m: 45.0,
            power_coefficient: 0.45,
            air_density: AIR_DENSITY_SEA_LEVEL,
            wind_speed_m_s: 10.0,
            rotor_rpm: 15.0,
            cut_in_m_s: 3.0,
            rated_m_s: 12.0,
            cut_out_m_s: 25.0,
            rated_power_mw: 3.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Wind Turbine Workbench right-side panel. A no-op when the
/// `show_windturbine_workbench` toggle is off.
pub fn draw_windturbine_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_windturbine_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_windturbine_workbench",
        "Wind Turbine",
        |app, ui| {
            ui.label(
                egui::RichText::new("native actuator-disc wind-turbine power · valenx-windturbine")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.windturbine;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Rotor + air").strong());
                    ui.horizontal(|ui| {
                        ui.label("rotor radius R (m)");
                        ui.add(egui::DragValue::new(&mut s.rotor_radius_m).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("power coeff Cp");
                        ui.add(egui::DragValue::new(&mut s.power_coefficient).speed(0.005));
                    });
                    ui.horizontal(|ui| {
                        ui.label("air ρ (kg/m³)");
                        ui.add(egui::DragValue::new(&mut s.air_density).speed(0.01));
                    });
                    ui.horizontal(|ui| {
                        ui.label("rotor speed (rpm)");
                        ui.add(egui::DragValue::new(&mut s.rotor_rpm).speed(0.2));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Power curve").strong());
                    ui.horizontal(|ui| {
                        ui.label("cut-in (m/s)");
                        ui.add(egui::DragValue::new(&mut s.cut_in_m_s).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("rated (m/s)");
                        ui.add(egui::DragValue::new(&mut s.rated_m_s).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("cut-out (m/s)");
                        ui.add(egui::DragValue::new(&mut s.cut_out_m_s).speed(0.2));
                    });
                    ui.horizontal(|ui| {
                        ui.label("rated power (MW)");
                        ui.add(egui::DragValue::new(&mut s.rated_power_mw).speed(0.1));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Operating point").strong());
                    ui.horizontal(|ui| {
                        ui.label("wind speed v (m/s)");
                        ui.add(egui::DragValue::new(&mut s.wind_speed_m_s).speed(0.2));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_turbine(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D turbine").strong())
                        .on_hover_text(
                            "Build the tower + nacelle + three-blade rotor as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Performance").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_windturbine_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.windturbine` borrow is
    // released here): build the turbine's 3-D solid and load it.
    if app.windturbine.show_3d_request {
        app.windturbine.show_3d_request = false;
        load_turbine_3d(app);
    }
}

/// Compute the full performance readout, mapping any domain error to a
/// display string. Extracted so it is unit-testable.
fn compute(s: &WindTurbineWorkbenchState) -> Result<String, String> {
    let area = swept_area(s.rotor_radius_m).map_err(|e| e.to_string())?;
    let v = s.wind_speed_m_s;
    let avail = available_power(s.air_density, area, v).map_err(|e| e.to_string())?;
    let betz = betz_power(s.air_density, area, v).map_err(|e| e.to_string())?;
    let extracted =
        extracted_power(s.air_density, area, v, s.power_coefficient).map_err(|e| e.to_string())?;
    let rated_w = s.rated_power_mw * 1.0e6;
    let curve = PowerCurve::new(s.cut_in_m_s, s.rated_m_s, s.cut_out_m_s, rated_w)
        .map_err(|e| e.to_string())?;
    let p_out = curve.power(v);
    let region = match curve.region(v) {
        Region::Idle => "idle (out of band)",
        Region::Ramp => "ramp (cube-law)",
        Region::Rated => "rated plateau",
    };
    let omega = rpm_to_rad_per_s(s.rotor_rpm).map_err(|e| e.to_string())?;
    let tip = tip_speed(omega, s.rotor_radius_m).map_err(|e| e.to_string())?;
    let tsr = tip_speed_ratio(omega, s.rotor_radius_m, v).map_err(|e| e.to_string())?;
    let cf = if rated_w > 0.0 { p_out / rated_w } else { 0.0 };
    let annual_gwh = p_out * 8760.0 / 1.0e9;
    Ok(format!(
        "rotor radius R : {:.1} m  (diameter {:.1} m)\n\
         swept area A   : {:.0} m²\n\
         power coeff Cp : {:.3}  (Betz max {:.3})\n\
         air ρ          : {:.3} kg/m³\n\n\
         at wind v = {:.1} m/s:\n\
         wind power     : {:.2} MW  (kinetic, in the wind)\n\
         Betz max 16/27 : {:.2} MW\n\
         extracted (Cp) : {:.2} MW\n\
         power curve    : {:.2} MW  [{}]\n\
         capacity factor: {:.2}\n\
         tip speed      : {:.1} m/s\n\
         tip-speed ratio: {:.2}\n\n\
         annual energy  : {:.1} GWh  (steady-wind approx.)",
        s.rotor_radius_m,
        2.0 * s.rotor_radius_m,
        area,
        s.power_coefficient,
        BETZ_LIMIT,
        s.air_density,
        v,
        avail / 1.0e6,
        betz / 1.0e6,
        extracted / 1.0e6,
        p_out / 1.0e6,
        region,
        cf,
        tip,
        tsr,
        annual_gwh,
    ))
}

/// Validate the form, compute the performance and format the readout.
fn run_turbine(s: &mut WindTurbineWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
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

/// Append a flat (double-sided) quad `a-b-c-d` (used for the thin blades).
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

/// Append a vertical (z-axis) cylinder side wall, double-sided, of radius `r`
/// rising `height` from the base centre `base`, with `seg` segments. Used for
/// the tower.
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

/// Build the turbine as a triangle [`Mesh`] — a tapered cylindrical tower, a
/// nacelle box at hub height and a three-blade rotor (thin tapered blades in
/// the rotor plane, perpendicular to the wind along +x). The blade length
/// follows the rotor radius. `None` for an invalid configuration.
/// Presentation spin rate of the rotor (hub + blades), rad/s (~0.2 rev/s) — a
/// slow, readable turbine spin keyed to the rotor's large size.
const ROTOR_RAD_PER_S: f32 = 1.2;

/// Build the turbine as a triangle [`Mesh`] together with the
/// [`crate::RigidPart`] for the rotating rotor (the hub + three blades), so the
/// tower and nacelle stay put while the rotor spins about the wind (+x) axis
/// through the hub. `None` for an invalid configuration.
fn turbine_solid_mesh_parts(
    s: &WindTurbineWorkbenchState,
) -> Option<(Mesh, Vec<crate::RigidPart>)> {
    // Gate the 3-D build on a physically valid configuration.
    let r = s.rotor_radius_m;
    if !(r.is_finite() && r > 0.0) {
        return None;
    }
    validate_cp(s.power_coefficient).ok()?;
    PowerCurve::new(
        s.cut_in_m_s,
        s.rated_m_s,
        s.cut_out_m_s,
        s.rated_power_mw * 1.0e6,
    )
    .ok()?;

    let hub_height = r * 1.6;
    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Tower (vertical cylinder from the ground to the nacelle).
    push_cyl_z(
        &mut nodes,
        &mut tris,
        Vector3::zeros(),
        hub_height,
        r * 0.035,
        16,
    );
    // Nacelle (machine housing atop the tower).
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, hub_height),
        Vector3::new(r * 0.12, r * 0.05, r * 0.05),
    );
    // Rotating rotor (hub + three blades) is built last; record its node span so
    // the animation rotates exactly the rotor and leaves the tower/nacelle still.
    let rotor_start = nodes.len();
    // Hub (at the upwind front of the nacelle, +x).
    let hub = Vector3::new(r * 0.15, 0.0, hub_height);
    push_box(
        &mut nodes,
        &mut tris,
        hub,
        Vector3::new(r * 0.04, r * 0.04, r * 0.04),
    );
    // Three blades, 120° apart, in the rotor plane (the y-z plane at the hub).
    for k in 0..3 {
        let th = k as f64 / 3.0 * TAU;
        let (c, sn) = (th.cos(), th.sin());
        let radial = Vector3::new(0.0, c, sn); // outward in the rotor plane
        let chord = Vector3::new(0.0, -sn, c); // perpendicular (chord) in-plane
        let root = hub + radial * (r * 0.05);
        let tip = hub + radial * r;
        let root_hc = r * 0.04;
        let tip_hc = r * 0.012;
        push_quad(
            &mut nodes,
            &mut tris,
            root + chord * root_hc,
            root - chord * root_hc,
            tip - chord * tip_hc,
            tip + chord * tip_hc,
        );
    }
    let rotor_end = nodes.len();

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-windturbine");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    // The rotor spins about the wind axis (+x) through the hub centre.
    let parts = vec![crate::RigidPart {
        node_range: rotor_start..rotor_end,
        axis: [1.0, 0.0, 0.0],
        pivot: [hub.x as f32, hub.y as f32, hub.z as f32],
        rad_per_s: ROTOR_RAD_PER_S,
    }];
    Some((mesh, parts))
}

/// Build the turbine as a triangle [`Mesh`] (without the rotor part metadata)
/// for the central viewport. See [`turbine_solid_mesh_parts`].
fn turbine_solid_mesh(s: &WindTurbineWorkbenchState) -> Option<Mesh> {
    turbine_solid_mesh_parts(s).map(|(mesh, _parts)| mesh)
}

/// Build the 3-D turbine solid and load it into the central viewport.
fn load_turbine_3d(app: &mut ValenxApp) {
    let Some(mesh) = turbine_solid_mesh(&app.windturbine) else {
        app.windturbine.error =
            Some("turbine parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<turbine>/valenx-windturbine"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"windturbine"}`** product: the canonical
/// tower + nacelle + three-blade rotor solid (the panel's "Show 3-D turbine"
/// geometry) paired with the workbench's own actuator-disc power headline
/// numbers, at a fixed 3/4 camera. Registered in
/// [`crate::products_registry`]; the per-tool builder the registry dispatches
/// to. Pure — driven off [`WindTurbineWorkbenchState::default`].
///
/// The readout rows mirror the panel's `compute()` performance readout.
pub(crate) fn windturbine_product() -> crate::WorkspaceProduct {
    let s = WindTurbineWorkbenchState::default();
    let (mesh, parts) = turbine_solid_mesh_parts(&s).expect("default 3 MW turbine ⇒ a 3-D solid");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<turbine>/valenx-windturbine");
    let readout = compute(&s).expect("default 3 MW turbine ⇒ a valid readout");
    let lines = crate::products_registry::lines_from_readout(&readout);
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Wind Turbine".into(),
        lines,
        mesh: Some(loaded),
        vertex_colors: None,
        camera,
        kind2d: None,
        last_export: None,
        image: None,
        image_texture: None,
        // Animated: the rotor (hub + three blades) spins about the wind axis
        // through the hub while the tower and nacelle stay put. Paused at t = 0.
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
        let s = WindTurbineWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_power_curve_and_tsr() {
        let mut s = WindTurbineWorkbenchState::default();
        run_turbine(&mut s);
        assert!(s.error.is_none());
        assert!(s.result.contains("power curve"));
        assert!(s.result.contains("tip-speed ratio"));
        assert!(s.result.contains("Betz"));
    }

    #[test]
    fn analyze_rejects_super_betz_cp() {
        let mut s = WindTurbineWorkbenchState {
            power_coefficient: 0.7, // above the Betz limit 16/27
            ..Default::default()
        };
        run_turbine(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn analyze_rejects_inconsistent_power_curve() {
        let mut s = WindTurbineWorkbenchState {
            rated_m_s: 2.0, // rated < cut_in -> inconsistent
            ..Default::default()
        };
        run_turbine(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn turbine_mesh_for_default_is_nonempty_and_in_range() {
        let s = WindTurbineWorkbenchState::default();
        let mesh = turbine_solid_mesh(&s).expect("default turbine yields a solid");
        assert!(mesh.nodes.len() > 8, "expected tower + nacelle + blades");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn turbine_mesh_none_for_invalid() {
        let s = WindTurbineWorkbenchState {
            rotor_radius_m: 0.0,
            ..Default::default()
        };
        assert!(turbine_solid_mesh(&s).is_none());
    }

    #[test]
    fn windturbine_product_spins_the_rotor_only() {
        // The product carries a RigidParts animation: the rotor (hub + blades, a
        // non-empty trailing node range) spins about +x at a non-zero rate; the
        // leading tower + nacelle nodes are left static.
        let product = windturbine_product();
        let loaded = product.mesh.as_ref().expect("turbine product has a mesh");
        let node_count = loaded.mesh.nodes.len();
        let anim = product.animation.expect("turbine product is animated");
        assert!(!anim.playing, "starts paused");
        match anim.motion {
            crate::ProductMotion::RigidParts(parts) => {
                assert_eq!(parts.len(), 1, "one rotating part: the rotor");
                let p = &parts[0];
                assert!(
                    p.node_range.start < p.node_range.end,
                    "non-empty rotor range"
                );
                assert!(
                    p.node_range.start > 0,
                    "the tower + nacelle precede the rotor (static)"
                );
                assert_eq!(p.node_range.end, node_count, "rotor reaches the final node");
                assert_eq!(p.axis, [1.0, 0.0, 0.0], "spins about the wind axis");
                assert!(p.rad_per_s.abs() > 0.0, "non-zero spin rate");
            }
            crate::ProductMotion::Turntable { .. } => {
                panic!("turbine must use per-part rigid motion")
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
            draw_windturbine_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_windturbine_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_windturbine_workbench = true;
        run_turbine(&mut app.windturbine);
        draw_workbench(&mut app);
    }
}
