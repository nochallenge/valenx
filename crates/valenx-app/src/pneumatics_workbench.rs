//! The right-side **Pneumatics Workbench** panel — native cylinder force,
//! compression-ratio, choked-flow and free-air-consumption analysis over
//! `valenx-pneumatics`.
//!
//! Mirrors the Heat Transfer / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_pneumatics_workbench`,
//! toggled from the View menu. The form sets a cylinder (single- or
//! double-acting bore / rod), a gauge supply pressure over an atmosphere, a
//! stroke length and a cycle count; "Analyze" reports the theoretical
//! extend / retract thrust (`F = p·A`), the compression ratio, whether the
//! exhaust to atmosphere is choked (the `~0.528` critical ratio for air),
//! the swept volume per cycle and the free-air consumption; "Show 3-D"
//! loads a representative actuator solid (barrel + protruding rod) into the
//! central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use valenx_pneumatics::compression::{
    absolute_pressure, compression_ratio, STANDARD_ATMOSPHERE_PA,
};
use valenx_pneumatics::consumption::{free_air_consumption, swept_volume_per_cycle, Action};
use valenx_pneumatics::cylinder::Cylinder;
use valenx_pneumatics::flow::{is_choked_air, pressure_ratio, CRITICAL_PRESSURE_RATIO_AIR};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Tau (a full turn in radians), for the swept-tube cylinder mesh.
const TAU: f64 = std::f64::consts::TAU;

/// Persistent form + result state for the Pneumatics Workbench.
pub struct PneumaticsWorkbenchState {
    /// Whether the cylinder is single- or double-acting (governs the rod
    /// term and whether a cycle counts one or both strokes).
    action: Action,
    /// Bore (piston) diameter (m).
    bore_m: f64,
    /// Rod diameter (m) — used only on a double-acting cylinder.
    rod_m: f64,
    /// Gauge supply pressure above atmospheric (Pa).
    gauge_pressure_pa: f64,
    /// Local atmospheric pressure (Pa).
    atmospheric_pa: f64,
    /// Stroke length (m).
    stroke_m: f64,
    /// Number of full cycles for the consumption figure.
    cycles: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D actuator solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for PneumaticsWorkbenchState {
    fn default() -> Self {
        // A 50 mm bore / 20 mm rod double-acting cylinder at 6 bar gauge
        // over a standard atmosphere, 100 mm stroke, 100 cycles:
        //   extend  F = 6e5 * pi/4 * 0.05^2          = 1178.10 N
        //   retract F = 6e5 * pi/4 * (0.05^2-0.02^2) =  989.60 N
        //   compression ratio r = 1 + 6e5/101325    = 6.922
        //   exhaust to atmosphere: ratio 0.1445 < 0.528 -> choked
        //   free air / 100 cycles                    = 250.064 L
        Self {
            action: Action::DoubleActing,
            bore_m: 0.05,
            rod_m: 0.02,
            gauge_pressure_pa: 600_000.0,
            atmospheric_pa: STANDARD_ATMOSPHERE_PA,
            stroke_m: 0.1,
            cycles: 100.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Pneumatics Workbench right-side panel. A no-op when the
/// `show_pneumatics_workbench` toggle is off.
pub fn draw_pneumatics_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_pneumatics_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_pneumatics_workbench",
        "Pneumatics",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native cylinder force / consumption / choked-flow · valenx-pneumatics",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.pneumatics;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Cylinder").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.action, Action::SingleActing, "single-acting");
                        ui.radio_value(&mut s.action, Action::DoubleActing, "double-acting");
                    });
                    ui.horizontal(|ui| {
                        ui.label("bore Ø (m)");
                        ui.add(egui::DragValue::new(&mut s.bore_m).speed(0.005));
                    });
                    ui.horizontal(|ui| {
                        ui.label("rod Ø (m)");
                        ui.add(egui::DragValue::new(&mut s.rod_m).speed(0.002));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Supply").strong());
                    ui.horizontal(|ui| {
                        ui.label("gauge pressure (Pa)");
                        ui.add(egui::DragValue::new(&mut s.gauge_pressure_pa).speed(5000.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("atmosphere (Pa)");
                        ui.add(egui::DragValue::new(&mut s.atmospheric_pa).speed(500.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Duty").strong());
                    ui.horizontal(|ui| {
                        ui.label("stroke (m)");
                        ui.add(egui::DragValue::new(&mut s.stroke_m).speed(0.01));
                    });
                    ui.horizontal(|ui| {
                        ui.label("cycles");
                        ui.add(egui::DragValue::new(&mut s.cycles).speed(1.0));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_pneumatics(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative pneumatic actuator (cylinder barrel with a protruding rod) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Analysis").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_pneumatics_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.pneumatics` borrow is
    // released here): build the actuator's 3-D solid and load it.
    if app.pneumatics.show_3d_request {
        app.pneumatics.show_3d_request = false;
        load_cylinder_3d(app);
    }
}

/// Validate the form, evaluate the cylinder and format the readout.
fn run_pneumatics(s: &mut PneumaticsWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build the [`Cylinder`] for the current form (single- or double-acting),
/// mapping any domain error to a display string. Shared by the readout and
/// the 3-D gate.
fn build_cylinder(s: &PneumaticsWorkbenchState) -> Result<Cylinder, String> {
    match s.action {
        Action::SingleActing => Cylinder::single_acting(s.bore_m).map_err(|e| e.to_string()),
        Action::DoubleActing => {
            Cylinder::double_acting(s.bore_m, s.rod_m).map_err(|e| e.to_string())
        }
    }
}

/// Evaluate the cylinder and format the full readout, mapping any domain
/// error to a display string. Extracted so it is unit-testable.
fn compute(s: &PneumaticsWorkbenchState) -> Result<String, String> {
    let cyl = build_cylinder(s)?;

    let extend_f = cyl
        .extend_force(s.gauge_pressure_pa)
        .map_err(|e| e.to_string())?;
    let retract_f = cyl
        .retract_force(s.gauge_pressure_pa)
        .map_err(|e| e.to_string())?;

    let ratio =
        compression_ratio(s.gauge_pressure_pa, s.atmospheric_pa).map_err(|e| e.to_string())?;

    // Exhaust check: the supply (absolute) venting down to the atmosphere.
    let p_up =
        absolute_pressure(s.gauge_pressure_pa, s.atmospheric_pa).map_err(|e| e.to_string())?;
    let vent_ratio = pressure_ratio(p_up, s.atmospheric_pa).map_err(|e| e.to_string())?;
    let choked = is_choked_air(p_up, s.atmospheric_pa).map_err(|e| e.to_string())?;
    let choke_state = if choked { "choked (sonic)" } else { "subsonic" };

    let swept_cycle =
        swept_volume_per_cycle(&cyl, s.stroke_m, s.action).map_err(|e| e.to_string())?;
    let free_air = free_air_consumption(
        &cyl,
        s.stroke_m,
        s.action,
        s.cycles,
        s.gauge_pressure_pa,
        s.atmospheric_pa,
    )
    .map_err(|e| e.to_string())?;

    let acting = match s.action {
        Action::SingleActing => "single-acting",
        Action::DoubleActing => "double-acting",
    };
    // Bind every value to a local matching its placeholder so the format
    // string can capture them inline (no trailing named-arg list).
    let bore = s.bore_m;
    let rod = s.rod_m;
    let bore_area_cm2 = cyl.bore_area() * 1.0e4;
    let pg = s.gauge_pressure_pa;
    let pa = s.atmospheric_pa;
    let stroke = s.stroke_m;
    let cyc = s.cycles;
    let crit = CRITICAL_PRESSURE_RATIO_AIR;
    let swept_l = swept_cycle * 1000.0;
    let free_l = free_air * 1000.0;

    Ok(format!(
        "cylinder        : {acting}\n\
         bore / rod      : {bore:.3} / {rod:.3} m\n\
         bore area       : {bore_area_cm2:.3} cm²\n\
         gauge / atm     : {pg:.0} / {pa:.0} Pa\n\
         stroke / cycles : {stroke:.3} m / {cyc:.0}\n\n\
         extend force    : {extend_f:.2} N\n\
         retract force   : {retract_f:.2} N\n\
         compression r   : {ratio:.3}\n\
         vent ratio      : {vent_ratio:.4} (crit {crit:.4})\n\
         exhaust flow    : {choke_state}\n\
         swept / cycle   : {swept_l:.4} L\n\
         free air        : {free_l:.3} L"
    ))
}

/// Append a solid (capped) cylinder whose axis runs along `+x`, spanning
/// `base.x ..= base.x + length` with circle centre `(base.y, base.z)`:
/// a swept-tube side wall plus a triangle-fan cap at each end.
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
    // Side wall (double-sided so it reads as a solid from either face).
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
    // End caps via a centre node + triangle fan at each face.
    let c0 = nodes.len();
    nodes.push(Vector3::new(x0, base.y, base.z));
    let c1 = nodes.len();
    nodes.push(Vector3::new(x1, base.y, base.z));
    for j in 0..seg {
        let jn = (j + 1) % seg;
        tris.extend_from_slice(&[c0, lo + jn, lo + j, c1, hi + j, hi + jn]);
    }
}

/// Build the actuator as a triangle [`Mesh`] — a cylinder barrel with a
/// thinner rod protruding from one end. Representative geometry (not to
/// scale; the forces / consumption are the `valenx-pneumatics` result).
/// `None` for an invalid configuration.
fn cylinder_solid_mesh(s: &PneumaticsWorkbenchState) -> Option<Mesh> {
    // Gate the 3-D build on a cylinder that actually constructs.
    let cyl = build_cylinder(s).ok()?;

    // Scale the rod relative to the bore so a fat-rod cylinder still looks
    // right; clamp to keep the rod strictly thinner than the barrel.
    let rod_frac = (cyl.rod() / cyl.bore()).clamp(0.15, 0.6);

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();
    let seg = 32;

    // Barrel: radius 0.4, length 1.2, centred on the x axis.
    push_cyl_x(
        &mut nodes,
        &mut tris,
        Vector3::new(-0.6, 0.0, 0.6),
        1.2,
        0.4,
        seg,
    );
    // Rod: thinner, protruding from the +x end of the barrel.
    push_cyl_x(
        &mut nodes,
        &mut tris,
        Vector3::new(0.6, 0.0, 0.6),
        0.9,
        0.4 * rod_frac,
        seg,
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-pneumatics");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D actuator solid and load it into the central viewport.
fn load_cylinder_3d(app: &mut ValenxApp) {
    let Some(mesh) = cylinder_solid_mesh(&app.pneumatics) else {
        app.pneumatics.error =
            Some("cylinder parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<cylinder>/valenx-pneumatics"),
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
        let s = PneumaticsWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_force_ratio_and_consumption() {
        let mut s = PneumaticsWorkbenchState::default();
        run_pneumatics(&mut s);
        assert!(
            s.error.is_none(),
            "default cylinder should analyze: {:?}",
            s.error
        );
        // 50 mm bore at 6 bar gauge -> extend F = 1178.10 N.
        assert!(s.result.contains("1178.10"), "got {}", s.result);
        assert!(s.result.contains("compression r"));
        // 6 bar gauge over a standard atmosphere -> r = 6.922.
        assert!(s.result.contains("6.922"));
        // Venting to atmosphere is well below the critical ratio.
        assert!(s.result.contains("choked (sonic)"));
        assert!(s.result.contains("free air"));
    }

    #[test]
    fn analyze_rejects_zero_bore() {
        let mut s = PneumaticsWorkbenchState {
            bore_m: 0.0,
            ..Default::default()
        };
        run_pneumatics(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn analyze_rejects_rod_at_or_above_bore() {
        // Double-acting with rod >= bore has no valid annular area.
        let mut s = PneumaticsWorkbenchState {
            bore_m: 0.05,
            rod_m: 0.05,
            ..Default::default()
        };
        run_pneumatics(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn force_is_gauge_pressure_times_bore_area() {
        // Ground truth, computed by hand: a 50 mm bore single-acting
        // cylinder at 6 bar (600 000 Pa) gauge gives an extend thrust of
        //   F = p * pi/4 * d^2 = 6e5 * pi/4 * 0.05^2 = 1178.0972... N.
        let pi: f64 = std::f64::consts::PI;
        let expected = 600_000.0 * 0.25 * pi * 0.05 * 0.05;
        let s = PneumaticsWorkbenchState {
            action: Action::SingleActing,
            bore_m: 0.05,
            gauge_pressure_pa: 600_000.0,
            ..Default::default()
        };
        let cyl = build_cylinder(&s).unwrap();
        let f = cyl.extend_force(s.gauge_pressure_pa).unwrap();
        assert!((f - expected).abs() < 1e-6);
        assert!((f - 1_178.097_245_096_172).abs() < 1e-6);
    }

    #[test]
    fn cylinder_mesh_for_default_is_nonempty_and_in_range() {
        let s = PneumaticsWorkbenchState::default();
        let mesh = cylinder_solid_mesh(&s).expect("default cylinder yields a solid");
        assert!(mesh.nodes.len() > 8, "expected barrel + rod");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn cylinder_mesh_none_for_invalid() {
        let s = PneumaticsWorkbenchState {
            bore_m: 0.0,
            ..Default::default()
        };
        assert!(cylinder_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_pneumatics_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_pneumatics_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_pneumatics_workbench = true;
        run_pneumatics(&mut app.pneumatics);
        draw_workbench(&mut app);
    }
}
