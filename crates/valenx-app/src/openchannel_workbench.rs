//! The right-side **Open Channel Workbench** panel — native steady-uniform
//! free-surface hydraulics over `valenx-openchannel`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_openchannel_workbench`,
//! toggled from the View menu. The form sets a prismatic channel
//! (rectangular or trapezoidal), a flow depth, Manning roughness and bed
//! slope; "Analyze" reports the section geometry (area, wetted perimeter,
//! hydraulic radius), the Manning velocity and discharge, the Froude number
//! with the flow regime, and — for the rectangular channel — the
//! hydraulic-jump sequent depth; "Show 3-D" loads a representative open
//! prismatic trough into the central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_openchannel::{
    classify_regime, discharge, froude_for_discharge, sequent_depth, velocity, Channel, FlowRegime,
    GRAVITY_M_S2,
};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Cross-section shape selector. A plain `PartialEq` enum (distinct from
/// the data-carrying [`Channel`]) so it can drive a `radio_value` group;
/// the validated [`Channel`] is built from it at analyze time.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ChannelShape {
    /// Flat bottom, vertical walls.
    Rectangular,
    /// Flat bottom, symmetric sloping walls.
    Trapezoidal,
}

/// Persistent form + result state for the Open Channel Workbench.
pub struct OpenChannelWorkbenchState {
    /// Cross-section shape.
    shape: ChannelShape,
    /// Bottom width `b` (m).
    bottom_width_m: f64,
    /// Side slope `z` = horizontal run per unit rise (trapezoidal only).
    side_slope: f64,
    /// Flow depth `y` (m), measured from the invert to the free surface.
    depth_m: f64,
    /// Manning roughness coefficient `n` (SI form, dimensionless).
    manning_n: f64,
    /// Channel-bed (energy-line) slope `S` (m/m).
    slope: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D trough solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for OpenChannelWorkbenchState {
    fn default() -> Self {
        // A concrete-lined trapezoidal canal: b = 3 m, side slope z = 1.5,
        // flowing 1.2 m deep at n = 0.013, S = 0.001. R ~ 0.78 m, the flow
        // is well subcritical (Fr ~ 0.3).
        Self {
            shape: ChannelShape::Trapezoidal,
            bottom_width_m: 3.0,
            side_slope: 1.5,
            depth_m: 1.2,
            manning_n: 0.013,
            slope: 0.001,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Open Channel Workbench right-side panel. A no-op when the
/// `show_openchannel_workbench` toggle is off.
pub fn draw_openchannel_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_openchannel_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_openchannel_workbench",
        "Open Channel",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native steady-uniform free-surface hydraulics · valenx-openchannel",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.openchannel;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Cross-section").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.shape, ChannelShape::Rectangular, "rectangular");
                        ui.radio_value(&mut s.shape, ChannelShape::Trapezoidal, "trapezoidal");
                    });
                    ui.horizontal(|ui| {
                        ui.label("bottom width b (m)");
                        ui.add(egui::DragValue::new(&mut s.bottom_width_m).speed(0.1));
                    });
                    if s.shape == ChannelShape::Trapezoidal {
                        ui.horizontal(|ui| {
                            ui.label("side slope z (H:V)");
                            ui.add(egui::DragValue::new(&mut s.side_slope).speed(0.1));
                        });
                    }
                    ui.horizontal(|ui| {
                        ui.label("flow depth y (m)");
                        ui.add(egui::DragValue::new(&mut s.depth_m).speed(0.05));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Flow").strong());
                    ui.horizontal(|ui| {
                        ui.label("Manning n");
                        ui.add(egui::DragValue::new(&mut s.manning_n).speed(0.001));
                    });
                    ui.horizontal(|ui| {
                        ui.label("bed slope S (m/m)");
                        ui.add(egui::DragValue::new(&mut s.slope).speed(0.0005));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_openchannel(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative open prismatic trough (channel walls + bed) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Hydraulics").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_openchannel_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.openchannel` borrow is
    // released here): build the trough's 3-D solid and load it.
    if app.openchannel.show_3d_request {
        app.openchannel.show_3d_request = false;
        load_channel_3d(app);
    }
}

/// Validate the form, evaluate the channel and format the readout.
fn run_openchannel(s: &mut OpenChannelWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build the validated [`Channel`] from the form's shape selector and
/// dimensions. Extracted so it is shared by the readout and the 3-D gate.
fn build_channel(s: &OpenChannelWorkbenchState) -> Result<Channel, String> {
    match s.shape {
        ChannelShape::Rectangular => Channel::rectangular(s.bottom_width_m),
        ChannelShape::Trapezoidal => Channel::trapezoidal(s.bottom_width_m, s.side_slope),
    }
    .map_err(|e| e.to_string())
}

/// Evaluate the channel and format the full readout, mapping any domain
/// error to a display string. Extracted so it is unit-testable.
fn compute(s: &OpenChannelWorkbenchState) -> Result<String, String> {
    let channel = build_channel(s)?;
    let area = channel.area_m2(s.depth_m).map_err(|e| e.to_string())?;
    let perimeter = channel
        .wetted_perimeter_m(s.depth_m)
        .map_err(|e| e.to_string())?;
    let radius = channel
        .hydraulic_radius_m(s.depth_m)
        .map_err(|e| e.to_string())?;
    let v = velocity(&channel, s.depth_m, s.manning_n, s.slope).map_err(|e| e.to_string())?;
    let q = discharge(&channel, s.depth_m, s.manning_n, s.slope).map_err(|e| e.to_string())?;

    // Froude / regime need a positive discharge (a zero slope stalls the
    // flow); report the regime only when there is flow to classify.
    let (froude_line, regime_line) = if q > 0.0 {
        let fr = froude_for_discharge(&channel, q, s.depth_m, GRAVITY_M_S2)
            .map_err(|e| e.to_string())?;
        let regime = match classify_regime(fr, 1e-6).map_err(|e| e.to_string())? {
            FlowRegime::Subcritical => "subcritical (Fr < 1)",
            FlowRegime::Critical => "critical (Fr = 1)",
            FlowRegime::Supercritical => "supercritical (Fr > 1)",
        };
        (
            format!("Froude Fr      : {fr:.3}"),
            format!("regime         : {regime}"),
        )
    } else {
        (
            "Froude Fr      : — (no flow at S = 0)".to_string(),
            "regime         : —".to_string(),
        )
    };

    // Sequent depth of a hydraulic jump is the rectangular Bélanger closed
    // form; only meaningful for a supercritical (Fr > 1) rectangular state.
    let jump_line = match s.shape {
        ChannelShape::Rectangular if q > 0.0 => {
            let fr = froude_for_discharge(&channel, q, s.depth_m, GRAVITY_M_S2)
                .map_err(|e| e.to_string())?;
            if fr > 1.0 {
                let y2 = sequent_depth(s.depth_m, fr).map_err(|e| e.to_string())?;
                format!("\nsequent depth  : {y2:.3} m (hydraulic jump)")
            } else {
                "\nsequent depth  : — (subcritical, no jump)".to_string()
            }
        }
        ChannelShape::Rectangular => "\nsequent depth  : — (no flow)".to_string(),
        ChannelShape::Trapezoidal => {
            "\nsequent depth  : — (rectangular closed form only)".to_string()
        }
    };

    let shape_name = match s.shape {
        ChannelShape::Rectangular => "rectangular",
        ChannelShape::Trapezoidal => "trapezoidal",
    };

    Ok(format!(
        "shape          : {shape_name}\n\
         bottom width b : {:.3} m\n\
         flow depth y   : {:.3} m\n\
         Manning n / S  : {:.3} / {:.4}\n\n\
         flow area A    : {area:.3} m²\n\
         wetted perim P : {perimeter:.3} m\n\
         hydraulic R    : {radius:.3} m\n\
         velocity v     : {v:.3} m/s\n\
         discharge Q    : {q:.3} m³/s\n\
         {froude_line}\n\
         {regime_line}{jump_line}",
        s.bottom_width_m, s.depth_m, s.manning_n, s.slope,
    ))
}

/// Build the open channel as a triangle [`Mesh`] — a prismatic trough run
/// along `y`: a flat bed plus the two side walls (vertical for the
/// rectangular section, battered by the side slope `z` for the
/// trapezoidal one). The free surface is left open (no top face), which is
/// the defining feature of an open channel. Representative geometry; the
/// hydraulic numbers are the `valenx-openchannel` result. `None` for an
/// invalid configuration.
fn channel_trough_mesh(s: &OpenChannelWorkbenchState) -> Option<Mesh> {
    let channel = build_channel(s).ok()?;
    // Validate the flow depth through the geometry (rejects non-positive /
    // non-finite); on success draw the trough at that depth.
    channel.area_m2(s.depth_m).ok()?;
    let b = channel.bottom_width_m();
    let z = channel.side_slope();
    let y = s.depth_m;
    let length = 6.0_f64;

    // Cross-section corners in the x (transverse) / z (vertical) plane,
    // extruded along y. Half the bottom width, and the top half-width
    // widened by the side slope.
    let half_b = 0.5 * b;
    let half_top = half_b + z * y;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Eight corners: bed (z = 0) and rim (z = y), at the two ends y = 0
    // and y = length. Index layout per end: 0 bed-left, 1 bed-right,
    // 2 rim-right, 3 rim-left.
    for &yy in &[0.0_f64, length] {
        nodes.push(Vector3::new(-half_b, yy, 0.0));
        nodes.push(Vector3::new(half_b, yy, 0.0));
        nodes.push(Vector3::new(half_top, yy, y));
        nodes.push(Vector3::new(-half_top, yy, y));
    }

    // Wetted faces (bed + two walls) as the swept trough, plus the two end
    // caps so the solid is closed apart from the open free surface.
    let quads = [
        [0, 1, 5, 4], // bed
        [1, 2, 6, 5], // right wall
        [3, 0, 4, 7], // left wall
    ];
    for q in quads {
        tris.extend_from_slice(&[q[0], q[1], q[2], q[0], q[2], q[3]]);
    }
    // End caps (U-shaped cross-section, open at the top): two triangles per
    // end spanning bed-left, bed-right, rim-right, rim-left.
    let caps = [[0, 1, 2, 3], [4, 7, 6, 5]];
    for c in caps {
        tris.extend_from_slice(&[c[0], c[1], c[2], c[0], c[2], c[3]]);
    }

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-openchannel");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D trough solid and load it into the central viewport.
fn load_channel_3d(app: &mut ValenxApp) {
    let Some(mesh) = channel_trough_mesh(&app.openchannel) else {
        app.openchannel.error =
            Some("channel parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<channel>/valenx-openchannel"),
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
        let s = OpenChannelWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_geometry_velocity_and_regime() {
        let mut s = OpenChannelWorkbenchState::default();
        run_openchannel(&mut s);
        assert!(
            s.error.is_none(),
            "default channel should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("hydraulic R"));
        assert!(s.result.contains("discharge Q"));
        assert!(s.result.contains("Froude Fr"));
        // The default canal flows deep and slow: subcritical.
        assert!(s.result.contains("subcritical"));
    }

    #[test]
    fn analyze_rejects_zero_width() {
        let mut s = OpenChannelWorkbenchState {
            bottom_width_m: 0.0,
            ..Default::default()
        };
        run_openchannel(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn rectangular_hydraulic_radius_is_area_over_perimeter() {
        // GROUND TRUTH: for a rectangular channel R = A/P = b*y/(b + 2y).
        // b = 3, y = 1  ->  A = 3, P = 5, R = 3/5 = 0.6 m, and Manning
        // gives v = (1/n) R^(2/3) S^(1/2). Pin the literals before .sqrt().
        let s = OpenChannelWorkbenchState {
            shape: ChannelShape::Rectangular,
            bottom_width_m: 3.0,
            side_slope: 0.0,
            depth_m: 1.0,
            manning_n: 0.013,
            slope: 0.001,
            ..Default::default()
        };
        let channel = build_channel(&s).unwrap();
        let r = channel.hydraulic_radius_m(1.0).unwrap();
        let r_hand = (3.0 * 1.0) / (3.0 + 2.0 * 1.0);
        assert!((r - r_hand).abs() < 1e-12);
        assert!((r - 0.6).abs() < 1e-12);
        let one = 1.0_f64;
        let v = velocity(&channel, 1.0, 0.013, 0.001).unwrap();
        let v_hand = (one / 0.013) * r_hand.powf(2.0 / 3.0) * 0.001_f64.sqrt();
        assert!((v - v_hand).abs() < 1e-9);
    }

    #[test]
    fn trough_mesh_for_default_is_nonempty_and_in_range() {
        let s = OpenChannelWorkbenchState::default();
        let mesh = channel_trough_mesh(&s).expect("default channel yields a trough");
        assert!(mesh.nodes.len() >= 8, "expected the eight trough corners");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn trough_mesh_none_for_invalid() {
        let s = OpenChannelWorkbenchState {
            bottom_width_m: 0.0,
            ..Default::default()
        };
        assert!(channel_trough_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_openchannel_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_openchannel_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_openchannel_workbench = true;
        run_openchannel(&mut app.openchannel);
        draw_workbench(&mut app);
    }
}
