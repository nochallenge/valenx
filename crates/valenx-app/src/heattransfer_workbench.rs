//! The right-side **Heat Transfer Workbench** panel — native composite-wall
//! 1-D heat-loss analysis over `valenx-heat-transfer`.
//!
//! Mirrors the Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_heattransfer_workbench`,
//! toggled from the View menu. The form sets a plane wall (area, thickness,
//! conductivity) with an inside and outside convective film; "Analyze"
//! sums the series thermal resistances and reports the total resistance,
//! heat-loss rate and the overall U-value, and "Show 3-D wall" loads a
//! representative wall-slab solid into the central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_heat_transfer::conduction::PlaneWall;
use valenx_heat_transfer::convection::ConvectiveSurface;
use valenx_heat_transfer::network::series_resistance;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Heat Transfer Workbench.
pub struct HeatTransferWorkbenchState {
    /// Wall area `A` (m^2).
    area_m2: f64,
    /// Wall thickness `L` (m).
    thickness_m: f64,
    /// Wall thermal conductivity `k` (W/m·K).
    conductivity_w_per_mk: f64,
    /// Inside-film convection coefficient `h_in` (W/m^2·K).
    h_in: f64,
    /// Outside-film convection coefficient `h_out` (W/m^2·K).
    h_out: f64,
    /// Inside (warm) temperature (deg C).
    t_in_c: f64,
    /// Outside (cold) temperature (deg C).
    t_out_c: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D wall solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for HeatTransferWorkbenchState {
    fn default() -> Self {
        // A 10 m^2 insulated building wall (200 mm of k=0.04 insulation)
        // with still inside air (h=8) and windy outside (h=25), 21 -> 0 C:
        // U ~ 0.19 W/m^2K, ~41 W of loss.
        Self {
            area_m2: 10.0,
            thickness_m: 0.2,
            conductivity_w_per_mk: 0.04,
            h_in: 8.0,
            h_out: 25.0,
            t_in_c: 21.0,
            t_out_c: 0.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Heat Transfer Workbench right-side panel. A no-op when the
/// `show_heattransfer_workbench` toggle is off.
pub fn draw_heattransfer_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_heattransfer_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_heattransfer_workbench",
        "Heat Transfer",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native composite-wall 1-D heat loss · valenx-heat-transfer",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.heattransfer;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Wall").strong());
                    ui.horizontal(|ui| {
                        ui.label("area (m²)");
                        ui.add(egui::DragValue::new(&mut s.area_m2).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("thickness (m)");
                        ui.add(egui::DragValue::new(&mut s.thickness_m).speed(0.01));
                    });
                    ui.horizontal(|ui| {
                        ui.label("conductivity k (W/m·K)");
                        ui.add(egui::DragValue::new(&mut s.conductivity_w_per_mk).speed(0.005));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Convective films").strong());
                    ui.horizontal(|ui| {
                        ui.label("inside h (W/m²K)");
                        ui.add(egui::DragValue::new(&mut s.h_in).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("outside h (W/m²K)");
                        ui.add(egui::DragValue::new(&mut s.h_out).speed(0.5));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Temperatures").strong());
                    ui.horizontal(|ui| {
                        ui.label("inside (°C)");
                        ui.add(egui::DragValue::new(&mut s.t_in_c).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("outside (°C)");
                        ui.add(egui::DragValue::new(&mut s.t_out_c).speed(0.5));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_wall(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D wall").strong())
                        .on_hover_text(
                            "Build a representative composite wall slab (with inside / outside surface films) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Heat loss").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_heattransfer_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.heattransfer` borrow is
    // released here): build the wall's 3-D solid and load it.
    if app.heattransfer.show_3d_request {
        app.heattransfer.show_3d_request = false;
        load_wall_3d(app);
    }
}

/// Validate the form, evaluate the wall and format the readout.
fn run_wall(s: &mut HeatTransferWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// The three series thermal resistances `(R_inside_film, R_conduction,
/// R_outside_film)` in K/W, the quantities both the readout and the 3-D
/// gate need. Extracted so it is unit-testable and shared.
fn resistances(s: &HeatTransferWorkbenchState) -> Result<(f64, f64, f64), String> {
    let wall = PlaneWall::new(s.thickness_m, s.area_m2, s.conductivity_w_per_mk)
        .map_err(|e| e.to_string())?;
    let surf_in = ConvectiveSurface::new(s.area_m2, s.h_in).map_err(|e| e.to_string())?;
    let surf_out = ConvectiveSurface::new(s.area_m2, s.h_out).map_err(|e| e.to_string())?;
    Ok((
        surf_in.resistance(),
        wall.resistance(),
        surf_out.resistance(),
    ))
}

/// Evaluate the wall and format the full readout, mapping any domain error
/// to a display string. Extracted so it is unit-testable.
fn compute(s: &HeatTransferWorkbenchState) -> Result<String, String> {
    let (r_in, r_cond, r_out) = resistances(s)?;
    let r_total = series_resistance(&[r_in, r_cond, r_out]).map_err(|e| e.to_string())?;
    let dt = s.t_in_c - s.t_out_c;
    let q = dt / r_total;
    let u = 1.0 / (r_total * s.area_m2);

    Ok(format!(
        "wall area       : {:.2} m²\n\
         thickness / k   : {:.3} m / {:.3} W/m·K\n\
         films h in/out  : {:.1} / {:.1} W/m²K\n\
         inside / outside: {:.1} / {:.1} °C\n\n\
         R inside film   : {:.4} K/W\n\
         R conduction    : {:.4} K/W\n\
         R outside film  : {:.4} K/W\n\
         R total         : {:.4} K/W\n\
         heat loss Q     : {:.2} W\n\
         U-value         : {:.3} W/m²K",
        s.area_m2,
        s.thickness_m,
        s.conductivity_w_per_mk,
        s.h_in,
        s.h_out,
        s.t_in_c,
        s.t_out_c,
        r_in,
        r_cond,
        r_out,
        r_total,
        q,
        u,
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

/// Build the composite wall as a triangle [`Mesh`] — a wall slab (thin in
/// the through-thickness `x` direction) with thin inside / outside surface
/// films and a base. Representative geometry (not to scale; the heat-loss
/// numbers are the `valenx-heat-transfer` result). `None` for an invalid
/// configuration.
fn wall_solid_mesh(s: &HeatTransferWorkbenchState) -> Option<Mesh> {
    resistances(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Main wall slab.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.6),
        Vector3::new(0.07, 0.7, 0.5),
    );
    // Inside surface film (-x face).
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(-0.09, 0.0, 0.6),
        Vector3::new(0.012, 0.66, 0.46),
    );
    // Outside surface film (+x face).
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.09, 0.0, 0.6),
        Vector3::new(0.012, 0.66, 0.46),
    );
    // Base.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.05),
        Vector3::new(0.25, 0.7, 0.05),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-heat-transfer");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D wall solid and load it into the central viewport.
fn load_wall_3d(app: &mut ValenxApp) {
    let Some(mesh) = wall_solid_mesh(&app.heattransfer) else {
        app.heattransfer.error =
            Some("wall parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<wall>/valenx-heat-transfer"),
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
        let s = HeatTransferWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_heat_loss_and_uvalue() {
        let mut s = HeatTransferWorkbenchState::default();
        run_wall(&mut s);
        assert!(
            s.error.is_none(),
            "default wall should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("R total"));
        assert!(s.result.contains("heat loss Q"));
        assert!(s.result.contains("U-value"));
        // 200 mm of k=0.04 dominates: U ~ 0.19 W/m^2K.
        assert!(s.result.contains("0.19"));
    }

    #[test]
    fn analyze_rejects_zero_area() {
        let mut s = HeatTransferWorkbenchState {
            area_m2: 0.0,
            ..Default::default()
        };
        run_wall(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn series_resistances_sum_and_q_times_r_equals_delta_t() {
        // Ground truth: series thermal resistances add, and Q * R_total
        // equals the temperature drop across the assembly.
        let area = 10.0;
        let r_in = 1.0 / (8.0 * area);
        let r_cond = 0.2 / (0.04 * area);
        let r_out = 1.0 / (25.0 * area);
        let r_total = series_resistance(&[r_in, r_cond, r_out]).unwrap();
        assert!((r_total - (r_in + r_cond + r_out)).abs() < 1e-12);
        let dt = 21.0;
        let q = dt / r_total;
        assert!((q * r_total - dt).abs() < 1e-9);
    }

    #[test]
    fn wall_mesh_for_default_is_nonempty_and_in_range() {
        let s = HeatTransferWorkbenchState::default();
        let mesh = wall_solid_mesh(&s).expect("default wall yields a solid");
        assert!(mesh.nodes.len() > 8, "expected slab + films + base");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn wall_mesh_none_for_invalid() {
        let s = HeatTransferWorkbenchState {
            area_m2: 0.0,
            ..Default::default()
        };
        assert!(wall_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_heattransfer_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_heattransfer_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_heattransfer_workbench = true;
        run_wall(&mut app.heattransfer);
        draw_workbench(&mut app);
    }
}
