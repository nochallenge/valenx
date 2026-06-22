//! The right-side **Heat Pump Workbench** panel — native Carnot COP +
//! second-law (Carnot-fraction) derating over `valenx-heatpump`.
//!
//! Mirrors the Battery Pack / Pipe Flow workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_heatpump_workbench`,
//! toggled from the View menu. The form sets the cold and hot reservoir
//! temperatures and a Carnot fraction (second-law efficiency); "Analyze"
//! reports the temperature lift, the ideal heating / cooling COP and the
//! real-machine derated COPs, and "Show 3-D unit" loads a representative
//! outdoor-unit solid into the central viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_heatpump::cop::CarnotCop;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Absolute-zero offset for the Celsius ↔ Kelvin conversion.
const KELVIN_OFFSET: f64 = 273.15;

/// Persistent form + result state for the Heat Pump Workbench.
pub struct HeatPumpWorkbenchState {
    /// Cold-reservoir temperature `T_c` (deg C) — the source when heating.
    t_cold_c: f64,
    /// Hot-reservoir temperature `T_h` (deg C) — the sink when heating.
    t_hot_c: f64,
    /// Carnot fraction (second-law efficiency) in (0, 1].
    carnot_fraction: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D unit solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for HeatPumpWorkbenchState {
    fn default() -> Self {
        // A typical air-source heat pump heating from 0 deg C outdoor air
        // to a 35 deg C low-temperature supply at ~50% of the Carnot
        // limit: lift 35 K, ideal COP_heat ~8.8, real COP_heat ~4.4.
        Self {
            t_cold_c: 0.0,
            t_hot_c: 35.0,
            carnot_fraction: 0.5,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Heat Pump Workbench right-side panel. A no-op when the
/// `show_heatpump_workbench` toggle is off.
pub fn draw_heatpump_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_heatpump_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_heatpump_workbench",
        "Heat Pump",
        |app, ui| {
            ui.label(
                egui::RichText::new("native Carnot COP + second-law derating · valenx-heatpump")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.heatpump;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Reservoir temperatures").strong());
                    ui.horizontal(|ui| {
                        ui.label("cold T_c (°C)");
                        ui.add(egui::DragValue::new(&mut s.t_cold_c).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("hot  T_h (°C)");
                        ui.add(egui::DragValue::new(&mut s.t_hot_c).speed(0.5));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Real-machine efficiency").strong());
                    ui.horizontal(|ui| {
                        ui.label("Carnot fraction");
                        ui.add(egui::DragValue::new(&mut s.carnot_fraction).speed(0.01));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_heatpump(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D unit").strong())
                        .on_hover_text(
                            "Build a representative outdoor heat-pump unit (enclosure, fan shroud, compressor) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Coefficient of performance").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_heatpump_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.heatpump` borrow is
    // released here): build the unit's 3-D solid and load it.
    if app.heatpump.show_3d_request {
        app.heatpump.show_3d_request = false;
        load_unit_3d(app);
    }
}

/// Validate the form, evaluate the COPs and format the readout.
fn run_heatpump(s: &mut HeatPumpWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build a validated [`CarnotCop`] from the form, mapping the domain error
/// to a display string.
fn build_cop(s: &HeatPumpWorkbenchState) -> Result<CarnotCop, String> {
    CarnotCop::new(s.t_cold_c + KELVIN_OFFSET, s.t_hot_c + KELVIN_OFFSET).map_err(|e| e.to_string())
}

/// Evaluate the COPs and format the full readout, mapping any domain error
/// to a display string. Extracted so it is unit-testable.
fn compute(s: &HeatPumpWorkbenchState) -> Result<String, String> {
    let cop = build_cop(s)?;
    let derated = cop.derated(s.carnot_fraction).map_err(|e| e.to_string())?;
    Ok(format!(
        "reservoirs      : {:.1} °C  ↔  {:.1} °C\n\
         temperature lift: {:.1} K\n\n\
         ideal COP heat  : {:.2}\n\
         ideal COP cool  : {:.2}   (heat = cool + 1)\n\n\
         Carnot fraction : {:.2}\n\
         real COP heat   : {:.2}\n\
         real COP cool   : {:.2}",
        s.t_cold_c,
        s.t_hot_c,
        cop.lift_k(),
        cop.cop_heat,
        cop.cop_cool,
        derated.fraction,
        derated.cop_heat,
        derated.cop_cool,
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

/// Append a (double-sided) cylinder whose axis runs along `+x`, spanning
/// `base.x ..= base.x + length` with circle centre `(base.y, base.z)`.
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
}

/// Append a (double-sided) cylinder whose axis runs along `+z`, spanning
/// `base.z ..= base.z + height` with circle centre `(base.x, base.y)`.
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

/// Build the heat-pump outdoor unit as a triangle [`Mesh`] — an enclosure
/// box with a circular fan shroud on the front face, an internal
/// compressor cylinder and a base. Representative geometry (the COP
/// physics is the `valenx-heatpump` result). `None` for an inconsistent
/// configuration (non-positive lift).
fn unit_solid_mesh(s: &HeatPumpWorkbenchState) -> Option<Mesh> {
    build_cop(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Enclosure cabinet.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.55),
        Vector3::new(0.6, 0.35, 0.45),
    );
    // Fan shroud — a short fat cylinder on the +x (front) face.
    push_cyl_x(
        &mut nodes,
        &mut tris,
        Vector3::new(0.6, 0.0, 0.6),
        0.08,
        0.3,
        28,
    );
    // Compressor — an upright cylinder beside the fan inside the cabinet.
    push_cyl_z(
        &mut nodes,
        &mut tris,
        Vector3::new(-0.25, 0.0, 0.18),
        0.5,
        0.16,
        20,
    );
    // Base plate / feet.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.05),
        Vector3::new(0.62, 0.37, 0.05),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-heatpump");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D unit solid and load it into the central viewport.
fn load_unit_3d(app: &mut ValenxApp) {
    let Some(mesh) = unit_solid_mesh(&app.heatpump) else {
        app.heatpump.error =
            Some("hot reservoir must exceed the cold one — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<unit>/valenx-heatpump"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// Agent-bridge product: the canonical heat-pump workbench as a 3-D solid plus
/// its `compute()` readout rows (see [`crate::products_registry`]).
pub(crate) fn heatpump_product() -> crate::WorkspaceProduct {
    let s = HeatPumpWorkbenchState::default();
    let mesh = unit_solid_mesh(&s).expect("canonical heat pump ⇒ unit solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<heatpump>/valenx-unit");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical heat pump ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Heat pump (COP/capacity)".into(),
        lines,
        mesh: Some(loaded),
        vertex_colors: None,
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
        let s = HeatPumpWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_cop_and_lift() {
        let mut s = HeatPumpWorkbenchState::default();
        run_heatpump(&mut s);
        assert!(
            s.error.is_none(),
            "default heat pump should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("temperature lift"));
        assert!(s.result.contains("ideal COP heat"));
        assert!(s.result.contains("real COP heat"));
        // 308.15 / 35 = 8.80 ideal heating COP.
        assert!(s.result.contains("8.80"));
    }

    #[test]
    fn analyze_rejects_nonpositive_lift() {
        let mut s = HeatPumpWorkbenchState {
            t_cold_c: 40.0,
            t_hot_c: 35.0,
            ..Default::default()
        };
        run_heatpump(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn ideal_cop_heat_exceeds_cool_by_exactly_one() {
        // Ground truth: COP_heat - COP_cool = T_h/(T_h-T_c) - T_c/(T_h-T_c)
        // = (T_h - T_c)/(T_h - T_c) = 1, exactly, for any valid lift.
        let cop = CarnotCop::new(273.15, 308.15).unwrap();
        assert!((cop.cop_heat - cop.cop_cool - 1.0).abs() < 1e-9);
    }

    #[test]
    fn unit_mesh_for_default_is_nonempty_and_in_range() {
        let s = HeatPumpWorkbenchState::default();
        let mesh = unit_solid_mesh(&s).expect("default unit yields a solid");
        assert!(mesh.nodes.len() > 8, "expected cabinet + fan + compressor");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn unit_mesh_none_for_invalid() {
        let s = HeatPumpWorkbenchState {
            t_cold_c: 40.0,
            t_hot_c: 35.0,
            ..Default::default()
        };
        assert!(unit_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_heatpump_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_heatpump_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_heatpump_workbench = true;
        run_heatpump(&mut app.heatpump);
        draw_workbench(&mut app);
    }
}
