//! The right-side **Battery Pack Workbench** panel — native series /
//! parallel pack sizing over `valenx-batterypack`.
//!
//! Mirrors the Pipe Flow / Heat Exchanger workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_batterypack_workbench`,
//! toggled from the View menu. The form sets one repeated [`Cell`] and the
//! `nSmP` topology plus a discharge C-rate; "Analyze" reports the pack
//! voltage, capacity, energy, cell count, discharge current and runtime,
//! and "Show 3-D pack" loads a cell-array enclosure solid into the central
//! viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_batterypack::{current_from_c_rate, runtime_hours_at_c_rate, Cell, PackConfig};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Largest cell grid drawn in the 3-D enclosure (representative — the
/// real cell count can be far larger).
const MAX_COLS: u32 = 10;
const MAX_ROWS: u32 = 6;

/// Persistent form + result state for the Battery Pack Workbench.
pub struct BatteryPackWorkbenchState {
    /// Cell nominal voltage `V_cell` (V).
    cell_voltage_v: f64,
    /// Cell capacity `Q_cell` (Ah).
    cell_capacity_ah: f64,
    /// Cells in series per string `n` (sets the voltage).
    series: u32,
    /// Parallel strings `m` (sets the capacity).
    parallel: u32,
    /// Discharge rate `C` (multiples of the pack capacity per hour).
    c_rate: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D pack solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for BatteryPackWorkbenchState {
    fn default() -> Self {
        // A 13S4P pack of 3.7 V / 3.0 Ah Li-ion 18650 cells — 48.1 V,
        // 12 Ah, ~577 Wh (a small e-bike pack).
        Self {
            cell_voltage_v: 3.7,
            cell_capacity_ah: 3.0,
            series: 13,
            parallel: 4,
            c_rate: 1.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Battery Pack Workbench right-side panel. A no-op when the
/// `show_batterypack_workbench` toggle is off.
pub fn draw_batterypack_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_batterypack_workbench {
        return;
    }

    egui::SidePanel::right("valenx_batterypack_workbench")
        .resizable(true)
        .default_width(360.0)
        .width_range(300.0..=560.0)
        .show(ctx, |ui| {
            if crate::workbench_ui::header(
                ui,
                "Battery Pack",
                "native series / parallel pack sizing · valenx-batterypack",
            ) {
                app.show_batterypack_workbench = false;
            }

            let s = &mut app.batterypack;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Cell").strong());
                    ui.horizontal(|ui| {
                        ui.label("nominal voltage (V)");
                        ui.add(egui::DragValue::new(&mut s.cell_voltage_v).speed(0.05));
                    });
                    ui.horizontal(|ui| {
                        ui.label("capacity (Ah)");
                        ui.add(egui::DragValue::new(&mut s.cell_capacity_ah).speed(0.1));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Topology  nSmP").strong());
                    ui.horizontal(|ui| {
                        ui.label("series n");
                        ui.add(egui::DragValue::new(&mut s.series).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("parallel m");
                        ui.add(egui::DragValue::new(&mut s.parallel).speed(1.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Discharge").strong());
                    ui.horizontal(|ui| {
                        ui.label("C-rate (1/h)");
                        ui.add(egui::DragValue::new(&mut s.c_rate).speed(0.05));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_pack(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D pack").strong())
                        .on_hover_text(
                            "Build a cell-array enclosure as a 3-D solid (a representative grid of cylindrical cells, capped in size) and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Pack").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        });

    // Serviced after the panel draws (the `&mut app.batterypack` borrow is
    // released here): build the pack's 3-D solid and load it.
    if app.batterypack.show_3d_request {
        app.batterypack.show_3d_request = false;
        load_pack_3d(app);
    }
}

/// Validate the form, size the pack and format the readout.
fn run_pack(s: &mut BatteryPackWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build a validated [`PackConfig`] from the form, mapping the domain
/// error to a display string.
fn build_pack(s: &BatteryPackWorkbenchState) -> Result<PackConfig, String> {
    let cell = Cell::new(s.cell_voltage_v, s.cell_capacity_ah).map_err(|e| e.to_string())?;
    PackConfig::new(cell, s.series, s.parallel).map_err(|e| e.to_string())
}

/// Size the pack and format the full readout, mapping any domain error to
/// a display string. Extracted so it is unit-testable.
fn compute(s: &BatteryPackWorkbenchState) -> Result<String, String> {
    let pack = build_pack(s)?;
    let current_a =
        current_from_c_rate(s.c_rate, pack.pack_capacity_ah()).map_err(|e| e.to_string())?;
    let runtime_h = runtime_hours_at_c_rate(s.c_rate).map_err(|e| e.to_string())?;

    Ok(format!(
        "cell             : {:.2} V / {:.2} Ah  ({:.1} Wh)\n\
         topology         : {}S{}P  ({} cells)\n\n\
         pack voltage     : {:.2} V\n\
         pack capacity    : {:.2} Ah\n\
         pack energy      : {:.1} Wh  ({:.3} kWh)\n\n\
         discharge C-rate : {:.2} C\n\
         discharge current: {:.2} A\n\
         runtime at C     : {:.2} h  ({:.0} min)",
        pack.cell.nominal_voltage_v,
        pack.cell.capacity_ah,
        pack.cell.energy_wh(),
        pack.series,
        pack.parallel,
        pack.total_cells(),
        pack.pack_voltage_v(),
        pack.pack_capacity_ah(),
        pack.pack_energy_wh(),
        pack.pack_energy_wh() / 1000.0,
        s.c_rate,
        current_a,
        runtime_h,
        runtime_h * 60.0,
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

/// Build the pack as a triangle [`Mesh`] — a grid of upright cylindrical
/// cells (`series` columns by `parallel` rows, capped at
/// [`MAX_COLS`]x[`MAX_ROWS`] for the drawing) sitting in an enclosure
/// tray. Representative geometry (the electrical sizing is the
/// `valenx-batterypack` result). `None` for an invalid pack.
fn pack_solid_mesh(s: &BatteryPackWorkbenchState) -> Option<Mesh> {
    let pack = build_pack(s).ok()?;

    let cols = pack.series.clamp(1, MAX_COLS);
    let rows = pack.parallel.clamp(1, MAX_ROWS);
    let pitch = 0.22;
    let r = 0.08;
    let height = 0.65;
    let w = (cols as f64) * pitch;
    let d = (rows as f64) * pitch;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Enclosure tray (a shallow open box under and around the cells).
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.06),
        Vector3::new(0.5 * w + 0.08, 0.5 * d + 0.08, 0.06),
    );
    // Cell grid, centred on the origin.
    let x0 = -0.5 * w + 0.5 * pitch;
    let y0 = -0.5 * d + 0.5 * pitch;
    for i in 0..cols {
        for j in 0..rows {
            let cx = x0 + (i as f64) * pitch;
            let cy = y0 + (j as f64) * pitch;
            push_cyl_z(
                &mut nodes,
                &mut tris,
                Vector3::new(cx, cy, 0.12),
                height,
                r,
                14,
            );
        }
    }

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-batterypack");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D pack solid and load it into the central viewport.
fn load_pack_3d(app: &mut ValenxApp) {
    let Some(mesh) = pack_solid_mesh(&app.batterypack) else {
        app.batterypack.error =
            Some("pack configuration is invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<pack>/valenx-batterypack"),
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
        let s = BatteryPackWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_voltage_and_energy() {
        let mut s = BatteryPackWorkbenchState::default();
        run_pack(&mut s);
        assert!(
            s.error.is_none(),
            "default pack should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("pack voltage"));
        assert!(s.result.contains("pack energy"));
        assert!(s.result.contains("runtime"));
        // 13S of 3.7 V = 48.1 V; 13*4 = 52 cells.
        assert!(s.result.contains("48.10"));
        assert!(s.result.contains("52 cells"));
    }

    #[test]
    fn analyze_rejects_zero_series() {
        let mut s = BatteryPackWorkbenchState {
            series: 0,
            ..Default::default()
        };
        run_pack(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn pack_energy_matches_cell_count_times_cell_energy() {
        // Ground truth: Wh_pack = n * m * V_cell * Q_cell, and the
        // cell-by-cell sum must agree with the topology formula.
        let cell = Cell::new(3.7, 3.0).unwrap();
        let pack = PackConfig::new(cell, 13, 4).unwrap();
        let expected = 13.0 * 4.0 * 3.7 * 3.0;
        assert!((pack.pack_energy_wh() - expected).abs() < 1e-9);
        assert!((pack.total_cell_energy_wh() - pack.pack_energy_wh()).abs() < 1e-9);
    }

    #[test]
    fn pack_mesh_for_default_is_nonempty_and_in_range() {
        let s = BatteryPackWorkbenchState::default();
        let mesh = pack_solid_mesh(&s).expect("default pack yields a solid");
        assert!(mesh.nodes.len() > 8, "expected tray + cell grid");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn pack_mesh_none_for_invalid() {
        let s = BatteryPackWorkbenchState {
            parallel: 0,
            ..Default::default()
        };
        assert!(pack_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_batterypack_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_batterypack_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_batterypack_workbench = true;
        run_pack(&mut app.batterypack);
        draw_workbench(&mut app);
    }
}
