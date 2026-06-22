//! The right-side **Psychrometrics Workbench** panel — native moist-air
//! state analysis over `valenx-psychrometrics`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on
//! `crate::ValenxApp::show_psychrometrics_workbench`, toggled from the View
//! menu. The form sets a moist-air parcel (dry-bulb temperature, relative
//! humidity and barometric pressure); "Analyze" resolves the consistent
//! psychrometric state — saturation pressure (Magnus), partial vapour
//! pressure, humidity ratio `W = 0.622 pv / (p - pv)`, dew point, degree of
//! saturation and specific enthalpy `h = 1.006 T + W (2501 + 1.86 T)` — and
//! "Show 3-D" loads a representative air-volume control box into the central
//! viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_psychrometrics::moist_air::degree_of_saturation;
use valenx_psychrometrics::state::MoistAirState;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Which derived quantity headlines the 3-D control box request. Cosmetic
/// only — the box geometry is identical; the choice records what the user is
/// inspecting and is echoed in the readout.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ParcelView {
    /// Highlight the air volume itself (the control box).
    Volume,
    /// Highlight the moisture content (humidity ratio).
    Moisture,
}

/// Persistent form + result state for the Psychrometrics Workbench.
pub struct PsychrometricsWorkbenchState {
    /// Dry-bulb temperature (deg C).
    dry_bulb_c: f64,
    /// Relative humidity (percent, `0..=100`).
    relative_humidity_pct: f64,
    /// Total (barometric) pressure (kPa).
    pressure_kpa: f64,
    /// Which quantity the 3-D request highlights.
    parcel_view: ParcelView,
    /// Formatted psychrometric readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D control box (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for PsychrometricsWorkbenchState {
    fn default() -> Self {
        // A comfortable indoor parcel at one standard atmosphere:
        // 25 degC, 50% RH, 101.325 kPa -> W ~ 0.0099 kg/kg, dew point
        // ~13.8 degC, enthalpy ~50.4 kJ/kg.
        Self {
            dry_bulb_c: 25.0,
            relative_humidity_pct: 50.0,
            pressure_kpa: 101.325,
            parcel_view: ParcelView::Volume,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Psychrometrics Workbench right-side panel. A no-op when the
/// `show_psychrometrics_workbench` toggle is off.
pub fn draw_psychrometrics_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_psychrometrics_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_psychrometrics_workbench",
        "Psychrometrics",
        |app, ui| {
            ui.label(
                egui::RichText::new("native moist-air state (Magnus) · valenx-psychrometrics")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.psychrometrics;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Moist-air parcel").strong());
                    ui.horizontal(|ui| {
                        ui.label("dry-bulb (°C)");
                        ui.add(egui::DragValue::new(&mut s.dry_bulb_c).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("relative humidity (%)");
                        ui.add(
                            egui::DragValue::new(&mut s.relative_humidity_pct)
                                .speed(1.0)
                                .range(0.0..=100.0),
                        );
                    });
                    ui.horizontal(|ui| {
                        ui.label("pressure (kPa)");
                        ui.add(egui::DragValue::new(&mut s.pressure_kpa).speed(0.5));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("3-D highlight").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.parcel_view, ParcelView::Volume, "air volume");
                        ui.radio_value(&mut s.parcel_view, ParcelView::Moisture, "moisture");
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_psychrometrics(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative air-volume control box for the parcel as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Psychrometric state").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_psychrometrics_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.psychrometrics` borrow is
    // released here): build the parcel's 3-D control box and load it.
    if app.psychrometrics.show_3d_request {
        app.psychrometrics.show_3d_request = false;
        load_parcel_3d(app);
    }
}

/// Validate the form, resolve the moist-air state and format the readout.
fn run_psychrometrics(s: &mut PsychrometricsWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Resolve the consistent moist-air state and format the full readout,
/// mapping any domain error to a display string. Extracted so it is
/// unit-testable.
fn compute(s: &PsychrometricsWorkbenchState) -> Result<String, String> {
    let rh = s.relative_humidity_pct / 100.0;
    let p_pa = s.pressure_kpa * 1000.0;
    let state =
        MoistAirState::from_relative_humidity(s.dry_bulb_c, rh, p_pa).map_err(|e| e.to_string())?;
    // Degree of saturation is not stored on the state — derive it from the
    // resolved humidity ratio at the same dry-bulb temperature and pressure.
    let mu = degree_of_saturation(state.humidity_ratio, state.dry_bulb_c, state.pressure_pa)
        .map_err(|e| e.to_string())?;

    Ok(format!(
        "dry-bulb        : {:.2} °C\n\
         pressure        : {:.3} kPa\n\
         relative humid. : {:.1} %\n\n\
         saturation pres.: {:.1} Pa\n\
         vapour pressure : {:.1} Pa\n\
         humidity ratio W: {:.5} kg/kg\n\
         dew point       : {:.2} °C\n\
         degree of sat.  : {:.4}\n\
         enthalpy h      : {:.2} kJ/kg",
        state.dry_bulb_c,
        state.pressure_pa / 1000.0,
        state.relative_humidity * 100.0,
        state.saturation_pressure_pa,
        state.vapour_pressure_pa,
        state.humidity_ratio,
        state.dew_point_c,
        mu,
        state.enthalpy_kj_per_kg,
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

/// Build the air parcel as a triangle [`Mesh`] — a unit air-volume control
/// box on a thin base. Representative geometry only (the psychrometric
/// numbers are the `valenx-psychrometrics` result); `None` for an invalid
/// configuration (e.g. relative humidity outside `[0, 100]`, a non-positive
/// pressure, or a temperature below absolute zero).
fn parcel_solid_mesh(s: &PsychrometricsWorkbenchState) -> Option<Mesh> {
    let rh = s.relative_humidity_pct / 100.0;
    let p_pa = s.pressure_kpa * 1000.0;
    MoistAirState::from_relative_humidity(s.dry_bulb_c, rh, p_pa).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Air-volume control box.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.6),
        Vector3::new(0.5, 0.5, 0.5),
    );
    // Thin base.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.05),
        Vector3::new(0.6, 0.6, 0.05),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-psychrometrics");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D control box and load it into the central viewport.
fn load_parcel_3d(app: &mut ValenxApp) {
    let Some(mesh) = parcel_solid_mesh(&app.psychrometrics) else {
        app.psychrometrics.error =
            Some("parcel parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<parcel>/valenx-psychrometrics"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// Agent-bridge product: a DATA-ONLY text card of the psychrometrics
/// workbench's `compute()` readout rows (see [`crate::products_registry`]). A
/// moist-air state has no characteristic shape — the panel's air-parcel
/// control box is a generic stand-in, not a real object — so the bridge
/// product is right-sized to a card (`mesh: None`) carrying just the readout
/// (the confidence badge is appended centrally). The panel's "Show 3-D" button
/// still builds that representative box into the central viewport.
pub(crate) fn psychrometrics_product() -> crate::WorkspaceProduct {
    let s = PsychrometricsWorkbenchState::default();
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical psychrometrics ⇒ readout computes"),
    );
    crate::WorkspaceProduct {
        title: "Psychrometrics (humidity/enthalpy)".into(),
        lines,
        mesh: None,
        vertex_colors: None,
        camera: Default::default(),
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
        let s = PsychrometricsWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_full_state() {
        let mut s = PsychrometricsWorkbenchState::default();
        run_psychrometrics(&mut s);
        assert!(
            s.error.is_none(),
            "default parcel should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("humidity ratio W"));
        assert!(s.result.contains("dew point"));
        assert!(s.result.contains("enthalpy h"));
        // 25 degC / 50% RH: dew point ~13.8 degC.
        assert!(s.result.contains("13.8"));
    }

    #[test]
    fn analyze_rejects_relative_humidity_above_full() {
        let mut s = PsychrometricsWorkbenchState {
            relative_humidity_pct: 150.0,
            ..Default::default()
        };
        run_psychrometrics(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn humidity_ratio_matches_hand_computed_formula() {
        // Ground truth: W = 0.622 * pv / (p - pv) for the default parcel.
        // Resolve the state, then recompute W by hand from its own vapour
        // and total pressures and confirm they agree.
        let s = PsychrometricsWorkbenchState::default();
        let rh = s.relative_humidity_pct / 100.0;
        let p_pa = s.pressure_kpa * 1000.0;
        let state = MoistAirState::from_relative_humidity(s.dry_bulb_c, rh, p_pa).unwrap();
        let pv = state.vapour_pressure_pa;
        let expect = 0.622 * pv / (p_pa - pv);
        let diff = (state.humidity_ratio - expect).abs();
        assert!(diff < 1e-12, "W {} vs {expect}", state.humidity_ratio);
    }

    #[test]
    fn parcel_mesh_for_default_is_nonempty_and_in_range() {
        let s = PsychrometricsWorkbenchState::default();
        let mesh = parcel_solid_mesh(&s).expect("default parcel yields a solid");
        assert!(mesh.nodes.len() > 8, "expected control box + base");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn parcel_mesh_none_for_invalid() {
        let s = PsychrometricsWorkbenchState {
            relative_humidity_pct: 150.0,
            ..Default::default()
        };
        assert!(parcel_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_psychrometrics_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_psychrometrics_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_psychrometrics_workbench = true;
        run_psychrometrics(&mut app.psychrometrics);
        draw_workbench(&mut app);
    }
}
