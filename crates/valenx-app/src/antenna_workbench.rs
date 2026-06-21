//! The right-side **Antenna Workbench** panel — native parabolic-dish
//! gain / beamwidth analysis over `valenx-antenna`.
//!
//! Mirrors the Gearbox / Induction Motor workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_antenna_workbench`,
//! toggled from the View menu. The form sets the operating frequency, the
//! dish diameter and the aperture efficiency; "Analyze" reports the
//! wavelength, the effective aperture, the gain (linear and dBi) and a
//! half-power beamwidth estimate, and "Show 3-D dish" loads a
//! representative parabolic reflector solid into the central viewport.

use std::f64::consts::{PI, TAU};
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_antenna::beamwidth::{beamwidth_k_rad, radians_to_degrees, K_UNIFORM_LINE_HPBW};
use valenx_antenna::gain::{effective_aperture, gain_from_aperture, to_dbi};
use valenx_antenna::wave::wavelength_from_frequency;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Antenna Workbench.
pub struct AntennaWorkbenchState {
    /// Operating frequency `f` (GHz).
    frequency_ghz: f64,
    /// Dish diameter `D` (m).
    diameter_m: f64,
    /// Aperture efficiency `eta` in (0, 1].
    efficiency: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D dish solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for AntennaWorkbenchState {
    fn default() -> Self {
        // A 1.2 m Ku-band (12 GHz) dish at 60% aperture efficiency —
        // ~41 dBi gain, ~1 deg half-power beamwidth.
        Self {
            frequency_ghz: 12.0,
            diameter_m: 1.2,
            efficiency: 0.6,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Antenna Workbench right-side panel. A no-op when the
/// `show_antenna_workbench` toggle is off.
pub fn draw_antenna_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_antenna_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_antenna_workbench",
        "Antenna",
        |app, ui| {
            ui.label(
                egui::RichText::new("native parabolic-dish gain / beamwidth · valenx-antenna")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.antenna;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Parabolic dish").strong());
                    ui.horizontal(|ui| {
                        ui.label("frequency (GHz)");
                        ui.add(egui::DragValue::new(&mut s.frequency_ghz).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("diameter D (m)");
                        ui.add(egui::DragValue::new(&mut s.diameter_m).speed(0.05));
                    });
                    ui.horizontal(|ui| {
                        ui.label("aperture efficiency");
                        ui.add(egui::DragValue::new(&mut s.efficiency).speed(0.01));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_antenna(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D dish").strong())
                        .on_hover_text(
                            "Build a representative parabolic dish (reflector, feed boom, feed and mount) as a 3-D solid and load it into the central viewport to orbit",
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
        app.show_antenna_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.antenna` borrow is
    // released here): build the dish's 3-D solid and load it.
    if app.antenna.show_3d_request {
        app.antenna.show_3d_request = false;
        load_dish_3d(app);
    }
}

/// Validate the form, evaluate the dish and format the readout.
fn run_antenna(s: &mut AntennaWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// The wavelength and effective aperture for the form, the two quantities
/// the gain and the 3-D gate both need. Extracted so it is unit-testable
/// and shared.
fn wavelength_and_aperture(s: &AntennaWorkbenchState) -> Result<(f64, f64), String> {
    let lambda = wavelength_from_frequency(s.frequency_ghz * 1.0e9).map_err(|e| e.to_string())?;
    let radius = s.diameter_m / 2.0;
    let physical_area = PI * radius * radius;
    let aeff = effective_aperture(physical_area, s.efficiency).map_err(|e| e.to_string())?;
    Ok((lambda, aeff))
}

/// Evaluate the dish and format the full readout, mapping any domain error
/// to a display string. Extracted so it is unit-testable.
fn compute(s: &AntennaWorkbenchState) -> Result<String, String> {
    let (lambda, aeff) = wavelength_and_aperture(s)?;
    let gain_lin = gain_from_aperture(aeff, lambda).map_err(|e| e.to_string())?;
    let gain_dbi = to_dbi(gain_lin).map_err(|e| e.to_string())?;
    let hpbw_rad =
        beamwidth_k_rad(K_UNIFORM_LINE_HPBW, lambda, s.diameter_m).map_err(|e| e.to_string())?;
    let hpbw_deg = radians_to_degrees(hpbw_rad).map_err(|e| e.to_string())?;
    let physical_area = PI * (s.diameter_m / 2.0).powi(2);

    Ok(format!(
        "frequency       : {:.2} GHz\n\
         wavelength      : {:.1} mm\n\
         dish diameter   : {:.2} m\n\
         aperture eff.   : {:.0} %\n\n\
         physical area   : {:.3} m²\n\
         effective aper. : {:.3} m²\n\
         gain            : {:.2} dBi  (×{:.0})\n\
         HPBW estimate   : {:.2} °  (0.886·λ/D)",
        s.frequency_ghz,
        lambda * 1000.0,
        s.diameter_m,
        s.efficiency * 100.0,
        physical_area,
        aeff,
        gain_dbi,
        gain_lin,
        hpbw_deg,
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

/// Build the dish antenna as a triangle [`Mesh`] — a circular reflector
/// disc facing `+x`, a feed boom out to the focus carrying a feed, and a
/// rear mount on a base. Representative geometry (the reflector is a flat
/// disc, not a true paraboloid; the gain / beamwidth are the
/// `valenx-antenna` result). `None` for an invalid configuration.
fn dish_solid_mesh(s: &AntennaWorkbenchState) -> Option<Mesh> {
    wavelength_and_aperture(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();
    let axis_z = 0.7;

    // Reflector disc.
    push_cyl_x(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, axis_z),
        0.06,
        0.5,
        32,
    );
    // Feed boom from the dish out to the focus (+x).
    push_cyl_x(
        &mut nodes,
        &mut tris,
        Vector3::new(0.06, 0.0, axis_z),
        0.45,
        0.025,
        12,
    );
    // Feed horn at the focus.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.5, 0.0, axis_z),
        Vector3::new(0.05, 0.05, 0.05),
    );
    // Rear mount column.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(-0.15, 0.0, 0.42),
        Vector3::new(0.08, 0.08, 0.22),
    );
    // Base.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(-0.15, 0.0, 0.06),
        Vector3::new(0.3, 0.3, 0.06),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-antenna");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D dish solid and load it into the central viewport.
fn load_dish_3d(app: &mut ValenxApp) {
    let Some(mesh) = dish_solid_mesh(&app.antenna) else {
        app.antenna.error = Some("dish parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<dish>/valenx-antenna"),
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
        let s = AntennaWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_gain_and_beamwidth() {
        let mut s = AntennaWorkbenchState::default();
        run_antenna(&mut s);
        assert!(
            s.error.is_none(),
            "default dish should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("wavelength"));
        assert!(s.result.contains("dBi"));
        assert!(s.result.contains("HPBW"));
        // A 1.2 m, 60% dish at 12 GHz sits around 41 dBi.
        assert!(s.result.contains("41."));
    }

    #[test]
    fn analyze_rejects_zero_diameter() {
        let mut s = AntennaWorkbenchState {
            diameter_m: 0.0,
            ..Default::default()
        };
        run_antenna(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn gain_and_aperture_are_inverse() {
        // Ground truth: gain_from_aperture and aperture_from_gain are
        // exact inverses — G = 4*pi*Ae/lambda^2 and Ae = G*lambda^2/(4*pi).
        let lambda = 0.025;
        let aeff = 0.68;
        let g = gain_from_aperture(aeff, lambda).unwrap();
        let back = valenx_antenna::gain::aperture_from_gain(g, lambda).unwrap();
        assert!((back - aeff).abs() < 1e-9, "round-trip {back} != {aeff}");
        // And the closed form itself.
        let expected = 4.0 * PI * aeff / (lambda * lambda);
        assert!((g - expected).abs() < 1e-6);
    }

    #[test]
    fn dish_mesh_for_default_is_nonempty_and_in_range() {
        let s = AntennaWorkbenchState::default();
        let mesh = dish_solid_mesh(&s).expect("default dish yields a solid");
        assert!(
            mesh.nodes.len() > 8,
            "expected reflector + boom + feed + mount"
        );
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn dish_mesh_none_for_invalid() {
        let s = AntennaWorkbenchState {
            diameter_m: 0.0,
            ..Default::default()
        };
        assert!(dish_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_antenna_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_antenna_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_antenna_workbench = true;
        run_antenna(&mut app.antenna);
        draw_workbench(&mut app);
    }
}
