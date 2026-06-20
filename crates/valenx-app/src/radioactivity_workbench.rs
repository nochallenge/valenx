//! The right-side **Radioactive Decay Workbench** panel — native
//! single-nuclide exponential-decay analysis over `valenx-radioactivity`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_radioactivity_workbench`,
//! toggled from the View menu. The form sets a nuclide by its half-life
//! plus an initial population `N0` and an elapsed time `t`; "Analyze"
//! builds the [`Nuclide`] (`lambda = ln(2) / t_half`) and reports the decay
//! constant, mean life, the surviving population `N(t) = N0 * exp(-lambda *
//! t)`, the remaining fraction, the initial and current activity
//! `A = lambda * N`, and the number of decays over `[0, t]`. "Show 3-D"
//! loads a representative swept decay-curve ribbon (height tracing the
//! `N(t)` exponential) into the central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_radioactivity::Nuclide;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Radioactive Decay Workbench.
pub struct RadioactivityWorkbenchState {
    /// Half-life `t_half` (in the same time unit as the elapsed time;
    /// these defaults are in days).
    half_life: f64,
    /// Initial number of un-decayed nuclei `N0`.
    n0: f64,
    /// Elapsed time `t` at which to evaluate the decay.
    elapsed: f64,
    /// Formatted decay readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D decay-curve solid (serviced
    /// after the panel draws).
    show_3d_request: bool,
}

impl Default for RadioactivityWorkbenchState {
    fn default() -> Self {
        // Iodine-131 (half-life ~ 8.02 days), the crate's own worked
        // example: a million-nucleus sample observed for one half-life, so
        // the surviving fraction is ~ 0.5 and ~ half the sample decays.
        Self {
            half_life: 8.02,
            n0: 1.0e6,
            elapsed: 8.02,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Radioactive Decay Workbench right-side panel. A no-op when the
/// `show_radioactivity_workbench` toggle is off.
pub fn draw_radioactivity_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_radioactivity_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_radioactivity_workbench",
        "Radioactive Decay",
        |app, ui| {
            ui.label(
                egui::RichText::new("native single-nuclide exponential decay · valenx-radioactivity")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.radioactivity;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Nuclide").strong());
                    ui.horizontal(|ui| {
                        ui.label("half-life t½ (time)");
                        ui.add(egui::DragValue::new(&mut s.half_life).speed(0.1));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Sample").strong());
                    ui.horizontal(|ui| {
                        ui.label("initial N₀ (nuclei)");
                        ui.add(egui::DragValue::new(&mut s.n0).speed(1.0e4));
                    });
                    ui.horizontal(|ui| {
                        ui.label("elapsed t (time)");
                        ui.add(egui::DragValue::new(&mut s.elapsed).speed(0.1));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_radioactivity(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative swept ribbon tracing the N(t) = N0·exp(−λt) decay curve as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Decay").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_radioactivity_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.radioactivity` borrow
    // is released here): build the decay curve's 3-D solid and load it.
    if app.radioactivity.show_3d_request {
        app.radioactivity.show_3d_request = false;
        load_decay_3d(app);
    }
}

/// Validate the form, evaluate the decay and format the readout.
fn run_radioactivity(s: &mut RadioactivityWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build the [`Nuclide`] from the form's half-life — the single value the
/// readout and the 3-D gate both need. Extracted so it is unit-testable
/// and shared.
fn nuclide(s: &RadioactivityWorkbenchState) -> Result<Nuclide, String> {
    Nuclide::from_half_life(s.half_life).map_err(|e| e.to_string())
}

/// Evaluate the single-nuclide decay and format the full readout, mapping
/// any domain error to a display string. Extracted so it is unit-testable.
fn compute(s: &RadioactivityWorkbenchState) -> Result<String, String> {
    let nuc = nuclide(s)?;
    let lambda = nuc.decay_constant();
    let mean_life = nuc.mean_life();
    let n_t = nuc.remaining(s.n0, s.elapsed).map_err(|e| e.to_string())?;
    let fraction = nuc
        .remaining_fraction(s.elapsed)
        .map_err(|e| e.to_string())?;
    let a0 = nuc.activity(s.n0).map_err(|e| e.to_string())?;
    let a_t = nuc
        .activity_at(s.n0, s.elapsed)
        .map_err(|e| e.to_string())?;
    let decayed = nuc
        .decays_in_interval(s.n0, 0.0, s.elapsed)
        .map_err(|e| e.to_string())?;

    Ok(format!(
        "half-life t½    : {:.4} time\n\
         decay const λ   : {:.6} /time\n\
         mean life τ     : {:.4} time\n\
         initial N₀      : {:.4e} nuclei\n\
         elapsed t       : {:.4} time\n\n\
         remaining N(t)  : {:.4e} nuclei\n\
         fraction N/N₀   : {:.6}\n\
         decayed in [0,t]: {:.4e} nuclei\n\
         activity A₀     : {:.4e} /time\n\
         activity A(t)   : {:.4e} /time",
        s.half_life, lambda, mean_life, s.n0, s.elapsed, n_t, fraction, decayed, a0, a_t,
    ))
}

/// Build the decay curve `N(t) = N0 * exp(-lambda * t)` as a swept-ribbon
/// triangle [`Mesh`]: a normalised exponential profile sampled out to a
/// few half-lives, extruded a little in the depth (`y`) direction into a
/// flat ribbon. Representative geometry (the decay numbers themselves are
/// the `valenx-radioactivity` result). `None` for an invalid nuclide.
fn decay_curve_mesh(s: &RadioactivityWorkbenchState) -> Option<Mesh> {
    let nuc = nuclide(s).ok()?;

    // Sample the normalised curve N(t)/N0 = exp(-lambda*t) across five
    // half-lives so the full grow-down shape is visible regardless of the
    // chosen elapsed time.
    let segments = 48usize;
    let span = 5.0 * nuc.half_life();
    let width = 4.0_f64; // total horizontal extent of the plotted axis
    let height = 2.0_f64; // height of the curve at t = 0 (fraction == 1)
    let depth = 0.18_f64; // half-thickness of the ribbon in y

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    for i in 0..=segments {
        let frac = i as f64 / segments as f64;
        let t = frac * span;
        // remaining_fraction only fails on a non-finite/negative t, which
        // cannot happen here; clamp into [0, 1] for geometric safety.
        let y_curve = nuc.remaining_fraction(t).unwrap_or(0.0).clamp(0.0, 1.0);
        let x = frac * width;
        let z = y_curve * height;
        // Front (y = -depth) and back (y = +depth) rails at this station.
        nodes.push(Vector3::new(x, -depth, z));
        nodes.push(Vector3::new(x, depth, z));
        // Baseline rails (z = 0) so the ribbon reads as a filled area.
        nodes.push(Vector3::new(x, -depth, 0.0));
        nodes.push(Vector3::new(x, depth, 0.0));
    }

    // Stitch consecutive stations: the curved top rail and the flat base,
    // each as a pair of triangles, giving a swept filled ribbon.
    for i in 0..segments {
        let a = 4 * i; // this station's four nodes: top-front/back, base-front/back
        let b = 4 * (i + 1);
        // Top surface (curve rail).
        tris.extend_from_slice(&[a, b, b + 1, a, b + 1, a + 1]);
        // Base surface.
        tris.extend_from_slice(&[a + 2, a + 3, b + 3, a + 2, b + 3, b + 2]);
        // Front wall (y = -depth) closing curve-top to base.
        tris.extend_from_slice(&[a, a + 2, b + 2, a, b + 2, b]);
        // Back wall (y = +depth).
        tris.extend_from_slice(&[a + 1, b + 1, b + 3, a + 1, b + 3, a + 3]);
    }

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-radioactivity");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D decay-curve solid and load it into the central viewport.
fn load_decay_3d(app: &mut ValenxApp) {
    let Some(mesh) = decay_curve_mesh(&app.radioactivity) else {
        app.radioactivity.error =
            Some("nuclide parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<decay>/valenx-radioactivity"),
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
        let s = RadioactivityWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_decay_quantities() {
        let mut s = RadioactivityWorkbenchState::default();
        run_radioactivity(&mut s);
        assert!(
            s.error.is_none(),
            "default nuclide should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("decay const λ"));
        assert!(s.result.contains("remaining N(t)"));
        assert!(s.result.contains("activity A(t)"));
        // The default elapsed time is exactly one half-life, so the
        // surviving fraction is 0.500000 to compute()'s {:.6} precision.
        assert!(s.result.contains("0.500000"));
    }

    #[test]
    fn analyze_rejects_nonpositive_half_life() {
        let mut s = RadioactivityWorkbenchState {
            half_life: 0.0,
            ..Default::default()
        };
        run_radioactivity(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn population_halves_after_one_half_life_ground_truth() {
        // GROUND TRUTH: lambda = ln(2)/t_half, and N(t_half) = N0/2.
        let t_half = 8.0_f64;
        let nuc = Nuclide::from_half_life(t_half).unwrap();
        // Hand-computed decay constant.
        let lambda_hand = std::f64::consts::LN_2 / t_half;
        assert!((nuc.decay_constant() - lambda_hand).abs() < 1e-15);
        // After exactly one half-life half the sample remains.
        let n0 = 1.0e6;
        let n = nuc.remaining(n0, t_half).unwrap();
        assert!((n - n0 / 2.0).abs() < 1e-6, "N(t_half) = {n}");
    }

    #[test]
    fn decay_mesh_for_default_is_nonempty_and_in_range() {
        let s = RadioactivityWorkbenchState::default();
        let mesh = decay_curve_mesh(&s).expect("default nuclide yields a solid");
        assert!(mesh.nodes.len() > 8, "expected a multi-station ribbon");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn decay_mesh_none_for_invalid() {
        let s = RadioactivityWorkbenchState {
            half_life: 0.0,
            ..Default::default()
        };
        assert!(decay_curve_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_radioactivity_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_radioactivity_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_radioactivity_workbench = true;
        run_radioactivity(&mut app.radioactivity);
        draw_workbench(&mut app);
    }
}
