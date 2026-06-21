//! The right-side **Diffusion Workbench** panel — native 1-D Fickian
//! diffusion analysis over `valenx-diffusion`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_diffusion_workbench`,
//! toggled from the View menu. The form picks one of two closed-form
//! modes — **Fick-1 flux** (the steady through-slab flux `J = -D dC/dx`
//! from a fixed concentration drop across a slab) or **Gaussian spread**
//! (the instantaneous point-source kernel `C = M / sqrt(4 pi D t) exp(...)`
//! and its `var = 2 D t` spreading) — and "Analyze" reports the result.
//! "Show 3-D" loads a representative Gaussian-spread bell surface into the
//! central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_diffusion::{
    distance_for_concentration_fraction, first_law_flux, gaussian_point_source, gaussian_std,
    gaussian_variance, steady_flux, steady_gradient, time_to_reach_std,
};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Which closed-form diffusion model the workbench evaluates.
#[derive(Clone, Copy, PartialEq, Eq)]
enum DiffusionMode {
    /// Fick's first law: the steady through-slab flux from a fixed
    /// concentration drop across a slab of finite length.
    Flux,
    /// The instantaneous point-source spreading: the Gaussian kernel and
    /// its `var = 2 D t` growth.
    Spread,
}

/// Persistent form + result state for the Diffusion Workbench.
pub struct DiffusionWorkbenchState {
    /// Which closed-form model to evaluate.
    mode: DiffusionMode,
    /// Diffusion coefficient `D` (m^2/s).
    diffusivity_m2_s: f64,
    /// Concentration at the left wall (Flux mode), arbitrary units.
    c_left: f64,
    /// Concentration at the right wall (Flux mode), arbitrary units.
    c_right: f64,
    /// Slab length `L` (m), the span the drop occurs over (Flux mode).
    length_m: f64,
    /// Released mass per unit area `M` at the origin (Spread mode).
    mass: f64,
    /// Elapsed time `t` since release (s), used by the Spread mode.
    time_s: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D spread surface (serviced after
    /// the panel draws).
    show_3d_request: bool,
}

impl Default for DiffusionWorkbenchState {
    fn default() -> Self {
        // A small-molecule tracer in water (D ~ 1e-9 m^2/s). In Flux mode
        // a unit concentration drop over a 1 mm slab. In Spread mode a
        // unit mass spreading for 1 hour reaches sigma = sqrt(2 D t) ~
        // 2.7e-3 m.
        Self {
            mode: DiffusionMode::Spread,
            diffusivity_m2_s: 1.0e-9,
            c_left: 1.0,
            c_right: 0.0,
            length_m: 1.0e-3,
            mass: 1.0,
            time_s: 3600.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Diffusion Workbench right-side panel. A no-op when the
/// `show_diffusion_workbench` toggle is off.
pub fn draw_diffusion_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_diffusion_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_diffusion_workbench",
        "Diffusion",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native 1-D Fickian diffusion (Fick-1 flux / Gaussian spread) · valenx-diffusion",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.diffusion;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Model").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.mode, DiffusionMode::Flux, "Fick-1 flux");
                        ui.radio_value(&mut s.mode, DiffusionMode::Spread, "Gaussian spread");
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Medium").strong());
                    ui.horizontal(|ui| {
                        ui.label("diffusivity D (m²/s)");
                        ui.add(
                            egui::DragValue::new(&mut s.diffusivity_m2_s)
                                .speed(1.0e-10)
                                .max_decimals(12),
                        );
                    });

                    match s.mode {
                        DiffusionMode::Flux => {
                            ui.add_space(4.0);
                            ui.label(egui::RichText::new("Slab").strong());
                            ui.horizontal(|ui| {
                                ui.label("C left");
                                ui.add(egui::DragValue::new(&mut s.c_left).speed(0.05));
                            });
                            ui.horizontal(|ui| {
                                ui.label("C right");
                                ui.add(egui::DragValue::new(&mut s.c_right).speed(0.05));
                            });
                            ui.horizontal(|ui| {
                                ui.label("length L (m)");
                                ui.add(
                                    egui::DragValue::new(&mut s.length_m)
                                        .speed(1.0e-4)
                                        .max_decimals(6),
                                );
                            });
                        }
                        DiffusionMode::Spread => {
                            ui.add_space(4.0);
                            ui.label(egui::RichText::new("Point source").strong());
                            ui.horizontal(|ui| {
                                ui.label("mass M (per area)");
                                ui.add(egui::DragValue::new(&mut s.mass).speed(0.05));
                            });
                            ui.horizontal(|ui| {
                                ui.label("time t (s)");
                                ui.add(egui::DragValue::new(&mut s.time_s).speed(60.0));
                            });
                        }
                    }

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_diffusion(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build the representative Gaussian point-source concentration bell C(x, t) as a 3-D surface and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Result").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_diffusion_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.diffusion` borrow is
    // released here): build the spread surface and load it.
    if app.diffusion.show_3d_request {
        app.diffusion.show_3d_request = false;
        load_diffusion_3d(app);
    }
}

/// Validate the form, evaluate the selected model and format the readout.
fn run_diffusion(s: &mut DiffusionWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Evaluate the selected model and format the full readout, mapping any
/// domain error to a display string. Extracted so it is unit-testable.
fn compute(s: &DiffusionWorkbenchState) -> Result<String, String> {
    match s.mode {
        DiffusionMode::Flux => {
            let grad =
                steady_gradient(s.c_left, s.c_right, s.length_m).map_err(|e| e.to_string())?;
            let flux = steady_flux(s.diffusivity_m2_s, s.c_left, s.c_right, s.length_m)
                .map_err(|e| e.to_string())?;
            // Cross-check via the explicit first-law form J = -D * grad.
            let flux_check = first_law_flux(s.diffusivity_m2_s, grad).map_err(|e| e.to_string())?;
            Ok(format!(
                "mode            : Fick-1 flux (steady slab)\n\
                 diffusivity D   : {:.3e} m²/s\n\
                 C left / right  : {:.3} / {:.3}\n\
                 slab length L   : {:.3e} m\n\n\
                 gradient dC/dx  : {:.3e} /m\n\
                 flux J = -D∇C   : {:.3e} (per area·s)\n\
                 first-law check : {:.3e}",
                s.diffusivity_m2_s, s.c_left, s.c_right, s.length_m, grad, flux, flux_check,
            ))
        }
        DiffusionMode::Spread => {
            let var = gaussian_variance(s.diffusivity_m2_s, s.time_s).map_err(|e| e.to_string())?;
            let sigma = gaussian_std(s.diffusivity_m2_s, s.time_s).map_err(|e| e.to_string())?;
            let peak = gaussian_point_source(s.mass, s.diffusivity_m2_s, 0.0, s.time_s)
                .map_err(|e| e.to_string())?;
            // Distance out to the half-peak (1/2) concentration contour.
            let x_half = distance_for_concentration_fraction(s.diffusivity_m2_s, s.time_s, 0.5)
                .map_err(|e| e.to_string())?;
            // Inverse: time for the RMS spread to reach one sigma again
            // (a self-consistency readout of the std law).
            let t_back = time_to_reach_std(s.diffusivity_m2_s, sigma).map_err(|e| e.to_string())?;
            Ok(format!(
                "mode            : Gaussian point source\n\
                 diffusivity D   : {:.3e} m²/s\n\
                 mass M          : {:.3}\n\
                 time t          : {:.1} s\n\n\
                 variance 2 D t  : {:.3e} m²\n\
                 RMS spread σ    : {:.3e} m\n\
                 peak C(0, t)    : {:.3e} (per area)\n\
                 half-peak x@½   : {:.3e} m\n\
                 t to reach σ    : {:.1} s",
                s.diffusivity_m2_s, s.mass, s.time_s, var, sigma, peak, x_half, t_back,
            ))
        }
    }
}

/// Build the representative Gaussian point-source bell `C(x, t)` as a
/// triangle [`Mesh`] surface: a height field `z = C(x, t)` over a square
/// `(x, y)` patch, the `y` axis carrying the same profile so the dome is
/// radially suggestive. Heights are normalised to the peak and scaled for
/// a viewable solid; the spread `sigma = sqrt(2 D t)` sets the footprint.
/// `None` for an invalid configuration.
fn spread_surface_mesh(s: &DiffusionWorkbenchState) -> Option<Mesh> {
    let sigma = gaussian_std(s.diffusivity_m2_s, s.time_s).ok()?;
    let peak = gaussian_point_source(s.mass, s.diffusivity_m2_s, 0.0, s.time_s).ok()?;
    // A degenerate footprint (t = 0 -> sigma = 0) or a flat (zero-mass)
    // bell carries no surface worth building.
    if !(sigma.is_finite() && sigma > 0.0 && peak.is_finite() && peak.abs() > 0.0) {
        return None;
    }

    // A grid spanning +/- 3 sigma, normalised so the patch is order-unity
    // regardless of the physical scale (D, t can be tiny).
    let half = 3.0 * sigma;
    let n = 33_usize; // odd, so a node sits on the peak
    let span = 2.0 * half;
    let height = 0.8; // viewable dome height

    let mut nodes: Vec<Vector3<f64>> = Vec::with_capacity(n * n);
    for j in 0..n {
        let y = -half + span * j as f64 / (n as f64 - 1.0);
        for i in 0..n {
            let x = -half + span * i as f64 / (n as f64 - 1.0);
            // Radial Gaussian bell normalised to the peak, in [0, 1].
            let r2 = x * x + y * y;
            let bell = (-r2 / (2.0 * sigma * sigma)).exp();
            let nx = (x / half).clamp(-1.0, 1.0);
            let ny = (y / half).clamp(-1.0, 1.0);
            nodes.push(Vector3::new(nx, ny, height * bell));
        }
    }

    let mut tris: Vec<usize> = Vec::with_capacity((n - 1) * (n - 1) * 6);
    for j in 0..n - 1 {
        for i in 0..n - 1 {
            let a = j * n + i;
            let b = j * n + i + 1;
            let c = (j + 1) * n + i;
            let d = (j + 1) * n + i + 1;
            tris.extend_from_slice(&[a, b, d, a, d, c]);
        }
    }

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-diffusion");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D spread surface and load it into the central viewport.
fn load_diffusion_3d(app: &mut ValenxApp) {
    let Some(mesh) = spread_surface_mesh(&app.diffusion) else {
        app.diffusion.error =
            Some("diffusion parameters are invalid — cannot build the 3-D surface".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<spread>/valenx-diffusion"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"diffusion"}`** product: the canonical
/// 1-D Fickian diffusion concentration profile built as a 3-D surface, paired
/// with the workbench's own `compute()` readout rows, at a fixed 3/4 camera.
/// Registered in [`crate::products_registry`]; the per-tool builder the
/// registry dispatches to. Pure — driven off
/// [`DiffusionWorkbenchState::default`].
pub(crate) fn diffusion_product() -> crate::WorkspaceProduct {
    let s = DiffusionWorkbenchState::default();
    let mesh = spread_surface_mesh(&s).expect("canonical diffusion ⇒ spread surface builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<diffusion>/valenx-diffusion");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical diffusion ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Diffusion (Fickian 1-D)".into(),
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
        let s = DiffusionWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_spread_reports_variance_and_spread() {
        let mut s = DiffusionWorkbenchState::default();
        run_diffusion(&mut s);
        assert!(
            s.error.is_none(),
            "default spread should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("Gaussian point source"));
        assert!(s.result.contains("RMS spread"));
        assert!(s.result.contains("variance 2 D t"));
        // D=1e-9, t=3600: var = 2*1e-9*3600 = 7.2e-6 m^2.
        assert!(s.result.contains("7.200e-6"));
    }

    #[test]
    fn analyze_flux_reports_flux_and_first_law_check() {
        let mut s = DiffusionWorkbenchState {
            mode: DiffusionMode::Flux,
            ..Default::default()
        };
        run_diffusion(&mut s);
        assert!(s.error.is_none(), "flux mode should analyze: {:?}", s.error);
        assert!(s.result.contains("Fick-1 flux"));
        assert!(s.result.contains("flux J"));
        assert!(s.result.contains("first-law check"));
    }

    #[test]
    fn analyze_rejects_nonpositive_diffusivity() {
        let mut s = DiffusionWorkbenchState {
            diffusivity_m2_s: 0.0,
            ..Default::default()
        };
        run_diffusion(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn spread_std_is_sqrt_two_d_t_ground_truth() {
        // Ground truth: the RMS spread of the point-source kernel is
        // sigma = sqrt(2 D t), and the steady flux is J = -D * gradient.
        let d = 2.0e-9;
        let t = 1800.0;
        let sigma = gaussian_std(d, t).unwrap();
        // Hand-computed: sqrt(2 * 2e-9 * 1800) = sqrt(7.2e-6).
        let expected_sigma = (7.2e-6_f64).sqrt();
        assert!((sigma - expected_sigma).abs() < 1e-15, "sigma = {sigma}");

        // Fick's first law: a 1.0 -> 0.0 drop over L = 1e-3 gives a
        // gradient of -1000 /m and a flux of +D*1000.
        let grad = steady_gradient(1.0, 0.0, 1.0e-3).unwrap();
        assert!((grad - (-1000.0)).abs() < 1e-9, "grad = {grad}");
        let flux = steady_flux(d, 1.0, 0.0, 1.0e-3).unwrap();
        assert!((flux - (-d * grad)).abs() < 1e-18, "flux = {flux}");
        assert!(flux > 0.0, "high-left wall drives +x flux");
    }

    #[test]
    fn spread_mesh_for_default_is_nonempty_and_in_range() {
        let s = DiffusionWorkbenchState::default();
        let mesh = spread_surface_mesh(&s).expect("default spread yields a surface");
        assert!(mesh.nodes.len() > 8, "expected a height-field grid");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn spread_mesh_none_for_invalid() {
        let s = DiffusionWorkbenchState {
            time_s: 0.0,
            ..Default::default()
        };
        assert!(spread_surface_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_diffusion_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_diffusion_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_diffusion_workbench = true;
        run_diffusion(&mut app.diffusion);
        draw_workbench(&mut app);
    }
}
