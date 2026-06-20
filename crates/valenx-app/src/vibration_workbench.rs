//! The right-side **Vibration Workbench** panel — native single-degree-of-
//! freedom (SDOF) forced-vibration analysis over `valenx-vibration`.
//!
//! Mirrors the Heat Transfer / Springs workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_vibration_workbench`,
//! toggled from the View menu. The form sets a mass-spring-damper
//! (`m`, `k`, `c`) and a harmonic forcing frequency ratio `r = w/wn`;
//! "Analyze" reports the natural frequency, damping ratio and regime, the
//! steady-state magnification factor and transmissibility at `r`, the
//! resonant peak, and the logarithmic decrement, and "Show 3-D" loads a
//! representative mass-spring-damper solid into the central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_vibration::{
    logarithmic_decrement, magnification_factor, peak_magnification, transmissibility,
    DampingRegime, SdofSystem,
};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// How the system is specified: by its three physical constants, or by the
/// modal pair `(natural frequency, damping ratio)`.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
enum SystemInput {
    /// Mass `m`, stiffness `k`, damping coefficient `c`.
    Physical,
    /// Natural frequency `wn` (rad/s) and damping ratio `zeta` (unit mass).
    Modal,
}

/// Persistent form + result state for the Vibration Workbench.
pub struct VibrationWorkbenchState {
    /// Whether the system is entered as physical constants or modal pair.
    input: SystemInput,
    /// Mass `m` (kg) — used in [`SystemInput::Physical`].
    mass_kg: f64,
    /// Spring stiffness `k` (N/m) — used in [`SystemInput::Physical`].
    stiffness_n_per_m: f64,
    /// Viscous damping coefficient `c` (N·s/m) — used in
    /// [`SystemInput::Physical`].
    damping_n_s_per_m: f64,
    /// Natural frequency `wn` (rad/s) — used in [`SystemInput::Modal`].
    natural_freq_rad_s: f64,
    /// Damping ratio `zeta` (dimensionless) — used in
    /// [`SystemInput::Modal`].
    damping_ratio: f64,
    /// Harmonic forcing frequency ratio `r = w/wn` (dimensionless).
    frequency_ratio: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D mass-spring-damper solid (serviced
    /// after the panel draws).
    show_3d_request: bool,
}

impl Default for VibrationWorkbenchState {
    fn default() -> Self {
        // m = 1 kg, k = 400 N/m, c = 4 N·s/m:
        //   wn = sqrt(400) = 20 rad/s, c_crit = 2*sqrt(400) = 40,
        //   zeta = 4/40 = 0.1 (lightly damped). Driven at resonance
        //   (r = 1) the magnification factor is 1/(2*zeta) = 5.
        Self {
            input: SystemInput::Physical,
            mass_kg: 1.0,
            stiffness_n_per_m: 400.0,
            damping_n_s_per_m: 4.0,
            natural_freq_rad_s: 20.0,
            damping_ratio: 0.1,
            frequency_ratio: 1.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Vibration Workbench right-side panel. A no-op when the
/// `show_vibration_workbench` toggle is off.
pub fn draw_vibration_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_vibration_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_vibration_workbench",
        "Vibration (SDOF)",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native single-DOF forced-vibration analysis · valenx-vibration",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.vibration;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("System input").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.input, SystemInput::Physical, "Physical (m, k, c)");
                        ui.radio_value(&mut s.input, SystemInput::Modal, "Modal (ωₙ, ζ)");
                    });

                    ui.add_space(4.0);
                    match s.input {
                        SystemInput::Physical => {
                            ui.label(egui::RichText::new("Mass-spring-damper").strong());
                            ui.horizontal(|ui| {
                                ui.label("mass m (kg)");
                                ui.add(egui::DragValue::new(&mut s.mass_kg).speed(0.1));
                            });
                            ui.horizontal(|ui| {
                                ui.label("stiffness k (N/m)");
                                ui.add(egui::DragValue::new(&mut s.stiffness_n_per_m).speed(5.0));
                            });
                            ui.horizontal(|ui| {
                                ui.label("damping c (N·s/m)");
                                ui.add(egui::DragValue::new(&mut s.damping_n_s_per_m).speed(0.2));
                            });
                        }
                        SystemInput::Modal => {
                            ui.label(egui::RichText::new("Modal pair").strong());
                            ui.horizontal(|ui| {
                                ui.label("natural freq ωₙ (rad/s)");
                                ui.add(egui::DragValue::new(&mut s.natural_freq_rad_s).speed(0.5));
                            });
                            ui.horizontal(|ui| {
                                ui.label("damping ratio ζ");
                                ui.add(egui::DragValue::new(&mut s.damping_ratio).speed(0.01));
                            });
                        }
                    }

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Harmonic forcing").strong());
                    ui.horizontal(|ui| {
                        ui.label("frequency ratio r = ω/ωₙ");
                        ui.add(egui::DragValue::new(&mut s.frequency_ratio).speed(0.02));
                    });
                    ui.label(
                        egui::RichText::new("r = 1 is resonance; isolation (T < 1) only for r > √2")
                            .weak()
                            .small(),
                    );

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_vibration(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative mass-spring-damper (a mass block on a coil spring and damper post) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Response").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_vibration_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.vibration` borrow is
    // released here): build the mass-spring-damper 3-D solid and load it.
    if app.vibration.show_3d_request {
        app.vibration.show_3d_request = false;
        load_sdof_3d(app);
    }
}

/// Validate the form, evaluate the system and format the readout.
fn run_vibration(s: &mut VibrationWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build the [`SdofSystem`] from whichever input mode is selected, mapping
/// any domain error to a display string. Shared by the readout and the 3-D
/// gate so both reject the same invalid forms.
fn build_system(s: &VibrationWorkbenchState) -> Result<SdofSystem, String> {
    match s.input {
        SystemInput::Physical => {
            SdofSystem::new(s.mass_kg, s.stiffness_n_per_m, s.damping_n_s_per_m)
                .map_err(|e| e.to_string())
        }
        SystemInput::Modal => {
            SdofSystem::from_modal(s.natural_freq_rad_s, s.damping_ratio).map_err(|e| e.to_string())
        }
    }
}

/// Human-readable name for a [`DampingRegime`].
fn regime_label(regime: DampingRegime) -> &'static str {
    match regime {
        DampingRegime::Undamped => "undamped (ζ = 0)",
        DampingRegime::Underdamped => "underdamped (ζ < 1)",
        DampingRegime::CriticallyDamped => "critically damped (ζ = 1)",
        DampingRegime::Overdamped => "overdamped (ζ > 1)",
    }
}

/// Evaluate the system and format the full readout, mapping any domain
/// error to a display string. Extracted so it is unit-testable.
fn compute(s: &VibrationWorkbenchState) -> Result<String, String> {
    let sys = build_system(s)?;
    let m = magnification_factor(&sys, s.frequency_ratio).map_err(|e| e.to_string())?;
    let t = transmissibility(&sys, s.frequency_ratio).map_err(|e| e.to_string())?;

    let zeta = sys.damping_ratio();
    let wn = sys.natural_freq_rad_s();

    // Resonant peak and logarithmic decrement only exist for some regimes;
    // report "n/a" rather than failing the whole readout.
    let peak = peak_magnification(&sys)
        .map(|p| format!("{p:.4}"))
        .unwrap_or_else(|_| "n/a".to_string());
    let decrement = logarithmic_decrement(&sys)
        .map(|d| format!("{d:.4}"))
        .unwrap_or_else(|_| "n/a".to_string());

    Ok(format!(
        "mass m          : {:.3} kg\n\
         stiffness k     : {:.2} N/m\n\
         damping c       : {:.3} N·s/m\n\n\
         natural freq ωₙ : {:.4} rad/s ({:.4} Hz)\n\
         damping ratio ζ : {:.4}\n\
         regime          : {}\n\n\
         frequency ratio : {:.4}\n\
         magnification M : {:.4}\n\
         transmissibility: {:.4}\n\
         peak magnif.    : {peak}\n\
         log decrement δ : {decrement}",
        sys.mass_kg(),
        sys.stiffness_n_per_m(),
        sys.damping_n_s_per_m(),
        wn,
        sys.natural_freq_hz(),
        zeta,
        regime_label(sys.regime()),
        s.frequency_ratio,
        m,
        t,
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

/// Append a helical coil spring (a swept circular wire) running axially in
/// `z` from `z0` to `z1`, centred on `(cx, cy)`, with helix radius `coil_r`,
/// wire radius `wire_r` and `turns` full turns. The cross-section is laid
/// out in the local tangent frame so the tube never self-intersects at the
/// resolution used here.
#[allow(clippy::too_many_arguments)]
fn push_coil(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    cx: f64,
    cy: f64,
    z0: f64,
    z1: f64,
    coil_r: f64,
    wire_r: f64,
    turns: f64,
) {
    // Centreline samples along the helix.
    let samples = (turns * 24.0).round().max(8.0) as usize;
    let mut path: Vec<Vector3<f64>> = Vec::with_capacity(samples + 1);
    for i in 0..=samples {
        let f = i as f64 / samples as f64;
        let theta = f * turns * std::f64::consts::TAU;
        path.push(Vector3::new(
            cx + coil_r * theta.cos(),
            cy + coil_r * theta.sin(),
            z0 + (z1 - z0) * f,
        ));
    }

    let seg = 8usize;
    let up = Vector3::new(0.0, 0.0, 1.0);
    let alt = Vector3::new(1.0, 0.0, 0.0);
    let ring_base = nodes.len();

    for i in 0..path.len() {
        let raw = if i + 1 < path.len() {
            path[i + 1] - path[i]
        } else {
            path[i] - path[i - 1]
        };
        let tan = if raw.norm() > 1e-12 {
            raw.normalize()
        } else {
            up
        };
        let reference = if tan.dot(&up).abs() < 0.9 { up } else { alt };
        let n1 = tan.cross(&reference).normalize();
        let n2 = tan.cross(&n1).normalize();
        for j in 0..seg {
            let phi = j as f64 / seg as f64 * std::f64::consts::TAU;
            nodes.push(path[i] + (n1 * phi.cos() + n2 * phi.sin()) * wire_r);
        }
    }

    for i in 0..path.len() - 1 {
        let a = ring_base + i * seg;
        let b = ring_base + (i + 1) * seg;
        for j in 0..seg {
            let jn = (j + 1) % seg;
            tris.extend_from_slice(&[a + j, b + j, b + jn, a + j, b + jn, a + jn]);
        }
    }
}

/// Build the representative mass-spring-damper as a triangle [`Mesh`]: a
/// mass block on top, a coil spring and a parallel damper post carrying it
/// down to a base plate. Representative geometry (not to scale; the
/// magnification / transmissibility numbers are the `valenx-vibration`
/// result). `None` for a configuration that does not yield a valid system.
fn sdof_solid_mesh(s: &VibrationWorkbenchState) -> Option<Mesh> {
    build_system(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Base plate.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.05),
        Vector3::new(0.6, 0.6, 0.05),
    );
    // Mass block on top.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 1.3),
        Vector3::new(0.45, 0.45, 0.3),
    );
    // Coil spring carrying the mass (left of centre).
    push_coil(&mut nodes, &mut tris, -0.25, 0.0, 0.1, 1.0, 0.18, 0.03, 5.0);
    // Damper post (a thin cylinder approximated as a slim box, right side).
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.25, 0.0, 0.55),
        Vector3::new(0.08, 0.08, 0.45),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-vibration");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D mass-spring-damper solid and load it into the central
/// viewport.
fn load_sdof_3d(app: &mut ValenxApp) {
    let Some(mesh) = sdof_solid_mesh(&app.vibration) else {
        app.vibration.error =
            Some("system parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<sdof>/valenx-vibration"),
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
        let s = VibrationWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_modal_and_response() {
        let mut s = VibrationWorkbenchState::default();
        run_vibration(&mut s);
        assert!(
            s.error.is_none(),
            "default system should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("natural freq ωₙ"));
        assert!(s.result.contains("magnification M"));
        assert!(s.result.contains("transmissibility"));
        // wn = sqrt(400) = 20 rad/s, zeta = 4/40 = 0.1.
        assert!(s.result.contains("20.0000 rad/s"));
        assert!(s.result.contains("damping ratio ζ : 0.1000"));
        // Driven at r = 1 with zeta = 0.1, so M = 1/(2*zeta) = 5.
        assert!(s.result.contains("magnification M : 5.0000"));
    }

    #[test]
    fn analyze_rejects_zero_mass() {
        let mut s = VibrationWorkbenchState {
            mass_kg: 0.0,
            ..Default::default()
        };
        run_vibration(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn magnification_at_resonance_equals_half_over_zeta() {
        // Ground truth: for a harmonically driven SDOF system at resonance
        // (r = 1), the magnification factor equals 1/(2*zeta). Hand value
        // for zeta = 0.1 is exactly 5.0.
        let sys = SdofSystem::from_modal(20.0, 0.1).expect("valid");
        let m = magnification_factor(&sys, 1.0).expect("ok");
        let hand = 1.0_f64 / (2.0 * 0.1);
        assert!((hand - 5.0).abs() < 1e-12);
        assert!((m - hand).abs() < 1e-9);
    }

    #[test]
    fn sdof_mesh_for_default_is_nonempty_and_in_range() {
        let s = VibrationWorkbenchState::default();
        let mesh = sdof_solid_mesh(&s).expect("default system yields a solid");
        assert!(mesh.nodes.len() > 8, "expected mass + coil + damper + base");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn sdof_mesh_none_for_invalid() {
        let s = VibrationWorkbenchState {
            mass_kg: 0.0,
            ..Default::default()
        };
        assert!(sdof_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_vibration_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_vibration_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_vibration_workbench = true;
        run_vibration(&mut app.vibration);
        draw_workbench(&mut app);
    }
}
