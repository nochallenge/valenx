//! The right-side **Acid-Base Workbench** panel — native aqueous pH /
//! buffer equilibrium evaluation over `valenx-acidbase`.
//!
//! Mirrors the Heat Transfer / Enzyme Kinetics workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_acidbase_workbench`,
//! toggled from the View menu. The form picks a solution mode — a strong or
//! weak monoprotic acid, a strong or weak monoprotic base, or a
//! Henderson-Hasselbalch buffer — and sets the relevant concentration(s)
//! and dissociation constant; "Analyze" reports the pH (and pOH, plus
//! degree of dissociation / buffer capacity where the model defines them),
//! and "Show 3-D beaker" loads a representative liquid-filled beaker solid
//! into the central viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_acidbase::{
    buffer_capacity, fraction_dissociated, fraction_protonated, henderson_hasselbalch,
    ph_strong_acid, ph_strong_base, ph_weak_acid, ph_weak_base, KW_25C, PKW_25C,
};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Which aqueous acid-base model the form evaluates.
///
/// Each variant maps to one `valenx-acidbase` model: a fully dissociating
/// strong acid / base, the partially dissociating weak acid / base
/// equilibrium, or a weak-acid / conjugate-base Henderson-Hasselbalch
/// buffer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SolutionMode {
    /// Strong monoprotic acid: `[H+] = C`, `pH = -log10(C)`.
    StrongAcid,
    /// Weak monoprotic acid: `[H+] ~= sqrt(Ka*C)`.
    WeakAcid,
    /// Strong monoprotic base: `[OH-] = C`, `pH = pKw + log10(C)`.
    StrongBase,
    /// Weak monoprotic base: `[OH-] ~= sqrt(Kb*C)`, `pH = pKw - pOH`.
    WeakBase,
    /// Weak-acid / conjugate-base buffer: `pH = pKa + log10([A-]/[HA])`.
    Buffer,
}

impl SolutionMode {
    /// Short human-readable label for the readout.
    fn label(self) -> &'static str {
        match self {
            SolutionMode::StrongAcid => "strong acid",
            SolutionMode::WeakAcid => "weak acid",
            SolutionMode::StrongBase => "strong base",
            SolutionMode::WeakBase => "weak base",
            SolutionMode::Buffer => "buffer (Henderson-Hasselbalch)",
        }
    }
}

/// Persistent form + result state for the Acid-Base Workbench.
pub struct AcidBaseWorkbenchState {
    /// Selected aqueous model.
    mode: SolutionMode,
    /// Formal concentration `C` (mol/L); the acid / base concentration for
    /// every mode except the buffer.
    concentration_m: f64,
    /// Dissociation constant `Ka` (weak acid) or `Kb` (weak base); used only
    /// for the weak-acid / weak-base modes.
    k_diss: f64,
    /// Buffer acid `pKa`; used only for the buffer mode.
    buffer_pka: f64,
    /// Buffer weak-acid concentration `[HA]` (mol/L); buffer mode only.
    buffer_acid_m: f64,
    /// Buffer conjugate-base concentration `[A-]` (mol/L); buffer mode only.
    buffer_base_m: f64,
    /// Formatted result readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D beaker solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for AcidBaseWorkbenchState {
    fn default() -> Self {
        // A canonical acetate buffer at equal acid / conjugate-base
        // concentrations (acetic acid, pKa 4.76, 0.1 M / 0.1 M): the log
        // term vanishes, so the Henderson-Hasselbalch pH equals the pKa
        // exactly — pH 4.760. The acid / base fields default to a 0.1 M
        // acetic acid (Ka 1.8e-5) for when the mode is switched.
        Self {
            mode: SolutionMode::Buffer,
            concentration_m: 0.1,
            k_diss: 1.8e-5,
            buffer_pka: 4.76,
            buffer_acid_m: 0.1,
            buffer_base_m: 0.1,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Acid-Base Workbench right-side panel. A no-op when the
/// `show_acidbase_workbench` toggle is off.
pub fn draw_acidbase_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_acidbase_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_acidbase_workbench",
        "Acid-Base",
        |app, ui| {
            ui.label(
                egui::RichText::new("native aqueous pH / buffer equilibria · valenx-acidbase")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.acidbase;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Solution").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.mode, SolutionMode::StrongAcid, "Strong acid");
                        ui.radio_value(&mut s.mode, SolutionMode::WeakAcid, "Weak acid");
                    });
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.mode, SolutionMode::StrongBase, "Strong base");
                        ui.radio_value(&mut s.mode, SolutionMode::WeakBase, "Weak base");
                    });
                    ui.radio_value(&mut s.mode, SolutionMode::Buffer, "Buffer");

                    ui.add_space(4.0);
                    if s.mode == SolutionMode::Buffer {
                        ui.label(egui::RichText::new("Buffer (HA / A⁻)").strong());
                        ui.horizontal(|ui| {
                            ui.label("pKa");
                            ui.add(egui::DragValue::new(&mut s.buffer_pka).speed(0.05));
                        });
                        ui.horizontal(|ui| {
                            ui.label("[HA] acid (mol/L)");
                            ui.add(egui::DragValue::new(&mut s.buffer_acid_m).speed(0.01));
                        });
                        ui.horizontal(|ui| {
                            ui.label("[A⁻] base (mol/L)");
                            ui.add(egui::DragValue::new(&mut s.buffer_base_m).speed(0.01));
                        });
                    } else {
                        ui.label(egui::RichText::new("Concentration").strong());
                        ui.horizontal(|ui| {
                            ui.label("C (mol/L)");
                            ui.add(egui::DragValue::new(&mut s.concentration_m).speed(0.01));
                        });
                        if s.mode == SolutionMode::WeakAcid || s.mode == SolutionMode::WeakBase {
                            let k_label = if s.mode == SolutionMode::WeakAcid {
                                "Ka"
                            } else {
                                "Kb"
                            };
                            ui.horizontal(|ui| {
                                ui.label(k_label);
                                ui.add(egui::DragValue::new(&mut s.k_diss).speed(1e-6));
                            });
                        }
                    }

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_acidbase(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D beaker").strong())
                        .on_hover_text(
                            "Build a representative liquid-filled beaker (glass wall, floor and a liquid level inside) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("pH").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_acidbase_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.acidbase` borrow is
    // released here): build the beaker's 3-D solid and load it.
    if app.acidbase.show_3d_request {
        app.acidbase.show_3d_request = false;
        load_beaker_3d(app);
    }
}

/// Validate the form, evaluate the equilibrium and format the readout.
fn run_acidbase(s: &mut AcidBaseWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// The solution pH for the current form — the single quantity both the 3-D
/// gate and the readout always need. Dispatches to the selected
/// `valenx-acidbase` model and maps any domain error to a display string.
/// Extracted so it is unit-testable and shared.
fn solution_ph(s: &AcidBaseWorkbenchState) -> Result<f64, String> {
    match s.mode {
        SolutionMode::StrongAcid => ph_strong_acid(s.concentration_m).map_err(|e| e.to_string()),
        SolutionMode::WeakAcid => {
            ph_weak_acid(s.k_diss, s.concentration_m).map_err(|e| e.to_string())
        }
        SolutionMode::StrongBase => {
            ph_strong_base(s.concentration_m, KW_25C).map_err(|e| e.to_string())
        }
        SolutionMode::WeakBase => {
            ph_weak_base(s.k_diss, s.concentration_m, KW_25C).map_err(|e| e.to_string())
        }
        SolutionMode::Buffer => {
            henderson_hasselbalch(s.buffer_pka, s.buffer_acid_m, s.buffer_base_m)
                .map_err(|e| e.to_string())
        }
    }
}

/// Evaluate the equilibrium and format the full readout, mapping any domain
/// error to a display string. Extracted so it is unit-testable.
fn compute(s: &AcidBaseWorkbenchState) -> Result<String, String> {
    let ph = solution_ph(s)?;
    // pOH closes via the water relation pH + pOH = pKw at 25 C.
    let poh = PKW_25C - ph;

    let inputs = match s.mode {
        SolutionMode::Buffer => format!(
            "pKa             : {:.3}\n\
             [HA] / [A⁻]     : {:.4} / {:.4} mol/L\n",
            s.buffer_pka, s.buffer_acid_m, s.buffer_base_m,
        ),
        SolutionMode::WeakAcid => format!(
            "concentration C : {:.4} mol/L\n\
             Ka              : {:.3e}\n",
            s.concentration_m, s.k_diss,
        ),
        SolutionMode::WeakBase => format!(
            "concentration C : {:.4} mol/L\n\
             Kb              : {:.3e}\n",
            s.concentration_m, s.k_diss,
        ),
        SolutionMode::StrongAcid | SolutionMode::StrongBase => {
            format!("concentration C : {:.4} mol/L\n", s.concentration_m)
        }
    };

    // Mode-specific extra line: degree of dissociation / protonation for the
    // weak models, Van Slyke buffer capacity for the buffer.
    let extra = match s.mode {
        SolutionMode::WeakAcid => {
            let alpha =
                fraction_dissociated(s.k_diss, s.concentration_m).map_err(|e| e.to_string())?;
            let pct = alpha * 100.0;
            format!("\ndissociated α   : {alpha:.4} ({pct:.2} %)")
        }
        SolutionMode::WeakBase => {
            let alpha =
                fraction_protonated(s.k_diss, s.concentration_m).map_err(|e| e.to_string())?;
            let pct = alpha * 100.0;
            format!("\nprotonated α    : {alpha:.4} ({pct:.2} %)")
        }
        SolutionMode::Buffer => {
            let beta = buffer_capacity(s.buffer_pka, s.buffer_acid_m, s.buffer_base_m)
                .map_err(|e| e.to_string())?;
            format!("\nbuffer capacity : {beta:.4} mol/L per pH")
        }
        SolutionMode::StrongAcid | SolutionMode::StrongBase => String::new(),
    };

    Ok(format!(
        "model           : {}\n\
         {inputs}\n\
         pH              : {ph:.3}\n\
         pOH             : {poh:.3}{extra}",
        s.mode.label(),
    ))
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
            bot + jn,
            top + jn,
        ]);
    }
}

/// Build the beaker as a triangle [`Mesh`] — an upright cylindrical glass
/// wall with a floor disk and an inner liquid column at a representative
/// fill level, on a wider base disk. Representative geometry (not to scale;
/// the pH numbers are the `valenx-acidbase` result). `None` for an invalid
/// configuration.
fn beaker_solid_mesh(s: &AcidBaseWorkbenchState) -> Option<Mesh> {
    // Gate the 3-D build on a successful equilibrium evaluation.
    solution_ph(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Glass wall — the tall upright cylinder.
    push_cyl_z(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.1),
        1.0,
        0.45,
        32,
    );
    // Glass floor disk.
    push_cyl_z(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.08),
        0.04,
        0.45,
        32,
    );
    // Liquid column inside, filling to a representative level.
    push_cyl_z(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.12),
        0.62,
        0.42,
        32,
    );
    // Wider base disk.
    push_cyl_z(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.04),
        0.04,
        0.6,
        32,
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-acidbase");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D beaker solid and load it into the central viewport.
fn load_beaker_3d(app: &mut ValenxApp) {
    let Some(mesh) = beaker_solid_mesh(&app.acidbase) else {
        app.acidbase.error =
            Some("solution parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<beaker>/valenx-acidbase"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"acidbase"}`** product: the canonical
/// beaker-of-solution solid (the panel's "Show 3-D beaker" geometry) paired
/// with the workbench's own buffer / pH headline numbers, at a fixed 3/4
/// camera. Registered in [`crate::products_registry`]; the per-tool builder
/// the registry dispatches to. Pure — driven off
/// [`AcidBaseWorkbenchState::default`].
///
/// The readout rows mirror the panel's `compute()` readout.
pub(crate) fn acidbase_product() -> crate::WorkspaceProduct {
    let s = AcidBaseWorkbenchState::default();
    let mesh = beaker_solid_mesh(&s).expect("default solution ⇒ a 3-D beaker");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<beaker>/valenx-acidbase");
    let readout = compute(&s).expect("default solution ⇒ a valid readout");
    let lines = crate::products_registry::lines_from_readout(&readout);
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Acid-Base".into(),
        lines,
        mesh: Some(loaded),
        vertex_colors: None,
        camera,
        kind2d: None,
        last_export: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_acidbase::poh_weak_base;

    #[test]
    fn default_state_is_idle() {
        let s = AcidBaseWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_ph_and_capacity() {
        let mut s = AcidBaseWorkbenchState::default();
        run_acidbase(&mut s);
        assert!(
            s.error.is_none(),
            "default buffer should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("pH"));
        assert!(s.result.contains("pOH"));
        assert!(s.result.contains("buffer capacity"));
        // Equal acid/base at pKa 4.76 → the log term vanishes → pH = pKa.
        assert!(s.result.contains("4.760"));
    }

    #[test]
    fn analyze_rejects_zero_concentration() {
        let mut s = AcidBaseWorkbenchState {
            mode: SolutionMode::StrongAcid,
            concentration_m: 0.0,
            ..Default::default()
        };
        run_acidbase(&mut s);
        assert!(s.error.is_some());
    }

    /// Ground truth: a strong monoprotic acid dissociates completely, so
    /// `[H+] = C` and `pH = -log10(C)` exactly — 0.1 M HCl is pH 1.
    #[test]
    fn strong_acid_ph_is_neg_log10_of_concentration() {
        let s = AcidBaseWorkbenchState {
            mode: SolutionMode::StrongAcid,
            concentration_m: 0.1,
            ..Default::default()
        };
        let ph = solution_ph(&s).expect("0.1 M strong acid is valid");
        assert!((ph - 1.0).abs() < 1e-12, "pH = {ph}, want 1");
        // And the general relation -log10(C) for another concentration.
        let s2 = AcidBaseWorkbenchState {
            concentration_m: 1.0e-3,
            ..s
        };
        let ph2 = solution_ph(&s2).expect("1 mM strong acid is valid");
        assert!((ph2 - 3.0).abs() < 1e-12, "pH = {ph2}, want 3");
    }

    /// Ground truth: Henderson-Hasselbalch gives `pH = pKa` exactly when the
    /// weak acid and its conjugate base are at equal concentration (the log
    /// term is zero), independent of that shared concentration.
    #[test]
    fn buffer_ph_equals_pka_at_equal_concentrations() {
        let s = AcidBaseWorkbenchState {
            mode: SolutionMode::Buffer,
            buffer_pka: 4.76,
            buffer_acid_m: 0.1,
            buffer_base_m: 0.1,
            ..Default::default()
        };
        let ph = solution_ph(&s).expect("valid buffer");
        assert!((ph - 4.76).abs() < 1e-12, "pH = {ph}, want pKa 4.76");
    }

    /// The strong-base mode is consistent with the strong-acid mode through
    /// the water relation: at equal concentration the two pH values sum to
    /// `pKw = 14` at 25 C.
    #[test]
    fn strong_acid_and_base_sum_to_pkw() {
        let acid = AcidBaseWorkbenchState {
            mode: SolutionMode::StrongAcid,
            concentration_m: 0.02,
            ..Default::default()
        };
        let base = AcidBaseWorkbenchState {
            mode: SolutionMode::StrongBase,
            concentration_m: 0.02,
            ..Default::default()
        };
        let pa = solution_ph(&acid).expect("valid");
        let pb = solution_ph(&base).expect("valid");
        assert!((pa + pb - 14.0).abs() < 1e-12, "{pa} + {pb} != 14");
    }

    #[test]
    fn weak_base_poh_matches_poh_helper() {
        // The weak-base pH equals pKw minus the model's own pOH (25 C).
        let s = AcidBaseWorkbenchState {
            mode: SolutionMode::WeakBase,
            concentration_m: 0.1,
            k_diss: 1.8e-5,
            ..Default::default()
        };
        let ph = solution_ph(&s).expect("valid weak base");
        let poh = poh_weak_base(1.8e-5, 0.1).expect("valid");
        assert!((ph - (PKW_25C - poh)).abs() < 1e-9, "pH = {ph}");
    }

    #[test]
    fn beaker_mesh_for_default_is_nonempty_and_in_range() {
        let s = AcidBaseWorkbenchState::default();
        let mesh = beaker_solid_mesh(&s).expect("default solution yields a solid");
        assert!(
            mesh.nodes.len() > 8,
            "expected wall + floor + liquid + base"
        );
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn beaker_mesh_none_for_invalid() {
        let s = AcidBaseWorkbenchState {
            mode: SolutionMode::StrongAcid,
            concentration_m: 0.0,
            ..Default::default()
        };
        assert!(beaker_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_acidbase_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_acidbase_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_acidbase_workbench = true;
        run_acidbase(&mut app.acidbase);
        draw_workbench(&mut app);
    }
}
