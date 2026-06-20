//! The right-side **Osmosis / Starling Workbench** panel — native
//! semipermeable-membrane water-balance analysis over `valenx-osmosis`.
//!
//! Mirrors the Acid-Base / Heat Transfer workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_osmosis_workbench`,
//! toggled from the View menu. The form picks a mode — the colligative
//! van't Hoff osmotic pressure of a single solution, or Starling
//! transcapillary filtration across a capillary wall — and sets the
//! relevant parameters; "Analyze" reports the osmotic pressure and
//! osmolarity (plus a tonicity verdict against a cell reference) for the
//! van't Hoff mode, or the net driving pressure, net filtration flux and
//! its direction for the Starling mode, and "Show 3-D" loads a
//! representative semipermeable membrane between two compartments into the
//! central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_osmosis::{classify_tonicity, FluxDirection, Solution, StarlingParams, Tonicity};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Which `valenx-osmosis` model the form evaluates.
///
/// Each variant maps to one model family: the colligative van't Hoff
/// osmotic pressure of a single dilute solution, or Starling's
/// transcapillary fluid-exchange equation across a capillary wall.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OsmosisMode {
    /// Van't Hoff colligative osmotic pressure, `Pi = i * c * R * T`.
    VantHoff,
    /// Starling transcapillary filtration,
    /// `Jv = Kf * ((Pc - Pi) - sigma * (pi_c - pi_i))`.
    Starling,
}

impl OsmosisMode {
    /// Short human-readable label for the readout.
    fn label(self) -> &'static str {
        match self {
            OsmosisMode::VantHoff => "van't Hoff osmotic pressure",
            OsmosisMode::Starling => "Starling filtration",
        }
    }
}

/// Persistent form + result state for the Osmosis / Starling Workbench.
pub struct OsmosisWorkbenchState {
    /// Selected osmosis model.
    mode: OsmosisMode,
    /// Solute formula-unit molarity `c` (mol/L); van't Hoff mode.
    concentration_mol_per_l: f64,
    /// Van't Hoff dissociation factor `i` (van't Hoff mode).
    vant_hoff_i: f64,
    /// Temperature (deg C); van't Hoff mode.
    temperature_c: f64,
    /// Cell-interior reference osmolarity (osmol/L) the solution's tonicity
    /// is judged against; van't Hoff mode.
    reference_osmolarity: f64,
    /// Capillary hydrostatic pressure `Pc` (mmHg); Starling mode.
    capillary_hydrostatic: f64,
    /// Interstitial hydrostatic pressure `Pi` (mmHg); Starling mode.
    interstitial_hydrostatic: f64,
    /// Capillary (plasma) oncotic pressure `pi_c` (mmHg); Starling mode.
    capillary_oncotic: f64,
    /// Interstitial oncotic pressure `pi_i` (mmHg); Starling mode.
    interstitial_oncotic: f64,
    /// Staverman reflection coefficient `sigma` (0..=1); Starling mode.
    reflection_sigma: f64,
    /// Filtration coefficient `Kf` (volume / time / mmHg); Starling mode.
    filtration_kf: f64,
    /// Formatted result readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D membrane solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for OsmosisWorkbenchState {
    fn default() -> Self {
        // Van't Hoff defaults: physiological saline (0.9% NaCl ~= 0.154
        // mol/L, i = 2 full dissociation) at body temperature 37 C, judged
        // against a 0.300 osmol/L cell reference. Osmolarity ~= 0.308
        // osmol/L (about isotonic), Pi ~= 7.84 atm.
        // Starling defaults: Guyton arteriolar-end capillary (mmHg) —
        // Pc 35, Pi -2, pi_c 28, pi_i 5, sigma 1, Kf 1 — NDP = +14 mmHg,
        // net outward filtration.
        Self {
            mode: OsmosisMode::VantHoff,
            concentration_mol_per_l: 0.154,
            vant_hoff_i: 2.0,
            temperature_c: 37.0,
            reference_osmolarity: 0.300,
            capillary_hydrostatic: 35.0,
            interstitial_hydrostatic: -2.0,
            capillary_oncotic: 28.0,
            interstitial_oncotic: 5.0,
            reflection_sigma: 1.0,
            filtration_kf: 1.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Osmosis / Starling Workbench right-side panel. A no-op when
/// the `show_osmosis_workbench` toggle is off.
pub fn draw_osmosis_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_osmosis_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_osmosis_workbench",
        "Osmosis / Starling",
        |app, ui| {
            ui.label(
                egui::RichText::new("native semipermeable-membrane water balance · valenx-osmosis")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.osmosis;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Model").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.mode, OsmosisMode::VantHoff, "Van't Hoff");
                        ui.radio_value(&mut s.mode, OsmosisMode::Starling, "Starling");
                    });

                    ui.add_space(4.0);
                    if s.mode == OsmosisMode::VantHoff {
                        ui.label(egui::RichText::new("Solution").strong());
                        ui.horizontal(|ui| {
                            ui.label("c (mol/L)");
                            ui.add(
                                egui::DragValue::new(&mut s.concentration_mol_per_l).speed(0.01),
                            );
                        });
                        ui.horizontal(|ui| {
                            ui.label("van't Hoff i");
                            ui.add(egui::DragValue::new(&mut s.vant_hoff_i).speed(0.05));
                        });
                        ui.horizontal(|ui| {
                            ui.label("temperature (°C)");
                            ui.add(egui::DragValue::new(&mut s.temperature_c).speed(0.5));
                        });
                        ui.horizontal(|ui| {
                            ui.label("cell ref (osmol/L)");
                            ui.add(
                                egui::DragValue::new(&mut s.reference_osmolarity).speed(0.01),
                            );
                        });
                    } else {
                        ui.label(egui::RichText::new("Hydrostatic (mmHg)").strong());
                        ui.horizontal(|ui| {
                            ui.label("capillary Pc");
                            ui.add(
                                egui::DragValue::new(&mut s.capillary_hydrostatic).speed(0.5),
                            );
                        });
                        ui.horizontal(|ui| {
                            ui.label("interstitial Pi");
                            ui.add(
                                egui::DragValue::new(&mut s.interstitial_hydrostatic).speed(0.5),
                            );
                        });

                        ui.add_space(4.0);
                        ui.label(egui::RichText::new("Oncotic (mmHg)").strong());
                        ui.horizontal(|ui| {
                            ui.label("capillary πc");
                            ui.add(egui::DragValue::new(&mut s.capillary_oncotic).speed(0.5));
                        });
                        ui.horizontal(|ui| {
                            ui.label("interstitial πi");
                            ui.add(
                                egui::DragValue::new(&mut s.interstitial_oncotic).speed(0.5),
                            );
                        });

                        ui.add_space(4.0);
                        ui.label(egui::RichText::new("Membrane").strong());
                        ui.horizontal(|ui| {
                            ui.label("reflection σ");
                            ui.add(
                                egui::DragValue::new(&mut s.reflection_sigma)
                                    .speed(0.02)
                                    .range(0.0..=1.0),
                            );
                        });
                        ui.horizontal(|ui| {
                            ui.label("filtration Kf");
                            ui.add(egui::DragValue::new(&mut s.filtration_kf).speed(0.1));
                        });
                    }

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_osmosis(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative semipermeable membrane (a thin slab) between two compartment boxes as a 3-D solid and load it into the central viewport to orbit",
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
        app.show_osmosis_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.osmosis` borrow is
    // released here): build the membrane's 3-D solid and load it.
    if app.osmosis.show_3d_request {
        app.osmosis.show_3d_request = false;
        load_membrane_3d(app);
    }
}

/// Validate the form, evaluate the selected model and format the readout.
fn run_osmosis(s: &mut OsmosisWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// A short human-readable label for a [`FluxDirection`].
fn flux_label(d: FluxDirection) -> &'static str {
    match d {
        FluxDirection::Filtration => "filtration (out of capillary)",
        FluxDirection::Reabsorption => "reabsorption (into capillary)",
        FluxDirection::Equilibrium => "equilibrium (no net flux)",
    }
}

/// A short human-readable label for a [`Tonicity`] verdict.
fn tonicity_label(t: Tonicity) -> &'static str {
    match t {
        Tonicity::Hypotonic => "hypotonic (water enters cell)",
        Tonicity::Isotonic => "isotonic (no net water)",
        Tonicity::Hypertonic => "hypertonic (water leaves cell)",
    }
}

/// Build the validated [`Solution`] for the current van't Hoff form — the
/// quantity both the readout and the 3-D gate need. Extracted so it is
/// unit-testable and shared.
fn solution(s: &OsmosisWorkbenchState) -> Result<Solution, String> {
    Solution::new_celsius(s.concentration_mol_per_l, s.vant_hoff_i, s.temperature_c)
        .map_err(|e| e.to_string())
}

/// Build the validated [`StarlingParams`] for the current Starling form —
/// the quantity both the readout and the 3-D gate need. Extracted so it is
/// unit-testable and shared.
fn starling(s: &OsmosisWorkbenchState) -> Result<StarlingParams, String> {
    StarlingParams::new(
        s.capillary_hydrostatic,
        s.interstitial_hydrostatic,
        s.capillary_oncotic,
        s.interstitial_oncotic,
        s.reflection_sigma,
        s.filtration_kf,
    )
    .map_err(|e| e.to_string())
}

/// Evaluate the selected model and format the full readout, mapping any
/// domain error to a display string. Extracted so it is unit-testable.
fn compute(s: &OsmosisWorkbenchState) -> Result<String, String> {
    match s.mode {
        OsmosisMode::VantHoff => {
            let sol = solution(s)?;
            let osmolarity = sol.osmolarity_osmol_per_l();
            let pi_atm = sol.osmotic_pressure_atm();
            let pi_kpa = sol.osmotic_pressure_pa() / 1000.0;
            let tonicity =
                classify_tonicity(osmolarity, s.reference_osmolarity).map_err(|e| e.to_string())?;
            Ok(format!(
                "model           : {}\n\
                 c / i / T       : {:.4} mol/L / {:.2} / {:.1} °C\n\
                 cell reference  : {:.4} osmol/L\n\n\
                 osmolarity i·c  : {osmolarity:.4} osmol/L\n\
                 osmotic Π       : {pi_atm:.3} atm\n\
                 osmotic Π       : {pi_kpa:.1} kPa\n\
                 tonicity        : {}",
                s.mode.label(),
                s.concentration_mol_per_l,
                s.vant_hoff_i,
                s.temperature_c,
                s.reference_osmolarity,
                tonicity_label(tonicity),
            ))
        }
        OsmosisMode::Starling => {
            let p = starling(s)?;
            let net_hydro = p.net_hydrostatic();
            let net_onco = p.net_oncotic();
            let ndp = p.net_driving_pressure();
            let jv = p.net_filtration();
            let dir = p.direction();
            Ok(format!(
                "model           : {}\n\
                 Pc / Pi         : {:.1} / {:.1} mmHg\n\
                 πc / πi         : {:.1} / {:.1} mmHg\n\
                 σ / Kf          : {:.2} / {:.3}\n\n\
                 net hydrostatic : {net_hydro:.2} mmHg\n\
                 net oncotic     : {net_onco:.2} mmHg\n\
                 net driving NDP : {ndp:.2} mmHg\n\
                 net flux Jv     : {jv:.3}\n\
                 direction       : {}",
                s.mode.label(),
                s.capillary_hydrostatic,
                s.interstitial_hydrostatic,
                s.capillary_oncotic,
                s.interstitial_oncotic,
                s.reflection_sigma,
                s.filtration_kf,
                flux_label(dir),
            ))
        }
    }
}

/// Whether the current form yields a valid model (either family),
/// gating the 3-D build.
fn is_valid(s: &OsmosisWorkbenchState) -> bool {
    match s.mode {
        OsmosisMode::VantHoff => solution(s).is_ok(),
        OsmosisMode::Starling => starling(s).is_ok(),
    }
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

/// Build the osmosis cell as a triangle [`Mesh`] — a thin semipermeable
/// membrane slab (thin in the through-thickness `x` direction) flanked by
/// two compartment boxes (left and right), on a base. Representative
/// geometry (not to scale; the osmosis / Starling numbers are the
/// `valenx-osmosis` result). `None` for an invalid configuration.
fn membrane_solid_mesh(s: &OsmosisWorkbenchState) -> Option<Mesh> {
    // Gate the 3-D build on a successful model evaluation.
    if !is_valid(s) {
        return None;
    }

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Thin semipermeable membrane at the centre (the x = 0 plane).
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.6),
        Vector3::new(0.015, 0.6, 0.5),
    );
    // Left compartment box.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(-0.35, 0.0, 0.6),
        Vector3::new(0.3, 0.55, 0.45),
    );
    // Right compartment box.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.35, 0.0, 0.6),
        Vector3::new(0.3, 0.55, 0.45),
    );
    // Base.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.05),
        Vector3::new(0.75, 0.65, 0.05),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-osmosis");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D membrane solid and load it into the central viewport.
fn load_membrane_3d(app: &mut ValenxApp) {
    let Some(mesh) = membrane_solid_mesh(&app.osmosis) else {
        app.osmosis.error =
            Some("model parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<membrane>/valenx-osmosis"),
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
    use valenx_osmosis::R_L_ATM;

    #[test]
    fn default_state_is_idle() {
        let s = OsmosisWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_pressure_and_tonicity() {
        let mut s = OsmosisWorkbenchState::default();
        run_osmosis(&mut s);
        assert!(
            s.error.is_none(),
            "default solution should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("osmolarity"));
        assert!(s.result.contains("osmotic Π"));
        assert!(s.result.contains("tonicity"));
        // 0.154 mol/L NaCl, i = 2 -> 0.308 osmol/L.
        assert!(s.result.contains("0.3080"));
    }

    #[test]
    fn analyze_starling_reports_flux_and_direction() {
        let mut s = OsmosisWorkbenchState {
            mode: OsmosisMode::Starling,
            ..Default::default()
        };
        run_osmosis(&mut s);
        assert!(s.error.is_none(), "default Starling should analyze");
        assert!(s.result.contains("net driving NDP"));
        assert!(s.result.contains("net flux Jv"));
        // Guyton arteriolar end: NDP = (35-(-2)) - 1*(28-5) = +14 mmHg.
        assert!(s.result.contains("14.00"));
        assert!(s.result.contains("filtration"));
    }

    #[test]
    fn analyze_rejects_below_unity_vant_hoff_factor() {
        let mut s = OsmosisWorkbenchState {
            mode: OsmosisMode::VantHoff,
            vant_hoff_i: 0.5,
            ..Default::default()
        };
        run_osmosis(&mut s);
        assert!(s.error.is_some());
    }

    /// Ground truth: van't Hoff `Pi = i * c * R * T`, hand-computed for a
    /// known molarity. A 1 mol/L ideal (i = 1) solute at 0 C = 273.15 K
    /// gives `Pi = R * T = 0.0820574 * 273.15 = 22.414 atm` — exactly the
    /// molar volume of an ideal gas at STP.
    #[test]
    fn vant_hoff_pressure_matches_hand_computed_value() {
        let s = OsmosisWorkbenchState {
            mode: OsmosisMode::VantHoff,
            concentration_mol_per_l: 1.0,
            vant_hoff_i: 1.0,
            temperature_c: 0.0,
            ..Default::default()
        };
        let sol = solution(&s).expect("1 mol/L ideal solute at 0 C is valid");
        let pi = sol.osmotic_pressure_atm();
        let expected = 1.0 * 1.0 * R_L_ATM * 273.15;
        assert!((pi - expected).abs() < 1e-9, "Pi = {pi} vs {expected}");
        assert!((pi - 22.414).abs() < 1e-3, "expected ~22.414 atm, got {pi}");
    }

    /// Ground truth: Starling net filtration flips sign with the gradient
    /// balance. The Guyton arteriolar end (NDP +14 mmHg) filters; the
    /// venular end with Pc dropped to 15 reabsorbs (NDP -6 mmHg).
    #[test]
    fn starling_direction_flips_between_capillary_ends() {
        let arterial = OsmosisWorkbenchState {
            mode: OsmosisMode::Starling,
            capillary_hydrostatic: 35.0,
            ..Default::default()
        };
        let p = starling(&arterial).expect("valid arterial end");
        assert!((p.net_driving_pressure() - 14.0).abs() < 1e-9);
        assert_eq!(p.direction(), FluxDirection::Filtration);

        let venular = OsmosisWorkbenchState {
            mode: OsmosisMode::Starling,
            capillary_hydrostatic: 15.0,
            ..Default::default()
        };
        let p = starling(&venular).expect("valid venular end");
        assert!((p.net_driving_pressure() - (-6.0)).abs() < 1e-9);
        assert_eq!(p.direction(), FluxDirection::Reabsorption);
    }

    #[test]
    fn membrane_mesh_for_default_is_nonempty_and_in_range() {
        let s = OsmosisWorkbenchState::default();
        let mesh = membrane_solid_mesh(&s).expect("default model yields a solid");
        assert!(
            mesh.nodes.len() > 8,
            "expected membrane + two compartments + base"
        );
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn membrane_mesh_none_for_invalid() {
        let s = OsmosisWorkbenchState {
            mode: OsmosisMode::Starling,
            reflection_sigma: 1.5,
            ..Default::default()
        };
        assert!(membrane_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_osmosis_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_osmosis_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_osmosis_workbench = true;
        run_osmosis(&mut app.osmosis);
        draw_workbench(&mut app);
    }
}
