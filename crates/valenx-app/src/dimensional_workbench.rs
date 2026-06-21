//! The right-side **Dimensionless Numbers Workbench** panel — native
//! similitude analysis over `valenx-dimensional`.
//!
//! Mirrors the Heat Transfer / Antenna workbenches: a resizable
//! [`egui::SidePanel`] gated on
//! `crate::ValenxApp::show_dimensional_workbench`, toggled from the View
//! menu. A selector picks one of the canonical dimensionless groups
//! (Biot, Reynolds, Nusselt, Prandtl, Mach, Froude, Peclet); only the
//! inputs that group needs are shown. "Analyze" builds the group with the
//! validating `valenx-dimensional` constructor and reports its value plus
//! the relevant regime / validity classifier, and "Show 3-D" loads a
//! representative labelled box into the central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_dimensional::biot::{Biot, LumpedCapacitance};
use valenx_dimensional::froude::{ChannelRegime, Froude};
use valenx_dimensional::mach::{Mach, SpeedRegime};
use valenx_dimensional::nusselt::Nusselt;
use valenx_dimensional::peclet::Peclet;
use valenx_dimensional::prandtl::Prandtl;
use valenx_dimensional::reynolds::{PipeRegime, Reynolds};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Which dimensionless group the workbench evaluates. Drives both the form
/// (only the inputs the selected group needs are shown) and `compute`.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum DimensionlessNumber {
    /// Biot number `Bi = h L / k` (solid conductivity) with the
    /// lumped-capacitance validity test.
    Biot,
    /// Reynolds number `Re = rho v L / mu` with a pipe-flow regime.
    Reynolds,
    /// Nusselt number `Nu = h L / k` (fluid conductivity).
    Nusselt,
    /// Prandtl number `Pr = cp mu / k`.
    Prandtl,
    /// Mach number `Ma = v / c` with a flow-speed regime.
    Mach,
    /// Froude number `Fr = v / sqrt(g L)` with an open-channel regime.
    Froude,
    /// Peclet number `Pe = rho cp v L / k` (and `= Re * Pr`).
    Peclet,
}

/// Persistent form + result state for the Dimensionless Numbers Workbench.
pub struct DimensionalWorkbenchState {
    /// Which group to evaluate.
    selected: DimensionlessNumber,
    /// Characteristic length `L` (m). Used by Biot / Reynolds / Nusselt /
    /// Froude / Peclet.
    length_m: f64,
    /// Convective heat-transfer coefficient `h` (W/m^2·K). Used by Biot
    /// and Nusselt.
    h_coeff: f64,
    /// Thermal conductivity `k` (W/m·K). Used by Biot / Nusselt / Prandtl
    /// / Peclet (solid for Biot, fluid for the others).
    conductivity_w_per_mk: f64,
    /// Fluid density `rho` (kg/m^3). Used by Reynolds and Peclet.
    density_kg_m3: f64,
    /// Characteristic velocity / speed `v` (m/s). Used by Reynolds / Mach
    /// / Froude / Peclet.
    velocity_m_s: f64,
    /// Dynamic viscosity `mu` (Pa·s). Used by Reynolds and Prandtl.
    dynamic_viscosity_pas: f64,
    /// Specific heat at constant pressure `cp` (J/kg·K). Used by Prandtl
    /// and Peclet.
    specific_heat_j_per_kgk: f64,
    /// Local speed of sound `c` (m/s). Used by Mach.
    speed_of_sound_m_s: f64,
    /// Gravitational acceleration `g` (m/s^2). Used by Froude.
    gravity_m_s2: f64,
    /// Formatted readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the representative 3-D box (serviced
    /// after the panel draws).
    show_3d_request: bool,
}

impl Default for DimensionalWorkbenchState {
    fn default() -> Self {
        // Realistic defaults that all yield a valid group:
        //   Biot (default): h=20, L=0.05, k=50 (steel) -> Bi = 0.02,
        //   well below 0.1, so lumped capacitance is valid.
        //   Reynolds: water in a 50 mm pipe at 2 m/s -> Re ~ 1e5 turbulent.
        //   Mach: 250 m/s at c=340 -> Ma ~ 0.735 subsonic.
        //   Froude: 2 m/s over 1 m depth -> Fr ~ 0.639 subcritical.
        Self {
            selected: DimensionlessNumber::Biot,
            length_m: 0.05,
            h_coeff: 20.0,
            conductivity_w_per_mk: 50.0,
            density_kg_m3: 998.0,
            velocity_m_s: 2.0,
            dynamic_viscosity_pas: 1.0e-3,
            specific_heat_j_per_kgk: 4182.0,
            speed_of_sound_m_s: 340.0,
            gravity_m_s2: 9.81,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Dimensionless Numbers Workbench right-side panel. A no-op when
/// the `show_dimensional_workbench` toggle is off.
pub fn draw_dimensional_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_dimensional_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_dimensional_workbench",
        "Dimensionless Numbers",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native similitude groups + regime classifiers · valenx-dimensional",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.dimensional;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Group").strong());
                    ui.horizontal_wrapped(|ui| {
                        ui.radio_value(&mut s.selected, DimensionlessNumber::Biot, "Biot");
                        ui.radio_value(&mut s.selected, DimensionlessNumber::Reynolds, "Reynolds");
                        ui.radio_value(&mut s.selected, DimensionlessNumber::Nusselt, "Nusselt");
                        ui.radio_value(&mut s.selected, DimensionlessNumber::Prandtl, "Prandtl");
                        ui.radio_value(&mut s.selected, DimensionlessNumber::Mach, "Mach");
                        ui.radio_value(&mut s.selected, DimensionlessNumber::Froude, "Froude");
                        ui.radio_value(&mut s.selected, DimensionlessNumber::Peclet, "Peclet");
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Inputs").strong());

                    // Only the inputs the selected group consumes are shown.
                    match s.selected {
                        DimensionlessNumber::Biot => {
                            drag(ui, "h (W/m²K)", &mut s.h_coeff, 0.5);
                            drag(ui, "length L (m)", &mut s.length_m, 0.005);
                            drag(ui, "k solid (W/m·K)", &mut s.conductivity_w_per_mk, 0.5);
                        }
                        DimensionlessNumber::Reynolds => {
                            drag(ui, "density ρ (kg/m³)", &mut s.density_kg_m3, 1.0);
                            drag(ui, "velocity v (m/s)", &mut s.velocity_m_s, 0.1);
                            drag(ui, "length L (m)", &mut s.length_m, 0.005);
                            drag(ui, "viscosity μ (Pa·s)", &mut s.dynamic_viscosity_pas, 1.0e-4);
                        }
                        DimensionlessNumber::Nusselt => {
                            drag(ui, "h (W/m²K)", &mut s.h_coeff, 0.5);
                            drag(ui, "length L (m)", &mut s.length_m, 0.005);
                            drag(ui, "k fluid (W/m·K)", &mut s.conductivity_w_per_mk, 0.01);
                        }
                        DimensionlessNumber::Prandtl => {
                            drag(ui, "cp (J/kg·K)", &mut s.specific_heat_j_per_kgk, 1.0);
                            drag(ui, "viscosity μ (Pa·s)", &mut s.dynamic_viscosity_pas, 1.0e-4);
                            drag(ui, "k fluid (W/m·K)", &mut s.conductivity_w_per_mk, 0.01);
                        }
                        DimensionlessNumber::Mach => {
                            drag(ui, "speed v (m/s)", &mut s.velocity_m_s, 1.0);
                            drag(ui, "speed of sound c (m/s)", &mut s.speed_of_sound_m_s, 1.0);
                        }
                        DimensionlessNumber::Froude => {
                            drag(ui, "velocity v (m/s)", &mut s.velocity_m_s, 0.1);
                            drag(ui, "gravity g (m/s²)", &mut s.gravity_m_s2, 0.01);
                            drag(ui, "length L (m)", &mut s.length_m, 0.005);
                        }
                        DimensionlessNumber::Peclet => {
                            drag(ui, "density ρ (kg/m³)", &mut s.density_kg_m3, 1.0);
                            drag(ui, "cp (J/kg·K)", &mut s.specific_heat_j_per_kgk, 1.0);
                            drag(ui, "velocity v (m/s)", &mut s.velocity_m_s, 0.1);
                            drag(ui, "length L (m)", &mut s.length_m, 0.005);
                            drag(ui, "k fluid (W/m·K)", &mut s.conductivity_w_per_mk, 0.01);
                        }
                    }

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_dimensional(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative labelled box for the selected group and load it into the central viewport to orbit (the number itself is the valenx-dimensional result, not the geometry)",
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
        app.show_dimensional_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.dimensional` borrow is
    // released here): build the representative box and load it.
    if app.dimensional.show_3d_request {
        app.dimensional.show_3d_request = false;
        load_regime_3d(app);
    }
}

/// One labelled [`egui::DragValue`] row. A closure capturing `&mut` over
/// each field would borrow `s` for the whole form, so this is a free
/// helper taking the field by reference instead.
fn drag(ui: &mut egui::Ui, label: &str, value: &mut f64, speed: f64) {
    ui.horizontal(|ui| {
        ui.label(label);
        ui.add(egui::DragValue::new(value).speed(speed));
    });
}

/// Validate the form, evaluate the selected group and format the readout.
fn run_dimensional(s: &mut DimensionalWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build the selected dimensionless group with its validating
/// `valenx-dimensional` constructor and format the readout, mapping any
/// domain error to a display string. Extracted so it is unit-testable.
fn compute(s: &DimensionalWorkbenchState) -> Result<String, String> {
    match s.selected {
        DimensionlessNumber::Biot => {
            let bi = Biot::new(s.h_coeff, s.length_m, s.conductivity_w_per_mk)
                .map_err(|e| e.to_string())?;
            let regime = match bi.lumped_capacitance() {
                LumpedCapacitance::Valid => "lumped capacitance VALID (Bi < 0.1)",
                LumpedCapacitance::Invalid => "lumped capacitance INVALID (Bi ≥ 0.1)",
            };
            Ok(format!(
                "group           : Biot  Bi = h L / k\n\
                 h / L / k       : {:.3} / {:.4} / {:.3}\n\n\
                 Bi              : {:.4}\n\
                 verdict         : {regime}",
                s.h_coeff,
                s.length_m,
                s.conductivity_w_per_mk,
                bi.value(),
            ))
        }
        DimensionlessNumber::Reynolds => {
            let re = Reynolds::new(
                s.density_kg_m3,
                s.velocity_m_s,
                s.length_m,
                s.dynamic_viscosity_pas,
            )
            .map_err(|e| e.to_string())?;
            let regime = match re.pipe_regime() {
                PipeRegime::Laminar => "laminar (Re < ~2300)",
                PipeRegime::Transitional => "transitional (~2300 ≤ Re < ~4000)",
                PipeRegime::Turbulent => "turbulent (Re ≥ ~4000)",
            };
            Ok(format!(
                "group           : Reynolds  Re = ρ v L / μ\n\
                 ρ / v / L / μ   : {:.1} / {:.3} / {:.4} / {:.2e}\n\n\
                 Re              : {:.1}\n\
                 pipe regime     : {regime}",
                s.density_kg_m3,
                s.velocity_m_s,
                s.length_m,
                s.dynamic_viscosity_pas,
                re.value(),
            ))
        }
        DimensionlessNumber::Nusselt => {
            let nu = Nusselt::new(s.h_coeff, s.length_m, s.conductivity_w_per_mk)
                .map_err(|e| e.to_string())?;
            Ok(format!(
                "group           : Nusselt  Nu = h L / k\n\
                 h / L / k       : {:.3} / {:.4} / {:.3}\n\n\
                 Nu              : {:.4}",
                s.h_coeff,
                s.length_m,
                s.conductivity_w_per_mk,
                nu.value(),
            ))
        }
        DimensionlessNumber::Prandtl => {
            let pr = Prandtl::new(
                s.specific_heat_j_per_kgk,
                s.dynamic_viscosity_pas,
                s.conductivity_w_per_mk,
            )
            .map_err(|e| e.to_string())?;
            Ok(format!(
                "group           : Prandtl  Pr = cp μ / k\n\
                 cp / μ / k      : {:.1} / {:.2e} / {:.3}\n\n\
                 Pr              : {:.4}",
                s.specific_heat_j_per_kgk,
                s.dynamic_viscosity_pas,
                s.conductivity_w_per_mk,
                pr.value(),
            ))
        }
        DimensionlessNumber::Mach => {
            let ma = Mach::new(s.velocity_m_s, s.speed_of_sound_m_s).map_err(|e| e.to_string())?;
            let regime = match ma.speed_regime() {
                SpeedRegime::Subsonic => "subsonic (Ma < 1)",
                SpeedRegime::Sonic => "sonic (Ma = 1)",
                SpeedRegime::Supersonic => "supersonic (Ma > 1)",
            };
            Ok(format!(
                "group           : Mach  Ma = v / c\n\
                 v / c           : {:.3} / {:.3}\n\n\
                 Ma              : {:.4}\n\
                 speed regime    : {regime}",
                s.velocity_m_s,
                s.speed_of_sound_m_s,
                ma.value(),
            ))
        }
        DimensionlessNumber::Froude => {
            let fr = Froude::new(s.velocity_m_s, s.gravity_m_s2, s.length_m)
                .map_err(|e| e.to_string())?;
            let regime = match fr.channel_regime() {
                ChannelRegime::Subcritical => "subcritical / tranquil (Fr < 1)",
                ChannelRegime::Critical => "critical (Fr = 1)",
                ChannelRegime::Supercritical => "supercritical / rapid (Fr > 1)",
            };
            Ok(format!(
                "group           : Froude  Fr = v / sqrt(g L)\n\
                 v / g / L       : {:.3} / {:.3} / {:.4}\n\n\
                 Fr              : {:.4}\n\
                 channel regime  : {regime}",
                s.velocity_m_s,
                s.gravity_m_s2,
                s.length_m,
                fr.value(),
            ))
        }
        DimensionlessNumber::Peclet => {
            let pe = Peclet::new(
                s.density_kg_m3,
                s.specific_heat_j_per_kgk,
                s.velocity_m_s,
                s.length_m,
                s.conductivity_w_per_mk,
            )
            .map_err(|e| e.to_string())?;
            Ok(format!(
                "group           : Peclet  Pe = ρ cp v L / k\n\
                 ρ / cp / v      : {:.1} / {:.1} / {:.3}\n\
                 L / k           : {:.4} / {:.3}\n\n\
                 Pe              : {:.1}",
                s.density_kg_m3,
                s.specific_heat_j_per_kgk,
                s.velocity_m_s,
                s.length_m,
                s.conductivity_w_per_mk,
                pe.value(),
            ))
        }
    }
}

/// Build a representative labelled box (centre `c`, half-extents `h`) as a
/// triangle [`Mesh`]. The geometry is a stand-in for the selected group;
/// the number itself is the `valenx-dimensional` result, not the box.
fn box_mesh(c: Vector3<f64>, h: Vector3<f64>) -> Mesh {
    let mut nodes: Vec<Vector3<f64>> = Vec::with_capacity(8);
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
    let mut tris: Vec<usize> = Vec::new();
    for f in faces {
        tris.extend_from_slice(&[f[0], f[1], f[2], f[0], f[2], f[3]]);
    }

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-dimensional");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    mesh
}

/// Build the representative box for the current configuration. `None` for
/// an invalid configuration (the same validation `compute` applies), so an
/// out-of-domain input does not load stray geometry.
fn regime_box_mesh(s: &DimensionalWorkbenchState) -> Option<Mesh> {
    compute(s).ok()?;
    Some(box_mesh(
        Vector3::new(0.0, 0.0, 0.5),
        Vector3::new(0.5, 0.5, 0.5),
    ))
}

/// Build the representative 3-D box and load it into the central viewport.
fn load_regime_3d(app: &mut ValenxApp) {
    let Some(mesh) = regime_box_mesh(&app.dimensional) else {
        app.dimensional.error = Some("inputs are invalid — cannot build the 3-D box".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<dimensional>/valenx-dimensional"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// Agent-bridge product: the canonical dimensional-analysis workbench as a 3-D
/// solid plus its `compute()` readout rows (see [`crate::products_registry`]).
pub(crate) fn dimensional_product() -> crate::WorkspaceProduct {
    let s = DimensionalWorkbenchState::default();
    let mesh = regime_box_mesh(&s).expect("canonical dimensional ⇒ regime box solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<dimensional>/valenx-regime");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical dimensional ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Dimensional analysis (flow regime)".into(),
        lines,
        mesh: Some(loaded),
        vertex_colors: None,
        camera,
        kind2d: None,
        last_export: None,
        image: None,
        image_texture: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_is_idle() {
        let s = DimensionalWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
        assert_eq!(s.selected, DimensionlessNumber::Biot);
    }

    #[test]
    fn analyze_default_biot_reports_value_and_lumped_valid() {
        // Ground truth: Bi = h L / k = 20 * 0.05 / 50 = 0.02, which is
        // below the 0.1 lumped-capacitance limit, so the verdict is VALID.
        let mut s = DimensionalWorkbenchState::default();
        run_dimensional(&mut s);
        assert!(
            s.error.is_none(),
            "default Biot should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("Biot"));
        assert!(s.result.contains("0.0200"));
        assert!(s.result.contains("lumped capacitance VALID"));
    }

    #[test]
    fn high_biot_is_lumped_invalid() {
        // h=500, L=0.1, k=0.2 -> Bi = 250, far above 0.1.
        let mut s = DimensionalWorkbenchState {
            h_coeff: 500.0,
            length_m: 0.1,
            conductivity_w_per_mk: 0.2,
            ..DimensionalWorkbenchState::default()
        };
        run_dimensional(&mut s);
        assert!(s.error.is_none());
        assert!(s.result.contains("lumped capacitance INVALID"));
    }

    #[test]
    fn analyze_reynolds_reports_turbulent_pipe_regime() {
        // Water 50 mm pipe at 2 m/s: Re = 998*2*0.05/1e-3 ~ 99800, turbulent.
        let mut s = DimensionalWorkbenchState {
            selected: DimensionlessNumber::Reynolds,
            ..DimensionalWorkbenchState::default()
        };
        run_dimensional(&mut s);
        assert!(s.error.is_none(), "{:?}", s.error);
        assert!(s.result.contains("Reynolds"));
        assert!(s.result.contains("turbulent"));
    }

    #[test]
    fn analyze_froude_default_is_subcritical() {
        // v=2, g=9.81, L=0.05 -> Fr = 2/sqrt(0.4905) = 2.856 -> use L=1.
        let mut s = DimensionalWorkbenchState {
            selected: DimensionlessNumber::Froude,
            length_m: 1.0,
            ..DimensionalWorkbenchState::default()
        };
        run_dimensional(&mut s);
        assert!(s.error.is_none(), "{:?}", s.error);
        assert!(s.result.contains("Froude"));
        assert!(s.result.contains("subcritical"));
    }

    #[test]
    fn analyze_rejects_non_positive_conductivity() {
        // Biot with k = 0 is out of domain (zero denominator).
        let mut s = DimensionalWorkbenchState {
            conductivity_w_per_mk: 0.0,
            ..DimensionalWorkbenchState::default()
        };
        run_dimensional(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn box_mesh_for_default_is_nonempty_and_in_range() {
        let s = DimensionalWorkbenchState::default();
        let mesh = regime_box_mesh(&s).expect("default config yields a box");
        assert_eq!(mesh.nodes.len(), 8);
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn box_mesh_none_for_invalid() {
        let s = DimensionalWorkbenchState {
            conductivity_w_per_mk: 0.0,
            ..DimensionalWorkbenchState::default()
        };
        assert!(regime_box_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_dimensional_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_dimensional_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_dimensional_workbench = true;
        run_dimensional(&mut app.dimensional);
        draw_workbench(&mut app);
    }
}
