//! The right-side **Combustion Workbench** panel — native hydrocarbon
//! `CxHy` air-fuel stoichiometry over `valenx-combustion`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_combustion_workbench`,
//! toggled from the View menu. The form picks a pure hydrocarbon fuel and
//! an equivalence ratio `phi` (lean / stoichiometric); "Analyze" reports
//! the stoichiometric air-fuel ratio (mass and molar), the percent excess
//! air, the per-mole product balance (CO2 / H2O / N2 / excess O2) and a
//! mean-cp adiabatic flame temperature, and "Show 3-D" loads a
//! representative combustor can (a cylinder) into the central viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_combustion::flame::{
    adiabatic_flame_temperature, CP_MOLAR_DEFAULT, LHV_METHANE, LHV_OCTANE, LHV_PROPANE, T_REF_K,
};
use valenx_combustion::fuel::Fuel;
use valenx_combustion::stoich::{
    afr_stoich_mass, afr_stoich_molar, percent_excess_air, product_moles,
};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Pure-hydrocarbon fuel choice for the form. Each maps to a
/// [`Fuel`] and its lower heating value.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum FuelChoice {
    /// Methane, CH4.
    Methane,
    /// Propane, C3H8.
    Propane,
    /// iso-Octane, C8H18 (gasoline surrogate).
    Octane,
}

impl FuelChoice {
    /// The `valenx-combustion` [`Fuel`] for this choice.
    fn fuel(self) -> Fuel {
        match self {
            FuelChoice::Methane => Fuel::methane(),
            FuelChoice::Propane => Fuel::propane(),
            FuelChoice::Octane => Fuel::octane(),
        }
    }

    /// Lower heating value of this fuel, J/kg.
    fn lhv(self) -> f64 {
        match self {
            FuelChoice::Methane => LHV_METHANE,
            FuelChoice::Propane => LHV_PROPANE,
            FuelChoice::Octane => LHV_OCTANE,
        }
    }

    /// Short formula label for the readout.
    fn formula(self) -> &'static str {
        match self {
            FuelChoice::Methane => "CH4",
            FuelChoice::Propane => "C3H8",
            FuelChoice::Octane => "C8H18",
        }
    }
}

/// Persistent form + result state for the Combustion Workbench.
pub struct CombustionWorkbenchState {
    /// Hydrocarbon fuel being burned.
    fuel: FuelChoice,
    /// Equivalence ratio `phi` (lean / stoichiometric, `0 < phi <= 1`).
    phi: f64,
    /// Mean molar heat capacity of the product gas, J/(mol·K).
    cp_molar: f64,
    /// Reactant inlet temperature, K.
    t_in_k: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D combustor solid (serviced after
    /// the panel draws).
    show_3d_request: bool,
}

impl Default for CombustionWorkbenchState {
    fn default() -> Self {
        // Methane (CH4) burned at stoichiometric (phi = 1) in air, with
        // the round mean-cp = 40 J/mol·K and a 298.15 K inlet: AFR ~ 17.12
        // by mass, 0% excess air, and a mean-cp adiabatic flame
        // temperature of ~2204 K.
        Self {
            fuel: FuelChoice::Methane,
            phi: 1.0,
            cp_molar: CP_MOLAR_DEFAULT,
            t_in_k: T_REF_K,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Combustion Workbench right-side panel. A no-op when the
/// `show_combustion_workbench` toggle is off.
pub fn draw_combustion_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_combustion_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_combustion_workbench",
        "Combustion",
        |app, ui| {
            ui.label(
                egui::RichText::new("native CxHy air-fuel stoichiometry · valenx-combustion")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.combustion;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Fuel").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.fuel, FuelChoice::Methane, "CH4");
                        ui.radio_value(&mut s.fuel, FuelChoice::Propane, "C3H8");
                        ui.radio_value(&mut s.fuel, FuelChoice::Octane, "C8H18");
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Mixture").strong());
                    ui.horizontal(|ui| {
                        ui.label("equivalence ratio φ");
                        ui.add(
                            egui::DragValue::new(&mut s.phi)
                                .speed(0.01)
                                .range(0.05..=1.0),
                        );
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Flame energy balance").strong());
                    ui.horizontal(|ui| {
                        ui.label("mean cp (J/mol·K)");
                        ui.add(egui::DragValue::new(&mut s.cp_molar).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("inlet T (K)");
                        ui.add(egui::DragValue::new(&mut s.t_in_k).speed(1.0));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_combustion(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative combustor can (a cylindrical burner) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Stoichiometry").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_combustion_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.combustion` borrow is
    // released here): build the combustor's 3-D solid and load it.
    if app.combustion.show_3d_request {
        app.combustion.show_3d_request = false;
        load_combustor_3d(app);
    }
}

/// Validate the form, evaluate the combustion and format the readout.
fn run_combustion(s: &mut CombustionWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Evaluate the air-fuel stoichiometry and flame balance, formatting the
/// full readout and mapping any domain error to a display string.
/// Extracted so it is unit-testable.
fn compute(s: &CombustionWorkbenchState) -> Result<String, String> {
    let fuel = s.fuel.fuel();
    let phi = s.phi;
    let afr_mass = afr_stoich_mass(&fuel);
    let afr_molar = afr_stoich_molar(&fuel);
    let excess = percent_excess_air(phi).map_err(|e| e.to_string())?;
    let products = product_moles(&fuel, phi).map_err(|e| e.to_string())?;
    let t_ad = adiabatic_flame_temperature(&fuel, phi, s.fuel.lhv(), s.cp_molar, s.t_in_k)
        .map_err(|e| e.to_string())?;

    Ok(format!(
        "fuel            : {formula} ({mm:.3} g/mol)\n\
         AFR stoich mass : {afr_mass:.2} kg air/kg fuel\n\
         AFR stoich molar: {afr_molar:.2} mol air/mol fuel\n\
         equivalence φ   : {phi:.3}\n\
         excess air      : {excess:.1} %\n\n\
         products (mol/mol fuel)\n\
         CO2             : {co2:.3}\n\
         H2O             : {h2o:.3}\n\
         N2              : {n2:.3}\n\
         excess O2       : {o2:.3}\n\
         total moles     : {total:.3}\n\n\
         adiabatic flame T: {t_ad:.2} K",
        formula = s.fuel.formula(),
        mm = fuel.molar_mass(),
        co2 = products.co2,
        h2o = products.h2o,
        n2 = products.n2,
        o2 = products.excess_o2,
        total = products.total(),
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
            top + jn,
            top + j,
            bot + j,
            bot + jn,
            top + jn,
        ]);
    }
}

/// Build the combustor as a triangle [`Mesh`] — a cylindrical burner can
/// (the main combustion chamber) with a narrower nozzle throat stub on top.
/// Representative geometry (not to scale; the stoichiometry numbers are the
/// `valenx-combustion` result). `None` for an invalid configuration.
fn combustor_solid_mesh(s: &CombustionWorkbenchState) -> Option<Mesh> {
    // Gate the 3-D build on a configuration the analysis itself accepts.
    product_moles(&s.fuel.fuel(), s.phi).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Main combustion can.
    push_cyl_z(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.0),
        1.0,
        0.35,
        48,
    );
    // Narrower nozzle / throat stub on top.
    push_cyl_z(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 1.0),
        0.3,
        0.18,
        48,
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-combustion");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D combustor solid and load it into the central viewport.
fn load_combustor_3d(app: &mut ValenxApp) {
    let Some(mesh) = combustor_solid_mesh(&app.combustion) else {
        app.combustion.error =
            Some("combustion parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<combustor>/valenx-combustion"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"combustion"}`** product: the canonical
/// combustor can built as a 3-D solid, paired with the workbench's own
/// `compute()` readout rows, at a fixed 3/4 camera. Registered in
/// [`crate::products_registry`]; the per-tool builder the registry dispatches
/// to. Pure — driven off [`CombustionWorkbenchState::default`].
pub(crate) fn combustion_product() -> crate::WorkspaceProduct {
    let s = CombustionWorkbenchState::default();
    let mesh = combustor_solid_mesh(&s).expect("canonical combustion ⇒ combustor solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<combustion>/valenx-combustor");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical combustion ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Combustion (air-fuel stoichiometry)".into(),
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

    #[test]
    fn default_state_is_idle() {
        let s = CombustionWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_afr_and_flame_temperature() {
        let mut s = CombustionWorkbenchState::default();
        run_combustion(&mut s);
        assert!(
            s.error.is_none(),
            "default methane mixture should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("AFR stoich mass"));
        assert!(s.result.contains("adiabatic flame T"));
        // Stoichiometric methane in air: AFR ~ 17.12 by mass, 0% excess.
        assert!(s.result.contains("17.12"));
        assert!(s.result.contains("0.0 %"));
        // Mean-cp adiabatic flame temperature for the default inputs.
        assert!(s.result.contains("2204.40"));
    }

    #[test]
    fn analyze_rejects_rich_mixture() {
        // phi > 1 (rich) is outside the closed-form complete-combustion
        // product balance, so compute() must surface the domain error.
        let mut s = CombustionWorkbenchState {
            phi: 1.5,
            ..Default::default()
        };
        run_combustion(&mut s);
        assert!(s.error.is_some());
        assert!(s.result.is_empty());
    }

    #[test]
    fn ground_truth_methane_stoich_afr_mass() {
        // Hand-computed stoichiometric AFR for CH4 in air (O2 + 3.76 N2):
        //   a = x + y/4 = 1 + 1 = 2
        //   AFR = a*(M_O2 + 3.76*M_N2)/M_fuel
        //       = 2*(31.998 + 3.76*28.013)/16.043 = 17.1199...
        let afr: f64 = afr_stoich_mass(&Fuel::methane());
        let expected: f64 = 2.0 * (31.998 + 3.76 * 28.013) / (12.011 + 4.0 * 1.008);
        assert!((afr - expected).abs() < 1e-9);
        assert!((afr - 17.1199).abs() < 1e-3);
    }

    #[test]
    fn lean_mixture_has_positive_excess_air_and_o2() {
        // phi = 0.5 => 100% excess air, and unburned O2 in the products
        // equal to the stoichiometric O2 demand (a_supplied - a = 2a - a).
        let s = CombustionWorkbenchState {
            phi: 0.5,
            ..Default::default()
        };
        let excess: f64 = percent_excess_air(s.phi).unwrap();
        assert!((excess - 100.0).abs() < 1e-9);
        let products = product_moles(&s.fuel.fuel(), s.phi).unwrap();
        assert!(products.excess_o2 > 0.0);
        // CH4: a = 2, supplied = a/phi = 4, excess O2 = 4 - 2 = 2.
        assert!((products.excess_o2 - 2.0).abs() < 1e-9);
    }

    #[test]
    fn combustor_mesh_for_default_is_nonempty_and_in_range() {
        let s = CombustionWorkbenchState::default();
        let mesh = combustor_solid_mesh(&s).expect("default config yields a solid");
        assert!(mesh.nodes.len() > 8, "expected can + nozzle stub");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_combustion_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_combustion_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_combustion_workbench = true;
        run_combustion(&mut app.combustion);
        draw_workbench(&mut app);
    }
}
