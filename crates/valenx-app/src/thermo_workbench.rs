//! The right-side **Thermodynamics (EoS)** workbench — an in-house
//! industrial fluid-thermodynamics panel over the [`valenx_thermo`] crate
//! (two-parameter cubic equations of state + vapor–liquid phase behavior).
//!
//! The user picks a pure fluid (CO₂ / N₂ / CH₄ / H₂O / C₃H₈), an EoS model
//! ([`valenx_thermo::EosModel::PengRobinson`] / `Srk`), and sets a
//! temperature `T` [K] and pressure `P` [Pa] with labelled, AI-settable
//! DragValues. **Compute** binds the EoS to the fluid and reports the
//! liquid & vapor compressibility roots `Z` ([`Eos::z_roots`]) plus the
//! pure-component saturation pressure `Psat(T)`
//! ([`valenx_thermo::saturation_pressure`]).
//!
//! It mirrors the other real-time workbenches (`brep_workbench`): a
//! [`crate::workbench_chrome::workbench_shell`] panel gated on
//! [`crate::ValenxApp::show_thermo_workbench`], toggled from the View menu
//! and openable by the agent bridge under the workbench id `"thermo"`. The
//! bridge can set the controls (`agent_set` / `agent_control_names`), read
//! a status line (`agent_readout`), and fire the compute via the
//! RunCommand id `thermo.compute`.

use eframe::egui;

use valenx_thermo::{saturation_pressure, Eos, EosModel, Fluid};

use crate::ValenxApp;

// ---------------------------------------------------------------------------
// Workbench state
// ---------------------------------------------------------------------------

/// Which pure fluid the EoS is parametrized for. A plain `Copy` selector
/// enum mapped to the [`valenx_thermo::Fluid`] library constructors.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum FluidKind {
    /// Carbon dioxide (CO₂).
    Co2,
    /// Nitrogen (N₂).
    N2,
    /// Methane (CH₄).
    Ch4,
    /// Water (H₂O).
    H2o,
    /// Propane (C₃H₈).
    C3h8,
}

impl FluidKind {
    /// Menu label.
    fn label(self) -> &'static str {
        match self {
            FluidKind::Co2 => "Carbon dioxide (CO\u{2082})",
            FluidKind::N2 => "Nitrogen (N\u{2082})",
            FluidKind::Ch4 => "Methane (CH\u{2084})",
            FluidKind::H2o => "Water (H\u{2082}O)",
            FluidKind::C3h8 => "Propane (C\u{2083}H\u{2088})",
        }
    }

    /// Construct the matching library [`Fluid`].
    fn fluid(self) -> Fluid {
        match self {
            FluidKind::Co2 => Fluid::carbon_dioxide(),
            FluidKind::N2 => Fluid::nitrogen(),
            FluidKind::Ch4 => Fluid::methane(),
            FluidKind::H2o => Fluid::water(),
            FluidKind::C3h8 => Fluid::propane(),
        }
    }

    /// Lowercase id / alias parse for the agent bridge.
    fn from_id(s: &str) -> Option<FluidKind> {
        match s.trim().to_ascii_lowercase().as_str() {
            "co2" | "carbondioxide" | "carbon_dioxide" => Some(FluidKind::Co2),
            "n2" | "nitrogen" => Some(FluidKind::N2),
            "ch4" | "methane" => Some(FluidKind::Ch4),
            "h2o" | "water" => Some(FluidKind::H2o),
            "c3h8" | "propane" => Some(FluidKind::C3h8),
            _ => None,
        }
    }
}

/// EoS-model selector aliases for the agent bridge (the crate enum has no
/// id parser of its own).
fn eos_model_from_id(s: &str) -> Option<EosModel> {
    match s.trim().to_ascii_lowercase().as_str() {
        "pengrobinson" | "peng-robinson" | "pr" => Some(EosModel::PengRobinson),
        "srk" | "soave" | "soave-redlich-kwong" => Some(EosModel::Srk),
        _ => None,
    }
}

fn eos_model_label(m: EosModel) -> &'static str {
    match m {
        EosModel::PengRobinson => "Peng\u{2013}Robinson",
        EosModel::Srk => "Soave\u{2013}Redlich\u{2013}Kwong (SRK)",
    }
}

fn eos_model_id(m: EosModel) -> &'static str {
    match m {
        EosModel::PengRobinson => "PengRobinson",
        EosModel::Srk => "Srk",
    }
}

/// The result of a compute: the compressibility roots at `(T, P)` and the
/// saturation pressure at `T`.
struct ThermoResult {
    /// Liquid-like (smallest) compressibility root `Z`.
    z_liquid: f64,
    /// Vapor-like (largest) compressibility root `Z`.
    z_vapor: f64,
    /// Whether the cubic had three distinct real roots (two-phase possible).
    three_real: bool,
    /// Saturation pressure `Psat(T)` [Pa], or `None` if `T` is out of the
    /// subcritical range where a vapor pressure is defined.
    psat: Option<f64>,
}

/// Persistent state for the Thermodynamics workbench: the fluid + EoS
/// model selectors, the `(T, P)` state point, and the latest result.
pub struct ThermoWorkbenchState {
    /// Selected pure fluid.
    pub fluid: FluidKind,
    /// Selected cubic EoS model.
    pub model: EosModel,
    /// Temperature `T` [K].
    pub t_k: f64,
    /// Pressure `P` [Pa].
    pub p_pa: f64,

    /// Latest compute result, or `None` before the first compute.
    result: Option<ThermoResult>,
    /// Last error string (shown in the panel), cleared on a good compute.
    error: Option<String>,
}

impl Default for ThermoWorkbenchState {
    fn default() -> Self {
        // A known-good seed: CO₂ with Peng–Robinson at 350 K, 5 MPa (the
        // crate doc example — vapor Z in (0.7, 1.0)). Compute is NOT run
        // until the user (or the bridge) presses Compute.
        Self {
            fluid: FluidKind::Co2,
            model: EosModel::PengRobinson,
            t_k: 350.0,
            p_pa: 5.0e6,
            result: None,
            error: None,
        }
    }
}

impl ThermoWorkbenchState {
    /// Bind the EoS to the fluid, solve the cubic for `Z` at `(T, P)` and
    /// the saturation pressure at `T`, and store the result (or an error).
    /// Factored out so the in-panel **Compute** button and the
    /// `thermo.compute` bridge id share one path.
    fn compute_now(&mut self) {
        match self.try_compute() {
            Ok(res) => {
                self.result = Some(res);
                self.error = None;
            }
            Err(e) => {
                self.error = Some(e);
                // Keep any previous result on screen so a failed retune
                // doesn't blank the readout; the error line explains why.
            }
        }
    }

    /// The fallible compute pipeline. Separated so `compute_now` can map a
    /// typed error into the panel's error line.
    fn try_compute(&self) -> Result<ThermoResult, String> {
        let eos = Eos::new(self.fluid.fluid(), self.model);
        let z = eos
            .z_roots(self.t_k, self.p_pa)
            .map_err(|e| e.to_string())?;
        // Saturation pressure is only defined subcritically; a supercritical
        // T simply yields no Psat (reported as "n/a") rather than an error.
        let psat = saturation_pressure(&eos, self.t_k).ok();
        Ok(ThermoResult {
            z_liquid: z.liquid,
            z_vapor: z.vapor,
            three_real: z.three_real,
            psat,
        })
    }

    /// The user-visible captions of every control the agent bridge can set
    /// via `SetControl` (returned by `ListControls`). Order follows the form.
    pub fn agent_control_names() -> &'static [&'static str] {
        &["Fluid", "EoS model", "Temperature [K]", "Pressure [Pa]"]
    }

    /// Set one labelled control by its caption, for the agent `SetControl`
    /// bridge. Fail-loud on an unknown caption / wrong type / out-of-range;
    /// no state is written on error and nothing panics. `Temperature [K]`
    /// reads a positive `f64`; `Pressure [Pa]` reads a positive `f64`;
    /// `Fluid` reads a fluid id (`co2`/`n2`/`ch4`/`h2o`/`c3h8`); `EoS model`
    /// reads a model id (`pengrobinson`/`srk`).
    pub fn agent_set(
        &mut self,
        name: &str,
        value: &crate::agent_commands::AgentValue,
    ) -> Result<(), String> {
        fn positive(v: f64, what: &str, max: f64) -> Result<f64, String> {
            if !v.is_finite() || v <= 0.0 {
                return Err(format!("{what} must be finite and > 0, got {v}"));
            }
            if v > max {
                return Err(format!("{what} must be <= {max:e}, got {v}"));
            }
            Ok(v)
        }
        match name {
            "Fluid" => {
                let s = value.as_str()?;
                self.fluid = FluidKind::from_id(s)
                    .ok_or_else(|| format!("Fluid must be co2/n2/ch4/h2o/c3h8, got {s:?}"))?;
            }
            "EoS model" => {
                let s = value.as_str()?;
                self.model = eos_model_from_id(s)
                    .ok_or_else(|| format!("EoS model must be pengrobinson/srk, got {s:?}"))?;
            }
            "Temperature [K]" => self.t_k = positive(value.as_f64()?, "Temperature [K]", 5.0e3)?,
            "Pressure [Pa]" => self.p_pa = positive(value.as_f64()?, "Pressure [Pa]", 1.0e12)?,
            other => return Err(format!("unknown thermo control: {other:?}")),
        }
        Ok(())
    }

    /// The current readout text for the agent `ReadReadout` bridge: the
    /// inputs plus the compute outcome (Z roots + Psat) once a compute
    /// exists. `Some` once computed, `None` before the first compute (or
    /// after a compute error with no prior result).
    pub fn agent_readout(&self) -> Option<String> {
        if let Some(err) = &self.error {
            if self.result.is_none() {
                return Some(format!("Thermodynamics compute failed: {err}"));
            }
        }
        let r = self.result.as_ref()?;
        let psat = match r.psat {
            Some(p) => format!("{p:.4e} Pa"),
            None => "n/a (supercritical)".to_string(),
        };
        Some(format!(
            "Thermo \u{00B7} {} \u{00B7} {} \u{00B7} T={:.2} K P={:.4e} Pa \u{00B7} \
             Z_vap={:.5} Z_liq={:.5} \u{00B7} {} \u{00B7} Psat(T)={}",
            self.fluid.label(),
            eos_model_id(self.model),
            self.t_k,
            self.p_pa,
            r.z_vapor,
            r.z_liquid,
            if r.three_real {
                "three real roots"
            } else {
                "single real root"
            },
            psat,
        ))
    }
}

// ---------------------------------------------------------------------------
// Bridge run action (compute)
// ---------------------------------------------------------------------------

/// Run the compute (the in-panel **Compute** action). Factored out so the
/// button and the `thermo.compute` bridge id share one path.
pub(crate) fn run(app: &mut ValenxApp) {
    app.thermo.compute_now();
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Draw the Thermodynamics workbench. A no-op unless toggled on via
/// View → Thermodynamics (EoS).
pub fn draw_thermo_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_thermo_workbench {
        return;
    }
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_thermo_workbench",
        "Thermodynamics (EoS) \u{2014} cubic equation of state + phase behavior",
        thermo_workbench_body,
    );
    if close {
        app.show_thermo_workbench = false;
    }
}

// ---------------------------------------------------------------------------
// Workbench body
// ---------------------------------------------------------------------------

fn thermo_workbench_body(app: &mut ValenxApp, ui: &mut egui::Ui) {
    egui::ScrollArea::vertical()
        .auto_shrink([false, false])
        .show(ui, |ui| thermo_workbench_body_inner(app, ui));
}

fn thermo_workbench_body_inner(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.label(
        egui::RichText::new(
            "In-house industrial fluid thermodynamics \u{00B7} two-parameter cubic equations of \
             state (Peng\u{2013}Robinson / SRK) + vapor\u{2013}liquid phase behavior [pick a fluid, \
             a model, and a (T, P) state; Compute solves the cubic for the liquid & vapor \
             compressibility Z and the pure-component saturation pressure Psat(T)].",
        )
        .weak()
        .small(),
    );
    ui.separator();

    let s = &mut app.thermo;

    // --- Fluid + EoS model --------------------------------------------------
    ui.horizontal(|ui| {
        let lbl = ui.label("Fluid");
        let mut chosen: Option<FluidKind> = None;
        egui::ComboBox::from_id_source("thermo_fluid")
            .selected_text(s.fluid.label())
            .show_ui(ui, |ui| {
                for k in [
                    FluidKind::Co2,
                    FluidKind::N2,
                    FluidKind::Ch4,
                    FluidKind::H2o,
                    FluidKind::C3h8,
                ] {
                    if ui.selectable_label(s.fluid == k, k.label()).clicked() {
                        chosen = Some(k);
                    }
                }
            })
            .response
            .labelled_by(lbl.id);
        if let Some(k) = chosen {
            s.fluid = k;
        }
    });

    ui.horizontal(|ui| {
        let lbl = ui.label("EoS model");
        let mut chosen: Option<EosModel> = None;
        egui::ComboBox::from_id_source("thermo_eos_model")
            .selected_text(eos_model_label(s.model))
            .show_ui(ui, |ui| {
                for m in [EosModel::PengRobinson, EosModel::Srk] {
                    if ui
                        .selectable_label(s.model == m, eos_model_label(m))
                        .clicked()
                    {
                        chosen = Some(m);
                    }
                }
            })
            .response
            .labelled_by(lbl.id);
        if let Some(m) = chosen {
            s.model = m;
        }
    });

    ui.add_space(6.0);

    // --- State point (T, P) -------------------------------------------------
    egui::Grid::new("thermo_state")
        .num_columns(2)
        .striped(true)
        .show(ui, |ui| {
            let lbl = ui.label("Temperature [K]");
            ui.add(
                egui::DragValue::new(&mut s.t_k)
                    .speed(1.0)
                    .range(0.1..=5.0e3)
                    .max_decimals(3),
            )
            .labelled_by(lbl.id)
            .on_hover_text("Absolute temperature of the state point, in kelvin.");
            ui.end_row();

            let lbl = ui.label("Pressure [Pa]");
            ui.add(
                egui::DragValue::new(&mut s.p_pa)
                    .speed(1.0e4)
                    .range(1.0..=1.0e12)
                    .max_decimals(1),
            )
            .labelled_by(lbl.id)
            .on_hover_text("Absolute pressure of the state point, in pascal.");
            ui.end_row();
        });

    // --- Compute ------------------------------------------------------------
    ui.add_space(6.0);
    if ui
        .button("\u{25B6} Compute")
        .on_hover_text(
            "Bind the EoS to the fluid, solve the cubic for the liquid & vapor compressibility Z \
             at (T, P), and compute the saturation pressure Psat(T).",
        )
        .clicked()
    {
        s.compute_now();
    }

    ui.add_space(6.0);
    ui.separator();
    draw_result(s, ui);
}

// ---------------------------------------------------------------------------
// Result render
// ---------------------------------------------------------------------------

fn draw_result(s: &ThermoWorkbenchState, ui: &mut egui::Ui) {
    if let Some(err) = &s.error {
        ui.label(
            egui::RichText::new(format!("Compute error: {err}"))
                .color(egui::Color32::from_rgb(220, 110, 90))
                .strong(),
        );
        ui.add_space(4.0);
    }

    let Some(r) = s.result.as_ref() else {
        ui.label(
            egui::RichText::new(
                "No result yet \u{2014} pick a fluid + model + (T, P), then press Compute.",
            )
            .italics()
            .weak(),
        );
        return;
    };

    ui.label(
        egui::RichText::new(format!(
            "Compressibility Z: vapor {:.5} \u{00B7} liquid {:.5} \u{00B7} {}",
            r.z_vapor,
            r.z_liquid,
            if r.three_real {
                "three real roots (two-phase possible)"
            } else {
                "single real root"
            },
        ))
        .strong()
        .color(egui::Color32::from_rgb(150, 200, 230)),
    );
    match r.psat {
        Some(p) => ui.label(format!(
            "Saturation pressure Psat(T) = {p:.4e} Pa ({:.4} MPa)",
            p / 1.0e6
        )),
        None => {
            ui.label(egui::RichText::new("Saturation pressure: n/a (T is supercritical)").italics())
        }
    };
}

// ---------------------------------------------------------------------------
// Tests (unit)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_seed_is_co2_pr() {
        let s = ThermoWorkbenchState::default();
        assert_eq!(s.fluid, FluidKind::Co2);
        assert_eq!(s.model, EosModel::PengRobinson);
        assert!(s.result.is_none());
    }

    #[test]
    fn compute_produces_a_result_and_readout() {
        let mut s = ThermoWorkbenchState::default();
        s.compute_now();
        assert!(
            s.error.is_none(),
            "default seed should compute: {:?}",
            s.error
        );
        let r = s.result.as_ref().expect("a result after compute");
        // CO₂ vapor Z at 350 K / 5 MPa is in (0.7, 1.0) (crate doc example).
        assert!(
            r.z_vapor > 0.7 && r.z_vapor < 1.0,
            "vapor Z out of expected range: {}",
            r.z_vapor
        );
        let readout = s.agent_readout().expect("readout after compute");
        assert!(readout.contains("Thermo"), "readout: {readout}");
        assert!(readout.contains("Z_vap"), "readout names Z: {readout}");
    }

    #[test]
    fn saturation_pressure_co2_matches_experiment() {
        // CO₂ Psat at 273.15 K ≈ 3.49 MPa (the crate's validated example).
        let mut s = ThermoWorkbenchState::default();
        s.t_k = 273.15;
        s.compute_now();
        let r = s.result.as_ref().expect("result");
        let psat = r.psat.expect("subcritical T has a Psat");
        assert!((psat - 3.49e6).abs() / 3.49e6 < 0.05, "Psat off: {psat:e}");
    }

    #[test]
    fn agent_set_numeric_and_enum_controls() {
        use crate::agent_commands::AgentValue;
        let mut s = ThermoWorkbenchState::default();
        s.agent_set("Temperature [K]", &AgentValue::Float(300.0))
            .expect("T");
        assert!((s.t_k - 300.0).abs() < 1e-9);
        s.agent_set("Pressure [Pa]", &AgentValue::Float(1.0e6))
            .expect("P");
        assert!((s.p_pa - 1.0e6).abs() < 1e-3);
        s.agent_set("Fluid", &AgentValue::Str("methane".into()))
            .expect("fluid");
        assert_eq!(s.fluid, FluidKind::Ch4);
        s.agent_set("EoS model", &AgentValue::Str("srk".into()))
            .expect("model");
        assert_eq!(s.model, EosModel::Srk);
    }

    #[test]
    fn agent_set_rejects_bad_values() {
        use crate::agent_commands::AgentValue;
        let mut s = ThermoWorkbenchState::default();
        assert!(s
            .agent_set("Temperature [K]", &AgentValue::Float(-1.0))
            .is_err());
        assert!(s
            .agent_set("Pressure [Pa]", &AgentValue::Float(f64::NAN))
            .is_err());
        assert!(s
            .agent_set("Fluid", &AgentValue::Str("xenon".into()))
            .is_err());
        assert!(s
            .agent_set("EoS model", &AgentValue::Str("nope".into()))
            .is_err());
        assert!(s.agent_set("bogus", &AgentValue::Float(1.0)).is_err());
    }

    #[test]
    fn readout_is_none_before_compute() {
        let s = ThermoWorkbenchState::default();
        assert!(s.agent_readout().is_none(), "no readout before a compute");
    }

    #[test]
    fn run_bridge_helper_computes_through_app() {
        let mut app = ValenxApp::default();
        run(&mut app);
        assert!(
            app.thermo.result.is_some(),
            "the thermo.compute bridge helper should produce a result"
        );
    }

    #[test]
    fn control_names_are_listed() {
        let names = ThermoWorkbenchState::agent_control_names();
        for c in ["Fluid", "EoS model", "Temperature [K]", "Pressure [Pa]"] {
            assert!(names.contains(&c), "missing control name {c}");
        }
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use egui::accesskit::{Node, NodeId, Role};

    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_thermo_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    fn has_named_node(nodes: &[(NodeId, Node)], name: &str) -> bool {
        nodes.iter().any(|(_, n)| n.name() == Some(name))
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_thermo_workbench);
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_thermo_workbench(&mut app, ctx);
        });
        // No panic = pass.
    }

    #[test]
    fn workbench_draws_when_shown() {
        let mut app = ValenxApp::default();
        app.show_thermo_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);
        assert!(!nodes.is_empty(), "a shown workbench produces a11y nodes");
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // Every numeric DragValue is a SpinButton and must be `labelled_by`
        // its caption so an AI / screen reader can find it by caption text.
        let mut app = ValenxApp::default();
        app.show_thermo_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let by_id: std::collections::HashMap<NodeId, &Node> =
            nodes.iter().map(|(id, n)| (*id, n)).collect();

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        // T and P are the two numeric DragValues.
        assert!(
            spin_buttons.len() >= 2,
            "expected the T/P numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every DragValue must be labelled_by a caption (AI-drivable name)"
        );
        assert!(
            spin_buttons.iter().all(|n| {
                n.labelled_by()
                    .iter()
                    .any(|id| by_id.get(id).is_some_and(|t| t.name().is_some()))
            }),
            "every DragValue's labelled_by must point at a named caption node"
        );

        for caption in ["Temperature [K]", "Pressure [Pa]"] {
            assert!(
                has_named_node(&nodes, caption),
                "caption '{caption}' must be a named node in the a11y tree"
            );
        }
    }
}
