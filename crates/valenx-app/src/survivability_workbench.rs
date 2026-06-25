//! The right-side **Survivability / Protection Workbench** panel — a native
//! front-end over the in-house `valenx-survivability` crate.
//!
//! This workbench is the **M7 track** of the Valenx defense modeling-&-simulation
//! roadmap and the **protective / defensive** side of the shared impact/blast
//! physics. It answers exactly one question, and only that question:
//!
//! > **Given a threat — a blast at a stand-off, or an impacting projectile —
//! > what is the minimum protection that lets a structure and its occupants
//! > *survive*?**
//!
//! It is the same physics as **civil blast-resistant building design**
//! (UFC 3-340-02, ASCE 59) and **automotive crash safety**. The user describes a
//! threat and a protective element and reads back four protective screens:
//!
//! - **Blast loading** — the free-field peak side-on overpressure `Pso` and
//!   positive-phase impulse `i_s` versus the Hopkinson–Cranz scaled distance
//!   `Z = R / W^(1/3)`, assembled into a Friedlander pressure pulse `P(t)`.
//! - **Protective structural response** — the single-degree-of-freedom (SDOF)
//!   dynamic response of a protective wall / plate / panel to that pulse: peak
//!   deflection, ductility ratio `μ = x_max / x_yield`, and the dynamic
//!   amplification factor (DAF).
//! - **Pressure–impulse (P–I) diagram** — the canonical iso-damage chart with
//!   its impulsive (vertical) and quasi-static (horizontal) asymptotes, and the
//!   current `(Pso, i_s)` design point plotted so the user can read "safe" vs
//!   "exceeds the protective limit" at a glance.
//! - **Armor / impact protection sizing** — the minimum protective plate
//!   thickness and areal density that just defeats a projectile of a given
//!   kinetic energy, plus an occupant acceleration-injury margin.
//!
//! ## Dual-use boundary (hard line — protective only)
//!
//! This workbench, like the crate it drives, is **protective / defensive only**.
//! It models how to *protect* against a threat, **never** how to penetrate,
//! defeat, or optimize a weapon against protection. There are **no** warhead,
//! lethality, fragmentation, shaped-charge, or penetration-optimization models
//! anywhere. The blast *source* appears only as the standard idealized far-field
//! overpressure curve used as a design *load*; the armor model sizes the
//! *minimum plate that keeps a threat out*, never a means of getting through one.
//! Every output is framed as "minimum protection to survive threat X". The UI
//! states explicitly that **no penetration / lethality is modeled**.
//!
//! Mirrors the other workbenches (`missionsim_workbench`, `uas_workbench`): a
//! [`crate::workbench_chrome::workbench_shell`] panel gated on
//! [`crate::ValenxApp::show_survivability_workbench`], toggled from the View menu
//! and openable by the agent bridge under the workbench id `"survivability"`
//! (aliases `"protection"` / `"blast"`; see [`crate::project_tabs::TabKind`] /
//! [`crate::agent_commands`]).
//!
//! Honesty: `valenx-survivability` is **research / educational grade,
//! validation-pending** — every model is a well-established open-literature
//! protective-design equation (Brode 1955 overpressure, Newmark–Hansen 1961
//! impulse, the Baker/Smith–Hetherington P–I hyperbola, a Recht–Ipson plug-shear
//! ballistic-limit, an Eiband-style occupant tolerance screen), but absolute
//! numbers depend on uncertain dynamic material data and the assumed failure
//! mechanism, so the outputs are fast trade-study *screens*, not certified blast
//! or ballistic ratings. Every error from `valenx-survivability` surfaces
//! verbatim in-panel; degenerate parameters (a zero / negative charge, stand-off,
//! mass, stiffness, area, diameter, strength, or tolerance, or a scaled distance
//! outside the validated band) show an in-panel error, NOT a panic.

use eframe::egui;
use valenx_survivability::{
    assess_occupant, minimum_protection, sdof_response, BlastLoad, OccupantAssessment, PiPoint,
    PressureImpulseDiagram, SdofResponse, Threat, TransientControls,
};

use crate::ValenxApp;

// ---------------------------------------------------------------------------
// Parameters
// ---------------------------------------------------------------------------

/// Editable survivability inputs: the blast threat (charge / stand-off), the
/// protective element (effective mass / stiffness / yield deflection / loaded
/// area → the SDOF system + P–I diagram), the armor threat (projectile mass /
/// velocity / diameter and the plate's dynamic shear strength + density), and
/// the occupant whole-body tolerance.
#[derive(Clone, Debug)]
pub struct SurvivabilityParams {
    // -- Blast threat (the design load) --
    /// TNT-equivalent charge mass `W` (kg), finite and positive.
    pub charge_kg: f64,
    /// Stand-off distance `R` (m) from charge centre to the element, finite and
    /// positive. (Together `W`, `R` set the scaled distance `Z = R/W^(1/3)`,
    /// which must land in the validated band.)
    pub standoff_m: f64,

    // -- Protective element (SDOF + P–I) --
    /// Effective (lumped) mass `m` (kg) of the protective wall / plate / panel.
    pub element_mass_kg: f64,
    /// Effective stiffness `k` (N/m) of the protective element.
    pub element_stiffness_n_m: f64,
    /// Effective damping `c` (N·s/m) of the protective element (`>= 0`).
    pub element_damping_n_s_m: f64,
    /// Yield deflection `x_yield` (m) — the deflection at first yield, sets the
    /// ductility ratio `μ = x_max / x_yield`. Also used as the elastic
    /// deflection limit for the P–I iso-damage curve.
    pub yield_deflection_m: f64,
    /// Loaded tributary area `A` (m²) that converts the side-on overpressure
    /// `Pso` (Pa) into the effective force on the single DOF (`F = Pso · A`),
    /// per the reused SDOF-blast convention.
    pub loaded_area_m2: f64,

    // -- Armor / impact threat --
    /// Projectile mass `m_p` (kg).
    pub projectile_mass_kg: f64,
    /// Projectile impact speed `v` (m/s).
    pub projectile_velocity_m_s: f64,
    /// Projectile (contact) diameter `d` (m) — sets the sheared plug area.
    pub projectile_diameter_m: f64,
    /// Protective plate dynamic shear strength `τ_d` (Pa).
    pub plate_shear_strength_pa: f64,
    /// Protective plate material density `ρ` (kg/m³).
    pub plate_density_kg_m3: f64,

    // -- Occupant tolerance --
    /// Whole-body acceleration tolerance limit `g_tol` (g) for the relevant
    /// pulse duration / direction. The occupant margin is `g_tol / g_peak`.
    pub occupant_tolerance_g: f64,
}

impl Default for SurvivabilityParams {
    fn default() -> Self {
        Self {
            // 100 kg TNT at 20 m: Z ≈ 4.31 m/kg^(1/3) — comfortably inside the
            // validated [0.1, 40] band; a realistic protective-design screen.
            charge_kg: 100.0,
            standoff_m: 20.0,
            // A wall panel reduced to an effective SDOF (mirrors the crate's
            // end-to-end protective screen).
            element_mass_kg: 500.0,
            element_stiffness_n_m: 5.0e6,
            element_damping_n_s_m: 100.0,
            yield_deflection_m: 0.02,
            loaded_area_m2: 1.0,
            // 10 g compact projectile at 800 m/s onto a steel plate.
            projectile_mass_kg: 0.01,
            projectile_velocity_m_s: 800.0,
            projectile_diameter_m: 0.008,
            plate_shear_strength_pa: 5.0e8, // ~500 MPa dynamic shear (steel)
            plate_density_kg_m3: 7850.0,    // steel
            // A representative short-pulse whole-body tolerance (caller's choice).
            occupant_tolerance_g: 50.0,
        }
    }
}

// ---------------------------------------------------------------------------
// Result
// ---------------------------------------------------------------------------

/// One sample of the Friedlander pressure pulse `P(t)` for the painter.
#[derive(Clone, Copy, Debug)]
pub struct PressurePoint {
    /// Time `t` (s) from shock-front arrival.
    pub t: f64,
    /// Overpressure `P(t)` (Pa) at that time.
    pub p: f64,
}

/// Cached survivability output for the painter + readouts.
#[derive(Clone)]
pub struct SurvivabilityResult {
    /// The assembled free-field blast load (descriptors + Friedlander pulse).
    pub blast: BlastLoad,
    /// Sampled Friedlander pressure-time curve `P(t)` for the plot.
    pub pressure_curve: Vec<PressurePoint>,
    /// The SDOF protective-element response (peak deflection, ductility, DAF).
    pub response: SdofResponse,
    /// The pressure–impulse iso-damage diagram (the two asymptotes).
    pub pi: PressureImpulseDiagram,
    /// Sampled P–I iso-damage curve for the plot.
    pub pi_curve: Vec<PiPoint>,
    /// Whether the current `(Pso, i_s)` design point is safe (at / below the
    /// iso-damage curve).
    pub design_point_safe: bool,
    /// The minimum-protection armor sizing (thickness + areal density).
    pub armor: valenx_survivability::ArmorSizing,
    /// The occupant acceleration-injury screen.
    pub occupant: OccupantAssessment,
}

// ---------------------------------------------------------------------------
// Workbench state
// ---------------------------------------------------------------------------

/// Persistent state for the survivability workbench.
#[derive(Default)]
pub struct SurvivabilityWorkbenchState {
    /// User-editable parameters.
    pub params: SurvivabilityParams,
    /// Last successful result (populated after a successful run).
    pub result: Option<SurvivabilityResult>,
    /// Status / error line shown below the controls.
    pub status: String,
}

impl SurvivabilityWorkbenchState {
    /// The user-visible captions of every control the agent bridge can set via
    /// `SetControl` (see [`crate::agent_commands`]). Returned by `ListControls`
    /// so an agent can discover the name space. Every control is an `f64` drag
    /// value; order follows the form. The two non-ASCII captions carry the same
    /// `·` / superscript characters the UI draws, so the agent must match them
    /// exactly (or read them from this list).
    pub fn agent_control_names() -> &'static [&'static str] {
        &[
            "charge mass W (kg)",
            "stand-off R (m)",
            "element mass (kg)",
            "element stiffness (N/m)",
            "element damping (N\u{00B7}s/m)",
            "yield deflection (m)",
            "loaded area (m\u{00B2})",
            "projectile mass (kg)",
            "projectile velocity (m/s)",
            "projectile diameter (m)",
            "plate shear strength (Pa)",
            "plate density (kg/m\u{00B3})",
            "occupant tolerance (g)",
        ]
    }

    /// Set one labelled control by its user-visible caption, for the agent
    /// `SetControl` bridge. Captions match exactly what the form draws (incl. the
    /// two non-ASCII captions). Fail-loud on an unknown caption / wrong type (the
    /// bridge posts a `warn` note); no field is written on error and nothing
    /// panics. Every control is an `f64` drag value ([`AgentValue::as_f64`]); the
    /// downstream `valenx-survivability` constructors re-validate the physical
    /// ranges on the next `run`, so this setter only type-checks.
    pub fn agent_set(
        &mut self,
        name: &str,
        value: &crate::agent_commands::AgentValue,
    ) -> Result<(), String> {
        let p = &mut self.params;
        match name {
            // -- Blast threat --
            "charge mass W (kg)" => p.charge_kg = value.as_f64()?,
            "stand-off R (m)" => p.standoff_m = value.as_f64()?,
            // -- Protective element (SDOF) --
            "element mass (kg)" => p.element_mass_kg = value.as_f64()?,
            "element stiffness (N/m)" => p.element_stiffness_n_m = value.as_f64()?,
            "element damping (N\u{00B7}s/m)" => p.element_damping_n_s_m = value.as_f64()?,
            "yield deflection (m)" => p.yield_deflection_m = value.as_f64()?,
            "loaded area (m\u{00B2})" => p.loaded_area_m2 = value.as_f64()?,
            // -- Armor / impact threat --
            "projectile mass (kg)" => p.projectile_mass_kg = value.as_f64()?,
            "projectile velocity (m/s)" => p.projectile_velocity_m_s = value.as_f64()?,
            "projectile diameter (m)" => p.projectile_diameter_m = value.as_f64()?,
            "plate shear strength (Pa)" => p.plate_shear_strength_pa = value.as_f64()?,
            "plate density (kg/m\u{00B3})" => p.plate_density_kg_m3 = value.as_f64()?,
            // -- Occupant tolerance --
            "occupant tolerance (g)" => p.occupant_tolerance_g = value.as_f64()?,
            other => return Err(format!("unknown survivability control: {other:?}")),
        }
        Ok(())
    }

    /// The current computed-result text for the agent `ReadReadout` bridge (see
    /// [`crate::agent_commands`]). This workbench keeps its result as a structured
    /// [`SurvivabilityResult`] and renders a one-line `status` summary (a `✔ …`
    /// line on success, a `⚠ …` line on error) — that same `status` string is
    /// returned here. `None` when it is empty, i.e. the pipeline has not been run
    /// yet. Read-only — lets an agent read the answer back after driving a run,
    /// closing the live-driving loop.
    pub fn agent_readout(&self) -> Option<String> {
        if self.status.is_empty() {
            None
        } else {
            Some(self.status.clone())
        }
    }

    /// Run the full survivability pipeline: assemble the free-field blast load
    /// (`Pso`, `i_s`, Friedlander pulse), march the SDOF protective element under
    /// the pulse (peak deflection / ductility / DAF), build the P–I iso-damage
    /// diagram and test the current design point, size the minimum protective
    /// plate for the projectile threat, and screen the occupant against the
    /// whole-body tolerance.
    ///
    /// Every failure is returned as an `Err(String)` — no panics, no invented
    /// numbers. Degenerate inputs (a non-positive charge / stand-off / mass /
    /// stiffness / area / diameter / strength / density / tolerance, a scaled
    /// distance outside the validated band, or a non-positive yield deflection)
    /// surface `valenx-survivability`'s own error verbatim.
    pub fn run(&self) -> Result<SurvivabilityResult, String> {
        let p = &self.params;

        // --- 1. Free-field blast load (the design load) ---------------------
        let blast =
            BlastLoad::tnt_free_air(p.charge_kg, p.standoff_m).map_err(|e| e.to_string())?;
        let pressure_pulse = blast.friedlander().map_err(|e| e.to_string())?;

        // Sample the positive phase (arrival -> first zero crossing) plus a
        // little of the suction tail for the P(t) plot.
        let pressure_curve = sample_pressure_curve(&pressure_pulse, blast.positive_duration_s);

        // --- 2. SDOF protective-element response ----------------------------
        // The reused SDOF-blast solver treats the pulse peak as the *force* on
        // the single DOF, so convert the side-on overpressure into a force by
        // the loaded tributary area: F_peak = Pso · A. A must be > 0.
        let area = require_pos("loaded area A", p.loaded_area_m2)?;
        let force_pulse = valenx_survivability::FriedlanderPulse::new(
            blast.peak_overpressure_pa * area,
            blast.positive_duration_s,
            blast.decay_b,
        )
        .map_err(|e| e.to_string())?;
        // Resolve the positive phase and several structural periods: pick a step
        // that resolves the load rise, then march well past the first peak.
        let controls = sdof_controls(
            p.element_mass_kg,
            p.element_stiffness_n_m,
            blast.positive_duration_s,
        )?;
        let response = sdof_response(
            p.element_mass_kg,
            p.element_damping_n_s_m,
            p.element_stiffness_n_m,
            Some(p.yield_deflection_m),
            &force_pulse,
            &controls,
        )
        .map_err(|e| e.to_string())?;

        // --- 3. Pressure–impulse iso-damage diagram + design point ----------
        // The elastic deflection limit for the P–I curve is the yield
        // deflection. The asymptotes are P_min = k·x_lim (per loaded area: a
        // pressure) and i_min = x_lim·√(km). The current design point is the
        // free-field (Pso, i_s); a per-area resistance comparison is consistent
        // because both P_min and Pso are pressures here (unit loaded area).
        let pi = PressureImpulseDiagram::elastic(
            p.element_mass_kg,
            p.element_stiffness_n_m,
            p.yield_deflection_m,
        )
        .map_err(|e| e.to_string())?;
        let pi_curve = pi.curve(160, 100.0);
        let design_point_safe = pi.is_safe(blast.peak_overpressure_pa, blast.impulse_pa_s);

        // --- 4. Minimum armor protection for the projectile threat ----------
        let threat = Threat::new(
            p.projectile_mass_kg,
            p.projectile_velocity_m_s,
            p.projectile_diameter_m,
        )
        .map_err(|e| e.to_string())?;
        let armor = minimum_protection(&threat, p.plate_density_kg_m3, p.plate_shear_strength_pa)
            .map_err(|e| e.to_string())?;

        // --- 5. Occupant acceleration-injury screen -------------------------
        // The peak acceleration comes straight from the SDOF response.
        let occupant = assess_occupant(response.peak_acceleration_m_s2, p.occupant_tolerance_g)
            .map_err(|e| e.to_string())?;

        Ok(SurvivabilityResult {
            blast,
            pressure_curve,
            response,
            pi,
            pi_curve,
            design_point_safe,
            armor,
            occupant,
        })
    }
}

/// Validate a finite, strictly-positive quantity (used for the loaded area,
/// which the workbench owns rather than the crate), fail-loud.
fn require_pos(name: &str, v: f64) -> Result<f64, String> {
    if v.is_finite() && v > 0.0 {
        Ok(v)
    } else {
        Err(format!("{name} must be finite and > 0, got {v}"))
    }
}

/// Choose [`TransientControls`] that resolve the blast load rise *and* march the
/// SDOF system well past its first deflection peak (≈ a quarter period), so the
/// peak response is captured. Fail-loud on a non-positive mass / stiffness /
/// duration (the crate also validates, but we need `ω` here for the step).
fn sdof_controls(
    mass: f64,
    stiffness: f64,
    positive_duration_s: f64,
) -> Result<TransientControls, String> {
    let m = require_pos("effective mass m", mass)?;
    let k = require_pos("effective stiffness k", stiffness)?;
    let td = require_pos("positive duration t_d", positive_duration_s)?;
    let omega = (k / m).sqrt();
    let period = 2.0 * std::f64::consts::PI / omega;
    // Step resolves the faster of (load rise, structural period); 1/50 of it.
    let dt = (td.min(period) / 50.0).max(1.0e-9);
    // March across several periods plus several positive durations so the peak
    // (which can occur after the pulse for a stiff element) is always seen.
    let span = 4.0 * period + 4.0 * td;
    let n_steps = ((span / dt).ceil() as usize).clamp(50, 200_000);
    Ok(TransientControls {
        dt,
        n_steps,
        newmark: valenx_survivability::NewmarkBeta::average_acceleration(),
    })
}

/// Sample the Friedlander pulse `P(t)` from arrival to `1.4 · t_d` (the positive
/// phase plus a little of the suction tail) at a fixed number of points.
fn sample_pressure_curve(
    pulse: &valenx_survivability::FriedlanderPulse,
    positive_duration_s: f64,
) -> Vec<PressurePoint> {
    const SAMPLES: usize = 96;
    let t_end = positive_duration_s * 1.4;
    (0..=SAMPLES)
        .map(|i| {
            let t = t_end * i as f64 / SAMPLES as f64;
            PressurePoint {
                t,
                p: pulse.pressure_at(t),
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Draw the survivability workbench. A no-op unless toggled on via
/// View → Survivability / protection.
///
/// Mirrors [`crate::missionsim_workbench::draw_missionsim_workbench`].
pub fn draw_survivability_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_survivability_workbench {
        return;
    }
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_survivability_workbench",
        "Survivability / protection",
        survivability_workbench_body,
    );
    if close {
        app.show_survivability_workbench = false;
    }
}

// ---------------------------------------------------------------------------
// Workbench body
// ---------------------------------------------------------------------------

fn survivability_workbench_body(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.label(
        egui::RichText::new(
            "Defensive survivability & protection \u{00B7} valenx-survivability  \
             [research / educational \u{2014} minimum protection to survive threat X; \
             NO penetration / lethality is modeled]",
        )
        .weak()
        .small(),
    );
    ui.separator();

    let mut do_run = false;

    {
        let s = &mut app.survivability;
        let p = &mut s.params;

        // --- Blast threat (design load) -------------------------------------
        ui.label(egui::RichText::new("Blast threat (design load)").strong());
        ui.label(
            egui::RichText::new(
                "Idealized free-field overpressure used as a design LOAD (Brode 1955 / \
                 Newmark-Hansen 1961). Not a warhead model.",
            )
            .weak()
            .small(),
        );
        egui::Grid::new("survivability_blast_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                drag_row(
                    ui,
                    "charge mass W (kg)",
                    &mut p.charge_kg,
                    1.0,
                    "TNT-equivalent charge mass. Finite and > 0. With stand-off sets Z = R/W^(1/3).",
                );
                drag_row(
                    ui,
                    "stand-off R (m)",
                    &mut p.standoff_m,
                    1.0,
                    "Stand-off from charge centre to the protected element. Finite and > 0. \
                     Z = R/W^(1/3) must be in [0.1, 40].",
                );
            });

        // --- Protective element (SDOF + P–I) --------------------------------
        ui.add_space(4.0);
        ui.label(egui::RichText::new("Protective element (SDOF)").strong());
        ui.label(
            egui::RichText::new(
                "A wall / plate / panel reduced to one effective degree of freedom. Yield \
                 deflection is also the P-I elastic deflection limit.",
            )
            .weak()
            .small(),
        );
        egui::Grid::new("survivability_element_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                drag_row(
                    ui,
                    "element mass (kg)",
                    &mut p.element_mass_kg,
                    10.0,
                    "Effective lumped mass m of the protective element. Finite and > 0.",
                );
                drag_row(
                    ui,
                    "element stiffness (N/m)",
                    &mut p.element_stiffness_n_m,
                    1.0e5,
                    "Effective stiffness k of the protective element. Finite and > 0.",
                );
                drag_row(
                    ui,
                    "element damping (N\u{00B7}s/m)",
                    &mut p.element_damping_n_s_m,
                    10.0,
                    "Effective damping c of the protective element. Finite and >= 0.",
                );
                drag_row(
                    ui,
                    "yield deflection (m)",
                    &mut p.yield_deflection_m,
                    0.001,
                    "Deflection at first yield x_yield. Sets ductility ratio and the P-I \
                     elastic deflection limit. Finite and > 0.",
                );
                drag_row(
                    ui,
                    "loaded area (m\u{00B2})",
                    &mut p.loaded_area_m2,
                    0.1,
                    "Tributary area A loaded by the blast; converts overpressure to the SDOF \
                     force (F = Pso\u{00B7}A). Finite and > 0.",
                );
            });

        // --- Armor / impact threat ------------------------------------------
        ui.add_space(4.0);
        ui.label(egui::RichText::new("Armor / impact protection").strong());
        ui.label(
            egui::RichText::new(
                "Sizes the MINIMUM plate that keeps the threat OUT (plug-shear energy balance). \
                 There is no model for getting through a plate.",
            )
            .weak()
            .small(),
        );
        egui::Grid::new("survivability_armor_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                drag_row(
                    ui,
                    "projectile mass (kg)",
                    &mut p.projectile_mass_kg,
                    0.001,
                    "Projectile mass m_p. Finite and > 0. With velocity sets the threat KE = 1/2 m v^2.",
                );
                drag_row(
                    ui,
                    "projectile velocity (m/s)",
                    &mut p.projectile_velocity_m_s,
                    10.0,
                    "Projectile impact speed v. Finite and > 0.",
                );
                drag_row(
                    ui,
                    "projectile diameter (m)",
                    &mut p.projectile_diameter_m,
                    0.001,
                    "Projectile / contact diameter d (sets the sheared plug area). Finite and > 0.",
                );
                drag_row(
                    ui,
                    "plate shear strength (Pa)",
                    &mut p.plate_shear_strength_pa,
                    1.0e7,
                    "Plate dynamic shear strength tau_d (~0.5-0.6x dynamic UTS). Finite and > 0.",
                );
                drag_row(
                    ui,
                    "plate density (kg/m\u{00B3})",
                    &mut p.plate_density_kg_m3,
                    50.0,
                    "Plate material density rho (sets areal density). Finite and > 0.",
                );
            });

        // --- Occupant tolerance ---------------------------------------------
        ui.add_space(4.0);
        ui.label(egui::RichText::new("Occupant tolerance").strong());
        egui::Grid::new("survivability_occupant_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                drag_row(
                    ui,
                    "occupant tolerance (g)",
                    &mut p.occupant_tolerance_g,
                    1.0,
                    "Whole-body acceleration tolerance g_tol for this pulse duration / direction. \
                     Margin = g_tol / g_peak. Finite and > 0.",
                );
            });

        ui.add_space(6.0);
        ui.horizontal(|ui| {
            if ui
                .button(egui::RichText::new("Run").strong())
                .on_hover_text(
                    "Assemble the blast load, march the SDOF protective response, build the \
                     P-I diagram + design point, size the minimum armor, and screen the occupant.",
                )
                .clicked()
            {
                do_run = true;
            }
        });
    }

    // --- Execute (outside borrow) -------------------------------------------
    if do_run {
        run_and_store(app);
    }

    // --- Status line ---------------------------------------------------------
    let s = &app.survivability;
    if !s.status.is_empty() {
        ui.add_space(6.0);
        let color = if s.status.starts_with('\u{26A0}') {
            egui::Color32::from_rgb(220, 120, 60)
        } else {
            egui::Color32::from_rgb(90, 180, 110)
        };
        ui.label(egui::RichText::new(&s.status).color(color).strong());
    }

    // --- Visualisation -------------------------------------------------------
    ui.add_space(6.0);
    ui.separator();
    draw_survivability_viz(s, ui);
}

/// A labelled `DragValue` row in a 2-column grid: a caption cell, the drag value
/// (`labelled_by` the caption so it carries an accessible name), then `end_row`.
/// Mirrors the missionsim workbench's per-row caption pattern.
fn drag_row(ui: &mut egui::Ui, caption: &str, value: &mut f64, speed: f64, hint: &str) {
    let lbl = ui.label(caption);
    ui.add(egui::DragValue::new(value).speed(speed))
        .labelled_by(lbl.id)
        .on_hover_text(hint);
    ui.end_row();
}

/// Run the pipeline and fold the result (or error) into the workbench status.
/// Factored out so the Run button (and tests) can share it.
pub(crate) fn run_and_store(app: &mut ValenxApp) {
    let s = &mut app.survivability;
    match s.run() {
        Ok(res) => {
            let duct = res
                .response
                .ductility_ratio
                .map(|m| format!("{m:.2}"))
                .unwrap_or_else(|| "\u{2014}".to_string());
            let pi_state = if res.design_point_safe {
                "SAFE"
            } else {
                "EXCEEDS limit"
            };
            s.status = format!(
                "\u{2714} Pso {:.1} kPa \u{00B7} \u{03BC} {} \u{00B7} P-I {} \u{00B7} \
                 armor {:.1} mm \u{00B7} occupant {}",
                res.blast.peak_overpressure_pa / 1.0e3,
                duct,
                pi_state,
                res.armor.min_thickness_m * 1.0e3,
                if res.occupant.survivable {
                    "OK"
                } else {
                    "OVER tolerance"
                },
            );
            s.result = Some(res);
        }
        Err(e) => {
            s.status = format!("\u{26A0} {e}");
            s.result = None;
        }
    }
}

// ---------------------------------------------------------------------------
// 2-D visualisation (pressure-time + P–I diagram + armor / occupant readouts)
// ---------------------------------------------------------------------------

fn draw_survivability_viz(s: &SurvivabilityWorkbenchState, ui: &mut egui::Ui) {
    let Some(res) = &s.result else {
        ui.label(
            egui::RichText::new(
                "press \"Run\" to compute the blast load P(t), the P-I iso-damage diagram, \
                 the minimum armor sizing, and the occupant margin",
            )
            .weak(),
        );
        return;
    };

    draw_pressure_curve(res, ui);
    ui.add_space(8.0);
    draw_pi_diagram(res, ui);
    ui.add_space(8.0);
    draw_armor_readout(res, ui);
    ui.add_space(8.0);
    draw_occupant_readout(res, ui);
}

/// View (a): the **Friedlander pressure-time curve** `P(t)` — the design blast
/// load on the protective element over the positive phase and a little of the
/// suction tail.
fn draw_pressure_curve(res: &SurvivabilityResult, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Blast pressure-time P(t)").strong());
    ui.label(
        egui::RichText::new(
            "cyan = side-on overpressure P(t) \u{00B7} dashed = zero \u{00B7} x = time since \
             arrival (ms)",
        )
        .weak()
        .small(),
    );

    let available = ui.available_size();
    let (rect, _) = ui.allocate_exact_size(
        egui::vec2(available.x.min(460.0), 180.0),
        egui::Sense::hover(),
    );
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 0.0, egui::Color32::from_rgb(12, 22, 32));

    if res.pressure_curve.len() < 2 {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "too few pressure points",
            egui::FontId::monospace(12.0),
            egui::Color32::from_gray(120),
        );
        return;
    }

    let mut t_hi = f64::NEG_INFINITY;
    let mut p_lo = f64::INFINITY;
    let mut p_hi = f64::NEG_INFINITY;
    for pt in &res.pressure_curve {
        if pt.t.is_finite() {
            t_hi = t_hi.max(pt.t);
        }
        if pt.p.is_finite() {
            p_lo = p_lo.min(pt.p);
            p_hi = p_hi.max(pt.p);
        }
    }
    // Always include zero on the y-axis (the suction tail dips below).
    p_lo = p_lo.min(0.0);
    p_hi = p_hi.max(0.0);
    if !(t_hi.is_finite() && p_lo.is_finite() && p_hi.is_finite()) || t_hi <= 0.0 {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "non-finite pressure curve",
            egui::FontId::monospace(12.0),
            egui::Color32::from_gray(120),
        );
        return;
    }
    if (p_hi - p_lo).abs() < 1e-12 {
        p_hi += 1.0;
    }

    let margin = 30.0_f32;
    let inner = rect.shrink(margin);
    let to_px = |t: f64, p: f64| -> egui::Pos2 {
        let fx = ((t - 0.0) / (t_hi - 0.0)).clamp(0.0, 1.0) as f32;
        let fy = ((p - p_lo) / (p_hi - p_lo)).clamp(0.0, 1.0) as f32;
        egui::pos2(
            inner.left() + fx * inner.width(),
            inner.bottom() - fy * inner.height(),
        )
    };

    // Axes + zero line.
    let axis = egui::Stroke::new(1.0, egui::Color32::from_gray(70));
    painter.line_segment(
        [
            egui::pos2(inner.left(), inner.bottom()),
            egui::pos2(inner.left(), inner.top()),
        ],
        axis,
    );
    let zero_y = to_px(0.0, 0.0).y;
    let mut x = inner.left();
    while x < inner.right() {
        painter.line_segment(
            [
                egui::pos2(x, zero_y),
                egui::pos2((x + 6.0).min(inner.right()), zero_y),
            ],
            egui::Stroke::new(1.0, egui::Color32::from_gray(60)),
        );
        x += 12.0;
    }

    // P(t) polyline.
    let pts: Vec<egui::Pos2> = res
        .pressure_curve
        .iter()
        .map(|pt| to_px(pt.t, pt.p))
        .collect();
    for w in pts.windows(2) {
        painter.line_segment(
            [w[0], w[1]],
            egui::Stroke::new(1.8, egui::Color32::from_rgb(90, 200, 220)),
        );
    }

    // Axis labels + the peak readout.
    painter.text(
        egui::pos2(inner.center().x, rect.bottom() - 2.0),
        egui::Align2::CENTER_BOTTOM,
        "time (ms)",
        egui::FontId::monospace(10.0),
        egui::Color32::from_gray(150),
    );
    painter.text(
        egui::pos2(rect.left() + 2.0, inner.top()),
        egui::Align2::LEFT_TOP,
        "P (Pa)",
        egui::FontId::monospace(10.0),
        egui::Color32::from_gray(150),
    );
    painter.text(
        egui::pos2(inner.right(), inner.top()),
        egui::Align2::RIGHT_TOP,
        format!(
            "Pso {:.1} kPa \u{00B7} i_s {:.0} Pa\u{00B7}s \u{00B7} t_d {:.2} ms",
            res.blast.peak_overpressure_pa / 1.0e3,
            res.blast.impulse_pa_s,
            res.blast.positive_duration_s * 1.0e3,
        ),
        egui::FontId::monospace(10.0),
        egui::Color32::from_rgb(150, 200, 210),
    );
}

/// View (b): the **pressure–impulse (P–I) diagram** — the iso-damage hyperbola
/// for the protective element with both asymptotes (impulsive = vertical,
/// quasi-static = horizontal) and the current `(Pso, i_s)` design point plotted
/// (green = safe, red = exceeds the protective limit). Drawn in log-log so the
/// hyperbola knee and both asymptotes are visible across decades.
fn draw_pi_diagram(res: &SurvivabilityResult, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Pressure-impulse (P-I) iso-damage diagram").strong());
    ui.label(
        egui::RichText::new(
            "white = iso-damage curve \u{00B7} cyan = impulsive (i_min) + quasi-static (P_min) \
             asymptotes \u{00B7} dot = design point (green safe / red exceeds) \u{00B7} log-log",
        )
        .weak()
        .small(),
    );

    let available = ui.available_size();
    let (rect, _) = ui.allocate_exact_size(
        egui::vec2(available.x.min(460.0), 240.0),
        egui::Sense::hover(),
    );
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 0.0, egui::Color32::from_rgb(16, 18, 26));

    let p_min = res.pi.pressure_asymptote_pa;
    let i_min = res.pi.impulse_asymptote_pa_s;
    if !(p_min.is_finite() && p_min > 0.0 && i_min.is_finite() && i_min > 0.0) {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "non-finite P-I asymptotes",
            egui::FontId::monospace(12.0),
            egui::Color32::from_gray(120),
        );
        return;
    }

    // Build the log-log bounds to span both asymptotes, the curve, and the
    // design point with a margin.
    let dp_p = res.blast.peak_overpressure_pa;
    let dp_i = res.blast.impulse_pa_s;
    let mut p_lo = p_min;
    let mut p_hi = p_min;
    let mut i_lo = i_min;
    let mut i_hi = i_min;
    let mut grow = |pp: f64, ii: f64| {
        if pp.is_finite() && pp > 0.0 {
            p_lo = p_lo.min(pp);
            p_hi = p_hi.max(pp);
        }
        if ii.is_finite() && ii > 0.0 {
            i_lo = i_lo.min(ii);
            i_hi = i_hi.max(ii);
        }
    };
    for pt in &res.pi_curve {
        grow(pt.pressure_pa, pt.impulse_pa_s);
    }
    grow(dp_p, dp_i);
    // Pad by a decade factor on each side so the asymptotes are readable.
    p_lo *= 0.3;
    p_hi *= 3.0;
    i_lo *= 0.3;
    i_hi *= 3.0;
    let (lpx0, lpx1) = (i_lo.ln(), i_hi.ln()); // x-axis = impulse
    let (lpy0, lpy1) = (p_lo.ln(), p_hi.ln()); // y-axis = pressure
    if !(lpx1 > lpx0 && lpy1 > lpy0) {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "degenerate P-I range",
            egui::FontId::monospace(12.0),
            egui::Color32::from_gray(120),
        );
        return;
    }

    let margin = 34.0_f32;
    let inner = rect.shrink(margin);
    let to_px = |impulse: f64, pressure: f64| -> egui::Pos2 {
        let fx = ((impulse.ln() - lpx0) / (lpx1 - lpx0)).clamp(0.0, 1.0) as f32;
        let fy = ((pressure.ln() - lpy0) / (lpy1 - lpy0)).clamp(0.0, 1.0) as f32;
        egui::pos2(
            inner.left() + fx * inner.width(),
            inner.bottom() - fy * inner.height(),
        )
    };

    // Axes.
    let axis = egui::Stroke::new(1.0, egui::Color32::from_gray(70));
    painter.line_segment(
        [
            egui::pos2(inner.left(), inner.bottom()),
            egui::pos2(inner.right(), inner.bottom()),
        ],
        axis,
    );
    painter.line_segment(
        [
            egui::pos2(inner.left(), inner.top()),
            egui::pos2(inner.left(), inner.bottom()),
        ],
        axis,
    );

    // Asymptotes (cyan): vertical at i_min (impulsive), horizontal at P_min
    // (quasi-static).
    let cyan = egui::Color32::from_rgb(70, 170, 190);
    if (i_lo..=i_hi).contains(&i_min) {
        let x = to_px(i_min, (p_lo * p_hi).sqrt()).x;
        painter.line_segment(
            [egui::pos2(x, inner.top()), egui::pos2(x, inner.bottom())],
            egui::Stroke::new(1.2, cyan),
        );
    }
    if (p_lo..=p_hi).contains(&p_min) {
        let y = to_px((i_lo * i_hi).sqrt(), p_min).y;
        painter.line_segment(
            [egui::pos2(inner.left(), y), egui::pos2(inner.right(), y)],
            egui::Stroke::new(1.2, cyan),
        );
    }

    // The iso-damage curve (white), plotted as (impulse, pressure).
    let curve_pts: Vec<egui::Pos2> = res
        .pi_curve
        .iter()
        .filter(|pt| pt.pressure_pa > 0.0 && pt.impulse_pa_s > 0.0)
        .map(|pt| to_px(pt.impulse_pa_s, pt.pressure_pa))
        .collect();
    for w in curve_pts.windows(2) {
        painter.line_segment(
            [w[0], w[1]],
            egui::Stroke::new(1.8, egui::Color32::from_rgb(225, 225, 235)),
        );
    }

    // The design point (Pso, i_s): green if safe, red if it exceeds the limit.
    if dp_p.is_finite() && dp_p > 0.0 && dp_i.is_finite() && dp_i > 0.0 {
        let c = to_px(dp_i, dp_p);
        let col = if res.design_point_safe {
            egui::Color32::from_rgb(110, 210, 130)
        } else {
            egui::Color32::from_rgb(225, 90, 80)
        };
        painter.circle_filled(c, 4.5, col);
        painter.circle_stroke(c, 6.5, egui::Stroke::new(1.2, col));
    }

    // Axis labels.
    painter.text(
        egui::pos2(inner.center().x, rect.bottom() - 2.0),
        egui::Align2::CENTER_BOTTOM,
        "impulse i (Pa\u{00B7}s, log)",
        egui::FontId::monospace(10.0),
        egui::Color32::from_gray(150),
    );
    painter.text(
        egui::pos2(rect.left() + 2.0, inner.center().y),
        egui::Align2::LEFT_CENTER,
        "P (Pa, log)",
        egui::FontId::monospace(10.0),
        egui::Color32::from_gray(150),
    );
}

/// View (c): the **armor-sizing readout** — the threat kinetic energy and the
/// minimum protective plate thickness + areal density that just defeats it.
fn draw_armor_readout(res: &SurvivabilityResult, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Minimum armor protection").strong());
    ui.label(
        egui::RichText::new(
            "minimum plate to keep the threat OUT (plug-shear energy balance) \u{2014} \
             no penetration is modeled",
        )
        .weak()
        .small(),
    );
    egui::Grid::new("survivability_armor_readout")
        .num_columns(2)
        .striped(true)
        .show(ui, |ui| {
            let row = |ui: &mut egui::Ui, k: &str, v: String| {
                ui.label(k);
                ui.label(egui::RichText::new(v).monospace());
                ui.end_row();
            };
            row(
                ui,
                "threat kinetic energy",
                format!("{:.0} J", res.armor.threat_energy_j),
            );
            row(
                ui,
                "min plate thickness",
                format!(
                    "{:.2} mm ({:.4} m)",
                    res.armor.min_thickness_m * 1.0e3,
                    res.armor.min_thickness_m
                ),
            );
            row(
                ui,
                "areal density",
                format!("{:.1} kg/m\u{00B2}", res.armor.areal_density_kg_m2),
            );
        });
}

/// View (d): the **occupant-margin readout** — peak g, the supplied whole-body
/// tolerance, the margin, and the survivable verdict, plus the SDOF response
/// summary that feeds it.
fn draw_occupant_readout(res: &SurvivabilityResult, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Protective response & occupant margin").strong());
    egui::Grid::new("survivability_occupant_readout")
        .num_columns(2)
        .striped(true)
        .show(ui, |ui| {
            let row = |ui: &mut egui::Ui, k: &str, v: String| {
                ui.label(k);
                ui.label(egui::RichText::new(v).monospace());
                ui.end_row();
            };
            row(
                ui,
                "peak deflection",
                format!("{:.4} m", res.response.peak_deflection_m),
            );
            row(
                ui,
                "ductility ratio \u{03BC}",
                res.response
                    .ductility_ratio
                    .map(|m| format!("{m:.3}"))
                    .unwrap_or_else(|| "\u{2014}".to_string()),
            );
            row(
                ui,
                "dynamic amplification",
                format!("{:.3}", res.response.dynamic_amplification),
            );
            row(
                ui,
                "peak acceleration",
                format!("{:.1} m/s\u{00B2}", res.response.peak_acceleration_m_s2),
            );
            row(ui, "peak g", format!("{:.2} g", res.occupant.peak_g));
            row(
                ui,
                "tolerance",
                format!("{:.1} g", res.occupant.tolerance_g),
            );
            row(
                ui,
                "margin (g_tol/g_peak)",
                format!("{:.2}", res.occupant.margin),
            );
            ui.label("verdict");
            let (txt, col) = if res.occupant.survivable {
                ("within tolerance", egui::Color32::from_rgb(110, 210, 130))
            } else {
                ("EXCEEDS tolerance", egui::Color32::from_rgb(225, 90, 80))
            };
            ui.label(egui::RichText::new(txt).color(col).strong());
            ui.end_row();
        });
}

// ---------------------------------------------------------------------------
// Tests (unit + headless_ui_tests, mirroring missionsim_workbench)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_run_succeeds_and_is_populated() {
        let s = SurvivabilityWorkbenchState::default();
        let res = s.run().expect("default survivability run should succeed");
        // Blast load is physical.
        assert!(res.blast.peak_overpressure_pa > 0.0);
        assert!(res.blast.impulse_pa_s > 0.0);
        // Pressure curve has samples + 1 points and starts at the peak.
        assert!(res.pressure_curve.len() >= 2);
        // SDOF response is physical.
        assert!(res.response.peak_deflection_m > 0.0);
        assert!(res.response.dynamic_amplification > 0.0);
        assert!(res.response.ductility_ratio.unwrap() > 0.0);
        // P–I asymptotes are physical.
        assert!(res.pi.pressure_asymptote_pa > 0.0);
        assert!(res.pi.impulse_asymptote_pa_s > 0.0);
        assert!(res.pi_curve.len() >= 2);
        // Armor sizing is physical.
        assert!(res.armor.min_thickness_m > 0.0);
        assert!(res.armor.areal_density_kg_m2 > 0.0);
        // Occupant screen produced a finite peak g.
        assert!(res.occupant.peak_g.is_finite());
    }

    #[test]
    fn pi_asymptotes_match_closed_form() {
        // P_min = k·x_lim, i_min = x_lim·√(km), per the crate's pinned form.
        let s = SurvivabilityWorkbenchState::default();
        let res = s.run().expect("run");
        let p = &s.params;
        let expect_p = p.element_stiffness_n_m * p.yield_deflection_m;
        let expect_i = p.yield_deflection_m * (p.element_stiffness_n_m * p.element_mass_kg).sqrt();
        assert!((res.pi.pressure_asymptote_pa - expect_p).abs() / expect_p < 1e-9);
        assert!((res.pi.impulse_asymptote_pa_s - expect_i).abs() / expect_i < 1e-9);
    }

    #[test]
    fn larger_standoff_lowers_overpressure() {
        // PIN: more stand-off ⇒ smaller design load (monotone).
        let mut near = SurvivabilityWorkbenchState::default();
        near.params.standoff_m = 15.0;
        let mut far = SurvivabilityWorkbenchState::default();
        far.params.standoff_m = 35.0;
        let pn = near.run().expect("near").blast.peak_overpressure_pa;
        let pf = far.run().expect("far").blast.peak_overpressure_pa;
        assert!(
            pf < pn,
            "overpressure should drop with stand-off: {pf} !< {pn}"
        );
    }

    #[test]
    fn faster_projectile_needs_more_armor() {
        // PIN: more threat KE ⇒ more minimum protection.
        let mut slow = SurvivabilityWorkbenchState::default();
        slow.params.projectile_velocity_m_s = 400.0;
        let mut fast = SurvivabilityWorkbenchState::default();
        fast.params.projectile_velocity_m_s = 1000.0;
        let ts = slow.run().expect("slow").armor.min_thickness_m;
        let tf = fast.run().expect("fast").armor.min_thickness_m;
        assert!(tf > ts, "faster threat should need a thicker plate");
    }

    #[test]
    fn occupant_over_tolerance_flags_unsurvivable() {
        // A tiny tolerance forces an over-tolerance verdict (still no panic).
        let mut s = SurvivabilityWorkbenchState::default();
        s.params.occupant_tolerance_g = 1e-6;
        let res = s.run().expect("run");
        assert!(!res.occupant.survivable);
        assert!(res.occupant.margin < 1.0);
    }

    // ---- degenerate-param tests — must return Err, NOT panic ----

    #[test]
    fn zero_charge_returns_err() {
        let mut s = SurvivabilityWorkbenchState::default();
        s.params.charge_kg = 0.0;
        assert!(s.run().is_err(), "zero charge must return Err, not panic");
    }

    #[test]
    fn zero_standoff_returns_err() {
        let mut s = SurvivabilityWorkbenchState::default();
        s.params.standoff_m = 0.0;
        assert!(
            s.run().is_err(),
            "zero stand-off must return Err, not panic"
        );
    }

    #[test]
    fn scaled_distance_out_of_band_returns_err() {
        // A huge charge at tiny stand-off drives Z below Z_MIN.
        let mut s = SurvivabilityWorkbenchState::default();
        s.params.charge_kg = 1.0e6;
        s.params.standoff_m = 0.5;
        assert!(
            s.run().is_err(),
            "a scaled distance outside the validated band must return Err"
        );
    }

    #[test]
    fn zero_element_mass_returns_err() {
        let mut s = SurvivabilityWorkbenchState::default();
        s.params.element_mass_kg = 0.0;
        assert!(s.run().is_err(), "zero element mass must return Err");
    }

    #[test]
    fn zero_stiffness_returns_err() {
        let mut s = SurvivabilityWorkbenchState::default();
        s.params.element_stiffness_n_m = 0.0;
        assert!(s.run().is_err(), "zero stiffness must return Err");
    }

    #[test]
    fn negative_damping_returns_err() {
        let mut s = SurvivabilityWorkbenchState::default();
        s.params.element_damping_n_s_m = -1.0;
        assert!(s.run().is_err(), "negative damping must return Err");
    }

    #[test]
    fn zero_yield_deflection_returns_err() {
        let mut s = SurvivabilityWorkbenchState::default();
        s.params.yield_deflection_m = 0.0;
        assert!(s.run().is_err(), "zero yield deflection must return Err");
    }

    #[test]
    fn zero_loaded_area_returns_err() {
        let mut s = SurvivabilityWorkbenchState::default();
        s.params.loaded_area_m2 = 0.0;
        assert!(s.run().is_err(), "zero loaded area must return Err");
    }

    #[test]
    fn zero_projectile_diameter_returns_err() {
        let mut s = SurvivabilityWorkbenchState::default();
        s.params.projectile_diameter_m = 0.0;
        assert!(s.run().is_err(), "zero projectile diameter must return Err");
    }

    #[test]
    fn zero_plate_strength_returns_err() {
        let mut s = SurvivabilityWorkbenchState::default();
        s.params.plate_shear_strength_pa = 0.0;
        assert!(
            s.run().is_err(),
            "zero plate shear strength must return Err"
        );
    }

    #[test]
    fn zero_tolerance_returns_err() {
        let mut s = SurvivabilityWorkbenchState::default();
        s.params.occupant_tolerance_g = 0.0;
        assert!(s.run().is_err(), "zero occupant tolerance must return Err");
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
            draw_survivability_workbench(app, ctx);
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
        assert!(!app.show_survivability_workbench);
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_survivability_workbench(&mut app, ctx);
        });
        // No panic = pass.
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_survivability_workbench = true;
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn workbench_draws_with_populated_result_without_panic() {
        let mut app = ValenxApp::default();
        app.show_survivability_workbench = true;
        let res = app.survivability.run().expect("run should succeed");
        app.survivability.result = Some(res);
        app.survivability.status = "\u{2714} test result".to_string();
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn workbench_draws_with_error_status_without_panic() {
        let mut app = ValenxApp::default();
        app.show_survivability_workbench = true;
        // Trigger an error state (zero charge is fail-loud in run()).
        app.survivability.params.charge_kg = 0.0;
        let result = app.survivability.run();
        app.survivability.status = match result {
            Err(e) => format!("\u{26A0} {e}"),
            Ok(_) => "\u{26A0} simulated error for testing".to_string(),
        };
        app.survivability.result = None;
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // Every numeric DragValue is a SpinButton and must be `labelled_by` its
        // caption (egui clears a DragValue's own Name), so an AI / screen reader
        // can find the control by caption text. Each `labelled_by` target must
        // RESOLVE to a real named caption node, not a dangling id.
        let mut app = ValenxApp::default();
        app.show_survivability_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let by_id: std::collections::HashMap<NodeId, &Node> =
            nodes.iter().map(|(id, n)| (*id, n)).collect();

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        // Many numeric controls across blast / element / armor / occupant; a
        // conservative lower bound that all are present and named.
        assert!(
            spin_buttons.len() >= 13,
            "expected many numeric controls as spin buttons, got {}",
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

        // Representative captions present as named accessibility nodes.
        for caption in [
            "charge mass W (kg)",
            "stand-off R (m)",
            "element mass (kg)",
            "element stiffness (N/m)",
            "yield deflection (m)",
            "loaded area (m\u{00B2})",
            "projectile mass (kg)",
            "projectile velocity (m/s)",
            "projectile diameter (m)",
            "plate shear strength (Pa)",
            "plate density (kg/m\u{00B3})",
            "occupant tolerance (g)",
        ] {
            assert!(
                has_named_node(&nodes, caption),
                "caption '{caption}' must be a named node in the a11y tree"
            );
        }

        // The Run button must be a named, invokable node.
        assert!(
            nodes.iter().any(|(_, n)| {
                n.role() == Role::Button && n.name().is_some_and(|s| s.contains("Run"))
            }),
            "the Run button must be a named, invokable node"
        );
    }

    #[test]
    fn degenerate_params_show_error_not_panic() {
        // When the charge is 0 or the scaled distance is out of band the
        // workbench must surface the error in-panel, not panic.
        let mut state = SurvivabilityWorkbenchState::default();
        state.params.charge_kg = 0.0;
        assert!(
            state.run().is_err(),
            "zero charge must produce Err, not panic"
        );
        state.params.charge_kg = 1.0e6;
        state.params.standoff_m = 0.5;
        assert!(
            state.run().is_err(),
            "an out-of-band scaled distance must produce Err, not panic"
        );
    }

    #[test]
    fn agent_bridge_survivability_id_resolves_and_sets_flag() {
        // Verify the two mechanisms the agent bridge uses for
        //   `OpenWorkbench { id: "survivability" }` (and the aliases):
        //   1. TabKind::from_id("survivability") -> Some(TabKind::Survivability)
        //   2. set_workbench_flag(app, "survivability", true) -> flag set
        use crate::project_tabs::{set_workbench_flag, TabKind};

        // 1. Lookup + aliases + case/whitespace tolerance.
        assert_eq!(
            TabKind::from_id("survivability"),
            Some(TabKind::Survivability),
            "\"survivability\" must resolve to TabKind::Survivability"
        );
        assert_eq!(
            TabKind::from_id("SURVIVABILITY"),
            Some(TabKind::Survivability)
        );
        assert_eq!(
            TabKind::from_id("  survivability  "),
            Some(TabKind::Survivability)
        );
        assert_eq!(TabKind::from_id("protection"), Some(TabKind::Survivability));
        assert_eq!(TabKind::from_id("blast"), Some(TabKind::Survivability));

        // 2. Flag toggle.
        let mut app = ValenxApp::default();
        assert!(!app.show_survivability_workbench);
        set_workbench_flag(&mut app, "survivability", true);
        assert!(
            app.show_survivability_workbench,
            "set_workbench_flag(\"survivability\", true) must set show_survivability_workbench"
        );
        set_workbench_flag(&mut app, "survivability", false);
        assert!(!app.show_survivability_workbench);
    }
}
