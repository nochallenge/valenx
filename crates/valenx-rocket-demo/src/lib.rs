//! # valenx-rocket-demo
//!
//! A worked, end-to-end **rocket design → simulate** example that ties two
//! valenx engines into one pipeline — the concrete answer to "can I design a
//! rocket *and* test it under simulation, all in valenx?":
//!
//! - **Trajectory** — [`valenx_astro`] flies a multi-stage vehicle from the pad
//!   to orbit (US-Standard-Atmosphere-1976, staged thrust/Isp with mass flow,
//!   Mach-dependent drag, gravity-turn guidance), reporting the orbit reached,
//!   the Δv budget, max-Q and the peak axial g-load.
//! - **Structure** — [`valenx_fem::beam`] sizes a load-bearing member against
//!   that flight: the **interstage thrust structure** must carry the upper
//!   stage's weight at peak acceleration without yielding.
//!
//! The two legs are **coupled**: the structural load is not assumed — it is
//! *derived from the trajectory result* via `F = m · a_max` (with
//! `a_max = g_max · g₀`). That makes this a genuine preliminary
//! design-and-check loop rather than two unrelated calculators.
//!
//! ## Honest scope
//!
//! Both engines are **research / preliminary-design grade** (see their own
//! module docs): the ascent is a planar 3-DOF point-mass flight, and the
//! structural check is a closed-form axial-member sizing (the
//! [`valenx_fem`] beam library — which also offers a full 3-D beam-element FE
//! solver for higher fidelity). This demonstrates that the design +
//! multi-physics-analysis *workflow* lives in one place; it is **not** a
//! flight-certified design tool.

use valenx_astro::constants::G0;
use valenx_astro::{
    presets, simulate_ascent, AscentConfig, AstroError, DragModel, GuidanceMode, GuidanceProgram,
    Stage, Vehicle, WindModel,
};
use valenx_fem::beam::{axial_force_capacity, axial_stress};

/// A minimal rocket design. The trajectory flies the medium-lift two-stage
/// preset; these parameters size the **interstage thrust structure** — the
/// ring of axial struts that carries the upper stage during powered flight.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RocketDesign {
    /// Mass carried by the interstage (upper stage + payload), kg.
    pub supported_mass_kg: f64,
    /// Number of axial load-bearing struts in the interstage.
    pub strut_count: usize,
    /// Cross-sectional area of each strut, m².
    pub strut_area_m2: f64,
    /// Yield strength of the strut material, Pa (≈ 324 MPa for Al-2024-T3).
    pub material_yield_pa: f64,
}

impl Default for RocketDesign {
    fn default() -> Self {
        // A plausible LEO upper stage on eight 10 cm² aluminium struts.
        RocketDesign {
            supported_mass_kg: 20_000.0,
            strut_count: 8,
            strut_area_m2: 1.0e-3,
            material_yield_pa: 324.0e6,
        }
    }
}

/// The combined result of the trajectory + structural analysis.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RocketReport {
    /// Did the vehicle reach a bound orbit with periapsis above the atmosphere?
    pub reached_orbit: bool,
    /// Apoapsis altitude (km).
    pub apoapsis_km: f64,
    /// Periapsis altitude (km).
    pub periapsis_km: f64,
    /// Ideal Δv budget (m/s).
    pub delta_v_budget_ms: f64,
    /// Peak dynamic pressure, max-Q (Pa).
    pub max_q_pa: f64,
    /// Peak axial acceleration during flight (g).
    pub max_acceleration_g: f64,
    /// Peak axial load the interstage carries, `F = m · g_max · g₀` (N) —
    /// derived from the trajectory's peak g-load.
    pub peak_axial_load_n: f64,
    /// Axial stress in each strut, `σ = F / (N · A)` (Pa).
    pub strut_stress_pa: f64,
    /// Structural safety factor against yield, `σ_yield / σ`.
    pub structural_safety_factor: f64,
    /// Whether the structure survives the peak load (safety factor ≥ 1).
    pub structurally_safe: bool,
}

/// Couple an ascent — `vehicle` flown under `config` — with an interstage
/// **strut check**: the trajectory's peak g-load becomes the inertial load
/// the struts in `design` must carry (`F = m · a_max`), yielding the
/// per-strut stress `σ = F / (N · A)` and the safety factor `σ_yield / σ`.
/// Shared by [`design_and_simulate`] (the medium-lift preset) and
/// [`fly_valenx_lv1`] (the from-scratch LV-1 launcher).
///
/// # Errors
///
/// Propagates [`AstroError`] if the trajectory integration rejects the
/// vehicle / configuration.
pub fn simulate_coupled(
    vehicle: &Vehicle,
    config: &AscentConfig,
    design: &RocketDesign,
) -> Result<RocketReport, AstroError> {
    // --- Trajectory leg (valenx-astro) ---
    let flight = simulate_ascent(vehicle, config)?;

    // --- Structural leg (valenx-fem), loaded BY the trajectory result ---
    // Peak axial acceleration → inertial load the interstage must react.
    let a_max = flight.max_acceleration_g * G0; // m/s²
    let peak_axial_load_n = design.supported_mass_kg * a_max; // F = m·a
    let total_area = design.strut_count.max(1) as f64 * design.strut_area_m2;
    // Per-strut axial stress σ = F/A (the struts share the load equally).
    let strut_stress_pa = axial_stress(peak_axial_load_n, total_area);
    // Capacity-based safety factor — an independent valenx-fem::beam path: the
    // struts can carry σ_yield·A_total before yielding. (Equals σ_yield/σ.)
    let capacity_n = axial_force_capacity(design.material_yield_pa, total_area);
    let structural_safety_factor = if peak_axial_load_n > 0.0 {
        capacity_n / peak_axial_load_n
    } else {
        f64::INFINITY
    };

    Ok(RocketReport {
        reached_orbit: flight.reached_orbit,
        apoapsis_km: flight.apoapsis_km(),
        periapsis_km: flight.periapsis_km(),
        delta_v_budget_ms: flight.ideal_delta_v,
        max_q_pa: flight.max_dynamic_pressure,
        max_acceleration_g: flight.max_acceleration_g,
        peak_axial_load_n,
        strut_stress_pa,
        structural_safety_factor,
        structurally_safe: structural_safety_factor >= 1.0,
    })
}

/// Run the end-to-end **design → simulate** pipeline for `design` on the
/// generic medium-lift two-stage preset (see [`simulate_coupled`]).
///
/// # Errors
///
/// Propagates [`AstroError`] if the trajectory integration rejects the inputs.
pub fn design_and_simulate(design: &RocketDesign) -> Result<RocketReport, AstroError> {
    simulate_coupled(
        &presets::two_stage_medium_lift(),
        &presets::leo_ascent_config(),
        design,
    )
}

/// **Valenx LV-1** — a from-scratch two-stage kerolox small-lift launch
/// vehicle, designed in valenx to lift ~2 t to low Earth orbit. The stage
/// masses are sized to a Tsiolkovsky `Δv` budget near 9.9 km/s (LEO orbital
/// velocity ~7.8 km/s plus ~2 km/s of gravity + drag + steering losses),
/// with a ~1.37 liftoff thrust-to-weight. Illustrative / research-grade —
/// not flight data for any real vehicle.
pub fn valenx_lv1() -> Vehicle {
    Vehicle {
        stages: vec![
            Stage {
                name: "LV-1 first stage (kerolox)".into(),
                dry_mass: 6_000.0,
                propellant_mass: 90_000.0,
                thrust_sl: 1_500_000.0,
                thrust_vac: 1_650_000.0,
                isp_sl: 283.0,
                isp_vac: 311.0,
            },
            Stage {
                name: "LV-1 second stage (kerolox, vacuum)".into(),
                dry_mass: 1_500.0,
                propellant_mass: 12_000.0,
                thrust_sl: 180_000.0,
                thrust_vac: 180_000.0,
                isp_sl: 345.0,
                isp_vac: 345.0,
            },
        ],
        payload_mass: 2_000.0,
        // ~1.8 m diameter -> ~2.5 m² frontal area.
        reference_area: 2.5,
        drag: DragModel::generic_launch_vehicle(),
    }
}

/// The Valenx LV-1's ascent profile — an open-loop gravity turn tuned to
/// fly the vehicle into a bound orbit.
pub fn valenx_lv1_ascent() -> AscentConfig {
    AscentConfig {
        launch_altitude_m: 0.0,
        guidance: GuidanceProgram {
            vertical_rise_time: 20.0,
            pitch_kick_deg: 12.0,
            kick_duration: 5.0,
        },
        time_step: 0.1,
        max_time: 1_800.0,
        sample_interval: 2.0,
        mode: GuidanceMode::OpenLoopGravityTurn,
        wind: WindModel::None,
    }
}

/// The Valenx LV-1's interstage thrust structure: 8 aluminium (Al-2024-T3)
/// struts carrying the wet upper stage + payload during first-stage burn.
pub fn valenx_lv1_interstage() -> RocketDesign {
    let v = valenx_lv1();
    let upper = &v.stages[1];
    RocketDesign {
        supported_mass_kg: upper.dry_mass + upper.propellant_mass + v.payload_mass,
        strut_count: 8,
        strut_area_m2: 1.5e-3, // 15 cm² each
        material_yield_pa: 324.0e6,
    }
}

/// Fly the **Valenx LV-1** end to end: design → ascent → structural check.
///
/// # Errors
///
/// Propagates [`AstroError`] if the trajectory integration rejects the design.
pub fn fly_valenx_lv1() -> Result<RocketReport, AstroError> {
    simulate_coupled(
        &valenx_lv1(),
        &valenx_lv1_ascent(),
        &valenx_lv1_interstage(),
    )
}

impl std::fmt::Display for RocketReport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "── valenx rocket design → simulate ──")?;
        writeln!(f, "Trajectory (valenx-astro):")?;
        writeln!(f, "  reached orbit:   {}", self.reached_orbit)?;
        writeln!(f, "  apoapsis:        {:.0} km", self.apoapsis_km)?;
        writeln!(f, "  periapsis:       {:.0} km", self.periapsis_km)?;
        writeln!(f, "  Δv budget:       {:.0} m/s", self.delta_v_budget_ms)?;
        writeln!(f, "  max-Q:           {:.1} kPa", self.max_q_pa / 1000.0)?;
        writeln!(f, "  peak g-load:     {:.1} g", self.max_acceleration_g)?;
        writeln!(f, "Structure (valenx-fem, loaded by peak g):")?;
        writeln!(
            f,
            "  axial load:      {:.0} kN",
            self.peak_axial_load_n / 1000.0
        )?;
        writeln!(
            f,
            "  strut stress:    {:.0} MPa",
            self.strut_stress_pa / 1.0e6
        )?;
        writeln!(f, "  safety factor:   {:.2}", self.structural_safety_factor)?;
        write!(
            f,
            "  verdict:         {}",
            if self.structurally_safe {
                "SAFE"
            } else {
                "OVER-STRESSED"
            }
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_design_reaches_orbit_and_is_structurally_safe() {
        let r = design_and_simulate(&RocketDesign::default()).unwrap();
        // Trajectory: the medium-lift preset reaches a bound LEO (astro-validated).
        assert!(r.reached_orbit, "should reach orbit");
        assert!(r.periapsis_km > 100.0, "periapsis {} km", r.periapsis_km);
        assert!(r.delta_v_budget_ms > 9_400.0, "Δv {}", r.delta_v_budget_ms);
        assert!(r.max_q_pa > 0.0 && r.max_acceleration_g > 0.0);
        // Structure: a reasonable interstage on 8×10 cm² Al struts survives.
        assert!(r.structurally_safe, "SF {}", r.structural_safety_factor);
        assert!(r.structural_safety_factor > 1.0);
    }

    #[test]
    fn structural_load_is_coupled_to_the_trajectory_g_load() {
        // GROUND TRUTH: the axial load is derived from the SIM's peak g
        // (F = m·a_max), and the per-strut stress + safety factor follow the
        // closed forms exactly — and the capacity-based safety factor
        // (valenx-fem axial_force_capacity) agrees with σ_yield/σ.
        let d = RocketDesign::default();
        let r = design_and_simulate(&d).unwrap();

        let expected_load = d.supported_mass_kg * r.max_acceleration_g * G0;
        assert!(
            (r.peak_axial_load_n - expected_load).abs() / expected_load < 1e-12,
            "F = m·a_max"
        );
        let total_area = d.strut_count as f64 * d.strut_area_m2;
        assert!(
            (r.strut_stress_pa - r.peak_axial_load_n / total_area).abs() / r.strut_stress_pa
                < 1e-12,
            "σ = F/(N·A)"
        );
        assert!(
            (r.structural_safety_factor - d.material_yield_pa / r.strut_stress_pa).abs()
                / r.structural_safety_factor
                < 1e-12,
            "SF = σ_yield/σ (two independent fem::beam paths agree)"
        );
    }

    #[test]
    fn undersized_struts_are_flagged_over_stressed() {
        // Shrink the struts 100× → stress 100× → safety factor below 1 → flagged.
        let weak = RocketDesign {
            strut_area_m2: 1.0e-5,
            ..RocketDesign::default()
        };
        let r = design_and_simulate(&weak).unwrap();
        assert!(
            !r.structurally_safe,
            "tiny struts should be over-stressed: SF {}",
            r.structural_safety_factor
        );
        assert!(r.structural_safety_factor < 1.0);
        // Doubling the strut count halves the stress / doubles the safety factor.
        let stronger = RocketDesign {
            strut_count: weak.strut_count * 2,
            ..weak
        };
        let r2 = design_and_simulate(&stronger).unwrap();
        assert!(
            (r2.structural_safety_factor / r.structural_safety_factor - 2.0).abs() < 1e-9,
            "2× struts → 2× safety factor"
        );
    }

    #[test]
    fn report_display_renders_both_legs() {
        let text = format!("{}", design_and_simulate(&RocketDesign::default()).unwrap());
        assert!(text.contains("Trajectory (valenx-astro)"));
        assert!(text.contains("Structure (valenx-fem"));
        assert!(text.contains("safety factor"));
    }
}

#[cfg(test)]
mod lv1_tests {
    use super::*;

    #[test]
    fn lv1_reaches_orbit_and_is_structurally_sound() {
        // The from-scratch LV-1 design flies to a bound orbit on its tuned
        // gravity-turn profile, and its 8×15 cm² aluminium interstage carries
        // the trajectory's peak g-load with margin.
        let r = fly_valenx_lv1().unwrap();
        assert!(
            r.reached_orbit,
            "LV-1 should reach orbit; periapsis {} km",
            r.periapsis_km
        );
        assert!(
            r.periapsis_km > 100.0,
            "periapsis {} km above the Kármán line",
            r.periapsis_km
        );
        assert!(
            r.structurally_safe,
            "interstage SF {}",
            r.structural_safety_factor
        );
        assert!(r.structural_safety_factor > 1.0);
        // The vehicle was sized to a ~9.9–10.1 km/s Tsiolkovsky LEO budget.
        assert!(
            (9_900.0..10_200.0).contains(&r.delta_v_budget_ms),
            "Δv budget {} m/s",
            r.delta_v_budget_ms
        );
    }

    #[test]
    fn lv1_interstage_supports_the_upper_stack() {
        // The interstage carries the wet upper stage + payload during the
        // first-stage burn — sized straight off the vehicle masses.
        let v = valenx_lv1();
        let upper = &v.stages[1];
        let d = valenx_lv1_interstage();
        assert!(
            (d.supported_mass_kg - (upper.dry_mass + upper.propellant_mass + v.payload_mass)).abs()
                < 1e-9
        );
    }

    #[test]
    #[ignore = "prints the LV-1 flight profile — run with --ignored --nocapture"]
    fn dump_lv1_flight() {
        let v = valenx_lv1();
        let r = simulate_ascent(&v, &valenx_lv1_ascent()).unwrap();
        println!(
            "SUMMARY apo={:.0}km peri={:.0}km dv={:.0} maxQ={:.0}Pa peakG={:.2} outcome={:?}",
            r.apoapsis_km(),
            r.periapsis_km(),
            r.ideal_delta_v,
            r.max_dynamic_pressure,
            r.max_acceleration_g,
            r.outcome
        );
        for e in &r.events {
            println!(
                "EVENT t={:.1} alt={:.0}m v={:.0} {}",
                e.time, e.altitude_m, e.speed, e.kind
            );
        }
        for s in r.samples.iter().step_by(4) {
            println!(
                "S {:.0} {:.2} {:.0} {:.0} {:.2}",
                s.time,
                s.altitude_m / 1000.0,
                s.speed_inertial,
                s.dynamic_pressure,
                s.acceleration_g
            );
        }
    }
}
