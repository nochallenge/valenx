//! Automated end-to-end vehicle design — the AI finds the best design with
//! no manual tuning.
//!
//! [`auto_design`] couples two optimizations so a single call produces a
//! finished design:
//!
//! 1. **Engine** — [`valenx_astro::optimize_engine`] searches chamber
//!    pressure × expansion ratio for the highest sea-level specific impulse
//!    whose first-order Bartz cooling margin clears a target. That optimized
//!    engine's `Isp` then powers the first stage.
//! 2. **Ascent** — a seeded search over payload × pitch-kick × vertical-rise
//!    flies each candidate to orbit ([`valenx_astro::simulate_ascent`]) and
//!    keeps the heaviest payload that reaches a bound orbit while the
//!    interstage thrust structure stays at safety factor ≥ 1.5.
//!
//! Honest scope: first-order preliminary design (ideal-nozzle engine, a
//! point-mass ascent, a closed-form strut check) — not multidisciplinary
//! design optimization with full CFD / FEA in the loop.

use valenx_astro::{
    constants::G0, optimize_engine, simulate_ascent, AscentConfig, CoolingInputs, DragModel,
    EngineDesign, EngineOptimum, GuidanceMode, GuidanceProgram, Stage, Vehicle, WindModel,
};

/// The AI-found best LV-1-class design — the optimized engine plus the
/// best-payload ascent that uses it.
#[derive(Debug, Clone)]
pub struct BestDesign {
    /// The optimized first-stage engine (and its performance + cooling).
    pub engine: EngineOptimum,
    /// Heaviest payload delivered to a bound orbit (kg).
    pub payload_kg: f64,
    /// Pitch-kick angle of the winning ascent (deg).
    pub pitch_kick_deg: f64,
    /// Vertical-rise time of the winning ascent (s).
    pub vertical_rise_s: f64,
    /// Apoapsis of the achieved orbit (km).
    pub apoapsis_km: f64,
    /// Periapsis of the achieved orbit (km).
    pub periapsis_km: f64,
    /// Peak axial g-load over the ascent.
    pub peak_g: f64,
    /// Max dynamic pressure (kPa).
    pub max_q_kpa: f64,
    /// Interstage structural safety factor at that peak g-load.
    pub structural_sf: f64,
    /// Engine candidates evaluated.
    pub engine_evals: usize,
    /// Ascent candidates evaluated.
    pub ascent_evals: usize,
}

/// A kerolox gas-generator-class engine design point (the optimizer's seed).
fn kerolox_base() -> EngineDesign {
    EngineDesign {
        chamber_pressure: 9.7e6,
        chamber_temperature: 3_500.0,
        gamma: 1.2,
        molar_mass: 22.0,
        throat_area: 0.05,
        expansion_ratio: 16.0,
    }
}

/// Run the full automated design search and return the best design found, or
/// `None` if nothing in range reaches orbit within the structural limit.
///
/// `ascent_evals` is the number of seeded ascent candidates to fly (the
/// engine grid is fixed); a few hundred is plenty since each sim is sub-ms.
/// Deterministic (seeded) so the answer is reproducible.
pub fn auto_design(ascent_evals: usize) -> Option<BestDesign> {
    // ── 1. AI-optimize the engine ────────────────────────────────────────
    let opt = optimize_engine(
        &kerolox_base(),
        &CoolingInputs::default(),
        1.5,
        (4.0e6, 20.0e6),
        (5.0, 60.0),
        28,
    )?;
    // The optimized engine's Isp powers the first stage; thrust is held at the
    // LV-1's nominal ~1.5 MN class so the airframe stays sane.
    let vac = opt.vacuum;
    let sl = opt.sea_level;

    // Interstage capacity: 8 Al-2024-T3 struts of 15 cm² each.
    let capacity_n = 324.0e6 * 8.0 * 1.5e-3;
    let (s1_dry, s1_prop) = (6_000.0, 90_000.0);
    let (s2_dry, s2_prop) = (1_500.0, 12_000.0);

    // ── 2. AI-optimize the ascent on the engine-powered vehicle ──────────
    // A tiny deterministic LCG (no rng dependency).
    let mut seed: u64 = 0x9E37_79B9_7F4A_7C15;
    let mut unit = || {
        seed = seed
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (seed >> 33) as f64 / ((1u64 << 31) as f64)
    };

    let mut best: Option<BestDesign> = None;
    let mut best_payload = 0.0_f64;

    for _ in 0..ascent_evals.max(1) {
        let payload = 500.0 + unit() * 7_500.0;
        let pitch = 8.0 + unit() * 10.0;
        let rise = 14.0 + unit() * 31.0;

        let stage1 = Stage {
            name: "LV-1 stage 1 (AI-optimized engine)".into(),
            dry_mass: s1_dry,
            propellant_mass: s1_prop,
            thrust_vac: 1_650_000.0,
            thrust_sl: 1_500_000.0,
            isp_vac: vac.isp,
            isp_sl: sl.isp,
        };
        let stage2 = Stage {
            name: "LV-1 stage 2".into(),
            dry_mass: s2_dry,
            propellant_mass: s2_prop,
            thrust_vac: 180_000.0,
            thrust_sl: 180_000.0,
            isp_vac: 345.0,
            isp_sl: 345.0,
        };
        let vehicle = Vehicle {
            stages: vec![stage1, stage2],
            payload_mass: payload,
            reference_area: 2.5,
            drag: DragModel::generic_launch_vehicle(),
        };
        let config = AscentConfig {
            launch_altitude_m: 0.0,
            guidance: GuidanceProgram {
                vertical_rise_time: rise,
                pitch_kick_deg: pitch,
                kick_duration: 5.0,
            },
            time_step: 0.1,
            max_time: 900.0,
            sample_interval: 5.0,
            mode: GuidanceMode::OpenLoopGravityTurn,
            wind: WindModel::None,
        };

        if let Ok(r) = simulate_ascent(&vehicle, &config) {
            if r.reached_orbit && r.periapsis_km() > 100.0 {
                let load = (s2_dry + s2_prop + payload) * r.max_acceleration_g * G0;
                let sf = if load > 0.0 {
                    capacity_n / load
                } else {
                    f64::INFINITY
                };
                if sf >= 1.5 && payload > best_payload {
                    best_payload = payload;
                    best = Some(BestDesign {
                        engine: opt,
                        payload_kg: payload,
                        pitch_kick_deg: pitch,
                        vertical_rise_s: rise,
                        apoapsis_km: r.apoapsis_km(),
                        periapsis_km: r.periapsis_km(),
                        peak_g: r.max_acceleration_g,
                        max_q_kpa: r.max_dynamic_pressure / 1000.0,
                        structural_sf: sf,
                        engine_evals: opt.evaluations,
                        ascent_evals,
                    });
                }
            }
        }
    }

    best
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn auto_design_finds_a_complete_feasible_design() {
        let d = auto_design(400).expect("a feasible design is found");
        // Engine: cooled within the target margin, sensible kerolox Isp.
        assert!(d.engine.cooling.cooling_margin >= 1.5 - 1e-9);
        assert!((250.0..380.0).contains(&d.engine.vacuum.isp));
        // Ascent: a real bound orbit, a positive payload, structurally sound.
        assert!(d.payload_kg > 500.0 && d.payload_kg <= 8000.0);
        assert!(d.periapsis_km > 100.0, "periapsis {}", d.periapsis_km);
        assert!(d.structural_sf >= 1.5, "SF {}", d.structural_sf);
        assert!(d.peak_g > 0.0 && d.peak_g < 20.0);
    }

    #[test]
    fn auto_design_is_deterministic() {
        let a = auto_design(200).expect("feasible");
        let b = auto_design(200).expect("feasible");
        assert_eq!(a.payload_kg as u64, b.payload_kg as u64);
        assert_eq!(
            a.engine.design.chamber_pressure as u64,
            b.engine.design.chamber_pressure as u64
        );
    }

    /// Prints the AI-found best design — run with
    /// `cargo test -p valenx-rocket-demo --lib auto_design::tests::dump -- --ignored --nocapture`.
    #[test]
    #[ignore = "prints the auto-designed rocket"]
    fn dump_best_design() {
        let d = auto_design(800).expect("feasible");
        println!("\n=== AI-AUTO-DESIGNED LV-1 ===");
        println!(
            "engine : Pc {:.0} bar · ε {:.1} → SL Isp {:.0} s · vac Isp {:.0} s",
            d.engine.design.chamber_pressure / 1.0e5,
            d.engine.design.expansion_ratio,
            d.engine.sea_level.isp,
            d.engine.vacuum.isp,
        );
        println!(
            "cooling: throat flux {:.1} MW/m² · margin {:.2}",
            d.engine.cooling.throat_heat_flux / 1.0e6,
            d.engine.cooling.cooling_margin,
        );
        println!(
            "ascent : {:.0} kg → {:.0} × {:.0} km · max-Q {:.0} kPa · peak {:.1} g",
            d.payload_kg, d.periapsis_km, d.apoapsis_km, d.max_q_kpa, d.peak_g,
        );
        println!(
            "guidance: pitch {:.1}° · vertical rise {:.0} s · interstage SF {:.2}",
            d.pitch_kick_deg, d.vertical_rise_s, d.structural_sf,
        );
        println!(
            "searched: {} engine candidates × {} ascent sims\n",
            d.engine_evals, d.ascent_evals,
        );
    }
}
