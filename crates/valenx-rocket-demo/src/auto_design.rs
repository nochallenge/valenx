//! Automated end-to-end vehicle design — the AI iterates through *thousands*
//! of virtual designs to converge on the best, with no manual tuning.
//!
//! [`auto_design`] runs a single coupled optimization over the whole design
//! vector — engine chamber pressure × expansion ratio, plus payload,
//! pitch-kick and vertical-rise — evaluating each candidate end to end:
//!
//! 1. build the engine, check its first-order Bartz cooling margin;
//! 2. power the first stage with that engine's `Isp`;
//! 3. fly the vehicle to orbit ([`valenx_astro::simulate_ascent`]);
//! 4. size the interstage thrust structure against the peak g-load.
//!
//! The search is **iterative** — an explore/exploit loop that perturbs the
//! best design so far with a shrinking radius — so over thousands of
//! evaluations it converges toward the heaviest payload that reaches a bound
//! orbit while staying cooled (margin ≥ 1.5) and structurally sound (SF ≥
//! 1.5). Deterministic (seeded) so the answer is reproducible.
//!
//! Honest scope: this is the same *iterate-through-thousands-of-designs*
//! methodology a computational-engineering model like LEAP 71's Noyron uses,
//! but at **preliminary-design fidelity** — ideal-nozzle thermodynamics, a
//! Bartz heat-flux estimate, and a point-mass ascent. It is NOT a
//! combustion / regen-cooling multiphysics solver, it does not generate
//! print-ready hardware, and nothing here is hot-fire validated.

use valenx_astro::{
    constants::G0, simulate_ascent, AscentConfig, CoolingInputs, CoolingPerformance, DragModel,
    EngineDesign, EnginePerformance, GuidanceMode, GuidanceProgram, Stage, Vehicle, WindModel,
};

/// The AI-found best LV-1-class design after iterating thousands of
/// candidates.
#[derive(Debug, Clone)]
pub struct BestDesign {
    /// The optimized first-stage engine design point.
    pub engine: EngineDesign,
    /// Its vacuum performance.
    pub engine_vacuum: EnginePerformance,
    /// Its sea-level performance.
    pub engine_sea_level: EnginePerformance,
    /// Its first-order cooling balance (margin ≥ 1.5).
    pub engine_cooling: CoolingPerformance,
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
    /// Total virtual designs evaluated.
    pub designs_evaluated: usize,
    /// How many were feasible (cooled, in orbit, structurally sound).
    pub feasible_count: usize,
    /// `[design_index, best_payload_kg_so_far]` convergence series.
    pub convergence: Vec<[f64; 2]>,
}

/// A kerolox gas-generator-class engine; the optimizer varies chamber
/// pressure and expansion ratio, the rest of the gas properties + throat are
/// held fixed.
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

/// Run the automated design search over `n_designs` virtual candidates and
/// return the best, or `None` if nothing in range was feasible.
pub fn auto_design(n_designs: usize) -> Option<BestDesign> {
    auto_design_with(n_designs, |_, _| {})
}

/// As [`auto_design`], but calls `on_progress(designs_done, best_payload_kg)`
/// each iteration so a background run can drive a live counter.
pub fn auto_design_with(
    n_designs: usize,
    mut on_progress: impl FnMut(usize, f64),
) -> Option<BestDesign> {
    let n = n_designs.max(1);
    let base = kerolox_base();
    let cooling = CoolingInputs::default();
    // Interstage capacity: 8 Al-2024-T3 struts of 15 cm² each.
    let capacity_n = 324.0e6 * 8.0 * 1.5e-3;
    let (s1_dry, s1_prop) = (6_000.0, 90_000.0);
    let (s2_dry, s2_prop) = (1_500.0, 12_000.0);

    // Design vector ranges: Pc, ε, payload, pitch, rise.
    let ranges = [
        (4.0e6, 20.0e6),
        (5.0, 60.0),
        (500.0, 8_000.0),
        (8.0, 18.0),
        (14.0, 45.0),
    ];
    let denorm = |v: f64, (lo, hi): (f64, f64)| lo + v * (hi - lo);

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
    let mut best_vars = [0.0_f64; 5];
    let mut feasible_count = 0usize;
    let mut convergence = Vec::with_capacity(n);

    for i in 0..n {
        // Explore (random) ~40 % of the time, otherwise exploit — perturb the
        // best-so-far with a radius that shrinks as the search progresses.
        let vars: [f64; 5] = if best.is_none() || unit() < 0.4 {
            [unit(), unit(), unit(), unit(), unit()]
        } else {
            let radius = 0.25 * (1.0 - i as f64 / n as f64).max(0.05);
            let mut v = best_vars;
            for slot in v.iter_mut() {
                *slot = (*slot + (unit() - 0.5) * 2.0 * radius).clamp(0.0, 1.0);
            }
            v
        };
        let pc = denorm(vars[0], ranges[0]);
        let eps = denorm(vars[1], ranges[1]);
        let payload = denorm(vars[2], ranges[2]);
        let pitch = denorm(vars[3], ranges[3]);
        let rise = denorm(vars[4], ranges[4]);

        // 1. Engine + cooling gate.
        let engine = EngineDesign {
            chamber_pressure: pc,
            expansion_ratio: eps,
            ..base
        };
        let (Ok(vac), Ok(sl), Ok(cool)) = (
            engine.vacuum(),
            engine.sea_level(),
            engine.cooling(&cooling),
        ) else {
            convergence.push([i as f64, best_payload]);
            on_progress(i + 1, best_payload);
            continue;
        };
        if cool.cooling_margin >= 1.5 {
            // 2-3. Fly the engine-powered vehicle.
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
                    // 4. Interstage structural check.
                    let load = (s2_dry + s2_prop + payload) * r.max_acceleration_g * G0;
                    let sf = if load > 0.0 {
                        capacity_n / load
                    } else {
                        f64::INFINITY
                    };
                    if sf >= 1.5 {
                        feasible_count += 1;
                        if payload > best_payload {
                            best_payload = payload;
                            best_vars = vars;
                            best = Some(BestDesign {
                                engine,
                                engine_vacuum: vac,
                                engine_sea_level: sl,
                                engine_cooling: cool,
                                payload_kg: payload,
                                pitch_kick_deg: pitch,
                                vertical_rise_s: rise,
                                apoapsis_km: r.apoapsis_km(),
                                periapsis_km: r.periapsis_km(),
                                peak_g: r.max_acceleration_g,
                                max_q_kpa: r.max_dynamic_pressure / 1000.0,
                                structural_sf: sf,
                                designs_evaluated: 0, // patched in below
                                feasible_count: 0,    // patched in below
                                convergence: Vec::new(),
                            });
                        }
                    }
                }
            }
        }
        convergence.push([i as f64, best_payload]);
        on_progress(i + 1, best_payload);
    }

    best.map(|mut b| {
        b.designs_evaluated = n;
        b.feasible_count = feasible_count;
        b.convergence = convergence;
        b
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn auto_design_finds_a_complete_feasible_design() {
        let d = auto_design(3_000).expect("a feasible design is found");
        assert_eq!(d.designs_evaluated, 3_000);
        assert!(d.feasible_count > 0);
        // Engine cooled within target, sensible kerolox Isp.
        assert!(d.engine_cooling.cooling_margin >= 1.5 - 1e-9);
        assert!((250.0..380.0).contains(&d.engine_vacuum.isp));
        // Ascent: a real bound orbit, structurally sound.
        assert!(d.payload_kg > 500.0 && d.payload_kg <= 8_000.0);
        assert!(d.periapsis_km > 100.0, "periapsis {}", d.periapsis_km);
        assert!(d.structural_sf >= 1.5, "SF {}", d.structural_sf);
        assert!(d.peak_g > 0.0 && d.peak_g < 20.0);
    }

    #[test]
    fn convergence_is_monotone_and_iterative() {
        let d = auto_design(2_000).expect("feasible");
        assert_eq!(d.convergence.len(), 2_000);
        // Best-so-far never decreases — the hallmark of an iterative search.
        for w in d.convergence.windows(2) {
            assert!(w[1][1] >= w[0][1] - 1e-9, "convergence must not regress");
        }
        // The exploit phase improves on the first feasible hit (the final best
        // is at least as good as the first non-zero convergence value).
        let first_hit = d.convergence.iter().find(|p| p[1] > 0.0).map(|p| p[1]);
        if let Some(first) = first_hit {
            assert!(d.payload_kg >= first - 1e-9);
        }
    }

    #[test]
    fn auto_design_is_deterministic() {
        let a = auto_design(1_000).expect("feasible");
        let b = auto_design(1_000).expect("feasible");
        assert_eq!(a.payload_kg as u64, b.payload_kg as u64);
        assert_eq!(
            a.engine.chamber_pressure as u64,
            b.engine.chamber_pressure as u64
        );
    }

    #[test]
    fn progress_callback_fires_each_design() {
        let mut count = 0usize;
        let mut last = 0usize;
        auto_design_with(500, |done, _best| {
            count += 1;
            last = done;
        });
        assert_eq!(count, 500);
        assert_eq!(last, 500);
    }

    /// Prints the AI-found best design — run with
    /// `cargo test -p valenx-rocket-demo --lib auto_design::tests::dump -- --ignored --nocapture`.
    #[test]
    #[ignore = "prints the auto-designed rocket"]
    fn dump_best_design() {
        let d = auto_design(8_000).expect("feasible");
        println!(
            "\n=== AI-AUTO-DESIGNED LV-1 ({} virtual designs) ===",
            d.designs_evaluated
        );
        println!(
            "engine : Pc {:.0} bar · ε {:.1} → SL Isp {:.0} s · vac Isp {:.0} s",
            d.engine.chamber_pressure / 1.0e5,
            d.engine.expansion_ratio,
            d.engine_sea_level.isp,
            d.engine_vacuum.isp,
        );
        println!(
            "cooling: throat flux {:.1} MW/m² · margin {:.2}",
            d.engine_cooling.throat_heat_flux / 1.0e6,
            d.engine_cooling.cooling_margin,
        );
        println!(
            "ascent : {:.0} kg → {:.0} × {:.0} km · max-Q {:.0} kPa · peak {:.1} g",
            d.payload_kg, d.periapsis_km, d.apoapsis_km, d.max_q_kpa, d.peak_g,
        );
        println!(
            "guidance: pitch {:.1}° · rise {:.0} s · interstage SF {:.2}",
            d.pitch_kick_deg, d.vertical_rise_s, d.structural_sf,
        );
        println!(
            "search : {} feasible of {} designs evaluated\n",
            d.feasible_count, d.designs_evaluated,
        );
    }
}
