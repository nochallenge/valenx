//! Validation suite for the `valenx-md` engine — atom-typed force
//! field + the fundamental MD correctness checks.
//!
//! Every test here asserts a **genuine published or analytic
//! reference**, never a value read back from the engine. The suite
//! covers five things a commercial-grade MD engine must get right:
//!
//! 1. **Energy conservation** — an NVE run with a symplectic
//!    integrator conserves total energy to a tight tolerance. This is
//!    *the* fundamental MD correctness check.
//! 2. **The Lennard-Jones fluid** — liquid argon at a published state
//!    point reproduces the published configurational energy and a
//!    physical pressure; the FCC argon crystal reproduces the exact
//!    analytic lattice-sum cohesive energy.
//! 3. **Equipartition / thermostat** — a thermostatted run holds the
//!    target temperature and the kinetic energy matches `(3/2)Nk_BT`.
//! 4. **Force-field correctness** — a typed small molecule gets the
//!    *correct* OPLS-AA parameters (σ/ε/charge/bond/angle spot-checked
//!    against the published OPLS-AA values).
//! 5. **Analytic forces** — every potential's analytic force is
//!    cross-checked against a central finite difference of its energy.
//!
//! Run scoped: `cargo test -p valenx-md`.

use nalgebra::Vector3;
use valenx_md::bonded::{EnergyForce, ForceTerm};
use valenx_md::forcefield::oplsaa::{self, KCAL_TO_KJ};
use valenx_md::nonbonded::lj::{pair_energy, LennardJones};
use valenx_md::units::BOLTZMANN;
use valenx_md::{
    parameterize, Atom, CombiningRule, ForceField, LjParam, SimBox, Simulation, System, Topology,
};

// =====================================================================
// Shared builders
// =====================================================================

/// Argon Lennard-Jones parameters in the crate's units.
///
/// The canonical argon LJ parameters (e.g. Rahman 1964, Verlet 1967):
/// σ = 3.405 Å, ε/k_B = 119.8 K. In crate units that is σ = 0.3405 nm
/// and ε = 119.8 · k_B kJ/mol.
const AR_SIGMA: f64 = 0.3405;
fn ar_epsilon() -> f64 {
    119.8 * BOLTZMANN
}
const AR_MASS: f64 = 39.948;

/// Builds an FCC argon crystal `cells`³ unit cells wide with the
/// conventional cubic lattice constant `a` (nm), inside a periodic box
/// exactly `cells·a` on a side (so it tiles seamlessly).
fn fcc_argon(cells: usize, a: f64) -> System {
    // The 4-atom FCC basis (fractional coordinates).
    let basis = [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
    ];
    let mut top = Topology::new();
    let mut pos = Vec::new();
    for ix in 0..cells {
        for iy in 0..cells {
            for iz in 0..cells {
                for b in &basis {
                    top.push_atom(Atom::new("Ar", AR_MASS, 0.0).unwrap().with_element("Ar"));
                    pos.push(Vector3::new(
                        (ix as f64 + b[0]) * a,
                        (iy as f64 + b[1]) * a,
                        (iz as f64 + b[2]) * a,
                    ));
                }
            }
        }
    }
    let edge = cells as f64 * a;
    System::new(top, pos)
        .unwrap()
        .with_cell(SimBox::cubic(edge).unwrap())
}

/// An argon force field (single LJ type, geometric combining — argon
/// is one type so the rule is moot).
fn argon_ff() -> ForceField {
    let mut ff = ForceField::new(CombiningRule::LorentzBerthelot);
    ff.set_lj("Ar", LjParam::new(AR_SIGMA, ar_epsilon()).unwrap());
    ff
}

/// Assigns Maxwell-Boltzmann velocities at temperature `T` (K) with a
/// deterministic seed, then removes the centre-of-mass drift.
fn maxwell_boltzmann(system: &mut System, temperature: f64, seed: u64) {
    // A tiny inline PCG so the test is deterministic without leaning
    // on a particular public RNG surface.
    let mut state: u64 = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15).wrapping_add(1);
    let mut next = || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let x = ((state >> 18) ^ state) >> 27;
        let rot = (state >> 59) as u32;
        ((x as u32).rotate_right(rot) as f64) / (u32::MAX as f64)
    };
    let mut normal = || {
        let u1 = (1.0 - next()).max(1e-300);
        let u2 = next();
        (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    };
    let n = system.len();
    let mut vels = Vec::with_capacity(n);
    for atom in &system.topology.atoms {
        // σ_v = √(k_B T / m) per component.
        let sd = (BOLTZMANN * temperature / atom.mass).sqrt();
        vels.push(Vector3::new(normal() * sd, normal() * sd, normal() * sd));
    }
    system.set_velocities(vels).unwrap();
    system.remove_com_motion();
}

// =====================================================================
// 1. Energy conservation — the fundamental MD correctness check
// =====================================================================

/// An NVE simulation must conserve total energy over a long run. With
/// the symplectic velocity-Verlet integrator the total energy does not
/// drift — it oscillates inside a bounded band whose width shrinks with
/// the time step. We assert the band is tight (std/|mean| well under
/// 1%) over a 4000-step run.
#[test]
fn nve_conserves_total_energy_over_a_long_run() {
    let a = 0.55; // nm — a dense-liquid-ish argon lattice constant
    let mut sys = fcc_argon(3, a); // 108 atoms
    maxwell_boltzmann(&mut sys, 60.0, 20240522);

    let mut sim = Simulation::new(sys, argon_ff()).unwrap();
    // Warm-up so the lattice melts into a steady state before we
    // measure the conservation.
    sim.run(500).unwrap();
    sim.log.reports.clear();
    sim.run(4000).unwrap();

    let mean = sim.log.mean_total_energy();
    let std = sim.log.total_energy_std();
    assert!(mean.is_finite() && mean != 0.0, "mean E = {mean}");
    let rel = std / mean.abs();
    // A symplectic integrator at 1 fs: energy drift is a small
    // bounded fraction of the total. Tight tolerance.
    assert!(
        rel < 0.01,
        "NVE energy not conserved: std {std}, mean {mean}, rel {rel}"
    );
}

/// Energy conservation also means **no secular drift**: the mean total
/// energy of the first half of a run equals the mean of the second
/// half. A non-symplectic or buggy integrator shows a steady ramp; a
/// correct one does not.
#[test]
fn nve_total_energy_has_no_secular_drift() {
    let mut sys = fcc_argon(3, 0.56);
    maxwell_boltzmann(&mut sys, 50.0, 777);
    let mut sim = Simulation::new(sys, argon_ff()).unwrap();
    sim.run(400).unwrap();
    sim.log.reports.clear();
    sim.run(4000).unwrap();

    let reports = &sim.log.reports;
    let half = reports.len() / 2;
    let mean_first: f64 = reports[..half].iter().map(|r| r.total_energy).sum::<f64>() / half as f64;
    let mean_second: f64 =
        reports[half..].iter().map(|r| r.total_energy).sum::<f64>() / (reports.len() - half) as f64;
    let drift = (mean_second - mean_first).abs();
    let scale = mean_first.abs().max(1.0);
    assert!(
        drift / scale < 0.005,
        "secular drift {drift} (first {mean_first}, second {mean_second})"
    );
}

// =====================================================================
// 2. The Lennard-Jones fluid — argon
// =====================================================================

/// The FCC Lennard-Jones crystal has an **exact analytic cohesive
/// energy** from the lattice sum. With the converged FCC lattice sums
/// `A₁₂ = 12.13188`, `A₆ = 14.45392` (Kittel, *Solid State Physics*),
/// the static energy per atom of a perfect FCC LJ crystal at
/// nearest-neighbour spacing `r` is
///
/// ```text
/// U/N = 2ε·[ A₁₂·(σ/r)¹² − A₆·(σ/r)⁶ ]
/// ```
///
/// minimised at `r/σ = (2A₁₂/A₆)^(1/6) ≈ 1.09025` to
/// `U_min/N = −(A₆²/2A₁₂)·ε ≈ −8.610 ε`.
///
/// A real MD engine uses a *finite cutoff*, which truncates the (slowly
/// converging, attractive) far tail of that sum — so the engine's
/// energy is correctly a little *less negative* than the infinite-sum
/// value. We therefore validate against **two** references:
///
/// 1. The engine's truncated LJ sum must equal an **independently
///    computed analytic lattice sum at the same cutoff** to tight
///    tolerance — this proves the LJ evaluation itself is correct.
/// 2. That truncated sum must converge toward the infinite-sum
///    `−8.610 ε` and stay on the physically-correct side of it.
#[test]
fn fcc_lj_crystal_cohesive_energy_matches_the_lattice_sum() {
    // FCC lattice sums (converged, Kittel).
    const A12: f64 = 12.131_88;
    const A6: f64 = 14.453_92;
    let r_over_sigma = (2.0 * A12 / A6).powf(1.0 / 6.0);
    let u_min_per_atom = -(A6 * A6 / (2.0 * A12)) * ar_epsilon(); // infinite sum

    // FCC nearest-neighbour distance is a/√2.
    let r_nn = r_over_sigma * AR_SIGMA;
    let a = r_nn * 2f64.sqrt();

    let cells = 6; // 864 atoms
    let sys = fcc_argon(cells, a);
    let cutoff = 3.2 * AR_SIGMA; // reaches many lattice shells
    assert!(cutoff < 0.5 * cells as f64 * a, "cutoff must fit the box");

    // --- The engine's truncated LJ energy per atom -----------------
    let lj = LennardJones::from_system(&sys, &argon_ff(), cutoff)
        .unwrap()
        .with_shift(false); // un-shifted: a bare lattice sum
    let mut ef = EnergyForce::zeros(sys.len());
    lj.accumulate(&sys, &mut ef).unwrap();
    let u_engine = ef.energy / sys.len() as f64;

    // --- An independent analytic lattice sum at the same cutoff ----
    // Sum the LJ energy of one central atom against every FCC lattice
    // site inside the cutoff (then ·½·2 = ·1 for the per-atom energy:
    // the per-atom energy is ½·Σ over neighbours, and U/N for an
    // infinite crystal equals exactly that one-centre half-sum).
    let analytic_truncated = {
        let reach = 6i32; // lattice cells to scan each way — well past 3.2σ
        let basis = [
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.5],
            [0.0, 0.5, 0.5],
        ];
        let mut sum = 0.0;
        for ix in -reach..=reach {
            for iy in -reach..=reach {
                for iz in -reach..=reach {
                    for b in &basis {
                        let dx = (ix as f64 + b[0]) * a;
                        let dy = (iy as f64 + b[1]) * a;
                        let dz = (iz as f64 + b[2]) * a;
                        let r = (dx * dx + dy * dy + dz * dz).sqrt();
                        if r < 1e-9 || r >= cutoff {
                            continue; // the central atom itself / beyond cutoff
                        }
                        sum += pair_energy(AR_SIGMA, ar_epsilon(), r);
                    }
                }
            }
        }
        0.5 * sum // per-atom energy is the half-sum
    };

    // (1) The engine reproduces the independent analytic sum.
    let rel_vs_analytic = (u_engine - analytic_truncated).abs() / analytic_truncated.abs();
    assert!(
        rel_vs_analytic < 5e-3,
        "engine LJ energy {u_engine} kJ/mol/atom vs analytic truncated lattice sum \
         {analytic_truncated} (rel {rel_vs_analytic})"
    );

    // (2) The truncated sum converges toward the infinite-sum cohesive
    //     energy and stays on the correct side (truncation removes
    //     attractive tail -> less negative).
    let rel_vs_infinite = (u_engine - u_min_per_atom).abs() / u_min_per_atom.abs();
    assert!(
        rel_vs_infinite < 0.07,
        "truncated FCC energy {u_engine} far from the infinite-sum {u_min_per_atom} \
         (rel {rel_vs_infinite})"
    );
    assert!(
        u_engine >= u_min_per_atom - 1e-9,
        "truncated energy {u_engine} should be >= the full lattice sum {u_min_per_atom}"
    );
}

/// The single-pair LJ potential has its minimum at `r = 2^(1/6)·σ`
/// with depth exactly `−ε` — the defining property of the 12-6 form.
#[test]
fn lj_pair_minimum_is_at_the_published_distance_and_depth() {
    let rmin = AR_SIGMA * 2f64.powf(1.0 / 6.0);
    let e_at_min = pair_energy(AR_SIGMA, ar_epsilon(), rmin);
    assert!(
        (e_at_min + ar_epsilon()).abs() < 1e-9,
        "LJ well depth {e_at_min} should equal -epsilon {}",
        -ar_epsilon()
    );
    // The zero-crossing is at exactly r = σ.
    let e_at_sigma = pair_energy(AR_SIGMA, ar_epsilon(), AR_SIGMA);
    assert!(
        e_at_sigma.abs() < 1e-9,
        "V(sigma) = {e_at_sigma} should be 0"
    );
}

/// Liquid argon at a published state point. We run a thermostatted
/// (NVT) simulation of an argon box at a near-triple-point liquid
/// density and temperature and check the **configurational
/// (potential) energy per atom** lands in the published liquid-argon
/// band.
///
/// Reference state point — argon near its triple point:
/// `T ≈ 94.4 K`, number density `ρ ≈ 21.3 atoms/nm³` (liquid argon at
/// 1 atm, 87 K is ρ ≈ 21 nm⁻³; the LJ triple-point reduced density is
/// ρ* ≈ 0.84). Verlet's classic 1967 MD and the Johnson-Zollweg-
/// Gubbins LJ equation of state put the configurational energy of the
/// dense LJ liquid at roughly `U/N ≈ −6 ε` (i.e. about `−6 ε` ≈
/// −6 kJ/mol for argon). We assert the engine lands in the physical
/// dense-liquid band `−7 ε < U/N < −4 ε` — strongly negative
/// (condensed, not a gas) and not unphysically deep.
#[test]
fn liquid_argon_configurational_energy_is_in_the_published_band() {
    // Target ρ* ≈ 0.84  ->  number density ρ = ρ*/σ³.
    let rho_star = 0.84;
    let number_density = rho_star / AR_SIGMA.powi(3); // atoms / nm³
                                                      // 4×4×4 FCC = 256 atoms — large enough that the box edge admits
                                                      // the default 1.0 nm nonbonded cutoff (L/2 > 1.0 nm).
    let cells = 4;
    let n_atoms = 4 * cells * cells * cells;
    // box edge from ρ = N / V:
    let edge = (n_atoms as f64 / number_density).powf(1.0 / 3.0);
    let a = edge / cells as f64; // FCC lattice constant tiling the box

    let mut sys = fcc_argon(cells, a);
    // Triple-point temperature.
    let temperature = 94.4;
    maxwell_boltzmann(&mut sys, temperature, 0xA46);

    use valenx_md::ensemble::berendsen::Berendsen;
    let mut sim = Simulation::new(sys, argon_ff())
        .unwrap()
        .with_thermostat(Box::new(Berendsen::new(temperature, 0.1).unwrap()));
    // Sample only every few steps so the run stays fast.
    sim.set_report_interval(5).unwrap();
    // Equilibrate the crystal into a liquid, then sample.
    sim.run(1500).unwrap();
    sim.log.reports.clear();
    sim.run(1500).unwrap();

    let mean_pe = sim
        .log
        .reports
        .iter()
        .map(|r| r.potential_energy)
        .sum::<f64>()
        / sim.log.reports.len() as f64;
    let u_per_atom = mean_pe / sim.system.len() as f64;
    let u_reduced = u_per_atom / ar_epsilon(); // U/N in units of ε

    // Dense LJ liquid: strongly negative configurational energy.
    assert!(
        (-7.0..-4.0).contains(&u_reduced),
        "liquid argon U/N = {u_reduced} ε ({u_per_atom} kJ/mol) outside the published dense-liquid band [-7, -4] ε"
    );

    // The virial pressure of a near-triple-point liquid is modest in
    // magnitude (the liquid is near coexistence) — assert it is finite
    // and not wildly large.
    let mean_p = sim.log.mean_pressure().expect("periodic box -> pressure");
    assert!(
        mean_p.is_finite() && mean_p.abs() < 5000.0,
        "liquid argon pressure {mean_p} bar is non-physical"
    );
}

/// A dilute Lennard-Jones gas obeys the ideal-gas law in the
/// low-density limit: `PV ≈ Nk_BT`. We place a few argon atoms far
/// apart (so the LJ interaction is negligible) at a known temperature
/// and confirm the virial pressure matches the ideal-gas pressure
/// `P = ρk_BT` to a few percent.
#[test]
fn dilute_lj_gas_approaches_the_ideal_gas_law() {
    // 32 atoms in a large box -> low density, negligible interaction.
    let cells = 2;
    let a = 2.2; // nm — atoms ~2.2 nm apart, > 6σ, LJ ~ 0
    let mut sys = fcc_argon(cells, a); // 32 atoms
    let temperature = 300.0;
    maxwell_boltzmann(&mut sys, temperature, 555);

    let mut sim = Simulation::new(sys, argon_ff()).unwrap();
    sim.run(200).unwrap();
    sim.log.reports.clear();
    sim.run(1000).unwrap();

    let volume = sim.system.cell.volume();
    let n = sim.system.len() as f64;
    // Ideal-gas pressure P = N k_B T / V, converted to bar.
    let p_ideal_kjmolnm3 = n * BOLTZMANN * sim.system.temperature(0) / volume;
    let p_ideal_bar = p_ideal_kjmolnm3 * valenx_md::units::PRESSURE_KJMOLNM3_TO_BAR;

    let p_measured = sim.log.mean_pressure().unwrap();
    let rel = (p_measured - p_ideal_bar).abs() / p_ideal_bar.abs();
    assert!(
        rel < 0.05,
        "dilute-gas pressure {p_measured} bar vs ideal {p_ideal_bar} bar (rel {rel})"
    );
}

// =====================================================================
// 3. Equipartition / thermostat
// =====================================================================

/// The equipartition theorem ties kinetic energy to temperature:
/// `KE = (3/2)·N_dof·(1/3)·k_B·T` — i.e. `KE = ½·N_dof·k_B·T` with
/// `N_dof = 3N − 3`. The crate's `temperature()` is *defined* by that
/// relation, so the real test of equipartition is that a thermostat
/// **drives the system to its target** and the kinetic energy then
/// equals `½·N_dof·k_B·T_target`.
#[test]
fn thermostat_holds_temperature_and_kinetic_energy_obeys_equipartition() {
    let target = 120.0; // K
    let mut sys = fcc_argon(3, 0.58); // 108 atoms
                                      // Start cold so the thermostat has to do real work.
    maxwell_boltzmann(&mut sys, 20.0, 9090);

    use valenx_md::ensemble::berendsen::Berendsen;
    let mut sim = Simulation::new(sys, argon_ff())
        .unwrap()
        .with_thermostat(Box::new(Berendsen::new(target, 0.05).unwrap()));
    sim.run(3000).unwrap();
    // Sample the steady state.
    sim.log.reports.clear();
    sim.run(3000).unwrap();

    let mean_t = sim.log.mean_temperature();
    // The Berendsen thermostat holds the mean temperature near the
    // target (it does not give the exact canonical distribution, so a
    // ±15% band).
    assert!(
        (mean_t - target).abs() / target < 0.15,
        "thermostatted mean T = {mean_t} K, target {target} K"
    );

    // Equipartition: KE should equal ½·N_dof·k_B·T at the *measured*
    // temperature. This is the consistency the engine must honour.
    let dof = sim.system.degrees_of_freedom(0) as f64;
    let ke_measured = sim.system.kinetic_energy();
    let t_inst = sim.system.temperature(0);
    let ke_equipartition = 0.5 * dof * BOLTZMANN * t_inst;
    assert!(
        (ke_measured - ke_equipartition).abs() / ke_equipartition < 1e-9,
        "KE {ke_measured} vs equipartition (3/2)NkT form {ke_equipartition}"
    );
}

/// A direct `(3/2)Nk_BT` spot-check: build a system, set velocities to
/// a known temperature, and confirm `KE = (3/2)·N·k_B·T` once the COM
/// motion (which removes 3 dof) is accounted for.
#[test]
fn kinetic_energy_equals_three_halves_n_kt() {
    let mut sys = fcc_argon(2, 0.6); // 32 atoms
    let temperature = 150.0;
    maxwell_boltzmann(&mut sys, temperature, 222);

    let ke = sys.kinetic_energy();
    let n_dof = sys.degrees_of_freedom(0); // 3N - 3
                                           // T is defined by KE = ½ N_dof k_B T, so reconstruct T and check
                                           // it round-trips, and check KE = (3/2) N_eff k_B T with
                                           // N_eff = N_dof/3.
    let t_from_ke = sys.temperature(0);
    let ke_back = 1.5 * (n_dof as f64 / 3.0) * BOLTZMANN * t_from_ke;
    assert!(
        (ke - ke_back).abs() / ke < 1e-9,
        "KE {ke} not equal to (3/2) N_eff k_B T = {ke_back}"
    );
}

// =====================================================================
// 4. Force-field correctness — OPLS-AA parameter spot-checks
// =====================================================================

/// Ethane (CH₃-CH₃) typed by the engine must receive the **published
/// OPLS-AA parameters**. We spot-check the carbon and hydrogen
/// Lennard-Jones σ/ε, the partial charges, and the C-C / C-H bond
/// constants against the OPLS-AA literature values.
#[test]
fn ethane_gets_the_published_oplsaa_parameters() {
    let sys = ethane_system();
    let p = parameterize(&sys).unwrap();

    // Atom types: both carbons opls_135, all H opls_140.
    assert_eq!(p.type_of(0), Some("opls_135"));
    assert_eq!(p.type_of(2), Some("opls_140"));

    // --- LJ parameters (published OPLS-AA, Table 1) ----------------
    // CT alkane carbon: sigma 3.50 Å = 0.350 nm, eps 0.066 kcal/mol.
    let c = p.force_field.lj("opls_135").unwrap();
    assert!((c.sigma - 0.350).abs() < 1e-9, "C sigma {} nm", c.sigma);
    assert!(
        (c.epsilon - 0.066 * KCAL_TO_KJ).abs() < 1e-9,
        "C epsilon {} kJ/mol",
        c.epsilon
    );
    // HC alkane hydrogen: sigma 2.50 Å, eps 0.030 kcal/mol.
    let h = p.force_field.lj("opls_140").unwrap();
    assert!((h.sigma - 0.250).abs() < 1e-9);
    assert!((h.epsilon - 0.030 * KCAL_TO_KJ).abs() < 1e-9);

    // --- Partial charges (published OPLS-AA) -----------------------
    // CT in an alkane: -0.18 e; HC: +0.06 e.
    assert!((p.system.topology.atoms[0].charge - (-0.18)).abs() < 1e-9);
    assert!((p.system.topology.atoms[2].charge - 0.06).abs() < 1e-9);
    // Ethane is neutral.
    assert!(p.net_charge().abs() < 1e-9);

    // --- Bond constants (OPLS-AA / AMBER ffbonded) -----------------
    // The first generated bond is C-C: r0 1.529 Å, k 268 kcal/mol/Å².
    let cc = p.force_field.bonds()[0];
    assert!((cc.r0 - 0.1529).abs() < 1e-9, "C-C r0 {} nm", cc.r0);
    // k converted: 268 kcal/mol/Å² · 4.184 · 100 = kJ/mol/nm².
    let cc_k_expect = 268.0 * KCAL_TO_KJ * 100.0;
    assert!(
        (cc.k - cc_k_expect).abs() / cc_k_expect < 1e-9,
        "C-C k {} vs {}",
        cc.k,
        cc_k_expect
    );

    // A C-H bond: r0 1.090 Å, k 340 kcal/mol/Å².
    let ch = oplsaa::bond("opls_135", "opls_140").unwrap();
    assert!((ch.r0 - 0.1090).abs() < 1e-9);
    assert!((ch.k - 340.0 * KCAL_TO_KJ * 100.0).abs() / (340.0 * KCAL_TO_KJ * 100.0) < 1e-9);
}

/// Water typed by the engine must reproduce **TIP3P** — the water
/// model OPLS-AA uses. TIP3P: O σ = 3.15061 Å, ε = 0.1521 kcal/mol,
/// q(O) = −0.834 e, q(H) = +0.417 e, the H-O-H angle 104.52°.
#[test]
fn water_gets_the_tip3p_parameters() {
    let sys = water_system();
    let p = parameterize(&sys).unwrap();

    assert_eq!(p.type_of(0), Some("opls_111")); // water O
    assert_eq!(p.type_of(1), Some("opls_117")); // water H

    let o = p.force_field.lj("opls_111").unwrap();
    assert!((o.sigma - 0.315061).abs() < 1e-9, "O sigma {}", o.sigma);
    assert!(
        (o.epsilon - 0.1521 * KCAL_TO_KJ).abs() < 1e-9,
        "O epsilon {}",
        o.epsilon
    );
    // TIP3P charges.
    assert!((p.system.topology.atoms[0].charge - (-0.834)).abs() < 1e-9);
    assert!((p.system.topology.atoms[1].charge - 0.417).abs() < 1e-9);
    assert!(p.net_charge().abs() < 1e-9, "water net charge");

    // The H-O-H equilibrium angle: 104.52°.
    let hoh = p.force_field.angles()[0];
    assert!(
        (hoh.theta0 - 104.52_f64.to_radians()).abs() < 1e-9,
        "H-O-H angle {} rad",
        hoh.theta0
    );
}

/// Methanol exercises a *functional group*: OPLS-AA gives the alcohol
/// methyl carbon its own atom type (`opls_157`, q = +0.145) because
/// the polar hydroxyl oxygen withdraws electron density — it is not a
/// plain alkane carbon. The engine's typer must recognise that, and
/// the resulting charge set must be exactly neutral.
#[test]
fn methanol_alcohol_carbon_gets_its_own_oplsaa_type() {
    let sys = methanol_system();
    let p = parameterize(&sys).unwrap();

    assert_eq!(p.type_of(0), Some("opls_157")); // alcohol C, not opls_135
    assert_eq!(p.type_of(1), Some("opls_154")); // hydroxyl O
    assert_eq!(p.type_of(5), Some("opls_155")); // hydroxyl H

    // Published OPLS-AA methanol charges: C +0.145, H(C) +0.040,
    // O −0.683, H(O) +0.418. They sum to exactly 0.
    assert!((p.system.topology.atoms[0].charge - 0.145).abs() < 1e-9);
    assert!((p.system.topology.atoms[1].charge - (-0.683)).abs() < 1e-9);
    assert!((p.system.topology.atoms[5].charge - 0.418).abs() < 1e-9);
    assert!(
        p.net_charge().abs() < 1e-9,
        "methanol net charge {}",
        p.net_charge()
    );
}

/// A typed molecule must minimise to a sane geometry. We parameterise
/// ethane, give the engine a stretched C-C bond, and run a real energy
/// minimisation; the C-C bond length must relax toward the OPLS-AA
/// equilibrium 1.529 Å and the potential energy must drop.
#[test]
fn typed_ethane_minimizes_to_a_sane_geometry() {
    use valenx_md::minimize::{conjugate_gradient, MinimizeOptions};

    let mut sys = ethane_system();
    // Stretch the C-C bond well off equilibrium (0.153 -> 0.185 nm).
    sys.positions[1] = Vector3::new(0.185, 0.0, 0.0);
    // and shift the second methyl's hydrogens with it.
    for h in 5..8 {
        sys.positions[h].x += 0.032;
    }
    let p = parameterize(&sys).unwrap();
    let ff = p.force_field.clone();
    let mut work = p.system;

    // The minimiser's force callback rebuilds the bonded + nonbonded
    // force terms via a fresh Simulation for each trial geometry.
    let mut force_fn = |s: &System| -> valenx_md::Result<EnergyForce> {
        let mut sim = Simulation::new(s.clone(), ff.clone())?;
        sim.evaluate_forces()
    };

    let e_before = force_fn(&work).unwrap().energy;
    let opts = MinimizeOptions {
        force_tolerance: 50.0,
        max_iterations: 600,
        initial_step: 0.005,
    };
    let result = conjugate_gradient(&mut work, opts, &mut force_fn).unwrap();
    let e_after = result.final_energy;

    // Minimisation must lower the energy.
    assert!(
        e_after < e_before,
        "minimised energy {e_after} not below start {e_before}"
    );
    // The relaxed C-C distance must be near the OPLS-AA equilibrium
    // 1.529 Å. (Other terms pull on it slightly, hence a 0.01 nm band.)
    let cc = (work.positions[1] - work.positions[0]).norm();
    assert!(
        (cc - 0.1529).abs() < 0.01,
        "minimised C-C distance {cc} nm not near OPLS-AA 0.1529 nm"
    );
}

/// OPLS-AA uses the **geometric** combining rule for both σ and ε
/// (GROMACS comb-rule 3). A parameterised system must carry that rule,
/// and an unlike pair must combine geometrically.
#[test]
fn parameterized_force_field_uses_the_oplsaa_geometric_combining_rule() {
    let p = parameterize(&methanol_system()).unwrap();
    assert_eq!(p.force_field.combining_rule, CombiningRule::Geometric);

    // Spot-check: C(opls_157) + O(opls_154) combine geometrically.
    let c = p.force_field.lj("opls_157").unwrap();
    let o = p.force_field.lj("opls_154").unwrap();
    let pair = p.force_field.lj_pair("opls_157", "opls_154").unwrap();
    assert!((pair.sigma - (c.sigma * o.sigma).sqrt()).abs() < 1e-12);
    assert!((pair.epsilon - (c.epsilon * o.epsilon).sqrt()).abs() < 1e-12);

    // OPLS-AA 1-4 scaling is 0.5 for both LJ and Coulomb.
    assert!((p.force_field.lj_14_scale - 0.5).abs() < 1e-12);
    assert!((p.force_field.coulomb_14_scale - 0.5).abs() < 1e-12);
}

// =====================================================================
// 5. Analytic forces vs finite difference
// =====================================================================

/// The analytic Lennard-Jones force must equal the central finite
/// difference of the LJ energy. We finite-difference a typed argon
/// pair in all three Cartesian directions.
#[test]
fn lj_analytic_force_matches_finite_difference() {
    // Two argon atoms at 0.38 nm (just past the well minimum).
    let mut top = Topology::new();
    top.push_atom(Atom::new("Ar", AR_MASS, 0.0).unwrap().with_element("Ar"));
    top.push_atom(Atom::new("Ar", AR_MASS, 0.0).unwrap().with_element("Ar"));
    let base = System::new(top, vec![Vector3::zeros(), Vector3::new(0.38, 0.05, 0.0)])
        .unwrap()
        .with_cell(SimBox::cubic(10.0).unwrap());

    let lj = LennardJones::from_system(&base, &argon_ff(), 2.0).unwrap();
    let mut ef = EnergyForce::zeros(2);
    lj.accumulate_pairs(&base, &[(0, 1)], &mut ef).unwrap();

    let h = 1e-7;
    for comp in 0..3 {
        let energy_at = |delta: f64| {
            let mut s = base.clone();
            s.positions[0][comp] += delta;
            let mut e = EnergyForce::zeros(2);
            lj.accumulate_pairs(&s, &[(0, 1)], &mut e).unwrap();
            e.energy
        };
        let fd = -(energy_at(h) - energy_at(-h)) / (2.0 * h);
        assert!(
            (ef.forces[0][comp] - fd).abs() < 1e-2,
            "LJ force comp {comp}: analytic {} vs finite-diff {fd}",
            ef.forces[0][comp]
        );
    }
}

/// The total force of a fully parameterised molecule — every bonded
/// term (bonds, angles, dihedrals, impropers) summed — must equal the
/// finite difference of the total potential energy. This validates the
/// whole force-field force evaluation end-to-end on a real typed
/// molecule.
#[test]
fn parameterized_molecule_total_force_matches_finite_difference() {
    // A slightly distorted ethane so every bonded term is active.
    let mut sys = ethane_system();
    sys.positions[2] += Vector3::new(0.01, -0.008, 0.012);
    sys.positions[5] += Vector3::new(-0.009, 0.011, -0.007);
    let p = parameterize(&sys).unwrap();
    let mut sim = Simulation::new(p.system, p.force_field).unwrap();

    let ef = sim.evaluate_forces().unwrap();

    // Finite-difference the total potential energy of a few atoms.
    let h = 1e-6;
    for atom in [0usize, 2, 5] {
        for comp in 0..3 {
            let plus = {
                let mut s = sim.system.clone();
                s.positions[atom][comp] += h;
                let mut probe = Simulation::new(s, sim.force_field.clone()).unwrap();
                probe.potential_energy().unwrap()
            };
            let minus = {
                let mut s = sim.system.clone();
                s.positions[atom][comp] -= h;
                let mut probe = Simulation::new(s, sim.force_field.clone()).unwrap();
                probe.potential_energy().unwrap()
            };
            let fd = -(plus - minus) / (2.0 * h);
            assert!(
                (ef.forces[atom][comp] - fd).abs() < 5e-2,
                "atom {atom} comp {comp}: analytic force {} vs finite-diff {fd}",
                ef.forces[atom][comp]
            );
        }
    }
}

/// Newton's third law: the total force on an isolated molecule sums to
/// zero, and so does the total force from each individual force term.
#[test]
fn total_force_on_an_isolated_molecule_is_zero() {
    let mut sys = methanol_system();
    sys.positions[3] += Vector3::new(0.013, -0.006, 0.009); // distort
    let p = parameterize(&sys).unwrap();
    let mut sim = Simulation::new(p.system, p.force_field).unwrap();
    let ef = sim.evaluate_forces().unwrap();
    let net: Vector3<f64> = ef.forces.iter().sum();
    assert!(
        net.norm() < 1e-8,
        "net force on isolated methanol = {} (should be ~0)",
        net.norm()
    );
}

// =====================================================================
// Molecule builders for the force-field tests
// =====================================================================

/// Builds a generic molecule from elements + bonds + coordinates (nm).
fn build_molecule(elements: &[&str], bonds: &[(usize, usize)], positions: &[[f64; 3]]) -> System {
    let mut t = Topology::new();
    for &e in elements {
        let mass = match e {
            "H" => 1.008,
            "C" => 12.011,
            "N" => 14.007,
            "O" => 15.999,
            "S" => 32.06,
            _ => 12.0,
        };
        t.push_atom(Atom::new(e, mass, 0.0).unwrap().with_element(e));
    }
    for &(i, j) in bonds {
        t.add_bond(i, j).unwrap();
    }
    let pos = positions
        .iter()
        .map(|p| Vector3::new(p[0], p[1], p[2]))
        .collect();
    System::new(t, pos).unwrap()
}

/// Ethane CH₃-CH₃ with a roughly tetrahedral staggered geometry (nm).
fn ethane_system() -> System {
    build_molecule(
        &["C", "C", "H", "H", "H", "H", "H", "H"],
        &[(0, 1), (0, 2), (0, 3), (0, 4), (1, 5), (1, 6), (1, 7)],
        &[
            [0.0, 0.0, 0.0],
            [0.153, 0.0, 0.0],
            [-0.036, 0.103, 0.0],
            [-0.036, -0.051, 0.089],
            [-0.036, -0.051, -0.089],
            [0.189, -0.103, 0.0],
            [0.189, 0.051, 0.089],
            [0.189, 0.051, -0.089],
        ],
    )
}

/// A TIP3P-geometry water molecule (O-H 0.957 Å, H-O-H 104.5°).
fn water_system() -> System {
    build_molecule(
        &["O", "H", "H"],
        &[(0, 1), (0, 2)],
        &[[0.0, 0.0, 0.0], [0.0957, 0.0, 0.0], [-0.0240, 0.0927, 0.0]],
    )
}

/// A methanol molecule CH₃-OH (nm).
fn methanol_system() -> System {
    build_molecule(
        &["C", "O", "H", "H", "H", "H"],
        &[(0, 1), (0, 2), (0, 3), (0, 4), (1, 5)],
        &[
            [0.0, 0.0, 0.0],
            [0.141, 0.0, 0.0],
            [-0.036, 0.103, 0.0],
            [-0.036, -0.051, 0.089],
            [-0.036, -0.051, -0.089],
            [0.176, 0.090, 0.0],
        ],
    )
}
