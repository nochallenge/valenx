//! The AIMD engine and the engine-agnostic data model.
//!
//! [`AimdEngine`] is the Phase-1 [`ReactionEngine`]: a Born-Oppenheimer
//! ab-initio MD run — velocity-Verlet ([`crate::integrator`]) over
//! numerical qchem forces ([`crate::forces`]). The [`System`],
//! [`Controls`], [`Frame`] and [`Trajectory`] types — and the
//! [`ReactionEngine`] trait — are deliberately engine-agnostic so the
//! future QM/MM and ReaxFF backends share them unchanged.
//!
//! State is in atomic units (bohr, hartree, electron-masses, atomic
//! time); the timestep crosses the boundary in femtoseconds.

use valenx_qchem::element::Element;
use valenx_qchem::geometry::BOHR_PER_ANGSTROM;

use crate::error::{ReactDynError, Result};
use crate::forces::{numerical_forces, single_point_energy, Method};
use crate::integrator::{kinetic_energy, velocity_verlet_step};
use crate::units::{amu_to_au_mass, fs_to_au};

/// Boltzmann constant in hartree per kelvin (CODATA 2018).
const BOLTZMANN_HARTREE_PER_K: f64 = 3.166_811_563e-6;

/// A molecular system to simulate. Positions are in **bohr**.
#[derive(Clone, Debug, PartialEq)]
pub struct System {
    /// One element per atom.
    pub elements: Vec<Element>,
    /// Atomic positions in bohr, one per element.
    pub pos_bohr: Vec<[f64; 3]>,
    /// Total molecular charge (units of e).
    pub charge: i32,
    /// Spin multiplicity `2S + 1` (always `>= 1`).
    pub multiplicity: u32,
}

impl System {
    /// Number of atoms.
    pub fn n_atoms(&self) -> usize {
        self.elements.len()
    }
}

/// Temperature-control scheme.
#[derive(Copy, Clone, Debug, PartialEq, Default)]
pub enum Thermostat {
    /// Microcanonical — no coupling; total energy is conserved (this is
    /// the built-in correctness check).
    #[default]
    Nve,
    /// Berendsen weak coupling to a target temperature.
    Berendsen {
        /// Target temperature in kelvin.
        target_kelvin: f64,
        /// Coupling time constant in femtoseconds.
        tau_fs: f64,
    },
}

/// Run controls for an AIMD trajectory.
#[derive(Clone, Debug, PartialEq)]
pub struct Controls {
    /// Electronic-structure method used for the forces.
    pub method: Method,
    /// Basis-set name (e.g. `"STO-3G"`).
    pub basis: String,
    /// Timestep in femtoseconds.
    pub dt_fs: f64,
    /// Number of integration steps.
    pub n_steps: usize,
    /// Finite-difference displacement for the forces, in bohr.
    pub fd_delta_bohr: f64,
    /// Temperature-control scheme.
    pub thermostat: Thermostat,
    /// Safety cap on `atoms × steps`, bounding the (expensive)
    /// numerical-gradient compute so a run can't lock up the machine.
    pub max_cost_guard: usize,
}

impl Default for Controls {
    fn default() -> Self {
        Controls {
            method: Method::Rhf,
            basis: "STO-3G".to_string(),
            dt_fs: 0.5,
            n_steps: 50,
            fd_delta_bohr: 0.01,
            thermostat: Thermostat::Nve,
            max_cost_guard: 4000, // e.g. 8 atoms × 500 steps
        }
    }
}

/// One recorded frame. Positions/velocities are atomic units (bohr,
/// bohr/atomic-time); energies are hartree.
#[derive(Clone, Debug, PartialEq)]
pub struct Frame {
    /// Simulation time in femtoseconds.
    pub time_fs: f64,
    /// Atomic positions in bohr.
    pub pos_bohr: Vec<[f64; 3]>,
    /// Atomic velocities (bohr / atomic-time).
    pub vel: Vec<[f64; 3]>,
    /// Electronic (potential) energy in hartree.
    pub potential_hartree: f64,
    /// Kinetic energy in hartree.
    pub kinetic_hartree: f64,
}

impl Frame {
    /// Total (potential + kinetic) energy in hartree.
    pub fn total_hartree(&self) -> f64 {
        self.potential_hartree + self.kinetic_hartree
    }

    /// Positions converted to ångström — for display / export.
    pub fn pos_angstrom(&self) -> Vec<[f64; 3]> {
        self.pos_bohr
            .iter()
            .map(|p| {
                [
                    p[0] / BOHR_PER_ANGSTROM,
                    p[1] / BOHR_PER_ANGSTROM,
                    p[2] / BOHR_PER_ANGSTROM,
                ]
            })
            .collect()
    }
}

/// A computed trajectory.
#[derive(Clone, Debug, PartialEq)]
pub struct Trajectory {
    /// The system that was simulated.
    pub system: System,
    /// Frames in time order — `n_steps + 1` of them (the initial frame
    /// plus one per step).
    pub frames: Vec<Frame>,
}

/// A reaction-dynamics engine. [`AimdEngine`] is the Phase-1 implementor;
/// QM/MM and ReaxFF backends will implement the same trait so the
/// workbench drives them identically.
pub trait ReactionEngine {
    /// Run a trajectory for `system` under `controls`. `progress(step)`
    /// is invoked after each completed step (0-based) for UI feedback.
    fn run(
        &self,
        system: &System,
        controls: &Controls,
        progress: &mut dyn FnMut(usize),
    ) -> Result<Trajectory>;
}

/// Born-Oppenheimer ab-initio MD: velocity-Verlet over numerical qchem
/// forces.
pub struct AimdEngine;

impl ReactionEngine for AimdEngine {
    fn run(
        &self,
        system: &System,
        controls: &Controls,
        progress: &mut dyn FnMut(usize),
    ) -> Result<Trajectory> {
        let n = system.n_atoms();
        if n == 0 {
            return Err(ReactDynError::Invalid {
                reason: "empty system (no atoms)".into(),
            });
        }
        if system.pos_bohr.len() != n {
            return Err(ReactDynError::Invalid {
                reason: "positions and elements differ in length".into(),
            });
        }
        if controls.n_steps == 0 {
            return Err(ReactDynError::Invalid {
                reason: "n_steps must be at least 1".into(),
            });
        }
        if !(controls.dt_fs > 0.0 && controls.dt_fs.is_finite()) {
            return Err(ReactDynError::Invalid {
                reason: format!("timestep must be positive (got {} fs)", controls.dt_fs),
            });
        }
        let cost = n.saturating_mul(controls.n_steps);
        if cost > controls.max_cost_guard {
            return Err(ReactDynError::GuardExceeded {
                atoms: n,
                steps: controls.n_steps,
                cap: controls.max_cost_guard,
            });
        }

        let masses: Vec<f64> = system
            .elements
            .iter()
            .map(|e| amu_to_au_mass(e.atomic_mass()))
            .collect();
        let dt = fs_to_au(controls.dt_fs);
        let basis = controls.basis.as_str();
        let method = controls.method;
        let (charge, mult) = (system.charge, system.multiplicity);
        let elements = &system.elements;

        // Fallible force evaluation, reused by the integrator each step.
        let force_fn = |p: &[[f64; 3]]| -> Result<Vec<[f64; 3]>> {
            numerical_forces(
                elements,
                p,
                charge,
                mult,
                method,
                basis,
                controls.fd_delta_bohr,
            )
        };

        let mut pos = system.pos_bohr.clone();
        let mut vel = init_velocities(&masses, controls.thermostat);
        let mut forces = force_fn(&pos)?;

        let mut frames = Vec::with_capacity(controls.n_steps + 1);
        let pe0 = single_point_energy(elements, &pos, charge, mult, method, basis)?;
        frames.push(Frame {
            time_fs: 0.0,
            pos_bohr: pos.clone(),
            vel: vel.clone(),
            potential_hartree: pe0,
            kinetic_hartree: kinetic_energy(&vel, &masses),
        });

        for step in 0..controls.n_steps {
            forces = velocity_verlet_step(&mut pos, &mut vel, &forces, &masses, dt, force_fn)?;
            if let Thermostat::Berendsen {
                target_kelvin,
                tau_fs,
            } = controls.thermostat
            {
                berendsen_rescale(&mut vel, &masses, target_kelvin, controls.dt_fs, tau_fs);
            }
            let pe = single_point_energy(elements, &pos, charge, mult, method, basis)?;
            frames.push(Frame {
                time_fs: (step + 1) as f64 * controls.dt_fs,
                pos_bohr: pos.clone(),
                vel: vel.clone(),
                potential_hartree: pe,
                kinetic_hartree: kinetic_energy(&vel, &masses),
            });
            progress(step);
        }

        Ok(Trajectory {
            system: system.clone(),
            frames,
        })
    }
}

/// Initial velocities for the chosen thermostat. Shared with the QM/MM
/// engine.
pub(crate) fn init_velocities(masses: &[f64], thermostat: Thermostat) -> Vec<[f64; 3]> {
    match thermostat {
        Thermostat::Nve => vec![[0.0; 3]; masses.len()],
        Thermostat::Berendsen { target_kelvin, .. } => seed_velocities(masses, target_kelvin),
    }
}

/// Deterministic initial velocities at ≈ `target_kelvin`: each atom gets
/// speed `√(3kT/m)` along a cycling axis with alternating sign, then the
/// net momentum is removed so the centre of mass stays put. No RNG — the
/// run is fully reproducible.
fn seed_velocities(masses: &[f64], target_kelvin: f64) -> Vec<[f64; 3]> {
    let mut vel = vec![[0.0; 3]; masses.len()];
    if target_kelvin <= 0.0 {
        return vel;
    }
    let kt = BOLTZMANN_HARTREE_PER_K * target_kelvin;
    for (i, v) in vel.iter_mut().enumerate() {
        let speed = (3.0 * kt / masses[i]).sqrt();
        let axis = i % 3;
        let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
        v[axis] = sign * speed;
    }
    remove_net_momentum(&mut vel, masses);
    vel
}

/// Subtract the centre-of-mass velocity so the system doesn't drift.
fn remove_net_momentum(vel: &mut [[f64; 3]], masses: &[f64]) {
    let total_m: f64 = masses.iter().sum();
    if total_m <= 0.0 {
        return;
    }
    let mut p = [0.0; 3];
    for (v, &m) in vel.iter().zip(masses) {
        for d in 0..3 {
            p[d] += m * v[d];
        }
    }
    for d in 0..3 {
        let v_com = p[d] / total_m;
        for v in vel.iter_mut() {
            v[d] -= v_com;
        }
    }
}

/// Berendsen weak-coupling velocity rescale toward `target_kelvin`.
/// Shared with the QM/MM engine.
pub(crate) fn berendsen_rescale(
    vel: &mut [[f64; 3]],
    masses: &[f64],
    target_kelvin: f64,
    dt_fs: f64,
    tau_fs: f64,
) {
    // Degrees of freedom: 3N minus the 3 removed COM components.
    let n_dof = (3 * vel.len()).saturating_sub(3).max(1);
    let ke = kinetic_energy(vel, masses);
    let current_t = 2.0 * ke / (n_dof as f64 * BOLTZMANN_HARTREE_PER_K);
    if current_t <= 0.0 {
        return;
    }
    let tau = tau_fs.max(dt_fs); // never let dt/tau exceed 1
    let lambda = (1.0 + (dt_fs / tau) * (target_kelvin / current_t - 1.0))
        .max(0.0)
        .sqrt();
    for v in vel.iter_mut() {
        for x in v.iter_mut() {
            *x *= lambda;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn h2(separation_angstrom: f64) -> System {
        let h = Element::from_symbol("H").unwrap();
        let half = 0.5 * separation_angstrom * BOHR_PER_ANGSTROM;
        System {
            elements: vec![h, h],
            pos_bohr: vec![[0.0, 0.0, half], [0.0, 0.0, -half]],
            charge: 0,
            multiplicity: 1,
        }
    }

    #[test]
    fn h2_nve_conserves_total_energy() {
        let system = h2(0.9); // slightly stretched → it will vibrate
        let controls = Controls {
            method: Method::Rhf,
            basis: "STO-3G".to_string(),
            dt_fs: 0.5,
            n_steps: 20,
            fd_delta_bohr: 0.01,
            thermostat: Thermostat::Nve,
            max_cost_guard: 4000,
        };
        let mut steps_seen = 0usize;
        let traj = AimdEngine
            .run(&system, &controls, &mut |_| steps_seen += 1)
            .expect("AIMD run should succeed");

        assert_eq!(traj.frames.len(), controls.n_steps + 1);
        assert_eq!(steps_seen, controls.n_steps);

        let e0 = traj.frames[0].total_hartree();
        let max_drift = traj
            .frames
            .iter()
            .map(|f| (f.total_hartree() - e0).abs())
            .fold(0.0_f64, f64::max);
        assert!(
            max_drift < 5e-3,
            "NVE total-energy drift too large: {max_drift} hartree"
        );
    }

    #[test]
    fn empty_system_fails_loud() {
        let system = System {
            elements: vec![],
            pos_bohr: vec![],
            charge: 0,
            multiplicity: 1,
        };
        let r = AimdEngine.run(&system, &Controls::default(), &mut |_| {});
        assert!(matches!(r, Err(ReactDynError::Invalid { .. })));
    }

    #[test]
    fn cost_guard_refuses_oversized_runs() {
        let system = h2(0.9);
        let controls = Controls {
            n_steps: 10_000,
            max_cost_guard: 100,
            ..Default::default()
        };
        let r = AimdEngine.run(&system, &controls, &mut |_| {});
        assert!(matches!(r, Err(ReactDynError::GuardExceeded { .. })));
    }
}
