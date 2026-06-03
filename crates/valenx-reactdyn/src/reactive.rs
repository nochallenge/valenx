//! A compact **reactive** classical potential for many-atom / materials
//! dynamics: **Morse pair bonds with a smooth distance cutoff**.
//!
//! It is reactive in the simplest honest sense — a bond (a Morse well)
//! forms when two atoms come within the cutoff and breaks when they are
//! pulled past it, so bonds visibly appear and disappear during the run.
//! It is purely classical (no SCF), so it is fast and scales to many
//! atoms; the energy and forces are analytic (atomic units throughout).
//!
//! Honest scope: this is a **pair** potential — it has no valence
//! saturation, so atoms can over-coordinate. A Tersoff-class angular
//! bond-order potential and a full ReaxFF (bond order + over/under-
//! coordination + charge equilibration) are the documented upgrades.

use valenx_qchem::element::Element;

use crate::engine::{
    berendsen_rescale, init_velocities, Controls, Frame, System, Thermostat, Trajectory,
};
use crate::error::{ReactDynError, Result};
use crate::integrator::{kinetic_energy, velocity_verlet_step};
use crate::units::{amu_to_au_mass, fs_to_au, BOHR_PER_ANGSTROM};

/// eV → hartree (CODATA).
const EV_TO_HARTREE: f64 = 0.036_749_322;
/// Cutoff geometry (bohr): a bond tapers to zero between `r_e + MARGIN -
/// HALFWIDTH` and `r_e + MARGIN + HALFWIDTH`.
const CUTOFF_MARGIN_BOHR: f64 = 1.4;
const CUTOFF_HALFWIDTH_BOHR: f64 = 0.7;

/// Morse parameters in atomic units: well depth `d_e` (hartree), width
/// `a` (1/bohr), equilibrium bond length `r_e` (bohr).
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct MorseParam {
    /// Well depth (hartree).
    pub d_e: f64,
    /// Width parameter (1/bohr).
    pub a: f64,
    /// Equilibrium bond length (bohr).
    pub r_e: f64,
}

/// Per-element Morse parameters (converted from eV / Å). A coarse default
/// covers anything not in the small table.
pub fn morse_param(symbol: &str) -> MorseParam {
    // (D_e in eV, a in 1/Å, r_e in Å).
    let (d_ev, a_inv_ang, re_ang) = match symbol {
        "H" => (4.50, 1.90, 0.74),
        "C" => (3.60, 2.00, 1.54),
        "N" => (9.80, 2.70, 1.10),
        "O" => (5.20, 2.30, 1.21),
        "Si" => (3.20, 1.50, 2.35),
        _ => (3.00, 2.00, 1.50),
    };
    MorseParam {
        d_e: d_ev * EV_TO_HARTREE,
        a: a_inv_ang / BOHR_PER_ANGSTROM,
        r_e: re_ang * BOHR_PER_ANGSTROM,
    }
}

/// Lorentz-Berthelot-style combining of two atoms' Morse parameters.
fn combine(pi: MorseParam, pj: MorseParam) -> MorseParam {
    MorseParam {
        d_e: (pi.d_e * pj.d_e).sqrt(),
        a: 0.5 * (pi.a + pj.a),
        r_e: 0.5 * (pi.r_e + pj.r_e),
    }
}

/// Smooth cutoff `f_C(r)` and its derivative for cutoff centre `rc`,
/// half-width `d`: 1 below `rc-d`, 0 above `rc+d`, a cosine taper between.
fn cutoff(r: f64, rc: f64, d: f64) -> (f64, f64) {
    if r <= rc - d {
        (1.0, 0.0)
    } else if r >= rc + d {
        (0.0, 0.0)
    } else {
        let x = std::f64::consts::PI * (r - rc) / (2.0 * d);
        (
            0.5 - 0.5 * x.sin(),
            -0.25 * std::f64::consts::PI / d * x.cos(),
        )
    }
}

/// Morse value `V(r)` and derivative `V'(r)`.
fn morse(r: f64, p: MorseParam) -> (f64, f64) {
    let ex = (-p.a * (r - p.r_e)).exp();
    let one_minus = 1.0 - ex;
    let v = p.d_e * (one_minus * one_minus - 1.0);
    let dv = 2.0 * p.d_e * p.a * ex * one_minus;
    (v, dv)
}

/// Total energy (hartree) + forces (hartree/bohr) of the Morse+cutoff
/// reactive potential over all atom pairs. `params[i]` is atom `i`'s
/// Morse parameters; positions are in bohr.
pub fn reactive_energy_forces(
    pos_bohr: &[[f64; 3]],
    params: &[MorseParam],
) -> (f64, Vec<[f64; 3]>) {
    let n = pos_bohr.len();
    let mut energy = 0.0;
    let mut forces = vec![[0.0; 3]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let p = combine(params[i], params[j]);
            let rc = p.r_e + CUTOFF_MARGIN_BOHR;
            let dx = pos_bohr[j][0] - pos_bohr[i][0];
            let dy = pos_bohr[j][1] - pos_bohr[i][1];
            let dz = pos_bohr[j][2] - pos_bohr[i][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r < 1e-9 || r >= rc + CUTOFF_HALFWIDTH_BOHR {
                continue;
            }
            let (fc, dfc) = cutoff(r, rc, CUTOFF_HALFWIDTH_BOHR);
            let (vm, dvm) = morse(r, p);
            energy += fc * vm;
            // dV_pair/dr; force on i = (dV/dr)·r̂ with r̂ = (j - i)/r.
            let dvdr = dfc * vm + fc * dvm;
            let inv_r = 1.0 / r;
            let fx = dvdr * dx * inv_r;
            let fy = dvdr * dy * inv_r;
            let fz = dvdr * dz * inv_r;
            forces[i][0] += fx;
            forces[i][1] += fy;
            forces[i][2] += fz;
            forces[j][0] -= fx;
            forces[j][1] -= fy;
            forces[j][2] -= fz;
        }
    }
    (energy, forces)
}

/// A many-atom system for the reactive force field. No charge/multiplicity
/// — the potential is classical.
#[derive(Clone, Debug)]
pub struct ReactiveSystem {
    /// One element per atom (sets mass + Morse parameters + viewport symbol).
    pub elements: Vec<Element>,
    /// Atomic positions in bohr.
    pub pos_bohr: Vec<[f64; 3]>,
}

impl ReactiveSystem {
    /// Number of atoms.
    pub fn n_atoms(&self) -> usize {
        self.elements.len()
    }
}

/// The reactive-force-field MD engine (Morse + cutoff, classical, fast).
/// Produces the same [`Trajectory`] as the other engines, so the
/// workbench's 3-D playback + energy plot are reused unchanged.
pub struct ReactiveEngine;

impl ReactiveEngine {
    /// Run a reactive-MD trajectory. `progress(step)` fires after each step.
    pub fn run(
        &self,
        sys: &ReactiveSystem,
        controls: &Controls,
        progress: &mut dyn FnMut(usize),
    ) -> Result<Trajectory> {
        let n = sys.n_atoms();
        if n == 0 {
            return Err(ReactDynError::Invalid {
                reason: "system has no atoms".into(),
            });
        }
        if sys.pos_bohr.len() != n {
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

        let params: Vec<MorseParam> = sys.elements.iter().map(|e| morse_param(e.symbol())).collect();
        let masses: Vec<f64> = sys.elements.iter().map(|e| amu_to_au_mass(e.atomic_mass())).collect();
        let dt = fs_to_au(controls.dt_fs);

        let mut pos = sys.pos_bohr.clone();
        let force_fn = |p: &[[f64; 3]]| -> Result<Vec<[f64; 3]>> {
            Ok(reactive_energy_forces(p, &params).1)
        };
        let mut vel = init_velocities(&masses, controls.thermostat);
        let mut forces = force_fn(&pos)?;

        let traj_system = System {
            elements: sys.elements.clone(),
            pos_bohr: pos.clone(),
            charge: 0,
            multiplicity: 1,
        };

        let mut frames = Vec::with_capacity(controls.n_steps + 1);
        let (pe0, _) = reactive_energy_forces(&pos, &params);
        frames.push(Frame {
            time_fs: 0.0,
            pos_bohr: pos.clone(),
            vel: vel.clone(),
            potential_hartree: pe0,
            kinetic_hartree: kinetic_energy(&vel, &masses),
        });

        for step in 0..controls.n_steps {
            forces = velocity_verlet_step(&mut pos, &mut vel, &forces, &masses, dt, force_fn)?;
            if let Thermostat::Berendsen { target_kelvin, tau_fs } = controls.thermostat {
                berendsen_rescale(&mut vel, &masses, target_kelvin, controls.dt_fs, tau_fs);
            }
            let (pe, _) = reactive_energy_forces(&pos, &params);
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
            system: traj_system,
            frames,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cc() -> [MorseParam; 2] {
        [morse_param("C"), morse_param("C")]
    }

    fn energy_at(sep_bohr: f64) -> f64 {
        reactive_energy_forces(&[[0.0; 3], [0.0, 0.0, sep_bohr]], &cc()).0
    }

    #[test]
    fn dimer_bonds_near_equilibrium_and_dissociates_far() {
        let re = morse_param("C").r_e;
        let e_eq = energy_at(re);
        assert!(e_eq < 0.0, "should be bonded (E<0) at r_e: {e_eq}");
        // Compressed → higher energy than the equilibrium well.
        assert!(energy_at(re * 0.7) > e_eq, "compressed should cost energy");
        // Past the cutoff → exactly zero (dissociated).
        let far = re + CUTOFF_MARGIN_BOHR + CUTOFF_HALFWIDTH_BOHR + 1.0;
        assert!(energy_at(far).abs() < 1e-12, "dissociated beyond cutoff");
    }

    #[test]
    fn forces_are_restoring_and_balanced() {
        let re = morse_param("C").r_e;
        // Stretched (within cutoff) → attractive: atom 0 pulled +z toward atom 1.
        let (_, f) = reactive_energy_forces(&[[0.0; 3], [0.0, 0.0, re + 0.5]], &cc());
        assert!(f[0][2] > 0.0, "stretched: atom 0 should be pulled +z, got {}", f[0][2]);
        assert!((f[0][2] + f[1][2]).abs() < 1e-12, "Newton's third law");
        // Compressed → repulsive: atom 0 pushed -z.
        let (_, fc) = reactive_energy_forces(&[[0.0; 3], [0.0, 0.0, re - 0.5]], &cc());
        assert!(fc[0][2] < 0.0, "compressed: atom 0 should be pushed -z, got {}", fc[0][2]);
    }

    #[test]
    fn reactive_cluster_conserves_energy_nve() {
        let c = Element::from_symbol("C").unwrap();
        let re = morse_param("C").r_e;
        // A short carbon chain at the equilibrium bond spacing.
        let sys = ReactiveSystem {
            elements: vec![c, c, c, c],
            pos_bohr: (0..4).map(|i| [0.0, 0.0, i as f64 * re]).collect(),
        };
        let controls = Controls {
            method: crate::forces::Method::Rhf, // unused by the classical engine
            basis: "STO-3G".to_string(),
            dt_fs: 0.2,
            n_steps: 30,
            fd_delta_bohr: 0.01,
            thermostat: Thermostat::Nve,
            max_cost_guard: 100_000,
        };
        let traj = ReactiveEngine.run(&sys, &controls, &mut |_| {}).expect("reactive run");
        assert_eq!(traj.frames.len(), 31);
        assert_eq!(traj.system.n_atoms(), 4);
        let e0 = traj.frames[0].total_hartree();
        let drift = traj
            .frames
            .iter()
            .map(|f| (f.total_hartree() - e0).abs())
            .fold(0.0_f64, f64::max);
        assert!(drift < 1e-3, "reactive NVE drift too large: {drift} Ha");
    }
}
