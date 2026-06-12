//! The **QM/MM engine** — a quantum reacting region embedded in an
//! explicit classical environment ([`crate::mm`]).
//!
//! Two embedding schemes:
//! - **Mechanical** ([`Embedding::Mechanical`]): the QM region is solved
//!   in isolation; the QM–MM interaction is purely classical LJ + Coulomb.
//!   Cheap — QM forces by FD over the QM region only, MM forces analytic.
//! - **Electrostatic** ([`Embedding::Electrostatic`]): the MM point
//!   charges enter the QM SCF ([`single_point_energy_embedded`]) and
//!   **polarize the density** — the accurate scheme. The total energy is
//!   the embedded QM energy + the classical MM–MM and QM–MM-LJ terms +
//!   the classical QM-nuclei–MM-charge term; forces are the finite-
//!   difference gradient of that total over *all* atoms (so they are the
//!   true gradient and the dynamics conserve energy). v1 is RHF-only.
//!
//! Both produce the same [`Trajectory`] (QM elements then MM elements),
//! so the workbench's 3-D playback + energy plot are reused unchanged.
//! Everything is in atomic units.

use valenx_qchem::element::Element;

use crate::engine::{
    berendsen_rescale, init_velocities, Controls, Frame, System, Thermostat, Trajectory,
};
use crate::error::{ReactDynError, Result};
use crate::forces::{numerical_forces, single_point_energy, single_point_energy_embedded, Method};
use crate::integrator::{kinetic_energy, velocity_verlet_step};
use crate::mm::{classical_forces, Particle};
use crate::units::{amu_to_au_mass, fs_to_au};

/// QM/MM coupling scheme.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub enum Embedding {
    /// Classical LJ + Coulomb coupling; the MM environment does not
    /// polarize the QM density.
    #[default]
    Mechanical,
    /// MM charges enter the QM SCF and polarize the density (RHF, v1).
    Electrostatic,
}

/// One explicit MM environment atom (atomic units; `element` supplies the
/// mass + the viewport symbol).
#[derive(Clone, Debug)]
pub struct MmAtom {
    /// Chemical element (mass + display symbol).
    pub element: Element,
    /// Position in bohr.
    pub pos_bohr: [f64; 3],
    /// Partial charge in e.
    pub charge: f64,
    /// Lennard-Jones σ in bohr.
    pub sigma_bohr: f64,
    /// Lennard-Jones ε in hartree.
    pub epsilon_hartree: f64,
    /// Mass in amu.
    pub mass_amu: f64,
}

/// A QM/MM system: a quantum reacting region + an explicit classical
/// environment + the embedding scheme.
#[derive(Clone, Debug)]
pub struct QmMmSystem {
    /// QM-region elements.
    pub qm_elements: Vec<Element>,
    /// QM-region positions (bohr).
    pub qm_pos_bohr: Vec<[f64; 3]>,
    /// QM-region total charge.
    pub qm_charge: i32,
    /// QM-region spin multiplicity.
    pub qm_mult: u32,
    /// Classical params per QM atom for the coupling: `(charge, σ_bohr, ε_hartree)`.
    pub qm_classical: Vec<(f64, f64, f64)>,
    /// The explicit MM environment.
    pub mm: Vec<MmAtom>,
    /// The embedding scheme.
    pub embedding: Embedding,
}

impl QmMmSystem {
    /// Number of QM atoms.
    pub fn n_qm(&self) -> usize {
        self.qm_elements.len()
    }
    /// Number of MM atoms.
    pub fn n_mm(&self) -> usize {
        self.mm.len()
    }
}

/// Build the classical [`Particle`] lists for the QM and MM atoms.
/// `include_qm_charge` is `true` for mechanical embedding (QM atoms carry
/// their coupling charge) and `false` for electrostatic embedding (the
/// QM–MM Coulomb is handled inside the SCF, so the classical QM charge is
/// zeroed — only the QM–MM LJ remains classical).
fn build_particles(
    qm_pos: &[[f64; 3]],
    mm_pos: &[[f64; 3]],
    sys: &QmMmSystem,
    include_qm_charge: bool,
) -> (Vec<Particle>, Vec<Particle>) {
    let qm = qm_pos
        .iter()
        .zip(&sys.qm_classical)
        .map(|(p, &(charge, sigma_bohr, epsilon_hartree))| Particle {
            pos_bohr: *p,
            charge: if include_qm_charge { charge } else { 0.0 },
            sigma_bohr,
            epsilon_hartree,
        })
        .collect();
    let mm = mm_pos
        .iter()
        .zip(&sys.mm)
        .map(|(p, m)| Particle {
            pos_bohr: *p,
            charge: m.charge,
            sigma_bohr: m.sigma_bohr,
            epsilon_hartree: m.epsilon_hartree,
        })
        .collect();
    (qm, mm)
}

/// The classical QM-nuclei ↔ MM-charge Coulomb energy (hartree). Used in
/// electrostatic embedding, where the QM-electron ↔ MM-charge part lives
/// in the SCF but the nuclei part stays classical.
fn nuclei_charge_energy(
    qm_elements: &[Element],
    qm_pos: &[[f64; 3]],
    mm_pos: &[[f64; 3]],
    mm: &[MmAtom],
) -> f64 {
    let mut e = 0.0;
    for (a, ra) in qm_pos.iter().enumerate() {
        let z = qm_elements[a].nuclear_charge();
        for (m, rm) in mm_pos.iter().enumerate() {
            let dx = ra[0] - rm[0];
            let dy = ra[1] - rm[1];
            let dz = ra[2] - rm[2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r > 1e-9 {
                e += z * mm[m].charge / r;
            }
        }
    }
    e
}

/// Total potential energy of a QM/MM configuration (`all` = QM positions
/// then MM positions), branching on the embedding scheme.
fn qmmm_total_energy(
    sys: &QmMmSystem,
    all: &[[f64; 3]],
    method: Method,
    basis: &str,
) -> Result<f64> {
    let n_qm = sys.n_qm();
    let qm_pos = &all[..n_qm];
    let mm_pos = &all[n_qm..];
    match sys.embedding {
        Embedding::Mechanical => {
            let e_qm = single_point_energy(
                &sys.qm_elements,
                qm_pos,
                sys.qm_charge,
                sys.qm_mult,
                method,
                basis,
            )?;
            let (qm_p, mm_p) = build_particles(qm_pos, mm_pos, sys, true);
            let (e_cl, _, _) = classical_forces(&qm_p, &mm_p);
            Ok(e_qm + e_cl)
        }
        Embedding::Electrostatic => {
            let ext: Vec<(f64, [f64; 3])> = mm_pos
                .iter()
                .zip(&sys.mm)
                .map(|(p, m)| (m.charge, *p))
                .collect();
            let e_qm = single_point_energy_embedded(
                &sys.qm_elements,
                qm_pos,
                sys.qm_charge,
                sys.qm_mult,
                basis,
                &ext,
            )?;
            // Classical: MM–MM (LJ + Coulomb) + QM–MM LJ only (QM charge
            // zeroed — the QM–MM Coulomb is in the SCF).
            let (qm_p, mm_p) = build_particles(qm_pos, mm_pos, sys, false);
            let (e_cl, _, _) = classical_forces(&qm_p, &mm_p);
            let e_nuc = nuclei_charge_energy(&sys.qm_elements, qm_pos, mm_pos, &sys.mm);
            Ok(e_qm + e_cl + e_nuc)
        }
    }
}

/// Central finite-difference gradient of `energy` over every coordinate
/// of `pos`, returned as forces (= −∇E). Used for electrostatic
/// embedding, where the MM forces depend on the QM density.
fn numerical_gradient(
    pos: &[[f64; 3]],
    delta: f64,
    energy: impl Fn(&[[f64; 3]]) -> Result<f64>,
) -> Result<Vec<[f64; 3]>> {
    let n = pos.len();
    let mut forces = vec![[0.0; 3]; n];
    let mut p = pos.to_vec();
    for i in 0..n {
        for d in 0..3 {
            let orig = p[i][d];
            p[i][d] = orig + delta;
            let ep = energy(&p)?;
            p[i][d] = orig - delta;
            let em = energy(&p)?;
            p[i][d] = orig;
            forces[i][d] = -(ep - em) / (2.0 * delta);
        }
    }
    Ok(forces)
}

/// Forces on all atoms (QM then MM), branching on the embedding scheme.
fn qmmm_forces(
    sys: &QmMmSystem,
    all: &[[f64; 3]],
    method: Method,
    basis: &str,
    delta: f64,
) -> Result<Vec<[f64; 3]>> {
    let n_qm = sys.n_qm();
    match sys.embedding {
        Embedding::Mechanical => {
            // Cheap: qchem forces on the QM region + classical coupling.
            let qm_pos = &all[..n_qm];
            let mm_pos = &all[n_qm..];
            let mut f = numerical_forces(
                &sys.qm_elements,
                qm_pos,
                sys.qm_charge,
                sys.qm_mult,
                method,
                basis,
                delta,
            )?;
            let (qm_p, mm_p) = build_particles(qm_pos, mm_pos, sys, true);
            let (_e, f_qm_cl, f_mm_cl) = classical_forces(&qm_p, &mm_p);
            for (i, c) in f_qm_cl.iter().enumerate() {
                f[i][0] += c[0];
                f[i][1] += c[1];
                f[i][2] += c[2];
            }
            f.extend(f_mm_cl);
            Ok(f)
        }
        Embedding::Electrostatic => {
            // The true gradient of the total (polarized) energy over all
            // atoms — so the MM atoms feel the QM density. Expensive.
            numerical_gradient(all, delta, |p| qmmm_total_energy(sys, p, method, basis))
        }
    }
}

/// QM/MM dynamics engine (mechanical or electrostatic embedding).
pub struct QmMmEngine;

impl QmMmEngine {
    /// Run a QM/MM trajectory. `progress(step)` fires after each step.
    pub fn run(
        &self,
        sys: &QmMmSystem,
        controls: &Controls,
        progress: &mut dyn FnMut(usize),
    ) -> Result<Trajectory> {
        let n_qm = sys.n_qm();
        let n_mm = sys.n_mm();
        let n = n_qm + n_mm;
        if n_qm == 0 {
            return Err(ReactDynError::Invalid {
                reason: "QM region has no atoms".into(),
            });
        }
        if sys.qm_pos_bohr.len() != n_qm || sys.qm_classical.len() != n_qm {
            return Err(ReactDynError::Invalid {
                reason: "QM elements / positions / classical-params length mismatch".into(),
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
        if sys.embedding == Embedding::Electrostatic && controls.method != Method::Rhf {
            return Err(ReactDynError::Invalid {
                reason: "electrostatic embedding supports RHF only in v1".into(),
            });
        }

        // Cost guard: mechanical does FD over the QM region only;
        // electrostatic does FD over all atoms, so it scales with `n`.
        let fd_atoms = match sys.embedding {
            Embedding::Mechanical => n_qm,
            Embedding::Electrostatic => n,
        };
        let cost = fd_atoms.saturating_mul(controls.n_steps);
        if cost > controls.max_cost_guard {
            return Err(ReactDynError::GuardExceeded {
                atoms: fd_atoms,
                steps: controls.n_steps,
                cap: controls.max_cost_guard,
            });
        }

        let mut masses = Vec::with_capacity(n);
        for e in &sys.qm_elements {
            masses.push(amu_to_au_mass(e.atomic_mass()));
        }
        for m in &sys.mm {
            masses.push(amu_to_au_mass(m.mass_amu));
        }

        let dt = fs_to_au(controls.dt_fs);
        let method = controls.method;
        let basis = controls.basis.as_str();
        let delta = controls.fd_delta_bohr;

        let mut pos: Vec<[f64; 3]> = Vec::with_capacity(n);
        pos.extend_from_slice(&sys.qm_pos_bohr);
        for m in &sys.mm {
            pos.push(m.pos_bohr);
        }

        let force_fn = |all: &[[f64; 3]]| -> Result<Vec<[f64; 3]>> {
            qmmm_forces(sys, all, method, basis, delta)
        };

        let mut vel = init_velocities(&masses, controls.thermostat);
        let mut forces = force_fn(&pos)?;

        let combined_elements: Vec<Element> = sys
            .qm_elements
            .iter()
            .copied()
            .chain(sys.mm.iter().map(|m| m.element))
            .collect();
        let traj_system = System {
            elements: combined_elements,
            pos_bohr: pos.clone(),
            charge: sys.qm_charge,
            multiplicity: sys.qm_mult,
        };

        let mut frames = Vec::with_capacity(controls.n_steps + 1);
        let pe0 = qmmm_total_energy(sys, &pos, method, basis)?;
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
            let pe = qmmm_total_energy(sys, &pos, method, basis)?;
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

    fn h() -> Element {
        Element::from_symbol("H").unwrap()
    }

    fn h2_with_mm(embedding: Embedding, mm_charge: f64) -> QmMmSystem {
        QmMmSystem {
            qm_elements: vec![h(), h()],
            qm_pos_bohr: vec![[0.0, 0.0, 0.0], [0.0, 0.0, 1.4]],
            qm_charge: 0,
            qm_mult: 1,
            qm_classical: vec![(0.0, 2.0, 0.005), (0.0, 2.0, 0.005)],
            mm: vec![
                MmAtom {
                    element: h(),
                    pos_bohr: [5.0, 0.0, 0.0],
                    charge: mm_charge,
                    sigma_bohr: 3.0,
                    epsilon_hartree: 0.001,
                    mass_amu: 20.0,
                },
                MmAtom {
                    element: h(),
                    pos_bohr: [-5.0, 0.0, 0.7],
                    charge: -mm_charge,
                    sigma_bohr: 3.0,
                    epsilon_hartree: 0.001,
                    mass_amu: 20.0,
                },
            ],
            embedding,
        }
    }

    fn controls(n_steps: usize) -> Controls {
        Controls {
            method: Method::Rhf,
            basis: "STO-3G".to_string(),
            dt_fs: 0.5,
            n_steps,
            fd_delta_bohr: 0.01,
            thermostat: Thermostat::Nve,
            max_cost_guard: 8000,
        }
    }

    #[test]
    fn mechanical_runs_and_builds_combined_system() {
        let sys = h2_with_mm(Embedding::Mechanical, 0.0);
        let traj = QmMmEngine
            .run(&sys, &controls(6), &mut |_| {})
            .expect("run");
        assert_eq!(traj.frames.len(), 7);
        assert_eq!(traj.system.n_atoms(), 4);
    }

    #[test]
    fn mechanical_nve_conserves_total_energy() {
        let sys = h2_with_mm(Embedding::Mechanical, 0.0);
        let traj = QmMmEngine.run(&sys, &controls(10), &mut |_| {}).unwrap();
        let e0 = traj.frames[0].total_hartree();
        let drift = traj
            .frames
            .iter()
            .map(|f| (f.total_hartree() - e0).abs())
            .fold(0.0_f64, f64::max);
        assert!(drift < 1e-2, "mechanical NVE drift too large: {drift} Ha");
    }

    #[test]
    fn electrostatic_runs_and_conserves_energy() {
        // Charged MM atoms so the embedding electrostatics matter.
        let sys = h2_with_mm(Embedding::Electrostatic, 0.4);
        let traj = QmMmEngine
            .run(&sys, &controls(6), &mut |_| {})
            .expect("electrostatic run");
        assert_eq!(traj.system.n_atoms(), 4);
        let e0 = traj.frames[0].total_hartree();
        let drift = traj
            .frames
            .iter()
            .map(|f| (f.total_hartree() - e0).abs())
            .fold(0.0_f64, f64::max);
        // FD-of-total-energy forces are the true gradient → energy conserves.
        assert!(
            drift < 2e-2,
            "electrostatic NVE drift too large: {drift} Ha"
        );
    }

    #[test]
    fn electrostatic_rejects_non_rhf() {
        let sys = h2_with_mm(Embedding::Electrostatic, 0.4);
        let c = Controls {
            method: Method::Uhf,
            ..controls(4)
        };
        assert!(matches!(
            QmMmEngine.run(&sys, &c, &mut |_| {}),
            Err(ReactDynError::Invalid { .. })
        ));
    }

    #[test]
    fn empty_qm_region_fails_loud() {
        let sys = QmMmSystem {
            qm_elements: vec![],
            qm_pos_bohr: vec![],
            qm_charge: 0,
            qm_mult: 1,
            qm_classical: vec![],
            mm: vec![],
            embedding: Embedding::Mechanical,
        };
        assert!(QmMmEngine.run(&sys, &controls(4), &mut |_| {}).is_err());
    }
}
