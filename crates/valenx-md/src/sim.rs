//! The simulation driver — wiring a run together.
//!
//! Everything else in the crate is a building block; [`Simulation`] is
//! the type that assembles them into a runnable molecular-dynamics
//! loop. It owns:
//!
//! - a [`System`] (topology + coordinates + box),
//! - a [`ForceField`] and the [`ForceTerm`]s built from it (bonded +
//!   nonbonded),
//! - an [`Integrator`] (default velocity-Verlet),
//! - an optional [`Thermostat`] and [`Barostat`],
//! - an optional [`Constraints`] set (SHAKE/RATTLE),
//! - an [`ObservableLog`] that records the run.
//!
//! Each step the driver: integrates one time step (the integrator
//! calls back into the summed force terms), applies the constraints,
//! then applies the thermostat and barostat. [`Simulation::run`]
//! repeats that for `n` steps and returns a [`SimulationReport`].
//!
//! The nonbonded force loop is accelerated by a Verlet
//! [`NeighborList`] that is rebuilt only when an atom has drifted past
//! half the skin (see [`crate::nonbonded::neighbor`]).
//!
//! ## Building a simulation
//!
//! [`Simulation::new`] builds a sensible default from a system + force
//! field: harmonic bonded terms (if the force field supplies bonded
//! parameters), cut-off Lennard-Jones + reaction-field Coulomb for a
//! periodic system, velocity-Verlet at 1 fs, no thermostat. The
//! builder methods then customise it.

use crate::bonded::angle::HarmonicAngles;
use crate::bonded::bond::HarmonicBonds;
use crate::bonded::dihedral::ProperDihedrals;
use crate::bonded::improper::ImproperDihedrals;
use crate::bonded::{EnergyForce, ForceTerm};
use crate::constrain::Constraints;
use crate::ensemble::{Barostat, Thermostat};
use crate::error::{MdError, Result};
use crate::forcefield::ForceField;
use crate::integrate::velocity_verlet::VelocityVerlet;
use crate::integrate::Integrator;
use crate::nonbonded::coulomb::{Coulomb, CoulombMethod};
use crate::nonbonded::lj::LennardJones;
use crate::nonbonded::neighbor::NeighborList;
use crate::analysis::reporters::{state_report, ObservableLog, StateReport};
use crate::system::System;

/// Which nonbonded interactions a simulation evaluates pairwise
/// through the neighbour list.
struct NonbondedTerms {
    lj: Option<LennardJones>,
    coulomb: Option<Coulomb>,
    /// Cutoff used to size the neighbour list (nm).
    cutoff: f64,
    /// Verlet skin (nm).
    skin: f64,
}

/// A runnable molecular-dynamics simulation.
pub struct Simulation {
    /// The system being simulated.
    pub system: System,
    /// The force field (kept for reference / re-derivation).
    pub force_field: ForceField,
    /// Bonded force terms summed every step.
    bonded: Vec<Box<dyn ForceTerm>>,
    /// Nonbonded terms (use the neighbour list).
    nonbonded: Option<NonbondedTerms>,
    /// The cached neighbour list, rebuilt on drift.
    neighbor_list: Option<NeighborList>,
    /// The time integrator.
    integrator: Box<dyn Integrator>,
    /// Optional temperature coupling.
    thermostat: Option<Box<dyn Thermostat>>,
    /// Optional pressure coupling.
    barostat: Option<Box<dyn Barostat>>,
    /// Optional holonomic constraints.
    constraints: Option<Constraints>,
    /// The thermodynamic log.
    pub log: ObservableLog,
    /// Steps completed so far.
    step_count: usize,
    /// Simulation time so far (ps).
    time: f64,
    /// Record an observable every `report_interval` steps.
    report_interval: usize,
}

/// The summary returned by [`Simulation::run`].
#[derive(Clone, Debug, PartialEq)]
pub struct SimulationReport {
    /// Steps performed in this `run` call.
    pub steps: usize,
    /// Potential energy at the end (kJ/mol).
    pub final_potential_energy: f64,
    /// Kinetic energy at the end (kJ/mol).
    pub final_kinetic_energy: f64,
    /// Total energy at the end (kJ/mol).
    pub final_total_energy: f64,
    /// Temperature at the end (K).
    pub final_temperature: f64,
    /// Simulation time at the end (ps).
    pub final_time: f64,
}

impl Simulation {
    /// Builds a default simulation from a system and a force field.
    ///
    /// Bonded terms are added for whichever bonded parameter lists the
    /// force field supplies. For a periodic system, cut-off
    /// Lennard-Jones and conductor-reaction-field Coulomb are added
    /// with a 1.0 nm cutoff and a 0.1 nm Verlet skin. The integrator
    /// is velocity-Verlet at a 1 fs (0.001 ps) step.
    ///
    /// # Errors
    /// [`MdError::DimensionMismatch`] if the force field does not
    /// cover the system (propagated from
    /// [`ForceField::validate_against`]); propagates term-construction
    /// errors.
    pub fn new(system: System, force_field: ForceField) -> Result<Self> {
        force_field.validate_against(&system.topology)?;

        // --- Bonded terms ------------------------------------------
        let mut bonded: Vec<Box<dyn ForceTerm>> = Vec::new();
        if !force_field.bonds().is_empty() {
            bonded.push(Box::new(HarmonicBonds::from_system(
                &system,
                force_field.bonds(),
            )?));
        }
        if !force_field.angles().is_empty() {
            bonded.push(Box::new(HarmonicAngles::from_system(
                &system,
                force_field.angles(),
            )?));
        }
        if !force_field.dihedrals().is_empty() {
            bonded.push(Box::new(ProperDihedrals::from_system(
                &system,
                force_field.dihedrals(),
            )?));
        }
        if !force_field.impropers().is_empty() {
            bonded.push(Box::new(ImproperDihedrals::from_system(
                &system,
                force_field.impropers(),
            )?));
        }

        // --- Nonbonded terms ---------------------------------------
        let cutoff = 1.0;
        let skin = 0.1;
        let nonbonded = if system.cell.is_periodic()
            && system.cell.max_cutoff() > cutoff
            && system.len() >= 2
        {
            let lj = LennardJones::from_system(&system, &force_field, cutoff)?;
            let coulomb = Coulomb::from_system(
                &system,
                cutoff,
                CoulombMethod::conductor_reaction_field(),
            )?;
            Some(NonbondedTerms {
                lj: Some(lj),
                coulomb: Some(coulomb),
                cutoff,
                skin,
            })
        } else {
            None
        };

        Ok(Simulation {
            system,
            force_field,
            bonded,
            nonbonded,
            neighbor_list: None,
            integrator: Box::new(VelocityVerlet::new(0.001)?),
            thermostat: None,
            barostat: None,
            constraints: None,
            log: ObservableLog::new(),
            step_count: 0,
            time: 0.0,
            report_interval: 1,
        })
    }

    /// Replaces the integrator.
    pub fn with_integrator(mut self, integrator: Box<dyn Integrator>) -> Self {
        self.integrator = integrator;
        self
    }

    /// Attaches a thermostat (NVT ensemble).
    pub fn with_thermostat(mut self, thermostat: Box<dyn Thermostat>) -> Self {
        self.thermostat = Some(thermostat);
        self
    }

    /// Attaches a barostat (NPT ensemble — usually paired with a
    /// thermostat).
    pub fn with_barostat(mut self, barostat: Box<dyn Barostat>) -> Self {
        self.barostat = Some(barostat);
        self
    }

    /// Attaches a SHAKE/RATTLE constraint set.
    pub fn with_constraints(mut self, constraints: Constraints) -> Self {
        self.constraints = Some(constraints);
        self
    }

    /// Sets how often (in steps) a [`StateReport`] is logged.
    ///
    /// # Errors
    /// [`MdError::Invalid`] if `interval` is zero.
    pub fn set_report_interval(&mut self, interval: usize) -> Result<()> {
        if interval == 0 {
            return Err(MdError::invalid("report_interval", "must be at least 1"));
        }
        self.report_interval = interval;
        Ok(())
    }

    /// Total steps completed across all `run` calls.
    pub fn step_count(&self) -> usize {
        self.step_count
    }

    /// Simulation time elapsed (ps).
    pub fn time(&self) -> f64 {
        self.time
    }

    /// The number of holonomic constraints (0 if none) — used for the
    /// degree-of-freedom count.
    fn constraint_count(&self) -> usize {
        self.constraints.as_ref().map(|c| c.len()).unwrap_or(0)
    }

    /// Evaluates the total energy + forces of the current system.
    ///
    /// Sums every bonded term and, if present, the nonbonded terms
    /// over the (lazily rebuilt) neighbour list.
    ///
    /// # Errors
    /// Propagates any term's evaluation error.
    pub fn evaluate_forces(&mut self) -> Result<EnergyForce> {
        // Rebuild the neighbour list if needed.
        if let Some(nb) = &self.nonbonded {
            let rebuild = match &self.neighbor_list {
                None => true,
                Some(nl) => nl.needs_rebuild(&self.system.positions, &self.system.cell),
            };
            if rebuild {
                self.neighbor_list = Some(NeighborList::build(
                    &self.system.positions,
                    &self.system.cell,
                    nb.cutoff,
                    nb.skin,
                )?);
            }
        }
        Self::sum_forces(
            &self.system,
            &self.bonded,
            self.nonbonded.as_ref(),
            self.neighbor_list.as_ref(),
        )
    }

    /// Sums all force terms for a system (a free function so the
    /// integrator closure can call it without borrowing `self`).
    fn sum_forces(
        system: &System,
        bonded: &[Box<dyn ForceTerm>],
        nonbonded: Option<&NonbondedTerms>,
        neighbor_list: Option<&NeighborList>,
    ) -> Result<EnergyForce> {
        let mut ef = EnergyForce::zeros(system.len());
        for term in bonded {
            term.accumulate(system, &mut ef)?;
        }
        if let (Some(nb), Some(nl)) = (nonbonded, neighbor_list) {
            if let Some(lj) = &nb.lj {
                lj.accumulate_pairs(system, nl.pairs(), &mut ef)?;
            }
            if let Some(coulomb) = &nb.coulomb {
                coulomb.accumulate_pairs(system, nl.pairs(), &mut ef)?;
            }
        }
        Ok(ef)
    }

    /// Advances the simulation by one step.
    ///
    /// # Errors
    /// Propagates integrator / force / constraint / thermostat errors.
    pub fn step(&mut self) -> Result<StateReport> {
        // Make sure the neighbour list exists / is current before the
        // integrator's force callbacks run.
        if let Some(nb) = &self.nonbonded {
            let rebuild = match &self.neighbor_list {
                None => true,
                Some(nl) => nl.needs_rebuild(&self.system.positions, &self.system.cell),
            };
            if rebuild {
                self.neighbor_list = Some(NeighborList::build(
                    &self.system.positions,
                    &self.system.cell,
                    nb.cutoff,
                    nb.skin,
                )?);
            }
        }
        let dt = self.integrator.dt();
        // Positions before the step — needed by SHAKE.
        let reference = self.system.positions.clone();

        // Integrate. The closure sums forces with the current
        // neighbour list (captured by reference, not via `self`).
        let bonded = &self.bonded;
        let nonbonded = self.nonbonded.as_ref();
        let neighbor_list = self.neighbor_list.as_ref();
        let mut force_fn = |s: &System| {
            Self::sum_forces(s, bonded, nonbonded, neighbor_list)
        };
        let ef = self.integrator.step(&mut self.system, &mut force_fn)?;

        // Apply constraints (SHAKE positions + RATTLE velocities).
        if let Some(constraints) = &self.constraints {
            constraints.shake(&mut self.system, &reference, Some(dt))?;
            constraints.rattle_velocities(&mut self.system)?;
        }

        let nconstraints = self.constraint_count();

        // Thermostat.
        if let Some(thermostat) = &mut self.thermostat {
            thermostat.apply(&mut self.system, dt, nconstraints)?;
        }

        // Barostat (needs the current virial pressure).
        if let Some(barostat) = &mut self.barostat {
            if let Some(pressure) =
                crate::analysis::reporters::pressure_bar(&self.system, &ef)
            {
                barostat.apply(&mut self.system, pressure, dt)?;
                // The box changed: force a neighbour-list rebuild.
                self.neighbor_list = None;
            }
        }

        self.step_count += 1;
        self.time += dt;
        Ok(state_report(
            &self.system,
            &ef,
            self.step_count,
            self.time,
            nconstraints,
        ))
    }

    /// Runs `n` steps, logging an observable every `report_interval`
    /// steps, and returns a [`SimulationReport`].
    ///
    /// # Errors
    /// Propagates any per-step error.
    pub fn run(&mut self, n: usize) -> Result<SimulationReport> {
        let mut last = state_report(
            &self.system,
            &EnergyForce::zeros(self.system.len()),
            self.step_count,
            self.time,
            self.constraint_count(),
        );
        for s in 0..n {
            last = self.step()?;
            if s % self.report_interval == 0 {
                self.log.record(last);
            }
        }
        Ok(SimulationReport {
            steps: n,
            final_potential_energy: last.potential_energy,
            final_kinetic_energy: last.kinetic_energy,
            final_total_energy: last.total_energy,
            final_temperature: last.temperature,
            final_time: last.time,
        })
    }

    /// Computes the current potential energy without advancing the
    /// simulation.
    ///
    /// # Errors
    /// Propagates force-evaluation errors.
    pub fn potential_energy(&mut self) -> Result<f64> {
        Ok(self.evaluate_forces()?.energy)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ensemble::berendsen::Berendsen;
    use crate::forcefield::{BondParam, CombiningRule, LjParam};
    use crate::integrate::langevin::Langevin;
    use crate::pbc::SimBox;
    use crate::system::{Atom, Topology};
    use nalgebra::Vector3;

    /// A small periodic argon-like system.
    fn argon_box(n_per_side: usize) -> (System, ForceField) {
        let spacing = 0.4;
        let edge = spacing * n_per_side as f64 + 0.4;
        let mut top = Topology::new();
        let mut pos = Vec::new();
        for i in 0..n_per_side {
            for j in 0..n_per_side {
                for k in 0..n_per_side {
                    top.push_atom(Atom::new("Ar", 39.95, 0.0).unwrap());
                    pos.push(Vector3::new(
                        i as f64 * spacing + 0.2,
                        j as f64 * spacing + 0.2,
                        k as f64 * spacing + 0.2,
                    ));
                }
            }
        }
        let sys = System::new(top, pos)
            .unwrap()
            .with_cell(SimBox::cubic(edge).unwrap());
        let mut ff = ForceField::new(CombiningRule::LorentzBerthelot);
        ff.set_lj("Ar", LjParam::new(0.34, 0.996).unwrap());
        (sys, ff)
    }

    #[test]
    fn new_rejects_incomplete_force_field() {
        let (sys, _) = argon_box(2);
        let empty_ff = ForceField::new(CombiningRule::LorentzBerthelot);
        // No LJ params for "Ar" -> error.
        assert!(Simulation::new(sys, empty_ff).is_err());
    }

    #[test]
    fn runs_an_nve_simulation() {
        let (sys, ff) = argon_box(3);
        let mut sim = Simulation::new(sys, ff).unwrap();
        let report = sim.run(50).unwrap();
        assert_eq!(report.steps, 50);
        assert_eq!(sim.step_count(), 50);
        assert!(report.final_total_energy.is_finite());
        assert!(sim.time() > 0.0);
    }

    #[test]
    fn energy_is_roughly_conserved_in_nve() {
        // A modest argon box with the velocity-Verlet integrator
        // should keep the total energy in a bounded band.
        let (sys, ff) = argon_box(3);
        let mut sim = Simulation::new(sys, ff).unwrap();
        sim.run(200).unwrap();
        let e_std = sim.log.total_energy_std();
        let e_mean = sim.log.mean_total_energy().abs().max(1.0);
        assert!(e_std < 0.2 * e_mean, "energy drift std {e_std} / mean {e_mean}");
    }

    #[test]
    fn thermostat_drives_temperature() {
        // Start cold, attach a Berendsen thermostat at 200 K.
        let (mut sys, ff) = argon_box(3);
        // Give a small initial temperature.
        sys.set_velocities(
            (0..sys.len())
                .map(|i| {
                    let s = if i % 2 == 0 { 0.1 } else { -0.1 };
                    Vector3::new(s, 0.0, 0.0)
                })
                .collect(),
        )
        .unwrap();
        let mut sim = Simulation::new(sys, ff)
            .unwrap()
            .with_thermostat(Box::new(Berendsen::new(200.0, 0.1).unwrap()));
        sim.run(2000).unwrap();
        let t = sim.system.temperature(0);
        // Should be pulled toward 200 K (loose tolerance — coupled to
        // the LJ dynamics).
        assert!(t > 80.0 && t < 350.0, "thermostatted T = {t}");
    }

    #[test]
    fn runs_with_a_langevin_integrator() {
        let (sys, ff) = argon_box(2);
        let langevin = Langevin::new(0.002, 2.0, 150.0, 42).unwrap();
        let mut sim = Simulation::new(sys, ff)
            .unwrap()
            .with_integrator(Box::new(langevin));
        let report = sim.run(100).unwrap();
        assert_eq!(report.steps, 100);
        assert!(report.final_temperature.is_finite());
    }

    #[test]
    fn bonded_only_system_runs_without_a_box() {
        // A non-periodic diatomic with only a harmonic bond.
        let mut top = Topology::new();
        top.push_atom(Atom::new("A", 12.0, 0.0).unwrap());
        top.push_atom(Atom::new("B", 12.0, 0.0).unwrap());
        top.add_bond(0, 1).unwrap();
        let sys = System::new(
            top,
            vec![Vector3::zeros(), Vector3::new(0.18, 0.0, 0.0)],
        )
        .unwrap();
        let mut ff = ForceField::new(CombiningRule::LorentzBerthelot);
        ff.set_lj("A", LjParam::new(0.3, 0.5).unwrap());
        ff.set_lj("B", LjParam::new(0.3, 0.5).unwrap());
        ff.push_bond(BondParam::new(0.15, 3000.0).unwrap());
        let mut sim = Simulation::new(sys, ff).unwrap();
        let report = sim.run(100).unwrap();
        assert!(report.final_potential_energy.is_finite());
    }

    #[test]
    fn report_interval_controls_log_density() {
        let (sys, ff) = argon_box(2);
        let mut sim = Simulation::new(sys, ff).unwrap();
        sim.set_report_interval(10).unwrap();
        sim.run(100).unwrap();
        // 100 steps, recorded every 10 -> 10 entries.
        assert_eq!(sim.log.len(), 10);
        assert!(sim.set_report_interval(0).is_err());
    }

    #[test]
    fn potential_energy_query_does_not_advance() {
        let (sys, ff) = argon_box(2);
        let mut sim = Simulation::new(sys, ff).unwrap();
        let _ = sim.potential_energy().unwrap();
        assert_eq!(sim.step_count(), 0);
    }
}
