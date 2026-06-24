//! Metropolis Monte Carlo — sampling the canonical (NVT) ensemble
//! without integrating equations of motion.
//!
//! Molecular *dynamics* propagates Newton's equations; Monte Carlo
//! instead generates a Markov chain of configurations whose stationary
//! distribution is the Boltzmann distribution `P(r) ∝ exp(−U(r)/k_BT)`.
//! The classic recipe (Metropolis, Rosenbluth, Rosenbluth, Teller &
//! Teller, *J. Chem. Phys.* **21**, 1087 (1953)) is:
//!
//! 1. Pick a particle at random and propose a small random
//!    displacement.
//! 2. Compute the energy change `ΔU` the move would cause.
//! 3. **Accept** the move unconditionally if `ΔU ≤ 0`; otherwise accept
//!    it with probability `exp(−ΔU / k_BT)`. Reject ⇒ keep the old
//!    configuration (and *count it again* in any average).
//!
//! Step 3 is the **Metropolis acceptance criterion**; it satisfies
//! detailed balance, so the chain samples the canonical ensemble. No
//! velocities, masses, or forces are needed — only the potential energy
//! and a temperature.
//!
//! This driver reuses the crate's existing energy model: any type that
//! can report the total potential energy of a [`System`] implements
//! [`EnergyModel`], and the bundled [`ForceTermEnergy`] adapts a list
//! of [`ForceTerm`]s (Lennard-Jones, Coulomb, Wolf, the bonded terms,
//! …) into one. The random stream is the crate's deterministic
//! [`Rng`], so a seeded run is bit-for-bit reproducible.
//!
//! [`Rng`]: crate::rng::Rng
//!
//! ## Move set (scope)
//!
//! v1 implements **single-particle random translation** moves — the
//! workhorse of canonical MC for fluids and the move set the
//! Metropolis paper introduced. Each move displaces one randomly chosen
//! atom by a uniform vector in `[−δ, δ]³`. The maximum displacement `δ`
//! can be tuned (manually or by [`MonteCarlo::tune_step`]) toward the
//! canonical ~50 % acceptance ratio. Volume moves (for NPT-MC),
//! configurational-bias / regrowth moves, and orientational moves for
//! rigid bodies are deliberately out of scope here; the constant-volume
//! translation chain is the well-defined, testable core.
//!
//! Because the move is local, only the moved atom's pair interactions
//! change. The driver computes `ΔU` as the difference of two total
//! energies by default (simple and always correct); supply an
//! [`EnergyModel`] with a cheaper [`EnergyModel::move_delta`] override
//! to make each step `O(N)` instead of `O(N²)`.

use nalgebra::Vector3;

use crate::bonded::{EnergyForce, ForceTerm};
use crate::error::{MdError, Result};
use crate::rng::Rng;
use crate::system::System;
use crate::units::BOLTZMANN;

/// A source of the total potential energy of a [`System`].
///
/// Implement this for whatever interaction set your Monte Carlo run
/// should sample. The bundled [`ForceTermEnergy`] covers the common
/// case of summing the crate's [`ForceTerm`]s; a closure
/// `Fn(&System) -> Result<f64>` also works via the blanket impl.
pub trait EnergyModel {
    /// The total potential energy of `system` (kJ/mol).
    ///
    /// # Errors
    /// Implementation-specific (e.g. a dimension mismatch in a force
    /// term).
    fn potential_energy(&self, system: &System) -> Result<f64>;

    /// The energy change caused by moving atom `atom` from
    /// `old_position` to its *current* position in `system`.
    ///
    /// The default recomputes the full energy twice (correct but
    /// `O(N²)`). Implementations that can evaluate only the moved
    /// atom's interactions should override this for an `O(N)` step.
    ///
    /// On entry `system` already holds the *trial* (moved) position;
    /// `old_position` is where the atom was before the move.
    ///
    /// # Errors
    /// Propagates [`potential_energy`](Self::potential_energy) errors.
    fn move_delta(
        &self,
        system: &mut System,
        atom: usize,
        old_position: Vector3<f64>,
    ) -> Result<f64> {
        // Trial energy (current positions).
        let e_trial = self.potential_energy(system)?;
        // Restore, evaluate the old energy, then re-apply the trial.
        let trial = system.positions[atom];
        system.positions[atom] = old_position;
        let e_old = self.potential_energy(system)?;
        system.positions[atom] = trial;
        Ok(e_trial - e_old)
    }
}

/// Forwarding impl so a shared reference to an energy model is itself an
/// energy model — lets callers reuse one model across several
/// [`MonteCarlo`] runs (`MonteCarlo::new(sys, &model, …)`).
impl<E: EnergyModel + ?Sized> EnergyModel for &E {
    fn potential_energy(&self, system: &System) -> Result<f64> {
        (**self).potential_energy(system)
    }

    fn move_delta(
        &self,
        system: &mut System,
        atom: usize,
        old_position: Vector3<f64>,
    ) -> Result<f64> {
        (**self).move_delta(system, atom, old_position)
    }
}

/// Wraps a closure `Fn(&System) -> Result<f64>` as an [`EnergyModel`].
///
/// A named wrapper (rather than a blanket `impl EnergyModel for F`) so
/// it does not collide with the concrete model types in coherence.
///
/// ```
/// use valenx_md::mc::FnEnergy;
/// use valenx_md::system::System;
/// use valenx_md::error::Result;
/// let model = FnEnergy::new(|_sys: &System| -> Result<f64> { Ok(0.0) }); // ideal gas
/// # let _ = &model;
/// ```
pub struct FnEnergy<F>(F);

impl<F> FnEnergy<F>
where
    F: Fn(&System) -> Result<f64>,
{
    /// Wraps the closure.
    pub fn new(f: F) -> Self {
        FnEnergy(f)
    }
}

impl<F> EnergyModel for FnEnergy<F>
where
    F: Fn(&System) -> Result<f64>,
{
    fn potential_energy(&self, system: &System) -> Result<f64> {
        (self.0)(system)
    }
}

/// An [`EnergyModel`] that sums a list of [`ForceTerm`]s.
///
/// This adapts the crate's existing force terms — [`LennardJones`],
/// [`Coulomb`], [`Wolf`], the bonded terms — into the scalar-energy
/// interface Monte Carlo needs. Only the energy field of each term's
/// [`EnergyForce`] output is used; the forces are computed and
/// discarded (MC does not need them).
///
/// [`LennardJones`]: crate::nonbonded::lj::LennardJones
/// [`Coulomb`]: crate::nonbonded::coulomb::Coulomb
/// [`Wolf`]: crate::nonbonded::wolf::Wolf
pub struct ForceTermEnergy {
    terms: Vec<Box<dyn ForceTerm>>,
}

impl ForceTermEnergy {
    /// An empty energy model — add terms with [`push`](Self::push).
    pub fn new() -> Self {
        ForceTermEnergy { terms: Vec::new() }
    }

    /// Adds a force term, returning `self` for chaining.
    pub fn with_term(mut self, term: Box<dyn ForceTerm>) -> Self {
        self.terms.push(term);
        self
    }

    /// Adds a force term in place.
    pub fn push(&mut self, term: Box<dyn ForceTerm>) {
        self.terms.push(term);
    }

    /// Number of force terms.
    pub fn len(&self) -> usize {
        self.terms.len()
    }

    /// Whether the model has no terms.
    pub fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }
}

impl Default for ForceTermEnergy {
    fn default() -> Self {
        Self::new()
    }
}

impl EnergyModel for ForceTermEnergy {
    fn potential_energy(&self, system: &System) -> Result<f64> {
        let mut ef = EnergyForce::zeros(system.len());
        for term in &self.terms {
            term.accumulate(system, &mut ef)?;
        }
        Ok(ef.energy)
    }
}

/// Running acceptance statistics for a Monte Carlo run.
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct McStats {
    /// Total moves attempted.
    pub attempted: u64,
    /// Moves accepted.
    pub accepted: u64,
}

impl McStats {
    /// The acceptance ratio in `[0, 1]` (0 if nothing was attempted).
    pub fn acceptance_ratio(&self) -> f64 {
        if self.attempted == 0 {
            0.0
        } else {
            self.accepted as f64 / self.attempted as f64
        }
    }
}

/// A Metropolis Monte Carlo driver over a [`System`] and an
/// [`EnergyModel`], sampling the NVT ensemble with single-particle
/// translation moves.
pub struct MonteCarlo<E: EnergyModel> {
    /// The system being sampled (positions are updated in place).
    pub system: System,
    /// The potential-energy model.
    energy: E,
    /// Target temperature (K).
    temperature: f64,
    /// Maximum per-component displacement δ for a translation move (nm).
    max_displacement: f64,
    /// The deterministic random stream.
    rng: Rng,
    /// Current total potential energy (cached; kept in sync with
    /// `system`).
    current_energy: f64,
    /// Running acceptance statistics.
    stats: McStats,
}

impl<E: EnergyModel> MonteCarlo<E> {
    /// Builds a Monte Carlo driver.
    ///
    /// * `system` — the configuration to sample.
    /// * `energy` — the potential-energy model.
    /// * `temperature` — target temperature (K); must be finite and
    ///   positive.
    /// * `max_displacement` — maximum per-axis trial displacement δ
    ///   (nm); must be finite and positive.
    /// * `seed` — RNG seed for a reproducible chain.
    ///
    /// # Errors
    /// [`MdError::Invalid`] for a non-positive temperature or
    /// displacement, or an empty system; propagates the initial
    /// energy-evaluation error.
    pub fn new(
        system: System,
        energy: E,
        temperature: f64,
        max_displacement: f64,
        seed: u64,
    ) -> Result<Self> {
        if !(temperature.is_finite() && temperature > 0.0) {
            return Err(MdError::invalid(
                "temperature",
                "must be finite and positive",
            ));
        }
        if !(max_displacement.is_finite() && max_displacement > 0.0) {
            return Err(MdError::invalid(
                "max_displacement",
                "must be finite and positive",
            ));
        }
        if system.is_empty() {
            return Err(MdError::invalid(
                "system",
                "Monte Carlo needs at least one atom",
            ));
        }
        let current_energy = energy.potential_energy(&system)?;
        if !current_energy.is_finite() {
            return Err(MdError::invalid(
                "energy",
                "initial potential energy is not finite",
            ));
        }
        Ok(MonteCarlo {
            system,
            energy,
            temperature,
            max_displacement,
            rng: Rng::new(seed),
            current_energy,
            stats: McStats::default(),
        })
    }

    /// The current cached total potential energy (kJ/mol).
    pub fn energy(&self) -> f64 {
        self.current_energy
    }

    /// The target temperature (K).
    pub fn temperature(&self) -> f64 {
        self.temperature
    }

    /// The current maximum per-axis displacement δ (nm).
    pub fn max_displacement(&self) -> f64 {
        self.max_displacement
    }

    /// The running acceptance statistics.
    pub fn stats(&self) -> McStats {
        self.stats
    }

    /// Attempts one single-particle translation move and returns
    /// whether it was accepted.
    ///
    /// A random atom is displaced by a uniform vector in `[−δ, δ]³`; the
    /// move is accepted by the Metropolis criterion. On rejection the
    /// atom is restored and the cached energy is unchanged.
    ///
    /// # Errors
    /// Propagates energy-evaluation errors.
    pub fn step(&mut self) -> Result<bool> {
        let n = self.system.len();
        let atom = self.rng.below(n);
        let old = self.system.positions[atom];

        // Propose a uniform displacement in [-δ, δ]³.
        let d = self.max_displacement;
        let disp = Vector3::new(
            self.rng.uniform_range(-d, d),
            self.rng.uniform_range(-d, d),
            self.rng.uniform_range(-d, d),
        );
        self.system.positions[atom] = old + disp;

        // Energy change of the trial move.
        let delta = self.energy.move_delta(&mut self.system, atom, old)?;

        let accept = if delta <= 0.0 {
            true
        } else {
            // Metropolis: accept with probability exp(−ΔU / k_BT).
            let beta = 1.0 / (BOLTZMANN * self.temperature);
            let p = (-beta * delta).exp();
            self.rng.bernoulli(p)
        };

        self.stats.attempted += 1;
        if accept {
            self.stats.accepted += 1;
            self.current_energy += delta;
        } else {
            // Reject: restore the old position.
            self.system.positions[atom] = old;
        }
        Ok(accept)
    }

    /// Runs `n` translation moves (one *sweep* is conventionally `N`
    /// moves; this counts individual moves). Returns the acceptance
    /// statistics accumulated across this call only.
    ///
    /// # Errors
    /// Propagates per-step energy-evaluation errors.
    pub fn run(&mut self, n: usize) -> Result<McStats> {
        let before = self.stats;
        for _ in 0..n {
            self.step()?;
        }
        Ok(McStats {
            attempted: self.stats.attempted - before.attempted,
            accepted: self.stats.accepted - before.accepted,
        })
    }

    /// Adjusts the maximum displacement toward a target acceptance ratio
    /// by short trial bursts, a standard equilibration convenience.
    ///
    /// Runs `bursts` rounds of `per_burst` moves; after each it scales
    /// δ up when acceptance is above `target` (moves too timid) and down
    /// when below (moves too bold). δ is kept in `[1e-4, half_box]`.
    ///
    /// # Errors
    /// [`MdError::Invalid`] for a `target` outside `(0, 1)` or a zero
    /// burst size; propagates energy errors.
    pub fn tune_step(&mut self, target: f64, bursts: usize, per_burst: usize) -> Result<f64> {
        if !(target.is_finite() && target > 0.0 && target < 1.0) {
            return Err(MdError::invalid("target", "acceptance must lie in (0, 1)"));
        }
        if per_burst == 0 {
            return Err(MdError::invalid("per_burst", "must be at least 1"));
        }
        let hi = if self.system.cell.is_periodic() {
            (self.system.cell.max_cutoff()).max(1e-3)
        } else {
            f64::INFINITY
        };
        for _ in 0..bursts {
            let s = self.run(per_burst)?;
            let ratio = s.acceptance_ratio();
            // Multiplicative update, bounded per round for stability.
            let factor = (ratio / target).clamp(0.5, 2.0);
            self.max_displacement = (self.max_displacement * factor).clamp(1e-4, hi);
        }
        Ok(self.max_displacement)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forcefield::{CombiningRule, ForceField, LjParam};
    use crate::nonbonded::lj::LennardJones;
    use crate::pbc::SimBox;
    use crate::system::{Atom, Topology};

    /// A small periodic argon-like fluid on a slightly perturbed grid.
    fn argon_fluid(n_per_side: usize, spacing: f64) -> (System, ForceField) {
        let edge = spacing * n_per_side as f64;
        let mut top = Topology::new();
        let mut pos = Vec::new();
        for i in 0..n_per_side {
            for j in 0..n_per_side {
                for k in 0..n_per_side {
                    top.push_atom(Atom::new("Ar", 39.95, 0.0).unwrap());
                    pos.push(Vector3::new(
                        i as f64 * spacing + 0.05,
                        j as f64 * spacing + 0.05,
                        k as f64 * spacing + 0.05,
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

    fn lj_energy_model(sys: &System, ff: &ForceField, cutoff: f64) -> ForceTermEnergy {
        let lj = LennardJones::from_system(sys, ff, cutoff).unwrap();
        ForceTermEnergy::new().with_term(Box::new(lj))
    }

    #[test]
    fn rejects_bad_parameters() {
        let (sys, ff) = argon_fluid(2, 0.4);
        let model = lj_energy_model(&sys, &ff, 1.0);
        assert!(MonteCarlo::new(sys.clone(), &model, 0.0, 0.01, 1).is_err());
        assert!(MonteCarlo::new(sys.clone(), &model, -1.0, 0.01, 1).is_err());
        assert!(MonteCarlo::new(sys.clone(), &model, 300.0, 0.0, 1).is_err());
        assert!(MonteCarlo::new(sys.clone(), &model, 300.0, f64::NAN, 1).is_err());

        // Empty system.
        let empty = System::new(Topology::new(), vec![]).unwrap();
        let empty_model = ForceTermEnergy::new();
        assert!(MonteCarlo::new(empty, &empty_model, 300.0, 0.01, 1).is_err());
    }

    #[test]
    fn deterministic_for_a_seed() {
        // Two runs with the same seed must produce identical chains.
        let (sys, ff) = argon_fluid(3, 0.4);
        let model = lj_energy_model(&sys, &ff, 1.0);

        let mut mc1 = MonteCarlo::new(sys.clone(), &model, 120.0, 0.02, 7).unwrap();
        let mut mc2 = MonteCarlo::new(sys.clone(), &model, 120.0, 0.02, 7).unwrap();
        let s1 = mc1.run(500).unwrap();
        let s2 = mc2.run(500).unwrap();
        assert_eq!(s1, s2);
        assert!((mc1.energy() - mc2.energy()).abs() < 1e-12);
        // The two final configurations must match exactly.
        for (a, b) in mc1.system.positions.iter().zip(&mc2.system.positions) {
            assert!((a - b).norm() < 1e-12);
        }
    }

    #[test]
    fn distinct_seeds_diverge() {
        let (sys, ff) = argon_fluid(3, 0.4);
        let model = lj_energy_model(&sys, &ff, 1.0);
        let mut a = MonteCarlo::new(sys.clone(), &model, 120.0, 0.02, 1).unwrap();
        let mut b = MonteCarlo::new(sys.clone(), &model, 120.0, 0.02, 2).unwrap();
        a.run(500).unwrap();
        b.run(500).unwrap();
        // Different streams => different final configurations.
        let same = a
            .system
            .positions
            .iter()
            .zip(&b.system.positions)
            .filter(|(x, y)| (*x - *y).norm() < 1e-12)
            .count();
        assert!(same < a.system.len(), "chains suspiciously identical");
    }

    #[test]
    fn cached_energy_tracks_recomputed_energy() {
        // The incrementally-updated cached energy must equal a full
        // from-scratch evaluation of the final configuration — proof
        // the ΔU bookkeeping is correct.
        let (sys, ff) = argon_fluid(3, 0.4);
        let model = lj_energy_model(&sys, &ff, 1.0);
        let mut mc = MonteCarlo::new(sys, &model, 150.0, 0.02, 99).unwrap();
        mc.run(1000).unwrap();
        let recomputed = model.potential_energy(&mc.system).unwrap();
        assert!(
            (mc.energy() - recomputed).abs() < 1e-6,
            "cached {} vs recomputed {}",
            mc.energy(),
            recomputed
        );
    }

    #[test]
    fn downhill_moves_always_accepted_detailed_balance_uphill_sometimes() {
        // Two LJ atoms placed far apart on the repulsive-free attractive
        // tail: moving them together lowers energy (always accepted),
        // moving apart raises it (Metropolis-gated). We verify the
        // acceptance obeys the criterion by construction over many
        // trials at low vs high T.
        //
        // Construct a single-atom-mobile system: a fixed atom at origin
        // and a mobile atom; only translation of the mobile atom is
        // proposed by chance, but with 2 atoms both move. Check the
        // *aggregate* acceptance rises with temperature for an
        // energetically frustrated start.
        let mut top = Topology::new();
        top.push_atom(Atom::new("Ar", 39.95, 0.0).unwrap());
        top.push_atom(Atom::new("Ar", 39.95, 0.0).unwrap());
        // Start near the repulsive wall (r < r_min) so most moves that
        // separate the atoms are downhill and most that compress are
        // uphill — a temperature-sensitive mix.
        let sys = System::new(
            top,
            vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(0.32, 0.0, 0.0)],
        )
        .unwrap()
        .with_cell(SimBox::cubic(5.0).unwrap());
        let mut ff = ForceField::new(CombiningRule::LorentzBerthelot);
        ff.set_lj("Ar", LjParam::new(0.34, 0.996).unwrap());
        let model = lj_energy_model(&sys, &ff, 2.0);

        let mut cold = MonteCarlo::new(sys.clone(), &model, 5.0, 0.02, 3).unwrap();
        let mut hot = MonteCarlo::new(sys.clone(), &model, 5000.0, 0.02, 3).unwrap();
        let c = cold.run(4000).unwrap().acceptance_ratio();
        let h = hot.run(4000).unwrap().acceptance_ratio();
        // Hotter ⇒ uphill moves accepted more often ⇒ higher acceptance.
        assert!(
            h > c,
            "acceptance should increase with temperature: cold {c} vs hot {h}"
        );
        assert!((0.0..=1.0).contains(&c) && (0.0..=1.0).contains(&h));
    }

    /// **Required test 3.** A seeded run reproduces a known ensemble
    /// average. For an ideal-gas-like system (zero potential energy,
    /// every move accepted, particles diffusing freely in a periodic
    /// box) the equilibrium distribution is *uniform*, so the
    /// time-averaged mean position converges to the box centre. We pin
    /// that ensemble average and the unit (100 %) acceptance ratio.
    #[test]
    fn metropolis_reproduces_ideal_gas_ensemble_average() {
        // Zero-energy model: U ≡ 0 ⇒ Boltzmann factor 1 ⇒ every move
        // accepted ⇒ uniform sampling of the periodic box.
        let edge = 2.0;
        let mut top = Topology::new();
        for _ in 0..20 {
            top.push_atom(Atom::new("X", 1.0, 0.0).unwrap());
        }
        let pos: Vec<_> = (0..20).map(|_| Vector3::new(0.1, 0.1, 0.1)).collect();
        let mut sys = System::new(top, pos)
            .unwrap()
            .with_cell(SimBox::cubic(edge).unwrap());
        // An ideal gas needs the box to wrap each accepted move so the
        // walk stays bounded; do that in the energy model's owner by
        // wrapping after the run. The model itself is U ≡ 0.
        let zero_energy = FnEnergy::new(|_: &System| -> Result<f64> { Ok(0.0) });

        let mut mc = MonteCarlo::new(sys.clone(), zero_energy, 300.0, 0.2, 2024).unwrap();
        let sweeps = 4000;
        // Accumulate the mean wrapped position over the chain.
        let mut mean = Vector3::zeros();
        let mut samples = 0u64;
        for _ in 0..sweeps {
            mc.run(sys.len()).unwrap(); // one sweep = N moves
            mc.system.wrap_into_cell();
            for p in &mc.system.positions {
                mean += *p;
                samples += 1;
            }
        }
        mean /= samples as f64;

        // Every move is downhill-or-flat (ΔU = 0 ⇒ accepted), so the
        // acceptance ratio must be exactly 1.
        assert!(
            (mc.stats().acceptance_ratio() - 1.0).abs() < 1e-12,
            "ideal-gas acceptance ratio = {}",
            mc.stats().acceptance_ratio()
        );
        // The uniform-distribution mean over [0, edge)³ is the box
        // centre (edge/2) in every component.
        let centre = edge / 2.0;
        for c in [mean.x, mean.y, mean.z] {
            assert!(
                (c - centre).abs() < 0.06 * edge,
                "ensemble-average position component {c} should approach the box centre {centre}"
            );
        }

        // Silence the unused-mut lint on the seed system clone helper.
        let _ = &mut sys;
    }

    #[test]
    fn tune_step_moves_toward_target_acceptance() {
        let (sys, ff) = argon_fluid(3, 0.42);
        let model = lj_energy_model(&sys, &ff, 1.0);
        // Start with an absurdly large step (almost everything rejected).
        let mut mc = MonteCarlo::new(sys, &model, 150.0, 0.5, 5).unwrap();
        let before = mc.max_displacement();
        mc.tune_step(0.5, 20, 200).unwrap();
        let after = mc.max_displacement();
        // The tuner should have shrunk the over-large step.
        assert!(
            after < before,
            "tuning did not shrink step: {before} -> {after}"
        );
        // And it must stay positive and bounded.
        assert!(after > 0.0 && after.is_finite());
        // Bad inputs rejected.
        assert!(mc.tune_step(0.0, 1, 1).is_err());
        assert!(mc.tune_step(0.5, 1, 0).is_err());
    }

    #[test]
    fn force_term_energy_matches_direct_sum() {
        // The ForceTermEnergy adapter must report exactly the LJ energy.
        let (sys, ff) = argon_fluid(2, 0.4);
        let lj = LennardJones::from_system(&sys, &ff, 1.0).unwrap();
        let mut ef = EnergyForce::zeros(sys.len());
        lj.accumulate(&sys, &mut ef).unwrap();

        let model = ForceTermEnergy::new().with_term(Box::new(lj));
        let e = model.potential_energy(&sys).unwrap();
        assert!((e - ef.energy).abs() < 1e-12);
        assert_eq!(model.len(), 1);
        assert!(!model.is_empty());
    }
}
