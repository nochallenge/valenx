//! # valenx-md — molecular-dynamics engine
//!
//! Round 6 Block 6 of the Valenx roadmap. A native-Rust molecular-
//! dynamics engine — a from-scratch reimplementation of the *classical*
//! core shared by GROMACS, OpenMM, LAMMPS, NAMD, AmberTools, Tinker,
//! HOOMD-blue, MDAnalysis, MDTraj and oxDNA. Pure algorithms: no GPU
//! kernels, no external processes, no neural-network weights.
//!
//! Everything works in the GROMACS unit system (nm, ps, u, e, kJ/mol);
//! see [`units`].
//!
//! ## What it does
//!
//! - **Model & I/O** ([`system`], [`forcefield`], [`pbc`], [`io`]) —
//!   the [`System`] / [`Atom`] / [`Topology`] data model, a force-field
//!   parameter model with Lorentz-Berthelot and geometric combining
//!   rules, periodic [`SimBox`]es (orthorhombic + triclinic), PDB /
//!   XYZ / GRO readers and writers, and DCD-class + framed-text
//!   trajectory writers.
//! - **A real atom-typed force field**
//!   ([`forcefield::typing`], [`forcefield::oplsaa`],
//!   [`forcefield::parameterize`]) — atom-type perception (elements +
//!   bonded connectivity + perceived hybridization → force-field atom
//!   types), a faithful representative subset of the **OPLS-AA** force
//!   field (per-type Lennard-Jones σ/ε and partial charges, bonded
//!   constants keyed by atom-type tuples, the geometric combining
//!   rule), and a [`parameterize`] path that types a molecule, looks
//!   up the database and assigns real bonded + nonbonded parameters —
//!   the validated force field commercial MD (GROMACS / AMBER /
//!   OpenMM) is built on. The generic positional parameter path stays
//!   available for systems the force field cannot type.
//! - **Nonbonded** ([`nonbonded`]) — cutoff + shifted Lennard-Jones,
//!   direct Coulomb with reaction-field correction, an Ewald / PME v1,
//!   a cell list + Verlet neighbour list, and the minimum-image
//!   convention.
//! - **Bonded** ([`bonded`]) — harmonic bonds, harmonic angles,
//!   periodic and Ryckaert-Bellemans proper dihedrals, harmonic
//!   improper dihedrals — each returning energy *and* analytic forces.
//! - **Integrators** ([`integrate`]) — velocity-Verlet, leapfrog, and
//!   a Langevin / Brownian-dynamics integrator.
//! - **Ensembles** ([`ensemble`]) — Berendsen, Andersen,
//!   velocity-rescale (Bussi) and Nosé-Hoover-chain thermostats; the
//!   Berendsen and Parrinello-Rahman barostats.
//! - **Constraints & minimisation** ([`constrain`], [`minimize`]) —
//!   SHAKE and RATTLE constraint solvers, and steepest-descent,
//!   conjugate-gradient and L-BFGS energy minimisers.
//! - **Analysis** ([`analysis`]) — energy / temperature / virial-
//!   pressure reporters, the radial distribution function, mean-
//!   squared displacement + the Einstein diffusion coefficient, and
//!   Kabsch-superposed RMSD / RMSF.
//! - **Coarse-grained & implicit solvent** ([`coarse`], [`implicit`])
//!   — an oxDNA-class coarse-grained DNA bead model and a
//!   generalized-Born implicit-solvent model.
//! - **Driver** ([`sim`]) — the [`Simulation`] type wires a system,
//!   force field, integrator and thermostat into a runnable loop.
//!
//! ## Errors
//!
//! Every fallible function returns
//! [`Result<_, MdError>`](error::MdError). The error type carries
//! stable [`code`](error::MdError::code) and
//! [`category`](error::MdError::category) accessors for telemetry.
//!
//! ## v1 scope
//!
//! This is a real working v1, not production parity with the decades-
//! old reference engines. Each module documents its own
//! simplifications inline; the notable ones are:
//!
//! - **No GPU kernels and no domain decomposition** — the force loop
//!   is a single-threaded CPU evaluation. Correct, not GROMACS-fast.
//! - **PME is a textbook Ewald split, not a true mesh code.** The
//!   reciprocal sum is evaluated by a direct structure-factor sum over
//!   k-vectors (an `O(N·K)` Ewald, not the `O(N log N)` FFT-based
//!   smooth-PME). It is exact Ewald within the k-space cutoff; it is
//!   simply slower than a real mesh. See [`nonbonded::pme`].
//! - **The Nosé-Hoover thermostat** is implemented as a single
//!   thermostat variable by default with an optional fixed-length
//!   chain; the chain uses a first-order operator split, not the full
//!   higher-order Suzuki-Yoshida factorisation. See
//!   [`ensemble::nose_hoover`].
//! - **The oxDNA-class model** is a *simplified* coarse-grained DNA
//!   force field — backbone (FENE), excluded volume, and a
//!   distance/orientation-gated hydrogen-bond + stacking term. It
//!   reproduces duplex formation qualitatively; it is not the
//!   published oxDNA2 parameterisation. See [`coarse`].
//! - **The generalized-Born implicit solvent** uses the
//!   Hawkins-Cramer-Truhlar pairwise-descreening Born radii and the
//!   Still GB energy kernel — a real GB/HCT v1, without the
//!   surface-area nonpolar term or the GBn2 corrections. See
//!   [`implicit`].
//! - **The OPLS-AA force field is a faithful representative subset,
//!   not the full release.** [`forcefield::oplsaa`] encodes genuine
//!   published OPLS-AA parameters for the common organic chemistry —
//!   C/H/N/O/S and the halogens in their usual hybridizations, and the
//!   bonded terms connecting them — and the proper-torsion table
//!   returns the dominant Fourier term per torsion class rather than
//!   the full three-term series. A type or functional group outside
//!   that coverage returns an honest error so the caller can fall back
//!   to the generic parameter path. Validated biomolecular coverage
//!   (the full protein / nucleic-acid residue libraries) is the
//!   documented next step. See [`forcefield::parameterize`].
//! - Free-energy methods (FEP / TI / umbrella sampling), Drude /
//!   polarisable force fields, QM/MM and replica exchange are out of
//!   scope and surface as [`MdError::NotYetImplemented`] where a hook
//!   exists.

#![forbid(unsafe_code)]

pub mod analysis;
pub mod bonded;
pub mod coarse;
pub mod constrain;
pub mod ensemble;
pub mod error;
pub mod forcefield;
pub mod implicit;
pub mod integrate;
pub mod io;
pub mod minimize;
pub mod nonbonded;
pub mod pbc;
pub mod rng;
pub mod sim;
pub mod system;
pub mod units;

// --- Convenience re-exports of the most-used types --------------------

pub use error::{ErrorCategory, MdError, Result};
pub use forcefield::{
    AngleParam, BondParam, CombiningRule, DihedralKind, DihedralParam, ForceField, ImproperParam,
    LjParam,
};
pub use forcefield::parameterize::{parameterize, Parameterized, ParameterizeOptions};
pub use forcefield::typing::{AtomType, Hybridization};
pub use pbc::SimBox;
pub use rng::Rng;
pub use system::{Angle, Atom, Bond, Dihedral, Improper, System, Topology};

pub use bonded::{EnergyForce, ForceTerm};
pub use integrate::Integrator;
pub use nonbonded::neighbor::NeighborList;
pub use sim::{Simulation, SimulationReport};

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    /// End-to-end: build a tiny system, attach a force field, run a
    /// few velocity-Verlet steps, and confirm nothing diverges.
    #[test]
    fn build_and_run_end_to_end() {
        // Two argon-like atoms.
        let mut top = Topology::new();
        top.push_atom(Atom::new("Ar", 39.95, 0.0).unwrap());
        top.push_atom(Atom::new("Ar", 39.95, 0.0).unwrap());
        let positions = vec![Vector3::zeros(), Vector3::new(0.4, 0.0, 0.0)];
        let system = System::new(top, positions)
            .unwrap()
            .with_cell(SimBox::cubic(5.0).unwrap());

        let mut ff = ForceField::new(CombiningRule::LorentzBerthelot);
        ff.set_lj("Ar", LjParam::new(0.34, 0.996).unwrap());

        let mut sim = Simulation::new(system, ff).unwrap();
        let report = sim.run(20).unwrap();
        assert_eq!(report.steps, 20);
        assert!(report.final_total_energy.is_finite());
    }

    #[test]
    fn re_exports_are_wired() {
        let e = MdError::not_yet("x");
        assert_eq!(e.category_enum(), ErrorCategory::Capability);
        let _ = SimBox::cubic(2.0).unwrap();
        let _ = Rng::new(1);
    }
}
