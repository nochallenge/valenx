//! Top-level calculation drivers and the [`QchemReport`].
//!
//! This module is the one most callers touch: it bundles geometry +
//! basis setup, the SCF, and the property analysis behind four
//! one-call entry points —
//!
//! - [`run_rhf`] — restricted Hartree-Fock for a closed-shell molecule;
//! - [`run_uhf`] — unrestricted Hartree-Fock for an open-shell one;
//! - [`run_mp2`] — RHF followed by the MP2 correlation correction;
//! - [`run_dft`] — Kohn-Sham density-functional theory —
//!
//! each returning a [`QchemReport`] with the energy, atomic charges,
//! dipole, orbital energies, HOMO-LUMO gap and per-stage timings.
//!
//! ## The remaining honest stub — geometry optimisation
//!
//! One genuine subsystem is *not* implemented and says so plainly.
//! [`GeometryOptRequest::run`] is a typed request struct whose `run`
//! method returns
//! [`QchemError::NotYetImplemented`](crate::error::QchemError) with a
//! clear message: **geometry optimisation** needs analytic energy
//! gradients (the integral derivatives) and a step-control optimiser.
//! It is stubbed honestly rather than faked so the request API is
//! stable for when that subsystem lands.
//!
//! Density-functional theory, formerly a stub, is now real — see
//! [`DftRequest`] and the [`dft`](crate::dft) module.

use crate::basis::BasisSet;
use crate::dft::ks::run_ks;
use crate::dft::{Functional, GridQuality, KsResult};
use crate::error::{QchemError, Result};
use crate::geometry::MolecularGeometry;
use crate::integrals::IntegralSet;
use crate::post::mp2::{mp2_energy, Mp2Result};
use crate::properties::dipole::{dipole_moment, DipoleMoment};
use crate::properties::orbitals::{restricted_summary, unrestricted_spin_summary, OrbitalSummary};
use crate::properties::population::{mulliken, PopulationAnalysis};
use crate::scf::rhf::{run_rhf_scf, RhfResult, ScfSettings};
use crate::scf::uhf::{run_uhf_scf, UhfResult};
use std::time::{Duration, Instant};

/// Which method produced a report.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Method {
    /// Restricted Hartree-Fock (closed shell).
    Rhf,
    /// Unrestricted Hartree-Fock (open shell).
    Uhf,
    /// RHF reference plus an MP2 correlation correction.
    Mp2,
    /// Restricted Kohn-Sham density-functional theory.
    Dft,
}

impl Method {
    /// A short human-readable label.
    pub fn label(self) -> &'static str {
        match self {
            Method::Rhf => "RHF",
            Method::Uhf => "UHF",
            Method::Mp2 => "MP2",
            Method::Dft => "DFT",
        }
    }
}

/// Wall-clock timings for the stages of a calculation.
#[derive(Copy, Clone, Debug, Default)]
pub struct Timings {
    /// Time spent computing the molecular integrals.
    pub integrals: Duration,
    /// Time spent in the SCF iteration.
    pub scf: Duration,
    /// Time spent in post-SCF work (MP2, property analysis).
    pub post: Duration,
}

impl Timings {
    /// The total wall-clock time of the calculation.
    pub fn total(&self) -> Duration {
        self.integrals + self.scf + self.post
    }
}

/// A complete result report for a quantum-chemistry calculation.
#[derive(Clone, Debug)]
pub struct QchemReport {
    /// Which method produced this report.
    pub method: Method,
    /// The basis-set name used.
    pub basis_name: &'static str,
    /// The final total energy in Hartree — for MP2 this includes the
    /// correlation correction.
    pub total_energy: f64,
    /// The Hartree-Fock total energy (equals `total_energy` for plain
    /// RHF / UHF; the reference energy for MP2).
    pub hartree_fock_energy: f64,
    /// The MP2 correlation energy, when a correlation method was run.
    pub correlation_energy: Option<f64>,
    /// Nuclear-repulsion energy (Hartree).
    pub nuclear_repulsion: f64,
    /// Mulliken atomic partial charges, one per atom (`e`).
    pub partial_charges: Vec<f64>,
    /// The molecular dipole moment.
    pub dipole: DipoleMoment,
    /// The molecular-orbital summary (the α set for a UHF report).
    pub orbitals: OrbitalSummary,
    /// The HOMO-LUMO gap in Hartree (`None` when undefined).
    pub homo_lumo_gap: Option<f64>,
    /// Number of SCF cycles taken.
    pub scf_iterations: usize,
    /// `⟨S²⟩` for a UHF calculation (`None` for RHF / MP2 / DFT).
    pub s_squared: Option<f64>,
    /// DFT-specific results — `None` unless the method is
    /// [`Method::Dft`].
    pub dft: Option<DftInfo>,
    /// Per-stage wall-clock timings.
    pub timings: Timings,
}

/// The density-functional-theory–specific results of a [`QchemReport`].
#[derive(Clone, Debug)]
pub struct DftInfo {
    /// The exchange-correlation functional used.
    pub functional: Functional,
    /// The exchange-correlation energy `E_xc` (Hartree).
    pub xc_energy: f64,
    /// The number of electrons recovered by integrating the density on
    /// the molecular grid — a grid-quality diagnostic that should
    /// equal the true electron count.
    pub grid_electron_count: f64,
}

impl QchemReport {
    /// A compact multi-line human-readable summary.
    pub fn summary(&self) -> String {
        let mut s = format!(
            "{} / {}\n  total energy   {:.8} Ha\n",
            self.method.label(),
            self.basis_name,
            self.total_energy
        );
        if let Some(corr) = self.correlation_energy {
            s.push_str(&format!(
                "  HF energy      {:.8} Ha\n",
                self.hartree_fock_energy
            ));
            s.push_str(&format!("  correlation    {corr:.8} Ha\n"));
        }
        if let Some(dft) = &self.dft {
            s.push_str(&format!("  functional     {}\n", dft.functional.label()));
            s.push_str(&format!("  XC energy      {:.8} Ha\n", dft.xc_energy));
        }
        s.push_str(&format!(
            "  nuclear rep.   {:.8} Ha\n",
            self.nuclear_repulsion
        ));
        if let Some(gap) = self.homo_lumo_gap {
            s.push_str(&format!("  HOMO-LUMO gap  {gap:.6} Ha\n"));
        }
        if let Some(s2) = self.s_squared {
            s.push_str(&format!("  <S^2>          {s2:.4}\n"));
        }
        s.push_str(&format!("  SCF cycles     {}\n", self.scf_iterations));
        s.push_str(&format!(
            "  wall time      {:.3} s\n",
            self.timings.total().as_secs_f64()
        ));
        s
    }
}

/// Validate a molecule and build its basis set, the shared front end of
/// every driver.
fn setup(geometry: &MolecularGeometry, basis_name: &str) -> Result<BasisSet> {
    geometry.validate()?;
    BasisSet::build(basis_name, geometry)
}

/// Run a restricted-Hartree-Fock calculation and return a full report.
///
/// The molecule must be a closed-shell singlet.
///
/// # Errors
///
/// Propagates [`QchemError`] from geometry validation, basis lookup or
/// SCF convergence; returns [`QchemError::InvalidInput`] when the
/// molecule is not closed-shell.
pub fn run_rhf(
    geometry: &MolecularGeometry,
    basis_name: &str,
    settings: ScfSettings,
) -> Result<QchemReport> {
    run_rhf_embedded(geometry, basis_name, settings, &[])
}

/// [`run_rhf`] in the field of an external electrostatic environment —
/// a set of point charges `(q, position_bohr)` added to the one-electron
/// Hamiltonian (**electrostatic QM/MM embedding**). The MM charges enter
/// the SCF and polarize the density. The returned `total_energy` includes
/// the electron–charge interaction; the classical nuclei–charge term is
/// the caller's to add. With an empty `external_charges` this is exactly
/// [`run_rhf`].
pub fn run_rhf_embedded(
    geometry: &MolecularGeometry,
    basis_name: &str,
    settings: ScfSettings,
    external_charges: &[(f64, [f64; 3])],
) -> Result<QchemReport> {
    if !geometry.is_closed_shell() {
        return Err(QchemError::invalid(
            "run_rhf requires a closed-shell singlet — use run_uhf for \
             open-shell systems",
        ));
    }
    let basis = setup(geometry, basis_name)?;
    let n_electrons = geometry.n_electrons()?;

    let t0 = Instant::now();
    let mut integrals = IntegralSet::compute(geometry, &basis);
    if !external_charges.is_empty() {
        integrals.nuclear += crate::integrals::external_charge_potential(&basis, external_charges);
    }
    let t_integrals = t0.elapsed();

    let t1 = Instant::now();
    let rhf = run_rhf_scf(&integrals, n_electrons, settings)?;
    let t_scf = t1.elapsed();

    let t2 = Instant::now();
    let report = build_rhf_report(geometry, &basis, &integrals, &rhf, None);
    let t_post = t2.elapsed();

    Ok(QchemReport {
        timings: Timings {
            integrals: t_integrals,
            scf: t_scf,
            post: t_post,
        },
        ..report
    })
}

/// Run an unrestricted-Hartree-Fock calculation and return a full
/// report. Works for both open- and closed-shell molecules.
///
/// # Errors
///
/// Propagates [`QchemError`] from geometry validation, basis lookup or
/// SCF convergence.
pub fn run_uhf(
    geometry: &MolecularGeometry,
    basis_name: &str,
    settings: ScfSettings,
) -> Result<QchemReport> {
    let basis = setup(geometry, basis_name)?;
    let (n_alpha, n_beta) = geometry.alpha_beta_counts()?;

    let t0 = Instant::now();
    let integrals = IntegralSet::compute(geometry, &basis);
    let t_integrals = t0.elapsed();

    let t1 = Instant::now();
    let uhf = run_uhf_scf(&integrals, n_alpha, n_beta, settings)?;
    let t_scf = t1.elapsed();

    let t2 = Instant::now();
    let report = build_uhf_report(geometry, &basis, &integrals, &uhf);
    let t_post = t2.elapsed();

    Ok(QchemReport {
        timings: Timings {
            integrals: t_integrals,
            scf: t_scf,
            post: t_post,
        },
        ..report
    })
}

/// Run an RHF calculation followed by the MP2 correlation correction.
///
/// The molecule must be a closed-shell singlet (MP2 here uses an RHF
/// reference).
///
/// # Errors
///
/// Propagates [`QchemError`] from setup, SCF or the MP2 step.
pub fn run_mp2(
    geometry: &MolecularGeometry,
    basis_name: &str,
    settings: ScfSettings,
) -> Result<QchemReport> {
    if !geometry.is_closed_shell() {
        return Err(QchemError::invalid(
            "run_mp2 requires a closed-shell singlet (RHF-reference MP2)",
        ));
    }
    let basis = setup(geometry, basis_name)?;
    let n_electrons = geometry.n_electrons()?;

    let t0 = Instant::now();
    let integrals = IntegralSet::compute(geometry, &basis);
    let t_integrals = t0.elapsed();

    let t1 = Instant::now();
    let rhf = run_rhf_scf(&integrals, n_electrons, settings)?;
    let t_scf = t1.elapsed();

    let t2 = Instant::now();
    let mp2: Mp2Result = mp2_energy(&rhf, &integrals.eri)?;
    let report = build_rhf_report(geometry, &basis, &integrals, &rhf, Some(mp2));
    let t_post = t2.elapsed();

    Ok(QchemReport {
        timings: Timings {
            integrals: t_integrals,
            scf: t_scf,
            post: t_post,
        },
        ..report
    })
}

/// Run a restricted Kohn-Sham density-functional-theory calculation and
/// return a full report.
///
/// The molecule must be a closed-shell singlet (restricted Kohn-Sham).
/// `functional` selects the exchange-correlation functional —
/// [`Functional::Lda`], [`Functional::Pbe`] or [`Functional::B3lyp`];
/// `quality` selects the molecular integration grid.
///
/// # Errors
///
/// Propagates [`QchemError`] from geometry validation, basis lookup or
/// the Kohn-Sham SCF; returns [`QchemError::InvalidInput`] when the
/// molecule is not closed-shell.
pub fn run_dft(
    geometry: &MolecularGeometry,
    basis_name: &str,
    functional: Functional,
    quality: GridQuality,
    settings: ScfSettings,
) -> Result<QchemReport> {
    if !geometry.is_closed_shell() {
        return Err(QchemError::invalid(
            "run_dft requires a closed-shell singlet — restricted \
             Kohn-Sham only",
        ));
    }
    let basis = setup(geometry, basis_name)?;
    let n_electrons = geometry.n_electrons()?;

    let t0 = Instant::now();
    let integrals = IntegralSet::compute(geometry, &basis);
    let t_integrals = t0.elapsed();

    let t1 = Instant::now();
    let ks = run_ks(
        &integrals,
        &basis,
        geometry,
        n_electrons,
        functional,
        quality,
        settings,
    )?;
    let t_scf = t1.elapsed();

    let t2 = Instant::now();
    let report = build_dft_report(geometry, &basis, &integrals, &ks);
    let t_post = t2.elapsed();

    Ok(QchemReport {
        timings: Timings {
            integrals: t_integrals,
            scf: t_scf,
            post: t_post,
        },
        ..report
    })
}

/// Assemble a report from an RHF (optionally MP2-corrected) result.
fn build_rhf_report(
    geometry: &MolecularGeometry,
    basis: &BasisSet,
    integrals: &IntegralSet,
    rhf: &RhfResult,
    mp2: Option<Mp2Result>,
) -> QchemReport {
    let pop: PopulationAnalysis = mulliken(geometry, basis, &rhf.density, &integrals.overlap);
    let dipole = dipole_moment(geometry, integrals, &rhf.density);
    let orbitals = restricted_summary(&rhf.orbital_energies, rhf.n_occupied);
    let homo_lumo_gap = orbitals.homo_lumo_gap();

    let (method, total_energy, correlation_energy) = match mp2 {
        Some(m) => (Method::Mp2, m.total_energy(), Some(m.correlation_energy)),
        None => (Method::Rhf, rhf.total_energy, None),
    };

    QchemReport {
        method,
        basis_name: basis.name,
        total_energy,
        hartree_fock_energy: rhf.total_energy,
        correlation_energy,
        nuclear_repulsion: rhf.nuclear_repulsion,
        partial_charges: pop.partial_charges,
        dipole,
        orbitals,
        homo_lumo_gap,
        scf_iterations: rhf.iterations.len(),
        s_squared: None,
        dft: None,
        timings: Timings::default(),
    }
}

/// Assemble a report from a UHF result.
fn build_uhf_report(
    geometry: &MolecularGeometry,
    basis: &BasisSet,
    integrals: &IntegralSet,
    uhf: &UhfResult,
) -> QchemReport {
    let total_density = &uhf.alpha_density + &uhf.beta_density;
    let pop = mulliken(geometry, basis, &total_density, &integrals.overlap);
    let dipole = dipole_moment(geometry, integrals, &total_density);
    // Report the alpha orbital set.
    let orbitals = unrestricted_spin_summary(&uhf.alpha_orbital_energies, uhf.n_alpha);
    let homo_lumo_gap = orbitals.homo_lumo_gap();

    QchemReport {
        method: Method::Uhf,
        basis_name: basis.name,
        total_energy: uhf.total_energy,
        hartree_fock_energy: uhf.total_energy,
        correlation_energy: None,
        nuclear_repulsion: uhf.nuclear_repulsion,
        partial_charges: pop.partial_charges,
        dipole,
        orbitals,
        homo_lumo_gap,
        scf_iterations: uhf.iterations.len(),
        s_squared: Some(uhf.s_squared),
        dft: None,
        timings: Timings::default(),
    }
}

/// Assemble a report from a Kohn-Sham DFT result.
fn build_dft_report(
    geometry: &MolecularGeometry,
    basis: &BasisSet,
    integrals: &IntegralSet,
    ks: &KsResult,
) -> QchemReport {
    let pop: PopulationAnalysis = mulliken(geometry, basis, &ks.density, &integrals.overlap);
    let dipole = dipole_moment(geometry, integrals, &ks.density);
    let orbitals = restricted_summary(&ks.orbital_energies, ks.n_occupied);
    let homo_lumo_gap = orbitals.homo_lumo_gap();

    QchemReport {
        method: Method::Dft,
        basis_name: basis.name,
        total_energy: ks.total_energy,
        // For DFT, the "Hartree-Fock energy" slot carries the total
        // KS-DFT energy (there is no separate HF reference).
        hartree_fock_energy: ks.total_energy,
        correlation_energy: None,
        nuclear_repulsion: ks.nuclear_repulsion,
        partial_charges: pop.partial_charges,
        dipole,
        orbitals,
        homo_lumo_gap,
        scf_iterations: ks.iterations.len(),
        s_squared: None,
        dft: Some(DftInfo {
            functional: ks.functional,
            xc_energy: ks.xc_energy,
            grid_electron_count: ks.grid_electron_count,
        }),
        timings: Timings::default(),
    }
}

// =====================================================================
// DftRequest — the real Kohn-Sham DFT request API
// =====================================================================

/// A density-functional-theory calculation request.
///
/// A typed request struct over the real Kohn-Sham DFT subsystem (the
/// [`dft`](crate::dft) module). [`DftRequest::run`] builds the
/// molecular integration grid, runs the Kohn-Sham SCF and returns a
/// [`QchemReport`]. The functional is named by string — `"lda"` /
/// `"svwn"`, `"pbe"`, `"b3lyp"` — and parsed by
/// [`Functional::from_name`].
#[derive(Clone, Debug)]
pub struct DftRequest {
    /// The molecule to treat.
    pub geometry: MolecularGeometry,
    /// The basis-set name.
    pub basis_name: String,
    /// The requested exchange-correlation functional (`"lda"`,
    /// `"pbe"`, `"b3lyp"`).
    pub functional: String,
    /// The molecular integration-grid quality.
    pub grid_quality: GridQuality,
}

impl DftRequest {
    /// Construct a DFT request with the default (`Medium`) integration
    /// grid.
    pub fn new(
        geometry: MolecularGeometry,
        basis_name: impl Into<String>,
        functional: impl Into<String>,
    ) -> Self {
        DftRequest {
            geometry,
            basis_name: basis_name.into(),
            functional: functional.into(),
            grid_quality: GridQuality::default(),
        }
    }

    /// Set the molecular integration-grid quality.
    pub fn with_grid_quality(mut self, quality: GridQuality) -> Self {
        self.grid_quality = quality;
        self
    }

    /// Run the Kohn-Sham DFT calculation.
    ///
    /// # Errors
    ///
    /// Returns [`QchemError::InvalidInput`] when the functional name is
    /// unknown or the molecule is not closed-shell, and propagates
    /// every [`QchemError`] from geometry validation, basis lookup and
    /// the Kohn-Sham SCF.
    pub fn run(&self) -> Result<QchemReport> {
        let functional = Functional::from_name(&self.functional).ok_or_else(|| {
            QchemError::invalid(format!(
                "unknown exchange-correlation functional `{}` — \
                     supported: lda / svwn, pbe, b3lyp",
                self.functional
            ))
        })?;
        run_dft(
            &self.geometry,
            &self.basis_name,
            functional,
            self.grid_quality,
            ScfSettings::default(),
        )
    }
}

/// A geometry-optimisation request.
///
/// **This is an honest stub.** Optimising a molecular geometry needs
/// analytic energy gradients — the derivatives of the integrals with
/// respect to the nuclear coordinates — plus a step-controlled
/// optimiser. [`GeometryOptRequest::run`] returns
/// [`QchemError::NotYetImplemented`](crate::error::QchemError) rather
/// than a fabricated geometry.
#[derive(Clone, Debug)]
pub struct GeometryOptRequest {
    /// The starting geometry.
    pub geometry: MolecularGeometry,
    /// The basis-set name.
    pub basis_name: String,
    /// The energy method to optimise on (`"rhf"`, `"uhf"`).
    pub method: String,
    /// Maximum optimisation steps the caller would allow.
    pub max_steps: usize,
}

impl GeometryOptRequest {
    /// Construct a geometry-optimisation request.
    pub fn new(
        geometry: MolecularGeometry,
        basis_name: impl Into<String>,
        method: impl Into<String>,
    ) -> Self {
        GeometryOptRequest {
            geometry,
            basis_name: basis_name.into(),
            method: method.into(),
            max_steps: 50,
        }
    }

    /// Run the geometry optimisation.
    ///
    /// # Errors
    ///
    /// Always returns
    /// [`QchemError::NotYetImplemented`](crate::error::QchemError) —
    /// geometry optimisation needs analytic energy gradients and a
    /// step-control optimiser, a separate subsystem this v1 does not
    /// ship.
    pub fn run(&self) -> Result<MolecularGeometry> {
        Err(QchemError::not_yet(
            "geometry_optimization",
            "geometry optimisation needs analytic energy gradients \
             (integral derivatives) and a step-controlled optimiser — a \
             separate subsystem from the single-point Hartree-Fock core",
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Atom;

    fn h2() -> MolecularGeometry {
        MolecularGeometry::new(vec![
            Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.0]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.7414]).unwrap(),
        ])
    }

    fn water() -> MolecularGeometry {
        MolecularGeometry::new(vec![
            Atom::from_symbol_angstrom("O", [0.0, 0.0, 0.1173]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, 0.7572, -0.4692]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, -0.7572, -0.4692]).unwrap(),
        ])
    }

    #[test]
    fn run_rhf_produces_a_full_report() {
        let report = run_rhf(&h2(), "sto-3g", ScfSettings::default()).unwrap();
        assert_eq!(report.method, Method::Rhf);
        assert!((report.total_energy - (-1.1167)).abs() < 2.0e-3);
        assert_eq!(report.partial_charges.len(), 2);
        assert!(report.homo_lumo_gap.unwrap() > 0.0);
        assert!(report.s_squared.is_none());
        // The summary string mentions the method.
        assert!(report.summary().contains("RHF"));
    }

    #[test]
    fn external_point_charge_polarizes_and_shifts_energy() {
        let g = h2();
        let e0 = run_rhf(&g, "sto-3g", ScfSettings::default())
            .unwrap()
            .total_energy;
        // A +1 point charge ~4 bohr out attracts the electrons → the
        // electron–charge term lowers the energy.
        let near = run_rhf_embedded(
            &g,
            "sto-3g",
            ScfSettings::default(),
            &[(1.0, [0.0, 0.0, 4.0])],
        )
        .unwrap()
        .total_energy;
        assert!(
            near < e0,
            "a nearby + charge should lower the energy: {near} vs {e0}"
        );
        // The effect decays with distance: a charge farther out shifts
        // the energy less than the near one.
        let shift_near = (near - e0).abs();
        let farther = run_rhf_embedded(
            &g,
            "sto-3g",
            ScfSettings::default(),
            &[(1.0, [0.0, 0.0, 20.0])],
        )
        .unwrap()
        .total_energy;
        let shift_far = (farther - e0).abs();
        assert!(
            shift_far < shift_near,
            "the charge's effect should decay with distance: {shift_far} (d=20) vs {shift_near} (d=4)"
        );
        // Empty external charges is exactly plain run_rhf.
        let none = run_rhf_embedded(&g, "sto-3g", ScfSettings::default(), &[])
            .unwrap()
            .total_energy;
        assert!((none - e0).abs() < 1e-12);
    }

    #[test]
    fn run_rhf_rejects_open_shell() {
        let radical = MolecularGeometry::with_charge_multiplicity(
            vec![Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.0]).unwrap()],
            0,
            2,
        );
        assert!(run_rhf(&radical, "sto-3g", ScfSettings::default()).is_err());
    }

    #[test]
    fn run_uhf_on_closed_shell_matches_rhf() {
        let rhf = run_rhf(&h2(), "sto-3g", ScfSettings::default()).unwrap();
        let uhf = run_uhf(&h2(), "sto-3g", ScfSettings::default()).unwrap();
        assert!((uhf.total_energy - rhf.total_energy).abs() < 1.0e-6);
        assert_eq!(uhf.method, Method::Uhf);
        assert!(uhf.s_squared.is_some());
    }

    #[test]
    fn run_uhf_on_hydrogen_atom() {
        let h_atom = MolecularGeometry::with_charge_multiplicity(
            vec![Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.0]).unwrap()],
            0,
            2,
        );
        let report = run_uhf(&h_atom, "sto-3g", ScfSettings::default()).unwrap();
        assert!((report.total_energy - (-0.4666)).abs() < 1.0e-3);
    }

    #[test]
    fn run_mp2_lowers_energy_below_rhf() {
        let rhf = run_rhf(&water(), "sto-3g", ScfSettings::default()).unwrap();
        let mp2 = run_mp2(&water(), "sto-3g", ScfSettings::default()).unwrap();
        assert_eq!(mp2.method, Method::Mp2);
        assert!(mp2.total_energy < rhf.total_energy);
        assert!(mp2.correlation_energy.unwrap() < 0.0);
        assert!((mp2.hartree_fock_energy - rhf.total_energy).abs() < 1.0e-6);
        assert!(mp2.summary().contains("correlation"));
    }

    #[test]
    fn timings_are_recorded() {
        let report = run_rhf(&h2(), "sto-3g", ScfSettings::default()).unwrap();
        // Total time is the sum of the stages.
        let total = report.timings.total();
        assert_eq!(
            total,
            report.timings.integrals + report.timings.scf + report.timings.post
        );
    }

    #[test]
    fn dft_request_runs_real_kohn_sham() {
        // DftRequest is no longer a stub — it runs real Kohn-Sham DFT.
        let req = DftRequest::new(h2(), "sto-3g", "lda");
        let report = req.run().unwrap();
        assert_eq!(report.method, Method::Dft);
        assert!(report.total_energy.is_finite());
        // The DFT-specific info block is populated.
        let dft = report.dft.as_ref().unwrap();
        assert_eq!(dft.functional, crate::dft::Functional::Lda);
        assert!(dft.xc_energy < 0.0);
    }

    #[test]
    fn dft_request_rejects_unknown_functional() {
        let req = DftRequest::new(h2(), "sto-3g", "made-up-functional");
        let err = req.run().unwrap_err();
        assert_eq!(err.code(), "qchem.invalid_input");
        assert!(err.to_string().contains("functional"));
    }

    #[test]
    fn geometry_opt_request_is_an_honest_stub() {
        let req = GeometryOptRequest::new(h2(), "sto-3g", "rhf");
        let err = req.run().unwrap_err();
        assert_eq!(err.code(), "qchem.not_yet_implemented");
        assert!(err.to_string().contains("gradient"));
    }

    #[test]
    fn unknown_basis_propagates_from_driver() {
        assert!(run_rhf(&h2(), "cc-pvtz", ScfSettings::default()).is_err());
    }
}

/// Reference-value validation of restricted Hartree-Fock against
/// published STO-3G total energies. The references are the standard
/// textbook / literature numbers for the minimal STO-3G basis; the
/// tolerance (1 mHa) absorbs the small geometry differences between the
/// experimental geometries used here and the reference geometries.
#[cfg(test)]
mod validation {
    use super::*;
    use crate::geometry::Atom;

    /// H₂ / STO-3G near the experimental bond length (0.7414 Å).
    /// Szabo & Ostlund (Modern Quantum Chemistry) give the STO-3G
    /// minimal-basis RHF energy of H₂ as ≈ −1.1167 Ha.
    #[test]
    fn h2_sto3g_total_energy_matches_published() {
        let h2 = MolecularGeometry::new(vec![
            Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.0]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.7414]).unwrap(),
        ]);
        let r = run_rhf(&h2, "sto-3g", ScfSettings::default()).unwrap();
        assert!(
            (r.total_energy - (-1.1167)).abs() < 1.0e-3,
            "H2/STO-3G E = {:.6} Ha, expected ≈ -1.1167",
            r.total_energy
        );
    }

    /// HeH⁺ / STO-3G — the classic Szabo & Ostlund worked example.
    /// Their tabulated converged RHF energy for HeH⁺ in the STO-3G
    /// minimal basis is ≈ −2.8418 Ha.
    #[test]
    fn heh_cation_sto3g_total_energy_matches_szabo_ostlund() {
        let hehp = MolecularGeometry::with_charge_multiplicity(
            vec![
                Atom::from_symbol_angstrom("He", [0.0, 0.0, 0.0]).unwrap(),
                Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.7743]).unwrap(),
            ],
            1,
            1,
        );
        let r = run_rhf(&hehp, "sto-3g", ScfSettings::default()).unwrap();
        assert!(
            (r.total_energy - (-2.8418)).abs() < 2.0e-3,
            "HeH+/STO-3G E = {:.6} Ha, expected ≈ -2.8418",
            r.total_energy
        );
    }

    /// H₂O / STO-3G at the experimental geometry. The published STO-3G
    /// RHF energy of water sits at ≈ −74.96 Ha (−74.94 to −74.97 Ha
    /// depending on the precise geometry used).
    #[test]
    fn water_sto3g_total_energy_matches_published() {
        let h2o = MolecularGeometry::new(vec![
            Atom::from_symbol_angstrom("O", [0.0, 0.0, 0.1173]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, 0.7572, -0.4692]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, -0.7572, -0.4692]).unwrap(),
        ]);
        let r = run_rhf(&h2o, "sto-3g", ScfSettings::default()).unwrap();
        assert!(
            (r.total_energy - (-74.963)).abs() < 2.0e-2,
            "H2O/STO-3G E = {:.6} Ha, expected ≈ -74.96",
            r.total_energy
        );
        // The nuclear-repulsion energy of this geometry is a pure
        // geometric quantity ≈ 9.19 Ha.
        assert!(
            (r.nuclear_repulsion - 9.19).abs() < 0.05,
            "H2O nuclear repulsion = {:.4} Ha, expected ≈ 9.19",
            r.nuclear_repulsion
        );
    }

    /// MP2 must lower the energy below RHF and the STO-3G water MP2
    /// correlation energy is on the order of a few tens of mHa.
    #[test]
    fn water_sto3g_mp2_correlation_is_physical() {
        let h2o = MolecularGeometry::new(vec![
            Atom::from_symbol_angstrom("O", [0.0, 0.0, 0.1173]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, 0.7572, -0.4692]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, -0.7572, -0.4692]).unwrap(),
        ]);
        let r = run_mp2(&h2o, "sto-3g", ScfSettings::default()).unwrap();
        let corr = r.correlation_energy.unwrap();
        // MP2 correlation energy is always negative.
        assert!(corr < 0.0, "MP2 correlation = {corr:.6} Ha, must be < 0");
        // For STO-3G water it is a small but non-trivial correction.
        assert!(
            (-0.10..-0.005).contains(&corr),
            "MP2 correlation = {corr:.6} Ha, expected ~ -0.02..-0.05"
        );
        assert!(r.total_energy < r.hartree_fock_energy);
    }
}
