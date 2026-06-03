//! The force-field parameter model.
//!
//! **Roadmap feature 2.** A [`ForceField`] maps the *symbolic*
//! topology produced by [`crate::system`] to the *numeric* constants a
//! force calculation needs:
//!
//! - **Nonbonded** — per atom type a Lennard-Jones σ (nm) and ε
//!   (kJ/mol). Pair parameters for two different types follow a
//!   [`CombiningRule`]: Lorentz-Berthelot (σ arithmetic, ε geometric —
//!   AMBER / CHARMM / OPLS) or fully geometric (GROMACS `comb-rule 3`).
//! - **Bonded** — per [`Bond`](crate::system::Bond) /
//!   [`Angle`](crate::system::Angle) /
//!   [`Dihedral`](crate::system::Dihedral) /
//!   [`Improper`](crate::system::Improper) a [`BondParam`],
//!   [`AngleParam`], [`DihedralParam`] or [`ImproperParam`] — looked
//!   up positionally (the *n*-th bond uses the *n*-th [`BondParam`]).
//!
//! Charges already live on the [`Atom`](crate::system::Atom)s; the
//! force field carries the *Lennard-Jones* and *bonded* halves.
//!
//! Positional lookup keeps the generic path simple and unambiguous: a
//! topology builder pushes a bond and its parameters together.
//!
//! ## A real atom-typed force field — [`typing`], [`oplsaa`], [`parameterize`]
//!
//! The positional [`ForceField`] above is the *generic* path: a caller
//! supplies σ/ε and bonded constants directly. Commercial MD
//! (GROMACS, AMBER, OpenMM) instead uses a **validated, atom-typed
//! force field** — a parameter database keyed by chemically-meaningful
//! atom types, plus an *atom-typer* that assigns those types from a
//! molecule's elements and bonded connectivity.
//!
//! This crate ships a faithful representative subset of **OPLS-AA**
//! (Jorgensen, Maxwell & Tirado-Rives, *JACS* 1996):
//!
//! - [`typing`] — atom-type perception: elements + bonded
//!   connectivity + perceived hybridization → OPLS-AA atom types.
//! - [`oplsaa`] — the OPLS-AA parameter database: per-type LJ σ/ε and
//!   partial charge, and bond / angle / dihedral / improper
//!   parameters keyed by atom-type tuples.
//! - [`parameterize`] — [`parameterize`](parameterize::parameterize)
//!   runs the typer, looks up the database, and returns a fully
//!   populated [`ForceField`] + per-atom charges. This is the default
//!   path where the force field can type a molecule; the generic path
//!   stays available for everything it cannot.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

use crate::error::{MdError, Result};
use crate::system::Topology;

pub mod oplsaa;
pub mod parameterize;
pub mod typing;

/// How two unlike Lennard-Jones atom types are combined into a pair
/// interaction.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum CombiningRule {
    /// Lorentz-Berthelot: `σ_ij = (σ_i + σ_j)/2`, `ε_ij = √(ε_i·ε_j)`.
    /// The AMBER / CHARMM / OPLS-AA convention.
    LorentzBerthelot,
    /// Geometric: `σ_ij = √(σ_i·σ_j)`, `ε_ij = √(ε_i·ε_j)`. GROMACS
    /// `comb-rule 3`, used by the GROMOS and OPLS united-atom sets.
    Geometric,
}

impl CombiningRule {
    /// Combines two single-type LJ parameters into a pair `(σ, ε)`.
    pub fn combine(&self, a: LjParam, b: LjParam) -> LjParam {
        let epsilon = (a.epsilon * b.epsilon).max(0.0).sqrt();
        let sigma = match self {
            CombiningRule::LorentzBerthelot => 0.5 * (a.sigma + b.sigma),
            CombiningRule::Geometric => (a.sigma * b.sigma).max(0.0).sqrt(),
        };
        LjParam { sigma, epsilon }
    }
}

/// Lennard-Jones parameters for one atom type or one pair.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct LjParam {
    /// Finite-distance zero-crossing σ (nm).
    pub sigma: f64,
    /// Well depth ε (kJ/mol).
    pub epsilon: f64,
}

impl LjParam {
    /// Builds an LJ parameter pair.
    ///
    /// # Errors
    /// [`MdError::Invalid`] if σ is negative or ε is negative / not
    /// finite. σ = 0 (a charge-only "dummy") is allowed.
    pub fn new(sigma: f64, epsilon: f64) -> Result<Self> {
        if !(sigma.is_finite() && sigma >= 0.0) {
            return Err(MdError::invalid("sigma", "must be finite and non-negative"));
        }
        if !(epsilon.is_finite() && epsilon >= 0.0) {
            return Err(MdError::invalid(
                "epsilon",
                "must be finite and non-negative",
            ));
        }
        Ok(LjParam { sigma, epsilon })
    }

    /// The `c6`/`c12` form: `c6 = 4εσ⁶`, `c12 = 4εσ¹²`.
    pub fn c6_c12(&self) -> (f64, f64) {
        let s6 = self.sigma.powi(6);
        (4.0 * self.epsilon * s6, 4.0 * self.epsilon * s6 * s6)
    }
}

/// Harmonic-bond parameters: `V = ½·k·(r − r₀)²`.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BondParam {
    /// Equilibrium length r₀ (nm).
    pub r0: f64,
    /// Force constant k (kJ/(mol·nm²)).
    pub k: f64,
}

impl BondParam {
    /// Builds a harmonic-bond parameter set.
    ///
    /// # Errors
    /// [`MdError::Invalid`] on non-finite / negative values.
    pub fn new(r0: f64, k: f64) -> Result<Self> {
        check_nonneg("bond.r0", r0)?;
        check_nonneg("bond.k", k)?;
        Ok(BondParam { r0, k })
    }
}

/// Harmonic-angle parameters: `V = ½·k·(θ − θ₀)²`, θ in radians.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct AngleParam {
    /// Equilibrium angle θ₀ (radians).
    pub theta0: f64,
    /// Force constant k (kJ/(mol·rad²)).
    pub k: f64,
}

impl AngleParam {
    /// Builds a harmonic-angle parameter set. `theta0` is in radians.
    ///
    /// # Errors
    /// [`MdError::Invalid`] on non-finite values or a negative k.
    pub fn new(theta0: f64, k: f64) -> Result<Self> {
        if !theta0.is_finite() {
            return Err(MdError::invalid("angle.theta0", "must be finite"));
        }
        check_nonneg("angle.k", k)?;
        Ok(AngleParam { theta0, k })
    }
}

/// Proper-dihedral parameters. Two functional forms share the struct,
/// tagged by [`DihedralKind`].
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DihedralParam {
    /// Which functional form the coefficients describe.
    pub kind: DihedralKind,
}

/// The two supported proper-dihedral functional forms.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum DihedralKind {
    /// Periodic / cosine form: `V = k·(1 + cos(n·φ − φ₀))`. The
    /// standard AMBER / CHARMM proper torsion.
    Periodic {
        /// Barrier height k (kJ/mol).
        k: f64,
        /// Multiplicity n (an integer, ≥ 1).
        multiplicity: u32,
        /// Phase φ₀ (radians).
        phase: f64,
    },
    /// Ryckaert-Bellemans form:
    /// `V = Σ_{n=0}^{5} cₙ·cosⁿ(ψ)` with `ψ = φ − 180°`. The GROMOS
    /// alkane torsion.
    RyckaertBellemans {
        /// The six RB coefficients c₀..c₅ (kJ/mol).
        c: [f64; 6],
    },
}

impl DihedralParam {
    /// A periodic (cosine) proper dihedral.
    ///
    /// # Errors
    /// [`MdError::Invalid`] on a non-finite k / phase or a zero
    /// multiplicity.
    pub fn periodic(k: f64, multiplicity: u32, phase: f64) -> Result<Self> {
        if !k.is_finite() {
            return Err(MdError::invalid("dihedral.k", "must be finite"));
        }
        if multiplicity == 0 {
            return Err(MdError::invalid(
                "dihedral.multiplicity",
                "must be at least 1",
            ));
        }
        if !phase.is_finite() {
            return Err(MdError::invalid("dihedral.phase", "must be finite"));
        }
        Ok(DihedralParam {
            kind: DihedralKind::Periodic {
                k,
                multiplicity,
                phase,
            },
        })
    }

    /// A Ryckaert-Bellemans proper dihedral.
    ///
    /// # Errors
    /// [`MdError::Invalid`] if any coefficient is non-finite.
    pub fn ryckaert_bellemans(c: [f64; 6]) -> Result<Self> {
        if c.iter().any(|x| !x.is_finite()) {
            return Err(MdError::invalid(
                "dihedral.c",
                "all RB coefficients must be finite",
            ));
        }
        Ok(DihedralParam {
            kind: DihedralKind::RyckaertBellemans { c },
        })
    }
}

/// Harmonic improper-dihedral parameters: `V = ½·k·(ξ − ξ₀)²`.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ImproperParam {
    /// Equilibrium improper angle ξ₀ (radians) — usually 0.
    pub xi0: f64,
    /// Force constant k (kJ/(mol·rad²)).
    pub k: f64,
}

impl ImproperParam {
    /// Builds a harmonic improper-dihedral parameter set.
    ///
    /// # Errors
    /// [`MdError::Invalid`] on non-finite values or a negative k.
    pub fn new(xi0: f64, k: f64) -> Result<Self> {
        if !xi0.is_finite() {
            return Err(MdError::invalid("improper.xi0", "must be finite"));
        }
        check_nonneg("improper.k", k)?;
        Ok(ImproperParam { xi0, k })
    }
}

/// A complete force field: LJ parameters keyed by atom-type name,
/// bonded parameters held positionally, plus the global combining
/// rule and 1-4 scaling factors.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ForceField {
    /// LJ parameters per atom-type name.
    lj: HashMap<String, LjParam>,
    /// Bonded-term parameters, indexed parallel to the topology lists.
    bonds: Vec<BondParam>,
    angles: Vec<AngleParam>,
    dihedrals: Vec<DihedralParam>,
    impropers: Vec<ImproperParam>,
    /// Combining rule for unlike LJ pairs.
    pub combining_rule: CombiningRule,
    /// Scale factor applied to 1-4 (dihedral-end) Lennard-Jones
    /// interactions. AMBER uses 0.5.
    pub lj_14_scale: f64,
    /// Scale factor applied to 1-4 Coulomb interactions. AMBER uses
    /// 1/1.2 ≈ 0.8333.
    pub coulomb_14_scale: f64,
}

impl ForceField {
    /// An empty force field with the given combining rule and default
    /// AMBER-style 1-4 scaling (LJ 0.5, Coulomb 0.8333).
    pub fn new(combining_rule: CombiningRule) -> Self {
        ForceField {
            lj: HashMap::new(),
            bonds: Vec::new(),
            angles: Vec::new(),
            dihedrals: Vec::new(),
            impropers: Vec::new(),
            combining_rule,
            lj_14_scale: 0.5,
            coulomb_14_scale: 1.0 / 1.2,
        }
    }

    /// Registers (or overwrites) the LJ parameters for an atom type.
    pub fn set_lj(&mut self, type_name: impl Into<String>, param: LjParam) {
        self.lj.insert(type_name.into(), param);
    }

    /// Looks up the single-type LJ parameters for an atom type.
    ///
    /// # Errors
    /// [`MdError::Invalid`] if the type is not registered.
    pub fn lj(&self, type_name: &str) -> Result<LjParam> {
        self.lj.get(type_name).copied().ok_or_else(|| {
            MdError::invalid(
                "atom-type",
                format!("no Lennard-Jones parameters registered for type `{type_name}`"),
            )
        })
    }

    /// The combined LJ parameters for an *unlike* pair of types.
    ///
    /// # Errors
    /// [`MdError::Invalid`] if either type is not registered.
    pub fn lj_pair(&self, a: &str, b: &str) -> Result<LjParam> {
        Ok(self.combining_rule.combine(self.lj(a)?, self.lj(b)?))
    }

    /// Appends a bond parameter (positionally matched to the topology's
    /// next bond).
    pub fn push_bond(&mut self, param: BondParam) {
        self.bonds.push(param);
    }

    /// Appends an angle parameter.
    pub fn push_angle(&mut self, param: AngleParam) {
        self.angles.push(param);
    }

    /// Appends a proper-dihedral parameter.
    pub fn push_dihedral(&mut self, param: DihedralParam) {
        self.dihedrals.push(param);
    }

    /// Appends an improper-dihedral parameter.
    pub fn push_improper(&mut self, param: ImproperParam) {
        self.impropers.push(param);
    }

    /// The bond parameters, parallel to `topology.bonds`.
    pub fn bonds(&self) -> &[BondParam] {
        &self.bonds
    }

    /// The angle parameters, parallel to `topology.angles`.
    pub fn angles(&self) -> &[AngleParam] {
        &self.angles
    }

    /// The proper-dihedral parameters, parallel to `topology.dihedrals`.
    pub fn dihedrals(&self) -> &[DihedralParam] {
        &self.dihedrals
    }

    /// The improper-dihedral parameters, parallel to
    /// `topology.impropers`.
    pub fn impropers(&self) -> &[ImproperParam] {
        &self.impropers
    }

    /// Checks that this force field fully covers `topology`: every atom
    /// type has LJ parameters, and the four bonded-parameter lists are
    /// each the same length as the matching topology list.
    ///
    /// # Errors
    /// [`MdError::DimensionMismatch`] for a length disagreement,
    /// [`MdError::Invalid`] for a missing atom type.
    pub fn validate_against(&self, topology: &Topology) -> Result<()> {
        for (n, atom) in topology.atoms.iter().enumerate() {
            if !self.lj.contains_key(&atom.type_name) {
                return Err(MdError::invalid(
                    "atom-type",
                    format!(
                        "atom {n} has type `{}` with no LJ parameters",
                        atom.type_name
                    ),
                ));
            }
        }
        let pairs = [
            ("bonds", self.bonds.len(), topology.bonds.len()),
            ("angles", self.angles.len(), topology.angles.len()),
            ("dihedrals", self.dihedrals.len(), topology.dihedrals.len()),
            ("impropers", self.impropers.len(), topology.impropers.len()),
        ];
        for (name, have, need) in pairs {
            if have != need {
                return Err(MdError::dimension(format!(
                    "{name}: force field has {have} parameter sets but topology has {need}"
                )));
            }
        }
        Ok(())
    }
}

fn check_nonneg(what: &'static str, v: f64) -> Result<()> {
    if !(v.is_finite() && v >= 0.0) {
        Err(MdError::invalid(what, "must be finite and non-negative"))
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::system::Atom;

    #[test]
    fn lorentz_berthelot_combination() {
        let a = LjParam::new(0.3, 0.5).unwrap();
        let b = LjParam::new(0.5, 0.8).unwrap();
        let c = CombiningRule::LorentzBerthelot.combine(a, b);
        assert!((c.sigma - 0.4).abs() < 1e-12);
        assert!((c.epsilon - (0.5_f64 * 0.8).sqrt()).abs() < 1e-12);
    }

    #[test]
    fn geometric_combination() {
        let a = LjParam::new(0.4, 1.0).unwrap();
        let b = LjParam::new(0.9, 4.0).unwrap();
        let c = CombiningRule::Geometric.combine(a, b);
        assert!((c.sigma - (0.4_f64 * 0.9).sqrt()).abs() < 1e-12);
        assert!((c.epsilon - 2.0).abs() < 1e-12);
    }

    #[test]
    fn c6_c12_round_trip() {
        let p = LjParam::new(0.35, 0.65).unwrap();
        let (c6, c12) = p.c6_c12();
        // Recover sigma from c12/c6 = sigma^6.
        let s6 = c12 / c6;
        assert!((s6.powf(1.0 / 6.0) - 0.35).abs() < 1e-9);
    }

    #[test]
    fn lj_params_reject_bad_input() {
        assert!(LjParam::new(-0.1, 0.5).is_err());
        assert!(LjParam::new(0.3, -0.5).is_err());
        assert!(LjParam::new(0.0, 0.5).is_ok()); // dummy / charge-only
    }

    #[test]
    fn dihedral_constructors_validate() {
        assert!(DihedralParam::periodic(2.0, 0, 0.0).is_err());
        assert!(DihedralParam::periodic(2.0, 2, 0.0).is_ok());
        assert!(DihedralParam::ryckaert_bellemans([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]).is_ok());
        assert!(DihedralParam::ryckaert_bellemans([f64::NAN, 0.0, 0.0, 0.0, 0.0, 0.0]).is_err());
    }

    #[test]
    fn validate_against_detects_missing_type_and_count_mismatch() {
        let mut top = Topology::new();
        top.push_atom(Atom::new("CT", 12.0, 0.0).unwrap());
        top.push_atom(Atom::new("HC", 1.0, 0.0).unwrap());
        top.add_bond(0, 1).unwrap();

        let mut ff = ForceField::new(CombiningRule::LorentzBerthelot);
        ff.set_lj("CT", LjParam::new(0.34, 0.45).unwrap());
        // HC missing -> error.
        assert!(ff.validate_against(&top).is_err());

        ff.set_lj("HC", LjParam::new(0.26, 0.08).unwrap());
        // bond count mismatch (0 params, 1 bond).
        assert!(ff.validate_against(&top).is_err());

        ff.push_bond(BondParam::new(0.11, 3e5).unwrap());
        assert!(ff.validate_against(&top).is_ok());
    }

    #[test]
    fn lj_pair_errors_on_unknown_type() {
        let ff = ForceField::new(CombiningRule::Geometric);
        assert!(ff.lj_pair("A", "B").is_err());
    }
}
