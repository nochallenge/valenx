//! Exchange-correlation functionals.
//!
//! This module is the exchange-correlation half of the Kohn-Sham DFT
//! subsystem. It evaluates, at a point of given density (and density
//! gradient), the XC energy density and the XC potential — the two
//! quantities the Kohn-Sham energy and the Kohn-Sham matrix need.
//!
//! - [`lda`] — the local density approximation: Slater (Dirac)
//!   exchange + VWN5 correlation (the `"SVWN"` / `"LDA"` functional).
//! - [`gga`] — generalised-gradient functionals: **PBE** exchange and
//!   correlation.
//! - [`b3lyp`] — the B88 exchange and LYP correlation pieces and the
//!   B3LYP hybrid mixing.
//!
//! ## The unit of currency — [`XcContribution`]
//!
//! Every functional returns an [`XcContribution`] at a point: the
//! energy density per electron `ε_xc` (so `E_xc = ∫ ρ ε_xc dr`), the
//! density potential `∂(ρε_xc)/∂ρ`, and — for a GGA — the
//! gradient potential `∂(ρε_xc)/∂|∇ρ|`. The Kohn-Sham build in
//! [`crate::dft::ks`] consumes exactly these three numbers.

pub mod b3lyp;
pub mod gga;
pub mod lda;

/// The exchange-correlation contribution at a single grid point.
///
/// Splitting the XC into an energy density and the two potential
/// components lets the Kohn-Sham energy build (`E_xc = Σ w ρ ε_xc`) and
/// the Kohn-Sham matrix build (which needs the potentials) share one
/// functional evaluation.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct XcContribution {
    /// Exchange-correlation energy density *per electron* `ε_xc`. The
    /// XC energy is `E_xc = ∫ ρ ε_xc dr`.
    pub energy_density: f64,
    /// The density part of the XC potential — `∂(ρε_xc)/∂ρ` (an LDA's
    /// whole `v_xc`; a GGA's `∂f/∂ρ` term).
    pub potential: f64,
    /// The gradient part of the XC potential — `∂(ρε_xc)/∂|∇ρ|`. Zero
    /// for an LDA; the term a GGA's divergence contribution needs.
    pub gradient_potential: f64,
}

impl XcContribution {
    /// The all-zero contribution — returned for a vanishing density.
    pub const fn zero() -> Self {
        XcContribution {
            energy_density: 0.0,
            potential: 0.0,
            gradient_potential: 0.0,
        }
    }

    /// Add two contributions component-wise (e.g. exchange + correlation).
    pub fn add(&self, other: &XcContribution) -> XcContribution {
        XcContribution {
            energy_density: self.energy_density + other.energy_density,
            potential: self.potential + other.potential,
            gradient_potential: self.gradient_potential + other.gradient_potential,
        }
    }

    /// Subtract `other` from `self` component-wise.
    pub fn sub(&self, other: &XcContribution) -> XcContribution {
        XcContribution {
            energy_density: self.energy_density - other.energy_density,
            potential: self.potential - other.potential,
            gradient_potential: self.gradient_potential - other.gradient_potential,
        }
    }

    /// Scale every component by `factor` (the hybrid-mixing operation).
    pub fn scale(&self, factor: f64) -> XcContribution {
        XcContribution {
            energy_density: self.energy_density * factor,
            potential: self.potential * factor,
            gradient_potential: self.gradient_potential * factor,
        }
    }
}

/// The exchange-correlation functional a DFT calculation uses.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Functional {
    /// Local density approximation — Slater exchange + VWN5 correlation
    /// (`"SVWN"` / `"LDA"`). No density-gradient dependence.
    Lda,
    /// PBE — the Perdew-Burke-Ernzerhof generalised-gradient
    /// functional. Needs the density gradient.
    Pbe,
    /// B3LYP — the Becke three-parameter hybrid (B88 + LYP + 20 %
    /// exact exchange). Needs the density gradient *and* the
    /// Hartree-Fock exchange matrix.
    B3lyp,
}

impl Functional {
    /// Parse a functional name (case-insensitive). Accepts the common
    /// aliases — `"lda"` / `"svwn"` / `"slater"`, `"pbe"`,
    /// `"b3lyp"` / `"b3-lyp"`.
    pub fn from_name(name: &str) -> Option<Functional> {
        match name.trim().to_ascii_lowercase().as_str() {
            "lda" | "svwn" | "svwn5" | "slater" | "s-vwn" => Some(Functional::Lda),
            "pbe" | "pbepbe" => Some(Functional::Pbe),
            "b3lyp" | "b3-lyp" | "becke3lyp" => Some(Functional::B3lyp),
            _ => None,
        }
    }

    /// A short human-readable label.
    pub fn label(self) -> &'static str {
        match self {
            Functional::Lda => "LDA (SVWN5)",
            Functional::Pbe => "PBE",
            Functional::B3lyp => "B3LYP",
        }
    }

    /// `true` when the functional is a GGA / hybrid that needs the
    /// density gradient on the grid.
    pub fn needs_gradient(self) -> bool {
        match self {
            Functional::Lda => false,
            Functional::Pbe | Functional::B3lyp => true,
        }
    }

    /// The fraction of exact (Hartree-Fock) exchange the functional
    /// mixes in — `0` for a pure functional, `0.2` for B3LYP.
    pub fn exact_exchange_fraction(self) -> f64 {
        match self {
            Functional::Lda | Functional::Pbe => 0.0,
            Functional::B3lyp => b3lyp::B3LYP_A_EXACT,
        }
    }

    /// Evaluate the *density-functional* exchange-correlation
    /// contribution at a point of density `rho` and gradient norm
    /// `grad_norm`.
    ///
    /// For B3LYP this is everything *except* the exact-exchange term —
    /// the Kohn-Sham build adds `exact_exchange_fraction · K` itself.
    pub fn evaluate(self, rho: f64, grad_norm: f64) -> XcContribution {
        match self {
            Functional::Lda => lda::lda_xc(rho),
            Functional::Pbe => gga::pbe_xc(rho, grad_norm),
            Functional::B3lyp => b3lyp::b3lyp_dft_xc(rho, grad_norm),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn contribution_arithmetic() {
        let a = XcContribution {
            energy_density: 1.0,
            potential: 2.0,
            gradient_potential: 3.0,
        };
        let b = XcContribution {
            energy_density: 0.5,
            potential: 1.0,
            gradient_potential: 1.5,
        };
        let s = a.add(&b);
        assert_eq!(s.energy_density, 1.5);
        assert_eq!(s.potential, 3.0);
        let d = a.sub(&b);
        assert_eq!(d.energy_density, 0.5);
        let sc = a.scale(2.0);
        assert_eq!(sc.gradient_potential, 6.0);
    }

    #[test]
    fn functional_name_parsing() {
        assert_eq!(Functional::from_name("LDA"), Some(Functional::Lda));
        assert_eq!(Functional::from_name("svwn"), Some(Functional::Lda));
        assert_eq!(Functional::from_name("PBE"), Some(Functional::Pbe));
        assert_eq!(Functional::from_name("b3lyp"), Some(Functional::B3lyp));
        assert_eq!(Functional::from_name("B3-LYP"), Some(Functional::B3lyp));
        assert_eq!(Functional::from_name("ccsd"), None);
    }

    #[test]
    fn gradient_and_exchange_flags() {
        assert!(!Functional::Lda.needs_gradient());
        assert!(Functional::Pbe.needs_gradient());
        assert!(Functional::B3lyp.needs_gradient());
        assert_eq!(Functional::Lda.exact_exchange_fraction(), 0.0);
        assert_eq!(Functional::Pbe.exact_exchange_fraction(), 0.0);
        assert!((Functional::B3lyp.exact_exchange_fraction() - 0.20).abs() < 1.0e-12);
    }

    /// Dispatch through `evaluate` matches calling the functional
    /// directly.
    #[test]
    fn evaluate_dispatches_correctly() {
        let lda = Functional::Lda.evaluate(1.0, 0.0);
        assert_eq!(lda, lda::lda_xc(1.0));
        let pbe = Functional::Pbe.evaluate(1.0, 1.5);
        assert_eq!(pbe, gga::pbe_xc(1.0, 1.5));
    }
}
