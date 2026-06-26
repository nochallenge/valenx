//! Pure-component definitions: the critical constants and acentric factor that
//! parametrize a cubic equation of state.
//!
//! A [`Fluid`] carries the minimum data a Peng–Robinson or Soave–Redlich–Kwong
//! model needs: the critical temperature `tc` (K), critical pressure `pc` (Pa),
//! and Pitzer acentric factor `omega` (dimensionless). A small library of
//! common fluids with literature constants is provided for convenience and for
//! validation against known ground truth.

use crate::error::{require_positive, Result, ThermoError};

/// A pure fluid described by its critical constants and acentric factor.
///
/// These three numbers are sufficient to construct a two-parameter cubic
/// equation of state (Peng–Robinson or SRK). Construct one with [`Fluid::new`]
/// (which validates the inputs) or pick a predefined fluid from the library
/// functions on this type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Fluid {
    /// Component name (for display / diagnostics).
    pub name: &'static str,
    /// Critical temperature, in kelvin. Must be strictly positive.
    pub tc: f64,
    /// Critical pressure, in pascal. Must be strictly positive.
    pub pc: f64,
    /// Pitzer acentric factor (dimensionless). Physically in roughly
    /// `[-0.5, 1.5]`; enforced to lie in `[-1, 2]`.
    pub omega: f64,
}

impl Fluid {
    /// Construct a fluid from critical constants and acentric factor, validating
    /// that `tc` and `pc` are strictly positive and `omega` is within `[-1, 2]`.
    ///
    /// # Errors
    ///
    /// Returns [`ThermoError::NonPositive`] if `tc` or `pc` is not strictly
    /// positive, or [`ThermoError::OutOfRange`] if `omega` is outside `[-1, 2]`.
    pub fn new(name: &'static str, tc: f64, pc: f64, omega: f64) -> Result<Self> {
        require_positive("critical temperature (tc)", tc)?;
        require_positive("critical pressure (pc)", pc)?;
        if !omega.is_finite() || !(-1.0..=2.0).contains(&omega) {
            return Err(ThermoError::OutOfRange {
                name: "acentric factor (omega)",
                value: omega,
                expected: "in [-1, 2]",
            });
        }
        Ok(Fluid {
            name,
            tc,
            pc,
            omega,
        })
    }

    /// Reduced temperature `T / Tc` for a given temperature `t` (K).
    #[must_use]
    pub fn reduced_temperature(&self, t: f64) -> f64 {
        t / self.tc
    }

    /// Carbon dioxide (CO₂). Tc = 304.13 K, Pc = 7.3773 MPa, ω = 0.2236.
    ///
    /// Constants from the NIST/DIPPR compilations (Span–Wagner critical point).
    #[must_use]
    pub fn carbon_dioxide() -> Self {
        Fluid {
            name: "carbon dioxide",
            tc: 304.13,
            pc: 7.377_3e6,
            omega: 0.223_6,
        }
    }

    /// Nitrogen (N₂). Tc = 126.19 K, Pc = 3.3958 MPa, ω = 0.0372.
    #[must_use]
    pub fn nitrogen() -> Self {
        Fluid {
            name: "nitrogen",
            tc: 126.19,
            pc: 3.395_8e6,
            omega: 0.037_2,
        }
    }

    /// Methane (CH₄). Tc = 190.56 K, Pc = 4.5992 MPa, ω = 0.0114.
    #[must_use]
    pub fn methane() -> Self {
        Fluid {
            name: "methane",
            tc: 190.56,
            pc: 4.599_2e6,
            omega: 0.011_4,
        }
    }

    /// Water (H₂O). Tc = 647.096 K, Pc = 22.064 MPa, ω = 0.3443.
    ///
    /// Cubic EoS are poor for water's vapor pressure; included for completeness.
    #[must_use]
    pub fn water() -> Self {
        Fluid {
            name: "water",
            tc: 647.096,
            pc: 22.064e6,
            omega: 0.344_3,
        }
    }

    /// Propane (C₃H₈). Tc = 369.89 K, Pc = 4.2512 MPa, ω = 0.1521.
    #[must_use]
    pub fn propane() -> Self {
        Fluid {
            name: "propane",
            tc: 369.89,
            pc: 4.251_2e6,
            omega: 0.152_1,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn library_fluids_are_valid() {
        for f in [
            Fluid::carbon_dioxide(),
            Fluid::nitrogen(),
            Fluid::methane(),
            Fluid::water(),
            Fluid::propane(),
        ] {
            // Reconstructing through the validating constructor must succeed.
            Fluid::new(f.name, f.tc, f.pc, f.omega).unwrap();
            assert!(f.tc > 0.0 && f.pc > 0.0);
        }
    }

    #[test]
    fn rejects_bad_inputs() {
        assert!(matches!(
            Fluid::new("bad", -1.0, 1e6, 0.1),
            Err(ThermoError::NonPositive { .. })
        ));
        assert!(matches!(
            Fluid::new("bad", 300.0, 0.0, 0.1),
            Err(ThermoError::NonPositive { .. })
        ));
        assert!(matches!(
            Fluid::new("bad", 300.0, 1e6, 9.0),
            Err(ThermoError::OutOfRange { .. })
        ));
    }
}
