//! Air-fuel ratio, equivalence ratio, and product mole balance.
//!
//! ## Model
//!
//! Complete combustion of one mole of `CxHy` with `a = x + y/4` moles
//! of O2 (and its inert N2 carrier):
//!
//! ```text
//! CxHy + a (O2 + 3.76 N2) -> x CO2 + (y/2) H2O + 3.76 a N2
//! ```
//!
//! With excess air (lean, `phi < 1`) the unburned O2 appears in the
//! products; with `phi > 1` (rich) there is not enough air for complete
//! combustion, so this closed-form complete-combustion balance no longer
//! holds and [`product_moles`] reports the rich case as unsupported.

use serde::{Deserialize, Serialize};

use crate::error::CombustionError;
use crate::fuel::{Fuel, M_AIR_PER_MOLE, M_N2, M_O2, N2_PER_O2};

/// Combustion-regime classification from the equivalence ratio.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Mixture {
    /// `phi < 1`: excess air, fuel-lean.
    Lean,
    /// `phi == 1` (within tolerance): chemically correct.
    Stoichiometric,
    /// `phi > 1`: excess fuel, fuel-rich.
    Rich,
}

/// Stoichiometric air-fuel ratio of `fuel`, by **mass** (kg air / kg
/// fuel).
///
/// `AFR_stoich = a * (M_O2 + 3.76 M_N2) / M_fuel`, with `a = x + y/4`.
///
/// For methane (CH4) this is approximately `17.1`.
pub fn afr_stoich_mass(fuel: &Fuel) -> f64 {
    let a = fuel.stoich_o2_moles();
    let air_mass = a * (M_O2 + N2_PER_O2 * M_N2);
    air_mass / fuel.molar_mass()
}

/// Stoichiometric air-fuel ratio of `fuel`, by **mole** (mol air / mol
/// fuel), where one mole of "air" is `O2 + 3.76 N2`.
///
/// `AFR_molar = a * (1 + 3.76)`, with `a = x + y/4`.
pub fn afr_stoich_molar(fuel: &Fuel) -> f64 {
    fuel.stoich_o2_moles() * (1.0 + N2_PER_O2)
}

/// Equivalence ratio `phi = AFR_stoich / AFR_actual` (mass basis).
///
/// `phi = 1` at stoichiometric, `phi > 1` fuel-rich, `phi < 1`
/// fuel-lean. Because the ratio is taken on a consistent (mass) basis
/// the molar AFR would give the identical `phi`.
///
/// # Errors
///
/// Returns [`CombustionError::BadParameter`] when `afr_actual_mass` is
/// not strictly positive.
pub fn equivalence_ratio(fuel: &Fuel, afr_actual_mass: f64) -> Result<f64, CombustionError> {
    CombustionError::require_positive(
        "afr_actual_mass",
        afr_actual_mass,
        "actual air-fuel ratio must be strictly positive",
    )?;
    Ok(afr_stoich_mass(fuel) / afr_actual_mass)
}

/// Actual mass air-fuel ratio that produces a target equivalence ratio
/// `phi`.
///
/// Inverse of [`equivalence_ratio`]: `AFR_actual = AFR_stoich / phi`.
///
/// # Errors
///
/// Returns [`CombustionError::BadParameter`] when `phi` is not strictly
/// positive.
pub fn afr_from_phi(fuel: &Fuel, phi: f64) -> Result<f64, CombustionError> {
    CombustionError::require_positive("phi", phi, "equivalence ratio must be strictly positive")?;
    Ok(afr_stoich_mass(fuel) / phi)
}

/// Percent **theoretical air** supplied at equivalence ratio `phi`:
/// `100 / phi` (percent).
///
/// This is the combustion-engineering way of stating mixture strength on
/// the air side. `100%` is stoichiometric; a lean mixture (`phi < 1`)
/// supplies more than `100%` theoretical air, and a rich mixture
/// (`phi > 1`) less than `100%` (an air deficiency).
///
/// # Errors
///
/// Returns [`CombustionError::BadParameter`] when `phi` is not strictly
/// positive (or non-finite).
pub fn percent_theoretical_air(phi: f64) -> Result<f64, CombustionError> {
    CombustionError::require_positive("phi", phi, "equivalence ratio must be strictly positive")?;
    Ok(100.0 / phi)
}

/// Percent **excess air** supplied at equivalence ratio `phi`:
/// `100 * (1 - phi) / phi` (percent), i.e. [`percent_theoretical_air`]
/// minus `100`.
///
/// `0%` at stoichiometric, positive for a lean mixture (`phi < 1`), and
/// negative — an air deficiency — for a rich mixture (`phi > 1`). It ties
/// back to the air-fuel ratio through
/// `AFR_actual = AFR_stoich * (1 + excess_air/100)`.
///
/// # Errors
///
/// Returns [`CombustionError::BadParameter`] when `phi` is not strictly
/// positive (or non-finite).
pub fn percent_excess_air(phi: f64) -> Result<f64, CombustionError> {
    CombustionError::require_positive("phi", phi, "equivalence ratio must be strictly positive")?;
    Ok(100.0 * (1.0 - phi) / phi)
}

/// Classify a mixture from its equivalence ratio `phi`.
///
/// Values within `tol` of `1.0` are reported as
/// [`Mixture::Stoichiometric`].
///
/// # Errors
///
/// Returns [`CombustionError::BadParameter`] when `phi <= 0` or
/// `tol < 0`.
pub fn classify(phi: f64, tol: f64) -> Result<Mixture, CombustionError> {
    CombustionError::require_positive("phi", phi, "equivalence ratio must be strictly positive")?;
    if tol < 0.0 {
        return Err(CombustionError::BadParameter {
            name: "tol",
            value: tol,
            reason: "tolerance must be non-negative",
        });
    }
    if (phi - 1.0).abs() <= tol {
        Ok(Mixture::Stoichiometric)
    } else if phi > 1.0 {
        Ok(Mixture::Rich)
    } else {
        Ok(Mixture::Lean)
    }
}

/// Molar amounts of each combustion product per mole of fuel.
///
/// All quantities are moles per mole of fuel burned.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ProductMoles {
    /// Carbon dioxide, `x`.
    pub co2: f64,
    /// Water vapour, `y / 2`.
    pub h2o: f64,
    /// Inert nitrogen carried through, `3.76 * a_supplied`.
    pub n2: f64,
    /// Excess (unreacted) oxygen, `a_supplied - a_stoich` (>= 0 for the
    /// lean / stoichiometric cases this function supports).
    pub excess_o2: f64,
}

impl ProductMoles {
    /// Total moles of product gas per mole of fuel.
    pub fn total(&self) -> f64 {
        self.co2 + self.h2o + self.n2 + self.excess_o2
    }
}

/// Balanced product moles per mole of fuel at equivalence ratio `phi`,
/// for the **lean and stoichiometric** regimes (`phi <= 1`).
///
/// The supplied O2 is `a_stoich / phi`, so excess O2 is
/// `a_stoich (1/phi - 1)` and the N2 carried is `3.76 * a_stoich / phi`.
/// Carbon and hydrogen burn completely to CO2 and H2O.
///
/// # Errors
///
/// Returns [`CombustionError::BadParameter`] when `phi <= 0`, or when
/// `phi > 1` (the rich regime is outside this closed-form complete-
/// combustion balance, which cannot allocate the products without
/// equilibrium chemistry for CO / H2 / soot).
pub fn product_moles(fuel: &Fuel, phi: f64) -> Result<ProductMoles, CombustionError> {
    CombustionError::require_positive("phi", phi, "equivalence ratio must be strictly positive")?;
    if phi > 1.0 {
        return Err(CombustionError::BadParameter {
            name: "phi",
            value: phi,
            reason: "rich combustion (phi > 1) is outside the complete-combustion product balance",
        });
    }
    let a_stoich = fuel.stoich_o2_moles();
    let a_supplied = a_stoich / phi;
    Ok(ProductMoles {
        co2: fuel.carbon as f64,
        h2o: fuel.hydrogen as f64 / 2.0,
        n2: N2_PER_O2 * a_supplied,
        excess_o2: a_supplied - a_stoich,
    })
}

/// Mean molar mass of the air mixture (`O2 + 3.76 N2`), g/mol — a
/// re-export of the fuel-module constant so callers can size air masses
/// without reaching across modules.
pub fn air_mean_molar_mass() -> f64 {
    M_AIR_PER_MOLE
}
