//! Constant-pressure adiabatic flame temperature (mean-cp model).
//!
//! ## Model
//!
//! A first-law energy balance at constant pressure with **no heat loss**
//! and a **single constant mean molar heat capacity** for the product
//! gas. All chemical energy released — the fuel lower heating value
//! (LHV) times the fuel mass — goes into raising the product mixture
//! from the reactant inlet temperature:
//!
//! ```text
//! Q_released = n_products * cp_molar * (T_ad - T_in)
//! =>  T_ad = T_in + Q_released / (n_products * cp_molar)
//! ```
//!
//! where `Q_released = LHV * M_fuel` (J per mole of fuel) and
//! `n_products` is the total product moles per mole of fuel from
//! [`crate::stoich::product_moles`].
//!
//! ## Honest scope
//!
//! This is a textbook teaching estimate. It deliberately ignores
//! dissociation (CO2 <-> CO + ½O2, H2O <-> H2 + ½O2), the temperature
//! dependence of cp (real cp climbs from ~29 to ~60 J/mol/K across the
//! range), and finite-rate kinetics. Real adiabatic flame temperatures
//! are several hundred kelvin lower than the no-dissociation balance.
//! Treat the result as a comparative, order-of-magnitude figure only.

use crate::error::CombustionError;
use crate::fuel::Fuel;
use crate::stoich::product_moles;

/// Standard reference inlet temperature, 298.15 K (25 degrees C).
pub const T_REF_K: f64 = 298.15;

/// A representative mean molar heat capacity of burned-gas products,
/// J/(mol*K). A round value between the cold (~29) and hot (~55) limits;
/// callers should pass their own value via
/// [`adiabatic_flame_temperature`] when accuracy matters.
pub const CP_MOLAR_DEFAULT: f64 = 40.0;

/// Lower heating value of methane, J/kg (LHV, ~50.0 MJ/kg).
pub const LHV_METHANE: f64 = 50.0e6;
/// Lower heating value of propane, J/kg (~46.4 MJ/kg).
pub const LHV_PROPANE: f64 = 46.4e6;
/// Lower heating value of iso-octane, J/kg (~44.3 MJ/kg).
pub const LHV_OCTANE: f64 = 44.3e6;

/// Heat released per mole of fuel burned, joules.
///
/// `Q = LHV [J/kg] * M_fuel [kg/mol] = LHV * (molar_mass_g / 1000)`.
///
/// # Errors
///
/// Returns [`CombustionError::BadParameter`] when `lhv_j_per_kg` is not
/// strictly positive.
pub fn heat_release_per_mole_fuel(fuel: &Fuel, lhv_j_per_kg: f64) -> Result<f64, CombustionError> {
    CombustionError::require_positive(
        "lhv_j_per_kg",
        lhv_j_per_kg,
        "lower heating value must be strictly positive",
    )?;
    let molar_mass_kg = fuel.molar_mass() / 1000.0;
    Ok(lhv_j_per_kg * molar_mass_kg)
}

/// Adiabatic flame temperature, kelvin, from the constant-pressure
/// mean-cp energy balance.
///
/// Inputs:
/// - `fuel`: the hydrocarbon being burned.
/// - `phi`: equivalence ratio (`phi <= 1`; lean / stoichiometric only,
///   matching [`product_moles`]).
/// - `lhv_j_per_kg`: fuel lower heating value, J/kg.
/// - `cp_molar`: mean molar heat capacity of the product gas,
///   J/(mol*K).
/// - `t_in_k`: reactant inlet temperature, K.
///
/// # Errors
///
/// Propagates [`CombustionError`] from [`product_moles`] /
/// [`heat_release_per_mole_fuel`], and returns
/// [`CombustionError::BadParameter`] when `cp_molar <= 0` or
/// `t_in_k <= 0`.
pub fn adiabatic_flame_temperature(
    fuel: &Fuel,
    phi: f64,
    lhv_j_per_kg: f64,
    cp_molar: f64,
    t_in_k: f64,
) -> Result<f64, CombustionError> {
    CombustionError::require_positive(
        "cp_molar",
        cp_molar,
        "mean molar heat capacity must be strictly positive",
    )?;
    CombustionError::require_positive(
        "t_in_k",
        t_in_k,
        "inlet temperature must be strictly positive",
    )?;
    let q = heat_release_per_mole_fuel(fuel, lhv_j_per_kg)?;
    let products = product_moles(fuel, phi)?;
    let n_products = products.total();
    Ok(t_in_k + q / (n_products * cp_molar))
}
