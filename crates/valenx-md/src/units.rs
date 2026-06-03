//! Physical constants and the crate's internal unit system.
//!
//! `valenx-md` works in the **GROMACS unit system**, which keeps every
//! number close to order 1 for typical biomolecular input:
//!
//! | Quantity    | Unit                          |
//! |-------------|-------------------------------|
//! | length      | nanometre (nm)                |
//! | time        | picosecond (ps)               |
//! | mass        | atomic mass unit (u = g/mol)  |
//! | charge      | elementary charge (e)         |
//! | energy      | kilojoule per mole (kJ/mol)   |
//! | temperature | kelvin (K)                    |
//! | force       | kJ/(mol·nm)                   |
//! | pressure    | bar                           |
//!
//! In this system velocity is nm/ps, and the relation
//! `KE = 1/2 m v^2` already comes out in kJ/mol with no conversion
//! factor — that consistency is the whole reason for choosing it.
//!
//! The constants below are the CODATA-derived values GROMACS uses.
//! They are deliberately exposed so callers and tests can reason about
//! magnitudes explicitly.

/// Boltzmann constant in kJ/(mol·K).
///
/// `k_B = 1.380649e-23 J/K * 6.02214076e23 /mol / 1000 = 0.00831446...`
pub const BOLTZMANN: f64 = 0.008_314_462_618;

/// Coulomb prefactor `1 / (4·π·ε₀)` in kJ·nm/(mol·e²).
///
/// This is the factor `f` such that the Coulomb energy of two unit
/// charges 1 nm apart is `f` kJ/mol. GROMACS calls it `ONE_4PI_EPS0`.
pub const COULOMB: f64 = 138.935_457_64;

/// Avogadro's number (per mole).
pub const AVOGADRO: f64 = 6.022_140_76e23;

/// Pressure-unit conversion: kJ/(mol·nm³) expressed in bar.
///
/// `1 kJ/(mol·nm³) = 1e30 / 6.02214076e23 / 1e5 bar ≈ 16.6054 bar`.
/// The virial pressure is naturally computed in kJ/(mol·nm³); multiply
/// by this to report bar.
pub const PRESSURE_KJMOLNM3_TO_BAR: f64 = 16.605_390_666;

/// Picoseconds per femtosecond (`1 fs = 0.001 ps`). A typical MD time
/// step is 1–2 fs, i.e. 0.001–0.002 in the crate's internal ps unit.
pub const FS_TO_PS: f64 = 0.001;

/// Converts a temperature to the thermal energy `k_B·T` (kJ/mol).
#[inline]
pub fn thermal_energy(temperature_k: f64) -> f64 {
    BOLTZMANN * temperature_k
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn boltzmann_thermal_energy_at_room_temperature() {
        // k_B·T at 300 K is the famous ~2.49 kJ/mol (~0.6 kcal/mol).
        let kt = thermal_energy(300.0);
        assert!((kt - 2.494).abs() < 0.01, "kT(300) = {kt}");
    }

    #[test]
    fn coulomb_prefactor_is_gromacs_value() {
        // Two unit charges 1 nm apart: ~139 kJ/mol.
        assert!((COULOMB - 138.935).abs() < 0.01);
    }

    #[test]
    fn pressure_conversion_magnitude() {
        // Sanity: the factor is ~16.6.
        assert!((PRESSURE_KJMOLNM3_TO_BAR - 16.6054).abs() < 0.01);
    }
}
