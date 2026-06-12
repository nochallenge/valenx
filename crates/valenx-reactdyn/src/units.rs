//! Unit conversions between user-facing units and the atomic units the
//! integrator runs in.
//!
//! The whole AIMD loop runs in **atomic units** (a.u.) — the units
//! `valenx-qchem` already uses internally: length in bohr, energy in
//! hartree, and therefore force in hartree/bohr. Mass is in electron
//! masses and time in atomic time units. Only two quantities cross the
//! boundary from human units — the timestep (femtoseconds) and atomic
//! masses (amu) — so this module stays small.

/// Bohr radii per ångström (CODATA 2018). Multiply ångström by this to
/// get bohr. Mirrors `valenx_qchem::geometry::BOHR_PER_ANGSTROM`.
pub const BOHR_PER_ANGSTROM: f64 = 1.889_726_124_625_770_2;

/// Electron masses per unified atomic mass unit (amu). In atomic units
/// `mₑ = 1`, and 1 amu = 1822.888486 mₑ (CODATA: mₑ = 5.485799×10⁻⁴ u).
/// Multiply an amu mass by this to get the mass the integrator uses.
pub const AMU_TO_AU_MASS: f64 = 1_822.888_486;

/// Atomic time units per femtosecond. One a.u. of time = 0.024188843 fs,
/// so 1 fs = 41.341375 a.u. Multiply a femtosecond timestep by this to
/// get the integrator's `dt`.
pub const FS_TO_AU_TIME: f64 = 41.341_374_575;

/// Convert a timestep in femtoseconds to atomic time units.
pub fn fs_to_au(fs: f64) -> f64 {
    fs * FS_TO_AU_TIME
}

/// Convert a mass in amu to electron masses (atomic units of mass).
pub fn amu_to_au_mass(amu: f64) -> f64 {
    amu * AMU_TO_AU_MASS
}

// --- GROMACS (valenx-md) <-> atomic-unit conversions (for QM/MM) ------
// valenx-md works in the GROMACS unit system (nm, kJ·mol⁻¹, u, e); the
// QM/MM integrator runs in atomic units, so these reconcile the MM side.
// Charge needs no conversion: e == 1 in atomic units. Mass reuses
// `AMU_TO_AU_MASS` (u and amu are the same unit).

/// Bohr per nanometre. (1 nm = 10 Å = 18.897261 bohr.)
pub const BOHR_PER_NM: f64 = 18.897_261_246_257_7;

/// Hartree per kJ·mol⁻¹. (1 hartree = 2625.499639 kJ·mol⁻¹.)
pub const HARTREE_PER_KJ_MOL: f64 = 1.0 / 2_625.499_639_479_9;

/// Atomic-unit force (hartree·bohr⁻¹) per GROMACS force (kJ·mol⁻¹·nm⁻¹).
/// Force is energy / length, so this is `HARTREE_PER_KJ_MOL / BOHR_PER_NM`.
pub const FORCE_AU_PER_KJ_MOL_NM: f64 = HARTREE_PER_KJ_MOL / BOHR_PER_NM;

/// Convert a length from nanometres to bohr.
pub fn nm_to_bohr(nm: f64) -> f64 {
    nm * BOHR_PER_NM
}

/// Convert an energy from kJ·mol⁻¹ to hartree.
pub fn kj_mol_to_hartree(e: f64) -> f64 {
    e * HARTREE_PER_KJ_MOL
}

/// Convert a force from kJ·mol⁻¹·nm⁻¹ to hartree·bohr⁻¹.
pub fn force_kj_mol_nm_to_au(f: f64) -> f64 {
    f * FORCE_AU_PER_KJ_MOL_NM
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn one_femtosecond_in_atomic_time() {
        assert!((fs_to_au(1.0) - 41.341_374_575).abs() < 1e-6);
    }

    #[test]
    fn carbon_mass_in_atomic_units() {
        // 12.011 amu × 1822.888486 ≈ 21894.7 electron masses.
        let m = amu_to_au_mass(12.011);
        assert!((m - 21894.7).abs() < 1.0, "got {m}");
    }

    #[test]
    fn angstrom_to_bohr_matches_codata() {
        assert!((BOHR_PER_ANGSTROM - 1.889_726_124_625_770_2).abs() < 1e-12);
    }

    #[test]
    fn one_nanometre_is_ten_angstrom_in_bohr() {
        // 1 nm = 10 Å = 18.897261 bohr; 0.1 nm = 1 Å.
        assert!((nm_to_bohr(1.0) - 18.897_261_246_257_7).abs() < 1e-9);
        assert!((nm_to_bohr(0.1) - BOHR_PER_ANGSTROM).abs() < 1e-9);
    }

    #[test]
    fn one_hartree_is_2625_kj_per_mol() {
        assert!((kj_mol_to_hartree(2_625.499_639_479_9) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn force_factor_is_energy_over_length() {
        assert!((force_kj_mol_nm_to_au(1.0) - HARTREE_PER_KJ_MOL / BOHR_PER_NM).abs() < 1e-15);
    }
}
