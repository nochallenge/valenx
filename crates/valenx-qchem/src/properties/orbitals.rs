//! Molecular-orbital energies, occupations and the HOMO-LUMO gap.
//!
//! After the SCF converges, the orbital energies `ε` and the occupation
//! pattern describe the molecule's electronic structure. This module
//! packages them into an [`OrbitalSummary`]: the energy and occupation
//! of every molecular orbital, the HOMO and LUMO indices, and the
//! frontier-orbital gap.
//!
//! Energies are reported in Hartree and, via [`HARTREE_TO_EV`], in
//! electron-volts (1 Ha = 27.211386 eV).

use nalgebra::DVector;

/// Electron-volts per Hartree (CODATA 2018).
pub const HARTREE_TO_EV: f64 = 27.211_386_245_988;

/// One molecular orbital — its energy and occupation.
#[derive(Copy, Clone, Debug)]
pub struct MolecularOrbital {
    /// 0-based orbital index, ascending in energy.
    pub index: usize,
    /// Orbital energy `ε` in Hartree.
    pub energy: f64,
    /// Occupation number — `2.0` for a doubly-occupied RHF orbital,
    /// `1.0` for a singly-occupied spin orbital, `0.0` for a virtual.
    pub occupation: f64,
}

impl MolecularOrbital {
    /// `true` when the orbital carries any electrons.
    pub fn is_occupied(&self) -> bool {
        self.occupation > 0.0
    }

    /// Orbital energy in electron-volts.
    pub fn energy_ev(&self) -> f64 {
        self.energy * HARTREE_TO_EV
    }
}

/// A summary of a molecule's molecular orbitals.
#[derive(Clone, Debug)]
pub struct OrbitalSummary {
    /// Every molecular orbital, ascending in energy.
    pub orbitals: Vec<MolecularOrbital>,
    /// Index of the highest occupied molecular orbital (`None` for a
    /// system with no electrons).
    pub homo_index: Option<usize>,
    /// Index of the lowest unoccupied molecular orbital (`None` when
    /// every orbital is occupied).
    pub lumo_index: Option<usize>,
}

impl OrbitalSummary {
    /// The HOMO energy in Hartree.
    pub fn homo_energy(&self) -> Option<f64> {
        self.homo_index.map(|i| self.orbitals[i].energy)
    }

    /// The LUMO energy in Hartree.
    pub fn lumo_energy(&self) -> Option<f64> {
        self.lumo_index.map(|i| self.orbitals[i].energy)
    }

    /// The HOMO-LUMO gap in Hartree.
    pub fn homo_lumo_gap(&self) -> Option<f64> {
        Some(self.lumo_energy()? - self.homo_energy()?)
    }

    /// The HOMO-LUMO gap in electron-volts.
    pub fn homo_lumo_gap_ev(&self) -> Option<f64> {
        Some(self.homo_lumo_gap()? * HARTREE_TO_EV)
    }
}

/// Build the orbital summary for a restricted (closed-shell)
/// calculation. The lowest `n_occupied` orbitals are doubly occupied.
pub fn restricted_summary(
    orbital_energies: &DVector<f64>,
    n_occupied: usize,
) -> OrbitalSummary {
    let orbitals: Vec<MolecularOrbital> = orbital_energies
        .iter()
        .enumerate()
        .map(|(i, &e)| MolecularOrbital {
            index: i,
            energy: e,
            occupation: if i < n_occupied { 2.0 } else { 0.0 },
        })
        .collect();
    let homo_index = if n_occupied == 0 {
        None
    } else {
        Some(n_occupied - 1)
    };
    let lumo_index = if n_occupied < orbitals.len() {
        Some(n_occupied)
    } else {
        None
    };
    OrbitalSummary {
        orbitals,
        homo_index,
        lumo_index,
    }
}

/// Build the orbital summary for one spin of an unrestricted
/// calculation. The lowest `n_spin` orbitals are singly occupied.
pub fn unrestricted_spin_summary(
    orbital_energies: &DVector<f64>,
    n_spin: usize,
) -> OrbitalSummary {
    let orbitals: Vec<MolecularOrbital> = orbital_energies
        .iter()
        .enumerate()
        .map(|(i, &e)| MolecularOrbital {
            index: i,
            energy: e,
            occupation: if i < n_spin { 1.0 } else { 0.0 },
        })
        .collect();
    let homo_index = if n_spin == 0 { None } else { Some(n_spin - 1) };
    let lumo_index = if n_spin < orbitals.len() {
        Some(n_spin)
    } else {
        None
    };
    OrbitalSummary {
        orbitals,
        homo_index,
        lumo_index,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn restricted_occupations_and_frontier() {
        let eps = DVector::from_vec(vec![-1.0, -0.5, 0.3, 0.8]);
        let summary = restricted_summary(&eps, 2);
        assert_eq!(summary.orbitals[0].occupation, 2.0);
        assert_eq!(summary.orbitals[2].occupation, 0.0);
        assert_eq!(summary.homo_index, Some(1));
        assert_eq!(summary.lumo_index, Some(2));
        assert!((summary.homo_lumo_gap().unwrap() - 0.8).abs() < 1.0e-12);
    }

    #[test]
    fn ev_conversion() {
        let eps = DVector::from_vec(vec![-1.0, 1.0]);
        let summary = restricted_summary(&eps, 1);
        // Gap is 2 Ha = ~54.4 eV.
        let gap = summary.homo_lumo_gap_ev().unwrap();
        assert!((gap - 2.0 * HARTREE_TO_EV).abs() < 1.0e-9);
    }

    #[test]
    fn fully_occupied_has_no_lumo() {
        let eps = DVector::from_vec(vec![-1.0, -0.5]);
        let summary = restricted_summary(&eps, 2);
        assert_eq!(summary.lumo_index, None);
        assert!(summary.homo_lumo_gap().is_none());
    }

    #[test]
    fn empty_system_has_no_homo() {
        let eps = DVector::from_vec(vec![0.5, 1.0]);
        let summary = restricted_summary(&eps, 0);
        assert_eq!(summary.homo_index, None);
    }

    #[test]
    fn unrestricted_spin_is_singly_occupied() {
        let eps = DVector::from_vec(vec![-1.0, -0.4, 0.5]);
        let summary = unrestricted_spin_summary(&eps, 2);
        assert_eq!(summary.orbitals[0].occupation, 1.0);
        assert_eq!(summary.orbitals[1].occupation, 1.0);
        assert_eq!(summary.orbitals[2].occupation, 0.0);
        assert_eq!(summary.homo_index, Some(1));
    }
}
