//! MP2 — second-order Møller-Plesset correlation energy.
//!
//! Hartree-Fock leaves out *electron correlation* — the instantaneous
//! avoidance of electrons beyond the mean field. Second-order
//! Møller-Plesset perturbation theory (MP2) is the cheapest systematic
//! correction. From a converged closed-shell RHF reference its
//! correlation energy is
//!
//! ```text
//! E_MP2 = Σ_{ijab} (ia|jb) [ 2(ia|jb) − (ib|ja) ]
//!                  / (εᵢ + εⱼ − εₐ − ε_b)
//! ```
//!
//! where `i, j` run over occupied and `a, b` over virtual orbitals,
//! `(ia|jb)` are electron-repulsion integrals in the *molecular*-orbital
//! basis, and `ε` are the orbital energies.
//!
//! ## The AO→MO transformation
//!
//! The ERIs come out of the integral engine in the atomic-orbital
//! basis. [`ao_to_mo_transform`] contracts them with the orbital
//! coefficients,
//! `(pq|rs) = Σ_{μνλσ} C_{μp} C_{νq} C_{λr} C_{σs} (μν|λσ)`,
//! done as four successive quarter-transforms.
//!
//! ## v1 scope
//!
//! Closed-shell (RHF-reference) MP2 only — no UHF-reference (UMP2),
//! no frozen-core option, no spin-component scaling. The naive
//! quarter-transform is `O(n⁵)`; correct and adequate for the
//! small-molecule regime, not a production transformation.

use crate::error::{QchemError, Result};
use crate::integrals::two_electron::EriTensor;
use crate::scf::RhfResult;
use nalgebra::DMatrix;

/// The result of an MP2 calculation.
#[derive(Copy, Clone, Debug)]
pub struct Mp2Result {
    /// The MP2 correlation energy `E_MP2` (Hartree, always `≤ 0`).
    pub correlation_energy: f64,
    /// The reference RHF total energy (Hartree).
    pub reference_energy: f64,
    /// The same-spin contribution to the correlation energy.
    pub same_spin: f64,
    /// The opposite-spin contribution to the correlation energy.
    pub opposite_spin: f64,
}

impl Mp2Result {
    /// The total MP2 energy — the RHF reference plus the correlation
    /// correction.
    pub fn total_energy(&self) -> f64 {
        self.reference_energy + self.correlation_energy
    }
}

/// Transform the four-index ERI tensor from the atomic-orbital basis
/// into the molecular-orbital basis.
///
/// `coefficients` is the `n_basis × n_mo` matrix `C`. The transform is
/// the four quarter-contractions
/// `(pq|rs) = Σ C_{μp} C_{νq} C_{λr} C_{σs} (μν|λσ)`.
pub fn ao_to_mo_transform(ao: &EriTensor, coefficients: &DMatrix<f64>) -> EriTensor {
    let n = ao.n;
    let n_mo = coefficients.ncols();

    // Quarter 1: μ → p. step1[p][ν][λ][σ].
    let mut step1 = vec![0.0; n_mo * n * n * n];
    for p in 0..n_mo {
        for nu in 0..n {
            for la in 0..n {
                for si in 0..n {
                    let mut acc = 0.0;
                    for mu in 0..n {
                        acc += coefficients[(mu, p)] * ao.get(mu, nu, la, si);
                    }
                    step1[((p * n + nu) * n + la) * n + si] = acc;
                }
            }
        }
    }
    // Quarter 2: ν → q. step2[p][q][λ][σ].
    let mut step2 = vec![0.0; n_mo * n_mo * n * n];
    for p in 0..n_mo {
        for q in 0..n_mo {
            for la in 0..n {
                for si in 0..n {
                    let mut acc = 0.0;
                    for nu in 0..n {
                        acc += coefficients[(nu, q)]
                            * step1[((p * n + nu) * n + la) * n + si];
                    }
                    step2[((p * n_mo + q) * n + la) * n + si] = acc;
                }
            }
        }
    }
    // Quarter 3: λ → r. step3[p][q][r][σ].
    let mut step3 = vec![0.0; n_mo * n_mo * n_mo * n];
    for p in 0..n_mo {
        for q in 0..n_mo {
            for r in 0..n_mo {
                for si in 0..n {
                    let mut acc = 0.0;
                    for la in 0..n {
                        acc += coefficients[(la, r)]
                            * step2[((p * n_mo + q) * n + la) * n + si];
                    }
                    step3[((p * n_mo + q) * n_mo + r) * n + si] = acc;
                }
            }
        }
    }
    // Quarter 4: σ → s. mo[p][q][r][s].
    let mut mo = EriTensor::zeros(n_mo);
    for p in 0..n_mo {
        for q in 0..n_mo {
            for r in 0..n_mo {
                for s in 0..n_mo {
                    let mut acc = 0.0;
                    for si in 0..n {
                        acc += coefficients[(si, s)]
                            * step3[((p * n_mo + q) * n_mo + r) * n + si];
                    }
                    mo.set(p, q, r, s, acc);
                }
            }
        }
    }
    mo
}

/// Compute the MP2 correlation energy from a converged RHF reference.
///
/// `ao_eri` is the atomic-orbital ERI tensor used for the RHF
/// calculation.
///
/// # Errors
///
/// Returns [`QchemError::InvalidInput`] when the reference has no
/// virtual orbitals (MP2 has nothing to correlate into).
pub fn mp2_energy(rhf: &RhfResult, ao_eri: &EriTensor) -> Result<Mp2Result> {
    let n_occ = rhf.n_occupied;
    let n_mo = rhf.orbital_coefficients.ncols();
    let n_virt = n_mo - n_occ;
    if n_virt == 0 {
        return Err(QchemError::invalid(
            "MP2 needs virtual orbitals — the reference basis is fully \
             occupied",
        ));
    }

    let mo = ao_to_mo_transform(ao_eri, &rhf.orbital_coefficients);
    let eps = &rhf.orbital_energies;

    let mut same_spin = 0.0;
    let mut opposite_spin = 0.0;
    // i, j occupied; a, b virtual (offset by n_occ).
    for i in 0..n_occ {
        for j in 0..n_occ {
            for a in n_occ..n_mo {
                for b in n_occ..n_mo {
                    let iajb = mo.get(i, a, j, b);
                    let ibja = mo.get(i, b, j, a);
                    let denom = eps[i] + eps[j] - eps[a] - eps[b];
                    // Opposite-spin: (ia|jb)²/denom.
                    // Same-spin: (ia|jb)[(ia|jb) − (ib|ja)]/denom.
                    opposite_spin += iajb * iajb / denom;
                    same_spin += iajb * (iajb - ibja) / denom;
                }
            }
        }
    }
    // E_MP2 = Σ (ia|jb)[2(ia|jb) − (ib|ja)] / denom
    //       = opposite_spin + same_spin  (with the decomposition above).
    let correlation_energy = opposite_spin + same_spin;

    Ok(Mp2Result {
        correlation_energy,
        reference_energy: rhf.total_energy,
        same_spin,
        opposite_spin,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basis::BasisSet;
    use crate::geometry::{Atom, MolecularGeometry};
    use crate::integrals::IntegralSet;
    use crate::scf::rhf::{run_rhf_scf, ScfSettings};

    fn h2_rhf() -> (RhfResult, IntegralSet) {
        let geom = MolecularGeometry::new(vec![
            Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.0]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.7414]).unwrap(),
        ]);
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        let ints = IntegralSet::compute(&geom, &basis);
        let rhf = run_rhf_scf(&ints, 2, ScfSettings::default()).unwrap();
        (rhf, ints)
    }

    #[test]
    fn ao_to_mo_preserves_a_diagonal_transform() {
        // With C = identity the MO tensor equals the AO tensor.
        let geom = MolecularGeometry::new(vec![
            Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.0]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.7414]).unwrap(),
        ]);
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        let ao = EriTensor::build(&basis);
        let id = DMatrix::<f64>::identity(ao.n, ao.n);
        let mo = ao_to_mo_transform(&ao, &id);
        for p in 0..ao.n {
            for q in 0..ao.n {
                assert!((mo.get(p, q, p, q) - ao.get(p, q, p, q)).abs() < 1.0e-12);
            }
        }
    }

    #[test]
    fn mp2_correlation_is_negative() {
        // The MP2 correlation energy is always ≤ 0.
        let (rhf, ints) = h2_rhf();
        let mp2 = mp2_energy(&rhf, &ints.eri).unwrap();
        assert!(
            mp2.correlation_energy <= 0.0,
            "E_corr = {}",
            mp2.correlation_energy
        );
        // STO-3G H2 MP2 correlation is small but nonzero (~ -0.013 Ha).
        assert!(mp2.correlation_energy < -1.0e-4);
    }

    #[test]
    fn mp2_lowers_the_total_energy() {
        let (rhf, ints) = h2_rhf();
        let mp2 = mp2_energy(&rhf, &ints.eri).unwrap();
        assert!(mp2.total_energy() < rhf.total_energy);
        // The total = reference + correlation.
        assert!(
            (mp2.total_energy()
                - (mp2.reference_energy + mp2.correlation_energy))
                .abs()
                < 1.0e-14
        );
    }

    #[test]
    fn same_and_opposite_spin_sum_to_total() {
        let (rhf, ints) = h2_rhf();
        let mp2 = mp2_energy(&rhf, &ints.eri).unwrap();
        assert!(
            (mp2.same_spin + mp2.opposite_spin - mp2.correlation_energy).abs()
                < 1.0e-12
        );
    }

    #[test]
    fn mp2_water_correlation_is_substantial() {
        // STO-3G water has a clearly nonzero MP2 correlation energy.
        let geom = MolecularGeometry::new(vec![
            Atom::from_symbol_angstrom("O", [0.0, 0.0, 0.1173]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, 0.7572, -0.4692]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, -0.7572, -0.4692]).unwrap(),
        ]);
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        let ints = IntegralSet::compute(&geom, &basis);
        let rhf = run_rhf_scf(&ints, 10, ScfSettings::default()).unwrap();
        let mp2 = mp2_energy(&rhf, &ints.eri).unwrap();
        // STO-3G water MP2 correlation is roughly -0.03 to -0.05 Ha.
        assert!(
            mp2.correlation_energy < -0.01 && mp2.correlation_energy > -0.2,
            "water E_corr = {}",
            mp2.correlation_energy
        );
    }
}
