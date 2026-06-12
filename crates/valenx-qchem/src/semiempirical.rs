//! Extended-Hückel theory — a fast, no-SCF molecular-orbital preview.
//!
//! Extended Hückel (Hoffmann 1963) is the simplest all-valence-electron
//! molecular-orbital method. It needs no two-electron integrals and no
//! self-consistent iteration, so it is *much* cheaper than Hartree-Fock
//! — useful for a quick qualitative look at a molecule's orbitals
//! before committing to a real SCF.
//!
//! ## The method
//!
//! A minimal valence basis of Slater-type orbitals is placed on each
//! atom (1s for hydrogen, 2s + 2p for the first row). The single
//! one-electron Hamiltonian matrix is built directly:
//!
//! ```text
//! H_μμ = VSIP_μ                                 (diagonal — ionisation potentials)
//! H_μν = K · S_μν · (H_μμ + H_νν) / 2           (off-diagonal — Wolfsberg-Helmholz)
//! ```
//!
//! with `K = 1.75` the standard Wolfsberg-Helmholz constant and `S` the
//! overlap matrix over the valence Slater orbitals. Diagonalising the
//! single generalized eigenproblem `H C = S C ε` gives the orbitals.
//!
//! ## v1 scope
//!
//! Valence STO basis as a single Gaussian fit per orbital (an STO-1G
//! contraction); the standard VSIP parameter table for H–Ne; no charge
//! self-consistency (this is plain extended Hückel, not the iterated
//! EHT). It is a *preview* tool — the energies are qualitative.

use crate::basis::{AngularMomentum, BasisFunction, BasisSet, Primitive, Shell};
use crate::error::{QchemError, Result};
use crate::geometry::MolecularGeometry;
use crate::integrals::one_electron::overlap;
use crate::scf::linalg::{solve_roothaan, symmetric_orthogonalizer};
use nalgebra::{DMatrix, DVector};

/// The Wolfsberg-Helmholz proportionality constant for the off-diagonal
/// Hamiltonian elements.
pub const WOLFSBERG_HELMHOLZ_K: f64 = 1.75;

/// Valence-shell ionisation potentials in Hartree for H–Ne, used as the
/// diagonal Hamiltonian elements `H_μμ`.
///
/// Each entry is `(Z, [VSIP_s, VSIP_p])` — the `s` and `p` valence-shell
/// ionisation potentials. Values are the standard Hoffmann / Anderson
/// parameters converted from eV (1 eV = 0.0367493 Ha).
struct VsipEntry {
    /// `s`-orbital valence ionisation potential (Hartree, negative).
    s: f64,
    /// `p`-orbital valence ionisation potential (Hartree, negative).
    p: f64,
}

/// Look up the valence ionisation potentials for element `z`.
fn vsip(z: u8) -> Option<VsipEntry> {
    // eV → Hartree.
    const EV: f64 = 0.036_749_322;
    let (s_ev, p_ev) = match z {
        1 => (-13.60, 0.0),
        2 => (-23.40, 0.0),
        3 => (-5.40, -3.50),
        4 => (-10.00, -6.00),
        5 => (-15.20, -8.50),
        6 => (-21.40, -11.40),
        7 => (-26.00, -13.40),
        8 => (-32.30, -14.80),
        9 => (-40.00, -18.10),
        10 => (-43.20, -20.00),
        _ => return None,
    };
    Some(VsipEntry {
        s: s_ev * EV,
        p: p_ev * EV,
    })
}

/// Single-Gaussian (STO-1G) exponents approximating a valence Slater
/// orbital. `(Z, [α_s, α_p])`. The values are STO-1G fits of the
/// valence shell.
fn valence_exponents(z: u8) -> Option<(f64, f64)> {
    // STO-1G valence exponents (least-squares Gaussian fit of the
    // Slater valence orbital).
    Some(match z {
        1 => (0.282_94, 0.0),
        2 => (0.480_00, 0.0),
        3 => (0.082_18, 0.090_00),
        4 => (0.150_50, 0.160_00),
        5 => (0.231_00, 0.230_00),
        6 => (0.314_70, 0.310_00),
        7 => (0.412_90, 0.400_00),
        8 => (0.524_30, 0.500_00),
        9 => (0.648_00, 0.620_00),
        10 => (0.785_00, 0.750_00),
        _ => return None,
    })
}

/// The result of an extended-Hückel calculation.
#[derive(Clone, Debug)]
pub struct HuckelResult {
    /// Valence molecular-orbital energies, ascending (Hartree).
    pub orbital_energies: DVector<f64>,
    /// Valence molecular-orbital coefficients.
    pub orbital_coefficients: DMatrix<f64>,
    /// The sum of the occupied orbital energies — the extended-Hückel
    /// "total energy" (Hartree). A qualitative quantity only.
    pub orbital_energy_sum: f64,
    /// Number of doubly-occupied valence orbitals.
    pub n_occupied: usize,
}

impl HuckelResult {
    /// The HOMO energy (`None` for a system with no valence electrons).
    pub fn homo_energy(&self) -> Option<f64> {
        if self.n_occupied == 0 {
            None
        } else {
            Some(self.orbital_energies[self.n_occupied - 1])
        }
    }

    /// The LUMO energy (`None` when every valence orbital is occupied).
    pub fn lumo_energy(&self) -> Option<f64> {
        if self.n_occupied >= self.orbital_energies.len() {
            None
        } else {
            Some(self.orbital_energies[self.n_occupied])
        }
    }

    /// The HOMO-LUMO gap in Hartree (`None` when undefined).
    pub fn homo_lumo_gap(&self) -> Option<f64> {
        Some(self.lumo_energy()? - self.homo_energy()?)
    }
}

/// Build the minimal valence basis for the extended-Hückel calculation:
/// a 1s shell on hydrogen / helium, a 2s + 2p pair on the first row,
/// each shell a single normalised Gaussian.
fn valence_basis(geometry: &MolecularGeometry) -> Result<BasisSet> {
    let mut shells = Vec::new();
    for (atom_index, atom) in geometry.atoms.iter().enumerate() {
        let z = atom.element.atomic_number();
        let (alpha_s, alpha_p) = valence_exponents(z).ok_or_else(|| {
            QchemError::invalid(format!(
                "extended Hückel has no valence parameters for {}",
                atom.element.symbol()
            ))
        })?;
        // Valence s shell.
        shells.push(
            Shell {
                atom_index,
                centre: atom.position,
                angular: AngularMomentum::S,
                primitives: vec![Primitive {
                    exponent: alpha_s,
                    coefficient: 1.0,
                }],
            }
            .normalised(),
        );
        // Valence p shell only for the first row (Z >= 3).
        if z >= 3 {
            shells.push(
                Shell {
                    atom_index,
                    centre: atom.position,
                    angular: AngularMomentum::P,
                    primitives: vec![Primitive {
                        exponent: alpha_p,
                        coefficient: 1.0,
                    }],
                }
                .normalised(),
            );
        }
    }
    // Expand into basis functions.
    let mut functions = Vec::new();
    for shell in &shells {
        for cart in shell.angular.cartesian_components() {
            functions.push(BasisFunction {
                atom_index: shell.atom_index,
                centre: shell.centre,
                cart,
                primitives: shell.primitives.clone(),
            });
        }
    }
    Ok(BasisSet {
        name: "eht-valence",
        shells,
        functions,
    })
}

/// The valence-electron count — total electrons minus the core
/// electrons (2 per first-row atom, 0 for H / He).
fn valence_electron_count(geometry: &MolecularGeometry) -> Result<u32> {
    let total = geometry.n_electrons()?;
    let core: u32 = geometry
        .atoms
        .iter()
        .map(|a| if a.element.atomic_number() >= 3 { 2 } else { 0 })
        .sum();
    Ok(total.saturating_sub(core))
}

/// The diagonal Hamiltonian element for a basis function — its valence
/// ionisation potential, chosen by angular momentum.
fn diagonal_h(f: &BasisFunction, z: u8) -> Result<f64> {
    let entry = vsip(z)
        .ok_or_else(|| QchemError::invalid(format!("extended Hückel has no VSIP for Z={z}")))?;
    Ok(if f.l() == 0 { entry.s } else { entry.p })
}

/// Run an extended-Hückel calculation on a molecule.
///
/// Builds the minimal valence basis, the overlap matrix and the
/// Wolfsberg-Helmholz Hamiltonian, and solves the single generalized
/// eigenproblem `H C = S C ε`.
///
/// # Errors
///
/// Returns [`QchemError::InvalidInput`] when an element lacks
/// extended-Hückel parameters or the molecule is empty.
pub fn run_extended_huckel(geometry: &MolecularGeometry) -> Result<HuckelResult> {
    let basis = valence_basis(geometry)?;
    let n = basis.n_functions();
    if n == 0 {
        return Err(QchemError::invalid("molecule has no valence orbitals"));
    }

    // Overlap and diagonal Hamiltonian.
    let mut s = DMatrix::<f64>::zeros(n, n);
    let mut h_diag = vec![0.0; n];
    for mu in 0..n {
        let z = geometry.atoms[basis.functions[mu].atom_index]
            .element
            .atomic_number();
        h_diag[mu] = diagonal_h(&basis.functions[mu], z)?;
        for nu in 0..=mu {
            let s_v = overlap(&basis.functions[mu], &basis.functions[nu]);
            s[(mu, nu)] = s_v;
            s[(nu, mu)] = s_v;
        }
    }

    // Wolfsberg-Helmholz off-diagonal Hamiltonian.
    let mut h = DMatrix::<f64>::zeros(n, n);
    for mu in 0..n {
        h[(mu, mu)] = h_diag[mu];
        for nu in 0..mu {
            let h_off = WOLFSBERG_HELMHOLZ_K * s[(mu, nu)] * 0.5 * (h_diag[mu] + h_diag[nu]);
            h[(mu, nu)] = h_off;
            h[(nu, mu)] = h_off;
        }
    }

    // Solve H C = S C ε.
    let ortho = symmetric_orthogonalizer(&s)?;
    let (eps, c) = solve_roothaan(&h, &ortho);

    let n_valence = valence_electron_count(geometry)? as usize;
    let n_occupied = n_valence / 2;
    let orbital_energy_sum: f64 = (0..n_occupied.min(eps.len())).map(|i| 2.0 * eps[i]).sum();

    Ok(HuckelResult {
        orbital_energies: eps,
        orbital_coefficients: c,
        orbital_energy_sum,
        n_occupied,
    })
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

    #[test]
    fn h2_has_two_valence_orbitals() {
        let res = run_extended_huckel(&h2()).unwrap();
        assert_eq!(res.orbital_energies.len(), 2);
        assert_eq!(res.n_occupied, 1);
    }

    #[test]
    fn h2_bonding_orbital_below_antibonding() {
        let res = run_extended_huckel(&h2()).unwrap();
        // The occupied bonding MO must lie below the virtual one.
        assert!(res.orbital_energies[0] < res.orbital_energies[1]);
        assert!(res.homo_lumo_gap().unwrap() > 0.0);
    }

    #[test]
    fn h2_bonding_orbital_is_stabilised() {
        // The bonding MO must drop below the bare hydrogen 1s VSIP.
        let res = run_extended_huckel(&h2()).unwrap();
        let h_1s = vsip(1).unwrap().s;
        assert!(
            res.orbital_energies[0] < h_1s,
            "bonding MO {} vs H 1s {}",
            res.orbital_energies[0],
            h_1s
        );
    }

    #[test]
    fn water_valence_basis_size() {
        // Water valence basis: O (2s+2p = 4) + 2×H (1s) = 6 functions.
        let geom = MolecularGeometry::new(vec![
            Atom::from_symbol_angstrom("O", [0.0, 0.0, 0.1173]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, 0.7572, -0.4692]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, -0.7572, -0.4692]).unwrap(),
        ]);
        let res = run_extended_huckel(&geom).unwrap();
        assert_eq!(res.orbital_energies.len(), 6);
        // Water has 8 valence electrons → 4 occupied valence MOs.
        assert_eq!(res.n_occupied, 4);
    }

    #[test]
    fn methane_orbital_ordering() {
        // CH4 — a quick qualitative check that EHT produces a sensible
        // bound set with a positive HOMO-LUMO gap.
        let d = 0.629; // C-H, ångström along the cube diagonals.
        let geom = MolecularGeometry::new(vec![
            Atom::from_symbol_angstrom("C", [0.0, 0.0, 0.0]).unwrap(),
            Atom::from_symbol_angstrom("H", [d, d, d]).unwrap(),
            Atom::from_symbol_angstrom("H", [-d, -d, d]).unwrap(),
            Atom::from_symbol_angstrom("H", [-d, d, -d]).unwrap(),
            Atom::from_symbol_angstrom("H", [d, -d, -d]).unwrap(),
        ]);
        let res = run_extended_huckel(&geom).unwrap();
        // C (2s+2p) + 4 H = 8 valence functions; 8 valence electrons.
        assert_eq!(res.orbital_energies.len(), 8);
        assert_eq!(res.n_occupied, 4);
        assert!(res.homo_lumo_gap().unwrap() > 0.0);
        // Energies are sorted ascending.
        for w in res.orbital_energies.as_slice().windows(2) {
            assert!(w[0] <= w[1] + 1.0e-12);
        }
    }
}
