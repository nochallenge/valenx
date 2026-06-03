//! Gasteiger-Marsili partial atomic charges (the PEOE method).
//!
//! Partial Equalization of Orbital Electronegativity assigns a partial
//! charge to every atom by an iterative scheme:
//!
//! 1. each atom type has three coefficients `a, b, c` defining its
//!    electronegativity as a quadratic in its current charge `q`:
//!    `χ(q) = a + b·q + c·q²`;
//! 2. across every bond, charge flows from the less to the more
//!    electronegative atom, scaled by the difference and by a damping
//!    factor `0.5^k` that shrinks each iteration `k`;
//! 3. after ~6 iterations the orbital electronegativities are
//!    partially equalised and the charges have converged.
//!
//! The `(a, b, c)` parameters here are the published Gasteiger-Marsili
//! values for the common hybridisation states of H, C, N, O, F, Cl,
//! Br, I, S and P. Atoms outside the table get a neutral fallback
//! (no charge flow) — the documented v1 limitation.

use crate::element;
use crate::molecule::{BondOrder, Molecule};

/// PEOE electronegativity coefficients for one atom type.
#[derive(Copy, Clone, Debug)]
struct Peoe {
    a: f64,
    b: f64,
    c: f64,
}

impl Peoe {
    /// Orbital electronegativity at charge `q`: `a + b·q + c·q²`.
    fn chi(&self, q: f64) -> f64 {
        self.a + self.b * q + self.c * q * q
    }
    /// `χ⁺` — the electronegativity of the cation, used to normalise
    /// the charge transferred across a bond.
    fn chi_plus(&self) -> f64 {
        self.a + self.b + self.c
    }
}

/// Choose PEOE coefficients for an atom from its element and a simple
/// hybridisation guess (number of π bonds / aromatic membership).
fn peoe_params(mol: &Molecule, atom: usize) -> Peoe {
    let a = &mol.atoms[atom];
    let z = a.atomic_number;
    // count multiple bonds to classify sp3 / sp2 / sp
    let mut doubles = 0;
    let mut triples = 0;
    let mut aromatic = a.aromatic;
    for bi in mol.bonds_on(atom) {
        match mol.bonds[bi].order {
            BondOrder::Double => doubles += 1,
            BondOrder::Triple => triples += 1,
            BondOrder::Aromatic => aromatic = true,
            _ => {}
        }
    }
    // Gasteiger-Marsili (1980) coefficients, in volts.
    match z {
        1 => Peoe {
            a: 7.17,
            b: 6.24,
            c: -0.56,
        },
        6 => {
            if triples > 0 {
                Peoe {
                    a: 10.39,
                    b: 9.45,
                    c: 0.73,
                } // sp
            } else if doubles > 0 || aromatic {
                Peoe {
                    a: 8.79,
                    b: 9.32,
                    c: 1.51,
                } // sp2
            } else {
                Peoe {
                    a: 7.98,
                    b: 9.18,
                    c: 1.88,
                } // sp3
            }
        }
        7 => {
            if triples > 0 {
                Peoe {
                    a: 15.68,
                    b: 11.70,
                    c: -0.27,
                }
            } else if doubles > 0 || aromatic {
                Peoe {
                    a: 12.87,
                    b: 11.15,
                    c: 0.85,
                }
            } else {
                Peoe {
                    a: 11.54,
                    b: 10.82,
                    c: 1.36,
                }
            }
        }
        8 => {
            if doubles > 0 || aromatic {
                Peoe {
                    a: 17.07,
                    b: 13.79,
                    c: 0.47,
                }
            } else {
                Peoe {
                    a: 14.18,
                    b: 12.92,
                    c: 1.39,
                }
            }
        }
        9 => Peoe {
            a: 14.66,
            b: 13.85,
            c: 2.31,
        },
        15 => Peoe {
            a: 8.90,
            b: 8.32,
            c: 1.58,
        },
        16 => Peoe {
            a: 10.14,
            b: 9.13,
            c: 1.38,
        },
        17 => Peoe {
            a: 11.00,
            b: 9.69,
            c: 1.35,
        },
        35 => Peoe {
            a: 10.08,
            b: 8.47,
            c: 1.16,
        },
        53 => Peoe {
            a: 9.90,
            b: 7.96,
            c: 0.96,
        },
        _ => {
            // fallback from raw Pauling electronegativity — yields no
            // charge flow against tabulated atoms beyond a tiny bias
            let chi = element::electronegativity(z) * 4.0;
            Peoe {
                a: chi,
                b: chi,
                c: 0.0,
            }
        }
    }
}

/// Compute Gasteiger-Marsili partial charges for every atom, returned
/// parallel to [`Molecule::atoms`]. Formal charges seed the iteration
/// (`q₀ = formal_charge`); hydrogens attached implicitly are *not*
/// charged separately — only atoms present as nodes get a value.
pub fn gasteiger_charges(mol: &Molecule) -> Vec<f64> {
    let n = mol.atoms.len();
    let mut q: Vec<f64> = mol
        .atoms
        .iter()
        .map(|a| f64::from(a.formal_charge))
        .collect();
    if n == 0 {
        return q;
    }
    let params: Vec<Peoe> = (0..n).map(|i| peoe_params(mol, i)).collect();

    // 6 PEOE iterations with geometric damping
    let iterations = 6;
    let mut damping = 1.0;
    for _ in 0..iterations {
        damping *= 0.5;
        let chi: Vec<f64> = (0..n).map(|i| params[i].chi(q[i])).collect();
        let mut dq = vec![0.0f64; n];
        for bond in &mol.bonds {
            let (i, j) = (bond.a, bond.b);
            if i >= n || j >= n {
                continue;
            }
            // Electron density flows toward the higher electronegativity.
            // `acceptor` is the more electronegative atom; it GAINS
            // electron density, so its partial charge `q` becomes more
            // NEGATIVE (q is positive for an electron-deficient atom).
            // The donor loses density and becomes more positive.
            let (donor, acceptor, diff) = if chi[i] > chi[j] {
                (j, i, chi[i] - chi[j])
            } else {
                (i, j, chi[j] - chi[i])
            };
            // normalise by the cation electronegativity of the donor
            let scale = params[donor].chi_plus().max(1.0);
            let transfer = (diff / scale) * damping;
            // The acceptor's charge goes negative, the donor's positive.
            dq[acceptor] -= transfer;
            dq[donor] += transfer;
        }
        for i in 0..n {
            q[i] += dq[i];
        }
    }
    q
}

/// Total Gasteiger charge — should equal the molecule's net formal
/// charge (the method conserves charge). Exposed as a self-check.
pub fn total_gasteiger_charge(mol: &Molecule) -> f64 {
    gasteiger_charges(mol).iter().sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol_from_smiles;

    #[test]
    fn charge_is_conserved() {
        let m = mol_from_smiles("CCO").unwrap();
        let total: f64 = gasteiger_charges(&m).iter().sum();
        assert!(total.abs() < 1e-6, "neutral molecule total = {total}");
    }

    #[test]
    fn electronegative_atom_is_negative() {
        // in ethanol the oxygen must carry a negative partial charge
        let m = mol_from_smiles("CCO").unwrap();
        let q = gasteiger_charges(&m);
        let o = m.atoms.iter().position(|a| a.atomic_number == 8).unwrap();
        assert!(q[o] < 0.0, "oxygen charge = {}", q[o]);
        // the carbon bonded to oxygen is pulled positive
        let c_next = m.neighbors(o)[0];
        assert!(q[c_next] > 0.0, "alpha carbon = {}", q[c_next]);
    }

    #[test]
    fn carbonyl_carbon_positive() {
        let m = mol_from_smiles("CC(=O)C").unwrap(); // acetone
        let q = gasteiger_charges(&m);
        let carbonyl_c = m
            .atoms
            .iter()
            .enumerate()
            .find(|(i, a)| {
                a.atomic_number == 6
                    && m.bonds_on(*i)
                        .iter()
                        .any(|&b| m.bonds[b].order == BondOrder::Double)
            })
            .map(|(i, _)| i)
            .unwrap();
        assert!(q[carbonyl_c] > 0.0);
    }

    #[test]
    fn cation_total_matches_formal() {
        let m = mol_from_smiles("[NH4+]").unwrap();
        let total = total_gasteiger_charge(&m);
        assert!((total - 1.0).abs() < 1e-6, "ammonium total = {total}");
    }

    #[test]
    fn empty_molecule_has_no_charges() {
        // The `n == 0` early-return path in `gasteiger_charges`.
        let m = Molecule::new();
        assert!(gasteiger_charges(&m).is_empty());
        assert_eq!(total_gasteiger_charge(&m), 0.0);
    }

    #[test]
    fn alkyne_triple_bond_carbon_path() {
        // Propyne CC#C exercises the sp-carbon (`triples > 0`) PEOE
        // arm and the sp-nitrogen-style triple-bond branch detection.
        // Charge must still be conserved for the neutral molecule.
        let m = mol_from_smiles("CC#C").unwrap();
        let q = gasteiger_charges(&m);
        let total: f64 = q.iter().sum();
        assert!(total.abs() < 1e-6, "propyne total = {total}");
        assert_eq!(q.len(), m.atoms.len());
    }

    #[test]
    fn nitrile_triple_bond_nitrogen_path() {
        // Acetonitrile CC#N exercises the triple-bond nitrogen PEOE
        // arm. The nitrile nitrogen (electronegative) is negative.
        let m = mol_from_smiles("CC#N").unwrap();
        let q = gasteiger_charges(&m);
        let n = m.atoms.iter().position(|a| a.atomic_number == 7).unwrap();
        assert!(q[n] < 0.0, "nitrile N charge = {}", q[n]);
        let total: f64 = q.iter().sum();
        assert!(total.abs() < 1e-6, "acetonitrile total = {total}");
    }

    #[test]
    fn amine_single_bond_nitrogen_path() {
        // Methylamine CN — the sp3 single-bond nitrogen arm (no
        // doubles, no triples, not aromatic).
        let m = mol_from_smiles("CN").unwrap();
        let q = gasteiger_charges(&m);
        let n = m.atoms.iter().position(|a| a.atomic_number == 7).unwrap();
        assert!(q[n] < 0.0, "amine N pulls electron density: {}", q[n]);
    }

    #[test]
    fn ether_oxygen_single_bond_path() {
        // Dimethyl ether COC — the single-bond oxygen arm (no double,
        // not aromatic). The ether O is negative.
        let m = mol_from_smiles("COC").unwrap();
        let q = gasteiger_charges(&m);
        let o = m.atoms.iter().position(|a| a.atomic_number == 8).unwrap();
        assert!(q[o] < 0.0, "ether O charge = {}", q[o]);
    }

    #[test]
    fn halogen_phosphorus_sulfur_peoe_arms() {
        // One molecule per remaining tabulated PEOE element so each
        // explicit `match z` arm is exercised. Each is neutral, so
        // the total Gasteiger charge must vanish.
        for (smiles, z) in [
            ("CF", 9u8),    // fluoromethane — F arm
            ("CCl", 17),    // chloromethane — Cl arm
            ("CBr", 35),    // bromomethane — Br arm
            ("CI", 53),     // iodomethane — I arm
            ("CS", 16),     // methanethiol — S arm
            ("CP", 15),     // methylphosphine — P arm
        ] {
            let m = mol_from_smiles(smiles).unwrap();
            assert!(
                m.atoms.iter().any(|a| a.atomic_number == z),
                "{smiles} should contain Z={z}",
            );
            let total: f64 = gasteiger_charges(&m).iter().sum();
            assert!(
                total.abs() < 1e-6,
                "{smiles} charge not conserved: total = {total}",
            );
            // The halogen / heteroatom is more electronegative than
            // carbon, so it carries a non-positive partial charge.
            let q = gasteiger_charges(&m);
            let hetero = m
                .atoms
                .iter()
                .position(|a| a.atomic_number == z)
                .unwrap();
            assert!(
                q[hetero] <= 1e-9,
                "{smiles}: Z={z} should not be positive, got {}",
                q[hetero],
            );
        }
    }

    #[test]
    fn untabulated_element_uses_fallback_arm() {
        // Boron is not in the PEOE table, so `peoe_params` takes the
        // `_` fallback arm (electronegativity-seeded, no charge flow).
        // Trimethylborane B(C)(C)C stays charge-conserved.
        let m = mol_from_smiles("B(C)(C)C").unwrap();
        assert!(m.atoms.iter().any(|a| a.atomic_number == 5));
        let total: f64 = gasteiger_charges(&m).iter().sum();
        assert!(total.abs() < 1e-6, "borane total = {total}");
    }
}
