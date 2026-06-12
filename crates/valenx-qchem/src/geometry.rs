//! Molecular geometry — atoms, coordinates, charge and multiplicity.
//!
//! [`MolecularGeometry`] is the entry point for every calculation: a
//! list of [`Atom`]s (an [`Element`] plus a Cartesian position), an
//! overall molecular charge and a spin multiplicity `2S + 1`.
//!
//! ## Units
//!
//! All integrals and SCF code work in **atomic units** (bohr). The xyz
//! file format conventionally stores ångström, so [`MolecularGeometry`]
//! keeps its coordinates in bohr internally and the xyz reader / writer
//! convert. [`BOHR_PER_ANGSTROM`] is the conversion factor (CODATA
//! 2018: 1 Å = 1.8897261246257702 a₀).
//!
//! ## The xyz format
//!
//! ```text
//! 3
//! water, charge 0 mult 1
//! O   0.000000   0.000000   0.117300
//! H   0.000000   0.757200  -0.469200
//! H   0.000000  -0.757200  -0.469200
//! ```
//!
//! Line 1 is the atom count, line 2 a free-form comment (the reader
//! also scans it for an optional `charge N` / `mult M` hint), and the
//! remaining lines are `symbol x y z` in ångström.

use crate::element::Element;
use crate::error::{QchemError, Result};
use serde::{Deserialize, Serialize};

/// Bohr radii per ångström (CODATA 2018). Multiply ångström by this to
/// get bohr; divide bohr by this to get ångström.
pub const BOHR_PER_ANGSTROM: f64 = 1.889_726_124_625_770_2;

/// A single atom: an element and a Cartesian position in **bohr**.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Atom {
    /// The chemical element.
    pub element: Element,
    /// Cartesian coordinates `[x, y, z]` in bohr (atomic units).
    pub position: [f64; 3],
}

impl Atom {
    /// Construct an atom from an element and a bohr-unit position.
    pub fn new(element: Element, position: [f64; 3]) -> Self {
        Atom { element, position }
    }

    /// Construct an atom from a symbol and an ångström-unit position.
    ///
    /// # Errors
    ///
    /// Returns [`QchemError::Parse`] when the symbol is unknown.
    pub fn from_symbol_angstrom(symbol: &str, xyz_angstrom: [f64; 3]) -> Result<Self> {
        let element = Element::from_symbol(symbol)?;
        Ok(Atom {
            element,
            position: [
                xyz_angstrom[0] * BOHR_PER_ANGSTROM,
                xyz_angstrom[1] * BOHR_PER_ANGSTROM,
                xyz_angstrom[2] * BOHR_PER_ANGSTROM,
            ],
        })
    }

    /// The atom's position converted to ångström.
    pub fn position_angstrom(&self) -> [f64; 3] {
        [
            self.position[0] / BOHR_PER_ANGSTROM,
            self.position[1] / BOHR_PER_ANGSTROM,
            self.position[2] / BOHR_PER_ANGSTROM,
        ]
    }
}

/// A complete molecular specification.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct MolecularGeometry {
    /// The atoms, in input order. Basis functions are emitted in this
    /// order, so the ordering is stable and meaningful.
    pub atoms: Vec<Atom>,
    /// Total molecular charge in units of `e` (`0` for a neutral
    /// molecule, `+1` for a cation, `-1` for an anion).
    pub charge: i32,
    /// Spin multiplicity `2S + 1` (`1` = singlet, `2` = doublet,
    /// `3` = triplet, …). Always `>= 1`.
    pub multiplicity: u32,
}

impl MolecularGeometry {
    /// Build a neutral closed-shell singlet from a list of atoms.
    pub fn new(atoms: Vec<Atom>) -> Self {
        MolecularGeometry {
            atoms,
            charge: 0,
            multiplicity: 1,
        }
    }

    /// Build a geometry with an explicit charge and multiplicity.
    pub fn with_charge_multiplicity(atoms: Vec<Atom>, charge: i32, multiplicity: u32) -> Self {
        MolecularGeometry {
            atoms,
            charge,
            multiplicity,
        }
    }

    /// Number of atoms.
    #[inline]
    pub fn n_atoms(&self) -> usize {
        self.atoms.len()
    }

    /// Sum of the nuclear charges `Z` across all atoms.
    pub fn total_nuclear_charge(&self) -> i64 {
        self.atoms
            .iter()
            .map(|a| i64::from(a.element.atomic_number()))
            .sum()
    }

    /// Total number of electrons: the summed nuclear charge minus the
    /// molecular charge.
    ///
    /// # Errors
    ///
    /// Returns [`QchemError::InvalidInput`] when the molecule is empty
    /// or when the charge implies a negative electron count.
    pub fn n_electrons(&self) -> Result<u32> {
        if self.atoms.is_empty() {
            return Err(QchemError::invalid("molecule has no atoms"));
        }
        let signed = self.total_nuclear_charge() - i64::from(self.charge);
        if signed < 0 {
            return Err(QchemError::invalid(format!(
                "charge {} implies a negative electron count",
                self.charge
            )));
        }
        Ok(signed as u32)
    }

    /// The `(n_alpha, n_beta)` electron counts implied by the electron
    /// count and the multiplicity.
    ///
    /// With `N` electrons and multiplicity `2S + 1 = m`, the number of
    /// unpaired electrons is `m - 1`, so
    /// `n_alpha = (N + m - 1) / 2` and `n_beta = (N - m + 1) / 2`.
    ///
    /// # Errors
    ///
    /// Returns [`QchemError::InvalidInput`] when the multiplicity is
    /// zero, when `N` and `m` have incompatible parity (no valid
    /// determinant), or when the implied `n_beta` is negative.
    pub fn alpha_beta_counts(&self) -> Result<(u32, u32)> {
        let n = self.n_electrons()?;
        let m = self.multiplicity;
        if m == 0 {
            return Err(QchemError::invalid("multiplicity must be >= 1"));
        }
        let unpaired = m - 1;
        // N and (m-1) must have the same parity for N+m-1 to be even.
        if (n + unpaired) % 2 != 0 {
            return Err(QchemError::invalid(format!(
                "electron count {n} is incompatible with multiplicity {m} \
                 (parity mismatch — no valid determinant)"
            )));
        }
        if unpaired > n {
            return Err(QchemError::invalid(format!(
                "multiplicity {m} requires {unpaired} unpaired electrons \
                 but the molecule only has {n}"
            )));
        }
        let n_alpha = (n + unpaired) / 2;
        let n_beta = (n - unpaired) / 2;
        Ok((n_alpha, n_beta))
    }

    /// `true` when the molecule is a closed-shell singlet — an even
    /// electron count and multiplicity 1. Restricted Hartree-Fock
    /// applies exactly to this case.
    pub fn is_closed_shell(&self) -> bool {
        match (self.n_electrons(), self.multiplicity) {
            (Ok(n), 1) => n % 2 == 0,
            _ => false,
        }
    }

    /// Centre of mass in bohr, weighting each atom by its atomic mass.
    pub fn centre_of_mass(&self) -> [f64; 3] {
        let mut total = 0.0;
        let mut com = [0.0; 3];
        for a in &self.atoms {
            let m = a.element.atomic_mass();
            total += m;
            for (k, c) in com.iter_mut().enumerate() {
                *c += m * a.position[k];
            }
        }
        if total > 0.0 {
            for c in &mut com {
                *c /= total;
            }
        }
        com
    }

    /// Parse a molecule from xyz-format text. Coordinates are read in
    /// ångström and converted to bohr. The comment line is scanned for
    /// optional `charge <i>` and `mult <u>` tokens (any case); when
    /// absent the molecule defaults to neutral and singlet.
    ///
    /// # Errors
    ///
    /// Returns [`QchemError::Parse`] for a malformed header, a bad atom
    /// line, an unknown element or a count mismatch.
    pub fn from_xyz_str(text: &str) -> Result<Self> {
        let mut lines = text.lines();
        let count_line = lines
            .next()
            .ok_or_else(|| QchemError::parse("xyz", "empty input"))?;
        let declared: usize = count_line.trim().parse().map_err(|_| {
            QchemError::parse(
                "xyz",
                format!("first line `{count_line}` is not an atom count"),
            )
        })?;
        let comment = lines
            .next()
            .ok_or_else(|| QchemError::parse("xyz", "missing comment line"))?;
        let (charge, multiplicity) = parse_charge_mult(comment);

        // Do NOT pre-size from the declared count: it is caller-
        // controlled, so `Vec::with_capacity(declared)` on an
        // adversarial header (e.g. `99999999999`) would attempt a
        // multi-hundred-gigabyte allocation and abort the process. Grow
        // the vector on demand instead; the `atoms.len() != declared`
        // check below still rejects a mismatched header.
        let mut atoms: Vec<Atom> = Vec::new();
        for (i, line) in lines.enumerate() {
            if line.trim().is_empty() {
                continue;
            }
            let mut tok = line.split_whitespace();
            let sym = tok
                .next()
                .ok_or_else(|| QchemError::parse("xyz", format!("atom line {} is empty", i + 1)))?;
            let mut coord = [0.0f64; 3];
            for (k, c) in coord.iter_mut().enumerate() {
                let raw = tok.next().ok_or_else(|| {
                    QchemError::parse(
                        "xyz",
                        format!("atom line {} is missing coordinate {}", i + 1, k + 1),
                    )
                })?;
                *c = raw.parse().map_err(|_| {
                    QchemError::parse(
                        "xyz",
                        format!("atom line {}: `{raw}` is not a number", i + 1),
                    )
                })?;
            }
            atoms.push(Atom::from_symbol_angstrom(sym, coord)?);
        }

        if atoms.len() != declared {
            return Err(QchemError::parse(
                "xyz",
                format!(
                    "header declared {declared} atoms but {} were found",
                    atoms.len()
                ),
            ));
        }

        Ok(MolecularGeometry {
            atoms,
            charge,
            multiplicity,
        })
    }

    /// Render the molecule as xyz-format text. Coordinates are written
    /// in ångström; the comment line records the charge and
    /// multiplicity so a round trip through [`from_xyz_str`] is exact.
    ///
    /// [`from_xyz_str`]: MolecularGeometry::from_xyz_str
    pub fn to_xyz_str(&self) -> String {
        let mut out = String::new();
        out.push_str(&self.atoms.len().to_string());
        out.push('\n');
        out.push_str(&format!(
            "valenx-qchem geometry charge {} mult {}\n",
            self.charge, self.multiplicity
        ));
        for a in &self.atoms {
            let p = a.position_angstrom();
            out.push_str(&format!(
                "{:<3} {:>14.8} {:>14.8} {:>14.8}\n",
                a.element.symbol(),
                p[0],
                p[1],
                p[2]
            ));
        }
        out
    }

    /// Validate the geometry: at least one atom, a physically possible
    /// `(charge, multiplicity)` pair, and no two atoms coincident.
    ///
    /// # Errors
    ///
    /// Returns [`QchemError::InvalidInput`] describing the first
    /// problem found.
    pub fn validate(&self) -> Result<()> {
        // alpha_beta_counts performs the empty / charge / parity checks.
        self.alpha_beta_counts()?;
        for i in 0..self.atoms.len() {
            for j in (i + 1)..self.atoms.len() {
                let d2: f64 = (0..3)
                    .map(|k| {
                        let d = self.atoms[i].position[k] - self.atoms[j].position[k];
                        d * d
                    })
                    .sum();
                if d2 < 1.0e-8 {
                    return Err(QchemError::invalid(format!(
                        "atoms {i} and {j} are at the same point"
                    )));
                }
            }
        }
        Ok(())
    }
}

/// Scan a free-form xyz comment line for `charge <i>` and `mult <u>`
/// tokens (case-insensitive). Returns `(0, 1)` when neither is present.
fn parse_charge_mult(comment: &str) -> (i32, u32) {
    let mut charge = 0i32;
    let mut mult = 1u32;
    let lower = comment.to_ascii_lowercase();
    let toks: Vec<&str> = lower.split_whitespace().collect();
    for w in toks.windows(2) {
        match w[0] {
            "charge" => {
                if let Ok(v) = w[1].parse() {
                    charge = v;
                }
            }
            "mult" | "multiplicity" => {
                if let Ok(v) = w[1].parse() {
                    mult = v;
                }
            }
            _ => {}
        }
    }
    (charge, mult)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn h2() -> MolecularGeometry {
        // H2 at the STO-3G equilibrium, ~0.7414 Å bond.
        let atoms = vec![
            Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.0]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.7414]).unwrap(),
        ];
        MolecularGeometry::new(atoms)
    }

    #[test]
    fn angstrom_bohr_conversion() {
        let a = Atom::from_symbol_angstrom("H", [1.0, 0.0, 0.0]).unwrap();
        assert!((a.position[0] - BOHR_PER_ANGSTROM).abs() < 1.0e-10);
        let back = a.position_angstrom();
        assert!((back[0] - 1.0).abs() < 1.0e-10);
    }

    #[test]
    fn electron_counts() {
        let m = h2();
        assert_eq!(m.n_electrons().unwrap(), 2);
        assert_eq!(m.alpha_beta_counts().unwrap(), (1, 1));
        assert!(m.is_closed_shell());
    }

    #[test]
    fn cation_loses_an_electron() {
        let mut m = h2();
        m.charge = 1;
        m.multiplicity = 2;
        assert_eq!(m.n_electrons().unwrap(), 1);
        assert_eq!(m.alpha_beta_counts().unwrap(), (1, 0));
        assert!(!m.is_closed_shell());
    }

    #[test]
    fn parity_mismatch_is_rejected() {
        let mut m = h2();
        m.multiplicity = 2; // 2 electrons cannot be a doublet
        assert!(m.alpha_beta_counts().is_err());
    }

    #[test]
    fn triplet_oxygen_molecule() {
        // O2 ground state is a triplet: 16 electrons, 9 alpha 7 beta.
        let atoms = vec![
            Atom::from_symbol_angstrom("O", [0.0, 0.0, 0.0]).unwrap(),
            Atom::from_symbol_angstrom("O", [0.0, 0.0, 1.21]).unwrap(),
        ];
        let m = MolecularGeometry::with_charge_multiplicity(atoms, 0, 3);
        assert_eq!(m.n_electrons().unwrap(), 16);
        assert_eq!(m.alpha_beta_counts().unwrap(), (9, 7));
    }

    #[test]
    fn xyz_round_trip() {
        let m = h2();
        let text = m.to_xyz_str();
        let back = MolecularGeometry::from_xyz_str(&text).unwrap();
        assert_eq!(back.atoms.len(), 2);
        assert_eq!(back.charge, 0);
        assert_eq!(back.multiplicity, 1);
        for (a, b) in m.atoms.iter().zip(&back.atoms) {
            assert_eq!(a.element, b.element);
            for k in 0..3 {
                assert!((a.position[k] - b.position[k]).abs() < 1.0e-6);
            }
        }
    }

    #[test]
    fn xyz_reads_charge_mult_hint() {
        let text = "1\nradical cation charge 1 mult 2\nH 0.0 0.0 0.0\n";
        let m = MolecularGeometry::from_xyz_str(text).unwrap();
        assert_eq!(m.charge, 1);
        assert_eq!(m.multiplicity, 2);
    }

    #[test]
    fn xyz_count_mismatch_errors() {
        let text = "3\ncomment\nH 0 0 0\n";
        assert!(MolecularGeometry::from_xyz_str(text).is_err());
    }

    // ---- Input hardening: adversarial xyz input ----------------------

    #[test]
    fn xyz_absurd_atom_count_is_rejected_without_a_huge_allocation() {
        // A header claiming ~1e11 atoms must NOT drive a multi-hundred-
        // gigabyte `Vec::with_capacity` (which aborts the process).
        // The parser grows on demand; with only one real atom line the
        // count-mismatch check rejects it with a typed error.
        let text = "99999999999\ncomment\nH 0.0 0.0 0.0\n";
        let err = MolecularGeometry::from_xyz_str(text).expect_err("must reject");
        assert_eq!(err.category(), "parse");
    }

    #[test]
    fn xyz_garbage_and_truncated_input_never_panics() {
        // None of these may panic — each must be a typed error.
        for bad in [
            "",                        // empty
            "notanumber\n",            // bad count
            "2\n",                     // missing comment + atoms
            "2\ncomment\n",            // declared 2, zero atoms
            "1\ncomment\nH 0.0 0.0\n", // atom missing a coordinate
            "1\ncomment\nH x y z\n",   // non-numeric coordinates
            "1\ncomment\nXx 0 0 0\n",  // unknown element
            "\n\n\n",                  // only blank lines
        ] {
            let r = MolecularGeometry::from_xyz_str(bad);
            assert!(r.is_err(), "input {bad:?} should be a typed error");
        }
    }

    #[test]
    fn coincident_atoms_rejected() {
        let atoms = vec![
            Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.0]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.0]).unwrap(),
        ];
        let m = MolecularGeometry::new(atoms);
        assert!(m.validate().is_err());
    }

    #[test]
    fn centre_of_mass_of_symmetric_h2() {
        let com = h2().centre_of_mass();
        let mid = 0.7414 * BOHR_PER_ANGSTROM / 2.0;
        assert!((com[2] - mid).abs() < 1.0e-8);
    }
}
