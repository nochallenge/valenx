//! AutoDock-4 atom type taxonomy + physical properties.
//!
//! Vina uses AD4's 22 atom-type symbols rather than just the element.
//! The distinction matters because, e.g., a polar `OA` (acceptor)
//! oxygen participates in the H-bond term but an aliphatic `C` does
//! not. Vina's source (`atom_constants.h`) is the ground truth.

use std::str::FromStr;

/// AutoDock-4 atom type. Variants follow the AD4 single/double letter
/// codes; aliases like `Cl` -> `Cl` are mapped in [`FromStr`].
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
#[allow(missing_docs)]
pub enum Ad4AtomType {
    /// Aliphatic / sp3 carbon.
    C,
    /// Aromatic carbon.
    A,
    /// Nitrogen, non-acceptor.
    N,
    /// Nitrogen, H-bond acceptor (lone pair).
    NA,
    /// Nitrogen, sp2 acceptor — Vina folds this into NA for scoring.
    NS,
    /// Oxygen, H-bond acceptor.
    OA,
    /// Oxygen, sp2 acceptor (Vina folds into OA for scoring).
    OS,
    /// Sulfur.
    S,
    /// Sulfur, acceptor.
    SA,
    /// Polar (donor) hydrogen.
    HD,
    /// Non-polar hydrogen (skipped in Vina — usually merged into heavy).
    H,
    /// Phosphorus.
    P,
    /// Fluorine.
    F,
    /// Chlorine.
    Cl,
    /// Bromine.
    Br,
    /// Iodine.
    I,
    /// Metal (Zn, Mg, Ca, Mn, Fe, Cu, K, Na — Vina treats them all the same).
    Metal,
}

/// Physical properties used by the scoring function.
#[derive(Copy, Clone, Debug)]
pub struct AtomProps {
    /// Van der Waals radius in Å (Vina's `xs_radius`).
    pub vdw_radius: f64,
    /// True if hydrophobic (participates in [`crate::score::hydrophobic_pair`]).
    pub is_hydrophobic: bool,
    /// True if H-bond donor (Vina maps via `xs_h_bond_donor`).
    pub is_donor: bool,
    /// True if H-bond acceptor.
    pub is_acceptor: bool,
}

impl Ad4AtomType {
    /// Look up [`AtomProps`] for this type. Values from Vina's
    /// `atom_constants.h` — DO NOT edit without cross-checking upstream.
    pub fn props(self) -> AtomProps {
        use Ad4AtomType::*;
        match self {
            C => AtomProps {
                vdw_radius: 1.9,
                is_hydrophobic: true,
                is_donor: false,
                is_acceptor: false,
            },
            A => AtomProps {
                vdw_radius: 1.9,
                is_hydrophobic: true,
                is_donor: false,
                is_acceptor: false,
            },
            N => AtomProps {
                vdw_radius: 1.8,
                is_hydrophobic: false,
                is_donor: false,
                is_acceptor: false,
            },
            NA => AtomProps {
                vdw_radius: 1.8,
                is_hydrophobic: false,
                is_donor: false,
                is_acceptor: true,
            },
            NS => AtomProps {
                vdw_radius: 1.8,
                is_hydrophobic: false,
                is_donor: false,
                is_acceptor: true,
            },
            OA => AtomProps {
                vdw_radius: 1.7,
                is_hydrophobic: false,
                is_donor: false,
                is_acceptor: true,
            },
            OS => AtomProps {
                vdw_radius: 1.7,
                is_hydrophobic: false,
                is_donor: false,
                is_acceptor: true,
            },
            S => AtomProps {
                vdw_radius: 2.0,
                is_hydrophobic: false,
                is_donor: false,
                is_acceptor: false,
            },
            SA => AtomProps {
                vdw_radius: 2.0,
                is_hydrophobic: false,
                is_donor: false,
                is_acceptor: true,
            },
            HD => AtomProps {
                vdw_radius: 0.0,
                is_hydrophobic: false,
                is_donor: true,
                is_acceptor: false,
            },
            H => AtomProps {
                vdw_radius: 0.0,
                is_hydrophobic: false,
                is_donor: false,
                is_acceptor: false,
            },
            P => AtomProps {
                vdw_radius: 2.1,
                is_hydrophobic: false,
                is_donor: false,
                is_acceptor: false,
            },
            F => AtomProps {
                vdw_radius: 1.5,
                is_hydrophobic: true,
                is_donor: false,
                is_acceptor: false,
            },
            Cl => AtomProps {
                vdw_radius: 1.8,
                is_hydrophobic: true,
                is_donor: false,
                is_acceptor: false,
            },
            Br => AtomProps {
                vdw_radius: 2.0,
                is_hydrophobic: true,
                is_donor: false,
                is_acceptor: false,
            },
            I => AtomProps {
                vdw_radius: 2.2,
                is_hydrophobic: true,
                is_donor: false,
                is_acceptor: false,
            },
            Metal => AtomProps {
                vdw_radius: 1.2,
                is_hydrophobic: false,
                is_donor: true,
                is_acceptor: false,
            },
        }
    }
}

/// Parse the two-character AD4 atom-type code that appears in column
/// 78-79 of a PDBQT ATOM/HETATM record.
///
/// **AD4 sodium / nitrogen-acceptor collision.** The two-letter code
/// `NA` is ambiguous between sodium (Na) and nitrogen-acceptor (NA).
/// valenx-dock follows Vina's case-insensitive interpretation where
/// `NA` always means nitrogen-acceptor (the field is upper-cased by
/// PDBQT writers in practice). Sodium atoms in receptor files should
/// be typed as `Na` (mixed case) by the upstream preparation tool,
/// which we parse as Metal. AutoDockTools 4.2 emits sodium as `Na`
/// by default, so real-world PDBQTs round-trip correctly.
impl FromStr for Ad4AtomType {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        use Ad4AtomType::*;
        match s.trim() {
            "C" => Ok(C),
            "A" => Ok(A),
            "N" => Ok(N),
            "NA" => Ok(NA),
            "NS" => Ok(NS),
            "OA" => Ok(OA),
            "OS" => Ok(OS),
            "S" => Ok(S),
            "SA" => Ok(SA),
            "HD" => Ok(HD),
            "H" => Ok(H),
            "P" => Ok(P),
            "F" => Ok(F),
            "Cl" | "CL" => Ok(Cl),
            "Br" | "BR" => Ok(Br),
            "I" => Ok(I),
            // Vina folds every metal into one slot. Sodium is `Na`
            // (mixed case); `NA` always means nitrogen-acceptor.
            "Mg" | "MG" | "Zn" | "ZN" | "Ca" | "CA" | "Mn" | "MN" | "Fe" | "FE" | "Cu" | "CU"
            | "K" | "Na" => Ok(Metal),
            other => Err(format!("unknown AD4 atom type `{other}`")),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_canonical_types() {
        assert_eq!("C".parse::<Ad4AtomType>().unwrap(), Ad4AtomType::C);
        assert_eq!("OA".parse::<Ad4AtomType>().unwrap(), Ad4AtomType::OA);
        assert_eq!("HD".parse::<Ad4AtomType>().unwrap(), Ad4AtomType::HD);
        assert_eq!("Cl".parse::<Ad4AtomType>().unwrap(), Ad4AtomType::Cl);
    }

    #[test]
    fn rejects_unknown() {
        let err = "ZZ".parse::<Ad4AtomType>().unwrap_err();
        assert!(err.contains("ZZ"));
    }

    #[test]
    fn carbon_is_hydrophobic_oxygen_acceptor_is_not() {
        assert!(Ad4AtomType::C.props().is_hydrophobic);
        assert!(!Ad4AtomType::OA.props().is_hydrophobic);
        assert!(Ad4AtomType::OA.props().is_acceptor);
        assert!(!Ad4AtomType::C.props().is_acceptor);
    }

    #[test]
    fn donor_hydrogen_has_zero_vdw() {
        // Vina collapses HD's VDW to zero so the donor-acceptor surface
        // distance reduces to the heavy-acceptor radius alone.
        assert_eq!(Ad4AtomType::HD.props().vdw_radius, 0.0);
        assert!(Ad4AtomType::HD.props().is_donor);
    }

    #[test]
    fn na_means_nitrogen_acceptor_not_sodium() {
        // The `NA` two-letter code is reserved for nitrogen-acceptor;
        // sodium must be spelled `Na` (mixed case) per the AD4
        // convention documented above the FromStr impl.
        assert_eq!("NA".parse::<Ad4AtomType>().unwrap(), Ad4AtomType::NA);
        assert_eq!("Na".parse::<Ad4AtomType>().unwrap(), Ad4AtomType::Metal);
    }
}
