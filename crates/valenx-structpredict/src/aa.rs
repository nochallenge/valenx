//! Amino-acid constants shared across the crate.
//!
//! A small, dependency-free table of the 20 standard amino acids:
//! one-letter ↔ three-letter codes, a hydropathy scale, an
//! approximate sidechain volume, and the per-residue rotamer-count
//! used by the design and repacking search. Every other module keys
//! into this table rather than re-deriving it.
//!
//! The numbers are the textbook published values: the hydropathy
//! scale is Kyte-Doolittle (1982); sidechain volumes are the
//! Zamyatnin (1972) residue volumes minus the backbone contribution
//! (Å³, rounded). They are used as *features* of a knowledge-based
//! score, not as a force field — exact parameterisation is not the
//! point of a classical v1.

/// The 20 standard amino acids, in a fixed canonical order. Indexing
/// into this slice gives a stable per-amino-acid integer id.
pub const AMINO_ACIDS: [char; 20] = [
    'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y',
    'V',
];

/// Number of standard amino acids.
pub const N_AA: usize = 20;

/// Maps a one-letter amino-acid code to its index in [`AMINO_ACIDS`]
/// (`0..20`). `None` for any character that is not one of the 20
/// standard codes (lower-case included — pass an upper-cased char).
pub fn aa_index(c: char) -> Option<usize> {
    AMINO_ACIDS.iter().position(|&a| a == c)
}

/// Maps a three-letter residue name (`"ALA"`, `"GLY"`, …) to its
/// one-letter code. Recognises the 20 standard residues plus the
/// common modified residue `MSE` (selenomethionine → `M`). `None`
/// otherwise.
pub fn three_to_one(name: &str) -> Option<char> {
    Some(match name.trim().to_ascii_uppercase().as_str() {
        "ALA" => 'A',
        "ARG" => 'R',
        "ASN" => 'N',
        "ASP" => 'D',
        "CYS" => 'C',
        "GLN" => 'Q',
        "GLU" => 'E',
        "GLY" => 'G',
        "HIS" => 'H',
        "ILE" => 'I',
        "LEU" => 'L',
        "LYS" => 'K',
        "MET" | "MSE" => 'M',
        "PHE" => 'F',
        "PRO" => 'P',
        "SER" => 'S',
        "THR" => 'T',
        "TRP" => 'W',
        "TYR" => 'Y',
        "VAL" => 'V',
        _ => return None,
    })
}

/// Maps a one-letter amino-acid code to its standard three-letter
/// residue name. `None` for a non-standard code.
pub fn one_to_three(c: char) -> Option<&'static str> {
    Some(match c {
        'A' => "ALA",
        'R' => "ARG",
        'N' => "ASN",
        'D' => "ASP",
        'C' => "CYS",
        'Q' => "GLN",
        'E' => "GLU",
        'G' => "GLY",
        'H' => "HIS",
        'I' => "ILE",
        'L' => "LEU",
        'K' => "LYS",
        'M' => "MET",
        'F' => "PHE",
        'P' => "PRO",
        'S' => "SER",
        'T' => "THR",
        'W' => "TRP",
        'Y' => "TYR",
        'V' => "VAL",
        _ => return None,
    })
}

/// Kyte-Doolittle hydropathy of an amino acid: positive is
/// hydrophobic, negative hydrophilic. Used as a burial-preference
/// feature by the design score. Returns `0.0` for an unknown code.
pub fn hydropathy(c: char) -> f64 {
    match c {
        'I' => 4.5,
        'V' => 4.2,
        'L' => 3.8,
        'F' => 2.8,
        'C' => 2.5,
        'M' => 1.9,
        'A' => 1.8,
        'G' => -0.4,
        'T' => -0.7,
        'S' => -0.8,
        'W' => -0.9,
        'Y' => -1.3,
        'P' => -1.6,
        'H' => -3.2,
        'E' => -3.5,
        'Q' => -3.5,
        'D' => -3.5,
        'N' => -3.5,
        'K' => -3.9,
        'R' => -4.5,
        _ => 0.0,
    }
}

/// Approximate sidechain heavy-atom volume in Å³. Used by the
/// packing-quality score and the design search's cavity term. `0.0`
/// for glycine (no sidechain) and for an unknown code.
pub fn sidechain_volume(c: char) -> f64 {
    match c {
        'G' => 0.0,
        'A' => 25.0,
        'S' => 36.0,
        'C' => 49.0,
        'D' => 60.0,
        'P' => 56.0,
        'N' => 64.0,
        'T' => 61.0,
        'E' => 85.0,
        'V' => 76.0,
        'Q' => 89.0,
        'H' => 92.0,
        'M' => 95.0,
        'I' => 102.0,
        'L' => 102.0,
        'K' => 102.0,
        'R' => 124.0,
        'F' => 122.0,
        'Y' => 130.0,
        'W' => 163.0,
        _ => 0.0,
    }
}

/// Number of χ (chi) dihedral angles in the sidechain of an amino
/// acid. `0` for `A` and `G` (no rotatable sidechain dihedral).
/// Drives the size of the per-residue rotamer set in the repacking
/// and design search.
pub fn chi_count(c: char) -> usize {
    match c {
        'A' | 'G' => 0,
        'S' | 'T' | 'C' | 'V' | 'P' => 1,
        'I' | 'L' | 'D' | 'N' | 'F' | 'Y' | 'H' | 'W' => 2,
        'M' | 'E' | 'Q' => 3,
        'K' | 'R' => 4,
        _ => 0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip_codes() {
        for &c in &AMINO_ACIDS {
            let three = one_to_three(c).expect("known one-letter");
            assert_eq!(three_to_one(three), Some(c));
            let idx = aa_index(c).expect("indexed");
            assert!(idx < N_AA);
        }
    }

    #[test]
    fn mse_maps_to_met() {
        assert_eq!(three_to_one("MSE"), Some('M'));
        assert_eq!(three_to_one("mse"), Some('M'));
    }

    #[test]
    fn hydropathy_sign_is_right() {
        assert!(hydropathy('I') > 0.0, "Ile is hydrophobic");
        assert!(hydropathy('R') < 0.0, "Arg is hydrophilic");
    }

    #[test]
    fn glycine_has_no_sidechain() {
        assert_eq!(sidechain_volume('G'), 0.0);
        assert_eq!(chi_count('G'), 0);
        assert!(sidechain_volume('W') > sidechain_volume('A'));
    }
}
