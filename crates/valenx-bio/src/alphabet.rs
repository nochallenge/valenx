//! IUPAC alphabet validation for nucleotide + protein sequences.
//!
//! Three alphabets ship today:
//! - [`Alphabet::Dna`] — `ACGT` plus IUPAC ambiguity codes
//!   (`N` = any, `R` = puRine A/G, `Y` = pYrimidine C/T, …).
//! - [`Alphabet::Rna`] — same as DNA with `U` substituted for `T`.
//! - [`Alphabet::Protein`] — 20 canonical amino-acid one-letter
//!   codes + `X` (any), `B` (Asx = D/N), `Z` (Glx = E/Q),
//!   `J` (Xle = L/I), `U` (selenocysteine), `O` (pyrrolysine), `*` (stop).

use serde::{Deserialize, Serialize};

/// The three sequence alphabets the crate accepts.
///
/// IUPAC ambiguity codes (`N`, `R`, `Y`, …) and the gap character `-`
/// are valid in every alphabet; see [`Alphabet::is_valid`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum Alphabet {
    /// DNA alphabet: `A`, `C`, `G`, `T` + IUPAC ambiguity codes.
    Dna,
    /// RNA alphabet: `A`, `C`, `G`, `U` + IUPAC ambiguity codes.
    Rna,
    /// Protein alphabet: 20 canonical amino-acid one-letter codes plus
    /// `X` (any), `B`/`Z`/`J` (ambiguity), `U` (selenocysteine), `O`
    /// (pyrrolysine), and `*` (stop).
    Protein,
}

impl Alphabet {
    /// `true` if `byte` is a valid character for this alphabet.
    /// Comparisons are case-insensitive.
    pub fn is_valid(self, byte: u8) -> bool {
        let upper = byte.to_ascii_uppercase();
        match self {
            Self::Dna => matches!(
                upper,
                b'A' | b'C'
                    | b'G'
                    | b'T'
                    | b'N'
                    | b'R'
                    | b'Y'
                    | b'M'
                    | b'K'
                    | b'S'
                    | b'W'
                    | b'B'
                    | b'D'
                    | b'H'
                    | b'V'
                    | b'-'
            ),
            Self::Rna => matches!(
                upper,
                b'A' | b'C'
                    | b'G'
                    | b'U'
                    | b'N'
                    | b'R'
                    | b'Y'
                    | b'M'
                    | b'K'
                    | b'S'
                    | b'W'
                    | b'B'
                    | b'D'
                    | b'H'
                    | b'V'
                    | b'-'
            ),
            Self::Protein => matches!(
                upper,
                b'A' | b'C'
                    | b'D'
                    | b'E'
                    | b'F'
                    | b'G'
                    | b'H'
                    | b'I'
                    | b'K'
                    | b'L'
                    | b'M'
                    | b'N'
                    | b'P'
                    | b'Q'
                    | b'R'
                    | b'S'
                    | b'T'
                    | b'V'
                    | b'W'
                    | b'Y'
                    | b'X'
                    | b'B'
                    | b'Z'
                    | b'J'
                    | b'U'
                    | b'O'
                    | b'*'
                    | b'-'
            ),
        }
    }

    /// Stable string identifier used by serde + tracing fields.
    pub fn id(self) -> &'static str {
        match self {
            Self::Dna => "dna",
            Self::Rna => "rna",
            Self::Protein => "protein",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dna_accepts_iupac_unambiguous_codes() {
        let alpha = Alphabet::Dna;
        for c in b"ACGT" {
            assert!(alpha.is_valid(*c), "rejected canonical {}", *c as char);
        }
    }

    #[test]
    fn dna_rejects_protein_residues() {
        let alpha = Alphabet::Dna;
        assert!(!alpha.is_valid(b'E')); // glutamate — protein only
        assert!(!alpha.is_valid(b'X')); // ambiguous protein only
    }

    #[test]
    fn protein_accepts_20_canonical() {
        let alpha = Alphabet::Protein;
        for c in b"ACDEFGHIKLMNPQRSTVWY" {
            assert!(alpha.is_valid(*c), "rejected {}", *c as char);
        }
    }

    #[test]
    fn rna_accepts_u_not_t() {
        assert!(Alphabet::Rna.is_valid(b'U'));
        assert!(!Alphabet::Rna.is_valid(b'T'));
    }

    #[test]
    fn gap_is_valid_in_all_alphabets() {
        assert!(Alphabet::Dna.is_valid(b'-'));
        assert!(Alphabet::Rna.is_valid(b'-'));
        assert!(Alphabet::Protein.is_valid(b'-'));
    }
}
