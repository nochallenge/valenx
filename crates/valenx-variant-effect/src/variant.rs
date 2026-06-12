//! The variant model and its parser.
//!
//! A [`Variant`] is a single genetic substitution in one of two
//! coordinate systems:
//!
//! - a **protein** substitution, e.g. `p.R273H` (1-letter) or
//!   `p.Arg273His` (3-letter), naming a wild-type amino acid, a 1-based
//!   residue position, and a mutant amino acid;
//! - a **coding** (CDS) nucleotide substitution, e.g. `c.817C>T`, naming
//!   a 1-based coding position, the wild-type base and the mutant base.
//!
//! [`parse`] is a pure function — no I/O — that turns the HGVS-style
//! string into a validated [`Variant`], rejecting anything malformed with
//! [`VariantError::Parse`].

use serde::Serialize;

use crate::error::VariantError;

/// A single amino-acid residue (one of the 20 standard amino acids),
/// stored as an uppercase ASCII 1-letter code.
///
/// Serializes as its 1-letter code string (e.g. `"R"`).
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize)]
#[serde(into = "String")]
pub struct AminoAcid(u8);

impl From<AminoAcid> for String {
    fn from(a: AminoAcid) -> String {
        a.as_char().to_string()
    }
}

impl AminoAcid {
    /// The 20 standard amino acids, as uppercase 1-letter codes.
    const STANDARD: &'static [u8] = b"ACDEFGHIKLMNPQRSTVWY";

    /// Builds an [`AminoAcid`] from a 1-letter code (case-insensitive).
    ///
    /// Returns [`VariantError::InvalidResidue`] for any character outside
    /// the 20 standard amino acids (so `B`, `Z`, `X`, `*`, digits, etc.
    /// are rejected).
    pub fn from_one_letter(c: char) -> Result<Self, VariantError> {
        let upper = c.to_ascii_uppercase();
        if upper.is_ascii() && Self::STANDARD.contains(&(upper as u8)) {
            Ok(AminoAcid(upper as u8))
        } else {
            Err(VariantError::InvalidResidue { residue: c })
        }
    }

    /// Builds an [`AminoAcid`] from a 3-letter code (case-insensitive),
    /// e.g. `"Arg"`, `"his"`, `"GLY"`.
    ///
    /// Returns [`VariantError::InvalidResidue`] (carrying the first
    /// character of `code`, or `'?'` if empty) for any code outside the
    /// standard 20.
    pub fn from_three_letter(code: &str) -> Result<Self, VariantError> {
        match three_to_one(code) {
            Some(one) => Ok(AminoAcid(one)),
            None => Err(VariantError::InvalidResidue {
                residue: code.chars().next().unwrap_or('?'),
            }),
        }
    }

    /// The uppercase 1-letter code as a `char`.
    pub fn as_char(self) -> char {
        self.0 as char
    }

    /// The uppercase 1-letter code as a byte.
    pub fn as_byte(self) -> u8 {
        self.0
    }
}

/// A single nucleotide base (`A`, `C`, `G` or `T`), stored as an
/// uppercase ASCII byte. Coding-sequence variants are expressed on DNA,
/// so `U` is not accepted here.
///
/// Serializes as its base-letter string (e.g. `"C"`).
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize)]
#[serde(into = "String")]
pub struct Nucleotide(u8);

impl From<Nucleotide> for String {
    fn from(n: Nucleotide) -> String {
        n.as_char().to_string()
    }
}

impl Nucleotide {
    /// Builds a [`Nucleotide`] from a base character (case-insensitive).
    ///
    /// Accepts only the four canonical DNA bases `A C G T`; returns
    /// [`VariantError::InvalidResidue`] otherwise (including `U` and IUPAC
    /// ambiguity codes — a substitution names a single concrete base).
    pub fn from_char(c: char) -> Result<Self, VariantError> {
        match c.to_ascii_uppercase() {
            u @ ('A' | 'C' | 'G' | 'T') => Ok(Nucleotide(u as u8)),
            _ => Err(VariantError::InvalidResidue { residue: c }),
        }
    }

    /// The uppercase base as a `char`.
    pub fn as_char(self) -> char {
        self.0 as char
    }

    /// The uppercase base as a byte.
    pub fn as_byte(self) -> u8 {
        self.0
    }
}

/// A single genetic substitution.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize)]
pub enum Variant {
    /// A protein (amino-acid) substitution, e.g. `p.R273H`.
    ProteinSub {
        /// Wild-type amino acid.
        wt: AminoAcid,
        /// 1-based residue position in the protein.
        pos: usize,
        /// Mutant amino acid.
        mt: AminoAcid,
    },
    /// A coding-sequence (CDS) nucleotide substitution, e.g. `c.817C>T`.
    CodingSub {
        /// 1-based position in the coding sequence.
        pos: usize,
        /// Wild-type base.
        wt: Nucleotide,
        /// Mutant base.
        mt: Nucleotide,
    },
}

/// Maps a 3-letter amino-acid code (case-insensitive) to its uppercase
/// 1-letter code byte. Returns `None` outside the standard 20.
fn three_to_one(code: &str) -> Option<u8> {
    // Normalise to Titlecase-insensitive by comparing uppercased.
    let up = code.to_ascii_uppercase();
    let one: u8 = match up.as_str() {
        "ALA" => b'A',
        "ARG" => b'R',
        "ASN" => b'N',
        "ASP" => b'D',
        "CYS" => b'C',
        "GLN" => b'Q',
        "GLU" => b'E',
        "GLY" => b'G',
        "HIS" => b'H',
        "ILE" => b'I',
        "LEU" => b'L',
        "LYS" => b'K',
        "MET" => b'M',
        "PHE" => b'F',
        "PRO" => b'P',
        "SER" => b'S',
        "THR" => b'T',
        "TRP" => b'W',
        "TYR" => b'Y',
        "VAL" => b'V',
        _ => return None,
    };
    Some(one)
}

/// Splits a leading run of ASCII letters from the rest of `s`.
fn take_letters(s: &str) -> (&str, &str) {
    let end = s
        .find(|c: char| !c.is_ascii_alphabetic())
        .unwrap_or(s.len());
    s.split_at(end)
}

/// Splits a leading run of ASCII digits from the rest of `s`.
fn take_digits(s: &str) -> (&str, &str) {
    let end = s.find(|c: char| !c.is_ascii_digit()).unwrap_or(s.len());
    s.split_at(end)
}

/// Parses a 1-based position from a digit string, rejecting empty input
/// and `0` (HGVS positions are 1-based).
fn parse_pos(digits: &str, input: &str) -> Result<usize, VariantError> {
    if digits.is_empty() {
        return Err(VariantError::parse(input, "missing position"));
    }
    let pos: usize = digits
        .parse()
        .map_err(|_| VariantError::parse(input, format!("position `{digits}` is not a number")))?;
    if pos == 0 {
        return Err(VariantError::parse(
            input,
            "position must be 1-based (got 0)",
        ));
    }
    Ok(pos)
}

/// Builds a [`VariantError::Parse`] for a non-standard amino-acid token,
/// carrying the full `input` so callers see *which variant string* failed
/// (not just a bare residue). Used to remap an [`AminoAcid`]-construction
/// failure into a parse error with full context.
fn bad_aa(input: &str, residue: &str) -> VariantError {
    VariantError::parse(input, format!("`{residue}` is not a standard amino acid"))
}

/// Parses a single variant string in HGVS-style syntax.
///
/// Supported forms (surrounding whitespace is trimmed):
///
/// - protein, 1-letter: `p.R273H`
/// - protein, 3-letter: `p.Arg273His`
/// - coding (CDS): `c.817C>T`
///
/// Amino-acid codes must name one of the 20 standard residues; bases must
/// be one of `A C G T`. Anything else yields [`VariantError::Parse`].
pub fn parse(s: &str) -> Result<Variant, VariantError> {
    let trimmed = s.trim();
    if let Some(rest) = trimmed.strip_prefix("p.") {
        parse_protein(rest, trimmed)
    } else if let Some(rest) = trimmed.strip_prefix("c.") {
        parse_coding(rest, trimmed)
    } else {
        Err(VariantError::parse(
            trimmed,
            "expected a `p.` (protein) or `c.` (coding) prefix",
        ))
    }
}

/// Parses the body of a protein substitution (the part after `p.`).
fn parse_protein(body: &str, input: &str) -> Result<Variant, VariantError> {
    // Layout: <wt-letters><digits><mt-letters>
    let (wt_letters, rest) = take_letters(body);
    let (digits, mt_letters) = take_digits(rest);

    if wt_letters.is_empty() || mt_letters.is_empty() {
        return Err(VariantError::parse(
            input,
            "expected wild-type and mutant amino acids around the position",
        ));
    }

    let pos = parse_pos(digits, input)?;

    // The wild-type and mutant codes must use the same width; both are
    // either 1-letter or 3-letter. (A 1-letter wt with a 3-letter mt is
    // malformed.)
    let (wt, mt) = match (wt_letters.len(), mt_letters.len()) {
        (1, 1) => {
            let wt = AminoAcid::from_one_letter(wt_letters.chars().next().unwrap())
                .map_err(|_| bad_aa(input, wt_letters))?;
            let mt = AminoAcid::from_one_letter(mt_letters.chars().next().unwrap())
                .map_err(|_| bad_aa(input, mt_letters))?;
            (wt, mt)
        }
        (3, 3) => {
            let wt =
                AminoAcid::from_three_letter(wt_letters).map_err(|_| bad_aa(input, wt_letters))?;
            let mt =
                AminoAcid::from_three_letter(mt_letters).map_err(|_| bad_aa(input, mt_letters))?;
            (wt, mt)
        }
        _ => {
            return Err(VariantError::parse(
                input,
                "amino-acid codes must both be 1-letter or both be 3-letter",
            ));
        }
    };

    Ok(Variant::ProteinSub { wt, pos, mt })
}

/// Parses the body of a coding substitution (the part after `c.`),
/// e.g. `817C>T`.
fn parse_coding(body: &str, input: &str) -> Result<Variant, VariantError> {
    let (digits, rest) = take_digits(body);
    let pos = parse_pos(digits, input)?;

    // `rest` must be exactly `<base>><base>`.
    let (wt_str, after_wt) = take_letters(rest);
    let gt = after_wt;
    let Some(mt_str) = gt.strip_prefix('>') else {
        return Err(VariantError::parse(
            input,
            "coding substitution must be of the form `<pos><ref>><alt>` (e.g. `817C>T`)",
        ));
    };

    if wt_str.len() != 1 || mt_str.len() != 1 {
        return Err(VariantError::parse(
            input,
            "coding substitution names exactly one reference and one alternate base",
        ));
    }

    let wt = Nucleotide::from_char(wt_str.chars().next().unwrap()).map_err(|_| {
        VariantError::parse(input, format!("`{wt_str}` is not a DNA base (A/C/G/T)"))
    })?;
    let mt = Nucleotide::from_char(mt_str.chars().next().unwrap()).map_err(|_| {
        VariantError::parse(input, format!("`{mt_str}` is not a DNA base (A/C/G/T)"))
    })?;

    Ok(Variant::CodingSub { pos, wt, mt })
}

/// Transition/transversion (Ti/Tv) ratio over a slice of variants. A **transition** is a
/// purine↔purine (A↔G) or pyrimidine↔pyrimidine (C↔T) single-base substitution; a
/// **transversion** is a purine↔pyrimidine substitution. Only [`Variant::CodingSub`] with
/// differing bases are counted (protein and silent substitutions are skipped). Returns the
/// ratio Ti/Tv, or `0.0` when there are no transversions (otherwise the ratio is undefined).
pub fn transition_transversion_ratio(variants: &[Variant]) -> f64 {
    let mut transitions = 0u32;
    let mut transversions = 0u32;
    for variant in variants {
        if let Variant::CodingSub { wt, mt, .. } = variant {
            if wt == mt {
                continue;
            }
            let is_purine = |n: &Nucleotide| matches!(n.as_char(), 'A' | 'G');
            if is_purine(wt) == is_purine(mt) {
                transitions += 1;
            } else {
                transversions += 1;
            }
        }
    }
    if transversions == 0 {
        0.0
    } else {
        f64::from(transitions) / f64::from(transversions)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn aa(c: char) -> AminoAcid {
        AminoAcid::from_one_letter(c).unwrap()
    }
    fn nt(c: char) -> Nucleotide {
        Nucleotide::from_char(c).unwrap()
    }

    // ---- amino-acid / nucleotide newtypes --------------------------------

    #[test]
    fn amino_acid_one_letter_roundtrip() {
        for &b in AminoAcid::STANDARD {
            let a = AminoAcid::from_one_letter(b as char).unwrap();
            assert_eq!(a.as_byte(), b);
            assert_eq!(a.as_char(), b as char);
        }
    }

    #[test]
    fn amino_acid_one_letter_is_case_insensitive() {
        assert_eq!(AminoAcid::from_one_letter('r').unwrap(), aa('R'));
    }

    #[test]
    fn amino_acid_rejects_non_standard() {
        // B, Z, J, X, O, U, the stop `*`, and digits are not among the 20.
        for c in ['B', 'Z', 'J', 'X', 'O', 'U', '*', '1', '-'] {
            assert!(
                AminoAcid::from_one_letter(c).is_err(),
                "expected `{c}` to be rejected"
            );
        }
    }

    #[test]
    fn three_letter_map_covers_all_twenty() {
        let pairs = [
            ("Ala", 'A'),
            ("Arg", 'R'),
            ("Asn", 'N'),
            ("Asp", 'D'),
            ("Cys", 'C'),
            ("Gln", 'Q'),
            ("Glu", 'E'),
            ("Gly", 'G'),
            ("His", 'H'),
            ("Ile", 'I'),
            ("Leu", 'L'),
            ("Lys", 'K'),
            ("Met", 'M'),
            ("Phe", 'F'),
            ("Pro", 'P'),
            ("Ser", 'S'),
            ("Thr", 'T'),
            ("Trp", 'W'),
            ("Tyr", 'Y'),
            ("Val", 'V'),
        ];
        for (three, one) in pairs {
            assert_eq!(
                AminoAcid::from_three_letter(three).unwrap(),
                aa(one),
                "{three}"
            );
            // case-insensitive
            assert_eq!(
                AminoAcid::from_three_letter(&three.to_uppercase()).unwrap(),
                aa(one)
            );
            assert_eq!(
                AminoAcid::from_three_letter(&three.to_lowercase()).unwrap(),
                aa(one)
            );
        }
    }

    #[test]
    fn three_letter_rejects_unknown() {
        assert!(AminoAcid::from_three_letter("Xaa").is_err());
        assert!(AminoAcid::from_three_letter("Sec").is_err()); // selenocysteine
        assert!(AminoAcid::from_three_letter("").is_err());
    }

    #[test]
    fn nucleotide_accepts_acgt_only() {
        for c in ['A', 'C', 'G', 'T', 'a', 'c', 'g', 't'] {
            assert!(Nucleotide::from_char(c).is_ok(), "{c}");
        }
        for c in ['U', 'u', 'N', 'R', 'X', '1'] {
            assert!(Nucleotide::from_char(c).is_err(), "{c}");
        }
    }

    // ---- protein parsing -------------------------------------------------

    #[test]
    fn parse_protein_one_letter() {
        let v = parse("p.R273H").unwrap();
        assert_eq!(
            v,
            Variant::ProteinSub {
                wt: aa('R'),
                pos: 273,
                mt: aa('H'),
            }
        );
    }

    #[test]
    fn parse_protein_three_letter_equals_one_letter() {
        let three = parse("p.Arg273His").unwrap();
        let one = parse("p.R273H").unwrap();
        assert_eq!(three, one);
    }

    #[test]
    fn parse_protein_three_letter_case_insensitive() {
        assert_eq!(parse("p.ARG273HIS").unwrap(), parse("p.R273H").unwrap());
        assert_eq!(parse("p.arg273his").unwrap(), parse("p.R273H").unwrap());
    }

    #[test]
    fn parse_protein_synonymous_like_same_residue() {
        // wt == mt is syntactically valid (a "silent" protein notation).
        let v = parse("p.A100A").unwrap();
        assert_eq!(
            v,
            Variant::ProteinSub {
                wt: aa('A'),
                pos: 100,
                mt: aa('A')
            }
        );
    }

    #[test]
    fn parse_protein_rejects_non_standard_aa() {
        // B is not one of the 20.
        assert!(matches!(parse("p.B100A"), Err(VariantError::Parse { .. })));
        assert!(matches!(parse("p.A100B"), Err(VariantError::Parse { .. })));
        // 'X' (any) is not accepted by this strict parser.
        assert!(parse("p.X100A").is_err());
    }

    #[test]
    fn parse_protein_rejects_mixed_letter_widths() {
        // 1-letter wt with 3-letter mt (or vice-versa) is malformed.
        assert!(parse("p.R273His").is_err());
        assert!(parse("p.Arg273H").is_err());
    }

    #[test]
    fn parse_protein_rejects_missing_or_zero_position() {
        assert!(parse("p.RH").is_err()); // no digits
        assert!(parse("p.R0H").is_err()); // 1-based, 0 illegal
    }

    // ---- coding parsing --------------------------------------------------

    #[test]
    fn parse_coding_basic() {
        let v = parse("c.817C>T").unwrap();
        assert_eq!(
            v,
            Variant::CodingSub {
                pos: 817,
                wt: nt('C'),
                mt: nt('T')
            }
        );
    }

    #[test]
    fn parse_coding_is_case_insensitive() {
        assert_eq!(parse("c.817c>t").unwrap(), parse("c.817C>T").unwrap());
    }

    #[test]
    fn parse_coding_position_is_decimal_not_octal() {
        // Leading zeros must be parsed as decimal. `007` is 7 (not an octal
        // literal), and the decimal-vs-octal distinction is made explicit
        // by `010`, which is 10 decimal (it would be 8 if read as octal).
        assert_eq!(
            parse("c.007C>T").unwrap(),
            Variant::CodingSub {
                pos: 7,
                wt: nt('C'),
                mt: nt('T')
            }
        );
        assert_eq!(
            parse("c.010C>T").unwrap(),
            Variant::CodingSub {
                pos: 10,
                wt: nt('C'),
                mt: nt('T')
            }
        );
    }

    #[test]
    fn parse_coding_rejects_bad_shape() {
        assert!(parse("c.817CT").is_err()); // missing '>'
        assert!(parse("c.817C>").is_err()); // missing alt
        assert!(parse("c.817>T").is_err()); // missing ref
        assert!(parse("c.C>T").is_err()); // missing position
        assert!(parse("c.817CC>T").is_err()); // ref is 2 bases
        assert!(parse("c.817C>U").is_err()); // U is not DNA
        assert!(parse("c.0C>T").is_err()); // 1-based
    }

    // ---- prefix / garbage ------------------------------------------------

    #[test]
    fn parse_rejects_missing_prefix_and_garbage() {
        assert!(matches!(parse("R273H"), Err(VariantError::Parse { .. })));
        assert!(parse("").is_err());
        assert!(parse("g.123A>T").is_err()); // genomic prefix unsupported
        assert!(parse("p.").is_err());
        assert!(parse("c.").is_err());
        assert!(parse("hello world").is_err());
    }

    #[test]
    fn parse_trims_surrounding_whitespace() {
        assert_eq!(parse("  p.R273H \n").unwrap(), parse("p.R273H").unwrap());
    }

    #[test]
    fn parse_error_carries_the_input() {
        match parse("p.B100A") {
            Err(VariantError::Parse { input, .. }) => assert_eq!(input, "p.B100A"),
            other => panic!("expected Parse error, got {other:?}"),
        }
    }

    #[test]
    fn transition_transversion_ratio_classifies_snvs() {
        // 2 transitions (A>G, C>T) : 1 transversion (A>C) = 2.0.
        let v = [
            parse("c.1A>G").unwrap(),
            parse("c.2C>T").unwrap(),
            parse("c.3A>C").unwrap(),
        ];
        assert!((transition_transversion_ratio(&v) - 2.0).abs() < 1e-9);
        // No transversions → 0.0 sentinel (ratio undefined).
        let ti_only = [parse("c.1A>G").unwrap(), parse("c.2C>T").unwrap()];
        assert_eq!(transition_transversion_ratio(&ti_only), 0.0);
        // Empty → 0.0; a protein substitution is ignored (no coding SNV present).
        assert_eq!(transition_transversion_ratio(&[]), 0.0);
        assert_eq!(
            transition_transversion_ratio(&[parse("p.R273H").unwrap()]),
            0.0
        );
        // Mixed batch: A>G, A>T, C>T, G>C, G>A = 3 Ti : 2 Tv = 1.5.
        let mixed = [
            parse("c.1A>G").unwrap(),
            parse("c.2A>T").unwrap(),
            parse("c.3C>T").unwrap(),
            parse("c.4G>C").unwrap(),
            parse("c.5G>A").unwrap(),
        ];
        assert!((transition_transversion_ratio(&mixed) - 1.5).abs() < 1e-9);
    }
}
