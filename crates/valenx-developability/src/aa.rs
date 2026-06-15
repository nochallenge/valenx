//! Amino-acid constants: the Kyte-Doolittle hydropathy scale and a standard
//! pKa set for the ionizable groups.

/// The 20 standard amino acids in canonical column order.
pub const ALPHABET: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";

/// Kyte-Doolittle (1982) hydropathy index, in [`ALPHABET`] order. Positive is
/// hydrophobic.
const HYDROPATHY: [f64; 20] = [
    1.8,  // A
    2.5,  // C
    -3.5, // D
    -3.5, // E
    2.8,  // F
    -0.4, // G
    -3.2, // H
    4.5,  // I
    -3.9, // K
    3.8,  // L
    1.9,  // M
    -3.5, // N
    -1.6, // P
    -3.5, // Q
    -4.5, // R
    -0.8, // S
    -0.7, // T
    4.2,  // V
    -0.9, // W
    -1.3, // Y
];

/// pKa of the N-terminal alpha-amino group (EMBOSS set).
pub const PKA_N_TERM: f64 = 8.6;
/// pKa of the C-terminal alpha-carboxyl group (EMBOSS set).
pub const PKA_C_TERM: f64 = 3.6;

/// Column index of `residue` (either case) in `0..20`, or `None` if not a
/// standard amino acid.
pub fn aa_index(residue: u8) -> Option<usize> {
    let up = residue.to_ascii_uppercase();
    ALPHABET.iter().position(|&a| a == up)
}

/// Kyte-Doolittle hydropathy of `residue`, or `None` if non-standard.
pub fn hydropathy(residue: u8) -> Option<f64> {
    aa_index(residue).map(|i| HYDROPATHY[i])
}

/// pKa of a positively-ionizable side chain (His/Lys/Arg), or `None`.
/// These groups carry `+1` when protonated (pH below pKa).
pub fn positive_sidechain_pka(residue: u8) -> Option<f64> {
    match residue.to_ascii_uppercase() {
        b'H' => Some(6.5),
        b'K' => Some(10.8),
        b'R' => Some(12.5),
        _ => None,
    }
}

/// pKa of a negatively-ionizable side chain (Cys/Asp/Glu/Tyr), or `None`.
/// These groups carry `-1` when deprotonated (pH above pKa).
pub fn negative_sidechain_pka(residue: u8) -> Option<f64> {
    match residue.to_ascii_uppercase() {
        b'C' => Some(8.5),
        b'D' => Some(3.9),
        b'E' => Some(4.1),
        b'Y' => Some(10.1),
        _ => None,
    }
}

/// The first non-standard residue as `(position, char)`, if any.
pub fn first_invalid(seq: &str) -> Option<(usize, char)> {
    seq.bytes()
        .enumerate()
        .find(|&(_, b)| aa_index(b).is_none())
        .map(|(i, b)| (i, char::from(b)))
}
