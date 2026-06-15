//! Minimal amino-acid validation shared by the similarity routines.

/// The 20 standard amino acids.
pub const ALPHABET: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";

/// Whether `residue` (either case) is one of the 20 standard amino acids.
pub fn is_standard(residue: u8) -> bool {
    ALPHABET.contains(&residue.to_ascii_uppercase())
}

/// Validate that `seq` is non-empty and all-standard, returning the offending
/// `(position, char)` on the first non-standard residue.
pub fn first_invalid(seq: &str) -> Option<(usize, char)> {
    seq.bytes()
        .enumerate()
        .find(|&(_, b)| !is_standard(b))
        .map(|(i, b)| (i, char::from(b)))
}
