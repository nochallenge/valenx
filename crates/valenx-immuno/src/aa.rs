//! The 20-letter amino-acid alphabet and residue indexing.
//!
//! [`Pssm`](crate::matrix::Pssm) stores one weight per position per residue;
//! [`aa_index`] maps a residue byte to that residue's column.

/// The 20 standard amino acids in a fixed canonical column order.
///
/// Column `i` of a [`Pssm`](crate::matrix::Pssm) row corresponds to
/// `ALPHABET[i]`.
pub const ALPHABET: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";

/// The number of standard amino acids (20).
pub const N_AA: usize = 20;

/// Map a residue byte to its column index in `0..20`, accepting either case.
///
/// Returns `None` for any byte that is not one of the 20 standard amino acids
/// (so ambiguity codes `B`, `Z`, `X`, the stop `*`, and gaps are rejected).
///
/// ```
/// use valenx_immuno::aa::aa_index;
/// assert_eq!(aa_index(b'A'), Some(0));
/// assert_eq!(aa_index(b'y'), Some(19));
/// assert_eq!(aa_index(b'X'), None);
/// ```
pub fn aa_index(residue: u8) -> Option<usize> {
    let upper = residue.to_ascii_uppercase();
    ALPHABET.iter().position(|&a| a == upper)
}
