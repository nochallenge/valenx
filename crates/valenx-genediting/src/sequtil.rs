//! Small DNA / RNA sequence helpers shared across the editing modules.
//!
//! Every gene-editing and mRNA module needs the same handful of
//! primitive operations — complement, reverse complement, transcription
//! to RNA, GC content, ACGT validation. They live here once rather
//! than being re-derived in `crispr`, `base_edit`, `prime_edit`,
//! `mrna` and `therapy`.
//!
//! These are intentionally tiny and `pub(crate)`: the public crate
//! surface exposes the *workflow* types, not a second sequence
//! library — [`valenx_bioseq`] already owns that role.

/// Watson-Crick complement of a single DNA base (`A↔T`, `C↔G`).
/// Any non-ACGT byte (including IUPAC ambiguity codes) maps to `N`.
pub(crate) fn complement(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' | b'U' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => b'N',
    }
}

/// Reverse complement of a DNA slice, uppercased.
pub(crate) fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

/// Transcribes a DNA slice to RNA: uppercase, `T → U`. Non-ACGTU bytes
/// pass through uppercased (callers that need strict validation check
/// [`is_acgt`] first).
pub(crate) fn transcribe(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| {
            let u = b.to_ascii_uppercase();
            if u == b'T' {
                b'U'
            } else {
                u
            }
        })
        .collect()
}

/// Reverse-transcribes an RNA slice back to DNA: uppercase, `U → T`.
pub(crate) fn reverse_transcribe(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| {
            let u = b.to_ascii_uppercase();
            if u == b'U' {
                b'T'
            } else {
                u
            }
        })
        .collect()
}

/// `true` when every byte of `seq` is an unambiguous DNA base
/// (`A C G T`, case-insensitive). An empty slice is `false`.
pub(crate) fn is_acgt(seq: &[u8]) -> bool {
    !seq.is_empty()
        && seq
            .iter()
            .all(|&b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
}

/// `true` when every byte of `seq` is an unambiguous RNA base
/// (`A C G U`, case-insensitive). An empty slice is `false`.
pub(crate) fn is_acgu(seq: &[u8]) -> bool {
    !seq.is_empty()
        && seq
            .iter()
            .all(|&b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'U'))
}

/// Uppercases a slice into an owned `Vec`.
pub(crate) fn upper(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|b| b.to_ascii_uppercase()).collect()
}

/// Longest run of identical bases anywhere in `seq` (homopolymer
/// length). An empty slice has run length `0`.
pub(crate) fn max_homopolymer(seq: &[u8]) -> usize {
    let mut best = 0usize;
    let mut run = 0usize;
    let mut prev = 0u8;
    for &b in seq {
        let u = b.to_ascii_uppercase();
        if u == prev {
            run += 1;
        } else {
            run = 1;
            prev = u;
        }
        best = best.max(run);
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complement_and_revcomp() {
        assert_eq!(complement(b'A'), b'T');
        assert_eq!(complement(b'g'), b'C');
        assert_eq!(complement(b'N'), b'N');
        assert_eq!(revcomp(b"ATGC"), b"GCAT");
    }

    #[test]
    fn transcription_round_trip() {
        let dna = b"ATGCTA";
        let rna = transcribe(dna);
        assert_eq!(rna, b"AUGCUA");
        assert_eq!(reverse_transcribe(&rna), dna);
    }

    #[test]
    fn validation_helpers() {
        assert!(is_acgt(b"ACGT"));
        assert!(!is_acgt(b"ACGN"));
        assert!(!is_acgt(b""));
        assert!(is_acgu(b"ACGU"));
        assert!(!is_acgu(b"ACGT"));
    }

    #[test]
    fn homopolymer_runs() {
        assert_eq!(max_homopolymer(b"AACCCGT"), 3);
        assert_eq!(max_homopolymer(b""), 0);
        assert_eq!(max_homopolymer(b"ATTTTACG"), 4);
    }
}
