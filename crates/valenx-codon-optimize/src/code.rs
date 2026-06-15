//! The standard genetic code (NCBI translation table 1).

/// Amino acids for all 64 codons, indexed in `TCAG × TCAG × TCAG` base order
/// (the canonical NCBI table-1 encoding). `*` denotes a stop codon.
const AAS: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

/// Bases in the canonical index order (`T=0, C=1, A=2, G=3`).
const BASES: [u8; 4] = [b'T', b'C', b'A', b'G'];

fn base_index(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'T' | b'U' => Some(0),
        b'C' => Some(1),
        b'A' => Some(2),
        b'G' => Some(3),
        _ => None,
    }
}

/// Translate a 3-letter DNA codon to its amino acid (`*` for a stop codon), or
/// `None` if the codon is not exactly three valid bases (`T/C/A/G`, with `U`
/// accepted as `T`).
pub fn translate_codon(codon: &str) -> Option<char> {
    let b = codon.as_bytes();
    if b.len() != 3 {
        return None;
    }
    let idx = base_index(b[0])? * 16 + base_index(b[1])? * 4 + base_index(b[2])?;
    Some(AAS[idx] as char)
}

/// All synonymous codons for an amino acid (use `*` for stop), in canonical
/// order. Empty if `aa` is not a coded residue.
pub fn synonymous_codons(aa: char) -> Vec<String> {
    let target = aa.to_ascii_uppercase() as u8;
    let mut out = Vec::new();
    for (idx, &a) in AAS.iter().enumerate() {
        if a == target {
            let codon: String = [BASES[idx / 16], BASES[(idx / 4) % 4], BASES[idx % 4]]
                .iter()
                .map(|&b| b as char)
                .collect();
            out.push(codon);
        }
    }
    out
}

/// Whether `aa` codes for exactly one codon (Met, Trp) — excluded from CAI.
pub fn is_single_codon(aa: char) -> bool {
    matches!(aa.to_ascii_uppercase(), 'M' | 'W')
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_codons_translate() {
        assert_eq!(translate_codon("ATG"), Some('M'));
        assert_eq!(translate_codon("TGG"), Some('W'));
        assert_eq!(translate_codon("TTT"), Some('F'));
        assert_eq!(translate_codon("GGG"), Some('G'));
        assert_eq!(translate_codon("TAA"), Some('*'));
        assert_eq!(translate_codon("TAG"), Some('*'));
        assert_eq!(translate_codon("TGA"), Some('*'));
        assert_eq!(translate_codon("AUG"), Some('M')); // U accepted as T
    }

    #[test]
    fn invalid_codons_are_none() {
        assert_eq!(translate_codon("AT"), None);
        assert_eq!(translate_codon("ATGC"), None);
        assert_eq!(translate_codon("ATX"), None);
    }

    #[test]
    fn synonymous_counts_match_the_code() {
        assert_eq!(synonymous_codons('M'), vec!["ATG"]);
        assert_eq!(synonymous_codons('W'), vec!["TGG"]);
        assert_eq!(synonymous_codons('L').len(), 6); // Leu: 6 codons
        assert_eq!(synonymous_codons('R').len(), 6); // Arg: 6 codons
        assert_eq!(synonymous_codons('S').len(), 6); // Ser: 6 codons
        assert_eq!(synonymous_codons('*').len(), 3); // 3 stop codons
                                                     // every synonymous codon really translates back to the amino acid
        for c in synonymous_codons('V') {
            assert_eq!(translate_codon(&c), Some('V'));
        }
    }

    #[test]
    fn all_61_sense_codons_plus_3_stops() {
        let mut sense = 0;
        let mut stop = 0;
        for &b1 in &BASES {
            for &b2 in &BASES {
                for &b3 in &BASES {
                    let c: String = [b1, b2, b3].iter().map(|&b| b as char).collect();
                    match translate_codon(&c) {
                        Some('*') => stop += 1,
                        Some(_) => sense += 1,
                        None => unreachable!(),
                    }
                }
            }
        }
        assert_eq!(sense, 61);
        assert_eq!(stop, 3);
    }
}
