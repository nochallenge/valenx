//! Translation: codon → amino acid via the NCBI genetic-code tables.
//!
//! A [`GeneticCode`] is built from the NCBI "Genetic Codes" data
//! (`AAs` / `Starts` strings). All published NCBI translation tables
//! (1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 16, 21, 22, 23, 24, 25,
//! 26, 27, 28, 29, 30, 31, 33) ship built-in. Tables 7, 8, 15, 17, 18,
//! 19, 20, 32 are reserved / withdrawn in the NCBI source and are
//! intentionally absent.
//!
//! Source: NCBI Taxonomy "The Genetic Codes",
//! `ftp.ncbi.nlm.nih.gov/entrez/misc/data/gc.prt`,
//! <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>.
//!
//! Translation handles `*` stop codons, partial trailing codons
//! (dropped), and ambiguity codons: an ambiguous codon translates to a
//! definite residue only if every concrete codon it expands to yields
//! the same amino acid, otherwise `X`.

use crate::alphabet::{self, SeqKind};
use crate::error::{BioseqError, Result};
use crate::seq::{Seq, Topology};

/// Codon order used by every NCBI `AAs` / `Starts` string: the third
/// base cycles fastest, with base order `T C A G`.
const NCBI_BASES: [u8; 4] = [b'T', b'C', b'A', b'G'];

/// A genetic code — an amino-acid assignment for all 64 codons plus a
/// set of start codons.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GeneticCode {
    /// NCBI translation-table id (1 = Standard).
    pub id: u8,
    /// Human-readable name.
    pub name: &'static str,
    /// `aa[i]` is the amino acid for codon `i` in NCBI order; `*` =
    /// stop.
    aa: [u8; 64],
    /// `is_start[i]` marks codon `i` as a valid start codon.
    is_start: [bool; 64],
}

/// Maps a base to its NCBI index (`T C A G` → `0 1 2 3`). `U` folds to
/// `T`. Returns `None` for non-canonical bases.
fn base_index(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'T' | b'U' => Some(0),
        b'C' => Some(1),
        b'A' => Some(2),
        b'G' => Some(3),
        _ => None,
    }
}

/// Computes the NCBI codon index `0..64` for three concrete bases.
fn codon_index(c0: u8, c1: u8, c2: u8) -> Option<usize> {
    Some(base_index(c0)? * 16 + base_index(c1)? * 4 + base_index(c2)?)
}

impl GeneticCode {
    /// Builds a code from the NCBI `aas` and `starts` strings (each 64
    /// chars, in NCBI codon order). `starts` uses `M` for a start
    /// codon and `-` otherwise.
    fn from_ncbi(id: u8, name: &'static str, aas: &str, starts: &str) -> GeneticCode {
        let ab = aas.as_bytes();
        let sb = starts.as_bytes();
        assert_eq!(ab.len(), 64, "NCBI aas string must be 64 chars");
        assert_eq!(sb.len(), 64, "NCBI starts string must be 64 chars");
        let mut aa = [0u8; 64];
        let mut is_start = [false; 64];
        for i in 0..64 {
            aa[i] = ab[i];
            is_start[i] = sb[i] == b'M';
        }
        GeneticCode {
            id,
            name,
            aa,
            is_start,
        }
    }

    /// Returns the built-in code for an NCBI table id, or
    /// [`BioseqError::Invalid`] for an unsupported id.
    ///
    /// All published NCBI tables are encoded:
    /// 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 16, 21, 22, 23, 24, 25,
    /// 26, 27, 28, 29, 30, 31, 33. Tables 7, 8, 15, 17, 18, 19, 20, 32
    /// are reserved / withdrawn in the NCBI source and therefore not
    /// defined.
    pub fn by_id(id: u8) -> Result<GeneticCode> {
        // The NCBI codon order for these strings is third-base-fastest
        // with base order TCAG. Source: NCBI Taxonomy "The Genetic
        // Codes", https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
        // and ftp.ncbi.nlm.nih.gov/entrez/misc/data/gc.prt.
        let code = match id {
            1 => GeneticCode::from_ncbi(
                1,
                "Standard",
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "---M------**--*----M---------------M----------------------------",
            ),
            2 => GeneticCode::from_ncbi(
                2,
                "Vertebrate Mitochondrial",
                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
                "----------**--------------------MMMM----------**----------------",
            ),
            3 => GeneticCode::from_ncbi(
                3,
                "Yeast Mitochondrial",
                "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "----------**----------------------MM---------------M------------",
            ),
            4 => GeneticCode::from_ncbi(
                4,
                "Mold/Protozoan/Coelenterate Mitochondrial; Mycoplasma/Spiroplasma",
                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "--MM------**-------M------------MMMM---------------M------------",
            ),
            5 => GeneticCode::from_ncbi(
                5,
                "Invertebrate Mitochondrial",
                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
                "---M------**--------------------MMMM---------------M------------",
            ),
            6 => GeneticCode::from_ncbi(
                6,
                "Ciliate, Dasycladacean and Hexamita Nuclear",
                "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "--------------*--------------------M----------------------------",
            ),
            9 => GeneticCode::from_ncbi(
                9,
                "Echinoderm and Flatworm Mitochondrial",
                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
                "----------**-----------------------M---------------M------------",
            ),
            10 => GeneticCode::from_ncbi(
                10,
                "Euplotid Nuclear",
                "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "----------**-----------------------M----------------------------",
            ),
            11 => GeneticCode::from_ncbi(
                11,
                "Bacterial, Archaeal and Plant Plastid",
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "---M------**--*----M------------MMMM---------------M------------",
            ),
            12 => GeneticCode::from_ncbi(
                12,
                "Alternative Yeast Nuclear",
                "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "----------**--*----M---------------M----------------------------",
            ),
            13 => GeneticCode::from_ncbi(
                13,
                "Ascidian Mitochondrial",
                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
                "---M------**----------------------MM---------------M------------",
            ),
            14 => GeneticCode::from_ncbi(
                14,
                "Alternative Flatworm Mitochondrial",
                "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
                "-----------*-----------------------M----------------------------",
            ),
            16 => GeneticCode::from_ncbi(
                16,
                "Chlorophycean Mitochondrial",
                "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "----------*---*--------------------M----------------------------",
            ),
            21 => GeneticCode::from_ncbi(
                21,
                "Trematode Mitochondrial",
                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
                "----------**-----------------------M---------------M------------",
            ),
            22 => GeneticCode::from_ncbi(
                22,
                "Scenedesmus obliquus Mitochondrial",
                "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "------*---*---*--------------------M----------------------------",
            ),
            23 => GeneticCode::from_ncbi(
                23,
                "Thraustochytrium Mitochondrial",
                "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "--*-------**--*-----------------M--M---------------M------------",
            ),
            24 => GeneticCode::from_ncbi(
                24,
                "Pterobranchia Mitochondrial",
                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
                "---M------**-------M---------------M---------------M------------",
            ),
            25 => GeneticCode::from_ncbi(
                25,
                "Candidate Division SR1 and Gracilibacteria",
                "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "---M------**-----------------------M---------------M------------",
            ),
            26 => GeneticCode::from_ncbi(
                26,
                "Pachysolen tannophilus Nuclear",
                "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "----------**--*----M---------------M----------------------------",
            ),
            27 => GeneticCode::from_ncbi(
                27,
                "Karyorelict Nuclear",
                "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "--------------*--------------------M----------------------------",
            ),
            28 => GeneticCode::from_ncbi(
                28,
                "Condylostoma Nuclear",
                "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "----------**--*--------------------M----------------------------",
            ),
            29 => GeneticCode::from_ncbi(
                29,
                "Mesodinium Nuclear",
                "FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "--------------*--------------------M----------------------------",
            ),
            30 => GeneticCode::from_ncbi(
                30,
                "Peritrich Nuclear",
                "FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "--------------*--------------------M----------------------------",
            ),
            31 => GeneticCode::from_ncbi(
                31,
                "Blastocrithidia Nuclear",
                "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "----------**-----------------------M----------------------------",
            ),
            33 => GeneticCode::from_ncbi(
                33,
                "Cephalodiscidae Mitochondrial",
                "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSKVVVVAAAADDEEGGGG",
                "---M-------*-------M---------------M---------------M------------",
            ),
            other => {
                return Err(BioseqError::invalid(
                    "genetic_code",
                    format!(
                        "unsupported NCBI table {other}; supported: 1,2,3,4,5,6,9,10,11,12,13,14,16,21,22,23,24,25,26,27,28,29,30,31,33"
                    ),
                ))
            }
        };
        Ok(code)
    }

    /// The standard (NCBI table 1) genetic code.
    pub fn standard() -> GeneticCode {
        GeneticCode::by_id(1).expect("table 1 is built-in")
    }

    /// The list of NCBI table ids this crate ships, in ascending order.
    pub fn supported_ids() -> &'static [u8] {
        &[
            1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
            31, 33,
        ]
    }

    /// Translates one codon (three bytes) to an amino acid.
    ///
    /// For concrete codons this is a table lookup. For codons
    /// containing IUPAC ambiguity codes, every concrete expansion is
    /// translated; the result is the common amino acid if unique, else
    /// `X`. A codon with a non-nucleotide byte yields `X`.
    pub fn translate_codon(&self, codon: &[u8]) -> u8 {
        if codon.len() != 3 {
            return b'X';
        }
        // Fast path: fully concrete codon.
        if let Some(i) = codon_index(codon[0], codon[1], codon[2]) {
            return self.aa[i];
        }
        // Ambiguous path: expand each position.
        let e0 = match alphabet::expand_iupac(codon[0]) {
            Some(s) => s,
            None => return b'X',
        };
        let e1 = match alphabet::expand_iupac(codon[1]) {
            Some(s) => s,
            None => return b'X',
        };
        let e2 = match alphabet::expand_iupac(codon[2]) {
            Some(s) => s,
            None => return b'X',
        };
        let mut seen: Option<u8> = None;
        for &b0 in e0 {
            for &b1 in e1 {
                for &b2 in e2 {
                    let i = codon_index(b0, b1, b2).expect("expanded bases are concrete");
                    let a = self.aa[i];
                    match seen {
                        None => seen = Some(a),
                        Some(prev) if prev != a => return b'X',
                        Some(_) => {}
                    }
                }
            }
        }
        seen.unwrap_or(b'X')
    }

    /// `true` if `codon` is a start codon under this code.
    pub fn is_start_codon(&self, codon: &[u8]) -> bool {
        if codon.len() != 3 {
            return false;
        }
        match codon_index(codon[0], codon[1], codon[2]) {
            Some(i) => self.is_start[i],
            None => false,
        }
    }

    /// `true` if `codon` is a stop codon under this code.
    pub fn is_stop_codon(&self, codon: &[u8]) -> bool {
        self.translate_codon(codon) == b'*'
    }

    /// All start codons in this code, as 3-byte arrays.
    pub fn start_codons(&self) -> Vec<[u8; 3]> {
        let mut out = Vec::new();
        for (i, &start) in self.is_start.iter().enumerate() {
            if start {
                let b2 = NCBI_BASES[i % 4];
                let b1 = NCBI_BASES[(i / 4) % 4];
                let b0 = NCBI_BASES[(i / 16) % 4];
                out.push([b0, b1, b2]);
            }
        }
        out
    }
}

/// Options controlling [`translate`].
#[derive(Copy, Clone, Debug, Default)]
pub struct TranslateOptions {
    /// Stop translating at the first stop codon (Biopython
    /// `to_stop=True`). Default `false`.
    pub to_stop: bool,
    /// Translate the first codon as `M` if it is a start codon under
    /// the code (Biopython `cds`-style behavior). Default `false`.
    pub start_as_met: bool,
}

/// Translates a nucleotide sequence into a protein [`Seq`].
///
/// Reads codons left to right starting at index 0; a partial trailing
/// codon (1 or 2 leftover bases) is dropped. RNA is accepted directly;
/// DNA is translated as if transcribed. Returns
/// [`BioseqError::Invalid`] for a protein input.
pub fn translate(seq: &Seq, code: &GeneticCode, opts: TranslateOptions) -> Result<Seq> {
    if seq.kind() == SeqKind::Protein {
        return Err(BioseqError::invalid(
            "kind",
            "cannot translate a protein sequence",
        ));
    }
    let bytes = seq.as_bytes();
    let mut protein: Vec<u8> = Vec::with_capacity(bytes.len() / 3);
    let n_codons = bytes.len() / 3;
    for c in 0..n_codons {
        let codon = &bytes[c * 3..c * 3 + 3];
        let mut aa = code.translate_codon(codon);
        if c == 0 && opts.start_as_met && code.is_start_codon(codon) {
            aa = b'M';
        }
        if aa == b'*' && opts.to_stop {
            break;
        }
        protein.push(aa);
    }
    Ok(Seq::new_unchecked(
        SeqKind::Protein,
        protein,
        Topology::Linear,
    ))
}

/// Translates with default options (full read-through, no
/// start-as-Met).
pub fn translate_default(seq: &Seq, code: &GeneticCode) -> Result<Seq> {
    translate(seq, code, TranslateOptions::default())
}

/// One reading frame's translation result.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FrameTranslation {
    /// Frame index: `0..3` forward (offset 0/1/2), `3..6` reverse.
    pub frame: u8,
    /// `true` for the three reverse-strand frames.
    pub reverse: bool,
    /// Nucleotide offset (0, 1 or 2) the frame starts at.
    pub offset: usize,
    /// The translated protein.
    pub protein: Seq,
}

/// Six-frame translation — three forward frames (offsets 0/1/2) and
/// three reverse frames (the reverse complement at offsets 0/1/2).
///
/// Returns [`BioseqError::Invalid`] for a protein input.
pub fn six_frame_translation(seq: &Seq, code: &GeneticCode) -> Result<Vec<FrameTranslation>> {
    if seq.kind() == SeqKind::Protein {
        return Err(BioseqError::invalid(
            "kind",
            "six-frame translation needs a nucleotide sequence",
        ));
    }
    use crate::ops::revcomp::reverse_complement;
    let mut frames = Vec::with_capacity(6);
    for offset in 0..3 {
        let sub = seq.slice(offset.min(seq.len()), seq.len())?;
        let protein = translate_default(&sub, code)?;
        frames.push(FrameTranslation {
            frame: offset as u8,
            reverse: false,
            offset,
            protein,
        });
    }
    let rc = reverse_complement(seq)?;
    for offset in 0..3 {
        let sub = rc.slice(offset.min(rc.len()), rc.len())?;
        let protein = translate_default(&sub, code)?;
        frames.push(FrameTranslation {
            frame: 3 + offset as u8,
            reverse: true,
            offset,
            protein,
        });
    }
    Ok(frames)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn standard_code_known_codons() {
        let c = GeneticCode::standard();
        assert_eq!(c.translate_codon(b"ATG"), b'M');
        assert_eq!(c.translate_codon(b"TTT"), b'F');
        assert_eq!(c.translate_codon(b"GGG"), b'G');
        assert_eq!(c.translate_codon(b"TAA"), b'*');
        assert_eq!(c.translate_codon(b"TAG"), b'*');
        assert_eq!(c.translate_codon(b"TGA"), b'*');
        assert_eq!(c.translate_codon(b"TGG"), b'W');
    }

    #[test]
    fn rna_codons_translate() {
        let c = GeneticCode::standard();
        assert_eq!(c.translate_codon(b"AUG"), b'M');
        assert_eq!(c.translate_codon(b"UUU"), b'F');
    }

    #[test]
    fn vertebrate_mito_differences() {
        // Table 2: AGA/AGG are stops, ATA is Met, TGA is Trp.
        let c = GeneticCode::by_id(2).unwrap();
        assert_eq!(c.translate_codon(b"AGA"), b'*');
        assert_eq!(c.translate_codon(b"AGG"), b'*');
        assert_eq!(c.translate_codon(b"ATA"), b'M');
        assert_eq!(c.translate_codon(b"TGA"), b'W');
    }

    #[test]
    fn ciliate_code_differences() {
        // Table 6: TAA and TAG code for Gln, not stop.
        let c = GeneticCode::by_id(6).unwrap();
        assert_eq!(c.translate_codon(b"TAA"), b'Q');
        assert_eq!(c.translate_codon(b"TAG"), b'Q');
        assert_eq!(c.translate_codon(b"TGA"), b'*'); // still stop
    }

    #[test]
    fn bacterial_code_extra_starts() {
        let c = GeneticCode::by_id(11).unwrap();
        assert!(c.is_start_codon(b"ATG"));
        assert!(c.is_start_codon(b"GTG"));
        assert!(c.is_start_codon(b"TTG"));
        assert!(!c.is_start_codon(b"AAA"));
    }

    #[test]
    fn unsupported_table_errs() {
        assert!(GeneticCode::by_id(99).is_err());
        // Withdrawn / reserved tables are also reported as unsupported.
        assert!(GeneticCode::by_id(7).is_err());
        assert!(GeneticCode::by_id(8).is_err());
        assert!(GeneticCode::by_id(15).is_err());
        assert!(GeneticCode::by_id(17).is_err());
        assert!(GeneticCode::by_id(18).is_err());
        assert!(GeneticCode::by_id(19).is_err());
        assert!(GeneticCode::by_id(20).is_err());
        assert!(GeneticCode::by_id(32).is_err());
    }

    #[test]
    fn translate_orf() {
        let c = GeneticCode::standard();
        // ATG AAA GGG TAA  -> M K G *
        let dna = Seq::new(SeqKind::Dna, "ATGAAAGGGTAA").unwrap();
        assert_eq!(translate_default(&dna, &c).unwrap().as_str(), "MKG*");
    }

    #[test]
    fn translate_to_stop() {
        let c = GeneticCode::standard();
        let dna = Seq::new(SeqKind::Dna, "ATGAAATAAGGG").unwrap();
        let opts = TranslateOptions {
            to_stop: true,
            ..Default::default()
        };
        assert_eq!(translate(&dna, &c, opts).unwrap().as_str(), "MK");
    }

    #[test]
    fn partial_trailing_codon_dropped() {
        let c = GeneticCode::standard();
        let dna = Seq::new(SeqKind::Dna, "ATGAA").unwrap(); // 5 nt
        assert_eq!(translate_default(&dna, &c).unwrap().as_str(), "M");
    }

    #[test]
    fn ambiguous_codon_resolves_or_x() {
        let c = GeneticCode::standard();
        // CTN — all four expansions are Leu -> definite L.
        assert_eq!(c.translate_codon(b"CTN"), b'L');
        // ATN — ATT/ATC/ATA = I, ATG = M -> ambiguous -> X.
        assert_eq!(c.translate_codon(b"ATN"), b'X');
        // GGN — all Gly.
        assert_eq!(c.translate_codon(b"GGN"), b'G');
    }

    #[test]
    fn start_as_met_option() {
        let c = GeneticCode::by_id(11).unwrap();
        // GTG normally translates to Val; as a start it becomes M.
        let dna = Seq::new(SeqKind::Dna, "GTGAAA").unwrap();
        let opts = TranslateOptions {
            start_as_met: true,
            ..Default::default()
        };
        assert_eq!(translate(&dna, &c, opts).unwrap().as_str(), "MK");
        // Without the option it stays V.
        assert_eq!(translate_default(&dna, &c).unwrap().as_str(), "VK");
    }

    #[test]
    fn six_frames_count_and_shape() {
        let c = GeneticCode::standard();
        let dna = Seq::new(SeqKind::Dna, "ATGAAAGGGTTTCCC").unwrap();
        let frames = six_frame_translation(&dna, &c).unwrap();
        assert_eq!(frames.len(), 6);
        assert_eq!(frames[0].protein.as_str(), "MKGFP");
        assert!(!frames[0].reverse);
        assert!(frames[3].reverse);
    }

    #[test]
    fn protein_input_rejected() {
        let c = GeneticCode::standard();
        let p = Seq::new(SeqKind::Protein, "MKVL").unwrap();
        assert!(translate_default(&p, &c).is_err());
        assert!(six_frame_translation(&p, &c).is_err());
    }

    // -------------------------------------------------------------------
    // All-tables spot-checks: the diagnostic codons for each NCBI table.
    // These are the exact divergences from the standard code listed in
    // the NCBI gc.prt file — i.e. the reason each non-standard table
    // exists. Asserting them all catches any AAs-string typo.
    // -------------------------------------------------------------------

    #[test]
    fn table_1_standard_landmarks() {
        let c = GeneticCode::by_id(1).unwrap();
        assert_eq!(c.translate_codon(b"ATG"), b'M');
        assert!(c.is_start_codon(b"ATG"));
        assert!(c.is_start_codon(b"CTG")); // alternative start in T1
        assert!(c.is_start_codon(b"TTG"));
        assert_eq!(c.translate_codon(b"TGA"), b'*');
    }

    #[test]
    fn table_3_yeast_mito_landmarks() {
        // CTN -> Thr (the four leucine CTN codons become T in yeast mito).
        let c = GeneticCode::by_id(3).unwrap();
        for codon in [b"CTT", b"CTC", b"CTA", b"CTG"] {
            assert_eq!(
                c.translate_codon(codon),
                b'T',
                "{}",
                std::str::from_utf8(codon).unwrap()
            );
        }
        assert_eq!(c.translate_codon(b"TGA"), b'W'); // UGA -> Trp
        assert_eq!(c.translate_codon(b"ATA"), b'M'); // AUA -> Met
    }

    #[test]
    fn table_4_mold_mito_landmarks() {
        // TGA -> Trp; alternative start TTA / TTG / CTG / GTG / ATT / ATC / ATA.
        let c = GeneticCode::by_id(4).unwrap();
        assert_eq!(c.translate_codon(b"TGA"), b'W');
        assert!(c.is_start_codon(b"TTA"));
        assert!(c.is_start_codon(b"TTG"));
        assert!(c.is_start_codon(b"CTG"));
        assert!(c.is_start_codon(b"GTG"));
    }

    #[test]
    fn table_5_invertebrate_mito_landmarks() {
        let c = GeneticCode::by_id(5).unwrap();
        assert_eq!(c.translate_codon(b"AGA"), b'S'); // S, not R
        assert_eq!(c.translate_codon(b"AGG"), b'S');
        assert_eq!(c.translate_codon(b"ATA"), b'M');
        assert_eq!(c.translate_codon(b"TGA"), b'W');
    }

    #[test]
    fn table_9_echinoderm_mito_landmarks() {
        let c = GeneticCode::by_id(9).unwrap();
        assert_eq!(c.translate_codon(b"AAA"), b'N'); // N, not K
        assert_eq!(c.translate_codon(b"AGA"), b'S');
        assert_eq!(c.translate_codon(b"AGG"), b'S');
        assert_eq!(c.translate_codon(b"TGA"), b'W');
    }

    #[test]
    fn table_10_euplotid_landmarks() {
        let c = GeneticCode::by_id(10).unwrap();
        assert_eq!(c.translate_codon(b"TGA"), b'C'); // TGA -> Cys
    }

    #[test]
    fn table_12_alt_yeast_landmarks() {
        let c = GeneticCode::by_id(12).unwrap();
        assert_eq!(c.translate_codon(b"CTG"), b'S'); // CTG -> Ser
        assert!(c.is_start_codon(b"CTG"));
    }

    #[test]
    fn table_13_ascidian_mito_landmarks() {
        let c = GeneticCode::by_id(13).unwrap();
        assert_eq!(c.translate_codon(b"AGA"), b'G'); // AGA -> Gly
        assert_eq!(c.translate_codon(b"AGG"), b'G');
        assert_eq!(c.translate_codon(b"ATA"), b'M');
        assert_eq!(c.translate_codon(b"TGA"), b'W');
    }

    #[test]
    fn table_14_alt_flatworm_landmarks() {
        let c = GeneticCode::by_id(14).unwrap();
        assert_eq!(c.translate_codon(b"AAA"), b'N');
        assert_eq!(c.translate_codon(b"AGA"), b'S');
        assert_eq!(c.translate_codon(b"AGG"), b'S');
        assert_eq!(c.translate_codon(b"TAA"), b'Y'); // UAA -> Tyr in T14
        assert_eq!(c.translate_codon(b"TGA"), b'W');
    }

    #[test]
    fn table_16_chlorophycean_landmarks() {
        let c = GeneticCode::by_id(16).unwrap();
        assert_eq!(c.translate_codon(b"TAG"), b'L'); // TAG -> Leu
    }

    #[test]
    fn table_21_trematode_mito_landmarks() {
        let c = GeneticCode::by_id(21).unwrap();
        assert_eq!(c.translate_codon(b"TGA"), b'W');
        assert_eq!(c.translate_codon(b"ATA"), b'M');
        assert_eq!(c.translate_codon(b"AGA"), b'S');
        assert_eq!(c.translate_codon(b"AAA"), b'N');
    }

    #[test]
    fn table_22_scenedesmus_landmarks() {
        let c = GeneticCode::by_id(22).unwrap();
        assert_eq!(c.translate_codon(b"TCA"), b'*'); // TCA -> stop
        assert_eq!(c.translate_codon(b"TAG"), b'L'); // TAG -> Leu
    }

    #[test]
    fn table_23_thraustochytrium_landmarks() {
        let c = GeneticCode::by_id(23).unwrap();
        assert_eq!(c.translate_codon(b"TTA"), b'*'); // TTA -> stop
    }

    #[test]
    fn table_24_pterobranchia_landmarks() {
        let c = GeneticCode::by_id(24).unwrap();
        assert_eq!(c.translate_codon(b"AGA"), b'S');
        assert_eq!(c.translate_codon(b"AGG"), b'K'); // AGG -> Lys
        assert_eq!(c.translate_codon(b"TGA"), b'W');
    }

    #[test]
    fn table_25_sr1_landmarks() {
        let c = GeneticCode::by_id(25).unwrap();
        assert_eq!(c.translate_codon(b"TGA"), b'G'); // TGA -> Gly
    }

    #[test]
    fn table_33_cephalodiscidae_landmarks() {
        let c = GeneticCode::by_id(33).unwrap();
        // T33 keeps T9-style AAA->N, AGA->S; adds AGG->K, TGA->W; and
        // (the diagnostic of T33) TAA -> Tyr.
        assert_eq!(c.translate_codon(b"TAA"), b'Y');
        assert_eq!(c.translate_codon(b"TGA"), b'W');
        assert_eq!(c.translate_codon(b"AGA"), b'S');
        assert_eq!(c.translate_codon(b"AGG"), b'K');
    }

    #[test]
    fn supported_ids_round_trip_all_tables() {
        // Every supported id must actually load.
        for &id in GeneticCode::supported_ids() {
            let c = GeneticCode::by_id(id).unwrap();
            assert_eq!(c.id, id);
            assert!(!c.name.is_empty());
            // At least one start codon is defined.
            assert!(!c.start_codons().is_empty(), "table {id} has no starts");
        }
        // We ship 25 tables (all non-withdrawn NCBI tables).
        assert_eq!(GeneticCode::supported_ids().len(), 25);
    }
}
