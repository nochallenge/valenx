//! Standard part library, annotation and codon-context helpers —
//! feature 29.
//!
//! A synthetic-biology design tool needs a *parts library* — a set of
//! well-characterised promoters, RBSs, CDSs and terminators a designer
//! draws from — and the ability to *annotate* an arbitrary sequence
//! against that library (the SBOLCanvas / pLannotate workflow).
//!
//! This module provides:
//!
//! - [`standard_part_library`] — a small built-in library of canonical
//!   iGEM / Anderson-collection-style parts (the Anderson constitutive
//!   promoters, the B0034 RBS, a GFP CDS, the B0015 terminator).
//! - [`PartLibrary`] — a searchable collection with lookup by id and
//!   by role.
//! - [`annotate_sequence`] — scans a target sequence for exact matches
//!   to any library part (both strands) and returns typed
//!   [`Annotation`]s — feeding the [`crate::synbio::sbol`] data model.
//! - [`codon_adaptation_index`] and [`relative_codon_frequency`] — the
//!   codon-context helpers used when designing a CDS for a given host:
//!   the CAI scores how well a coding sequence matches a host's
//!   preferred codon usage.

use std::collections::HashMap;

use valenx_bioseq::ops::revcomp::reverse_complement_dna_bytes;

use crate::error::{Result, SysbioError};
use crate::synbio::sbol::{Part, PartRole};

/// A searchable collection of genetic parts.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct PartLibrary {
    /// Library name.
    pub name: String,
    /// The parts.
    pub parts: Vec<Part>,
}

impl PartLibrary {
    /// An empty library with the given name.
    pub fn new(name: impl Into<String>) -> Self {
        PartLibrary {
            name: name.into(),
            parts: Vec::new(),
        }
    }

    /// Append a part.
    pub fn add(&mut self, part: Part) {
        self.parts.push(part);
    }

    /// Look up a part by its id.
    pub fn get(&self, id: &str) -> Option<&Part> {
        self.parts.iter().find(|p| p.id == id)
    }

    /// Every part playing the given role.
    pub fn by_role(&self, role: PartRole) -> Vec<&Part> {
        self.parts.iter().filter(|p| p.role == role).collect()
    }

    /// Number of parts in the library.
    pub fn len(&self) -> usize {
        self.parts.len()
    }

    /// Whether the library is empty.
    pub fn is_empty(&self) -> bool {
        self.parts.is_empty()
    }
}

/// Build the built-in standard part library (feature 29).
///
/// A compact, real set of canonical parts: three Anderson-collection
/// constitutive promoters of graded strength, the consensus B0034 RBS,
/// a GFP coding sequence, and the B0015 double terminator. Enough to
/// design a working expression cassette out of the box.
pub fn standard_part_library() -> PartLibrary {
    let mut lib = PartLibrary::new("valenx-standard-parts");
    // Anderson constitutive promoters (real iGEM registry sequences).
    for (id, seq) in [
        ("J23100", "TTGACGGCTAGCTCAGTCCTAGGTACAGTGCTAGC"),
        ("J23106", "TTTACGGCTAGCTCAGTCCTAGGTATAGTGCTAGC"),
        ("J23114", "TTTATGGCTAGCTCAGTCCTAGGTACAATGCTAGC"),
    ] {
        if let Ok(p) = Part::new(id, PartRole::Promoter, seq) {
            lib.add(p);
        }
    }
    // Consensus RBS.
    if let Ok(p) = Part::new("B0034", PartRole::Rbs, "AAAGAGGAGAAA") {
        lib.add(p);
    }
    // GFP coding sequence (truncated representative ORF).
    if let Ok(p) = Part::new(
        "E0040_gfp",
        PartRole::Cds,
        "ATGCGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATTAA",
    ) {
        lib.add(p);
    }
    // B0015 double terminator (representative sequence).
    if let Ok(p) = Part::new(
        "B0015",
        PartRole::Terminator,
        "CCAGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTG",
    ) {
        lib.add(p);
    }
    lib
}

/// A typed match of a library part against a target sequence.
#[derive(Debug, Clone, PartialEq)]
pub struct Annotation {
    /// Id of the matched library part.
    pub part_id: String,
    /// Role of the matched part.
    pub role: PartRole,
    /// 0-based inclusive start in the target.
    pub start: usize,
    /// 0-based exclusive end in the target.
    pub end: usize,
    /// `true` if the match is on the reverse strand.
    pub reverse_strand: bool,
}

/// Annotate `target` (a DNA byte slice) against `library` (feature
/// 29).
///
/// Scans for every exact occurrence of each library part on both
/// strands and returns the matches sorted by start position. Empty
/// parts are skipped. The scan is case-insensitive.
pub fn annotate_sequence(target: &[u8], library: &PartLibrary) -> Vec<Annotation> {
    let hay: Vec<u8> = target.iter().map(|b| b.to_ascii_uppercase()).collect();
    let mut annotations = Vec::new();
    for part in &library.parts {
        if part.is_empty() {
            continue;
        }
        let fwd: Vec<u8> = part
            .sequence
            .as_bytes()
            .iter()
            .map(|b| b.to_ascii_uppercase())
            .collect();
        for start in find_all(&hay, &fwd) {
            annotations.push(Annotation {
                part_id: part.id.clone(),
                role: part.role,
                start,
                end: start + fwd.len(),
                reverse_strand: false,
            });
        }
        let rc = reverse_complement_dna_bytes(&fwd);
        if rc != fwd {
            for start in find_all(&hay, &rc) {
                annotations.push(Annotation {
                    part_id: part.id.clone(),
                    role: part.role,
                    start,
                    end: start + rc.len(),
                    reverse_strand: true,
                });
            }
        }
    }
    annotations.sort_by_key(|a| (a.start, a.end));
    annotations
}

/// Every start offset of `needle` in `haystack` (overlapping matches).
fn find_all(haystack: &[u8], needle: &[u8]) -> Vec<usize> {
    let mut out = Vec::new();
    if needle.is_empty() || haystack.len() < needle.len() {
        return out;
    }
    for i in 0..=haystack.len() - needle.len() {
        if &haystack[i..i + needle.len()] == needle {
            out.push(i);
        }
    }
    out
}

/// The relative codon frequency table of a coding sequence.
///
/// Splits `cds` into codons and returns, for every codon, its
/// frequency *relative to the most-used synonymous codon for the same
/// amino acid* — the `w_ij` weights of the codon-adaptation index.
/// `cds` must be a multiple of three nucleotides long.
pub fn relative_codon_frequency(cds: &[u8]) -> Result<HashMap<String, f64>> {
    if cds.is_empty() || cds.len() % 3 != 0 {
        return Err(SysbioError::invalid(
            "cds",
            "coding sequence length must be a positive multiple of three",
        ));
    }
    // Count each codon.
    let mut counts: HashMap<String, usize> = HashMap::new();
    for chunk in cds.chunks(3) {
        let codon = String::from_utf8_lossy(chunk).to_uppercase();
        *counts.entry(codon).or_insert(0) += 1;
    }
    // Group by amino acid; the per-AA max count normalises the group.
    let mut by_aa: HashMap<char, usize> = HashMap::new();
    for (codon, &c) in &counts {
        let aa = translate_codon(codon);
        let e = by_aa.entry(aa).or_insert(0);
        *e = (*e).max(c);
    }
    let mut w = HashMap::new();
    for (codon, &c) in &counts {
        let aa = translate_codon(codon);
        let max = by_aa[&aa].max(1);
        w.insert(codon.clone(), c as f64 / max as f64);
    }
    Ok(w)
}

/// The codon-adaptation index (CAI) of `cds` against a host weight
/// table (feature 29).
///
/// CAI is the geometric mean of the per-codon relative-adaptiveness
/// weights `w` from the host: `CAI = (∏ w_codon)^{1/L}`. A value near
/// `1.0` means the coding sequence uses the host's preferred codons
/// throughout; a low value flags codons that may translate poorly.
/// Codons absent from the weight table are skipped (treated as
/// neutral) rather than zeroing the whole product.
pub fn codon_adaptation_index(cds: &[u8], host_weights: &HashMap<String, f64>) -> Result<f64> {
    if cds.is_empty() || cds.len() % 3 != 0 {
        return Err(SysbioError::invalid(
            "cds",
            "coding sequence length must be a positive multiple of three",
        ));
    }
    let mut log_sum = 0.0;
    let mut counted = 0usize;
    for chunk in cds.chunks(3) {
        let codon = String::from_utf8_lossy(chunk).to_uppercase();
        if let Some(&w) = host_weights.get(&codon) {
            if w > 0.0 {
                log_sum += w.ln();
                counted += 1;
            }
        }
    }
    if counted == 0 {
        return Err(SysbioError::invalid(
            "cds",
            "no codons matched the host weight table",
        ));
    }
    Ok((log_sum / counted as f64).exp())
}

/// Translate one DNA codon to its amino-acid letter (standard genetic
/// code; `*` for a stop, `X` for an unrecognised codon). A compact
/// table sufficient for the codon-context helpers.
fn translate_codon(codon: &str) -> char {
    match codon {
        "TTT" | "TTC" => 'F',
        "TTA" | "TTG" | "CTT" | "CTC" | "CTA" | "CTG" => 'L',
        "ATT" | "ATC" | "ATA" => 'I',
        "ATG" => 'M',
        "GTT" | "GTC" | "GTA" | "GTG" => 'V',
        "TCT" | "TCC" | "TCA" | "TCG" | "AGT" | "AGC" => 'S',
        "CCT" | "CCC" | "CCA" | "CCG" => 'P',
        "ACT" | "ACC" | "ACA" | "ACG" => 'T',
        "GCT" | "GCC" | "GCA" | "GCG" => 'A',
        "TAT" | "TAC" => 'Y',
        "TAA" | "TAG" | "TGA" => '*',
        "CAT" | "CAC" => 'H',
        "CAA" | "CAG" => 'Q',
        "AAT" | "AAC" => 'N',
        "AAA" | "AAG" => 'K',
        "GAT" | "GAC" => 'D',
        "GAA" | "GAG" => 'E',
        "TGT" | "TGC" => 'C',
        "TGG" => 'W',
        "CGT" | "CGC" | "CGA" | "CGG" | "AGA" | "AGG" => 'R',
        "GGT" | "GGC" | "GGA" | "GGG" => 'G',
        _ => 'X',
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_bioseq::{Seq, SeqKind};

    #[test]
    fn standard_library_has_each_role() {
        let lib = standard_part_library();
        assert!(!lib.by_role(PartRole::Promoter).is_empty());
        assert!(!lib.by_role(PartRole::Rbs).is_empty());
        assert!(!lib.by_role(PartRole::Cds).is_empty());
        assert!(!lib.by_role(PartRole::Terminator).is_empty());
        assert!(lib.get("J23100").is_some());
    }

    #[test]
    fn annotation_finds_a_planted_part() {
        let lib = standard_part_library();
        let promoter = lib.get("J23100").unwrap().sequence.as_bytes().to_vec();
        // Embed the promoter inside a longer sequence.
        let mut target = b"AAAAAAAAAA".to_vec();
        target.extend_from_slice(&promoter);
        target.extend_from_slice(b"TTTTTTTTTT");
        let annos = annotate_sequence(&target, &lib);
        let hit = annos
            .iter()
            .find(|a| a.part_id == "J23100")
            .expect("promoter annotated");
        assert_eq!(hit.start, 10);
        assert_eq!(hit.end, 10 + promoter.len());
        assert!(!hit.reverse_strand);
    }

    #[test]
    fn annotation_finds_reverse_strand_match() {
        let lib = standard_part_library();
        let rbs = lib.get("B0034").unwrap().sequence.as_bytes().to_vec();
        let rc = reverse_complement_dna_bytes(&rbs);
        let mut target = b"GGGGGG".to_vec();
        target.extend_from_slice(&rc);
        let annos = annotate_sequence(&target, &lib);
        let hit = annos
            .iter()
            .find(|a| a.part_id == "B0034")
            .expect("RBS annotated");
        assert!(hit.reverse_strand);
    }

    #[test]
    fn relative_codon_frequency_normalises_per_amino_acid() {
        // Two leucine codons, CTG used 3x and CTA used 1x.
        // CTG -> w 1.0, CTA -> w 1/3.
        let cds = b"CTGCTGCTGCTA";
        let w = relative_codon_frequency(cds).unwrap();
        assert!((w["CTG"] - 1.0).abs() < 1e-12);
        assert!((w["CTA"] - 1.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn cai_is_one_for_perfectly_adapted_cds() {
        // Host strongly prefers CTG, GAA.
        let mut host = HashMap::new();
        host.insert("CTG".to_string(), 1.0);
        host.insert("GAA".to_string(), 1.0);
        // A CDS using only the preferred codons scores CAI = 1.
        let cds = b"CTGGAACTGGAA";
        let cai = codon_adaptation_index(cds, &host).unwrap();
        assert!((cai - 1.0).abs() < 1e-9, "cai {cai}");
    }

    #[test]
    fn cai_drops_for_rare_codon_usage() {
        let mut host = HashMap::new();
        host.insert("CTG".to_string(), 1.0); // preferred Leu
        host.insert("CTA".to_string(), 0.1); // rare Leu
        let good = codon_adaptation_index(b"CTGCTGCTG", &host).unwrap();
        let poor = codon_adaptation_index(b"CTACTACTA", &host).unwrap();
        assert!(good > poor);
        assert!((poor - 0.1).abs() < 1e-9);
    }

    #[test]
    fn cai_rejects_non_triplet_length() {
        let host = HashMap::new();
        assert!(codon_adaptation_index(b"ATGC", &host).is_err());
    }

    #[test]
    fn library_lookup_by_role_and_id() {
        let mut lib = PartLibrary::new("tiny");
        lib.add(Part::new("p1", PartRole::Promoter, "ACGT").unwrap());
        let _ = Seq::new(SeqKind::Dna, "ACGT").unwrap();
        assert_eq!(lib.len(), 1);
        assert!(lib.get("p1").is_some());
        assert_eq!(lib.by_role(PartRole::Promoter).len(), 1);
        assert!(lib.by_role(PartRole::Cds).is_empty());
    }
}
