//! VCF normalisation — left-align indels, decompose multiallelics, trim.
//!
//! The same variant can be written many ways. `bcftools norm` and
//! `vt normalize` put every record into a single canonical form so two
//! callers' output can be compared. This module implements the three
//! core transforms:
//!
//! 1. **Decompose** ([`decompose_multiallelic`]) — split a record with
//!    *N* ALT alleles into *N* biallelic records.
//! 2. **Trim** ([`trim_alleles`]) — remove shared leading and trailing
//!    bases common to REF and ALT, keeping one anchor base when the
//!    trim would otherwise empty an allele.
//! 3. **Left-align** ([`left_align`]) — shift an indel as far left
//!    (5′) as the reference allows, the way `vt normalize` does, so a
//!    deletion in a homopolymer run always reports the same position.
//!
//! [`normalize_record`] runs all three; [`normalize_vcf`] applies it
//! across a whole file. A [`crate::format::pileup::Reference`] supplies
//! the reference bases needed for left-alignment.

use crate::format::pileup::Reference;
use crate::format::vcf::{VcfFile, VcfRecord};

/// Splits a multiallelic record into one biallelic record per ALT
/// allele. A record with 0 or 1 ALT alleles is returned unchanged
/// (wrapped in a one-element `Vec`).
///
/// INFO and per-sample fields are copied verbatim to every split
/// record — a `bcftools norm -m-` style decomposition that does **not**
/// re-apportion `Number=A` INFO fields (a documented v1 simplification;
/// re-apportionment needs the `##INFO` Number metadata).
pub fn decompose_multiallelic(rec: &VcfRecord) -> Vec<VcfRecord> {
    if rec.alt.len() <= 1 {
        return vec![rec.clone()];
    }
    rec.alt
        .iter()
        .map(|alt| {
            let mut r = rec.clone();
            r.alt = vec![alt.clone()];
            r
        })
        .collect()
}

/// Trims a single REF/ALT allele pair.
///
/// Shared trailing bases are removed first, then shared leading bases;
/// `pos` is advanced by the number of leading bases removed. At least
/// one base is always kept in each allele (the anchor). Returns the
/// new `(pos, ref, alt)`.
pub fn trim_alleles(pos: i64, reference: &str, alt: &str) -> (i64, String, String) {
    let mut r: Vec<u8> = reference.bytes().collect();
    let mut a: Vec<u8> = alt.bytes().collect();
    let mut new_pos = pos;

    // Trim shared suffix while both keep > 1 base.
    while r.len() > 1 && a.len() > 1 && r.last() == a.last() {
        r.pop();
        a.pop();
    }
    // Trim shared prefix while both keep > 1 base.
    while r.len() > 1 && a.len() > 1 && r.first() == a.first() {
        r.remove(0);
        a.remove(0);
        new_pos += 1;
    }
    (
        new_pos,
        String::from_utf8_lossy(&r).into_owned(),
        String::from_utf8_lossy(&a).into_owned(),
    )
}

/// Left-aligns a biallelic indel against the reference.
///
/// For a pure insertion or deletion (one allele a strict prefix-anchor
/// of the other after trimming), the routine repeatedly shifts the
/// variant one base 5′ while the base it would roll over equals the
/// last base of the indel — the standard left-alignment loop. SNVs and
/// MNVs (equal-length alleles) are returned unchanged. `contig` /
/// `reference` provide the bases; `pos` is 1-based.
///
/// Returns the left-aligned `(pos, ref, alt)`.
pub fn left_align(
    contig: &str,
    reference_seq: &Reference,
    pos: i64,
    reference: &str,
    alt: &str,
) -> (i64, String, String) {
    // Only indels (length difference) can be shifted.
    if reference.len() == alt.len() {
        return (pos, reference.to_string(), alt.to_string());
    }
    let mut r: Vec<u8> = reference.bytes().map(|b| b.to_ascii_uppercase()).collect();
    let mut a: Vec<u8> = alt.bytes().map(|b| b.to_ascii_uppercase()).collect();
    let mut p = pos;

    // Shift left while there is room and the rolled-over base matches.
    // Guard against an unbounded loop on a degenerate reference.
    let mut guard = 0usize;
    while p > 1 && guard < 10_000 {
        guard += 1;
        // The reference base immediately 5′ of the current anchor.
        let prev = reference_seq.base_at(contig, (p - 2) as usize);
        if prev == b'N' {
            break;
        }
        // Both alleles must currently end with the same base for a
        // shift to be valid (the indel is "rotatable").
        if r.last() != a.last() {
            break;
        }
        // Roll: drop the shared last base, prepend `prev` to both.
        r.pop();
        a.pop();
        r.insert(0, prev);
        a.insert(0, prev);
        p -= 1;
    }
    (
        p,
        String::from_utf8_lossy(&r).into_owned(),
        String::from_utf8_lossy(&a).into_owned(),
    )
}

/// Fully normalises one VCF record: decompose → trim → left-align each
/// resulting biallelic record. Returns one record per ALT allele.
pub fn normalize_record(rec: &VcfRecord, reference: &Reference) -> Vec<VcfRecord> {
    let mut out = Vec::new();
    for biallelic in decompose_multiallelic(rec) {
        let alt = biallelic.alt.first().cloned().unwrap_or_default();
        // Trim first.
        let (p1, r1, a1) = trim_alleles(biallelic.pos, &biallelic.reference, &alt);
        // Left-align (only shifts indels; uses the reference if present).
        let (p2, r2, a2) = if reference.has(&biallelic.chrom) {
            let (lp, lr, la) = left_align(&biallelic.chrom, reference, p1, &r1, &a1);
            // Re-trim after the shift in case it exposed a shared base.
            trim_alleles(lp, &lr, &la)
        } else {
            (p1, r1, a1)
        };
        let mut r = biallelic;
        r.pos = p2;
        r.reference = r2;
        r.alt = vec![a2];
        out.push(r);
    }
    out
}

/// Normalises every record of a VCF file, returning a new file. The
/// record count grows when multiallelic sites are decomposed. Records
/// are re-sorted by `(chrom, pos)`.
pub fn normalize_vcf(vcf: &VcfFile, reference: &Reference) -> VcfFile {
    let mut records: Vec<VcfRecord> = Vec::new();
    for rec in &vcf.records {
        records.extend(normalize_record(rec, reference));
    }
    records.sort_by(|a, b| (&a.chrom, a.pos).cmp(&(&b.chrom, b.pos)));
    VcfFile {
        header: vcf.header.clone(),
        records,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decompose_splits_alt_alleles() {
        let mut rec = VcfRecord::snv("chr1", 100, "C", "A");
        rec.alt = vec!["A".to_string(), "G".to_string(), "T".to_string()];
        let split = decompose_multiallelic(&rec);
        assert_eq!(split.len(), 3);
        assert_eq!(split[0].alt, vec!["A"]);
        assert_eq!(split[1].alt, vec!["G"]);
        assert_eq!(split[2].alt, vec!["T"]);
    }

    #[test]
    fn decompose_leaves_biallelic_alone() {
        let rec = VcfRecord::snv("chr1", 100, "C", "A");
        assert_eq!(decompose_multiallelic(&rec).len(), 1);
    }

    #[test]
    fn trim_shared_suffix() {
        // REF=ATCC ALT=AGCC: trailing CC trims, then the leading A
        // trims too — `trim_alleles` removes shared bases at *both*
        // ends (the minimal `vt normalize` representation). What is
        // left is the genuine SNV T->G, shifted to pos 101.
        let (pos, r, a) = trim_alleles(100, "ATCC", "AGCC");
        assert_eq!(pos, 101);
        assert_eq!(r, "T");
        assert_eq!(a, "G");
    }

    #[test]
    fn trim_shared_prefix_advances_pos() {
        // REF=GGA ALT=GGT -> leading GG trimmed -> A / T, pos +2.
        let (pos, r, a) = trim_alleles(100, "GGA", "GGT");
        assert_eq!(pos, 102);
        assert_eq!(r, "A");
        assert_eq!(a, "T");
    }

    #[test]
    fn trim_keeps_one_anchor_base() {
        // A simple deletion REF=AT ALT=A — already minimal.
        let (pos, r, a) = trim_alleles(100, "AT", "A");
        assert_eq!((pos, r.as_str(), a.as_str()), (100, "AT", "A"));
    }

    #[test]
    fn trim_indel_with_redundant_padding() {
        // REF=CTT ALT=CT is a 1bp deletion; trailing T shared.
        let (pos, r, a) = trim_alleles(100, "CTT", "CT");
        assert_eq!(r, "CT");
        assert_eq!(a, "C");
        assert_eq!(pos, 100);
    }

    #[test]
    fn left_align_shifts_deletion_in_homopolymer() {
        // Reference: positions 1.. = "GAAAAAT".
        // A deletion of one A written at pos 5 (REF="AA" ALT="A")
        // left-aligns past the whole A-run and anchors on the G at
        // pos 1 — the standard `vt`/`bcftools` result is REF="GA"
        // ALT="G" at pos 1, since the anchor base must be the base
        // immediately 5' of the homopolymer.
        let mut refr = Reference::new();
        refr.add("chr1", "GAAAAAT");
        let (pos, r, a) = left_align("chr1", &refr, 5, "AA", "A");
        assert_eq!(pos, 1, "expected shift to the anchor 5' of the A-run");
        assert_eq!(r, "GA");
        assert_eq!(a, "G");
    }

    #[test]
    fn left_align_leaves_snv_unchanged() {
        let mut refr = Reference::new();
        refr.add("chr1", "ACGTACGT");
        let (pos, r, a) = left_align("chr1", &refr, 4, "T", "G");
        assert_eq!((pos, r.as_str(), a.as_str()), (4, "T", "G"));
    }

    #[test]
    fn normalize_record_decomposes_and_trims() {
        let mut refr = Reference::new();
        refr.add("chr1", "GAAAAAT");
        let mut rec = VcfRecord::snv("chr1", 5, "AA", "A");
        rec.alt = vec!["A".to_string(), "AAA".to_string()];
        let norm = normalize_record(&rec, &refr);
        // Two ALT alleles -> two records.
        assert_eq!(norm.len(), 2);
        // The deletion allele (AA -> A) left-aligns past the A-run and
        // anchors on the G at pos 1.
        assert_eq!(norm[0].pos, 1);
    }

    #[test]
    fn normalize_vcf_grows_on_multiallelic() {
        let header = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t100\t.\tC\tA,G\t.\tPASS\t.\n";
        let vcf = VcfFile::parse(header).unwrap();
        let norm = normalize_vcf(&vcf, &Reference::new());
        assert_eq!(norm.records.len(), 2);
    }
}
