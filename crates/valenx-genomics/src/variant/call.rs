//! Pileup-based SNV and short-indel variant caller (v1).
//!
//! Given a [`crate::format::pileup::PileupColumn`] stream
//! built from aligned reads, this caller proposes variant sites the
//! way VarScan / bcftools-`call` do at heart: tally the alleles at each
//! position, pick the strongest non-reference allele, and decide
//! whether the evidence clears a depth + allele-fraction + quality
//! threshold. A passing site is genotyped with the diploid
//! likelihood model in [`crate::variant::genotype`].
//!
//! Two variant classes are emitted:
//!
//! - **SNV / MNV** — a non-reference *base* with enough supporting
//!   reads and allele fraction at a column.
//! - **Short indel** — an insertion (read bases attached to a pileup
//!   base) or a deletion (`*` placeholders) whose count clears the
//!   thresholds.
//!
//! ## v1 scope
//!
//! This is a real per-site model — pileup tally, threshold, Bayesian
//! genotype likelihood — but **not** GATK HaplotypeCaller's local
//! de-novo assembly of candidate haplotypes, nor a joint multi-sample
//! model. Each column is genotyped independently; indel
//! representation is the strongest single indel allele per anchor
//! position. Strand bias is computed and exposed (see
//! [`crate::variant::filter`]) but not used as a hard gate here.

use crate::error::{GenomicsError, Result};
use crate::format::pileup::PileupColumn;
use crate::variant::genotype::{default_priors, genotype_site, AlleleObs, Genotype, GenotypeCall};
use std::collections::HashMap;

/// The kind of a called variant.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum VariantKind {
    /// A single-nucleotide variant — REF base → ALT base.
    Snv,
    /// An insertion — ALT carries inserted bases after the anchor.
    Insertion,
    /// A deletion — REF carries the deleted bases, ALT is the anchor.
    Deletion,
}

/// One called variant site.
#[derive(Clone, Debug, PartialEq)]
pub struct Variant {
    /// The reference contig.
    pub chrom: String,
    /// The 1-based reference position (the anchor for an indel — the
    /// base **before** the inserted / deleted bases, VCF-style).
    pub pos: i64,
    /// The reference allele.
    pub reference: String,
    /// The alternate allele.
    pub alt: String,
    /// The variant class.
    pub kind: VariantKind,
    /// Total read depth at the site.
    pub depth: usize,
    /// Number of reads supporting the ALT allele.
    pub alt_count: usize,
    /// ALT allele fraction (`alt_count / depth`).
    pub alt_fraction: f64,
    /// Phred-scaled site quality (derived from the genotype call).
    pub qual: f64,
    /// The diploid genotype call.
    pub genotype: GenotypeCall,
    /// Forward / reverse ALT-supporting read counts (strand bias
    /// evidence).
    pub strand: (usize, usize),
}

impl Variant {
    /// `true` for a pure SNV.
    pub fn is_snv(&self) -> bool {
        self.kind == VariantKind::Snv
    }

    /// `true` for an insertion or deletion.
    pub fn is_indel(&self) -> bool {
        matches!(self.kind, VariantKind::Insertion | VariantKind::Deletion)
    }
}

/// Tunable thresholds for the caller.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct CallParams {
    /// Minimum total depth for a site to be considered.
    pub min_depth: usize,
    /// Minimum reads supporting the ALT allele.
    pub min_alt_count: usize,
    /// Minimum ALT allele fraction.
    pub min_alt_fraction: f64,
    /// Minimum mean base quality of the ALT-supporting reads.
    pub min_base_quality: f64,
    /// Minimum phred-scaled site quality to emit a call.
    pub min_qual: f64,
    /// Genotype priors (see [`crate::variant::genotype`]).
    pub priors: [f64; 3],
}

impl Default for CallParams {
    /// Conservative defaults suitable for ~30× germline data.
    fn default() -> Self {
        CallParams {
            min_depth: 8,
            min_alt_count: 3,
            min_alt_fraction: 0.20,
            min_base_quality: 15.0,
            min_qual: 20.0,
            priors: default_priors(),
        }
    }
}

/// Calls variants from a slice of pileup columns.
///
/// Columns are processed independently. For each column the routine
/// finds the strongest SNV allele and the strongest indel allele,
/// genotypes whichever clears the [`CallParams`] gates, and emits a
/// [`Variant`]. Returns the calls sorted by `(chrom, pos)`.
pub fn call_variants(columns: &[PileupColumn], params: &CallParams) -> Result<Vec<Variant>> {
    if params.min_depth == 0 {
        return Err(GenomicsError::invalid("min_depth", "must be positive"));
    }
    if !(0.0..=1.0).contains(&params.min_alt_fraction) {
        return Err(GenomicsError::invalid(
            "min_alt_fraction",
            "must be in [0, 1]",
        ));
    }
    let mut out = Vec::new();
    for col in columns {
        if let Some(v) = call_snv_at(col, params) {
            out.push(v);
        }
        if let Some(v) = call_indel_at(col, params) {
            out.push(v);
        }
    }
    out.sort_by(|a, b| (&a.chrom, a.pos).cmp(&(&b.chrom, b.pos)));
    Ok(out)
}

/// Attempts an SNV call at one column.
fn call_snv_at(col: &PileupColumn, params: &CallParams) -> Option<Variant> {
    let ref_base = col.ref_base.to_ascii_uppercase();
    if !matches!(ref_base, b'A' | b'C' | b'G' | b'T') {
        return None; // an N reference base cannot anchor an SNV
    }
    let base_depth = col.base_depth();
    if base_depth < params.min_depth {
        return None;
    }

    // Tally non-deletion, non-reference bases.
    let mut alt_counts: HashMap<u8, (usize, f64, usize, usize)> = HashMap::new();
    // value = (count, quality sum, fwd, rev)
    for b in &col.bases {
        if b.is_deletion() {
            continue;
        }
        let base = b.base.to_ascii_uppercase();
        if base == ref_base || !matches!(base, b'A' | b'C' | b'G' | b'T') {
            continue;
        }
        let e = alt_counts.entry(base).or_insert((0, 0.0, 0, 0));
        e.0 += 1;
        e.1 += b.quality as f64;
        if b.reverse {
            e.3 += 1;
        } else {
            e.2 += 1;
        }
    }

    // Pick the strongest ALT base.
    let (&alt_base, &(alt_count, qual_sum, fwd, rev)) =
        alt_counts.iter().max_by_key(|(_, v)| v.0)?;
    if alt_count < params.min_alt_count {
        return None;
    }
    let alt_fraction = alt_count as f64 / base_depth as f64;
    if alt_fraction < params.min_alt_fraction {
        return None;
    }
    let mean_alt_q = qual_sum / alt_count as f64;
    if mean_alt_q < params.min_base_quality {
        return None;
    }

    // Genotype: each non-deletion base is ref/alt for this allele.
    let obs: Vec<AlleleObs> = col
        .bases
        .iter()
        .filter(|b| !b.is_deletion())
        .filter_map(|b| {
            let base = b.base.to_ascii_uppercase();
            if base == ref_base {
                Some(AlleleObs {
                    is_ref: true,
                    quality: b.quality,
                })
            } else if base == alt_base {
                Some(AlleleObs {
                    is_ref: false,
                    quality: b.quality,
                })
            } else {
                None // a third allele — ignore for this biallelic call
            }
        })
        .collect();
    let gt = genotype_site(&obs, params.priors);
    if gt.best == Genotype::HomRef {
        return None;
    }
    let qual = site_qual(&gt);
    if qual < params.min_qual {
        return None;
    }

    Some(Variant {
        chrom: col.chrom.clone(),
        pos: col.pos,
        reference: (ref_base as char).to_string(),
        alt: (alt_base as char).to_string(),
        kind: VariantKind::Snv,
        depth: base_depth,
        alt_count,
        alt_fraction,
        qual,
        genotype: gt,
        strand: (fwd, rev),
    })
}

/// Attempts a short-indel call at one column.
fn call_indel_at(col: &PileupColumn, params: &CallParams) -> Option<Variant> {
    let depth = col.depth();
    if depth < params.min_depth {
        return None;
    }

    // Insertion tally: keyed by the inserted-base string.
    let mut ins: HashMap<Vec<u8>, (usize, usize, usize)> = HashMap::new();
    // Deletion: a `*` placeholder means this column is *inside* a
    // deletion; the anchor is one base earlier, so we instead count
    // insertions here and deletions are detected at the anchor base by
    // looking at whether the next position holds `*`. To keep the v1
    // single-column, we treat the deletion as anchored *at this
    // column* when the reads carry a `*` and call REF = ref_base + the
    // deleted run is unknown — so we represent a 1bp deletion. A full
    // multi-base deletion needs the reference; see the module note.
    let mut del_count = 0usize;
    let mut del_fwd = 0usize;
    let mut del_rev = 0usize;
    for b in &col.bases {
        if b.is_deletion() {
            del_count += 1;
            if b.reverse {
                del_rev += 1;
            } else {
                del_fwd += 1;
            }
        }
        if !b.insertion.is_empty() {
            let e = ins.entry(b.insertion.clone()).or_insert((0, 0, 0));
            e.0 += 1;
            if b.reverse {
                e.2 += 1;
            } else {
                e.1 += 1;
            }
        }
    }

    // Strongest insertion vs the deletion — pick the higher count.
    let best_ins = ins.iter().max_by_key(|(_, v)| v.0);
    let ins_count = best_ins.map(|(_, v)| v.0).unwrap_or(0);

    let ref_base = col.ref_base.to_ascii_uppercase();
    let anchor = if matches!(ref_base, b'A' | b'C' | b'G' | b'T') {
        ref_base as char
    } else {
        'N'
    };

    if ins_count >= del_count && ins_count >= params.min_alt_count {
        let (seq, &(count, fwd, rev)) = best_ins.unwrap();
        let frac = count as f64 / depth as f64;
        if frac < params.min_alt_fraction {
            return None;
        }
        // VCF insertion: REF = anchor, ALT = anchor + inserted bases.
        let alt: String = std::iter::once(anchor)
            .chain(seq.iter().map(|&b| b as char))
            .collect();
        let obs = indel_obs(col, count, depth);
        let gt = genotype_site(&obs, params.priors);
        if gt.best == Genotype::HomRef {
            return None;
        }
        let qual = site_qual(&gt);
        if qual < params.min_qual {
            return None;
        }
        Some(Variant {
            chrom: col.chrom.clone(),
            pos: col.pos,
            reference: anchor.to_string(),
            alt,
            kind: VariantKind::Insertion,
            depth,
            alt_count: count,
            alt_fraction: frac,
            qual,
            genotype: gt,
            strand: (fwd, rev),
        })
    } else if del_count >= params.min_alt_count {
        let frac = del_count as f64 / depth as f64;
        if frac < params.min_alt_fraction {
            return None;
        }
        // VCF 1bp deletion: REF = anchor + deleted base (here unknown
        // beyond a single base), ALT = anchor. We anchor one base
        // before this column; without the reference we model a 1bp
        // deletion of `N`.
        let reference = format!("{anchor}N");
        let obs = indel_obs(col, del_count, depth);
        let gt = genotype_site(&obs, params.priors);
        if gt.best == Genotype::HomRef {
            return None;
        }
        let qual = site_qual(&gt);
        if qual < params.min_qual {
            return None;
        }
        Some(Variant {
            chrom: col.chrom.clone(),
            // anchor is the base before this column
            pos: (col.pos - 1).max(1),
            reference,
            alt: anchor.to_string(),
            kind: VariantKind::Deletion,
            depth,
            alt_count: del_count,
            alt_fraction: frac,
            qual,
            genotype: gt,
            strand: (del_fwd, del_rev),
        })
    } else {
        None
    }
}

/// Builds genotype observations for an indel: `alt_count` ALT-supporting
/// observations and `depth - alt_count` REF observations, all at a
/// nominal indel quality (indels carry no per-base Phred score).
fn indel_obs(_col: &PileupColumn, alt_count: usize, depth: usize) -> Vec<AlleleObs> {
    const INDEL_Q: u8 = 30;
    let mut v = Vec::with_capacity(depth);
    for _ in 0..alt_count {
        v.push(AlleleObs {
            is_ref: false,
            quality: INDEL_Q,
        });
    }
    for _ in 0..depth.saturating_sub(alt_count) {
        v.push(AlleleObs {
            is_ref: true,
            quality: INDEL_Q,
        });
    }
    v
}

/// Site QUAL from a genotype call: the phred-scaled probability that
/// the site is *not* homozygous-reference.
fn site_qual(gt: &GenotypeCall) -> f64 {
    // P(variant) = 1 - P(0/0); QUAL = -10 log10(P(0/0)).
    let p_homref = 10f64.powf(gt.log10_posteriors[0]).clamp(0.0, 1.0);
    if p_homref <= 0.0 {
        return 99.0;
    }
    (-10.0 * p_homref.log10()).clamp(0.0, 99.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::format::pileup::{PileupBase, PileupColumn};

    fn pbase(base: u8, q: u8, rev: bool, ins: &[u8]) -> PileupBase {
        PileupBase {
            base,
            quality: q,
            reverse: rev,
            mapq: 60,
            read_pos: 0,
            insertion: ins.to_vec(),
        }
    }

    fn column(ref_base: u8, bases: Vec<PileupBase>) -> PileupColumn {
        PileupColumn {
            chrom: "chr1".to_string(),
            pos: 100,
            ref_base,
            bases,
        }
    }

    #[test]
    fn het_snv_called() {
        // 10 ref A, 10 alt G at quality 35.
        let mut bases = Vec::new();
        for _ in 0..10 {
            bases.push(pbase(b'A', 35, false, &[]));
        }
        for _ in 0..10 {
            bases.push(pbase(b'G', 35, false, &[]));
        }
        let col = column(b'A', bases);
        let vars = call_variants(&[col], &CallParams::default()).unwrap();
        assert_eq!(vars.len(), 1);
        assert_eq!(vars[0].kind, VariantKind::Snv);
        assert_eq!(vars[0].reference, "A");
        assert_eq!(vars[0].alt, "G");
        assert_eq!(vars[0].genotype.best, Genotype::Het);
    }

    #[test]
    fn hom_alt_snv_called() {
        let mut bases = Vec::new();
        for _ in 0..20 {
            bases.push(pbase(b'T', 35, false, &[]));
        }
        let col = column(b'C', bases);
        let vars = call_variants(&[col], &CallParams::default()).unwrap();
        assert_eq!(vars.len(), 1);
        assert_eq!(vars[0].genotype.best, Genotype::HomAlt);
        assert_eq!(vars[0].alt, "T");
    }

    #[test]
    fn pure_reference_makes_no_call() {
        let mut bases = Vec::new();
        for _ in 0..30 {
            bases.push(pbase(b'A', 35, false, &[]));
        }
        let col = column(b'A', bases);
        let vars = call_variants(&[col], &CallParams::default()).unwrap();
        assert!(vars.is_empty());
    }

    #[test]
    fn low_depth_skipped() {
        let mut bases = Vec::new();
        for _ in 0..2 {
            bases.push(pbase(b'G', 35, false, &[]));
        }
        let col = column(b'A', bases);
        let vars = call_variants(&[col], &CallParams::default()).unwrap();
        assert!(vars.is_empty());
    }

    #[test]
    fn low_allele_fraction_skipped() {
        // 1 alt out of 30 -> ~3% AF, well below 20%.
        let mut bases = Vec::new();
        for _ in 0..29 {
            bases.push(pbase(b'A', 35, false, &[]));
        }
        bases.push(pbase(b'G', 35, false, &[]));
        let col = column(b'A', bases);
        let vars = call_variants(&[col], &CallParams::default()).unwrap();
        assert!(vars.is_empty());
    }

    #[test]
    fn insertion_called() {
        // 12 reads carry an insertion of "GG"; 8 are plain reference.
        let mut bases = Vec::new();
        for _ in 0..12 {
            bases.push(pbase(b'A', 35, false, b"GG"));
        }
        for _ in 0..8 {
            bases.push(pbase(b'A', 35, false, &[]));
        }
        let col = column(b'A', bases);
        let vars = call_variants(&[col], &CallParams::default()).unwrap();
        let ins: Vec<_> = vars.iter().filter(|v| v.is_indel()).collect();
        assert_eq!(ins.len(), 1);
        assert_eq!(ins[0].kind, VariantKind::Insertion);
        assert_eq!(ins[0].reference, "A");
        assert_eq!(ins[0].alt, "AGG");
    }

    #[test]
    fn deletion_called() {
        // 14 reads show a `*` deletion placeholder; 6 reference.
        let mut bases = Vec::new();
        for _ in 0..14 {
            bases.push(pbase(b'*', 0, false, &[]));
        }
        for _ in 0..6 {
            bases.push(pbase(b'A', 35, false, &[]));
        }
        let col = column(b'A', bases);
        let vars = call_variants(&[col], &CallParams::default()).unwrap();
        let del: Vec<_> = vars
            .iter()
            .filter(|v| v.kind == VariantKind::Deletion)
            .collect();
        assert_eq!(del.len(), 1);
        assert_eq!(del[0].alt, "A");
    }

    #[test]
    fn strand_counts_recorded() {
        let mut bases = Vec::new();
        for _ in 0..10 {
            bases.push(pbase(b'A', 35, false, &[]));
        }
        for i in 0..10 {
            bases.push(pbase(b'G', 35, i % 2 == 0, &[]));
        }
        let col = column(b'A', bases);
        let vars = call_variants(&[col], &CallParams::default()).unwrap();
        assert_eq!(vars[0].strand, (5, 5));
    }

    #[test]
    fn rejects_bad_params() {
        let p = CallParams {
            min_depth: 0,
            ..CallParams::default()
        };
        assert!(call_variants(&[], &p).is_err());
    }
}
