//! Variant filtering and basic annotation.
//!
//! After calling, a variant set is filtered against quality gates and
//! annotated with derived metrics — exactly what GATK
//! `VariantFiltration` and bcftools `filter` do. This module operates
//! on the [`crate::variant::call::Variant`] values produced by
//! [`crate::variant::call`] (or on [`VcfRecord`]s) and stamps each
//! with a [`FilterStatus`] plus a list of failed-filter tags.
//!
//! The hard gates are depth, site quality, allele fraction, genotype
//! quality and a **strand-bias** test (a symmetric Fisher-style
//! odds-ratio proxy — a real bias signal without pulling a stats
//! crate).

use crate::format::vcf::VcfRecord;
use crate::variant::call::Variant;

/// Whether a variant passed every filter.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum FilterStatus {
    /// Cleared every gate.
    Pass,
    /// Failed one or more gates (see the tag list).
    Fail,
}

/// Quality gates applied to each variant.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct VariantFilter {
    /// Minimum total depth.
    pub min_depth: usize,
    /// Maximum total depth (guards against collapsed-repeat pileups;
    /// `None` disables).
    pub max_depth: Option<usize>,
    /// Minimum phred-scaled site quality.
    pub min_qual: f64,
    /// Minimum ALT allele fraction.
    pub min_alt_fraction: f64,
    /// Minimum genotype quality (`GQ`).
    pub min_gq: u8,
    /// Maximum tolerated strand-bias score (see [`strand_bias_score`]);
    /// `0` is perfectly balanced, larger is worse.
    pub max_strand_bias: f64,
}

impl Default for VariantFilter {
    /// Sensible germline defaults.
    fn default() -> Self {
        VariantFilter {
            min_depth: 8,
            max_depth: None,
            min_qual: 20.0,
            min_alt_fraction: 0.20,
            min_gq: 20,
            max_strand_bias: 3.0,
        }
    }
}

/// A variant annotated with its filter verdict.
#[derive(Clone, Debug, PartialEq)]
pub struct AnnotatedVariant {
    /// The variant.
    pub variant: Variant,
    /// The overall verdict.
    pub status: FilterStatus,
    /// The tags of every failed gate (empty when [`FilterStatus::Pass`]).
    pub failed: Vec<&'static str>,
    /// The computed strand-bias score.
    pub strand_bias: f64,
}

impl AnnotatedVariant {
    /// `true` when the variant passed every filter.
    pub fn passed(&self) -> bool {
        self.status == FilterStatus::Pass
    }

    /// The VCF `FILTER` field value: `"PASS"` or the failed tags joined
    /// by `;`.
    pub fn filter_field(&self) -> String {
        if self.failed.is_empty() {
            "PASS".to_string()
        } else {
            self.failed.join(";")
        }
    }
}

/// A symmetric strand-bias score.
///
/// Given the forward/reverse counts of the reference allele
/// `(ref_fwd, ref_rev)` and the alternate allele `(alt_fwd, alt_rev)`,
/// the score is the absolute log of the odds ratio
///
/// ```text
/// |ln( (alt_fwd · ref_rev) / (alt_rev · ref_fwd) )|
/// ```
///
/// with a `+0.5` Haldane-Anscombe correction on every cell so a zero
/// count cannot blow the ratio up. A balanced site scores near `0`; an
/// allele seen on only one strand scores high.
pub fn strand_bias_score(ref_strand: (usize, usize), alt_strand: (usize, usize)) -> f64 {
    let a = alt_strand.0 as f64 + 0.5;
    let b = alt_strand.1 as f64 + 0.5;
    let c = ref_strand.0 as f64 + 0.5;
    let d = ref_strand.1 as f64 + 0.5;
    let odds = (a * d) / (b * c);
    odds.ln().abs()
}

/// Filters and annotates a slice of called variants.
///
/// Each variant is tested against the [`VariantFilter`] gates; the
/// strand-bias score uses the variant's recorded ALT strand counts and
/// derives the REF strand counts from `depth - alt_count` split evenly
/// (the caller does not retain per-strand REF counts in v1, so this is
/// a conservative even-split assumption — documented honestly).
pub fn filter_variants(variants: &[Variant], filter: &VariantFilter) -> Vec<AnnotatedVariant> {
    variants.iter().map(|v| annotate_one(v, filter)).collect()
}

fn annotate_one(v: &Variant, filter: &VariantFilter) -> AnnotatedVariant {
    let mut failed: Vec<&'static str> = Vec::new();

    if v.depth < filter.min_depth {
        failed.push("LowDepth");
    }
    if let Some(max) = filter.max_depth {
        if v.depth > max {
            failed.push("HighDepth");
        }
    }
    if v.qual < filter.min_qual {
        failed.push("LowQual");
    }
    if v.alt_fraction < filter.min_alt_fraction {
        failed.push("LowAltFraction");
    }
    if v.genotype.gq < filter.min_gq {
        failed.push("LowGQ");
    }

    // REF strand counts: even-split of the non-ALT depth (v1 caveat).
    let ref_total = v.depth.saturating_sub(v.alt_count);
    let ref_strand = (ref_total / 2, ref_total - ref_total / 2);
    let sb = strand_bias_score(ref_strand, v.strand);
    if sb > filter.max_strand_bias {
        failed.push("StrandBias");
    }

    let status = if failed.is_empty() {
        FilterStatus::Pass
    } else {
        FilterStatus::Fail
    };
    AnnotatedVariant {
        variant: v.clone(),
        status,
        failed,
        strand_bias: sb,
    }
}

/// Converts an [`AnnotatedVariant`] into a [`VcfRecord`] — REF / ALT,
/// QUAL, the FILTER field, an `INFO` block (`DP`, `AF`, `SB`) and a
/// single-sample `GT:GQ:DP` genotype column.
pub fn variant_to_vcf(av: &AnnotatedVariant, sample: &str) -> VcfRecord {
    let v = &av.variant;
    let mut rec = VcfRecord::snv(&v.chrom, v.pos, &v.reference, &v.alt);
    rec.qual = Some(v.qual);
    rec.filter = if av.failed.is_empty() {
        vec!["PASS".to_string()]
    } else {
        av.failed.iter().map(|s| s.to_string()).collect()
    };
    rec.info.insert("DP".to_string(), v.depth.to_string());
    rec.info
        .insert("AF".to_string(), format!("{:.4}", v.alt_fraction));
    rec.info
        .insert("SB".to_string(), format!("{:.3}", av.strand_bias));
    if v.is_indel() {
        rec.info.insert("INDEL".to_string(), String::new());
    }
    if !sample.is_empty() {
        rec.format = vec!["GT".to_string(), "GQ".to_string(), "DP".to_string()];
        rec.samples = vec![vec![
            v.genotype.best.gt_string().to_string(),
            v.genotype.gq.to_string(),
            v.depth.to_string(),
        ]];
    }
    rec
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::call::VariantKind;
    use crate::variant::genotype::{default_priors, genotype_site, AlleleObs};

    fn variant(depth: usize, alt: usize, qual: f64, strand: (usize, usize)) -> Variant {
        let obs: Vec<AlleleObs> = (0..depth)
            .map(|i| AlleleObs {
                is_ref: i >= alt,
                quality: 35,
            })
            .collect();
        let gt = genotype_site(&obs, default_priors());
        Variant {
            chrom: "chr1".to_string(),
            pos: 100,
            reference: "A".to_string(),
            alt: "G".to_string(),
            kind: VariantKind::Snv,
            depth,
            alt_count: alt,
            alt_fraction: alt as f64 / depth as f64,
            qual,
            genotype: gt,
            strand,
        }
    }

    #[test]
    fn strand_bias_balanced_is_low() {
        let sb = strand_bias_score((20, 20), (20, 20));
        assert!(sb < 0.1, "sb = {sb}");
    }

    #[test]
    fn strand_bias_one_sided_is_high() {
        // ALT only on the forward strand.
        let sb = strand_bias_score((20, 20), (40, 0));
        assert!(sb > 2.0, "sb = {sb}");
    }

    #[test]
    fn passing_variant() {
        let v = variant(30, 15, 60.0, (8, 7));
        let av = annotate_one(&v, &VariantFilter::default());
        assert!(av.passed());
        assert_eq!(av.filter_field(), "PASS");
    }

    #[test]
    fn low_depth_fails() {
        let v = variant(4, 2, 60.0, (1, 1));
        let av = annotate_one(&v, &VariantFilter::default());
        assert!(!av.passed());
        assert!(av.failed.contains(&"LowDepth"));
    }

    #[test]
    fn low_qual_fails() {
        let v = variant(30, 15, 5.0, (8, 7));
        let av = annotate_one(&v, &VariantFilter::default());
        assert!(av.failed.contains(&"LowQual"));
    }

    #[test]
    fn strand_bias_fails() {
        // ALT entirely on the forward strand.
        let v = variant(30, 15, 60.0, (15, 0));
        let av = annotate_one(&v, &VariantFilter::default());
        assert!(av.failed.contains(&"StrandBias"));
    }

    #[test]
    fn high_depth_gate() {
        let f = VariantFilter {
            max_depth: Some(100),
            ..VariantFilter::default()
        };
        let v = variant(500, 250, 60.0, (125, 125));
        let av = annotate_one(&v, &f);
        assert!(av.failed.contains(&"HighDepth"));
    }

    #[test]
    fn vcf_conversion_carries_info() {
        let v = variant(30, 15, 60.0, (8, 7));
        let av = annotate_one(&v, &VariantFilter::default());
        let rec = variant_to_vcf(&av, "sampleA");
        assert_eq!(rec.info_get("DP"), Some("30"));
        assert!(rec.info_get("AF").is_some());
        assert_eq!(
            rec.sample_field(0, "GT"),
            Some(av.variant.genotype.best.gt_string())
        );
        assert_eq!(rec.filter, vec!["PASS"]);
    }

    #[test]
    fn filter_variants_batch() {
        let vs = vec![
            variant(30, 15, 60.0, (8, 7)), // pass
            variant(3, 2, 60.0, (1, 1)),   // low depth
        ];
        let out = filter_variants(&vs, &VariantFilter::default());
        assert_eq!(out.len(), 2);
        assert!(out[0].passed());
        assert!(!out[1].passed());
    }
}
