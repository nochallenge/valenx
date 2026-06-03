//! Allele-frequency and variant-set statistics.
//!
//! Summary statistics over a VCF — the numbers VCFtools `--freq` and
//! bcftools `stats` report: per-site allele frequency from the
//! genotype columns, the transition/transversion ratio, the
//! SNV/indel breakdown, the genotype-class counts and a Hardy-Weinberg
//! observed-vs-expected heterozygosity check.

use crate::format::vcf::{VcfFile, VcfRecord};

/// Per-site allele-frequency result.
#[derive(Clone, Debug, PartialEq)]
pub struct SiteFrequency {
    /// The contig.
    pub chrom: String,
    /// The 1-based position.
    pub pos: i64,
    /// The reference allele.
    pub reference: String,
    /// The alternate alleles.
    pub alt: Vec<String>,
    /// Allele counts; index `0` is REF, `1..` are the ALT alleles.
    pub allele_counts: Vec<usize>,
    /// Total called alleles (`2 × non-missing genotypes`).
    pub total_alleles: usize,
}

impl SiteFrequency {
    /// Allele frequency of allele `index` (`0` = REF). Returns `0.0`
    /// when no alleles were called.
    pub fn frequency(&self, index: usize) -> f64 {
        if self.total_alleles == 0 {
            return 0.0;
        }
        self.allele_counts
            .get(index)
            .map(|&c| c as f64 / self.total_alleles as f64)
            .unwrap_or(0.0)
    }

    /// The frequency of the first ALT allele.
    pub fn alt_frequency(&self) -> f64 {
        self.frequency(1)
    }

    /// The minor-allele frequency — the smaller of the REF and the
    /// summed-ALT frequency.
    pub fn minor_allele_frequency(&self) -> f64 {
        let ref_f = self.frequency(0);
        (ref_f).min(1.0 - ref_f)
    }
}

/// Parses the integer allele indices from a VCF `GT` value such as
/// `0/1`, `1|2` or `./.`. Missing (`.`) alleles are skipped.
fn parse_gt(gt: &str) -> Vec<usize> {
    gt.split(['/', '|'])
        .filter_map(|a| a.parse::<usize>().ok())
        .collect()
}

/// Computes the [`SiteFrequency`] of one VCF record from its genotype
/// columns.
pub fn site_frequency(rec: &VcfRecord) -> SiteFrequency {
    let n_alleles = 1 + rec.alt.len();
    let mut counts = vec![0usize; n_alleles];
    let mut total = 0usize;
    for s_idx in 0..rec.samples.len() {
        if let Some(gt) = rec.sample_field(s_idx, "GT") {
            for a in parse_gt(gt) {
                if a < n_alleles {
                    counts[a] += 1;
                    total += 1;
                }
            }
        }
    }
    SiteFrequency {
        chrom: rec.chrom.clone(),
        pos: rec.pos,
        reference: rec.reference.clone(),
        alt: rec.alt.clone(),
        allele_counts: counts,
        total_alleles: total,
    }
}

/// `true` when REF→ALT is a transition (A↔G or C↔T), `false` for a
/// transversion. Only meaningful for single-base alleles.
pub fn is_transition(reference: &str, alt: &str) -> bool {
    if reference.len() != 1 || alt.len() != 1 {
        return false;
    }
    let r = reference.as_bytes()[0].to_ascii_uppercase();
    let a = alt.as_bytes()[0].to_ascii_uppercase();
    matches!(
        (r, a),
        (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C')
    )
}

/// Genotype-class counts at a biallelic site across all samples.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub struct GenotypeCounts {
    /// Count of homozygous-reference (`0/0`) genotypes.
    pub hom_ref: usize,
    /// Count of heterozygous genotypes.
    pub het: usize,
    /// Count of homozygous-alternate genotypes.
    pub hom_alt: usize,
    /// Count of missing (`./.`) genotypes.
    pub missing: usize,
}

impl GenotypeCounts {
    /// Total non-missing genotypes.
    pub fn called(&self) -> usize {
        self.hom_ref + self.het + self.hom_alt
    }

    /// Observed heterozygosity (`het / called`).
    pub fn observed_heterozygosity(&self) -> f64 {
        let n = self.called();
        if n == 0 {
            0.0
        } else {
            self.het as f64 / n as f64
        }
    }

    /// Hardy-Weinberg **expected** heterozygosity `2pq` from the
    /// observed allele frequencies.
    pub fn expected_heterozygosity(&self) -> f64 {
        let n = self.called();
        if n == 0 {
            return 0.0;
        }
        let p = (2 * self.hom_ref + self.het) as f64 / (2 * n) as f64;
        let q = 1.0 - p;
        2.0 * p * q
    }

    /// The inbreeding coefficient `F = 1 − H_obs / H_exp` — a quick
    /// Hardy-Weinberg deviation signal. Positive `F` means a
    /// heterozygote deficit; `0.0` is HWE. Returns `0.0` when `H_exp`
    /// is zero.
    pub fn inbreeding_coefficient(&self) -> f64 {
        let he = self.expected_heterozygosity();
        if he <= 0.0 {
            0.0
        } else {
            1.0 - self.observed_heterozygosity() / he
        }
    }
}

/// Counts the genotype classes of a biallelic VCF record.
pub fn genotype_counts(rec: &VcfRecord) -> GenotypeCounts {
    let mut c = GenotypeCounts::default();
    for s_idx in 0..rec.samples.len() {
        match rec.sample_field(s_idx, "GT") {
            None => c.missing += 1,
            Some(gt) => {
                let alleles = parse_gt(gt);
                if alleles.len() < 2 {
                    c.missing += 1;
                } else {
                    let a = alleles[0];
                    let b = alleles[1];
                    if a == 0 && b == 0 {
                        c.hom_ref += 1;
                    } else if a == b {
                        c.hom_alt += 1;
                    } else {
                        c.het += 1;
                    }
                }
            }
        }
    }
    c
}

/// A whole-VCF summary.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct VcfStats {
    /// Total records.
    pub total: usize,
    /// SNV records.
    pub snvs: usize,
    /// Indel records.
    pub indels: usize,
    /// Multiallelic records.
    pub multiallelic: usize,
    /// Transition count (single-base sites).
    pub transitions: usize,
    /// Transversion count (single-base sites).
    pub transversions: usize,
    /// Records with `FILTER == PASS`.
    pub passing: usize,
}

impl VcfStats {
    /// The transition/transversion ratio (`0.0` when no transversions).
    pub fn ts_tv_ratio(&self) -> f64 {
        if self.transversions == 0 {
            0.0
        } else {
            self.transitions as f64 / self.transversions as f64
        }
    }
}

/// Computes the [`VcfStats`] of a whole VCF file.
pub fn vcf_stats(vcf: &VcfFile) -> VcfStats {
    let mut s = VcfStats {
        total: vcf.records.len(),
        ..VcfStats::default()
    };
    for rec in &vcf.records {
        if rec.is_snv() {
            s.snvs += 1;
        }
        if rec.is_indel() {
            s.indels += 1;
        }
        if rec.is_multiallelic() {
            s.multiallelic += 1;
        }
        if rec.filter.iter().any(|f| f == "PASS") {
            s.passing += 1;
        }
        for alt in &rec.alt {
            if rec.reference.len() == 1 && alt.len() == 1 {
                if is_transition(&rec.reference, alt) {
                    s.transitions += 1;
                } else {
                    s.transversions += 1;
                }
            }
        }
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    fn vcf(body: &str) -> VcfFile {
        let header = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4\n";
        VcfFile::parse(&format!("{header}{body}")).unwrap()
    }

    #[test]
    fn transition_vs_transversion() {
        assert!(is_transition("A", "G"));
        assert!(is_transition("C", "T"));
        assert!(!is_transition("A", "C"));
        assert!(!is_transition("A", "T"));
    }

    #[test]
    fn site_frequency_from_genotypes() {
        // 4 samples: 0/0, 0/1, 1/1, 0/1.
        //   REF alleles: 2 + 1 + 0 + 1 = 4
        //   ALT alleles: 0 + 1 + 2 + 1 = 4
        let f = vcf("chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/1\n");
        let sf = site_frequency(&f.records[0]);
        assert_eq!(sf.total_alleles, 8);
        assert_eq!(sf.allele_counts, vec![4, 4]);
        assert!((sf.alt_frequency() - 4.0 / 8.0).abs() < 1e-9);
    }

    #[test]
    fn minor_allele_frequency() {
        // REF freq 4/8 = 0.5 -> MAF = min(0.5, 0.5) = 0.5.
        let f = vcf("chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/1\n");
        let sf = site_frequency(&f.records[0]);
        assert!((sf.minor_allele_frequency() - 0.5).abs() < 1e-9);
    }

    #[test]
    fn genotype_counts_classified() {
        let f = vcf("chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\t./.\n");
        let c = genotype_counts(&f.records[0]);
        assert_eq!(c.hom_ref, 1);
        assert_eq!(c.het, 1);
        assert_eq!(c.hom_alt, 1);
        assert_eq!(c.missing, 1);
        assert_eq!(c.called(), 3);
    }

    #[test]
    fn hardy_weinberg_at_equilibrium() {
        // p = q = 0.5: expected 2pq = 0.5. Build 4 samples with
        // het frequency exactly 0.5: 0/0, 0/1, 0/1, 1/1.
        let f = vcf("chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t0/1\t1/1\n");
        let c = genotype_counts(&f.records[0]);
        assert!((c.observed_heterozygosity() - 0.5).abs() < 1e-9);
        assert!((c.expected_heterozygosity() - 0.5).abs() < 1e-9);
        assert!(c.inbreeding_coefficient().abs() < 1e-9);
    }

    #[test]
    fn heterozygote_deficit_positive_f() {
        // All homozygous, allele freq still 0.5 -> H_obs 0, F = 1.
        let f = vcf("chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/0\t1/1\t1/1\n");
        let c = genotype_counts(&f.records[0]);
        assert!(c.inbreeding_coefficient() > 0.9);
    }

    #[test]
    fn whole_file_stats() {
        let body = "chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/1\t0/1\t0/1\n\
chr1\t200\t.\tAT\tA\t.\tq10\t.\tGT\t0/1\t0/1\t0/1\t0/1\n\
chr1\t300\t.\tC\tA\t.\tPASS\t.\tGT\t0/1\t0/1\t0/1\t0/1\n";
        let f = vcf(body);
        let s = vcf_stats(&f);
        assert_eq!(s.total, 3);
        assert_eq!(s.snvs, 2);
        assert_eq!(s.indels, 1);
        assert_eq!(s.transitions, 1); // A->G
        assert_eq!(s.transversions, 1); // C->A
        assert_eq!(s.passing, 2);
        assert!((s.ts_tv_ratio() - 1.0).abs() < 1e-9);
    }
}
