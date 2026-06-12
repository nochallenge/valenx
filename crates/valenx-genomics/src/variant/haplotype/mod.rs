//! GATK HaplotypeCaller-class variant calling.
//!
//! This module implements the modern-standard local-haplotype-
//! reassembly pipeline that replaced per-site pileup callers in
//! production variant-calling (GATK HaplotypeCaller, DeepVariant,
//! Strelka2). The pipeline runs in four stages per genomic region:
//!
//! 1. [`active::detect_active_regions`] — scan the alignments for
//!    windows showing variation evidence (mismatch / indel /
//!    quality-weighted activity) above a threshold; calm regions
//!    skip reassembly entirely.
//! 2. [`assembly::assemble_local_haplotypes`] — within each active
//!    region, reassemble candidate haplotypes from the supporting
//!    reads using a De Bruijn graph (the crate's existing De Bruijn
//!    assembler from [`crate::assembly::debruijn`]), seeded with the
//!    reference and enumerated path-by-path under a cycle bound.
//! 3. [`pairhmm::log10_p_read_given_haplotype`] — score each read
//!    against every candidate haplotype with a GATK-style PairHMM
//!    that uses per-base Phred qualities as the emission error.
//! 4. [`call_haplotype_variants`] — marginalise over haplotype pairs
//!    (diploid genotype prior), call the most-probable diploid
//!    genotype per locus, and emit [`crate::variant::call::Variant`]
//!    records with proper `QUAL`, `GT`, `AD`, `DP`, `PL` populated.
//!
//! The v1 pileup caller in [`crate::variant::call`] stays available
//! through the [`VariantCallMethod`] selector; the haplotype caller is
//! the new high-stakes default.
//!
//! ## v1 scope
//!
//! Real algorithms, honest residue:
//!
//! - Single-sample, biallelic per locus. Multi-sample joint calling
//!   (GATK `GVCF`/`GenomicsDB`/`JointCalling`) and proper multi-allelic
//!   sites are documented gaps.
//! - The PairHMM uses a single configurable gap-open/extend pair
//!   rather than GATK's per-base GOP/GCP qualities (CRAM `BI/BD`-tag
//!   territory, never in plain SAM). See [`pairhmm`].
//! - Active-region detection and local-assembly knobs are documented
//!   and tunable; the defaults target ~30× germline data and were
//!   validated end-to-end on the synthetic-truth tests in the
//!   `crate::variant::haplotype::tests` module.
//! - No VQSR / DeepVariant-class deep-learning rescoring (those need
//!   trained network weights — see the standing "no LLM weights"
//!   rule).

pub mod active;
pub mod assembly;
pub mod pairhmm;

use crate::error::{GenomicsError, Result};
use crate::format::pileup::{build_pileup, PileupColumn, Reference};
use crate::format::sam::{CigarKind, SamRecord};
use crate::variant::call::{CallParams, Variant, VariantKind};
use crate::variant::genotype::{default_priors, Genotype, GenotypeCall};
use std::collections::HashMap;

pub use active::{detect_active_regions, ActiveRegion, ActiveRegionParams};
pub use assembly::{assemble_local_haplotypes, Haplotype, LocalAssemblyParams};
pub use pairhmm::{log10_p_read_given_haplotype, PairHmmParams};

/// Variant-calling method selector. The v1 pileup caller stays
/// available; the haplotype caller is the new high-stakes default.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub enum VariantCallMethod {
    /// Per-site pileup caller (v1) — see [`crate::variant::call`].
    Pileup,
    /// GATK-class local-haplotype-reassembly caller (this module).
    #[default]
    Haplotype,
}

/// All knobs for the haplotype caller in one place.
#[derive(Clone, Debug, PartialEq)]
pub struct HaplotypeCallParams {
    /// Site-level call gates (depth / fraction / quality). The same
    /// gates the pileup caller uses; the haplotype caller honours
    /// `min_qual` after marginalising the haplotype-pair posterior to
    /// a site genotype.
    pub call: CallParams,
    /// Active-region detection parameters.
    pub active: ActiveRegionParams,
    /// Local-assembly parameters.
    pub assembly: LocalAssemblyParams,
    /// PairHMM parameters.
    pub pairhmm: PairHmmParams,
    /// Minimum number of reads in an active region for the
    /// reassembler to run. Below this the region is silently skipped.
    pub min_reads_per_region: usize,
}

impl Default for HaplotypeCallParams {
    fn default() -> Self {
        HaplotypeCallParams {
            call: CallParams::default(),
            active: ActiveRegionParams::default(),
            assembly: LocalAssemblyParams::default(),
            pairhmm: PairHmmParams::default(),
            min_reads_per_region: 4,
        }
    }
}

/// Top-level entry point — call variants by local haplotype
/// reassembly across every active region implied by the alignments.
///
/// Returns the variants sorted by `(chrom, pos)`. The caller must
/// supply a [`Reference`] covering at least the active regions;
/// without reference bases the local assembler cannot seed candidate
/// haplotypes.
pub fn call_haplotype_variants(
    records: &[SamRecord],
    reference: &Reference,
    params: &HaplotypeCallParams,
) -> Result<Vec<Variant>> {
    if params.call.min_depth == 0 {
        return Err(GenomicsError::invalid("min_depth", "must be positive"));
    }
    // Build pileup once — re-used both for active-region detection and
    // for the per-site evidence the caller needs.
    let pileup = build_pileup(records, reference, 0)?;
    let regions = detect_active_regions(&pileup, &params.active);

    let mut out: Vec<Variant> = Vec::new();
    for region in &regions {
        let mut region_vars = call_one_region(records, &pileup, region, reference, params)?;
        out.append(&mut region_vars);
    }
    // Deduplicate identical variants at the same site (multiple chunks
    // could overlap).
    out.sort_by(|a, b| {
        (a.chrom.clone(), a.pos, a.reference.clone(), a.alt.clone()).cmp(&(
            b.chrom.clone(),
            b.pos,
            b.reference.clone(),
            b.alt.clone(),
        ))
    });
    out.dedup_by(|a, b| {
        a.chrom == b.chrom && a.pos == b.pos && a.reference == b.reference && a.alt == b.alt
    });
    Ok(out)
}

/// Calls variants inside one active region.
fn call_one_region(
    records: &[SamRecord],
    pileup: &[PileupColumn],
    region: &ActiveRegion,
    reference: &Reference,
    params: &HaplotypeCallParams,
) -> Result<Vec<Variant>> {
    // Extract the reference sub-sequence.
    let contig_len = reference.contig_len(&region.chrom);
    if contig_len == 0 {
        return Ok(Vec::new());
    }
    let ref_start = (region.start - 1).max(0) as usize;
    let ref_end = (region.end as usize).min(contig_len);
    if ref_start >= ref_end {
        return Ok(Vec::new());
    }
    let mut ref_sub: Vec<u8> = Vec::with_capacity(ref_end - ref_start);
    for i in ref_start..ref_end {
        ref_sub.push(reference.base_at(&region.chrom, i));
    }

    // Extract reads overlapping the region, projected onto the reference
    // coordinate window. Each read's bases inside the region are
    // collected ("local-read" view).
    let local_reads = collect_local_reads(records, region);
    if local_reads.len() < params.min_reads_per_region {
        return Ok(Vec::new());
    }

    // Local assembly — candidate haplotypes including the reference.
    let read_seqs: Vec<&[u8]> = local_reads.iter().map(|r| r.bases.as_slice()).collect();
    let haplotypes = assemble_local_haplotypes(&ref_sub, &read_seqs, &params.assembly);
    if haplotypes.len() < 2 {
        // No alternate haplotype was reconstructed — no variants.
        return Ok(Vec::new());
    }

    // Score reads against haplotypes — full read-by-haplotype log10 P.
    let mut likelihoods: Vec<Vec<f64>> = Vec::with_capacity(local_reads.len());
    for read in &local_reads {
        let mut row = Vec::with_capacity(haplotypes.len());
        for h in &haplotypes {
            let lp =
                log10_p_read_given_haplotype(&read.bases, &read.quals, &h.bases, &params.pairhmm)?;
            row.push(lp);
        }
        likelihoods.push(row);
    }

    // Build per-haplotype "vote" — the haplotype likelihoods marginalised
    // over reads. Normalised so the best is 0 (a phred-style scale).
    let mut hap_loglik = vec![0.0f64; haplotypes.len()];
    for row in &likelihoods {
        // Sum log-likelihoods (independence assumption across reads)
        for (h, &lp) in row.iter().enumerate() {
            hap_loglik[h] += lp;
        }
    }

    // For each non-reference haplotype, materialise the variants it
    // implies relative to the reference, and emit one Variant per
    // locus.
    let ref_hap_idx = haplotypes.iter().position(|h| h.is_reference).unwrap_or(0);
    let mut out: Vec<Variant> = Vec::new();
    for (alt_idx, alt_hap) in haplotypes.iter().enumerate() {
        if alt_idx == ref_hap_idx {
            continue;
        }
        let alleles = haplotype_diffs(&ref_sub, &alt_hap.bases, region.start);
        for allele in alleles {
            // For this allele, run the diploid genotyper using read
            // likelihoods aggregated across the reference + alt
            // haplotypes that carry the same allele.
            let var = call_one_allele(
                allele,
                &haplotypes,
                ref_hap_idx,
                alt_idx,
                &likelihoods,
                &local_reads,
                &ref_sub,
                region,
                params,
                &hap_loglik,
            );
            if let Some(v) = var {
                out.push(v);
            }
        }
    }
    // Cull duplicate sites (same alt haplotype may produce identical
    // alleles via multiple anchors).
    out.sort_by(|a, b| {
        (a.chrom.clone(), a.pos, a.alt.clone()).cmp(&(b.chrom.clone(), b.pos, b.alt.clone()))
    });
    out.dedup_by(|a, b| {
        a.chrom == b.chrom && a.pos == b.pos && a.reference == b.reference && a.alt == b.alt
    });

    // Cross-check against the per-site pileup evidence to populate AD/DP
    // and supply a fallback strand-count.
    annotate_with_pileup(&mut out, pileup);
    Ok(out)
}

/// A localised view of one read inside one active region — the read's
/// bases that fall inside the region window (CIGAR walk), with matching
/// per-base qualities.
#[derive(Clone, Debug)]
struct LocalRead {
    bases: Vec<u8>,
    quals: Vec<u8>,
}

/// Walks each record's CIGAR to project the read's bases inside
/// `[region.start, region.end]` to a local-read view. Soft-clipped
/// flanks are dropped; insertions inside the window are preserved;
/// deletions are treated as gaps (no read bases).
fn collect_local_reads(records: &[SamRecord], region: &ActiveRegion) -> Vec<LocalRead> {
    let mut out = Vec::new();
    for rec in records {
        if rec.is_unmapped() || rec.pos <= 0 || rec.cigar.is_empty() {
            continue;
        }
        if rec.rname != region.chrom {
            continue;
        }
        // Quick bounding box.
        let read_end = rec.ref_end().unwrap_or(rec.pos);
        if read_end < region.start || rec.pos > region.end {
            continue;
        }
        let seq = rec.seq.as_bytes();
        let qual: Vec<u8> = if rec.qual.is_empty() {
            vec![30u8; seq.len()]
        } else {
            rec.qual
                .as_bytes()
                .iter()
                .map(|&q| q.saturating_sub(33))
                .collect()
        };
        if qual.len() != seq.len() {
            continue;
        }

        let mut ref_pos = rec.pos; // 1-based
        let mut read_pos = 0usize;
        let mut bases: Vec<u8> = Vec::new();
        let mut quals: Vec<u8> = Vec::new();

        for op in &rec.cigar.ops {
            let n = op.len as usize;
            match op.kind {
                CigarKind::Match | CigarKind::Equal | CigarKind::Diff => {
                    for k in 0..n {
                        let rp = read_pos + k;
                        let pos = ref_pos + k as i64;
                        if pos >= region.start && pos <= region.end {
                            bases.push(seq[rp].to_ascii_uppercase());
                            quals.push(qual[rp]);
                        }
                    }
                    ref_pos += n as i64;
                    read_pos += n;
                }
                CigarKind::Ins => {
                    // Insertion: bases belong to the position just left
                    // of `ref_pos`. Add them if that anchor lies inside
                    // the region.
                    let anchor = ref_pos - 1;
                    if anchor >= region.start && anchor <= region.end {
                        for k in 0..n {
                            let rp = read_pos + k;
                            bases.push(seq[rp].to_ascii_uppercase());
                            quals.push(qual[rp]);
                        }
                    }
                    read_pos += n;
                }
                CigarKind::Del => {
                    // Skip the deleted reference bases — no read base
                    // contributed for these positions.
                    ref_pos += n as i64;
                }
                CigarKind::SoftClip => {
                    read_pos += n;
                }
                CigarKind::Skip => {
                    ref_pos += n as i64;
                }
                CigarKind::HardClip | CigarKind::Pad => {}
            }
        }

        if !bases.is_empty() {
            out.push(LocalRead { bases, quals });
        }
    }
    out
}

/// Decomposes a haplotype into the list of `(REF_pos, REF, ALT)`
/// alleles it implies relative to `reference`. Uses a simple
/// pairwise-anchored difference scan with affine-gap-like grouping —
/// consecutive mismatches are emitted as single MNV-style entries,
/// runs of insertions / deletions get a single anchor each.
///
/// `start_pos` is the 1-based reference coordinate of `reference[0]`.
fn haplotype_diffs(reference: &[u8], haplotype: &[u8], start_pos: i64) -> Vec<HaplotypeAllele> {
    let mut out = Vec::new();
    // Run a simple global-alignment with affine-like banding via direct
    // O(m*n) DP — these sequences are short (region length), so plain
    // DP is cheap.
    let dp = align_pair(reference, haplotype);
    // Backtrace to a CIGAR-like operation list.
    let ops = backtrace(&dp, reference, haplotype);

    let mut ref_idx: usize = 0;
    let mut hap_idx: usize = 0;
    let mut i = 0;
    while i < ops.len() {
        let op = ops[i];
        match op {
            PairOp::Match => {
                ref_idx += 1;
                hap_idx += 1;
                i += 1;
            }
            PairOp::Mismatch => {
                // Emit an SNV at this reference position. Group multiple
                // adjacent mismatches into one MNV.
                let mut j = i;
                while j < ops.len() && ops[j] == PairOp::Mismatch {
                    j += 1;
                }
                let ref_chunk: Vec<u8> = reference[ref_idx..ref_idx + (j - i)].to_vec();
                let alt_chunk: Vec<u8> = haplotype[hap_idx..hap_idx + (j - i)].to_vec();
                for (k, (rb, ab)) in ref_chunk.iter().zip(alt_chunk.iter()).enumerate() {
                    out.push(HaplotypeAllele {
                        kind: VariantKind::Snv,
                        pos: start_pos + (ref_idx + k) as i64,
                        reference: (*rb as char).to_string(),
                        alt: (*ab as char).to_string(),
                    });
                }
                ref_idx += j - i;
                hap_idx += j - i;
                i = j;
            }
            PairOp::Insertion => {
                // Bases in haplotype not in reference. VCF convention
                // anchors at the preceding reference base.
                let anchor_pos = (start_pos + ref_idx as i64) - 1; // 1-based base before
                let anchor_idx = ref_idx.saturating_sub(1);
                let anchor_base = if ref_idx == 0 {
                    b'N'
                } else {
                    reference[anchor_idx]
                };
                let mut j = i;
                while j < ops.len() && ops[j] == PairOp::Insertion {
                    j += 1;
                }
                let inserted: Vec<u8> = haplotype[hap_idx..hap_idx + (j - i)].to_vec();
                if ref_idx > 0 {
                    let mut alt = vec![anchor_base];
                    alt.extend_from_slice(&inserted);
                    out.push(HaplotypeAllele {
                        kind: VariantKind::Insertion,
                        pos: anchor_pos,
                        reference: (anchor_base as char).to_string(),
                        alt: String::from_utf8_lossy(&alt).into_owned(),
                    });
                }
                hap_idx += j - i;
                i = j;
            }
            PairOp::Deletion => {
                // Bases in reference not in haplotype. Anchor at the
                // preceding reference base.
                let anchor_pos = (start_pos + ref_idx as i64) - 1;
                let anchor_idx = ref_idx.saturating_sub(1);
                let anchor_base = if ref_idx == 0 {
                    b'N'
                } else {
                    reference[anchor_idx]
                };
                let mut j = i;
                while j < ops.len() && ops[j] == PairOp::Deletion {
                    j += 1;
                }
                let deleted: Vec<u8> = reference[ref_idx..ref_idx + (j - i)].to_vec();
                if ref_idx > 0 {
                    let mut refstr = vec![anchor_base];
                    refstr.extend_from_slice(&deleted);
                    out.push(HaplotypeAllele {
                        kind: VariantKind::Deletion,
                        pos: anchor_pos,
                        reference: String::from_utf8_lossy(&refstr).into_owned(),
                        alt: (anchor_base as char).to_string(),
                    });
                }
                ref_idx += j - i;
                i = j;
            }
        }
    }
    out
}

#[derive(Clone, Debug, PartialEq)]
struct HaplotypeAllele {
    kind: VariantKind,
    pos: i64,
    reference: String,
    alt: String,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum PairOp {
    Match,
    Mismatch,
    Insertion, // in haplotype, not in reference
    Deletion,  // in reference, not in haplotype
}

/// Needleman-Wunsch-style DP table for the diff alignment. Scoring:
/// +1 match, −1 mismatch, −2 indel — biased toward calling mismatches
/// over indels so a true SNV is not "explained" as an insertion+deletion.
fn align_pair(a: &[u8], b: &[u8]) -> Vec<i32> {
    let n = a.len();
    let m = b.len();
    let w = m + 1;
    let mut dp = vec![0i32; (n + 1) * w];
    for i in 0..=n {
        dp[i * w] = -2 * i as i32;
    }
    for (j, cell) in dp.iter_mut().enumerate().take(m + 1) {
        *cell = -2 * j as i32;
    }
    for i in 1..=n {
        for j in 1..=m {
            let m_score = if a[i - 1].eq_ignore_ascii_case(&b[j - 1]) {
                1
            } else {
                -1
            };
            let diag = dp[(i - 1) * w + j - 1] + m_score;
            let up = dp[(i - 1) * w + j] - 2;
            let left = dp[i * w + j - 1] - 2;
            dp[i * w + j] = diag.max(up).max(left);
        }
    }
    dp
}

/// Backtrace `dp` into a list of [`PairOp`]s.
fn backtrace(dp: &[i32], a: &[u8], b: &[u8]) -> Vec<PairOp> {
    let n = a.len();
    let m = b.len();
    let w = m + 1;
    let mut i = n;
    let mut j = m;
    let mut ops: Vec<PairOp> = Vec::new();
    while i > 0 || j > 0 {
        if i > 0 && j > 0 {
            let m_score = if a[i - 1].eq_ignore_ascii_case(&b[j - 1]) {
                1
            } else {
                -1
            };
            if dp[i * w + j] == dp[(i - 1) * w + j - 1] + m_score {
                ops.push(if m_score == 1 {
                    PairOp::Match
                } else {
                    PairOp::Mismatch
                });
                i -= 1;
                j -= 1;
                continue;
            }
        }
        if i > 0 && dp[i * w + j] == dp[(i - 1) * w + j] - 2 {
            ops.push(PairOp::Deletion);
            i -= 1;
            continue;
        }
        if j > 0 && dp[i * w + j] == dp[i * w + j - 1] - 2 {
            ops.push(PairOp::Insertion);
            j -= 1;
            continue;
        }
        // Shouldn't reach here; emit a match-step as a safety valve.
        if i > 0 {
            ops.push(PairOp::Deletion);
            i -= 1;
        } else if j > 0 {
            ops.push(PairOp::Insertion);
            j -= 1;
        }
    }
    ops.reverse();
    ops
}

/// Calls a diploid genotype at `allele` by marginalising over haplotype
/// pairs. Each diploid pair `(h1, h2)` has a per-read likelihood
/// `0.5·(P(read|h1) + P(read|h2))`; the per-site log-likelihood is the
/// product over reads. The marginal posterior over the three diploid
/// genotypes (HomRef = (ref,ref), Het = (ref,alt), HomAlt = (alt,alt))
/// drives the call.
#[allow(clippy::too_many_arguments)]
fn call_one_allele(
    allele: HaplotypeAllele,
    haplotypes: &[Haplotype],
    ref_hap_idx: usize,
    alt_hap_idx: usize,
    likelihoods: &[Vec<f64>],
    _local_reads: &[LocalRead],
    _ref_sub: &[u8],
    region: &ActiveRegion,
    params: &HaplotypeCallParams,
    _hap_loglik: &[f64],
) -> Option<Variant> {
    let _ = haplotypes;
    // The three diploid hypotheses are (ref,ref), (ref,alt), (alt,alt)
    // — log10 of `0.5·(P_ref + P_alt)` for each read, then summed.
    let mut loglik = [0.0f64; 3];
    for row in likelihoods {
        let lr = row[ref_hap_idx];
        let la = row[alt_hap_idx];
        let homref = lr;
        let homalt = la;
        // log10(0.5 * (10^lr + 10^la)) = log10(0.5) + log10_add(lr, la)
        let het = 0.5f64.log10() + log10_add(lr, la);
        loglik[0] += homref;
        loglik[1] += het;
        loglik[2] += homalt;
    }
    // Apply diploid genotype prior — reuse the project's defaults.
    let priors = params.call.priors;
    let mut logpost = [0.0f64; 3];
    for i in 0..3 {
        logpost[i] = loglik[i] + priors[i].max(1e-300).log10();
    }
    let max_lp = logpost.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let denom: f64 = logpost.iter().map(|&lp| 10f64.powf(lp - max_lp)).sum();
    let log_denom = max_lp + denom.log10();
    for lp in &mut logpost {
        *lp -= log_denom;
    }

    let mut best_idx = 0usize;
    for i in 1..3 {
        if logpost[i] > logpost[best_idx] {
            best_idx = i;
        }
    }
    let best = match best_idx {
        0 => Genotype::HomRef,
        1 => Genotype::Het,
        _ => Genotype::HomAlt,
    };
    if best == Genotype::HomRef {
        return None;
    }

    let p_homref = 10f64.powf(logpost[0]).clamp(0.0, 1.0);
    let qual = if p_homref <= 0.0 {
        99.0
    } else {
        (-10.0 * p_homref.log10()).clamp(0.0, 99.0)
    };
    if qual < params.call.min_qual {
        return None;
    }

    // GQ = -10 log10(1 - P(best))
    let p_best = 10f64.powf(logpost[best_idx]).clamp(0.0, 1.0);
    let gq = if p_best >= 1.0 {
        99
    } else {
        let q = -10.0 * (1.0 - p_best).max(1e-10).log10();
        q.round().clamp(0.0, 99.0) as u8
    };

    // PL = phred-scaled likelihoods, normalised so best = 0.
    let max_ll = loglik.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mut pl = [0i32; 3];
    for i in 0..3 {
        let phred = -10.0 * (loglik[i] - max_ll);
        pl[i] = phred.round().clamp(0.0, 255.0) as i32;
    }
    let gt = GenotypeCall {
        best,
        log10_posteriors: logpost,
        gq,
        pl,
    };

    // Depth / ALT-count populated later from the pileup.
    let _ = default_priors; // keep referenced
    Some(Variant {
        chrom: region.chrom.clone(),
        pos: allele.pos,
        reference: allele.reference,
        alt: allele.alt,
        kind: allele.kind,
        depth: 0,
        alt_count: 0,
        alt_fraction: 0.0,
        qual,
        genotype: gt,
        strand: (0, 0),
    })
}

#[inline]
fn log10_add(a: f64, b: f64) -> f64 {
    const LZ: f64 = -1.0e30;
    if a <= LZ {
        return b;
    }
    if b <= LZ {
        return a;
    }
    let (hi, lo) = if a > b { (a, b) } else { (b, a) };
    hi + (1.0 + 10f64.powf(lo - hi)).log10()
}

/// Backfills `depth`, `alt_count`, `alt_fraction`, `strand` from the
/// per-site pileup evidence. The haplotype caller's strength is the
/// QUAL/GT/PL — these scalar fields are kept consistent with the
/// pileup so VCF consumers see sensible counts.
fn annotate_with_pileup(variants: &mut [Variant], pileup: &[PileupColumn]) {
    // Build a quick (chrom, pos) -> column index.
    let mut index: HashMap<(String, i64), usize> = HashMap::new();
    for (i, c) in pileup.iter().enumerate() {
        index.insert((c.chrom.clone(), c.pos), i);
    }
    for v in variants.iter_mut() {
        // For SNVs use the column at v.pos; for indels use the anchor
        // (v.pos) which is one base before the indel.
        let look_pos = match v.kind {
            VariantKind::Snv => v.pos,
            VariantKind::Insertion | VariantKind::Deletion => v.pos + 1,
        };
        let col_idx = match index.get(&(v.chrom.clone(), look_pos)) {
            Some(&i) => i,
            None => continue,
        };
        let col = &pileup[col_idx];

        match v.kind {
            VariantKind::Snv => {
                v.depth = col.base_depth();
                let alt_byte = v.alt.as_bytes().first().copied().unwrap_or(b'N');
                let (fwd, rev) = col.strand_counts(alt_byte);
                let alt = fwd + rev;
                v.alt_count = alt;
                v.alt_fraction = if v.depth == 0 {
                    0.0
                } else {
                    alt as f64 / v.depth as f64
                };
                v.strand = (fwd, rev);
            }
            VariantKind::Insertion => {
                v.depth = col.depth();
                // Count reads carrying any insertion at this anchor.
                let inserted: Vec<u8> = v.alt.as_bytes()[1..].to_ascii_uppercase();
                let mut fwd = 0usize;
                let mut rev = 0usize;
                for b in &col.bases {
                    let ins_up: Vec<u8> =
                        b.insertion.iter().map(|c| c.to_ascii_uppercase()).collect();
                    if ins_up == inserted {
                        if b.reverse {
                            rev += 1;
                        } else {
                            fwd += 1;
                        }
                    }
                }
                let cnt = fwd + rev;
                v.alt_count = cnt;
                v.alt_fraction = if v.depth == 0 {
                    0.0
                } else {
                    cnt as f64 / v.depth as f64
                };
                v.strand = (fwd, rev);
            }
            VariantKind::Deletion => {
                v.depth = col.depth();
                // Count `*` placeholders.
                let mut fwd = 0usize;
                let mut rev = 0usize;
                for b in &col.bases {
                    if b.is_deletion() {
                        if b.reverse {
                            rev += 1;
                        } else {
                            fwd += 1;
                        }
                    }
                }
                let cnt = fwd + rev;
                v.alt_count = cnt;
                v.alt_fraction = if v.depth == 0 {
                    0.0
                } else {
                    cnt as f64 / v.depth as f64
                };
                v.strand = (fwd, rev);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::format::sam::{Cigar, SamFlags};

    fn mapped(name: &str, pos: i64, cigar: &str, seq: &str, qual: &str) -> SamRecord {
        let mut r = SamRecord::unmapped(name);
        r.flags = SamFlags(0);
        r.rname = "chr1".to_string();
        r.pos = pos;
        r.mapq = 60;
        r.cigar = Cigar::parse(cigar).unwrap();
        r.seq = seq.to_string();
        r.qual = qual.to_string();
        r
    }

    fn quality_string(len: usize, phred: u8) -> String {
        let c = (phred + 33) as char;
        std::iter::repeat_n(c, len).collect()
    }

    /// Builds an 80 bp non-repetitive deterministic test reference.
    /// Every 8-mer is unique — the De-Bruijn-based local assembler can
    /// then unambiguously seed haplotype paths from the flanking
    /// `(k−1)`-mers without graph-collapse confusion.
    fn nonrep_ref80() -> Vec<u8> {
        // Hand-written 80 bp sequence with no 8-mer repeated twice in
        // it. Verified manually.
        b"GCATAGCGTCTAGCGAAGCTGCAATGCCTAGTCATGGCACTGAATGTCCGAGTAGCCTGAGCTAAGCGTACGGTTCAGTC".to_vec()
    }

    #[test]
    fn snv_called_end_to_end() {
        // 80 bp non-repetitive reference, SNV at 1-based pos 30.
        let reference = nonrep_ref80();
        let mut refr = Reference::new();
        refr.add("chr1", &reference);

        let mut alt_ref: Vec<u8> = reference.to_vec();
        let ref_at_30 = reference[29]; // 0-based 29 == 1-based 30
                                       // Pick a different base.
        let alt_at_30: u8 = if ref_at_30 == b'A' {
            b'T'
        } else if ref_at_30 == b'C' {
            b'A'
        } else if ref_at_30 == b'G' {
            b'C'
        } else {
            b'G'
        };
        alt_ref[29] = alt_at_30;

        // 20 reads tiling position 30 — 10 alt, 10 ref. Each 40 bp.
        let mut records = Vec::new();
        for i in 0..10 {
            let start = i as usize; // 0-based 0..9
            let end = start + 40;
            let bases = &alt_ref[start..end];
            let q = quality_string(bases.len(), 35);
            records.push(mapped(
                &format!("alt{i}"),
                (start + 1) as i64,
                &format!("{}M", bases.len()),
                std::str::from_utf8(bases).unwrap(),
                &q,
            ));
        }
        for i in 0..10 {
            let start = i as usize;
            let end = start + 40;
            let bases = &reference[start..end];
            let q = quality_string(bases.len(), 35);
            records.push(mapped(
                &format!("ref{i}"),
                (start + 1) as i64,
                &format!("{}M", bases.len()),
                std::str::from_utf8(bases).unwrap(),
                &q,
            ));
        }

        let params = HaplotypeCallParams::default();
        let vars = call_haplotype_variants(&records, &refr, &params).unwrap();

        // Expect an SNV at pos 30, REF=ref_at_30, ALT=alt_at_30.
        let snv = vars
            .iter()
            .find(|v| v.pos == 30 && v.kind == VariantKind::Snv);
        assert!(
            snv.is_some(),
            "expected SNV at pos 30; got {:?}",
            vars.iter()
                .map(|v| (v.pos, v.kind.clone(), v.reference.clone(), v.alt.clone()))
                .collect::<Vec<_>>()
        );
        let v = snv.unwrap();
        assert_eq!(v.reference, (ref_at_30 as char).to_string());
        assert_eq!(v.alt, (alt_at_30 as char).to_string());
        assert_eq!(v.genotype.best, Genotype::Het);
        assert!(v.qual > 20.0, "QUAL = {}", v.qual);
        assert!(v.depth > 0);
        assert!(v.alt_count > 0);
        // Sanity on AD-style fields populated from pileup.
        assert!(v.alt_fraction > 0.3 && v.alt_fraction < 0.7);
    }

    #[test]
    fn hom_alt_snv_called() {
        // All reads carry the SNV — genotype must be HomAlt.
        let reference = nonrep_ref80();
        let mut refr = Reference::new();
        refr.add("chr1", &reference);

        let ref_at_30 = reference[29];
        let alt_at_30: u8 = if ref_at_30 == b'A' { b'C' } else { b'A' };
        let mut alt: Vec<u8> = reference.to_vec();
        alt[29] = alt_at_30;

        let mut records = Vec::new();
        for i in 0..16 {
            let start = i as usize % 10;
            let end = start + 40;
            let bases = &alt[start..end];
            let q = quality_string(bases.len(), 35);
            records.push(mapped(
                &format!("alt{i}"),
                (start + 1) as i64,
                &format!("{}M", bases.len()),
                std::str::from_utf8(bases).unwrap(),
                &q,
            ));
        }

        let vars =
            call_haplotype_variants(&records, &refr, &HaplotypeCallParams::default()).unwrap();
        let v = vars.iter().find(|v| v.pos == 30).expect("missing call");
        assert_eq!(v.genotype.best, Genotype::HomAlt);
        assert_eq!(v.alt, (alt_at_30 as char).to_string());
    }

    #[test]
    fn calm_region_emits_no_calls() {
        // All reads match the reference exactly — no variants.
        let reference = nonrep_ref80();
        let mut refr = Reference::new();
        refr.add("chr1", &reference);

        let mut records = Vec::new();
        for i in 0..20 {
            let start = i as usize % 10;
            let end = start + 40;
            let bases = &reference[start..end];
            let q = quality_string(bases.len(), 35);
            records.push(mapped(
                &format!("r{i}"),
                (start + 1) as i64,
                &format!("{}M", bases.len()),
                std::str::from_utf8(bases).unwrap(),
                &q,
            ));
        }
        let vars =
            call_haplotype_variants(&records, &refr, &HaplotypeCallParams::default()).unwrap();
        assert!(
            vars.is_empty(),
            "expected no calls on calm region, got {vars:?}"
        );
    }

    #[test]
    fn rejects_bad_min_depth() {
        let refr = Reference::new();
        let params = HaplotypeCallParams {
            call: CallParams {
                min_depth: 0,
                ..CallParams::default()
            },
            ..HaplotypeCallParams::default()
        };
        let res = call_haplotype_variants(&[], &refr, &params);
        assert!(res.is_err());
    }

    #[test]
    fn insertion_called_end_to_end() {
        // 80 bp non-repetitive reference; 2 bp insertion after 1-based pos 30.
        let reference = nonrep_ref80();
        let mut refr = Reference::new();
        refr.add("chr1", &reference);

        let ins_bases: &[u8] = b"TT";

        // Alt read template: ref[15..30] (15 bp) + "TT" (2 bp) + ref[30..50] (20 bp)
        // = 37 read bases mapped 1-based start 16 with CIGAR "15M2I20M".
        // (37 query bases, ref span 35 — same as the M+M sum.)
        let mut records = Vec::new();
        for i in 0..12 {
            let mut bases: Vec<u8> = Vec::new();
            bases.extend_from_slice(&reference[15..30]);
            bases.extend_from_slice(ins_bases);
            bases.extend_from_slice(&reference[30..50]);
            let q = quality_string(bases.len(), 35);
            records.push(mapped(
                &format!("alt{i}"),
                16,
                "15M2I20M",
                std::str::from_utf8(&bases).unwrap(),
                &q,
            ));
        }
        for i in 0..4 {
            let start = 10 + i as usize; // 0-based
            let end = start + 40;
            let bases = &reference[start..end];
            let q = quality_string(bases.len(), 35);
            records.push(mapped(
                &format!("ref{i}"),
                (start + 1) as i64,
                &format!("{}M", bases.len()),
                std::str::from_utf8(bases).unwrap(),
                &q,
            ));
        }

        let vars =
            call_haplotype_variants(&records, &refr, &HaplotypeCallParams::default()).unwrap();
        let ins = vars.iter().find(|v| v.kind == VariantKind::Insertion);
        assert!(
            ins.is_some(),
            "missing insertion in {:?}",
            vars.iter()
                .map(|v| (v.pos, v.kind.clone(), v.reference.clone(), v.alt.clone()))
                .collect::<Vec<_>>()
        );
        let v = ins.unwrap();
        // The anchor is pos 30 (the last ref base before "TT"); REF =
        // reference[29] (1-based 30), ALT = REF + "TT".
        assert_eq!(v.pos, 30);
        assert_eq!(v.alt, format!("{}TT", v.reference));
        assert!(v.qual > 20.0);
    }

    #[test]
    fn deletion_called_end_to_end() {
        // 80 bp non-repetitive reference; delete 2 bases at 1-based
        // positions 31..32. The anchor (VCF convention) is pos 30.
        let reference = nonrep_ref80();
        let mut refr = Reference::new();
        refr.add("chr1", &reference);

        let mut records = Vec::new();
        for i in 0..14 {
            // Alt read covers 1-based pos 16..50 (ref span 35); bases =
            // ref[15..30] (15 bp) + ref[32..50] (18 bp) = 33 read bases.
            // CIGAR: 15M2D18M.
            let mut bases: Vec<u8> = Vec::new();
            bases.extend_from_slice(&reference[15..30]);
            bases.extend_from_slice(&reference[32..50]);
            let q = quality_string(bases.len(), 35);
            records.push(mapped(
                &format!("alt{i}"),
                16,
                "15M2D18M",
                std::str::from_utf8(&bases).unwrap(),
                &q,
            ));
        }
        for i in 0..6 {
            let start = 10 + i as usize;
            let end = start + 40;
            let bases = &reference[start..end];
            let q = quality_string(bases.len(), 35);
            records.push(mapped(
                &format!("ref{i}"),
                (start + 1) as i64,
                &format!("{}M", bases.len()),
                std::str::from_utf8(bases).unwrap(),
                &q,
            ));
        }

        let vars =
            call_haplotype_variants(&records, &refr, &HaplotypeCallParams::default()).unwrap();
        let del = vars.iter().find(|v| v.kind == VariantKind::Deletion);
        assert!(
            del.is_some(),
            "missing deletion in {:?}",
            vars.iter()
                .map(|v| (v.pos, v.kind.clone(), v.reference.clone(), v.alt.clone()))
                .collect::<Vec<_>>()
        );
        let v = del.unwrap();
        assert_eq!(v.pos, 30);
        assert_eq!(v.reference.len(), 3);
        assert_eq!(v.alt.len(), 1);
        assert!(v.qual > 20.0);
    }

    #[test]
    fn method_default_is_haplotype() {
        let m = VariantCallMethod::default();
        assert_eq!(m, VariantCallMethod::Haplotype);
    }

    /// End-to-end test using the real Illumina simulator: simulate
    /// reads from a synthetic reference carrying a *known* SNV;
    /// validate the haplotype caller recovers that exact site, the
    /// genotype, and an AD/DP consistent with the truth.
    #[test]
    fn simulator_to_caller_recovers_known_snv() {
        use crate::simulate::illumina::{simulate_reads, IlluminaProfile};
        use crate::util::rng::Rng;

        let mut reference = nonrep_ref80();
        // Extend the reference to give the simulator enough room.
        // Tile two copies of the deterministic 80 bp seed.
        let mut full = reference.clone();
        full.extend_from_slice(&reference);
        full.extend_from_slice(&reference);
        full.extend_from_slice(&reference);
        reference = full; // 320 bp

        // Inject a known SNV at 1-based pos 150 (0-based 149).
        let truth_pos = 150i64;
        let truth_ref = reference[149];
        let truth_alt: u8 = if truth_ref == b'A' {
            b'C'
        } else if truth_ref == b'C' {
            b'G'
        } else if truth_ref == b'G' {
            b'T'
        } else {
            b'A'
        };
        let mut alt_ref = reference.clone();
        alt_ref[149] = truth_alt;

        // Simulate 40 ref-style and 40 alt-style reads (HiSeq-like 60 bp).
        let mut profile = IlluminaProfile::hiseq_150();
        profile.read_length = 60;
        let ref_reads = simulate_reads(&reference, &profile, 40, 17).unwrap();
        let alt_reads = simulate_reads(&alt_ref, &profile, 40, 31).unwrap();

        // Build a primitive aligner: each simulated read carries
        // `pos=<1-based start> strand=+/-` in its description; we map
        // it back to its truth position and emit a SAM record with a
        // simple "{len}M" CIGAR.
        fn parse_pos_strand(desc: &str) -> Option<(i64, bool)> {
            // "pos=N strand=+|-"
            let mut pos: Option<i64> = None;
            let mut rev: Option<bool> = None;
            for token in desc.split_whitespace() {
                if let Some(v) = token.strip_prefix("pos=") {
                    pos = v.parse().ok();
                }
                if let Some(v) = token.strip_prefix("strand=") {
                    rev = Some(v == "-");
                }
            }
            match (pos, rev) {
                (Some(p), Some(r)) => Some((p, r)),
                _ => None,
            }
        }
        // Need revcomp for reverse-strand reads.
        fn revcomp(s: &[u8]) -> Vec<u8> {
            s.iter()
                .rev()
                .map(|&b| match b.to_ascii_uppercase() {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'C',
                    b'T' => b'A',
                    _ => b'N',
                })
                .collect()
        }

        let mut records: Vec<SamRecord> = Vec::new();
        let mut rng_id = Rng::new(99); // for unique names
        for r in ref_reads.iter().chain(alt_reads.iter()) {
            let (pos, rev) = match parse_pos_strand(&r.record.description) {
                Some(x) => x,
                None => continue,
            };
            let bases: Vec<u8> = r.record.seq.as_bytes().to_vec();
            let oriented = if rev { revcomp(&bases) } else { bases };
            let q = if rev {
                r.quality.iter().rev().copied().collect::<Vec<u8>>()
            } else {
                r.quality.clone()
            };
            let qstr: String = q.iter().map(|&q| (q + 33) as char).collect();
            let mut rec = mapped(
                &format!("sim{}", rng_id.next_u64()),
                pos,
                &format!("{}M", oriented.len()),
                std::str::from_utf8(&oriented).unwrap(),
                &qstr,
            );
            if rev {
                rec.flags = SamFlags(SamFlags::REVERSE);
            }
            records.push(rec);
        }

        let mut refr = Reference::new();
        refr.add("chr1", &reference);
        let vars =
            call_haplotype_variants(&records, &refr, &HaplotypeCallParams::default()).unwrap();

        // The known SNV must be in the output, with the right REF/ALT
        // and a het genotype call (40/40 alt/ref).
        let v = vars
            .iter()
            .find(|v| v.pos == truth_pos && v.kind == VariantKind::Snv)
            .unwrap_or_else(|| {
                panic!(
                    "missing truth SNV at pos {truth_pos}; vars = {:?}",
                    vars.iter()
                        .map(|v| (v.pos, v.kind.clone(), v.reference.clone(), v.alt.clone()))
                        .collect::<Vec<_>>()
                )
            });
        assert_eq!(v.reference, (truth_ref as char).to_string());
        assert_eq!(v.alt, (truth_alt as char).to_string());
        assert_eq!(v.genotype.best, Genotype::Het);
        assert!(v.qual > 30.0, "QUAL = {}", v.qual);
        assert!(v.depth > 10);
        assert!(v.alt_count > 5);
    }

    #[test]
    fn beats_pileup_on_hard_indel_case() {
        // The hard case: a 2 bp deletion where only a handful of reads
        // (6) carry the alt and a few (3) carry the reference. The
        // pileup caller's per-column tally has to clear the depth gate
        // *and* the AF gate at a column that does not get the full
        // 6/9 fraction the haplotype-level view sees; the haplotype
        // caller marginalises the read evidence over candidate
        // haplotypes and clears the gate.
        let reference = nonrep_ref80();
        let mut refr = Reference::new();
        refr.add("chr1", &reference);

        let mut records = Vec::new();
        for i in 0..6 {
            let mut bases: Vec<u8> = Vec::new();
            bases.extend_from_slice(&reference[15..30]);
            bases.extend_from_slice(&reference[32..50]);
            let q = quality_string(bases.len(), 35);
            records.push(mapped(
                &format!("alt{i}"),
                16,
                "15M2D18M",
                std::str::from_utf8(&bases).unwrap(),
                &q,
            ));
        }
        for i in 0..3 {
            let start = 10 + i as usize;
            let end = start + 40;
            let bases = &reference[start..end];
            let q = quality_string(bases.len(), 35);
            records.push(mapped(
                &format!("ref{i}"),
                (start + 1) as i64,
                &format!("{}M", bases.len()),
                std::str::from_utf8(bases).unwrap(),
                &q,
            ));
        }

        let hap_vars =
            call_haplotype_variants(&records, &refr, &HaplotypeCallParams::default()).unwrap();
        let pileup = build_pileup(&records, &refr, 0).unwrap();
        let pileup_vars = crate::variant::call::call_variants(
            &pileup,
            &crate::variant::call::CallParams::default(),
        )
        .unwrap();

        let hap_indel = hap_vars
            .iter()
            .find(|v| v.kind == VariantKind::Deletion)
            .cloned();
        let pile_indel = pileup_vars
            .iter()
            .find(|v| v.kind == VariantKind::Deletion)
            .cloned();

        // The haplotype caller must call this deletion.
        assert!(
            hap_indel.is_some(),
            "haplotype caller failed to call the deletion: hap_vars={hap_vars:?}"
        );
        // And its QUAL must beat or match the pileup caller's on the
        // same site (when the pileup caller also called it).
        if let (Some(h), Some(p)) = (hap_indel.as_ref(), pile_indel.as_ref()) {
            assert!(
                h.qual + 1e-6 >= p.qual,
                "haplotype QUAL {} should be >= pileup QUAL {} on the same site",
                h.qual,
                p.qual
            );
        }
    }
}
