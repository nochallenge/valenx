//! Genome-wide off-target search via the production FM-index.
//!
//! The default off-target enumerator in `valenx-genomics`
//! ([`valenx_genomics::crispr::offtarget::enumerate_off_targets`]) is a
//! direct `O(|genome| × |guide|)` two-strand sweep with a per-window
//! mismatch count. It is correct, but it does not scale to a real
//! mammalian genome: at 3·10⁹ bp × 20 nt × 2 strands the constant is
//! enormous, and the whole `valenx-genediting` v1 only used it through
//! the `off_target_genome: Vec<(String, Vec<u8>)>` field on a guide
//! design request — a single user-supplied window.
//!
//! This module closes the genome-wide gap. It builds an FM-index over
//! each contig of the reference (reusing
//! [`valenx_align::search::FmIndex`] — a real SA-IS index with a
//! block-sampled rank + sampled SA, the same layout BWA uses), then
//! runs the classical **seed-and-extend mismatch-tolerant search**
//! BWA / Cas-OFFinder use:
//!
//! 1. Split the guide into `k + 1` non-overlapping seeds of roughly
//!    equal length. The Pigeonhole Principle says that for any `k`
//!    mismatches at least one of the `k + 1` seeds must be an exact
//!    match.
//! 2. For every seed, every contig and every strand: backward-search
//!    the seed in the FM-index — `O(|seed|)` per query, independent of
//!    the genome size — and resolve the matches to text positions.
//! 3. For each seed hit: extend left and right to the full guide
//!    footprint, counting mismatches over the whole guide and bailing
//!    once the budget is blown.
//! 4. Filter by PAM. Score with the existing CFD heuristic from
//!    `valenx-genomics`. Deduplicate seeds that hit the same
//!    forward-strand position.
//!
//! The reverse strand is handled by reverse-complementing the guide
//! and running the same procedure: the FM-index is built on the
//! forward strand only, but the guide rotates.
//!
//! ## v1 scope
//!
//! The FM-index from `valenx-align` is the production SA-IS / sampled-
//! rank layout; the seed-and-extend loop is the standard mismatch-
//! tolerant search. The "k+1 seeds" split is the pigeonhole-bounded
//! recipe — no heuristic fall-back is required for the proven mismatch
//! budget. The CFD score is the same transparent heuristic the
//! existing crate uses (no trained model — the project's "no trained-
//! weights" rule). Bulges (RNA/DNA insertions) are not modelled.

use crate::error::{GeneditingError, Result};
use crate::sequtil::{is_acgt, revcomp};
use std::collections::HashSet;
use valenx_align::search::FmIndex;
use valenx_genomics::crispr::guide::{iupac_match, PamSide, PamSpec};
use valenx_genomics::crispr::offtarget::{cfd_score, OffTarget};

/// One indexed contig of the search genome.
///
/// The FM-index is built once per contig and reused across every
/// guide query — that is the whole point of the data structure.
#[derive(Clone, Debug)]
pub struct IndexedContig {
    /// Contig name (e.g. `"chr1"`).
    pub name: String,
    /// Forward-strand sequence (uppercased).
    forward: Vec<u8>,
    /// FM-index over the forward strand.
    fm_forward: FmIndex,
}

impl IndexedContig {
    /// Builds a single-contig FM-index.
    ///
    /// # Errors
    /// - [`GeneditingError::InvalidTarget`] for a non-ACGT contig (the
    ///   FM-index can index any byte string with no NUL byte, but the
    ///   off-target scan is over A/C/G/T only).
    /// - [`GeneditingError::Invalid`] if the FM-index builder rejects
    ///   the contig (empty input).
    pub fn build(name: impl Into<String>, sequence: &[u8]) -> Result<Self> {
        if !is_acgt(sequence) {
            return Err(GeneditingError::invalid_target(
                "region",
                "contig must be a non-empty ACGT sequence",
            ));
        }
        let forward: Vec<u8> = sequence.iter().map(|b| b.to_ascii_uppercase()).collect();
        let fm_forward = FmIndex::build(&forward).map_err(|e| {
            GeneditingError::invalid("contig", format!("FM-index build failed: {e}"))
        })?;
        Ok(IndexedContig {
            name: name.into(),
            forward,
            fm_forward,
        })
    }

    /// Length of the indexed forward strand.
    pub fn len(&self) -> usize {
        self.forward.len()
    }

    /// `true` when the contig is empty (never produced — present for
    /// clippy hygiene).
    pub fn is_empty(&self) -> bool {
        self.forward.is_empty()
    }

    /// Read-only access to the forward-strand bases (uppercased).
    pub fn forward_bases(&self) -> &[u8] {
        &self.forward
    }
}

/// A genome FM-index — one contig at a time, named.
///
/// Build it **once** per reference and share it across every guide
/// query. The expensive step (SA-IS suffix-array + BWT + rank tables)
/// happens inside [`IndexedContig::build`]; off-target queries reuse
/// the index for every guide.
#[derive(Clone, Debug, Default)]
pub struct GenomeIndex {
    contigs: Vec<IndexedContig>,
}

impl GenomeIndex {
    /// An empty index. Add contigs with [`Self::add_contig`].
    pub fn new() -> Self {
        GenomeIndex {
            contigs: Vec::new(),
        }
    }

    /// Builds the index from a list of named contigs in one shot.
    ///
    /// # Errors
    /// Forwards [`IndexedContig::build`] errors.
    pub fn build(contigs: &[(String, Vec<u8>)]) -> Result<Self> {
        let mut g = GenomeIndex::new();
        for (name, seq) in contigs {
            g.add_contig(name.clone(), seq)?;
        }
        Ok(g)
    }

    /// Adds (and builds the FM-index for) one contig.
    pub fn add_contig(&mut self, name: impl Into<String>, sequence: &[u8]) -> Result<()> {
        self.contigs.push(IndexedContig::build(name, sequence)?);
        Ok(())
    }

    /// Number of contigs in the index.
    pub fn contig_count(&self) -> usize {
        self.contigs.len()
    }

    /// Total indexed length (sum of contig lengths in bp).
    pub fn total_length(&self) -> usize {
        self.contigs.iter().map(|c| c.len()).sum()
    }

    /// Read-only view of the indexed contigs.
    pub fn contigs(&self) -> &[IndexedContig] {
        &self.contigs
    }
}

/// Finds every off-target site for `guide` across the genome FM-index,
/// within `max_mismatches`.
///
/// The standard BWA-style **seed-and-extend** loop:
///
/// 1. Split the guide into `k + 1` non-overlapping seeds. By the
///    Pigeonhole Principle, at least one seed is an exact match at any
///    real off-target.
/// 2. Exact-match every seed against every contig + strand FM-index.
/// 3. For each seed hit, extend left and right to the full guide
///    footprint and verify the mismatch count is `<= k`.
/// 4. Filter by the PAM (read from the genome at the predicted PAM
///    span on the matching strand).
/// 5. Deduplicate (same forward-strand start + strand can be reached
///    via several seeds).
/// 6. Score with the existing CFD heuristic and sort by descending
///    score (most-dangerous-first).
///
/// # Errors
/// - [`GeneditingError::InvalidTarget`] for a non-ACGT guide.
/// - [`GeneditingError::Invalid`] for a guide whose length disagrees
///   with the [`PamSpec`].
pub fn find_off_targets_genome(
    guide: &[u8],
    genome: &GenomeIndex,
    pam: &PamSpec,
    max_mismatches: usize,
) -> Result<Vec<OffTarget>> {
    if !is_acgt(guide) {
        return Err(GeneditingError::invalid_target(
            "locus",
            "guide must be a non-empty ACGT sequence",
        ));
    }
    if guide.len() != pam.protospacer_len {
        return Err(GeneditingError::invalid(
            "guide",
            format!(
                "guide length {} != PAM-spec protospacer length {}",
                guide.len(),
                pam.protospacer_len
            ),
        ));
    }
    let guide_u: Vec<u8> = guide.iter().map(|b| b.to_ascii_uppercase()).collect();
    let guide_rc = revcomp(&guide_u);

    // Pigeonhole seeding: k+1 non-overlapping seeds of length
    // floor(plen / (k+1)). If max_mismatches == 0, a single full-length
    // exact seed is correct.
    let plen = guide_u.len();
    let n_seeds = max_mismatches + 1;
    let seeds_fwd = pigeonhole_seeds(&guide_u, n_seeds);
    let seeds_rev = pigeonhole_seeds(&guide_rc, n_seeds);

    let motif = pam.motif.as_bytes();
    let pamlen = pam.pam_len();

    // Deduplicate by (contig_idx, fwd_start, reverse).
    let mut seen: HashSet<(usize, usize, bool)> = HashSet::new();
    let mut hits: Vec<OffTarget> = Vec::new();

    for (ci, contig) in genome.contigs.iter().enumerate() {
        let n = contig.len();
        if n < plen + pamlen {
            continue;
        }
        // Forward strand: search the guide directly, candidates are
        // forward-strand protospacer windows.
        scan_seeds_one_strand(
            ci,
            contig,
            &guide_u,
            &guide_u,
            &seeds_fwd,
            motif,
            pam.side,
            pamlen,
            plen,
            max_mismatches,
            false,
            &mut seen,
            &mut hits,
        );
        // Reverse strand: an off-target on the reverse strand of the
        // forward axis is an exact match of `guide_rc` against the
        // forward FM-index — those text positions are the protospacer's
        // **forward-strand** coordinates already. The mismatch check
        // compares the guide-orientation match (guide_rc vs
        // forward-strand bases), but the *reported* protospacer is
        // read on the reverse strand (so we revcomp it back).
        scan_seeds_one_strand(
            ci,
            contig,
            &guide_u,
            &guide_rc,
            &seeds_rev,
            motif,
            pam.side,
            pamlen,
            plen,
            max_mismatches,
            true,
            &mut seen,
            &mut hits,
        );
    }

    hits.sort_by(|a, b| {
        b.cfd_score
            .partial_cmp(&a.cfd_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    Ok(hits)
}

/// Splits `guide` into `n_seeds` non-overlapping seeds and records each
/// `(offset, seed_bytes)`. The first `plen % n_seeds` seeds get one
/// extra base so the seeds collectively cover the guide.
fn pigeonhole_seeds(guide: &[u8], n_seeds: usize) -> Vec<(usize, Vec<u8>)> {
    let plen = guide.len();
    let n_seeds = n_seeds.max(1);
    let base = plen / n_seeds;
    let rem = plen % n_seeds;
    let mut out = Vec::with_capacity(n_seeds);
    let mut off = 0usize;
    for i in 0..n_seeds {
        let len = base + if i < rem { 1 } else { 0 };
        if len == 0 {
            continue;
        }
        out.push((off, guide[off..off + len].to_vec()));
        off += len;
    }
    out
}

/// Scans the FM-index of one contig for off-targets of one orientation
/// of the guide (forward or reverse-complement).
///
/// `guide_canonical` is the original guide as the user supplied it
/// (5′→3′ on its own strand) — what gets reported and CFD-scored.
/// `guide_oriented` is what was actually backward-searched: equal to
/// `guide_canonical` on the forward sweep, equal to `revcomp(guide)`
/// on the reverse sweep.
#[allow(clippy::too_many_arguments)]
fn scan_seeds_one_strand(
    contig_idx: usize,
    contig: &IndexedContig,
    guide_canonical: &[u8],
    guide_oriented: &[u8],
    seeds: &[(usize, Vec<u8>)],
    pam_motif: &[u8],
    pam_side: PamSide,
    pamlen: usize,
    plen: usize,
    max_mismatches: usize,
    reverse: bool,
    seen: &mut HashSet<(usize, usize, bool)>,
    out: &mut Vec<OffTarget>,
) {
    let n = contig.len();
    let total = plen + pamlen;
    for (seed_offset, seed) in seeds {
        // Exact backward-search every position of this seed in the
        // FM-index. The hits are text positions where the seed starts
        // on the contig (forward-strand coordinates).
        let positions = contig.fm_forward.locate(seed);
        for &p in &positions {
            // Extend back to the guide's hypothetical start on the
            // contig's forward strand.
            if p < *seed_offset {
                continue;
            }
            let proto_start = p - seed_offset;
            // The protospacer window must fit AND the PAM window
            // (on the same forward-axis layout used by genomics).
            //
            // For a 3' PAM on the forward strand → PAM at
            // `proto_end..proto_end+pamlen`. For a 3' PAM scanned via
            // the reverse-complemented guide on the forward FM-index,
            // the "protospacer" we matched is at proto_start..proto_end
            // on the forward strand, but its PAM (on the reverse
            // strand of the forward axis) sits *5'* of proto_start
            // on the forward axis (i.e. at proto_start - pamlen).
            //
            // For a 5' PAM (Cas12a) the geometry mirrors.
            let (window_start, window_end) = if !reverse {
                match pam_side {
                    PamSide::ThreePrime => (proto_start, proto_start + total),
                    PamSide::FivePrime => {
                        if proto_start < pamlen {
                            continue;
                        }
                        (proto_start - pamlen, proto_start + plen)
                    }
                }
            } else {
                match pam_side {
                    PamSide::ThreePrime => {
                        if proto_start < pamlen {
                            continue;
                        }
                        (proto_start - pamlen, proto_start + plen)
                    }
                    PamSide::FivePrime => (proto_start, proto_start + total),
                }
            };
            if window_end > n {
                continue;
            }
            let proto_end = proto_start + plen;
            // Read the PAM at the predicted forward-strand position;
            // for a reverse-strand hit we still read from the forward
            // strand but flip it to the reverse-complement to compare.
            let pam_window_fwd: &[u8] = if !reverse {
                match pam_side {
                    PamSide::ThreePrime => &contig.forward[proto_end..proto_end + pamlen],
                    PamSide::FivePrime => &contig.forward[window_start..window_start + pamlen],
                }
            } else {
                match pam_side {
                    PamSide::ThreePrime => &contig.forward[window_start..window_start + pamlen],
                    PamSide::FivePrime => &contig.forward[proto_end..proto_end + pamlen],
                }
            };
            let pam_strand: Vec<u8> = if reverse {
                // PAM on the reverse strand: revcomp of the forward-
                // strand window at that location.
                revcomp(pam_window_fwd)
            } else {
                pam_window_fwd.to_vec()
            };
            if !motif_matches(pam_motif, &pam_strand) {
                continue;
            }
            // Full-guide mismatch verification. The pigeonhole seed
            // is *some* slice of `guide_oriented` that matched exactly
            // at the corresponding forward-strand position; the rest
            // of the guide window may contain mismatches up to `k`.
            //
            // Compare in the *backward-search orientation*: guide_oriented
            // vs the forward-strand bases of the contig at proto_start.
            // Mismatch positions in that frame map directly to the
            // guide's 5′→3′ axis for the forward sweep; for the reverse
            // sweep, position `i` in `guide_rc` corresponds to position
            // `plen - 1 - i` in the original guide (because revcomp
            // reverses the order).
            let proto_window_fwd = &contig.forward[proto_start..proto_end];
            let (mm, positions_in_oriented) = mismatches(guide_oriented, proto_window_fwd);
            if mm > max_mismatches {
                continue;
            }
            let key = (contig_idx, proto_start, reverse);
            if !seen.insert(key) {
                continue;
            }
            // Report the protospacer in the guide's canonical 5′→3′
            // orientation (revcomp of the forward bases for a reverse
            // hit), so CFD scoring lines up with the canonical guide.
            let proto_canonical: Vec<u8> = if reverse {
                revcomp(proto_window_fwd)
            } else {
                proto_window_fwd.to_vec()
            };
            // Map mismatch positions into the canonical guide axis.
            let positions_mm: Vec<usize> = if reverse {
                let mut p: Vec<usize> = positions_in_oriented
                    .iter()
                    .map(|&i| plen - 1 - i)
                    .collect();
                p.sort_unstable();
                p
            } else {
                positions_in_oriented
            };
            let cfd = cfd_score(guide_canonical, &proto_canonical, &pam_strand);
            out.push(OffTarget {
                chrom: contig.name.clone(),
                start: proto_start,
                reverse,
                protospacer: String::from_utf8_lossy(&proto_canonical).into_owned(),
                pam: String::from_utf8_lossy(&pam_strand).into_owned(),
                mismatches: mm,
                mismatch_positions: positions_mm,
                cfd_score: cfd,
            });
        }
    }
}

/// Returns `(mismatch_count, mismatch_positions)` for two equal-length
/// sequences. Positions are 0-based against the guide's 5′→3′ axis.
fn mismatches(guide: &[u8], proto: &[u8]) -> (usize, Vec<usize>) {
    let mut count = 0usize;
    let mut positions = Vec::new();
    for (i, (&g, &p)) in guide.iter().zip(proto).enumerate() {
        if !g.eq_ignore_ascii_case(&p) {
            count += 1;
            positions.push(i);
        }
    }
    (count, positions)
}

/// `true` when `window` matches IUPAC `motif` base-for-base.
fn motif_matches(motif: &[u8], window: &[u8]) -> bool {
    motif.len() == window.len() && motif.iter().zip(window).all(|(&c, &b)| iupac_match(c, b))
}

/// A compact summary of an FM-index-backed off-target scan: where the
/// hits are, plus the CRISPOR-style aggregate specificity score.
///
/// Not `serde`-derived: it carries [`valenx_genomics`] `OffTarget`
/// values which are transient analysis outputs (see the matching note
/// on [`crate::therapy::safety::SafetyScreenInput`]).
#[derive(Clone, Debug, PartialEq)]
pub struct OffTargetReport {
    /// All off-target hits (including the perfect on-target site if it
    /// is in the indexed genome), sorted by descending CFD-style score.
    pub hits: Vec<OffTarget>,
    /// CRISPOR-style aggregate guide specificity in `(0, 1]`:
    /// `1 / (1 + Σ cfd_off)` over the non-perfect hits. A guide with no
    /// off-targets scores `1.0`; many active off-targets drive it
    /// toward `0`.
    pub specificity: f64,
    /// Number of non-perfect off-target hits considered.
    pub off_target_count: usize,
    /// Number of perfect (zero-mismatch) hits — the on-target site
    /// itself in the indexed reference.
    pub perfect_hits: usize,
}

/// Wraps [`find_off_targets_genome`] and adds the CRISPOR-style
/// specificity aggregate the existing guide-design module also uses.
pub fn off_target_report(
    guide: &[u8],
    genome: &GenomeIndex,
    pam: &PamSpec,
    max_mismatches: usize,
) -> Result<OffTargetReport> {
    let hits = find_off_targets_genome(guide, genome, pam, max_mismatches)?;
    let mut off_target_count = 0usize;
    let mut perfect_hits = 0usize;
    let mut cfd_sum = 0.0f64;
    for h in &hits {
        if h.is_perfect() {
            perfect_hits += 1;
        } else {
            off_target_count += 1;
            cfd_sum += h.cfd_score;
        }
    }
    let specificity = 1.0 / (1.0 + cfd_sum);
    Ok(OffTargetReport {
        hits,
        specificity: specificity.clamp(0.0, 1.0),
        off_target_count,
        perfect_hits,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pam_spcas9() -> PamSpec {
        PamSpec::spcas9()
    }

    #[test]
    fn pigeonhole_split_covers_guide() {
        let g = b"ACGTACGTACGTACGTACGT";
        // k = 2 → 3 seeds; 20 / 3 = 6 r 2 → seeds of 7,7,6.
        let seeds = pigeonhole_seeds(g, 3);
        assert_eq!(seeds.len(), 3);
        assert_eq!(seeds[0].1.len(), 7);
        assert_eq!(seeds[1].1.len(), 7);
        assert_eq!(seeds[2].1.len(), 6);
        let mut reconstructed = Vec::new();
        for (off, s) in &seeds {
            assert_eq!(*off, reconstructed.len());
            reconstructed.extend_from_slice(s);
        }
        assert_eq!(reconstructed, g);
    }

    #[test]
    fn finds_perfect_on_target_hit() {
        // A 1 kb contig with the on-target site planted at position 200.
        let proto = b"ACGTACGTACGTACGTACGT";
        let pam = b"AGG";
        let mut chrom = vec![b'T'; 1000];
        chrom[200..220].copy_from_slice(proto);
        chrom[220..223].copy_from_slice(pam);
        let genome = GenomeIndex::build(&[("chr1".to_string(), chrom)]).unwrap();
        let hits = find_off_targets_genome(proto, &genome, &pam_spcas9(), 0).unwrap();
        assert!(!hits.is_empty());
        let h = hits.iter().find(|h| !h.reverse && h.start == 200).unwrap();
        assert_eq!(h.mismatches, 0);
        assert!((h.cfd_score - 1.0).abs() < 1e-9);
        assert_eq!(h.chrom, "chr1");
    }

    #[test]
    fn finds_known_off_targets_with_mismatches() {
        // Plant three sites: perfect at 200, 1-mismatch at 500, 3-mismatch
        // at 800. All have a valid PAM. With k=3, all three must be
        // recovered.
        //
        // Guide: A C G T A C G T A C G T A C G T A C G T
        //   pos: 0 1 2 3 4 5 6 7 8 9 ...               19
        let mut chrom = vec![b'T'; 2000];
        // Perfect.
        chrom[200..220].copy_from_slice(b"ACGTACGTACGTACGTACGT");
        chrom[220..223].copy_from_slice(b"AGG");
        // 1 mismatch at position 19 (T->C).
        chrom[500..520].copy_from_slice(b"ACGTACGTACGTACGTACGC");
        chrom[520..523].copy_from_slice(b"TGG");
        // 3 mismatches at positions 14, 17 and 19:
        //   target: A C G T A C G T A C G T A C [A] T A [G] G [A]
        //   guide : A C G T A C G T A C G T A C [G] T A [C] G [T]
        chrom[800..820].copy_from_slice(b"ACGTACGTACGTACATAGGA");
        chrom[820..823].copy_from_slice(b"CGG");

        let genome = GenomeIndex::build(&[("chr1".to_string(), chrom)]).unwrap();
        let hits =
            find_off_targets_genome(b"ACGTACGTACGTACGTACGT", &genome, &pam_spcas9(), 3).unwrap();
        assert!(hits.iter().any(|h| h.start == 200 && h.mismatches == 0));
        assert!(hits.iter().any(|h| h.start == 500 && h.mismatches == 1));
        assert!(hits.iter().any(|h| h.start == 800 && h.mismatches == 3));
        // CFD scores are ranked: perfect 1.0 strictly dominates.
        let perfect_score = hits
            .iter()
            .find(|h| h.start == 200)
            .map(|h| h.cfd_score)
            .unwrap();
        let three_mm_score = hits
            .iter()
            .find(|h| h.start == 800)
            .map(|h| h.cfd_score)
            .unwrap();
        assert!(perfect_score > three_mm_score);
    }

    #[test]
    fn respects_mismatch_budget() {
        // A 4-mismatch off-target with k=3 must be rejected, k=4 must
        // catch it.
        // guide:  A C G T A C G T A C G T A C [G] T A [C] [G] [T]
        // target: A C G T A C G T A C G T A C [A] T A [A] [A] [A]
        // → mismatches at positions 14, 17, 18, 19 = 4 total.
        let mut chrom = vec![b'T'; 1000];
        chrom[500..520].copy_from_slice(b"ACGTACGTACGTACATAAAA");
        chrom[520..523].copy_from_slice(b"AGG");
        let genome = GenomeIndex::build(&[("chr1".to_string(), chrom)]).unwrap();
        let guide = b"ACGTACGTACGTACGTACGT";
        let pam = pam_spcas9();
        let k3 = find_off_targets_genome(guide, &genome, &pam, 3).unwrap();
        assert!(k3.iter().all(|h| h.start != 500));
        let k4 = find_off_targets_genome(guide, &genome, &pam, 4).unwrap();
        assert!(k4.iter().any(|h| h.start == 500 && h.mismatches == 4));
    }

    #[test]
    fn rejects_hit_without_pam() {
        // A protospacer at 200 but no valid NGG PAM downstream — must
        // not be reported.
        let mut chrom = vec![b'T'; 1000];
        chrom[200..220].copy_from_slice(b"ACGTACGTACGTACGTACGT");
        chrom[220..223].copy_from_slice(b"AAA"); // no NGG
        let genome = GenomeIndex::build(&[("chr1".to_string(), chrom)]).unwrap();
        let hits =
            find_off_targets_genome(b"ACGTACGTACGTACGTACGT", &genome, &pam_spcas9(), 0).unwrap();
        assert!(hits.iter().all(|h| h.start != 200 || h.reverse));
    }

    #[test]
    fn finds_reverse_strand_hit() {
        // Plant the reverse-complement of the protospacer on the
        // forward strand so the guide matches the reverse strand.
        // Guide: ACGTACGTACGTACGTACGT, revcomp = ACGTACGTACGTACGTACGT.
        // That guide is its own revcomp, so use a non-palindromic one.
        let guide = b"AAAACCCCGGGGTTTTACGT";
        let guide_rc = revcomp(guide);
        // On the reverse strand, the PAM (3') is "downstream" of the
        // protospacer's 3' end. Reading on the forward axis the PAM
        // sits *5' of* the protospacer's forward position, i.e. at
        // proto_start - 3. So place the revcomp of "NGG" PAM at the
        // forward position immediately before the guide_rc window.
        let mut chrom = vec![b'A'; 1000];
        // PAM at 197..200 on forward strand = "CCT" (revcomp of "AGG").
        chrom[197..200].copy_from_slice(b"CCT");
        // Revcomp-guide at 200..220 on forward strand.
        chrom[200..220].copy_from_slice(&guide_rc);
        let genome = GenomeIndex::build(&[("chr1".to_string(), chrom)]).unwrap();
        let hits = find_off_targets_genome(guide, &genome, &pam_spcas9(), 0).unwrap();
        let h = hits
            .iter()
            .find(|h| h.reverse && h.start == 200)
            .expect("expected reverse-strand hit at 200");
        assert_eq!(h.mismatches, 0);
        assert_eq!(h.pam, "AGG");
    }

    #[test]
    fn deduplicates_multi_seed_hits() {
        // A perfect on-target hit will be found by all (k+1) seeds.
        // The reported list must contain exactly one entry per
        // (contig, start, strand).
        let mut chrom = vec![b'T'; 1000];
        chrom[200..220].copy_from_slice(b"ACGTACGTACGTACGTACGT");
        chrom[220..223].copy_from_slice(b"AGG");
        let genome = GenomeIndex::build(&[("chr1".to_string(), chrom)]).unwrap();
        let hits =
            find_off_targets_genome(b"ACGTACGTACGTACGTACGT", &genome, &pam_spcas9(), 3).unwrap();
        let n_at_200 = hits.iter().filter(|h| h.start == 200 && !h.reverse).count();
        assert_eq!(n_at_200, 1, "expected exactly one perfect hit");
    }

    #[test]
    fn aggregate_report_specificity() {
        // One perfect site + one weak off-target.
        let mut chrom = vec![b'T'; 1000];
        chrom[200..220].copy_from_slice(b"ACGTACGTACGTACGTACGT");
        chrom[220..223].copy_from_slice(b"AGG");
        chrom[500..520].copy_from_slice(b"ACGTACGTACGTACGTACGC"); // 1mm
        chrom[520..523].copy_from_slice(b"TGG");
        let genome = GenomeIndex::build(&[("chr1".to_string(), chrom)]).unwrap();
        let report = off_target_report(b"ACGTACGTACGTACGTACGT", &genome, &pam_spcas9(), 3).unwrap();
        assert_eq!(report.perfect_hits, 1);
        assert!(report.off_target_count >= 1);
        assert!(report.specificity > 0.0 && report.specificity < 1.0);
    }

    #[test]
    fn matches_legacy_enumerator_on_small_genome() {
        // Cross-check: the FM-index search must find the same set of
        // (chrom, start, reverse, mismatches) as the legacy O(N×L)
        // enumerator in `valenx-genomics` on a small genome with planted
        // off-targets.
        let mut chrom = vec![b'A'; 800];
        let guide = b"GCGTACGTACGTACGTACGT";
        chrom[100..120].copy_from_slice(guide);
        chrom[120..123].copy_from_slice(b"CGG");
        // 1 mm site.
        let mm1 = b"GCATACGTACGTACGTACGT";
        chrom[300..320].copy_from_slice(mm1);
        chrom[320..323].copy_from_slice(b"AGG");
        // 2 mm site (different positions).
        let mm2 = b"GCATATGTACGTACGTACGT";
        chrom[500..520].copy_from_slice(mm2);
        chrom[520..523].copy_from_slice(b"TGG");
        let genome_vec = vec![("chr1".to_string(), chrom.clone())];
        let genome = GenomeIndex::build(&genome_vec).unwrap();
        let pam = pam_spcas9();

        let fm = find_off_targets_genome(guide, &genome, &pam, 3).unwrap();
        let legacy =
            valenx_genomics::crispr::offtarget::enumerate_off_targets(guide, &genome_vec, &pam, 3)
                .unwrap();
        // Compare as sets of (chrom, start, reverse, mismatches).
        let fm_set: HashSet<(String, usize, bool, usize)> = fm
            .iter()
            .map(|h| (h.chrom.clone(), h.start, h.reverse, h.mismatches))
            .collect();
        let legacy_set: HashSet<(String, usize, bool, usize)> = legacy
            .iter()
            .map(|h| (h.chrom.clone(), h.start, h.reverse, h.mismatches))
            .collect();
        assert_eq!(
            fm_set, legacy_set,
            "FM-index search must agree with the legacy enumerator on a small genome"
        );
    }

    #[test]
    fn empty_genome_returns_no_hits() {
        let g = GenomeIndex::new();
        let hits = find_off_targets_genome(b"ACGTACGTACGTACGTACGT", &g, &pam_spcas9(), 3).unwrap();
        assert!(hits.is_empty());
    }

    #[test]
    fn rejects_non_acgt_guide() {
        let g = GenomeIndex::build(&[("chr1".to_string(), b"ACGT".repeat(50))]).unwrap();
        let err =
            find_off_targets_genome(b"ACGTACGTNNNNNNNNACGT", &g, &pam_spcas9(), 3).unwrap_err();
        assert_eq!(err.code(), "genediting.invalid_target");
    }

    #[test]
    fn rejects_wrong_length_guide() {
        let g = GenomeIndex::build(&[("chr1".to_string(), b"ACGT".repeat(50))]).unwrap();
        let err = find_off_targets_genome(b"ACGTACGT", &g, &pam_spcas9(), 1).unwrap_err();
        assert_eq!(err.code(), "genediting.invalid");
    }

    #[test]
    fn genome_index_book_keeping() {
        let g = GenomeIndex::build(&[
            ("chr1".to_string(), b"ACGT".repeat(50)),
            ("chr2".to_string(), b"GCTA".repeat(30)),
        ])
        .unwrap();
        assert_eq!(g.contig_count(), 2);
        assert_eq!(g.total_length(), 200 + 120);
        assert_eq!(g.contigs()[0].name, "chr1");
        assert_eq!(g.contigs()[1].len(), 120);
    }

    #[test]
    fn rejects_non_acgt_contig() {
        let err = IndexedContig::build("chr1", b"ACGTNNNN").unwrap_err();
        assert_eq!(err.code(), "genediting.invalid_target");
    }
}
