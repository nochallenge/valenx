//! Read mapper — index a reference, map short reads, emit SAM.
//!
//! Closes the loop on the alignment crate by tying together its
//! pieces into a working short-read aligner in the BWA-MEM /
//! minimap2 class:
//!
//! 1. **Index** — for each reference sequence, build an
//!    [`crate::search::FmIndex`] *and* a minimizer index. Both halves
//!    of the seeding pipeline (SMEMs from the FM-index, minimizer
//!    anchors against the minimizer index) are stored.
//! 2. **Seed** — for each read, collect SMEMs (super-maximal exact
//!    matches via FM-index backward search) and minimizer hits;
//!    convert each occurrence into an [`Anchor`] tagged with its
//!    reference id and strand.
//! 3. **Chain** — group anchors per `(reference, strand)` and run a
//!    minimap2-style colinear DP chainer ([`chain_anchors`]) to
//!    produce candidate chains.
//! 4. **Extend** — for each top-scoring chain, run banded affine-gap
//!    DP (Gotoh, banded) around the chained region to recover an
//!    exact CIGAR. Smith-Waterman is used to trim local clips so the
//!    reported alignment is genuinely the best local hit, not the
//!    full extension window.
//! 5. **MAPQ** — derive mapping quality from the best vs.
//!    second-best chain/alignment scores (the BWA-MEM rule):
//!    `MAPQ ≈ 60 · (1 − S2/S1)`, clamped to `[0, 60]`. Unique
//!    placements get high MAPQ, ambiguous placements get low.
//! 6. **Paired-end** — for a pair of reads, place each mate
//!    independently, then *rescue* the partner: among the candidate
//!    placements of each mate, pick the pair that maximises
//!    `S₁ + S₂ + insert_size_bonus(d)` where `d` is the absolute
//!    insert size and the bonus is the negative-square-distance log-
//!    likelihood of `d` under a Normal `(mean, sd)` insert-size model.
//!
//! The output is a stream of [`crate::util::sam::SamRecord`] records;
//! `map_to_sam` renders a SAM file body with `@HD`/`@SQ` headers and
//! one alignment line per read (or per mate, for paired ends).

use crate::error::Result;
use crate::matrix::ScoringScheme;
use crate::pairwise::banded::banded_affine;
use crate::pairwise::local::smith_waterman;
use crate::pairwise::result::Cigar;
use crate::search::chain::{chain_anchors, Anchor, Chain, ChainParams};
use crate::search::fmindex::FmIndex;
use crate::search::minimizer::minimizer_sketch;
use crate::util::sam::SamRecord;
use std::collections::HashMap;
use valenx_bioseq::ops::revcomp::reverse_complement_dna_bytes;

/// SAM flag bit: read paired.
const FLAG_PAIRED: u16 = 0x1;
/// SAM flag bit: each segment properly aligned (i.e. pair within
/// insert-size expectations).
const FLAG_PROPER_PAIR: u16 = 0x2;
/// SAM flag bit: mate is unmapped.
const FLAG_MATE_UNMAPPED: u16 = 0x8;
/// SAM flag bit: query mapped to reverse strand.
const FLAG_REVERSE: u16 = 0x10;
/// SAM flag bit: mate mapped to reverse strand.
const FLAG_MATE_REVERSE: u16 = 0x20;
/// SAM flag bit: first segment in template.
const FLAG_FIRST: u16 = 0x40;
/// SAM flag bit: last segment in template.
const FLAG_LAST: u16 = 0x80;

/// Strand of a placement.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Strand {
    /// Forward strand (read aligns directly to the reference).
    Forward,
    /// Reverse strand (read is reverse-complemented before aligning).
    Reverse,
}

/// A reference sequence registered with a [`ReadMapper`].
#[derive(Clone, Debug)]
struct Reference {
    name: String,
    seq: Vec<u8>,
    fm: FmIndex,
}

/// Tunable parameters for [`ReadMapper`].
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct MapperParams {
    /// Minimum SMEM length to keep a seed (BWA-MEM default ~ 19).
    pub min_seed_len: usize,
    /// Minimizer k-mer length used for the secondary seed index.
    pub minimizer_k: usize,
    /// Minimizer window size.
    pub minimizer_w: usize,
    /// Anchor length used when synthesising minimizer-hit anchors.
    pub minimizer_anchor_len: usize,
    /// Half-bandwidth of the extension banded-affine DP, in residues.
    pub extension_band: usize,
    /// Maximum chains per `(reference, strand)` to extend before
    /// settling on a best alignment.
    pub max_chains_per_orientation: usize,
    /// Maximum number of hits to keep per SMEM (drops SMEMs that
    /// occur too often, the BWA `-c` heuristic).
    pub max_occurrences_per_seed: usize,
    /// Anchor weight gap penalty for [`chain_anchors`].
    pub chain_gap_weight: f64,
    /// Maximum residue gap between consecutive anchors in a chain.
    pub chain_max_gap: usize,
    /// Minimum *raw alignment* score to report a placement.
    pub min_score: i32,
}

impl Default for MapperParams {
    fn default() -> Self {
        MapperParams {
            min_seed_len: 15,
            minimizer_k: 11,
            minimizer_w: 5,
            minimizer_anchor_len: 11,
            extension_band: 25,
            max_chains_per_orientation: 4,
            max_occurrences_per_seed: 200,
            chain_gap_weight: 0.5,
            chain_max_gap: 5_000,
            min_score: 20,
        }
    }
}

/// A short-read mapper over a fixed set of reference sequences.
#[derive(Clone, Debug)]
pub struct ReadMapper {
    references: Vec<Reference>,
    /// A minimizer index keyed by `(seq_id, hash) -> Vec<positions>`.
    /// Built once at construction; minimap2's `mm_idx_t` analogue.
    minimizer_index: HashMap<u64, Vec<(usize, usize)>>,
    scheme: ScoringScheme,
    params: MapperParams,
}

/// The outcome of mapping one read.
#[derive(Clone, Debug, PartialEq)]
pub struct MappingResult {
    /// The SAM record (mapped or unmapped).
    pub record: SamRecord,
    /// The alignment score of the best placement (`0` if unmapped).
    pub score: i32,
    /// Strand the read was placed on.
    pub strand: Strand,
}

/// Insert-size model for paired-end mapping. A Normal distribution
/// of insert sizes scored as `−((d − mean)² / (2 sd²))` (the
/// log-density up to additive constants), then scaled to a small
/// integer bonus comparable to alignment scores.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct InsertSizeModel {
    /// Expected insert size in base pairs.
    pub mean: f64,
    /// Standard deviation in base pairs.
    pub sd: f64,
    /// Multiplier converting the log-density into score units.
    pub scale: f64,
    /// Maximum bonus (subtracted from the absolute score floor of `0`)
    /// awarded to a perfectly-sized pair.
    pub max_bonus: i32,
}

impl Default for InsertSizeModel {
    fn default() -> Self {
        InsertSizeModel {
            mean: 300.0,
            sd: 100.0,
            scale: 1.0,
            max_bonus: 20,
        }
    }
}

impl InsertSizeModel {
    /// Score bonus (integer, in alignment-score units) for an absolute
    /// insert size `d`. Capped at `max_bonus`; never negative.
    pub fn bonus(&self, d: f64) -> i32 {
        let z = (d - self.mean) / self.sd;
        // log density of normal, dropping the constant: −z²/2.
        let lp = -0.5 * z * z;
        // lp <= 0; scale to [0, max_bonus].
        let scaled = (self.scale * (self.max_bonus as f64) * lp.exp()).round() as i32;
        scaled.clamp(0, self.max_bonus)
    }
}

/// A paired-end mapping outcome.
#[derive(Clone, Debug, PartialEq)]
pub struct PairedMappingResult {
    /// Mate 1's record.
    pub mate1: SamRecord,
    /// Mate 2's record.
    pub mate2: SamRecord,
    /// Raw insert size (mate2 POS − mate1 POS + mate2 length), 0 if
    /// either is unmapped or they aren't on the same reference.
    pub insert_size: i64,
}

/// Internal candidate placement (one chain extended to base level).
#[derive(Clone, Debug)]
struct Candidate {
    seq_id: usize,
    strand: Strand,
    ref_pos: usize, // 0-based on the forward reference
    score: i32,
    cigar: Cigar,
}

impl ReadMapper {
    /// Builds a mapper indexing `references`. Each reference gets its
    /// own FM-index; minimizers are pooled into a single hash table
    /// keyed by `(seq_id, position)`.
    pub fn new(
        references: &[(&str, &[u8])],
        scheme: ScoringScheme,
        params: MapperParams,
    ) -> Result<Self> {
        let mut refs: Vec<Reference> = Vec::with_capacity(references.len());
        for (name, seq) in references {
            let fm = FmIndex::build(seq)?;
            refs.push(Reference {
                name: (*name).to_string(),
                seq: seq.to_vec(),
                fm,
            });
        }
        // Build a pooled minimizer index over the forward strand of
        // each reference. Long reads / repetitive minimizers create
        // chains; we use both this and SMEMs as seeds.
        let mut minimizer_index: HashMap<u64, Vec<(usize, usize)>> = HashMap::new();
        for (i, r) in refs.iter().enumerate() {
            if let Ok(sk) = minimizer_sketch(&r.seq, params.minimizer_k, params.minimizer_w) {
                for m in sk {
                    minimizer_index.entry(m.hash).or_default().push((i, m.pos));
                }
            }
        }
        Ok(ReadMapper {
            references: refs,
            minimizer_index,
            scheme,
            params,
        })
    }

    /// Construct with default [`MapperParams`].
    pub fn with_defaults(references: &[(&str, &[u8])], scheme: ScoringScheme) -> Result<Self> {
        Self::new(references, scheme, MapperParams::default())
    }

    /// The number of indexed reference sequences.
    pub fn reference_count(&self) -> usize {
        self.references.len()
    }

    /// Name of the i-th reference.
    pub fn reference_name(&self, i: usize) -> Option<&str> {
        self.references.get(i).map(|r| r.name.as_str())
    }

    /// Length of the i-th reference.
    pub fn reference_len(&self, i: usize) -> Option<usize> {
        self.references.get(i).map(|r| r.seq.len())
    }

    /// Maps one read against the reference set on both strands,
    /// returning its best placement.
    pub fn map_read(&self, read_name: &str, read: &[u8]) -> MappingResult {
        let read_rc = reverse_complement_dna_bytes(read);

        // Collect candidate placements from both strands.
        let mut candidates: Vec<Candidate> = Vec::new();
        candidates.extend(self.candidates_for_strand(read, Strand::Forward));
        candidates.extend(self.candidates_for_strand(&read_rc, Strand::Reverse));
        if candidates.is_empty() {
            return MappingResult {
                record: SamRecord::unmapped(read_name, read.to_vec()),
                score: 0,
                strand: Strand::Forward,
            };
        }

        candidates.sort_by(|a, b| b.score.cmp(&a.score));
        let best = candidates[0].clone();
        // The "second best" must be at a *different* placement —
        // ignore candidates that overlap the best one's reference
        // span by more than half. Two chains over the same window
        // are not competing placements.
        let second = second_best_distinct(&candidates, &best);

        if best.score < self.params.min_score {
            return MappingResult {
                record: SamRecord::unmapped(read_name, read.to_vec()),
                score: 0,
                strand: Strand::Forward,
            };
        }

        let mapq = mapping_quality(best.score, second);
        let mut record = SamRecord::mapped(
            read_name,
            &self.references[best.seq_id].name,
            best.ref_pos,
            mapq,
            best.cigar,
            // SAM SEQ is on the forward strand of the read; for a
            // reverse placement we report the reverse-complement and
            // set the FLAG bit. minimap2/BWA convention.
            if best.strand == Strand::Reverse {
                read_rc.clone()
            } else {
                read.to_vec()
            },
        );
        if best.strand == Strand::Reverse {
            record.flag |= FLAG_REVERSE;
        }
        MappingResult {
            record,
            score: best.score,
            strand: best.strand,
        }
    }

    /// Maps a batch of `(name, read)` pairs, returning one
    /// [`MappingResult`] each in input order.
    pub fn map_reads(&self, reads: &[(&str, &[u8])]) -> Vec<MappingResult> {
        reads
            .iter()
            .map(|(name, read)| self.map_read(name, read))
            .collect()
    }

    /// Map a pair of reads. Both mates are mapped independently, then
    /// the best *consistent* pair (same reference, opposite strands,
    /// insert size within `insert.mean ± several · sd`) is selected;
    /// if no consistent pair exists each mate keeps its best single-
    /// end placement.
    pub fn map_pair(
        &self,
        pair_name: &str,
        read1: &[u8],
        read2: &[u8],
        insert: InsertSizeModel,
    ) -> PairedMappingResult {
        let read1_rc = reverse_complement_dna_bytes(read1);
        let read2_rc = reverse_complement_dna_bytes(read2);

        let cands1_fwd: Vec<Candidate> = self.candidates_for_strand(read1, Strand::Forward);
        let cands1_rev: Vec<Candidate> = self.candidates_for_strand(&read1_rc, Strand::Reverse);
        let cands2_fwd: Vec<Candidate> = self.candidates_for_strand(read2, Strand::Forward);
        let cands2_rev: Vec<Candidate> = self.candidates_for_strand(&read2_rc, Strand::Reverse);

        let mut all1: Vec<Candidate> = cands1_fwd.clone();
        all1.extend(cands1_rev.clone());
        all1.sort_by(|a, b| b.score.cmp(&a.score));
        let mut all2: Vec<Candidate> = cands2_fwd.clone();
        all2.extend(cands2_rev.clone());
        all2.sort_by(|a, b| b.score.cmp(&a.score));

        // Try to find a consistent pair: same reference, opposite
        // strands, insert size positive and reasonable. Score the pair
        // by sum of mate scores plus the insert-size bonus.
        let mut best_pair: Option<(Candidate, Candidate, i32, i64)> = None;
        let take = self.params.max_chains_per_orientation * 4;
        for c1 in all1.iter().take(take) {
            for c2 in all2.iter().take(take) {
                if c1.seq_id != c2.seq_id || c1.strand == c2.strand {
                    continue;
                }
                // Insert size: distance between the two outer endpoints
                // on the reference. Standard convention: signed,
                // positive when mate2 lies to the right of mate1.
                let r1_start = c1.ref_pos as i64;
                let r1_end = (c1.ref_pos + c1.cigar.ref_len()) as i64;
                let r2_start = c2.ref_pos as i64;
                let r2_end = (c2.ref_pos + c2.cigar.ref_len()) as i64;
                let isize_signed = if c1.strand == Strand::Forward {
                    // Mate1 forward => mate2 should be downstream and reverse.
                    r2_end - r1_start
                } else {
                    -(r1_end - r2_start)
                };
                let isize_abs = isize_signed.unsigned_abs() as f64;
                let bonus = insert.bonus(isize_abs);
                let combined = c1.score + c2.score + bonus;
                if best_pair.as_ref().is_none_or(|p| combined > p.2) {
                    best_pair = Some((c1.clone(), c2.clone(), combined, isize_signed));
                }
            }
        }

        let make_record = |cand: Option<&Candidate>,
                           read_fwd: &[u8],
                           read_rc: &[u8],
                           name: &str,
                           base_flag: u16,
                           mate_unmapped: bool,
                           mate_reverse: bool|
         -> (SamRecord, i32, Strand) {
            match cand {
                Some(c) => {
                    let mut flag = base_flag | FLAG_PAIRED;
                    if c.strand == Strand::Reverse {
                        flag |= FLAG_REVERSE;
                    }
                    if mate_unmapped {
                        flag |= FLAG_MATE_UNMAPPED;
                    }
                    if mate_reverse {
                        flag |= FLAG_MATE_REVERSE;
                    }
                    let seq = if c.strand == Strand::Reverse {
                        read_rc.to_vec()
                    } else {
                        read_fwd.to_vec()
                    };
                    let mut rec = SamRecord::mapped(
                        name,
                        &self.references[c.seq_id].name,
                        c.ref_pos,
                        // MAPQ filled in below from best/second
                        0,
                        c.cigar.clone(),
                        seq,
                    );
                    rec.flag |= flag;
                    (rec, c.score, c.strand)
                }
                None => {
                    let mut rec = SamRecord::unmapped(name, read_fwd.to_vec());
                    rec.flag |= base_flag | FLAG_PAIRED;
                    if mate_reverse {
                        rec.flag |= FLAG_MATE_REVERSE;
                    }
                    (rec, 0, Strand::Forward)
                }
            }
        };

        match best_pair {
            Some((c1, c2, _combined, isize_signed)) => {
                let second1 = second_best_distinct(&all1, &c1);
                let second2 = second_best_distinct(&all2, &c2);
                let mapq1 = mapping_quality(c1.score, second1);
                let mapq2 = mapping_quality(c2.score, second2);
                let (mut r1, _, _) = make_record(
                    Some(&c1),
                    read1,
                    &read1_rc,
                    pair_name,
                    FLAG_FIRST,
                    false,
                    c2.strand == Strand::Reverse,
                );
                let (mut r2, _, _) = make_record(
                    Some(&c2),
                    read2,
                    &read2_rc,
                    pair_name,
                    FLAG_LAST,
                    false,
                    c1.strand == Strand::Reverse,
                );
                r1.mapq = mapq1;
                r2.mapq = mapq2;
                r1.flag |= FLAG_PROPER_PAIR;
                r2.flag |= FLAG_PROPER_PAIR;
                PairedMappingResult {
                    mate1: r1,
                    mate2: r2,
                    insert_size: isize_signed,
                }
            }
            None => {
                // No consistent pair — fall back to independent best
                // single-end placements.
                let best1 = all1.first().cloned();
                let best2 = all2.first().cloned();
                let second1 = best1
                    .as_ref()
                    .map(|b| second_best_distinct(&all1, b))
                    .unwrap_or(i32::MIN);
                let second2 = best2
                    .as_ref()
                    .map(|b| second_best_distinct(&all2, b))
                    .unwrap_or(i32::MIN);
                let (mut r1, _, _) = make_record(
                    best1.as_ref(),
                    read1,
                    &read1_rc,
                    pair_name,
                    FLAG_FIRST,
                    best2.is_none(),
                    best2
                        .as_ref()
                        .map(|c| c.strand == Strand::Reverse)
                        .unwrap_or(false),
                );
                let (mut r2, _, _) = make_record(
                    best2.as_ref(),
                    read2,
                    &read2_rc,
                    pair_name,
                    FLAG_LAST,
                    best1.is_none(),
                    best1
                        .as_ref()
                        .map(|c| c.strand == Strand::Reverse)
                        .unwrap_or(false),
                );
                if let Some(ref c) = best1 {
                    r1.mapq = mapping_quality(c.score, second1);
                }
                if let Some(ref c) = best2 {
                    r2.mapq = mapping_quality(c.score, second2);
                }
                PairedMappingResult {
                    mate1: r1,
                    mate2: r2,
                    insert_size: 0,
                }
            }
        }
    }

    /// Produce candidate placements for the read on a single strand.
    fn candidates_for_strand(&self, read: &[u8], strand: Strand) -> Vec<Candidate> {
        let mut out: Vec<Candidate> = Vec::new();

        // Per-reference anchor cloud. For each reference we collect
        // FM-index SMEM hits + minimizer hits, then chain.
        let mut anchors_by_ref: HashMap<usize, Vec<Anchor>> = HashMap::new();

        // --- SMEM anchors (FM-index per reference) ----------------
        for (seq_id, r) in self.references.iter().enumerate() {
            let smems = r.fm.smems(read, self.params.min_seed_len);
            for s in smems {
                if s.count > self.params.max_occurrences_per_seed {
                    continue;
                }
                for ref_pos in r.fm.smem_positions(&s) {
                    anchors_by_ref.entry(seq_id).or_default().push(Anchor::new(
                        s.query_start,
                        ref_pos,
                        s.len(),
                    ));
                }
            }
        }

        // --- Minimizer anchors -----------------------------------
        if let Ok(sk) = minimizer_sketch(read, self.params.minimizer_k, self.params.minimizer_w) {
            for m in sk {
                if let Some(refs) = self.minimizer_index.get(&m.hash) {
                    if refs.len() > self.params.max_occurrences_per_seed {
                        continue;
                    }
                    for &(seq_id, ref_pos) in refs {
                        anchors_by_ref.entry(seq_id).or_default().push(Anchor::new(
                            m.pos,
                            ref_pos,
                            self.params.minimizer_anchor_len,
                        ));
                    }
                }
            }
        }

        let chain_params = ChainParams {
            gap_weight: self.params.chain_gap_weight,
            max_gap: self.params.chain_max_gap,
        };

        for (seq_id, anchors) in anchors_by_ref {
            if anchors.is_empty() {
                continue;
            }
            // Run repeated chaining: take the best chain, drop its
            // anchors, run again. This gives `max_chains_per_orientation`
            // distinct chains per reference (minimap2's "chain
            // pruning" — simpler than the score-decay heuristic).
            let mut remaining = anchors;
            for _ in 0..self.params.max_chains_per_orientation {
                if remaining.is_empty() {
                    break;
                }
                let chain = chain_anchors(&remaining, chain_params);
                if chain.is_empty() {
                    break;
                }
                if let Some(cand) = self.extend_chain(read, seq_id, strand, &chain) {
                    out.push(cand);
                }
                // Drop the anchors that participated in this chain so
                // the next round finds a different one.
                let used: std::collections::HashSet<Anchor> =
                    chain.anchors.iter().copied().collect();
                remaining.retain(|a| !used.contains(a));
            }
        }

        out
    }

    /// Base-level extension of one chain into a [`Candidate`].
    /// Runs banded affine alignment over a windowed slice of the
    /// reference around the chain's target span, then trims the local
    /// best-scoring substring via Smith-Waterman to drop poor flanks.
    fn extend_chain(
        &self,
        read: &[u8],
        seq_id: usize,
        strand: Strand,
        chain: &Chain,
    ) -> Option<Candidate> {
        let r = &self.references[seq_id];
        let (q_lo, q_hi) = chain.query_span();
        let (t_lo, t_hi) = chain.target_span();
        // Diagonal of the first anchor — use it to anchor the window
        // even when q_lo > 0.
        let diag = if let Some(first) = chain.anchors.first() {
            first.target_pos as isize - first.query_pos as isize
        } else {
            (t_lo as isize) - (q_lo as isize)
        };

        // Predict the reference region that should contain the whole
        // read: extend the chain's target span by enough to cover the
        // unaligned read prefix/suffix plus a margin.
        let pad = self.params.extension_band.max(16) + read.len() / 4;
        let ext_start = (diag - pad as isize).max(0) as usize;
        let ext_end =
            (((read.len() as isize + diag) + pad as isize).max(0) as usize).min(r.seq.len());
        let _ = (q_hi, t_hi);
        if ext_start >= ext_end {
            return None;
        }
        let window = &r.seq[ext_start..ext_end];

        // Run local Smith-Waterman to find the best-scoring trimmed
        // alignment of the read within this window. SW returns
        // `span1` (in the read) and `span2` (in the window).
        let al = match smith_waterman(read, window, &self.scheme) {
            Ok(a) => a,
            Err(_) => return None,
        };
        if al.is_empty() || al.score < self.params.min_score {
            return None;
        }
        // Refine with a banded affine pass over the aligned core to
        // get a globally-optimal alignment in the band — exactness for
        // BWA-MEM-class CIGARs.
        let read_lo = al.span1.0;
        let read_hi = al.span1.1;
        let win_lo = al.span2.0;
        let win_hi = al.span2.1;
        let ref_pos = ext_start + win_lo;
        // Re-align in a wider band if possible.
        let trimmed_read = &read[read_lo..read_hi];
        let trimmed_ref = &window[win_lo..win_hi];
        let band = self
            .params
            .extension_band
            .max((trimmed_read.len() as isize - trimmed_ref.len() as isize).unsigned_abs());
        let cigar = match banded_affine(trimmed_read, trimmed_ref, &self.scheme, band) {
            Ok(refined) => refined.cigar(),
            Err(_) => al.cigar(),
        };
        Some(Candidate {
            seq_id,
            strand,
            ref_pos,
            score: al.score,
            cigar,
        })
    }

    /// Maps a batch and renders the results as a SAM file body
    /// (`@SQ` header lines for each reference, then one alignment line
    /// per read).
    pub fn map_to_sam(&self, reads: &[(&str, &[u8])]) -> String {
        let mut out = String::from("@HD\tVN:1.6\tSO:unsorted\n");
        for r in &self.references {
            out.push_str(&format!("@SQ\tSN:{}\tLN:{}\n", r.name, r.seq.len()));
        }
        for res in self.map_reads(reads) {
            out.push_str(&res.record.to_sam_line());
            out.push('\n');
        }
        out
    }

    /// Renders paired-end results as a SAM file body. Each pair
    /// contributes two alignment lines (one per mate).
    pub fn map_pairs_to_sam(
        &self,
        pairs: &[(&str, &[u8], &[u8])],
        insert: InsertSizeModel,
    ) -> String {
        let mut out = String::from("@HD\tVN:1.6\tSO:unsorted\n");
        for r in &self.references {
            out.push_str(&format!("@SQ\tSN:{}\tLN:{}\n", r.name, r.seq.len()));
        }
        for (name, r1, r2) in pairs {
            let p = self.map_pair(name, r1, r2, insert);
            out.push_str(&p.mate1.to_sam_line());
            out.push('\n');
            out.push_str(&p.mate2.to_sam_line());
            out.push('\n');
        }
        out
    }
}

/// `true` when two [`Candidate`]s refer to the same placement (same
/// reference, strand, and reference position) — used to skip the
/// best when scanning for the "second best".
fn same_placement(a: &Candidate, b: &Candidate) -> bool {
    a.seq_id == b.seq_id && a.strand == b.strand && a.ref_pos == b.ref_pos
}

/// Return the highest-scoring candidate that is *not* the same
/// placement as `best` (different reference, or strand, or reference
/// position differing by more than half of `best`'s ref span). When
/// no such candidate exists, returns `i32::MIN` so [`mapping_quality`]
/// treats `best` as uncontested.
fn second_best_distinct(candidates: &[Candidate], best: &Candidate) -> i32 {
    let half = (best.cigar.ref_len().max(1) / 2) as i64;
    for c in candidates.iter() {
        if same_placement(c, best) {
            continue;
        }
        if c.seq_id == best.seq_id && c.strand == best.strand {
            let d = (c.ref_pos as i64) - (best.ref_pos as i64);
            if d.abs() <= half {
                // Same placement, just a slightly different chain
                // recovered the same alignment — skip.
                continue;
            }
        }
        return c.score;
    }
    i32::MIN
}

/// A heuristic mapping quality in `[0, 60]` — the BWA-MEM rule.
///
/// `MAPQ = 60 · (1 − S2/S1)` clamped to `[0, 60]`, with the special
/// case that a placement with no competitor (`second == i32::MIN`)
/// gets the maximum `60` and any negative best score gets `0`.
pub fn mapping_quality(best: i32, second: i32) -> u8 {
    if best <= 0 {
        return 0;
    }
    if second == i32::MIN {
        return 60;
    }
    let s1 = best as f64;
    let s2 = second.max(0) as f64;
    let frac = 1.0 - (s2 / s1);
    let mapq = (60.0 * frac).round().clamp(0.0, 60.0);
    mapq as u8
}

// =====================================================================
// Back-compat re-exports — historically the mapper was a v1 keyed on a
// KmerIndex. The new implementation keeps no API for that v1 (it was
// always heuristic), but a few callers used the constructor signature
// `ReadMapper::new(refs, k, scheme)`; we *replace* that with the new
// signature. The single-end `map_read(name, read, min_score)` lives on
// via `map_read_with_floor` to support callers that pre-supply a
// minimum-score floor; the default `map_read` uses `params.min_score`.
// =====================================================================

impl ReadMapper {
    /// Map one read with an explicit minimum-score floor (overrides
    /// the configured `params.min_score`). Otherwise identical to
    /// [`map_read`](ReadMapper::map_read).
    pub fn map_read_with_floor(
        &self,
        read_name: &str,
        read: &[u8],
        min_score: i32,
    ) -> MappingResult {
        let saved = self.params.min_score;
        let mut mapper = self.clone();
        mapper.params.min_score = min_score;
        let mut res = mapper.map_read(read_name, read);
        mapper.params.min_score = saved;
        // Belt-and-braces: if the floor demoted it to unmapped at the
        // very end, mark as such.
        if res.score < min_score {
            res.record = SamRecord::unmapped(read_name, read.to_vec());
            res.score = 0;
        }
        res
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::{GapCost, ScoringScheme, SubstitutionMatrix};

    /// A 200 bp pseudo-random DNA reference, deterministic.
    fn reference_200() -> Vec<u8> {
        let mut out = Vec::with_capacity(200);
        let alphabet = b"ACGT";
        let mut h: u32 = 0xdeadbeef;
        for _ in 0..200 {
            h ^= h << 13;
            h ^= h >> 17;
            h ^= h << 5;
            out.push(alphabet[(h as usize) % 4]);
        }
        out
    }

    fn scheme() -> ScoringScheme {
        ScoringScheme::new(SubstitutionMatrix::dna_simple(2, -3), GapCost::new(5, 2))
    }

    // (No factory helper — each test builds its own owned reference
    // and mapper; ReadMapper owns its references via Vec<u8>, so once
    // built it does not borrow caller-owned buffers.)

    #[test]
    fn maps_exact_read_to_reference() {
        let r = reference_200();
        let refs: &[(&str, &[u8])] = &[("chr1", &r)];
        let m = ReadMapper::with_defaults(refs, scheme()).unwrap();
        // Take a 50 bp window at offset 40 — should map back exactly.
        let read = &r[40..90];
        let res = m.map_read("r1", read);
        assert!(!res.record.is_unmapped(), "exact read must map");
        assert_eq!(res.record.rname, "chr1");
        assert_eq!(res.record.pos, 41, "0-based 40 should map to 1-based 41");
        assert_eq!(res.strand, Strand::Forward);
        // CIGAR should be perfect 50M.
        let cig = res.record.cigar.to_string();
        assert!(
            cig.contains("50M") || cig == "50M",
            "expected 50M, got {cig}"
        );
    }

    #[test]
    fn maps_read_with_substitutions() {
        let r = reference_200();
        let refs: &[(&str, &[u8])] = &[("chr1", &r)];
        let m = ReadMapper::with_defaults(refs, scheme()).unwrap();
        // 60 bp window at offset 70, mutate a couple of positions.
        let mut read = r[70..130].to_vec();
        for pos in [10usize, 30] {
            read[pos] = if read[pos] == b'A' { b'C' } else { b'A' };
        }
        let res = m.map_read("r1", &read);
        assert!(!res.record.is_unmapped(), "2-mismatch read should map");
        assert_eq!(res.record.rname, "chr1");
        // The placement should be at or very near offset 70.
        let pos = res.record.pos as i64;
        assert!((pos - 71).abs() <= 1, "expected 1-based ~71, got {pos}");
    }

    #[test]
    fn maps_read_with_indel() {
        let r = reference_200();
        let refs: &[(&str, &[u8])] = &[("chr1", &r)];
        let m = ReadMapper::with_defaults(refs, scheme()).unwrap();
        // 60 bp window, delete 2 bases from the middle.
        let mut read = r[50..110].to_vec();
        read.drain(28..30);
        let res = m.map_read("r1", &read);
        assert!(!res.record.is_unmapped(), "1-indel read should map");
        assert_eq!(res.record.rname, "chr1");
        // CIGAR must consume more reference than query — i.e. a
        // deletion appears.
        let cig = &res.record.cigar;
        assert!(
            cig.ref_len() > cig.query_len(),
            "CIGAR must contain a deletion"
        );
    }

    #[test]
    fn reverse_strand_read_maps_reverse() {
        let r = reference_200();
        let refs: &[(&str, &[u8])] = &[("chr1", &r)];
        let m = ReadMapper::with_defaults(refs, scheme()).unwrap();
        let window = &r[80..140];
        let read_rc = reverse_complement_dna_bytes(window);
        let res = m.map_read("r1", &read_rc);
        assert!(!res.record.is_unmapped(), "reverse-strand read should map");
        assert_eq!(res.strand, Strand::Reverse);
        assert_ne!(
            res.record.flag & FLAG_REVERSE,
            0,
            "FLAG must have REVERSE bit"
        );
        // POS reports the forward-strand reference position; for our
        // 60 bp window at offset 80, that's 1-based 81.
        let pos = res.record.pos as i64;
        assert!((pos - 81).abs() <= 1, "expected 1-based ~81, got {pos}");
    }

    #[test]
    fn unmappable_read_reported_unmapped() {
        let refs: &[(&str, &[u8])] =
            &[("chr1", b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")];
        let m = ReadMapper::with_defaults(refs, scheme()).unwrap();
        let res = m.map_read("r1", b"CCCCCCCCCCCCCCCCCCCCCCCC");
        assert!(res.record.is_unmapped());
        assert_eq!(res.score, 0);
    }

    #[test]
    fn unique_placement_gets_high_mapq() {
        let r = reference_200();
        let refs: &[(&str, &[u8])] = &[("chr1", &r)];
        let m = ReadMapper::with_defaults(refs, scheme()).unwrap();
        // A long read from a unique part of the reference.
        let read = &r[20..100];
        let res = m.map_read("r1", read);
        assert!(!res.record.is_unmapped());
        assert!(
            res.record.mapq >= 40,
            "unique placement expected MAPQ>=40, got {}",
            res.record.mapq
        );
    }

    #[test]
    fn repeat_placement_gets_lower_mapq() {
        // Two identical 30 bp windows in different references — any
        // 30 bp read from that window has two equally good placements,
        // so MAPQ must drop.
        let core = b"ACGTACGTACGTACGTACGTACGTACGTAC";
        let mut chr1 = vec![b'A'; 60];
        chr1.extend_from_slice(core);
        chr1.extend(std::iter::repeat_n(b'A', 60));
        let mut chr2 = vec![b'C'; 60];
        chr2.extend_from_slice(core);
        chr2.extend(std::iter::repeat_n(b'C', 60));
        let refs: &[(&str, &[u8])] = &[("chr1", &chr1), ("chr2", &chr2)];
        let m = ReadMapper::with_defaults(refs, scheme()).unwrap();
        let res = m.map_read("r1", core);
        assert!(!res.record.is_unmapped());
        assert!(
            res.record.mapq < 40,
            "ambiguous (2-copy) placement should have MAPQ < 40, got {}",
            res.record.mapq,
        );
    }

    #[test]
    fn paired_end_consistent_pair() {
        let r = reference_200();
        let refs: &[(&str, &[u8])] = &[("chr1", &r)];
        let m = ReadMapper::with_defaults(refs, scheme()).unwrap();
        // Mate1 forward at offset 30, mate2 reverse-complement at
        // offset 130 — insert size ~ 130 - 30 + 50 = 150.
        let mate1 = r[30..80].to_vec();
        let mate2_fwd = r[130..180].to_vec();
        let mate2 = reverse_complement_dna_bytes(&mate2_fwd);
        let insert = InsertSizeModel {
            mean: 150.0,
            sd: 30.0,
            ..InsertSizeModel::default()
        };
        let pair = m.map_pair("p1", &mate1, &mate2, insert);
        assert!(!pair.mate1.is_unmapped(), "mate1 must map");
        assert!(!pair.mate2.is_unmapped(), "mate2 must map");
        assert_eq!(pair.mate1.rname, "chr1");
        assert_eq!(pair.mate2.rname, "chr1");
        assert_ne!(pair.mate1.flag & FLAG_PAIRED, 0);
        assert_ne!(pair.mate1.flag & FLAG_FIRST, 0);
        assert_ne!(pair.mate2.flag & FLAG_LAST, 0);
        assert_ne!(
            pair.mate2.flag & FLAG_REVERSE,
            0,
            "mate2 should be on reverse"
        );
        assert_ne!(
            pair.mate1.flag & FLAG_PROPER_PAIR,
            0,
            "should be a proper pair"
        );
        assert!(
            pair.insert_size.abs() > 100 && pair.insert_size.abs() < 250,
            "insert size {} should be ~150",
            pair.insert_size,
        );
    }

    #[test]
    fn paired_end_unmappable_fallback() {
        let r = reference_200();
        let refs: &[(&str, &[u8])] = &[("chr1", &r)];
        let m = ReadMapper::with_defaults(refs, scheme()).unwrap();
        let mate1 = r[30..80].to_vec();
        let mate2 = b"NNNNNNNNNNNNNNNNNNNNNNNNNN".to_vec();
        let pair = m.map_pair("p1", &mate1, &mate2, InsertSizeModel::default());
        assert!(!pair.mate1.is_unmapped(), "mate1 maps");
        assert!(pair.mate2.is_unmapped(), "mate2 unmappable");
        assert_ne!(pair.mate1.flag & FLAG_MATE_UNMAPPED, 0);
    }

    #[test]
    fn cross_check_against_smith_waterman_on_small_case() {
        // Build a tiny reference with one obvious 20-bp homology
        // region; check that the mapper's reported alignment score
        // equals the score of a plain Smith-Waterman on the same
        // window. Verifies the mapper's scoring isn't off.
        let r = b"AAAAAAAAAACTTGTTAACGGTCTAAACAAAAAAAAAA";
        let refs: &[(&str, &[u8])] = &[("chr1", r)];
        let m = ReadMapper::with_defaults(refs, scheme()).unwrap();
        let read = &r[10..30];
        let res = m.map_read("r1", read);
        assert!(!res.record.is_unmapped());
        let sw = smith_waterman(read, r, &scheme()).unwrap();
        // The mapper picks the best Smith-Waterman placement of the
        // read in the window; on an exact substring the scores must
        // match the unconstrained Smith-Waterman score.
        assert_eq!(res.score, sw.score);
    }

    #[test]
    fn batch_to_sam_has_header_and_records() {
        let r = reference_200();
        let refs: &[(&str, &[u8])] = &[("chr1", &r)];
        let m = ReadMapper::with_defaults(refs, scheme()).unwrap();
        let read = r[20..80].to_vec();
        let reads: &[(&str, &[u8])] = &[("r1", &read), ("r2", b"ZZZZZZZZ".as_slice())];
        let sam = m.map_to_sam(reads);
        assert!(sam.contains("@HD"));
        assert!(sam.contains("@SQ\tSN:chr1"));
        assert!(sam.contains("r1\t"));
        assert!(sam.contains("r2\t"));
    }

    #[test]
    fn mapping_quality_bounds() {
        assert_eq!(mapping_quality(100, i32::MIN), 60); // unique
        assert_eq!(mapping_quality(100, 100), 0); // tied -> ambiguous
        assert_eq!(mapping_quality(0, 0), 0); // floor
                                              // 100 vs 50 => 30. With our BWA-MEM rule, 60*(1 - 0.5) = 30.
        assert_eq!(mapping_quality(100, 50), 30);
    }

    #[test]
    fn picks_the_right_reference() {
        let r = reference_200();
        // chr2 is reference_200 reshuffled to break local matches.
        let mut r2 = r.clone();
        r2.reverse();
        let refs: &[(&str, &[u8])] = &[
            ("chrA", b"AAAAAAAAAAAAAAAAAAAA"),
            ("chrB", &r),
            ("chrC", &r2),
        ];
        let m = ReadMapper::with_defaults(refs, scheme()).unwrap();
        let read = &r[10..70];
        let res = m.map_read("r1", read);
        assert!(!res.record.is_unmapped());
        assert_eq!(res.record.rname, "chrB");
    }

    #[test]
    fn insert_size_bonus_is_reasonable() {
        let m = InsertSizeModel {
            mean: 300.0,
            sd: 50.0,
            scale: 1.0,
            max_bonus: 20,
        };
        // At the mean: full bonus.
        assert_eq!(m.bonus(300.0), 20);
        // 4 sd away: essentially zero.
        assert!(m.bonus(500.0) <= 1);
    }
}
