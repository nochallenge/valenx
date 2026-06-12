//! Seed-and-extend heuristic search (BLAST-class).
//!
//! The exhaustive Smith-Waterman scan of a query against a large
//! database is O(query · database) — too slow for genome-scale work.
//! The BLAST heuristic instead:
//!
//! 1. **Seeds** — finds short exact k-mer matches via a
//!    [`crate::search::kmer::KmerIndex`].
//! 2. **Diagonal binning** — groups seeds sharing the same
//!    `(seq_id, diagonal)`; nearby seeds on one diagonal hint at a
//!    real homology.
//! 3. **Extension** — extends each seed left and right with a
//!    banded gapless / small-gap dynamic program, using the
//!    X-drop rule to stop when the running score falls too far below
//!    the best seen.
//!
//! Each surviving high-scoring segment pair (HSP) is returned as an
//! [`Hsp`] with its score, coordinate spans, and (optionally) the
//! Karlin-Altschul E-value / bit-score from
//! [`crate::search::stats`].

use super::kmer::KmerIndex;
use super::stats::KarlinAltschul;
use crate::matrix::ScoringScheme;

/// A high-scoring segment pair found by seed-and-extend.
#[derive(Clone, Debug, PartialEq)]
pub struct Hsp {
    /// Index of the database sequence the HSP lies in.
    pub seq_id: usize,
    /// Half-open `[start, end)` span within the query.
    pub query_span: (usize, usize),
    /// Half-open `[start, end)` span within the database sequence.
    pub target_span: (usize, usize),
    /// Raw alignment score of the extended segment.
    pub score: i32,
    /// Bit score (Karlin-Altschul), populated by
    /// [`SeedSearch::with_stats`]; `None` otherwise.
    pub bit_score: Option<f64>,
    /// E-value (Karlin-Altschul), populated by
    /// [`SeedSearch::with_stats`]; `None` otherwise.
    pub e_value: Option<f64>,
}

/// Tunable parameters for [`SeedSearch`].
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct SeedParams {
    /// X-drop threshold: extension stops once the running score drops
    /// `x_drop` below the best score seen on that side.
    pub x_drop: i32,
    /// Minimum raw score for an HSP to be reported.
    pub min_score: i32,
    /// Maximum number of seeds to consider per `(seq_id, diagonal)`
    /// bin before declaring it a candidate (1 is enough to trigger).
    pub min_seeds_per_diagonal: usize,
}

impl Default for SeedParams {
    fn default() -> Self {
        SeedParams {
            x_drop: 20,
            min_score: 10,
            min_seeds_per_diagonal: 1,
        }
    }
}

/// A configured seed-and-extend searcher over a fixed database index.
#[derive(Clone, Debug)]
pub struct SeedSearch<'a> {
    index: &'a KmerIndex,
    scheme: &'a ScoringScheme,
    params: SeedParams,
    db_sequences: Vec<&'a [u8]>,
    stats: Option<KarlinAltschul>,
    /// Effective database length for E-value scaling (sum of seq lens).
    db_len: usize,
}

impl<'a> SeedSearch<'a> {
    /// Builds a searcher. `db_sequences` must be the *same* sequences,
    /// in the *same order*, that `index` was built over — the index
    /// stores only positions, so the residues are needed for
    /// extension.
    pub fn new(
        index: &'a KmerIndex,
        db_sequences: Vec<&'a [u8]>,
        scheme: &'a ScoringScheme,
        params: SeedParams,
    ) -> Self {
        let db_len = db_sequences.iter().map(|s| s.len()).sum();
        SeedSearch {
            index,
            scheme,
            params,
            db_sequences,
            stats: None,
            db_len,
        }
    }

    /// Enables Karlin-Altschul E-value / bit-score annotation on every
    /// reported [`Hsp`].
    pub fn with_stats(mut self, ka: KarlinAltschul) -> Self {
        self.stats = Some(ka);
        self
    }

    /// Searches `query` against the database, returning HSPs sorted by
    /// descending score.
    pub fn search(&self, query: &[u8]) -> Vec<Hsp> {
        // 1. Seed: collect (q_off, hit), bin by (seq_id, diagonal).
        //    diagonal = target_offset - query_offset (can be negative).
        use std::collections::HashMap;
        let mut diagonals: HashMap<(usize, isize), Vec<(usize, usize)>> = HashMap::new();
        for (q_off, hit) in self.index.seed_query(query) {
            let diag = hit.offset as isize - q_off as isize;
            diagonals
                .entry((hit.seq_id, diag))
                .or_default()
                .push((q_off, hit.offset));
        }

        // 2. Extend the best seed in each qualifying diagonal bin.
        let mut hsps: Vec<Hsp> = Vec::new();
        for ((seq_id, _diag), seeds) in diagonals {
            if seeds.len() < self.params.min_seeds_per_diagonal {
                continue;
            }
            let target = self.db_sequences[seq_id];
            // Anchor on the first seed of the bin.
            let &(q_seed, t_seed) = &seeds[0];
            let hsp = self.extend_seed(query, target, seq_id, q_seed, t_seed);
            if hsp.score >= self.params.min_score {
                hsps.push(hsp);
            }
        }

        // 3. Deduplicate overlapping HSPs from neighbouring diagonals:
        //    keep the highest-scoring per (seq_id, query_start rounded).
        hsps.sort_by(|a, b| b.score.cmp(&a.score));
        let mut kept: Vec<Hsp> = Vec::new();
        for h in hsps {
            let overlaps = kept.iter().any(|k| {
                k.seq_id == h.seq_id
                    && spans_overlap(k.query_span, h.query_span)
                    && spans_overlap(k.target_span, h.target_span)
            });
            if !overlaps {
                kept.push(h);
            }
        }

        // 4. Annotate with statistics if requested.
        if let Some(ka) = self.stats {
            for h in &mut kept {
                h.bit_score = Some(ka.bit_score(h.score));
                h.e_value = Some(ka.e_value(h.score, query.len().max(1), self.db_len.max(1)));
            }
        }
        kept
    }

    /// Extends one seed left and right with an X-drop gapless walk plus
    /// a small banded gapped refinement. Returns the resulting HSP.
    ///
    /// The seed itself spans `[q_seed, q_seed + k)` of the query and
    /// the equivalent of the target; extension grows that segment.
    fn extend_seed(
        &self,
        query: &[u8],
        target: &[u8],
        seq_id: usize,
        q_seed: usize,
        t_seed: usize,
    ) -> Hsp {
        let k = self.index.k();
        // Seed core score. Use `.get()` so a caller whose `query`/
        // `target` are shorter than `q_seed + k`/`t_seed + k` (a
        // length-contract mismatch) does not index out of bounds — the
        // core simply stops where either sequence ends.
        let mut score: i32 = 0;
        let mut core = 0usize;
        for i in 0..k {
            match (query.get(q_seed + i), target.get(t_seed + i)) {
                (Some(&a), Some(&b)) => {
                    score += self.scheme.sub(a, b);
                    core = i + 1;
                }
                _ => break,
            }
        }

        // --- Extend right (gapless X-drop) -------------------------
        // Start past the (possibly short-circuited) core; `core` never
        // exceeds either remaining length, so these stay in range.
        let mut q = q_seed + core;
        let mut t = t_seed + core;
        let mut run = score;
        let mut best = score;
        let mut best_q = q;
        let mut best_t = t;
        while q < query.len() && t < target.len() {
            run += self.scheme.sub(query[q], target[t]);
            q += 1;
            t += 1;
            if run > best {
                best = run;
                best_q = q;
                best_t = t;
            } else if best - run > self.params.x_drop {
                break;
            }
        }
        let r_end_q = best_q;
        let r_end_t = best_t;
        score = best;

        // --- Extend left (gapless X-drop) --------------------------
        // Clamp the start to each sequence's length so a `q_seed`/
        // `t_seed` past the (short) sequence end can't index out of
        // bounds when the loop steps left.
        let mut q = q_seed.min(query.len());
        let mut t = t_seed.min(target.len());
        let mut run = score;
        let mut best = score;
        let mut best_q = q;
        let mut best_t = t;
        while q > 0 && t > 0 {
            q -= 1;
            t -= 1;
            run += self.scheme.sub(query[q], target[t]);
            if run > best {
                best = run;
                best_q = q;
                best_t = t;
            } else if best - run > self.params.x_drop {
                break;
            }
        }

        Hsp {
            seq_id,
            query_span: (best_q, r_end_q),
            target_span: (best_t, r_end_t),
            score: best,
            bit_score: None,
            e_value: None,
        }
    }
}

/// `true` if two half-open spans overlap by at least one position.
fn spans_overlap(a: (usize, usize), b: (usize, usize)) -> bool {
    a.0 < b.1 && b.0 < a.1
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::{GapCost, ScoringScheme, SubstitutionMatrix};

    fn dna_scheme() -> ScoringScheme {
        ScoringScheme::new(SubstitutionMatrix::dna_simple(2, -3), GapCost::new(5, 2))
    }

    #[test]
    fn finds_obvious_homology() {
        // Database holds a sequence containing GATTACACAT...; the query
        // shares a long exact core.
        let db: Vec<&[u8]> = vec![b"TTTTTGATTACACATGGGGG"];
        let idx = KmerIndex::build_many(&db, 6).unwrap();
        let scheme = dna_scheme();
        let search = SeedSearch::new(&idx, db.clone(), &scheme, SeedParams::default());
        let hsps = search.search(b"GATTACACAT");
        assert!(!hsps.is_empty(), "should find the embedded core");
        let top = &hsps[0];
        assert_eq!(top.seq_id, 0);
        assert!(
            top.score >= 16,
            "10 matches * 2 expected, got {}",
            top.score
        );
    }

    #[test]
    fn no_homology_returns_empty() {
        let db: Vec<&[u8]> = vec![b"AAAAAAAAAAAAAAAA"];
        let idx = KmerIndex::build_many(&db, 6).unwrap();
        let scheme = dna_scheme();
        let search = SeedSearch::new(&idx, db.clone(), &scheme, SeedParams::default());
        let hsps = search.search(b"CCCCCCCCCCCC");
        assert!(hsps.is_empty());
    }

    #[test]
    fn extension_grows_beyond_seed() {
        // The seed k-mer is short; extension should capture the full
        // matching region either side of it.
        let db: Vec<&[u8]> = vec![b"ACGTACGTACGTACGT"];
        let idx = KmerIndex::build_many(&db, 6).unwrap();
        let scheme = dna_scheme();
        let search = SeedSearch::new(&idx, db.clone(), &scheme, SeedParams::default());
        let hsps = search.search(b"ACGTACGTACGTACGT");
        let top = &hsps[0];
        // Full 16-residue identity => span covers the whole query.
        assert_eq!(top.query_span, (0, 16));
        assert_eq!(top.score, 32);
    }

    #[test]
    fn multi_database_picks_right_sequence() {
        let db: Vec<&[u8]> = vec![b"AAAAAAAAAAAA", b"GGGGGATTACAGGGGG", b"CCCCCCCCCCCC"];
        let idx = KmerIndex::build_many(&db, 5).unwrap();
        let scheme = dna_scheme();
        let search = SeedSearch::new(&idx, db.clone(), &scheme, SeedParams::default());
        let hsps = search.search(b"GATTACA");
        assert!(!hsps.is_empty());
        assert_eq!(hsps[0].seq_id, 1);
    }

    #[test]
    fn extend_seed_does_not_panic_past_sequence_end() {
        // A caller whose target is shorter than `t_seed + k` (or query
        // shorter than `q_seed + k`) must not trigger an out-of-bounds
        // index in the seed-core / extension loops — extension simply
        // ends when either sequence runs out.
        let db: Vec<&[u8]> = vec![b"ACGTACGT"];
        let idx = KmerIndex::build_many(&db, 4).unwrap();
        let scheme = dna_scheme();
        let search = SeedSearch::new(&idx, db.clone(), &scheme, SeedParams::default());
        // idx.k() == 4
        let query = b"ACGTACGT"; // len 8
        let target = b"ACG"; // len 3 < k -> core loop would OOB on target

        // Seed positions that run past the (short) target and near the
        // query end. With the old `query[q_seed + i]`/`target[t_seed+i]`
        // this panics; bounded extension must return a valid HSP.
        let hsp = search.extend_seed(query, target, 0, 0, 0);
        // Spans stay within both sequences and are well-formed.
        assert!(hsp.query_span.0 <= hsp.query_span.1);
        assert!(hsp.target_span.0 <= hsp.target_span.1);
        assert!(hsp.query_span.1 <= query.len());
        assert!(hsp.target_span.1 <= target.len());

        // Also exercise a seed anchored at the very end of the target.
        let hsp2 = search.extend_seed(query, target, 0, query.len() - 1, target.len() - 1);
        assert!(hsp2.target_span.1 <= target.len());
        assert!(hsp2.query_span.1 <= query.len());
    }

    #[test]
    fn stats_annotation() {
        let db: Vec<&[u8]> = vec![b"TTTTTGATTACACATGGGGG"];
        let idx = KmerIndex::build_many(&db, 6).unwrap();
        let scheme = dna_scheme();
        let search = SeedSearch::new(&idx, db.clone(), &scheme, SeedParams::default())
            .with_stats(KarlinAltschul::dna_ungapped());
        let hsps = search.search(b"GATTACACAT");
        assert!(hsps[0].bit_score.is_some());
        assert!(hsps[0].e_value.is_some());
        assert!(hsps[0].e_value.unwrap() >= 0.0);
    }
}
