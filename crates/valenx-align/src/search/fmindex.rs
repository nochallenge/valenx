//! FM-index — Burrows-Wheeler transform + backward search.
//!
//! The FM-index is the core data structure of BWA and Bowtie: it
//! supports *exact substring search* in time proportional to the
//! pattern length, independent of the reference length, in space close
//! to the reference itself.
//!
//! This module ships the production layout:
//!
//! - the **suffix array** is built by **SA-IS** — the linear-time
//!   `O(n)` induced-sorting algorithm of Nong, Zhang & Chan (2009) —
//!   rather than the textbook `O(n log²n)` prefix-doubling sort. SA-IS
//!   is the algorithm BWA's `bwtsw2` and minimap2's index builder use.
//! - the **BWT** is `BWT[i] = T[SA[i] − 1]`.
//! - the **C array** (`C[c]` = number of reference characters
//!   alphabetically smaller than `c`).
//! - the **rank structure** uses a sampled / block-rank layout: for
//!   each character we store *only every `BLOCK_SIZE`-th* cumulative
//!   count (a sample), and fill in the remaining ranks by scanning
//!   inside the block. Memory is `O((σ·n)/BLOCK_SIZE + n)` rather than
//!   the `O(σ·n)` of a dense Occ table; rank queries are O(1)
//!   amortised (one direct sample read + ≤ `BLOCK_SIZE` byte
//!   comparisons in the BWT block). This is the same trade-off BWA
//!   makes (its block size is `1<<7`; we use `1<<6` = 64).
//! - the **suffix-array samples** keep one in every `SA_SAMPLE_RATE`
//!   suffix positions; `locate()` walks the LF-mapping back to the
//!   nearest sample to recover an unstored position, the BWA / Bowtie
//!   trick to cut SA memory from `n·log n` bits to a fraction.
//!
//! [`FmIndex::count`] returns how many times a pattern occurs;
//! [`FmIndex::locate`] returns every 0-based start position. Both run
//! the standard LF-mapping backward-search loop.

use crate::error::{AlignError, Result};

/// The sentinel terminator appended to the text. It must be
/// lexicographically smaller than every real character; `0x00` works
/// because biological sequences are printable ASCII.
const SENTINEL: u8 = 0;

/// Block size for the rank structure: every `BLOCK_SIZE`-th cumulative
/// count is sampled; a rank query reads the sample plus scans within
/// the block. 64 keeps the worst-case scan at one cache line.
const BLOCK_SIZE: usize = 64;

/// Default suffix-array sample rate. Every `SA_SAMPLE_RATE`-th SA
/// entry is kept; the rest are recovered by LF-mapping. A value of 32
/// trades 1/32 the SA memory for at most `SA_SAMPLE_RATE` rank steps
/// per `locate()` per occurrence.
const SA_SAMPLE_RATE: usize = 32;

/// A Burrows-Wheeler FM-index over a byte text.
#[derive(Clone, Debug)]
pub struct FmIndex {
    /// The original text length (without the sentinel).
    text_len: usize,
    /// The BWT string (same length as the sentinel-terminated text).
    bwt: Vec<u8>,
    /// Sampled suffix array: `sa_sampled[k]` = `SA[k * SA_SAMPLE_RATE]`
    /// when sampled. A parallel `sa_sample_present` bitmap flags
    /// which `i` carry a stored value.
    sa_sampled: Vec<usize>,
    /// `sa_sample_present[i]` = true iff `SA[i]` is stored directly in
    /// `sa_sampled[sa_sample_rank[i]]`.
    sa_sample_present: Vec<bool>,
    /// Prefix-sum over `sa_sample_present` so we can index into
    /// `sa_sampled` in O(1).
    sa_sample_rank: Vec<u32>,
    /// SA sample rate actually used.
    sa_sample_rate: usize,
    /// Sorted list of distinct characters present (incl. sentinel).
    alphabet: Vec<u8>,
    /// `c_array[idx]` = count of characters smaller than
    /// `alphabet[idx]`. Indexed by alphabet position.
    c_array: Vec<usize>,
    /// Block-sampled occurrence table: `occ_blocks[idx]` is a length
    /// `(n / BLOCK_SIZE + 1)` array giving the count of
    /// `alphabet[idx]` in `bwt[0..i * BLOCK_SIZE]`.
    occ_blocks: Vec<Vec<u32>>,
    /// Dense map ASCII byte -> alphabet index, or `usize::MAX`.
    rank_of: [usize; 256],
}

impl FmIndex {
    /// Builds an FM-index over `text`.
    ///
    /// Returns [`AlignError::Invalid`] if `text` is empty or contains a
    /// `0x00` byte (reserved as the BWT sentinel).
    pub fn build(text: &[u8]) -> Result<Self> {
        Self::build_with(text, SA_SAMPLE_RATE)
    }

    /// Builds an FM-index with a custom SA sample rate. `sample_rate`
    /// must be `>= 1`; `1` stores the full SA (no sampling).
    pub fn build_with(text: &[u8], sample_rate: usize) -> Result<Self> {
        if text.is_empty() {
            return Err(AlignError::invalid("text", "cannot index an empty text"));
        }
        if text.contains(&SENTINEL) {
            return Err(AlignError::invalid(
                "text",
                "text must not contain a NUL (0x00) byte (reserved sentinel)",
            ));
        }
        if sample_rate == 0 {
            return Err(AlignError::invalid("sample_rate", "must be >= 1"));
        }

        // Sentinel-terminated text.
        let mut t = text.to_vec();
        t.push(SENTINEL);
        let n = t.len();

        // SA-IS suffix array (linear-time).
        let sa = sa_is(&t);

        // BWT[i] = T[SA[i] - 1], wrapping the sentinel at SA[i] == 0.
        let bwt: Vec<u8> = sa
            .iter()
            .map(|&p| if p == 0 { t[n - 1] } else { t[p - 1] })
            .collect();

        // Alphabet = sorted distinct bytes of t.
        let mut alphabet: Vec<u8> = t.clone();
        alphabet.sort_unstable();
        alphabet.dedup();

        let mut rank_of = [usize::MAX; 256];
        for (idx, &c) in alphabet.iter().enumerate() {
            rank_of[c as usize] = idx;
        }

        // C array: number of characters strictly smaller than each.
        let mut counts = vec![0usize; alphabet.len()];
        for &c in &t {
            counts[rank_of[c as usize]] += 1;
        }
        let mut c_array = vec![0usize; alphabet.len()];
        let mut running = 0;
        for (idx, &cnt) in counts.iter().enumerate() {
            c_array[idx] = running;
            running += cnt;
        }

        // Block-sampled Occ: `occ_blocks[idx][b]` = count of
        // `alphabet[idx]` in `bwt[0..b * BLOCK_SIZE]`. The largest `b`
        // a rank query may need is `n / BLOCK_SIZE` (for `i = n`), so
        // we allocate `n / BLOCK_SIZE + 1` blocks.
        let num_blocks = n / BLOCK_SIZE + 1;
        let mut occ_blocks: Vec<Vec<u32>> = vec![vec![0u32; num_blocks]; alphabet.len()];
        let mut running_per_char = vec![0u32; alphabet.len()];
        for i in 0..n {
            // bwt[i] gets consumed into running_per_char; before doing
            // so, if `i+1` lands on a block boundary then the block
            // *after* consuming this character is the cumulative
            // count up to (and including) i.
            let cr = rank_of[bwt[i] as usize];
            running_per_char[cr] += 1;
            if (i + 1) % BLOCK_SIZE == 0 {
                let blk = (i + 1) / BLOCK_SIZE;
                for (idx, row) in occ_blocks.iter_mut().enumerate() {
                    row[blk] = running_per_char[idx];
                }
            }
        }
        // Block 0 stays at 0 (count over empty prefix).

        // Sample the suffix array: keep entries whose *text* position
        // is divisible by `sample_rate`. We also always keep SA[0]
        // (the sentinel position) so the LF-walk has a known anchor.
        let mut sa_sample_present = vec![false; n];
        let mut sa_sampled = Vec::new();
        let mut sa_sample_rank = vec![0u32; n + 1];
        for i in 0..n {
            sa_sample_rank[i] = sa_sampled.len() as u32;
            let keep = sample_rate == 1 || sa[i] % sample_rate == 0 || i == 0;
            if keep {
                sa_sample_present[i] = true;
                sa_sampled.push(sa[i]);
            }
        }
        sa_sample_rank[n] = sa_sampled.len() as u32;

        Ok(FmIndex {
            text_len: text.len(),
            bwt,
            sa_sampled,
            sa_sample_present,
            sa_sample_rank,
            sa_sample_rate: sample_rate,
            alphabet,
            c_array,
            occ_blocks,
            rank_of,
        })
    }

    /// Length of the indexed text (excluding the sentinel).
    pub fn text_len(&self) -> usize {
        self.text_len
    }

    /// The Burrows-Wheeler transform string (with the sentinel).
    pub fn bwt(&self) -> &[u8] {
        &self.bwt
    }

    /// The full suffix array, reconstructed if it was sampled. O(n) in
    /// the index size — provided for tests / debugging only.
    pub fn suffix_array(&self) -> Vec<usize> {
        let n = self.bwt.len();
        (0..n).map(|i| self.locate_sa(i)).collect()
    }

    /// The sorted distinct characters present in the indexed text
    /// (including the `0x00` sentinel as the first entry).
    pub fn alphabet(&self) -> &[u8] {
        &self.alphabet
    }

    /// The SA sample rate this index was built with.
    pub fn sa_sample_rate(&self) -> usize {
        self.sa_sample_rate
    }

    /// Inverse BWT — reconstructs the original text from this index.
    /// Returns the text **without** the sentinel byte.
    pub fn inverse_bwt(&self) -> Vec<u8> {
        let n = self.bwt.len();
        let mut out = vec![0u8; n - 1];
        // Walk LF-mapping from position 0 (the sentinel row) backwards.
        // T[n - 1 - step] = self.bwt[i]; i = LF(i).
        let mut i = 0usize;
        for step in 0..(n - 1) {
            let c = self.bwt[i];
            // SA position of the next row in inverse order.
            let cr = self.rank_of[c as usize];
            let rank_before = self.rank(cr, i);
            i = self.c_array[cr] + rank_before;
            // T[n - 2 - step] = c (we skip the sentinel itself).
            let pos = n - 2 - step;
            out[pos] = c;
        }
        out
    }

    /// Number of occurrences of `pattern` in the text.
    ///
    /// An empty pattern matches `text_len + 1` times (every position
    /// including the end); an unknown character yields `0`.
    pub fn count(&self, pattern: &[u8]) -> usize {
        match self.backward_search(pattern) {
            Some((lo, hi)) => hi - lo,
            None => 0,
        }
    }

    /// `true` if `pattern` occurs at least once.
    pub fn contains(&self, pattern: &[u8]) -> bool {
        self.count(pattern) > 0
    }

    /// Every 0-based start position of `pattern`, sorted ascending.
    pub fn locate(&self, pattern: &[u8]) -> Vec<usize> {
        let mut hits = match self.backward_search(pattern) {
            Some((lo, hi)) => (lo..hi).map(|i| self.locate_sa(i)).collect::<Vec<_>>(),
            None => Vec::new(),
        };
        hits.sort_unstable();
        hits
    }

    /// Recovers `SA[i]` — either reading the sample directly or
    /// walking LF-mapping until we hit a sampled position.
    fn locate_sa(&self, mut i: usize) -> usize {
        let mut steps = 0usize;
        while !self.sa_sample_present[i] {
            // LF-step: i <- C[BWT[i]] + Occ(BWT[i], i).
            let c = self.bwt[i];
            let cr = self.rank_of[c as usize];
            let rank_before = self.rank(cr, i);
            i = self.c_array[cr] + rank_before;
            steps += 1;
            // Safety net: the walk must terminate inside `sample_rate`
            // steps because every `sample_rate`-th text position is
            // stored.
            if steps > self.sa_sample_rate + self.bwt.len() {
                // unreachable on a well-formed index
                break;
            }
        }
        let sample_idx = self.sa_sample_rank[i] as usize;
        let stored = self.sa_sampled[sample_idx];
        // `stored` is the text position of the *current* row; we have
        // walked `steps` LF-mapping moves backwards, so the row we
        // *started* at corresponds to text position `stored + steps`
        // (modulo n).
        let n = self.bwt.len();
        (stored + steps) % n
    }

    /// `rank(c, i)` = number of `alphabet[cr]` in `bwt[0..i]`.
    fn rank(&self, cr: usize, i: usize) -> usize {
        let block = i / BLOCK_SIZE;
        let block_start = block * BLOCK_SIZE;
        let base = self.occ_blocks[cr][block] as usize;
        let c = self.alphabet[cr];
        let mut extra = 0usize;
        // Scan within the block.
        for j in block_start..i {
            if self.bwt[j] == c {
                extra += 1;
            }
        }
        base + extra
    }

    /// The LF-mapping backward search. Returns the half-open suffix-
    /// array interval `[lo, hi)` matching `pattern`, or `None` when the
    /// pattern does not occur.
    fn backward_search(&self, pattern: &[u8]) -> Option<(usize, usize)> {
        let n = self.bwt.len();
        let mut lo = 0usize;
        let mut hi = n;
        // Process the pattern right-to-left.
        for &ch in pattern.iter().rev() {
            let cr = *self.rank_of.get(ch as usize)?;
            if cr == usize::MAX {
                return None; // character absent from the reference
            }
            let c = self.c_array[cr];
            lo = c + self.rank(cr, lo);
            hi = c + self.rank(cr, hi);
            if lo >= hi {
                return None;
            }
        }
        Some((lo, hi))
    }

    /// Find super-maximal exact matches (SMEMs) of `query` against
    /// the indexed text, of length `>= min_len`. An SMEM ending at
    /// query position `i` is the *longest* substring of the query
    /// ending at `i` that matches the reference at least once, and
    /// that is **not strictly contained** in any longer match.
    ///
    /// This is the BWA-MEM seeding primitive: a single scan over the
    /// query that, for each starting position `i`, grows the match
    /// rightward via backward search on the reverse of the growing
    /// substring until the SA interval becomes empty, then records
    /// the longest substring whose interval was still alive. Total
    /// work is `O(|query| · L̄ · σ_lookup)` where L̄ is the average
    /// SMEM length — fine for typical short reads.
    ///
    /// Returned matches are filtered to be maximal: SMEMs strictly
    /// contained in another are dropped.
    pub fn smems(&self, query: &[u8], min_len: usize) -> Vec<Smem> {
        let qn = query.len();
        let mut raw: Vec<Smem> = Vec::new();
        let n = self.bwt.len();

        // For each start `i`, extend as far right as possible. We
        // recompute the backward search from scratch each step; this
        // costs `O(L̄^2)` total but is simple and correct.
        let mut i = 0usize;
        while i < qn {
            let mut best_end: usize = i;
            let mut best_lo = 0usize;
            let mut best_hi = n;
            // Grow rightward j = i+1, i+2, ...; for each j compute the
            // SA interval of query[i..j] via backward search.
            let mut j = i + 1;
            while j <= qn {
                // Backward search on query[i..j].
                let mut lo = 0usize;
                let mut hi = n;
                let mut alive = true;
                for &ch in query[i..j].iter().rev() {
                    let cr = self.rank_of[ch as usize];
                    if cr == usize::MAX {
                        alive = false;
                        break;
                    }
                    let c = self.c_array[cr];
                    lo = c + self.rank(cr, lo);
                    hi = c + self.rank(cr, hi);
                    if lo >= hi {
                        alive = false;
                        break;
                    }
                }
                if !alive {
                    break;
                }
                best_end = j;
                best_lo = lo;
                best_hi = hi;
                j += 1;
            }
            if best_end > i && best_end - i >= min_len {
                raw.push(Smem {
                    query_start: i,
                    query_end: best_end,
                    count: best_hi - best_lo,
                    sa_lo: best_lo,
                    sa_hi: best_hi,
                });
                // Advance to the first position past the start of the
                // recorded SMEM. Stepping by 1 keeps every position's
                // longest-match candidate so the maximality filter
                // below can pick the genuine SMEMs.
                i += 1;
            } else {
                i += 1;
            }
        }

        // Maximality filter: drop SMEMs strictly contained in another
        // (same-or-larger span on at least one side, with at least one
        // strict). Equal spans are kept once.
        raw.sort_by(|a, b| {
            a.query_start
                .cmp(&b.query_start)
                .then(b.query_end.cmp(&a.query_end))
        });
        let mut maximal: Vec<Smem> = Vec::new();
        for s in raw.iter() {
            let strictly_contained = maximal.iter().any(|m| {
                m.query_start <= s.query_start
                    && m.query_end >= s.query_end
                    && (m.query_start < s.query_start || m.query_end > s.query_end)
            });
            let equal_already = maximal
                .iter()
                .any(|m| m.query_start == s.query_start && m.query_end == s.query_end);
            if !strictly_contained && !equal_already {
                maximal.push(*s);
            }
        }
        maximal
    }

    /// Resolve every text position of an SMEM into a sorted vector of
    /// 0-based reference offsets — `locate()` over the SA interval.
    pub fn smem_positions(&self, smem: &Smem) -> Vec<usize> {
        let mut hits: Vec<usize> = (smem.sa_lo..smem.sa_hi)
            .map(|i| self.locate_sa(i))
            .collect();
        hits.sort_unstable();
        hits
    }
}

/// A super-maximal exact match against the FM-index reference.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Smem {
    /// 0-based start of the match in the query.
    pub query_start: usize,
    /// Half-open end of the match in the query.
    pub query_end: usize,
    /// Number of occurrences in the reference.
    pub count: usize,
    /// Suffix-array interval lower bound.
    pub sa_lo: usize,
    /// Suffix-array interval upper bound (exclusive).
    pub sa_hi: usize,
}

impl Smem {
    /// Length of the match (always positive when the SMEM is real).
    pub fn len(&self) -> usize {
        self.query_end - self.query_start
    }

    /// Always `false` for a real SMEM (zero-length matches never
    /// appear). Provided for clippy hygiene.
    pub fn is_empty(&self) -> bool {
        self.query_end <= self.query_start
    }
}

// =====================================================================
// SA-IS — linear-time suffix-array construction
// =====================================================================
//
// Nong, Zhang & Chan, "Two efficient algorithms for linear time suffix
// array construction" (2009). The algorithm classifies each suffix as
// `L` (larger than its right neighbour) or `S` (smaller); an `S` suffix
// whose left neighbour is `L` is an `LMS` ("leftmost S") suffix. The
// trick: sort the LMS suffixes recursively (their reduced text is at
// most half the original length, giving the `T(n) = T(n/2) + O(n)`
// recurrence), then *induce* the sort of all other suffixes from them
// in two linear scans of the text.

/// Public entry point for SA-IS. Returns the suffix array over
/// `t` (which must end with a unique sentinel that is the strict
/// minimum, e.g. `0u8` over a non-zero text).
fn sa_is(t: &[u8]) -> Vec<usize> {
    let n = t.len();
    // Promote bytes to `i32` and run the recursive algorithm. We use
    // signed integers so that recursive levels can encode "undefined"
    // as `-1`.
    let t_i32: Vec<i32> = t.iter().map(|&b| b as i32).collect();
    // Sigma = 256 at the top level (bytes).
    let mut sa = vec![0i32; n];
    sais_core(&t_i32, &mut sa, 256);
    sa.into_iter().map(|x| x as usize).collect()
}

/// Recursive SA-IS over a generic alphabet.
fn sais_core(t: &[i32], sa: &mut [i32], sigma: usize) {
    let n = t.len();
    debug_assert_eq!(sa.len(), n);
    if n == 0 {
        return;
    }
    if n == 1 {
        sa[0] = 0;
        return;
    }
    if n == 2 {
        if t[0] < t[1] {
            sa[0] = 0;
            sa[1] = 1;
        } else {
            sa[0] = 1;
            sa[1] = 0;
        }
        return;
    }

    // 1. Type classification. `is_s[i]` = true iff suffix i is S-type.
    let mut is_s = vec![false; n];
    is_s[n - 1] = true; // sentinel by convention
    for i in (0..n - 1).rev() {
        if t[i] < t[i + 1] {
            is_s[i] = true;
        } else if t[i] == t[i + 1] {
            is_s[i] = is_s[i + 1];
        } // else L, leave false
    }

    // LMS positions: S that is preceded by L.
    let is_lms = |i: usize| i > 0 && is_s[i] && !is_s[i - 1];

    // 2. Bucket sizes.
    let mut bkt = vec![0usize; sigma];
    for &c in t {
        bkt[c as usize] += 1;
    }
    let bucket_ends = |bkt: &Vec<usize>| -> Vec<usize> {
        let mut e = vec![0usize; sigma];
        let mut sum = 0usize;
        for c in 0..sigma {
            sum += bkt[c];
            e[c] = sum;
        }
        e
    };
    let bucket_starts = |bkt: &Vec<usize>| -> Vec<usize> {
        let mut s = vec![0usize; sigma];
        let mut sum = 0usize;
        for c in 0..sigma {
            s[c] = sum;
            sum += bkt[c];
        }
        s
    };

    // 3. First induced sort: place LMS suffixes at the END of their
    // buckets in arbitrary (text) order, then induce L then S.
    for x in sa.iter_mut() {
        *x = -1;
    }
    let mut ends = bucket_ends(&bkt);
    for (i, &ch) in t.iter().enumerate().take(n) {
        if is_lms(i) {
            let c = ch as usize;
            ends[c] -= 1;
            sa[ends[c]] = i as i32;
        }
    }
    induce_l(t, sa, &is_s, &bkt, sigma);
    induce_s(t, sa, &is_s, &bkt, sigma);

    // 4. Name LMS substrings.
    let mut lms_names = vec![-1i32; n];
    let mut name = 0i32;
    let mut prev_lms: Option<usize> = None;
    for &sa_v in sa.iter().take(n) {
        let p = sa_v as usize;
        if !is_lms(p) {
            continue;
        }
        let mut diff = true;
        if let Some(q) = prev_lms {
            diff = !lms_equal(t, &is_s, p, q);
        }
        if diff && prev_lms.is_some() {
            name += 1;
        }
        lms_names[p] = name;
        prev_lms = Some(p);
    }

    // 5. Compact LMS names into a reduced string.
    let mut reduced: Vec<i32> = Vec::new();
    let mut reduced_pos: Vec<usize> = Vec::new();
    for (i, &nm) in lms_names.iter().enumerate().take(n) {
        if is_lms(i) {
            reduced_pos.push(i);
            reduced.push(nm);
        }
    }
    let reduced_sigma = (name + 1) as usize;

    // 6. Solve the reduced problem.
    let mut reduced_sa = vec![0i32; reduced.len()];
    if reduced_sigma == reduced.len() {
        // All names unique => sort directly.
        for (i, &v) in reduced.iter().enumerate() {
            reduced_sa[v as usize] = i as i32;
        }
    } else {
        sais_core(&reduced, &mut reduced_sa, reduced_sigma);
    }

    // 7. Final induced sort using the LMS order from the reduction.
    for x in sa.iter_mut() {
        *x = -1;
    }
    let mut ends = bucket_ends(&bkt);
    // Process LMS suffixes in reduced-SA order so they go to the END
    // of their buckets in *correct* sorted order.
    for k in (0..reduced_sa.len()).rev() {
        let idx = reduced_pos[reduced_sa[k] as usize];
        let c = t[idx] as usize;
        ends[c] -= 1;
        sa[ends[c]] = idx as i32;
    }
    induce_l(t, sa, &is_s, &bkt, sigma);
    induce_s(t, sa, &is_s, &bkt, sigma);

    let _ = bucket_starts; // (kept for symmetry; unused here)
}

/// Induced sort of L suffixes (left-to-right scan).
fn induce_l(t: &[i32], sa: &mut [i32], is_s: &[bool], bkt: &[usize], sigma: usize) {
    let mut starts = vec![0usize; sigma];
    let mut sum = 0usize;
    for c in 0..sigma {
        starts[c] = sum;
        sum += bkt[c];
    }
    let n = t.len();
    for i in 0..n {
        let v = sa[i];
        if v <= 0 {
            continue;
        }
        let j = v as usize - 1;
        if !is_s[j] {
            let c = t[j] as usize;
            sa[starts[c]] = j as i32;
            starts[c] += 1;
        }
    }
}

/// Induced sort of S suffixes (right-to-left scan).
fn induce_s(t: &[i32], sa: &mut [i32], is_s: &[bool], bkt: &[usize], sigma: usize) {
    let mut ends = vec![0usize; sigma];
    let mut sum = 0usize;
    for c in 0..sigma {
        sum += bkt[c];
        ends[c] = sum;
    }
    let n = t.len();
    for i in (0..n).rev() {
        let v = sa[i];
        if v <= 0 {
            continue;
        }
        let j = v as usize - 1;
        if is_s[j] {
            let c = t[j] as usize;
            ends[c] -= 1;
            sa[ends[c]] = j as i32;
        }
    }
}

/// Compare two LMS substrings starting at `p` and `q` for equality.
/// Two LMS substrings are equal iff they agree character-by-character
/// and have the same type sequence up to and including the next LMS.
fn lms_equal(t: &[i32], is_s: &[bool], p: usize, q: usize) -> bool {
    let n = t.len();
    if p == q {
        return true;
    }
    let mut i = 0usize;
    loop {
        let pi = p + i;
        let qi = q + i;
        let p_at_end = pi >= n;
        let q_at_end = qi >= n;
        if p_at_end || q_at_end {
            return p_at_end == q_at_end;
        }
        let p_lms = i > 0 && is_s[pi] && !is_s[pi - 1];
        let q_lms = i > 0 && is_s[qi] && !is_s[qi - 1];
        if i > 0 && (p_lms || q_lms) {
            return p_lms == q_lms && t[pi] == t[qi];
        }
        if t[pi] != t[qi] || is_s[pi] != is_s[qi] {
            return false;
        }
        i += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A brute-force suffix array over an arbitrary byte text.
    fn brute_sa(t: &[u8]) -> Vec<usize> {
        let mut sa: Vec<usize> = (0..t.len()).collect();
        sa.sort_by(|&a, &b| t[a..].cmp(&t[b..]));
        sa
    }

    #[test]
    fn rejects_empty_and_nul() {
        assert!(FmIndex::build(b"").is_err());
        assert!(FmIndex::build(&[b'A', 0u8, b'C']).is_err());
    }

    #[test]
    fn sais_matches_brute_force_on_small_strings() {
        for &text in &[
            b"banana".as_slice(),
            b"mississippi",
            b"abracadabra",
            b"GATTACAGATTACA",
            b"AAAAAA",
            b"abcdefghij",
            b"jihgfedcba",
            b"ababababab",
            b"the quick brown fox jumps",
        ] {
            let mut t = text.to_vec();
            t.push(SENTINEL);
            let mine = sa_is(&t);
            let brute = brute_sa(&t);
            assert_eq!(mine, brute, "SA-IS mismatch on {:?}", std::str::from_utf8(text));
        }
    }

    #[test]
    fn sais_handles_repeated_lms_substrings() {
        // A pathological string for SA-IS: many equal LMS substrings.
        let mut t = b"ABABABABABABABAB".to_vec();
        t.push(SENTINEL);
        let mine = sa_is(&t);
        let brute = brute_sa(&t);
        assert_eq!(mine, brute);
    }

    #[test]
    fn sais_random_strings_match_brute_force() {
        // Pseudo-random strings: cycle a small alphabet so we exercise
        // many tie-breaking paths. Deterministic — no rand crate dep.
        let alphabet = b"ACGT";
        for seed in 0..32u32 {
            let mut h: u32 = seed.wrapping_mul(0x9E3779B1).wrapping_add(1);
            let len = 5 + (seed as usize % 60);
            let mut text = Vec::with_capacity(len);
            for _ in 0..len {
                h ^= h << 13;
                h ^= h >> 17;
                h ^= h << 5;
                text.push(alphabet[(h as usize) % alphabet.len()]);
            }
            let mut t = text.clone();
            t.push(SENTINEL);
            let mine = sa_is(&t);
            let brute = brute_sa(&t);
            assert_eq!(mine, brute, "SA-IS mismatch on seed {seed}");
        }
    }

    #[test]
    fn suffix_array_is_a_permutation() {
        let fm = FmIndex::build(b"BANANA").unwrap();
        let mut sorted = fm.suffix_array().to_vec();
        sorted.sort_unstable();
        assert_eq!(sorted, (0..7).collect::<Vec<_>>()); // 6 + sentinel
    }

    #[test]
    fn banana_bwt_known() {
        // The textbook "banana$" BWT is "annb$aa".
        let fm = FmIndex::build(b"banana").unwrap();
        assert_eq!(fm.bwt(), b"annb\0aa");
    }

    #[test]
    fn exact_count_and_locate() {
        let fm = FmIndex::build(b"GATTACAGATTACA").unwrap();
        // "GATTACA" occurs twice, at 0 and 7.
        assert_eq!(fm.count(b"GATTACA"), 2);
        assert_eq!(fm.locate(b"GATTACA"), vec![0, 7]);
        // "ATTACA" occurs twice as well.
        assert_eq!(fm.count(b"ATTACA"), 2);
        // "TT" occurs twice.
        assert_eq!(fm.count(b"TT"), 2);
    }

    #[test]
    fn absent_pattern() {
        let fm = FmIndex::build(b"ACGTACGT").unwrap();
        assert_eq!(fm.count(b"TTTT"), 0);
        assert!(!fm.contains(b"GGGG"));
        // A character entirely absent from the reference.
        assert_eq!(fm.count(b"N"), 0);
        assert!(fm.locate(b"NN").is_empty());
    }

    #[test]
    fn single_character_text() {
        let fm = FmIndex::build(b"AAAA").unwrap();
        assert_eq!(fm.count(b"A"), 4);
        assert_eq!(fm.count(b"AA"), 3);
        assert_eq!(fm.count(b"AAAA"), 1);
        assert_eq!(fm.locate(b"AA"), vec![0, 1, 2]);
    }

    #[test]
    fn alphabet_lists_distinct_characters() {
        let fm = FmIndex::build(b"ACGTACGT").unwrap();
        // Sentinel 0x00 plus the four bases, sorted.
        assert_eq!(fm.alphabet(), &[0u8, b'A', b'C', b'G', b'T']);
    }

    #[test]
    fn locate_matches_naive_scan() {
        let text = b"ACGTACGATACGATACGTACG";
        let fm = FmIndex::build(text).unwrap();
        for pat in [b"ACG".as_slice(), b"GAT", b"TACG", b"ACGT"] {
            let mut naive = Vec::new();
            for start in 0..=text.len().saturating_sub(pat.len()) {
                if &text[start..start + pat.len()] == pat {
                    naive.push(start);
                }
            }
            assert_eq!(fm.locate(pat), naive, "mismatch for {:?}", std::str::from_utf8(pat));
        }
    }

    #[test]
    fn locate_with_sampled_sa() {
        // Aggressive sampling: 1-in-8. The LF-walk must recover every
        // unsampled position correctly.
        let text = b"GATTACAGATTACAGATTACA";
        let fm = FmIndex::build_with(text, 8).unwrap();
        assert_eq!(fm.locate(b"GATTACA"), vec![0, 7, 14]);
        assert_eq!(fm.count(b"ATT"), 3);
    }

    #[test]
    fn inverse_bwt_recovers_text() {
        for &text in &[
            b"banana".as_slice(),
            b"mississippi",
            b"GATTACAGATTACA",
            b"ACGTACGTACGT",
            b"AAAAAA",
        ] {
            let fm = FmIndex::build(text).unwrap();
            let recovered = fm.inverse_bwt();
            assert_eq!(recovered, text, "inverse BWT failed on {:?}", std::str::from_utf8(text));
        }
    }

    #[test]
    fn rank_matches_naive_scan() {
        // The rank function is the workhorse of backward search; spot-
        // check it directly against a naive byte count.
        let fm = FmIndex::build(b"GATTACAGATTACA").unwrap();
        for (idx, &c) in fm.alphabet.iter().enumerate() {
            for i in 0..=fm.bwt.len() {
                let expected = fm.bwt[..i].iter().filter(|&&b| b == c).count();
                let got = fm.rank(idx, i);
                assert_eq!(got, expected, "rank mismatch at c={c}, i={i}");
            }
        }
    }

    #[test]
    fn smems_recover_unique_match() {
        // Reference: TTTTACGTACGGGGG; SMEM ACGTACG of length 7 unique.
        let text = b"TTTTACGTACGGGGGG";
        let fm = FmIndex::build(text).unwrap();
        let smems = fm.smems(b"ACGTACG", 4);
        assert!(!smems.is_empty(), "SMEM should exist for ACGTACG");
        let best = smems.iter().max_by_key(|s| s.len()).unwrap();
        assert_eq!(best.query_start, 0);
        assert_eq!(best.query_end, 7);
        assert_eq!(best.count, 1);
        let positions = fm.smem_positions(best);
        assert_eq!(positions, vec![4]);
    }

    #[test]
    fn smems_handle_repeats() {
        // ACGT repeats — long substring matches many places.
        let text = b"ACGTACGTACGTACGT";
        let fm = FmIndex::build(text).unwrap();
        let smems = fm.smems(b"ACGTACGT", 4);
        // Longest SMEM is the whole query, exact match at offsets
        // 0, 4, 8.
        let best = smems.iter().max_by_key(|s| s.len()).unwrap();
        assert_eq!(best.len(), 8);
        let pos = fm.smem_positions(best);
        assert_eq!(pos, vec![0, 4, 8]);
    }
}
