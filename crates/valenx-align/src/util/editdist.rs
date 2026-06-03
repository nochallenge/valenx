//! Edit distance — Levenshtein DP and Myers' bit-parallel algorithm.
//!
//! The **Levenshtein distance** of two strings is the minimum number
//! of single-character insertions, deletions and substitutions that
//! turn one into the other. Two routines:
//!
//! - [`levenshtein`] — the textbook O(nm) two-row dynamic program.
//! - [`myers_bit_parallel`] — Myers' (1999) bit-vector algorithm,
//!   which packs one DP column into machine words and advances it with
//!   a handful of bitwise operations, giving O(n · m / w) time
//!   (`w = 64`). It computes the *same* distance, just faster.
//!
//! [`levenshtein_bounded`] adds an early-exit threshold for the common
//! "are these within `k` edits?" question.

/// Levenshtein edit distance via the classic two-row DP.
///
/// Runs in O(`a.len()` · `b.len()`) time and O(min(n, m)) space.
pub fn levenshtein(a: &[u8], b: &[u8]) -> usize {
    // Iterate with the shorter string as the inner (row) dimension.
    let (a, b) = if a.len() <= b.len() { (a, b) } else { (b, a) };
    let n = a.len();
    let mut prev: Vec<usize> = (0..=n).collect();
    let mut cur = vec![0usize; n + 1];
    for (i, &cb) in b.iter().enumerate() {
        cur[0] = i + 1;
        for (j, &ca) in a.iter().enumerate() {
            let cost = usize::from(!ca.eq_ignore_ascii_case(&cb));
            cur[j + 1] = (prev[j] + cost)
                .min(prev[j + 1] + 1)
                .min(cur[j] + 1);
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    prev[n]
}

/// Levenshtein distance with an early-exit threshold `k`.
///
/// Returns `Some(distance)` if the true distance is `<= k`, otherwise
/// `None`. A length difference greater than `k` short-circuits
/// immediately; each DP row is abandoned once its minimum exceeds `k`.
pub fn levenshtein_bounded(a: &[u8], b: &[u8], k: usize) -> Option<usize> {
    let (a, b) = if a.len() <= b.len() { (a, b) } else { (b, a) };
    let n = a.len();
    let m = b.len();
    if m - n > k {
        return None;
    }
    let mut prev: Vec<usize> = (0..=n).collect();
    let mut cur = vec![0usize; n + 1];
    for (i, &cb) in b.iter().enumerate() {
        cur[0] = i + 1;
        let mut row_min = cur[0];
        for (j, &ca) in a.iter().enumerate() {
            let cost = usize::from(!ca.eq_ignore_ascii_case(&cb));
            let v = (prev[j] + cost).min(prev[j + 1] + 1).min(cur[j] + 1);
            cur[j + 1] = v;
            row_min = row_min.min(v);
        }
        if row_min > k {
            return None; // every cell already exceeds the budget
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    let d = prev[n];
    if d <= k {
        Some(d)
    } else {
        None
    }
}

/// Myers' bit-parallel edit distance.
///
/// Computes the Levenshtein distance of `pattern` against `text` with
/// the bit-vector recurrence of Myers (1999), "A fast bit-vector
/// algorithm for approximate string matching based on dynamic
/// programming". Time is O(`text.len()` · `ceil(pattern.len()/64)`).
///
/// The result is identical to [`levenshtein`]; this routine is the
/// fast path when one string is long.
pub fn myers_bit_parallel(pattern: &[u8], text: &[u8]) -> usize {
    let m = pattern.len();
    if m == 0 {
        return text.len();
    }
    if m <= 64 {
        myers_single_word(pattern, text)
    } else {
        // Fall back to the DP for patterns wider than one machine word
        // — a multi-word Myers implementation is a documented future
        // optimisation; correctness is unaffected.
        levenshtein(pattern, text)
    }
}

/// Single-word (`m <= 64`) Myers bit-parallel kernel.
fn myers_single_word(pattern: &[u8], text: &[u8]) -> usize {
    let m = pattern.len();
    debug_assert!(m <= 64);

    // Peq[c] has bit j set iff pattern[j] == c. Keyed by byte value.
    let mut peq = [0u64; 256];
    for (j, &c) in pattern.iter().enumerate() {
        peq[c.to_ascii_uppercase() as usize] |= 1u64 << j;
    }

    let mut pv: u64 = !0u64; // positive vertical delta = all ones
    let mut mv: u64 = 0; // negative vertical delta
    let top_bit = 1u64 << (m - 1);
    let mut score = m; // distance at column 0 == |pattern|

    for &c in text {
        let eq = peq[c.to_ascii_uppercase() as usize];
        // Myers' recurrence.
        let xv = eq | mv;
        let xh = (((eq & pv).wrapping_add(pv)) ^ pv) | eq;
        let mut ph = mv | !(xh | pv);
        let mut mh = pv & xh;
        // Update the running score from the high bit of the
        // horizontal deltas.
        if ph & top_bit != 0 {
            score += 1;
        } else if mh & top_bit != 0 {
            score -= 1;
        }
        // Shift the horizontal deltas into the vertical ones.
        ph = (ph << 1) | 1;
        mh <<= 1;
        pv = mh | !(xv | ph);
        mv = ph & xv;
    }
    score
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn levenshtein_known_cases() {
        assert_eq!(levenshtein(b"", b""), 0);
        assert_eq!(levenshtein(b"abc", b"abc"), 0);
        assert_eq!(levenshtein(b"abc", b""), 3);
        // "kitten" -> "sitting": classic distance 3.
        assert_eq!(levenshtein(b"kitten", b"sitting"), 3);
        // One substitution.
        assert_eq!(levenshtein(b"ACGT", b"ACTT"), 1);
        // One insertion.
        assert_eq!(levenshtein(b"ACGT", b"ACGGT"), 1);
    }

    #[test]
    fn levenshtein_is_symmetric() {
        let a = b"GATTACA";
        let b = b"GCATGCU";
        assert_eq!(levenshtein(a, b), levenshtein(b, a));
    }

    #[test]
    fn levenshtein_case_insensitive() {
        assert_eq!(levenshtein(b"ACGT", b"acgt"), 0);
    }

    #[test]
    fn bounded_early_exit() {
        // Within budget.
        assert_eq!(levenshtein_bounded(b"kitten", b"sitting", 3), Some(3));
        assert_eq!(levenshtein_bounded(b"kitten", b"sitting", 5), Some(3));
        // Over budget.
        assert_eq!(levenshtein_bounded(b"kitten", b"sitting", 2), None);
        // Length difference alone exceeds the budget.
        assert_eq!(levenshtein_bounded(b"AAAA", b"AAAAAAAAAA", 2), None);
    }

    #[test]
    fn myers_matches_levenshtein_short() {
        let cases: &[(&[u8], &[u8])] = &[
            (b"", b"ACGT"),
            (b"ACGT", b"ACGT"),
            (b"ACGT", b"ACTT"),
            (b"kitten", b"sitting"),
            (b"GATTACA", b"GCATGCU"),
            (b"AAAA", b"AAAAAAAA"),
        ];
        for &(p, t) in cases {
            assert_eq!(
                myers_bit_parallel(p, t),
                levenshtein(p, t),
                "Myers != Levenshtein on {:?}/{:?}",
                std::str::from_utf8(p),
                std::str::from_utf8(t),
            );
        }
    }

    #[test]
    fn myers_matches_levenshtein_long() {
        // A pattern longer than one machine word exercises the
        // documented DP fallback; still must agree.
        let pattern: Vec<u8> = b"ACGT".iter().cycle().take(100).copied().collect();
        let mut text = pattern.clone();
        text[50] = b'A'; // single edit somewhere in the middle
        text[51] = b'A';
        assert_eq!(
            myers_bit_parallel(&pattern, &text),
            levenshtein(&pattern, &text)
        );
    }

    #[test]
    fn myers_single_word_boundary() {
        // Exactly 64-residue pattern — the single-word kernel's edge.
        let pattern: Vec<u8> = b"AC".iter().cycle().take(64).copied().collect();
        let text = pattern.clone();
        assert_eq!(myers_bit_parallel(&pattern, &text), 0);
        let mut text2 = pattern.clone();
        text2[10] = b'G';
        assert_eq!(
            myers_bit_parallel(&pattern, &text2),
            levenshtein(&pattern, &text2)
        );
    }
}
