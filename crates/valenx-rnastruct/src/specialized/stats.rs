//! Structure statistics — loop and element counts.
//!
//! A secondary structure decomposes into a small set of element
//! types; counting them characterises the structure at a glance. This
//! module walks a (nested) structure and reports:
//!
//! - the number of **hairpin loops**, **bulge loops**, **internal
//!   loops** and **multiloops**;
//! - the number of **stems** (maximal stacks of consecutive pairs)
//!   and the total **base pairs**;
//! - the count of **unpaired bases** and the size of the largest
//!   hairpin loop.
//!
//! The decomposition is the same loop-classification rule
//! [`crate::fold::eval`] uses: every pair closes one loop, classified
//! by the number of pairs it directly encloses.

use crate::structure::Structure;

/// A tally of the structural elements of a secondary structure.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub struct StructureStats {
    /// Total base pairs.
    pub base_pairs: usize,
    /// Hairpin loops (a closing pair enclosing 0 pairs).
    pub hairpins: usize,
    /// Bulge loops (1 enclosed pair, unpaired bases on exactly one
    /// side).
    pub bulges: usize,
    /// Internal loops (1 enclosed pair, unpaired bases on both sides).
    pub internal_loops: usize,
    /// Multiloops (a closing pair enclosing ≥ 2 pairs).
    pub multiloops: usize,
    /// Stems — maximal runs of stacked (consecutive) base pairs.
    pub stems: usize,
    /// Bases not in any pair.
    pub unpaired_bases: usize,
    /// The size (unpaired bases) of the largest hairpin loop.
    pub largest_hairpin_loop: usize,
}

impl StructureStats {
    /// The total number of loops of all kinds.
    pub fn total_loops(&self) -> usize {
        self.hairpins + self.bulges + self.internal_loops + self.multiloops
    }
}

/// Computes the [`StructureStats`] of a (nested) structure.
///
/// A pseudoknotted structure is decomposed by treating its crossing
/// pairs as ordinary pairs for the stem count; loop classification
/// then describes the nested skeleton — documented here because loop
/// types are only strictly defined for a nested structure.
pub fn structure_stats(s: &Structure) -> StructureStats {
    let n = s.len();
    let partner = s.partner_array();
    let mut st = StructureStats {
        base_pairs: s.n_pairs(),
        ..Default::default()
    };

    st.unpaired_bases = (0..n).filter(|&i| partner[i].is_none()).count();

    // Loop classification: for each pair (i, j) count enclosed pairs.
    for i in 0..n {
        let Some(j) = partner[i] else { continue };
        if i >= j {
            continue;
        }
        let mut enclosed: Vec<(usize, usize)> = Vec::new();
        let mut k = i + 1;
        while k < j {
            match partner[k] {
                Some(p) if p > k => {
                    enclosed.push((k, p));
                    k = p + 1;
                }
                _ => k += 1,
            }
        }
        match enclosed.len() {
            0 => {
                st.hairpins += 1;
                let loop_size = j - i - 1;
                st.largest_hairpin_loop = st.largest_hairpin_loop.max(loop_size);
            }
            1 => {
                let (k, l) = enclosed[0];
                let left = k - i - 1;
                let right = j - l - 1;
                if left == 0 && right == 0 {
                    // a stacked pair — not a loop, not counted here
                } else if left == 0 || right == 0 {
                    st.bulges += 1;
                } else {
                    st.internal_loops += 1;
                }
            }
            _ => st.multiloops += 1,
        }
    }

    // Stem count: a stem is a maximal run of pairs (i, j), (i+1, j-1),
    // ... A pair starts a new stem if (i-1, j+1) is not also a pair.
    let mut stems = 0;
    for i in 0..n {
        let Some(j) = partner[i] else { continue };
        if i >= j {
            continue;
        }
        let extends_outward = i > 0 && j + 1 < n && partner[i - 1] == Some(j + 1);
        if !extends_outward {
            stems += 1;
        }
    }
    st.stems = stems;

    st
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_structure_has_no_elements() {
        let s = Structure::empty(10);
        let st = structure_stats(&s);
        assert_eq!(st.base_pairs, 0);
        assert_eq!(st.total_loops(), 0);
        assert_eq!(st.unpaired_bases, 10);
        assert_eq!(st.stems, 0);
    }

    #[test]
    fn simple_hairpin() {
        // (((....))) — one stem, one hairpin, 3 pairs
        let s = Structure::from_dot_bracket("(((....)))").unwrap();
        let st = structure_stats(&s);
        assert_eq!(st.base_pairs, 3);
        assert_eq!(st.hairpins, 1);
        assert_eq!(st.stems, 1);
        assert_eq!(st.unpaired_bases, 4);
        assert_eq!(st.largest_hairpin_loop, 4);
        assert_eq!(st.bulges, 0);
        assert_eq!(st.internal_loops, 0);
        assert_eq!(st.multiloops, 0);
    }

    #[test]
    fn bulge_loop() {
        // ((.((....)))) — outer 2-pair stem, a 1-nt bulge, inner stem
        let s = Structure::from_dot_bracket("((.((....))))").unwrap();
        let st = structure_stats(&s);
        assert_eq!(st.bulges, 1, "expected one bulge");
        assert_eq!(st.hairpins, 1);
        assert_eq!(st.stems, 2, "the bulge breaks the helix into 2 stems");
    }

    #[test]
    fn internal_loop() {
        // ((.((....)).)) — unpaired on both sides of the inner stem
        let s = Structure::from_dot_bracket("((.((....)).))").unwrap();
        let st = structure_stats(&s);
        assert_eq!(st.internal_loops, 1);
        assert_eq!(st.hairpins, 1);
    }

    #[test]
    fn multiloop() {
        // (((....))((....))) wrapped: an outer pair enclosing two
        // hairpins -> one multiloop
        let s = Structure::from_dot_bracket("((((....))((....))))").unwrap();
        let st = structure_stats(&s);
        assert_eq!(st.multiloops, 1, "expected one multiloop");
        assert_eq!(st.hairpins, 2, "two hairpin arms");
    }

    #[test]
    fn two_separate_hairpins() {
        // (((...)))(((...))) — two stems, two hairpins, no multiloop
        let s = Structure::from_dot_bracket("(((...)))(((...)))").unwrap();
        let st = structure_stats(&s);
        assert_eq!(st.hairpins, 2);
        assert_eq!(st.stems, 2);
        assert_eq!(st.multiloops, 0);
    }

    #[test]
    fn total_loops_sums_correctly() {
        let s = Structure::from_dot_bracket("((((....))((....))))").unwrap();
        let st = structure_stats(&s);
        assert_eq!(
            st.total_loops(),
            st.hairpins + st.bulges + st.internal_loops + st.multiloops
        );
    }
}
