//! Mountain-plot data generation.
//!
//! A *mountain plot* is a 1-D representation of a secondary
//! structure: at each sequence position `k` the "height" is the
//! number of base pairs `(i, j)` that *enclose* `k` — that is, with
//! `i < k < j` (the partners of a pair are not enclosed by their own
//! pair, so a helix forms a symmetric mountain whose stem ends sit at
//! ground level). Plotted against position the structure becomes a
//! profile of peaks (helices climbing then descending) and plateaus
//! (loops).
//! Mountain plots make two structures easy to compare by eye and a
//! probability-weighted mountain plot summarises an ensemble.
//!
//! This module produces the height vector — the data a plotting
//! layer renders. Both the single-structure mountain plot
//! ([`mountain_plot`]) and the ensemble-averaged mountain plot
//! ([`ensemble_mountain_plot`], weighted by base-pair probabilities)
//! are provided.

use crate::ensemble::partition::PartitionFunction;
use crate::structure::Structure;

/// A mountain plot — one height per sequence position.
#[derive(Clone, Debug, PartialEq)]
pub struct MountainPlot {
    /// `heights[k]` = number of pairs enclosing position `k`
    /// (or the expected number, for an ensemble plot).
    pub heights: Vec<f64>,
}

impl MountainPlot {
    /// The sequence length the plot spans.
    pub fn len(&self) -> usize {
        self.heights.len()
    }

    /// `true` if the plot is empty.
    pub fn is_empty(&self) -> bool {
        self.heights.is_empty()
    }

    /// The maximum height — the depth of the deepest nesting.
    pub fn peak(&self) -> f64 {
        self.heights.iter().cloned().fold(0.0, f64::max)
    }

    /// The area under the mountain (sum of heights) — a scalar
    /// measure of how much structure the molecule has.
    pub fn area(&self) -> f64 {
        self.heights.iter().sum()
    }
}

/// Builds the mountain plot of a single secondary structure.
///
/// `heights[k]` counts the base pairs `(i, j)` that enclose `k`,
/// i.e. with `i < k < j`. A pair's own partners are not counted, so a
/// hairpin stem climbs from ground level at `i` to the loop plateau
/// and back to ground at `j`.
pub fn mountain_plot(s: &Structure) -> MountainPlot {
    let n = s.len();
    // A pair (i, j) raises every interior position in (i, j). Use a
    // difference array for an O(n) accumulation.
    let mut delta = vec![0.0_f64; n + 1];
    for bp in s.pairs() {
        // raise positions i+1 ..= j-1
        delta[bp.i + 1] += 1.0;
        delta[bp.j] -= 1.0;
    }
    let mut heights = vec![0.0; n];
    let mut acc = 0.0;
    for (k, h) in heights.iter_mut().enumerate() {
        acc += delta[k];
        *h = acc;
    }
    MountainPlot { heights }
}

/// Builds the ensemble-averaged mountain plot from a partition
/// function: each pair contributes its base-pair *probability*
/// instead of 1, so the height at `k` is the *expected* number of
/// pairs enclosing `k`.
pub fn ensemble_mountain_plot(pf: &PartitionFunction) -> MountainPlot {
    let n = pf.len();
    let mut delta = vec![0.0_f64; n + 1];
    for i in 0..n {
        for j in (i + 1)..n {
            let p = pf.pair_probability(i, j);
            if p > 0.0 {
                // raise interior positions i+1 ..= j-1 (see mountain_plot)
                delta[i + 1] += p;
                delta[j] -= p;
            }
        }
    }
    let mut heights = vec![0.0; n];
    let mut acc = 0.0;
    for (k, h) in heights.iter_mut().enumerate() {
        acc += delta[k];
        *h = acc.max(0.0);
    }
    MountainPlot { heights }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ensemble::partition::partition_function;
    use crate::rna::RnaSeq;

    #[test]
    fn empty_structure_is_flat() {
        let s = Structure::empty(8);
        let m = mountain_plot(&s);
        assert_eq!(m.len(), 8);
        assert_eq!(m.peak(), 0.0);
        assert_eq!(m.area(), 0.0);
    }

    #[test]
    fn hairpin_makes_a_single_peak() {
        // (((....))) — height climbs to 3 in the loop, back to 0
        let s = Structure::from_dot_bracket("(((....)))").unwrap();
        let m = mountain_plot(&s);
        assert_eq!(m.peak(), 3.0, "three nested pairs -> peak height 3");
        // the ends are at ground level
        assert_eq!(m.heights[0], 0.0);
        assert_eq!(*m.heights.last().unwrap(), 0.0);
        // the loop region sits at the peak
        assert_eq!(m.heights[5], 3.0);
    }

    #[test]
    fn two_hairpins_make_two_peaks() {
        let s = Structure::from_dot_bracket("(((...)))(((...)))").unwrap();
        let m = mountain_plot(&s);
        assert_eq!(m.peak(), 3.0);
        // it returns to ground between the hairpins
        assert_eq!(m.heights[9], 0.0);
    }

    #[test]
    fn deeper_nesting_gives_higher_peak() {
        let shallow = mountain_plot(&Structure::from_dot_bracket("((....))").unwrap());
        let deep = mountain_plot(&Structure::from_dot_bracket("(((((....)))))").unwrap());
        assert!(deep.peak() > shallow.peak());
    }

    #[test]
    fn ensemble_mountain_plot_runs() {
        let seq = RnaSeq::parse("GGGGGAAAACCCCC").unwrap();
        let pf = partition_function(&seq).unwrap();
        let m = ensemble_mountain_plot(&pf);
        assert_eq!(m.len(), seq.len());
        // every height is non-negative
        for &h in &m.heights {
            assert!(h >= 0.0);
        }
    }

    #[test]
    fn ensemble_unpairable_is_flat() {
        let seq = RnaSeq::parse("AAAAAAAAAAAA").unwrap();
        let pf = partition_function(&seq).unwrap();
        let m = ensemble_mountain_plot(&pf);
        assert!(m.peak() < 1e-6, "an unpairable RNA has a flat mountain");
    }
}
