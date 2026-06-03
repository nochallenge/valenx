//! Energy dot plot (mfold-class).
//!
//! A *dot plot* is the classic 2-D summary of an RNA folding
//! landscape: an upper-triangular matrix where cell `(i, j)` carries
//! a value for the pair `i·j`. mfold renders two kinds:
//!
//! - a **probability dot plot** — cell `(i, j)` is the base-pair
//!   probability `p(i,j)` from the McCaskill ensemble;
//! - an **MFE dot plot** — cell `(i, j)` marks the pairs of the
//!   single minimum-free-energy structure.
//!
//! Both are produced here as a [`DotPlot`] matrix, the data a heat-map
//! / dot-plot rendering layer consumes.

use crate::ensemble::partition::PartitionFunction;
use crate::structure::Structure;

/// An upper-triangular dot-plot matrix.
#[derive(Clone, Debug, PartialEq)]
pub struct DotPlot {
    n: usize,
    /// Flat `n×n` matrix; only `cells[i*n+j]` for `i < j` is
    /// meaningful, the rest is 0.
    cells: Vec<f64>,
}

impl DotPlot {
    /// Sequence length (the matrix is `n × n`).
    pub fn len(&self) -> usize {
        self.n
    }

    /// `true` if the plot spans zero positions.
    pub fn is_empty(&self) -> bool {
        self.n == 0
    }

    /// The value of cell `(i, j)` (order-insensitive; 0 for `i == j`
    /// or out of range).
    pub fn value(&self, i: usize, j: usize) -> f64 {
        if i >= self.n || j >= self.n || i == j {
            return 0.0;
        }
        let (a, b) = if i < j { (i, j) } else { (j, i) };
        self.cells[a * self.n + b]
    }

    /// The flat `n×n` matrix backing the plot.
    pub fn matrix(&self) -> &[f64] {
        &self.cells
    }

    /// All non-zero cells as `(i, j, value)` triples — the "dots".
    pub fn dots(&self) -> Vec<(usize, usize, f64)> {
        let mut out = Vec::new();
        for i in 0..self.n {
            for j in (i + 1)..self.n {
                let v = self.cells[i * self.n + j];
                if v != 0.0 {
                    out.push((i, j, v));
                }
            }
        }
        out
    }

    /// The largest cell value in the plot.
    pub fn max_value(&self) -> f64 {
        self.cells.iter().cloned().fold(0.0, f64::max)
    }
}

/// Builds a *probability* dot plot from a McCaskill partition
/// function: cell `(i, j)` is the base-pair probability `p(i,j)`.
pub fn probability_dot_plot(pf: &PartitionFunction) -> DotPlot {
    let n = pf.len();
    let mut cells = vec![0.0_f64; n * n];
    for i in 0..n {
        for j in (i + 1)..n {
            cells[i * n + j] = pf.pair_probability(i, j);
        }
    }
    DotPlot { n, cells }
}

/// Builds an *MFE* dot plot from a single structure: cell `(i, j)` is
/// `1.0` if `i` pairs `j` in the structure, `0.0` otherwise.
pub fn mfe_dot_plot(s: &Structure) -> DotPlot {
    let n = s.len();
    let mut cells = vec![0.0_f64; n * n];
    for bp in s.pairs() {
        cells[bp.i * n + bp.j] = 1.0;
    }
    DotPlot { n, cells }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ensemble::partition::partition_function;
    use crate::rna::RnaSeq;

    #[test]
    fn mfe_dot_plot_marks_the_pairs() {
        let s = Structure::from_dot_bracket("(((...)))").unwrap();
        let dp = mfe_dot_plot(&s);
        assert_eq!(dp.len(), 9);
        assert_eq!(dp.value(0, 8), 1.0);
        assert_eq!(dp.value(1, 7), 1.0);
        assert_eq!(dp.value(0, 5), 0.0); // not a pair
        assert_eq!(dp.dots().len(), 3);
    }

    #[test]
    fn dot_plot_is_order_insensitive() {
        let s = Structure::from_dot_bracket("(((...)))").unwrap();
        let dp = mfe_dot_plot(&s);
        assert_eq!(dp.value(0, 8), dp.value(8, 0));
    }

    #[test]
    fn probability_dot_plot_values_in_range() {
        let seq = RnaSeq::parse("GGGGGAAAACCCCC").unwrap();
        let pf = partition_function(&seq).unwrap();
        let dp = probability_dot_plot(&pf);
        for (_, _, v) in dp.dots() {
            assert!((0.0..=1.0).contains(&v));
        }
        assert!(dp.max_value() <= 1.0 + 1e-9);
    }

    #[test]
    fn empty_structure_has_no_dots() {
        let dp = mfe_dot_plot(&Structure::empty(6));
        assert!(dp.dots().is_empty());
        assert_eq!(dp.max_value(), 0.0);
    }

    #[test]
    fn probability_plot_of_stable_stem_has_dots() {
        let seq = RnaSeq::parse("GGGGGGAAAACCCCCC").unwrap();
        let pf = partition_function(&seq).unwrap();
        let dp = probability_dot_plot(&pf);
        assert!(!dp.dots().is_empty(), "a stable stem should leave dots");
    }
}
