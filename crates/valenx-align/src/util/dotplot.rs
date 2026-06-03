//! Dot-plot / self-alignment matrix — windowed identity grid.
//!
//! A **dot-plot** compares two sequences (or one against itself) by
//! placing a dot at `(i, j)` wherever a window starting at `x[i]` and
//! `y[j]` is sufficiently similar. Diagonal runs of dots reveal
//! conserved segments; an off-diagonal run reveals a repeat or an
//! inversion. It is the oldest sequence-comparison visualisation and
//! still the quickest way to *see* repeats and rearrangements.
//!
//! [`dot_plot`] builds the boolean matrix for a `(window, stringency)`
//! filter; [`DotPlot::diagonals`] extracts the maximal diagonal runs
//! (the "lines" a renderer would draw); [`DotPlot::to_ascii`] renders
//! it as text.

use crate::error::{AlignError, Result};
use crate::limits::{check_dp_size_with, MAX_DP_CELLS};

/// A windowed-identity dot-plot matrix.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DotPlot {
    /// Width — number of window positions along sequence X.
    rows: usize,
    /// Height — number of window positions along sequence Y.
    cols: usize,
    /// `cell[i*cols + j]` is `true` if window `(i, j)` passed the
    /// similarity filter.
    cells: Vec<bool>,
}

impl DotPlot {
    /// Number of X-axis window positions.
    pub fn rows(&self) -> usize {
        self.rows
    }

    /// Number of Y-axis window positions.
    pub fn cols(&self) -> usize {
        self.cols
    }

    /// The dot at `(i, j)`; `false` if out of range.
    pub fn get(&self, i: usize, j: usize) -> bool {
        if i < self.rows && j < self.cols {
            self.cells[i * self.cols + j]
        } else {
            false
        }
    }

    /// Total number of set dots.
    pub fn dot_count(&self) -> usize {
        self.cells.iter().filter(|&&d| d).count()
    }

    /// The fraction of cells that are dots, in `[0, 1]`.
    pub fn density(&self) -> f64 {
        let total = self.rows * self.cols;
        if total == 0 {
            0.0
        } else {
            self.dot_count() as f64 / total as f64
        }
    }

    /// The maximal *diagonal* runs of dots — each a
    /// [`DiagonalRun`] of length at least `min_len`. A renderer draws
    /// one line segment per run. Runs are found on every diagonal,
    /// including off-diagonal ones (repeats).
    pub fn diagonals(&self, min_len: usize) -> Vec<DiagonalRun> {
        let min_len = min_len.max(1);
        let mut runs = Vec::new();
        // Each diagonal is identified by d = j - i.
        for d in -(self.rows as isize - 1)..(self.cols as isize) {
            let mut run_start: Option<(usize, usize)> = None;
            let mut run_len = 0usize;
            // Walk the diagonal.
            let mut i = if d < 0 { (-d) as usize } else { 0 };
            let mut j = if d < 0 { 0 } else { d as usize };
            while i < self.rows && j < self.cols {
                if self.get(i, j) {
                    if run_start.is_none() {
                        run_start = Some((i, j));
                    }
                    run_len += 1;
                } else {
                    if run_len >= min_len {
                        let (si, sj) = run_start.unwrap();
                        runs.push(DiagonalRun {
                            start: (si, sj),
                            len: run_len,
                        });
                    }
                    run_start = None;
                    run_len = 0;
                }
                i += 1;
                j += 1;
            }
            if run_len >= min_len {
                let (si, sj) = run_start.unwrap();
                runs.push(DiagonalRun {
                    start: (si, sj),
                    len: run_len,
                });
            }
        }
        runs
    }

    /// Renders the dot-plot as ASCII — `#` for a dot, `.` for empty.
    /// Row 0 is the top line.
    pub fn to_ascii(&self) -> String {
        let mut out = String::with_capacity((self.cols + 1) * self.rows);
        for i in 0..self.rows {
            for j in 0..self.cols {
                out.push(if self.get(i, j) { '#' } else { '.' });
            }
            out.push('\n');
        }
        out
    }
}

/// A maximal run of dots along one diagonal of a [`DotPlot`].
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct DiagonalRun {
    /// `(i, j)` of the first dot in the run.
    pub start: (usize, usize),
    /// Number of consecutive dots.
    pub len: usize,
}

impl DiagonalRun {
    /// `true` if the run is on the main diagonal (`i == j`).
    pub fn is_main_diagonal(&self) -> bool {
        self.start.0 == self.start.1
    }
}

/// Builds a windowed-identity dot-plot of `x` against `y`.
///
/// For every pair of `window`-length windows `(x[i..i+w], y[j..j+w])`
/// the number of matching positions is counted; the cell is set when
/// that count is `>= stringency`. A `window` of `1` with `stringency`
/// `1` is the raw single-residue dot-plot.
///
/// Returns [`AlignError::Invalid`] for `window == 0` or
/// `stringency > window`, and
/// [`AlignError::TooLarge`] when the
/// `rows × cols` cell matrix would exceed
/// [`MAX_DP_CELLS`].
pub fn dot_plot(x: &[u8], y: &[u8], window: usize, stringency: usize) -> Result<DotPlot> {
    dot_plot_capped(x, y, window, stringency, MAX_DP_CELLS)
}

/// [`dot_plot`] with an explicit cell cap (test seam).
fn dot_plot_capped(
    x: &[u8],
    y: &[u8],
    window: usize,
    stringency: usize,
    max_cells: usize,
) -> Result<DotPlot> {
    if window == 0 {
        return Err(AlignError::invalid("window", "window must be >= 1"));
    }
    if stringency > window {
        return Err(AlignError::invalid(
            "stringency",
            format!("stringency {stringency} exceeds window {window}"),
        ));
    }
    let rows = x.len().saturating_sub(window) + 1;
    let cols = y.len().saturating_sub(window) + 1;
    // If either sequence is shorter than the window the plot is empty.
    let (rows, cols) = if x.len() < window || y.len() < window {
        (0, 0)
    } else {
        (rows, cols)
    };

    // Bound the boolean-matrix allocation (DoS guard for two long
    // sequences).
    check_dp_size_with(rows, cols, max_cells)?;

    let mut cells = vec![false; rows * cols];
    for i in 0..rows {
        for j in 0..cols {
            let matches = (0..window)
                .filter(|&k| x[i + k].eq_ignore_ascii_case(&y[j + k]))
                .count();
            if matches >= stringency {
                cells[i * cols + j] = true;
            }
        }
    }
    Ok(DotPlot { rows, cols, cells })
}

/// Builds a *self* dot-plot — `seq` against itself. The main diagonal
/// is always fully set; off-diagonal runs are internal repeats.
pub fn self_dot_plot(seq: &[u8], window: usize, stringency: usize) -> Result<DotPlot> {
    dot_plot(seq, seq, window, stringency)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_bad_params() {
        assert!(dot_plot(b"ACGT", b"ACGT", 0, 0).is_err());
        assert!(dot_plot(b"ACGT", b"ACGT", 3, 5).is_err());
    }

    #[test]
    fn identical_sequences_fill_main_diagonal() {
        let dp = dot_plot(b"ACGTACGT", b"ACGTACGT", 1, 1).unwrap();
        assert_eq!(dp.rows(), 8);
        assert_eq!(dp.cols(), 8);
        // Every main-diagonal cell is a dot.
        for i in 0..8 {
            assert!(dp.get(i, i), "main diagonal cell {i} should be set");
        }
    }

    #[test]
    fn unrelated_sequences_sparse() {
        let dp = dot_plot(b"AAAAAAAA", b"CCCCCCCC", 1, 1).unwrap();
        assert_eq!(dp.dot_count(), 0);
        assert_eq!(dp.density(), 0.0);
    }

    #[test]
    fn window_filtering_reduces_noise() {
        // Single-residue plot of two similar-but-noisy sequences is
        // busy; a wider window with full stringency is sparser.
        let x = b"ACGTACGTACGT";
        let y = b"ACGTACGTACGT";
        let raw = dot_plot(x, y, 1, 1).unwrap();
        let filtered = dot_plot(x, y, 4, 4).unwrap();
        assert!(filtered.dot_count() <= raw.dot_count());
    }

    #[test]
    fn diagonal_runs_detected() {
        // Identical 10-mers: one long main-diagonal run.
        let dp = dot_plot(b"ACGTACGTAC", b"ACGTACGTAC", 1, 1).unwrap();
        let runs = dp.diagonals(5);
        // The main-diagonal run of length 10 must be present.
        assert!(runs.iter().any(|r| r.is_main_diagonal() && r.len == 10));
    }

    #[test]
    fn self_dotplot_reveals_repeat() {
        // "ABCABC"-style sequence: an off-diagonal run marks the repeat.
        let dp = self_dot_plot(b"ACGTACGT", 1, 1).unwrap();
        let runs = dp.diagonals(3);
        // There is at least one off-main-diagonal run (the repeat).
        assert!(runs.iter().any(|r| !r.is_main_diagonal()));
    }

    #[test]
    fn short_sequence_empty_plot() {
        // Window wider than the sequence -> empty plot.
        let dp = dot_plot(b"ACG", b"ACG", 5, 5).unwrap();
        assert_eq!(dp.rows(), 0);
        assert_eq!(dp.dot_count(), 0);
    }

    #[test]
    fn ascii_rendering() {
        let dp = dot_plot(b"AC", b"AC", 1, 1).unwrap();
        let art = dp.to_ascii();
        // 2x2 grid: diagonal dots.
        assert!(art.contains('#'));
        assert_eq!(art.lines().count(), 2);
    }

    #[test]
    fn dotplot_over_cap_errors() {
        use crate::error::AlignError;
        // window 1 over 8x8 sequences -> 8x8 = 64 cells; cap of 8 rejects
        // without allocating the boolean matrix.
        let err = dot_plot_capped(b"ACGTACGT", b"ACGTACGT", 1, 1, 8).unwrap_err();
        assert!(matches!(err, AlignError::TooLarge { .. }), "got {err:?}");
        // Generous cap computes normally.
        assert!(dot_plot_capped(b"ACGTACGT", b"ACGTACGT", 1, 1, usize::MAX).is_ok());
    }
}
