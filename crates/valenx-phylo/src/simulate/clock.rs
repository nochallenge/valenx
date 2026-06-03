//! Root-to-tip molecular-clock dating.
//!
//! When tips are sampled at *known, different* times (heterochronous
//! data — e.g. virus genomes from dated outbreaks), a tree's branch
//! lengths can be calibrated to absolute time by a simple regression
//! (the TempEst / "root-to-tip" method, Rambaut et al. 2016):
//!
//! - Under a strict molecular clock, a tip's **root-to-tip genetic
//!   distance** grows linearly with its **sampling time**.
//! - Regressing root-to-tip distance (`y`) on sampling time (`x`) over
//!   all tips gives a line `y = rate · x + intercept`.
//! - The **slope** is the substitution rate (substitutions per site per
//!   time unit); the **x-intercept** (`-intercept / rate`) estimates
//!   the age of the root (the time of the most recent common
//!   ancestor).
//! - The regression `R²` measures how clock-like the data are.
//!
//! This is a fast point estimate, not a relaxed-clock Bayesian analysis
//! (BEAST's domain) — but it is exactly the standard first-pass dating
//! diagnostic.

use crate::error::{PhyloError, Result};
use crate::tree::Tree;

/// The fitted molecular-clock regression.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClockEstimate {
    /// Substitution rate — the regression slope (substitutions per site
    /// per unit time).
    pub rate: f64,
    /// Regression intercept (root-to-tip distance at time 0).
    pub intercept: f64,
    /// Estimated age (sampling-time value) of the root — the line's
    /// x-intercept, `-intercept / rate`.
    pub root_time: f64,
    /// Coefficient of determination `R²` in `[0, 1]` — how clock-like
    /// the data are.
    pub r_squared: f64,
    /// Number of tips used in the regression.
    pub n_tips: usize,
}

/// Fits a root-to-tip molecular-clock regression.
///
/// `tree` supplies the branch lengths; `sampling_times` maps each leaf
/// label to its sampling time (any consistent unit — calendar years,
/// days since an epoch, …). Every leaf of the tree must have an entry.
///
/// # Errors
/// - [`PhyloError::Invalid`] if a leaf has no sampling time, fewer than
///   two distinct sampling times are present (the slope is then
///   undefined), or the tree has fewer than two leaves.
pub fn root_to_tip_regression(
    tree: &Tree,
    sampling_times: &[(String, f64)],
) -> Result<ClockEstimate> {
    let leaves = tree.leaves();
    if leaves.len() < 2 {
        return Err(PhyloError::invalid("tree", "need at least two tips"));
    }
    let root = tree.root();

    // Build the (sampling_time, root_to_tip_distance) point cloud.
    let mut xs = Vec::with_capacity(leaves.len());
    let mut ys = Vec::with_capacity(leaves.len());
    for &leaf in &leaves {
        let label = tree
            .node(leaf)
            .label
            .as_deref()
            .ok_or_else(|| PhyloError::invalid("tree", "leaf without a label"))?;
        let time = sampling_times
            .iter()
            .find(|(name, _)| name == label)
            .map(|(_, t)| *t)
            .ok_or_else(|| {
                PhyloError::invalid(
                    "sampling_times",
                    format!("no sampling time for `{label}`"),
                )
            })?;
        xs.push(time);
        ys.push(tree.patristic_distance(root, leaf));
    }

    // Need spread in x for a defined slope.
    let x_min = xs.iter().cloned().fold(f64::INFINITY, f64::min);
    let x_max = xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if (x_max - x_min).abs() < 1e-12 {
        return Err(PhyloError::invalid(
            "sampling_times",
            "all tips share one sampling time — the clock rate is undefined",
        ));
    }

    let (rate, intercept, r_squared) = ordinary_least_squares(&xs, &ys);
    // x-intercept: where the fitted line crosses y = 0.
    let root_time = if rate.abs() > 1e-15 {
        -intercept / rate
    } else {
        f64::NAN
    };

    Ok(ClockEstimate {
        rate,
        intercept,
        root_time,
        r_squared,
        n_tips: leaves.len(),
    })
}

/// Ordinary least-squares fit of `y = slope·x + intercept`.
///
/// Returns `(slope, intercept, r_squared)`. `r_squared` is 1.0 when `y`
/// has zero variance (a degenerate but well-defined perfect fit).
fn ordinary_least_squares(xs: &[f64], ys: &[f64]) -> (f64, f64, f64) {
    let n = xs.len() as f64;
    let mean_x = xs.iter().sum::<f64>() / n;
    let mean_y = ys.iter().sum::<f64>() / n;

    let mut sxx = 0.0;
    let mut sxy = 0.0;
    let mut syy = 0.0;
    for (&x, &y) in xs.iter().zip(ys) {
        let dx = x - mean_x;
        let dy = y - mean_y;
        sxx += dx * dx;
        sxy += dx * dy;
        syy += dy * dy;
    }
    let slope = if sxx.abs() > 1e-15 { sxy / sxx } else { 0.0 };
    let intercept = mean_y - slope * mean_x;
    // R² = (explained variance) / (total variance).
    let r_squared = if syy.abs() > 1e-15 {
        (sxy * sxy / (sxx * syy)).clamp(0.0, 1.0)
    } else {
        1.0
    };
    (slope, intercept, r_squared)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::newick::read_newick;

    #[test]
    fn perfect_clock_recovers_the_rate() {
        // An ultrametric tree where all tips are at distance 1.0 from
        // the root, but sampled at different times — a contrived but
        // exact clock test using an explicit dated tree.
        // Tip A sampled at t=10 with root-to-tip 0.5; B at t=20 with
        // 1.0 — slope 0.05.
        let tree = read_newick("(A:0.5,B:1.0);").unwrap();
        let times = vec![("A".to_string(), 10.0), ("B".to_string(), 20.0)];
        let est = root_to_tip_regression(&tree, &times).unwrap();
        // Two points => a perfect line.
        assert!((est.rate - 0.05).abs() < 1e-9, "rate = {}", est.rate);
        assert!((est.r_squared - 1.0).abs() < 1e-9);
        assert_eq!(est.n_tips, 2);
    }

    #[test]
    fn root_time_is_the_x_intercept() {
        // Line through (10, 0.5) and (20, 1.0): y = 0.05x + 0.0, so
        // the x-intercept (root age) is 0.
        let tree = read_newick("(A:0.5,B:1.0);").unwrap();
        let times = vec![("A".to_string(), 10.0), ("B".to_string(), 20.0)];
        let est = root_to_tip_regression(&tree, &times).unwrap();
        assert!(est.root_time.abs() < 1e-6, "root_time = {}", est.root_time);
    }

    #[test]
    fn positive_rate_for_increasing_distance_with_time() {
        let tree =
            read_newick("(((A:0.1,B:0.2):0.1,C:0.4):0.1,D:0.6);").unwrap();
        let times = vec![
            ("A".to_string(), 2000.0),
            ("B".to_string(), 2001.0),
            ("C".to_string(), 2003.0),
            ("D".to_string(), 2005.0),
        ];
        let est = root_to_tip_regression(&tree, &times).unwrap();
        assert!(est.rate > 0.0, "expected a positive rate");
        assert!((0.0..=1.0).contains(&est.r_squared));
    }

    #[test]
    fn rejects_isochronous_sampling() {
        // Every tip sampled at the same time => undefined slope.
        let tree = read_newick("((A:0.1,B:0.1):0.1,C:0.2);").unwrap();
        let times = vec![
            ("A".to_string(), 2020.0),
            ("B".to_string(), 2020.0),
            ("C".to_string(), 2020.0),
        ];
        assert!(root_to_tip_regression(&tree, &times).is_err());
    }

    #[test]
    fn rejects_a_missing_sampling_time() {
        let tree = read_newick("((A,B),C);").unwrap();
        let times = vec![("A".to_string(), 1.0), ("B".to_string(), 2.0)];
        assert!(root_to_tip_regression(&tree, &times).is_err());
    }

    #[test]
    fn ols_fits_a_known_line() {
        // y = 3x + 1 exactly.
        let xs = [0.0, 1.0, 2.0, 3.0, 4.0];
        let ys = [1.0, 4.0, 7.0, 10.0, 13.0];
        let (slope, intercept, r2) = ordinary_least_squares(&xs, &ys);
        assert!((slope - 3.0).abs() < 1e-9);
        assert!((intercept - 1.0).abs() < 1e-9);
        assert!((r2 - 1.0).abs() < 1e-9);
    }
}
