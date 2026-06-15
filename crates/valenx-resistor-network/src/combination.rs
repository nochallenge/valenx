//! Series and parallel equivalent-resistance reductions.
//!
//! ## Model
//!
//! For ideal lumped resistors `R1, R2, ... Rn` (each a real,
//! strictly-positive ohmic value):
//!
//! - Series: `R_eq = R1 + R2 + ... + Rn`.
//! - Parallel: `1 / R_eq = 1/R1 + 1/R2 + ... + 1/Rn`, equivalently
//!   `R_eq = 1 / (sum of conductances)`.
//!
//! Two named consequences the unit tests pin down:
//!
//! - Series resistance is at least the largest element, so adding a
//!   resistor in series never lowers the total.
//! - Parallel resistance is at most the smallest element, so adding
//!   a resistor in parallel never raises the total; two equal
//!   resistors `R` in parallel give exactly `R / 2`.

use crate::error::{check_positive, ResistorError};

/// Equivalent resistance of resistors wired in series.
///
/// Computes `R_eq = sum(Ri)` over `resistors`. Each entry must be a
/// finite, strictly-positive resistance.
///
/// # Errors
///
/// - [`ResistorError::EmptyNetwork`] if `resistors` is empty.
/// - [`ResistorError::NonPositive`] / [`ResistorError::NonFinite`]
///   if any element is out of domain.
///
/// # Examples
///
/// ```
/// use valenx_resistor_network::combination::series;
/// let r = series(&[100.0, 220.0, 330.0]).unwrap();
/// assert!((r - 650.0).abs() < 1e-9);
/// ```
pub fn series(resistors: &[f64]) -> Result<f64, ResistorError> {
    if resistors.is_empty() {
        return Err(ResistorError::empty_network());
    }
    let mut total = 0.0_f64;
    for &r in resistors {
        total += check_positive("resistance", r)?;
    }
    Ok(total)
}

/// Equivalent resistance of resistors wired in parallel.
///
/// Computes `R_eq = 1 / sum(1/Ri)` over `resistors`. Each entry
/// must be a finite, strictly-positive resistance.
///
/// # Errors
///
/// - [`ResistorError::EmptyNetwork`] if `resistors` is empty.
/// - [`ResistorError::NonPositive`] / [`ResistorError::NonFinite`]
///   if any element is out of domain.
///
/// # Examples
///
/// ```
/// use valenx_resistor_network::combination::parallel;
/// // Two equal resistors in parallel halve the resistance.
/// let r = parallel(&[1000.0, 1000.0]).unwrap();
/// assert!((r - 500.0).abs() < 1e-9);
/// ```
pub fn parallel(resistors: &[f64]) -> Result<f64, ResistorError> {
    if resistors.is_empty() {
        return Err(ResistorError::empty_network());
    }
    let mut conductance = 0.0_f64;
    for &r in resistors {
        conductance += 1.0 / check_positive("resistance", r)?;
    }
    Ok(1.0 / conductance)
}

/// Equivalent resistance of exactly two resistors in parallel.
///
/// A convenience wrapper that uses the product-over-sum form
/// `R_eq = (R1 * R2) / (R1 + R2)`, which is numerically identical to
/// [`parallel`] for two elements and is the form most often quoted in
/// textbooks. Both inputs must be finite and strictly positive.
///
/// # Errors
///
/// [`ResistorError::NonPositive`] / [`ResistorError::NonFinite`] if
/// either input is out of domain.
///
/// # Examples
///
/// ```
/// use valenx_resistor_network::combination::parallel_pair;
/// // 2 kohm || 3 kohm = 1.2 kohm.
/// let r = parallel_pair(2000.0, 3000.0).unwrap();
/// assert!((r - 1200.0).abs() < 1e-9);
/// ```
pub fn parallel_pair(r1: f64, r2: f64) -> Result<f64, ResistorError> {
    let r1 = check_positive("resistance", r1)?;
    let r2 = check_positive("resistance", r2)?;
    Ok((r1 * r2) / (r1 + r2))
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn series_adds() {
        // Ground truth: 100 + 220 + 330 = 650.
        let r = series(&[100.0, 220.0, 330.0]).expect("valid");
        assert!((r - 650.0).abs() < EPS, "got {r}");
    }

    #[test]
    fn series_single_resistor_is_itself() {
        let r = series(&[470.0]).expect("valid");
        assert!((r - 470.0).abs() < EPS, "got {r}");
    }

    #[test]
    fn series_is_at_least_the_largest_element() {
        let parts = [10.0, 47.0, 5.0, 22.0];
        let max = 47.0_f64;
        let r = series(&parts).expect("valid");
        assert!(r >= max, "series {r} should be >= largest {max}");
        // And strictly greater once a second resistor is present.
        assert!(r > max);
    }

    #[test]
    fn parallel_two_equal_halves() {
        // Ground truth: R || R = R / 2.
        let r = parallel(&[1000.0, 1000.0]).expect("valid");
        assert!((r - 500.0).abs() < EPS, "got {r}");
    }

    #[test]
    fn parallel_classic_two_resistor() {
        // Ground truth: 2k || 3k = 6e6 / 5e3 = 1200.
        let r = parallel(&[2000.0, 3000.0]).expect("valid");
        assert!((r - 1200.0).abs() < EPS, "got {r}");
    }

    #[test]
    fn parallel_three_equal_is_third() {
        // Ground truth: R || R || R = R / 3.
        let r = parallel(&[300.0, 300.0, 300.0]).expect("valid");
        assert!((r - 100.0).abs() < EPS, "got {r}");
    }

    #[test]
    fn parallel_is_below_the_smallest_element() {
        let parts = [100.0, 220.0, 47.0, 330.0];
        let min = 47.0_f64;
        let r = parallel(&parts).expect("valid");
        assert!(r < min, "parallel {r} should be < smallest {min}");
    }

    #[test]
    fn parallel_single_resistor_is_itself() {
        let r = parallel(&[820.0]).expect("valid");
        assert!((r - 820.0).abs() < EPS, "got {r}");
    }

    #[test]
    fn parallel_pair_matches_slice_parallel() {
        let a = parallel_pair(2000.0, 3000.0).expect("valid");
        let b = parallel(&[2000.0, 3000.0]).expect("valid");
        assert!((a - b).abs() < EPS, "pair {a} vs slice {b}");
    }

    #[test]
    fn parallel_pair_equal_halves() {
        let r = parallel_pair(680.0, 680.0).expect("valid");
        assert!((r - 340.0).abs() < EPS, "got {r}");
    }

    #[test]
    fn empty_is_rejected() {
        assert_eq!(series(&[]), Err(ResistorError::empty_network()));
        assert_eq!(parallel(&[]), Err(ResistorError::empty_network()));
    }

    #[test]
    fn non_positive_is_rejected() {
        assert_eq!(
            series(&[100.0, 0.0]),
            Err(ResistorError::non_positive("resistance", 0.0))
        );
        assert_eq!(
            parallel(&[100.0, -10.0]),
            Err(ResistorError::non_positive("resistance", -10.0))
        );
        assert_eq!(
            parallel_pair(-1.0, 5.0),
            Err(ResistorError::non_positive("resistance", -1.0))
        );
    }

    #[test]
    fn non_finite_is_rejected() {
        match series(&[f64::INFINITY]) {
            Err(ResistorError::NonFinite { name, .. }) => assert_eq!(name, "resistance"),
            other => panic!("expected NonFinite, got {other:?}"),
        }
        match parallel(&[f64::NAN, 100.0]) {
            Err(ResistorError::NonFinite { name, .. }) => assert_eq!(name, "resistance"),
            other => panic!("expected NonFinite, got {other:?}"),
        }
    }

    #[test]
    fn series_of_parallel_pairs_known_value() {
        // Two 100-ohm in parallel (= 50) in series with two 200-ohm
        // in parallel (= 100) gives 150 ohm.
        let left = parallel(&[100.0, 100.0]).expect("valid");
        let right = parallel(&[200.0, 200.0]).expect("valid");
        let total = series(&[left, right]).expect("valid");
        assert!((total - 150.0).abs() < EPS, "got {total}");
    }
}
