//! System reliability from component reliabilities (reliability block
//! diagrams).
//!
//! Given the reliabilities `R_i in [0, 1]` of independent components, the
//! reliability of a *system* built from them depends on its topology:
//!
//! - **Series** — the system works only if *every* component works, so
//!   `R_series = prod(R_i)`. Adding a component can never increase
//!   reliability, hence `R_series <= min(R_i)`.
//! - **Parallel** (active redundancy) — the system works if *at least
//!   one* component works, so `R_parallel = 1 - prod(1 - R_i)`. Adding a
//!   redundant path can never decrease reliability, hence
//!   `R_parallel >= max(R_i)`.
//! - **`k`-out-of-`n`** — the system works if at least `k` of its `n`
//!   identical components work; series is the `n`-of-`n` case and
//!   parallel is the `1`-of-`n` case.
//!
//! All functions assume the components fail *independently*.

use crate::error::{require_probability, ReliabilityError};

/// The reliability of a **series** system, `R = prod(R_i)`.
///
/// Every component must survive for the system to survive, so the result
/// is bounded above by the weakest component: `R <= min(R_i)`.
///
/// # Errors
///
/// Returns [`ReliabilityError::EmptySystem`] if `components` is empty, or
/// [`ReliabilityError::ProbabilityOutOfRange`] / [`ReliabilityError::NotFinite`]
/// if any supplied reliability is not a finite probability in `[0, 1]`.
pub fn series(components: &[f64]) -> Result<f64, ReliabilityError> {
    if components.is_empty() {
        return Err(ReliabilityError::EmptySystem { kind: "series" });
    }
    let mut product = 1.0;
    for &r in components {
        product *= require_probability(r)?;
    }
    Ok(product)
}

/// The reliability of a **parallel** (active-redundant) system,
/// `R = 1 - prod(1 - R_i)`.
///
/// The system survives as long as at least one component survives, so the
/// result is bounded below by the strongest component: `R >= max(R_i)`.
///
/// # Errors
///
/// Returns [`ReliabilityError::EmptySystem`] if `components` is empty, or
/// [`ReliabilityError::ProbabilityOutOfRange`] / [`ReliabilityError::NotFinite`]
/// if any supplied reliability is not a finite probability in `[0, 1]`.
pub fn parallel(components: &[f64]) -> Result<f64, ReliabilityError> {
    if components.is_empty() {
        return Err(ReliabilityError::EmptySystem { kind: "parallel" });
    }
    let mut product_of_failures = 1.0;
    for &r in components {
        product_of_failures *= 1.0 - require_probability(r)?;
    }
    Ok(1.0 - product_of_failures)
}

/// The reliability of a **`k`-out-of-`n`** system of `n` *identical*
/// components, each with reliability `r`, that works iff at least `k`
/// components work.
///
/// This is the upper tail of a binomial distribution:
///
/// `R = sum_{j=k}^{n} C(n, j) r^j (1 - r)^(n - j)`.
///
/// The `1`-out-of-`n` case equals [`parallel`] of `n` equal components
/// and the `n`-out-of-`n` case equals [`series`] of `n` equal components.
///
/// # Errors
///
/// Returns [`ReliabilityError::InvalidKofN`] unless `1 <= k <= n`, or
/// [`ReliabilityError::ProbabilityOutOfRange`] / [`ReliabilityError::NotFinite`]
/// if `r` is not a finite probability in `[0, 1]`.
pub fn k_out_of_n(k: usize, n: usize, r: f64) -> Result<f64, ReliabilityError> {
    if k == 0 || k > n {
        return Err(ReliabilityError::InvalidKofN { k, n });
    }
    let r = require_probability(r)?;
    let q = 1.0 - r;
    let mut total = 0.0;
    for j in k..=n {
        total += binomial_coefficient(n, j) * r.powi(j as i32) * q.powi((n - j) as i32);
    }
    Ok(total)
}

/// The binomial coefficient `C(n, k)` ("n choose k") as an `f64`.
///
/// Computed multiplicatively to avoid forming intermediate factorials,
/// which keeps it exact for the small `n` typical of reliability block
/// diagrams. Returns `0.0` when `k > n`.
#[must_use]
pub fn binomial_coefficient(n: usize, k: usize) -> f64 {
    if k > n {
        return 0.0;
    }
    // C(n, k) == C(n, n - k); pick the smaller to minimise iterations.
    let k = k.min(n - k);
    let mut result = 1.0;
    for i in 0..k {
        result *= (n - i) as f64;
        result /= (i + 1) as f64;
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tolerance for floating-point ground-truth comparisons.
    const EPS: f64 = 1e-12;

    fn close(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    #[test]
    fn series_is_product() {
        // GROUND TRUTH: R = prod(R_i).
        let r = series(&[0.9, 0.8, 0.95]).unwrap();
        assert!(close(r, 0.9 * 0.8 * 0.95, EPS), "series = {r}");
    }

    #[test]
    fn series_single_component_is_itself() {
        let r = series(&[0.73]).unwrap();
        assert!(close(r, 0.73, EPS), "series([0.73]) = {r}");
    }

    #[test]
    fn series_is_at_most_the_weakest_component() {
        // GROUND TRUTH: R_series <= min(R_i).
        let comps = [0.99, 0.6, 0.85, 0.7];
        let r = series(&comps).unwrap();
        let min = comps.iter().cloned().fold(f64::INFINITY, f64::min);
        assert!(r <= min + EPS, "series {r} must be <= min {min}");
        // With at least two non-unit components it is strictly less.
        assert!(r < min, "series {r} should be strictly below min {min}");
    }

    #[test]
    fn series_with_a_perfect_component_unchanged() {
        // A component with R = 1 contributes a factor of 1.
        let with = series(&[1.0, 0.8, 0.9]).unwrap();
        let without = series(&[0.8, 0.9]).unwrap();
        assert!(close(with, without, EPS), "{with} vs {without}");
    }

    #[test]
    fn parallel_is_one_minus_product_of_failures() {
        // GROUND TRUTH: R = 1 - prod(1 - R_i).
        let r = parallel(&[0.9, 0.8]).unwrap();
        let expected = 1.0 - (1.0 - 0.9) * (1.0 - 0.8);
        assert!(close(r, expected, EPS), "parallel = {r} expected {expected}");
    }

    #[test]
    fn parallel_single_component_is_itself() {
        let r = parallel(&[0.42]).unwrap();
        assert!(close(r, 0.42, EPS), "parallel([0.42]) = {r}");
    }

    #[test]
    fn parallel_is_at_least_the_strongest_component() {
        // GROUND TRUTH: R_parallel >= max(R_i).
        let comps = [0.5, 0.7, 0.6];
        let r = parallel(&comps).unwrap();
        let max = comps.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(r >= max - EPS, "parallel {r} must be >= max {max}");
        assert!(r > max, "parallel {r} should be strictly above max {max}");
    }

    #[test]
    fn redundancy_improves_reliability() {
        // Two identical 0.9 units in parallel beat a single one.
        let single = 0.9;
        let dual = parallel(&[0.9, 0.9]).unwrap();
        // GROUND TRUTH: 1 - 0.1^2 = 0.99.
        assert!(close(dual, 0.99, EPS), "dual = {dual}");
        assert!(dual > single);
    }

    #[test]
    fn series_and_parallel_bracket_each_other() {
        // For the same components, series <= any single <= parallel.
        let comps = [0.8, 0.85, 0.9];
        let s = series(&comps).unwrap();
        let p = parallel(&comps).unwrap();
        assert!(s < p, "series {s} should be below parallel {p}");
        for &c in &comps {
            assert!(s <= c + EPS && c <= p + EPS, "component {c} not bracketed");
        }
    }

    #[test]
    fn empty_systems_are_rejected() {
        assert!(matches!(
            series(&[]),
            Err(ReliabilityError::EmptySystem { kind: "series" })
        ));
        assert!(matches!(
            parallel(&[]),
            Err(ReliabilityError::EmptySystem { kind: "parallel" })
        ));
    }

    #[test]
    fn out_of_range_component_is_rejected() {
        assert!(matches!(
            series(&[0.5, 1.2]),
            Err(ReliabilityError::ProbabilityOutOfRange { .. })
        ));
        assert!(matches!(
            parallel(&[-0.1, 0.5]),
            Err(ReliabilityError::ProbabilityOutOfRange { .. })
        ));
        assert!(matches!(
            series(&[f64::NAN]),
            Err(ReliabilityError::NotFinite { .. })
        ));
    }

    #[test]
    fn binomial_coefficients_are_correct() {
        // GROUND TRUTH: Pascal's-triangle values.
        assert!(close(binomial_coefficient(0, 0), 1.0, EPS));
        assert!(close(binomial_coefficient(5, 0), 1.0, EPS));
        assert!(close(binomial_coefficient(5, 5), 1.0, EPS));
        assert!(close(binomial_coefficient(5, 2), 10.0, EPS));
        assert!(close(binomial_coefficient(6, 3), 20.0, EPS));
        assert!(close(binomial_coefficient(10, 4), 210.0, EPS));
        // k > n yields 0.
        assert!(close(binomial_coefficient(3, 4), 0.0, EPS));
    }

    #[test]
    fn k_of_n_n_of_n_equals_series_of_equal_components() {
        // GROUND TRUTH: n-out-of-n is a pure series of equal parts.
        let r = 0.92;
        let n = 4;
        let k_of_n = k_out_of_n(n, n, r).unwrap();
        let series_eq = series(&[r; 4]).unwrap();
        assert!(close(k_of_n, series_eq, EPS), "{k_of_n} vs {series_eq}");
        // And that equals r^n directly.
        assert!(close(k_of_n, r.powi(4), EPS));
    }

    #[test]
    fn k_of_n_1_of_n_equals_parallel_of_equal_components() {
        // GROUND TRUTH: 1-out-of-n is active parallel of equal parts.
        let r = 0.8;
        let n = 3;
        let one_of_n = k_out_of_n(1, n, r).unwrap();
        let parallel_eq = parallel(&[r; 3]).unwrap();
        assert!(close(one_of_n, parallel_eq, EPS), "{one_of_n} vs {parallel_eq}");
        // And that equals 1 - (1 - r)^n.
        assert!(close(one_of_n, 1.0 - (1.0 - r).powi(3), EPS));
    }

    #[test]
    fn k_of_n_two_of_three_closed_form() {
        // GROUND TRUTH (2-out-of-3 voting):
        //   R = 3 r^2 (1 - r) + r^3 = r^2 (3 - 2 r).
        let r = 0.9;
        let got = k_out_of_n(2, 3, r).unwrap();
        let expected = r * r * (3.0 - 2.0 * r);
        assert!(close(got, expected, EPS), "2oo3 = {got} expected {expected}");
    }

    #[test]
    fn k_of_n_is_monotone_decreasing_in_k() {
        // More required components -> lower reliability (for r < 1).
        let r = 0.85;
        let n = 5;
        let mut prev = f64::INFINITY;
        for k in 1..=n {
            let val = k_out_of_n(k, n, r).unwrap();
            assert!(val < prev, "k_out_of_n must decrease in k: k = {k}, val = {val}");
            prev = val;
        }
    }

    #[test]
    fn k_of_n_validates_k_range() {
        assert!(matches!(
            k_out_of_n(0, 3, 0.9),
            Err(ReliabilityError::InvalidKofN { k: 0, n: 3 })
        ));
        assert!(matches!(
            k_out_of_n(4, 3, 0.9),
            Err(ReliabilityError::InvalidKofN { k: 4, n: 3 })
        ));
        assert!(matches!(
            k_out_of_n(2, 3, 1.5),
            Err(ReliabilityError::ProbabilityOutOfRange { .. })
        ));
    }

    #[test]
    fn perfect_components_give_perfect_systems() {
        // R_i = 1 everywhere -> system reliability 1 in all topologies.
        assert!(close(series(&[1.0, 1.0, 1.0]).unwrap(), 1.0, EPS));
        assert!(close(parallel(&[1.0, 1.0]).unwrap(), 1.0, EPS));
        assert!(close(k_out_of_n(2, 3, 1.0).unwrap(), 1.0, EPS));
    }

    #[test]
    fn certain_failure_components() {
        // A series with a dead component is dead; a parallel of dead
        // components is dead.
        assert!(close(series(&[0.9, 0.0]).unwrap(), 0.0, EPS));
        assert!(close(parallel(&[0.0, 0.0]).unwrap(), 0.0, EPS));
    }
}
