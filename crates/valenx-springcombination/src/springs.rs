//! Linear (Hookean) springs and their parallel / series combinations.
//!
//! ## Model
//!
//! A linear spring is characterised by a single positive *rate* (or
//! *stiffness*) `k`, in newtons per metre (N/m). For a deflection `x`
//! (metres) from its free length it produces a restoring force and stores
//! elastic potential energy:
//!
//! ```text
//! F = k * x                 (Hooke's law,  N)
//! U = 0.5 * k * x^2         (stored energy, J)
//! ```
//!
//! Combining springs replaces several rates with a single *equivalent*
//! rate that behaves identically for the whole assembly:
//!
//! ```text
//! parallel:  k_eq = sum(k_i)                 (rates add — stiffer)
//! series:    1 / k_eq = sum(1 / k_i)         (compliances add — softer)
//! ```
//!
//! The parallel rule follows from every spring sharing the same
//! displacement and the forces adding; the series rule follows from every
//! spring carrying the same force and the displacements adding. For two
//! equal springs of rate `k` these give `2k` (parallel) and `k/2`
//! (series) respectively — the canonical sanity check.
//!
//! ## Honest scope
//!
//! These are the ideal-spring formulas from a first mechanics course.
//! They assume a single linear regime with no preload, no end effects, no
//! buckling, no fatigue, no hysteresis and no manufacturing tolerance,
//! and they say nothing about the geometry or material of a *physical*
//! coil. This is a research / educational tool, **not** a clinical,
//! medical, or production engineering one.

use crate::error::SpringError;
use serde::{Deserialize, Serialize};

/// A single ideal linear spring, defined solely by its rate.
///
/// Construct one with [`Spring::new`], which validates that the rate is
/// finite and strictly positive. The stored [`rate`](Spring::rate) is
/// then always a usable N/m value.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Spring {
    /// Spring rate (stiffness) `k`, in newtons per metre (N/m).
    rate: f64,
}

impl Spring {
    /// Create a spring of the given rate `k` (N/m).
    ///
    /// # Errors
    ///
    /// Returns [`SpringError::NonPositiveRate`] if `k` is not finite or is
    /// not strictly greater than zero.
    ///
    /// ```
    /// use valenx_springcombination::Spring;
    /// let s = Spring::new(100.0).unwrap();
    /// assert!(Spring::new(0.0).is_err());
    /// assert!(Spring::new(-5.0).is_err());
    /// ```
    pub fn new(rate: f64) -> Result<Self, SpringError> {
        validate_rate(rate, "Spring::new")?;
        Ok(Self { rate })
    }

    /// The spring rate `k`, in newtons per metre (N/m).
    ///
    /// Always finite and strictly positive by construction.
    pub fn rate(&self) -> f64 {
        self.rate
    }

    /// The restoring force `F = k * x` at deflection `x` (metres).
    ///
    /// The result is in newtons (N). The sign of the force follows the
    /// sign of `x` (a positive deflection yields a positive force on this
    /// convention); the caller chooses whether "positive" means
    /// compression or extension.
    ///
    /// # Errors
    ///
    /// Returns [`SpringError::NonFiniteDisplacement`] if `x` is NaN or
    /// infinite.
    ///
    /// ```
    /// use valenx_springcombination::Spring;
    /// let s = Spring::new(200.0).unwrap();
    /// // 200 N/m * 0.05 m = 10 N
    /// assert!((s.force(0.05).unwrap() - 10.0).abs() < 1e-12);
    /// ```
    pub fn force(&self, x: f64) -> Result<f64, SpringError> {
        validate_displacement(x)?;
        Ok(self.rate * x)
    }

    /// The stored elastic potential energy `U = 0.5 * k * x^2` at
    /// deflection `x` (metres).
    ///
    /// The result is in joules (J) and is always non-negative because the
    /// deflection is squared.
    ///
    /// # Errors
    ///
    /// Returns [`SpringError::NonFiniteDisplacement`] if `x` is NaN or
    /// infinite.
    ///
    /// ```
    /// use valenx_springcombination::Spring;
    /// let s = Spring::new(200.0).unwrap();
    /// // 0.5 * 200 N/m * (0.05 m)^2 = 0.25 J
    /// assert!((s.energy(0.05).unwrap() - 0.25).abs() < 1e-12);
    /// ```
    pub fn energy(&self, x: f64) -> Result<f64, SpringError> {
        validate_displacement(x)?;
        Ok(0.5 * self.rate * x * x)
    }
}

/// How a set of springs is arranged.
///
/// This is the discriminant consumed by [`combine`]; the two arms map to
/// the [`parallel_rate`] and [`series_rate`] reductions.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub enum Combination {
    /// Springs share the same displacement and their forces add; the
    /// equivalent rate is the *sum* of the member rates (stiffer than any
    /// single member).
    Parallel,
    /// Springs carry the same force and their displacements add; the
    /// equivalent *compliance* (reciprocal rate) is the sum of the member
    /// compliances, so the assembly is *softer* than any single member.
    Series,
}

/// The equivalent rate of springs combined in parallel: `k = sum(k_i)`.
///
/// Because each individual rate is positive, the parallel rate is greater
/// than or equal to the largest member rate — adding a spring in parallel
/// always stiffens the assembly.
///
/// # Errors
///
/// Returns [`SpringError::EmptyCombination`] if `springs` is empty.
///
/// ```
/// use valenx_springcombination::{parallel_rate, Spring};
/// let a = Spring::new(100.0).unwrap();
/// let b = Spring::new(150.0).unwrap();
/// // 100 + 150 = 250 N/m
/// assert!((parallel_rate(&[a, b]).unwrap() - 250.0).abs() < 1e-12);
/// ```
pub fn parallel_rate(springs: &[Spring]) -> Result<f64, SpringError> {
    if springs.is_empty() {
        return Err(SpringError::EmptyCombination {
            combination: "parallel",
        });
    }
    Ok(springs.iter().map(Spring::rate).sum())
}

/// The equivalent rate of springs combined in series:
/// `1 / k = sum(1 / k_i)`.
///
/// Because each individual rate is positive, the series rate is less than
/// or equal to the smallest member rate — adding a spring in series always
/// softens the assembly.
///
/// # Errors
///
/// Returns [`SpringError::EmptyCombination`] if `springs` is empty. (Each
/// member rate is already guaranteed positive by [`Spring::new`], so the
/// reciprocal sum can neither divide by zero nor be empty-of-finite-terms
/// here.)
///
/// ```
/// use valenx_springcombination::{series_rate, Spring};
/// let a = Spring::new(100.0).unwrap();
/// let b = Spring::new(100.0).unwrap();
/// // two equal springs in series -> k/2 = 50 N/m
/// assert!((series_rate(&[a, b]).unwrap() - 50.0).abs() < 1e-12);
/// ```
pub fn series_rate(springs: &[Spring]) -> Result<f64, SpringError> {
    if springs.is_empty() {
        return Err(SpringError::EmptyCombination {
            combination: "series",
        });
    }
    let reciprocal_sum: f64 = springs.iter().map(|s| 1.0 / s.rate()).sum();
    Ok(1.0 / reciprocal_sum)
}

/// The equivalent rate for a set of springs in the given [`Combination`].
///
/// A thin dispatcher over [`parallel_rate`] and [`series_rate`]; see those
/// for the per-mode formulas and error conditions.
///
/// # Errors
///
/// Propagates [`SpringError::EmptyCombination`] from the underlying
/// reduction when `springs` is empty.
///
/// ```
/// use valenx_springcombination::{combine, Combination, Spring};
/// let s = [Spring::new(100.0).unwrap(), Spring::new(100.0).unwrap()];
/// assert!((combine(Combination::Parallel, &s).unwrap() - 200.0).abs() < 1e-12);
/// assert!((combine(Combination::Series, &s).unwrap() - 50.0).abs() < 1e-12);
/// ```
pub fn combine(mode: Combination, springs: &[Spring]) -> Result<f64, SpringError> {
    match mode {
        Combination::Parallel => parallel_rate(springs),
        Combination::Series => series_rate(springs),
    }
}

/// The restoring force `F = k * x` for a bare rate `k` (N/m) at deflection
/// `x` (metres), returning newtons.
///
/// This is the free-function form of [`Spring::force`], convenient when
/// you already hold an equivalent rate from [`combine`] rather than a
/// [`Spring`].
///
/// # Errors
///
/// Returns [`SpringError::NonPositiveRate`] if `k` is not finite and
/// positive, or [`SpringError::NonFiniteDisplacement`] if `x` is not
/// finite.
///
/// ```
/// use valenx_springcombination::force_from_rate;
/// // 250 N/m * 0.02 m = 5 N
/// assert!((force_from_rate(250.0, 0.02).unwrap() - 5.0).abs() < 1e-12);
/// ```
pub fn force_from_rate(k: f64, x: f64) -> Result<f64, SpringError> {
    validate_rate(k, "force_from_rate")?;
    validate_displacement(x)?;
    Ok(k * x)
}

/// The stored energy `U = 0.5 * k * x^2` for a bare rate `k` (N/m) at
/// deflection `x` (metres), returning joules.
///
/// This is the free-function form of [`Spring::energy`], convenient when
/// you already hold an equivalent rate from [`combine`] rather than a
/// [`Spring`].
///
/// # Errors
///
/// Returns [`SpringError::NonPositiveRate`] if `k` is not finite and
/// positive, or [`SpringError::NonFiniteDisplacement`] if `x` is not
/// finite.
///
/// ```
/// use valenx_springcombination::energy_from_rate;
/// // 0.5 * 250 N/m * (0.02 m)^2 = 0.05 J
/// assert!((energy_from_rate(250.0, 0.02).unwrap() - 0.05).abs() < 1e-12);
/// ```
pub fn energy_from_rate(k: f64, x: f64) -> Result<f64, SpringError> {
    validate_rate(k, "energy_from_rate")?;
    validate_displacement(x)?;
    Ok(0.5 * k * x * x)
}

/// Reject a rate that is not finite and strictly positive.
fn validate_rate(rate: f64, context: &'static str) -> Result<(), SpringError> {
    if !rate.is_finite() || rate <= 0.0 {
        return Err(SpringError::NonPositiveRate {
            value: rate,
            context,
        });
    }
    Ok(())
}

/// Reject a displacement that is not finite.
fn validate_displacement(x: f64) -> Result<(), SpringError> {
    if !x.is_finite() {
        return Err(SpringError::NonFiniteDisplacement { value: x });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::ErrorCategory;

    /// Tight tolerance for analytic float comparisons.
    const EPS: f64 = 1e-12;

    // --- construction / validation -------------------------------------

    #[test]
    fn new_accepts_positive_rate() {
        let s = Spring::new(42.5).unwrap();
        assert!((s.rate() - 42.5).abs() < EPS);
    }

    #[test]
    fn new_rejects_zero_negative_and_non_finite() {
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let err = Spring::new(bad).unwrap_err();
            assert_eq!(err.code(), "spring.non-positive-rate");
            assert_eq!(err.category(), ErrorCategory::Input);
        }
    }

    // --- Hooke's law F = k x -------------------------------------------

    #[test]
    fn force_is_k_times_x() {
        // 200 N/m * 0.05 m = 10 N (worked by hand).
        let s = Spring::new(200.0).unwrap();
        assert!((s.force(0.05).unwrap() - 10.0).abs() < EPS);
    }

    #[test]
    fn force_is_linear_and_odd_in_x() {
        let s = Spring::new(73.0).unwrap();
        // Linearity: F(2x) = 2 F(x).
        let f1 = s.force(0.011).unwrap();
        let f2 = s.force(0.022).unwrap();
        assert!((f2 - 2.0 * f1).abs() < EPS);
        // Oddness: F(-x) = -F(x); zero deflection -> zero force.
        assert!((s.force(-0.011).unwrap() + f1).abs() < EPS);
        assert!(s.force(0.0).unwrap().abs() < EPS);
    }

    #[test]
    fn force_rejects_non_finite_displacement() {
        let s = Spring::new(10.0).unwrap();
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let err = s.force(bad).unwrap_err();
            assert_eq!(err.code(), "spring.non-finite-displacement");
        }
    }

    // --- energy U = 1/2 k x^2 ------------------------------------------

    #[test]
    fn energy_is_half_k_x_squared() {
        // 0.5 * 200 N/m * (0.05 m)^2 = 0.25 J (worked by hand).
        let s = Spring::new(200.0).unwrap();
        assert!((s.energy(0.05).unwrap() - 0.25).abs() < EPS);
    }

    #[test]
    fn energy_is_even_in_x_and_nonnegative() {
        let s = Spring::new(123.4).unwrap();
        let up = s.energy(0.03).unwrap();
        let down = s.energy(-0.03).unwrap();
        // Even function: U(x) == U(-x).
        assert!((up - down).abs() < EPS);
        // Non-negative, and zero only at zero deflection.
        assert!(up > 0.0);
        assert!(s.energy(0.0).unwrap().abs() < EPS);
    }

    #[test]
    fn energy_quadruples_when_displacement_doubles() {
        let s = Spring::new(500.0).unwrap();
        let u1 = s.energy(0.004).unwrap();
        let u2 = s.energy(0.008).unwrap();
        assert!((u2 - 4.0 * u1).abs() < EPS);
    }

    #[test]
    fn energy_equals_half_force_times_x() {
        // Work-energy consistency: U = 1/2 F x for a linear spring.
        let s = Spring::new(880.0).unwrap();
        let x = 0.017;
        let u = s.energy(x).unwrap();
        let f = s.force(x).unwrap();
        assert!((u - 0.5 * f * x).abs() < EPS);
    }

    // --- parallel: rates add -------------------------------------------

    #[test]
    fn parallel_adds_rates() {
        // 100 + 150 + 250 = 500 N/m (worked by hand).
        let springs = [
            Spring::new(100.0).unwrap(),
            Spring::new(150.0).unwrap(),
            Spring::new(250.0).unwrap(),
        ];
        assert!((parallel_rate(&springs).unwrap() - 500.0).abs() < EPS);
    }

    #[test]
    fn two_equal_in_parallel_is_2k() {
        let k = 320.0;
        let springs = [Spring::new(k).unwrap(), Spring::new(k).unwrap()];
        assert!((parallel_rate(&springs).unwrap() - 2.0 * k).abs() < EPS);
    }

    #[test]
    fn parallel_single_spring_is_identity() {
        let s = [Spring::new(77.0).unwrap()];
        assert!((parallel_rate(&s).unwrap() - 77.0).abs() < EPS);
    }

    #[test]
    fn parallel_is_at_least_the_max_member() {
        // Adding a spring in parallel can only stiffen the assembly.
        let springs = [Spring::new(10.0).unwrap(), Spring::new(90.0).unwrap()];
        let k_eq = parallel_rate(&springs).unwrap();
        assert!(k_eq >= 90.0 - EPS);
        assert!((k_eq - 100.0).abs() < EPS);
    }

    #[test]
    fn parallel_rejects_empty() {
        let err = parallel_rate(&[]).unwrap_err();
        assert_eq!(err.code(), "spring.empty-combination");
        assert_eq!(err.category(), ErrorCategory::Domain);
    }

    // --- series: softer -------------------------------------------------

    #[test]
    fn two_equal_in_series_is_k_over_2() {
        let k = 320.0;
        let springs = [Spring::new(k).unwrap(), Spring::new(k).unwrap()];
        assert!((series_rate(&springs).unwrap() - k / 2.0).abs() < EPS);
    }

    #[test]
    fn series_three_equal_is_k_over_3() {
        let k = 90.0;
        let springs = [
            Spring::new(k).unwrap(),
            Spring::new(k).unwrap(),
            Spring::new(k).unwrap(),
        ];
        assert!((series_rate(&springs).unwrap() - k / 3.0).abs() < EPS);
    }

    #[test]
    fn series_unequal_matches_reciprocal_formula() {
        // 1/k = 1/200 + 1/300 = 5/600  ->  k = 120 N/m (worked by hand).
        let springs = [Spring::new(200.0).unwrap(), Spring::new(300.0).unwrap()];
        assert!((series_rate(&springs).unwrap() - 120.0).abs() < EPS);
    }

    #[test]
    fn series_single_spring_is_identity() {
        let s = [Spring::new(77.0).unwrap()];
        assert!((series_rate(&s).unwrap() - 77.0).abs() < EPS);
    }

    #[test]
    fn series_is_at_most_the_min_member() {
        // Adding a spring in series can only soften the assembly.
        let springs = [Spring::new(50.0).unwrap(), Spring::new(200.0).unwrap()];
        let k_eq = series_rate(&springs).unwrap();
        assert!(k_eq <= 50.0 + EPS);
        // 1/40 = 1/50 + 1/200 -> 40 N/m.
        assert!((k_eq - 40.0).abs() < EPS);
    }

    #[test]
    fn series_rejects_empty() {
        let err = series_rate(&[]).unwrap_err();
        assert_eq!(err.code(), "spring.empty-combination");
    }

    // --- series is softer than parallel for the same set --------------

    #[test]
    fn series_is_softer_than_parallel() {
        let springs = [
            Spring::new(120.0).unwrap(),
            Spring::new(240.0).unwrap(),
            Spring::new(360.0).unwrap(),
        ];
        let kp = parallel_rate(&springs).unwrap();
        let ks = series_rate(&springs).unwrap();
        assert!(ks < kp);
        // For two equal springs the ratio is exactly 4 (2k vs k/2); for a
        // general positive set series is strictly the softest, parallel
        // the stiffest, of all wirings.
        for s in &springs {
            assert!(ks <= s.rate() + EPS);
            assert!(kp >= s.rate() - EPS);
        }
    }

    // --- combine dispatcher --------------------------------------------

    #[test]
    fn combine_matches_dedicated_reductions() {
        let springs = [Spring::new(100.0).unwrap(), Spring::new(100.0).unwrap()];
        assert!((combine(Combination::Parallel, &springs).unwrap() - 200.0).abs() < EPS);
        assert!((combine(Combination::Series, &springs).unwrap() - 50.0).abs() < EPS);
    }

    // --- free-function force / energy ----------------------------------

    #[test]
    fn free_force_and_energy_agree_with_methods() {
        let k = 654.0;
        let x = 0.013;
        let s = Spring::new(k).unwrap();
        assert!((force_from_rate(k, x).unwrap() - s.force(x).unwrap()).abs() < EPS);
        assert!((energy_from_rate(k, x).unwrap() - s.energy(x).unwrap()).abs() < EPS);
    }

    #[test]
    fn free_functions_validate_inputs() {
        assert_eq!(
            force_from_rate(-1.0, 0.0).unwrap_err().code(),
            "spring.non-positive-rate"
        );
        assert_eq!(
            energy_from_rate(0.0, 0.0).unwrap_err().code(),
            "spring.non-positive-rate"
        );
        assert_eq!(
            force_from_rate(10.0, f64::NAN).unwrap_err().code(),
            "spring.non-finite-displacement"
        );
        assert_eq!(
            energy_from_rate(10.0, f64::INFINITY).unwrap_err().code(),
            "spring.non-finite-displacement"
        );
    }

    // --- end-to-end physics: equivalent rate carries the dynamics ------

    #[test]
    fn equivalent_rate_reproduces_assembly_force() {
        // Two springs in parallel share displacement; total force is the
        // sum of member forces and equals k_eq * x.
        let a = Spring::new(150.0).unwrap();
        let b = Spring::new(250.0).unwrap();
        let x = 0.04;
        let summed = a.force(x).unwrap() + b.force(x).unwrap();
        let k_eq = parallel_rate(&[a, b]).unwrap();
        let via_eq = force_from_rate(k_eq, x).unwrap();
        assert!((summed - via_eq).abs() < EPS);
        // 150*0.04 + 250*0.04 = 16 N.
        assert!((summed - 16.0).abs() < EPS);
    }

    #[test]
    fn series_equivalent_energy_matches_member_sum_under_shared_force() {
        // Springs in series carry the same force F; each stretches by
        // x_i = F / k_i and stores 0.5 * F * x_i. Total stored energy
        // equals 0.5 * F * x_total = 0.5 * F^2 / k_eq.
        let a = Spring::new(200.0).unwrap();
        let b = Spring::new(300.0).unwrap();
        let force = 6.0; // N, common to both.
        let xa = force / a.rate();
        let xb = force / b.rate();
        let member_energy = a.energy(xa).unwrap() + b.energy(xb).unwrap();
        let k_eq = series_rate(&[a, b]).unwrap();
        let assembly_energy = energy_from_rate(k_eq, xa + xb).unwrap();
        assert!((member_energy - assembly_energy).abs() < EPS);
    }
}
