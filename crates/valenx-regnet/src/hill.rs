//! Hill kinetics — the sigmoidal building blocks of regulatory inputs.
//!
//! Transcriptional regulation is modelled with the classic **Hill
//! functions**, which describe the fractional saturation of a promoter by
//! a regulator at concentration `x`, given a half-saturation threshold `k`
//! and a cooperativity (Hill) coefficient `n`:
//!
//! - **Activation** (an activator turns a gene *on* as it accumulates):
//!   `hill_activate(x, k, n) = x^n / (k^n + x^n)`. It rises monotonically
//!   from `0` at `x = 0` toward `1` as `x → ∞`, passing through `0.5` at
//!   `x = k`.
//! - **Repression** (a repressor turns a gene *off* as it accumulates):
//!   `hill_repress(x, k, n) = k^n / (k^n + x^n) = 1 - hill_activate(...)`.
//!   It falls monotonically from `1` at `x = 0` toward `0` as `x → ∞`,
//!   also crossing `0.5` at `x = k`.
//!
//! `n = 1` gives the (non-cooperative) Michaelis-Menten shape; `n > 1`
//! sharpens the switch toward a step at `x = k`.
//!
//! These are the textbook forms used throughout systems biology (e.g.
//! Alon, *An Introduction to Systems Biology*). Inputs are validated by
//! [`hill_activate_checked`] / [`hill_repress_checked`]; the bare
//! [`hill_activate`] / [`hill_repress`] assume already-valid `k > 0`,
//! `n > 0` and a finite, non-negative `x`.

use crate::error::{RegnetError, Result};

/// Validate a `(k, n)` Hill parameter pair: both must be strictly
/// positive and finite.
fn check_kn(k: f64, n: f64) -> Result<()> {
    if !(k.is_finite() && k > 0.0) {
        return Err(RegnetError::InvalidHill {
            what: "threshold k",
            value: k,
        });
    }
    if !(n.is_finite() && n > 0.0) {
        return Err(RegnetError::InvalidHill {
            what: "coefficient n",
            value: n,
        });
    }
    Ok(())
}

/// Hill **activation** function `x^n / (k^n + x^n)`.
///
/// Returns the fraction (in `[0, 1]`) of promoters bound by an activator
/// present at concentration `x`, with half-saturation threshold `k` and
/// Hill coefficient `n`. Strictly increasing in `x`; equals `0.5` exactly
/// at `x == k`.
///
/// A negative `x` is clamped to `0` (concentrations cannot be negative).
/// Assumes `k > 0` and `n > 0`; use [`hill_activate_checked`] to validate
/// those at a boundary where they may be user-supplied.
///
/// The computation is arranged as `r^n / (1 + r^n)` with `r = x / k` so
/// that it stays finite (and tends to `1`) even when `x` is enormous.
#[must_use]
pub fn hill_activate(x: f64, k: f64, n: f64) -> f64 {
    let x = if x < 0.0 { 0.0 } else { x };
    if x == 0.0 {
        return 0.0;
    }
    // r = x / k > 0; rn = (x/k)^n. As x -> inf, rn -> inf and the ratio
    // rn / (1 + rn) -> 1 without overflowing to NaN.
    let rn = (x / k).powf(n);
    if rn.is_infinite() {
        1.0
    } else {
        rn / (1.0 + rn)
    }
}

/// Hill **repression** function `k^n / (k^n + x^n)`.
///
/// Returns the fraction (in `[0, 1]`) of promoters left *free* of a
/// repressor present at concentration `x`, with half-saturation threshold
/// `k` and Hill coefficient `n`. This is the de-repressed (still-active)
/// fraction: strictly decreasing in `x`, equal to `0.5` exactly at
/// `x == k`, and identically `1 - hill_activate(x, k, n)`.
///
/// A negative `x` is clamped to `0`. Assumes `k > 0` and `n > 0`; use
/// [`hill_repress_checked`] to validate those at an untrusted boundary.
#[must_use]
pub fn hill_repress(x: f64, k: f64, n: f64) -> f64 {
    let x = if x < 0.0 { 0.0 } else { x };
    if x == 0.0 {
        return 1.0;
    }
    let rn = (x / k).powf(n);
    if rn.is_infinite() {
        0.0
    } else {
        1.0 / (1.0 + rn)
    }
}

/// Validated [`hill_activate`]: errors with [`RegnetError::InvalidHill`]
/// if `k <= 0` or `n <= 0` (or either is non-finite).
///
/// # Errors
///
/// Returns [`RegnetError::InvalidHill`] when the threshold or coefficient
/// is not strictly positive and finite.
pub fn hill_activate_checked(x: f64, k: f64, n: f64) -> Result<f64> {
    check_kn(k, n)?;
    Ok(hill_activate(x, k, n))
}

/// Validated [`hill_repress`]: errors with [`RegnetError::InvalidHill`]
/// if `k <= 0` or `n <= 0` (or either is non-finite).
///
/// # Errors
///
/// Returns [`RegnetError::InvalidHill`] when the threshold or coefficient
/// is not strictly positive and finite.
pub fn hill_repress_checked(x: f64, k: f64, n: f64) -> Result<f64> {
    check_kn(k, n)?;
    Ok(hill_repress(x, k, n))
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    #[test]
    fn activate_is_half_at_threshold() {
        // hill_activate(k, k, n) = k^n / (k^n + k^n) = 1/2 for any n.
        for &n in &[1.0, 2.0, 4.0, 0.5] {
            let v = hill_activate(2.5, 2.5, n);
            assert!((v - 0.5).abs() < EPS, "n={n}: got {v}");
        }
    }

    #[test]
    fn repress_is_half_at_threshold() {
        for &n in &[1.0, 2.0, 4.0, 0.5] {
            let v = hill_repress(2.5, 2.5, n);
            assert!((v - 0.5).abs() < EPS, "n={n}: got {v}");
        }
    }

    #[test]
    fn activate_strictly_increasing_in_x() {
        let (k, n) = (1.0, 2.0);
        let mut prev = hill_activate(0.0, k, n);
        assert!((prev - 0.0).abs() < EPS, "at x=0 should be 0, got {prev}");
        let mut x = 0.05;
        while x <= 10.0 {
            let cur = hill_activate(x, k, n);
            assert!(cur > prev, "not increasing at x={x}: {cur} <= {prev}");
            prev = cur;
            x += 0.05;
        }
        // Approaches 1 from below for large x.
        let big = hill_activate(1.0e6, k, n);
        assert!(big > 0.999_999 && big <= 1.0, "tail not -> 1: {big}");
    }

    #[test]
    fn repress_strictly_decreasing_in_x() {
        let (k, n) = (1.0, 2.0);
        let mut prev = hill_repress(0.0, k, n);
        assert!((prev - 1.0).abs() < EPS, "at x=0 should be 1, got {prev}");
        let mut x = 0.05;
        while x <= 10.0 {
            let cur = hill_repress(x, k, n);
            assert!(cur < prev, "not decreasing at x={x}: {cur} >= {prev}");
            prev = cur;
            x += 0.05;
        }
        let big = hill_repress(1.0e6, k, n);
        assert!((0.0..1.0e-6).contains(&big), "tail not -> 0: {big}");
    }

    #[test]
    fn activate_and_repress_are_complementary() {
        // hill_repress = 1 - hill_activate for any (x, k, n).
        for &x in &[0.0, 0.3, 1.0, 2.7, 50.0] {
            let a = hill_activate(x, 1.4, 3.0);
            let r = hill_repress(x, 1.4, 3.0);
            assert!((a + r - 1.0).abs() < EPS, "x={x}: a+r={}", a + r);
        }
    }

    #[test]
    fn activate_matches_closed_form_michaelis_menten() {
        // n = 1: x / (k + x).
        let (x, k) = (3.0, 2.0);
        let expected = x / (k + x);
        let got = hill_activate(x, k, 1.0);
        assert!((got - expected).abs() < EPS, "got {got}, want {expected}");
    }

    #[test]
    fn activate_matches_closed_form_n2() {
        // n = 2: x^2 / (k^2 + x^2).
        let (x, k) = (3.0, 2.0);
        let expected = x * x / (k * k + x * x);
        let got = hill_activate(x, k, 2.0);
        assert!((got - expected).abs() < EPS, "got {got}, want {expected}");
    }

    #[test]
    fn negative_x_is_clamped_to_zero() {
        assert!((hill_activate(-5.0, 1.0, 2.0) - 0.0).abs() < EPS);
        assert!((hill_repress(-5.0, 1.0, 2.0) - 1.0).abs() < EPS);
    }

    #[test]
    fn checked_rejects_bad_parameters() {
        assert!(hill_activate_checked(1.0, 0.0, 2.0).is_err());
        assert!(hill_activate_checked(1.0, -1.0, 2.0).is_err());
        assert!(hill_activate_checked(1.0, 1.0, 0.0).is_err());
        assert!(hill_repress_checked(1.0, 1.0, -2.0).is_err());
        assert!(hill_activate_checked(1.0, f64::NAN, 2.0).is_err());
        // A valid pair passes through to the bare function.
        let ok = hill_activate_checked(2.0, 2.0, 2.0).unwrap();
        assert!((ok - 0.5).abs() < EPS, "got {ok}");
    }

    #[test]
    fn higher_n_is_sharper_switch() {
        // For x just below k, a larger n drives activation lower (steeper).
        let below = hill_activate(0.8, 1.0, 8.0);
        let below_soft = hill_activate(0.8, 1.0, 1.0);
        assert!(below < below_soft, "{below} !< {below_soft}");
        // For x just above k, a larger n drives activation higher.
        let above = hill_activate(1.25, 1.0, 8.0);
        let above_soft = hill_activate(1.25, 1.0, 1.0);
        assert!(above > above_soft, "{above} !> {above_soft}");
    }
}
