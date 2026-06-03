//! The Boys function `F_n(x)`.
//!
//! The Boys function
//!
//! ```text
//! F_n(x) = ∫₀¹ t^{2n} exp(-x t²) dt
//! ```
//!
//! is the special function at the heart of every Coulomb integral in
//! Gaussian-basis quantum chemistry — the nuclear-attraction integrals
//! and the two-electron repulsion integrals both reduce to a sum of
//! `F_n` values.
//!
//! ## Evaluation strategy
//!
//! Two regimes, switched at `x = 25`:
//!
//! - **Small / moderate `x`** — the rapidly convergent ascending series
//!   `F_n(x) = exp(-x) Σ_{k≥0} (2x)^k (2n-1)!! / (2n+2k+1)!!`. The
//!   highest order is summed directly and the lower orders are obtained
//!   by the exact stable downward recurrence
//!   `F_{n-1}(x) = (2x F_n(x) + exp(-x)) / (2n-1)`.
//! - **Large `x`** — the asymptotic limit `F_0(x) → ½√(π/x)` followed
//!   by the upward recurrence
//!   `F_{n+1}(x) = ((2n+1) F_n(x) - exp(-x)) / (2x)`, which is stable
//!   when `exp(-x)` is negligible.
//!
//! [`boys_array`] returns `F_0 … F_nmax` in one call — the form the
//! McMurchie-Davidson Coulomb routines consume.

/// Crossover between the series and the asymptotic evaluation.
const ASYMPTOTIC_X: f64 = 25.0;

/// Evaluate `F_0(x) … F_{n_max}(x)` and return them as a vector of
/// length `n_max + 1`.
pub fn boys_array(n_max: usize, x: f64) -> Vec<f64> {
    let mut out = vec![0.0; n_max + 1];
    if x < ASYMPTOTIC_X {
        // Sum the top order from the ascending series, then recur down.
        out[n_max] = boys_series(n_max, x);
        let ex = (-x).exp();
        for n in (0..n_max).rev() {
            // F_n = (2x F_{n+1} + e^{-x}) / (2n+1)
            out[n] = (2.0 * x * out[n + 1] + ex) / (2.0 * n as f64 + 1.0);
        }
    } else {
        // Asymptotic F_0, then recur up (e^{-x} ~ 0 so it is stable).
        out[0] = 0.5 * (std::f64::consts::PI / x).sqrt();
        let ex = (-x).exp();
        for n in 0..n_max {
            // F_{n+1} = ((2n+1) F_n - e^{-x}) / (2x)
            out[n + 1] = ((2.0 * n as f64 + 1.0) * out[n] - ex) / (2.0 * x);
        }
    }
    out
}

/// Single-order Boys function `F_n(x)`.
pub fn boys(n: usize, x: f64) -> f64 {
    boys_array(n, x)[n]
}

/// `F_n(x)` from the ascending power series. Converges for every `x`
/// but is only used below [`ASYMPTOTIC_X`], where it needs few terms.
fn boys_series(n: usize, x: f64) -> f64 {
    // F_n(x) = e^{-x} Σ_{k≥0} (2x)^k · (2n-1)!! / (2n+2k+1)!!
    //        = e^{-x} Σ_{k≥0} term_k  with
    // term_0 = 1/(2n+1), term_{k+1} = term_k · 2x / (2n+2k+3).
    let ex = (-x).exp();
    let mut term = 1.0 / (2.0 * n as f64 + 1.0);
    let mut sum = term;
    let mut k = 0usize;
    loop {
        let denom = 2.0 * n as f64 + 2.0 * k as f64 + 3.0;
        term *= 2.0 * x / denom;
        sum += term;
        k += 1;
        if term < 1.0e-17 * sum || k > 400 {
            break;
        }
    }
    ex * sum
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn boys_zero_argument() {
        // F_n(0) = 1/(2n+1).
        for n in 0..6 {
            let f = boys(n, 0.0);
            assert!(
                (f - 1.0 / (2.0 * n as f64 + 1.0)).abs() < 1.0e-12,
                "F_{n}(0) = {f}"
            );
        }
    }

    #[test]
    fn boys_zero_limit_large_x() {
        // F_0(x) → ½√(π/x) as x → ∞.
        let x = 60.0;
        let f0 = boys(0, x);
        let asym = 0.5 * (std::f64::consts::PI / x).sqrt();
        assert!((f0 - asym).abs() < 1.0e-6, "F_0({x}) = {f0}");
    }

    #[test]
    fn boys_known_value() {
        // F_0(x) = ½ √(π/x) erf(√x). At x = 1: erf(1) = 0.8427007929...
        let f0 = boys(0, 1.0);
        let expect = 0.5 * std::f64::consts::PI.sqrt() * 0.842_700_792_949_715;
        assert!((f0 - expect).abs() < 1.0e-10, "F_0(1) = {f0}");
    }

    #[test]
    fn downward_recurrence_consistent() {
        // The relation 2x F_{n+1} + e^{-x} = (2n+1) F_n must hold.
        let x = 3.7;
        let f = boys_array(5, x);
        let ex = (-x).exp();
        for n in 0..5 {
            let lhs = 2.0 * x * f[n + 1] + ex;
            let rhs = (2.0 * n as f64 + 1.0) * f[n];
            assert!((lhs - rhs).abs() < 1.0e-12, "recurrence broke at n={n}");
        }
    }

    #[test]
    fn series_and_asymptotic_agree_near_crossover() {
        // Both regimes must give close values just below / above 25.
        let lo = boys(3, 24.9);
        let hi = boys(3, 25.1);
        // Monotone decreasing and continuous-ish.
        assert!(lo > hi);
        assert!((lo - hi).abs() < 1.0e-3);
    }

    #[test]
    fn array_length_and_monotonicity() {
        let f = boys_array(4, 2.0);
        assert_eq!(f.len(), 5);
        // F_n decreases with n for fixed positive x.
        for w in f.windows(2) {
            assert!(w[0] > w[1]);
        }
    }
}
