//! Discrete-gamma rate heterogeneity (Yang 1994).
//!
//! Real sequences evolve at site-specific rates: some columns are
//! near-invariant, others hypervariable. The standard model draws each
//! site's relative rate from a Gamma distribution with mean 1 and shape
//! `α` (small `α` ⇒ strong heterogeneity, large `α` ⇒ near-uniform).
//!
//! Integrating the likelihood over a continuous Gamma is intractable,
//! so Yang's discrete-gamma approximation splits the distribution into
//! `k` equal-probability categories and represents each by a single
//! representative rate — the **mean of the Gamma density over that
//! category's quantile interval**. The per-site likelihood is then the
//! average of the likelihoods under the `k` category rates.
//!
//! This module provides [`DiscreteGamma`]: the category rates for a
//! given `α` and `k`. The category rates are normalised so their mean
//! is exactly 1, keeping branch lengths interpretable.

use crate::error::{PhyloError, Result};

/// A discrete-gamma rate-category set.
#[derive(Debug, Clone, PartialEq)]
pub struct DiscreteGamma {
    /// Shape parameter `α`.
    alpha: f64,
    /// Representative relative rate of each category (mean = 1).
    rates: Vec<f64>,
}

impl DiscreteGamma {
    /// Builds `k` equal-probability gamma rate categories for shape
    /// `alpha`.
    ///
    /// # Errors
    /// [`PhyloError::Invalid`] if `alpha <= 0` or `k == 0`.
    pub fn new(alpha: f64, k: usize) -> Result<Self> {
        if alpha <= 0.0 {
            return Err(PhyloError::invalid("alpha", "gamma shape must be > 0"));
        }
        if k == 0 {
            return Err(PhyloError::invalid("k", "need at least one category"));
        }
        if k == 1 {
            return Ok(DiscreteGamma {
                alpha,
                rates: vec![1.0],
            });
        }
        // A mean-1 Gamma has shape α and scale 1/α.
        let beta = 1.0 / alpha;
        let mut rates = Vec::with_capacity(k);
        // Category boundaries at the i/k quantiles; the representative
        // rate of a category is the conditional mean of the Gamma over
        // its quantile interval, which for the mean-1 Gamma equals
        // [I(b_{i+1}; α+1) - I(b_i; α+1)] · k, where I is the
        // regularised lower incomplete gamma. (Yang 1994, eq. 10.)
        let mut prev_bound = 0.0;
        let mut prev_partial = 0.0; // I(0; α+1) = 0
        for i in 1..=k {
            let cum = i as f64 / k as f64;
            let bound = if i == k {
                f64::INFINITY
            } else {
                gamma_quantile(cum, alpha, beta)
            };
            // Partial expectation up to `bound`.
            let partial = if bound.is_infinite() {
                1.0 // I(∞; α+1) = 1
            } else {
                regularized_lower_gamma(alpha + 1.0, bound / beta)
            };
            let cat_rate = (partial - prev_partial) * k as f64;
            rates.push(cat_rate);
            prev_bound = bound;
            prev_partial = partial;
        }
        let _ = prev_bound;
        // Normalise so the category mean is exactly 1.
        let mean: f64 = rates.iter().sum::<f64>() / k as f64;
        if mean > 0.0 {
            for r in &mut rates {
                *r /= mean;
            }
        }
        Ok(DiscreteGamma { alpha, rates })
    }

    /// The shape parameter `α`.
    pub fn alpha(&self) -> f64 {
        self.alpha
    }

    /// Number of rate categories.
    pub fn n_categories(&self) -> usize {
        self.rates.len()
    }

    /// The representative relative rates, one per category.
    pub fn rates(&self) -> &[f64] {
        &self.rates
    }

    /// Each category has equal prior probability `1/k`.
    pub fn category_probability(&self) -> f64 {
        1.0 / self.rates.len() as f64
    }
}

/// Regularised lower incomplete gamma function `P(s, x) = γ(s,x)/Γ(s)`.
///
/// Uses the series expansion for `x < s + 1` and the continued fraction
/// otherwise — the standard "Numerical Recipes" split. Accurate to
/// roughly `1e-12`.
pub(crate) fn regularized_lower_gamma(s: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x < s + 1.0 {
        // Series: P(s,x) = x^s e^{-x} / Γ(s) · Σ x^n / (s)(s+1)…(s+n).
        let mut term = 1.0 / s;
        let mut sum = term;
        let mut n = s;
        for _ in 0..1000 {
            n += 1.0;
            term *= x / n;
            sum += term;
            if term.abs() < sum.abs() * 1e-15 {
                break;
            }
        }
        (sum * (-x + s * x.ln() - ln_gamma(s)).exp()).clamp(0.0, 1.0)
    } else {
        // Continued fraction for the upper Q, then P = 1 - Q.
        let mut b = x + 1.0 - s;
        let mut c = 1e300;
        let mut d = 1.0 / b;
        let mut h = d;
        for i in 1..1000 {
            let an = -(i as f64) * (i as f64 - s);
            b += 2.0;
            d = an * d + b;
            if d.abs() < 1e-300 {
                d = 1e-300;
            }
            c = b + an / c;
            if c.abs() < 1e-300 {
                c = 1e-300;
            }
            d = 1.0 / d;
            let del = d * c;
            h *= del;
            if (del - 1.0).abs() < 1e-15 {
                break;
            }
        }
        let q = (-x + s * x.ln() - ln_gamma(s)).exp() * h;
        (1.0 - q).clamp(0.0, 1.0)
    }
}

/// Natural log of the gamma function (Lanczos approximation, g = 7).
pub(crate) fn ln_gamma(x: f64) -> f64 {
    // Lanczos coefficients (g = 7, n = 9).
    const C: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_403,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_1,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        // Reflection formula.
        let pi = std::f64::consts::PI;
        (pi / (pi * x).sin()).ln() - ln_gamma(1.0 - x)
    } else {
        let x = x - 1.0;
        let mut a = C[0];
        let t = x + 7.5;
        for (i, &c) in C.iter().enumerate().skip(1) {
            a += c / (x + i as f64);
        }
        0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
    }
}

/// Inverse CDF (quantile) of a Gamma distribution with the given
/// `shape` and `scale`, found by bisection on
/// [`regularized_lower_gamma`].
fn gamma_quantile(p: f64, shape: f64, scale: f64) -> f64 {
    if p <= 0.0 {
        return 0.0;
    }
    if p >= 1.0 {
        return f64::INFINITY;
    }
    // Bracket the root in x/scale space.
    let mut lo = 0.0;
    let mut hi = 1.0;
    while regularized_lower_gamma(shape, hi) < p {
        hi *= 2.0;
        if hi > 1e12 {
            break;
        }
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if regularized_lower_gamma(shape, mid) < p {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi) * scale
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_bad_parameters() {
        assert!(DiscreteGamma::new(-1.0, 4).is_err());
        assert!(DiscreteGamma::new(0.5, 0).is_err());
    }

    #[test]
    fn category_rates_have_mean_one() {
        for alpha in [0.1, 0.5, 1.0, 2.0, 5.0] {
            for k in [2, 4, 8] {
                let g = DiscreteGamma::new(alpha, k).unwrap();
                let mean: f64 = g.rates().iter().sum::<f64>() / k as f64;
                assert!((mean - 1.0).abs() < 1e-9, "alpha={alpha} k={k} mean={mean}");
            }
        }
    }

    #[test]
    fn rates_are_sorted_and_positive() {
        let g = DiscreteGamma::new(0.5, 6).unwrap();
        for &r in g.rates() {
            assert!(r >= 0.0);
        }
        // The category means are monotonically increasing.
        for w in g.rates().windows(2) {
            assert!(w[1] >= w[0] - 1e-12, "rates not sorted: {:?}", g.rates());
        }
    }

    #[test]
    fn large_alpha_approaches_uniform_rates() {
        // As α -> ∞ the Gamma collapses to a point mass at 1: its
        // coefficient of variation is 1/√α, so the discrete-category
        // spread shrinks as √α. α = 100 (CV = 0.1) still leaves the
        // slowest of four categories near 0.87 — a real, correct ~13 %
        // spread — so a meaningful "approaches 1" test needs a genuinely
        // large α. At α = 10000 (CV = 0.01) every category is within
        // ~2 % of 1.
        let g = DiscreteGamma::new(10_000.0, 4).unwrap();
        for &r in g.rates() {
            assert!((r - 1.0).abs() < 0.05, "rate {r} far from 1");
        }
        // The spread genuinely narrows with α: α=10000 is far tighter
        // than α=100.
        let wide = DiscreteGamma::new(100.0, 4).unwrap();
        let spread = |g: &DiscreteGamma| {
            let rs = g.rates();
            rs.iter().cloned().fold(f64::MIN, f64::max)
                - rs.iter().cloned().fold(f64::MAX, f64::min)
        };
        assert!(spread(&g) < spread(&wide));
    }

    #[test]
    fn small_alpha_gives_strong_heterogeneity() {
        // Small α => a wide spread between the slowest and fastest
        // category.
        let g = DiscreteGamma::new(0.1, 4).unwrap();
        let min = g.rates().first().copied().unwrap();
        let max = g.rates().last().copied().unwrap();
        assert!(max / min.max(1e-9) > 10.0, "expected a wide spread");
    }

    #[test]
    fn single_category_is_rate_one() {
        let g = DiscreteGamma::new(0.5, 1).unwrap();
        assert_eq!(g.rates(), &[1.0]);
        assert!((g.category_probability() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn ln_gamma_matches_known_values() {
        // Γ(1) = 1, Γ(5) = 24, Γ(0.5) = √π.
        assert!(ln_gamma(1.0).abs() < 1e-9);
        assert!((ln_gamma(5.0) - 24.0_f64.ln()).abs() < 1e-9);
        assert!((ln_gamma(0.5) - std::f64::consts::PI.sqrt().ln()).abs() < 1e-9);
    }

    #[test]
    fn regularized_gamma_is_a_cdf() {
        // P(s, 0) = 0, P(s, ∞) -> 1, monotone increasing.
        assert!(regularized_lower_gamma(2.0, 0.0).abs() < 1e-12);
        assert!((regularized_lower_gamma(2.0, 1e6) - 1.0).abs() < 1e-9);
        let mut prev = 0.0;
        for i in 1..50 {
            let v = regularized_lower_gamma(2.0, i as f64 * 0.2);
            assert!(v >= prev - 1e-12);
            prev = v;
        }
    }
}
