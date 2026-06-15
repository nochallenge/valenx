//! Inferential statistics — the one-sample t-statistic, the inverse
//! standard-normal quantile, and a normal (known-sigma) confidence interval
//! for the mean.
//!
//! ## Model
//!
//! - **One-sample t-statistic** for a hypothesised mean `μ₀`:
//!   `t = (x̄ − μ₀) / (s / √n)`, where `x̄` is the sample mean, `s` the
//!   sample standard deviation and `n` the sample size. The denominator
//!   `s / √n` is the standard error of the mean.
//! - **Inverse normal CDF** ([`probit`]) `Φ⁻¹(p)` by the Acklam rational
//!   approximation, used to turn a confidence level into a critical value.
//! - **Confidence interval** for the mean with *known* population standard
//!   deviation `σ`: `x̄ ± z · σ / √n`, where the critical value `z` is the
//!   `(1 + c) / 2` quantile of `N(0, 1)` for confidence level `c`.
//!
//! ## Honest scope
//!
//! The interval here is the **z-interval** (population spread assumed
//! known). A small-sample interval with *estimated* spread would use the
//! Student-t critical value instead; the t-distribution quantile is not
//! implemented, so for unknown `σ` and small `n` this interval is slightly
//! too narrow. The [`probit`] critical value carries the Acklam
//! approximation's `≈ 1.15 × 10⁻⁹` relative error, ample for teaching but
//! not a certified inverse-CDF.

use crate::descriptive::{mean, sample_std};
use crate::error::{Result, StatsError};
use crate::validate;

/// The one-sample t-statistic for a hypothesised population mean `mu0`.
///
/// `t = (x̄ − μ₀) / (s / √n)`. A `t` near zero is consistent with the null
/// hypothesis `μ = μ₀`; large `|t|` is evidence against it. The sign
/// follows whether the sample mean sits above or below `μ₀`.
///
/// ## Errors
///
/// Returns [`StatsError::TooFewObservations`] unless `n >= 2` (the sample
/// standard deviation needs at least two points), [`StatsError::NonFinite`]
/// if `mu0` or any element is not finite, and [`StatsError::NonPositiveScale`]
/// if the sample is constant (zero standard error).
///
/// ## Examples
///
/// ```
/// use valenx_statistics::inference::t_statistic;
/// // Sample mean equal to the hypothesised mean gives t = 0.
/// let data = [4.0, 5.0, 6.0];
/// assert!(t_statistic(&data, 5.0).unwrap().abs() < 1e-12);
/// ```
pub fn t_statistic(data: &[f64], mu0: f64) -> Result<f64> {
    validate::at_least(data, 2, "t_statistic")?;
    validate::finite(mu0, "mu0")?;
    let xbar = mean(data)?;
    let s = sample_std(data)?;
    if s == 0.0 {
        return Err(StatsError::NonPositiveScale {
            name: "std",
            value: 0.0,
        });
    }
    let n = data.len() as f64;
    let standard_error = s / n.sqrt();
    Ok((xbar - mu0) / standard_error)
}

/// The inverse of the standard-normal CDF, `Φ⁻¹(p)` (the *probit* / normal
/// quantile), by Peter Acklam's rational approximation.
///
/// Returns the `z` such that `Φ(z) = p`. The approximation has relative
/// error below `≈ 1.15 × 10⁻⁹` across `p ∈ (0, 1)`. It is **odd about the
/// median**: `Φ⁻¹(1 − p) = −Φ⁻¹(p)`, and `Φ⁻¹(0.5) = 0`.
///
/// ## Errors
///
/// Returns [`StatsError::OutOfRange`] if `p` is not in the open interval
/// `(0, 1)` (the quantile diverges to `∓∞` at the endpoints) and
/// [`StatsError::NonFinite`] if `p` is not finite.
///
/// ## Examples
///
/// ```
/// use valenx_statistics::inference::probit;
/// // Median of the standard normal.
/// assert!(probit(0.5).unwrap().abs() < 1e-9);
/// // The classic 97.5% quantile ≈ 1.959964 (the "1.96" of a 95% interval).
/// assert!((probit(0.975).unwrap() - 1.959_963_985).abs() < 1e-6);
/// ```
pub fn probit(p: f64) -> Result<f64> {
    validate::finite(p, "p")?;
    if !(p > 0.0 && p < 1.0) {
        return Err(StatsError::OutOfRange {
            name: "p",
            value: p,
            expected: "(0, 1)",
        });
    }

    // Acklam's coefficients.
    const A: [f64; 6] = [
        -3.969_683_028_665_376e1,
        2.209_460_984_245_205e2,
        -2.759_285_104_469_687e2,
        1.383_577_518_672_69e2,
        -3.066_479_806_614_716e1,
        2.506_628_277_459_239e0,
    ];
    const B: [f64; 5] = [
        -5.447_609_879_822_406e1,
        1.615_858_368_580_409e2,
        -1.556_989_798_598_866e2,
        6.680_131_188_771_972e1,
        -1.328_068_155_288_572e1,
    ];
    const C: [f64; 6] = [
        -7.784_894_002_430_293e-3,
        -3.223_964_580_411_365e-1,
        -2.400_758_277_161_838e0,
        -2.549_732_539_343_734e0,
        4.374_664_141_464_968e0,
        2.938_163_982_698_783e0,
    ];
    const D: [f64; 4] = [
        7.784_695_709_041_462e-3,
        3.224_671_290_700_398e-1,
        2.445_134_137_142_996e0,
        3.754_408_661_907_416e0,
    ];

    // Break points between the central rational region and the tails.
    const P_LOW: f64 = 0.024_25;
    const P_HIGH: f64 = 1.0 - P_LOW;

    if p < P_LOW {
        // Lower tail: rational approximation in q = √(−2 ln p).
        let q = (-2.0 * p.ln()).sqrt();
        Ok(
            (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
                / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0),
        )
    } else if p <= P_HIGH {
        // Central region: rational approximation in q = p − 0.5.
        let q = p - 0.5;
        let r = q * q;
        Ok(
            (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q
                / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0),
        )
    } else {
        // Upper tail: mirror of the lower-tail expansion.
        let q = (-2.0 * (1.0 - p).ln()).sqrt();
        Ok(
            -(((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
                / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0),
        )
    }
}

/// A closed interval `[lower, upper]` estimate for an unknown mean, at a
/// stated confidence level.
#[derive(Copy, Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ConfidenceInterval {
    /// The point estimate the interval is centred on (the sample mean).
    pub point: f64,
    /// The lower endpoint `point − margin`.
    pub lower: f64,
    /// The upper endpoint `point + margin`.
    pub upper: f64,
    /// The half-width `z · σ / √n` (the *margin of error*); the interval is
    /// `point ± margin`.
    pub margin: f64,
    /// The confidence level the interval was built for, in `(0, 1)`.
    pub confidence: f64,
}

impl ConfidenceInterval {
    /// The full width of the interval, `upper − lower = 2 · margin`.
    pub fn width(&self) -> f64 {
        self.upper - self.lower
    }

    /// Whether a candidate value lies inside the closed interval.
    ///
    /// ## Examples
    ///
    /// ```
    /// use valenx_statistics::inference::normal_ci;
    /// let ci = normal_ci(&[10.0, 10.0, 10.0, 10.0], 2.0, 0.95).unwrap();
    /// assert!(ci.contains(10.0));
    /// ```
    pub fn contains(&self, value: f64) -> bool {
        value >= self.lower && value <= self.upper
    }
}

/// A two-sided confidence interval for the mean assuming a **known**
/// population standard deviation `sigma` (the z-interval).
///
/// The interval is `x̄ ± z · σ / √n`, where `x̄` is the sample mean, `n` the
/// sample size, and `z = Φ⁻¹((1 + c) / 2)` is the critical value for
/// confidence level `c` (so `c = 0.95` uses `z ≈ 1.96`). The margin
/// `z · σ / √n` **grows with `σ`** and **shrinks as `n` increases** (like
/// `1 / √n`), and **widens with the confidence level**.
///
/// ## Errors
///
/// Returns [`StatsError::EmptySample`] for empty `data`,
/// [`StatsError::NonPositiveScale`] if `sigma <= 0`,
/// [`StatsError::OutOfRange`] if `confidence ∉ (0, 1)`, and
/// [`StatsError::NonFinite`] for any non-finite input.
///
/// ## Examples
///
/// ```
/// use valenx_statistics::inference::normal_ci;
/// // n = 4, sigma = 2, 95% level: margin = 1.96 · 2 / 2 = 1.96 (approx).
/// let ci = normal_ci(&[10.0, 10.0, 10.0, 10.0], 2.0, 0.95).unwrap();
/// assert!((ci.point - 10.0).abs() < 1e-12);
/// assert!((ci.margin - 1.959_963_985 * 2.0 / 2.0).abs() < 1e-6);
/// ```
pub fn normal_ci(data: &[f64], sigma: f64, confidence: f64) -> Result<ConfidenceInterval> {
    validate::non_empty(data, "normal_ci")?;
    validate::positive_scale(sigma, "sigma")?;
    validate::unit_open(confidence, "confidence")?;

    let xbar = mean(data)?;
    let n = data.len() as f64;
    // Two-sided critical value: the upper `(1 + c)/2` quantile of N(0,1).
    let z = probit((1.0 + confidence) / 2.0)?;
    let margin = z * sigma / n.sqrt();
    Ok(ConfidenceInterval {
        point: xbar,
        lower: xbar - margin,
        upper: xbar + margin,
        margin,
        confidence,
    })
}
