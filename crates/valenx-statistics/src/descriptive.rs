//! Descriptive statistics — location, dispersion and order statistics.
//!
//! Everything here is a textbook closed form over an `f64` slice:
//!
//! - **Location:** the arithmetic [`mean`].
//! - **Dispersion:** [`sample_variance`] (Bessel-corrected, `n - 1`) and
//!   [`population_variance`] (`n`), with their square-root standard
//!   deviations [`sample_std`] / [`population_std`].
//! - **Order statistics:** the [`median`] and the [`quartiles`] (and the
//!   general [`quantile`]) using the **linear-interpolation type-7** rule —
//!   the default of NumPy `percentile`, R `quantile`, and most spreadsheets.
//! - **Standardisation:** the [`z_score`] `(x − μ) / σ`.
//!
//! ## Sample vs population
//!
//! The *sample* variance divides the sum of squared deviations by `n − 1`
//! (Bessel's correction), giving an unbiased estimator of the variance of
//! the population the sample was drawn from. The *population* variance
//! divides by `n` and is the second central moment of the data treated as
//! the whole population. Pick the sample form when generalising beyond the
//! data in hand, the population form when the data *is* the population.

use crate::error::{Result, StatsError};
use crate::validate;

/// Arithmetic mean (average) of a sample.
///
/// `mean = (1 / n) · Σ xᵢ`.
///
/// ## Errors
///
/// Returns [`StatsError::EmptySample`] if `data` is empty and
/// [`StatsError::NonFinite`] if any element is not finite.
///
/// ## Examples
///
/// The mean of the integers `1..=n` is the well-known `(n + 1) / 2`:
///
/// ```
/// use valenx_statistics::descriptive::mean;
/// let data: Vec<f64> = (1..=100).map(|i| i as f64).collect();
/// let m = mean(&data).unwrap();
/// assert!((m - 50.5).abs() < 1e-12);
/// ```
pub fn mean(data: &[f64]) -> Result<f64> {
    validate::non_empty(data, "mean")?;
    validate::all_finite(data, "sample")?;
    let n = data.len() as f64;
    let sum: f64 = data.iter().sum();
    Ok(sum / n)
}

/// Sum of squared deviations from the mean, `Σ (xᵢ − μ)²`.
///
/// This is the shared numerator of both variance forms; computing it once
/// keeps [`sample_variance`] and [`population_variance`] consistent. A
/// two-pass formula is used (mean first, then deviations) for numerical
/// stability over the naïve `Σx² − (Σx)²/n` shortcut.
fn sum_squared_deviations(data: &[f64]) -> Result<f64> {
    let m = mean(data)?;
    Ok(data.iter().map(|&x| (x - m) * (x - m)).sum())
}

/// Sample variance with Bessel's correction (divisor `n − 1`).
///
/// `s² = (1 / (n − 1)) · Σ (xᵢ − μ)²`. This is the unbiased estimator of
/// the population variance from a finite sample.
///
/// ## Errors
///
/// Returns [`StatsError::TooFewObservations`] unless `n >= 2` (the divisor
/// `n − 1` would otherwise be zero), and [`StatsError::NonFinite`] for any
/// non-finite element.
///
/// ## Examples
///
/// ```
/// use valenx_statistics::descriptive::sample_variance;
/// // Population {2, 4, 4, 4, 5, 5, 7, 9}: sample variance = 32/7.
/// let data = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
/// let v = sample_variance(&data).unwrap();
/// assert!((v - 32.0 / 7.0).abs() < 1e-12);
/// ```
pub fn sample_variance(data: &[f64]) -> Result<f64> {
    validate::at_least(data, 2, "sample_variance")?;
    let ssd = sum_squared_deviations(data)?;
    Ok(ssd / (data.len() as f64 - 1.0))
}

/// Population variance (divisor `n`).
///
/// `σ² = (1 / n) · Σ (xᵢ − μ)²` — the second central moment, treating the
/// data as the entire population.
///
/// ## Errors
///
/// Returns [`StatsError::EmptySample`] if `data` is empty and
/// [`StatsError::NonFinite`] for any non-finite element.
///
/// ## Examples
///
/// ```
/// use valenx_statistics::descriptive::population_variance;
/// // Population {2, 4, 4, 4, 5, 5, 7, 9}: population variance = 4.
/// let data = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
/// let v = population_variance(&data).unwrap();
/// assert!((v - 4.0).abs() < 1e-12);
/// ```
pub fn population_variance(data: &[f64]) -> Result<f64> {
    validate::non_empty(data, "population_variance")?;
    let ssd = sum_squared_deviations(data)?;
    Ok(ssd / data.len() as f64)
}

/// Sample standard deviation — the square root of [`sample_variance`].
///
/// `s = √s²`. Carries the same `n >= 2` and finiteness preconditions.
///
/// ## Examples
///
/// ```
/// use valenx_statistics::descriptive::{sample_std, sample_variance};
/// let data = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
/// let s = sample_std(&data).unwrap();
/// assert!((s - sample_variance(&data).unwrap().sqrt()).abs() < 1e-12);
/// ```
pub fn sample_std(data: &[f64]) -> Result<f64> {
    Ok(sample_variance(data)?.sqrt())
}

/// Population standard deviation — the square root of
/// [`population_variance`].
///
/// `σ = √σ²`. Carries the same non-empty and finiteness preconditions.
///
/// ## Examples
///
/// ```
/// use valenx_statistics::descriptive::population_std;
/// // Population {2, 4, 4, 4, 5, 5, 7, 9}: σ² = 4, so σ = 2.
/// let data = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
/// let s = population_std(&data).unwrap();
/// assert!((s - 2.0).abs() < 1e-12);
/// ```
pub fn population_std(data: &[f64]) -> Result<f64> {
    Ok(population_variance(data)?.sqrt())
}

/// Sort a validated copy of `data` ascending, ready for order statistics.
///
/// Uses [`f64::total_cmp`] so the sort is total even though `f64` is only
/// partially ordered; finiteness has already been checked by the caller,
/// so no `NaN` can actually reach the comparator.
fn sorted_copy(data: &[f64]) -> Vec<f64> {
    let mut v = data.to_vec();
    v.sort_by(f64::total_cmp);
    v
}

/// The `q`-quantile of a sample by **linear interpolation (type 7)**.
///
/// On the already-sorted order statistics `x₍₀₎ ≤ … ≤ x₍ₙ₋₁₎`, the position
/// is `h = q · (n − 1)`; the result interpolates linearly between the two
/// neighbours bracketing `h`: `x₍⌊h⌋₎ + (h − ⌊h⌋) · (x₍⌊h⌋₊₁₎ − x₍⌊h⌋₎)`.
/// `q = 0` gives the minimum, `q = 1` the maximum, `q = 0.5` the [`median`].
///
/// This is the NumPy / R type-7 default. With `q ∈ {0.25, 0.5, 0.75}` it
/// agrees with [`quartiles`].
///
/// ## Errors
///
/// Returns [`StatsError::EmptySample`] for empty input,
/// [`StatsError::OutOfRange`] if `q ∉ [0, 1]`, and [`StatsError::NonFinite`]
/// for a non-finite `q` or element.
///
/// ## Examples
///
/// ```
/// use valenx_statistics::descriptive::quantile;
/// let data = [1.0, 2.0, 3.0, 4.0];
/// assert!((quantile(&data, 0.0).unwrap() - 1.0).abs() < 1e-12);
/// assert!((quantile(&data, 1.0).unwrap() - 4.0).abs() < 1e-12);
/// // h = 0.5·3 = 1.5 → between x[1]=2 and x[2]=3 → 2.5.
/// assert!((quantile(&data, 0.5).unwrap() - 2.5).abs() < 1e-12);
/// ```
pub fn quantile(data: &[f64], q: f64) -> Result<f64> {
    validate::non_empty(data, "quantile")?;
    validate::all_finite(data, "sample")?;
    validate::unit_closed(q, "q")?;

    let sorted = sorted_copy(data);
    let n = sorted.len();
    if n == 1 {
        return Ok(sorted[0]);
    }

    let h = q * (n as f64 - 1.0);
    let lo = h.floor();
    let lo_idx = lo as usize;
    // `lo_idx` is in `0..=n-1`; when it is the last index, `frac` is 0 so
    // the upper neighbour is never indexed out of bounds.
    let frac = h - lo;
    if lo_idx + 1 >= n {
        return Ok(sorted[n - 1]);
    }
    let lower = sorted[lo_idx];
    let upper = sorted[lo_idx + 1];
    Ok(lower + frac * (upper - lower))
}

/// The median — the `0.5`-quantile (type 7).
///
/// For an odd count this is the middle order statistic; for an even count
/// it is the average of the two central order statistics.
///
/// ## Errors
///
/// Returns [`StatsError::EmptySample`] for empty input and
/// [`StatsError::NonFinite`] for any non-finite element.
///
/// ## Examples
///
/// ```
/// use valenx_statistics::descriptive::median;
/// assert!((median(&[3.0, 1.0, 2.0]).unwrap() - 2.0).abs() < 1e-12);
/// assert!((median(&[1.0, 2.0, 3.0, 4.0]).unwrap() - 2.5).abs() < 1e-12);
/// ```
pub fn median(data: &[f64]) -> Result<f64> {
    quantile(data, 0.5)
}

/// The first quartile (`Q1`), median (`Q2`) and third quartile (`Q3`) of a
/// sample, computed by the type-7 [`quantile`] rule.
#[derive(Copy, Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Quartiles {
    /// First quartile — the `0.25`-quantile (25th percentile).
    pub q1: f64,
    /// Second quartile — the `0.50`-quantile, i.e. the [`median`].
    pub q2: f64,
    /// Third quartile — the `0.75`-quantile (75th percentile).
    pub q3: f64,
}

impl Quartiles {
    /// The interquartile range `IQR = Q3 − Q1`, the spread of the central
    /// half of the data. Always non-negative because `Q3 >= Q1`.
    pub fn iqr(&self) -> f64 {
        self.q3 - self.q1
    }
}

/// Compute the [`Quartiles`] of a sample.
///
/// ## Errors
///
/// Returns [`StatsError::EmptySample`] for empty input and
/// [`StatsError::NonFinite`] for any non-finite element.
///
/// ## Examples
///
/// ```
/// use valenx_statistics::descriptive::quartiles;
/// let data = [1.0, 2.0, 3.0, 4.0, 5.0];
/// let q = quartiles(&data).unwrap();
/// assert!((q.q1 - 2.0).abs() < 1e-12);
/// assert!((q.q2 - 3.0).abs() < 1e-12);
/// assert!((q.q3 - 4.0).abs() < 1e-12);
/// assert!((q.iqr() - 2.0).abs() < 1e-12);
/// ```
pub fn quartiles(data: &[f64]) -> Result<Quartiles> {
    Ok(Quartiles {
        q1: quantile(data, 0.25)?,
        q2: quantile(data, 0.50)?,
        q3: quantile(data, 0.75)?,
    })
}

/// The z-score (standard score) of a single observation.
///
/// `z = (x − μ) / σ`. It expresses how many standard deviations `x` lies
/// above (`z > 0`) or below (`z < 0`) the mean `μ`.
///
/// ## Errors
///
/// Returns [`StatsError::NonFinite`] if `x` or `mean` is not finite, and
/// [`StatsError::NonPositiveScale`] if `std <= 0` (a zero or negative
/// spread cannot standardise).
///
/// ## Examples
///
/// ```
/// use valenx_statistics::descriptive::z_score;
/// // 130 on a μ=100, σ=15 scale is exactly two standard deviations up.
/// let z = z_score(130.0, 100.0, 15.0).unwrap();
/// assert!((z - 2.0).abs() < 1e-12);
/// ```
pub fn z_score(x: f64, mean: f64, std: f64) -> Result<f64> {
    validate::finite(x, "x")?;
    validate::finite(mean, "mean")?;
    validate::positive_scale(std, "std")?;
    Ok((x - mean) / std)
}

/// Standardise an entire sample against its **own** sample mean and sample
/// standard deviation, returning the z-score of each element.
///
/// Equivalent to calling [`z_score`] with `μ = mean(data)` and
/// `σ = sample_std(data)` for every element. The returned scores have a
/// (sample) mean of zero and a sample standard deviation of one.
///
/// ## Errors
///
/// Returns [`StatsError::TooFewObservations`] unless `n >= 2` (the sample
/// standard deviation needs at least two points), and
/// [`StatsError::NonPositiveScale`] if the sample is constant (zero spread).
///
/// ## Examples
///
/// ```
/// use valenx_statistics::descriptive::{standardize, mean};
/// let data = [1.0, 2.0, 3.0, 4.0, 5.0];
/// let z = standardize(&data).unwrap();
/// // Standardised data has (sample) mean 0.
/// assert!(mean(&z).unwrap().abs() < 1e-12);
/// ```
pub fn standardize(data: &[f64]) -> Result<Vec<f64>> {
    let m = mean(data)?;
    let s = sample_std(data)?;
    if s == 0.0 {
        return Err(StatsError::NonPositiveScale {
            name: "std",
            value: 0.0,
        });
    }
    Ok(data.iter().map(|&x| (x - m) / s).collect())
}
