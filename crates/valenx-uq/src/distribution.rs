//! Input probability distributions.
//!
//! Each uncertain model input is described by one [`Distribution`]. A
//! distribution can be **sampled** (draw one realisation, using the crate's
//! deterministic [`SplitMix64`]) and can report its **mean** in closed form.
//! Constructors validate their parameters and return [`UqError`] on bad input,
//! so an impossible distribution can never be built.

use crate::error::UqError;
use crate::rng::SplitMix64;

/// A univariate probability distribution for a single model input.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum Distribution {
    /// Continuous uniform on `[lo, hi]`.
    Uniform {
        /// Lower bound (inclusive).
        lo: f64,
        /// Upper bound. Required `> lo`.
        hi: f64,
    },
    /// Gaussian / normal with the given mean and standard deviation.
    Normal {
        /// Mean.
        mean: f64,
        /// Standard deviation. Required `> 0`.
        std: f64,
    },
    /// Triangular distribution on `[lo, hi]` peaking at `mode`.
    Triangular {
        /// Lower bound.
        lo: f64,
        /// Mode (peak). Required `lo <= mode <= hi`.
        mode: f64,
        /// Upper bound. Required `> lo`.
        hi: f64,
    },
}

impl Distribution {
    /// Construct a validated [`Distribution::Uniform`] (`lo < hi`).
    ///
    /// # Errors
    /// Returns [`UqError::InvalidDistribution`] if `lo >= hi` or either bound
    /// is non-finite.
    pub fn uniform(lo: f64, hi: f64) -> Result<Self, UqError> {
        if !lo.is_finite() || !hi.is_finite() {
            return Err(UqError::InvalidDistribution(format!(
                "uniform bounds must be finite (got lo={lo}, hi={hi})"
            )));
        }
        if lo >= hi {
            return Err(UqError::InvalidDistribution(format!(
                "uniform requires lo < hi (got lo={lo}, hi={hi})"
            )));
        }
        Ok(Self::Uniform { lo, hi })
    }

    /// Construct a validated [`Distribution::Normal`] (`std > 0`).
    ///
    /// # Errors
    /// Returns [`UqError::InvalidDistribution`] if `std <= 0` or any parameter
    /// is non-finite.
    pub fn normal(mean: f64, std: f64) -> Result<Self, UqError> {
        if !mean.is_finite() || !std.is_finite() {
            return Err(UqError::InvalidDistribution(format!(
                "normal parameters must be finite (got mean={mean}, std={std})"
            )));
        }
        if std <= 0.0 {
            return Err(UqError::InvalidDistribution(format!(
                "normal requires std > 0 (got std={std})"
            )));
        }
        Ok(Self::Normal { mean, std })
    }

    /// Construct a validated [`Distribution::Triangular`]
    /// (`lo < hi` and `lo <= mode <= hi`).
    ///
    /// # Errors
    /// Returns [`UqError::InvalidDistribution`] if `lo >= hi`, the mode lies
    /// outside `[lo, hi]`, or any parameter is non-finite.
    pub fn triangular(lo: f64, mode: f64, hi: f64) -> Result<Self, UqError> {
        if !lo.is_finite() || !mode.is_finite() || !hi.is_finite() {
            return Err(UqError::InvalidDistribution(format!(
                "triangular parameters must be finite (got lo={lo}, mode={mode}, hi={hi})"
            )));
        }
        if lo >= hi {
            return Err(UqError::InvalidDistribution(format!(
                "triangular requires lo < hi (got lo={lo}, hi={hi})"
            )));
        }
        if mode < lo || mode > hi {
            return Err(UqError::InvalidDistribution(format!(
                "triangular requires lo <= mode <= hi (got lo={lo}, mode={mode}, hi={hi})"
            )));
        }
        Ok(Self::Triangular { lo, mode, hi })
    }

    /// The closed-form mean of the distribution.
    #[must_use]
    pub fn mean(&self) -> f64 {
        match *self {
            Distribution::Uniform { lo, hi } => 0.5 * (lo + hi),
            Distribution::Normal { mean, .. } => mean,
            // Mean of a triangular distribution is (lo + mode + hi) / 3.
            Distribution::Triangular { lo, mode, hi } => (lo + mode + hi) / 3.0,
        }
    }

    /// Draw one sample using the supplied deterministic PRNG.
    pub fn sample(&self, rng: &mut SplitMix64) -> f64 {
        match *self {
            Distribution::Uniform { lo, hi } => rng.next_range(lo, hi),
            Distribution::Normal { mean, std } => mean + std * rng.next_standard_normal(),
            Distribution::Triangular { lo, mode, hi } => {
                // Inverse-CDF (Smirnov) transform of a U(0,1) draw.
                let u = rng.next_f64();
                inverse_triangular_cdf(u, lo, mode, hi)
            }
        }
    }

    /// Map a probability `u ∈ [0, 1]` to the value at that quantile (the
    /// inverse CDF). Used internally by Latin-hypercube sampling, which needs
    /// to place a draw within a chosen probability stratum.
    ///
    /// `u` is clamped into `[0, 1]`.
    pub(crate) fn quantile(&self, u: f64) -> f64 {
        let u = u.clamp(0.0, 1.0);
        match *self {
            Distribution::Uniform { lo, hi } => lo + (hi - lo) * u,
            Distribution::Normal { mean, std } => mean + std * standard_normal_quantile(u),
            Distribution::Triangular { lo, mode, hi } => inverse_triangular_cdf(u, lo, mode, hi),
        }
    }
}

/// Inverse CDF of the triangular distribution on `[lo, hi]` with peak `mode`.
fn inverse_triangular_cdf(u: f64, lo: f64, mode: f64, hi: f64) -> f64 {
    let span = hi - lo;
    // Degenerate guard: a zero-width support collapses to the point.
    if span <= 0.0 {
        return lo;
    }
    // Fraction of the support to the left of the mode; this equals F(mode).
    let fc = (mode - lo) / span;
    if u < fc {
        lo + (u * span * (mode - lo)).sqrt()
    } else {
        hi - ((1.0 - u) * span * (hi - mode)).sqrt()
    }
}

/// Inverse CDF (quantile function) of the standard normal `N(0, 1)`.
///
/// Uses Acklam's rational approximation, which is accurate to roughly
/// `1.15e-9` in absolute value over the whole open interval `(0, 1)` — far
/// tighter than the sampling error of any practical UQ run. The two endpoints
/// are clamped to a large finite magnitude so the result is never infinite.
fn standard_normal_quantile(p: f64) -> f64 {
    // Coefficients for Acklam's algorithm (each is the exact shortest
    // round-trip `f64` literal, so no precision is lost vs. the published
    // constants).
    const A: [f64; 6] = [
        -39.696_830_286_653_76,
        220.946_098_424_520_5,
        -275.928_510_446_968_7,
        138.357_751_867_269,
        -30.664_798_066_147_16,
        2.506_628_277_459_239,
    ];
    const B: [f64; 5] = [
        -54.476_098_798_224_06,
        161.585_836_858_040_9,
        -155.698_979_859_886_6,
        66.801_311_887_719_72,
        -13.280_681_552_885_72,
    ];
    const C: [f64; 6] = [
        -0.007_784_894_002_430_293,
        -0.322_396_458_041_136_5,
        -2.400_758_277_161_838,
        -2.549_732_539_343_734,
        4.374_664_141_464_968,
        2.938_163_982_698_783,
    ];
    const D: [f64; 4] = [
        0.007_784_695_709_041_462,
        0.322_467_129_070_039_8,
        2.445_134_137_142_996,
        3.754_408_661_907_416,
    ];
    // Break-points between the central and tail rational approximations.
    const P_LOW: f64 = 0.024_25;
    const P_HIGH: f64 = 1.0 - P_LOW;

    if p <= 0.0 {
        return -1e10;
    }
    if p >= 1.0 {
        return 1e10;
    }

    if p < P_LOW {
        // Lower tail.
        let q = (-2.0 * p.ln()).sqrt();
        (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    } else if p <= P_HIGH {
        // Central region.
        let q = p - 0.5;
        let r = q * q;
        (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q
            / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0)
    } else {
        // Upper tail.
        let q = (-2.0 * (1.0 - p).ln()).sqrt();
        -(((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    }
}
