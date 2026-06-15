//! Error taxonomy for the thermodynamic-cycle calculators.
//!
//! Every fallible constructor or efficiency routine returns
//! [`Result<T, CycleError>`](Result). The error carries a stable
//! [`code`](CycleError::code) and a coarse [`category`](CycleError::category)
//! for telemetry and UI grouping.

use thiserror::Error;

/// Shorthand for `Result<T, CycleError>`.
pub type Result<T> = core::result::Result<T, CycleError>;

/// Anything that can go wrong validating cycle inputs or computing an
/// efficiency.
///
/// This enum is `#[non_exhaustive]`: new variants may be added without it
/// being a breaking change, so downstream `match` arms must include a
/// wildcard.
#[derive(Debug, Error, Clone, PartialEq)]
#[non_exhaustive]
pub enum CycleError {
    /// A named scalar parameter was non-finite (`NaN` / `±∞`).
    #[error("parameter `{name}` is not finite (got {value})")]
    NotFinite {
        /// Parameter name (e.g. `"t_hot"`, `"gamma"`).
        name: &'static str,
        /// The offending value.
        value: f64,
    },

    /// A parameter that must be strictly positive was zero or negative.
    ///
    /// Used for absolute temperatures (kelvin), the heat-capacity ratio
    /// `γ`, compression / pressure ratios, and the cutoff ratio.
    #[error("parameter `{name}` must be > {floor} (got {value})")]
    NotPositive {
        /// Parameter name.
        name: &'static str,
        /// The exclusive lower bound that was violated.
        floor: f64,
        /// The offending value.
        value: f64,
    },

    /// A pair of temperatures was ordered wrongly: the cold-reservoir
    /// temperature was not strictly below the hot-reservoir temperature,
    /// so the heat engine has no temperature drop to exploit.
    #[error("cold reservoir T_c = {t_cold} K must be < hot reservoir T_h = {t_hot} K")]
    TemperatureOrder {
        /// Cold-reservoir absolute temperature in kelvin.
        t_cold: f64,
        /// Hot-reservoir absolute temperature in kelvin.
        t_hot: f64,
    },

    /// The heat-capacity ratio `γ` was out of the physically meaningful
    /// open interval `(1, ∞)`. A real gas has `γ = c_p / c_v > 1`; the
    /// air-standard efficiency formulas degenerate (zero efficiency) at
    /// `γ = 1` and are non-physical below it.
    #[error("heat-capacity ratio gamma = {value} must be > 1")]
    GammaTooLow {
        /// The offending `γ`.
        value: f64,
    },

    /// A compression / pressure ratio was not greater than one. With a
    /// ratio of exactly `1` the cycle does no net work (zero efficiency),
    /// and below `1` the air-standard formula is non-physical.
    #[error("ratio `{name}` = {value} must be > 1 (no compression otherwise)")]
    RatioTooLow {
        /// Which ratio (`"compression_ratio"`, `"pressure_ratio"`).
        name: &'static str,
        /// The offending value.
        value: f64,
    },

    /// A Rankine-cycle enthalpy set was inconsistent: the heat added in
    /// the boiler, `h3 - h2`, was not strictly positive, so the cycle
    /// absorbs no heat and its efficiency is undefined (division by a
    /// non-positive heat input).
    #[error(
        "boiler heat input h3 - h2 = {q_in} kJ/kg must be > 0 \
         (turbine inlet enthalpy must exceed pump outlet enthalpy)"
    )]
    NoHeatInput {
        /// The computed boiler heat input `h3 - h2` in kJ/kg.
        q_in: f64,
    },
}

/// Coarse error category for grouping in a UI or telemetry pipeline.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// A user-supplied input value was invalid (out of range, non-finite,
    /// or inconsistent with another input).
    Input,
    /// The requested computation has no valid result for these inputs
    /// (a domain / degeneracy condition of the model itself).
    Domain,
}

impl CycleError {
    /// Stable, kebab-cased identifier suitable for logs and tests.
    ///
    /// The string is part of the crate's contract and will not change for
    /// an existing variant.
    pub fn code(&self) -> &'static str {
        match self {
            CycleError::NotFinite { .. } => "thermocycle.not-finite",
            CycleError::NotPositive { .. } => "thermocycle.not-positive",
            CycleError::TemperatureOrder { .. } => "thermocycle.temperature-order",
            CycleError::GammaTooLow { .. } => "thermocycle.gamma-too-low",
            CycleError::RatioTooLow { .. } => "thermocycle.ratio-too-low",
            CycleError::NoHeatInput { .. } => "thermocycle.no-heat-input",
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            CycleError::NotFinite { .. }
            | CycleError::NotPositive { .. }
            | CycleError::TemperatureOrder { .. }
            | CycleError::GammaTooLow { .. }
            | CycleError::RatioTooLow { .. } => ErrorCategory::Input,
            CycleError::NoHeatInput { .. } => ErrorCategory::Domain,
        }
    }
}

/// Validate that `value` is finite, returning [`CycleError::NotFinite`]
/// otherwise. Internal helper shared by the cycle constructors.
pub(crate) fn finite(name: &'static str, value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(CycleError::NotFinite { name, value })
    }
}

/// Validate that `value` is finite and strictly positive, returning the
/// appropriate [`CycleError`] otherwise. Internal helper.
pub(crate) fn positive(name: &'static str, value: f64) -> Result<f64> {
    let v = finite(name, value)?;
    if v > 0.0 {
        Ok(v)
    } else {
        Err(CycleError::NotPositive {
            name,
            floor: 0.0,
            value: v,
        })
    }
}
