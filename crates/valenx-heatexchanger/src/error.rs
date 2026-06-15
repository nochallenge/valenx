//! Heat-exchanger error taxonomy.

use thiserror::Error;

/// Errors raised by heat-exchanger thermal analysis.
#[derive(Debug, Error)]
pub enum HeatExchangerError {
    /// A scalar input fell outside its valid domain (non-positive area,
    /// negative coefficient, etc.).
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Offending parameter name.
        name: &'static str,
        /// Human-readable reason the value was rejected.
        reason: String,
    },

    /// The supplied inlet temperatures are thermodynamically
    /// inconsistent — the hot stream is not hotter than the cold
    /// stream, so no heat can flow in the assumed direction.
    #[error("inconsistent temperatures: {0}")]
    InconsistentTemperatures(String),

    /// A computed log-mean / effectiveness expression hit a degenerate
    /// branch that the caller's inputs cannot represent (e.g. a
    /// temperature approach crossing zero in a way that inverts the
    /// LMTD logarithm).
    #[error("degenerate configuration: {0}")]
    Degenerate(String),
}

/// Coarse category for routing / metrics.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// Caller supplied an invalid input value.
    Input,
    /// The physical configuration is degenerate / unrepresentable.
    Algorithm,
}

impl HeatExchangerError {
    /// Stable kebab-cased identifier for logs and tests.
    pub fn code(&self) -> &'static str {
        match self {
            HeatExchangerError::BadParameter { .. } => "heatexchanger.bad_parameter",
            HeatExchangerError::InconsistentTemperatures(_) => {
                "heatexchanger.inconsistent_temperatures"
            }
            HeatExchangerError::Degenerate(_) => "heatexchanger.degenerate",
        }
    }

    /// Coarse error category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            HeatExchangerError::BadParameter { .. } => ErrorCategory::Input,
            HeatExchangerError::InconsistentTemperatures(_) => ErrorCategory::Input,
            HeatExchangerError::Degenerate(_) => ErrorCategory::Algorithm,
        }
    }
}
