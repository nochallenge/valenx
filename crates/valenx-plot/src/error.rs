//! Plot workbench error taxonomy.

use thiserror::Error;

/// Errors raised by plotting.
#[derive(Debug, Error)]
pub enum PlotError {
    /// Plot has no series.
    #[error("plot is empty")]
    Empty,

    /// Bad parameter (width <= 0, etc).
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// PNG rendering deferred to a follow-up phase.
    #[error("PNG rendering deferred to Phase 31.5 (would need the `image` crate)")]
    PngDeferred,

    /// IO error.
    #[error("io: {0}")]
    Io(#[from] std::io::Error),

    /// RON error.
    #[error("ron: {0}")]
    Ron(String),
}

/// Coarse error category.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User input.
    Input,
    /// Tunable knob.
    Config,
    /// Not implemented in v1.
    NotImplemented,
    /// Transient / IO / parse.
    Runtime,
}

impl PlotError {
    /// Stable kebab-cased identifier.
    pub fn code(&self) -> &'static str {
        match self {
            PlotError::Empty => "plot.empty",
            PlotError::BadParameter { .. } => "plot.bad_parameter",
            PlotError::PngDeferred => "plot.png_deferred",
            PlotError::Io(_) => "plot.io",
            PlotError::Ron(_) => "plot.ron",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            PlotError::Empty => ErrorCategory::Input,
            PlotError::BadParameter { .. } => ErrorCategory::Config,
            PlotError::PngDeferred => ErrorCategory::NotImplemented,
            PlotError::Io(_) | PlotError::Ron(_) => ErrorCategory::Runtime,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codes_and_cats() {
        assert_eq!(PlotError::Empty.code(), "plot.empty");
        assert_eq!(PlotError::PngDeferred.category(), ErrorCategory::NotImplemented);
        assert_eq!(
            PlotError::BadParameter {
                name: "width",
                reason: "zero".into()
            }
            .category(),
            ErrorCategory::Config
        );
    }
}
