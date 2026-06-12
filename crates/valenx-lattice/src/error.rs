//! Lattice2 workbench error taxonomy.

use thiserror::Error;

/// Errors raised by lattice generation.
#[derive(Debug, Error)]
pub enum LatticeError {
    /// Caller supplied a degenerate lattice (zero count, empty curve).
    #[error("degenerate lattice: {0}")]
    Degenerate(String),

    /// Parameter out of range.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Offending parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// Round-10 fix: a hostile / mistyped lattice file (e.g.
    /// `Grid { rows: usize::MAX, cols: usize::MAX, levels: 1, ... }`)
    /// could ask `Vec::with_capacity(rows * cols * levels)` for an
    /// allocation that wraps `usize` and then back-fills with billions
    /// of placements. The generator now checked-muls the placement
    /// count and rejects anything past
    /// [`crate::generate::MAX_LATTICE_PLACEMENTS`].
    #[error("lattice would place {count} instances; cap is {max}")]
    TooManyPlacements {
        /// Placement count the recipe would produce.
        count: u128,
        /// Hard cap enforced before allocation.
        max: usize,
    },

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
    /// User-tunable knob.
    Config,
    /// Transient / IO.
    Runtime,
}

impl LatticeError {
    /// Stable code string.
    pub fn code(&self) -> &'static str {
        match self {
            LatticeError::Degenerate(_) => "lattice.degenerate",
            LatticeError::BadParameter { .. } => "lattice.bad_parameter",
            LatticeError::TooManyPlacements { .. } => "lattice.too_many_placements",
            LatticeError::Io(_) => "lattice.io",
            LatticeError::Ron(_) => "lattice.ron",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            LatticeError::Degenerate(_) => ErrorCategory::Input,
            LatticeError::BadParameter { .. } | LatticeError::TooManyPlacements { .. } => {
                ErrorCategory::Config
            }
            LatticeError::Io(_) | LatticeError::Ron(_) => ErrorCategory::Runtime,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codes_and_cats() {
        assert_eq!(
            LatticeError::Degenerate("x".into()).code(),
            "lattice.degenerate"
        );
        assert_eq!(
            LatticeError::BadParameter {
                name: "n",
                reason: "zero".into()
            }
            .category(),
            ErrorCategory::Config
        );
    }
}
