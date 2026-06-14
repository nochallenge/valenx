//! Error taxonomy for `valenx-regnet`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, RegnetError>`](crate::Result). The variants separate the
//! failure modes a caller cares about:
//!
//! 1. A model parameter was physically meaningless — a non-positive Hill
//!    threshold or coefficient ([`RegnetError::InvalidHill`]), a negative
//!    production or degradation rate ([`RegnetError::InvalidRate`]).
//! 2. The network was assembled inconsistently — a state / parameter
//!    length that does not match the declared gene count
//!    ([`RegnetError::DimensionMismatch`]), or a regulator index that
//!    points past the last gene ([`RegnetError::GeneIndexOutOfRange`]).
//! 3. The integrator was asked to do something impossible — a
//!    non-positive timestep or zero steps ([`RegnetError::InvalidStep`]).

use thiserror::Error;

/// Errors produced by `valenx-regnet`.
///
/// Derives [`thiserror::Error`]; each variant carries a human-readable
/// `Display` message via its `#[error(...)]` attribute.
///
/// Marked `#[non_exhaustive]`: more failure modes may be added as the
/// model grows, so downstream matches must include a wildcard arm.
#[derive(Debug, Clone, PartialEq, Error)]
#[non_exhaustive]
pub enum RegnetError {
    /// A Hill-kinetics parameter was out of range. Both the threshold
    /// `k` and the coefficient `n` must be strictly positive: `k` is a
    /// concentration scale (so `k > 0`) and `n` is a cooperativity
    /// exponent (so `n > 0`).
    #[error("invalid Hill parameter ({what}): got {value}, but it must be strictly positive")]
    InvalidHill {
        /// Which parameter was bad — `"threshold k"` or `"coefficient n"`.
        what: &'static str,
        /// The offending value.
        value: f64,
    },

    /// A production or degradation rate was negative. Rates are
    /// non-negative by construction (a negative production or decay rate
    /// is not physical), so the validated constructors reject them.
    #[error("invalid {what} rate for gene {gene}: got {value}, but it must be >= 0")]
    InvalidRate {
        /// Which rate was bad — `"production"` or `"degradation"`.
        what: &'static str,
        /// Index of the gene whose rate was rejected.
        gene: usize,
        /// The offending value.
        value: f64,
    },

    /// A vector length did not match the network's gene count. `expected`
    /// is the declared number of genes `N`; `found` is the length of the
    /// state or parameter vector that was supplied.
    #[error("dimension mismatch ({what}): expected length {expected}, found {found}")]
    DimensionMismatch {
        /// What was being checked — e.g. `"initial state"`.
        what: &'static str,
        /// The declared gene count `N`.
        expected: usize,
        /// The length actually supplied.
        found: usize,
    },

    /// A regulator referenced a gene index `index` that is past the last
    /// gene of a network with `count` genes (valid indices are
    /// `0..count`).
    #[error("regulator references gene {index}, but the network has only {count} genes")]
    GeneIndexOutOfRange {
        /// The out-of-range gene index.
        index: usize,
        /// The number of genes in the network.
        count: usize,
    },

    /// The integrator was given a non-usable step. The timestep `dt` must
    /// be strictly positive and the requested number of steps must be
    /// non-zero.
    #[error("invalid integration step: {reason}")]
    InvalidStep {
        /// Human-readable reason (e.g. `"dt must be > 0, got 0"`).
        reason: String,
    },
}

impl RegnetError {
    /// Convenience constructor for [`RegnetError::InvalidStep`].
    pub fn invalid_step(reason: impl Into<String>) -> Self {
        RegnetError::InvalidStep {
            reason: reason.into(),
        }
    }
}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, RegnetError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn invalid_hill_display_names_parameter_and_value() {
        let e = RegnetError::InvalidHill {
            what: "threshold k",
            value: -1.0,
        };
        let msg = e.to_string();
        assert!(msg.contains("threshold k"), "got: {msg}");
        assert!(msg.contains("-1"), "got: {msg}");
    }

    #[test]
    fn invalid_rate_display_names_gene() {
        let e = RegnetError::InvalidRate {
            what: "degradation",
            gene: 3,
            value: -0.5,
        };
        let msg = e.to_string();
        assert!(msg.contains("degradation"), "got: {msg}");
        assert!(msg.contains('3'), "got: {msg}");
    }

    #[test]
    fn dimension_mismatch_display() {
        let e = RegnetError::DimensionMismatch {
            what: "initial state",
            expected: 2,
            found: 5,
        };
        let msg = e.to_string();
        assert!(msg.contains('2'), "got: {msg}");
        assert!(msg.contains('5'), "got: {msg}");
    }

    #[test]
    fn gene_index_out_of_range_display() {
        let e = RegnetError::GeneIndexOutOfRange { index: 7, count: 3 };
        let msg = e.to_string();
        assert!(msg.contains('7'), "got: {msg}");
        assert!(msg.contains('3'), "got: {msg}");
    }

    #[test]
    fn invalid_step_constructor_and_display() {
        let e = RegnetError::invalid_step("dt must be > 0, got 0");
        assert!(e.to_string().contains("dt must be > 0"), "got: {e}");
    }

    #[test]
    fn error_is_a_std_error_trait_object() {
        let err: Box<dyn std::error::Error> =
            Box::new(RegnetError::GeneIndexOutOfRange { index: 1, count: 0 });
        assert!(err.to_string().contains('1'));
    }
}
