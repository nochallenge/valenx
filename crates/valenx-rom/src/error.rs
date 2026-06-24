//! Error taxonomy for reduced-order modelling.
//!
//! Every fallible routine in this crate returns [`RomError`] rather than
//! panicking or emitting a plausible-but-wrong number. The variants are
//! deliberately specific so a caller (or a test) can branch on *why* a fit
//! refused to run — empty data, a rank-deficient snapshot set, a dimension
//! mismatch, a non-finite entry, or an ill-posed reduced problem.

use thiserror::Error;

/// Errors raised while fitting or evaluating a reduced-order model.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum RomError {
    /// The snapshot matrix had zero rows or zero columns.
    #[error("empty snapshot matrix ({rows}x{cols}); need at least 1 row and 1 column")]
    Empty {
        /// Number of rows (state dimension) supplied.
        rows: usize,
        /// Number of columns (time samples) supplied.
        cols: usize,
    },

    /// A finite-valued container held a `NaN` or an infinity.
    #[error("non-finite value in {what}")]
    NonFinite {
        /// Which container was non-finite.
        what: &'static str,
    },

    /// Two operands had incompatible shapes.
    #[error("dimension mismatch in {what}: expected {expected}, got {got}")]
    DimensionMismatch {
        /// Where the mismatch occurred.
        what: &'static str,
        /// The required size.
        expected: usize,
        /// The supplied size.
        got: usize,
    },

    /// A requested truncation rank was zero or exceeded the data rank.
    #[error("invalid rank {requested}: must be in 1..={max}")]
    InvalidRank {
        /// The requested rank.
        requested: usize,
        /// The largest admissible rank (min dimension of the data).
        max: usize,
    },

    /// The data (or a reduced operand) was numerically rank-deficient: every
    /// singular value was at or below the working tolerance, so there is no
    /// well-defined basis / inverse to return.
    #[error("rank-deficient {what}: no singular value exceeds tolerance {tol:e}")]
    RankDeficient {
        /// What was rank-deficient.
        what: &'static str,
        /// The tolerance the singular values failed to exceed.
        tol: f64,
    },

    /// An energy tolerance was outside the half-open `(0, 1]` interval.
    #[error("energy tolerance {value} must be in (0, 1]")]
    BadEnergyTol {
        /// The offending tolerance.
        value: f64,
    },

    /// Fewer time samples than required for the requested operation (e.g. DMD
    /// needs at least two snapshots to form the `X / X'` pair).
    #[error("need at least {needed} time samples for {what}, got {got}")]
    NotEnoughSamples {
        /// What needed the samples.
        what: &'static str,
        /// How many samples are required.
        needed: usize,
        /// How many were supplied.
        got: usize,
    },

    /// A required numerical factorisation failed to converge.
    #[error("{what} factorisation failed to converge")]
    NotConverged {
        /// Which factorisation failed.
        what: &'static str,
    },

    /// A time step `dt` was non-positive or non-finite.
    #[error("dt {value} must be finite and > 0")]
    BadTimeStep {
        /// The offending time step.
        value: f64,
    },
}

impl RomError {
    /// A short, stable machine-readable code for this error.
    ///
    /// Useful for matching in tests and for logging without formatting the
    /// human-readable message.
    pub fn code(&self) -> &'static str {
        match self {
            RomError::Empty { .. } => "empty",
            RomError::NonFinite { .. } => "non_finite",
            RomError::DimensionMismatch { .. } => "dimension_mismatch",
            RomError::InvalidRank { .. } => "invalid_rank",
            RomError::RankDeficient { .. } => "rank_deficient",
            RomError::BadEnergyTol { .. } => "bad_energy_tol",
            RomError::NotEnoughSamples { .. } => "not_enough_samples",
            RomError::NotConverged { .. } => "not_converged",
            RomError::BadTimeStep { .. } => "bad_time_step",
        }
    }
}
