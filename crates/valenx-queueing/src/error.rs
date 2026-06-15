//! Error taxonomy for `valenx-queueing`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, QueueingError>`]. The variants are intentionally coarse —
//! a queueing-theory caller usually only cares about three things:
//!
//! 1. Did the caller pass a non-positive or non-finite rate — a zero or
//!    negative arrival rate `lambda`, a zero or negative service rate
//!    `mu`, a `NaN`/infinite input ([`QueueingError::Invalid`])?
//! 2. Is the queue unstable — the offered load `rho = lambda / mu` is
//!    at or above 1, so no finite steady state exists and the closed
//!    forms diverge ([`QueueingError::Unstable`])?
//! 3. Was an out-of-domain index requested — a negative-count state
//!    probability `P(n)` ([`QueueingError::Domain`])?
//!
//! Use [`QueueingError::code`] for stable log / telemetry tagging and
//! [`QueueingError::category`] to bucket failures into Input / Stability
//! / Domain without matching every variant. The pattern mirrors
//! `valenx-springs`'s `SpringsError` and `valenx-popgen`'s `PopgenError`.

use thiserror::Error;

/// Errors produced by `valenx-queueing`.
///
/// Construct these through the validated helpers
/// ([`QueueingError::invalid`], [`QueueingError::unstable`],
/// [`QueueingError::domain`]) rather than building the variants by hand,
/// so the `Display` strings stay uniform.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum QueueingError {
    /// Caller passed a rate the model cannot accept: a non-positive
    /// arrival rate `lambda`, a non-positive service rate `mu`, or any
    /// non-finite (`NaN`/`±inf`) value. A property of the *call*, not
    /// of the queue's stability.
    #[error("invalid `{what}`: {reason}")]
    Invalid {
        /// Logical parameter name (e.g. `"lambda"`, `"mu"`).
        what: &'static str,
        /// Human-readable reason.
        reason: String,
    },

    /// The queue is unstable: the traffic intensity
    /// `rho = lambda / mu` is at or above 1, so arrivals meet or exceed
    /// the server's capacity, the backlog grows without bound, and the
    /// steady-state metrics (`L`, `Lq`, `W`, `Wq`) diverge to infinity.
    /// Stability requires `rho < 1`, i.e. `lambda < mu`.
    #[error(
        "unstable queue: rho = lambda / mu = {rho} >= 1 (lambda = {lambda}, mu = {mu}); \
         the M/M/1 steady state requires rho < 1"
    )]
    Unstable {
        /// Arrival rate that produced the instability.
        lambda: f64,
        /// Service rate that produced the instability.
        mu: f64,
        /// The offending traffic intensity `lambda / mu` (`>= 1`).
        rho: f64,
    },

    /// An out-of-domain quantity was requested — currently only a
    /// stationary state probability `P(n)` for an index `n` that cannot
    /// be a valid customer count. A property of the *query*, not of the
    /// queue parameters.
    #[error("out-of-domain {what}: {reason}")]
    Domain {
        /// Logical quantity requested (e.g. `"state_probability"`).
        what: &'static str,
        /// Human-readable reason.
        reason: String,
    },
}

/// Coarse category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this rather
/// than on the individual error variants.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// A supplied rate is non-positive or non-finite (bad argument).
    Input,
    /// The queue is saturated (`rho >= 1`); no finite steady state.
    Stability,
    /// A requested quantity is outside its valid domain.
    Domain,
}

impl QueueingError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"queueing.<sub_id>"`. Codes never change across
    /// minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            QueueingError::Invalid { .. } => "queueing.invalid",
            QueueingError::Unstable { .. } => "queueing.unstable",
            QueueingError::Domain { .. } => "queueing.domain",
        }
    }

    /// Coarse category string — see [`ErrorCategory`].
    pub fn category(&self) -> &'static str {
        match self {
            QueueingError::Invalid { .. } => "input",
            QueueingError::Unstable { .. } => "stability",
            QueueingError::Domain { .. } => "domain",
        }
    }

    /// Typed category enum (for callers that want to `match` instead of
    /// comparing the [`category`](Self::category) string).
    pub fn category_enum(&self) -> ErrorCategory {
        match self {
            QueueingError::Invalid { .. } => ErrorCategory::Input,
            QueueingError::Unstable { .. } => ErrorCategory::Stability,
            QueueingError::Domain { .. } => ErrorCategory::Domain,
        }
    }

    /// Validated constructor for [`QueueingError::Invalid`].
    pub fn invalid(what: &'static str, reason: impl Into<String>) -> Self {
        QueueingError::Invalid {
            what,
            reason: reason.into(),
        }
    }

    /// Validated constructor for [`QueueingError::Unstable`].
    ///
    /// Stores the offending `lambda`, `mu` and their ratio so the
    /// message can report all three.
    pub fn unstable(lambda: f64, mu: f64, rho: f64) -> Self {
        QueueingError::Unstable { lambda, mu, rho }
    }

    /// Validated constructor for [`QueueingError::Domain`].
    pub fn domain(what: &'static str, reason: impl Into<String>) -> Self {
        QueueingError::Domain {
            what,
            reason: reason.into(),
        }
    }
}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, QueueingError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = QueueingError::invalid("lambda", "must be positive");
        assert_eq!(err.code(), "queueing.invalid");
        assert_eq!(err.category(), "input");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = QueueingError::unstable(5.0, 4.0, 1.25);
        assert_eq!(err.code(), "queueing.unstable");
        assert_eq!(err.category(), "stability");
        assert_eq!(err.category_enum(), ErrorCategory::Stability);

        let err = QueueingError::domain("state_probability", "negative index");
        assert_eq!(err.code(), "queueing.domain");
        assert_eq!(err.category(), "domain");
        assert_eq!(err.category_enum(), ErrorCategory::Domain);
    }

    #[test]
    fn display_is_informative() {
        let msg = QueueingError::invalid("mu", "must be finite").to_string();
        assert!(msg.contains("mu"), "got: {msg}");
        assert!(msg.contains("finite"), "got: {msg}");

        // The unstable message names lambda, mu and rho.
        let msg = QueueingError::unstable(5.0, 4.0, 1.25).to_string();
        assert!(msg.contains("1.25"), "got: {msg}");
        assert!(msg.contains('5') && msg.contains('4'), "got: {msg}");
        assert!(msg.contains("rho"), "got: {msg}");

        let msg = QueueingError::domain("state_probability", "n out of range").to_string();
        assert!(msg.contains("state_probability"), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(QueueingError::invalid("x", "y"));
        assert!(err.to_string().contains('x'));
    }
}
