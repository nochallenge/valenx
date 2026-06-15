//! Error taxonomy for `valenx-rootfind`.
//!
//! Every fallible public function in this crate returns
//! [`Result<T>`](crate::Result), i.e. `Result<T, RootError>`. The
//! variants are intentionally coarse — a root-finding caller usually
//! only cares about four failure modes:
//!
//! 1. Did the caller pass a nonsense argument — a non-positive
//!    tolerance, a zero iteration budget, two equal secant seeds
//!    ([`RootError::Invalid`])?
//! 2. Is the supplied interval not a valid bracket — does `f` fail to
//!    change sign across `[a, b]` so bisection cannot start
//!    ([`RootError::NoSignChange`])?
//! 3. Did an open method stall because the local slope (Newton's
//!    derivative, or the secant's finite difference) is effectively
//!    zero ([`RootError::ZeroDerivative`])?
//! 4. Did `f` (or `f'`) return a non-finite value (`NaN` / `±∞`)
//!    ([`RootError::NotFinite`]), or did the method exhaust its
//!    iteration budget without meeting the tolerance
//!    ([`RootError::MaxIterations`])?
//!
//! Use [`RootError::code`] for stable log / telemetry tagging and
//! [`RootError::category`] to bucket failures into Input / Convergence
//! without matching every variant. The shape mirrors `valenx-popgen`'s
//! `PopgenError` and the other Valenx numerical crates.

use thiserror::Error;

/// Errors produced by the `valenx-rootfind` solvers.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum RootError {
    /// Caller passed an argument the solver cannot accept: a tolerance
    /// that is not strictly positive (or not finite), a zero iteration
    /// budget, or two equal seed points for the secant method. A
    /// property of the *call*, not of the function being solved.
    #[error("invalid `{what}`: {reason}")]
    Invalid {
        /// Logical parameter name (e.g. `"tol"`, `"max_iter"`, `"x0/x1"`).
        what: &'static str,
        /// Human-readable reason.
        reason: String,
    },

    /// The interval handed to a bracketing method is not a valid
    /// bracket: `f(a)` and `f(b)` have the same sign (and neither is
    /// zero), so no root is guaranteed to lie between them and bisection
    /// cannot begin.
    #[error("no sign change on the bracket: f({a}) = {fa} and f({b}) = {fb} have the same sign")]
    NoSignChange {
        /// Lower bracket endpoint.
        a: f64,
        /// Upper bracket endpoint.
        b: f64,
        /// Function value at `a`.
        fa: f64,
        /// Function value at `b`.
        fb: f64,
    },

    /// An open method stalled: the slope it divides by — Newton's
    /// supplied derivative `f'(x)`, or the secant's finite difference
    /// `f(x_n) - f(x_{n-1})` — is smaller in magnitude than the guard
    /// threshold, so the update step is undefined or would explode.
    #[error("derivative effectively zero at x = {x} (|slope| = {slope} below {threshold})")]
    ZeroDerivative {
        /// Iterate at which the slope collapsed.
        x: f64,
        /// Magnitude of the offending slope.
        slope: f64,
        /// Guard threshold the slope fell below.
        threshold: f64,
    },

    /// `f` (or the supplied derivative `f'`) returned a value that is
    /// not finite (`NaN`, `+∞` or `-∞`) at the named iterate. The
    /// iteration cannot continue from a non-finite state.
    #[error("non-finite {what} at x = {x}")]
    NotFinite {
        /// What was non-finite (e.g. `"f(x)"`, `"f'(x)"`, `"iterate"`).
        what: &'static str,
        /// Iterate at which the non-finite value appeared.
        x: f64,
    },

    /// The method ran for `max_iter` iterations without the residual /
    /// step meeting the tolerance. `last` is the best iterate reached
    /// and `residual` its `|f(last)|`, so the caller can decide whether
    /// the near-miss is good enough.
    #[error(
        "did not converge to tol {tol} within {max_iter} iterations \
         (last x = {last}, |f| = {residual})"
    )]
    MaxIterations {
        /// Iteration budget that was exhausted.
        max_iter: usize,
        /// Tolerance that was being targeted.
        tol: f64,
        /// Best iterate reached before giving up.
        last: f64,
        /// Residual `|f(last)|` at that iterate.
        residual: f64,
    },
}

/// Coarse category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this
/// rather than on the full set of error variants.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// The caller supplied a bad argument or a non-bracketing interval
    /// before any iteration ran.
    Input,
    /// The iteration started but could not converge — a vanishing
    /// slope, a non-finite evaluation, or an exhausted budget.
    Convergence,
}

impl RootError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"rootfind.<sub_id>"`. Codes never change across
    /// minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            RootError::Invalid { .. } => "rootfind.invalid",
            RootError::NoSignChange { .. } => "rootfind.no_sign_change",
            RootError::ZeroDerivative { .. } => "rootfind.zero_derivative",
            RootError::NotFinite { .. } => "rootfind.not_finite",
            RootError::MaxIterations { .. } => "rootfind.max_iterations",
        }
    }

    /// Coarse category string — see [`ErrorCategory`].
    pub fn category(&self) -> &'static str {
        match self {
            RootError::Invalid { .. } | RootError::NoSignChange { .. } => "input",
            RootError::ZeroDerivative { .. }
            | RootError::NotFinite { .. }
            | RootError::MaxIterations { .. } => "convergence",
        }
    }

    /// Typed category enum (for callers that want to `match` instead of
    /// comparing the [`category`](Self::category) string).
    pub fn category_enum(&self) -> ErrorCategory {
        match self {
            RootError::Invalid { .. } | RootError::NoSignChange { .. } => ErrorCategory::Input,
            RootError::ZeroDerivative { .. }
            | RootError::NotFinite { .. }
            | RootError::MaxIterations { .. } => ErrorCategory::Convergence,
        }
    }

    /// Convenience constructor for [`RootError::Invalid`].
    pub fn invalid(what: &'static str, reason: impl Into<String>) -> Self {
        RootError::Invalid {
            what,
            reason: reason.into(),
        }
    }
}

/// Crate-wide result alias: `Result<T, RootError>`.
pub type Result<T> = std::result::Result<T, RootError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = RootError::invalid("tol", "must be positive");
        assert_eq!(err.code(), "rootfind.invalid");
        assert_eq!(err.category(), "input");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = RootError::NoSignChange {
            a: 0.0,
            b: 1.0,
            fa: -1.0,
            fb: -2.0,
        };
        assert_eq!(err.code(), "rootfind.no_sign_change");
        assert_eq!(err.category(), "input");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = RootError::ZeroDerivative {
            x: 0.0,
            slope: 0.0,
            threshold: 1e-12,
        };
        assert_eq!(err.code(), "rootfind.zero_derivative");
        assert_eq!(err.category(), "convergence");
        assert_eq!(err.category_enum(), ErrorCategory::Convergence);

        let err = RootError::NotFinite {
            what: "f(x)",
            x: 2.0,
        };
        assert_eq!(err.code(), "rootfind.not_finite");
        assert_eq!(err.category(), "convergence");

        let err = RootError::MaxIterations {
            max_iter: 10,
            tol: 1e-9,
            last: 1.4,
            residual: 0.2,
        };
        assert_eq!(err.code(), "rootfind.max_iterations");
        assert_eq!(err.category(), "convergence");
        assert_eq!(err.category_enum(), ErrorCategory::Convergence);
    }

    #[test]
    fn display_is_informative() {
        let msg = RootError::NoSignChange {
            a: 0.0,
            b: 1.0,
            fa: -1.0,
            fb: -2.0,
        }
        .to_string();
        assert!(msg.contains("sign change"), "got: {msg}");

        let msg = RootError::MaxIterations {
            max_iter: 7,
            tol: 1e-9,
            last: 1.4,
            residual: 0.2,
        }
        .to_string();
        assert!(msg.contains('7'), "got: {msg}");
        assert!(msg.contains("converge"), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> =
            Box::new(RootError::invalid("max_iter", "must be > 0"));
        assert!(err.to_string().contains("max_iter"));
    }
}
