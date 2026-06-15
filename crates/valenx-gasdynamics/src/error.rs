//! Error taxonomy for `valenx-gasdynamics`.
//!
//! Every fallible public function returns [`Result<_, GasError>`]. The
//! variants are deliberately coarse — a compressible-flow caller usually
//! only needs to know whether an argument fell outside its physical
//! domain (a non-positive Mach number, a specific-heat ratio that is not
//! strictly greater than one) or whether a relation that is only defined
//! supersonically was asked to evaluate a subsonic state (the normal
//! shock requires `M1 >= 1`).
//!
//! Validated constructors ([`GasError::bad_mach`], [`GasError::bad_gamma`],
//! [`GasError::subsonic_shock`]) build the variants so call sites stay
//! terse and the messages stay uniform. [`GasError::code`] gives a stable
//! kebab-cased identifier for logging / telemetry; [`GasError::category`]
//! buckets failures without matching every variant.

use thiserror::Error;

/// Errors produced by `valenx-gasdynamics`.
///
/// The relations in this crate are pure algebra over two inputs — a Mach
/// number `M` and a specific-heat ratio `gamma` — so the only ways a call
/// can fail are a Mach number outside its physical domain, a `gamma`
/// outside `(1, infinity)`, or a normal-shock evaluation requested for a
/// subsonic upstream state.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum GasError {
    /// A Mach number was outside the domain required by the relation.
    ///
    /// Stagnation ratios and the area-Mach relation accept any
    /// `M >= 0`; the value here is rejected for being negative,
    /// non-finite (NaN / infinity), or — where the relation needs it —
    /// not strictly positive.
    #[error("invalid Mach number for `{context}`: {value} ({reason})")]
    BadMach {
        /// The offending Mach number, surfaced verbatim.
        value: f64,
        /// Short label for the relation that rejected it
        /// (e.g. `"area_mach_ratio"`, `"normal_shock"`).
        context: &'static str,
        /// Human-readable reason (e.g. `"must be finite and >= 0"`).
        reason: &'static str,
    },

    /// The specific-heat ratio `gamma` was outside its physical domain.
    ///
    /// For a calorically-perfect ideal gas `gamma = cp / cv` is strictly
    /// greater than one (air is `~1.4`, a monatomic gas `5/3`, a steam
    /// approximation `~1.33`). Values `<= 1` make the relations singular
    /// or sign-flipped, so they are rejected here.
    #[error("invalid specific-heat ratio gamma = {value}: {reason}")]
    BadGamma {
        /// The offending `gamma`, surfaced verbatim.
        value: f64,
        /// Human-readable reason (e.g. `"must be finite and > 1"`).
        reason: &'static str,
    },

    /// A normal-shock relation was asked to evaluate a subsonic (or
    /// non-finite) upstream Mach number.
    ///
    /// A stationary normal shock only exists for a supersonic approach
    /// flow, so the jump relations require `M1 >= 1`. At exactly
    /// `M1 = 1` every ratio is unity (the degenerate, infinitesimally
    /// weak shock), which is accepted.
    #[error("normal shock requires a supersonic upstream Mach number M1 >= 1, got M1 = {value}")]
    SubsonicShock {
        /// The offending upstream Mach number `M1`.
        value: f64,
    },
}

/// Coarse error category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this rather
/// than on every [`GasError`] variant.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// A caller-supplied argument was outside its physical domain.
    Input,
    /// A relation was evaluated outside the flow regime it is defined for
    /// (e.g. a subsonic normal shock).
    Domain,
}

impl GasError {
    /// Build a [`GasError::BadMach`] for an out-of-domain Mach number.
    ///
    /// `context` names the relation that rejected the value and `reason`
    /// states the constraint it violated.
    pub fn bad_mach(value: f64, context: &'static str, reason: &'static str) -> Self {
        GasError::BadMach {
            value,
            context,
            reason,
        }
    }

    /// Build a [`GasError::BadGamma`] for an out-of-domain specific-heat
    /// ratio.
    pub fn bad_gamma(value: f64, reason: &'static str) -> Self {
        GasError::BadGamma { value, reason }
    }

    /// Build a [`GasError::SubsonicShock`] for a subsonic upstream Mach
    /// number passed to a normal-shock relation.
    pub fn subsonic_shock(value: f64) -> Self {
        GasError::SubsonicShock { value }
    }

    /// Stable kebab-cased identifier suitable for log / telemetry
    /// tagging. Format: `"gasdynamics.<sub_id>"`. Codes never change
    /// across minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            GasError::BadMach { .. } => "gasdynamics.bad_mach",
            GasError::BadGamma { .. } => "gasdynamics.bad_gamma",
            GasError::SubsonicShock { .. } => "gasdynamics.subsonic_shock",
        }
    }

    /// Coarse [`ErrorCategory`] for the variant.
    pub fn category(&self) -> ErrorCategory {
        match self {
            GasError::BadMach { .. } | GasError::BadGamma { .. } => ErrorCategory::Input,
            GasError::SubsonicShock { .. } => ErrorCategory::Domain,
        }
    }
}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, GasError>;

/// Validate a specific-heat ratio: it must be finite and strictly greater
/// than one. Returns `gamma` unchanged on success.
///
/// Shared by every relation so the `gamma` domain is enforced in exactly
/// one place.
pub(crate) fn check_gamma(gamma: f64) -> Result<f64> {
    if !gamma.is_finite() || gamma <= 1.0 {
        return Err(GasError::bad_gamma(gamma, "must be finite and > 1"));
    }
    Ok(gamma)
}

/// Validate a Mach number that must be finite and non-negative
/// (`M >= 0`). Returns `mach` unchanged on success.
pub(crate) fn check_mach_nonneg(mach: f64, context: &'static str) -> Result<f64> {
    if !mach.is_finite() || mach < 0.0 {
        return Err(GasError::bad_mach(mach, context, "must be finite and >= 0"));
    }
    Ok(mach)
}

/// Validate a Mach number that must be finite and strictly positive
/// (`M > 0`, e.g. the area-Mach relation, which divides by `M`). Returns
/// `mach` unchanged on success.
pub(crate) fn check_mach_pos(mach: f64, context: &'static str) -> Result<f64> {
    if !mach.is_finite() || mach <= 0.0 {
        return Err(GasError::bad_mach(mach, context, "must be finite and > 0"));
    }
    Ok(mach)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = GasError::bad_mach(-1.0, "area_mach_ratio", "must be finite and > 0");
        assert_eq!(err.code(), "gasdynamics.bad_mach");
        assert_eq!(err.category(), ErrorCategory::Input);

        let err = GasError::bad_gamma(0.5, "must be finite and > 1");
        assert_eq!(err.code(), "gasdynamics.bad_gamma");
        assert_eq!(err.category(), ErrorCategory::Input);

        let err = GasError::subsonic_shock(0.8);
        assert_eq!(err.code(), "gasdynamics.subsonic_shock");
        assert_eq!(err.category(), ErrorCategory::Domain);
    }

    #[test]
    fn display_is_informative() {
        let msg = GasError::bad_mach(-2.0, "normal_shock", "must be finite and >= 0").to_string();
        assert!(msg.contains("normal_shock"), "got: {msg}");
        assert!(msg.contains("-2"), "got: {msg}");

        let msg = GasError::subsonic_shock(0.5).to_string();
        assert!(msg.contains("0.5"), "got: {msg}");
        assert!(msg.contains("supersonic"), "got: {msg}");

        let msg = GasError::bad_gamma(1.0, "must be finite and > 1").to_string();
        assert!(msg.contains('1'), "got: {msg}");
    }

    #[test]
    fn check_gamma_accepts_air_and_rejects_unit() {
        assert_eq!(check_gamma(1.4).unwrap(), 1.4);
        assert!(check_gamma(1.0).is_err());
        assert!(check_gamma(0.9).is_err());
        assert!(check_gamma(f64::NAN).is_err());
        assert!(check_gamma(f64::INFINITY).is_err());
    }

    #[test]
    fn check_mach_domains() {
        assert_eq!(check_mach_nonneg(0.0, "ctx").unwrap(), 0.0);
        assert!(check_mach_nonneg(-0.1, "ctx").is_err());
        assert!(check_mach_nonneg(f64::NAN, "ctx").is_err());

        assert_eq!(check_mach_pos(0.5, "ctx").unwrap(), 0.5);
        assert!(check_mach_pos(0.0, "ctx").is_err());
        assert!(check_mach_pos(f64::INFINITY, "ctx").is_err());
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(GasError::subsonic_shock(0.3));
        assert!(err.to_string().contains("0.3"));
    }
}
