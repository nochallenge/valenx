//! Error type for the acoustics calculators.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, AcousticsError>`]. The variants are deliberately coarse:
//! a caller almost always only wants to know *which named input* was out
//! of its physical domain (a non-positive pressure, a negative absolute
//! temperature, a zero room dimension) so it can point the user at the
//! offending field.
//!
//! Use [`AcousticsError::code`] for stable log / telemetry tagging and
//! [`AcousticsError::category`] (or [`AcousticsError::category_enum`]) to
//! bucket failures without matching every variant. The pattern mirrors
//! `valenx-astro`'s `AstroError` and `valenx-springs`'s `SpringsError`.

use thiserror::Error;

/// Shorthand for `Result<T, AcousticsError>`.
pub type Result<T> = core::result::Result<T, AcousticsError>;

/// Anything that can go wrong evaluating an acoustics formula.
///
/// This enum is `#[non_exhaustive]`: new variants may be added in future
/// releases without it being a breaking change, so downstream `match`
/// arms must include a wildcard.
#[derive(Debug, Error, Clone, PartialEq)]
#[non_exhaustive]
pub enum AcousticsError {
    /// A pressure argument was non-positive or non-finite. SPL is the
    /// logarithm of a pressure *ratio*, so both the sound pressure and
    /// the reference pressure must be strictly positive and finite.
    #[error("invalid pressure `{name}` = {value}: must be finite and > 0")]
    InvalidPressure {
        /// Logical parameter name (e.g. `"pressure"`, `"p_ref"`).
        name: &'static str,
        /// The offending value.
        value: f64,
    },

    /// An absolute (kelvin) temperature was at or below absolute zero, or
    /// was non-finite. The speed-of-sound model
    /// `c = 331.3*sqrt(1 + T/273.15)` requires `1 + T/273.15 > 0`, i.e.
    /// `T > -273.15` degrees Celsius (`T > 0` kelvin).
    #[error("invalid temperature `{name}` = {value}: must be finite and above absolute zero")]
    InvalidTemperature {
        /// Logical parameter name (e.g. `"temperature_c"`).
        name: &'static str,
        /// The offending value.
        value: f64,
    },

    /// A speed of sound was non-positive or non-finite. The Doppler shift
    /// and the room-mode formula both divide by `c`, so it must be a
    /// strictly positive, finite propagation speed.
    #[error("invalid speed of sound `{name}` = {value}: must be finite and > 0")]
    InvalidSpeedOfSound {
        /// Logical parameter name (e.g. `"speed_of_sound"`).
        name: &'static str,
        /// The offending value.
        value: f64,
    },

    /// A frequency argument was negative or non-finite. A source emits at
    /// a non-negative frequency.
    #[error("invalid frequency `{name}` = {value}: must be finite and >= 0")]
    InvalidFrequency {
        /// Logical parameter name (e.g. `"source_hz"`).
        name: &'static str,
        /// The offending value.
        value: f64,
    },

    /// A velocity argument was non-finite (NaN / infinity). Velocities may
    /// be negative (receding) but must be a real, finite speed.
    #[error("invalid velocity `{name}` = {value}: must be finite")]
    InvalidVelocity {
        /// Logical parameter name (e.g. `"source_velocity"`).
        name: &'static str,
        /// The offending value.
        value: f64,
    },

    /// The Doppler geometry is singular: the source approaches (or is
    /// modelled as approaching) the observer at or above the speed of
    /// sound, so `c - v_s <= 0` and the classical (subsonic) shift has no
    /// finite value. Carries the offending closing speed and `c`.
    #[error(
        "Doppler singularity: source closing speed {source_velocity} >= speed of sound {speed_of_sound} (c - v_s <= 0)"
    )]
    DopplerSingularity {
        /// The source velocity component toward the observer.
        source_velocity: f64,
        /// The speed of sound used.
        speed_of_sound: f64,
    },

    /// A room dimension (edge length) was non-positive or non-finite. Each
    /// of `Lx`, `Ly`, `Lz` must be a strictly positive, finite length.
    #[error("invalid room dimension `{name}` = {value}: must be finite and > 0")]
    InvalidDimension {
        /// Which edge (`"length_x"`, `"length_y"`, `"length_z"`).
        name: &'static str,
        /// The offending value.
        value: f64,
    },

    /// The requested room mode `(0, 0, 0)` is the trivial zero-frequency
    /// DC mode and is rejected: at least one modal index must be positive
    /// for a physical standing wave.
    #[error("mode (0, 0, 0) is the trivial DC mode; at least one index must be > 0")]
    TrivialMode,
}

/// Coarse category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this rather
/// than on every error variant.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// A caller-supplied argument is outside its physical domain.
    Input,
    /// A formula is singular at the requested point (e.g. a sonic Doppler
    /// source, or the trivial zero mode).
    Domain,
}

impl AcousticsError {
    /// Stable snake-cased error code suitable for log / telemetry tagging.
    /// Format: `"acoustics.<sub_id>"`. Codes never change across minor
    /// versions.
    pub fn code(&self) -> &'static str {
        match self {
            AcousticsError::InvalidPressure { .. } => "acoustics.invalid_pressure",
            AcousticsError::InvalidTemperature { .. } => "acoustics.invalid_temperature",
            AcousticsError::InvalidSpeedOfSound { .. } => "acoustics.invalid_speed_of_sound",
            AcousticsError::InvalidFrequency { .. } => "acoustics.invalid_frequency",
            AcousticsError::InvalidVelocity { .. } => "acoustics.invalid_velocity",
            AcousticsError::DopplerSingularity { .. } => "acoustics.doppler_singularity",
            AcousticsError::InvalidDimension { .. } => "acoustics.invalid_dimension",
            AcousticsError::TrivialMode => "acoustics.trivial_mode",
        }
    }

    /// Coarse category string — see [`ErrorCategory`].
    pub fn category(&self) -> &'static str {
        match self.category_enum() {
            ErrorCategory::Input => "input",
            ErrorCategory::Domain => "domain",
        }
    }

    /// Typed category enum (for callers that want to `match` instead of
    /// comparing the [`category`](Self::category) string).
    pub fn category_enum(&self) -> ErrorCategory {
        match self {
            AcousticsError::InvalidPressure { .. }
            | AcousticsError::InvalidTemperature { .. }
            | AcousticsError::InvalidSpeedOfSound { .. }
            | AcousticsError::InvalidFrequency { .. }
            | AcousticsError::InvalidVelocity { .. }
            | AcousticsError::InvalidDimension { .. } => ErrorCategory::Input,
            AcousticsError::DopplerSingularity { .. } | AcousticsError::TrivialMode => {
                ErrorCategory::Domain
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = AcousticsError::InvalidPressure {
            name: "p_ref",
            value: 0.0,
        };
        assert_eq!(err.code(), "acoustics.invalid_pressure");
        assert_eq!(err.category(), "input");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = AcousticsError::DopplerSingularity {
            source_velocity: 343.0,
            speed_of_sound: 343.0,
        };
        assert_eq!(err.code(), "acoustics.doppler_singularity");
        assert_eq!(err.category(), "domain");
        assert_eq!(err.category_enum(), ErrorCategory::Domain);

        assert_eq!(AcousticsError::TrivialMode.code(), "acoustics.trivial_mode");
        assert_eq!(AcousticsError::TrivialMode.category(), "domain");
    }

    #[test]
    fn display_is_informative() {
        let msg = AcousticsError::InvalidTemperature {
            name: "temperature_c",
            value: -300.0,
        }
        .to_string();
        assert!(msg.contains("temperature_c"), "got: {msg}");
        assert!(msg.contains("-300"), "got: {msg}");

        let msg = AcousticsError::InvalidDimension {
            name: "length_x",
            value: 0.0,
        }
        .to_string();
        assert!(msg.contains("length_x"), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(AcousticsError::TrivialMode);
        assert!(err.to_string().contains("DC mode"));
    }
}
