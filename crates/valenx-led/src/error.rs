//! Error taxonomy for the LED circuit calculator.
//!
//! Every fallible constructor in this crate validates its inputs and
//! returns a [`LedError`] on violation, so an out-of-domain circuit can
//! never be represented. The variants are deliberately fine-grained so a
//! caller can branch on the precise failure (for example, distinguishing
//! "supply too low" from "non-positive current") without string parsing.

use thiserror::Error;

/// Errors raised when building or evaluating an LED circuit.
///
/// Construction helpers ([`crate::circuit::LedCircuit::new`],
/// [`crate::circuit::LedString::new`]) reject physically meaningless
/// inputs up front and surface the reason through one of these variants.
#[derive(Debug, Clone, PartialEq, Error)]
#[non_exhaustive]
pub enum LedError {
    /// A quantity that must be strictly positive was zero or negative.
    ///
    /// Carries the offending parameter `name` (a stable, machine-readable
    /// identifier such as `"current_a"`) and the `value` that was supplied.
    #[error("parameter `{name}` must be strictly positive, got {value}")]
    NonPositive {
        /// Stable identifier of the rejected parameter.
        name: &'static str,
        /// The supplied (invalid) value.
        value: f64,
    },

    /// A quantity that must be non-negative was negative.
    ///
    /// Used for values for which exactly zero is physically admissible
    /// (for example a forward voltage of `0.0`) but a negative magnitude
    /// is not.
    #[error("parameter `{name}` must be non-negative, got {value}")]
    Negative {
        /// Stable identifier of the rejected parameter.
        name: &'static str,
        /// The supplied (invalid) value.
        value: f64,
    },

    /// A supplied floating-point quantity was NaN or infinite.
    ///
    /// All inputs must be finite real numbers; this guards against `NaN`
    /// or `±∞` propagating silently through the closed-form expressions.
    #[error("parameter `{name}` must be finite, got {value}")]
    NotFinite {
        /// Stable identifier of the rejected parameter.
        name: &'static str,
        /// The supplied (non-finite) value.
        value: f64,
    },

    /// The supply voltage does not exceed the total forward voltage drop.
    ///
    /// The current-limiting resistor can only develop a voltage across it,
    /// and therefore conduct current, when `supply_v > total_forward_v`.
    /// At or below that threshold the LED(s) cannot be forward biased into
    /// conduction and no operating point exists, so the calculator refuses
    /// the circuit rather than reporting a zero or negative resistor.
    #[error(
        "supply voltage {supply_v} V must exceed total forward voltage \
         {forward_v} V (headroom {headroom} V must be positive)"
    )]
    InsufficientHeadroom {
        /// The supply (source) voltage in volts.
        supply_v: f64,
        /// The summed LED forward voltage in volts.
        forward_v: f64,
        /// `supply_v - forward_v`; non-positive here, hence the rejection.
        headroom: f64,
    },

    /// An LED string was constructed with zero LEDs.
    ///
    /// A series string must contain at least one LED for the forward-voltage
    /// sum and the operating point to be defined.
    #[error("LED string must contain at least one LED, got {count}")]
    EmptyString {
        /// The requested (zero) LED count.
        count: usize,
    },
}

impl LedError {
    /// Stable, kebab-cased identifier for the error variant.
    ///
    /// The returned string is part of the crate's public contract and is
    /// suitable for logging, metrics labels, or mapping to UI messages; it
    /// will not change for an existing variant across patch releases.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_led::error::LedError;
    ///
    /// let err = LedError::EmptyString { count: 0 };
    /// assert_eq!(err.code(), "led.empty-string");
    /// ```
    pub fn code(&self) -> &'static str {
        match self {
            LedError::NonPositive { .. } => "led.non-positive",
            LedError::Negative { .. } => "led.negative",
            LedError::NotFinite { .. } => "led.not-finite",
            LedError::InsufficientHeadroom { .. } => "led.insufficient-headroom",
            LedError::EmptyString { .. } => "led.empty-string",
        }
    }
}

/// Validate that `value` is a finite, strictly positive number.
///
/// Returns the value unchanged on success, or the appropriate
/// [`LedError`] (`NotFinite` then `NonPositive`) tagged with `name`.
pub(crate) fn require_positive(name: &'static str, value: f64) -> Result<f64, LedError> {
    if !value.is_finite() {
        return Err(LedError::NotFinite { name, value });
    }
    if value <= 0.0 {
        return Err(LedError::NonPositive { name, value });
    }
    Ok(value)
}

/// Validate that `value` is a finite, non-negative number.
///
/// Returns the value unchanged on success, or the appropriate
/// [`LedError`] (`NotFinite` then `Negative`) tagged with `name`.
pub(crate) fn require_non_negative(name: &'static str, value: f64) -> Result<f64, LedError> {
    if !value.is_finite() {
        return Err(LedError::NotFinite { name, value });
    }
    if value < 0.0 {
        return Err(LedError::Negative { name, value });
    }
    Ok(value)
}
