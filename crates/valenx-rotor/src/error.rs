//! Rotor / BEMT error taxonomy.
//!
//! Every fallible constructor and entry point in this crate returns
//! [`RotorError`]. The variants distinguish a plainly invalid physical
//! input (a non-positive radius, a non-finite density) from a structurally
//! inconsistent geometry (hub radius not below the tip, an empty or
//! out-of-order station table) and from a numerical failure of the
//! blade-element-momentum solver (the inflow-angle root-finder failed to
//! bracket or converge within the iteration cap). A solver failure is a
//! returned error, never a silently-emitted NaN or infinity.

use thiserror::Error;

/// Errors raised by rotor blade-element-momentum calculations.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum RotorError {
    /// A scalar parameter that must be finite and strictly positive was not.
    ///
    /// `name` is the offending parameter; `value` is what was supplied.
    #[error("`{name}` must be finite and positive, got {value}")]
    NonPositive {
        /// Parameter name (e.g. `"tip_radius"`, `"chord"`, `"air_density"`).
        name: &'static str,
        /// The offending value.
        value: f64,
    },

    /// A value that must merely be finite (it may be zero or negative, e.g.
    /// the freestream speed or a twist angle) was `NaN` or infinite.
    #[error("`{name}` must be finite, got {value}")]
    NotFinite {
        /// Parameter name.
        name: &'static str,
        /// The offending value.
        value: f64,
    },

    /// The blade count was zero.
    #[error("blade count must be at least 1")]
    NoBlades,

    /// The hub radius was not strictly below the tip radius.
    ///
    /// A blade must have positive span: `0 <= hub_radius < tip_radius`.
    #[error("hub radius {hub} must be < tip radius {tip}")]
    HubNotBelowTip {
        /// The supplied hub radius (m).
        hub: f64,
        /// The supplied tip radius (m).
        tip: f64,
    },

    /// The radial station table was empty or had fewer than two stations.
    ///
    /// Integrating loads over the span by the trapezoid rule needs at
    /// least two stations.
    #[error("need at least 2 radial stations, got {count}")]
    TooFewStations {
        /// Number of stations supplied.
        count: usize,
    },

    /// The radial stations were not strictly increasing in radius, or fell
    /// outside `[hub_radius, tip_radius]`.
    #[error("radial stations must be strictly increasing within [hub, tip]: {reason}")]
    BadStations {
        /// Human-readable description of the violated ordering / bound.
        reason: String,
    },

    /// A caller-supplied tabulated airfoil polar was malformed.
    ///
    /// A usable table needs at least two angle-of-attack samples that are
    /// strictly increasing in `alpha`, all finite.
    #[error("invalid airfoil polar table: {reason}")]
    BadPolar {
        /// Human-readable description of the problem.
        reason: String,
    },

    /// The inflow-angle root-finder could not solve a blade element.
    ///
    /// The momentum/blade-element consistency residual
    /// `tan(phi) - V(1-a)/(Omega r (1+a'))` either could not be bracketed
    /// on `(eps, pi/2 - eps)` or did not converge within the iteration cap.
    /// Reported per element so the caller can see which station failed.
    #[error("BEMT did not converge at station {station} (r = {radius} m): {reason}")]
    NoConvergence {
        /// Index of the offending radial station.
        station: usize,
        /// Radius of the offending station (m).
        radius: f64,
        /// Why the solve failed (no bracket / iteration cap hit).
        reason: &'static str,
    },
}

/// Coarse category for a [`RotorError`], for UI grouping / logging.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// The caller supplied an invalid physical input value.
    Input,
    /// The geometry / table structure is inconsistent.
    Config,
    /// The numerical solver failed to converge.
    Algorithm,
}

impl RotorError {
    /// Stable kebab-cased identifier, suitable for matching in logs or
    /// machine-readable error reporting.
    ///
    /// ```
    /// use valenx_rotor::RotorError;
    /// let e = RotorError::NoBlades;
    /// assert_eq!(e.code(), "rotor.no_blades");
    /// ```
    pub fn code(&self) -> &'static str {
        match self {
            RotorError::NonPositive { .. } => "rotor.non_positive",
            RotorError::NotFinite { .. } => "rotor.not_finite",
            RotorError::NoBlades => "rotor.no_blades",
            RotorError::HubNotBelowTip { .. } => "rotor.hub_not_below_tip",
            RotorError::TooFewStations { .. } => "rotor.too_few_stations",
            RotorError::BadStations { .. } => "rotor.bad_stations",
            RotorError::BadPolar { .. } => "rotor.bad_polar",
            RotorError::NoConvergence { .. } => "rotor.no_convergence",
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            RotorError::NonPositive { .. } | RotorError::NotFinite { .. } => ErrorCategory::Input,
            RotorError::NoBlades
            | RotorError::HubNotBelowTip { .. }
            | RotorError::TooFewStations { .. }
            | RotorError::BadStations { .. }
            | RotorError::BadPolar { .. } => ErrorCategory::Config,
            RotorError::NoConvergence { .. } => ErrorCategory::Algorithm,
        }
    }
}

/// Return `value` when it is finite and strictly positive, else
/// [`RotorError::NonPositive`].
pub(crate) fn require_positive(name: &'static str, value: f64) -> Result<f64, RotorError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(RotorError::NonPositive { name, value })
    }
}

/// Return `value` when it is finite (any sign / zero allowed), else
/// [`RotorError::NotFinite`].
pub(crate) fn require_finite(name: &'static str, value: f64) -> Result<f64, RotorError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(RotorError::NotFinite { name, value })
    }
}
