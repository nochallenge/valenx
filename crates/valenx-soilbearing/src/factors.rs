//! Bearing-capacity factors `Nc`, `Nq`, and `Ngamma`.
//!
//! The three dimensionless factors that weight the cohesion,
//! surcharge, and self-weight terms of the Terzaghi bearing-capacity
//! equation are each a closed-form function of the soil drained
//! friction angle `phi` alone.
//!
//! # Formulae
//!
//! With `phi` in radians:
//!
//! `Nq = exp(pi * tan(phi)) * tan^2(pi/4 + phi/2)`
//!
//! `Nc = (Nq - 1) * cot(phi)`  (the Prandtl/Reissner form)
//!
//! `Ngamma = 2 * (Nq + 1) * tan(phi)`  (the Vesic 1973 form)
//!
//! These are evaluated by [`BearingFactors::from_friction_angle`].
//!
//! # Limiting case `phi = 0` (undrained / cohesion-only)
//!
//! As `phi -> 0` the expressions above are taken in the limit:
//!
//! `Nq = 1`, `Ngamma = 0`, and `Nc = pi + 2 = 5.14159...`
//!
//! (the classical Prandtl value). The implementation special-cases a
//! tiny neighbourhood of `phi = 0` to avoid the `0 * cot(0)`
//! indeterminate form, returning these exact limits.

use serde::{Deserialize, Serialize};

use crate::soil::SoilProperties;

/// Below this friction angle (radians) the factors are evaluated at the
/// `phi -> 0` limit rather than through the general formulae, to dodge
/// the `cot(phi)` singularity in `Nc`. `1e-9 rad` is far smaller than
/// any physically meaningful friction angle.
const PHI_ZERO_TOL: f64 = 1e-9;

/// The classical Prandtl cohesion factor at `phi = 0`, `Nc = pi + 2`.
const NC_PHI_ZERO: f64 = std::f64::consts::PI + 2.0;

/// The three dimensionless Terzaghi bearing-capacity factors.
///
/// All three are pure functions of the soil friction angle; see the
/// [module documentation](self) for the closed forms.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct BearingFactors {
    /// Cohesion factor `Nc`, weighting the `c * Nc` term.
    pub nc: f64,
    /// Surcharge factor `Nq`, weighting the `q * Nq` term.
    pub nq: f64,
    /// Self-weight factor `Ngamma`, weighting the `0.5 * gamma * B * Ngamma` term.
    pub ngamma: f64,
}

impl BearingFactors {
    /// Evaluate `Nc`, `Nq`, and `Ngamma` from the soil friction angle.
    ///
    /// The friction angle is read from `soil` (already validated to lie
    /// in `[0, 90)` degrees by [`SoilProperties::new`]), so this is
    /// infallible.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_soilbearing::{BearingFactors, SoilProperties};
    ///
    /// // Purely cohesive soil (phi = 0): Nq = 1, Ngamma = 0, Nc = pi + 2.
    /// let clay = SoilProperties::new(0.0, 25.0, 18.0).unwrap();
    /// let f = BearingFactors::from_friction_angle(&clay);
    /// assert!((f.nq - 1.0).abs() < 1e-12);
    /// assert!(f.ngamma.abs() < 1e-12);
    /// assert!((f.nc - (std::f64::consts::PI + 2.0)).abs() < 1e-12);
    /// ```
    pub fn from_friction_angle(soil: &SoilProperties) -> Self {
        let phi = soil.friction_angle_rad();

        if phi <= PHI_ZERO_TOL {
            return BearingFactors {
                nc: NC_PHI_ZERO,
                nq: 1.0,
                ngamma: 0.0,
            };
        }

        let tan_phi = phi.tan();
        // Nq = e^(pi tan phi) * tan^2(45 deg + phi/2).
        let tan_half = (std::f64::consts::FRAC_PI_4 + phi / 2.0).tan();
        let nq = (std::f64::consts::PI * tan_phi).exp() * tan_half * tan_half;
        // Nc = (Nq - 1) cot(phi) = (Nq - 1) / tan(phi).
        let nc = (nq - 1.0) / tan_phi;
        // Ngamma = 2 (Nq + 1) tan(phi)  (Vesic 1973).
        let ngamma = 2.0 * (nq + 1.0) * tan_phi;

        BearingFactors { nc, nq, ngamma }
    }
}
