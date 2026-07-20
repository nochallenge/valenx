//! The Willis equation and the three named single-input speed ratios.
//!
//! ## Willis relation
//!
//! For a simple sun-planet-ring set the relative angular velocities,
//! measured in the rotating frame of the carrier, obey the ordinary
//! (carrier-fixed) gear ratio `e = omega_ring / omega_sun = -S/R`:
//!
//! ```text
//! (omega_ring - omega_carrier) / (omega_sun - omega_carrier) = e
//! ```
//!
//! Clearing the denominator gives the linear constraint this module
//! solves:
//!
//! ```text
//! omega_ring - e*omega_sun + (e - 1)*omega_carrier = 0
//! ```
//!
//! Given any two of `omega_sun`, `omega_ring`, `omega_carrier`, the
//! third follows directly. All angular velocities here share a single
//! arbitrary unit (rad/s, rev/min, ...); the relations are homogeneous,
//! so the output carries whatever unit the inputs were given in.
//!
//! ## Named single-input cases (one member grounded)
//!
//! With one central member held stationary the set becomes a fixed
//! reduction. The closed forms below are the textbook results and are
//! re-derived in the unit tests against the general solver:
//!
//! ```text
//! ring fixed:    omega_carrier / omega_sun    =  S / (S + R)
//! sun  fixed:    omega_carrier / omega_ring   =  R / (S + R)
//! carrier fixed: omega_ring    / omega_sun    = -S / R
//! ```

use crate::error::{require_finite, EpicyclicError};
use crate::teeth::TeethSet;

/// Resolved angular velocities of every member of the set, in the same
/// unit the inputs were supplied in.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TrainState {
    /// Sun angular velocity.
    pub omega_sun: f64,
    /// Ring (annulus) angular velocity.
    pub omega_ring: f64,
    /// Carrier (planet-arm) angular velocity.
    pub omega_carrier: f64,
}

impl TrainState {
    /// Absolute spin of a single planet about its own axis, in the same
    /// unit as the member velocities.
    ///
    /// In the carrier frame the planet turns at
    /// `(omega_sun - omega_carrier) * (-S/P)`; adding the carrier's own
    /// rotation back gives the planet's absolute spin.
    ///
    /// # Example
    ///
    /// ```
    /// use valenx_epicyclic::{TeethSet, willis};
    /// let t = TeethSet::new(24, 24, 72).unwrap();
    /// // Ring fixed, sun at 96, carrier comes out at 24.
    /// let s = willis::solve_carrier_ring_fixed(&t, 96.0).unwrap();
    /// let wp = s.planet_spin(&t);
    /// // carrier-frame: (96-24)*(-1) = -72, plus carrier 24 -> -48.
    /// assert!((wp - (-48.0)).abs() < 1e-9);
    /// ```
    pub fn planet_spin(&self, teeth: &TeethSet) -> f64 {
        let rel = (self.omega_sun - self.omega_carrier) * teeth.planet_sun_ratio_carrier_fixed();
        rel + self.omega_carrier
    }
}

/// Solve for the carrier velocity given the sun and ring velocities.
///
/// From `omega_ring - e*omega_sun + (e - 1)*omega_carrier = 0`:
/// `omega_carrier = (e*omega_sun - omega_ring) / (e - 1)`.
///
/// # Errors
///
/// Returns [`EpicyclicError::NonFinite`] if either input is not finite.
/// The coefficient `e - 1` is `-(S + R)/R`, which is strictly negative
/// for a real set, so this solve is never singular.
pub fn solve_carrier(
    teeth: &TeethSet,
    omega_sun: f64,
    omega_ring: f64,
) -> Result<f64, EpicyclicError> {
    let omega_sun = require_finite("omega_sun", omega_sun)?;
    let omega_ring = require_finite("omega_ring", omega_ring)?;
    let e = teeth.basic_train_ratio();
    Ok((e * omega_sun - omega_ring) / (e - 1.0))
}

/// Solve for the sun velocity given the ring and carrier velocities.
///
/// `omega_sun = (omega_ring + (e - 1)*omega_carrier) / e`.
///
/// # Errors
///
/// Returns [`EpicyclicError::NonFinite`] if either input is not finite,
/// or [`EpicyclicError::Singular`] if `e == 0` (which cannot occur for a
/// validated [`TeethSet`], whose `e = -S/R` is always non-zero, but is
/// guarded for completeness).
pub fn solve_sun(
    teeth: &TeethSet,
    omega_ring: f64,
    omega_carrier: f64,
) -> Result<f64, EpicyclicError> {
    let omega_ring = require_finite("omega_ring", omega_ring)?;
    let omega_carrier = require_finite("omega_carrier", omega_carrier)?;
    let e = teeth.basic_train_ratio();
    if e == 0.0 {
        return Err(EpicyclicError::Singular {
            reason: "basic train ratio e = 0; sun velocity is unconstrained",
        });
    }
    Ok((omega_ring + (e - 1.0) * omega_carrier) / e)
}

/// Solve for the ring velocity given the sun and carrier velocities.
///
/// `omega_ring = e*omega_sun - (e - 1)*omega_carrier`.
///
/// # Errors
///
/// Returns [`EpicyclicError::NonFinite`] if either input is not finite.
pub fn solve_ring(
    teeth: &TeethSet,
    omega_sun: f64,
    omega_carrier: f64,
) -> Result<f64, EpicyclicError> {
    let omega_sun = require_finite("omega_sun", omega_sun)?;
    let omega_carrier = require_finite("omega_carrier", omega_carrier)?;
    let e = teeth.basic_train_ratio();
    Ok(e * omega_sun - (e - 1.0) * omega_carrier)
}

/// Ring fixed, sun input: solve the carrier (output) velocity.
///
/// Closed form `omega_carrier = omega_sun * S / (S + R)`, the classic
/// planetary reduction with the annulus grounded.
///
/// # Errors
///
/// Returns [`EpicyclicError::NonFinite`] if `omega_sun` is not finite.
///
/// # Example
///
/// ```
/// use valenx_epicyclic::{TeethSet, willis};
/// // S = 24, R = 72: reduction carrier/sun = 24/96 = 1/4.
/// let t = TeethSet::new(24, 24, 72).unwrap();
/// let s = willis::solve_carrier_ring_fixed(&t, 100.0).unwrap();
/// assert!((s.omega_carrier - 25.0).abs() < 1e-9);
/// assert!((s.omega_ring).abs() < 1e-12);
/// ```
pub fn solve_carrier_ring_fixed(
    teeth: &TeethSet,
    omega_sun: f64,
) -> Result<TrainState, EpicyclicError> {
    let omega_sun = require_finite("omega_sun", omega_sun)?;
    let omega_carrier = solve_carrier(teeth, omega_sun, 0.0)?;
    Ok(TrainState {
        omega_sun,
        omega_ring: 0.0,
        omega_carrier,
    })
}

/// Sun fixed, ring input: solve the carrier (output) velocity.
///
/// Closed form `omega_carrier = omega_ring * R / (S + R)`.
///
/// # Errors
///
/// Returns [`EpicyclicError::NonFinite`] if `omega_ring` is not finite.
///
/// # Example
///
/// ```
/// use valenx_epicyclic::{TeethSet, willis};
/// // S = 24, R = 72: carrier/ring = 72/96 = 3/4.
/// let t = TeethSet::new(24, 24, 72).unwrap();
/// let s = willis::solve_carrier_sun_fixed(&t, 100.0).unwrap();
/// assert!((s.omega_carrier - 75.0).abs() < 1e-9);
/// assert!((s.omega_sun).abs() < 1e-12);
/// ```
pub fn solve_carrier_sun_fixed(
    teeth: &TeethSet,
    omega_ring: f64,
) -> Result<TrainState, EpicyclicError> {
    let omega_ring = require_finite("omega_ring", omega_ring)?;
    let omega_carrier = solve_carrier(teeth, 0.0, omega_ring)?;
    Ok(TrainState {
        omega_sun: 0.0,
        omega_ring,
        omega_carrier,
    })
}

/// Carrier fixed, sun input: solve the ring (output) velocity.
///
/// Closed form `omega_ring = -omega_sun * S / R`, the basic train ratio
/// applied directly. The ring turns opposite to the sun.
///
/// # Errors
///
/// Returns [`EpicyclicError::NonFinite`] if `omega_sun` is not finite.
///
/// # Example
///
/// ```
/// use valenx_epicyclic::{TeethSet, willis};
/// // S = 24, R = 72: ring/sun = -24/72 = -1/3.
/// let t = TeethSet::new(24, 24, 72).unwrap();
/// let s = willis::solve_ring_carrier_fixed(&t, 90.0).unwrap();
/// assert!((s.omega_ring - (-30.0)).abs() < 1e-9);
/// assert!((s.omega_carrier).abs() < 1e-12);
/// ```
pub fn solve_ring_carrier_fixed(
    teeth: &TeethSet,
    omega_sun: f64,
) -> Result<TrainState, EpicyclicError> {
    let omega_sun = require_finite("omega_sun", omega_sun)?;
    let omega_ring = solve_ring(teeth, omega_sun, 0.0)?;
    Ok(TrainState {
        omega_sun,
        omega_ring,
        omega_carrier: 0.0,
    })
}
