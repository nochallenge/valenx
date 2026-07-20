//! The tabular (superposition) method for resolving every member speed.
//!
//! ## The two-row recipe
//!
//! The classic textbook bookkeeping (Shigley, Norton) decomposes any
//! planetary motion into two rigid-body steps that are then summed:
//!
//! ```text
//! Row 1  lock the whole train, rotate everything by  x:
//!        sun = x,   ring = x,        carrier = x,   planet = x
//!
//! Row 2  hold the carrier, rotate the sun by  y:
//!        sun = y,   ring = y*e,      carrier = 0,   planet = y*(-S/P)
//!
//! Sum    sun = x + y
//!        ring = x + y*e
//!        carrier = x
//!        planet = x + y*(-S/P)
//! ```
//!
//! where `e = -S/R` is the carrier-fixed basic train ratio. The carrier
//! speed is simply `x`, which is what makes the method so direct: pick
//! the two known speeds, solve the resulting pair of linear equations
//! for `x` and `y`, then read every member off the sum row.
//!
//! This module exposes the raw [`TableRows`] for inspection and a
//! [`resolve`] entry point that takes the carrier speed plus one driving
//! member speed (sun or ring) — the dominant practical case — and
//! returns the full resolved state. The results are identical to the
//! [`crate::willis`] closed forms; the tests assert that equivalence.

use crate::error::{require_finite, EpicyclicError};
use crate::teeth::TeethSet;
use crate::willis::TrainState;

/// The two superposition rows plus their column sum, all in the input
/// speed unit. Each field is the contribution of that member.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TableRows {
    /// `x`: the whole-train rotation applied in row 1.
    pub x: f64,
    /// `y`: the carrier-fixed sun rotation applied in row 2.
    pub y: f64,
    /// Resolved member speeds (the summed bottom row).
    pub state: TrainState,
    /// Resolved absolute planet spin (the summed bottom row, planet
    /// column).
    pub omega_planet: f64,
}

/// Which member, alongside the carrier, supplies the second known speed
/// for [`resolve`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Driver {
    /// The sun speed is the second known.
    Sun,
    /// The ring speed is the second known.
    Ring,
}

/// Resolve the whole train by the tabular method from the carrier speed
/// and one driving member speed.
///
/// With `carrier = x` known directly, row 2's amount `y` is recovered
/// from the chosen driver:
///
/// ```text
/// driver = Sun:   sun  = x + y       -> y = sun  - x
/// driver = Ring:  ring = x + y*e     -> y = (ring - x) / e
/// ```
///
/// and then every member (including the planet) is read from the sum
/// row.
///
/// # Errors
///
/// Returns [`EpicyclicError::NonFinite`] if either speed is not finite.
/// When the driver is [`Driver::Ring`] the divisor is `e = -S/R`, which
/// is non-zero for a validated [`TeethSet`]; the case is nonetheless
/// guarded and reported as [`EpicyclicError::Singular`].
///
/// # Example
///
/// ```
/// use valenx_epicyclic::{tabular, tabular::Driver, TeethSet};
/// // S=24, R=72. Ring fixed (carrier-unknown style restated): drive the
/// // sun at 96 while the carrier turns at 24 -> ring must be 0.
/// let t = TeethSet::new(24, 24, 72).unwrap();
/// let rows = tabular::resolve(&t, 24.0, Driver::Sun, 96.0).unwrap();
/// assert!((rows.state.omega_ring).abs() < 1e-9);
/// assert!((rows.x - 24.0).abs() < 1e-12);
/// assert!((rows.y - 72.0).abs() < 1e-12);
/// ```
pub fn resolve(
    teeth: &TeethSet,
    omega_carrier: f64,
    driver: Driver,
    omega_driver: f64,
) -> Result<TableRows, EpicyclicError> {
    let x = require_finite("omega_carrier", omega_carrier)?;
    let driver_val = require_finite("omega_driver", omega_driver)?;
    let e = teeth.basic_train_ratio();

    let y = match driver {
        Driver::Sun => driver_val - x,
        Driver::Ring => {
            if e == 0.0 {
                return Err(EpicyclicError::Singular {
                    reason: "basic train ratio e = 0; ring cannot pin the table",
                });
            }
            (driver_val - x) / e
        }
    };

    let omega_sun = x + y;
    let omega_ring = x + y * e;
    let omega_planet = x + y * teeth.planet_sun_ratio_carrier_fixed();

    Ok(TableRows {
        x,
        y,
        state: TrainState {
            omega_sun,
            omega_ring,
            omega_carrier: x,
        },
        omega_planet,
    })
}
