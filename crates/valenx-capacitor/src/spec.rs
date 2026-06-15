//! Serialisable parallel-plate descriptor.
//!
//! ## Model
//!
//! [`ParallelPlate`] bundles the three geometric / material inputs of a
//! parallel-plate capacitor (relative permittivity, plate area, gap) into
//! one validated, [`serde`]-serialisable record so a design can round-trip
//! to RON / JSON and re-derive its capacitance, energy and charge on
//! demand. It is a thin convenience wrapper over the free functions in
//! [`crate::parallel_plate`].
//!
//! ## Honest scope
//!
//! The same idealisations as [`crate::parallel_plate`] apply: no fringing,
//! no dielectric loss, no breakdown limit.

use serde::{Deserialize, Serialize};

use crate::error::Result;
use crate::parallel_plate;

/// A validated parallel-plate capacitor geometry, in SI units.
///
/// Construct one with [`ParallelPlate::new`], which validates the inputs;
/// the fields are then guaranteed to be in their physical domain
/// (`eps_r >= 1`, `area_m2 > 0`, `gap_m > 0`).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ParallelPlate {
    /// Relative permittivity of the dielectric (dimensionless, `>= 1`).
    pub eps_r: f64,
    /// Overlapping plate area `A`, in square metres (`> 0`).
    pub area_m2: f64,
    /// Plate separation `d`, in metres (`> 0`).
    pub gap_m: f64,
}

impl ParallelPlate {
    /// Build a validated descriptor.
    ///
    /// # Errors
    ///
    /// Returns [`crate::CapacitorError::InvalidParameter`] if `eps_r < 1`
    /// or if `area_m2` / `gap_m` is not strictly positive and finite. The
    /// check is performed by re-using [`parallel_plate::capacitance`], so
    /// a successfully constructed value always yields a finite
    /// capacitance.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_capacitor::spec::ParallelPlate;
    ///
    /// let p = ParallelPlate::new(2.2, 5.0e-4, 25.0e-6).unwrap();
    /// assert!(p.capacitance().unwrap() > 0.0);
    /// assert!(ParallelPlate::new(1.0, -1.0, 1.0e-3).is_err());
    /// ```
    pub fn new(eps_r: f64, area_m2: f64, gap_m: f64) -> Result<Self> {
        // Validate by attempting the derived computation; discard the
        // result and keep the inputs only once they are known-good.
        parallel_plate::capacitance(eps_r, area_m2, gap_m)?;
        Ok(Self {
            eps_r,
            area_m2,
            gap_m,
        })
    }

    /// Capacitance of this geometry, in farads.
    ///
    /// Delegates to [`parallel_plate::capacitance`].
    ///
    /// # Errors
    ///
    /// Propagates any validation error, though a value built through
    /// [`ParallelPlate::new`] never fails here.
    pub fn capacitance(&self) -> Result<f64> {
        parallel_plate::capacitance(self.eps_r, self.area_m2, self.gap_m)
    }

    /// Energy stored at terminal voltage `voltage_v`, in joules.
    ///
    /// Computes the capacitance, then `E = 1/2 C V^2` via
    /// [`parallel_plate::stored_energy`].
    ///
    /// # Errors
    ///
    /// Returns [`crate::CapacitorError::InvalidParameter`] if `voltage_v`
    /// is non-finite.
    pub fn stored_energy(&self, voltage_v: f64) -> Result<f64> {
        let c = self.capacitance()?;
        parallel_plate::stored_energy(c, voltage_v)
    }

    /// Charge stored at terminal voltage `voltage_v`, in coulombs.
    ///
    /// Computes the capacitance, then `Q = C V` via
    /// [`parallel_plate::charge`].
    ///
    /// # Errors
    ///
    /// Returns [`crate::CapacitorError::InvalidParameter`] if `voltage_v`
    /// is non-finite.
    pub fn charge(&self, voltage_v: f64) -> Result<f64> {
        let c = self.capacitance()?;
        parallel_plate::charge(c, voltage_v)
    }
}
