//! The single-cell model: a nominal voltage and a charge capacity.
//!
//! A [`Cell`] is the atom every pack is built from. It is described by
//! two textbook quantities:
//!
//! - **Nominal voltage** `V_cell` (volts) — the chemistry's average
//!   working voltage (e.g. ~3.7 V for a Li-ion cell, ~1.2 V for NiMH).
//! - **Capacity** `Q_cell` (ampere-hours) — the charge the cell can
//!   deliver, the integral of current over a full discharge.
//!
//! Their product is the cell's energy in watt-hours, `Wh = V·Ah`.

use serde::{Deserialize, Serialize};

use crate::error::{require_positive, BatteryPackError};

/// A single electrochemical cell, the building block of a pack.
///
/// Construct one with [`Cell::new`], which validates that both the
/// voltage and the capacity are finite and strictly positive. The
/// fields are public for read access but are guaranteed valid by every
/// constructor path.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Cell {
    /// Nominal cell voltage `V_cell`, in volts.
    pub nominal_voltage_v: f64,
    /// Cell charge capacity `Q_cell`, in ampere-hours (Ah).
    pub capacity_ah: f64,
}

impl Cell {
    /// Build a [`Cell`] from a nominal voltage (volts) and a capacity
    /// (ampere-hours).
    ///
    /// # Errors
    ///
    /// Returns [`BatteryPackError::BadParameter`] if either argument is
    /// non-finite or not strictly positive.
    pub fn new(nominal_voltage_v: f64, capacity_ah: f64) -> Result<Self, BatteryPackError> {
        let nominal_voltage_v = require_positive("nominal_voltage_v", nominal_voltage_v)?;
        let capacity_ah = require_positive("capacity_ah", capacity_ah)?;
        Ok(Self {
            nominal_voltage_v,
            capacity_ah,
        })
    }

    /// A representative single Li-ion cell: 3.7 V nominal, 3.0 Ah.
    ///
    /// Handy as a default / example; the exact numbers are illustrative,
    /// not a spec for any particular part.
    pub fn li_ion_18650() -> Self {
        // Constructed from known-valid constants; the unwrap cannot fail.
        Self::new(3.7, 3.0).expect("3.7 V / 3.0 Ah is a valid cell")
    }

    /// Energy stored by this single cell, in watt-hours: `Wh = V·Ah`.
    pub fn energy_wh(&self) -> f64 {
        self.nominal_voltage_v * self.capacity_ah
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_accepts_valid_cell() {
        let c = Cell::new(3.7, 3.0).unwrap();
        assert!((c.nominal_voltage_v - 3.7).abs() < 1e-12);
        assert!((c.capacity_ah - 3.0).abs() < 1e-12);
    }

    #[test]
    fn new_rejects_non_positive_voltage() {
        assert!(matches!(
            Cell::new(0.0, 3.0),
            Err(BatteryPackError::BadParameter {
                name: "nominal_voltage_v",
                ..
            })
        ));
        assert!(matches!(
            Cell::new(-3.7, 3.0),
            Err(BatteryPackError::BadParameter {
                name: "nominal_voltage_v",
                ..
            })
        ));
    }

    #[test]
    fn new_rejects_non_positive_capacity() {
        assert!(matches!(
            Cell::new(3.7, 0.0),
            Err(BatteryPackError::BadParameter {
                name: "capacity_ah",
                ..
            })
        ));
    }

    #[test]
    fn new_rejects_non_finite() {
        assert!(Cell::new(f64::NAN, 3.0).is_err());
        assert!(Cell::new(3.7, f64::INFINITY).is_err());
    }

    #[test]
    fn energy_is_voltage_times_capacity() {
        // 3.7 V * 3.0 Ah = 11.1 Wh.
        let c = Cell::new(3.7, 3.0).unwrap();
        assert!((c.energy_wh() - 11.1).abs() < 1e-9);
    }

    #[test]
    fn reference_cell_is_valid() {
        let c = Cell::li_ion_18650();
        assert!((c.energy_wh() - 11.1).abs() < 1e-9);
    }
}
