//! Series / parallel (S×P) pack topology and the derived pack
//! quantities.
//!
//! A pack is `series` identical cells wired in series to form one
//! *string*, with `parallel` such strings wired in parallel. Writing
//! `n = series` and `m = parallel`, with a cell of voltage `V_cell` and
//! capacity `Q_cell`, the textbook closed-form relations are:
//!
//! - Pack voltage: `V_pack = n · V_cell` (series cells add voltage;
//!   capacity is unchanged along a string).
//! - Pack capacity: `Ah_pack = m · Q_cell` (parallel strings add
//!   capacity; voltage is unchanged across the parallel set).
//! - Pack energy: `Wh_pack = V_pack · Ah_pack = n · m · V_cell · Q_cell`.
//! - Total cells: `n · m`.
//!
//! These are *nominal* quantities. A physical pack also has
//! cell-to-cell imbalance, internal resistance, voltage sag under load,
//! ageing, and a protection / balancing circuit — none of which this
//! model captures (see the crate-level honest-scope note).

use serde::{Deserialize, Serialize};

use crate::cell::Cell;
use crate::error::{require_at_least_one, BatteryPackError};

/// A series / parallel battery pack built from one repeated [`Cell`].
///
/// Construct one with [`PackConfig::new`], which validates that both
/// multiplicities are at least one. The naming convention `nSmP` (e.g.
/// `13S4P`) reads as `series` cells in series, `parallel` strings in
/// parallel.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PackConfig {
    /// The repeated cell.
    pub cell: Cell,
    /// Number of cells in series per string, `n`. Sets the voltage.
    pub series: u32,
    /// Number of parallel strings, `m`. Sets the capacity.
    pub parallel: u32,
}

impl PackConfig {
    /// Build a pack from a [`Cell`], a series count `series` (`n`) and a
    /// parallel count `parallel` (`m`).
    ///
    /// # Errors
    ///
    /// Returns [`BatteryPackError::BadCount`] if either `series` or
    /// `parallel` is zero — a pack must have at least one cell in each
    /// dimension.
    pub fn new(cell: Cell, series: u32, parallel: u32) -> Result<Self, BatteryPackError> {
        let series = require_at_least_one("series", series)?;
        let parallel = require_at_least_one("parallel", parallel)?;
        Ok(Self {
            cell,
            series,
            parallel,
        })
    }

    /// Pack nominal voltage `V_pack = n · V_cell`, in volts.
    ///
    /// Series cells add their voltages; the string capacity equals a
    /// single cell's capacity.
    pub fn pack_voltage_v(&self) -> f64 {
        f64::from(self.series) * self.cell.nominal_voltage_v
    }

    /// Pack capacity `Ah_pack = m · Q_cell`, in ampere-hours.
    ///
    /// Parallel strings add their capacities; the pack voltage equals a
    /// single string's voltage.
    pub fn pack_capacity_ah(&self) -> f64 {
        f64::from(self.parallel) * self.cell.capacity_ah
    }

    /// Pack energy `Wh_pack = V_pack · Ah_pack`, in watt-hours.
    ///
    /// Equivalently `n · m · V_cell · Q_cell`, i.e. the single-cell
    /// energy scaled by the total cell count.
    pub fn pack_energy_wh(&self) -> f64 {
        self.pack_voltage_v() * self.pack_capacity_ah()
    }

    /// Total number of cells in the pack, `n · m`.
    ///
    /// Returned as [`u64`] so large packs cannot overflow `u32` when
    /// `series` and `parallel` are both near the `u32` maximum.
    pub fn total_cells(&self) -> u64 {
        u64::from(self.series) * u64::from(self.parallel)
    }

    /// Total energy summed cell-by-cell, `n · m · Wh_cell`.
    ///
    /// Mathematically identical to [`pack_energy_wh`](Self::pack_energy_wh)
    /// (both equal `n · m · V_cell · Q_cell`); provided as an
    /// independent cross-check of the topology arithmetic.
    pub fn total_cell_energy_wh(&self) -> f64 {
        self.total_cells() as f64 * self.cell.energy_wh()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A `13S4P` pack of 3.7 V / 3.0 Ah cells — the canonical worked
    /// example used across these tests.
    fn pack_13s4p() -> PackConfig {
        let cell = Cell::new(3.7, 3.0).unwrap();
        PackConfig::new(cell, 13, 4).unwrap()
    }

    #[test]
    fn rejects_zero_series_or_parallel() {
        let cell = Cell::new(3.7, 3.0).unwrap();
        assert!(matches!(
            PackConfig::new(cell, 0, 4),
            Err(BatteryPackError::BadCount { name: "series", .. })
        ));
        assert!(matches!(
            PackConfig::new(cell, 13, 0),
            Err(BatteryPackError::BadCount {
                name: "parallel",
                ..
            })
        ));
    }

    #[test]
    fn series_sets_voltage_capacity_unchanged() {
        // 13S1P: V = 13 * 3.7 = 48.1 V; Ah = 1 * 3.0 = 3.0 Ah.
        let cell = Cell::new(3.7, 3.0).unwrap();
        let pack = PackConfig::new(cell, 13, 1).unwrap();
        assert!((pack.pack_voltage_v() - 48.1).abs() < 1e-9);
        assert!((pack.pack_capacity_ah() - 3.0).abs() < 1e-12);
    }

    #[test]
    fn parallel_sets_capacity_voltage_unchanged() {
        // 1S4P: V = 1 * 3.7 = 3.7 V; Ah = 4 * 3.0 = 12.0 Ah.
        let cell = Cell::new(3.7, 3.0).unwrap();
        let pack = PackConfig::new(cell, 1, 4).unwrap();
        assert!((pack.pack_voltage_v() - 3.7).abs() < 1e-12);
        assert!((pack.pack_capacity_ah() - 12.0).abs() < 1e-9);
    }

    #[test]
    fn sxp_voltage_and_capacity() {
        // 13S4P: V = 13 * 3.7 = 48.1 V; Ah = 4 * 3.0 = 12.0 Ah.
        let pack = pack_13s4p();
        assert!((pack.pack_voltage_v() - 48.1).abs() < 1e-9);
        assert!((pack.pack_capacity_ah() - 12.0).abs() < 1e-9);
    }

    #[test]
    fn energy_equals_voltage_times_capacity() {
        // 48.1 V * 12.0 Ah = 577.2 Wh.
        let pack = pack_13s4p();
        assert!((pack.pack_energy_wh() - 577.2).abs() < 1e-7);
    }

    #[test]
    fn total_cells_is_series_times_parallel() {
        let pack = pack_13s4p();
        assert_eq!(pack.total_cells(), 52);
    }

    #[test]
    fn pack_energy_equals_per_cell_energy_sum() {
        // Two independent routes to the same number must agree:
        //   V_pack * Ah_pack   ==   (n*m) * (V_cell * Ah_cell)
        let pack = pack_13s4p();
        assert!((pack.pack_energy_wh() - pack.total_cell_energy_wh()).abs() < 1e-7);
    }

    #[test]
    fn total_cells_does_not_overflow_u32() {
        // 100_000 * 100_000 = 1e10, which overflows u32 (~4.29e9) but
        // fits comfortably in the u64 return type.
        let cell = Cell::new(3.7, 3.0).unwrap();
        let pack = PackConfig::new(cell, 100_000, 100_000).unwrap();
        assert_eq!(pack.total_cells(), 10_000_000_000);
    }
}
