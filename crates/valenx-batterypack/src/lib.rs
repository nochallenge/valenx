//! # valenx-batterypack — series / parallel battery-pack sizing
//!
//! A small, dependency-light calculator that turns one cell plus a
//! series / parallel (S×P) topology into the nominal electrical figures
//! of the assembled pack: voltage, capacity, energy, cell count, and
//! the C-rate ↔ current relations. It is the same software *category*
//! as the pack-sizing calculators bundled with battery configurators
//! and EV / pack-design spreadsheets — the closed-form arithmetic, made
//! explicit and validated.
//!
//! ## What
//!
//! - [`Cell`] — a single cell: nominal voltage `V_cell` (volts) and
//!   capacity `Q_cell` (ampere-hours), with a validated constructor and
//!   the per-cell energy `Wh = V·Ah`.
//! - [`PackConfig`] — `series` cells in series, `parallel` strings in
//!   parallel, exposing [`pack_voltage_v`](PackConfig::pack_voltage_v),
//!   [`pack_capacity_ah`](PackConfig::pack_capacity_ah),
//!   [`pack_energy_wh`](PackConfig::pack_energy_wh) and
//!   [`total_cells`](PackConfig::total_cells).
//! - [`current_from_c_rate`] / [`c_rate_from_current`] /
//!   [`runtime_hours_at_c_rate`] — the constant-current C-rate
//!   relations, in module [`rate`].
//! - [`BatteryPackError`] — the validated-constructor error type, in
//!   module [`error`], with stable
//!   [`code`](BatteryPackError::code) /
//!   [`category`](BatteryPackError::category) accessors.
//!
//! ## Model
//!
//! Write `n = series`, `m = parallel`, with a cell of voltage `V_cell`
//! and capacity `Q_cell`. The textbook closed-form relations are:
//!
//! - Series adds voltage at constant capacity: `V_pack = n · V_cell`.
//! - Parallel adds capacity at constant voltage: `Ah_pack = m · Q_cell`.
//! - Energy is the product: `Wh_pack = V_pack · Ah_pack
//!   = n · m · V_cell · Q_cell`.
//! - C-rate to current against a capacity `Q`: `I = C · Q` (amperes);
//!   the inverse is `C = I / Q`, and the nominal run-time is
//!   `t = 1 / C` hours.
//! - Total cell count: `n · m`.
//!
//! Worked example — a `13S4P` pack of 3.7 V / 3.0 Ah cells:
//! `V_pack = 13 · 3.7 = 48.1` V, `Ah_pack = 4 · 3.0 = 12.0` Ah,
//! `Wh_pack = 48.1 · 12.0 = 577.2` Wh, across `13 · 4 = 52` cells; a
//! `2C` discharge of the pack draws `2 · 12.0 = 24.0` A for `0.5` h.
//!
//! ## Honest scope
//!
//! **Research / educational grade.** This crate implements the standard
//! textbook closed-form / numerical models only. It is NOT a clinical,
//! medical, or production battery-engineering / safety tool, and must
//! not be used to design, qualify, or certify a real battery pack.
//!
//! Every quantity here is *nominal* and *idealised*. The model treats
//! all cells as identical and lossless and assumes perfect series /
//! parallel addition. A physical pack instead has cell-to-cell
//! imbalance, internal resistance and voltage sag under load,
//! temperature dependence, ageing / capacity fade, the Peukert
//! reduction in usable capacity at high C-rates, and a protection /
//! balancing circuit — none of which this crate models. Treat its
//! output as a first-order sizing estimate, never as a safety or
//! certification figure.
//!
//! ## Errors
//!
//! Every fallible function returns
//! [`Result<_, BatteryPackError>`](error::BatteryPackError). Inputs are
//! validated up front: voltages and capacities must be finite and
//! strictly positive; series and parallel counts must be at least one.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod cell;
pub mod error;
pub mod pack;
pub mod rate;

pub use cell::Cell;
pub use error::{BatteryPackError, ErrorCategory};
pub use pack::PackConfig;
pub use rate::{c_rate_from_current, current_from_c_rate, runtime_hours_at_c_rate};

#[cfg(test)]
mod tests {
    use super::*;

    /// End-to-end ground-truth check of the full 13S4P worked example
    /// quoted in the crate docs, tying every module together.
    #[test]
    fn worked_example_13s4p_end_to_end() {
        let cell = Cell::new(3.7, 3.0).unwrap();
        let pack = PackConfig::new(cell, 13, 4).unwrap();

        // Topology figures.
        assert!((pack.pack_voltage_v() - 48.1).abs() < 1e-9);
        assert!((pack.pack_capacity_ah() - 12.0).abs() < 1e-9);
        assert!((pack.pack_energy_wh() - 577.2).abs() < 1e-7);
        assert_eq!(pack.total_cells(), 52);

        // A 2C discharge of the pack: 24 A for half an hour.
        let i = current_from_c_rate(2.0, pack.pack_capacity_ah()).unwrap();
        assert!((i - 24.0).abs() < 1e-9);
        assert!((runtime_hours_at_c_rate(2.0).unwrap() - 0.5).abs() < 1e-12);

        // Charge / discharge over an hour at that current recovers the
        // pack capacity: I * t = 24 A * 0.5 h = 12 Ah.
        let delivered_ah = i * runtime_hours_at_c_rate(2.0).unwrap();
        assert!((delivered_ah - pack.pack_capacity_ah()).abs() < 1e-9);
    }

    /// A `1S1P` pack of a single cell must reproduce that cell exactly:
    /// the degenerate base case of the topology.
    #[test]
    fn single_cell_pack_matches_cell() {
        let cell = Cell::new(3.7, 3.0).unwrap();
        let pack = PackConfig::new(cell, 1, 1).unwrap();
        assert!((pack.pack_voltage_v() - cell.nominal_voltage_v).abs() < 1e-12);
        assert!((pack.pack_capacity_ah() - cell.capacity_ah).abs() < 1e-12);
        assert!((pack.pack_energy_wh() - cell.energy_wh()).abs() < 1e-12);
        assert_eq!(pack.total_cells(), 1);
    }

    /// Energy scales linearly and independently in each dimension:
    /// doubling `series` doubles energy; doubling `parallel` doubles
    /// energy; doing both quadruples it.
    #[test]
    fn energy_scales_with_each_dimension() {
        let cell = Cell::new(3.7, 3.0).unwrap();
        let base = PackConfig::new(cell, 2, 3).unwrap().pack_energy_wh();

        let double_s = PackConfig::new(cell, 4, 3).unwrap().pack_energy_wh();
        assert!((double_s - 2.0 * base).abs() < 1e-9);

        let double_p = PackConfig::new(cell, 2, 6).unwrap().pack_energy_wh();
        assert!((double_p - 2.0 * base).abs() < 1e-9);

        let double_both = PackConfig::new(cell, 4, 6).unwrap().pack_energy_wh();
        assert!((double_both - 4.0 * base).abs() < 1e-9);
    }
}
