//! # valenx-pneumatics — compressed-air sizing helpers
//!
//! Closed-form textbook models for the four everyday pneumatics
//! calculations: cylinder thrust, air consumption, the compression ratio,
//! and whether discharge flow is choked.
//!
//! ## What
//!
//! - **Cylinder force** ([`cylinder`]) — the theoretical thrust
//!   `F = p_gauge * A` of a [`Cylinder`], with single- and double-acting
//!   geometry (the rod steals area on the retract [`Stroke`]).
//! - **Air consumption** ([`consumption`]) — the free-air volume a
//!   reciprocating cylinder draws, swept volume `A * L` per stroke,
//!   normalised to atmosphere by the compression ratio and scaled by the
//!   cycle count or rate.
//! - **Compression ratio** ([`compression`]) — absolute-vs-gauge pressure
//!   bookkeeping and the ratio `r = p_abs / p_atm` that links compressed
//!   and free-air volumes.
//! - **Choked flow** ([`flow`]) — the isentropic critical pressure ratio
//!   `(2/(k+1))^(k/(k-1))` (about `0.528` for air) and an [`is_choked`]
//!   predicate over absolute upstream / downstream pressures.
//!
//! ## Model
//!
//! Everything is the standard ideal-gas / isentropic textbook treatment in
//! SI units (metres, pascals, square metres, cubic metres, newtons):
//!
//! ```text
//! F = p_gauge * A                          cylinder thrust
//! A = pi/4 * d_bore^2   (extend)           full-bore area
//! A = pi/4 * (d_bore^2 - d_rod^2) (retract) annular area
//! r = p_abs / p_atm = 1 + p_gauge / p_atm  compression ratio
//! Q = A * L * r * N                        free-air consumption (N cycles)
//! (p_down/p_up)_crit = (2/(k+1))^(k/(k-1)) choked-flow threshold
//! ```
//!
//! Gauge pressure drives cylinder force (the atmosphere balances the
//! opposite face); absolute pressure drives the gas-law conversions and the
//! choked-flow ratio.
//!
//! ## Honest scope
//!
//! Research / educational grade. These are textbook closed-form and
//! numerical models under ideal-gas and isentropic assumptions. They
//! deliberately ignore seal and rod friction, breakaway stiction,
//! exhaust back-pressure, port / hose dead volume, line and fitting
//! losses, valve flow coefficients (`Cv`/`Kv`), leakage, heat transfer,
//! humidity and real-gas compressibility. Cylinder thrust is an ideal
//! upper bound (real delivered force is lower); consumption is the
//! irreducible stroke demand (real demand is higher); the choked-flow
//! routines decide *whether* flow chokes and at *what* ratio, not the
//! actual mass-flow rate. This is **not** a clinical, medical, or
//! production-engineering tool — do not size safety-critical or
//! life-support pneumatics with it.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod compression;
pub mod consumption;
pub mod cylinder;
pub mod error;
pub mod flow;

// --- Convenience re-exports of the most-used items --------------------

pub use error::{PneumaticsError, Result};

pub use cylinder::{Cylinder, Stroke};

pub use compression::{
    absolute_pressure, compression_ratio, compression_ratio_from_absolute, STANDARD_ATMOSPHERE_PA,
};

pub use consumption::{
    free_air_consumption, free_air_flow_demand, free_air_per_stroke, swept_volume,
    swept_volume_per_cycle, Action,
};

pub use flow::{
    critical_pressure_ratio, is_choked, is_choked_air, pressure_ratio, CRITICAL_PRESSURE_RATIO_AIR,
    GAMMA_AIR,
};

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute-difference tolerance for floating comparisons.
    const EPS: f64 = 1e-9;

    /// End-to-end: a 50 mm-bore single-acting cylinder running at 6 bar
    /// gauge over a standard atmosphere. Walk force -> compression ratio
    /// -> consumption -> choked-flow check and confirm each step.
    #[test]
    fn cylinder_sizing_end_to_end() {
        let p_gauge = 600_000.0; // 6 bar gauge
        let p_atm = 100_000.0; // round atmosphere -> exact r = 7
        let bore = 0.05; // 50 mm
        let stroke = 0.10; // 100 mm
        let cyl = Cylinder::single_acting(bore).unwrap();

        // 1. Force F = p_gauge * A.
        let f = cyl.extend_force(p_gauge).unwrap();
        assert!((f - p_gauge * cyl.bore_area()).abs() < 1e-6);

        // 2. Compression ratio r = 1 + p_gauge/p_atm = 7.
        let r = compression_ratio(p_gauge, p_atm).unwrap();
        assert!((r - 7.0).abs() < EPS);

        // 3. Consumption over 200 cycles = A*L*r*N.
        let q = free_air_consumption(&cyl, stroke, Action::SingleActing, 200.0, p_gauge, p_atm)
            .unwrap();
        assert!((q - cyl.bore_area() * stroke * r * 200.0).abs() < 1e-9);

        // 4. Discharging the supply (7 bar abs) to atmosphere (1 bar abs)
        //    is choked: ratio 1/7 = 0.143 < 0.528.
        let up = absolute_pressure(p_gauge, p_atm).unwrap();
        assert!(is_choked_air(up, p_atm).unwrap());
    }

    /// The re-exported symbols resolve at the crate root (smoke test for
    /// the public surface).
    #[test]
    fn public_reexports_are_reachable() {
        assert!((GAMMA_AIR - 1.4).abs() < EPS);
        assert!((STANDARD_ATMOSPHERE_PA - 101_325.0).abs() < EPS);
        assert!(
            (CRITICAL_PRESSURE_RATIO_AIR - critical_pressure_ratio(GAMMA_AIR).unwrap()).abs() < EPS
        );
        let _: Stroke = Stroke::Extend;
        let _: Action = Action::DoubleActing;
        let _: PneumaticsError = PneumaticsError::Geometry("x");
        let _r: Result<f64> = PneumaticsError::positive("v", 1.0);
    }
}
