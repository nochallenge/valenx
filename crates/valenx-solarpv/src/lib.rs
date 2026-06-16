//! # valenx-solarpv — solar photovoltaic single-diode model
//!
//! A small, self-contained library for the electrical behaviour of a
//! solar photovoltaic (PV) cell or module, built on the canonical
//! *single-diode* equivalent circuit. Pure scalar numerics: no external
//! processes, no datasets, no platform dependencies.
//!
//! ## What
//!
//! Given the single-diode parameters of a cell — photocurrent `Iph`,
//! diode saturation current `I0`, ideality factor `n`, cell temperature
//! `T`, and optional series `Rs` / shunt `Rsh` resistances — this crate
//! computes:
//!
//! - the **terminal current** `I(V)` and **power** `P(V) = V*I(V)` at any
//!   operating voltage ([`SingleDiode::current_at`],
//!   [`SingleDiode::power_at`]);
//! - the **short-circuit current** `Isc = I(0)` ([`SingleDiode::isc`]);
//! - the **open-circuit voltage** `Voc` where `I(Voc) = 0`
//!   ([`SingleDiode::voc`]);
//! - the **maximum power point** `(Vmp, Imp, Pmax)`
//!   ([`max_power_point`]);
//! - the **fill factor** `FF = Pmax / (Voc * Isc)` ([`fill_factor`]);
//! - the **module efficiency** `eta = Pmax / (irradiance * area)`
//!   ([`efficiency`]);
//! - the **load-line operating point** `(V, I)` where a resistive load
//!   `I = V / R` meets the I-V curve ([`operating_point_at_load`]) — the
//!   maximum power point is its matched-load (`R = Vmp / Imp`) special case.
//!
//! ## Model
//!
//! The cell obeys the five-parameter single-diode equation derived from
//! Kirchhoff's current law on a photo-current source in parallel with a
//! diode and a shunt resistance, in series with `Rs`:
//!
//! ```text
//! I = Iph - I0 * (exp(q*(V + I*Rs) / (n*k*T)) - 1) - (V + I*Rs) / Rsh
//! ```
//!
//! with thermal voltage `Vt = k*T/q` (see
//! [`constants::thermal_voltage`]). In the ideal limit `Rs = 0`,
//! `Rsh = +inf` this is explicit in `I` and yields the textbook closed
//! forms `Isc = Iph` and `Voc = n*Vt*ln(Iph/I0 + 1)`. With `Rs > 0` (or
//! finite `Rsh`) the equation is implicit, so [`SingleDiode::current_at`]
//! solves it with a damped Newton iteration whose derivative is known
//! analytically; `Voc` is then found by bisection and the maximum power
//! point by a coarse scan plus golden-section refinement.
//!
//! ## Honest scope
//!
//! This is research/educational grade. It implements textbook
//! closed-form and numerical models (the single-diode equivalent
//! circuit, an ideal-diode closed form, and a Newton-solved implicit
//! I-V curve), validated against analytic ground truth in the unit
//! tests. It is NOT a clinical/medical/production engineering tool: it
//! does not model spectral response, angle-of-incidence and soiling
//! losses, temperature-dependent parameter drift, partial shading and
//! bypass diodes, module mismatch, degradation, or maximum-power-point
//! tracker dynamics, and it must not be used for grid-interconnection
//! certification, financial energy-yield guarantees, or safety-critical
//! array sizing. Treat its numbers as illustrative single-operating-point
//! physics, not as a substitute for measured datasheets or qualified
//! engineering software.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod constants;
pub mod diode;
pub mod error;
pub mod performance;

pub use constants::{
    celsius_to_kelvin, thermal_voltage, BOLTZMANN_J_PER_K, CELSIUS_ZERO_K, ELEMENTARY_CHARGE_C,
    STC_IRRADIANCE_W_PER_M2, STC_TEMPERATURE_K,
};
pub use diode::SingleDiode;
pub use error::{ErrorCategory, Result, SolarPvError};
pub use performance::{
    efficiency, fill_factor, max_power_point, operating_point_at_load, LoadPoint, MaxPowerPoint,
};

#[cfg(test)]
mod tests {
    use super::*;

    /// End-to-end ground-truth check on a representative crystalline-
    /// silicon cell at Standard Test Conditions: the whole pipeline
    /// (Isc, Voc, MPP, FF, efficiency) produces physically consistent,
    /// mutually compatible numbers.
    #[test]
    fn datasheet_pipeline_is_self_consistent() {
        let cell = SingleDiode::new(3.8, 1.0e-9, 1.2, STC_TEMPERATURE_K, 0.01, 200.0).unwrap();

        let isc = cell.isc().unwrap();
        let voc = cell.voc().unwrap();
        let mpp = max_power_point(&cell, 600).unwrap();
        let ff = fill_factor(&cell, 600).unwrap();

        // VALIDATE: I(0) ~ Isc, near Iph.
        assert!((isc - cell.iph).abs() < 1e-2, "Isc = {isc}");
        // VALIDATE: I(Voc) = 0.
        assert!(cell.current_at(voc).unwrap().abs() < 1e-9);
        // VALIDATE: Pmax <= Voc * Isc.
        assert!(mpp.p_max <= voc * isc + 1e-9);
        // VALIDATE: FF in (0,1), typically 0.7-0.85.
        assert!((0.70..=0.86).contains(&ff), "FF = {ff}");

        // The fill factor closes the loop: Pmax == FF * Voc * Isc.
        let pmax_from_ff = ff * voc * isc;
        assert!(
            (mpp.p_max - pmax_from_ff).abs() < 1e-6,
            "Pmax={} vs FF*Voc*Isc={}",
            mpp.p_max,
            pmax_from_ff
        );

        // Efficiency formula closes too: eta == Pmax / (G * A).
        let area = 0.0243;
        let eta = efficiency(&cell, STC_IRRADIANCE_W_PER_M2, area, 600).unwrap();
        assert!(
            (eta - mpp.p_max / (STC_IRRADIANCE_W_PER_M2 * area)).abs() < 1e-12,
            "eta = {eta}"
        );
    }

    /// Raising the irradiance (modelled here as a proportional rise in
    /// `Iph`) lifts both Isc and Pmax — the expected qualitative
    /// behaviour of a PV cell.
    #[test]
    fn higher_irradiance_lifts_isc_and_pmax() {
        let dim = SingleDiode::ideal(1.9, 1.0e-9, 1.2, STC_TEMPERATURE_K).unwrap();
        let bright = SingleDiode::ideal(3.8, 1.0e-9, 1.2, STC_TEMPERATURE_K).unwrap();
        assert!(bright.isc().unwrap() > dim.isc().unwrap());
        let p_dim = max_power_point(&dim, 400).unwrap().p_max;
        let p_bright = max_power_point(&bright, 400).unwrap().p_max;
        assert!(p_bright > p_dim, "p_dim={p_dim}, p_bright={p_bright}");
    }
}
