//! # valenx-capacitor
//!
//! Closed-form capacitor electrostatics and RC-circuit calculators.
//!
//! ## What
//!
//! A small, dependency-light library of the textbook capacitor formulas:
//!
//! - [`parallel_plate`] — parallel-plate [`capacitance`] `C = eps0 eps_r A / d`,
//!   stored [energy](parallel_plate::stored_energy) `E = 1/2 C V^2`, and stored
//!   [`charge`] `Q = C V`.
//! - [`mod@reactance`] — sinusoidal-steady-state capacitive [`reactance()`]
//!   `X_C = 1 / (2 pi f C)` (and an angular-frequency variant).
//! - [`transient`] — first-order RC [charge](transient::charging_voltage) and
//!   [discharge](transient::discharging_voltage) responses
//!   `V(t) = V0 (1 - exp(-t/RC))` / `V0 exp(-t/RC)` with the
//!   [time constant](transient::time_constant) `tau = R C`, plus their
//!   inverses [time-to-charge](transient::time_to_charge) /
//!   [time-to-discharge](transient::time_to_discharge) `t = R C ln(...)`
//!   for RC-timing design.
//! - [`network`] — ideal [`series`] and [`parallel`] combination.
//! - [`spec`] — a small serde-serialisable [`spec::ParallelPlate`]
//!   descriptor that ties the geometry to its derived quantities.
//!
//! ## Model
//!
//! Every quantity is the idealised lumped-element / infinite-parallel-plate
//! result, evaluated in SI base units (metres, farads, volts, hertz,
//! ohms, seconds, joules, coulombs). The governing equations are stated in
//! each module's documentation. All inputs are validated against their
//! physical domain (areas, gaps, frequencies, capacitances and
//! resistances must be positive; times must be non-negative) and an
//! out-of-domain input yields a typed [`CapacitorError`] rather than a
//! `NaN` or `panic`.
//!
//! ## Honest scope
//!
//! Research / educational grade. These are textbook closed-form and
//! ideal-circuit models only. They deliberately ignore fringing fields,
//! dielectric loss / non-linearity / breakdown, leakage and equivalent
//! series resistance / inductance (ESR / ESL), self-resonance, tolerance
//! spread and voltage-rating limits. The library is suitable for teaching,
//! intuition-building and order-of-magnitude estimation; it is **not** a
//! clinical / medical tool and **not** a production electrical-engineering
//! design tool.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod network;
pub mod parallel_plate;
pub mod reactance;
pub mod spec;
pub mod transient;

pub use error::{CapacitorError, ErrorCategory, Result};
pub use network::{parallel, series};
pub use parallel_plate::{capacitance, charge, stored_energy, VACUUM_PERMITTIVITY};
pub use reactance::{reactance, reactance_omega};
pub use spec::ParallelPlate;
pub use transient::{
    charging_voltage, discharging_voltage, time_constant, time_to_charge, time_to_discharge,
    CHARGE_FRACTION_ONE_TAU,
};

#[cfg(test)]
mod tests {
    use super::*;

    /// Tolerance for comparing `f64` results to closed-form references.
    const EPS: f64 = 1e-12;

    /// Capacitance scales **linearly** with plate area `A`: doubling `A`
    /// at fixed gap doubles `C`.
    #[test]
    fn capacitance_scales_linearly_with_area() {
        let c1 = capacitance(1.0, 1.0e-4, 1.0e-3).unwrap();
        let c2 = capacitance(1.0, 2.0e-4, 1.0e-3).unwrap();
        assert!((c2 - 2.0 * c1).abs() < EPS * c1.max(1.0));
    }

    /// Capacitance scales as `1/d`: doubling the gap halves `C`, and the
    /// product `C * d` is invariant.
    #[test]
    fn capacitance_scales_inversely_with_gap() {
        let c1 = capacitance(1.0, 1.0e-4, 1.0e-3).unwrap();
        let c2 = capacitance(1.0, 1.0e-4, 2.0e-3).unwrap();
        assert!((c2 - 0.5 * c1).abs() < EPS * c1.max(1.0));
        // C * d is constant = eps0 * eps_r * A.
        let cd1 = c1 * 1.0e-3;
        let cd2 = c2 * 2.0e-3;
        assert!((cd1 - cd2).abs() < EPS * cd1.max(1.0));
    }

    /// Dielectric scaling: filling with `eps_r` multiplies `C` by exactly
    /// `eps_r` relative to vacuum.
    #[test]
    fn capacitance_scales_with_relative_permittivity() {
        let vac = capacitance(1.0, 1.0e-4, 1.0e-3).unwrap();
        let diel = capacitance(4.7, 1.0e-4, 1.0e-3).unwrap();
        assert!((diel - 4.7 * vac).abs() < EPS * vac.max(1.0));
    }

    /// Absolute ground truth: vacuum, 1 m^2 plates, 1 m apart gives
    /// exactly `eps0` farads.
    #[test]
    fn capacitance_unit_geometry_equals_eps0() {
        let c = capacitance(1.0, 1.0, 1.0).unwrap();
        assert!((c - VACUUM_PERMITTIVITY).abs() < 1e-24);
    }

    /// Energy obeys `E = 1/2 C V^2` and therefore scales with `V^2`:
    /// tripling the voltage multiplies the energy by nine.
    #[test]
    fn energy_matches_half_c_v_squared() {
        let c = 100.0e-6;
        let e = stored_energy(c, 10.0).unwrap();
        assert!((e - 0.5 * c * 10.0 * 10.0).abs() < EPS);
        assert!((e - 5.0e-3).abs() < EPS);

        let e3 = stored_energy(c, 30.0).unwrap();
        assert!((e3 - 9.0 * e).abs() < EPS * e.max(1.0));
    }

    /// `Q = C V` and the energy can be recovered as `E = Q^2 / (2 C)` —
    /// an independent cross-check of the charge and energy relations.
    #[test]
    fn charge_and_energy_are_consistent() {
        let c = 2.2e-6;
        let v = 12.0;
        let q = charge(c, v).unwrap();
        assert!((q - c * v).abs() < EPS * (c * v).max(1.0));

        let e_from_v = stored_energy(c, v).unwrap();
        let e_from_q = q * q / (2.0 * c);
        assert!((e_from_v - e_from_q).abs() < EPS * e_from_v.max(1.0));
    }

    /// Reactance ground truth: `X_C = 1/(2 pi f C)`. 1 uF at 1 kHz is
    /// ~159.155 ohm, and the closed form matches to machine precision.
    #[test]
    fn reactance_matches_closed_form() {
        let f = 1.0e3;
        let c = 1.0e-6;
        let xc = reactance(f, c).unwrap();
        let expected = 1.0 / (2.0 * core::f64::consts::PI * f * c);
        assert!((xc - expected).abs() < EPS * expected);
        assert!((xc - 159.154_943_091_895_34).abs() < 1e-9);
    }

    /// Reactance is inversely proportional to frequency: a tenfold rise in
    /// `f` drops `X_C` by exactly a factor of ten.
    #[test]
    fn reactance_inverse_in_frequency() {
        let c = 1.0e-6;
        let lo = reactance(1.0e3, c).unwrap();
        let hi = reactance(1.0e4, c).unwrap();
        assert!((hi - lo / 10.0).abs() < EPS * lo);
    }

    /// The angular-frequency form agrees with the ordinary-frequency form
    /// when `omega = 2 pi f`.
    #[test]
    fn reactance_omega_agrees_with_frequency_form() {
        let f = 2.5e3;
        let c = 4.7e-9;
        let by_f = reactance(f, c).unwrap();
        let by_w = reactance_omega(2.0 * core::f64::consts::PI * f, c).unwrap();
        assert!((by_f - by_w).abs() < EPS * by_f);
    }

    /// The canonical RC fact: at `t = tau = RC` the charging voltage has
    /// reached `1 - 1/e ~= 63.2 %` of the final value.
    #[test]
    fn charge_reaches_63_percent_at_one_tau() {
        let v0 = 5.0;
        let r = 1.0e3;
        let c = 1.0e-6;
        let tau = time_constant(r, c).unwrap();
        assert!((tau - 1.0e-3).abs() < EPS);

        let v = charging_voltage(v0, r, c, tau).unwrap();
        let frac = v / v0;
        // 1 - 1/e.
        let expected = 1.0 - (-1.0f64).exp();
        assert!((frac - expected).abs() < EPS);
        // The published "63 %" figure, to three significant figures.
        assert!((frac - 0.632).abs() < 1e-3);
        // And the exported constant matches.
        assert!((frac - CHARGE_FRACTION_ONE_TAU).abs() < EPS);
    }

    /// Charging boundary conditions: `V(0) = 0` and `V(t) -> V0` for large
    /// `t`; at `5 tau` the response is within ~0.7 % of `V0`.
    #[test]
    fn charging_boundary_conditions() {
        let v0 = 5.0;
        let r = 1.0e3;
        let c = 1.0e-6;
        let tau = time_constant(r, c).unwrap();

        let v_start = charging_voltage(v0, r, c, 0.0).unwrap();
        assert!(v_start.abs() < EPS);

        let v_5tau = charging_voltage(v0, r, c, 5.0 * tau).unwrap();
        assert!((v0 - v_5tau) / v0 < 0.01);

        let v_far = charging_voltage(v0, r, c, 1000.0 * tau).unwrap();
        assert!((v_far - v0).abs() < EPS * v0);
    }

    /// Discharge ground truth: starts at `V0`, falls to `V0/e` at one time
    /// constant, and charge + discharge sum to `V0` at every instant
    /// (complementary first-order responses).
    #[test]
    fn discharge_decays_and_complements_charge() {
        let v0 = 5.0;
        let r = 4.7e3;
        let c = 2.2e-6;
        let tau = time_constant(r, c).unwrap();

        let v_start = discharging_voltage(v0, r, c, 0.0).unwrap();
        assert!((v_start - v0).abs() < EPS * v0);

        let v_tau = discharging_voltage(v0, r, c, tau).unwrap();
        assert!((v_tau - v0 * (-1.0f64).exp()).abs() < EPS * v0);

        // Charging + discharging from the same V0 sum to V0 for all t.
        for &t in &[0.0, 0.3e-3, tau, 2.5 * tau, 7.0 * tau] {
            let up = charging_voltage(v0, r, c, t).unwrap();
            let down = discharging_voltage(v0, r, c, t).unwrap();
            assert!((up + down - v0).abs() < EPS * v0);
        }
    }

    /// Parallel capacitances add; series capacitances combine reciprocally
    /// (and the series total is below the smallest branch).
    #[test]
    fn series_and_parallel_combination() {
        let caps = [1.0e-6, 2.0e-6, 3.0e-6];

        let par = parallel(&caps).unwrap();
        assert!((par - 6.0e-6).abs() < EPS * 6.0e-6);

        let ser = series(&caps).unwrap();
        let expected = 1.0 / (1.0 / 1.0e-6 + 1.0 / 2.0e-6 + 1.0 / 3.0e-6);
        assert!((ser - expected).abs() < EPS * expected);
        // Series total is strictly less than the smallest branch.
        assert!(ser < 1.0e-6);
    }

    /// Two equal capacitors: parallel doubles, series halves — the
    /// textbook special case.
    #[test]
    fn equal_pair_doubles_in_parallel_and_halves_in_series() {
        let c = 2.0e-6;
        let par = parallel(&[c, c]).unwrap();
        let ser = series(&[c, c]).unwrap();
        assert!((par - 2.0 * c).abs() < EPS * (2.0 * c));
        assert!((ser - c / 2.0).abs() < EPS * c);
    }

    /// A single-element network is the identity for both combinators.
    #[test]
    fn single_element_network_is_identity() {
        let c = 3.3e-6;
        assert!((parallel(&[c]).unwrap() - c).abs() < EPS * c);
        assert!((series(&[c]).unwrap() - c).abs() < EPS * c);
    }

    /// Domain validation: out-of-range inputs are rejected with the right
    /// error code and category rather than producing `NaN` / `inf`.
    #[test]
    fn invalid_inputs_are_rejected() {
        // Non-positive geometry / capacitance.
        assert!(capacitance(1.0, 0.0, 1.0e-3).is_err());
        assert!(capacitance(1.0, 1.0e-4, 0.0).is_err());
        // eps_r below 1 is unphysical for a passive dielectric.
        assert!(capacitance(0.5, 1.0e-4, 1.0e-3).is_err());
        // DC reactance diverges -> rejected.
        assert!(reactance(0.0, 1.0e-6).is_err());
        // Negative time is rejected; zero time is allowed.
        assert!(charging_voltage(5.0, 1.0e3, 1.0e-6, -1.0).is_err());
        assert!(charging_voltage(5.0, 1.0e3, 1.0e-6, 0.0).is_ok());
        // Empty networks are rejected with the dedicated variant.
        let err = parallel(&[]).unwrap_err();
        assert_eq!(err.code(), "capacitor.empty_network");
        assert_eq!(err.category(), ErrorCategory::Input);

        // A bad branch reports an InvalidParameter with the Input category.
        let err2 = series(&[1.0e-6, -2.0e-6]).unwrap_err();
        assert_eq!(err2.code(), "capacitor.invalid_parameter");
        assert_eq!(err2.category(), ErrorCategory::Input);
    }
}
