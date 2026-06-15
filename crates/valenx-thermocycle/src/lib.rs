//! # valenx-thermocycle
//!
//! Closed-form thermal efficiencies of the classic thermodynamic power
//! cycles: the **Carnot** bound, the air-standard **Otto**, **Diesel**
//! and **Brayton / Joule** cycles, and a basic ideal **Rankine** steam
//! cycle from supplied enthalpies.
//!
//! ## What
//!
//! Each cycle is a small validated value type with an `efficiency()`
//! method (plus a free `*_efficiency` convenience function):
//!
//! - [`Carnot`] — the reversible upper bound, `η = 1 - T_c / T_h`, from
//!   two reservoir temperatures in kelvin. Also exposes refrigerator and
//!   heat-pump coefficients of performance.
//! - [`Otto`] — the ideal spark-ignition engine,
//!   `η = 1 - 1 / r^(γ-1)`, rising with compression ratio `r`.
//! - [`Diesel`] — the ideal compression-ignition engine, the Otto form
//!   times a cutoff-ratio penalty factor.
//! - [`Brayton`] — the ideal gas turbine / jet,
//!   `η = 1 - 1 / r_p^((γ-1)/γ)`, rising with pressure ratio `r_p`.
//! - [`Rankine`] — the ideal steam cycle, `η = w_net / q_in` from the
//!   four corner-state enthalpies.
//!
//! The working-fluid heat-capacity ratio `γ = c_p / c_v` is carried by a
//! validated [`HeatCapacityRatio`] newtype with [air][HeatCapacityRatio::air],
//! [monatomic][HeatCapacityRatio::monatomic] and
//! [diatomic][HeatCapacityRatio::diatomic] presets.
//!
//! ```
//! use valenx_thermocycle::{Carnot, Otto, Brayton, HeatCapacityRatio};
//!
//! // A reversible engine between 300 K and 900 K.
//! let carnot = Carnot::new(300.0, 900.0).unwrap();
//! assert!((carnot.efficiency() - (1.0 - 300.0 / 900.0)).abs() < 1e-12);
//!
//! // A petrol engine at compression ratio 8 on air.
//! let otto = Otto::with_air(8.0).unwrap();
//! assert!(otto.efficiency() > 0.0 && otto.efficiency() < 1.0);
//!
//! // A higher compression ratio is always more efficient.
//! let otto_hi = Otto::with_air(10.0).unwrap();
//! assert!(otto_hi.efficiency() > otto.efficiency());
//!
//! // A gas turbine improves with pressure ratio.
//! let g = HeatCapacityRatio::air();
//! let lo = Brayton::new(8.0, g).unwrap().efficiency();
//! let hi = Brayton::new(16.0, g).unwrap().efficiency();
//! assert!(hi > lo);
//! ```
//!
//! ## Model
//!
//! The gas cycles use the standard **air-standard** assumptions: the
//! working fluid is a fixed mass of ideal gas with *constant* specific
//! heats, every process is internally reversible, and the combustion and
//! exhaust strokes are replaced by external constant-volume or
//! constant-pressure heat transfer. Under those assumptions the Otto,
//! Diesel and Brayton efficiencies reduce to the closed forms above,
//! functions only of a compression / pressure ratio and `γ`. The Rankine
//! cycle's working fluid changes phase, so there is no single-ratio
//! closed form; this crate computes its efficiency directly from the four
//! state enthalpies the caller supplies (typically steam-table lookups),
//! `η = ((h3 - h4) - (h2 - h1)) / (h3 - h2)`.
//!
//! All constructors validate their inputs and return
//! [`Result<_, CycleError>`](error::CycleError): temperatures must be
//! positive and correctly ordered, `γ` must exceed one, compression /
//! pressure ratios must exceed one, and the Rankine boiler heat input
//! must be positive. The error type carries stable
//! [`code`](error::CycleError::code) and
//! [`category`](error::CycleError::category) accessors.
//!
//! Two cross-cycle helpers are provided:
//! [`carnot_upper_bound_holds`] checks the Carnot-bound inequality for a
//! cycle running between a given temperature pair, and [`CycleKind`]
//! labels a result for reporting.
//!
//! ## Honest scope
//!
//! Research / educational grade. These are **textbook closed-form
//! air-standard and ideal-cycle models** — the genuine equations from any
//! introductory thermodynamics course — and they reproduce the standard
//! worked examples exactly. They are deliberately idealised and are **not
//! a clinical/medical or production engineering tool**:
//!
//! - **Constant specific heats** (cold-air-standard `γ`): real `c_p`,
//!   `c_v` and hence `γ` vary with temperature, which shifts the
//!   efficiency of a real engine by a few points. There is no
//!   variable-specific-heat or finite-rate combustion model here.
//! - **No irreversibilities**: friction, turbulence, heat leakage,
//!   incomplete combustion and finite-time effects are all ignored, so
//!   these are reversible-limit efficiencies, not measured ones.
//! - **No component sub-models**: compressor / turbine isentropic
//!   efficiencies, regeneration, reheat, intercooling and feedwater
//!   heating are out of scope — only the basic cycles are modelled.
//! - **Rankine is enthalpy-in, efficiency-out**: this crate does not
//!   bundle a steam-table / IAPWS-IF97 property library; the caller
//!   supplies the four state enthalpies.
//!
//! None of that makes the numbers meaningless: the Carnot bound, the
//! air-standard efficiencies and the Rankine energy balance are all real,
//! well-established results, and each omission above is a documented,
//! well-understood refinement on the path toward a fuller cycle-analysis
//! suite.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod brayton;
pub mod carnot;
pub mod diesel;
pub mod error;
pub mod gas;
pub mod otto;
pub mod rankine;

pub use brayton::{brayton_efficiency, Brayton};
pub use carnot::{carnot_efficiency, Carnot};
pub use diesel::{diesel_efficiency, Diesel};
pub use error::{CycleError, ErrorCategory, Result};
pub use gas::HeatCapacityRatio;
pub use otto::{otto_efficiency, Otto};
pub use rankine::{rankine_efficiency, Rankine};

/// A label identifying which thermodynamic cycle an efficiency came from.
///
/// Useful when collecting heterogeneous results into a comparison table.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum CycleKind {
    /// The reversible Carnot cycle (the upper bound).
    Carnot,
    /// The air-standard Otto (spark-ignition) cycle.
    Otto,
    /// The air-standard Diesel (compression-ignition) cycle.
    Diesel,
    /// The air-standard Brayton / Joule (gas-turbine) cycle.
    Brayton,
    /// The basic ideal Rankine (steam) cycle.
    Rankine,
}

impl CycleKind {
    /// A short human-readable name for the cycle.
    pub fn name(self) -> &'static str {
        match self {
            CycleKind::Carnot => "Carnot",
            CycleKind::Otto => "Otto",
            CycleKind::Diesel => "Diesel",
            CycleKind::Brayton => "Brayton",
            CycleKind::Rankine => "Rankine",
        }
    }
}

/// Check that a candidate cycle efficiency respects the Carnot bound for
/// an engine operating between absolute temperatures `t_cold` and
/// `t_hot`.
///
/// Carnot's theorem says no engine between those reservoirs can beat the
/// reversible efficiency `1 - T_c / T_h`. This returns `true` iff
/// `0 <= efficiency <= η_carnot` (within a small tolerance to absorb
/// floating-point round-off). A reversible engine sits exactly on the
/// bound; every irreversible one sits strictly below it.
///
/// # Errors
///
/// Propagates the temperature-validation errors of [`Carnot::new`]
/// (positive, finite, correctly ordered temperatures).
pub fn carnot_upper_bound_holds(efficiency: f64, t_cold: f64, t_hot: f64) -> Result<bool> {
    let bound = Carnot::new(t_cold, t_hot)?.efficiency();
    // A tiny absolute slack so an exactly-reversible cycle (whose
    // computed efficiency may differ from the bound only in the last
    // floating-point bit) is still counted as satisfying the inequality.
    const SLACK: f64 = 1e-12;
    Ok(efficiency >= -SLACK && efficiency <= bound + SLACK)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Loose epsilon for comparisons that go through several `powf`s.
    const EPS: f64 = 1e-9;

    // ----- Carnot ------------------------------------------------------

    #[test]
    fn carnot_matches_closed_form_and_is_in_unit_interval() {
        // Worked example: T_c = 300 K, T_h = 600 K -> η = 0.5 exactly.
        let c = Carnot::new(300.0, 600.0).unwrap();
        assert!((c.efficiency() - 0.5).abs() < EPS);
        assert!(c.efficiency() > 0.0 && c.efficiency() < 1.0);

        // A second independent point, 290 K / 1160 K -> 1 - 1/4 = 0.75.
        let c2 = Carnot::new(290.0, 1160.0).unwrap();
        assert!((c2.efficiency() - 0.75).abs() < EPS);
    }

    #[test]
    fn carnot_efficiency_rises_as_hot_reservoir_rises() {
        let lo = Carnot::new(300.0, 600.0).unwrap().efficiency();
        let hi = Carnot::new(300.0, 1200.0).unwrap().efficiency();
        assert!(hi > lo, "hotter source must raise the Carnot limit");
    }

    #[test]
    fn carnot_cops_satisfy_energy_conservation() {
        // COP_heat_pump = COP_refrigerator + 1 (since Q_h = Q_c + W).
        let c = Carnot::new(280.0, 320.0).unwrap();
        let cop_ref = c.cop_refrigerator();
        let cop_hp = c.cop_heat_pump();
        assert!((cop_hp - (cop_ref + 1.0)).abs() < EPS);

        // Known values for 280 K / 320 K: ref = 280/40 = 7, hp = 320/40 = 8.
        assert!((cop_ref - 7.0).abs() < EPS);
        assert!((cop_hp - 8.0).abs() < EPS);
    }

    #[test]
    fn carnot_rejects_bad_temperatures() {
        // Equal reservoirs: no temperature drop.
        let e = Carnot::new(500.0, 500.0).unwrap_err();
        assert_eq!(e.code(), "thermocycle.temperature-order");
        assert_eq!(e.category(), ErrorCategory::Input);

        // Cold above hot.
        assert!(matches!(
            Carnot::new(600.0, 500.0),
            Err(CycleError::TemperatureOrder { .. })
        ));

        // Non-positive kelvin.
        assert!(matches!(
            Carnot::new(0.0, 500.0),
            Err(CycleError::NotPositive { .. })
        ));
        assert!(matches!(
            Carnot::new(-10.0, 500.0),
            Err(CycleError::NotPositive { .. })
        ));

        // NaN.
        assert!(matches!(
            Carnot::new(f64::NAN, 500.0),
            Err(CycleError::NotFinite { .. })
        ));
    }

    // ----- Otto --------------------------------------------------------

    #[test]
    fn otto_matches_closed_form() {
        // r = 8, γ = 1.4 -> η = 1 - 8^(-0.4) = 0.5647... (classic value).
        let o = Otto::with_air(8.0).unwrap();
        let expected = 1.0 - 8.0_f64.powf(-0.4);
        assert!((o.efficiency() - expected).abs() < EPS);
        assert!((o.efficiency() - 0.564_724_7).abs() < 1e-6);
    }

    #[test]
    fn otto_efficiency_strictly_increases_with_compression_ratio() {
        let mut prev = f64::NEG_INFINITY;
        for r in [2.0_f64, 4.0, 6.0, 8.0, 10.0, 12.0, 16.0, 20.0] {
            let eta = Otto::with_air(r).unwrap().efficiency();
            assert!(eta > 0.0 && eta < 1.0, "η out of (0,1) at r = {r}");
            assert!(eta > prev, "η must rise with r; failed at r = {r}");
            prev = eta;
        }
    }

    #[test]
    fn otto_higher_gamma_is_more_efficient_at_fixed_ratio() {
        // Monatomic γ = 5/3 beats diatomic γ = 1.4 at the same r.
        let r = 9.0;
        let dia = Otto::new(r, HeatCapacityRatio::diatomic())
            .unwrap()
            .efficiency();
        let mon = Otto::new(r, HeatCapacityRatio::monatomic())
            .unwrap()
            .efficiency();
        assert!(mon > dia);
    }

    #[test]
    fn otto_rejects_non_compressing_ratio() {
        let e = Otto::with_air(1.0).unwrap_err();
        assert_eq!(e.code(), "thermocycle.ratio-too-low");
        assert!(matches!(
            Otto::with_air(0.5),
            Err(CycleError::RatioTooLow { .. })
        ));
    }

    // ----- Diesel ------------------------------------------------------

    #[test]
    fn diesel_matches_closed_form_worked_example() {
        // Cengel-style example: r = 18, r_c = 2, γ = 1.4.
        // η = 1 - (1/18^0.4) * (2^1.4 - 1) / (1.4 * (2 - 1)) = 0.6315...
        let d = Diesel::with_air(18.0, 2.0).unwrap();
        let base = 18.0_f64.powf(-0.4);
        let factor = (2.0_f64.powf(1.4) - 1.0) / (1.4 * (2.0 - 1.0));
        let expected = 1.0 - base * factor;
        assert!((d.efficiency() - expected).abs() < EPS);
        assert!((d.efficiency() - 0.631_577_5).abs() < 1e-6);
        assert!(d.efficiency() > 0.0 && d.efficiency() < 1.0);
    }

    #[test]
    fn diesel_reduces_to_otto_as_cutoff_ratio_tends_to_one() {
        // At r_c -> 1 the cutoff penalty -> 1, so Diesel == Otto.
        let r = 18.0;
        let otto = Otto::with_air(r).unwrap().efficiency();

        // Exact unity is handled analytically.
        let at_one = Diesel::with_air(r, 1.0).unwrap();
        assert!((at_one.cutoff_penalty() - 1.0).abs() < EPS);
        assert!((at_one.efficiency() - otto).abs() < EPS);

        // And the limit is approached smoothly from above.
        let near = Diesel::with_air(r, 1.000_001).unwrap();
        assert!((near.efficiency() - otto).abs() < 1e-5);
    }

    #[test]
    fn diesel_is_less_efficient_than_otto_at_same_compression_ratio() {
        // The constant-pressure heat addition costs efficiency when
        // r_c > 1, for the same compression ratio.
        let r = 18.0;
        let otto = Otto::with_air(r).unwrap().efficiency();
        for rc in [1.5_f64, 2.0, 2.5, 3.0] {
            let diesel = Diesel::with_air(r, rc).unwrap();
            assert!(diesel.cutoff_penalty() > 1.0);
            assert!(
                diesel.efficiency() < otto,
                "Diesel should trail Otto at r_c = {rc}"
            );
        }
    }

    #[test]
    fn diesel_efficiency_falls_as_cutoff_ratio_grows() {
        // Longer constant-pressure burn -> lower efficiency at fixed r.
        let r = 20.0;
        let mut prev = f64::INFINITY;
        for rc in [1.5_f64, 2.0, 2.5, 3.0, 3.5] {
            let eta = Diesel::with_air(r, rc).unwrap().efficiency();
            assert!(eta < prev, "η must fall as r_c rises; failed at {rc}");
            prev = eta;
        }
    }

    #[test]
    fn diesel_rejects_cutoff_below_one() {
        let e = Diesel::with_air(18.0, 0.9).unwrap_err();
        assert_eq!(e.code(), "thermocycle.not-positive");
        assert!(matches!(
            Diesel::with_air(1.0, 2.0),
            Err(CycleError::RatioTooLow { .. })
        ));
    }

    // ----- Brayton -----------------------------------------------------

    #[test]
    fn brayton_matches_closed_form() {
        // r_p = 10, γ = 1.4 -> η = 1 - 10^(-0.4/1.4) = 0.4821... (classic).
        let b = Brayton::with_air(10.0).unwrap();
        let expected = 1.0 - 10.0_f64.powf(-(0.4 / 1.4));
        assert!((b.efficiency() - expected).abs() < EPS);
        assert!((b.efficiency() - 0.482_052_5).abs() < 1e-6);
        assert!(b.efficiency() > 0.0 && b.efficiency() < 1.0);
    }

    #[test]
    fn brayton_efficiency_strictly_increases_with_pressure_ratio() {
        let mut prev = f64::NEG_INFINITY;
        for rp in [2.0_f64, 4.0, 8.0, 12.0, 16.0, 24.0, 32.0] {
            let eta = Brayton::with_air(rp).unwrap().efficiency();
            assert!(eta > 0.0 && eta < 1.0, "η out of (0,1) at r_p = {rp}");
            assert!(eta > prev, "η must rise with r_p; failed at {rp}");
            prev = eta;
        }
    }

    #[test]
    fn brayton_rejects_non_compressing_ratio() {
        assert!(matches!(
            Brayton::with_air(1.0),
            Err(CycleError::RatioTooLow { .. })
        ));
        assert!(matches!(
            brayton_efficiency(8.0, 1.0),
            Err(CycleError::GammaTooLow { .. })
        ));
    }

    // ----- Rankine -----------------------------------------------------

    #[test]
    fn rankine_matches_textbook_energy_balance() {
        // Standard worked example (Moran/Cengel-class), enthalpies kJ/kg:
        //   h1 = 191.8  (sat. liquid, condenser exit)
        //   h2 = 199.7  (pump exit)
        //   h3 = 3247.6 (turbine inlet)
        //   h4 = 2007.5 (turbine exit, wet steam)
        let r = Rankine::new(191.8, 199.7, 3247.6, 2007.5).unwrap();

        let w_turbine = 3247.6 - 2007.5;
        let w_pump = 199.7 - 191.8;
        let q_in = 3247.6 - 199.7;
        let expected = (w_turbine - w_pump) / q_in;

        assert!((r.turbine_work() - w_turbine).abs() < EPS);
        assert!((r.pump_work() - w_pump).abs() < EPS);
        assert!((r.heat_in() - q_in).abs() < EPS);
        assert!((r.efficiency() - expected).abs() < EPS);

        // Efficiency in range and close to the expected ~0.404.
        assert!(r.efficiency() > 0.0 && r.efficiency() < 1.0);
        assert!((r.efficiency() - 0.4043).abs() < 1e-3);
    }

    #[test]
    fn rankine_net_work_equals_heat_in_minus_heat_out() {
        // First-law consistency: w_net = q_in - q_out for any valid state
        // set. Use round numbers so the identity is exact.
        let r = Rankine::new(200.0, 210.0, 3400.0, 2200.0).unwrap();
        let identity = r.heat_in() - r.heat_out();
        assert!((r.net_work() - identity).abs() < EPS);
    }

    #[test]
    fn rankine_back_work_ratio_is_small() {
        // For a steam cycle the pump work is a small fraction of turbine
        // work — a defining practical feature of the Rankine cycle.
        let r = Rankine::new(191.8, 199.7, 3247.6, 2007.5).unwrap();
        assert!(r.back_work_ratio() > 0.0);
        assert!(r.back_work_ratio() < 0.01);
    }

    #[test]
    fn rankine_rejects_non_positive_heat_input() {
        // h3 <= h2 means no boiler heat -> undefined efficiency.
        let e = Rankine::new(100.0, 3000.0, 2000.0, 1500.0).unwrap_err();
        assert_eq!(e.code(), "thermocycle.no-heat-input");
        assert_eq!(e.category(), ErrorCategory::Domain);

        // NaN enthalpy is rejected.
        assert!(matches!(
            Rankine::new(f64::NAN, 199.7, 3247.6, 2007.5),
            Err(CycleError::NotFinite { .. })
        ));
    }

    // ----- Cross-cycle: the Carnot bound -------------------------------

    #[test]
    fn carnot_is_at_least_any_cycle_between_the_same_temperatures() {
        // Build an Otto cycle and bound it by the Carnot efficiency of an
        // engine spanning its own temperature extremes.
        //
        // For an ideal Otto cycle with intake state (T1) and after
        // isentropic compression (T2 = T1 * r^(γ-1)), peak temperature T3
        // (after heat addition), the cycle's highest temperature is T3 and
        // lowest is T1. The Carnot engine between T1 and T3 must dominate.
        let r = 8.0;
        let gamma = HeatCapacityRatio::air().value();
        let t1 = 300.0_f64; // intake, K
        let t3 = 2000.0_f64; // peak after combustion, K

        let otto_eta = Otto::with_air(r).unwrap().efficiency();
        let carnot_eta = Carnot::new(t1, t3).unwrap().efficiency();
        assert!(
            otto_eta <= carnot_eta,
            "Otto η {otto_eta} exceeded Carnot bound {carnot_eta}"
        );

        // The helper agrees.
        assert!(carnot_upper_bound_holds(otto_eta, t1, t3).unwrap());

        // Sanity: T2 < T3 so the chosen extremes are self-consistent.
        let t2 = t1 * r.powf(gamma - 1.0);
        assert!(t2 < t3);

        // The same bound holds for Brayton and Diesel between T1 and T3.
        let brayton_eta = Brayton::with_air(10.0).unwrap().efficiency();
        assert!(carnot_upper_bound_holds(brayton_eta, t1, t3).unwrap());
        let diesel_eta = Diesel::with_air(18.0, 2.0).unwrap().efficiency();
        assert!(carnot_upper_bound_holds(diesel_eta, t1, t3).unwrap());
    }

    #[test]
    fn carnot_bound_helper_rejects_a_super_carnot_efficiency() {
        // An (impossible) efficiency above the Carnot limit must fail the
        // check rather than pass it.
        let bound = Carnot::new(300.0, 900.0).unwrap().efficiency();
        assert!(!carnot_upper_bound_holds(bound + 0.05, 300.0, 900.0).unwrap());
        // Negative efficiency is also outside the valid window.
        assert!(!carnot_upper_bound_holds(-0.01, 300.0, 900.0).unwrap());
        // Exactly on the bound passes (reversible limit).
        assert!(carnot_upper_bound_holds(bound, 300.0, 900.0).unwrap());
    }

    #[test]
    fn free_functions_agree_with_methods() {
        assert!(
            (carnot_efficiency(300.0, 600.0).unwrap()
                - Carnot::new(300.0, 600.0).unwrap().efficiency())
            .abs()
                < EPS
        );
        assert!(
            (otto_efficiency(8.0, 1.4).unwrap() - Otto::with_air(8.0).unwrap().efficiency()).abs()
                < EPS
        );
        assert!(
            (diesel_efficiency(18.0, 2.0, 1.4).unwrap()
                - Diesel::with_air(18.0, 2.0).unwrap().efficiency())
            .abs()
                < EPS
        );
        assert!(
            (brayton_efficiency(10.0, 1.4).unwrap()
                - Brayton::with_air(10.0).unwrap().efficiency())
            .abs()
                < EPS
        );
        assert!(
            (rankine_efficiency(191.8, 199.7, 3247.6, 2007.5).unwrap()
                - Rankine::new(191.8, 199.7, 3247.6, 2007.5)
                    .unwrap()
                    .efficiency())
            .abs()
                < EPS
        );
    }

    #[test]
    fn cycle_kind_names_are_stable() {
        assert_eq!(CycleKind::Carnot.name(), "Carnot");
        assert_eq!(CycleKind::Otto.name(), "Otto");
        assert_eq!(CycleKind::Diesel.name(), "Diesel");
        assert_eq!(CycleKind::Brayton.name(), "Brayton");
        assert_eq!(CycleKind::Rankine.name(), "Rankine");
    }

    #[test]
    fn gamma_presets_and_validation() {
        assert!((HeatCapacityRatio::air().value() - 1.4).abs() < EPS);
        assert!((HeatCapacityRatio::monatomic().value() - 5.0 / 3.0).abs() < EPS);
        assert!((HeatCapacityRatio::diatomic().value() - 1.4).abs() < EPS);
        assert!((HeatCapacityRatio::default().value() - 1.4).abs() < EPS);

        // γ <= 1 is rejected; NaN is rejected.
        assert!(matches!(
            HeatCapacityRatio::new(1.0),
            Err(CycleError::GammaTooLow { .. })
        ));
        assert!(matches!(
            HeatCapacityRatio::new(0.5),
            Err(CycleError::GammaTooLow { .. })
        ));
        assert!(matches!(
            HeatCapacityRatio::new(f64::INFINITY),
            Err(CycleError::NotFinite { .. })
        ));
    }

    #[test]
    fn structs_round_trip_through_json() {
        // Serde derives on the public value types survive a round-trip.
        let c = Carnot::new(300.0, 900.0).unwrap();
        let back: Carnot = serde_json::from_str(&serde_json::to_string(&c).unwrap()).unwrap();
        assert_eq!(c, back);

        let d = Diesel::with_air(18.0, 2.0).unwrap();
        let back: Diesel = serde_json::from_str(&serde_json::to_string(&d).unwrap()).unwrap();
        assert_eq!(d, back);

        let r = Rankine::new(191.8, 199.7, 3247.6, 2007.5).unwrap();
        let back: Rankine = serde_json::from_str(&serde_json::to_string(&r).unwrap()).unwrap();
        assert_eq!(r, back);
    }
}
