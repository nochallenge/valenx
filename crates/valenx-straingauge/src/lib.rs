//! # valenx-straingauge
//!
//! Closed-form **strain-gauge and Wheatstone-bridge** calculator: turn a
//! mechanical strain into the resistance change a bonded foil gauge
//! reports, the voltage a quarter / half / full bridge puts out, and the
//! uniaxial stress the strain implies — and run every relation backwards.
//!
//! ## What
//!
//! Three textbook relations, wired together.
//!
//! 1. **Gauge factor** ([`Gauge`]). The dimensionless gauge factor
//!    `GF = (ΔR/R)/ε` ties strain to fractional resistance change:
//!    `ΔR/R = GF · ε`. [`Gauge::resistance_change_ohm`] scales that to an
//!    absolute `ΔR` for a nominal 120 Ω / 350 Ω gauge.
//! 2. **Wheatstone bridge** ([`Bridge`] / [`BridgeConfig`]). The
//!    linearised bridge equation `Vout/Vin = (N/4) · GF · ε` for `N = 1`
//!    (quarter), `N = 2` (half) and `N = 4` (full) active gauges. The
//!    half bridge doubles, and the full bridge quadruples, the
//!    quarter-bridge output.
//! 3. **Hooke's law** ([`stress`]). The uniaxial stress `σ = E · ε` from
//!    the measured strain and the material's Young's modulus.
//!
//! Each forward relation has a matching inverse
//! ([`Gauge::strain_from_fractional_resistance_change`],
//! [`Bridge::strain_from_output_ratio`], [`strain_from_stress`]) so a
//! measured `ΔR/R`, `Vout/Vin` or `σ` can be turned back into strain.
//!
//! ```
//! use valenx_straingauge::{Bridge, BridgeConfig, Gauge, stress};
//!
//! let gauge = Gauge::new(2.0).expect("GF > 0");
//! let bridge = Bridge::new(gauge, BridgeConfig::Quarter);
//!
//! let strain = 1.0e-3; // 1000 µε in tension
//! let vout_ratio = bridge.output_ratio(strain).expect("finite strain");
//! let sigma_pa = stress(200.0e9, strain).expect("steel, E = 200 GPa");
//!
//! // Quarter bridge, GF = 2:  Vout/Vin = 2 · 1e-3 / 4 = 5e-4.
//! assert!((vout_ratio - 5.0e-4).abs() < 1e-15);
//! // σ = 200 GPa · 1e-3 = 200 MPa.
//! assert!((sigma_pa - 200.0e6).abs() < 1e-3);
//! ```
//!
//! ## Model
//!
//! Strain `ε` is a pure ratio (m/m), signed positive in tension and
//! negative in compression; [`microstrain`] converts from `µε`. The
//! gauge factor, Young's modulus and excitation voltage must be finite
//! and strictly positive; strain and the measured ratios may be any
//! finite value (including zero — the balanced bridge). Pressures
//! (`E`, `σ`) are unit-agnostic: keep `E` and `σ` in the same unit and
//! the relations hold. The bridge equation is taken to **first order in
//! the strain**, which is the regime where strain gauges are used.
//!
//! ## Honest scope
//!
//! This is a **research / educational-grade** calculator built on
//! textbook closed-form models. It is **not** a clinical, medical, or
//! production engineering tool and must not be used to design, qualify,
//! or certify any load-bearing structure or safety-critical system.
//! What it deliberately leaves out:
//!
//! 1. **Transverse sensitivity** — the gauge factor here is the
//!    manufacturer's axial value; the cross-axis correction is ignored.
//! 2. **Lead-wire resistance** and the three-wire desensitisation it
//!    forces on a quarter bridge.
//! 3. **Temperature effects** — apparent (thermal-output) strain,
//!    self-temperature-compensation, and the temperature dependence of
//!    `GF` and `E` are all absent; the model is isothermal.
//! 4. **Bridge non-linearity** — the exact quarter-bridge response is
//!    `Vout/Vin = (GF·ε/4)/(1 + GF·ε/2)`; here only its leading linear
//!    term is kept, which is accurate to well under a percent for the
//!    small strains gauges operate at but diverges for large `ΔR/R`.
//! 5. **Excitation self-heating**, gauge-factor and resistance
//!    tolerances, creep, hysteresis, fatigue, and rosette / principal-
//!    strain resolution.
//!
//! Treat the numbers as the ideal first-order answer for learning and
//! cross-checking, not as a substitute for a calibrated measurement
//! chain or a structural sign-off.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod bridge;
pub mod error;
pub mod gauge;

pub use bridge::{Bridge, BridgeConfig};
pub use error::{ErrorCategory, Result, StrainGaugeError};
pub use gauge::{microstrain, strain_from_stress, stress, Gauge};

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute tolerance for the dimensionless / small-magnitude
    /// quantities (resistance ratios, bridge output ratios).
    const EPS: f64 = 1.0e-15;

    // ----- gauge factor: ΔR/R = GF · ε --------------------------------

    #[test]
    fn fractional_resistance_change_is_gf_times_strain() {
        let g = Gauge::new(2.0).unwrap();
        // GF = 2, ε = 1e-3  ⇒  ΔR/R = 2e-3.
        let dr = g.fractional_resistance_change(1.0e-3).unwrap();
        assert!((dr - 2.0e-3).abs() < EPS, "got {dr}");
    }

    #[test]
    fn fractional_resistance_change_scales_with_gauge_factor() {
        // Sweep a range of gauge factors and strains; ΔR/R must equal
        // GF · ε exactly at every point.
        for &gf in &[0.5_f64, 1.0, 2.0, 2.13, 100.0] {
            let g = Gauge::new(gf).unwrap();
            for &eps in &[-2.0e-3_f64, -1.0e-4, 0.0, 7.5e-4, 5.0e-3] {
                let got = g.fractional_resistance_change(eps).unwrap();
                let want = gf * eps;
                assert!((got - want).abs() < EPS, "gf={gf} eps={eps} got={got}");
            }
        }
    }

    #[test]
    fn compression_flips_sign_of_resistance_change() {
        let g = Gauge::constantan();
        let tension = g.fractional_resistance_change(1.0e-3).unwrap();
        let compression = g.fractional_resistance_change(-1.0e-3).unwrap();
        assert!(
            (tension + compression).abs() < EPS,
            "should be equal and opposite"
        );
        assert!(tension > 0.0 && compression < 0.0);
    }

    #[test]
    fn absolute_resistance_change_scales_with_nominal_resistance() {
        let g = Gauge::constantan();
        // 350 Ω, GF = 2, ε = 1e-3  ⇒  ΔR = 350 · 2 · 1e-3 = 0.7 Ω.
        let dr = g.resistance_change_ohm(350.0, 1.0e-3).unwrap();
        assert!((dr - 0.7).abs() < 1e-12, "got {dr}");

        // 120 Ω gauge, same strain  ⇒  ΔR = 120 · 2 · 1e-3 = 0.24 Ω.
        let dr120 = g.resistance_change_ohm(120.0, 1.0e-3).unwrap();
        assert!((dr120 - 0.24).abs() < 1e-12, "got {dr120}");
    }

    #[test]
    fn resistance_inverse_round_trips() {
        let g = Gauge::new(2.05).unwrap();
        for &eps in &[-3.0e-3_f64, 0.0, 1.2e-3, 4.0e-3] {
            let frac = g.fractional_resistance_change(eps).unwrap();
            let back = g.strain_from_fractional_resistance_change(frac).unwrap();
            assert!((back - eps).abs() < EPS, "eps={eps} back={back}");
        }
    }

    // ----- quarter bridge: Vout/Vin = GF · ε / 4 ----------------------

    #[test]
    fn quarter_bridge_is_gf_strain_over_four() {
        let q = Bridge::new(Gauge::new(2.0).unwrap(), BridgeConfig::Quarter);
        // GF = 2, ε = 1e-3  ⇒  Vout/Vin = 2 · 1e-3 / 4 = 5e-4.
        let ratio = q.output_ratio(1.0e-3).unwrap();
        assert!((ratio - 5.0e-4).abs() < EPS, "got {ratio}");
    }

    #[test]
    fn quarter_bridge_matches_closed_form_over_sweep() {
        let gf = 2.11;
        let q = Bridge::new(Gauge::new(gf).unwrap(), BridgeConfig::Quarter);
        for &eps in &[-1.5e-3_f64, -2.0e-4, 0.0, 6.0e-4, 2.5e-3] {
            let got = q.output_ratio(eps).unwrap();
            let want = gf * eps / 4.0;
            assert!((got - want).abs() < EPS, "eps={eps} got={got} want={want}");
        }
    }

    // ----- full bridge = 4 × quarter; half = 2 × quarter --------------

    #[test]
    fn full_bridge_is_four_times_quarter() {
        let gauge = Gauge::new(2.07).unwrap();
        let quarter = Bridge::new(gauge, BridgeConfig::Quarter);
        let full = Bridge::new(gauge, BridgeConfig::Full);
        for &eps in &[-2.0e-3_f64, -5.0e-4, 1.0e-4, 1.8e-3] {
            let q = quarter.output_ratio(eps).unwrap();
            let f = full.output_ratio(eps).unwrap();
            assert!((f - 4.0 * q).abs() < EPS, "eps={eps} q={q} f={f}");
        }
    }

    #[test]
    fn half_bridge_is_two_times_quarter() {
        let gauge = Gauge::constantan();
        let quarter = Bridge::new(gauge, BridgeConfig::Quarter);
        let half = Bridge::new(gauge, BridgeConfig::Half);
        for &eps in &[-1.0e-3_f64, 3.0e-4, 1.1e-3] {
            let q = quarter.output_ratio(eps).unwrap();
            let h = half.output_ratio(eps).unwrap();
            assert!((h - 2.0 * q).abs() < EPS, "eps={eps} q={q} h={h}");
        }
    }

    #[test]
    fn bridge_gain_table_is_one_two_four_quarters() {
        assert!((BridgeConfig::Quarter.gain() - 0.25).abs() < EPS);
        assert!((BridgeConfig::Half.gain() - 0.5).abs() < EPS);
        assert!((BridgeConfig::Full.gain() - 1.0).abs() < EPS);
        assert_eq!(BridgeConfig::Quarter.active_arms(), 1);
        assert_eq!(BridgeConfig::Half.active_arms(), 2);
        assert_eq!(BridgeConfig::Full.active_arms(), 4);
    }

    #[test]
    fn output_voltage_is_vin_times_ratio() {
        let full = Bridge::new(Gauge::new(2.0).unwrap(), BridgeConfig::Full);
        // Vin = 5 V, GF = 2, ε = 1e-3, full bridge ⇒ Vout = 5·1·2·1e-3 = 10 mV.
        let v = full.output_voltage(5.0, 1.0e-3).unwrap();
        assert!((v - 0.01).abs() < 1e-12, "got {v}");

        // 10 V excitation doubles it.
        let v10 = full.output_voltage(10.0, 1.0e-3).unwrap();
        assert!((v10 - 0.02).abs() < 1e-12, "got {v10}");
    }

    #[test]
    fn bridge_inverse_round_trips() {
        for cfg in [
            BridgeConfig::Quarter,
            BridgeConfig::Half,
            BridgeConfig::Full,
        ] {
            let b = Bridge::new(Gauge::new(1.97).unwrap(), cfg);
            for &eps in &[-2.0e-3_f64, 0.0, 9.0e-4, 3.3e-3] {
                let ratio = b.output_ratio(eps).unwrap();
                let back = b.strain_from_output_ratio(ratio).unwrap();
                assert!(
                    (back - eps).abs() < EPS,
                    "cfg={cfg:?} eps={eps} back={back}"
                );
            }
        }
    }

    // ----- zero strain ⇒ balanced bridge (exactly 0) ------------------

    #[test]
    fn zero_strain_gives_balanced_zero_output() {
        for cfg in [
            BridgeConfig::Quarter,
            BridgeConfig::Half,
            BridgeConfig::Full,
        ] {
            for &gf in &[0.5_f64, 2.0, 50.0] {
                let b = Bridge::new(Gauge::new(gf).unwrap(), cfg);
                let ratio = b.output_ratio(0.0).unwrap();
                assert_eq!(ratio, 0.0, "balanced bridge must be exactly zero");
                // And the absolute voltage is zero for any excitation.
                let v = b.output_voltage(7.3, 0.0).unwrap();
                assert_eq!(v, 0.0);
            }
        }
    }

    #[test]
    fn zero_strain_gives_zero_resistance_change() {
        let g = Gauge::constantan();
        assert_eq!(g.fractional_resistance_change(0.0).unwrap(), 0.0);
        assert_eq!(g.resistance_change_ohm(350.0, 0.0).unwrap(), 0.0);
    }

    // ----- Hooke's law: σ = E · ε -------------------------------------

    #[test]
    fn stress_is_modulus_times_strain() {
        // Steel: E = 200 GPa, ε = 1e-3  ⇒  σ = 200 MPa.
        let sigma = stress(200.0e9, 1.0e-3).unwrap();
        assert!((sigma - 200.0e6).abs() < 1e-3, "got {sigma}");

        // Aluminium: E = 69 GPa, ε = 2e-3  ⇒  σ = 138 MPa.
        let al = stress(69.0e9, 2.0e-3).unwrap();
        assert!((al - 138.0e6).abs() < 1e-3, "got {al}");
    }

    #[test]
    fn stress_unit_agnostic_in_mpa() {
        // E in MPa, ε dimensionless  ⇒  σ in MPa.
        // 200_000 MPa · 1e-3 = 200 MPa.
        let sigma_mpa = stress(200_000.0, 1.0e-3).unwrap();
        assert!((sigma_mpa - 200.0).abs() < 1e-9, "got {sigma_mpa}");
    }

    #[test]
    fn stress_compression_is_negative() {
        let sigma = stress(200.0e9, -1.0e-3).unwrap();
        assert!((sigma + 200.0e6).abs() < 1e-3, "got {sigma}");
        assert!(sigma < 0.0);
    }

    #[test]
    fn stress_zero_strain_is_zero() {
        assert_eq!(stress(200.0e9, 0.0).unwrap(), 0.0);
    }

    #[test]
    fn stress_inverse_round_trips() {
        let e = 113.0e9; // titanium-ish
        for &eps in &[-4.0e-3_f64, 0.0, 1.0e-3, 5.0e-3] {
            let sigma = stress(e, eps).unwrap();
            let back = strain_from_stress(e, sigma).unwrap();
            assert!((back - eps).abs() < EPS, "eps={eps} back={back}");
        }
    }

    // ----- end-to-end chain: ε → ΔR/R → bridge → σ --------------------

    #[test]
    fn end_to_end_strain_to_outputs_is_consistent() {
        // A measured quarter-bridge ratio, inverted to strain, must give
        // back the same fractional resistance change and the same stress.
        let gauge = Gauge::new(2.0).unwrap();
        let bridge = Bridge::new(gauge, BridgeConfig::Quarter);
        let e = 200.0e9;

        let eps_true = 8.0e-4;
        let ratio = bridge.output_ratio(eps_true).unwrap();
        // Quarter bridge: ratio = 2 · 8e-4 / 4 = 4e-4.
        assert!((ratio - 4.0e-4).abs() < EPS, "got {ratio}");

        let eps_back = bridge.strain_from_output_ratio(ratio).unwrap();
        assert!((eps_back - eps_true).abs() < EPS);

        let dr = gauge.fractional_resistance_change(eps_back).unwrap();
        assert!((dr - 1.6e-3).abs() < EPS, "got {dr}"); // 2 · 8e-4
        let sigma = stress(e, eps_back).unwrap();
        assert!((sigma - 160.0e6).abs() < 1e-3, "got {sigma}"); // 200e9 · 8e-4
    }

    // ----- validation / error paths -----------------------------------

    #[test]
    fn non_positive_gauge_factor_is_rejected() {
        assert!(Gauge::new(0.0).is_err());
        assert!(Gauge::new(-2.0).is_err());
        match Gauge::new(-1.0) {
            Err(StrainGaugeError::NonPositive { name, .. }) => assert_eq!(name, "gauge_factor"),
            other => panic!("expected NonPositive, got {other:?}"),
        }
    }

    #[test]
    fn non_finite_inputs_are_rejected() {
        let g = Gauge::constantan();
        assert!(g.fractional_resistance_change(f64::NAN).is_err());
        assert!(g.fractional_resistance_change(f64::INFINITY).is_err());
        assert!(stress(f64::NAN, 1e-3).is_err());
        assert!(stress(200e9, f64::NEG_INFINITY).is_err());

        let b = Bridge::new(g, BridgeConfig::Full);
        assert!(b.output_ratio(f64::NAN).is_err());
        assert!(b.output_voltage(f64::INFINITY, 1e-3).is_err());
    }

    #[test]
    fn non_positive_voltage_and_resistance_are_rejected() {
        let g = Gauge::constantan();
        assert!(g.resistance_change_ohm(0.0, 1e-3).is_err());
        assert!(g.resistance_change_ohm(-350.0, 1e-3).is_err());

        let b = Bridge::new(g, BridgeConfig::Quarter);
        assert!(b.output_voltage(0.0, 1e-3).is_err());
        assert!(b.output_voltage(-5.0, 1e-3).is_err());
    }

    #[test]
    fn error_metadata_is_stable() {
        let e = Gauge::new(-1.0).unwrap_err();
        assert_eq!(e.code(), "straingauge.non-positive");
        assert_eq!(e.category(), ErrorCategory::Domain);

        let nf = stress(f64::NAN, 0.0).unwrap_err();
        assert_eq!(nf.code(), "straingauge.non-finite");
        assert_eq!(nf.category(), ErrorCategory::NotFinite);
    }

    // ----- serde round-trips ------------------------------------------

    #[test]
    fn structs_round_trip_through_json() {
        let bridge = Bridge::new(Gauge::new(2.03).unwrap(), BridgeConfig::Half);
        let json = serde_json::to_string(&bridge).unwrap();
        let back: Bridge = serde_json::from_str(&json).unwrap();
        assert_eq!(bridge, back);
    }

    #[test]
    fn microstrain_helper_scales_by_1e_minus_6() {
        assert!((microstrain(1000.0) - 1.0e-3).abs() < EPS);
        assert!((microstrain(-250.0) + 2.5e-4).abs() < EPS);
        assert_eq!(microstrain(0.0), 0.0);
    }
}
