//! # valenx-transmissionline
//!
//! A small, exact **RF transmission-line calculator** for ideal
//! lossless lines terminated in resistive loads.
//!
//! ## What
//!
//! Describe a line by its characteristic impedance `Z0` — either
//! directly ([`Line::from_z0`], the catalogued 50 Ω / 75 Ω case) or
//! from its distributed series inductance and shunt capacitance per
//! unit length ([`Line::from_lc`], `Z0 = sqrt(L / C)`). Terminate it in
//! a resistive [`Load`] and read back the reflection and the standard
//! standing-wave figures of merit:
//!
//! ```
//! use valenx_transmissionline::{Line, Load};
//!
//! let line = Line::from_z0(50.0).unwrap();
//! let r = line.reflection(Load::resistive(100.0).unwrap());
//!
//! assert!((r.gamma() - 1.0 / 3.0).abs() < 1e-12); // reflection coeff
//! assert!((r.vswr().unwrap() - 2.0).abs() < 1e-12); // standing-wave ratio
//! ```
//!
//! ## Model
//!
//! Everything is the textbook **lossless** line (series resistance
//! `R = 0`, shunt conductance `G = 0`) with a purely **real** load, so
//! every quantity is a real closed form:
//!
//! ```text
//! Z0          = sqrt(L / C)
//! gamma       = (ZL - Z0) / (ZL + Z0)            (-1 <= gamma <= 1)
//! VSWR        = (1 + |gamma|) / (1 - |gamma|)    (VSWR >= 1)
//! return loss = -20 * log10(|gamma|)             dB, >= 0
//! mismatch    = -10 * log10(1 - |gamma|^2)       dB, >= 0
//! ```
//!
//! The exact limiting cases drop straight out and are represented
//! faithfully rather than approximated:
//!
//! - **Matched** (`ZL = Z0`): `gamma = 0`, `VSWR = 1`, return loss
//!   `+infinity` (reported as `None`).
//! - **Short** (`ZL = 0`): `gamma = -1`, total reflection.
//! - **Open** (`ZL -> infinity`, via [`Load::Open`]): `gamma = +1`,
//!   total reflection; `VSWR` diverges to `+infinity` (reported as
//!   `None`).
//!
//! - [`error`] — the [`TlError`] taxonomy and validated constructors.
//! - [`mod@line`] — the [`Line`] and `Z0 = sqrt(L / C)`.
//! - [`reflection`] — [`Load`], [`Reflection`], and the figures of
//!   merit, plus the [`load_from_gamma`] inverse and the
//!   [`gamma_magnitude_from_vswr`] bench inverse `|gamma| = (S-1)/(S+1)`
//!   that recovers the reflection magnitude from a measured VSWR.
//!
//! ## Honest scope
//!
//! Research/educational grade. This is a textbook closed-form /
//! numerical model of an **ideal lossless line with a purely resistive
//! termination**; it is NOT a clinical/medical or production
//! engineering / EMC-compliance tool. In particular it deliberately
//! does **not** model:
//!
//! - **Loss** — no series `R` or shunt `G`, so no attenuation per unit
//!   length, no frequency-dependent dielectric / copper loss, and no
//!   complex `Z0`.
//! - **Complex / reactive loads** — loads are real resistances only.
//!   There is no complex reflection coefficient, no Smith-chart phase,
//!   and no frequency sweep.
//! - **Position along the line** — no input-impedance transformation
//!   `Z_in(l)`, no electrical length / phase rotation `e^{-j*beta*l}`,
//!   and no quarter-wave or stub matching.
//!
//! Within that scope the formulas are exact and validated against the
//! analytic ground truth (see the crate tests).

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod line;
pub mod reflection;

pub use error::TlError;
pub use line::Line;
pub use reflection::{gamma_magnitude_from_vswr, load_from_gamma, Load, Reflection};

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute tolerance for floating-point comparisons in tests.
    const EPS: f64 = 1e-12;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < EPS
    }

    // ----- Z0 = sqrt(L / C) -------------------------------------------------

    #[test]
    fn z0_is_sqrt_l_over_c() {
        // L = 250 nH/m, C = 100 pF/m ⇒ sqrt(2.5e3) = 50 Ω.
        let line = Line::from_lc(250e-9, 100e-12).unwrap();
        assert!(close(line.z0_ohms(), 50.0));
    }

    #[test]
    fn z0_matches_closed_form_for_random_lc() {
        // Cross-check the constructor against sqrt(L/C) directly for a
        // spread of constants spanning several decades.
        let cases = [
            (1.0e-9, 1.0e-12),
            (3.7e-7, 2.1e-11),
            (8.2e-6, 5.0e-13),
            (1.0, 1.0),
        ];
        for (l, c) in cases {
            let line = Line::from_lc(l, c).unwrap();
            assert!(close(line.z0_ohms(), (l / c).sqrt()), "L={l}, C={c}");
        }
    }

    #[test]
    fn equal_l_and_c_give_unit_impedance() {
        let line = Line::from_lc(4.2, 4.2).unwrap();
        assert!(close(line.z0_ohms(), 1.0));
    }

    #[test]
    fn from_lc_rejects_non_positive() {
        assert!(matches!(
            Line::from_lc(0.0, 1.0),
            Err(TlError::NonPositive {
                name: "inductance_per_m",
                ..
            })
        ));
        assert!(matches!(
            Line::from_lc(1.0, -1.0),
            Err(TlError::NonPositive {
                name: "capacitance_per_m",
                ..
            })
        ));
    }

    #[test]
    fn from_lc_rejects_non_finite() {
        assert!(matches!(
            Line::from_lc(f64::NAN, 1.0),
            Err(TlError::NotFinite {
                name: "inductance_per_m",
                ..
            })
        ));
        assert!(matches!(
            Line::from_lc(1.0, f64::INFINITY),
            Err(TlError::NotFinite {
                name: "capacitance_per_m",
                ..
            })
        ));
    }

    #[test]
    fn from_z0_rejects_non_positive() {
        assert!(matches!(
            Line::from_z0(0.0),
            Err(TlError::NonPositive {
                name: "z0_ohms",
                ..
            })
        ));
        assert!(matches!(
            Line::from_z0(-50.0),
            Err(TlError::NonPositive {
                name: "z0_ohms",
                ..
            })
        ));
    }

    // ----- matched load (ZL = Z0) -------------------------------------------

    #[test]
    fn matched_load_has_zero_reflection() {
        let line = Line::from_z0(50.0).unwrap();
        let r = line.reflection(Load::resistive(50.0).unwrap());
        assert!(close(r.gamma(), 0.0));
        assert!(r.is_matched());
        assert!(!r.is_total_reflection());
    }

    #[test]
    fn matched_load_has_unit_vswr() {
        let line = Line::from_z0(75.0).unwrap();
        let r = line.reflection(Load::resistive(75.0).unwrap());
        assert!(close(r.vswr().unwrap(), 1.0));
    }

    #[test]
    fn matched_load_has_infinite_return_loss() {
        // Perfect match ⇒ |gamma| = 0 ⇒ return loss = +infinity,
        // represented as None.
        let line = Line::from_z0(50.0).unwrap();
        let r = line.reflection(Load::resistive(50.0).unwrap());
        assert!(r.return_loss_db().is_none());
    }

    #[test]
    fn matched_load_has_zero_mismatch_loss() {
        let line = Line::from_z0(50.0).unwrap();
        let r = line.reflection(Load::resistive(50.0).unwrap());
        assert!(close(r.mismatch_loss_db().unwrap(), 0.0));
        assert!(close(r.power_transmitted_fraction(), 1.0));
        assert!(close(r.power_reflected_fraction(), 0.0));
    }

    // ----- open circuit (ZL -> infinity) ------------------------------------

    #[test]
    fn open_circuit_reflects_fully_in_phase() {
        let line = Line::from_z0(50.0).unwrap();
        let r = line.reflection(Load::Open);
        assert!(close(r.gamma(), 1.0));
        assert!(close(r.gamma_magnitude(), 1.0));
        assert!(r.is_total_reflection());
    }

    #[test]
    fn open_circuit_has_infinite_vswr() {
        // |gamma| = 1 ⇒ VSWR diverges ⇒ None.
        let line = Line::from_z0(50.0).unwrap();
        assert!(line.reflection(Load::Open).vswr().is_none());
    }

    #[test]
    fn open_circuit_has_zero_return_loss_and_no_load_impedance() {
        let line = Line::from_z0(50.0).unwrap();
        let r = line.reflection(Load::Open);
        // |gamma| = 1 ⇒ return loss = -20*log10(1) = 0 dB.
        assert!(close(r.return_loss_db().unwrap(), 0.0));
        // Open circuit has no finite load resistance.
        assert!(r.load_ohms().is_none());
        // All power reflected, none transmitted ⇒ mismatch loss diverges.
        assert!(close(r.power_transmitted_fraction(), 0.0));
        assert!(r.mismatch_loss_db().is_none());
    }

    // ----- short circuit (ZL = 0) -------------------------------------------

    #[test]
    fn short_circuit_reflects_fully_inverted() {
        let line = Line::from_z0(50.0).unwrap();
        let r = line.reflection(Load::short());
        // ZL = 0 ⇒ gamma = (0 - Z0)/(0 + Z0) = -1.
        assert!(close(r.gamma(), -1.0));
        assert!(close(r.gamma_magnitude(), 1.0));
        assert!(r.is_total_reflection());
    }

    #[test]
    fn short_circuit_has_infinite_vswr_and_zero_return_loss() {
        let line = Line::from_z0(50.0).unwrap();
        let r = line.reflection(Load::short());
        assert!(r.vswr().is_none());
        assert!(close(r.return_loss_db().unwrap(), 0.0));
    }

    #[test]
    fn short_via_resistive_zero_matches_short_helper() {
        let line = Line::from_z0(50.0).unwrap();
        let a = line.reflection(Load::resistive(0.0).unwrap());
        let b = line.reflection(Load::short());
        assert!(close(a.gamma(), b.gamma()));
        assert!(close(a.gamma(), -1.0));
    }

    // ----- known mismatch (ZL = 100 Ω on 50 Ω) ------------------------------

    #[test]
    fn double_impedance_gives_one_third_gamma_and_vswr_two() {
        let line = Line::from_z0(50.0).unwrap();
        let r = line.reflection(Load::resistive(100.0).unwrap());
        // gamma = (100 - 50)/(100 + 50) = 1/3.
        assert!(close(r.gamma(), 1.0 / 3.0));
        // VSWR = (1 + 1/3)/(1 - 1/3) = 2.
        assert!(close(r.vswr().unwrap(), 2.0));
    }

    #[test]
    fn half_impedance_gives_negative_one_third_gamma_and_vswr_two() {
        let line = Line::from_z0(50.0).unwrap();
        let r = line.reflection(Load::resistive(25.0).unwrap());
        // gamma = (25 - 50)/(25 + 50) = -1/3.
        assert!(close(r.gamma(), -1.0 / 3.0));
        // VSWR depends on |gamma|, so it is also 2.
        assert!(close(r.vswr().unwrap(), 2.0));
    }

    #[test]
    fn vswr_two_has_known_return_loss() {
        // |gamma| = 1/3 ⇒ return loss = -20*log10(1/3) ≈ 9.5424 dB.
        let line = Line::from_z0(50.0).unwrap();
        let r = line.reflection(Load::resistive(100.0).unwrap());
        let expected = -20.0 * (1.0_f64 / 3.0).log10();
        assert!(close(r.return_loss_db().unwrap(), expected));
        assert!((r.return_loss_db().unwrap() - 9.542_425_094_393_249).abs() < 1e-9);
    }

    #[test]
    fn mismatch_loss_matches_closed_form() {
        // |gamma| = 1/3 ⇒ mismatch = -10*log10(1 - 1/9) = -10*log10(8/9).
        let line = Line::from_z0(50.0).unwrap();
        let r = line.reflection(Load::resistive(100.0).unwrap());
        let expected = -10.0 * (8.0_f64 / 9.0).log10();
        assert!(close(r.mismatch_loss_db().unwrap(), expected));
    }

    // ----- invariants over a sweep ------------------------------------------

    #[test]
    fn vswr_is_always_at_least_one() {
        let line = Line::from_z0(50.0).unwrap();
        // Sweep resistive loads from near-short to far-open.
        for zl_milli in 1..=500_000u64 {
            let zl = zl_milli as f64 / 1000.0; // 0.001 .. 500 Ω
            if let Some(vswr) = line.reflection(Load::resistive(zl).unwrap()).vswr() {
                assert!(vswr >= 1.0 - EPS, "VSWR {vswr} < 1 at ZL={zl}");
            }
        }
    }

    #[test]
    fn gamma_is_bounded_by_unity_for_passive_loads() {
        let line = Line::from_z0(50.0).unwrap();
        for zl_milli in 0..=500_000u64 {
            let zl = zl_milli as f64 / 1000.0;
            let mag = line
                .reflection(Load::resistive(zl).unwrap())
                .gamma_magnitude();
            assert!(mag <= 1.0 + EPS, "|gamma| {mag} > 1 at ZL={zl}");
        }
    }

    #[test]
    fn power_fractions_sum_to_one() {
        let line = Line::from_z0(50.0).unwrap();
        for zl in [10.0, 25.0, 50.0, 75.0, 100.0, 300.0] {
            let r = line.reflection(Load::resistive(zl).unwrap());
            let total = r.power_reflected_fraction() + r.power_transmitted_fraction();
            assert!(close(total, 1.0), "power split != 1 at ZL={zl}");
        }
    }

    #[test]
    fn z0_independent_match_holds_for_75_ohm_line() {
        // A 75 Ω line matched by a 75 Ω load behaves identically to the
        // 50 Ω case — confirms Z0 is threaded correctly, not hard-coded.
        let line = Line::from_z0(75.0).unwrap();
        let r = line.reflection(Load::resistive(75.0).unwrap());
        assert!(close(r.gamma(), 0.0));
        assert!(close(r.vswr().unwrap(), 1.0));
        assert!(r.return_loss_db().is_none());
        assert_eq!(r.z0_ohms(), Some(75.0));
        assert_eq!(r.load_ohms(), Some(75.0));
    }

    // ----- inverse: load_from_gamma -----------------------------------------

    #[test]
    fn load_from_gamma_inverts_reflection() {
        let line = Line::from_z0(50.0).unwrap();
        for zl in [1.0, 25.0, 50.0, 100.0, 250.0] {
            let g = line.reflection(Load::resistive(zl).unwrap()).gamma();
            let recovered = load_from_gamma(50.0, g).unwrap();
            assert!(
                (recovered - zl).abs() < 1e-9,
                "round-trip failed at ZL={zl}"
            );
        }
    }

    #[test]
    fn load_from_gamma_known_value() {
        // gamma = 1/3 on 50 Ω ⇒ ZL = 100 Ω.
        let zl = load_from_gamma(50.0, 1.0 / 3.0).unwrap();
        assert!((zl - 100.0).abs() < 1e-9);
    }

    #[test]
    fn load_from_gamma_rejects_open_circuit_limit() {
        // gamma = 1 ⇒ infinite resistance ⇒ rejected.
        assert!(matches!(
            load_from_gamma(50.0, 1.0),
            Err(TlError::GammaOutOfRange { .. })
        ));
    }

    #[test]
    fn load_from_gamma_rejects_bad_z0() {
        assert!(matches!(
            load_from_gamma(0.0, 0.0),
            Err(TlError::NonPositive {
                name: "z0_ohms",
                ..
            })
        ));
    }

    // ----- inverse: gamma_magnitude_from_vswr -------------------------------

    #[test]
    fn gamma_magnitude_from_vswr_round_trips_with_vswr() {
        // inverse ∘ forward: |gamma| -> VSWR -> |gamma|.
        for g in [0.0_f64, 1.0 / 3.0, 0.5, 0.9] {
            let s = Reflection::from_gamma(g).unwrap().vswr().unwrap();
            let back = gamma_magnitude_from_vswr(s).unwrap();
            assert!(close(back, g), "g={g} -> S={s} -> {back}");
        }
        // forward ∘ inverse: VSWR -> |gamma| -> VSWR.
        for s in [1.0_f64, 1.5, 2.0, 10.0] {
            let mag = gamma_magnitude_from_vswr(s).unwrap();
            let s_back = Reflection::from_gamma(mag).unwrap().vswr().unwrap();
            assert!(close(s_back, s), "S={s} -> mag={mag} -> {s_back}");
        }
    }

    #[test]
    fn gamma_magnitude_from_vswr_known_value() {
        // VSWR = 2 ⇒ |gamma| = 1/3 (the 100 Ω on 50 Ω case); VSWR = 1 ⇒ 0.
        assert!(close(gamma_magnitude_from_vswr(2.0).unwrap(), 1.0 / 3.0));
        assert!(close(gamma_magnitude_from_vswr(1.0).unwrap(), 0.0));
    }

    #[test]
    fn gamma_magnitude_from_vswr_matches_line_measurement() {
        // Cross-check vs the full line->load->reflection path: 100 Ω on a
        // 50 Ω line gives VSWR 2, whose inverse recovers the same
        // |gamma| = 1/3 the forward reflection reports.
        let line = Line::from_z0(50.0).unwrap();
        let r = line.reflection(Load::resistive(100.0).unwrap());
        let mag = gamma_magnitude_from_vswr(r.vswr().unwrap()).unwrap();
        assert!(close(mag, r.gamma_magnitude()));
    }

    #[test]
    fn gamma_magnitude_from_vswr_rejects_below_one_and_non_finite() {
        assert!(matches!(
            gamma_magnitude_from_vswr(0.5),
            Err(TlError::VswrBelowOne { .. })
        ));
        assert!(matches!(
            gamma_magnitude_from_vswr(f64::NAN),
            Err(TlError::NotFinite { name: "vswr", .. })
        ));
    }

    // ----- Reflection::from_gamma -------------------------------------------

    #[test]
    fn reflection_from_measured_gamma() {
        let r = Reflection::from_gamma(0.5).unwrap();
        // VSWR = (1 + 0.5)/(1 - 0.5) = 3.
        assert!(close(r.vswr().unwrap(), 3.0));
        // No impedances are known when built from a bare gamma.
        assert!(r.z0_ohms().is_none());
        assert!(r.load_ohms().is_none());
    }

    #[test]
    fn reflection_from_gamma_rejects_out_of_range() {
        assert!(matches!(
            Reflection::from_gamma(1.5),
            Err(TlError::GammaOutOfRange { .. })
        ));
        assert!(matches!(
            Reflection::from_gamma(-1.5),
            Err(TlError::GammaOutOfRange { .. })
        ));
    }

    #[test]
    fn reflection_from_gamma_rejects_non_finite() {
        assert!(matches!(
            Reflection::from_gamma(f64::NAN),
            Err(TlError::NotFinite { name: "gamma", .. })
        ));
    }

    // ----- load constructor validation --------------------------------------

    #[test]
    fn resistive_load_rejects_negative() {
        assert!(matches!(
            Load::resistive(-1.0),
            Err(TlError::Negative {
                name: "load_ohms",
                ..
            })
        ));
    }

    #[test]
    fn resistive_load_rejects_non_finite() {
        assert!(matches!(
            Load::resistive(f64::INFINITY),
            Err(TlError::NotFinite {
                name: "load_ohms",
                ..
            })
        ));
    }

    // ----- error metadata ---------------------------------------------------

    #[test]
    fn error_codes_are_stable() {
        assert_eq!(
            Line::from_z0(-1.0).unwrap_err().code(),
            "transmissionline.non-positive"
        );
        assert_eq!(
            Load::resistive(-1.0).unwrap_err().code(),
            "transmissionline.negative"
        );
        assert_eq!(
            Reflection::from_gamma(2.0).unwrap_err().code(),
            "transmissionline.gamma-out-of-range"
        );
        assert_eq!(
            Line::from_z0(f64::NAN).unwrap_err().code(),
            "transmissionline.not-finite"
        );
    }

    // ----- serde round-trip -------------------------------------------------

    #[test]
    fn reflection_serde_round_trips() {
        let line = Line::from_z0(50.0).unwrap();
        let r = line.reflection(Load::resistive(100.0).unwrap());
        let json = serde_json::to_string(&r).unwrap();
        let back: Reflection = serde_json::from_str(&json).unwrap();
        assert_eq!(r, back);
        assert!(close(back.gamma(), 1.0 / 3.0));
    }

    #[test]
    fn line_and_load_serde_round_trip() {
        let line = Line::from_lc(250e-9, 100e-12).unwrap();
        let line_back: Line = serde_json::from_str(&serde_json::to_string(&line).unwrap()).unwrap();
        assert_eq!(line, line_back);

        let load = Load::Open;
        let load_back: Load = serde_json::from_str(&serde_json::to_string(&load).unwrap()).unwrap();
        assert_eq!(load, load_back);
    }
}
