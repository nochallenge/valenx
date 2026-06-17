//! # valenx-shellandtube
//!
//! Thermal **sizing of shell-and-tube heat exchangers** by the classic
//! log-mean-temperature-difference (LMTD) design method.
//!
//! ## What
//!
//! Give the crate a thermal duty `Q`, an overall heat-transfer
//! coefficient `U`, an LMTD correction factor `F`, and the exchanger's
//! terminal temperatures, and it returns the heat-transfer surface area
//! required to do the job — and, for a chosen tube diameter and length,
//! the number of tubes that area implies.
//!
//! The surface is split into small, composable, individually-validated
//! pieces:
//!
//! - [`TerminalDeltas`] — the two end-of-exchanger temperature gaps,
//!   built directly or from the four stream temperatures for either
//!   counter-current or co-current flow, with the
//!   [`TerminalDeltas::lmtd`] log-mean.
//! - [`CorrectionFactor`] — a validated `F` in `(0, 1]` and the
//!   `F`-corrected effective driving temperature.
//! - [`SizingInput`] / [`size`] / [`SizingResult`] — the area solve
//!   `A = Q / (U F LMTD)`, with the rating dual [`duty_from_area`]
//!   `Q = U A F LMTD` for "what duty does this exchanger deliver?".
//! - [`TubeGeometry`] — a tube `(d, L)` and the per-tube area
//!   `pi d L`, feeding the whole-bundle [`SizingResult::tube_count`].
//!
//! ```
//! use valenx_shellandtube::{size, SizingInput, TerminalDeltas, TubeGeometry};
//!
//! // 100 kW duty, U = 500 W/m^2K, F = 0.9, terminals 30 K and 10 K.
//! let deltas = TerminalDeltas::counter_current(150.0, 90.0, 60.0, 120.0).unwrap();
//! let input = SizingInput::new(100_000.0, 500.0, 0.9, deltas).unwrap();
//! let result = size(&input);
//!
//! // Size the bundle for 19 mm OD, 4 m long tubes.
//! let tube = TubeGeometry::new(0.019, 4.0).unwrap();
//! let n = result.tube_count(&tube);
//! println!(
//!     "A = {area:.2} m^2, effective LMTD = {dt:.2} K, tubes = {n}",
//!     area = result.area_m2,
//!     dt = result.effective_lmtd_k,
//! );
//! ```
//!
//! ## Model
//!
//! The single governing relation is the integrated steady-state energy
//! balance for a two-stream exchanger,
//!
//! ```text
//! Q = U * A * F * LMTD ,   LMTD = (dt1 - dt2) / ln(dt1 / dt2)
//! ```
//!
//! solved for area, `A = Q / (U F LMTD)`. `LMTD` is the correct mean
//! driving temperature for constant `U` and constant specific heats; `F`
//! (in `(0, 1]`) corrects the ideal counter-current LMTD down to the
//! effective value a multi-pass shell-and-tube geometry actually
//! delivers, so lower `F` demands more area. Tube count comes from
//! dividing the area by a single tube's cylindrical surface `pi d L` and
//! rounding up to a whole bundle. When the two terminal differences are
//! equal the LMTD reduces, by its analytic limit, to that common value
//! (the code uses the exact limit rather than evaluating `0/0`).
//!
//! ## Honest scope
//!
//! Research/educational grade. These are **textbook closed-form,
//! steady-state** models — the LMTD design method and plain geometric
//! tube counting — and nothing more. The crate is **not** a clinical,
//! medical, or production engineering tool, and a result here is a
//! first-pass sizing estimate, not a design of record. In particular it
//! deliberately does **not**:
//!
//! - compute `U` (you supply it) — no film-coefficient correlations
//!   (Dittus-Boelter, Kern, Bell-Delaware) and no wall/fouling
//!   resistances;
//! - look `F` up from `P`/`R` charts — `F` is a validated input, not
//!   derived from a shell-and-tube pass arrangement;
//! - model pressure drop, flow-induced tube vibration, two-phase /
//!   condensing / boiling duties, transient behaviour, or partial-load
//!   off-design rating;
//! - perform any TEMA / ASME mechanical, layout, baffle, nozzle, or
//!   tube-count-vs-shell-diameter packing checks.
//!
//! Any real exchanger must be verified by the rating method and the
//! applicable mechanical and process codes before it is built.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod correction;
pub mod error;
pub mod lmtd;
pub mod sizing;

pub use correction::CorrectionFactor;
pub use error::{ErrorCategory, HxError};
pub use lmtd::{lmtd, TerminalDeltas};
pub use sizing::{duty_from_area, size, SizingInput, SizingResult, TubeGeometry};

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// Generous-but-tight tolerance for analytic float comparisons.
    const EPS: f64 = 1e-9;

    // ---- LMTD ground truth -------------------------------------------

    #[test]
    fn lmtd_matches_closed_form_hand_value() {
        // dt1 = 30, dt2 = 10 -> 20 / ln(3) = 18.2047845325...
        let v = lmtd(30.0, 10.0);
        assert!((v - 18.204_784_532_536_746).abs() < EPS, "got {v}");
    }

    #[test]
    fn lmtd_is_symmetric_in_its_terminals() {
        // The log-mean is invariant under swapping dt1 and dt2.
        let a = lmtd(30.0, 10.0);
        let b = lmtd(10.0, 30.0);
        assert!((a - b).abs() < EPS, "{a} vs {b}");
    }

    #[test]
    fn lmtd_equal_terminals_returns_common_value_not_nan() {
        // dt1 == dt2: closed form is 0/0; analytic limit is the value.
        let v = lmtd(15.0, 15.0);
        assert!(v.is_finite());
        assert!((v - 15.0).abs() < EPS, "got {v}");
    }

    #[test]
    fn lmtd_near_equal_terminals_is_continuous() {
        // Approaching equality must converge to the common value, not blow
        // up: the limit of LMTD as dt2 -> dt1 is dt1.
        let near = lmtd(20.000_000_1, 20.0);
        assert!((near - 20.0).abs() < 1e-3, "got {near}");
    }

    #[test]
    fn lmtd_is_bounded_by_arithmetic_and_geometric_means() {
        // Classic inequality: GM <= LMTD <= AM for positive terminals.
        let (dt1, dt2) = (40.0_f64, 8.0_f64);
        let l = lmtd(dt1, dt2);
        let am = 0.5 * (dt1 + dt2);
        let gm = (dt1 * dt2).sqrt();
        assert!(l <= am + EPS, "LMTD {l} should be <= AM {am}");
        assert!(l >= gm - EPS, "LMTD {l} should be >= GM {gm}");
    }

    #[test]
    fn terminal_deltas_counter_current_geometry() {
        // Th: 150->90, Tc: 60->120 counter-current.
        // dt1 = Th_in - Tc_out = 150 - 120 = 30; dt2 = 90 - 60 = 30.
        let d = TerminalDeltas::counter_current(150.0, 90.0, 60.0, 120.0).unwrap();
        assert!((d.dt1 - 30.0).abs() < EPS);
        assert!((d.dt2 - 30.0).abs() < EPS);
        // Equal terminals -> LMTD is the common value 30.
        assert!((d.lmtd() - 30.0).abs() < EPS);
    }

    #[test]
    fn terminal_deltas_co_current_geometry() {
        // Same streams, parallel flow: dt1 = 150-60 = 90, dt2 = 90-120 = -30
        // would be a cross; use a feasible co-current case instead.
        // Th: 150->100, Tc: 40->80 co-current: dt1 = 110, dt2 = 20.
        let d = TerminalDeltas::co_current(150.0, 100.0, 40.0, 80.0).unwrap();
        assert!((d.dt1 - 110.0).abs() < EPS);
        assert!((d.dt2 - 20.0).abs() < EPS);
        // 90 / ln(5.5) = 52.789...
        let expected = 90.0 / (110.0_f64 / 20.0).ln();
        assert!((d.lmtd() - expected).abs() < EPS, "got {}", d.lmtd());
    }

    // ---- Area: A = Q / (U F LMTD) ------------------------------------

    #[test]
    fn area_equals_q_over_u_f_lmtd() {
        let deltas = TerminalDeltas::new(30.0, 10.0).unwrap();
        let input = SizingInput::new(100_000.0, 500.0, 0.9, deltas).unwrap();
        let r = size(&input);

        let lmtd_expected = 20.0 / 3.0_f64.ln();
        let area_expected = 100_000.0 / (500.0 * 0.9 * lmtd_expected);

        assert!((r.lmtd_k - lmtd_expected).abs() < EPS, "lmtd {}", r.lmtd_k);
        assert!(
            (r.effective_lmtd_k - 0.9 * lmtd_expected).abs() < EPS,
            "eff {}",
            r.effective_lmtd_k
        );
        assert!(
            (r.area_m2 - area_expected).abs() < 1e-6,
            "area {}",
            r.area_m2
        );
    }

    #[test]
    fn area_back_substitutes_to_recover_duty() {
        // Round-trip: Q_recovered = U * A * F * LMTD must equal input Q.
        let deltas = TerminalDeltas::new(45.0, 12.0).unwrap();
        let input = SizingInput::new(250_000.0, 850.0, 0.78, deltas).unwrap();
        let r = size(&input);
        let q_recovered = input.u_w_per_m2k * r.area_m2 * input.correction.value() * r.lmtd_k;
        assert!(
            (q_recovered - input.duty_w).abs() < 1e-6,
            "recovered {q_recovered} vs {}",
            input.duty_w
        );
    }

    #[test]
    fn higher_duty_needs_more_area_all_else_equal() {
        let deltas = TerminalDeltas::new(30.0, 10.0).unwrap();
        let small = size(&SizingInput::new(50_000.0, 500.0, 0.9, deltas).unwrap());
        let large = size(&SizingInput::new(150_000.0, 500.0, 0.9, deltas).unwrap());
        assert!(
            large.area_m2 > small.area_m2,
            "more duty must need more area: {} !> {}",
            large.area_m2,
            small.area_m2
        );
        // Linear in Q: tripling duty triples area.
        assert!((large.area_m2 / small.area_m2 - 3.0).abs() < EPS);
    }

    #[test]
    fn higher_u_needs_less_area() {
        let deltas = TerminalDeltas::new(30.0, 10.0).unwrap();
        let low_u = size(&SizingInput::new(100_000.0, 300.0, 0.9, deltas).unwrap());
        let high_u = size(&SizingInput::new(100_000.0, 900.0, 0.9, deltas).unwrap());
        assert!(high_u.area_m2 < low_u.area_m2);
        // Inverse in U: tripling U thirds the area.
        assert!((low_u.area_m2 / high_u.area_m2 - 3.0).abs() < EPS);
    }

    // ---- Rating dual: duty_from_area ---------------------------------

    #[test]
    fn duty_from_area_inverts_size() {
        // size gives the area for a duty; duty_from_area gives the duty
        // back from that area at the same conditions.
        let deltas = TerminalDeltas::new(45.0, 12.0).unwrap();
        let input = SizingInput::new(250_000.0, 850.0, 0.78, deltas).unwrap();
        let r = size(&input);
        let q = duty_from_area(r.area_m2, input.u_w_per_m2k, input.correction, deltas).unwrap();
        assert!(
            (q - input.duty_w).abs() < 1e-6 * input.duty_w,
            "q {q} vs {}",
            input.duty_w
        );
    }

    #[test]
    fn duty_from_area_matches_energy_balance() {
        // U=500, A=10, F=0.9, deltas(30,10): LMTD=20/ln3, Q = U*A*F*LMTD.
        let deltas = TerminalDeltas::new(30.0, 10.0).unwrap();
        let f = CorrectionFactor::new(0.9).unwrap();
        let q = duty_from_area(10.0, 500.0, f, deltas).unwrap();
        let lmtd = 20.0 / 3.0_f64.ln();
        assert!((q - 500.0 * 10.0 * 0.9 * lmtd).abs() < 1e-6, "q = {q}");
    }

    #[test]
    fn duty_scales_with_area_and_u() {
        let deltas = TerminalDeltas::new(25.0, 15.0).unwrap();
        let f = CorrectionFactor::new(1.0).unwrap();
        let base = duty_from_area(5.0, 400.0, f, deltas).unwrap();
        let double_a = duty_from_area(10.0, 400.0, f, deltas).unwrap();
        let double_u = duty_from_area(5.0, 800.0, f, deltas).unwrap();
        assert!((double_a - 2.0 * base).abs() < EPS * base);
        assert!((double_u - 2.0 * base).abs() < EPS * base);
    }

    #[test]
    fn duty_from_area_rejects_bad_inputs() {
        let deltas = TerminalDeltas::new(30.0, 10.0).unwrap();
        let f = CorrectionFactor::new(0.9).unwrap();
        assert!(duty_from_area(0.0, 500.0, f, deltas).is_err());
        assert!(duty_from_area(10.0, 0.0, f, deltas).is_err());
        assert!(duty_from_area(f64::NAN, 500.0, f, deltas).is_err());
    }

    // ---- Correction factor F -----------------------------------------

    #[test]
    fn correction_factor_accepts_open_closed_unit_interval() {
        assert!(CorrectionFactor::new(0.5).is_ok());
        assert!(CorrectionFactor::new(1.0).is_ok());
        assert!(CorrectionFactor::new(1e-9).is_ok());
        assert!((CorrectionFactor::ideal().value() - 1.0).abs() < EPS);
    }

    #[test]
    fn correction_factor_rejects_out_of_range_values() {
        // Zero, negative, > 1, and non-finite are all rejected.
        for bad in [0.0, -0.1, 1.0001, f64::NAN, f64::INFINITY] {
            let err = CorrectionFactor::new(bad).unwrap_err();
            assert_eq!(err.code(), "shellandtube.correction_factor_out_of_range");
            assert_eq!(err.category(), ErrorCategory::Input);
        }
    }

    #[test]
    fn correction_factor_reduces_effective_lmtd() {
        // F in (0,1) must strictly reduce the effective LMTD below raw.
        let raw = lmtd(30.0, 10.0);
        let f = CorrectionFactor::new(0.85).unwrap();
        let eff = f.effective_lmtd(raw);
        assert!(eff < raw, "effective {eff} should be < raw {raw}");
        assert!((eff - 0.85 * raw).abs() < EPS);

        // F = 1 leaves it unchanged.
        let eff_ideal = CorrectionFactor::ideal().effective_lmtd(raw);
        assert!((eff_ideal - raw).abs() < EPS);
    }

    #[test]
    fn lower_correction_factor_needs_more_area() {
        // Since A ∝ 1/F, halving F doubles the required area.
        let deltas = TerminalDeltas::new(30.0, 10.0).unwrap();
        let good = size(&SizingInput::new(100_000.0, 500.0, 0.9, deltas).unwrap());
        let poor = size(&SizingInput::new(100_000.0, 500.0, 0.45, deltas).unwrap());
        assert!(poor.area_m2 > good.area_m2);
        assert!((poor.area_m2 / good.area_m2 - 2.0).abs() < EPS);
    }

    // ---- Tube count: n = A / (pi d L) --------------------------------

    #[test]
    fn area_per_tube_is_pi_d_l() {
        let tube = TubeGeometry::new(0.019, 4.0).unwrap();
        let expected = PI * 0.019 * 4.0;
        assert!((tube.area_per_tube_m2() - expected).abs() < EPS);
    }

    #[test]
    fn tube_count_from_area_matches_a_over_pi_d_l() {
        // Construct a case whose area is an exact multiple of pi d L so the
        // continuous tube count is a clean integer.
        let tube = TubeGeometry::new(0.020, 5.0).unwrap();
        // per_tube = pi * 0.02 * 5 = 0.1 pi m^2.
        let per_tube = tube.area_per_tube_m2();
        // Want area = 40 * per_tube exactly; choose U, F, LMTD then back out Q.
        let deltas = TerminalDeltas::new(20.0, 20.0).unwrap(); // LMTD = 20
        let f = 1.0;
        let u = 500.0;
        let target_area = 40.0 * per_tube;
        let q = target_area * u * f * deltas.lmtd();
        let r = size(&SizingInput::new(q, u, f, deltas).unwrap());

        assert!((r.area_m2 - target_area).abs() < 1e-9, "area {}", r.area_m2);
        let real = r.tubes_real(&tube);
        assert!((real - 40.0).abs() < 1e-6, "tubes_real {real}");
        assert_eq!(r.tube_count(&tube), 40);
    }

    #[test]
    fn tube_count_rounds_up_to_whole_bundle() {
        // 40.0 tubes-worth + a sliver must round up to 41 whole tubes.
        let tube = TubeGeometry::new(0.020, 5.0).unwrap();
        let per_tube = tube.area_per_tube_m2();
        let deltas = TerminalDeltas::new(20.0, 20.0).unwrap();
        let (u, f) = (500.0, 1.0);
        let target_area = 40.01 * per_tube;
        let q = target_area * u * f * deltas.lmtd();
        let r = size(&SizingInput::new(q, u, f, deltas).unwrap());
        let real = r.tubes_real(&tube);
        assert!(real > 40.0 && real < 41.0, "tubes_real {real}");
        assert_eq!(r.tube_count(&tube), 41);
    }

    #[test]
    fn longer_or_fatter_tubes_need_fewer_of_them() {
        let deltas = TerminalDeltas::new(30.0, 10.0).unwrap();
        let r = size(&SizingInput::new(500_000.0, 500.0, 0.9, deltas).unwrap());

        let short = TubeGeometry::new(0.019, 2.0).unwrap();
        let long = TubeGeometry::new(0.019, 6.0).unwrap();
        assert!(r.tubes_real(&long) < r.tubes_real(&short));
        // 3x the length -> 1/3 the tubes (continuous form).
        assert!((r.tubes_real(&short) / r.tubes_real(&long) - 3.0).abs() < 1e-6);

        let thin = TubeGeometry::new(0.010, 4.0).unwrap();
        let fat = TubeGeometry::new(0.030, 4.0).unwrap();
        assert!(r.tubes_real(&fat) < r.tubes_real(&thin));
    }

    // ---- Validation / error surface ----------------------------------

    #[test]
    fn rejects_non_positive_duty_and_u() {
        let deltas = TerminalDeltas::new(30.0, 10.0).unwrap();
        let bad_q = SizingInput::new(0.0, 500.0, 0.9, deltas).unwrap_err();
        assert_eq!(bad_q.code(), "shellandtube.bad_parameter");
        assert_eq!(bad_q.category(), ErrorCategory::Input);

        let bad_u = SizingInput::new(1000.0, -5.0, 0.9, deltas).unwrap_err();
        assert_eq!(bad_u.code(), "shellandtube.bad_parameter");
    }

    #[test]
    fn rejects_temperature_cross_as_infeasible() {
        // A negative terminal difference is a temperature cross.
        let err = TerminalDeltas::new(30.0, -5.0).unwrap_err();
        assert_eq!(err.code(), "shellandtube.infeasible_temperature_profile");
        assert_eq!(err.category(), ErrorCategory::Infeasible);
    }

    #[test]
    fn rejects_non_finite_terminal_difference() {
        let err = TerminalDeltas::new(f64::NAN, 10.0).unwrap_err();
        assert_eq!(err.code(), "shellandtube.bad_parameter");
    }

    #[test]
    fn rejects_non_positive_tube_dimensions() {
        assert_eq!(
            TubeGeometry::new(0.0, 4.0).unwrap_err().code(),
            "shellandtube.bad_parameter"
        );
        assert_eq!(
            TubeGeometry::new(0.019, -1.0).unwrap_err().code(),
            "shellandtube.bad_parameter"
        );
    }

    // ---- Serde round-trip --------------------------------------------

    #[test]
    fn result_serde_round_trips() {
        let deltas = TerminalDeltas::new(30.0, 10.0).unwrap();
        let r = size(&SizingInput::new(100_000.0, 500.0, 0.9, deltas).unwrap());
        let json = serde_json::to_string(&r).unwrap();
        let back: SizingResult = serde_json::from_str(&json).unwrap();
        assert!((back.area_m2 - r.area_m2).abs() < EPS);
        assert!((back.lmtd_k - r.lmtd_k).abs() < EPS);
        assert!((back.effective_lmtd_k - r.effective_lmtd_k).abs() < EPS);
    }
}
