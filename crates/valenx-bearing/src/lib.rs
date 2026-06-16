//! # valenx-bearing
//!
//! Rolling-element bearing **rating-life** calculator: from a bearing's
//! basic dynamic load rating and the load it carries, work out how long
//! it will last — in millions of revolutions and in operating hours.
//!
//! ## What
//!
//! Four textbook quantities, exactly as written in ISO 281 / ISO 76 and
//! the SKF General Catalogue, wired together with validated inputs:
//!
//! - The **dynamic equivalent load** [`EquivalentLoad`] collapses a
//!   combined radial + axial load into the single number the life
//!   formula needs: `P = X·Fr + Y·Fa`.
//! - The **basic rating life** [`RatingLife`] evaluates
//!   `L10 = (C / P)^p`, the life (in millions of revolutions) that
//!   90 % of bearings reach before fatigue, where the
//!   [exponent `p`](BearingType) is `3` for ball and `10/3` for roller
//!   bearings.
//! - That revolution life converts to **operating hours** at a shaft
//!   speed via [`RatingLife::life_hours`]:
//!   `L10h = L10 · 1e6 / (60 · n)`.
//! - The **static safety factor** [`StaticEquivalentLoad`] guards a slow
//!   or stationary bearing against brinelling (ISO 76):
//!   `s0 = C0 / P0` with `P0 = max(X0·Fr + Y0·Fa, Fr)`.
//! - The **required dynamic load rating** inverts the load-life relation
//!   for bearing *selection*: [`required_dynamic_load_rating`] gives the
//!   `C = P · L10^(1/p)` a bearing must have to reach a target life, and
//!   [`required_dynamic_load_rating_for_hours`] sizes straight from a
//!   target life in hours at a shaft speed.
//!
//! ```
//! use valenx_bearing::{BearingType, EquivalentLoad, RatingLife};
//!
//! // A deep-groove ball bearing: C = 50 kN dynamic rating, carrying
//! // Fr = 8 kN radial and Fa = 3 kN axial with X = 0.56, Y = 1.6.
//! let load = EquivalentLoad::new(8000.0, 3000.0, 0.56, 1.6).unwrap();
//! assert!((load.value() - 9280.0).abs() < 1e-9); // P = 9.28 kN
//!
//! let life = RatingLife::from_equivalent_load(50_000.0, &load, BearingType::Ball).unwrap();
//! let l10h = life.life_hours(1500.0).unwrap();
//! println!("L10 = {:.1} Mrev, {:.0} h at 1500 rpm", life.l10_million_revs(), l10h);
//! ```
//!
//! ## Model
//!
//! The crate implements the classic ISO 281 **basic** rating life and
//! nothing more:
//!
//! - `L10 = (C / P)^p` — the load-life relation. `C` is the basic
//!   dynamic load rating (read from the bearing's data sheet), `P` is
//!   the dynamic equivalent load, and `p` is fixed by the contact
//!   geometry (point-contact ball `p = 3`, line-contact roller
//!   `p = 10/3`).
//! - `P = X·Fr + Y·Fa` — the dynamic equivalent radial load. The `X`
//!   and `Y` factors are *inputs*: they come from the bearing
//!   manufacturer's table (they depend on the series and on the
//!   `Fa/Fr` ratio relative to the limit `e`), and this crate does not
//!   guess them for you.
//! - `L10h = L10 · 1e6 / (60 · n)` — hours from revolutions at a shaft
//!   speed `n` in rpm.
//!
//! Every public constructor validates its inputs and returns a
//! [`Result<_, BearingError>`](BearingError); the error carries stable
//! [`code`](BearingError::code) and [`category`](BearingError::category)
//! accessors.
//!
//! ## Honest scope
//!
//! This is **research/educational grade**. It evaluates the genuine,
//! textbook closed-form rating-life equations — the formulae are exact
//! and the arithmetic is checked against hand-computed ground truth in
//! the test suite — but it is deliberately only the *basic* model. It
//! is **NOT a clinical/medical tool and NOT a production engineering
//! tool**, and it does **not** replace a bearing manufacturer's rating
//! software or the judgement of a qualified engineer. In particular it
//! does **not** model:
//!
//! - the **modified rating life** `Lnm = a1 · aISO · L10` — no
//!   reliability factor `a1` for survival probabilities other than
//!   90 %, and no life-modification factor `aISO` for lubrication,
//!   contamination and the fatigue load limit `Cu`;
//! - **selection of the `X` / `Y` / `X0` / `Y0` / `e` factors** (the
//!   static safety factor `s0 = C0 / P0` itself *is* provided, but you
//!   still supply the `X0` / `Y0` table values), the limiting / reference
//!   speed, friction, heat, or lubricant film thickness;
//! - any **temperature, misalignment, preload, or mounting** effects.
//!
//! Use it to learn and to sanity-check, not to certify a design.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod bearing;
pub mod error;
pub mod life;
pub mod load;
pub mod static_load;

pub use bearing::BearingType;
pub use error::{BearingError, ErrorCategory};
pub use life::{
    l10_million_revs, life_hours_from_revs, required_dynamic_load_rating,
    required_dynamic_load_rating_for_hours, RatingLife,
};
pub use load::EquivalentLoad;
pub use static_load::{static_safety_factor, StaticEquivalentLoad};

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute tolerance for floating-point comparisons in tests.
    const EPS: f64 = 1e-9;

    // ----- L10 = (C / P)^p ground truth -------------------------------

    #[test]
    fn l10_ball_exact_cube() {
        // C/P = 5, ball (p = 3): L10 = 5^3 = 125 Mrev exactly.
        let life = RatingLife::new(50_000.0, 10_000.0, BearingType::Ball).unwrap();
        assert!((life.l10_million_revs() - 125.0).abs() < EPS);
    }

    #[test]
    fn l10_roller_exact_power() {
        // C/P = 2, roller (p = 10/3): L10 = 2^(10/3).
        let life = RatingLife::new(20_000.0, 10_000.0, BearingType::Roller).unwrap();
        let expected = 2.0_f64.powf(10.0 / 3.0);
        assert!((life.l10_million_revs() - expected).abs() < EPS);
    }

    #[test]
    fn l10_unit_ratio_is_one_million_rev() {
        // By definition C is the load giving L10 = 1 Mrev, so C == P
        // gives L10 = 1 for either bearing type.
        for bt in [BearingType::Ball, BearingType::Roller] {
            let life = RatingLife::new(7000.0, 7000.0, bt).unwrap();
            assert!((life.l10_million_revs() - 1.0).abs() < EPS);
        }
    }

    #[test]
    fn free_function_matches_struct() {
        let via_fn = l10_million_revs(40_000.0, 8000.0, BearingType::Ball).unwrap();
        let via_struct = RatingLife::new(40_000.0, 8000.0, BearingType::Ball)
            .unwrap()
            .l10_million_revs();
        assert!((via_fn - via_struct).abs() < EPS);
    }

    // ----- Doubling the load cuts ball life ~8x -----------------------

    #[test]
    fn doubling_load_cuts_ball_life_eightfold() {
        // Ball: L10 ∝ P^-3, so doubling P multiplies life by 2^-3 = 1/8.
        let base = RatingLife::new(60_000.0, 10_000.0, BearingType::Ball).unwrap();
        let doubled = RatingLife::new(60_000.0, 20_000.0, BearingType::Ball).unwrap();
        let ratio = base.l10_million_revs() / doubled.l10_million_revs();
        assert!((ratio - 8.0).abs() < 1e-6, "ratio was {ratio}");
    }

    #[test]
    fn doubling_load_cuts_roller_life_by_2_pow_10_3() {
        // Roller: L10 ∝ P^-(10/3), so doubling P divides life by
        // 2^(10/3) ≈ 10.0794, NOT 8 — the roller exponent is steeper.
        let base = RatingLife::new(60_000.0, 10_000.0, BearingType::Roller).unwrap();
        let doubled = RatingLife::new(60_000.0, 20_000.0, BearingType::Roller).unwrap();
        let ratio = base.l10_million_revs() / doubled.l10_million_revs();
        let expected = 2.0_f64.powf(10.0 / 3.0);
        assert!((ratio - expected).abs() < 1e-6, "ratio was {ratio}");
        // Sanity: the roller drop is larger than the ball's 8x.
        assert!(ratio > 8.0);
    }

    #[test]
    fn halving_load_octuples_ball_life() {
        // The other direction: halving P multiplies ball life by 8.
        let base = RatingLife::new(60_000.0, 10_000.0, BearingType::Ball).unwrap();
        let halved = RatingLife::new(60_000.0, 5_000.0, BearingType::Ball).unwrap();
        let ratio = halved.l10_million_revs() / base.l10_million_revs();
        assert!((ratio - 8.0).abs() < 1e-6, "ratio was {ratio}");
    }

    // ----- Life-hours formula -----------------------------------------

    #[test]
    fn life_hours_exact_value() {
        // 125 Mrev at 1500 rpm: 125e6 / (60 * 1500) = 1388.888... h.
        let life = RatingLife::new(50_000.0, 10_000.0, BearingType::Ball).unwrap();
        let hours = life.life_hours(1500.0).unwrap();
        let expected = 125.0 * 1.0e6 / (60.0 * 1500.0);
        assert!((hours - expected).abs() < 1e-6);
        assert!((expected - 1_388.888_888_888_889).abs() < 1e-6);
    }

    #[test]
    fn life_hours_free_function_matches_method() {
        let life = RatingLife::new(50_000.0, 10_000.0, BearingType::Roller).unwrap();
        let via_method = life.life_hours(2000.0).unwrap();
        let via_fn = life_hours_from_revs(life.l10_million_revs(), 2000.0).unwrap();
        assert!((via_method - via_fn).abs() < 1e-6);
    }

    // ----- Higher rpm => fewer hours for the same revolutions ---------

    #[test]
    fn doubling_rpm_halves_hours() {
        // Same revolution life; twice the speed reaches it in half the
        // time, so hours halve.
        let life = RatingLife::new(50_000.0, 10_000.0, BearingType::Ball).unwrap();
        let slow = life.life_hours(1000.0).unwrap();
        let fast = life.life_hours(2000.0).unwrap();
        assert!((slow / fast - 2.0).abs() < 1e-9, "ratio {}", slow / fast);
        assert!(fast < slow);
    }

    #[test]
    fn hours_inversely_proportional_to_rpm() {
        // L10h * rpm is constant for a fixed L10 (= L10 * 1e6 / 60).
        let life = RatingLife::new(30_000.0, 6_000.0, BearingType::Ball).unwrap();
        let product_a = life.life_hours(900.0).unwrap() * 900.0;
        let product_b = life.life_hours(3600.0).unwrap() * 3600.0;
        assert!((product_a - product_b).abs() < 1e-3);
        let expected = life.l10_million_revs() * 1.0e6 / 60.0;
        assert!((product_a - expected).abs() < 1e-3);
    }

    // ----- Equivalent load combines radial + axial --------------------

    #[test]
    fn equivalent_load_linear_combination() {
        // P = 0.56*8000 + 1.6*3000 = 4480 + 4800 = 9280 N.
        let p = EquivalentLoad::new(8000.0, 3000.0, 0.56, 1.6).unwrap();
        assert!((p.value() - 9280.0).abs() < EPS);
    }

    #[test]
    fn equivalent_load_radial_only_is_radial_force() {
        // X = 1, Y = 0, Fa = 0 => P = Fr.
        let p = EquivalentLoad::radial_only(5000.0).unwrap();
        assert!((p.value() - 5000.0).abs() < EPS);
        // And the explicit form agrees.
        let q = EquivalentLoad::new(5000.0, 1234.0, 1.0, 0.0).unwrap();
        assert!((q.value() - 5000.0).abs() < EPS); // axial ignored when Y = 0
    }

    #[test]
    fn equivalent_load_axial_contribution_adds() {
        // Increasing only the axial load raises P by Y * ΔFa.
        let low = EquivalentLoad::new(8000.0, 2000.0, 0.56, 1.6).unwrap();
        let high = EquivalentLoad::new(8000.0, 5000.0, 0.56, 1.6).unwrap();
        let delta = high.value() - low.value();
        assert!((delta - 1.6 * (5000.0 - 2000.0)).abs() < EPS);
    }

    #[test]
    fn equivalent_load_feeds_life() {
        // End-to-end: P from the factors, then L10 from C and P.
        let p = EquivalentLoad::new(8000.0, 3000.0, 0.56, 1.6).unwrap(); // 9280 N
        let life = RatingLife::from_equivalent_load(60_000.0, &p, BearingType::Ball).unwrap();
        let expected = (60_000.0_f64 / 9280.0).powf(3.0);
        assert!((life.l10_million_revs() - expected).abs() < 1e-6);
    }

    // ----- BearingType exponent ---------------------------------------

    #[test]
    fn exponents_are_iso_281() {
        assert!((BearingType::Ball.life_exponent() - 3.0).abs() < EPS);
        assert!((BearingType::Roller.life_exponent() - 10.0 / 3.0).abs() < EPS);
        // The roller exponent is the larger of the two.
        assert!(BearingType::Roller.life_exponent() > BearingType::Ball.life_exponent());
    }

    // ----- Required dynamic load rating (selection inverse) -----------

    #[test]
    fn required_rating_hand_values() {
        // Ball (p = 3): C = 10_000 · 125^(1/3) = 10_000 · 5 = 50_000 N.
        let cb = required_dynamic_load_rating(10_000.0, 125.0, BearingType::Ball).unwrap();
        assert!((cb - 50_000.0).abs() < 1e-6, "ball C = {cb}");
        // Roller (p = 10/3): target = 2^(10/3) -> C = 10_000 · 2 = 20_000 N.
        let target = 2.0_f64.powf(10.0 / 3.0);
        let cr = required_dynamic_load_rating(10_000.0, target, BearingType::Roller).unwrap();
        assert!((cr - 20_000.0).abs() < 1e-6, "roller C = {cr}");
    }

    #[test]
    fn required_rating_inverts_l10_both_directions() {
        for bt in [BearingType::Ball, BearingType::Roller] {
            // forward C -> L10 -> required C recovers C.
            let (c, p) = (48_000.0, 9_280.0);
            let l10 = l10_million_revs(c, p, bt).unwrap();
            let c_back = required_dynamic_load_rating(p, l10, bt).unwrap();
            assert!(
                (c_back / c - 1.0).abs() < 1e-9,
                "C round-trip {c_back} vs {c} for {bt:?}"
            );
            // inverse target -> required C -> L10 recovers the target.
            let target = 200.0;
            let c_req = required_dynamic_load_rating(p, target, bt).unwrap();
            let l10_back = l10_million_revs(c_req, p, bt).unwrap();
            assert!(
                (l10_back / target - 1.0).abs() < 1e-9,
                "L10 round-trip {l10_back} vs {target} for {bt:?}"
            );
        }
    }

    #[test]
    fn required_rating_sizing_closure() {
        // Size C for a target life, build a bearing with exactly that
        // rating, and confirm its life equals the target.
        let (p, target) = (7_500.0, 90.0);
        for bt in [BearingType::Ball, BearingType::Roller] {
            let c_req = required_dynamic_load_rating(p, target, bt).unwrap();
            let life = RatingLife::new(c_req, p, bt).unwrap();
            assert!(
                (life.l10_million_revs() / target - 1.0).abs() < 1e-9,
                "sized life {} vs target {target} for {bt:?}",
                life.l10_million_revs()
            );
        }
    }

    #[test]
    fn required_rating_for_hours_matches_revs_path() {
        // 1388.888... h at 1500 rpm is 125 Mrev; ball, P = 10 kN -> C = 50 kN.
        let c = required_dynamic_load_rating_for_hours(
            10_000.0,
            1_388.888_888_888_889,
            1500.0,
            BearingType::Ball,
        )
        .unwrap();
        assert!((c - 50_000.0).abs() < 1e-3, "C = {c}");
        // Must equal the revs route with the hand-converted target L10.
        let target_l10 = 1_388.888_888_888_889 * 60.0 * 1500.0 / 1.0e6;
        let via_revs =
            required_dynamic_load_rating(10_000.0, target_l10, BearingType::Ball).unwrap();
        assert!((c - via_revs).abs() < 1e-9);
        // Closure: a bearing with this C reaches ~1388.9 h at 1500 rpm.
        let life = RatingLife::new(c, 10_000.0, BearingType::Ball).unwrap();
        assert!((life.life_hours(1500.0).unwrap() - 1_388.888_888_888_889).abs() < 1e-3);
    }

    #[test]
    fn required_rating_grows_with_target_and_load() {
        let bt = BearingType::Ball;
        // More life under the same load needs a bigger bearing.
        let lo = required_dynamic_load_rating(10_000.0, 100.0, bt).unwrap();
        let hi = required_dynamic_load_rating(10_000.0, 200.0, bt).unwrap();
        assert!(hi > lo, "more life should need more C: {hi} vs {lo}");
        // A heavier load for the same life needs a proportionally bigger
        // bearing: C is linear in P.
        let light = required_dynamic_load_rating(5_000.0, 100.0, bt).unwrap();
        let heavy = required_dynamic_load_rating(10_000.0, 100.0, bt).unwrap();
        assert!(
            (heavy / light - 2.0).abs() < 1e-9,
            "C linear in P: {heavy} vs {light}"
        );
    }

    #[test]
    fn required_rating_rejects_bad_inputs() {
        assert!(required_dynamic_load_rating(0.0, 100.0, BearingType::Ball).is_err());
        assert!(required_dynamic_load_rating(-1.0, 100.0, BearingType::Ball).is_err());
        assert!(required_dynamic_load_rating(10_000.0, 0.0, BearingType::Ball).is_err());
        assert!(required_dynamic_load_rating(10_000.0, -5.0, BearingType::Ball).is_err());
        assert!(required_dynamic_load_rating(f64::NAN, 100.0, BearingType::Ball).is_err());
        assert!(required_dynamic_load_rating(10_000.0, f64::INFINITY, BearingType::Ball).is_err());
        // The hours variant also guards hours and rpm.
        assert!(
            required_dynamic_load_rating_for_hours(10_000.0, 0.0, 1500.0, BearingType::Ball)
                .is_err()
        );
        assert!(
            required_dynamic_load_rating_for_hours(10_000.0, 1000.0, 0.0, BearingType::Ball)
                .is_err()
        );
        assert!(required_dynamic_load_rating_for_hours(
            10_000.0,
            1000.0,
            f64::NAN,
            BearingType::Ball
        )
        .is_err());
    }

    // ----- Validation / error behaviour -------------------------------

    #[test]
    fn rejects_non_positive_rating_and_load() {
        assert!(RatingLife::new(0.0, 10_000.0, BearingType::Ball).is_err());
        assert!(RatingLife::new(-1.0, 10_000.0, BearingType::Ball).is_err());
        assert!(RatingLife::new(50_000.0, 0.0, BearingType::Ball).is_err());
        assert!(RatingLife::new(50_000.0, -5.0, BearingType::Ball).is_err());
    }

    #[test]
    fn rejects_non_finite_inputs() {
        assert!(RatingLife::new(f64::NAN, 1.0, BearingType::Ball).is_err());
        assert!(RatingLife::new(f64::INFINITY, 1.0, BearingType::Ball).is_err());
        let life = RatingLife::new(50_000.0, 10_000.0, BearingType::Ball).unwrap();
        assert!(life.life_hours(f64::NAN).is_err());
        assert!(life.life_hours(0.0).is_err());
        assert!(life.life_hours(-100.0).is_err());
    }

    #[test]
    fn equivalent_load_rejects_negative() {
        assert!(EquivalentLoad::new(-1.0, 0.0, 1.0, 0.0).is_err());
        assert!(EquivalentLoad::new(0.0, -1.0, 1.0, 0.0).is_err());
        assert!(EquivalentLoad::new(0.0, 0.0, -1.0, 0.0).is_err());
        assert!(EquivalentLoad::new(0.0, 0.0, 1.0, -1.0).is_err());
        // Zeros are allowed (purely radial, radial-only treatment).
        assert!(EquivalentLoad::new(0.0, 0.0, 0.0, 0.0).is_ok());
    }

    #[test]
    fn free_helpers_reject_bad_inputs() {
        assert!(l10_million_revs(0.0, 1.0, BearingType::Ball).is_err());
        assert!(life_hours_from_revs(0.0, 100.0).is_err());
        assert!(life_hours_from_revs(100.0, 0.0).is_err());
        assert!(life_hours_from_revs(f64::INFINITY, 100.0).is_err());
    }

    #[test]
    fn error_code_and_category_are_stable() {
        let err = RatingLife::new(-1.0, 10_000.0, BearingType::Ball).unwrap_err();
        assert_eq!(err.code(), "bearing.invalid-parameter");
        assert_eq!(err.category(), ErrorCategory::Input);

        let nf = RatingLife::new(f64::NAN, 1.0, BearingType::Ball).unwrap_err();
        assert_eq!(nf.code(), "bearing.not-finite");
        assert_eq!(nf.category(), ErrorCategory::Input);
    }

    // ----- Serde round-trips ------------------------------------------

    #[test]
    fn rating_life_serde_round_trip() {
        let life = RatingLife::new(50_000.0, 10_000.0, BearingType::Roller).unwrap();
        let json = serde_json::to_string(&life).unwrap();
        let back: RatingLife = serde_json::from_str(&json).unwrap();
        assert_eq!(life, back);
    }

    #[test]
    fn equivalent_load_serde_round_trip() {
        let p = EquivalentLoad::new(8000.0, 3000.0, 0.56, 1.6).unwrap();
        let json = serde_json::to_string(&p).unwrap();
        let back: EquivalentLoad = serde_json::from_str(&json).unwrap();
        assert_eq!(p, back);
    }
}
