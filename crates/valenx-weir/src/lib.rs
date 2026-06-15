//! # valenx-weir
//!
//! Sharp-crested **weir-flow hydraulics** for open channels: closed-form
//! discharge over a **rectangular** weir and a **triangular (V-notch)**
//! weir.
//!
//! ## What
//!
//! A weir is a notch or low dam across an open channel; the volumetric
//! flow rate `Q` (the *discharge*) over it is a known, monotone function
//! of the upstream **head** `H` — the height of the free surface above
//! the crest. That single-valued `Q(H)` relationship is what makes a
//! weir a *flow meter*: read the head, read off the flow.
//!
//! This crate provides the two textbook sharp-crested weirs:
//!
//! - [`rectangular::RectangularWeir`] — a horizontal crest of length
//!   `L`, with `Q ∝ H^(3/2)`.
//! - [`vnotch::VNotchWeir`] — a triangular notch of full vertex angle
//!   `θ`, with `Q ∝ H^(5/2)`.
//!
//! Each type is built through a validated constructor (see
//! [`error::WeirError`]) and evaluated with `.discharge(head_m)`.
//!
//! ```
//! use valenx_weir::{RectangularWeir, VNotchWeir};
//!
//! // A 2 m-wide rectangular weir, Cd = 0.62, at 0.30 m of head.
//! let rect = RectangularWeir::new(2.0, 0.62).unwrap();
//! let q_rect = rect.discharge(0.30).unwrap();
//!
//! // A 90-degree V-notch, Cd = 0.58, at 0.30 m of head.
//! let vee = VNotchWeir::ninety_degree(0.58).unwrap();
//! let q_vee = vee.discharge(0.30).unwrap();
//!
//! // At low head the V-notch passes far less flow — exactly why it is
//! // the gauging weir of choice for small discharges.
//! assert!(q_vee < q_rect);
//! ```
//!
//! ## Model
//!
//! Both formulae come from integrating the ideal spillover velocity
//! `v(z) = √(2 g z)` (Torricelli) across the weir opening and folding
//! every real-flow effect — viscosity, surface tension, the
//! contraction of the nappe, the approach velocity — into a single
//! empirical **discharge coefficient** `Cd`:
//!
//! ```text
//!   Rectangular:  Q = Cd · (2/3)  · √(2 g) · L          · H^(3/2)
//!   V-notch:      Q = Cd · (8/15) · √(2 g) · tan(θ/2)   · H^(5/2)
//! ```
//!
//! The constants `2/3` ([`rectangular::RECT_COEFFICIENT`]) and `8/15`
//! ([`vnotch::VNOTCH_COEFFICIENT`]) are the exact results of those
//! integrals over a rectangular and a triangular opening respectively.
//! All quantities are SI: `L` and `H` in metres, `g` in m·s⁻²,
//! `θ` in radians, and `Q` in m³·s⁻¹.
//!
//! The qualitative behaviour the formulae encode — and which the test
//! suite checks against hand-computed ground truth — is:
//!
//! - rectangular discharge grows as the **3/2 power** of head;
//! - V-notch discharge grows as the **5/2 power** of head;
//! - rectangular `Q` scales **linearly** with both `Cd` and `L`;
//! - more head always means more flow (strictly monotone);
//! - at the same small head the V-notch passes **less** flow, giving it
//!   finer resolution at low discharge.
//!
//! ## Honest scope
//!
//! This is a **research / educational-grade** library. It implements the
//! textbook sharp-crested closed forms only — nothing more:
//!
//! - It is **not** a clinical, medical, or production engineering tool,
//!   and is not a substitute for calibrated field gauging or a
//!   standards-validated hydraulic design code (ISO 1438, ASTM D5242,
//!   the USBR *Water Measurement Manual*, BS 3680, and the like).
//! - The discharge coefficient `Cd` is treated as a **single supplied
//!   constant**. Real weirs use head- and geometry-dependent `Cd`
//!   correlations (e.g. Kindsvater–Carter for rectangular weirs,
//!   Kindsvater–Shen for V-notches); none of those are modelled here.
//!   You must supply a `Cd` appropriate to your weir and head range.
//! - The formulae assume an **ideal sharp-crested, fully-ventilated,
//!   free-discharge** weir with negligible approach velocity. Broad-
//!   crested weirs, submerged / drowned flow, end-contraction
//!   corrections, the velocity-of-approach term, viscous and
//!   surface-tension effects at very low heads, and nappe aeration are
//!   **out of scope**.
//! - Results are deterministic evaluations of the equations above; they
//!   carry no uncertainty quantification.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod rectangular;
pub mod vnotch;

pub use error::{ErrorCategory, WeirError};
pub use rectangular::RectangularWeir;
pub use vnotch::VNotchWeir;

/// Standard gravitational acceleration `g₀ = 9.80665 m·s⁻²`.
///
/// The CGPM-defined standard gravity, used as the default `g` by the
/// [`RectangularWeir::new`] and [`VNotchWeir::new`] constructors. Supply
/// a site-specific value through the `with_gravity` constructors when a
/// different `g` is wanted.
pub const G_STANDARD: f64 = 9.806_65;

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// Absolute tolerance for float comparisons against hand-computed
    /// ground truth.
    const EPS: f64 = 1e-9;

    // ----------------------------------------------------------------
    // Rectangular weir — closed-form ground truth.
    // ----------------------------------------------------------------

    #[test]
    fn rectangular_matches_hand_computed_value() {
        // L = 2 m, Cd = 0.62, g = 9.80665, H = 0.30 m.
        // Q = 0.62 · (2/3) · sqrt(2·9.80665) · 2 · 0.30^1.5.
        let weir = RectangularWeir::new(2.0, 0.62).unwrap();
        let expected = 0.62 * (2.0 / 3.0) * (2.0 * G_STANDARD).sqrt() * 2.0 * 0.30_f64.powf(1.5);
        let got = weir.discharge(0.30).unwrap();
        assert!(
            (got - expected).abs() < EPS,
            "rectangular Q mismatch: got {got}, expected {expected}"
        );
        // Sanity on the magnitude: ~0.60 m^3/s for this 2 m-wide case.
        assert!((got - 0.601_572_0).abs() < 1e-6, "magnitude off: {got}");
    }

    #[test]
    fn rectangular_scales_as_head_to_the_three_halves() {
        // Doubling H must multiply Q by 2^1.5.
        let weir = RectangularWeir::new(1.5, 0.6).unwrap();
        let q1 = weir.discharge(0.2).unwrap();
        let q2 = weir.discharge(0.4).unwrap();
        let ratio = q2 / q1;
        assert!(
            (ratio - 2.0_f64.powf(1.5)).abs() < EPS,
            "head exponent wrong: ratio {ratio}, want {}",
            2.0_f64.powf(1.5)
        );
    }

    #[test]
    fn rectangular_scales_linearly_with_cd() {
        // Tripling Cd triples Q (everything else fixed).
        let a = RectangularWeir::new(1.0, 0.2).unwrap();
        let b = RectangularWeir::new(1.0, 0.6).unwrap();
        let qa = a.discharge(0.25).unwrap();
        let qb = b.discharge(0.25).unwrap();
        assert!(
            (qb - 3.0 * qa).abs() < EPS,
            "Cd not linear: qa {qa}, qb {qb}"
        );
    }

    #[test]
    fn rectangular_scales_linearly_with_length() {
        // Quadrupling L quadruples Q.
        let a = RectangularWeir::new(0.5, 0.61).unwrap();
        let b = RectangularWeir::new(2.0, 0.61).unwrap();
        let qa = a.discharge(0.18).unwrap();
        let qb = b.discharge(0.18).unwrap();
        assert!(
            (qb - 4.0 * qa).abs() < EPS,
            "L not linear: qa {qa}, qb {qb}"
        );
    }

    #[test]
    fn rectangular_higher_head_more_flow() {
        let weir = RectangularWeir::new(1.0, 0.62).unwrap();
        let mut prev = weir.discharge(0.05).unwrap();
        for h in [0.10, 0.20, 0.40, 0.80, 1.60] {
            let q = weir.discharge(h).unwrap();
            assert!(q > prev, "not monotone at H={h}: q {q}, prev {prev}");
            prev = q;
        }
    }

    // ----------------------------------------------------------------
    // V-notch weir — closed-form ground truth.
    // ----------------------------------------------------------------

    #[test]
    fn vnotch_matches_hand_computed_value() {
        // 90-degree notch (tan(45) = 1), Cd = 0.58, H = 0.20 m.
        // Q = 0.58 · (8/15) · sqrt(2·9.80665) · 1 · 0.20^2.5.
        let weir = VNotchWeir::ninety_degree(0.58).unwrap();
        let expected = 0.58 * (8.0 / 15.0) * (2.0 * G_STANDARD).sqrt() * 1.0 * 0.20_f64.powf(2.5);
        let got = weir.discharge(0.20).unwrap();
        assert!(
            (got - expected).abs() < EPS,
            "V-notch Q mismatch: got {got}, expected {expected}"
        );
        // 90-degree notch half-angle tangent is exactly 1.
        assert!((weir.half_angle_tangent() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn vnotch_scales_as_head_to_the_five_halves() {
        // Doubling H multiplies Q by 2^2.5.
        let weir = VNotchWeir::ninety_degree(0.6).unwrap();
        let q1 = weir.discharge(0.1).unwrap();
        let q2 = weir.discharge(0.2).unwrap();
        let ratio = q2 / q1;
        assert!(
            (ratio - 2.0_f64.powf(2.5)).abs() < EPS,
            "head exponent wrong: ratio {ratio}, want {}",
            2.0_f64.powf(2.5)
        );
    }

    #[test]
    fn vnotch_scales_with_half_angle_tangent() {
        // A 60-degree notch has tan(30) = 1/sqrt(3); a 90-degree notch
        // has tan(45) = 1. Their discharges at equal head, Cd and g
        // must be in the ratio of those tangents.
        let sixty = VNotchWeir::new(PI / 3.0, 0.6).unwrap();
        let ninety = VNotchWeir::new(PI / 2.0, 0.6).unwrap();
        let q60 = sixty.discharge(0.15).unwrap();
        let q90 = ninety.discharge(0.15).unwrap();
        let want = (PI / 6.0).tan() / (PI / 4.0).tan(); // tan(30)/tan(45)
        assert!(
            (q60 / q90 - want).abs() < EPS,
            "tan(θ/2) scaling wrong: ratio {}, want {want}",
            q60 / q90
        );
    }

    #[test]
    fn vnotch_higher_head_more_flow() {
        let weir = VNotchWeir::ninety_degree(0.58).unwrap();
        let mut prev = weir.discharge(0.02).unwrap();
        for h in [0.05, 0.10, 0.20, 0.40] {
            let q = weir.discharge(h).unwrap();
            assert!(q > prev, "not monotone at H={h}: q {q}, prev {prev}");
            prev = q;
        }
    }

    // ----------------------------------------------------------------
    // Cross-type behaviour: the V-notch's defining advantage.
    // ----------------------------------------------------------------

    #[test]
    fn vnotch_better_at_low_flow_than_rectangular() {
        // The steeper H^2.5 law makes a unit fractional change in head
        // produce a *bigger* fractional change in Q for the V-notch than
        // for the rectangular weir — the sensitivity that makes a
        // V-notch the better gauge at low flows.
        //
        // d(ln Q)/d(ln H) is exactly the head exponent: 2.5 vs 1.5.
        let rect = RectangularWeir::new(1.0, 0.62).unwrap();
        let vee = VNotchWeir::ninety_degree(0.58).unwrap();

        let h_lo = 0.10;
        let h_hi = 0.11; // +10% head.

        let rect_sens = rect.discharge(h_hi).unwrap() / rect.discharge(h_lo).unwrap();
        let vee_sens = vee.discharge(h_hi).unwrap() / vee.discharge(h_lo).unwrap();

        assert!(
            vee_sens > rect_sens,
            "V-notch should be more head-sensitive: vee {vee_sens}, rect {rect_sens}"
        );
        // And those sensitivities equal the power-law ratios.
        assert!((rect_sens - (h_hi / h_lo).powf(1.5)).abs() < EPS);
        assert!((vee_sens - (h_hi / h_lo).powf(2.5)).abs() < EPS);
    }

    #[test]
    fn vnotch_passes_less_flow_at_small_head() {
        // At a small, equal head a typical 90-degree V-notch passes much
        // less than a metre-wide rectangular weir — the small-flow regime
        // where the notch gives readable resolution.
        let rect = RectangularWeir::new(1.0, 0.62).unwrap();
        let vee = VNotchWeir::ninety_degree(0.58).unwrap();
        let h = 0.05;
        assert!(
            vee.discharge(h).unwrap() < rect.discharge(h).unwrap(),
            "expected V-notch < rectangular at small head"
        );
    }

    // ----------------------------------------------------------------
    // Validated constructors / error taxonomy.
    // ----------------------------------------------------------------

    #[test]
    fn rejects_non_positive_geometry_and_head() {
        assert_eq!(
            RectangularWeir::new(0.0, 0.6).unwrap_err(),
            WeirError::NonPositive {
                name: "crest_length",
                value: 0.0
            }
        );
        assert_eq!(
            RectangularWeir::new(1.0, -0.1).unwrap_err(),
            WeirError::NonPositive {
                name: "discharge_coefficient",
                value: -0.1
            }
        );
        let weir = RectangularWeir::new(1.0, 0.6).unwrap();
        assert_eq!(
            weir.discharge(0.0).unwrap_err(),
            WeirError::NonPositive {
                name: "head",
                value: 0.0
            }
        );
    }

    #[test]
    fn rejects_non_finite_inputs() {
        // NaN != NaN, so match structurally rather than with assert_eq!.
        match RectangularWeir::new(f64::NAN, 0.6).unwrap_err() {
            WeirError::NotFinite { name, value } => {
                assert_eq!(name, "crest_length");
                assert!(value.is_nan(), "expected NaN, got {value}");
            }
            other => panic!("wrong error: {other:?}"),
        }
        let weir = VNotchWeir::ninety_degree(0.58).unwrap();
        match weir.discharge(f64::INFINITY).unwrap_err() {
            WeirError::NotFinite { name, value } => {
                assert_eq!(name, "head");
                assert!(value.is_infinite(), "expected inf, got {value}");
            }
            other => panic!("wrong error: {other:?}"),
        }
    }

    #[test]
    fn rejects_out_of_range_notch_angle() {
        // Zero, negative, and >= pi are all invalid.
        for bad in [0.0, -0.5, PI, PI + 0.1] {
            match VNotchWeir::new(bad, 0.6).unwrap_err() {
                WeirError::NotchAngleOutOfRange { radians } => {
                    assert!((radians - bad).abs() < EPS, "echoed {radians}, sent {bad}");
                }
                other => panic!("wrong error for angle {bad}: {other:?}"),
            }
        }
    }

    #[test]
    fn error_codes_and_categories_are_stable() {
        let e = RectangularWeir::new(-1.0, 0.6).unwrap_err();
        assert_eq!(e.code(), "weir.non-positive");
        assert_eq!(e.category(), ErrorCategory::Input);

        let e = RectangularWeir::new(1.0, -1.0).unwrap_err();
        assert_eq!(e.category(), ErrorCategory::Config);

        let e = VNotchWeir::new(PI, 0.6).unwrap_err();
        assert_eq!(e.code(), "weir.notch-angle-out-of-range");
        assert_eq!(e.category(), ErrorCategory::Input);
    }

    #[test]
    fn accessors_round_trip_constructor_inputs() {
        let weir = RectangularWeir::with_gravity(1.7, 0.63, 9.81).unwrap();
        assert!((weir.crest_length_m() - 1.7).abs() < EPS);
        assert!((weir.discharge_coefficient() - 0.63).abs() < EPS);
        assert!((weir.gravity() - 9.81).abs() < EPS);

        let vee = VNotchWeir::with_gravity(PI / 2.0, 0.59, 9.78).unwrap();
        assert!((vee.vertex_angle_rad() - PI / 2.0).abs() < EPS);
        assert!((vee.discharge_coefficient() - 0.59).abs() < EPS);
        assert!((vee.gravity() - 9.78).abs() < EPS);
    }

    #[test]
    fn serde_round_trips_both_weirs() {
        let rect = RectangularWeir::new(2.0, 0.62).unwrap();
        let json = serde_json::to_string(&rect).unwrap();
        let back: RectangularWeir = serde_json::from_str(&json).unwrap();
        assert_eq!(rect, back);

        let vee = VNotchWeir::ninety_degree(0.58).unwrap();
        let json = serde_json::to_string(&vee).unwrap();
        let back: VNotchWeir = serde_json::from_str(&json).unwrap();
        assert_eq!(vee, back);
    }
}
