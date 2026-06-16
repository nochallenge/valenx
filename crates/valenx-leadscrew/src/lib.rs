//! # valenx-leadscrew
//!
//! Closed-form kinematics and statics of a stepper-driven **lead screw**
//! (or ball screw) linear stage.
//!
//! ## What
//!
//! Describe a screw by its [`LeadScrew::lead_mm`] (axial advance per
//! revolution) and [`LeadScrew::pitch_diameter_mm`], then read back the
//! four quantities that size a motion axis:
//!
//! - [`LeadScrew::linear_speed_mm_per_min`] /
//!   [`LeadScrew::linear_speed_mm_per_s`] — how fast the nut moves for a
//!   given screw RPM.
//! - [`LeadScrew::thrust_n`] — the axial force a given screw torque
//!   produces (and its inverse, [`LeadScrew::torque_for_thrust_n_mm`]).
//! - [`LeadScrew::resolution_mm`] — the smallest commandable linear step
//!   for a given microsteps-per-revolution count.
//! - [`LeadScrew::back_drive`] — whether the screw self-locks under a
//!   nut/thread friction coefficient, plus the
//!   [`LeadScrew::critical_friction`] boundary and the ideal raising
//!   [`LeadScrew::screw_efficiency`] `eta = tan(lambda)/tan(lambda + phi)`
//!   that friction implies (below `0.5` whenever the screw self-locks).
//!
//! ```
//! use valenx_leadscrew::LeadScrew;
//!
//! // T8x8: 8 mm lead, 8 mm pitch diameter.
//! let screw = LeadScrew::new(8.0, 8.0).unwrap();
//!
//! let feed = screw.linear_speed_mm_per_min(600.0).unwrap(); // 4800 mm/min
//! let thrust = screw.thrust_n(50.0, 0.4).unwrap();          // N from 50 N·mm @ eta 0.4
//! let step = screw.resolution_mm(3200).unwrap();            // 0.0025 mm/microstep
//! let lock = screw.back_drive(0.20).unwrap();               // self-locking?
//!
//! assert!((feed - 4800.0).abs() < 1e-9);
//! assert!(thrust > 0.0);
//! assert!((step - 0.0025).abs() < 1e-12);
//! assert!(lock.back_drivable()); // 8 mm lead on 8 mm screw is steep
//! ```
//!
//! ## Model
//!
//! The screw is an **ideal constant-lead helix** with no thread-form
//! detail and a single lumped efficiency:
//!
//! - **Kinematics.** One screw revolution advances the nut by exactly
//!   one lead, so linear speed is `v = lead * rpm` and resolution is
//!   `lead / microsteps`. These are exact, geometry-only relations.
//! - **Statics.** One revolution does `2 * pi * T` of input work and
//!   moves the load through `lead`, so the ideal axial force is
//!   `2 * pi * T / lead`. A dimensionless efficiency `eta in (0, 1]`
//!   scales it: `F = 2 * pi * eta * T / lead`. Larger lead trades
//!   force for speed; lower efficiency reduces thrust.
//! - **Back-drive.** With lead angle `lambda = atan(lead / (pi * d_m))`
//!   and friction angle `phi = atan(mu)`, the screw **self-locks** when
//!   `mu >= tan(lambda)` (equivalently `lambda <= phi`) and is
//!   **back-drivable** otherwise. The same two angles give the ideal
//!   raising efficiency `eta = tan(lambda) / tan(lambda + phi)` — the
//!   friction-derived value of the lumped `eta` above, necessarily below
//!   `0.5` whenever the screw self-locks.
//!
//! ## Honest scope
//!
//! Research/educational grade. These are **textbook closed-form
//! models** — the constant-lead screw mechanics and an ideal stepper —
//! and nothing more. In particular this crate does **not** model:
//!
//! - thread-form geometry (ACME / trapezoidal flank-angle corrections to
//!   the friction term), preload, lash/backlash, or wind-up;
//! - dynamic effects — inertia, acceleration limits, resonance, the
//!   stepper torque-speed curve, or missed steps;
//! - column buckling, critical screw speed, wear, lubrication regime, or
//!   thermal growth;
//! - the distinction between sliding (lead screw) and rolling (ball
//!   screw) contact beyond the single lumped `eta`.
//!
//! It is **NOT a clinical/medical or production engineering tool**. Do
//! not size load-bearing or safety-critical hardware from these numbers;
//! treat them as first-order estimates for learning and exploration.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod backdrive;
pub mod error;
pub mod screw;

pub use backdrive::BackDrive;
pub use error::LeadScrewError;
pub use screw::LeadScrew;

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// Tolerance for floating-point ground-truth comparisons.
    const EPS: f64 = 1e-9;

    // ---- speed = lead * rpm -------------------------------------------------

    #[test]
    fn speed_is_lead_times_rpm() {
        let screw = LeadScrew::new(2.0, 8.0).unwrap();
        // 2 mm/rev * 300 rev/min = 600 mm/min exactly.
        let v = screw.linear_speed_mm_per_min(300.0).unwrap();
        assert!((v - 600.0).abs() < EPS, "got {v}");
    }

    #[test]
    fn speed_scales_linearly_with_rpm() {
        let screw = LeadScrew::new(5.0, 10.0).unwrap();
        let v1 = screw.linear_speed_mm_per_min(100.0).unwrap();
        let v2 = screw.linear_speed_mm_per_min(250.0).unwrap();
        // 2.5x the rpm -> 2.5x the feed.
        assert!((v2 / v1 - 2.5).abs() < EPS, "v1={v1} v2={v2}");
    }

    #[test]
    fn speed_mm_per_s_is_per_min_over_60() {
        let screw = LeadScrew::new(4.0, 8.0).unwrap();
        let per_min = screw.linear_speed_mm_per_min(120.0).unwrap();
        let per_s = screw.linear_speed_mm_per_s(120.0).unwrap();
        assert!((per_s - per_min / 60.0).abs() < EPS, "per_s={per_s}");
        // 4 mm/rev * 120 rev/min = 480 mm/min = 8 mm/s.
        assert!((per_s - 8.0).abs() < EPS, "per_s={per_s}");
    }

    // ---- thrust = 2*pi*eta*T / lead -----------------------------------------

    #[test]
    fn thrust_matches_closed_form_ideal() {
        let screw = LeadScrew::new(2.0, 8.0).unwrap();
        // eta = 1, T = 1 N·mm -> F = 2*pi*1*1 / 2 = pi N.
        let f = screw.thrust_n(1.0, 1.0).unwrap();
        assert!((f - PI).abs() < EPS, "got {f}");
    }

    #[test]
    fn thrust_matches_closed_form_general() {
        let lead = 5.0;
        let torque = 40.0;
        let eta = 0.35;
        let screw = LeadScrew::new(lead, 12.0).unwrap();
        let expected = 2.0 * PI * eta * torque / lead;
        let f = screw.thrust_n(torque, eta).unwrap();
        assert!((f - expected).abs() < EPS, "got {f}, want {expected}");
    }

    #[test]
    fn efficiency_reduces_thrust() {
        let screw = LeadScrew::new(4.0, 10.0).unwrap();
        let ideal = screw.thrust_n(10.0, 1.0).unwrap();
        let real = screw.thrust_n(10.0, 0.4).unwrap();
        // Lower efficiency -> strictly less thrust for the same torque.
        assert!(real < ideal, "real={real} ideal={ideal}");
        // And exactly the efficiency ratio.
        assert!((real / ideal - 0.4).abs() < EPS, "ratio={}", real / ideal);
    }

    #[test]
    fn higher_lead_is_faster_but_less_force() {
        let torque = 20.0;
        let eta = 0.9;
        let rpm = 200.0;
        let coarse = LeadScrew::new(8.0, 10.0).unwrap();
        let fine = LeadScrew::new(2.0, 10.0).unwrap();

        // Faster: larger lead -> larger linear speed at the same rpm.
        let v_coarse = coarse.linear_speed_mm_per_min(rpm).unwrap();
        let v_fine = fine.linear_speed_mm_per_min(rpm).unwrap();
        assert!(v_coarse > v_fine, "v_coarse={v_coarse} v_fine={v_fine}");

        // Less force: larger lead -> smaller thrust at the same torque.
        let f_coarse = coarse.thrust_n(torque, eta).unwrap();
        let f_fine = fine.thrust_n(torque, eta).unwrap();
        assert!(f_coarse < f_fine, "f_coarse={f_coarse} f_fine={f_fine}");
    }

    #[test]
    fn torque_for_thrust_is_inverse_of_thrust() {
        let screw = LeadScrew::new(3.0, 9.0).unwrap();
        let eta = 0.5;
        let torque = 25.0;
        let f = screw.thrust_n(torque, eta).unwrap();
        let back = screw.torque_for_thrust_n_mm(f, eta).unwrap();
        // Round-trip recovers the original torque.
        assert!((back - torque).abs() < 1e-7, "back={back} torque={torque}");
    }

    // ---- resolution = lead / microsteps -------------------------------------

    #[test]
    fn resolution_is_lead_over_microsteps() {
        let screw = LeadScrew::new(8.0, 8.0).unwrap();
        // 200-step motor * 16 microstepping = 3200 microsteps/rev.
        let r = screw.resolution_mm(3200).unwrap();
        assert!((r - 8.0 / 3200.0).abs() < 1e-12, "got {r}");
        assert!((r - 0.0025).abs() < 1e-12, "got {r}");
    }

    #[test]
    fn finer_microstepping_gives_finer_resolution() {
        let screw = LeadScrew::new(2.0, 8.0).unwrap();
        let coarse = screw.resolution_mm(200).unwrap();
        let fine = screw.resolution_mm(3200).unwrap();
        // More steps per rev -> smaller increment.
        assert!(fine < coarse, "fine={fine} coarse={coarse}");
        // Exactly 16x finer (3200 / 200).
        assert!(
            (coarse / fine - 16.0).abs() < EPS,
            "ratio={}",
            coarse / fine
        );
    }

    // ---- lead angle ---------------------------------------------------------

    #[test]
    fn lead_angle_matches_atan_formula() {
        let lead = 2.0;
        let dm = 8.0;
        let screw = LeadScrew::new(lead, dm).unwrap();
        let expected = (lead / (PI * dm)).atan();
        assert!(
            (screw.lead_angle_rad() - expected).abs() < 1e-12,
            "got {}",
            screw.lead_angle_rad()
        );
        // Degrees wrapper agrees.
        assert!(
            (screw.lead_angle_deg() - expected.to_degrees()).abs() < 1e-9,
            "deg={}",
            screw.lead_angle_deg()
        );
    }

    // ---- back-drive / self-locking ------------------------------------------

    #[test]
    fn shallow_screw_self_locks() {
        // 1 mm lead on a 10 mm screw -> tan(lambda) = 1/(pi*10) ~= 0.0318.
        let screw = LeadScrew::new(1.0, 10.0).unwrap();
        // Typical steel-on-steel mu ~ 0.15 >> 0.0318 -> self-locking.
        let bd = screw.back_drive(0.15).unwrap();
        assert!(bd.self_locking, "expected self-locking");
        assert!(!bd.back_drivable());
        assert!(bd.locking_margin_rad() > 0.0, "margin should be positive");
    }

    #[test]
    fn steep_screw_back_drives() {
        // 12 mm lead on a 10 mm screw -> tan(lambda) = 12/(pi*10) ~= 0.382.
        let screw = LeadScrew::new(12.0, 10.0).unwrap();
        // mu = 0.15 < 0.382 -> back-drivable.
        let bd = screw.back_drive(0.15).unwrap();
        assert!(bd.back_drivable(), "expected back-drivable");
        assert!(!bd.self_locking);
        assert!(bd.locking_margin_rad() < 0.0, "margin should be negative");
    }

    #[test]
    fn critical_friction_is_tan_lead_angle() {
        let screw = LeadScrew::new(2.0, 8.0).unwrap();
        let mu_crit = screw.critical_friction();
        // tan(lambda) = lead / (pi * d_m).
        assert!(
            (mu_crit - 2.0 / (PI * 8.0)).abs() < 1e-12,
            "mu_crit={mu_crit}"
        );
        // At exactly mu_crit the conservative convention is self-locking.
        let at = screw.back_drive(mu_crit).unwrap();
        assert!(at.self_locking, "boundary should classify as self-locking");
        // Just below it back-drives; just above it locks.
        assert!(screw.back_drive(mu_crit * 0.999).unwrap().back_drivable());
        assert!(screw.back_drive(mu_crit * 1.001).unwrap().self_locking);
    }

    #[test]
    fn frictionless_screw_always_back_drives() {
        // mu = 0 is the ideal frictionless limit: tan(lambda) > 0 always,
        // so 0 >= tan(lambda) is false -> back-drivable for any real lead.
        for &(lead, dm) in &[(0.5, 20.0), (2.0, 8.0), (10.0, 10.0)] {
            let screw = LeadScrew::new(lead, dm).unwrap();
            assert!(
                screw.back_drive(0.0).unwrap().back_drivable(),
                "lead={lead} dm={dm} should back-drive at mu=0"
            );
        }
    }

    // ---- validation / error paths -------------------------------------------

    #[test]
    fn rejects_non_positive_geometry() {
        assert!(LeadScrew::new(0.0, 8.0).is_err());
        assert!(LeadScrew::new(-1.0, 8.0).is_err());
        assert!(LeadScrew::new(2.0, 0.0).is_err());
        assert!(LeadScrew::new(2.0, -3.0).is_err());
    }

    #[test]
    fn rejects_non_finite_geometry() {
        // NaN never equals NaN, so match on the variant + name and check
        // the carried value with `is_nan()` rather than `assert_eq!`.
        match LeadScrew::new(f64::NAN, 8.0).unwrap_err() {
            LeadScrewError::NotFinite { name, value } => {
                assert_eq!(name, "lead_mm");
                assert!(value.is_nan(), "value should be NaN, got {value}");
            }
            other => panic!("expected NotFinite, got {other:?}"),
        }
        assert!(matches!(
            LeadScrew::new(f64::INFINITY, 8.0).unwrap_err(),
            LeadScrewError::NotFinite { .. }
        ));
        // A non-finite pitch diameter is rejected too.
        assert!(matches!(
            LeadScrew::new(2.0, f64::NEG_INFINITY).unwrap_err(),
            LeadScrewError::NotFinite { .. }
        ));
    }

    #[test]
    fn rejects_bad_efficiency() {
        let screw = LeadScrew::new(2.0, 8.0).unwrap();
        assert!(matches!(
            screw.thrust_n(10.0, 0.0).unwrap_err(),
            LeadScrewError::EfficiencyOutOfRange { .. }
        ));
        assert!(matches!(
            screw.thrust_n(10.0, 1.5).unwrap_err(),
            LeadScrewError::EfficiencyOutOfRange { .. }
        ));
        assert!(matches!(
            screw.thrust_n(10.0, -0.1).unwrap_err(),
            LeadScrewError::EfficiencyOutOfRange { .. }
        ));
        // eta == 1.0 is the inclusive upper bound and is allowed.
        assert!(screw.thrust_n(10.0, 1.0).is_ok());
    }

    #[test]
    fn rejects_negative_friction() {
        let screw = LeadScrew::new(2.0, 8.0).unwrap();
        assert_eq!(
            screw.back_drive(-0.1).unwrap_err(),
            LeadScrewError::NegativeFriction {
                name: "mu",
                value: -0.1
            }
        );
    }

    #[test]
    fn rejects_zero_microsteps() {
        let screw = LeadScrew::new(2.0, 8.0).unwrap();
        assert_eq!(
            screw.resolution_mm(0).unwrap_err(),
            LeadScrewError::ZeroMicrosteps(0)
        );
    }

    #[test]
    fn rejects_non_positive_thrust_inputs() {
        let screw = LeadScrew::new(2.0, 8.0).unwrap();
        assert!(screw.thrust_n(0.0, 0.5).is_err());
        assert!(screw.thrust_n(-5.0, 0.5).is_err());
        assert!(screw.linear_speed_mm_per_min(0.0).is_err());
        assert!(screw.linear_speed_mm_per_min(-100.0).is_err());
        assert!(screw.torque_for_thrust_n_mm(0.0, 0.5).is_err());
    }

    // ---- error metadata + serde ---------------------------------------------

    #[test]
    fn error_codes_are_stable() {
        let screw = LeadScrew::new(2.0, 8.0).unwrap();
        assert_eq!(
            LeadScrew::new(0.0, 8.0).unwrap_err().code(),
            "leadscrew.not-positive"
        );
        assert_eq!(
            screw.thrust_n(1.0, 2.0).unwrap_err().code(),
            "leadscrew.efficiency-out-of-range"
        );
        assert_eq!(
            screw.back_drive(-1.0).unwrap_err().code(),
            "leadscrew.negative-friction"
        );
        assert_eq!(
            screw.resolution_mm(0).unwrap_err().code(),
            "leadscrew.zero-microsteps"
        );
    }

    #[test]
    fn leadscrew_roundtrips_through_serde_json() {
        let screw = LeadScrew::new(8.0, 8.0).unwrap();
        let json = serde_json::to_string(&screw).unwrap();
        let back: LeadScrew = serde_json::from_str(&json).unwrap();
        assert_eq!(screw, back);
    }
}
