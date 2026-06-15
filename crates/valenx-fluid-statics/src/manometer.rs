//! Manometers — reading pressure differences from liquid-column heights.
//!
//! A manometer balances an unknown pressure against the weight of a
//! column of liquid. The governing idea is the **manometer rule**:
//! walking along a connected body of static fluid, the pressure *rises*
//! by `rho * g * dz` for every metre you descend and *falls* by the same
//! for every metre you climb, and it is continuous across a fluid-fluid
//! interface. Summing those increments from one open end to the other
//! and setting the result equal to the known pressure gives the unknown.
//!
//! **Simple (open) U-tube.** A vessel of gas at gauge pressure `P` is
//! connected to one leg; the other leg is open to atmosphere; the
//! heavier gauge liquid (density `rho_m`) stands a height `R` higher in
//! the open leg. With the gas column's weight negligible,
//!
//! ```text
//! P_gauge = rho_m * g * R.
//! ```
//!
//! When the working fluid above the gauge liquid is itself a liquid of
//! density `rho_f` whose column height is `a`, its weight must be
//! subtracted:
//!
//! ```text
//! P_gauge = rho_m * g * R - rho_f * g * a.
//! ```
//!
//! **Differential U-tube.** To measure the pressure difference between
//! two points `A` and `B` in a working fluid of density `rho_f`,
//! connected through a gauge liquid of density `rho_m` reading `R`, with
//! `A` a height `z_a` and `B` a height `z_b` above an arbitrary datum,
//!
//! ```text
//! P_A - P_B = (rho_m - rho_f) * g * R + rho_f * g * (z_b - z_a).
//! ```
//!
//! When the two taps are at the same elevation (`z_a == z_b`) this
//! collapses to the familiar `(rho_m - rho_f) * g * R`.

use crate::error::{require_finite, require_non_negative, require_positive, Result};
use crate::fluid::{Fluid, STANDARD_GRAVITY};

/// Gauge pressure of a gas vessel read from a simple open U-tube
/// manometer, in pascals: `P_gauge = rho_m * g * R`.
///
/// `gauge_fluid` is the manometer (gauge) liquid and `reading_m` is the
/// height difference `R` of that liquid between the two legs. The weight
/// of the gas column is taken as negligible.
///
/// # Errors
///
/// Returns [`Invalid`](crate::FluidStaticsError::Invalid)
/// if `gravity` is not strictly positive or `reading_m` is negative /
/// non-finite.
pub fn open_manometer_gauge(gauge_fluid: &Fluid, gravity: f64, reading_m: f64) -> Result<f64> {
    let gravity = require_positive("gravity", gravity)?;
    let reading_m = require_non_negative("reading_m", reading_m)?;
    Ok(gauge_fluid.density() * gravity * reading_m)
}

/// Gauge pressure read from an open U-tube whose working fluid above the
/// gauge liquid is a liquid of non-negligible density, in pascals:
/// `P_gauge = rho_m * g * R - rho_f * g * a`.
///
/// `gauge_fluid` / `reading_m` describe the gauge liquid and its height
/// difference `R`; `working_fluid` / `working_column_m` describe the
/// liquid column of height `a` standing above the gauge liquid on the
/// vessel side. The result may be negative (a suction / partial vacuum).
///
/// # Errors
///
/// Returns [`Invalid`](crate::FluidStaticsError::Invalid)
/// if `gravity` is not strictly positive or either height is negative /
/// non-finite.
pub fn open_manometer_gauge_with_working_fluid(
    gauge_fluid: &Fluid,
    reading_m: f64,
    working_fluid: &Fluid,
    working_column_m: f64,
    gravity: f64,
) -> Result<f64> {
    let gravity = require_positive("gravity", gravity)?;
    let reading_m = require_non_negative("reading_m", reading_m)?;
    let working_column_m = require_non_negative("working_column_m", working_column_m)?;
    Ok(gauge_fluid.density() * gravity * reading_m
        - working_fluid.density() * gravity * working_column_m)
}

/// Geometry of a differential U-tube manometer measuring the pressure
/// difference `P_A - P_B` between two points in a working fluid.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct DifferentialManometer {
    /// Density of the working (process) fluid filling the connecting
    /// lines, in kilograms per cubic metre.
    pub working_fluid: Fluid,
    /// Density of the heavier gauge (manometer) liquid, in kilograms per
    /// cubic metre.
    pub gauge_fluid: Fluid,
    /// Manometer reading `R` — the height of the gauge-liquid column
    /// difference between the two legs, in metres. Non-negative.
    pub reading_m: f64,
    /// Elevation of tap `A` above an arbitrary common datum, in metres.
    pub elevation_a_m: f64,
    /// Elevation of tap `B` above the same datum, in metres.
    pub elevation_b_m: f64,
}

impl DifferentialManometer {
    /// The pressure difference `P_A - P_B`, in pascals:
    /// `(rho_m - rho_f) * g * R + rho_f * g * (z_b - z_a)`.
    ///
    /// A positive value means `A` is at the higher pressure. The sign of
    /// the result is meaningful and may be negative.
    ///
    /// # Errors
    ///
    /// Returns [`Invalid`](crate::FluidStaticsError::Invalid)
    /// if `gravity` is not strictly positive, `reading_m` is negative /
    /// non-finite, or either elevation is non-finite.
    pub fn pressure_difference(&self, gravity: f64) -> Result<f64> {
        let gravity = require_positive("gravity", gravity)?;
        let reading_m = require_non_negative("reading_m", self.reading_m)?;
        let z_a = require_finite("elevation_a_m", self.elevation_a_m)?;
        let z_b = require_finite("elevation_b_m", self.elevation_b_m)?;

        let rho_m = self.gauge_fluid.density();
        let rho_f = self.working_fluid.density();
        Ok((rho_m - rho_f) * gravity * reading_m + rho_f * gravity * (z_b - z_a))
    }
}

/// A single vertical segment of a manometer path: a fluid spanning a
/// signed vertical drop. Used by [`balance_path`] to add up the pressure
/// change along a hand-traced manometer circuit.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Segment {
    /// The fluid occupying this segment.
    pub fluid: Fluid,
    /// Signed vertical *descent* across the segment, in metres: positive
    /// when the path goes **down** (pressure rises by `rho * g * drop`),
    /// negative when it goes **up** (pressure falls).
    pub drop_m: f64,
}

impl Segment {
    /// Construct a downward segment of the given fluid and positive
    /// descent `drop_m` (pressure increases along it).
    pub fn down(fluid: Fluid, drop_m: f64) -> Self {
        Segment { fluid, drop_m }
    }

    /// Construct an upward segment of the given fluid and positive climb
    /// `rise_m` (pressure decreases along it).
    pub fn up(fluid: Fluid, rise_m: f64) -> Self {
        Segment {
            fluid,
            drop_m: -rise_m,
        }
    }
}

/// Sum the pressure change along a traced manometer path, in pascals.
///
/// Starting from a point at known pressure, each [`Segment`] contributes
/// `+rho * g * drop` (descending) or `-rho * g * rise` (ascending). The
/// returned total is `P_end - P_start`. For a *closed* manometer loop
/// that returns to its starting height through balanced columns the
/// total is zero — the manometer-balance condition this crate validates.
///
/// # Errors
///
/// Returns [`Invalid`](crate::FluidStaticsError::Invalid)
/// if `gravity` is not strictly positive or any segment drop is
/// non-finite.
pub fn balance_path(segments: &[Segment], gravity: f64) -> Result<f64> {
    let gravity = require_positive("gravity", gravity)?;
    let mut total = 0.0_f64;
    for (i, seg) in segments.iter().enumerate() {
        // Index the parameter name so a bad segment is identifiable.
        let drop = require_finite(
            SEGMENT_DROP_NAMES[i.min(SEGMENT_DROP_NAMES.len() - 1)],
            seg.drop_m,
        )?;
        total += seg.fluid.density() * gravity * drop;
    }
    Ok(total)
}

// A small pool of `&'static str` names so `balance_path` can report
// which segment had a non-finite drop without allocating.
const SEGMENT_DROP_NAMES: [&str; 4] = [
    "segment[0].drop_m",
    "segment[1].drop_m",
    "segment[2].drop_m",
    "segment[..].drop_m",
];

/// Gauge pressure from a simple open U-tube under [`STANDARD_GRAVITY`],
/// in pascals — a convenience wrapper around [`open_manometer_gauge`].
///
/// # Errors
///
/// Returns [`Invalid`](crate::FluidStaticsError::Invalid)
/// if `reading_m` is negative / non-finite.
pub fn open_manometer_gauge_standard(gauge_fluid: &Fluid, reading_m: f64) -> Result<f64> {
    open_manometer_gauge(gauge_fluid, STANDARD_GRAVITY, reading_m)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pressure::gauge_pressure;

    const EPS: f64 = 1e-6;

    #[test]
    fn open_manometer_reads_rho_m_g_r() {
        // 250 mm of mercury -> P = 13534 * 9.80665 * 0.25 = 33180.9... Pa.
        let p = open_manometer_gauge(&Fluid::mercury(), STANDARD_GRAVITY, 0.25).unwrap();
        let expected = 13_534.0 * STANDARD_GRAVITY * 0.25;
        assert!((p - expected).abs() < EPS, "got {p}");
        assert!((p - 33_180.94).abs() < 1.0, "got {p}");
    }

    #[test]
    fn open_manometer_matches_equivalent_head() {
        // The gauge pressure read by the manometer equals the hydrostatic
        // pressure of R metres of the gauge fluid: cross-check the two
        // independent code paths.
        let r = 0.42;
        let via_manometer = open_manometer_gauge(&Fluid::mercury(), STANDARD_GRAVITY, r).unwrap();
        let via_pressure = gauge_pressure(&Fluid::mercury(), STANDARD_GRAVITY, r).unwrap();
        assert!((via_manometer - via_pressure).abs() < EPS);
    }

    #[test]
    fn working_fluid_column_reduces_reading() {
        // A water column above the mercury lowers the inferred gauge
        // pressure: P = rho_m*g*R - rho_f*g*a.
        let gauge = open_manometer_gauge_with_working_fluid(
            &Fluid::mercury(),
            0.30,
            &Fluid::water(),
            0.10,
            STANDARD_GRAVITY,
        )
        .unwrap();
        let expected = 13_534.0 * STANDARD_GRAVITY * 0.30 - 1000.0 * STANDARD_GRAVITY * 0.10;
        assert!((gauge - expected).abs() < EPS, "got {gauge}");
        // It is strictly less than ignoring the water column.
        let bare = open_manometer_gauge(&Fluid::mercury(), STANDARD_GRAVITY, 0.30).unwrap();
        assert!(gauge < bare);
    }

    #[test]
    fn differential_same_elevation_is_density_contrast_times_reading() {
        // Two taps at the same height: P_A - P_B = (rho_m - rho_f)*g*R.
        let m = DifferentialManometer {
            working_fluid: Fluid::water(),
            gauge_fluid: Fluid::mercury(),
            reading_m: 0.15,
            elevation_a_m: 1.0,
            elevation_b_m: 1.0,
        };
        let dp = m.pressure_difference(STANDARD_GRAVITY).unwrap();
        let expected = (13_534.0 - 1000.0) * STANDARD_GRAVITY * 0.15;
        assert!((dp - expected).abs() < EPS, "got {dp}");
        assert!(dp > 0.0);
    }

    #[test]
    fn differential_accounts_for_tap_elevation() {
        // Raising tap B relative to A adds rho_f*g*(z_b - z_a) to P_A-P_B.
        let base = DifferentialManometer {
            working_fluid: Fluid::water(),
            gauge_fluid: Fluid::mercury(),
            reading_m: 0.15,
            elevation_a_m: 0.0,
            elevation_b_m: 0.0,
        };
        let raised = DifferentialManometer {
            elevation_b_m: 0.5,
            ..base
        };
        let dp_base = base.pressure_difference(STANDARD_GRAVITY).unwrap();
        let dp_raised = raised.pressure_difference(STANDARD_GRAVITY).unwrap();
        let delta = dp_raised - dp_base;
        let expected = 1000.0 * STANDARD_GRAVITY * 0.5;
        assert!((delta - expected).abs() < EPS, "delta {delta}");
    }

    #[test]
    fn elevation_contribution_is_antisymmetric() {
        // The gauge-deflection term (rho_m - rho_f)*g*R is fixed by the
        // instrument; the *elevation* term rho_f*g*(z_b - z_a) flips sign
        // when the two tap heights are exchanged. Subtracting out the
        // common reading term must therefore give equal-and-opposite
        // elevation contributions.
        let ab = DifferentialManometer {
            working_fluid: Fluid::water(),
            gauge_fluid: Fluid::mercury(),
            reading_m: 0.2,
            elevation_a_m: 0.3,
            elevation_b_m: 0.8,
        };
        let ba = DifferentialManometer {
            elevation_a_m: ab.elevation_b_m,
            elevation_b_m: ab.elevation_a_m,
            ..ab
        };
        // The reading-only term shared by both orientations.
        let reading_term = (ab.gauge_fluid.density() - ab.working_fluid.density())
            * STANDARD_GRAVITY
            * ab.reading_m;
        let elev_ab = ab.pressure_difference(STANDARD_GRAVITY).unwrap() - reading_term;
        let elev_ba = ba.pressure_difference(STANDARD_GRAVITY).unwrap() - reading_term;
        assert!(
            (elev_ab + elev_ba).abs() < EPS,
            "elev_ab {elev_ab} elev_ba {elev_ba}"
        );
        // And each equals rho_f*g*(z_b - z_a) with the expected sign.
        let expected_ab =
            ab.working_fluid.density() * STANDARD_GRAVITY * (ab.elevation_b_m - ab.elevation_a_m);
        assert!((elev_ab - expected_ab).abs() < EPS, "elev_ab {elev_ab}");
    }

    #[test]
    fn closed_balanced_loop_sums_to_zero() {
        // Down 0.5 m and back up 0.5 m through the SAME fluid returns to
        // the start pressure: the manometer-balance condition.
        let segs = [
            Segment::down(Fluid::water(), 0.5),
            Segment::up(Fluid::water(), 0.5),
        ];
        let total = balance_path(&segs, STANDARD_GRAVITY).unwrap();
        assert!(total.abs() < EPS, "got {total}");
    }

    #[test]
    fn balance_path_recovers_simple_manometer() {
        // Trace the open U-tube as a path: from the gas surface, descend
        // through the (massless) gas to the mercury, then climb R through
        // the mercury to the open atmosphere. P_atm - P_gas = -rho_m*g*R,
        // so the gas gauge pressure is +rho_m*g*R.
        let r = 0.18;
        let segs = [Segment::up(Fluid::mercury(), r)];
        let p_end_minus_start = balance_path(&segs, STANDARD_GRAVITY).unwrap();
        let gauge = -p_end_minus_start; // P_gas(gauge) = P_gas - P_atm
        let expected = open_manometer_gauge(&Fluid::mercury(), STANDARD_GRAVITY, r).unwrap();
        assert!((gauge - expected).abs() < EPS, "gauge {gauge}");
    }

    #[test]
    fn up_and_down_segments_have_opposite_sign() {
        let down = balance_path(&[Segment::down(Fluid::water(), 1.0)], STANDARD_GRAVITY).unwrap();
        let up = balance_path(&[Segment::up(Fluid::water(), 1.0)], STANDARD_GRAVITY).unwrap();
        assert!(down > 0.0 && up < 0.0, "down {down} up {up}");
        assert!((down + up).abs() < EPS);
    }

    #[test]
    fn rejects_bad_arguments() {
        assert!(open_manometer_gauge(&Fluid::mercury(), 0.0, 0.1).is_err());
        assert!(open_manometer_gauge(&Fluid::mercury(), STANDARD_GRAVITY, -0.1).is_err());
        let bad = DifferentialManometer {
            working_fluid: Fluid::water(),
            gauge_fluid: Fluid::mercury(),
            reading_m: -0.1,
            elevation_a_m: 0.0,
            elevation_b_m: 0.0,
        };
        assert!(bad.pressure_difference(STANDARD_GRAVITY).is_err());
        assert!(
            balance_path(&[Segment::down(Fluid::water(), f64::NAN)], STANDARD_GRAVITY).is_err()
        );
    }

    #[test]
    fn standard_gravity_wrapper_matches_explicit() {
        let a = open_manometer_gauge_standard(&Fluid::mercury(), 0.2).unwrap();
        let b = open_manometer_gauge(&Fluid::mercury(), STANDARD_GRAVITY, 0.2).unwrap();
        assert!((a - b).abs() < EPS);
    }
}
