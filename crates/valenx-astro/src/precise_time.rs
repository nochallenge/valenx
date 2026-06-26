//! Precise, multi-timescale epochs for the trajectory and frame code.
//!
//! The rest of the crate speaks **Julian Date** — [`crate::frames::gmst`]
//! takes a JD, [`crate::frames::J2000`] is one, and the on-orbit element
//! sets are referenced to an epoch. Historically those JDs were fed in
//! "raw", with `UTC` silently standing in for `UT1` and with no notion of
//! the **TAI − UTC** leap-second offset or the relativistic **TDB − TT**
//! difference. That is fine for visibility / coverage geometry but loses
//! the sub-second timescale bookkeeping that precise astrodynamics needs.
//!
//! This module closes that gap by wrapping the [`hifitime`] `Epoch` — a
//! femtosecond-resolution, leap-second-aware instant tagged with its
//! [time scale](TimeScale) — and exposing exactly the few conversions the
//! existing code consumes:
//!
//! - build an epoch from a civil **UTC** or **TAI** calendar date
//!   ([`utc`], [`tai`]),
//! - read it back as the **Julian Date** in the timescale a given
//!   calculation actually wants ([`Epoch::jd_utc`], [`Epoch::jd_tai`],
//!   [`Epoch::jd_tt`], [`Epoch::jd_tdb`]) — feed [`Epoch::jd_ut1_approx`]
//!   (≡ UTC) straight into [`crate::frames::gmst`],
//! - get the **leap-second** offset `TAI − UTC` in effect at that instant
//!   ([`Epoch::tai_minus_utc_seconds`]),
//! - get the relativistic **TDB − TT** offset ([`Epoch::tdb_minus_tt_seconds`]),
//!   the periodic term (≲ 1.7 ms) that separates barycentric dynamical
//!   time from terrestrial time, and
//! - advance an epoch by a number of **SI seconds** ([`Epoch::add_seconds`]),
//!   so a propagation can carry a precise wall-clock alongside its state.
//!
//! Time scales are real and distinct: at the 2017-01-01 leap second the
//! offset became **TAI − UTC = 37 s**, and **TT − TAI = 32.184 s** is the
//! fixed definitional constant. Getting these right is what lets a JD that
//! drives an ephemeris differ — correctly — from the JD that drives the
//! Earth-rotation angle.
//!
//! ## Honest scope
//!
//! This is the **time** half of precise astrodynamics. The leap-second
//! table and the TAI/UTC/TT/TDB transforms come from `hifitime` and are
//! authoritative. What it does *not* add is `UT1 − UTC` (length-of-day),
//! which needs an Earth-orientation-parameter feed — so
//! [`Epoch::jd_ut1_approx`] still equals the UTC JD, the same mean-rotation
//! approximation [`crate::frames`] already documents. The TDB model is the
//! standard analytic series (the dominant annual term), not a full
//! JPL-ephemeris evaluation.

use hifitime::{Epoch as HEpoch, TimeScale, Unit};

/// Fixed offset **TT − TAI = 32.184 s** (the definitional constant linking
/// Terrestrial Time to International Atomic Time).
pub const TT_MINUS_TAI_SECONDS: f64 = 32.184;

/// A precise instant in time, tagged with its time scale.
///
/// A thin, `Copy` wrapper over [`hifitime::Epoch`] exposing the
/// Julian-Date and offset conversions the rest of `valenx-astro` consumes.
/// Construct one with [`utc`] / [`tai`] / [`from_jde_utc`], read a JD with
/// the `jd_*` accessors, and step it with [`Epoch::add_seconds`].
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Epoch(HEpoch);

/// Build an epoch from a **UTC** civil date and time. `nanos` is the
/// sub-second part in nanoseconds (`0` for a whole second). Month and day
/// are one-indexed (January = 1).
///
/// ```
/// let e = valenx_astro::precise_time::utc(2017, 1, 1, 0, 0, 0, 0);
/// // At/after the 2017 leap second, TAI is 37 s ahead of UTC.
/// assert!((e.tai_minus_utc_seconds() - 37.0).abs() < 1e-9);
/// ```
pub fn utc(year: i32, month: u8, day: u8, hour: u8, minute: u8, second: u8, nanos: u32) -> Epoch {
    Epoch(HEpoch::from_gregorian_utc(
        year, month, day, hour, minute, second, nanos,
    ))
}

/// Build an epoch from a **TAI** civil date and time (International Atomic
/// Time — a continuous scale with no leap seconds). Arguments as for [`utc`].
pub fn tai(year: i32, month: u8, day: u8, hour: u8, minute: u8, second: u8, nanos: u32) -> Epoch {
    Epoch(HEpoch::from_gregorian_tai(
        year, month, day, hour, minute, second, nanos,
    ))
}

/// Build an epoch from a **UTC Julian Date** (days). The inverse of
/// [`Epoch::jd_utc`]; handy for round-tripping a JD that the frame / orbit
/// code is already working in back into a timescale-aware instant.
pub fn from_jde_utc(jde_utc_days: f64) -> Epoch {
    Epoch(HEpoch::from_jde_utc(jde_utc_days))
}

impl Epoch {
    /// The **J2000.0** epoch (2000-01-01 12:00:00 **TT**) as a precise
    /// instant — the reference the crate's element sets and
    /// [`crate::frames::J2000`] are tied to.
    pub fn j2000() -> Self {
        // hifitime's TT reference; equals JD 2451545.0 in the TT scale.
        Epoch(HEpoch::from_gregorian(
            2000,
            1,
            1,
            12,
            0,
            0,
            0,
            TimeScale::TT,
        ))
    }

    /// Julian Date in the **UTC** scale (days). This is the value the
    /// Earth-rotation / ground-geometry code wants: pass it (or the
    /// identical [`Epoch::jd_ut1_approx`]) to [`crate::frames::gmst`].
    pub fn jd_utc(&self) -> f64 {
        self.0.to_jde_utc_days()
    }

    /// Julian Date in the **TAI** scale (days) — continuous atomic time,
    /// `jd_utc + (TAI − UTC)/86400`.
    pub fn jd_tai(&self) -> f64 {
        self.0.to_jde_tai_days()
    }

    /// Julian Date in the **TT** (Terrestrial Time) scale (days), the usual
    /// independent variable for Earth-satellite dynamics:
    /// `jd_tai + 32.184/86400`.
    pub fn jd_tt(&self) -> f64 {
        // TT = TAI + 32.184 s, exactly.
        self.0.to_jde_tai_days() + TT_MINUS_TAI_SECONDS / 86_400.0
    }

    /// Julian Date in the **TDB** (Barycentric Dynamical Time) scale (days)
    /// — the timescale planetary/lunar ephemerides are tabulated in. Differs
    /// from TT only by a small (≲ 1.7 ms) periodic relativistic term.
    pub fn jd_tdb(&self) -> f64 {
        self.0.to_jde_tdb_days()
    }

    /// Julian Date suitable for the **UT1**-based sidereal-time series.
    ///
    /// `UT1 − UTC` (length-of-day) needs an Earth-orientation feed this crate
    /// does not carry, so this is the **UTC** JD — the same mean-rotation
    /// approximation [`crate::frames`] already uses. Provided as a named
    /// accessor so call sites read honestly ("UT1 ≈ UTC here") rather than
    /// silently reusing the UTC value.
    pub fn jd_ut1_approx(&self) -> f64 {
        self.0.to_jde_utc_days()
    }

    /// The leap-second offset **TAI − UTC** (seconds) in effect at this
    /// instant, from the IERS leap-second table.
    ///
    /// E.g. `37.0` for any instant on/after 2017-01-01, `36.0` for
    /// 2015-07-01 … 2016-12-31, and so on back through the table.
    pub fn tai_minus_utc_seconds(&self) -> f64 {
        // `false` => use the full (IERS + pre-1972 fractional) table.
        self.0.leap_seconds(false).unwrap_or(0.0)
    }

    /// The relativistic offset **TDB − TT** (seconds) at this instant — the
    /// small periodic term (annual, amplitude ≲ 1.7 ms) separating
    /// barycentric dynamical time from terrestrial time.
    pub fn tdb_minus_tt_seconds(&self) -> f64 {
        (self.jd_tdb() - self.jd_tt()) * 86_400.0
    }

    /// A new epoch advanced by `seconds` SI seconds (negative steps back).
    /// Leap seconds are handled by the underlying continuous time line, so
    /// adding `60.0` always advances exactly one minute of physical time.
    pub fn add_seconds(&self, seconds: f64) -> Epoch {
        Epoch(self.0 + seconds * Unit::Second)
    }

    /// Elapsed SI seconds from `self` to `later` (`later − self`); negative
    /// if `later` precedes `self`.
    pub fn seconds_until(&self, later: Epoch) -> f64 {
        (later.0 - self.0).to_seconds()
    }

    /// The civil **UTC** calendar breakdown `(year, month, day, hour,
    /// minute, second, nanosecond)` — the inverse of [`utc`], for display.
    pub fn to_utc_parts(&self) -> (i32, u8, u8, u8, u8, u8, u32) {
        self.0.to_gregorian_utc()
    }

    /// The underlying [`hifitime::Epoch`], for callers that need the full
    /// hifitime API (other time scales, formatting, durations).
    pub fn inner(&self) -> HEpoch {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::frames::{gmst, J2000};

    #[test]
    fn leap_second_offset_is_37_after_2017() {
        // GROUND TRUTH: the 2017-01-01 leap second brought TAI − UTC to 37 s,
        // and it has stood there since.
        let e = utc(2020, 6, 15, 12, 0, 0, 0);
        assert!(
            (e.tai_minus_utc_seconds() - 37.0).abs() < 1e-9,
            "TAI-UTC = {} s, expected 37",
            e.tai_minus_utc_seconds()
        );
        // Just after the leap second itself.
        let e2 = utc(2017, 1, 1, 0, 0, 1, 0);
        assert!((e2.tai_minus_utc_seconds() - 37.0).abs() < 1e-9);
    }

    #[test]
    fn leap_second_offset_steps_down_through_history() {
        // GROUND TRUTH (IERS): 36 s in the 2015-07 … 2016-12 window, and
        // 32 s across 1999-01 … 2005-12 (a long stretch with no inserted
        // leap second).
        let mid_2016 = utc(2016, 6, 1, 0, 0, 0, 0);
        assert!(
            (mid_2016.tai_minus_utc_seconds() - 36.0).abs() < 1e-9,
            "mid-2016 TAI-UTC = {}",
            mid_2016.tai_minus_utc_seconds()
        );
        let year_2000 = utc(2000, 1, 1, 0, 0, 0, 0);
        assert!(
            (year_2000.tai_minus_utc_seconds() - 32.0).abs() < 1e-9,
            "2000 TAI-UTC = {}",
            year_2000.tai_minus_utc_seconds()
        );
    }

    // A Julian Date near the present (~2.46e6 days) has an f64 ULP of
    // ~5e-10 days ≈ 4e-5 s, so any offset reconstructed as a *difference of
    // two JDs* carries that much floating-point granularity. The timescale
    // OFFSET tests below therefore check to ~1e-4 s (the JD float floor),
    // not to the nanosecond — the offsets themselves are exact, this is
    // honest about the JD representation, not slack in the model.
    const JD_DIFF_TOL_SECONDS: f64 = 1.0e-4;

    #[test]
    fn tt_is_tai_plus_32_184() {
        // TT − TAI is the fixed definitional 32.184 s, so the TT and TAI
        // Julian Dates differ by exactly that in days (to the JD float floor).
        let e = utc(2024, 3, 21, 6, 30, 0, 0);
        let diff_seconds = (e.jd_tt() - e.jd_tai()) * 86_400.0;
        assert!(
            (diff_seconds - 32.184).abs() < JD_DIFF_TOL_SECONDS,
            "TT-TAI = {diff_seconds} s, expected 32.184"
        );
    }

    #[test]
    fn tai_jd_leads_utc_jd_by_the_leap_offset() {
        // TAI is ahead of UTC by exactly the leap-second offset, so the JDs
        // differ by (TAI − UTC)/86400 days (to the JD float floor).
        let e = utc(2022, 9, 1, 0, 0, 0, 0);
        let diff_seconds = (e.jd_tai() - e.jd_utc()) * 86_400.0;
        assert!(
            (diff_seconds - e.tai_minus_utc_seconds()).abs() < JD_DIFF_TOL_SECONDS,
            "JD(TAI)-JD(UTC) = {diff_seconds} s vs leap {} s",
            e.tai_minus_utc_seconds()
        );
    }

    #[test]
    fn tdb_minus_tt_is_a_small_periodic_term() {
        // TDB − TT is bounded by ~1.7 ms (it never grows secularly). Sample
        // a few epochs and confirm the magnitude stays in band.
        for month in [1u8, 4, 7, 10] {
            let e = utc(2023, month, 1, 0, 0, 0, 0);
            let d = e.tdb_minus_tt_seconds().abs();
            assert!(
                d < 2.0e-3,
                "|TDB-TT| = {d} s at month {month} exceeds ~1.7 ms"
            );
        }
    }

    #[test]
    fn j2000_jd_is_2451545_in_tt() {
        // GROUND TRUTH: J2000.0 is JD 2451545.0 in TT, and matches the
        // crate's frames::J2000 constant.
        let e = Epoch::j2000();
        assert!(
            (e.jd_tt() - 2_451_545.0).abs() < 1e-6,
            "J2000 JD(TT) = {}",
            e.jd_tt()
        );
        assert!((e.jd_tt() - J2000).abs() < 1e-6);
    }

    #[test]
    fn utc_jd_feeds_gmst_consistently() {
        // The whole point of the bridge: the UTC/UT1-approx JD an Epoch hands
        // back is exactly what frames::gmst expects. Build the J2000 instant
        // as a UTC calendar date and confirm gmst on its UT1-approx JD lands
        // on the known 280.46° (within the UT1≈UTC leap offset, which at
        // J2000 is 32 s → a few arc-minutes; we check the value is sane and
        // matches calling gmst on the same number).
        let e = utc(2000, 1, 1, 12, 0, 0, 0);
        let theta = gmst(e.jd_ut1_approx());
        // gmst on the identical JD must agree bit-for-bit.
        assert_eq!(theta, gmst(e.jd_utc()));
        // And it is a valid wrapped angle near the J2000 GMST.
        assert!((0.0..std::f64::consts::TAU).contains(&theta));
        assert!(
            (theta.to_degrees() - 280.46).abs() < 0.5,
            "GMST {}° not near 280.46",
            theta.to_degrees()
        );
    }

    #[test]
    fn add_seconds_advances_physical_time_and_round_trips() {
        let e0 = utc(2021, 11, 5, 8, 15, 30, 0);
        let e1 = e0.add_seconds(3_600.0); // +1 hour
        assert!((e0.seconds_until(e1) - 3_600.0).abs() < 1e-6);
        // One hour later in JD is 1/24 day.
        assert!((e1.jd_utc() - e0.jd_utc() - 1.0 / 24.0).abs() < 1e-9);
        // Stepping back returns the original instant.
        let back = e1.add_seconds(-3_600.0);
        assert!((back.seconds_until(e0)).abs() < 1e-9);
    }

    #[test]
    fn utc_parts_round_trip() {
        let (y, mo, d, h, mi, s, ns) = (2019, 7, 22, 14, 3, 9, 250_000_000);
        let e = utc(y, mo, d, h, mi, s, ns);
        let (y2, mo2, d2, h2, mi2, s2, ns2) = e.to_utc_parts();
        assert_eq!((y, mo, d, h, mi, s), (y2, mo2, d2, h2, mi2, s2));
        // Sub-second part recovered to nanosecond precision.
        assert!((ns2 as i64 - ns as i64).abs() <= 1, "ns {ns} != {ns2}");
    }
}
