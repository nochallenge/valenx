//! Basic rating life `L10` and its conversion to operating hours.
//!
//! The ISO 281 **basic rating life** is the life, in millions of
//! revolutions, that 90 % of an apparently identical group of
//! bearings will reach or exceed under a given load before the first
//! evidence of rolling-contact fatigue. It is
//!
//! ```text
//! L10 = (C / P)^p   [millions of revolutions]
//! ```
//!
//! where `C` is the **basic dynamic load rating** (the constant load a
//! bearing can theoretically carry for one million revolutions), `P`
//! is the [dynamic equivalent load](crate::EquivalentLoad), and `p` is
//! the [load-life exponent](crate::BearingType) (`3` ball, `10/3`
//! roller).
//!
//! Because revolutions are awkward to plan maintenance around, the
//! life is usually re-expressed in **hours** at a fixed shaft speed:
//!
//! ```text
//! L10h = (L10 · 1e6) / (60 · n)   [hours]
//! ```
//!
//! where `n` is the rotational speed in revolutions per minute. The
//! `1e6` converts *millions* of revolutions to revolutions and the
//! `60` converts revolutions-per-minute to revolutions-per-hour.

use serde::{Deserialize, Serialize};

use crate::bearing::BearingType;
use crate::error::{require_positive, BearingError};
use crate::load::EquivalentLoad;

/// Revolutions per million revolution, used to turn `L10` (in millions
/// of revolutions) into raw revolutions.
const REVS_PER_MILLION: f64 = 1.0e6;

/// Minutes per hour, used to turn rpm into revolutions per hour.
const MINUTES_PER_HOUR: f64 = 60.0;

/// A complete basic-rating-life result: the `L10` in millions of
/// revolutions plus the inputs that produced it.
///
/// Build it with [`RatingLife::new`]; read the revolution life with
/// [`RatingLife::l10_million_revs`] and convert to operating hours at
/// a shaft speed with [`RatingLife::life_hours`].
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RatingLife {
    /// Basic dynamic load rating `C` (newtons).
    pub dynamic_load_rating: f64,
    /// Dynamic equivalent load `P` (newtons).
    pub equivalent_load: f64,
    /// Bearing element type, which fixes the exponent `p`.
    pub bearing_type: BearingType,
    /// Basic rating life `L10` in millions of revolutions.
    pub l10_million_revs: f64,
}

impl RatingLife {
    /// Evaluate `L10 = (C / P)^p` from the dynamic load rating `C`,
    /// the dynamic equivalent load `P`, and the bearing type (which
    /// supplies the exponent `p`).
    ///
    /// Both `C` and `P` must be finite and strictly positive — a
    /// zero-or-negative rating or load has no physical meaning and
    /// would make the ratio undefined or the life infinite.
    ///
    /// # Errors
    ///
    /// Returns [`BearingError`] when `dynamic_load_rating` or
    /// `equivalent_load` is `NaN`, infinite, or not greater than zero.
    ///
    /// ```
    /// use valenx_bearing::{BearingType, RatingLife};
    /// // C = 50 kN, P = 10 kN, ball (p = 3): L10 = 5^3 = 125 Mrev.
    /// let life = RatingLife::new(50_000.0, 10_000.0, BearingType::Ball).unwrap();
    /// assert!((life.l10_million_revs() - 125.0).abs() < 1e-9);
    /// ```
    pub fn new(
        dynamic_load_rating: f64,
        equivalent_load: f64,
        bearing_type: BearingType,
    ) -> Result<Self, BearingError> {
        let dynamic_load_rating = require_positive("dynamic_load_rating", dynamic_load_rating)?;
        let equivalent_load = require_positive("equivalent_load", equivalent_load)?;
        let p = bearing_type.life_exponent();
        let l10 = (dynamic_load_rating / equivalent_load).powf(p);
        Ok(Self {
            dynamic_load_rating,
            equivalent_load,
            bearing_type,
            l10_million_revs: l10,
        })
    }

    /// Convenience constructor taking a pre-built [`EquivalentLoad`]
    /// instead of a bare `P` value.
    ///
    /// Equivalent to calling [`RatingLife::new`] with
    /// `equivalent_load.value()`.
    ///
    /// # Errors
    ///
    /// Returns [`BearingError`] under the same conditions as
    /// [`RatingLife::new`].
    ///
    /// ```
    /// use valenx_bearing::{BearingType, EquivalentLoad, RatingLife};
    /// let p = EquivalentLoad::new(8000.0, 3000.0, 0.56, 1.6).unwrap(); // 9280 N
    /// let life = RatingLife::from_equivalent_load(60_000.0, &p, BearingType::Ball).unwrap();
    /// let expected = (60_000.0_f64 / 9280.0).powf(3.0);
    /// assert!((life.l10_million_revs() - expected).abs() < 1e-6);
    /// ```
    pub fn from_equivalent_load(
        dynamic_load_rating: f64,
        equivalent_load: &EquivalentLoad,
        bearing_type: BearingType,
    ) -> Result<Self, BearingError> {
        Self::new(dynamic_load_rating, equivalent_load.value(), bearing_type)
    }

    /// The basic rating life `L10` in millions of revolutions.
    #[must_use]
    pub fn l10_million_revs(&self) -> f64 {
        self.l10_million_revs
    }

    /// The basic rating life expressed in operating hours at a shaft
    /// speed of `rpm` revolutions per minute:
    /// `L10h = L10 · 1e6 / (60 · rpm)`.
    ///
    /// For a fixed revolution life, doubling the speed halves the
    /// hours, because the bearing reaches the same number of
    /// revolutions in half the time.
    ///
    /// # Errors
    ///
    /// Returns [`BearingError`] when `rpm` is `NaN`, infinite, or not
    /// greater than zero.
    ///
    /// ```
    /// use valenx_bearing::{BearingType, RatingLife};
    /// let life = RatingLife::new(50_000.0, 10_000.0, BearingType::Ball).unwrap(); // 125 Mrev
    /// // 125e6 revs / (60 * 1500 rpm) ≈ 1388.9 h
    /// let hours = life.life_hours(1500.0).unwrap();
    /// assert!((hours - 125.0 * 1.0e6 / (60.0 * 1500.0)).abs() < 1e-6);
    /// ```
    pub fn life_hours(&self, rpm: f64) -> Result<f64, BearingError> {
        let rpm = require_positive("rpm", rpm)?;
        Ok(self.l10_million_revs * REVS_PER_MILLION / (MINUTES_PER_HOUR * rpm))
    }
}

/// Evaluate the basic rating life `L10 = (C / P)^p` directly, without
/// building a [`RatingLife`] struct.
///
/// `dynamic_load_rating` is `C` and `equivalent_load` is `P`, both in
/// newtons; the exponent `p` comes from `bearing_type`. Returns the
/// life in millions of revolutions.
///
/// # Errors
///
/// Returns [`BearingError`] when `dynamic_load_rating` or
/// `equivalent_load` is `NaN`, infinite, or not greater than zero.
///
/// ```
/// use valenx_bearing::{l10_million_revs, BearingType};
/// // Roller (p = 10/3): C/P = 2 -> L10 = 2^(10/3) ≈ 10.0794.
/// let l10 = l10_million_revs(20_000.0, 10_000.0, BearingType::Roller).unwrap();
/// assert!((l10 - 2.0_f64.powf(10.0 / 3.0)).abs() < 1e-9);
/// ```
pub fn l10_million_revs(
    dynamic_load_rating: f64,
    equivalent_load: f64,
    bearing_type: BearingType,
) -> Result<f64, BearingError> {
    Ok(RatingLife::new(dynamic_load_rating, equivalent_load, bearing_type)?.l10_million_revs())
}

/// Convert a basic rating life in millions of revolutions to operating
/// hours at `rpm` revolutions per minute:
/// `L10h = L10 · 1e6 / (60 · rpm)`.
///
/// # Errors
///
/// Returns [`BearingError`] when `l10_million_revs` or `rpm` is `NaN`,
/// infinite, or not greater than zero.
///
/// ```
/// use valenx_bearing::life_hours_from_revs;
/// // 100 Mrev at 3000 rpm = 100e6 / (60*3000) ≈ 555.56 h.
/// let hours = life_hours_from_revs(100.0, 3000.0).unwrap();
/// assert!((hours - 100.0 * 1.0e6 / (60.0 * 3000.0)).abs() < 1e-9);
/// ```
pub fn life_hours_from_revs(l10_million_revs: f64, rpm: f64) -> Result<f64, BearingError> {
    let l10 = require_positive("l10_million_revs", l10_million_revs)?;
    let rpm = require_positive("rpm", rpm)?;
    Ok(l10 * REVS_PER_MILLION / (MINUTES_PER_HOUR * rpm))
}
