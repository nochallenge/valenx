//! Bearing element type and its life exponent.
//!
//! The basic-rating-life relation `L10 = (C / P)^p` uses a
//! load-life exponent `p` that depends only on the geometry of the
//! rolling-contact: it is `3` for point-contact (ball) bearings and
//! `10/3 ≈ 3.333…` for line-contact (roller) bearings. These values
//! are the ISO 281 standard; see Harris, *Rolling Bearing Analysis*,
//! and the SKF General Catalogue.

use serde::{Deserialize, Serialize};

/// The kind of rolling element, which fixes the load-life exponent
/// `p` in `L10 = (C / P)^p`.
///
/// The point-contact of a ball gives `p = 3`; the line-contact of a
/// roller gives `p = 10/3`. The larger roller exponent is why roller
/// bearings lose life *less* steeply than ball bearings as the load
/// rises (and conversely gain life faster as the load drops).
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub enum BearingType {
    /// Point-contact (ball) bearing: load-life exponent `p = 3`.
    Ball,
    /// Line-contact (roller) bearing: load-life exponent `p = 10/3`.
    Roller,
}

impl BearingType {
    /// The ISO 281 load-life exponent `p` for this bearing type.
    ///
    /// Returns `3.0` for [`BearingType::Ball`] and `10.0 / 3.0` for
    /// [`BearingType::Roller`].
    ///
    /// ```
    /// use valenx_bearing::BearingType;
    /// assert_eq!(BearingType::Ball.life_exponent(), 3.0);
    /// assert!((BearingType::Roller.life_exponent() - 10.0 / 3.0).abs() < 1e-12);
    /// ```
    #[must_use]
    pub fn life_exponent(self) -> f64 {
        match self {
            BearingType::Ball => 3.0,
            BearingType::Roller => 10.0 / 3.0,
        }
    }
}
