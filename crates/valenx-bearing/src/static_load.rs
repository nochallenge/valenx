//! Static load safety factor `s0 = C0 / P0` (ISO 76).
//!
//! The dynamic rating life ([`crate::life`]) governs a *rotating* bearing's
//! fatigue. A slowly-rotating, oscillating, or stationary bearing instead
//! fails by **permanent (brinelling) deformation** of the raceway, which
//! ISO 76 guards against with the **static safety factor**
//!
//! ```text
//! s0 = C0 / P0
//! ```
//!
//! where `C0` is the bearing's **basic static load rating** (the load that
//! produces a total permanent deformation of `0.0001 ×` the rolling-element
//! diameter, read from the data sheet) and `P0` is the **static equivalent
//! load** — the combined radial/axial load reduced to a single number,
//!
//! ```text
//! P0 = X0 · Fr + Y0 · Fa,   but never less than Fr.
//! ```
//!
//! The `X0` / `Y0` static load factors are *inputs* from the manufacturer's
//! table (just like the dynamic `X` / `Y`), and ISO 76 requires the result
//! to be taken as `Fr` whenever the formula falls below it. A larger `s0`
//! means a larger margin against indentation; typical required values run
//! from `~0.5` (smooth, low-demand) to `~2` and above (shock loads, high
//! running accuracy).
//!
//! Same honest, calorically-textbook scope as the rest of the crate: the
//! factors are supplied, not guessed.

use serde::{Deserialize, Serialize};

use crate::error::{require_non_negative, require_positive, BearingError};

/// The static equivalent load `P0` and the inputs it was built from.
///
/// Construct it with [`StaticEquivalentLoad::new`] (validated) and read the
/// combined value back with [`StaticEquivalentLoad::value`], which applies
/// the ISO 76 rule `P0 = max(X0·Fr + Y0·Fa, Fr)`.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct StaticEquivalentLoad {
    /// Radial load component `Fr` (newtons), `>= 0`.
    pub radial: f64,
    /// Axial (thrust) load component `Fa` (newtons), `>= 0`.
    pub axial: f64,
    /// Dimensionless static radial load factor `X0`, `>= 0`.
    pub x0_factor: f64,
    /// Dimensionless static axial load factor `Y0`, `>= 0`.
    pub y0_factor: f64,
}

impl StaticEquivalentLoad {
    /// Build a validated static-equivalent-load case from its radial and
    /// axial components and the static load factors `X0`, `Y0`.
    ///
    /// All four arguments must be finite and non-negative.
    ///
    /// # Errors
    ///
    /// Returns [`BearingError`] when any argument is `NaN`, infinite, or
    /// negative.
    ///
    /// ```
    /// use valenx_bearing::StaticEquivalentLoad;
    /// // Fr = 2000 N, Fa = 10000 N, X0 = 0.6, Y0 = 0.5
    /// let p0 = StaticEquivalentLoad::new(2000.0, 10000.0, 0.6, 0.5).unwrap();
    /// // 0.6*2000 + 0.5*10000 = 6200 N (> Fr, so used directly)
    /// assert!((p0.value() - 6200.0).abs() < 1e-9);
    /// ```
    pub fn new(
        radial: f64,
        axial: f64,
        x0_factor: f64,
        y0_factor: f64,
    ) -> Result<Self, BearingError> {
        let radial = require_non_negative("radial", radial)?;
        let axial = require_non_negative("axial", axial)?;
        let x0_factor = require_non_negative("x0_factor", x0_factor)?;
        let y0_factor = require_non_negative("y0_factor", y0_factor)?;
        Ok(Self {
            radial,
            axial,
            x0_factor,
            y0_factor,
        })
    }

    /// A purely radial static load: `P0 = Fr` (`X0 = 1`, `Y0 = 0`,
    /// `Fa = 0`).
    ///
    /// # Errors
    ///
    /// Returns [`BearingError`] when `radial` is `NaN`, infinite, or
    /// negative.
    pub fn radial_only(radial: f64) -> Result<Self, BearingError> {
        Self::new(radial, 0.0, 1.0, 0.0)
    }

    /// The static equivalent load `P0 = max(X0·Fr + Y0·Fa, Fr)` in newtons.
    ///
    /// The ISO 76 floor at `Fr` is what makes a purely radial case collapse
    /// to `P0 = Fr` even when `X0 < 1`.
    #[must_use]
    pub fn value(&self) -> f64 {
        (self.x0_factor * self.radial + self.y0_factor * self.axial).max(self.radial)
    }

    /// The static safety factor `s0 = C0 / P0` against this static
    /// equivalent load, given the bearing's basic static load rating `C0`.
    ///
    /// # Errors
    ///
    /// Returns [`BearingError`] when `basic_static_load_rating` is `NaN`,
    /// infinite, or not greater than zero, or when the static equivalent
    /// load `P0` is zero (an unloaded bearing has no defined safety factor).
    ///
    /// ```
    /// use valenx_bearing::StaticEquivalentLoad;
    /// let p0 = StaticEquivalentLoad::new(2000.0, 10000.0, 0.6, 0.5).unwrap(); // 6200 N
    /// // C0 = 31 kN -> s0 = 31000 / 6200 = 5.
    /// assert!((p0.safety_factor(31_000.0).unwrap() - 5.0).abs() < 1e-9);
    /// ```
    pub fn safety_factor(&self, basic_static_load_rating: f64) -> Result<f64, BearingError> {
        let c0 = require_positive("basic_static_load_rating", basic_static_load_rating)?;
        let p0 = require_positive("static_equivalent_load", self.value())?;
        Ok(c0 / p0)
    }
}

/// Evaluate the static safety factor `s0 = C0 / P0` directly from the basic
/// static load rating `C0` and a static equivalent load `P0` (both newtons),
/// without going through [`StaticEquivalentLoad`].
///
/// # Errors
///
/// Returns [`BearingError`] when either argument is `NaN`, infinite, or not
/// greater than zero.
///
/// ```
/// use valenx_bearing::static_safety_factor;
/// // C0 = 20 kN, P0 = 8 kN -> s0 = 2.5.
/// assert!((static_safety_factor(20_000.0, 8000.0).unwrap() - 2.5).abs() < 1e-9);
/// ```
pub fn static_safety_factor(
    basic_static_load_rating: f64,
    static_equivalent_load: f64,
) -> Result<f64, BearingError> {
    let c0 = require_positive("basic_static_load_rating", basic_static_load_rating)?;
    let p0 = require_positive("static_equivalent_load", static_equivalent_load)?;
    Ok(c0 / p0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) {
        assert!((a - b).abs() < 1e-9, "expected {a} ~= {b}");
    }

    #[test]
    fn static_equivalent_uses_the_formula_when_above_fr() {
        // Fr = 2000, Fa = 10000, X0 = 0.6, Y0 = 0.5 -> 1200 + 5000 = 6200 > Fr.
        let p0 = StaticEquivalentLoad::new(2000.0, 10000.0, 0.6, 0.5).unwrap();
        close(p0.value(), 6200.0);
    }

    #[test]
    fn static_equivalent_floors_at_fr() {
        // Fr = 8000, Fa = 1000, X0 = 0.6, Y0 = 0.5 -> 4800 + 500 = 5300 < Fr,
        // so the ISO 76 rule clamps P0 up to Fr = 8000.
        let p0 = StaticEquivalentLoad::new(8000.0, 1000.0, 0.6, 0.5).unwrap();
        close(p0.value(), 8000.0);
        // A purely radial load reduces to P0 = Fr.
        close(
            StaticEquivalentLoad::radial_only(5000.0).unwrap().value(),
            5000.0,
        );
    }

    #[test]
    fn safety_factor_is_c0_over_p0() {
        let p0 = StaticEquivalentLoad::new(2000.0, 10000.0, 0.6, 0.5).unwrap(); // 6200
        close(p0.safety_factor(31_000.0).unwrap(), 5.0);
        // The free function agrees with the method.
        close(static_safety_factor(31_000.0, p0.value()).unwrap(), 5.0);
        // A larger rating gives a proportionally larger margin.
        assert!(p0.safety_factor(62_000.0).unwrap() > p0.safety_factor(31_000.0).unwrap());
    }

    #[test]
    fn rejects_bad_inputs() {
        let p0 = StaticEquivalentLoad::new(2000.0, 10000.0, 0.6, 0.5).unwrap();
        // Non-positive rating.
        assert!(p0.safety_factor(0.0).is_err());
        assert!(p0.safety_factor(-1.0).is_err());
        assert!(static_safety_factor(f64::NAN, 6200.0).is_err());
        // Unloaded bearing -> P0 = 0 -> no defined safety factor.
        let unloaded = StaticEquivalentLoad::new(0.0, 0.0, 0.6, 0.5).unwrap();
        close(unloaded.value(), 0.0);
        assert!(unloaded.safety_factor(31_000.0).is_err());
        // Negative load components are rejected at construction.
        assert!(StaticEquivalentLoad::new(-1.0, 0.0, 0.6, 0.5).is_err());
    }
}
