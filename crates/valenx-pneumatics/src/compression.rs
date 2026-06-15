//! Absolute pressure bookkeeping and the compression ratio.
//!
//! ## Model
//!
//! Pneumatic gauges read *gauge* pressure (relative to the surrounding
//! atmosphere); the gas laws and free-air conversions need *absolute*
//! pressure. The two differ by the local atmospheric pressure:
//!
//! ```text
//! p_abs = p_gauge + p_atm
//! ```
//!
//! The **compression ratio** is how many times the air has been
//! compressed relative to free (atmospheric) air:
//!
//! ```text
//! r = p_abs / p_atm = (p_gauge + p_atm) / p_atm = 1 + p_gauge / p_atm
//! ```
//!
//! A volume `V` of compressed air at gauge pressure `p_gauge` expands to
//! `r * V` of free air when vented to atmosphere — this ratio is exactly
//! what converts a cylinder's swept volume into free-air consumption (see
//! [`crate::consumption`]).
//!
//! ## Honest scope
//!
//! This is the ideal-gas / isothermal volume ratio at a single reference
//! atmosphere. It uses whatever `p_atm` you supply (sea-level standard is
//! ~101 325 Pa; altitude, weather and temperature all shift it). It does
//! not model temperature changes during compression or expansion, real-gas
//! compressibility, or humidity.

use crate::error::{PneumaticsError, Result};

/// Standard sea-level atmospheric pressure, 101 325 Pa (1 atm).
///
/// A convenient default for `p_atm` when a more precise local value is not
/// available.
pub const STANDARD_ATMOSPHERE_PA: f64 = 101_325.0;

/// Convert a gauge pressure to an absolute pressure, `p_abs = p_gauge +
/// p_atm`, in pascals.
///
/// # Errors
///
/// - [`PneumaticsError::Negative`] if `gauge_pressure` is negative or
///   non-finite (gauge may be zero).
/// - [`PneumaticsError::NonPositive`] if `atmospheric_pressure` is not
///   finite and `> 0`.
///
/// # Examples
///
/// ```
/// use valenx_pneumatics::compression::{absolute_pressure, STANDARD_ATMOSPHERE_PA};
/// // 6 bar gauge at standard atmosphere -> 7.01325 bar absolute.
/// let p = absolute_pressure(600_000.0, STANDARD_ATMOSPHERE_PA).unwrap();
/// assert!((p - 701_325.0).abs() < 1e-9);
/// ```
pub fn absolute_pressure(gauge_pressure: f64, atmospheric_pressure: f64) -> Result<f64> {
    let pg = PneumaticsError::non_negative("gauge_pressure", gauge_pressure)?;
    let pa = PneumaticsError::positive("atmospheric_pressure", atmospheric_pressure)?;
    Ok(pg + pa)
}

/// The compression ratio `r = p_abs / p_atm` from a *gauge* supply
/// pressure, i.e. `1 + p_gauge / p_atm` (dimensionless).
///
/// This is the multiplier that turns a volume of compressed air into the
/// volume of free air it expands to at atmosphere.
///
/// # Errors
///
/// - [`PneumaticsError::Negative`] if `gauge_pressure` is negative or
///   non-finite.
/// - [`PneumaticsError::NonPositive`] if `atmospheric_pressure` is not
///   finite and `> 0`.
///
/// # Examples
///
/// ```
/// use valenx_pneumatics::compression::{compression_ratio, STANDARD_ATMOSPHERE_PA};
/// // 6 bar gauge -> r = 1 + 600000/101325 = 6.9214...
/// let r = compression_ratio(600_000.0, STANDARD_ATMOSPHERE_PA).unwrap();
/// assert!((r - 6.921_415_).abs() < 1e-3);
/// // At zero gauge the air is uncompressed: r = 1.
/// assert!((compression_ratio(0.0, STANDARD_ATMOSPHERE_PA).unwrap() - 1.0).abs() < 1e-12);
/// ```
pub fn compression_ratio(gauge_pressure: f64, atmospheric_pressure: f64) -> Result<f64> {
    let p_abs = absolute_pressure(gauge_pressure, atmospheric_pressure)?;
    // `atmospheric_pressure` is already validated positive inside
    // `absolute_pressure`, so this division is finite and non-zero.
    Ok(p_abs / atmospheric_pressure)
}

/// The compression ratio `r = p_abs / p_atm` directly from an *absolute*
/// supply pressure (dimensionless). Both arguments are absolute pressures
/// in pascals.
///
/// # Errors
///
/// - [`PneumaticsError::NonPositive`] if either pressure is not finite and
///   `> 0`.
/// - [`PneumaticsError::Geometry`] if `absolute_pressure < atmospheric`
///   (a compression ratio below 1 is not a compression).
///
/// # Examples
///
/// ```
/// use valenx_pneumatics::compression::compression_ratio_from_absolute;
/// // 700 kPa absolute over 100 kPa atmosphere -> r = 7.
/// let r = compression_ratio_from_absolute(700_000.0, 100_000.0).unwrap();
/// assert!((r - 7.0).abs() < 1e-12);
/// ```
pub fn compression_ratio_from_absolute(
    absolute_pressure: f64,
    atmospheric_pressure: f64,
) -> Result<f64> {
    let p_abs = PneumaticsError::positive("absolute_pressure", absolute_pressure)?;
    let p_atm = PneumaticsError::positive("atmospheric_pressure", atmospheric_pressure)?;
    if p_abs < p_atm {
        return Err(PneumaticsError::Geometry(
            "absolute pressure is below atmospheric (compression ratio < 1)",
        ));
    }
    Ok(p_abs / p_atm)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute-difference tolerance for floating comparisons.
    const EPS: f64 = 1e-9;

    #[test]
    fn absolute_is_gauge_plus_atmosphere() {
        let p = absolute_pressure(600_000.0, STANDARD_ATMOSPHERE_PA).unwrap();
        assert!((p - 701_325.0).abs() < EPS);
    }

    #[test]
    fn zero_gauge_is_one_atmosphere_absolute() {
        let p = absolute_pressure(0.0, STANDARD_ATMOSPHERE_PA).unwrap();
        assert!((p - STANDARD_ATMOSPHERE_PA).abs() < EPS);
    }

    #[test]
    fn compression_ratio_is_one_at_zero_gauge() {
        let r = compression_ratio(0.0, STANDARD_ATMOSPHERE_PA).unwrap();
        assert!((r - 1.0).abs() < EPS);
    }

    #[test]
    fn compression_ratio_round_number_with_round_atmosphere() {
        // Choose p_atm = 100 kPa so the arithmetic is exact: 6 bar gauge
        // -> p_abs = 700 kPa -> r = 7.
        let r = compression_ratio(600_000.0, 100_000.0).unwrap();
        assert!((r - 7.0).abs() < EPS);
    }

    #[test]
    fn gauge_and_absolute_paths_agree() {
        // r from gauge must equal r from the matching absolute pressure.
        let p_atm = STANDARD_ATMOSPHERE_PA;
        let pg = 450_000.0;
        let from_gauge = compression_ratio(pg, p_atm).unwrap();
        let from_abs = compression_ratio_from_absolute(pg + p_atm, p_atm).unwrap();
        assert!((from_gauge - from_abs).abs() < EPS);
    }

    #[test]
    fn compression_ratio_identity_one_plus_pg_over_patm() {
        // r == 1 + p_gauge / p_atm, computed two independent ways.
        let p_atm = STANDARD_ATMOSPHERE_PA;
        let pg = 800_000.0;
        let r = compression_ratio(pg, p_atm).unwrap();
        assert!((r - (1.0 + pg / p_atm)).abs() < EPS);
    }

    #[test]
    fn rejects_negative_gauge() {
        assert!(absolute_pressure(-1.0, STANDARD_ATMOSPHERE_PA).is_err());
        assert!(compression_ratio(-1.0, STANDARD_ATMOSPHERE_PA).is_err());
    }

    #[test]
    fn rejects_nonpositive_atmosphere() {
        assert!(absolute_pressure(100_000.0, 0.0).is_err());
        assert!(compression_ratio(100_000.0, -1.0).is_err());
    }

    #[test]
    fn absolute_path_rejects_subatmospheric() {
        // 90 kPa absolute under a 100 kPa atmosphere is not a compression.
        assert!(compression_ratio_from_absolute(90_000.0, 100_000.0).is_err());
    }
}
