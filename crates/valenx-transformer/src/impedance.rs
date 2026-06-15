//! Reflected (referred) impedance.
//!
//! A load impedance `Zs` connected across the secondary appears, when
//! viewed from the primary terminals, scaled by the square of the turns
//! ratio:
//!
//! ```text
//! Zp = a^2 * Zs,   a = Np / Ns
//! ```
//!
//! This follows directly from the [`crate::ratio`] relations. The
//! primary sees `Zp = Vp / Ip`; substitute `Vp = a Vs` and `Ip = Is / a`
//! to get `Zp = a^2 (Vs / Is) = a^2 Zs`. The same `a^2` law is why a
//! transformer is used for *impedance matching* — choosing the turns
//! ratio sets the reflected load resistance seen by a source.
//!
//! ## Honest scope
//!
//! Only the ideal reflection is modelled. A real transformer also adds
//! its own series leakage reactance and winding resistance, and a shunt
//! magnetising branch, none of which appear here. The impedance value
//! is treated as a single real magnitude `|Z|`; this module does not
//! carry phase, so it is the resistive / magnitude form of the law, not
//! a full complex two-port.

use crate::error::TransformerError;
use crate::ratio::TurnsRatio;

/// Reflect a secondary-side load impedance to the primary side:
/// `Zp = a^2 * Zs`.
///
/// The impedance is a real magnitude (ohms); a step-down ratio
/// (`a > 1`) makes the load look *larger* from the primary, and a
/// step-up ratio (`a < 1`) makes it look *smaller*.
///
/// # Errors
///
/// Returns [`TransformerError::Invalid`] if `load_secondary` is not
/// finite or is negative (a passive load impedance magnitude is
/// non-negative).
pub fn reflect_to_primary(
    ratio: &TurnsRatio,
    load_secondary: f64,
) -> Result<f64, TransformerError> {
    if !load_secondary.is_finite() || load_secondary < 0.0 {
        return Err(TransformerError::invalid(
            "load_secondary",
            format!(
                "secondary load impedance must be finite and non-negative, got {load_secondary}"
            ),
        ));
    }
    let a = ratio.ratio();
    Ok(a * a * load_secondary)
}

/// Refer a primary-side impedance to the secondary side:
/// `Zs = Zp / a^2`.
///
/// The inverse of [`reflect_to_primary`]; reflecting a load to the
/// primary and back to the secondary recovers the original value.
///
/// # Errors
///
/// Returns [`TransformerError::Invalid`] if `impedance_primary` is not
/// finite or is negative. The denominator `a^2` is guaranteed non-zero
/// because the turns ratio is validated strictly positive.
pub fn reflect_to_secondary(
    ratio: &TurnsRatio,
    impedance_primary: f64,
) -> Result<f64, TransformerError> {
    if !impedance_primary.is_finite() || impedance_primary < 0.0 {
        return Err(TransformerError::invalid(
            "impedance_primary",
            format!("primary impedance must be finite and non-negative, got {impedance_primary}"),
        ));
    }
    let a = ratio.ratio();
    Ok(impedance_primary / (a * a))
}

/// Choose the turns ratio that matches a source impedance to a load
/// impedance by reflection: `a = sqrt(Zsource / Zload)`.
///
/// At this ratio the load reflected to the primary equals the source
/// impedance, the maximum-power-transfer condition for a purely
/// resistive match.
///
/// # Errors
///
/// Returns [`TransformerError::Invalid`] if either argument is not
/// finite, if `source` is negative, or if `load` is not strictly
/// positive (it is the ratio denominator). The resulting ratio is then
/// validated by [`TurnsRatio::new`], which rejects a zero source
/// impedance.
pub fn matching_ratio(source: f64, load: f64) -> Result<TurnsRatio, TransformerError> {
    if !source.is_finite() || source < 0.0 {
        return Err(TransformerError::invalid(
            "source",
            format!("source impedance must be finite and non-negative, got {source}"),
        ));
    }
    if !load.is_finite() || load <= 0.0 {
        return Err(TransformerError::invalid(
            "load",
            format!("load impedance must be finite and positive, got {load}"),
        ));
    }
    TurnsRatio::new((source / load).sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute-difference tolerance for the analytic float checks.
    const EPS: f64 = 1e-9;

    #[test]
    fn reflected_impedance_is_a_squared_times_load() {
        // a = 10, Zs = 8 ohm. Zp = 100 * 8 = 800 ohm.
        let t = TurnsRatio::new(10.0).unwrap();
        let zp = reflect_to_primary(&t, 8.0).unwrap();
        assert!((zp - 800.0).abs() < EPS, "Zp got {zp}");
        // Cross-check against an explicit a^2 * Zs computation.
        let a = t.ratio();
        assert!((zp - a * a * 8.0).abs() < EPS, "Zp != a^2 Zs");
    }

    #[test]
    fn reflection_round_trips() {
        // Reflect to primary then back to secondary recovers Zs.
        let t = TurnsRatio::from_turns(120.0, 8.0).unwrap();
        let zs = 5.0;
        let zp = reflect_to_primary(&t, zs).unwrap();
        let zs_back = reflect_to_secondary(&t, zp).unwrap();
        assert!((zs_back - zs).abs() < EPS, "round-trip Zs got {zs_back}");
    }

    #[test]
    fn step_up_lowers_reflected_impedance() {
        // a = 0.5 (step up). Zp = 0.25 * Zload < Zload.
        let t = TurnsRatio::new(0.5).unwrap();
        let zload = 16.0;
        let zp = reflect_to_primary(&t, zload).unwrap();
        assert!(zp < zload, "step-up should lower reflected Z");
        assert!((zp - 4.0).abs() < EPS, "Zp got {zp}");
    }

    #[test]
    fn matching_ratio_reflects_load_to_source() {
        // Match an 8 ohm speaker to a 5000 ohm source. The reflected
        // load at the matching ratio must equal the source impedance.
        let source = 5000.0;
        let load = 8.0;
        let t = matching_ratio(source, load).unwrap();
        let expected_a = (source / load).sqrt();
        assert!((t.ratio() - expected_a).abs() < EPS, "a got {}", t.ratio());
        let reflected = reflect_to_primary(&t, load).unwrap();
        assert!(
            (reflected - source).abs() < 1e-6,
            "reflected load {reflected} should match source {source}"
        );
    }

    #[test]
    fn zero_load_reflects_to_zero() {
        // A short (Zs = 0) reflects to a short on the primary.
        let t = TurnsRatio::new(3.0).unwrap();
        assert!((reflect_to_primary(&t, 0.0).unwrap()).abs() < EPS);
    }

    #[test]
    fn rejects_negative_and_non_finite_impedances() {
        let t = TurnsRatio::new(2.0).unwrap();
        assert!(reflect_to_primary(&t, -1.0).is_err());
        assert!(reflect_to_primary(&t, f64::NAN).is_err());
        assert!(reflect_to_secondary(&t, -1.0).is_err());
        assert!(matching_ratio(50.0, 0.0).is_err());
        assert!(matching_ratio(-1.0, 8.0).is_err());
    }
}
