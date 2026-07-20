//! Receptor occupancy and the potency log-scale (`pEC50`).
//!
//! For a non-cooperative ligand (Hill slope `n = 1`) the Hill equation
//! reduces to the **Clark / Langmuir** occupancy isotherm: the fraction
//! of receptors bound at free concentration `C` is
//!
//! ```text
//! occupancy(C) = C / (K + C)
//! ```
//!
//! where `K` is the dissociation constant (numerically the `EC50` in the
//! occupancy-equals-effect idealisation). This is exactly the `n = 1`
//! case of [`HillCurve::fraction`](crate::hill::HillCurve::fraction);
//! it is provided as a standalone function because the `n = 1` Langmuir
//! form is the workhorse of binding pharmacology.
//!
//! Potency is conventionally reported on a negative-decadic-log scale,
//! the **`pEC50`**:
//!
//! ```text
//! pEC50 = -log10(EC50)
//! ```
//!
//! with `EC50` expressed in molar (mol/L). A higher `pEC50` means a more
//! potent drug. [`p_ec50`] and [`ec50_from_p`] convert between the two
//! representations and are exact inverses.

use crate::error::DoseResponseError;
use crate::Result;

/// Fractional receptor occupancy for a non-cooperative ligand,
/// `occupancy = C / (K + C)`.
///
/// This is the `n = 1` Langmuir isotherm: `0` at `C = 0`, exactly `1/2`
/// at `C = K`, and approaching `1` as `C -> infinity`.
///
/// # Errors
///
/// Returns [`DoseResponseError::NonPositive`] if `k <= 0` (a
/// dissociation constant is a concentration), [`DoseResponseError::Negative`]
/// if `concentration < 0`, or [`DoseResponseError::NonFinite`] for a
/// `NaN` / infinite argument.
pub fn occupancy(k: f64, concentration: f64) -> Result<f64> {
    let k = DoseResponseError::positive("k", k)?;
    let c = DoseResponseError::non_negative("concentration", concentration)?;
    if c == 0.0 {
        return Ok(0.0);
    }
    // r / (1 + r) with r = C / K is the n = 1 Hill fraction; identical to
    // C / (K + C) but avoids cancellation when C and K differ in scale.
    let r = c / k;
    if r.is_infinite() {
        return Ok(1.0);
    }
    Ok(r / (1.0 + r))
}

/// Convert an `EC50` (in molar) to its `pEC50 = -log10(EC50)`.
///
/// # Errors
///
/// Returns [`DoseResponseError::NonPositive`] if `ec50_molar <= 0` (the
/// logarithm is only defined for a positive concentration), or
/// [`DoseResponseError::NonFinite`] for a `NaN` / infinite argument.
pub fn p_ec50(ec50_molar: f64) -> Result<f64> {
    let ec50 = DoseResponseError::positive("ec50_molar", ec50_molar)?;
    Ok(-ec50.log10())
}

/// Convert a `pEC50` back to an `EC50` in molar: `EC50 = 10^(-pEC50)`.
///
/// Exact inverse of [`p_ec50`].
///
/// # Errors
///
/// Returns [`DoseResponseError::NonFinite`] if `p` is `NaN` / infinite.
pub fn ec50_from_p(p: f64) -> Result<f64> {
    let p = DoseResponseError::finite("p_ec50", p)?;
    Ok(10f64.powf(-p))
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    fn approx(a: f64, b: f64) {
        assert!((a - b).abs() < EPS, "expected {a} ≈ {b}");
    }

    #[test]
    fn occupancy_half_at_k() {
        // GROUND TRUTH: occupancy = 1/2 exactly at C = K.
        approx(occupancy(10.0, 10.0).unwrap(), 0.5);
    }

    #[test]
    fn occupancy_limits() {
        // C -> 0 gives 0; C >> K gives -> 1.
        approx(occupancy(2.0, 0.0).unwrap(), 0.0);
        let high = occupancy(1.0, 1.0e9).unwrap();
        assert!(high > 1.0 - 1e-8 && high <= 1.0, "got {high}");
    }

    #[test]
    fn occupancy_textbook_value() {
        // K = 5, C = 15: occ = 15 / (5 + 15) = 15/20 = 0.75.
        approx(occupancy(5.0, 15.0).unwrap(), 0.75);
        // K = 5, C = 5/3: occ = (5/3)/(5 + 5/3) = (5/3)/(20/3) = 0.25.
        approx(occupancy(5.0, 5.0 / 3.0).unwrap(), 0.25);
    }

    #[test]
    fn occupancy_matches_hill_n1() {
        // Must equal the n = 1 Hill fraction with EC50 = K, Emax = 1.
        let curve = crate::hill::HillCurve::new(1.0, 8.0, 1.0).unwrap();
        for &c in &[0.0, 1.0, 8.0, 50.0, 500.0] {
            approx(occupancy(8.0, c).unwrap(), curve.fraction(c).unwrap());
        }
    }

    #[test]
    fn occupancy_rejects_bad_input() {
        assert!(occupancy(0.0, 1.0).is_err()); // k = 0
        assert!(occupancy(-1.0, 1.0).is_err()); // k < 0
        assert!(occupancy(1.0, -1.0).is_err()); // C < 0
        assert!(occupancy(f64::NAN, 1.0).is_err());
    }

    #[test]
    fn p_ec50_known_values() {
        // EC50 = 1e-9 M (1 nM) -> pEC50 = 9.
        approx(p_ec50(1.0e-9).unwrap(), 9.0);
        // EC50 = 1e-6 M (1 µM) -> pEC50 = 6.
        approx(p_ec50(1.0e-6).unwrap(), 6.0);
        // EC50 = 1 M -> pEC50 = 0.
        approx(p_ec50(1.0).unwrap(), 0.0);
    }

    #[test]
    fn p_ec50_round_trip() {
        for &ec50 in &[1.0e-12, 3.3e-9, 1.0e-6, 2.5e-3, 1.0] {
            let p = p_ec50(ec50).unwrap();
            approx(ec50_from_p(p).unwrap(), ec50);
        }
    }

    #[test]
    fn ec50_from_p_known_values() {
        // pEC50 = 9 -> 1e-9 M.
        approx(ec50_from_p(9.0).unwrap(), 1.0e-9);
        approx(ec50_from_p(0.0).unwrap(), 1.0);
    }

    #[test]
    fn potency_log_rejects_bad_input() {
        assert!(p_ec50(0.0).is_err());
        assert!(p_ec50(-1.0e-9).is_err());
        assert!(ec50_from_p(f64::INFINITY).is_err());
    }
}
