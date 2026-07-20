//! Fin heat-dissipation rate for the insulated (adiabatic) tip.
//!
//! Integrating the conduction at the base of a constant-area fin with an
//! insulated tip gives the closed-form heat rate
//!
//! ```text
//! Q = sqrt(h * P * k * A) * theta_b * tanh(m * L)          [watts]
//! ```
//!
//! where `theta_b = T_base - T_inf` is the base temperature excess (K) and
//! `m = sqrt(h P / (k A))`. The same quantity can be written through the fin
//! efficiency as `Q = eta * h * (P * L) * theta_b`, since `h P L` is the
//! convecting surface area and `eta` the efficiency `tanh(mL)/(mL)`; the two
//! forms are algebraically identical and the tests below check they agree.
//!
//! ## Limiting behaviour
//!
//! For large `mL`, `tanh(mL) -> 1`, so the heat rate saturates at
//! `Q_max = sqrt(h P k A) * theta_b` no matter how much longer the fin is —
//! the analytic statement that "a very long fin's extra length is wasted".

use crate::error::{FinError, Result};

/// Fin tip heat rate `Q = sqrt(h P k A) * theta_b * tanh(mL)`, in watts.
///
/// Inputs: convection coefficient `h` (W/m^2 K), perimeter `perimeter`
/// (m), conductivity `k` (W/m K), area (m^2), base temperature excess
/// `theta_b = T_base - T_inf` (K), and the dimensionless length `mlength`
/// (`m * L`). The first four must be finite and strictly positive; `theta_b`
/// must be finite (it may be negative, meaning the base is colder than the
/// fluid, which simply flips the sign of `Q`); `mlength` must be finite and
/// non-negative.
///
/// # Errors
///
/// Propagates the positivity checks on `h`, `perimeter`, `k`, `area`, the
/// finiteness check on `theta_b`, and the non-negativity check on `mlength`.
///
/// # Example
///
/// ```
/// use valenx_fin::heat_rate;
/// // sqrt(hPkA) = sqrt(100 * 0.1 * 200 * 1e-4) = sqrt(0.2) = 0.4472136.
/// // theta_b = 80 K, mL = 1  ->  Q = 0.4472136 * 80 * tanh(1).
/// let q = heat_rate(100.0, 0.1, 200.0, 1.0e-4, 80.0, 1.0).unwrap();
/// let expected = 0.2_f64.sqrt() * 80.0 * 1.0_f64.tanh();
/// assert!((q - expected).abs() < 1e-9);
/// ```
#[allow(clippy::too_many_arguments)]
pub fn heat_rate(
    h: f64,
    perimeter: f64,
    k: f64,
    area: f64,
    theta_b: f64,
    mlength: f64,
) -> Result<f64> {
    let h = FinError::positive("h", h)?;
    let perimeter = FinError::positive("perimeter", perimeter)?;
    let k = FinError::positive("k", k)?;
    let area = FinError::positive("area", area)?;
    let theta_b = FinError::finite("theta_b", theta_b)?;
    let x = FinError::non_negative("mL", mlength)?;
    Ok((h * perimeter * k * area).sqrt() * theta_b * x.tanh())
}

/// Maximum (long-fin) heat rate `Q_max = sqrt(h P k A) * theta_b`, in watts.
///
/// This is the `tanh(mL) -> 1` ceiling of [`heat_rate`]: the most a fin of
/// the given cross-section and material can ever dissipate, approached as the
/// fin becomes long. Useful as a normaliser and as a sanity bound (`heat_rate`
/// is always `<= Q_max` in magnitude).
///
/// # Errors
///
/// Propagates the positivity checks on `h`, `perimeter`, `k`, `area` and the
/// finiteness check on `theta_b`.
///
/// # Example
///
/// ```
/// use valenx_fin::max_heat_rate;
/// let qmax = max_heat_rate(100.0, 0.1, 200.0, 1.0e-4, 80.0).unwrap();
/// assert!((qmax - 0.2_f64.sqrt() * 80.0).abs() < 1e-9);
/// ```
pub fn max_heat_rate(h: f64, perimeter: f64, k: f64, area: f64, theta_b: f64) -> Result<f64> {
    let h = FinError::positive("h", h)?;
    let perimeter = FinError::positive("perimeter", perimeter)?;
    let k = FinError::positive("k", k)?;
    let area = FinError::positive("area", area)?;
    let theta_b = FinError::finite("theta_b", theta_b)?;
    Ok((h * perimeter * k * area).sqrt() * theta_b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::efficiency::efficiency;

    const EPS: f64 = 1e-9;

    #[test]
    fn heat_rate_matches_hand_worked_value() {
        // sqrt(hPkA) = sqrt(0.2); theta_b = 80; tanh(1).
        let q = heat_rate(100.0, 0.1, 200.0, 1.0e-4, 80.0, 1.0).unwrap();
        let expected = 0.2_f64.sqrt() * 80.0 * 1.0_f64.tanh();
        assert!((q - expected).abs() < EPS);
    }

    #[test]
    fn heat_rate_zero_length_is_zero() {
        // tanh(0) = 0  ->  no length, no heat.
        let q = heat_rate(100.0, 0.1, 200.0, 1.0e-4, 80.0, 0.0).unwrap();
        assert!(q.abs() < EPS);
    }

    #[test]
    fn heat_rate_equals_efficiency_times_convective_form() {
        // Q = sqrt(hPkA) theta_b tanh(mL)  must equal  eta * h * P * L * theta_b
        // with eta = tanh(mL)/(mL) and L = mL / m, m = sqrt(hP/kA).
        let (h, p, k, a, theta_b) = (75.0_f64, 0.12, 180.0, 2.0e-4, 65.0);
        let m = (h * p / (k * a)).sqrt();
        let mlength = 1.7_f64;
        let length = mlength / m;
        let q_closed = heat_rate(h, p, k, a, theta_b, mlength).unwrap();
        let eta = efficiency(mlength).unwrap();
        let q_eff = eta * h * p * length * theta_b;
        assert!((q_closed - q_eff).abs() < 1e-7);
    }

    #[test]
    fn heat_rate_saturates_at_max_for_long_fin() {
        // Large mL: tanh -> 1, so Q -> Q_max.
        let q = heat_rate(100.0, 0.1, 200.0, 1.0e-4, 80.0, 40.0).unwrap();
        let qmax = max_heat_rate(100.0, 0.1, 200.0, 1.0e-4, 80.0).unwrap();
        assert!((q - qmax).abs() < 1e-9);
        // And Q never exceeds Q_max in magnitude.
        assert!(q <= qmax + EPS);
    }

    #[test]
    fn heat_rate_sign_follows_theta_b() {
        // A colder-than-fluid base (theta_b < 0) gives negative Q of equal
        // magnitude — heat flows into the fin instead of out.
        let q_pos = heat_rate(100.0, 0.1, 200.0, 1.0e-4, 50.0, 1.0).unwrap();
        let q_neg = heat_rate(100.0, 0.1, 200.0, 1.0e-4, -50.0, 1.0).unwrap();
        assert!((q_pos + q_neg).abs() < EPS);
    }

    #[test]
    fn max_heat_rate_matches_hand_worked_value() {
        let qmax = max_heat_rate(100.0, 0.1, 200.0, 1.0e-4, 80.0).unwrap();
        assert!((qmax - 0.2_f64.sqrt() * 80.0).abs() < EPS);
    }

    #[test]
    fn heat_rate_rejects_bad_inputs() {
        assert!(heat_rate(0.0, 0.1, 200.0, 1e-4, 80.0, 1.0).is_err());
        assert!(heat_rate(100.0, 0.1, 200.0, 1e-4, f64::NAN, 1.0).is_err());
        assert_eq!(
            heat_rate(100.0, 0.1, 200.0, 1e-4, 80.0, -1.0)
                .unwrap_err()
                .code(),
            "fin.negative"
        );
    }

    #[test]
    fn max_heat_rate_rejects_bad_inputs() {
        assert!(max_heat_rate(100.0, -0.1, 200.0, 1e-4, 80.0).is_err());
        assert!(max_heat_rate(100.0, 0.1, 200.0, 1e-4, f64::INFINITY).is_err());
    }
}
