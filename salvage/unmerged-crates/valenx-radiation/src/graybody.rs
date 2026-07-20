//! Gray-body radiative exchange between two surfaces.
//!
//! For a small gray surface of area `A` and emissivity `eps` at
//! temperature `T1` exchanging with large isothermal surroundings at
//! `T2` (the classic "small object in a large enclosure" limit), the
//! net rate of radiant heat *leaving* the surface is
//!
//! ```text
//! Q = eps * sigma * A * (T1^4 - T2^4)          [W]
//! ```
//!
//! When `T1 > T2` the surface loses heat (`Q > 0`); when `T1 < T2` it
//! gains heat (`Q < 0`). The sibling [`view_factor`](crate::view_factor)
//! module handles the general two-surface enclosure where neither
//! surface is large.
//!
//! Ground-truth checks: a black body (`eps = 1`) with `T2 = 0`
//! recovers `Q = sigma A T1^4 = A * Eb(T1)`; isothermal surfaces
//! (`T1 = T2`) exchange zero net power.

use crate::error::{check_area, check_temperature, check_unit_interval, Result};
use crate::stefan_boltzmann::SIGMA;

/// Net radiant power leaving a gray surface to large surroundings,
/// `Q = eps sigma A (T1^4 - T2^4)`.
///
/// # Arguments
///
/// `emissivity` — total hemispherical emissivity in `[0, 1]`;
/// `area_m2` — radiating area in `m^2` (`> 0`);
/// `surface_temp_k` — surface temperature `T1` in kelvin (`>= 0`);
/// `surroundings_temp_k` — surroundings temperature `T2` in kelvin
/// (`>= 0`).
///
/// # Returns
///
/// Net power in watts. Positive means net loss from the surface.
///
/// # Errors
///
/// Propagates the domain checks on every argument.
///
/// # Example
///
/// ```
/// use valenx_radiation::graybody::net_exchange_to_surroundings;
/// // Black 1 m^2 plate at 1000 K facing 0 K space: Q = sigma * 1e12.
/// let q = net_exchange_to_surroundings(1.0, 1.0, 1000.0, 0.0).unwrap();
/// assert!((q - 56_703.74).abs() < 1.0);
/// ```
pub fn net_exchange_to_surroundings(
    emissivity: f64,
    area_m2: f64,
    surface_temp_k: f64,
    surroundings_temp_k: f64,
) -> Result<f64> {
    let eps = check_unit_interval("emissivity", emissivity)?;
    let a = check_area("area_m2", area_m2)?;
    let t1 = check_temperature("surface_temp_k", surface_temp_k)?;
    let t2 = check_temperature("surroundings_temp_k", surroundings_temp_k)?;
    Ok(eps * SIGMA * a * (t1.powi(4) - t2.powi(4)))
}

/// Net radiant *flux* (per unit area) leaving a gray surface,
/// `q'' = eps sigma (T1^4 - T2^4)` in `W / m^2`.
///
/// The area-independent form of [`net_exchange_to_surroundings`].
///
/// # Errors
///
/// Propagates the domain checks on emissivity and both temperatures.
pub fn net_flux_to_surroundings(
    emissivity: f64,
    surface_temp_k: f64,
    surroundings_temp_k: f64,
) -> Result<f64> {
    let eps = check_unit_interval("emissivity", emissivity)?;
    let t1 = check_temperature("surface_temp_k", surface_temp_k)?;
    let t2 = check_temperature("surroundings_temp_k", surroundings_temp_k)?;
    Ok(eps * SIGMA * (t1.powi(4) - t2.powi(4)))
}

/// Linearised radiation heat-transfer coefficient
/// `h_r = eps sigma (T1 + T2) (T1^2 + T2^2)` in `W m^-2 K^-1`.
///
/// This is the coefficient such that the net flux equals
/// `q'' = h_r (T1 - T2)`; it follows from factoring the difference of
/// fourth powers, `T1^4 - T2^4 = (T1 - T2)(T1 + T2)(T1^2 + T2^2)`. It
/// is widely used to fold radiation into a combined convection /
/// radiation resistance network.
///
/// # Errors
///
/// Propagates the domain checks on emissivity and both temperatures.
pub fn radiation_coefficient(
    emissivity: f64,
    surface_temp_k: f64,
    surroundings_temp_k: f64,
) -> Result<f64> {
    let eps = check_unit_interval("emissivity", emissivity)?;
    let t1 = check_temperature("surface_temp_k", surface_temp_k)?;
    let t2 = check_temperature("surroundings_temp_k", surroundings_temp_k)?;
    Ok(eps * SIGMA * (t1 + t2) * (t1 * t1 + t2 * t2))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stefan_boltzmann::black_body_emissive_power;

    const EPS: f64 = 1e-9;

    #[test]
    fn isothermal_surfaces_exchange_nothing() {
        let q = net_exchange_to_surroundings(0.7, 2.5, 400.0, 400.0).unwrap();
        assert!(q.abs() < EPS, "got {q}");
        let qf = net_flux_to_surroundings(0.7, 400.0, 400.0).unwrap();
        assert!(qf.abs() < EPS, "got {qf}");
    }

    #[test]
    fn black_body_to_cold_space_recovers_stefan_boltzmann() {
        // eps = 1, T2 = 0  =>  Q = sigma * A * T1^4 = A * Eb(T1).
        let q = net_exchange_to_surroundings(1.0, 3.0, 800.0, 0.0).unwrap();
        let eb = black_body_emissive_power(800.0).unwrap();
        assert!((q - 3.0 * eb).abs() < 1e-6, "got {q}");
    }

    #[test]
    fn hand_worked_two_temperature_loss() {
        // eps=0.8, A=2, T1=500, T2=300.
        // 500^4 - 300^4 = 6.25e10 - 8.1e9 = 5.44e10.
        // Q = 0.8 * 5.670374419e-8 * 2 * 5.44e10 = 4935.494... W.
        let q = net_exchange_to_surroundings(0.8, 2.0, 500.0, 300.0).unwrap();
        let truth = 0.8 * 5.670_374_419e-8 * 2.0 * 5.44e10;
        assert!((q - truth).abs() < 1e-6, "got {q}");
        assert!((q - 4935.494).abs() < 1e-2, "got {q}");
    }

    #[test]
    fn sign_flips_when_surroundings_hotter() {
        let loss = net_exchange_to_surroundings(0.9, 1.0, 600.0, 300.0).unwrap();
        let gain = net_exchange_to_surroundings(0.9, 1.0, 300.0, 600.0).unwrap();
        assert!(loss > 0.0, "loss {loss}");
        assert!(gain < 0.0, "gain {gain}");
        // Antisymmetry: Q(T1,T2) = -Q(T2,T1).
        assert!((loss + gain).abs() < 1e-9, "sum {}", loss + gain);
    }

    #[test]
    fn flux_times_area_equals_power() {
        let q = net_exchange_to_surroundings(0.6, 4.0, 700.0, 350.0).unwrap();
        let qf = net_flux_to_surroundings(0.6, 700.0, 350.0).unwrap();
        assert!((q - 4.0 * qf).abs() < 1e-9, "q={q} qf={qf}");
    }

    #[test]
    fn linearised_coefficient_reproduces_net_flux() {
        // q'' = h_r * (T1 - T2) must equal the exact eps*sigma*(T1^4-T2^4).
        let (eps, t1, t2) = (0.85, 450.0, 300.0);
        let hr = radiation_coefficient(eps, t1, t2).unwrap();
        let q_lin = hr * (t1 - t2);
        let q_exact = net_flux_to_surroundings(eps, t1, t2).unwrap();
        assert!(
            (q_lin - q_exact).abs() < 1e-9,
            "lin={q_lin} exact={q_exact}"
        );
    }

    #[test]
    fn coefficient_equal_temps_limit() {
        // At T1 = T2 = T: h_r = eps*sigma*(2T)*(2T^2) = 4 eps sigma T^3.
        let t = 350.0;
        let hr = radiation_coefficient(0.5, t, t).unwrap();
        let truth = 4.0 * 0.5 * SIGMA * t.powi(3);
        assert!((hr - truth).abs() < 1e-12, "got {hr}");
    }

    #[test]
    fn rejects_bad_inputs() {
        assert!(net_exchange_to_surroundings(1.2, 1.0, 300.0, 200.0).is_err());
        assert!(net_exchange_to_surroundings(0.5, 0.0, 300.0, 200.0).is_err());
        assert!(net_exchange_to_surroundings(0.5, 1.0, -1.0, 200.0).is_err());
        assert!(net_exchange_to_surroundings(0.5, 1.0, 300.0, f64::NAN).is_err());
        assert!(radiation_coefficient(0.5, -1.0, 200.0).is_err());
    }
}
