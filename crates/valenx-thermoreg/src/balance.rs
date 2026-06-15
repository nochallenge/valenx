//! The whole-body heat balance and the core-temperature dynamics.
//!
//! ## The balance equation
//!
//! The conceptual human heat balance, in the standard physiology
//! sign convention, is
//!
//! ```text
//!   M - W = R + C + E + S
//! ```
//!
//! where every term is a heat *rate* in watts:
//!
//! - `M` metabolic heat production (always ≥ 0),
//! - `W` external mechanical work done *by* the body (≥ 0; positive
//!   work leaves the energy budget),
//! - `R` net heat *lost* by radiation,
//! - `C` net heat *lost* by convection,
//! - `E` net heat *lost* by evaporation,
//! - `S` heat *stored* in the body (the balancing term).
//!
//! Rearranging for the storage,
//!
//! ```text
//!   S = (M - W) - (R + C + E)
//! ```
//!
//! `S > 0` means the body is gaining heat (core temperature rising);
//! `S < 0` means it is losing heat; `S = 0` is thermal balance.
//!
//! ## Core-temperature dynamics
//!
//! Treating the body as a single lumped mass (see [`crate::body`]),
//! the stored heat raises the temperature through
//!
//! ```text
//!   dT = S * dt / (m * c)
//! ```
//!
//! the lumped-capacitance / "calorimeter" relation `Q = m c ΔT`. This
//! is the same first-order ODE that governs a stirred tank or a
//! one-node building model — appropriate for back-of-envelope work,
//! not for resolving the core-to-skin temperature gradient.

use crate::body::{Body, Sweat};
use crate::environment::Environment;
use crate::error::{finite, non_negative, ThermoregError};
use serde::{Deserialize, Serialize};

/// The metabolic side of the balance: heat produced minus external
/// work performed.
///
/// `metabolic_w` is the total metabolic rate `M` (e.g. ~100 W at rest,
/// hundreds during exercise); `work_w` is the external mechanical
/// power `W` delivered to the surroundings (e.g. pedalling a bike).
/// The *net* metabolic heat that must be dissipated or stored is
/// `M - W`.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Metabolism {
    /// Metabolic heat production `M` (W). Must be ≥ 0.
    pub metabolic_w: f64,
    /// External mechanical work rate `W` (W). Must be ≥ 0.
    pub work_w: f64,
}

impl Metabolism {
    /// Construct a validated [`Metabolism`].
    ///
    /// # Errors
    ///
    /// [`ThermoregError::NotFinite`] for non-finite inputs;
    /// [`ThermoregError::OutOfRange`] if either rate is negative.
    pub fn new(metabolic_w: f64, work_w: f64) -> Result<Self, ThermoregError> {
        Ok(Self {
            metabolic_w: non_negative("metabolic_w", metabolic_w)?,
            work_w: non_negative("work_w", work_w)?,
        })
    }

    /// At-rest metabolism with no external work, `M ≈ 100 W`
    /// (a typical basal-plus-resting value for an adult).
    pub fn resting() -> Self {
        Self {
            metabolic_w: 100.0,
            work_w: 0.0,
        }
    }

    /// The net metabolic heat `M - W` (W) entering the balance.
    pub fn net_heat_production(&self) -> f64 {
        self.metabolic_w - self.work_w
    }
}

/// A fully-resolved heat balance: every term of `M - W = R + C + E + S`
/// in watts.
///
/// Produced by [`heat_balance`]. The sign convention matches the
/// module documentation: `R`, `C`, `E` are heat *losses* (positive =
/// leaving the body) and `storage_w` is the residual `S`.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct HeatBalance {
    /// Net metabolic heat production, `M - W` (W).
    pub metabolic_net_w: f64,
    /// Radiative heat loss `R` (W).
    pub radiation_w: f64,
    /// Convective heat loss `C` (W).
    pub convection_w: f64,
    /// Evaporative heat loss `E` (W).
    pub evaporation_w: f64,
    /// Heat storage `S` (W); the residual that closes the balance.
    pub storage_w: f64,
}

impl HeatBalance {
    /// Total sensible-plus-latent heat *loss*, `R + C + E` (W).
    pub fn total_loss_w(&self) -> f64 {
        self.radiation_w + self.convection_w + self.evaporation_w
    }

    /// The net heat *retained* by the body per second (W); identical
    /// to [`HeatBalance::storage_w`], named for the dynamics.
    ///
    /// Positive heats the core, negative cools it.
    pub fn net_heat_w(&self) -> f64 {
        self.storage_w
    }

    /// `true` when the body is in thermal balance to within `tol`
    /// watts, i.e. the storage term is negligibly small.
    pub fn is_balanced(&self, tol_w: f64) -> bool {
        self.storage_w.abs() <= tol_w
    }

    /// Residual of the balance identity `M - W - (R + C + E) - S`,
    /// which is identically zero by construction.
    ///
    /// Exposed so callers / tests can assert the books are closed
    /// regardless of the input numbers.
    pub fn closure_residual_w(&self) -> f64 {
        self.metabolic_net_w - self.total_loss_w() - self.storage_w
    }
}

/// Resolve the heat balance `M - W = R + C + E + S` for a body in an
/// environment, with a given metabolism and sweat (evaporative) state.
///
/// The three loss terms are computed from first principles:
///
/// - `R` = [`Environment::radiative_power`] at the skin temperature,
/// - `C` = [`Environment::convective_power`] at the skin temperature,
/// - `E` = [`Sweat::evaporative_power`],
///
/// and the storage `S` is whatever is left of `M - W` after those
/// losses, so the returned [`HeatBalance`] always satisfies the
/// identity exactly.
pub fn heat_balance(
    body: &Body,
    env: &Environment,
    metabolism: &Metabolism,
    sweat: &Sweat,
) -> HeatBalance {
    let metabolic_net_w = metabolism.net_heat_production();
    let radiation_w = env.radiative_power(body.skin_temp_c, body.surface_area_m2);
    let convection_w = env.convective_power(body.skin_temp_c, body.surface_area_m2);
    let evaporation_w = sweat.evaporative_power();
    let storage_w = metabolic_net_w - (radiation_w + convection_w + evaporation_w);
    HeatBalance {
        metabolic_net_w,
        radiation_w,
        convection_w,
        evaporation_w,
        storage_w,
    }
}

/// Core-temperature change (K, equivalently °C) produced by storing
/// `net_heat_w` watts for `dt_s` seconds in a body, by the
/// lumped-capacitance relation `dT = Q / (m c) = (S * dt) / (m c)`.
///
/// # Errors
///
/// Returns [`ThermoregError::NotFinite`] if `net_heat_w` or `dt_s` is
/// non-finite, and [`ThermoregError::OutOfRange`] if `dt_s` is
/// negative. (`net_heat_w` may be negative — cooling.)
pub fn core_temp_change(body: &Body, net_heat_w: f64, dt_s: f64) -> Result<f64, ThermoregError> {
    let net_heat_w = finite("net_heat_w", net_heat_w)?;
    let dt_s = non_negative("dt_s", dt_s)?;
    Ok(net_heat_w * dt_s / body.heat_capacity())
}

/// Advance the body's core temperature by one explicit-Euler step of
/// `dt_s` seconds under the current balance, returning the updated
/// [`Body`].
///
/// Convenience composition of [`heat_balance`] →
/// [`core_temp_change`] → [`Body::with_core_shifted`]. Only the core
/// temperature is advanced; the skin temperature and the input rates
/// are held fixed across the step (the single-node assumption).
///
/// # Errors
///
/// Propagates the validation from [`core_temp_change`].
pub fn step_core_temp(
    body: &Body,
    env: &Environment,
    metabolism: &Metabolism,
    sweat: &Sweat,
    dt_s: f64,
) -> Result<Body, ThermoregError> {
    let balance = heat_balance(body, env, metabolism, sweat);
    let delta = core_temp_change(body, balance.net_heat_w(), dt_s)?;
    Ok(body.with_core_shifted(delta))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::body::{BODY_SPECIFIC_HEAT, LATENT_HEAT_SWEAT};

    fn approx(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    fn body() -> Body {
        Body::new(70.0, BODY_SPECIFIC_HEAT, 1.8, 37.0, 33.0).unwrap()
    }

    #[test]
    fn metabolism_rejects_negatives() {
        assert!(Metabolism::new(-1.0, 0.0).is_err());
        assert!(Metabolism::new(100.0, -1.0).is_err());
    }

    #[test]
    fn net_heat_production_is_m_minus_w() {
        let m = Metabolism::new(300.0, 80.0).unwrap();
        assert!(approx(m.net_heat_production(), 220.0, 1e-12));
    }

    #[test]
    fn balance_always_closes_exactly() {
        // Arbitrary, non-equilibrium inputs: the residual must still
        // be zero because S is defined as the closing term.
        let env = Environment::new(18.0, 18.0, 5.0, 0.95).unwrap();
        let met = Metabolism::new(250.0, 40.0).unwrap();
        let sweat = Sweat::from_rate(3.0e-5).unwrap();
        let bal = heat_balance(&body(), &env, &met, &sweat);
        assert!(approx(bal.closure_residual_w(), 0.0, 1e-9));
    }

    #[test]
    fn at_balance_storage_is_zero() {
        // Choose a sweat rate that makes losses exactly equal M - W,
        // so S = 0. First compute the sensible losses, then size E.
        let env = Environment::new(22.0, 22.0, 4.0, 0.95).unwrap();
        let b = body();
        let r = env.radiative_power(b.skin_temp_c, b.surface_area_m2);
        let c = env.convective_power(b.skin_temp_c, b.surface_area_m2);
        let met = Metabolism::new(200.0, 0.0).unwrap();
        // Required evaporative power to balance: E = (M - W) - (R + C).
        let needed_e = met.net_heat_production() - (r + c);
        assert!(needed_e > 0.0, "sanity: need positive sweat here");
        let rate = needed_e / LATENT_HEAT_SWEAT;
        let sweat = Sweat::from_rate(rate).unwrap();
        let bal = heat_balance(&b, &env, &met, &sweat);
        assert!(bal.is_balanced(1e-6), "storage = {}", bal.storage_w);
        assert!(approx(bal.storage_w, 0.0, 1e-6));
    }

    #[test]
    fn positive_net_heat_raises_core_temp() {
        // Hot, humid, exercising, no sweat evaporation: big positive S.
        let env = Environment::new(35.0, 35.0, 4.0, 0.95).unwrap();
        let met = Metabolism::new(500.0, 0.0).unwrap();
        let sweat = Sweat::none();
        let bal = heat_balance(&body(), &env, &met, &sweat);
        assert!(
            bal.storage_w > 0.0,
            "expected heat gain, got {}",
            bal.storage_w
        );
        let after = step_core_temp(&body(), &env, &met, &sweat, 60.0).unwrap();
        assert!(after.core_temp_c > body().core_temp_c);
    }

    #[test]
    fn negative_net_heat_lowers_core_temp() {
        // Cold room, resting, light sweat: strong net loss.
        let env = Environment::new(5.0, 5.0, 6.0, 0.95).unwrap();
        let met = Metabolism::resting();
        let sweat = Sweat::none();
        let bal = heat_balance(&body(), &env, &met, &sweat);
        assert!(
            bal.storage_w < 0.0,
            "expected heat loss, got {}",
            bal.storage_w
        );
        let after = step_core_temp(&body(), &env, &met, &sweat, 60.0).unwrap();
        assert!(after.core_temp_c < body().core_temp_c);
    }

    #[test]
    fn higher_sweat_rate_increases_evaporative_cooling_and_lowers_storage() {
        let env = Environment::new(30.0, 30.0, 4.0, 0.95).unwrap();
        let met = Metabolism::new(300.0, 0.0).unwrap();
        let lo = heat_balance(&body(), &env, &met, &Sweat::from_rate(2.0e-5).unwrap());
        let hi = heat_balance(&body(), &env, &met, &Sweat::from_rate(8.0e-5).unwrap());
        // More sweat -> more evaporative loss ...
        assert!(hi.evaporation_w > lo.evaporation_w);
        // ... and therefore less heat stored (the body is cooler).
        assert!(hi.storage_w < lo.storage_w);
    }

    #[test]
    fn core_temp_change_equals_q_over_mc() {
        // dT = Q / (m c). With Q = S * dt: pick S = 244.44 W, dt = 1000 s,
        // m c = 70 * 3492 = 244_440 J/K  ->  dT = 244.44*1000/244440 = 1.0 K.
        let b = body();
        let dt = core_temp_change(&b, 244.44, 1000.0).unwrap();
        assert!(approx(dt, 1.0, 1e-9), "got {dt}");
    }

    #[test]
    fn core_temp_change_hand_value() {
        // 100 W stored for 1 hour in a 70 kg, c=3492 body:
        // dT = 100 * 3600 / (70*3492) = 360000 / 244440 = 1.4727... K.
        let b = body();
        let dt = core_temp_change(&b, 100.0, 3600.0).unwrap();
        assert!(approx(dt, 360_000.0 / 244_440.0, 1e-9), "got {dt}");
        assert!(approx(dt, 1.472_754, 1e-5), "got {dt}");
    }

    #[test]
    fn core_temp_change_is_linear_in_both_heat_and_time() {
        let b = body();
        let base = core_temp_change(&b, 50.0, 600.0).unwrap();
        let dbl_q = core_temp_change(&b, 100.0, 600.0).unwrap();
        let dbl_t = core_temp_change(&b, 50.0, 1200.0).unwrap();
        assert!(approx(dbl_q, 2.0 * base, 1e-12));
        assert!(approx(dbl_t, 2.0 * base, 1e-12));
    }

    #[test]
    fn zero_step_leaves_temperature_unchanged() {
        let b = body();
        assert!(approx(
            core_temp_change(&b, 250.0, 0.0).unwrap(),
            0.0,
            1e-15
        ));
    }

    #[test]
    fn core_temp_change_rejects_negative_dt_and_nonfinite() {
        let b = body();
        assert!(core_temp_change(&b, 100.0, -1.0).is_err());
        assert!(core_temp_change(&b, f64::NAN, 1.0).is_err());
        assert!(core_temp_change(&b, f64::INFINITY, 1.0).is_err());
    }

    #[test]
    fn step_then_change_are_consistent() {
        // step_core_temp must move the core by exactly core_temp_change
        // of the balance's net heat.
        let env = Environment::new(28.0, 28.0, 4.0, 0.95).unwrap();
        let met = Metabolism::new(350.0, 20.0).unwrap();
        let sweat = Sweat::from_rate(4.0e-5).unwrap();
        let b = body();
        let bal = heat_balance(&b, &env, &met, &sweat);
        let expected_delta = core_temp_change(&b, bal.net_heat_w(), 120.0).unwrap();
        let after = step_core_temp(&b, &env, &met, &sweat, 120.0).unwrap();
        assert!(approx(
            after.core_temp_c,
            b.core_temp_c + expected_delta,
            1e-12
        ));
    }

    #[test]
    fn total_loss_sums_the_three_channels() {
        let bal = HeatBalance {
            metabolic_net_w: 200.0,
            radiation_w: 60.0,
            convection_w: 50.0,
            evaporation_w: 40.0,
            storage_w: 50.0,
        };
        assert!(approx(bal.total_loss_w(), 150.0, 1e-12));
        assert!(approx(bal.net_heat_w(), 50.0, 1e-12));
        assert!(approx(bal.closure_residual_w(), 0.0, 1e-12));
    }
}
