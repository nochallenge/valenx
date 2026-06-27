//! The pressure equation of state (EOS).
//!
//! Weakly-compressible SPH does not solve for pressure implicitly; it reads
//! pressure straight off the local density via a stiff equation of state. Two
//! standard choices are provided:
//!
//! * **Ideal-gas (linearised) form** — `p = k · (ρ − ρ₀)`, the "gas" pressure of
//!   Müller et al. 2003 (a Cole / Desbrun-style linear law). `k` is a gas
//!   stiffness constant. Simple and well-behaved; pressure goes negative when a
//!   particle is below rest density, which produces a (cohesive) attractive
//!   force.
//!
//! * **Tait form** — `p = B · ((ρ/ρ₀)^γ − 1)`, the Batchelor / Monaghan weakly-
//!   compressible law for liquids, with `γ = 7` and `B = ρ₀ c² / γ` where `c` is
//!   an artificial speed of sound. Stiffer near rest density, so it resists
//!   compression more sharply for the same sound speed.
//!
//! Both clamp negative pressure to zero **optionally** — left to the solver /
//! caller, since the cohesive negative-pressure behaviour is sometimes wanted.
//! Neither divides by anything that can be zero: `ρ₀ > 0` is validated at
//! construction and `ρ` is the (strictly positive) smoothed density the solver
//! always passes.

use crate::error::FluidError;

/// A pressure equation of state mapping density `ρ` to pressure `p`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EquationOfState {
    /// Linearised ideal-gas law `p = k · (ρ − ρ₀)`.
    IdealGas {
        /// Rest (reference) density `ρ₀` (kg/m³), strictly positive.
        rest_density: f64,
        /// Gas stiffness constant `k` (Pa·m³/kg), strictly positive.
        stiffness: f64,
    },
    /// Tait liquid law `p = B · ((ρ/ρ₀)^γ − 1)` with `B = ρ₀ c² / γ`.
    Tait {
        /// Rest density `ρ₀` (kg/m³), strictly positive.
        rest_density: f64,
        /// Artificial speed of sound `c` (m/s), strictly positive.
        sound_speed: f64,
        /// Adiabatic exponent `γ` (dimensionless), `> 0` (commonly `7`).
        gamma: f64,
    },
}

impl EquationOfState {
    /// A linearised ideal-gas EOS.
    ///
    /// # Errors
    /// [`FluidError::InvalidConfig`] if `rest_density` or `stiffness` is not
    /// finite and strictly positive.
    pub fn ideal_gas(rest_density: f64, stiffness: f64) -> Result<Self, FluidError> {
        check_positive(rest_density, "rest density")?;
        check_positive(stiffness, "gas stiffness")?;
        Ok(Self::IdealGas {
            rest_density,
            stiffness,
        })
    }

    /// A Tait (weakly-compressible liquid) EOS with the given exponent.
    ///
    /// # Errors
    /// [`FluidError::InvalidConfig`] if `rest_density`, `sound_speed`, or
    /// `gamma` is not finite and strictly positive.
    pub fn tait(rest_density: f64, sound_speed: f64, gamma: f64) -> Result<Self, FluidError> {
        check_positive(rest_density, "rest density")?;
        check_positive(sound_speed, "sound speed")?;
        check_positive(gamma, "gamma")?;
        Ok(Self::Tait {
            rest_density,
            sound_speed,
            gamma,
        })
    }

    /// A Tait EOS with the conventional `γ = 7` for water-like liquids.
    ///
    /// # Errors
    /// [`FluidError::InvalidConfig`] if `rest_density` or `sound_speed` is not
    /// finite and strictly positive.
    pub fn tait_water(rest_density: f64, sound_speed: f64) -> Result<Self, FluidError> {
        Self::tait(rest_density, sound_speed, 7.0)
    }

    /// The rest (reference) density `ρ₀`.
    #[must_use]
    pub fn rest_density(&self) -> f64 {
        match *self {
            Self::IdealGas { rest_density, .. } | Self::Tait { rest_density, .. } => rest_density,
        }
    }

    /// The pressure for a density `density` (`ρ`).
    ///
    /// May return a negative value when `ρ < ρ₀` (a cohesive pressure); use
    /// [`Self::pressure_clamped`] to floor it at zero. `density` is assumed to be
    /// the strictly positive smoothed density; the formula divides only by the
    /// validated, strictly positive `ρ₀`.
    #[must_use]
    pub fn pressure(&self, density: f64) -> f64 {
        match *self {
            Self::IdealGas {
                rest_density,
                stiffness,
            } => stiffness * (density - rest_density),
            Self::Tait {
                rest_density,
                sound_speed,
                gamma,
            } => {
                let b = rest_density * sound_speed * sound_speed / gamma;
                b * ((density / rest_density).powf(gamma) - 1.0)
            }
        }
    }

    /// The pressure for a density, floored at zero (drops the cohesive
    /// negative-pressure branch).
    #[must_use]
    pub fn pressure_clamped(&self, density: f64) -> f64 {
        self.pressure(density).max(0.0)
    }
}

fn check_positive(value: f64, what: &str) -> Result<(), FluidError> {
    if !(value.is_finite() && value > 0.0) {
        return Err(FluidError::InvalidConfig(format!(
            "{what} must be finite and > 0, got {value}"
        )));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_bad_parameters() {
        assert!(EquationOfState::ideal_gas(0.0, 1.0).is_err());
        assert!(EquationOfState::ideal_gas(1000.0, -1.0).is_err());
        assert!(EquationOfState::tait(1000.0, 0.0, 7.0).is_err());
        assert!(EquationOfState::tait(1000.0, 10.0, 0.0).is_err());
        assert!(EquationOfState::ideal_gas(1000.0, 100.0).is_ok());
        assert!(EquationOfState::tait_water(1000.0, 50.0).is_ok());
    }

    #[test]
    fn ideal_gas_is_zero_at_rest_density_and_signed_around_it() {
        let eos = EquationOfState::ideal_gas(1000.0, 3.0).unwrap();
        assert!((eos.pressure(1000.0)).abs() < 1e-12);
        assert!(eos.pressure(1010.0) > 0.0); // compressed ⇒ pushes out
        assert!(eos.pressure(990.0) < 0.0); // rarefied ⇒ cohesive pull
        assert_eq!(eos.pressure_clamped(990.0), 0.0);
    }

    #[test]
    fn tait_is_zero_at_rest_density_and_stiffer_than_ideal_when_compressed() {
        let rho0 = 1000.0;
        let c = 50.0;
        let tait = EquationOfState::tait_water(rho0, c).unwrap();
        assert!((tait.pressure(rho0)).abs() < 1e-9);

        // At ρ = ρ₀ the slope dp/dρ of Tait equals c² (= ρ₀ c² / γ · γ / ρ₀).
        // A linear gas with the same slope would have k = c². Past rest density
        // Tait's γ-th power makes it rise faster than that linear law.
        let k = c * c;
        let gas = EquationOfState::ideal_gas(rho0, k).unwrap();
        let rho = 1050.0;
        assert!(
            tait.pressure(rho) > gas.pressure(rho),
            "Tait should be stiffer than the slope-matched linear gas when compressed"
        );
    }

    #[test]
    fn pressure_increases_monotonically_with_density() {
        let eos = EquationOfState::tait_water(1000.0, 50.0).unwrap();
        let p1 = eos.pressure(1000.0);
        let p2 = eos.pressure(1010.0);
        let p3 = eos.pressure(1020.0);
        assert!(p1 < p2 && p2 < p3);
    }
}
