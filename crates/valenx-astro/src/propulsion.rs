//! Native rocket-engine performance from ideal-rocket (de Laval nozzle)
//! theory.
//!
//! Given the chamber conditions (pressure, temperature, combustion-gas
//! properties) and the nozzle geometry (throat area, expansion ratio),
//! this derives the engine's **thrust** and **specific impulse** at any
//! ambient pressure — the numbers a [`crate::vehicle::Stage`] otherwise
//! takes as hand-entered inputs. So a stage can now be built from an
//! *engine design* instead of guessed performance figures.
//!
//! The model is the textbook isentropic-nozzle / ideal-rocket theory:
//! characteristic velocity `c*`, the Vandenkerckhove function `Γ`,
//! exit-Mach from the area-ratio relation, exit pressure, and the
//! momentum + pressure thrust terms. It assumes a calorically-perfect
//! gas, frozen composition, full expansion (no separation), and 100 %
//! combustion / nozzle efficiency — a real first-order engine model,
//! not a finite-rate-chemistry combustion simulation.

use serde::{Deserialize, Serialize};

use crate::constants::G0;
use crate::error::AstroError;

/// Universal gas constant (J/(kmol·K)).
const R_UNIVERSAL: f64 = 8_314.462_618;

/// Sea-level reference pressure (Pa).
const SEA_LEVEL_PRESSURE: f64 = 101_325.0;

/// A liquid-rocket-engine design point.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct EngineDesign {
    /// Chamber stagnation pressure (Pa).
    pub chamber_pressure: f64,
    /// Chamber (combustion) temperature (K).
    pub chamber_temperature: f64,
    /// Ratio of specific heats of the combustion gas (dimensionless).
    pub gamma: f64,
    /// Mean molar mass of the combustion gas (kg/kmol, i.e. g/mol).
    pub molar_mass: f64,
    /// Nozzle throat area (m²).
    pub throat_area: f64,
    /// Nozzle area ratio `ε = A_exit / A_throat` (dimensionless).
    pub expansion_ratio: f64,
}

/// Engine performance at a particular ambient pressure.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct EnginePerformance {
    /// Thrust (N).
    pub thrust: f64,
    /// Specific impulse (s).
    pub isp: f64,
    /// Propellant mass-flow rate (kg/s).
    pub mass_flow: f64,
    /// Characteristic velocity `c*` (m/s).
    pub c_star: f64,
    /// Thrust coefficient `C_f` (dimensionless).
    pub thrust_coefficient: f64,
    /// Nozzle exit pressure (Pa).
    pub exit_pressure: f64,
    /// Nozzle exit Mach number.
    pub exit_mach: f64,
}

impl EngineDesign {
    /// Specific gas constant of the combustion gas (J/(kg·K)).
    fn r_specific(&self) -> f64 {
        R_UNIVERSAL / self.molar_mass
    }

    /// Vandenkerckhove function `Γ(γ)` — the choked-throat mass-flow
    /// coefficient.
    fn gamma_function(&self) -> f64 {
        let g = self.gamma;
        g.sqrt() * (2.0 / (g + 1.0)).powf((g + 1.0) / (2.0 * (g - 1.0)))
    }

    /// Characteristic velocity `c* = √(R·T_c) / Γ` (m/s) — depends only
    /// on the chamber, not the nozzle.
    ///
    /// # Errors
    ///
    /// Returns [`AstroError::InvalidPropulsion`] if the design point is
    /// non-physical (see [`EngineDesign::validate`]), so a bad `gamma` /
    /// `molar_mass` / temperature yields an error rather than a NaN.
    pub fn c_star(&self) -> Result<f64, AstroError> {
        self.validate()?;
        Ok((self.r_specific() * self.chamber_temperature).sqrt() / self.gamma_function())
    }

    /// Choked propellant mass-flow rate `ṁ = p_c · A_t / c*` (kg/s).
    ///
    /// # Errors
    ///
    /// Returns [`AstroError::InvalidPropulsion`] for a non-physical
    /// design point (see [`EngineDesign::validate`]).
    pub fn mass_flow(&self) -> Result<f64, AstroError> {
        Ok(self.chamber_pressure * self.throat_area / self.c_star()?)
    }

    /// Supersonic exit Mach number implied by the area ratio, found by
    /// bisection on the isentropic area–Mach relation.
    ///
    /// # Errors
    ///
    /// Returns [`AstroError::InvalidPropulsion`] for a non-physical design
    /// point (see [`EngineDesign::validate`]) — a hand-built `gamma ≤ 1`
    /// makes the `(γ+1)/(2(γ-1))` exponent `Inf`/`NaN` and the solved Mach
    /// a silent `NaN`/`Inf`.
    pub fn exit_mach(&self) -> Result<f64, AstroError> {
        self.validate()?;
        let g = self.gamma;
        let target = self.expansion_ratio;
        // ε(M) = (1/M)·[ (2/(γ+1))·(1 + (γ-1)/2·M²) ]^((γ+1)/(2(γ-1)))
        let area_ratio = |m: f64| -> f64 {
            let exp = (g + 1.0) / (2.0 * (g - 1.0));
            (1.0 / m) * ((2.0 / (g + 1.0)) * (1.0 + (g - 1.0) / 2.0 * m * m)).powf(exp)
        };
        // ε is monotonically increasing for M > 1; bisect.
        let (mut lo, mut hi) = (1.0_f64, 100.0_f64);
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if area_ratio(mid) < target {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Ok(0.5 * (lo + hi))
    }

    /// Engine performance at a given ambient pressure (Pa).
    ///
    /// # Errors
    ///
    /// Returns [`AstroError::InvalidPropulsion`] for a non-physical
    /// design point (see [`EngineDesign::validate`]), so an invalid
    /// engine yields an error instead of a silently-NaN performance.
    pub fn performance(&self, ambient_pressure: f64) -> Result<EnginePerformance, AstroError> {
        self.validate()?;
        let g = self.gamma;
        let c_star = self.c_star()?;
        let mass_flow = self.mass_flow()?;
        let exit_mach = self.exit_mach()?;

        // Isentropic exit pressure and temperature from the chamber.
        let stagnation_ratio = 1.0 + (g - 1.0) / 2.0 * exit_mach * exit_mach;
        let exit_pressure = self.chamber_pressure * stagnation_ratio.powf(-g / (g - 1.0));
        let exit_temperature = self.chamber_temperature / stagnation_ratio;

        // Exit velocity and areas.
        let exit_velocity = exit_mach * (g * self.r_specific() * exit_temperature).sqrt();
        let exit_area = self.throat_area * self.expansion_ratio;

        // Thrust = momentum + pressure terms.
        let thrust = mass_flow * exit_velocity + (exit_pressure - ambient_pressure) * exit_area;
        let thrust_coefficient = thrust / (self.chamber_pressure * self.throat_area);
        let isp = thrust / (mass_flow * G0);

        Ok(EnginePerformance {
            thrust,
            isp,
            mass_flow,
            c_star,
            thrust_coefficient,
            exit_pressure,
            exit_mach,
        })
    }

    /// Vacuum performance (ambient pressure = 0).
    ///
    /// # Errors
    ///
    /// See [`EngineDesign::performance`].
    pub fn vacuum(&self) -> Result<EnginePerformance, AstroError> {
        self.performance(0.0)
    }

    /// Sea-level performance (ambient pressure = 101.325 kPa).
    ///
    /// # Errors
    ///
    /// See [`EngineDesign::performance`].
    pub fn sea_level(&self) -> Result<EnginePerformance, AstroError> {
        self.performance(SEA_LEVEL_PRESSURE)
    }

    /// Validate the design point.
    pub fn validate(&self) -> Result<(), AstroError> {
        for (field, value) in [
            ("chamber_pressure", self.chamber_pressure),
            ("chamber_temperature", self.chamber_temperature),
            ("molar_mass", self.molar_mass),
            ("throat_area", self.throat_area),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(AstroError::InvalidPropulsion {
                    index: usize::MAX,
                    field,
                    value,
                });
            }
        }
        if !self.gamma.is_finite() || self.gamma <= 1.0 {
            return Err(AstroError::InvalidPropulsion {
                index: usize::MAX,
                field: "gamma",
                value: self.gamma,
            });
        }
        if !self.expansion_ratio.is_finite() || self.expansion_ratio < 1.0 {
            return Err(AstroError::InvalidPropulsion {
                index: usize::MAX,
                field: "expansion_ratio",
                value: self.expansion_ratio,
            });
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A kerolox-class gas-generator engine design point.
    fn kerolox() -> EngineDesign {
        EngineDesign {
            chamber_pressure: 9.7e6, // ~97 bar
            chamber_temperature: 3_500.0,
            gamma: 1.2,
            molar_mass: 22.0,
            throat_area: 0.05,
            expansion_ratio: 16.0,
        }
    }

    #[test]
    fn c_star_in_physical_range() {
        // Kerolox c* is ~1700–1900 m/s.
        let c = kerolox().c_star().expect("valid engine");
        assert!((1_600.0..2_000.0).contains(&c), "c* = {c}");
    }

    #[test]
    fn area_ratio_solver_recovers_known_mach() {
        // For γ = 1.4, an area ratio of 1.6875 corresponds to Mach 2.0
        // (a standard isentropic-table value).
        let d = EngineDesign {
            gamma: 1.4,
            expansion_ratio: 1.6875,
            ..kerolox()
        };
        let me = d.exit_mach().expect("valid engine");
        assert!((me - 2.0).abs() < 1e-3, "Me = {me}");
    }

    #[test]
    fn vacuum_isp_beats_sea_level_and_is_physical() {
        let d = kerolox();
        let vac = d.vacuum().expect("valid engine");
        let sl = d.sea_level().expect("valid engine");
        // Vacuum Isp exceeds sea-level Isp (no back-pressure loss).
        assert!(vac.isp > sl.isp, "vac {} <= sl {}", vac.isp, sl.isp);
        // Kerolox vacuum Isp lands in the ~300–360 s ballpark.
        assert!((280.0..360.0).contains(&vac.isp), "Isp_vac = {}", vac.isp);
        // Thrust is positive at both conditions.
        assert!(vac.thrust > 0.0 && sl.thrust > 0.0);
        // Vacuum thrust exceeds sea-level thrust by the pressure term.
        assert!(vac.thrust > sl.thrust);
    }

    #[test]
    fn thrust_equals_cf_times_pc_at() {
        let d = kerolox();
        let p = d.vacuum().expect("valid engine");
        let reconstructed = p.thrust_coefficient * d.chamber_pressure * d.throat_area;
        assert!((reconstructed - p.thrust).abs() / p.thrust < 1e-9);
    }

    #[test]
    fn mass_flow_matches_c_star_definition() {
        let d = kerolox();
        let p = d.vacuum().expect("valid engine");
        let expected = d.chamber_pressure * d.throat_area / d.c_star().expect("valid engine");
        assert!((p.mass_flow - expected).abs() / expected < 1e-9);
    }

    #[test]
    fn validation_rejects_bad_designs() {
        let mut d = kerolox();
        assert!(d.validate().is_ok());
        d.gamma = 1.0; // must be > 1
        assert!(d.validate().is_err());
    }

    #[test]
    fn performance_methods_error_on_invalid_design_instead_of_nan() {
        // gamma <= 1 and molar_mass = 0 used to flow straight into the
        // performance math and emit silent NaN/Inf. They must now error.
        let bad_gamma = EngineDesign {
            gamma: 1.0,
            ..kerolox()
        };
        assert!(bad_gamma.c_star().is_err());
        assert!(bad_gamma.mass_flow().is_err());
        // exit_mach: gamma = 1 -> (g+1)/(2(g-1)) exponent is Inf -> NaN Mach.
        assert!(bad_gamma.exit_mach().is_err());
        assert!(bad_gamma.performance(0.0).is_err());
        assert!(bad_gamma.vacuum().is_err());
        assert!(bad_gamma.sea_level().is_err());

        let zero_molar = EngineDesign {
            molar_mass: 0.0,
            ..kerolox()
        };
        assert!(zero_molar.c_star().is_err());
        assert!(zero_molar.vacuum().is_err());
    }
}
