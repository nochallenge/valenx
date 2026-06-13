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

/// Inputs for the first-order regenerative-cooling heat balance — the wall
/// and coolant parameters the Bartz model needs on top of the engine's own
/// gas properties.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct CoolingInputs {
    /// Gas-side wall temperature held by the cooling (K). The throat heat
    /// flux is driven by `T_aw − T_wall`.
    pub gas_wall_temperature: f64,
    /// Maximum throat heat flux the regen circuit can sustain (W/m²). The
    /// headline cooling margin is this divided by the actual throat flux, so
    /// it tightens as chamber pressure (and thus flux) rises.
    pub max_heat_flux: f64,
    /// Throat wall radius of curvature as a multiple of the throat diameter
    /// (the Bartz `D_t/R_c` term; ~0.5–1.5 typically).
    pub throat_curvature_ratio: f64,
    /// Coolant specific heat capacity (J/(kg·K)).
    pub coolant_specific_heat: f64,
    /// Fraction of the total propellant mass-flow routed as coolant (the
    /// fuel fraction in a regen circuit), 0–1.
    pub coolant_mass_fraction: f64,
    /// Heated chamber+throat length used to turn the throat flux into a
    /// total heat load (m) — a conservative first-order area.
    pub heated_length: f64,
}

impl Default for CoolingInputs {
    /// A representative RP-1 regenerative-cooling circuit.
    fn default() -> Self {
        Self {
            gas_wall_temperature: 800.0,
            max_heat_flux: 80.0e6, // advanced regen throat capability
            throat_curvature_ratio: 1.0,
            coolant_specific_heat: 2_000.0, // ~RP-1
            coolant_mass_fraction: 0.3,
            heated_length: 0.5,
        }
    }
}

/// First-order regenerative-cooling heat balance at the throat.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct CoolingPerformance {
    /// Bartz gas-side heat-transfer coefficient at the throat (W/(m²·K)).
    pub gas_heat_transfer_coeff: f64,
    /// Adiabatic (recovery) wall temperature at the throat (K).
    pub adiabatic_wall_temperature: f64,
    /// Throat heat flux `q = h_g·(T_aw − T_wall)` (W/m²).
    pub throat_heat_flux: f64,
    /// Total heat load `q · heated area` (W) — conservative (throat flux
    /// applied over the whole heated wall).
    pub heat_load: f64,
    /// Bulk coolant temperature rise `q·A / (ṁ_c·c_p)` (K).
    pub coolant_temperature_rise: f64,
    /// Cooling margin `max_heat_flux / throat_heat_flux` (> 1 ⇒ the throat
    /// flux is within the regen circuit's capability). Falls as chamber
    /// pressure rises — the real high-`p_c` cooling constraint.
    pub cooling_margin: f64,
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

    /// First-order regenerative-cooling heat balance at the throat, from a
    /// Bartz heat-transfer estimate.
    ///
    /// The gas-side film coefficient is the classic Bartz correlation
    /// `h_g = (0.026/D_t^0.2)·(μ^0.2·c_p/Pr^0.6)·(p_c/c*)^0.8·(D_t/R_c)^0.1`,
    /// the throat heat flux is `q = h_g·(T_aw − T_wall)` against the
    /// recovery temperature `T_aw`, and the cooling margin compares that flux
    /// to the regen circuit's sustainable maximum. It is a first-order design
    /// estimate — calorically-perfect gas, Eucken Prandtl number, a unity
    /// transport-property correction — not a conjugate heat-transfer solve.
    ///
    /// # Errors
    ///
    /// Returns [`AstroError::InvalidPropulsion`] for a non-physical engine
    /// (see [`EngineDesign::validate`]).
    pub fn cooling(&self, inputs: &CoolingInputs) -> Result<CoolingPerformance, AstroError> {
        self.validate()?;
        let g = self.gamma;
        let r = self.r_specific();
        // Gas transport properties (first-order estimates).
        let cp = g * r / (g - 1.0); // calorically-perfect c_p
        let pr = 4.0 * g / (9.0 * g - 5.0); // Eucken Prandtl number
                                            // Bartz viscosity correlation, SI-adapted: μ ≈ 1.184e-7·M^0.5·T^0.6.
        let mu = 1.184e-7 * self.molar_mass.sqrt() * self.chamber_temperature.powf(0.6);
        let d_t = (4.0 * self.throat_area / std::f64::consts::PI).sqrt();
        let c_star = self.c_star()?;
        let r_c = inputs.throat_curvature_ratio.max(1e-6) * d_t;

        // Bartz gas-side film coefficient (transport correction σ ≈ 1).
        let h_g = (0.026 / d_t.powf(0.2))
            * (mu.powf(0.2) * cp / pr.powf(0.6))
            * (self.chamber_pressure / c_star).powf(0.8)
            * (d_t / r_c).powf(0.1);

        // Recovery (adiabatic wall) temperature at the throat (M = 1),
        // turbulent recovery factor r = Pr^(1/3): static throat temp scaled
        // back up by the recovery of the kinetic energy.
        let recovery = pr.cbrt();
        let t_aw =
            self.chamber_temperature * (1.0 + recovery * (g - 1.0) / 2.0) / (1.0 + (g - 1.0) / 2.0);

        let throat_heat_flux = h_g * (t_aw - inputs.gas_wall_temperature);

        // Conservative heat load: peak throat flux over the heated wall area.
        let heated_area = std::f64::consts::PI * d_t * inputs.heated_length;
        let heat_load = throat_heat_flux * heated_area;

        // Bulk coolant temperature rise from that load.
        let m_dot = self.mass_flow()?;
        let coolant_flow = inputs.coolant_mass_fraction.clamp(0.0, 1.0) * m_dot;
        let coolant_temperature_rise = if coolant_flow > 0.0 && inputs.coolant_specific_heat > 0.0 {
            heat_load / (coolant_flow * inputs.coolant_specific_heat)
        } else {
            f64::INFINITY
        };

        // Headline margin: sustainable flux ÷ actual flux (falls with p_c).
        let cooling_margin = if throat_heat_flux > 0.0 {
            inputs.max_heat_flux / throat_heat_flux
        } else {
            f64::INFINITY
        };

        Ok(CoolingPerformance {
            gas_heat_transfer_coeff: h_g,
            adiabatic_wall_temperature: t_aw,
            throat_heat_flux,
            heat_load,
            coolant_temperature_rise,
            cooling_margin,
        })
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

/// The best engine found by [`optimize_engine`].
#[derive(Debug, Clone, Copy)]
pub struct EngineOptimum {
    /// The winning design point.
    pub design: EngineDesign,
    /// Its sea-level performance (the optimized objective is sea-level Isp).
    pub sea_level: EnginePerformance,
    /// Its vacuum performance.
    pub vacuum: EnginePerformance,
    /// Its cooling balance (margin ≥ the requested target).
    pub cooling: CoolingPerformance,
    /// Total candidates evaluated.
    pub evaluations: usize,
    /// How many candidates were feasible (cooling margin ≥ target).
    pub feasible_count: usize,
}

/// Search chamber-pressure × expansion-ratio space for the engine with the
/// highest **sea-level specific impulse** whose cooling margin stays at or
/// above `target_margin`. A deterministic grid over `pc_range` × `eps_range`
/// (`grid` points per axis), reusing `base` for the gas properties and throat
/// area. Returns `None` if no candidate is feasible.
///
/// This is the honest "let the optimizer design the engine" loop: sea-level
/// Isp rises with chamber pressure (more expansion before back-pressure
/// losses), but the Bartz cooling margin falls as `p_c^0.8`, so the optimum is
/// pushed up against the cooling limit — exactly the trade a preliminary-design
/// tool makes. It is NOT generative-geometry / combustion-CFD optimization.
pub fn optimize_engine(
    base: &EngineDesign,
    cooling: &CoolingInputs,
    target_margin: f64,
    pc_range: (f64, f64),
    eps_range: (f64, f64),
    grid: usize,
) -> Option<EngineOptimum> {
    let grid = grid.max(2);
    let (pc_lo, pc_hi) = pc_range;
    let (eps_lo, eps_hi) = eps_range;
    let mut best: Option<EngineOptimum> = None;
    let mut best_isp = f64::NEG_INFINITY;
    let mut evaluations = 0usize;
    let mut feasible_count = 0usize;

    for i in 0..grid {
        let pc = pc_lo + (i as f64 / (grid - 1) as f64) * (pc_hi - pc_lo);
        for j in 0..grid {
            let eps = eps_lo + (j as f64 / (grid - 1) as f64) * (eps_hi - eps_lo);
            evaluations += 1;
            let design = EngineDesign {
                chamber_pressure: pc,
                expansion_ratio: eps,
                ..*base
            };
            let (Ok(sl), Ok(vac), Ok(cool)) =
                (design.sea_level(), design.vacuum(), design.cooling(cooling))
            else {
                continue;
            };
            // Feasible = real operating point (positive thrust) whose cooling
            // margin clears the target.
            if !sl.isp.is_finite() || sl.thrust <= 0.0 {
                continue;
            }
            if !cool.cooling_margin.is_finite() || cool.cooling_margin < target_margin {
                continue;
            }
            feasible_count += 1;
            if sl.isp > best_isp {
                best_isp = sl.isp;
                best = Some(EngineOptimum {
                    design,
                    sea_level: sl,
                    vacuum: vac,
                    cooling: cool,
                    evaluations: 0,    // patched in after the search
                    feasible_count: 0, // patched in after the search
                });
            }
        }
    }

    best.map(|mut b| {
        b.evaluations = evaluations;
        b.feasible_count = feasible_count;
        b
    })
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

    #[test]
    fn bartz_throat_heat_flux_is_physical_for_kerolox() {
        let d = kerolox();
        let c = d.cooling(&CoolingInputs::default()).expect("valid engine");
        // Recovery temperature sits just below the chamber stagnation temp
        // and well above the cooled wall.
        assert!(c.adiabatic_wall_temperature < d.chamber_temperature);
        assert!(c.adiabatic_wall_temperature > 800.0);
        // Throat heat flux for a ~100-bar kerolox engine is tens of MW/m²
        // (the SSME throat is ~160 MW/m²); this broad band catches unit
        // errors without over-fitting the first-order model.
        assert!(
            (5.0e6..200.0e6).contains(&c.throat_heat_flux),
            "q_throat = {} W/m²",
            c.throat_heat_flux
        );
        assert!(c.gas_heat_transfer_coeff > 0.0);
        assert!(c.heat_load > 0.0);
        assert!(c.coolant_temperature_rise > 0.0 && c.coolant_temperature_rise.is_finite());
        assert!(c.cooling_margin > 0.0 && c.cooling_margin.is_finite());
    }

    #[test]
    fn raising_chamber_pressure_raises_flux_and_tightens_margin() {
        // Bartz h_g ∝ p_c^0.8, so the throat flux climbs with chamber
        // pressure while the sustainable-flux margin falls — the real
        // high-p_c cooling constraint that bounds how hard you can run.
        let inp = CoolingInputs::default();
        let lo = EngineDesign {
            chamber_pressure: 5.0e6,
            ..kerolox()
        }
        .cooling(&inp)
        .unwrap();
        let hi = EngineDesign {
            chamber_pressure: 15.0e6,
            ..kerolox()
        }
        .cooling(&inp)
        .unwrap();
        assert!(
            hi.throat_heat_flux > lo.throat_heat_flux,
            "flux should rise with p_c: {} -> {}",
            lo.throat_heat_flux,
            hi.throat_heat_flux
        );
        assert!(
            hi.cooling_margin < lo.cooling_margin,
            "margin should fall with p_c: {} -> {}",
            lo.cooling_margin,
            hi.cooling_margin
        );
    }

    #[test]
    fn more_coolant_flow_lowers_the_coolant_temperature_rise() {
        let d = kerolox();
        let lean = CoolingInputs {
            coolant_mass_fraction: 0.15,
            ..CoolingInputs::default()
        };
        let rich = CoolingInputs {
            coolant_mass_fraction: 0.45,
            ..CoolingInputs::default()
        };
        let dt_lean = d.cooling(&lean).unwrap().coolant_temperature_rise;
        let dt_rich = d.cooling(&rich).unwrap().coolant_temperature_rise;
        assert!(
            dt_rich < dt_lean,
            "more coolant ⇒ smaller ΔT: {dt_lean} -> {dt_rich}"
        );
    }

    #[test]
    fn cooling_rejects_invalid_design() {
        let bad = EngineDesign {
            gamma: 1.0,
            ..kerolox()
        };
        assert!(bad.cooling(&CoolingInputs::default()).is_err());
    }

    #[test]
    fn optimizer_finds_a_feasible_cooled_engine() {
        let opt = optimize_engine(
            &kerolox(),
            &CoolingInputs::default(),
            1.5,
            (4.0e6, 20.0e6),
            (8.0, 40.0),
            24,
        )
        .expect("a feasible engine exists in range");
        assert_eq!(opt.evaluations, 24 * 24);
        assert!(opt.feasible_count > 0);
        assert!(
            opt.cooling.cooling_margin >= 1.5 - 1e-9,
            "margin {}",
            opt.cooling.cooling_margin
        );
        assert!(opt.sea_level.isp > 0.0 && opt.sea_level.thrust > 0.0);
        assert!(opt.vacuum.isp > opt.sea_level.isp);
    }

    #[test]
    fn tighter_cooling_target_forces_a_lower_chamber_pressure() {
        // Sea-level Isp wants high p_c, but a stricter cooling margin caps how
        // high p_c can go — so the optimum p_c drops as the target rises.
        let base = kerolox();
        let inp = CoolingInputs::default();
        let lax = optimize_engine(&base, &inp, 1.2, (4.0e6, 20.0e6), (8.0, 40.0), 24).unwrap();
        let strict = optimize_engine(&base, &inp, 2.5, (4.0e6, 20.0e6), (8.0, 40.0), 24).unwrap();
        assert!(
            strict.design.chamber_pressure <= lax.design.chamber_pressure,
            "strict p_c {} should be ≤ lax p_c {}",
            strict.design.chamber_pressure,
            lax.design.chamber_pressure
        );
        assert!(strict.cooling.cooling_margin >= 2.5 - 1e-9);
    }

    #[test]
    fn optimizer_returns_none_when_nothing_meets_the_target() {
        // An impossibly high margin target has no feasible point.
        let none = optimize_engine(
            &kerolox(),
            &CoolingInputs::default(),
            1_000.0,
            (18.0e6, 20.0e6),
            (8.0, 40.0),
            8,
        );
        assert!(none.is_none());
    }
}
