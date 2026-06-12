//! Launch-vehicle description: stages, propulsion, mass and aerodynamics.
//!
//! A [`Vehicle`] is a stack of [`Stage`]s (index 0 burns first) plus a
//! payload and an aerodynamic model. All masses are in kilograms,
//! thrust in newtons, specific impulse in seconds. Thrust and `Isp`
//! vary linearly with ambient pressure between their sea-level and
//! vacuum values, which is the standard first-order nozzle model; the
//! propellant mass-flow rate is held constant (set by the vacuum
//! rating), so a stage's burn time is `propellant / ṁ`.

use serde::{Deserialize, Serialize};

use crate::constants::G0;
use crate::error::AstroError;
use crate::propulsion::EngineDesign;

/// Sea-level reference pressure (Pa) used to scale thrust / `Isp`.
const SEA_LEVEL_PRESSURE: f64 = 101_325.0;

/// A single propulsive stage.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Stage {
    /// Human-readable name (e.g. `"booster"`).
    pub name: String,
    /// Inert (structural) mass that is carried then jettisoned at
    /// staging (kg).
    pub dry_mass: f64,
    /// Usable propellant mass (kg).
    pub propellant_mass: f64,
    /// Vacuum thrust (N).
    pub thrust_vac: f64,
    /// Sea-level thrust (N). For an upper stage that never fires in the
    /// atmosphere set this equal to `thrust_vac`.
    pub thrust_sl: f64,
    /// Vacuum specific impulse (s).
    pub isp_vac: f64,
    /// Sea-level specific impulse (s).
    pub isp_sl: f64,
}

impl Stage {
    /// Build a stage from a native [`EngineDesign`] (and an engine
    /// count), deriving the vacuum / sea-level thrust and `Isp` from
    /// ideal-rocket nozzle theory instead of hand-entered figures.
    ///
    /// `num_engines` scales the thrust (mass flow) but not the `Isp`.
    ///
    /// # Errors
    ///
    /// Returns [`AstroError::InvalidPropulsion`] if the engine design is
    /// non-physical (see [`EngineDesign::validate`]), so an invalid
    /// engine can no longer silently produce a NaN-rated stage.
    pub fn from_engine(
        name: impl Into<String>,
        dry_mass: f64,
        propellant_mass: f64,
        engine: &EngineDesign,
        num_engines: u32,
    ) -> Result<Self, AstroError> {
        let n = num_engines as f64;
        let vac = engine.vacuum()?;
        let sl = engine.sea_level()?;
        Ok(Self {
            name: name.into(),
            dry_mass,
            propellant_mass,
            thrust_vac: vac.thrust * n,
            thrust_sl: sl.thrust * n,
            isp_vac: vac.isp,
            isp_sl: sl.isp,
        })
    }

    /// Constant propellant mass-flow rate (kg/s), fixed by the vacuum
    /// rating: `ṁ = F_vac / (Isp_vac · g₀)`.
    ///
    /// # Errors
    ///
    /// Returns [`AstroError::InvalidPropulsion`] if `isp_vac` or
    /// `thrust_vac` is non-finite or non-positive — a hand-built [`Stage`]
    /// that bypassed [`Stage::from_engine`]'s validation could otherwise
    /// have `Isp_vac · g₀ = 0`/NaN and yield a silent `NaN` (0/0) or
    /// `Inf` (F_vac/0) mass flow.
    pub fn mass_flow(&self) -> Result<f64, AstroError> {
        for (field, value) in [("thrust_vac", self.thrust_vac), ("isp_vac", self.isp_vac)] {
            if !value.is_finite() || value <= 0.0 {
                return Err(AstroError::InvalidPropulsion {
                    index: usize::MAX,
                    field,
                    value,
                });
            }
        }
        Ok(self.thrust_vac / (self.isp_vac * G0))
    }

    /// Full-throttle burn time of this stage in isolation (s).
    ///
    /// # Errors
    ///
    /// Returns [`AstroError::InvalidPropulsion`] if `thrust_vac` or
    /// `isp_vac` is non-finite or non-positive — a hand-built [`Stage`]
    /// that bypassed [`Stage::from_engine`]'s validation could otherwise
    /// have a zero/NaN mass flow, making the burn time silently `NaN`
    /// (0/0) or `Inf` (propellant/0).
    pub fn burn_time(&self) -> Result<f64, AstroError> {
        for (field, value) in [("thrust_vac", self.thrust_vac), ("isp_vac", self.isp_vac)] {
            if !value.is_finite() || value <= 0.0 {
                return Err(AstroError::InvalidPropulsion {
                    index: usize::MAX,
                    field,
                    value,
                });
            }
        }
        let mdot = self.mass_flow()?;
        if !mdot.is_finite() || mdot <= 0.0 {
            return Err(AstroError::InvalidPropulsion {
                index: usize::MAX,
                field: "mass_flow",
                value: mdot,
            });
        }
        Ok(self.propellant_mass / mdot)
    }

    /// Thrust (N) at a given ambient pressure (Pa), linearly
    /// interpolated between the sea-level and vacuum ratings and
    /// clamped to the `[sea-level, vacuum]` envelope.
    pub fn thrust(&self, ambient_pressure: f64) -> f64 {
        let frac = (ambient_pressure / SEA_LEVEL_PRESSURE).clamp(0.0, 1.0);
        self.thrust_vac - frac * (self.thrust_vac - self.thrust_sl)
    }

    /// Ideal velocity change this stage delivers (m/s) when it has
    /// `upper_mass` kilograms stacked above it, via Tsiolkovsky:
    /// `Δv = Isp_vac · g₀ · ln(m₀/m_f)`.
    ///
    /// # Errors
    ///
    /// Returns [`AstroError::InvalidPropulsion`] if `isp_vac` is
    /// non-finite or non-positive, or [`AstroError::InvalidMass`] if the
    /// masses are non-finite or the burnout mass `m_f = upper_mass +
    /// dry_mass` is not strictly positive — either would otherwise make
    /// the `ln(m₀/m_f)` term a silent `NaN`/`Inf`. A hand-built [`Stage`]
    /// or out-of-range `upper_mass` is the way to hit this.
    pub fn ideal_delta_v(&self, upper_mass: f64) -> Result<f64, AstroError> {
        if !self.isp_vac.is_finite() || self.isp_vac <= 0.0 {
            return Err(AstroError::InvalidPropulsion {
                index: usize::MAX,
                field: "isp_vac",
                value: self.isp_vac,
            });
        }
        let m0 = upper_mass + self.dry_mass + self.propellant_mass;
        let mf = upper_mass + self.dry_mass;
        if !m0.is_finite() {
            return Err(AstroError::InvalidMass {
                index: usize::MAX,
                field: "initial_mass",
                value: m0,
            });
        }
        if !mf.is_finite() || mf <= 0.0 {
            return Err(AstroError::InvalidMass {
                index: usize::MAX,
                field: "burnout_mass",
                value: mf,
            });
        }
        Ok(self.isp_vac * G0 * (m0 / mf).ln())
    }

    fn validate(&self, index: usize) -> Result<(), AstroError> {
        for (field, value) in [
            ("dry_mass", self.dry_mass),
            ("propellant_mass", self.propellant_mass),
        ] {
            if !value.is_finite() || value < 0.0 {
                return Err(AstroError::InvalidMass {
                    index,
                    field,
                    value,
                });
            }
        }
        if self.dry_mass <= 0.0 {
            return Err(AstroError::InvalidMass {
                index,
                field: "dry_mass",
                value: self.dry_mass,
            });
        }
        for (field, value) in [
            ("thrust_vac", self.thrust_vac),
            ("thrust_sl", self.thrust_sl),
            ("isp_vac", self.isp_vac),
            ("isp_sl", self.isp_sl),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(AstroError::InvalidPropulsion {
                    index,
                    field,
                    value,
                });
            }
        }
        Ok(())
    }
}

/// Drag coefficient as a piecewise-linear function of Mach number.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DragModel {
    /// `(Mach, Cd)` breakpoints, ascending in Mach. Lookups below the
    /// first / above the last point clamp to the end values.
    pub cd_curve: Vec<(f64, f64)>,
}

impl DragModel {
    /// A generic slender-body transonic drag curve: low subsonic drag,
    /// a transonic rise peaking near Mach 1.2, then a supersonic
    /// fall-off. Reasonable default for a launch vehicle silhouette.
    pub fn generic_launch_vehicle() -> Self {
        Self {
            cd_curve: vec![
                (0.0, 0.30),
                (0.8, 0.32),
                (1.0, 0.55),
                (1.2, 0.60),
                (2.0, 0.42),
                (4.0, 0.28),
                (8.0, 0.22),
            ],
        }
    }

    /// Drag coefficient at a given Mach number (clamped at the table
    /// ends, linearly interpolated between breakpoints).
    pub fn cd(&self, mach: f64) -> f64 {
        let curve = &self.cd_curve;
        if curve.is_empty() {
            return 0.0;
        }
        if mach <= curve[0].0 {
            return curve[0].1;
        }
        let last = curve[curve.len() - 1];
        if mach >= last.0 {
            return last.1;
        }
        // Find the bracketing interval.
        for w in curve.windows(2) {
            let (m0, c0) = w[0];
            let (m1, c1) = w[1];
            if mach >= m0 && mach <= m1 {
                let t = (mach - m0) / (m1 - m0);
                return c0 + t * (c1 - c0);
            }
        }
        last.1
    }

    fn validate(&self) -> Result<(), AstroError> {
        if self.cd_curve.is_empty() {
            return Err(AstroError::InvalidAero("drag curve is empty"));
        }
        for &(m, c) in &self.cd_curve {
            if !m.is_finite() || !c.is_finite() || m < 0.0 || c < 0.0 {
                return Err(AstroError::InvalidAero("drag curve has invalid point"));
            }
        }
        // Must be ascending in Mach so the bracket search is valid.
        if self.cd_curve.windows(2).any(|w| w[1].0 <= w[0].0) {
            return Err(AstroError::InvalidAero("drag curve Mach not ascending"));
        }
        Ok(())
    }
}

/// A complete launch vehicle: a stack of stages plus payload and
/// aerodynamics.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Vehicle {
    /// Stages in firing order: index 0 is the first stage (booster).
    pub stages: Vec<Stage>,
    /// Payload mass carried all the way to orbit (kg).
    pub payload_mass: f64,
    /// Aerodynamic reference (frontal) area (m²).
    pub reference_area: f64,
    /// Mach-dependent drag model.
    pub drag: DragModel,
}

impl Vehicle {
    /// Liftoff gross mass (kg): payload plus every stage's wet mass.
    pub fn initial_mass(&self) -> f64 {
        self.payload_mass
            + self
                .stages
                .iter()
                .map(|s| s.dry_mass + s.propellant_mass)
                .sum::<f64>()
    }

    /// Mass stacked above stage `i` at the moment it ignites (kg): the
    /// payload plus every higher stage's full wet mass.
    pub fn upper_mass_above(&self, i: usize) -> f64 {
        self.payload_mass
            + self.stages[i + 1..]
                .iter()
                .map(|s| s.dry_mass + s.propellant_mass)
                .sum::<f64>()
    }

    /// Ideal total `Δv` budget (m/s): the sum of each stage's
    /// Tsiolkovsky `Δv` given the live mass above it. This is the
    /// vehicle's theoretical performance ceiling, independent of
    /// gravity / drag / steering losses.
    ///
    /// # Errors
    ///
    /// Propagates [`Stage::ideal_delta_v`]'s error if any stage has a
    /// non-physical `isp_vac` or mass (which would otherwise make the
    /// summed budget a silent `NaN`/`Inf`).
    pub fn ideal_delta_v(&self) -> Result<f64, AstroError> {
        let mut total = 0.0;
        for i in 0..self.stages.len() {
            total += self.stages[i].ideal_delta_v(self.upper_mass_above(i))?;
        }
        Ok(total)
    }

    /// Validate the whole vehicle, returning the first problem found.
    pub fn validate(&self) -> Result<(), AstroError> {
        if self.stages.is_empty() {
            return Err(AstroError::NoStages);
        }
        for (i, s) in self.stages.iter().enumerate() {
            s.validate(i)?;
        }
        if !self.payload_mass.is_finite() || self.payload_mass < 0.0 {
            return Err(AstroError::InvalidMass {
                index: usize::MAX,
                field: "payload_mass",
                value: self.payload_mass,
            });
        }
        if !self.reference_area.is_finite() || self.reference_area <= 0.0 {
            return Err(AstroError::InvalidAero("reference area must be positive"));
        }
        self.drag.validate()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn demo_stage() -> Stage {
        Stage {
            name: "demo".into(),
            dry_mass: 1_000.0,
            propellant_mass: 9_000.0,
            thrust_vac: 200_000.0,
            thrust_sl: 180_000.0,
            isp_vac: 300.0,
            isp_sl: 270.0,
        }
    }

    #[test]
    fn from_engine_derives_sane_performance() {
        let engine = EngineDesign {
            chamber_pressure: 9.7e6,
            chamber_temperature: 3_500.0,
            gamma: 1.2,
            molar_mass: 22.0,
            throat_area: 0.05,
            expansion_ratio: 16.0,
        };
        let s =
            Stage::from_engine("booster", 20_000.0, 400_000.0, &engine, 9).expect("valid engine");
        // Nine engines -> roughly nine times one engine's thrust.
        let one = engine.vacuum().expect("valid engine").thrust;
        assert!((s.thrust_vac - 9.0 * one).abs() / s.thrust_vac < 1e-9);
        // Isp is per-engine (not scaled by count) and physical.
        assert!((280.0..360.0).contains(&s.isp_vac));
        assert!(s.isp_vac > s.isp_sl);
        assert!(s.thrust_vac > s.thrust_sl);
        // The derived stage is a valid building block.
        let v = Vehicle {
            stages: vec![s],
            payload_mass: 1_000.0,
            reference_area: 10.0,
            drag: DragModel::generic_launch_vehicle(),
        };
        assert!(v.validate().is_ok());
    }

    #[test]
    fn from_engine_rejects_invalid_engine_instead_of_nan_stage() {
        // gamma <= 1 used to flow into the nozzle math and build a stage
        // with NaN thrust/Isp. It must now be a clean Err.
        let bad = EngineDesign {
            chamber_pressure: 9.7e6,
            chamber_temperature: 3_500.0,
            gamma: 1.0, // invalid
            molar_mass: 22.0,
            throat_area: 0.05,
            expansion_ratio: 16.0,
        };
        let r = Stage::from_engine("booster", 20_000.0, 400_000.0, &bad, 9);
        assert!(
            matches!(r, Err(AstroError::InvalidPropulsion { .. })),
            "invalid engine must be rejected, got {r:?}"
        );
    }

    #[test]
    fn tsiolkovsky_matches_closed_form() {
        let s = demo_stage();
        // No upper mass: m0 = 10 000, mf = 1 000 -> ratio 10.
        let dv = s.ideal_delta_v(0.0).expect("valid stage");
        let expected = 300.0 * G0 * (10.0_f64).ln();
        assert!((dv - expected).abs() < 1e-6, "dv {dv} vs {expected}");
    }

    #[test]
    fn thrust_interpolates_with_pressure() {
        let s = demo_stage();
        // At sea level -> thrust_sl, in vacuum -> thrust_vac.
        assert!((s.thrust(SEA_LEVEL_PRESSURE) - 180_000.0).abs() < 1e-6);
        assert!((s.thrust(0.0) - 200_000.0).abs() < 1e-6);
        // Halfway pressure -> halfway thrust.
        let mid = s.thrust(SEA_LEVEL_PRESSURE / 2.0);
        assert!((mid - 190_000.0).abs() < 1.0);
    }

    #[test]
    fn mass_flow_and_burn_time_consistent() {
        let s = demo_stage();
        let mdot = s.mass_flow().expect("valid stage");
        assert!((mdot - 200_000.0 / (300.0 * G0)).abs() < 1e-9);
        assert!((s.burn_time().expect("valid stage") - 9_000.0 / mdot).abs() < 1e-6);
    }

    #[test]
    fn burn_time_rejects_non_physical_stage_instead_of_nan() {
        // A hand-built stage that bypassed `from_engine` validation with
        // isp_vac = 0 and thrust_vac = 0 used to yield mass_flow = NaN and
        // thus a NaN burn time, returned silently as a plain value.
        let mut s = demo_stage();
        s.isp_vac = 0.0;
        s.thrust_vac = 0.0;
        let r = s.burn_time();
        assert!(
            matches!(r, Err(AstroError::InvalidPropulsion { .. })),
            "zero isp/thrust must be rejected, got {r:?}"
        );
        // isp_vac = 0 alone -> mass_flow = Inf -> burn_time = 0 silently.
        let mut s2 = demo_stage();
        s2.isp_vac = 0.0;
        assert!(matches!(
            s2.burn_time(),
            Err(AstroError::InvalidPropulsion {
                field: "isp_vac",
                ..
            })
        ));
        // thrust_vac = 0 alone -> mass_flow = 0 -> burn_time = Inf silently.
        let mut s3 = demo_stage();
        s3.thrust_vac = 0.0;
        assert!(matches!(
            s3.burn_time(),
            Err(AstroError::InvalidPropulsion {
                field: "thrust_vac",
                ..
            })
        ));
    }

    #[test]
    fn mass_flow_rejects_non_physical_stage_instead_of_silent_nan_inf() {
        // The standalone public `mass_flow` is the sibling of `burn_time`:
        // a hand-built stage with isp_vac = 0 -> F/(0·g₀) = Inf; isp_vac =
        // NaN -> NaN. Both were returned silently as plain values.
        let mut s = demo_stage();
        s.isp_vac = 0.0;
        assert!(matches!(
            s.mass_flow(),
            Err(AstroError::InvalidPropulsion {
                field: "isp_vac",
                ..
            })
        ));
        let mut s2 = demo_stage();
        s2.isp_vac = f64::NAN;
        assert!(s2.mass_flow().is_err());
        let mut s3 = demo_stage();
        s3.thrust_vac = 0.0;
        assert!(matches!(
            s3.mass_flow(),
            Err(AstroError::InvalidPropulsion {
                field: "thrust_vac",
                ..
            })
        ));
    }

    #[test]
    fn ideal_delta_v_rejects_non_physical_stage_instead_of_silent_nan() {
        // isp_vac = NaN multiplies the ln term into NaN.
        let mut s = demo_stage();
        s.isp_vac = f64::NAN;
        assert!(s.ideal_delta_v(0.0).is_err());
        // isp_vac = 0 alone is non-physical and rejected (it would zero the
        // budget, masking a misconfigured stage).
        let mut s0 = demo_stage();
        s0.isp_vac = 0.0;
        assert!(matches!(
            s0.ideal_delta_v(0.0),
            Err(AstroError::InvalidPropulsion {
                field: "isp_vac",
                ..
            })
        ));
        // A burnout mass of zero (dry_mass = 0, no upper mass) -> m0/mf is
        // m0/0 = Inf -> ln(Inf) = Inf silently. Now an InvalidMass error.
        let mut s2 = demo_stage();
        s2.dry_mass = 0.0;
        assert!(matches!(
            s2.ideal_delta_v(0.0),
            Err(AstroError::InvalidMass {
                field: "burnout_mass",
                ..
            })
        ));
        // A non-finite upper_mass poisons both m0 and mf.
        assert!(demo_stage().ideal_delta_v(f64::INFINITY).is_err());

        // Vehicle::ideal_delta_v propagates a bad stage rather than
        // summing a NaN into the total.
        let v = Vehicle {
            stages: vec![{
                let mut bad = demo_stage();
                bad.isp_vac = f64::NAN;
                bad
            }],
            payload_mass: 500.0,
            reference_area: 10.0,
            drag: DragModel::generic_launch_vehicle(),
        };
        assert!(v.ideal_delta_v().is_err());
    }

    #[test]
    fn drag_curve_clamps_and_interpolates() {
        let d = DragModel::generic_launch_vehicle();
        assert!((d.cd(0.0) - 0.30).abs() < 1e-9); // clamp low
        assert!((d.cd(100.0) - 0.22).abs() < 1e-9); // clamp high
                                                    // Between Mach 1.0 (0.55) and 1.2 (0.60): midpoint ≈ 0.575.
        assert!((d.cd(1.1) - 0.575).abs() < 1e-6);
        // Transonic peak is the maximum.
        let peak = d.cd(1.2);
        assert!(peak >= d.cd(0.5) && peak >= d.cd(3.0));
    }

    #[test]
    fn two_stage_budget_sums_stages() {
        let v = Vehicle {
            stages: vec![demo_stage(), demo_stage()],
            payload_mass: 500.0,
            reference_area: 10.0,
            drag: DragModel::generic_launch_vehicle(),
        };
        let s = &v.stages[0];
        let dv0 = s.ideal_delta_v(v.upper_mass_above(0)).expect("valid stage");
        let dv1 = v.stages[1]
            .ideal_delta_v(v.upper_mass_above(1))
            .expect("valid stage");
        assert!((v.ideal_delta_v().expect("valid vehicle") - (dv0 + dv1)).abs() < 1e-6);
        // Upper stage above the booster = payload + full second stage.
        assert!((v.upper_mass_above(0) - (500.0 + 10_000.0)).abs() < 1e-9);
        assert!((v.upper_mass_above(1) - 500.0).abs() < 1e-9);
    }

    #[test]
    fn validation_rejects_bad_vehicles() {
        let mut v = Vehicle {
            stages: vec![demo_stage()],
            payload_mass: 100.0,
            reference_area: 5.0,
            drag: DragModel::generic_launch_vehicle(),
        };
        assert!(v.validate().is_ok());
        v.stages.clear();
        assert_eq!(v.validate(), Err(AstroError::NoStages));
    }
}
