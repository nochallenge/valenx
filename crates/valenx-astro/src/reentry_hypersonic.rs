//! Hypersonic / atmospheric-reentry analysis — the **Allen–Eggers**
//! closed-form anchor over an **exponential atmosphere**, with stagnation
//! heating correlations.
//!
//! This module is the *analytic-anchor* companion to [`crate::reentry`].
//! Where [`crate::reentry`] integrates a planar entry through the full **US
//! Standard Atmosphere 1976** (the general engineering case), this module
//! works in the idealised setting in which the classic Allen–Eggers (1958)
//! ballistic-entry theory is **exact** — a constant flight-path angle `γ`
//! and a strictly exponential density `ρ(h) = ρ₀·e^(−h/H)`. That lets the
//! numerically integrated trajectory here be cross-checked against the
//! closed form to machine-converging precision (the key benchmark), rather
//! than the few-percent agreement you get against the real atmosphere.
//!
//! ## The Allen–Eggers ballistic-entry results
//!
//! For an unpowered, non-lifting body diving at constant `γ` (below the
//! local horizontal, `γ > 0`) into an exponential atmosphere, with
//! ballistic coefficient `β = m/(C_d·A)` (kg/m²), the planar equation
//! `m dV/dt = −½ρV²C_dA` together with `dh/dt = −V·sin γ` integrates in
//! closed form to the **velocity–altitude profile**
//!
//! ```text
//!   V(h) = Vₑ · exp( −ρ(h)·H / (2·β·sin γ) ).
//! ```
//!
//! (The scale height `H` appears explicitly: it is required for
//! dimensional consistency once `β = m/(C_d·A)` carries units of kg/m².
//! Some texts fold `H` into a non-dimensional density and write
//! `exp(−ρ̃/(2β sin γ))`; this module keeps `H` explicit.)
//!
//! Differentiating the sensed drag deceleration `a = ρV²/(2β)` along this
//! profile and maximising gives the celebrated result that the **peak
//! deceleration is independent of the ballistic coefficient and the
//! mass**:
//!
//! ```text
//!   a_max = Vₑ²·sin γ / (2·e·H)            (e = Euler's number)
//! ```
//!
//! attained where the dimensionless group `u = ρH/(β sin γ) = 1`, i.e. at
//! the density `ρ* = β·sin γ / H` and hence the **peak-deceleration
//! altitude**
//!
//! ```text
//!   h* = H · ln( ρ₀·H / (β·sin γ) ).
//! ```
//!
//! Unlike `a_max`, `h*` *does* depend on `β`: a denser/heavier body (larger
//! `β`) penetrates to a **lower** altitude before peaking — a monotone
//! relationship the benchmarks pin.
//!
//! ## Stagnation-point heating
//!
//! Two engineering correlations for the convective stagnation-point heat
//! flux of a blunt body of nose radius `R_n`:
//!
//! - **Sutton–Graves** `q̇ = K·√(ρ/R_n)·V³`, `K ≈ 1.7415e-4` (SI, Earth
//!   air) — the workhorse optical-region correlation. (This is shared with
//!   [`crate::reentry::stagnation_heat_flux`]; re-exposed here for the
//!   self-contained exponential-atmosphere workflow.)
//! - **Fay–Riddell** equilibrium-catalytic flux, written in the common
//!   engineering form `q̇ = C·√(ρ/R_n)·V³` with `C ≈ 1.83e-4` (SI) for a
//!   fully-catalytic cold wall — the same `√(ρ/R_n)·V³` scaling, a slightly
//!   different constant.
//!
//! ## Honest scope
//!
//! This is **analytic / engineering-grade** reentry physics: a planar
//! point-mass body, a single exponential atmosphere layer, constant `γ`,
//! and optical-region stagnation-heating *correlations*. It is **not**
//! coupled CFD, not a real-gas / chemical-nonequilibrium flow solver, and
//! models **no ablation, no radiative heating, no shock-layer chemistry,
//! and no 3-D body shape** beyond the nose radius in the heating term. The
//! intent is **defensive survivability / heat-shield-sizing analysis** —
//! peak g-load, peak heat flux, integrated heat load — not weapons design.
//! For the general (real-atmosphere, lifting) case use [`crate::reentry`].

use serde::{Deserialize, Serialize};

use crate::error::{AstroError, Result};
use crate::reentry::{SCALE_HEIGHT, SUTTON_GRAVES_K};
use crate::sim::{check_step_count, MAX_SIM_STEPS};

/// Sea-level reference density `ρ₀` (kg/m³) of the exponential atmosphere,
/// matching the US Standard Atmosphere 1976 sea-level value used elsewhere
/// in the crate.
pub const RHO0_SEA_LEVEL: f64 = 1.225;

/// Fay–Riddell stagnation-point heat-flux constant `C` (SI, Earth air) for
/// a fully-catalytic cold wall, in the engineering form
/// `q̇ = C·√(ρ/R_n)·V³`. Slightly higher than the Sutton–Graves `K`
/// because the equilibrium-catalytic wall recombines dissociated species
/// and recovers their chemical enthalpy.
pub const FAY_RIDDELL_C: f64 = 1.83e-4;

/// A simple **exponential atmosphere**: `ρ(h) = ρ₀·e^(−h/H)`.
///
/// This is the density law under which the Allen–Eggers theory is exact.
/// It is deliberately separate from the seven-layer
/// [`crate::atmosphere`] US Standard model: here the *whole point* is an
/// analytically integrable density, so the integrated trajectory and the
/// closed form share an identical atmosphere and must converge.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ExponentialAtmosphere {
    /// Sea-level (h = 0) density `ρ₀` (kg/m³).
    pub rho0: f64,
    /// Scale height `H` (m).
    pub scale_height: f64,
}

impl Default for ExponentialAtmosphere {
    /// Earth: `ρ₀ = 1.225 kg/m³`, `H = 7200 m` (the [`SCALE_HEIGHT`] the
    /// Allen–Eggers oracle uses).
    fn default() -> Self {
        Self {
            rho0: RHO0_SEA_LEVEL,
            scale_height: SCALE_HEIGHT,
        }
    }
}

impl ExponentialAtmosphere {
    /// Build a validated exponential atmosphere.
    ///
    /// # Errors
    ///
    /// [`AstroError::InvalidParameter`] if `rho0` is non-finite or
    /// non-positive, or `scale_height` is non-finite or non-positive (it
    /// sits in a denominator and inside `exp(−h/H)`).
    pub fn new(rho0: f64, scale_height: f64) -> Result<Self> {
        if !rho0.is_finite() || rho0 <= 0.0 {
            return Err(AstroError::InvalidParameter("rho0 must be finite and > 0"));
        }
        if !scale_height.is_finite() || scale_height <= 0.0 {
            return Err(AstroError::InvalidParameter(
                "scale height must be finite and > 0",
            ));
        }
        Ok(Self { rho0, scale_height })
    }

    /// Density `ρ(h) = ρ₀·e^(−h/H)` (kg/m³) at geometric altitude `h` (m).
    ///
    /// Negative altitudes extrapolate the exponential (density rises),
    /// which is physically reasonable just below the datum; callers that
    /// must reject sub-surface states do so at the trajectory level.
    #[must_use]
    pub fn density(&self, altitude_m: f64) -> f64 {
        self.rho0 * (-altitude_m / self.scale_height).exp()
    }

    /// Invert the density law: the altitude (m) at which the density equals
    /// `rho`, i.e. `h = H·ln(ρ₀/ρ)`.
    ///
    /// # Errors
    ///
    /// [`AstroError::InvalidParameter`] if `rho` is non-finite or
    /// non-positive (the logarithm is undefined).
    pub fn altitude_of_density(&self, rho: f64) -> Result<f64> {
        if !rho.is_finite() || rho <= 0.0 {
            return Err(AstroError::InvalidParameter(
                "density must be finite and > 0",
            ));
        }
        Ok(self.scale_height * (self.rho0 / rho).ln())
    }
}

/// Sutton–Graves stagnation-point convective heat flux (W/m²):
/// `q̇ = K·√(ρ/R_n)·V³`, with `K =` [`SUTTON_GRAVES_K`].
///
/// Identical to [`crate::reentry::stagnation_heat_flux`]; provided here so
/// the exponential-atmosphere hypersonic workflow is self-contained.
///
/// # Errors
///
/// [`AstroError::InvalidParameter`] if `density` is non-finite or negative,
/// `nose_radius` is non-finite or non-positive, or `speed` is non-finite.
pub fn sutton_graves_heat_flux(density: f64, speed: f64, nose_radius: f64) -> Result<f64> {
    validate_heating_inputs(density, speed, nose_radius)?;
    Ok(SUTTON_GRAVES_K * (density / nose_radius).sqrt() * speed.powi(3))
}

/// Fay–Riddell equilibrium-catalytic stagnation-point heat flux (W/m²) in
/// the engineering form `q̇ = C·√(ρ/R_n)·V³`, with `C =`
/// [`FAY_RIDDELL_C`] (fully-catalytic cold wall, Earth air).
///
/// Shares the Sutton–Graves `√(ρ/R_n)·V³` scaling; the larger constant
/// reflects the chemical enthalpy recovered at a catalytic wall, so for the
/// same `(ρ, V, R_n)` this returns a flux a few percent above Sutton–Graves.
///
/// # Errors
///
/// [`AstroError::InvalidParameter`] under the same conditions as
/// [`sutton_graves_heat_flux`].
pub fn fay_riddell_heat_flux(density: f64, speed: f64, nose_radius: f64) -> Result<f64> {
    validate_heating_inputs(density, speed, nose_radius)?;
    Ok(FAY_RIDDELL_C * (density / nose_radius).sqrt() * speed.powi(3))
}

/// Shared fail-loud guard for the stagnation-heating correlations.
fn validate_heating_inputs(density: f64, speed: f64, nose_radius: f64) -> Result<()> {
    if !density.is_finite() || density < 0.0 {
        return Err(AstroError::InvalidParameter(
            "density must be finite and >= 0",
        ));
    }
    if !nose_radius.is_finite() || nose_radius <= 0.0 {
        return Err(AstroError::InvalidParameter(
            "nose_radius must be finite and > 0",
        ));
    }
    if !speed.is_finite() {
        return Err(AstroError::InvalidParameter("speed must be finite"));
    }
    Ok(())
}

/// The Allen–Eggers ballistic-entry problem: entry state, body ballistic
/// coefficient and nose radius, and the exponential atmosphere.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct BallisticEntry {
    /// Entry speed `Vₑ` (m/s) at the entry interface.
    pub entry_speed: f64,
    /// Constant flight-path angle `γ` (rad, below the local horizontal,
    /// `> 0` and `< π/2`). Allen–Eggers holds `γ` fixed along the dive.
    pub flight_path_angle: f64,
    /// Ballistic coefficient `β = m/(C_d·A)` (kg/m²).
    pub ballistic_coefficient: f64,
    /// Nose radius `R_n` (m), for the stagnation-heating correlations.
    pub nose_radius: f64,
    /// Entry altitude (m) where the dive begins.
    pub entry_altitude: f64,
    /// The exponential atmosphere the dive falls through.
    pub atmosphere: ExponentialAtmosphere,
}

impl Default for BallisticEntry {
    fn default() -> Self {
        Self {
            entry_speed: 7_000.0,
            flight_path_angle: 30.0_f64.to_radians(),
            ballistic_coefficient: 4_000.0,
            nose_radius: 0.5,
            entry_altitude: 120_000.0,
            atmosphere: ExponentialAtmosphere::default(),
        }
    }
}

impl BallisticEntry {
    /// Validate the entry configuration.
    ///
    /// # Errors
    ///
    /// [`AstroError::InvalidParameter`] for a non-finite or non-positive
    /// entry speed / ballistic coefficient / nose radius / entry altitude,
    /// or a `flight_path_angle` outside `(0, π/2)`. The atmosphere is
    /// validated via [`ExponentialAtmosphere::new`].
    pub fn validate(&self) -> Result<()> {
        if !self.entry_speed.is_finite() || self.entry_speed <= 0.0 {
            return Err(AstroError::InvalidParameter(
                "entry speed must be finite and > 0",
            ));
        }
        if !self.flight_path_angle.is_finite()
            || self.flight_path_angle <= 0.0
            || self.flight_path_angle >= std::f64::consts::FRAC_PI_2
        {
            return Err(AstroError::InvalidParameter(
                "flight-path angle must be in (0, pi/2)",
            ));
        }
        if !self.ballistic_coefficient.is_finite() || self.ballistic_coefficient <= 0.0 {
            return Err(AstroError::InvalidParameter(
                "ballistic coefficient must be finite and > 0",
            ));
        }
        if !self.nose_radius.is_finite() || self.nose_radius <= 0.0 {
            return Err(AstroError::InvalidParameter(
                "nose radius must be finite and > 0",
            ));
        }
        if !self.entry_altitude.is_finite() || self.entry_altitude <= 0.0 {
            return Err(AstroError::InvalidParameter(
                "entry altitude must be finite and > 0",
            ));
        }
        // Re-validate the atmosphere fields (a struct literal could carry a
        // bad rho0 / H that bypassed `ExponentialAtmosphere::new`).
        ExponentialAtmosphere::new(self.atmosphere.rho0, self.atmosphere.scale_height)?;
        Ok(())
    }

    /// Closed-form Allen–Eggers **velocity–altitude profile**
    /// `V(h) = Vₑ·exp(−ρ(h)·H/(2·β·sin γ))` (m/s) at altitude `h` (m).
    ///
    /// At the entry interface (`ρ → 0`) this returns `Vₑ`; deep in the
    /// atmosphere it decays toward zero.
    ///
    /// # Errors
    ///
    /// The validation errors of [`BallisticEntry::validate`].
    pub fn velocity_at_altitude(&self, altitude_m: f64) -> Result<f64> {
        self.validate()?;
        let rho = self.atmosphere.density(altitude_m);
        let exponent = -rho * self.atmosphere.scale_height
            / (2.0 * self.ballistic_coefficient * self.sin_gamma());
        Ok(self.entry_speed * exponent.exp())
    }

    /// Closed-form Allen–Eggers **peak deceleration**
    /// `a_max = Vₑ²·sin γ / (2·e·H)` (m/s²). Independent of the ballistic
    /// coefficient and the mass.
    ///
    /// # Errors
    ///
    /// The validation errors of [`BallisticEntry::validate`].
    pub fn analytic_peak_deceleration(&self) -> Result<f64> {
        self.validate()?;
        Ok(self.entry_speed * self.entry_speed * self.sin_gamma()
            / (2.0 * std::f64::consts::E * self.atmosphere.scale_height))
    }

    /// Closed-form **altitude of peak deceleration**
    /// `h* = H·ln(ρ₀·H/(β·sin γ))` (m), where the dimensionless group
    /// `ρH/(β sin γ)` reaches 1.
    ///
    /// A larger ballistic coefficient (denser/heavier body) gives a lower
    /// `h*` — it penetrates deeper before peaking.
    ///
    /// # Errors
    ///
    /// The validation errors of [`BallisticEntry::validate`].
    pub fn analytic_peak_deceleration_altitude(&self) -> Result<f64> {
        self.validate()?;
        // ρ* = β·sin γ / H (the density at which u = 1). Reuse the
        // atmosphere inverse so the value is consistent with `density`.
        let rho_star = self.ballistic_coefficient * self.sin_gamma() / self.atmosphere.scale_height;
        self.atmosphere.altitude_of_density(rho_star)
    }

    /// Closed-form **speed at peak deceleration**: along the Allen–Eggers
    /// profile the peak g occurs where `u = 1`, so
    /// `V(h*) = Vₑ·e^(−1/2) ≈ 0.6065·Vₑ`.
    ///
    /// # Errors
    ///
    /// The validation errors of [`BallisticEntry::validate`].
    pub fn analytic_speed_at_peak_deceleration(&self) -> Result<f64> {
        self.validate()?;
        Ok(self.entry_speed * (-0.5_f64).exp())
    }

    /// Sutton–Graves stagnation heat flux (W/m²) along the analytic profile
    /// at altitude `h`, using the closed-form `V(h)` and the exponential
    /// density.
    ///
    /// # Errors
    ///
    /// The validation errors of [`BallisticEntry::validate`].
    pub fn analytic_heat_flux_at_altitude(&self, altitude_m: f64) -> Result<f64> {
        let v = self.velocity_at_altitude(altitude_m)?;
        let rho = self.atmosphere.density(altitude_m);
        sutton_graves_heat_flux(rho, v, self.nose_radius)
    }

    /// Numerically integrate the ballistic dive (RK4 over the exponential
    /// atmosphere, constant-γ **Allen–Eggers idealisation**) and return the
    /// peaks. As the step shrinks the integrated peak deceleration
    /// converges to [`BallisticEntry::analytic_peak_deceleration`] to
    /// integrator precision — they solve the *same* equations of motion, so
    /// the only residual is RK4 truncation.
    ///
    /// **The dynamics are drag-only**, matching the Allen–Eggers analysis
    /// the closed forms come from: `m dV/dt = −½ρV²C_dA` with `dh/dt =
    /// −V·sin γ` and `γ` held fixed. Allen–Eggers neglects the gravity
    /// component along the path on the (steep, high-speed) entry where drag
    /// dominates; that is *why* the peak deceleration is independent of mass
    /// and ballistic coefficient. Including the `−g·sin γ` term shifts the
    /// integrated peak by a roughly constant few percent (gravity bleeds a
    /// little speed before the peak) and would *break* the clean
    /// convergence — so the **gravity-bearing, real-atmosphere, optionally
    /// lifting general case lives in [`crate::reentry`]**, and this method
    /// is deliberately the matching-idealisation numeric twin of the
    /// closed form.
    ///
    /// `time_step` is the RK4 step (s); the run stops at the surface, at
    /// near-stop (`V ≤ 1 m/s`), or at a generous time ceiling.
    ///
    /// # Errors
    ///
    /// The validation errors of [`BallisticEntry::validate`],
    /// [`AstroError::InvalidIntegration`] for a non-finite or non-positive
    /// `time_step`, or [`AstroError::StepBudgetExceeded`] if the dive does
    /// not terminate within the step ceiling.
    pub fn simulate(&self, time_step: f64) -> Result<BallisticEntryResult> {
        self.validate()?;
        if !time_step.is_finite() || time_step <= 0.0 {
            return Err(AstroError::InvalidIntegration(
                "time_step must be finite and > 0",
            ));
        }

        let dt = time_step;
        let beta = self.ballistic_coefficient;
        let rn = self.nose_radius;
        let atmos = self.atmosphere;
        // Constant flight-path angle (the Allen–Eggers idealisation): the
        // dive geometry is fixed, so only (h, V) evolve.
        let sin_g = self.sin_gamma();

        // Bound the loop. The dive is short; cap by a long ceiling and the
        // shared step budget.
        const MAX_TIME: f64 = 100_000.0;
        let max_steps_by_time = (MAX_TIME / dt).ceil();
        let max_steps = if max_steps_by_time > MAX_SIM_STEPS as f64 {
            MAX_SIM_STEPS
        } else {
            max_steps_by_time as u64
        };
        check_step_count(max_steps)?;

        let mut h = self.entry_altitude;
        let mut v = self.entry_speed;
        let mut t = 0.0f64;
        let mut steps = 0u64;

        let mut peak_decel = 0.0f64;
        let mut peak_decel_alt = h;
        let mut peak_q = 0.0f64;
        let mut peak_q_alt = h;
        let mut peak_q_speed = v;
        let mut heat_load = 0.0f64;

        // dh/dt, dV/dt as a function of (h, V) at fixed γ. Drag-only — the
        // Allen–Eggers idealisation the closed form solves, so the
        // integrated peak converges to `analytic_peak_deceleration` to
        // integrator precision. (The gravity-bearing, real-atmosphere,
        // optionally lifting general case is `crate::reentry`.)
        let deriv = |h: f64, v: f64| -> (f64, f64) {
            let rho = atmos.density(h);
            let drag = rho * v * v / (2.0 * beta);
            let dv = -drag;
            let dh = -v * sin_g;
            (dh, dv)
        };

        while h > 0.0 && v > 1.0 {
            if steps >= max_steps || t >= MAX_TIME {
                return Err(AstroError::StepBudgetExceeded(max_steps));
            }
            steps += 1;

            let rho = atmos.density(h);
            // Sensed (drag) deceleration — the Allen–Eggers quantity.
            let decel = rho * v * v / (2.0 * beta);
            if decel > peak_decel {
                peak_decel = decel;
                peak_decel_alt = h;
            }
            let q_dot = SUTTON_GRAVES_K * (rho / rn).sqrt() * v.powi(3);
            if q_dot > peak_q {
                peak_q = q_dot;
                peak_q_alt = h;
                peak_q_speed = v;
            }
            heat_load += q_dot * dt;

            // RK4 over (h, V).
            let k1 = deriv(h, v);
            let k2 = deriv(h + 0.5 * dt * k1.0, v + 0.5 * dt * k1.1);
            let k3 = deriv(h + 0.5 * dt * k2.0, v + 0.5 * dt * k2.1);
            let k4 = deriv(h + dt * k3.0, v + dt * k3.1);
            h += dt / 6.0 * (k1.0 + 2.0 * k2.0 + 2.0 * k3.0 + k4.0);
            v += dt / 6.0 * (k1.1 + 2.0 * k2.1 + 2.0 * k3.1 + k4.1);
            t += dt;
        }

        Ok(BallisticEntryResult {
            peak_deceleration: peak_decel,
            peak_deceleration_altitude: peak_decel_alt,
            peak_heat_flux: peak_q,
            peak_heating_altitude: peak_q_alt,
            peak_heating_speed: peak_q_speed,
            heat_load,
            final_speed: v,
            final_time: t,
        })
    }

    /// `sin γ`; `γ` is guaranteed in `(0, π/2)` by [`validate`], so this is
    /// strictly positive — the safe denominator for the velocity profile.
    fn sin_gamma(&self) -> f64 {
        self.flight_path_angle.sin()
    }
}

/// Peaks of a numerically integrated Allen–Eggers ballistic dive.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct BallisticEntryResult {
    /// Peak sensed (drag) deceleration (m/s²).
    pub peak_deceleration: f64,
    /// Altitude at peak deceleration (m).
    pub peak_deceleration_altitude: f64,
    /// Peak Sutton–Graves stagnation heat flux (W/m²).
    pub peak_heat_flux: f64,
    /// Altitude at peak heat flux (m).
    pub peak_heating_altitude: f64,
    /// Speed at peak heat flux (m/s).
    pub peak_heating_speed: f64,
    /// Integrated heat load `∫q̇ dt` (J/m²).
    pub heat_load: f64,
    /// Final speed when integration stopped (m/s).
    pub final_speed: f64,
    /// Elapsed time (s).
    pub final_time: f64,
}

impl BallisticEntryResult {
    /// Peak deceleration expressed in Earth g's (`/ g₀`).
    #[must_use]
    pub fn peak_deceleration_g(&self) -> f64 {
        self.peak_deceleration / crate::constants::G0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const E: f64 = std::f64::consts::E;

    fn rel_err(a: f64, b: f64) -> f64 {
        (a - b).abs() / b.abs()
    }

    // ------------------------------------------------------------------
    // Exponential atmosphere
    // ------------------------------------------------------------------

    #[test]
    fn exponential_atmosphere_is_the_textbook_law() {
        let atmos = ExponentialAtmosphere::default();
        // ρ(0) = ρ₀.
        assert!(rel_err(atmos.density(0.0), RHO0_SEA_LEVEL) < 1e-12);
        // One scale height up: ρ = ρ₀/e.
        let h = atmos.scale_height;
        assert!(rel_err(atmos.density(h), RHO0_SEA_LEVEL / E) < 1e-12);
        // Two scale heights: ρ = ρ₀/e².
        assert!(rel_err(atmos.density(2.0 * h), RHO0_SEA_LEVEL / (E * E)) < 1e-12);
        // Inverse round-trips.
        let rho = atmos.density(12_345.0);
        assert!(rel_err(atmos.altitude_of_density(rho).expect("ok"), 12_345.0) < 1e-9);
    }

    #[test]
    fn exponential_atmosphere_rejects_bad_config() {
        assert!(ExponentialAtmosphere::new(0.0, 7_200.0).is_err());
        assert!(ExponentialAtmosphere::new(-1.0, 7_200.0).is_err());
        assert!(ExponentialAtmosphere::new(1.225, 0.0).is_err());
        assert!(ExponentialAtmosphere::new(1.225, -1.0).is_err());
        assert!(ExponentialAtmosphere::new(f64::NAN, 7_200.0).is_err());
        let atmos = ExponentialAtmosphere::default();
        assert!(atmos.altitude_of_density(0.0).is_err());
        assert!(atmos.altitude_of_density(-1.0).is_err());
    }

    // ------------------------------------------------------------------
    // Stagnation heating — Sutton–Graves & Fay–Riddell
    // ------------------------------------------------------------------

    #[test]
    fn sutton_graves_matches_textbook_value() {
        // BENCHMARK (2): q̇ = K·√(ρ/R_n)·V³ at ρ=1e-3, R_n=0.5, V=6000.
        // √(1e-3/0.5)=√(2e-3)=0.04472135955, V³=2.16e11,
        // q̇ = 1.7415e-4 · 0.04472135955 · 2.16e11 = 1.681866...e6 W/m².
        let q = sutton_graves_heat_flux(1e-3, 6_000.0, 0.5).expect("ok");
        let expected = 1.7415e-4 * (1e-3_f64 / 0.5).sqrt() * 6_000.0_f64.powi(3);
        assert!(rel_err(q, expected) < 1e-12, "q = {q}");
        // Pin the absolute magnitude (~1.68 MW/m²) so a constant typo fails.
        assert!((1.60e6..1.76e6).contains(&q), "q = {q}");
        // V³ scaling: doubling V multiplies q̇ by 8.
        let q2 = sutton_graves_heat_flux(1e-3, 12_000.0, 0.5).expect("ok");
        assert!(rel_err(q2 / q, 8.0) < 1e-9);
    }

    #[test]
    fn fay_riddell_matches_form_and_exceeds_sutton_graves() {
        // q̇ = C·√(ρ/R_n)·V³ with C = 1.83e-4.
        let q = fay_riddell_heat_flux(1e-3, 6_000.0, 0.5).expect("ok");
        let expected = FAY_RIDDELL_C * (1e-3_f64 / 0.5).sqrt() * 6_000.0_f64.powi(3);
        assert!(rel_err(q, expected) < 1e-12, "q = {q}");
        // Catalytic wall recovers chemical enthalpy: slightly hotter than
        // Sutton–Graves, by exactly the constant ratio C/K.
        let q_sg = sutton_graves_heat_flux(1e-3, 6_000.0, 0.5).expect("ok");
        assert!(
            q > q_sg,
            "Fay-Riddell {q} should exceed Sutton-Graves {q_sg}"
        );
        assert!(rel_err(q / q_sg, FAY_RIDDELL_C / SUTTON_GRAVES_K) < 1e-12);
    }

    #[test]
    fn heating_correlations_reject_bad_inputs() {
        assert!(sutton_graves_heat_flux(-1.0, 6_000.0, 0.5).is_err()); // ρ<0
        assert!(sutton_graves_heat_flux(1e-3, 6_000.0, 0.0).is_err()); // R_n=0
        assert!(sutton_graves_heat_flux(1e-3, f64::NAN, 0.5).is_err()); // V NaN
        assert!(fay_riddell_heat_flux(-1.0, 6_000.0, 0.5).is_err());
        assert!(fay_riddell_heat_flux(1e-3, 6_000.0, -2.0).is_err());
    }

    // ------------------------------------------------------------------
    // Allen–Eggers closed forms
    // ------------------------------------------------------------------

    #[test]
    fn allen_eggers_amax_matches_closed_form() {
        // BENCHMARK (1a): a_max = Vₑ²·sin γ/(2·e·H), Vₑ=7000, γ=90°-ish
        // (use 90° via FRAC_PI_2 is rejected by validate, so use a steep
        // 80°), H=7200.
        let entry = BallisticEntry {
            entry_speed: 7_000.0,
            flight_path_angle: 80.0_f64.to_radians(),
            ..BallisticEntry::default()
        };
        let a = entry.analytic_peak_deceleration().expect("ok");
        let expected = 7_000.0_f64.powi(2) * 80.0_f64.to_radians().sin()
            / (2.0 * E * entry.atmosphere.scale_height);
        assert!(rel_err(a, expected) < 1e-12, "a = {a}");
    }

    #[test]
    fn allen_eggers_amax_is_beta_and_mass_invariant() {
        // BENCHMARK (1b): the elegant Allen–Eggers result — a_max does not
        // depend on the ballistic coefficient (and β = m/(C_d·A) carries
        // the mass, so it is mass-invariant too). Sweep β over 4 orders of
        // magnitude; a_max must not budge.
        let base = BallisticEntry {
            entry_speed: 7_500.0,
            flight_path_angle: 35.0_f64.to_radians(),
            ..BallisticEntry::default()
        };
        let a_ref = base.analytic_peak_deceleration().expect("ok");
        for beta in [10.0, 100.0, 1_000.0, 4_000.0, 50_000.0, 100_000.0] {
            let a = BallisticEntry {
                ballistic_coefficient: beta,
                ..base
            }
            .analytic_peak_deceleration()
            .expect("ok");
            assert!(
                rel_err(a, a_ref) < 1e-12,
                "a_max changed with β={beta}: {a} vs {a_ref}"
            );
        }
        // Steeper entry -> larger a_max (sin γ).
        let steep = BallisticEntry {
            flight_path_angle: 70.0_f64.to_radians(),
            ..base
        }
        .analytic_peak_deceleration()
        .expect("ok");
        assert!(steep > a_ref, "steeper {steep} should exceed {a_ref}");
    }

    #[test]
    fn allen_eggers_velocity_profile_matches_closed_form() {
        // BENCHMARK (companion to 1): V(h) = Vₑ·exp(−ρ(h)·H/(2β sin γ)).
        let entry = BallisticEntry::default();
        // At the entry interface the density is tiny -> V ≈ Vₑ.
        let v_top = entry
            .velocity_at_altitude(entry.entry_altitude)
            .expect("ok");
        assert!(rel_err(v_top, entry.entry_speed) < 1e-6, "v_top = {v_top}");
        // At a representative mid-altitude, check against the hand formula.
        let h = 20_000.0;
        let v = entry.velocity_at_altitude(h).expect("ok");
        let rho = entry.atmosphere.density(h);
        let expected = entry.entry_speed
            * (-rho * entry.atmosphere.scale_height
                / (2.0 * entry.ballistic_coefficient * entry.flight_path_angle.sin()))
            .exp();
        assert!(rel_err(v, expected) < 1e-12, "v = {v}");
        // Monotone: lower altitude (denser air) -> slower.
        let v_low = entry.velocity_at_altitude(5_000.0).expect("ok");
        assert!(v_low < v, "deeper should be slower: {v_low} vs {v}");
        // Speed at the analytic peak-decel altitude must equal Vₑ·e^(−1/2).
        let h_star = entry.analytic_peak_deceleration_altitude().expect("ok");
        let v_star = entry.velocity_at_altitude(h_star).expect("ok");
        let v_star_cf = entry.analytic_speed_at_peak_deceleration().expect("ok");
        assert!(rel_err(v_star, v_star_cf) < 1e-9, "{v_star} vs {v_star_cf}");
        assert!(rel_err(v_star, entry.entry_speed * (-0.5_f64).exp()) < 1e-9);
    }

    #[test]
    fn allen_eggers_peak_altitude_is_consistent() {
        // h* = H·ln(ρ₀H/(β sin γ)); equivalently the density there is
        // ρ* = β sin γ / H. Check both ways round-trip.
        let entry = BallisticEntry::default();
        let h_star = entry.analytic_peak_deceleration_altitude().expect("ok");
        let rho_star = entry.atmosphere.density(h_star);
        let expected_rho = entry.ballistic_coefficient * entry.flight_path_angle.sin()
            / entry.atmosphere.scale_height;
        assert!(rel_err(rho_star, expected_rho) < 1e-9, "ρ* mismatch");
    }

    #[test]
    fn higher_beta_penetrates_deeper_monotone() {
        // BENCHMARK (4): a larger ballistic coefficient peaks at a LOWER
        // altitude — monotone. Check both the closed form and the
        // integrated trajectory agree on the ordering.
        let base = BallisticEntry {
            entry_speed: 7_000.0,
            flight_path_angle: 30.0_f64.to_radians(),
            ..BallisticEntry::default()
        };
        let betas = [100.0, 500.0, 2_000.0, 8_000.0, 30_000.0];
        let mut prev_alt = f64::INFINITY;
        for beta in betas {
            let e = BallisticEntry {
                ballistic_coefficient: beta,
                ..base
            };
            let h_star = e.analytic_peak_deceleration_altitude().expect("ok");
            assert!(
                h_star < prev_alt,
                "β={beta} peak alt {h_star} should be below previous {prev_alt}"
            );
            prev_alt = h_star;
        }
        // The integrated solver shows the same ordering at its extremes.
        let low_beta = BallisticEntry {
            ballistic_coefficient: 100.0,
            ..base
        }
        .simulate(0.01)
        .expect("ok");
        let high_beta = BallisticEntry {
            ballistic_coefficient: 30_000.0,
            ..base
        }
        .simulate(0.01)
        .expect("ok");
        assert!(
            high_beta.peak_deceleration_altitude < low_beta.peak_deceleration_altitude,
            "integrated: high-β peak alt {} should be below low-β {}",
            high_beta.peak_deceleration_altitude,
            low_beta.peak_deceleration_altitude
        );
    }

    // ------------------------------------------------------------------
    // The key cross-check: integrated peak-g -> analytic a_max
    // ------------------------------------------------------------------

    #[test]
    fn integrated_peak_decel_converges_to_analytic_amax() {
        // BENCHMARK (3) — THE KEY CROSS-CHECK: the integrated trajectory's
        // peak deceleration converges to the Allen–Eggers closed form as
        // the step shrinks. The integrator solves the SAME drag-only
        // equations over the SAME exponential atmosphere the closed form
        // comes from, so the only residual is RK4 truncation and the
        // agreement is essentially exact at a fine step.
        for gamma_deg in [30.0_f64, 45.0, 60.0, 80.0] {
            let entry = BallisticEntry {
                entry_speed: 7_000.0,
                flight_path_angle: gamma_deg.to_radians(),
                ballistic_coefficient: 4_000.0,
                ..BallisticEntry::default()
            };
            let a_analytic = entry.analytic_peak_deceleration().expect("ok");

            let coarse = entry.simulate(1.0).expect("ok");
            let fine = entry.simulate(0.1).expect("ok");
            let finest = entry.simulate(0.005).expect("ok");

            let e_coarse = rel_err(coarse.peak_deceleration, a_analytic);
            let e_fine = rel_err(fine.peak_deceleration, a_analytic);
            let e_finest = rel_err(finest.peak_deceleration, a_analytic);

            // Refinement reduces the error (true convergence, not a fudged
            // tolerance masking a model gap).
            assert!(
                e_fine < e_coarse && e_finest <= e_fine,
                "γ={gamma_deg}: should converge: coarse {e_coarse:.2e}, \
                 fine {e_fine:.2e}, finest {e_finest:.2e}"
            );
            // The finest step lands on the closed form to <1e-4.
            assert!(
                e_finest < 1e-4,
                "γ={gamma_deg}: finest step within 1e-4 of analytic a_max: \
                 rel err {e_finest:.2e} (integrated {:.4} m/s², analytic {a_analytic:.4} m/s²)",
                finest.peak_deceleration
            );

            // The peak-decel altitude must also land on the closed form h*.
            let h_star = entry.analytic_peak_deceleration_altitude().expect("ok");
            assert!(
                rel_err(finest.peak_deceleration_altitude, h_star) < 0.02,
                "γ={gamma_deg}: integrated peak alt {} vs analytic h* {h_star}",
                finest.peak_deceleration_altitude
            );
        }
    }

    #[test]
    fn integrated_peak_heating_speed_near_allen_eggers_band() {
        // Sutton–Graves peak heating occurs higher and faster than peak
        // deceleration. Pin the qualitative ordering and positivity.
        let entry = BallisticEntry {
            entry_speed: 7_000.0,
            flight_path_angle: 45.0_f64.to_radians(),
            ..BallisticEntry::default()
        };
        let r = entry.simulate(0.02).expect("ok");
        assert!(r.peak_heat_flux > 0.0);
        assert!(r.heat_load > 0.0);
        assert!(
            r.peak_heating_altitude > r.peak_deceleration_altitude,
            "q-peak alt {} should be above g-peak alt {}",
            r.peak_heating_altitude,
            r.peak_deceleration_altitude
        );
        assert!(
            r.peak_heating_speed > entry.analytic_speed_at_peak_deceleration().expect("ok"),
            "heating happens earlier/faster than peak-g"
        );
    }

    // ------------------------------------------------------------------
    // Fail-loud configuration
    // ------------------------------------------------------------------

    #[test]
    fn rejects_non_physical_entry_config() {
        let base = BallisticEntry::default();
        assert!(base.validate().is_ok());

        assert!(BallisticEntry {
            entry_speed: 0.0,
            ..base
        }
        .validate()
        .is_err());
        assert!(BallisticEntry {
            entry_speed: -1.0,
            ..base
        }
        .validate()
        .is_err());
        // γ <= 0.
        assert!(BallisticEntry {
            flight_path_angle: 0.0,
            ..base
        }
        .validate()
        .is_err());
        // γ >= 90°.
        assert!(BallisticEntry {
            flight_path_angle: std::f64::consts::FRAC_PI_2,
            ..base
        }
        .validate()
        .is_err());
        assert!(BallisticEntry {
            flight_path_angle: 2.0,
            ..base
        }
        .validate()
        .is_err());
        // β <= 0.
        assert!(BallisticEntry {
            ballistic_coefficient: 0.0,
            ..base
        }
        .validate()
        .is_err());
        // R_n <= 0.
        assert!(BallisticEntry {
            nose_radius: 0.0,
            ..base
        }
        .validate()
        .is_err());
        assert!(BallisticEntry {
            nose_radius: -0.5,
            ..base
        }
        .validate()
        .is_err());
        // Negative entry altitude.
        assert!(BallisticEntry {
            entry_altitude: -100.0,
            ..base
        }
        .validate()
        .is_err());
        // Bad atmosphere embedded in the entry.
        assert!(BallisticEntry {
            atmosphere: ExponentialAtmosphere {
                rho0: -1.0,
                scale_height: 7_200.0,
            },
            ..base
        }
        .validate()
        .is_err());
        // Bad integration step.
        assert!(base.simulate(0.0).is_err());
        assert!(base.simulate(-0.1).is_err());
        assert!(base.simulate(f64::NAN).is_err());
    }

    #[test]
    fn types_round_trip_through_json() {
        // Pin the serde derives (front-end mission persistence relies on
        // them).
        let entry = BallisticEntry::default();
        let s = serde_json::to_string(&entry).expect("serialize");
        let back: BallisticEntry = serde_json::from_str(&s).expect("deserialize");
        assert_eq!(entry, back);

        let result = entry.simulate(0.05).expect("ok");
        let s = serde_json::to_string(&result).expect("serialize");
        let back: BallisticEntryResult = serde_json::from_str(&s).expect("deserialize");
        // Computed f64s: re-parsed shortest-decimal can land on an adjacent
        // f64, so compare with a tight relative tolerance.
        assert!(rel_err(back.peak_deceleration, result.peak_deceleration) < 1e-12);
        assert!(rel_err(back.peak_heat_flux, result.peak_heat_flux) < 1e-12);
        assert!(rel_err(back.heat_load, result.heat_load) < 1e-12);
    }
}
