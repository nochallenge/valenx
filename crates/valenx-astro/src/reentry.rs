//! Atmospheric **entry** trajectory and **aeroheating**.
//!
//! Integrates the planar entry equations of motion for a capsule or
//! booster diving back into the atmosphere — inverse-square gravity, drag
//! against the [`crate::atmosphere`] density, and an optional lift term
//! (a lift-to-drag ratio `L/D`) — and reports the engineering peaks that
//! size a heat shield and a structure: **peak stagnation heat flux**,
//! **peak deceleration**, the altitudes at which they occur, and the
//! **integrated heat load** (`∫q̇ dt`).
//!
//! The state is `(altitude h, speed V, flight-path angle γ)` with `γ`
//! measured **below the local horizontal** (so `γ > 0` is descending):
//!
//! ```text
//!   V̇ = −ρV²/(2β) − g·sin γ                     (drag + gravity)
//!   γ̇ = (−g/V + V/(R⊕+h))·cos γ − (L/D)·ρV/(2β)  (gravity, curvature, lift)
//!   ḣ = −V·sin γ
//! ```
//!
//! where `β = m/(C_d·A)` is the **ballistic coefficient** (kg/m²). The
//! integration is the shared fixed-step RK4, bounded by a step ceiling
//! and a minimum step.
//!
//! ## Stagnation heating — Sutton–Graves
//!
//! The convective stagnation-point heat flux uses the **Sutton–Graves**
//! correlation
//!
//! ```text
//!   q̇ = K · √(ρ / R_n) · V³,   K ≈ 1.7415e-4  (SI, Earth air)
//! ```
//!
//! with `R_n` the nose radius (m). The peak heat flux and its altitude,
//! and the time-integrated heat load, fall straight out of the
//! trajectory.
//!
//! ## Closed-form oracle — Allen–Eggers ballistic entry
//!
//! For a **ballistic** (`L/D = 0`) entry into an exponential atmosphere of
//! scale height `H`, the classic Allen–Eggers analysis gives the peak
//! deceleration in closed form — and, strikingly, **independent of the
//! ballistic coefficient**:
//!
//! ```text
//!   a_max = Vₑ²·sin(γₑ) / (2·e·H)        (e = Euler's number)
//! ```
//!
//! and the speed at peak Sutton–Graves heating
//!
//! ```text
//!   V(q̇_max) = Vₑ · e^(−1/6) ≈ 0.846·Vₑ.
//! ```
//!
//! The unit tests pin Sutton–Graves exactly, and validate that the fully
//! integrated trajectory (real US Standard atmosphere, gravity and
//! curvature included) reproduces both: a steep ballistic entry's peak-g
//! matches `a_max` to within a few percent, and the peak-heating speed
//! matches `0.846·Vₑ` to within a couple of percent (the small offset is
//! the gravity term the idealisation drops).

use serde::{Deserialize, Serialize};

use crate::atmosphere;
use crate::constants::{MU_EARTH, R_EARTH, G0};
use crate::error::{AstroError, Result};
use crate::sim::{check_step_count, MAX_SIM_STEPS};

/// Sutton–Graves stagnation-point heat-flux constant `K` (SI units, Earth
/// air): `q̇ = K·√(ρ/R_n)·V³`.
pub const SUTTON_GRAVES_K: f64 = 1.7415e-4;

/// Representative exponential-atmosphere **scale height** for Earth (m),
/// the value the Allen–Eggers oracle uses. The US Standard atmosphere is
/// close to exponential with this scale height through the entry-heating
/// band (~10–40 km).
pub const SCALE_HEIGHT: f64 = 7_200.0;

/// Smallest accepted entry-integration step (s).
pub const MIN_TIME_STEP: f64 = 1e-3;

/// Largest accepted entry-integration duration (s) — an entry that has
/// not terminated (surface or near-stop) within this is rejected.
pub const MAX_ENTRY_TIME: f64 = 100_000.0;

/// Sutton–Graves stagnation-point convective heat flux (W/m²):
/// `q̇ = K·√(ρ/R_n)·V³`.
///
/// # Errors
///
/// Returns [`AstroError::InvalidParameter`] if `density` is non-finite or
/// negative, `nose_radius` is non-finite or non-positive (it sits under a
/// square root in a denominator), or `speed` is non-finite.
pub fn stagnation_heat_flux(density: f64, speed: f64, nose_radius: f64) -> Result<f64> {
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
    Ok(SUTTON_GRAVES_K * (density / nose_radius).sqrt() * speed.powi(3))
}

/// Allen–Eggers closed-form **peak deceleration** (m/s²) of a ballistic
/// entry: `a_max = Vₑ²·sin(γₑ)/(2·e·H)`. Independent of the ballistic
/// coefficient.
///
/// `entry_speed` in m/s, `flight_path_angle` in **radians below the
/// horizontal** (`> 0`), `scale_height` in m.
///
/// # Errors
///
/// Returns [`AstroError::InvalidParameter`] for a non-finite or
/// non-positive entry speed / scale height, or a `flight_path_angle` not
/// in `(0, π/2]`.
pub fn allen_eggers_max_deceleration(
    entry_speed: f64,
    flight_path_angle: f64,
    scale_height: f64,
) -> Result<f64> {
    if !entry_speed.is_finite() || entry_speed <= 0.0 {
        return Err(AstroError::InvalidParameter(
            "entry speed must be finite and > 0",
        ));
    }
    if !scale_height.is_finite() || scale_height <= 0.0 {
        return Err(AstroError::InvalidParameter(
            "scale height must be finite and > 0",
        ));
    }
    if !flight_path_angle.is_finite()
        || flight_path_angle <= 0.0
        || flight_path_angle > std::f64::consts::FRAC_PI_2
    {
        return Err(AstroError::InvalidParameter(
            "flight-path angle must be in (0, pi/2]",
        ));
    }
    Ok(entry_speed * entry_speed * flight_path_angle.sin()
        / (2.0 * std::f64::consts::E * scale_height))
}

/// The entry initial conditions and vehicle aerodynamics.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct EntryConditions {
    /// Entry speed `Vₑ` (m/s).
    pub entry_speed: f64,
    /// Entry flight-path angle `γₑ` (rad, below the local horizontal,
    /// `> 0`).
    pub flight_path_angle: f64,
    /// Entry altitude (m) — the interface where integration begins.
    pub entry_altitude: f64,
    /// Ballistic coefficient `β = m/(C_d·A)` (kg/m²).
    pub ballistic_coefficient: f64,
    /// Nose radius `R_n` (m), for the Sutton–Graves heat flux.
    pub nose_radius: f64,
    /// Lift-to-drag ratio `L/D` (dimensionless); 0 for a ballistic entry.
    pub lift_to_drag: f64,
    /// Integration step (s).
    pub time_step: f64,
}

impl Default for EntryConditions {
    fn default() -> Self {
        Self {
            entry_speed: 7_000.0,
            flight_path_angle: std::f64::consts::FRAC_PI_2, // steep ballistic
            entry_altitude: 120_000.0,
            ballistic_coefficient: 4_000.0,
            nose_radius: 0.5,
            lift_to_drag: 0.0,
            time_step: 0.05,
        }
    }
}

/// The aeroheating / deceleration peaks of an integrated entry.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct EntryResult {
    /// Peak stagnation heat flux `q̇` (W/m²).
    pub peak_heat_flux: f64,
    /// Altitude at peak heat flux (m).
    pub peak_heating_altitude: f64,
    /// Speed at peak heat flux (m/s).
    pub peak_heating_speed: f64,
    /// Peak sensed (drag) deceleration, in Earth g's.
    pub peak_deceleration_g: f64,
    /// Altitude at peak deceleration (m).
    pub peak_deceleration_altitude: f64,
    /// Integrated heat load `∫q̇ dt` (J/m²).
    pub heat_load: f64,
    /// Final speed when integration stopped (m/s).
    pub final_speed: f64,
    /// Elapsed time (s).
    pub final_time: f64,
}

impl EntryConditions {
    /// Validate the entry conditions.
    ///
    /// # Errors
    ///
    /// [`AstroError::InvalidParameter`] for any non-physical aerodynamic
    /// or state field, [`AstroError::InvalidIntegration`] for a
    /// `time_step` below the [`MIN_TIME_STEP`] floor.
    pub fn validate(&self) -> Result<()> {
        if !self.entry_speed.is_finite() || self.entry_speed <= 0.0 {
            return Err(AstroError::InvalidParameter(
                "entry speed must be finite and > 0",
            ));
        }
        if !self.flight_path_angle.is_finite()
            || self.flight_path_angle <= 0.0
            || self.flight_path_angle > std::f64::consts::FRAC_PI_2
        {
            return Err(AstroError::InvalidParameter(
                "flight-path angle must be in (0, pi/2]",
            ));
        }
        if !self.entry_altitude.is_finite() || self.entry_altitude <= 0.0 {
            return Err(AstroError::InvalidParameter(
                "entry altitude must be finite and > 0",
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
        if !self.lift_to_drag.is_finite() {
            return Err(AstroError::InvalidParameter("lift-to-drag must be finite"));
        }
        if !self.time_step.is_finite() || self.time_step < MIN_TIME_STEP {
            return Err(AstroError::InvalidIntegration(
                "time_step below the minimum (1e-3 s)",
            ));
        }
        Ok(())
    }

    /// Integrate the entry and return its aeroheating / deceleration
    /// peaks.
    ///
    /// The run stops at the surface (`h ≤ 0`), when the vehicle has
    /// effectively stopped (speed below 1 m/s), on a **skip-out** (a
    /// lifting entry that climbs back above the entry interface re-enters
    /// vacuum on a ballistic arc — the heating phase is over), or at the
    /// time ceiling. The step count is bounded by [`MAX_SIM_STEPS`].
    ///
    /// # Errors
    ///
    /// The validation errors of [`EntryConditions::validate`], or
    /// [`AstroError::StepBudgetExceeded`] if the entry has not terminated
    /// within the step / time ceiling.
    pub fn simulate(&self) -> Result<EntryResult> {
        self.validate()?;

        let dt = self.time_step;
        let beta = self.ballistic_coefficient;
        let ld = self.lift_to_drag;
        let rn = self.nose_radius;

        // Bound the loop: a hard time ceiling and the shared step cap.
        let max_steps_by_time = (MAX_ENTRY_TIME / dt).ceil();
        let max_steps = if max_steps_by_time > MAX_SIM_STEPS as f64 {
            MAX_SIM_STEPS
        } else {
            max_steps_by_time as u64
        };
        check_step_count(max_steps)?;

        // State: altitude h (m), speed V (m/s), flight-path angle γ (rad).
        let mut h = self.entry_altitude;
        let mut v = self.entry_speed;
        let mut gamma = self.flight_path_angle;
        let mut t = 0.0f64;
        let mut steps = 0u64;

        let mut peak_q = 0.0f64;
        let mut peak_q_alt = h;
        let mut peak_q_speed = v;
        let mut peak_decel_g = 0.0f64;
        let mut peak_decel_alt = h;
        let mut heat_load = 0.0f64;
        // Track the deepest penetration so a lifting **skip-out** (climbing
        // back above the entry interface after dipping in) terminates the
        // run — past it the vehicle is in vacuum on a ballistic arc and
        // the aeroheating phase is over.
        let entry_altitude = h;
        let mut min_altitude = h;

        // dh/dt, dV/dt, dγ/dt as a function of (h, V, γ).
        let deriv = |h: f64, v: f64, g: f64| -> (f64, f64, f64) {
            let rho = atmosphere::sample(h).density;
            let g_local = MU_EARTH / (R_EARTH + h).powi(2);
            let drag_decel = rho * v * v / (2.0 * beta);
            let dv = -drag_decel - g_local * g.sin();
            // Guard the 1/V term: V is always > 1 m/s while the loop runs.
            // Lift acts to *reduce* the (below-horizontal) flight-path
            // angle — it pulls the nose up and flattens the dive — so the
            // L/D term enters with a negative sign.
            let dgamma = (-g_local / v + v / (R_EARTH + h)) * g.cos() - ld * drag_decel / v;
            let dh = -v * g.sin();
            (dh, dv, dgamma)
        };

        while h > 0.0 && v > 1.0 {
            if steps >= max_steps {
                return Err(AstroError::StepBudgetExceeded(max_steps));
            }
            if t >= MAX_ENTRY_TIME {
                return Err(AstroError::StepBudgetExceeded(max_steps));
            }
            steps += 1;

            min_altitude = min_altitude.min(h);
            // Skip-out: the vehicle has climbed back above the interface
            // having dipped well below it. The atmospheric pass is done.
            if h > entry_altitude && min_altitude < entry_altitude - 1_000.0 {
                break;
            }

            let rho = atmosphere::sample(h).density;
            // Sutton–Graves heat flux and sensed (drag) deceleration.
            let q_dot = SUTTON_GRAVES_K * (rho / rn).sqrt() * v.powi(3);
            if q_dot > peak_q {
                peak_q = q_dot;
                peak_q_alt = h;
                peak_q_speed = v;
            }
            let decel_g = (rho * v * v / (2.0 * beta)) / G0;
            if decel_g > peak_decel_g {
                peak_decel_g = decel_g;
                peak_decel_alt = h;
            }
            heat_load += q_dot * dt;

            // RK4 step over (h, V, γ).
            let k1 = deriv(h, v, gamma);
            let k2 = deriv(
                h + 0.5 * dt * k1.0,
                v + 0.5 * dt * k1.1,
                gamma + 0.5 * dt * k1.2,
            );
            let k3 = deriv(
                h + 0.5 * dt * k2.0,
                v + 0.5 * dt * k2.1,
                gamma + 0.5 * dt * k2.2,
            );
            let k4 = deriv(h + dt * k3.0, v + dt * k3.1, gamma + dt * k3.2);
            h += dt / 6.0 * (k1.0 + 2.0 * k2.0 + 2.0 * k3.0 + k4.0);
            v += dt / 6.0 * (k1.1 + 2.0 * k2.1 + 2.0 * k3.1 + k4.1);
            gamma += dt / 6.0 * (k1.2 + 2.0 * k2.2 + 2.0 * k3.2 + k4.2);
            t += dt;
        }

        Ok(EntryResult {
            peak_heat_flux: peak_q,
            peak_heating_altitude: peak_q_alt,
            peak_heating_speed: peak_q_speed,
            peak_deceleration_g: peak_decel_g,
            peak_deceleration_altitude: peak_decel_alt,
            heat_load,
            final_speed: v,
            final_time: t,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sutton_graves_is_exact() {
        // q̇ = K·√(ρ/R_n)·V³. ρ=1e-3, R_n=0.5, V=6000:
        // √(1e-3/0.5)=√(2e-3)=0.0447213595..., V³=2.16e11.
        // q̇ = 1.7415e-4 · 0.0447213595 · 2.16e11.
        let q = stagnation_heat_flux(1e-3, 6_000.0, 0.5).expect("ok");
        let expected = 1.7415e-4 * (1e-3_f64 / 0.5).sqrt() * 6_000.0_f64.powi(3);
        assert!((q - expected).abs() / expected < 1e-12, "q = {q}");
        // Heat flux scales as V³: doubling V multiplies q̇ by 8.
        let q2 = stagnation_heat_flux(1e-3, 12_000.0, 0.5).expect("ok");
        assert!((q2 / q - 8.0).abs() < 1e-9);
        // Smaller nose radius -> higher flux (√(1/R_n)).
        let q_blunt = stagnation_heat_flux(1e-3, 6_000.0, 2.0).expect("ok");
        assert!(q_blunt < q, "blunter nose should heat less");
    }

    #[test]
    fn allen_eggers_formula_is_exact_and_beta_independent() {
        // a_max = Vₑ²·sin(γₑ)/(2·e·H). Vₑ=7000, γ=90°, H=7200:
        // 7000²·1/(2·e·7200) = 4.9e7/(39150.8...) = 1251.6... m/s².
        let a = allen_eggers_max_deceleration(7_000.0, std::f64::consts::FRAC_PI_2, 7_200.0)
            .expect("ok");
        let expected =
            7_000.0_f64.powi(2) * 1.0 / (2.0 * std::f64::consts::E * 7_200.0);
        assert!((a - expected).abs() / expected < 1e-12, "a = {a}");
        // It does not depend on the ballistic coefficient (no β in it).
        // Steeper entry -> larger peak (sin γ).
        let shallow = allen_eggers_max_deceleration(7_000.0, 10.0_f64.to_radians(), 7_200.0)
            .expect("ok");
        assert!(shallow < a, "shallow {shallow} should be < steep {a}");
    }

    #[test]
    fn integrated_steep_entry_peak_g_matches_allen_eggers() {
        // ORACLE: a steep ballistic entry (γ=60°) integrated through the
        // real US Standard atmosphere with full gravity and curvature
        // reproduces the Allen–Eggers closed-form peak deceleration to
        // within a few percent.
        let gamma = 60.0_f64.to_radians();
        let ve = 7_000.0;
        let entry = EntryConditions {
            entry_speed: ve,
            flight_path_angle: gamma,
            entry_altitude: 120_000.0,
            ballistic_coefficient: 4_000.0,
            nose_radius: 0.5,
            lift_to_drag: 0.0,
            time_step: 0.02,
        };
        let r = entry.simulate().expect("valid entry");
        let a_max_ae = allen_eggers_max_deceleration(ve, gamma, SCALE_HEIGHT).expect("ok");
        let peak_g_ae = a_max_ae / G0;
        let ratio = r.peak_deceleration_g / peak_g_ae;
        assert!(
            (ratio - 1.0).abs() < 0.05,
            "integrated peak {:.2} g vs Allen-Eggers {:.2} g (ratio {ratio:.4})",
            r.peak_deceleration_g,
            peak_g_ae
        );
        // The peak deceleration occurs down in the lower atmosphere.
        assert!(
            r.peak_deceleration_altitude < 30_000.0,
            "peak-g altitude {} m",
            r.peak_deceleration_altitude
        );
    }

    #[test]
    fn integrated_peak_heating_speed_matches_allen_eggers() {
        // ORACLE: the speed at peak Sutton–Graves heating is ≈ Vₑ·e^(−1/6)
        // ≈ 0.846·Vₑ. The integrated value (with gravity) lands a couple
        // of percent below this; pin both the ratio and the band.
        let ve = 7_000.0;
        let entry = EntryConditions {
            entry_speed: ve,
            flight_path_angle: 60.0_f64.to_radians(),
            time_step: 0.02,
            ..EntryConditions::default()
        };
        let r = entry.simulate().expect("valid entry");
        let frac = r.peak_heating_speed / ve;
        let ae_frac = (-1.0_f64 / 6.0).exp(); // 0.8464817...
        assert!(
            (frac - ae_frac).abs() < 0.03,
            "peak-heating speed fraction {frac:.4} vs Allen-Eggers {ae_frac:.4}"
        );
        assert!((0.80..0.87).contains(&frac), "fraction {frac}");
        // Peak heating is high and positive; heat load is positive.
        assert!(r.peak_heat_flux > 0.0);
        assert!(r.heat_load > 0.0);
        // Peak heating occurs above peak deceleration (it happens earlier,
        // at higher speed and thinner air).
        assert!(
            r.peak_heating_altitude > r.peak_deceleration_altitude,
            "q-peak alt {} should be above g-peak alt {}",
            r.peak_heating_altitude,
            r.peak_deceleration_altitude
        );
    }

    #[test]
    fn lifting_entry_decelerates_less_than_ballistic() {
        // A lifting entry (L/D > 0) flattens the trajectory, spreading the
        // deceleration over a longer, gentler arc -> a lower peak-g than
        // the otherwise-identical ballistic entry.
        let base = EntryConditions {
            entry_speed: 7_500.0,
            flight_path_angle: 30.0_f64.to_radians(),
            time_step: 0.02,
            ..EntryConditions::default()
        };
        let ballistic = base.simulate().expect("ballistic");
        let lifting = EntryConditions {
            lift_to_drag: 0.5, // Apollo-class lifting entry
            ..base
        }
        .simulate()
        .expect("lifting");
        assert!(
            lifting.peak_deceleration_g < ballistic.peak_deceleration_g,
            "lifting {:.1} g should be < ballistic {:.1} g",
            lifting.peak_deceleration_g,
            ballistic.peak_deceleration_g
        );
    }

    #[test]
    fn rejects_non_physical_conditions() {
        assert!(stagnation_heat_flux(-1.0, 6_000.0, 0.5).is_err());
        assert!(stagnation_heat_flux(1e-3, 6_000.0, 0.0).is_err()); // R_n = 0
        assert!(allen_eggers_max_deceleration(0.0, 1.0, 7_200.0).is_err());
        assert!(allen_eggers_max_deceleration(7_000.0, 0.0, 7_200.0).is_err()); // γ = 0
        assert!(allen_eggers_max_deceleration(7_000.0, 2.0, 7_200.0).is_err()); // γ > 90°

        let base = EntryConditions::default();
        assert!(base.validate().is_ok());
        assert!(EntryConditions { entry_speed: 0.0, ..base }.validate().is_err());
        assert!(EntryConditions { ballistic_coefficient: 0.0, ..base }.validate().is_err());
        assert!(EntryConditions { nose_radius: -1.0, ..base }.validate().is_err());
        assert!(EntryConditions { flight_path_angle: 0.0, ..base }.validate().is_err());
        assert!(EntryConditions { time_step: 1e-9, ..base }.validate().is_err());
        assert!(EntryConditions { lift_to_drag: f64::NAN, ..base }.validate().is_err());
    }

    #[test]
    fn entry_types_round_trip_through_json() {
        // Pin the serde derives (front-end mission persistence relies on
        // them). EntryConditions holds exact input literals, so it must
        // round-trip byte-for-byte.
        let cond = EntryConditions::default();
        let s = serde_json::to_string(&cond).expect("serialize conditions");
        let back: EntryConditions = serde_json::from_str(&s).expect("deserialize conditions");
        assert_eq!(cond, back);

        // EntryResult carries computed f64s; JSON's shortest-decimal form
        // can land on an adjacent f64 on re-parse, so assert the re-parsed
        // fields match the originals within a tight relative tolerance.
        let result = cond.simulate().expect("valid entry");
        let s = serde_json::to_string(&result).expect("serialize result");
        let back: EntryResult = serde_json::from_str(&s).expect("deserialize result");
        assert!((back.peak_heat_flux - result.peak_heat_flux).abs() <= result.peak_heat_flux.abs() * 1e-12);
        assert!((back.heat_load - result.heat_load).abs() <= result.heat_load.abs() * 1e-12);
        assert!((back.final_speed - result.final_speed).abs() <= result.final_speed.abs() * 1e-12);
        assert!((back.final_time - result.final_time).abs() <= result.final_time.abs() * 1e-12);
    }
}
