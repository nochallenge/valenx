//! The Hill-type contractile model.
//!
//! A muscle's force is its peak isometric force `F₀ = σ · PCSA` (specific
//! tension × physiological cross-sectional area) scaled by three
//! dimensionless factors and projected onto the tendon through the pennation
//! angle:
//!
//! - [`active_force_length`] — peaks at the optimal fibre length
//!   (Gordon-Huxley-Julian 1966; here a Gaussian approximation),
//! - [`force_velocity`] — A. V. Hill's 1938 hyperbola for shortening, with a
//!   saturating rise above `F₀` for forced lengthening,
//! - [`passive_force_length`] — the connective-tissue spring resisting stretch
//!   beyond the optimal length.

use serde::{Deserialize, Serialize};

use crate::error::MuscleError;

/// Default maximum isometric stress of vertebrate skeletal muscle
/// ("specific tension"), in newtons per square centimetre. Reported values
/// cluster around 22–35 N/cm²; 25 is a standard mid-range figure.
pub const DEFAULT_SPECIFIC_TENSION_N_CM2: f64 = 25.0;

/// Hill force-velocity curvature constant `a / F₀` (dimensionless). Hill's
/// classic value is ≈ 0.25 for mixed muscle.
const HILL_A_REL: f64 = 0.25;

/// Eccentric (forced-lengthening) force plateau, as a multiple of `F₀`.
/// Lengthening muscle can bear roughly 1.4× its isometric maximum.
const ECCENTRIC_PLATEAU: f64 = 1.4;

/// Width of the Gaussian active force-length curve, in units of optimal
/// fibre length.
const FL_GAUSSIAN_WIDTH: f64 = 0.45;

/// Shape factor of the exponential passive force-length curve.
const PASSIVE_SHAPE_K: f64 = 5.0;

/// Active force-length factor `f_L ∈ [0, 1]` at a normalized fibre length
/// `L / L_opt`. Peaks at 1.0 at the optimal length and falls off either side.
pub fn active_force_length(normalized_length: f64) -> f64 {
    let z = (normalized_length - 1.0) / FL_GAUSSIAN_WIDTH;
    (-z * z).exp()
}

/// Force-velocity factor `f_V` at a normalized velocity `v / v_max`
/// (positive = shortening / concentric, negative = lengthening / eccentric).
/// Hill's hyperbola for shortening; a saturating rise toward the eccentric
/// plateau for lengthening.
pub fn force_velocity(normalized_velocity: f64) -> f64 {
    if normalized_velocity >= 0.0 {
        // Shortening: force falls to zero at the maximum shortening velocity.
        let v = normalized_velocity.min(1.0);
        (1.0 - v) / (1.0 + v / HILL_A_REL)
    } else {
        // Lengthening: force rises above F₀, saturating at the plateau.
        let s = -normalized_velocity;
        1.0 + (ECCENTRIC_PLATEAU - 1.0) * s / (s + HILL_A_REL)
    }
}

/// Passive (parallel-elastic) force-length factor, as a fraction of `F₀`.
/// Zero up to the optimal length, then an exponential rise; reaches ≈ `F₀` at
/// twice the optimal length.
pub fn passive_force_length(normalized_length: f64) -> f64 {
    if normalized_length <= 1.0 {
        0.0
    } else {
        ((PASSIVE_SHAPE_K * (normalized_length - 1.0)).exp() - 1.0) / (PASSIVE_SHAPE_K.exp() - 1.0)
    }
}

/// A skeletal muscle described by its force-relevant architecture.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Muscle {
    /// Physiological cross-sectional area (cm²). With specific tension this
    /// sets the peak isometric force `F₀`.
    pub pcsa_cm2: f64,
    /// Maximum isometric stress (N/cm²) — "muscle quality", force per area.
    pub specific_tension_n_cm2: f64,
    /// Optimal fibre length (m) — the length at which active force peaks.
    pub optimal_fiber_length_m: f64,
    /// Pennation angle at the optimal length (degrees, in `[0, 90)`).
    pub pennation_deg: f64,
    /// Maximum shortening velocity (optimal fibre lengths per second).
    pub max_shortening_velocity_lopt_per_s: f64,
    /// Voluntary activation capacity `[0, 1]` — the fraction of the muscle the
    /// nervous system can recruit. Untrained ≈ 0.90, well-trained ≈ 0.98.
    pub voluntary_activation: f64,
    /// Oxidative capacity `[0, 1]` — a fatigue-resistance proxy (mitochondrial
    /// / capillary density). Reported, not part of the peak-force calculation.
    pub oxidative_capacity: f64,
}

impl Muscle {
    /// Build a muscle from its core architecture, validating the parameters.
    /// `voluntary_activation` defaults to 0.90 and `oxidative_capacity` to
    /// 0.80; override them with [`Muscle::with_voluntary_activation`] /
    /// [`Muscle::with_oxidative_capacity`].
    pub fn new(
        pcsa_cm2: f64,
        specific_tension_n_cm2: f64,
        optimal_fiber_length_m: f64,
        pennation_deg: f64,
        max_shortening_velocity_lopt_per_s: f64,
    ) -> Result<Self, MuscleError> {
        check_positive("pcsa_cm2", pcsa_cm2)?;
        check_positive("specific_tension_n_cm2", specific_tension_n_cm2)?;
        check_positive("optimal_fiber_length_m", optimal_fiber_length_m)?;
        check_positive(
            "max_shortening_velocity_lopt_per_s",
            max_shortening_velocity_lopt_per_s,
        )?;
        if !(0.0..90.0).contains(&pennation_deg) {
            return Err(MuscleError::BadPennation(pennation_deg));
        }
        Ok(Self {
            pcsa_cm2,
            specific_tension_n_cm2,
            optimal_fiber_length_m,
            pennation_deg,
            max_shortening_velocity_lopt_per_s,
            voluntary_activation: 0.90,
            oxidative_capacity: 0.80,
        })
    }

    /// A textbook adult human quadriceps femoris (whole-muscle ballpark):
    /// large pennate muscle, ~140 cm² PCSA, 12° pennation, 95% voluntary
    /// activation.
    pub fn human_quadriceps() -> Self {
        Self::new(140.0, DEFAULT_SPECIFIC_TENSION_N_CM2, 0.10, 12.0, 8.0)
            .expect("valid human-quadriceps constants")
            .with_voluntary_activation(0.95)
    }

    /// Set the voluntary-activation fraction (clamped to `[0, 1]`).
    pub fn with_voluntary_activation(mut self, activation: f64) -> Self {
        self.voluntary_activation = activation.clamp(0.0, 1.0);
        self
    }

    /// Set the oxidative-capacity fraction (clamped to `[0, 1]`).
    pub fn with_oxidative_capacity(mut self, capacity: f64) -> Self {
        self.oxidative_capacity = capacity.clamp(0.0, 1.0);
        self
    }

    /// Peak isometric fibre force `F₀ = σ · PCSA` (newtons), before pennation.
    pub fn max_isometric_force_n(&self) -> f64 {
        self.specific_tension_n_cm2 * self.pcsa_cm2
    }

    /// Cosine of the pennation angle — the fraction of fibre force that acts
    /// along the tendon.
    pub fn pennation_factor(&self) -> f64 {
        self.pennation_deg.to_radians().cos()
    }

    /// Total force transmitted to the tendon (newtons) at a given activation,
    /// normalized length and normalized velocity. Combines the active
    /// (length × velocity) and passive contributions and projects through the
    /// pennation angle.
    pub fn tendon_force_n(
        &self,
        activation: f64,
        normalized_length: f64,
        normalized_velocity: f64,
    ) -> f64 {
        let a = activation.clamp(0.0, 1.0);
        let active =
            a * active_force_length(normalized_length) * force_velocity(normalized_velocity);
        let passive = passive_force_length(normalized_length);
        self.max_isometric_force_n() * (active + passive) * self.pennation_factor()
    }

    /// Maximum *usable* (voluntary) strength (newtons): the muscle's full
    /// available neural drive, at the optimal length and zero velocity.
    pub fn max_voluntary_strength_n(&self) -> f64 {
        self.tendon_force_n(self.voluntary_activation, 1.0, 0.0)
    }
}

fn check_positive(field: &'static str, value: f64) -> Result<(), MuscleError> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(MuscleError::NotPositive { field, value })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9
    }

    #[test]
    fn active_force_length_peaks_at_optimum() {
        assert!(close(active_force_length(1.0), 1.0));
        assert!(active_force_length(0.6) < 1.0);
        assert!(active_force_length(1.4) < 1.0);
        // Symmetric Gaussian: equal distance either side gives equal force.
        assert!(close(active_force_length(0.7), active_force_length(1.3)));
    }

    #[test]
    fn force_velocity_isometric_is_unity_and_falls_to_zero_at_vmax() {
        assert!(close(force_velocity(0.0), 1.0));
        assert!(close(force_velocity(1.0), 0.0));
        // Monotonic decreasing across the shortening range.
        assert!(force_velocity(0.25) > force_velocity(0.5));
        assert!(force_velocity(0.5) > force_velocity(0.75));
    }

    #[test]
    fn force_velocity_eccentric_exceeds_isometric_and_is_bounded() {
        let f = force_velocity(-0.5);
        assert!(f > 1.0);
        assert!(f < ECCENTRIC_PLATEAU + 1e-9);
        // Faster lengthening => closer to the plateau.
        assert!(force_velocity(-2.0) > force_velocity(-0.5));
    }

    #[test]
    fn passive_force_rises_only_beyond_optimum() {
        assert!(close(passive_force_length(1.0), 0.0));
        assert!(close(passive_force_length(0.8), 0.0));
        assert!(passive_force_length(1.3) > 0.0);
        // Reaches ≈ F₀ at twice the optimal length.
        assert!(close(passive_force_length(2.0), 1.0));
    }

    #[test]
    fn max_isometric_force_is_sigma_times_pcsa() {
        let m = Muscle::new(30.0, 25.0, 0.1, 0.0, 8.0).unwrap();
        assert!(close(m.max_isometric_force_n(), 750.0));
    }

    #[test]
    fn pennation_reduces_tendon_force() {
        let flat = Muscle::new(30.0, 25.0, 0.1, 0.0, 8.0).unwrap();
        let penn = Muscle::new(30.0, 25.0, 0.1, 25.0, 8.0).unwrap();
        assert!(penn.max_voluntary_strength_n() < flat.max_voluntary_strength_n());
    }

    #[test]
    fn human_quadriceps_is_in_a_physiological_range() {
        let q = Muscle::human_quadriceps();
        // Whole-quadriceps peak tendon force lands in the low thousands of N.
        let f = q.tendon_force_n(1.0, 1.0, 0.0);
        assert!((2500.0..5000.0).contains(&f), "got {f} N");
    }

    #[test]
    fn rejects_bad_parameters() {
        assert!(Muscle::new(0.0, 25.0, 0.1, 0.0, 8.0).is_err());
        assert!(Muscle::new(30.0, 25.0, 0.1, 95.0, 8.0).is_err());
        assert!(Muscle::new(30.0, f64::NAN, 0.1, 0.0, 8.0).is_err());
    }
}
