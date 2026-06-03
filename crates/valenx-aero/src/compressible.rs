//! Subsonic compressibility correction — Prandtl-Glauert.
//!
//! The core solver is **incompressible** — it assumes the air density
//! is constant. That is an excellent approximation up to roughly
//! Mach 0.3; above it, compressibility starts to matter (the air
//! noticeably speeds up and the pressures change). A full compressible
//! solver is a different beast, but for the subsonic range there is a
//! cheap, classic fix: the **Prandtl-Glauert correction**, which
//! rescales the incompressible pressure coefficient by `1/β` with
//! `β = √(1 − M²)`.
//!
//! This module applies that correction to an incompressible result so
//! the coefficients are usable up into the high-subsonic range —
//! roughly Mach 0.3–0.7, below the drag-divergence Mach number.
//!
//! # Honest scope
//!
//! The Prandtl-Glauert rule is a linearised small-perturbation
//! correction. It is real and standard, but it is a *correction*, not
//! a compressible solver: it has no shock waves, it breaks down as
//! the local flow approaches Mach 1, and it does not capture the
//! transonic drag rise. The correction is clamped to a maximum Mach
//! number well below 1 and reports a warning when the regime is
//! pushed. Karman-Tsien and Laitone are documented refinements.

/// The speed of sound in air at a given absolute temperature.
///
/// `a = √(γ·R·T)` with `γ = 1.4` and `R = 287 J·kg⁻¹·K⁻¹`. At the
/// 288 K standard temperature this gives `a ≈ 340 m·s⁻¹`.
pub fn speed_of_sound(temperature_k: f64) -> f64 {
    let gamma = 1.4;
    let r = 287.0;
    (gamma * r * temperature_k.max(1.0)).sqrt()
}

/// The free-stream Mach number for a given speed and speed of sound.
pub fn mach_number(speed: f64, speed_of_sound: f64) -> f64 {
    if speed_of_sound <= 0.0 {
        return 0.0;
    }
    speed / speed_of_sound
}

/// The flow regime a Mach number falls in — used to qualify whether
/// the incompressible result + the Prandtl-Glauert correction is
/// trustworthy.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FlowRegime {
    /// `M < 0.3` — the incompressible solver is directly valid; no
    /// correction needed.
    Incompressible,
    /// `0.3 ≤ M < 0.7` — subsonic; the Prandtl-Glauert correction is
    /// appropriate.
    Subsonic,
    /// `0.7 ≤ M < 1.0` — high subsonic / transonic onset; the
    /// correction is increasingly unreliable.
    HighSubsonic,
    /// `M ≥ 1.0` — supersonic; this engine does not model it at all.
    Supersonic,
}

impl FlowRegime {
    /// Classify a Mach number.
    pub fn classify(mach: f64) -> FlowRegime {
        if mach < 0.3 {
            FlowRegime::Incompressible
        } else if mach < 0.7 {
            FlowRegime::Subsonic
        } else if mach < 1.0 {
            FlowRegime::HighSubsonic
        } else {
            FlowRegime::Supersonic
        }
    }

    /// A human-readable caveat for the regime.
    pub fn caveat(self) -> &'static str {
        match self {
            FlowRegime::Incompressible => {
                "incompressible — no compressibility correction needed"
            }
            FlowRegime::Subsonic => {
                "subsonic — Prandtl-Glauert correction applied, reliable"
            }
            FlowRegime::HighSubsonic => {
                "high subsonic — Prandtl-Glauert correction applied but increasingly \
                 unreliable near drag-divergence; not a transonic solver"
            }
            FlowRegime::Supersonic => {
                "supersonic — outside this engine's scope; result is not valid"
            }
        }
    }
}

/// The Prandtl-Glauert factor `β = √(1 − M²)`.
///
/// The compressible pressure coefficient is the incompressible one
/// divided by `β`. The Mach number is clamped to `0.95` so the factor
/// stays finite — pushing the rule to `M → 1` is unphysical.
pub fn prandtl_glauert_beta(mach: f64) -> f64 {
    let m = mach.clamp(0.0, 0.95);
    (1.0 - m * m).sqrt().max(1e-3)
}

/// The compressibility-corrected aerodynamic coefficients.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CompressibleCoefficients {
    /// The free-stream Mach number.
    pub mach: f64,
    /// The flow regime.
    pub regime: FlowRegime,
    /// The Prandtl-Glauert factor `β` used.
    pub beta: f64,
    /// The corrected drag coefficient.
    pub cd: f64,
    /// The corrected lift coefficient.
    pub cl: f64,
}

/// Apply the Prandtl-Glauert correction to incompressible drag / lift
/// coefficients.
///
/// The pressure-driven part of a coefficient scales by `1/β`; the
/// correction is applied to both `Cd` and `Cl` (a v1 simplification —
/// strictly only the pressure component scales, but for a slender body
/// the pressure part dominates the lift and a sizeable fraction of the
/// drag). The result carries the regime so a caller can decide whether
/// to trust it.
pub fn correct_coefficients(
    cd_incompressible: f64,
    cl_incompressible: f64,
    mach: f64,
) -> CompressibleCoefficients {
    let beta = prandtl_glauert_beta(mach);
    let regime = FlowRegime::classify(mach);
    // Below M = 0.3 the correction is a no-op (β ≈ 1 anyway, and the
    // regime is "incompressible").
    let factor = if regime == FlowRegime::Incompressible {
        1.0
    } else {
        1.0 / beta
    };
    CompressibleCoefficients {
        mach,
        regime,
        beta,
        cd: cd_incompressible * factor,
        cl: cl_incompressible * factor,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn speed_of_sound_at_standard_temperature() {
        // a ≈ 340 m/s at 288 K.
        let a = speed_of_sound(288.0);
        assert!((a - 340.0).abs() < 5.0, "speed of sound {a} should be ~340");
    }

    #[test]
    fn mach_number_is_speed_over_sound_speed() {
        let m = mach_number(170.0, 340.0);
        assert!((m - 0.5).abs() < 1e-9);
    }

    #[test]
    fn regime_classification() {
        assert_eq!(FlowRegime::classify(0.1), FlowRegime::Incompressible);
        assert_eq!(FlowRegime::classify(0.5), FlowRegime::Subsonic);
        assert_eq!(FlowRegime::classify(0.85), FlowRegime::HighSubsonic);
        assert_eq!(FlowRegime::classify(1.5), FlowRegime::Supersonic);
    }

    #[test]
    fn prandtl_glauert_beta_decreases_with_mach() {
        // β = √(1−M²): 1 at M=0, smaller as M grows.
        assert!((prandtl_glauert_beta(0.0) - 1.0).abs() < 1e-9);
        // M = 0.6 → β = √(1−0.36) = 0.8.
        assert!((prandtl_glauert_beta(0.6) - 0.8).abs() < 1e-9);
        assert!(prandtl_glauert_beta(0.9) < prandtl_glauert_beta(0.6));
        // The factor stays finite even at M → 1.
        assert!(prandtl_glauert_beta(1.0).is_finite());
        assert!(prandtl_glauert_beta(1.0) > 0.0);
    }

    #[test]
    fn correction_is_a_noop_in_the_incompressible_regime() {
        // Below M = 0.3 the coefficients are returned unchanged.
        let c = correct_coefficients(0.30, 0.50, 0.2);
        assert!((c.cd - 0.30).abs() < 1e-12);
        assert!((c.cl - 0.50).abs() < 1e-12);
        assert_eq!(c.regime, FlowRegime::Incompressible);
    }

    #[test]
    fn correction_amplifies_coefficients_in_the_subsonic_regime() {
        // At M = 0.6, β = 0.8, so the coefficients scale by 1/0.8 =
        // 1.25.
        let c = correct_coefficients(0.40, 0.20, 0.6);
        assert_eq!(c.regime, FlowRegime::Subsonic);
        assert!((c.cd - 0.40 / 0.8).abs() < 1e-9);
        assert!((c.cl - 0.20 / 0.8).abs() < 1e-9);
        // The corrected coefficient is larger than the incompressible
        // one — compressibility raises the suction peaks.
        assert!(c.cd > 0.40);
    }

    #[test]
    fn high_subsonic_regime_carries_a_caveat() {
        let c = correct_coefficients(0.3, 0.1, 0.85);
        assert_eq!(c.regime, FlowRegime::HighSubsonic);
        assert!(c.regime.caveat().contains("unreliable"));
    }
}
