//! Myostatin loss-of-function parameterization.
//!
//! Myostatin (GDF-8) is a negative regulator of skeletal-muscle growth.
//! Losing its function makes muscle *bigger* (PCSA up) but, per the animal
//! literature (e.g. Amthor et al., *PNAS* 2007), also *lower quality*:
//! reduced specific tension, a shift toward fast / fatigable fibres (higher
//! shortening velocity, lower oxidative capacity). The defaults encode that
//! direction; tune them to explore.

use serde::{Deserialize, Serialize};

use crate::hill::Muscle;

/// How loss of myostatin function reshapes a muscle. Each effect is the
/// fractional change *at full knockout*, scaled linearly by [`severity`].
///
/// [`severity`]: MyostatinKnockout::severity
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MyostatinKnockout {
    /// Severity `[0, 1]`: 0 = wild-type, 1 = full loss of function (null).
    pub severity: f64,
    /// Fractional PCSA (size) gain at full knockout — e.g. 0.8 = +80%.
    pub mass_gain_fraction: f64,
    /// Fractional specific-tension (quality) loss at full knockout — e.g.
    /// 0.30 = −30% force per unit area.
    pub specific_tension_loss_fraction: f64,
    /// Fractional max-shortening-velocity gain at full knockout (fast-fibre
    /// shift).
    pub vmax_gain_fraction: f64,
    /// Fractional oxidative-capacity loss at full knockout (fatigue).
    pub oxidative_loss_fraction: f64,
}

impl Default for MyostatinKnockout {
    fn default() -> Self {
        Self {
            severity: 1.0,
            mass_gain_fraction: 0.80,
            specific_tension_loss_fraction: 0.30,
            vmax_gain_fraction: 0.20,
            oxidative_loss_fraction: 0.40,
        }
    }
}

impl MyostatinKnockout {
    /// A full (homozygous-null) knockout with literature-default effect sizes.
    pub fn null() -> Self {
        Self::default()
    }

    /// A partial knockout at the given severity `[0, 1]`, otherwise using the
    /// default effect sizes.
    pub fn with_severity(severity: f64) -> Self {
        Self {
            severity: severity.clamp(0.0, 1.0),
            ..Self::default()
        }
    }

    /// Apply the (severity-scaled) effects to a baseline muscle, returning the
    /// modified muscle. Size goes up; quality (specific tension), fibre type
    /// (faster) and oxidative capacity move the way a real knockout does.
    pub fn apply_to(&self, base: &Muscle) -> Muscle {
        let s = self.severity.clamp(0.0, 1.0);
        let mut ko = *base;
        ko.pcsa_cm2 = base.pcsa_cm2 * (1.0 + self.mass_gain_fraction * s);
        ko.specific_tension_n_cm2 =
            base.specific_tension_n_cm2 * (1.0 - self.specific_tension_loss_fraction * s);
        ko.max_shortening_velocity_lopt_per_s =
            base.max_shortening_velocity_lopt_per_s * (1.0 + self.vmax_gain_fraction * s);
        ko.oxidative_capacity =
            (base.oxidative_capacity * (1.0 - self.oxidative_loss_fraction * s)).clamp(0.0, 1.0);
        ko
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn full_knockout_grows_size_but_cuts_quality() {
        let base = Muscle::human_quadriceps();
        let ko = MyostatinKnockout::null().apply_to(&base);
        assert!(ko.pcsa_cm2 > base.pcsa_cm2);
        assert!(ko.specific_tension_n_cm2 < base.specific_tension_n_cm2);
        assert!(ko.oxidative_capacity < base.oxidative_capacity);
        assert!(ko.max_shortening_velocity_lopt_per_s > base.max_shortening_velocity_lopt_per_s);
    }

    #[test]
    fn zero_severity_is_a_no_op() {
        let base = Muscle::human_quadriceps();
        let ko = MyostatinKnockout::with_severity(0.0).apply_to(&base);
        assert_eq!(ko.pcsa_cm2, base.pcsa_cm2);
        assert_eq!(ko.specific_tension_n_cm2, base.specific_tension_n_cm2);
    }
}
