//! Joint strength and the size-vs-strength comparison.
//!
//! The headline result lives here: [`compare_wild_vs_knockout`] shows that a
//! myostatin knockout can add a great deal of muscle *size* while adding far
//! less *usable strength* — and, if quality drops far enough, a bigger muscle
//! can even be weaker. That is exactly why myostatin-blocking drugs raised
//! lean mass but disappointed on strength in human trials.

use serde::{Deserialize, Serialize};

use crate::hill::Muscle;
use crate::myostatin::MyostatinKnockout;

/// Joint torque (N·m) from a tendon force (N) acting at a moment arm (m).
pub fn joint_torque(tendon_force_n: f64, moment_arm_m: f64) -> f64 {
    tendon_force_n * moment_arm_m
}

/// A side-by-side comparison of a baseline muscle and its myostatin-knockout
/// variant — the "bigger ≠ stronger" result.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct StrengthComparison {
    /// Baseline physiological cross-sectional area (cm²).
    pub baseline_pcsa_cm2: f64,
    /// Knockout physiological cross-sectional area (cm²).
    pub knockout_pcsa_cm2: f64,
    /// Baseline maximum voluntary strength (N).
    pub baseline_strength_n: f64,
    /// Knockout maximum voluntary strength (N).
    pub knockout_strength_n: f64,
    /// Percentage change in size (PCSA).
    pub size_gain_pct: f64,
    /// Percentage change in usable strength.
    pub strength_gain_pct: f64,
}

impl StrengthComparison {
    /// How far size growth outruns strength growth (percentage points). A
    /// large positive number is the classic myostatin result: much bigger,
    /// only a little stronger.
    pub fn size_strength_gap_pct(&self) -> f64 {
        self.size_gain_pct - self.strength_gain_pct
    }
}

/// Compare a baseline muscle against the same muscle after a myostatin
/// knockout, reporting both the size change and the (usually much smaller)
/// strength change.
pub fn compare_wild_vs_knockout(base: &Muscle, knockout: &MyostatinKnockout) -> StrengthComparison {
    let ko = knockout.apply_to(base);
    let baseline_strength = base.max_voluntary_strength_n();
    let knockout_strength = ko.max_voluntary_strength_n();
    StrengthComparison {
        baseline_pcsa_cm2: base.pcsa_cm2,
        knockout_pcsa_cm2: ko.pcsa_cm2,
        baseline_strength_n: baseline_strength,
        knockout_strength_n: knockout_strength,
        size_gain_pct: pct_change(base.pcsa_cm2, ko.pcsa_cm2),
        strength_gain_pct: pct_change(baseline_strength, knockout_strength),
    }
}

fn pct_change(from: f64, to: f64) -> f64 {
    (to - from) / from * 100.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn joint_torque_is_force_times_arm() {
        assert!((joint_torque(3000.0, 0.045) - 135.0).abs() < 1e-9);
    }

    #[test]
    fn knockout_is_much_bigger_than_it_is_stronger() {
        let base = Muscle::human_quadriceps();
        let cmp = compare_wild_vs_knockout(&base, &MyostatinKnockout::null());
        // Size goes up a lot...
        assert!(cmp.size_gain_pct > 50.0, "size gain {}", cmp.size_gain_pct);
        // ...usable strength much less.
        assert!(cmp.strength_gain_pct < cmp.size_gain_pct);
        assert!(cmp.size_strength_gap_pct() > 20.0);
    }

    #[test]
    fn severe_quality_loss_can_make_a_bigger_muscle_weaker() {
        let base = Muscle::human_quadriceps();
        // A knockout that adds modest size but craters quality.
        let ko = MyostatinKnockout {
            severity: 1.0,
            mass_gain_fraction: 0.30,
            specific_tension_loss_fraction: 0.60,
            vmax_gain_fraction: 0.20,
            oxidative_loss_fraction: 0.50,
        };
        let cmp = compare_wild_vs_knockout(&base, &ko);
        assert!(cmp.size_gain_pct > 0.0);
        assert!(
            cmp.strength_gain_pct < 0.0,
            "strength gain {}",
            cmp.strength_gain_pct
        );
    }

    #[test]
    fn improving_quality_and_neural_drive_recovers_strength() {
        let base = Muscle::human_quadriceps();
        let ko = MyostatinKnockout::null().apply_to(&base);
        let ko_strength = ko.max_voluntary_strength_n();
        // Train the *whole chain*: restore specific tension and push neural
        // drive up. Strength recovers above the knockout's, exceeding baseline
        // because the extra size now sits on restored quality.
        let improved = Muscle {
            specific_tension_n_cm2: base.specific_tension_n_cm2,
            ..ko.with_voluntary_activation(0.99)
        };
        assert!(improved.max_voluntary_strength_n() > ko_strength);
        assert!(improved.max_voluntary_strength_n() > base.max_voluntary_strength_n());
    }
}
