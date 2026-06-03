//! Scoring functions and affinity grids.
//!
//! A docking score is the model's estimate of how favourably a pose
//! sits in the receptor. This crate ships two classical scoring
//! functions plus the grid machinery that makes them fast:
//!
//! - [`vina`] — a Vina-class empirical scoring function (gauss /
//!   repulsion / hydrophobic / H-bond inter-atomic terms plus the
//!   rotatable-bond entropy penalty). Built directly on
//!   [`valenx_dock`]'s frozen Vina weights.
//! - [`ad4`] — an AutoDock4-class force-field scoring function
//!   (12-6 van der Waals, 12-10 directional hydrogen bonding,
//!   screened-Coulomb electrostatics, and a volume-based desolvation
//!   term).
//! - [`gridmap`] — per-receptor-atom-type affinity-map precomputation
//!   over a [`crate::prep::GridBox`], plus trilinear interpolation for
//!   fast grid-based pose scoring.
//!
//! The common [`ScoringFunction`] enum lets the search and screening
//! layers stay generic over which function is in use.

pub mod ad4;
pub mod gridmap;
pub mod vina;

pub use ad4::Ad4Terms;
pub use gridmap::{AffinityMap, AffinityMapSet};
pub use vina::VinaTerms;

/// Which scoring function a docking / screening run uses.
///
/// The search and screening layers are generic over this so consensus
/// scoring can combine several functions without each caller branching
/// on the concrete type.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Default)]
pub enum ScoringFunction {
    /// The Vina-class empirical scoring function ([`vina`]).
    #[default]
    Vina,
    /// The AutoDock4-class force-field scoring function ([`ad4`]).
    AutoDock4,
}

impl ScoringFunction {
    /// A short stable label.
    pub fn label(self) -> &'static str {
        match self {
            ScoringFunction::Vina => "vina",
            ScoringFunction::AutoDock4 => "autodock4",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scoring_function_labels_are_stable() {
        assert_eq!(ScoringFunction::Vina.label(), "vina");
        assert_eq!(ScoringFunction::AutoDock4.label(), "autodock4");
        assert_eq!(ScoringFunction::default(), ScoringFunction::Vina);
    }
}
