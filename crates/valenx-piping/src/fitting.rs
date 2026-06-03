//! Pipe fittings + valves.

use serde::{Deserialize, Serialize};

/// Standard pipe fittings.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum PipeFitting {
    /// 90° elbow.
    Elbow90,
    /// 45° elbow.
    Elbow45,
    /// 3-way tee.
    Tee,
    /// Eccentric / concentric reducer.
    Reducer,
    /// End cap.
    Cap,
    /// Pipe coupling.
    Coupling,
    /// Demountable union.
    Union,
}

impl PipeFitting {
    /// True when the fitting has 3 connection points (Tee) rather than 2.
    pub fn is_branching(self) -> bool {
        matches!(self, PipeFitting::Tee)
    }

    /// Short label for UI display.
    pub fn label(self) -> &'static str {
        match self {
            PipeFitting::Elbow90 => "Elbow 90°",
            PipeFitting::Elbow45 => "Elbow 45°",
            PipeFitting::Tee => "Tee",
            PipeFitting::Reducer => "Reducer",
            PipeFitting::Cap => "Cap",
            PipeFitting::Coupling => "Coupling",
            PipeFitting::Union => "Union",
        }
    }
}

/// Standard valves.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum Valve {
    /// Linear-motion isolation valve.
    Gate,
    /// Throttling valve.
    Globe,
    /// Quarter-turn ball.
    Ball,
    /// Flap / lift check.
    Check,
    /// Quarter-turn disk valve.
    Butterfly,
    /// Fine throttling (small needle).
    Needle,
    /// Elastomer diaphragm.
    Diaphragm,
}

impl Valve {
    /// Short label for UI display.
    pub fn label(self) -> &'static str {
        match self {
            Valve::Gate => "Gate",
            Valve::Globe => "Globe",
            Valve::Ball => "Ball",
            Valve::Check => "Check",
            Valve::Butterfly => "Butterfly",
            Valve::Needle => "Needle",
            Valve::Diaphragm => "Diaphragm",
        }
    }

    /// True when the valve is one-way (Check valve only).
    pub fn is_one_way(self) -> bool {
        matches!(self, Valve::Check)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tee_is_only_branching_fitting() {
        for f in [
            PipeFitting::Elbow90,
            PipeFitting::Elbow45,
            PipeFitting::Tee,
            PipeFitting::Reducer,
            PipeFitting::Cap,
            PipeFitting::Coupling,
            PipeFitting::Union,
        ] {
            assert_eq!(f.is_branching(), f == PipeFitting::Tee);
        }
    }

    #[test]
    fn check_valve_is_only_one_way() {
        for v in [
            Valve::Gate,
            Valve::Globe,
            Valve::Ball,
            Valve::Check,
            Valve::Butterfly,
            Valve::Needle,
            Valve::Diaphragm,
        ] {
            assert_eq!(v.is_one_way(), v == Valve::Check);
        }
    }

    #[test]
    fn labels_are_non_empty() {
        for v in [Valve::Gate, Valve::Ball] {
            assert!(!v.label().is_empty());
        }
    }
}
