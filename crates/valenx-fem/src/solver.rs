//! FEM solver choice — which subprocess adapter to dispatch to.

use serde::{Deserialize, Serialize};

/// Pick the FEA backend.
#[derive(Copy, Clone, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
pub enum FemSolverChoice {
    /// CalculiX (CCX) — linear static / dynamic / thermal.
    #[default]
    CalculiX,
    /// Elmer — multiphysics FEA (Finnish CSC).
    Elmer,
    /// Code_Aster — EDF's large general-purpose FEA suite.
    CodeAster,
}

impl FemSolverChoice {
    /// Adapter ID string — matches the IDs registered in the desktop
    /// shell's adapter registry (see `valenx-app/src/lib.rs`).
    pub fn adapter_id(&self) -> &'static str {
        match self {
            FemSolverChoice::CalculiX => "calculix",
            FemSolverChoice::Elmer => "elmer",
            FemSolverChoice::CodeAster => "code-aster",
        }
    }

    /// User-facing label for the FEM panel's solver picker.
    pub fn label(&self) -> &'static str {
        match self {
            FemSolverChoice::CalculiX => "CalculiX (ccx)",
            FemSolverChoice::Elmer => "Elmer (ElmerSolver)",
            FemSolverChoice::CodeAster => "Code_Aster (aster)",
        }
    }

    /// Recommended input-file extension for this solver.
    pub fn input_extension(&self) -> &'static str {
        match self {
            FemSolverChoice::CalculiX => "inp",
            FemSolverChoice::Elmer => "sif",
            FemSolverChoice::CodeAster => "comm",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn adapter_ids_match_registry() {
        assert_eq!(FemSolverChoice::CalculiX.adapter_id(), "calculix");
        assert_eq!(FemSolverChoice::Elmer.adapter_id(), "elmer");
        assert_eq!(FemSolverChoice::CodeAster.adapter_id(), "code-aster");
    }

    #[test]
    fn input_extensions_are_solver_specific() {
        assert_eq!(FemSolverChoice::CalculiX.input_extension(), "inp");
        assert_eq!(FemSolverChoice::Elmer.input_extension(), "sif");
        assert_eq!(FemSolverChoice::CodeAster.input_extension(), "comm");
    }

    #[test]
    fn default_is_calculix() {
        assert_eq!(FemSolverChoice::default(), FemSolverChoice::CalculiX);
    }
}
