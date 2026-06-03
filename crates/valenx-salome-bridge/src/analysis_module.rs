//! Analysis module facade — bridges to `valenx-fem`. v1 ships a
//! planning struct that names the mesh + the load case + the
//! solver; the actual solve runs in `valenx-fem` when invoked from
//! `valenx-app`.

use serde::{Deserialize, Serialize};

use crate::error::SalomeError;

/// Supported solver kinds.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum Solver {
    /// Linear elastic statics.
    LinearElastic,
    /// Steady heat conduction.
    SteadyHeat,
    /// Modal analysis (eigenfrequencies).
    Modal,
}

/// Analysis specification.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Analysis {
    /// Underlying mesh-object name.
    pub mesh_name: String,
    /// Load-case name (free-form).
    pub load_case: String,
    /// Solver kind.
    pub solver: Solver,
}

/// Plan an analysis — returns the result-object name.
pub fn plan(analysis: &Analysis, result_name: &str) -> Result<String, SalomeError> {
    if analysis.mesh_name.is_empty() {
        return Err(SalomeError::BadParameter {
            name: "mesh_name",
            reason: "must not be empty".into(),
        });
    }
    if analysis.load_case.is_empty() {
        return Err(SalomeError::BadParameter {
            name: "load_case",
            reason: "must not be empty".into(),
        });
    }
    Ok(result_name.to_string())
}
