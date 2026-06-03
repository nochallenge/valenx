//! Mesh module facade — bridges to `valenx-adapter-gmsh` /
//! `valenx-adapter-netgen` when wired through the full app. v1 ships
//! a parameter struct + a planning function that returns an opaque
//! handle name; the actual mesher call is performed by the adapter
//! crates when invoked from `valenx-app`.

use serde::{Deserialize, Serialize};

use crate::error::SalomeError;

/// Meshing parameters.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MeshParams {
    /// Target characteristic length.
    pub h: f64,
    /// Element order (1 = linear, 2 = quadratic).
    pub order: u8,
}

impl Default for MeshParams {
    fn default() -> Self {
        Self { h: 1.0, order: 1 }
    }
}

/// Mesh `solid_name` with `params`. Returns the *name* of the mesh
/// object that the bridge expects the caller to look up via the
/// underlying adapter — the bridge itself is a planner, not the
/// runtime.
pub fn mesh_solid(solid_name: &str, params: &MeshParams) -> Result<String, SalomeError> {
    if !params.h.is_finite() || params.h <= 0.0 {
        return Err(SalomeError::BadParameter {
            name: "h",
            reason: format!("must be > 0 (got {})", params.h),
        });
    }
    if !(1..=2).contains(&params.order) {
        return Err(SalomeError::BadParameter {
            name: "order",
            reason: format!("must be 1 or 2 (got {})", params.order),
        });
    }
    Ok(format!("{solid_name}.mesh"))
}
