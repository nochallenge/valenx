//! Geom module facade — bridges to `valenx-cad`. v1 ships planner
//! functions that name an operation + its operands; the actual BRep
//! calls happen in `valenx-cad` when invoked from `valenx-app`.

use serde::{Deserialize, Serialize};

use crate::error::SalomeError;

/// Geometry operation tag.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum GeomOp {
    /// Box.
    Box {
        /// Length x.
        lx: f64,
        /// Length y.
        ly: f64,
        /// Length z.
        lz: f64,
    },
    /// Sphere.
    Sphere {
        /// Radius.
        r: f64,
    },
    /// Cylinder.
    Cylinder {
        /// Radius.
        r: f64,
        /// Height.
        h: f64,
    },
    /// Boolean union of two named operands.
    Union {
        /// Operand A name.
        a: String,
        /// Operand B name.
        b: String,
    },
    /// Boolean difference (a - b).
    Difference {
        /// Operand A name.
        a: String,
        /// Operand B name.
        b: String,
    },
    /// Boolean intersection.
    Intersection {
        /// Operand A name.
        a: String,
        /// Operand B name.
        b: String,
    },
}

/// Plan a geometry op — returns the result-object name on success.
pub fn plan(op: &GeomOp, name: &str) -> Result<String, SalomeError> {
    match op {
        GeomOp::Box { lx, ly, lz } => {
            if [lx, ly, lz].iter().any(|x| !x.is_finite() || **x <= 0.0) {
                return Err(SalomeError::BadParameter {
                    name: "box.lx/ly/lz",
                    reason: format!("must be > 0 (got {lx}, {ly}, {lz})"),
                });
            }
        }
        GeomOp::Sphere { r } => {
            if !r.is_finite() || *r <= 0.0 {
                return Err(SalomeError::BadParameter {
                    name: "sphere.r",
                    reason: format!("must be > 0 (got {r})"),
                });
            }
        }
        GeomOp::Cylinder { r, h } => {
            if [r, h].iter().any(|x| !x.is_finite() || **x <= 0.0) {
                return Err(SalomeError::BadParameter {
                    name: "cylinder.r/h",
                    reason: format!("must be > 0 (got r={r}, h={h})"),
                });
            }
        }
        GeomOp::Union { .. } | GeomOp::Difference { .. } | GeomOp::Intersection { .. } => {}
    }
    Ok(name.to_string())
}
