//! Canonical field type — the heart of `Results`.
//!
//! A `Field` is a named array of values living on a mesh region,
//! at a given location (node / cell / face / edge / region), with
//! units and a `TimeKey`. Rank is encoded in `FieldKind`.

use serde::{Deserialize, Serialize};

use crate::time::TimeKey;
use crate::units::Units;

/// Where on the mesh the values live.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Location {
    /// Stored at mesh vertices.
    OnNode,
    /// Stored per cell (one value per element).
    OnCell,
    /// Stored per face (shared face between cells or boundary face).
    OnFace,
    /// Stored per edge.
    OnEdge,
    /// A single value for the whole region (e.g. a reference
    /// temperature, a bulk density).
    RegionConstant,
}

/// Tensor rank of a field.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum FieldKind {
    /// Scalar value per sample point.
    Scalar,
    /// Vector of `dim` components (usually 2 or 3) per sample point.
    Vector { dim: u8 },
    /// Second-order tensor of `rows × cols` per sample point.
    Tensor { rows: u8, cols: u8 },
}

impl FieldKind {
    /// Number of f64 values per sample point.
    pub fn components(self) -> usize {
        match self {
            FieldKind::Scalar => 1,
            FieldKind::Vector { dim } => dim as usize,
            FieldKind::Tensor { rows, cols } => rows as usize * cols as usize,
        }
    }
}

/// Reference to a mesh region this field is defined on.
///
/// Regions come from `valenx-mesh`; we hold only a stable string ID
/// here so `valenx-fields` doesn't depend on the mesh crate.
#[derive(Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct RegionRef(pub String);

/// A named field of values over a mesh region.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Field {
    pub name: String,
    pub kind: FieldKind,
    pub location: Location,
    pub region: RegionRef,
    pub units: Units,
    pub time: TimeKey,
    pub data: Vec<f64>,
    /// Cached min/max across all components.
    pub range: Option<(f64, f64)>,
}

impl Field {
    /// Total sample-point count implied by `data.len()` and `kind`.
    pub fn samples(&self) -> usize {
        let c = self.kind.components();
        if c == 0 {
            0
        } else {
            self.data.len() / c
        }
    }

    /// Project the field onto a per-sample magnitude scalar.
    ///
    /// - `Scalar` returns a clone (round-trips through the same
    ///   "render-as-scalar" pipeline).
    /// - `Vector { dim }` returns `||v||_2` per sample as a new
    ///   scalar field; `name` becomes `"<original>_mag"`, `range`
    ///   is recomputed from the magnitudes.
    /// - `Tensor { rows, cols }` returns the Frobenius norm
    ///   `sqrt(sum_ij T_ij^2)` per sample as a scalar. Useful for
    ///   stress fields (3x3 Cauchy stress → equivalent-magnitude
    ///   surface plot). Other tensor norms (von Mises, principal
    ///   eigenvalue) are physics-specific and live in their own
    ///   helpers.
    pub fn magnitude_field(&self) -> Option<Field> {
        let components = self.kind.components();
        if components == 0 {
            return None;
        }
        match self.kind {
            FieldKind::Scalar => Some(self.clone()),
            FieldKind::Vector { .. } | FieldKind::Tensor { .. } => {
                let n = self.samples();
                if n * components != self.data.len() {
                    return None;
                }
                let mut data = Vec::with_capacity(n);
                for i in 0..n {
                    let mut sum_sq = 0.0_f64;
                    for c in 0..components {
                        let v = self.data[i * components + c];
                        sum_sq += v * v;
                    }
                    data.push(sum_sq.sqrt());
                }
                let mut out = Field {
                    name: format!("{}_mag", self.name),
                    kind: FieldKind::Scalar,
                    location: self.location,
                    region: self.region.clone(),
                    units: self.units,
                    time: self.time,
                    data,
                    range: None,
                };
                out.recompute_range();
                Some(out)
            }
        }
    }

    /// Recompute `range` from the current `data`. Safe to call on
    /// empty fields (leaves `range = None`).
    pub fn recompute_range(&mut self) {
        if self.data.is_empty() {
            self.range = None;
            return;
        }
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for &v in &self.data {
            if v.is_finite() {
                if v < lo {
                    lo = v;
                }
                if v > hi {
                    hi = v;
                }
            }
        }
        if lo.is_finite() && hi.is_finite() {
            self.range = Some((lo, hi));
        } else {
            self.range = None;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::PASCAL;

    fn sample_scalar() -> Field {
        Field {
            name: "pressure".into(),
            kind: FieldKind::Scalar,
            location: Location::OnCell,
            region: RegionRef("fluid".into()),
            units: PASCAL,
            time: TimeKey::Steady,
            data: vec![1.0, 2.0, 3.0, 4.0],
            range: None,
        }
    }

    #[test]
    fn components_per_kind() {
        assert_eq!(FieldKind::Scalar.components(), 1);
        assert_eq!(FieldKind::Vector { dim: 3 }.components(), 3);
        assert_eq!(FieldKind::Tensor { rows: 3, cols: 3 }.components(), 9);
    }

    #[test]
    fn samples_for_scalar() {
        let f = sample_scalar();
        assert_eq!(f.samples(), 4);
    }

    #[test]
    fn samples_for_vector() {
        let f = Field {
            name: "velocity".into(),
            kind: FieldKind::Vector { dim: 3 },
            location: Location::OnCell,
            region: RegionRef("fluid".into()),
            units: crate::units::METER_PER_SECOND,
            time: TimeKey::Steady,
            data: vec![1.0, 0.0, 0.0, 2.0, 0.0, 0.0],
            range: None,
        };
        assert_eq!(f.samples(), 2);
    }

    #[test]
    fn range_from_data() {
        let mut f = sample_scalar();
        f.recompute_range();
        assert_eq!(f.range, Some((1.0, 4.0)));
    }

    #[test]
    fn range_ignores_non_finite() {
        let mut f = sample_scalar();
        f.data = vec![1.0, f64::NAN, 3.0, f64::INFINITY];
        f.recompute_range();
        assert_eq!(f.range, Some((1.0, 3.0)));
    }

    #[test]
    fn magnitude_field_for_3d_vector_computes_per_sample_norm() {
        // Two velocity samples: (3,4,0) -> 5; (1,2,2) -> 3.
        let f = Field {
            name: "velocity".into(),
            kind: FieldKind::Vector { dim: 3 },
            location: Location::OnNode,
            region: RegionRef("fluid".into()),
            units: crate::units::METER_PER_SECOND,
            time: TimeKey::Steady,
            data: vec![3.0, 4.0, 0.0, 1.0, 2.0, 2.0],
            range: None,
        };
        let mag = f
            .magnitude_field()
            .expect("vector should produce magnitude");
        assert_eq!(mag.kind, FieldKind::Scalar);
        assert_eq!(mag.location, Location::OnNode);
        assert_eq!(mag.data.len(), 2);
        assert!((mag.data[0] - 5.0).abs() < 1e-12);
        assert!((mag.data[1] - 3.0).abs() < 1e-12);
        assert!(
            mag.name.contains("velocity"),
            "expected derived name, got {:?}",
            mag.name
        );
        // Range should be precomputed.
        assert_eq!(mag.range, Some((3.0, 5.0)));
    }

    #[test]
    fn magnitude_field_for_scalar_returns_clone() {
        // Calling magnitude on an already-scalar field returns a
        // copy of the field itself — convenient for code paths that
        // want "always a scalar with magnitude semantics".
        let f = sample_scalar();
        let mag = f.magnitude_field().expect("scalar should pass through");
        assert_eq!(mag.kind, FieldKind::Scalar);
        assert_eq!(mag.data, f.data);
        assert_eq!(mag.name, f.name);
    }

    #[test]
    fn magnitude_field_for_tensor_computes_frobenius_norm() {
        // 3x3 stress sample with components (1,2,2; 0,0,0; 0,0,0).
        // Frobenius = sqrt(1+4+4+0+...) = sqrt(9) = 3.
        let f = Field {
            name: "stress".into(),
            kind: FieldKind::Tensor { rows: 3, cols: 3 },
            location: Location::OnCell,
            region: RegionRef("solid".into()),
            units: PASCAL,
            time: TimeKey::Steady,
            data: vec![1.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            range: None,
        };
        let mag = f
            .magnitude_field()
            .expect("tensor should produce Frobenius norm");
        assert_eq!(mag.kind, FieldKind::Scalar);
        assert_eq!(mag.data.len(), 1);
        assert!((mag.data[0] - 3.0).abs() < 1e-12);
        assert!(mag.name.contains("stress"));
    }

    #[test]
    fn magnitude_field_for_voigt_stress_handles_six_component_vectors() {
        // Some adapters (CalculiX, Code_Aster) ship stress as a
        // Vector{dim:6} in Voigt notation rather than a 3x3 Tensor.
        // The magnitude is still |v|_2 — verify a known sample.
        // Voigt sample (sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_xz, sigma_yz)
        // = (3, 4, 0, 0, 0, 0) -> ||v|| = sqrt(9 + 16) = 5.
        let f = Field {
            name: "stress_voigt".into(),
            kind: FieldKind::Vector { dim: 6 },
            location: Location::OnCell,
            region: RegionRef("solid".into()),
            units: PASCAL,
            time: TimeKey::Steady,
            data: vec![3.0, 4.0, 0.0, 0.0, 0.0, 0.0],
            range: None,
        };
        let mag = f.magnitude_field().expect("voigt vector should be ok");
        assert_eq!(mag.data.len(), 1);
        assert!((mag.data[0] - 5.0).abs() < 1e-12);
    }

    #[test]
    fn magnitude_field_for_zero_dim_tensor_returns_none() {
        // Defensive: a 0x0 or 0xN tensor has no components — returning
        // None matches the Vector{dim:0} guard.
        let f = Field {
            name: "weird".into(),
            kind: FieldKind::Tensor { rows: 0, cols: 3 },
            location: Location::OnCell,
            region: RegionRef("x".into()),
            units: PASCAL,
            time: TimeKey::Steady,
            data: Vec::new(),
            range: None,
        };
        assert!(f.magnitude_field().is_none());
    }
}
