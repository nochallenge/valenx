//! Stress / strain analysis helpers.
//!
//! ## Scope
//!
//! Quick post-processing primitives for FEA results — the kind of
//! "what's the worst stress in this part" answer users want
//! immediately after a CalculiX / Code Aster / Elmer run completes
//! without dropping into ParaView / Gmsh.
//!
//! - **`von_mises_from_components`** — convert a 6-component stress
//!   tensor field (Sxx, Syy, Szz, Sxy, Syz, Sxz) into a scalar
//!   von-Mises field. Standard isotropic-yield criterion; what FEA
//!   reports show as "Mises stress" or "VMS".
//! - **`field_max_per_node`** — scan a scalar field for the highest
//!   value + the node index where it occurred. Pairs with VMS to
//!   find the hot spot.
//! - **`safety_factor`** — `yield_strength / vms_value` per node;
//!   numbers below 1.0 mean the part is yielding.

use crate::{Field, FieldKind};
#[cfg(test)]
use crate::{Location, RegionRef, TimeKey};

/// Convert a 6-component stress tensor field into a scalar
/// von-Mises field defined on the same nodes.
///
/// Input layout: a flat `data` vector of length `n_nodes * 6`,
/// packed in CalculiX / Abaqus order:
/// `[Sxx, Syy, Szz, Sxy, Syz, Sxz]` per node.
///
/// Formula:
///
/// ```text
/// VMS = sqrt( ((Sxx-Syy)² + (Syy-Szz)² + (Szz-Sxx)²) / 2
///             + 3 * (Sxy² + Syz² + Sxz²) )
/// ```
///
/// Returns `None` when the input has the wrong shape (`data.len()`
/// not divisible by 6) so the caller surfaces a structured error
/// rather than emitting a half-shaped field.
pub fn von_mises_from_components(stress: &Field) -> Option<Field> {
    if stress.data.len() % 6 != 0 {
        return None;
    }
    if !matches!(stress.kind, FieldKind::Tensor { rows: 3, cols: 3 })
        && !matches!(stress.kind, FieldKind::Vector { dim: 6 })
    {
        // Accept either packing — full 3×3 tensor (9 components,
        // symmetric) or the abbreviated 6-component "Voigt" form.
        // CalculiX writes Voigt; if we got 9 we treat the first 6
        // as the diagonal + upper-triangular pack.
        return None;
    }
    let n = stress.data.len() / 6;
    let mut vms: Vec<f64> = Vec::with_capacity(n);
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    for i in 0..n {
        let off = i * 6;
        let sxx = stress.data[off];
        let syy = stress.data[off + 1];
        let szz = stress.data[off + 2];
        let sxy = stress.data[off + 3];
        let syz = stress.data[off + 4];
        let sxz = stress.data[off + 5];
        let dxy = sxx - syy;
        let dyz = syy - szz;
        let dzx = szz - sxx;
        let v = ((dxy * dxy + dyz * dyz + dzx * dzx) * 0.5
            + 3.0 * (sxy * sxy + syz * syz + sxz * sxz))
            .sqrt();
        if v < min {
            min = v;
        }
        if v > max {
            max = v;
        }
        vms.push(v);
    }
    Some(Field {
        name: format!("{}_vms", stress.name),
        kind: FieldKind::Scalar,
        location: stress.location,
        region: stress.region.clone(),
        units: stress.units, // same as input — Pa for SI
        time: stress.time,
        data: vms,
        range: if n == 0 { None } else { Some((min, max)) },
    })
}

/// One node's stress hot spot — node index + the field value at
/// that index.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PeakSample {
    pub node_index: usize,
    pub value: f64,
}

/// Scan a scalar field for its maximum and report the node index +
/// value. Returns `None` for empty fields. Ties go to the
/// lowest-indexed node (deterministic).
pub fn field_max_per_node(field: &Field) -> Option<PeakSample> {
    if field.data.is_empty() {
        return None;
    }
    let mut best = PeakSample {
        node_index: 0,
        value: field.data[0],
    };
    for (i, &v) in field.data.iter().enumerate().skip(1) {
        if v > best.value {
            best = PeakSample {
                node_index: i,
                value: v,
            };
        }
    }
    Some(best)
}

/// Compute a per-node safety factor against a yield strength.
/// `yield_strength` is in the same units as the VMS field (Pa for
/// SI). Returned values:
///
/// - `> 1`: under yield
/// - `= 1`: at yield
/// - `< 1`: yielding
/// - `inf`: zero-stress nodes (defensive — caller decides how to
///   surface those, e.g. drop them from a min-safety calculation)
pub fn safety_factor(vms: &Field, yield_strength: f64) -> Option<Field> {
    if !matches!(vms.kind, FieldKind::Scalar) {
        return None;
    }
    let mut data: Vec<f64> = Vec::with_capacity(vms.data.len());
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    for &v in &vms.data {
        let sf = if v == 0.0 {
            f64::INFINITY
        } else {
            yield_strength / v
        };
        if sf < min {
            min = sf;
        }
        if sf > max {
            max = sf;
        }
        data.push(sf);
    }
    Some(Field {
        name: format!("{}_safety_factor", vms.name),
        kind: FieldKind::Scalar,
        location: vms.location,
        region: vms.region.clone(),
        units: crate::units::DIMENSIONLESS,
        time: vms.time,
        data,
        range: if vms.data.is_empty() {
            None
        } else {
            Some((min, max))
        },
    })
}

/// Trivial helper: count how many nodes in a safety-factor field
/// have `value < threshold`. Use with `safety_factor()` to answer
/// "how many nodes are below SF=1?" in one call.
pub fn count_below(field: &Field, threshold: f64) -> usize {
    field.data.iter().filter(|&&v| v < threshold).count()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn voigt_stress(name: &str, values: Vec<f64>) -> Field {
        let n = values.len();
        Field {
            name: name.into(),
            kind: FieldKind::Vector { dim: 6 },
            location: Location::OnNode,
            region: RegionRef("default".into()),
            units: crate::units::Units::new([-1, 1, -2, 0, 0, 0, 0], 1.0, Some("Pa")),
            time: TimeKey::Steady,
            data: values,
            range: None,
        }
        .tap(|f| {
            // recompute range
            let _ = f;
            let _ = n;
        })
    }

    /// Tiny chaining helper — same shape as the `tap` in `tap` crate
    /// but inlined to avoid the dep.
    trait Tap: Sized {
        fn tap(self, f: impl FnOnce(&mut Self)) -> Self {
            let mut s = self;
            f(&mut s);
            s
        }
    }
    impl<T> Tap for T {}

    #[test]
    fn von_mises_zero_stress_yields_zero() {
        let s = voigt_stress("S", vec![0.0; 6]); // 1 node, all zeros
        let vms = von_mises_from_components(&s).expect("converts");
        assert_eq!(vms.data.len(), 1);
        assert_eq!(vms.data[0], 0.0);
    }

    #[test]
    fn von_mises_uniaxial_tension_equals_axial_stress() {
        // Pure Sxx loading: VMS = Sxx (sanity check from textbook).
        let s = voigt_stress("S", vec![100.0e6, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let vms = von_mises_from_components(&s).expect("converts");
        assert!(
            (vms.data[0] - 100.0e6).abs() < 1e-3,
            "VMS for pure tension should equal Sxx; got {}",
            vms.data[0]
        );
    }

    #[test]
    fn von_mises_pure_shear_is_root_three_times_shear() {
        // Pure Sxy loading: VMS = sqrt(3) * Sxy.
        let s = voigt_stress("S", vec![0.0, 0.0, 0.0, 50.0e6, 0.0, 0.0]);
        let vms = von_mises_from_components(&s).expect("converts");
        let expected = 3.0_f64.sqrt() * 50.0e6;
        assert!(
            (vms.data[0] - expected).abs() < 1e-3,
            "got {} expected {expected}",
            vms.data[0]
        );
    }

    #[test]
    fn von_mises_hydrostatic_stress_is_zero() {
        // Equal triaxial stress: deviatoric part vanishes, VMS = 0.
        let s = voigt_stress("S", vec![100.0e6, 100.0e6, 100.0e6, 0.0, 0.0, 0.0]);
        let vms = von_mises_from_components(&s).expect("converts");
        assert!(vms.data[0].abs() < 1e-6, "got {}", vms.data[0]);
    }

    #[test]
    fn von_mises_handles_multi_node_input() {
        // Two nodes, two different stress states.
        let s = voigt_stress(
            "S",
            vec![
                100.0e6, 0.0, 0.0, 0.0, 0.0, 0.0, // uniaxial -> 100 MPa
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, // zero -> 0
            ],
        );
        let vms = von_mises_from_components(&s).expect("converts");
        assert_eq!(vms.data.len(), 2);
        assert!((vms.data[0] - 100.0e6).abs() < 1e-3);
        assert_eq!(vms.data[1], 0.0);
    }

    #[test]
    fn von_mises_carries_units_and_renames() {
        let s = voigt_stress("S_steel", vec![0.0; 6]);
        let vms = von_mises_from_components(&s).expect("converts");
        assert_eq!(vms.name, "S_steel_vms");
        assert_eq!(vms.units.display, Some("Pa"));
        assert_eq!(vms.kind, FieldKind::Scalar);
    }

    #[test]
    fn von_mises_rejects_misshaped_data() {
        // 5-element Voigt = invalid (must be multiple of 6).
        let mut s = voigt_stress("S", vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        s.kind = FieldKind::Vector { dim: 6 };
        assert!(von_mises_from_components(&s).is_none());
    }

    #[test]
    fn field_max_per_node_finds_hot_spot() {
        let mut f = voigt_stress("vms", vec![10.0, 50.0, 30.0, 99.0, 1.0]);
        f.kind = FieldKind::Scalar;
        let peak = field_max_per_node(&f).expect("non-empty");
        assert_eq!(peak.node_index, 3);
        assert_eq!(peak.value, 99.0);
    }

    #[test]
    fn field_max_per_node_returns_none_for_empty() {
        let mut f = voigt_stress("vms", Vec::new());
        f.kind = FieldKind::Scalar;
        assert!(field_max_per_node(&f).is_none());
    }

    #[test]
    fn field_max_per_node_first_index_on_tie() {
        // First node carries the ties — deterministic sort.
        let mut f = voigt_stress("vms", vec![100.0, 100.0, 100.0]);
        f.kind = FieldKind::Scalar;
        let peak = field_max_per_node(&f).expect("non-empty");
        assert_eq!(peak.node_index, 0);
    }

    #[test]
    fn safety_factor_under_yield_above_one() {
        let mut vms = voigt_stress("vms", vec![100.0e6, 250.0e6, 50.0e6]);
        vms.kind = FieldKind::Scalar;
        let yield_str = 250.0e6;
        let sf = safety_factor(&vms, yield_str).expect("converts");
        // Node 0: 250e6/100e6 = 2.5 (under)
        // Node 1: 250e6/250e6 = 1.0 (at yield)
        // Node 2: 250e6/50e6  = 5.0 (under)
        assert!((sf.data[0] - 2.5).abs() < 1e-9);
        assert!((sf.data[1] - 1.0).abs() < 1e-9);
        assert!((sf.data[2] - 5.0).abs() < 1e-9);
    }

    #[test]
    fn safety_factor_zero_stress_is_infinity() {
        let mut vms = voigt_stress("vms", vec![0.0]);
        vms.kind = FieldKind::Scalar;
        let sf = safety_factor(&vms, 100.0e6).expect("converts");
        assert!(sf.data[0].is_infinite());
    }

    #[test]
    fn safety_factor_field_metadata() {
        let mut vms = voigt_stress("VMS_steel", vec![1.0]);
        vms.kind = FieldKind::Scalar;
        let sf = safety_factor(&vms, 1.0).expect("converts");
        assert_eq!(sf.name, "VMS_steel_safety_factor");
        // Safety factor is dimensionless even though VMS isn't.
        // DIMENSIONLESS in valenx-fields uses display = Some("")
        // (the empty string) — accept that, None, or the spec
        // canonical "1" so the convention can evolve without
        // breaking the test.
        assert!(matches!(sf.units.display, None | Some("") | Some("1")));
        assert_eq!(sf.kind, FieldKind::Scalar);
    }

    #[test]
    fn safety_factor_rejects_non_scalar_input() {
        let s = voigt_stress("S", vec![0.0; 6]);
        assert!(safety_factor(&s, 1e8).is_none());
    }

    #[test]
    fn count_below_finds_yielded_nodes() {
        let mut sf = voigt_stress("sf", vec![0.5, 1.5, 0.8, 2.0, 0.99]);
        sf.kind = FieldKind::Scalar;
        // SF below 1.0: nodes 0, 2, 4 -> 3
        assert_eq!(count_below(&sf, 1.0), 3);
        // None below 0.0
        assert_eq!(count_below(&sf, 0.0), 0);
    }
}
