//! Field integral / statistics helpers for CFD + FEA post-
//! processing.
//!
//! Quick scalar-out-of-field-in primitives the user wants
//! immediately after a run: average pressure, peak velocity, total
//! kinetic energy, L2 residual norm. Each function takes a `Field`
//! reference + (where needed) per-node weights, and returns one
//! scalar.
//!
//! All helpers are pure / allocation-free apart from the borrow.
//! They handle empty fields by returning a sensible zero (so a
//! dashboard widget showing "average T = NaN" doesn't surface
//! when no data is loaded yet).

use crate::{Field, FieldKind};

/// Sum of every value in the field's `data` buffer. For scalar
/// fields this is the integral assuming unit weights — pair with
/// [`integrate_with_node_volumes`] when the per-node measure
/// matters.
///
/// Returns `0.0` for empty fields rather than NaN so dashboard
/// widgets stay readable.
pub fn field_sum(field: &Field) -> f64 {
    field.data.iter().copied().sum()
}

/// Arithmetic mean of the field's `data` buffer. Returns `0.0` for
/// empty fields. For per-cell fields this is the volume-weighted
/// mean only when cells are equal-volume; pair with
/// [`mean_with_weights`] otherwise.
pub fn field_mean(field: &Field) -> f64 {
    if field.data.is_empty() {
        return 0.0;
    }
    field_sum(field) / field.data.len() as f64
}

/// Population variance of the field. Uses N (not N-1) — stratifies
/// the field as the population, not a sample. Returns `0.0` for
/// empty / single-point fields.
pub fn field_variance(field: &Field) -> f64 {
    let n = field.data.len();
    if n < 2 {
        return 0.0;
    }
    let mean = field_mean(field);
    let sum_sq: f64 = field.data.iter().map(|v| (v - mean).powi(2)).sum();
    sum_sq / n as f64
}

/// Standard deviation = sqrt(variance).
pub fn field_std_dev(field: &Field) -> f64 {
    field_variance(field).sqrt()
}

/// L1 norm: sum of absolute values. Useful for scaled residual
/// magnitudes and total mass-flux balance checks.
pub fn field_l1_norm(field: &Field) -> f64 {
    field.data.iter().map(|v| v.abs()).sum()
}

/// L2 norm: sqrt(sum of squares). Common residual metric in
/// iterative solvers.
pub fn field_l2_norm(field: &Field) -> f64 {
    field.data.iter().map(|v| v * v).sum::<f64>().sqrt()
}

/// L-infinity norm: max absolute value. Pairs with
/// [`crate::stress::field_max_per_node`] when the user wants the
/// magnitude (not the signed value) of the worst cell.
pub fn field_linf_norm(field: &Field) -> f64 {
    field.data.iter().map(|v| v.abs()).fold(0.0_f64, f64::max)
}

/// (min, max) over the field's `data`. Returns `None` for empty
/// fields. Lifted out of the VTU parser's private helper so other
/// callers can use it without re-implementing the loop.
pub fn field_min_max(field: &Field) -> Option<(f64, f64)> {
    let mut iter = field.data.iter().copied();
    let first = iter.next()?;
    let mut min = first;
    let mut max = first;
    for v in iter {
        if v < min {
            min = v;
        }
        if v > max {
            max = v;
        }
    }
    Some((min, max))
}

/// Weighted mean over a scalar field. `weights` must be the same
/// length as `field.data`; mismatches return `None`. Use case:
/// volume-weighted average pressure across non-uniform cells.
pub fn mean_with_weights(field: &Field, weights: &[f64]) -> Option<f64> {
    if !matches!(field.kind, FieldKind::Scalar) {
        return None;
    }
    if field.data.len() != weights.len() {
        return None;
    }
    if weights.is_empty() {
        return Some(0.0);
    }
    let total_weight: f64 = weights.iter().sum();
    if total_weight == 0.0 {
        return Some(0.0);
    }
    let weighted: f64 = field
        .data
        .iter()
        .zip(weights.iter())
        .map(|(v, w)| v * w)
        .sum();
    Some(weighted / total_weight)
}

/// Integrate a field over per-node volumes: sum(value × volume).
/// Used for total mass / total energy when each node carries a
/// scalar value and a known control volume.
pub fn integrate_with_node_volumes(field: &Field, volumes: &[f64]) -> Option<f64> {
    if !matches!(field.kind, FieldKind::Scalar) {
        return None;
    }
    if field.data.len() != volumes.len() {
        return None;
    }
    Some(
        field
            .data
            .iter()
            .zip(volumes.iter())
            .map(|(v, vol)| v * vol)
            .sum(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Location, RegionRef, TimeKey};

    fn scalar(name: &str, data: Vec<f64>) -> Field {
        Field {
            name: name.into(),
            kind: FieldKind::Scalar,
            location: Location::OnNode,
            region: RegionRef("default".into()),
            units: crate::units::DIMENSIONLESS,
            time: TimeKey::Steady,
            data,
            range: None,
        }
    }

    #[test]
    fn field_sum_handles_empty_and_nonempty() {
        assert_eq!(field_sum(&scalar("e", vec![])), 0.0);
        assert_eq!(field_sum(&scalar("a", vec![1.0, 2.0, 3.0])), 6.0);
        assert_eq!(field_sum(&scalar("b", vec![-1.0, 1.0])), 0.0);
    }

    #[test]
    fn field_mean_handles_empty() {
        assert_eq!(field_mean(&scalar("e", vec![])), 0.0);
        assert_eq!(field_mean(&scalar("a", vec![2.0, 4.0])), 3.0);
        assert_eq!(field_mean(&scalar("b", vec![10.0, 10.0, 10.0])), 10.0);
    }

    #[test]
    fn field_variance_zero_for_constant_field() {
        assert_eq!(field_variance(&scalar("c", vec![5.0, 5.0, 5.0, 5.0])), 0.0);
        // Single-element + empty are zero by definition.
        assert_eq!(field_variance(&scalar("e", vec![])), 0.0);
        assert_eq!(field_variance(&scalar("s", vec![42.0])), 0.0);
    }

    #[test]
    fn field_variance_known_population() {
        // [1, 2, 3, 4, 5] -> mean=3, variance = (4+1+0+1+4)/5 = 2.0
        assert!((field_variance(&scalar("v", vec![1.0, 2.0, 3.0, 4.0, 5.0])) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn field_std_dev_is_sqrt_variance() {
        let f = scalar("v", vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        assert!((field_std_dev(&f) - 2.0_f64.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn field_l1_norm_sums_absolute_values() {
        assert_eq!(
            field_l1_norm(&scalar("a", vec![1.0, -2.0, 3.0, -4.0])),
            10.0
        );
        assert_eq!(field_l1_norm(&scalar("e", vec![])), 0.0);
    }

    #[test]
    fn field_l2_norm_sqrt_sum_of_squares() {
        // [3, 4] -> sqrt(9+16) = 5
        assert!((field_l2_norm(&scalar("v", vec![3.0, 4.0])) - 5.0).abs() < 1e-12);
        assert_eq!(field_l2_norm(&scalar("e", vec![])), 0.0);
    }

    #[test]
    fn field_linf_norm_max_abs() {
        assert_eq!(field_linf_norm(&scalar("a", vec![1.0, -7.0, 3.0])), 7.0);
        assert_eq!(field_linf_norm(&scalar("e", vec![])), 0.0);
    }

    #[test]
    fn field_min_max_empty_returns_none() {
        assert!(field_min_max(&scalar("e", vec![])).is_none());
        assert_eq!(
            field_min_max(&scalar("a", vec![3.0, 1.0, 2.0])),
            Some((1.0, 3.0))
        );
    }

    #[test]
    fn mean_with_weights_volume_weighted_average() {
        // Field: [10, 20], weights: [3, 7] -> (30 + 140) / 10 = 17
        let f = scalar("p", vec![10.0, 20.0]);
        let w = vec![3.0, 7.0];
        assert!((mean_with_weights(&f, &w).unwrap() - 17.0).abs() < 1e-12);
    }

    #[test]
    fn mean_with_weights_rejects_size_mismatch() {
        let f = scalar("a", vec![1.0, 2.0]);
        assert!(mean_with_weights(&f, &[1.0]).is_none());
    }

    #[test]
    fn mean_with_weights_zero_weight_returns_zero() {
        let f = scalar("a", vec![10.0, 20.0]);
        let w = vec![0.0, 0.0];
        assert_eq!(mean_with_weights(&f, &w), Some(0.0));
    }

    #[test]
    fn mean_with_weights_rejects_non_scalar() {
        let mut f = scalar("v", vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        f.kind = FieldKind::Vector { dim: 3 };
        assert!(mean_with_weights(&f, &[1.0; 6]).is_none());
    }

    #[test]
    fn integrate_with_node_volumes_sums_value_times_vol() {
        // Field: [10, 20, 30], volumes: [0.1, 0.2, 0.3]
        // -> 1 + 4 + 9 = 14
        let f = scalar("rho", vec![10.0, 20.0, 30.0]);
        let v = vec![0.1, 0.2, 0.3];
        assert!((integrate_with_node_volumes(&f, &v).unwrap() - 14.0).abs() < 1e-12);
    }

    #[test]
    fn integrate_with_node_volumes_rejects_size_mismatch() {
        let f = scalar("a", vec![1.0, 2.0]);
        assert!(integrate_with_node_volumes(&f, &[1.0]).is_none());
    }

    #[test]
    fn integrate_with_node_volumes_rejects_non_scalar() {
        let mut f = scalar("v", vec![1.0, 2.0, 3.0]);
        f.kind = FieldKind::Vector { dim: 3 };
        assert!(integrate_with_node_volumes(&f, &[1.0]).is_none());
    }

    /// Sanity: the `field_min_max` helper agrees with the
    /// `field_l1_norm`-based bound `|min|, |max| <= L1`. Catches
    /// regressions where one helper drifts from the other.
    #[test]
    fn helpers_agree_on_simple_fields() {
        let f = scalar("v", vec![-3.0, 2.0, -1.0, 4.0]);
        let (min, max) = field_min_max(&f).unwrap();
        let l1 = field_l1_norm(&f);
        let l2 = field_l2_norm(&f);
        let linf = field_linf_norm(&f);
        // L-inf == max(|min|, |max|)
        assert_eq!(linf, min.abs().max(max.abs()));
        // L2 <= L1 (Cauchy-Schwarz upper bound for finite sequences)
        assert!(l2 <= l1 + 1e-12);
    }
}
