//! Block constraint — pin every coordinate variable of an entity to
//! its current value, freezing the entity in place.
//!
//! Phase 12B Task 11.
//!
//! Residuals = `vars[i] - frozen_value[i]` for every variable the
//! entity references. For a point: 2 residuals. For a line: 4. For
//! a circle: 3 (center.x, center.y, radius). For an arc: 5.
//!
//! The "frozen value" is the value of the variable at the moment the
//! Block constraint was added (encoded in the constraint variant via
//! a `Vec<(var_idx, value)>` payload).

use crate::sketch::Sketch;

/// Returns the number of residuals for a frozen var list.
pub fn n_residuals(frozen: &[(usize, f64)]) -> usize {
    frozen.len()
}

/// Residual = `vars[idx] - target` for every `(idx, target)`.
pub fn residuals(sketch: &Sketch, frozen: &[(usize, f64)], out: &mut [f64]) {
    for (i, (var, target)) in frozen.iter().enumerate() {
        if let Some(v) = sketch.vars.get(*var) {
            out[i] = *v - *target;
        } else {
            out[i] = 0.0;
        }
    }
}

/// Jacobian: identity per row — one +1 entry on the variable.
pub fn jacobian(
    _sketch: &Sketch,
    frozen: &[(usize, f64)],
    triplets: &mut Vec<(usize, usize, f64)>,
) {
    for (i, (var, _)) in frozen.iter().enumerate() {
        triplets.push((i, *var, 1.0));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn block_residuals_zero_when_vars_unchanged() {
        let mut s = Sketch::new();
        let _ = s.add_point(3.0, 4.0);
        let frozen = vec![(0, 3.0), (1, 4.0)];
        let mut out = vec![0.0; 2];
        residuals(&s, &frozen, &mut out);
        assert!(out.iter().all(|r| r.abs() < 1e-12));
    }

    #[test]
    fn block_residuals_match_delta_when_changed() {
        let mut s = Sketch::new();
        let _ = s.add_point(3.0, 4.0);
        s.vars[0] = 5.0;
        let frozen = vec![(0, 3.0), (1, 4.0)];
        let mut out = vec![0.0; 2];
        residuals(&s, &frozen, &mut out);
        assert!((out[0] - 2.0).abs() < 1e-12);
        assert!(out[1].abs() < 1e-12);
    }
}
