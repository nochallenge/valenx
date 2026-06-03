//! Parameter scans — feature 18.
//!
//! A parameter scan sweeps one ([`scan_1d`]) or two ([`scan_2d`])
//! model parameters across a grid of values and records a scalar
//! **readout** at every grid point. The readout is whatever the caller
//! cares about — a steady-state concentration, an oscillation
//! amplitude, a final flux — supplied as a closure over the modified
//! [`Model`].
//!
//! This is the COPASI "parameter scan" task and the substrate of the
//! bifurcation and sensitivity modules. The scan itself is pure
//! bookkeeping; the science is in the readout closure.

use crate::analysis::param::ParamTarget;
use crate::error::{Result, SysbioError};
use crate::model::Model;

/// Result of a 1-D parameter scan.
#[derive(Debug, Clone, PartialEq)]
pub struct Scan1d {
    /// Label of the swept parameter.
    pub param: String,
    /// The grid of parameter values.
    pub values: Vec<f64>,
    /// The readout at each value (same length as `values`).
    pub readouts: Vec<f64>,
}

impl Scan1d {
    /// The `(value, readout)` pair with the largest readout.
    pub fn argmax(&self) -> Option<(f64, f64)> {
        self.values
            .iter()
            .zip(&self.readouts)
            .copied_pairs_max()
    }
}

/// Helper trait — `argmax` over a zipped iterator of `(f64, f64)`.
trait CopiedPairsMax {
    fn copied_pairs_max(self) -> Option<(f64, f64)>;
}

impl<'a, I> CopiedPairsMax for I
where
    I: Iterator<Item = (&'a f64, &'a f64)>,
{
    fn copied_pairs_max(self) -> Option<(f64, f64)> {
        let mut best: Option<(f64, f64)> = None;
        for (&v, &r) in self {
            if best.map(|(_, br)| r > br).unwrap_or(true) {
                best = Some((v, r));
            }
        }
        best
    }
}

/// Result of a 2-D parameter scan.
#[derive(Debug, Clone, PartialEq)]
pub struct Scan2d {
    /// Label of the first (row) parameter.
    pub param_a: String,
    /// Label of the second (column) parameter.
    pub param_b: String,
    /// Grid values of the first parameter.
    pub values_a: Vec<f64>,
    /// Grid values of the second parameter.
    pub values_b: Vec<f64>,
    /// `readouts[i][j]` for `values_a[i]`, `values_b[j]`.
    pub readouts: Vec<Vec<f64>>,
}

/// Build an inclusive linear grid of `n` points from `lo` to `hi`.
///
/// `n == 1` yields `[lo]`. Errors if `lo > hi` or `n == 0`.
pub fn linspace(lo: f64, hi: f64, n: usize) -> Result<Vec<f64>> {
    if n == 0 {
        return Err(SysbioError::invalid("n", "grid needs at least one point"));
    }
    if lo > hi {
        return Err(SysbioError::invalid("range", "lo must not exceed hi"));
    }
    if n == 1 {
        return Ok(vec![lo]);
    }
    Ok((0..n)
        .map(|k| lo + (hi - lo) * k as f64 / (n - 1) as f64)
        .collect())
}

/// Sweep one parameter across a linear grid, recording `readout` at
/// each point (feature 18, 1-D).
pub fn scan_1d<F>(
    model: &Model,
    target: &ParamTarget,
    lo: f64,
    hi: f64,
    n: usize,
    mut readout: F,
) -> Result<Scan1d>
where
    F: FnMut(&Model) -> Result<f64>,
{
    let values = linspace(lo, hi, n)?;
    let mut readouts = Vec::with_capacity(values.len());
    for &v in &values {
        let modified = target.apply(model, v)?;
        readouts.push(readout(&modified)?);
    }
    Ok(Scan1d {
        param: target.label(),
        values,
        readouts,
    })
}

/// Sweep two parameters across a rectangular grid (feature 18, 2-D).
#[allow(clippy::too_many_arguments)]
pub fn scan_2d<F>(
    model: &Model,
    target_a: &ParamTarget,
    a_lo: f64,
    a_hi: f64,
    a_n: usize,
    target_b: &ParamTarget,
    b_lo: f64,
    b_hi: f64,
    b_n: usize,
    mut readout: F,
) -> Result<Scan2d>
where
    F: FnMut(&Model) -> Result<f64>,
{
    let values_a = linspace(a_lo, a_hi, a_n)?;
    let values_b = linspace(b_lo, b_hi, b_n)?;
    let mut readouts = Vec::with_capacity(values_a.len());
    for &va in &values_a {
        let m_a = target_a.apply(model, va)?;
        let mut row = Vec::with_capacity(values_b.len());
        for &vb in &values_b {
            let m_ab = target_b.apply(&m_a, vb)?;
            row.push(readout(&m_ab)?);
        }
        readouts.push(row);
    }
    Ok(Scan2d {
        param_a: target_a.label(),
        param_b: target_b.label(),
        values_a,
        values_b,
        readouts,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{RateLaw, Reaction, Species};
    use crate::ode::{steady_state, OdeSystem};

    /// Source/decay model: 0 ->(rate s) A ->(k A) 0. Steady A* = s/k.
    fn source_decay() -> Model {
        let mut m = Model::new("sd");
        let a = m.add_species(Species::new("A", 0.0));
        m.add_reaction(Reaction {
            id: "src".into(),
            reactants: vec![],
            products: vec![(a, 1.0)],
            rate_law: RateLaw::Constant { rate: 1.0 },
            reversible: false,
        });
        m.add_reaction(Reaction {
            id: "dec".into(),
            reactants: vec![(a, 1.0)],
            products: vec![],
            rate_law: RateLaw::MassAction {
                k: 1.0,
                reactants: vec![(a, 1.0)],
            },
            reversible: false,
        });
        m
    }

    fn steady_a(m: &Model) -> Result<f64> {
        let sys = OdeSystem::from_model(m);
        Ok(steady_state(&sys, &m.initial_state(), 1e-9, 200)?.state[0])
    }

    #[test]
    fn linspace_endpoints() {
        let g = linspace(0.0, 10.0, 5).unwrap();
        assert_eq!(g, vec![0.0, 2.5, 5.0, 7.5, 10.0]);
        assert_eq!(linspace(3.0, 3.0, 1).unwrap(), vec![3.0]);
    }

    #[test]
    fn linspace_rejects_inverted_range() {
        assert!(linspace(5.0, 1.0, 3).is_err());
    }

    #[test]
    fn scan_1d_of_source_rate_is_linear() {
        // Sweeping the source rate s with k=1 gives A* = s.
        let m = source_decay();
        let scan = scan_1d(
            &m,
            &ParamTarget::ConstantRate { reaction: 0 },
            1.0,
            5.0,
            5,
            steady_a,
        )
        .unwrap();
        for (v, r) in scan.values.iter().zip(&scan.readouts) {
            assert!((v - r).abs() < 1e-4, "A*={r} should equal s={v}");
        }
        let (best_v, best_r) = scan.argmax().unwrap();
        assert!((best_v - 5.0).abs() < 1e-9);
        assert!((best_r - 5.0).abs() < 1e-4);
    }

    #[test]
    fn scan_2d_grid_shape_and_values() {
        let m = source_decay();
        let scan = scan_2d(
            &m,
            &ParamTarget::ConstantRate { reaction: 0 },
            2.0,
            4.0,
            3,
            &ParamTarget::MassActionK { reaction: 1 },
            1.0,
            2.0,
            2,
            steady_a,
        )
        .unwrap();
        assert_eq!(scan.readouts.len(), 3);
        assert_eq!(scan.readouts[0].len(), 2);
        // A* = s / k. Check one corner: s=4, k=2 -> 2.
        assert!((scan.readouts[2][1] - 2.0).abs() < 1e-4);
    }
}
