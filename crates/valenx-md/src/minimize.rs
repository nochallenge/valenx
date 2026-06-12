//! Energy minimisation — **roadmap features 24 & 25**.
//!
//! Before a dynamics run a structure is **energy-minimised** to remove
//! bad contacts (atoms too close, strained bonds) that would otherwise
//! make the integrator blow up. Minimisation walks the atoms downhill
//! on the potential-energy surface until the force is small.
//!
//! Three optimisers share a common driver. Each takes the same
//! energy+force callback the integrators use:
//!
//! - [`steepest_descent`] — **steepest descent** (feature 24). Step
//!   along the force (the negative gradient) with an adaptive step
//!   size: grow it after a successful (energy-lowering) step, shrink
//!   it after a failed one. Robust, never diverges, but slow in long
//!   narrow valleys.
//!
//! - [`conjugate_gradient`] — **conjugate gradient** (feature 25).
//!   Builds each search direction from the current force *plus* a
//!   Polak-Ribière-corrected fraction of the previous direction, so it
//!   does not undo its own progress. Much faster than steepest descent
//!   on quadratic-like surfaces.
//!
//! - [`lbfgs`] — **limited-memory BFGS** (feature 25). Reconstructs an
//!   implicit inverse-Hessian from the last `m` position/gradient
//!   changes (the two-loop recursion) and steps with it — quasi-Newton
//!   convergence at `O(m·N)` memory. The fastest of the three near a
//!   minimum.
//!
//! All three use a backtracking line search and stop when the maximum
//! force drops below a tolerance or the step count is exhausted.

use nalgebra::Vector3;

use crate::bonded::EnergyForce;
use crate::error::{MdError, Result};
use crate::system::System;

/// The outcome of a minimisation run.
#[derive(Clone, Debug, PartialEq)]
pub struct MinimizationResult {
    /// Potential energy at the minimised structure (kJ/mol).
    pub final_energy: f64,
    /// Maximum force on any atom at the end (kJ/(mol·nm)).
    pub final_max_force: f64,
    /// Iterations performed.
    pub iterations: usize,
    /// Whether the force tolerance was reached.
    pub converged: bool,
}

/// Settings shared by all three minimisers.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct MinimizeOptions {
    /// Stop when the maximum force is below this (kJ/(mol·nm)).
    pub force_tolerance: f64,
    /// Maximum iterations.
    pub max_iterations: usize,
    /// Initial trial step size (nm).
    pub initial_step: f64,
}

impl Default for MinimizeOptions {
    fn default() -> Self {
        MinimizeOptions {
            force_tolerance: 10.0,
            max_iterations: 1000,
            initial_step: 0.01,
        }
    }
}

impl MinimizeOptions {
    /// Validates the options.
    fn check(&self) -> Result<()> {
        if !(self.force_tolerance.is_finite() && self.force_tolerance > 0.0) {
            return Err(MdError::invalid(
                "force_tolerance",
                "must be finite and positive",
            ));
        }
        if self.max_iterations == 0 {
            return Err(MdError::invalid("max_iterations", "must be at least 1"));
        }
        if !(self.initial_step.is_finite() && self.initial_step > 0.0) {
            return Err(MdError::invalid(
                "initial_step",
                "must be finite and positive",
            ));
        }
        Ok(())
    }
}

/// Evaluates the energy of a trial configuration: `base` displaced by
/// `alpha` along `direction`.
fn energy_along(
    base: &System,
    direction: &[Vector3<f64>],
    alpha: f64,
    force_fn: &mut dyn FnMut(&System) -> Result<EnergyForce>,
) -> Result<(EnergyForce, System)> {
    let mut trial = base.clone();
    for (p, d) in trial.positions.iter_mut().zip(direction) {
        *p += alpha * d;
    }
    let ef = force_fn(&trial)?;
    Ok((ef, trial))
}

/// Steepest-descent minimisation — **feature 24**.
///
/// Steps along the force with an adaptive step size.
///
/// # Errors
/// [`MdError::Invalid`] for bad options; propagates `force_fn` errors.
pub fn steepest_descent(
    system: &mut System,
    options: MinimizeOptions,
    force_fn: &mut dyn FnMut(&System) -> Result<EnergyForce>,
) -> Result<MinimizationResult> {
    options.check()?;
    let mut ef = force_fn(system)?;
    let mut step = options.initial_step;
    let mut iterations = 0;

    for _ in 0..options.max_iterations {
        let max_f = ef.max_force();
        if max_f < options.force_tolerance {
            return Ok(MinimizationResult {
                final_energy: ef.energy,
                final_max_force: max_f,
                iterations,
                converged: true,
            });
        }
        // Search direction is the (normalised) force.
        let norm = max_f.max(1e-12);
        let direction: Vec<Vector3<f64>> = ef.forces.iter().map(|f| f / norm).collect();
        // Try a step; backtrack until the energy goes down.
        let (new_ef, new_sys) = line_search(system, &direction, step, ef.energy, force_fn)?;
        if new_ef.energy < ef.energy {
            *system = new_sys;
            ef = new_ef;
            step *= 1.2; // success: be bolder
        } else {
            step *= 0.5; // failure: be cautious
            if step < 1e-10 {
                break;
            }
        }
        iterations += 1;
    }
    Ok(MinimizationResult {
        final_energy: ef.energy,
        final_max_force: ef.max_force(),
        iterations,
        converged: ef.max_force() < options.force_tolerance,
    })
}

/// A backtracking line search: halve `alpha` until the energy along
/// `direction` improves on `e0`, or give up after a few tries (in
/// which case the best — possibly non-improving — trial is returned).
fn line_search(
    base: &System,
    direction: &[Vector3<f64>],
    alpha0: f64,
    e0: f64,
    force_fn: &mut dyn FnMut(&System) -> Result<EnergyForce>,
) -> Result<(EnergyForce, System)> {
    let mut alpha = alpha0;
    let mut best: Option<(EnergyForce, System)> = None;
    for _ in 0..20 {
        let (ef, sys) = energy_along(base, direction, alpha, force_fn)?;
        if ef.energy < e0 {
            return Ok((ef, sys));
        }
        if best
            .as_ref()
            .map(|(b, _)| ef.energy < b.energy)
            .unwrap_or(true)
        {
            best = Some((ef, sys));
        }
        alpha *= 0.5;
    }
    Ok(best.unwrap())
}

/// Conjugate-gradient minimisation — **feature 25**.
///
/// Polak-Ribière conjugate gradient with a restart whenever the
/// correction would point uphill.
///
/// # Errors
/// [`MdError::Invalid`] for bad options; propagates `force_fn` errors.
pub fn conjugate_gradient(
    system: &mut System,
    options: MinimizeOptions,
    force_fn: &mut dyn FnMut(&System) -> Result<EnergyForce>,
) -> Result<MinimizationResult> {
    options.check()?;
    let mut ef = force_fn(system)?;
    let n = system.len();
    // Search direction; initialised to the steepest-descent direction.
    let mut direction: Vec<Vector3<f64>> = ef.forces.clone();
    let mut prev_grad: Vec<Vector3<f64>> = ef.forces.clone();
    let mut iterations = 0;
    let mut step = options.initial_step;

    for _ in 0..options.max_iterations {
        let max_f = ef.max_force();
        if max_f < options.force_tolerance {
            return Ok(MinimizationResult {
                final_energy: ef.energy,
                final_max_force: max_f,
                iterations,
                converged: true,
            });
        }
        // Normalise the search direction for a scale-stable step.
        let dnorm = direction
            .iter()
            .map(|d| d.norm_squared())
            .sum::<f64>()
            .sqrt()
            .max(1e-12);
        let unit: Vec<Vector3<f64>> = direction.iter().map(|d| d / dnorm).collect();
        let (new_ef, new_sys) = line_search(system, &unit, step, ef.energy, force_fn)?;
        if new_ef.energy < ef.energy {
            *system = new_sys;
            step *= 1.1;
        } else {
            step *= 0.5;
            if step < 1e-10 {
                break;
            }
            iterations += 1;
            continue;
        }
        // Polak-Ribière beta. Gradient g = -force.
        let new_grad = &new_ef.forces;
        let mut numerator = 0.0;
        let mut denominator = 0.0;
        for i in 0..n {
            let g_new = -new_grad[i];
            let g_old = -prev_grad[i];
            numerator += g_new.dot(&(g_new - g_old));
            denominator += g_old.dot(&g_old);
        }
        let beta = if denominator > 1e-18 {
            (numerator / denominator).max(0.0) // PR+ (restart if < 0)
        } else {
            0.0
        };
        // New direction: force + beta * old direction.
        for i in 0..n {
            direction[i] = new_grad[i] + beta * direction[i];
        }
        prev_grad = new_ef.forces.clone();
        ef = new_ef;
        iterations += 1;
    }
    Ok(MinimizationResult {
        final_energy: ef.energy,
        final_max_force: ef.max_force(),
        iterations,
        converged: ef.max_force() < options.force_tolerance,
    })
}

/// Limited-memory BFGS minimisation — **feature 25**.
///
/// `history` is the number of (position-change, gradient-change) pairs
/// kept for the implicit inverse-Hessian (`5`–`10` is usual).
///
/// # Errors
/// [`MdError::Invalid`] for bad options or a zero history;
/// propagates `force_fn` errors.
pub fn lbfgs(
    system: &mut System,
    options: MinimizeOptions,
    history: usize,
    force_fn: &mut dyn FnMut(&System) -> Result<EnergyForce>,
) -> Result<MinimizationResult> {
    options.check()?;
    if history == 0 {
        return Err(MdError::invalid("history", "must be at least 1"));
    }
    let mut ef = force_fn(system)?;
    // Stored differences: s = Δposition, y = Δgradient.
    let mut s_hist: Vec<Vec<Vector3<f64>>> = Vec::new();
    let mut y_hist: Vec<Vec<Vector3<f64>>> = Vec::new();
    let mut rho_hist: Vec<f64> = Vec::new();
    let mut iterations = 0;
    let mut step = options.initial_step;

    for _ in 0..options.max_iterations {
        let max_f = ef.max_force();
        if max_f < options.force_tolerance {
            return Ok(MinimizationResult {
                final_energy: ef.energy,
                final_max_force: max_f,
                iterations,
                converged: true,
            });
        }
        // Gradient g = -force.
        let grad: Vec<Vector3<f64>> = ef.forces.iter().map(|f| -f).collect();
        // --- Two-loop recursion: q = -H·g (the descent direction) ---
        let mut q = grad.clone();
        let m = s_hist.len();
        let mut alpha = vec![0.0; m];
        for k in (0..m).rev() {
            let a = rho_hist[k] * dot(&s_hist[k], &q);
            alpha[k] = a;
            for (qi, yi) in q.iter_mut().zip(&y_hist[k]) {
                *qi -= a * yi;
            }
        }
        // Initial Hessian scaling.
        let gamma = if m > 0 {
            let last = m - 1;
            let sy = dot(&s_hist[last], &y_hist[last]);
            let yy = dot(&y_hist[last], &y_hist[last]).max(1e-18);
            (sy / yy).max(1e-8)
        } else {
            // First step: a small scaled steepest-descent.
            1.0
        };
        for qi in &mut q {
            *qi *= gamma;
        }
        for k in 0..m {
            let b = rho_hist[k] * dot(&y_hist[k], &q);
            for (qi, si) in q.iter_mut().zip(&s_hist[k]) {
                *qi += (alpha[k] - b) * si;
            }
        }
        // Descent direction is -q (q approximates H·g).
        let direction: Vec<Vector3<f64>> = q.iter().map(|v| -v).collect();
        let dnorm = direction
            .iter()
            .map(|d| d.norm_squared())
            .sum::<f64>()
            .sqrt()
            .max(1e-12);
        let unit: Vec<Vector3<f64>> = direction.iter().map(|d| d / dnorm).collect();

        let (new_ef, new_sys) = line_search(system, &unit, step, ef.energy, force_fn)?;
        if new_ef.energy >= ef.energy {
            step *= 0.5;
            if step < 1e-10 {
                break;
            }
            iterations += 1;
            continue;
        }
        // Record the curvature pair.
        let s: Vec<Vector3<f64>> = new_sys
            .positions
            .iter()
            .zip(&system.positions)
            .map(|(a, b)| a - b)
            .collect();
        let new_grad: Vec<Vector3<f64>> = new_ef.forces.iter().map(|f| -f).collect();
        let y: Vec<Vector3<f64>> = new_grad.iter().zip(&grad).map(|(a, b)| a - b).collect();
        let sy = dot(&s, &y);
        if sy > 1e-12 {
            s_hist.push(s);
            y_hist.push(y);
            rho_hist.push(1.0 / sy);
            if s_hist.len() > history {
                s_hist.remove(0);
                y_hist.remove(0);
                rho_hist.remove(0);
            }
        }
        *system = new_sys;
        ef = new_ef;
        step = (step * 1.1).min(options.initial_step * 10.0);
        iterations += 1;
    }
    Ok(MinimizationResult {
        final_energy: ef.energy,
        final_max_force: ef.max_force(),
        iterations,
        converged: ef.max_force() < options.force_tolerance,
    })
}

/// Flat dot product of two per-atom vector arrays.
fn dot(a: &[Vector3<f64>], b: &[Vector3<f64>]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x.dot(y)).sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bonded::bond::HarmonicBonds;
    use crate::bonded::ForceTerm;
    use crate::forcefield::BondParam;
    use crate::system::{Atom, Topology};

    /// A strained diatomic: the bond starts well off its equilibrium.
    fn strained_bond(start_len: f64, eq_len: f64) -> (System, HarmonicBonds) {
        let mut top = Topology::new();
        top.push_atom(Atom::new("A", 12.0, 0.0).unwrap());
        top.push_atom(Atom::new("B", 12.0, 0.0).unwrap());
        top.add_bond(0, 1).unwrap();
        let sys = System::new(
            top,
            vec![Vector3::zeros(), Vector3::new(start_len, 0.0, 0.0)],
        )
        .unwrap();
        let term =
            HarmonicBonds::from_system(&sys, &[BondParam::new(eq_len, 5000.0).unwrap()]).unwrap();
        (sys, term)
    }

    fn force_closure(term: &HarmonicBonds) -> impl FnMut(&System) -> Result<EnergyForce> + '_ {
        move |s: &System| {
            let mut ef = EnergyForce::zeros(s.len());
            term.accumulate(s, &mut ef)?;
            Ok(ef)
        }
    }

    #[test]
    fn steepest_descent_relaxes_a_strained_bond() {
        let (mut sys, term) = strained_bond(0.20, 0.15);
        let mut force = force_closure(&term);
        let opts = MinimizeOptions {
            force_tolerance: 1.0,
            max_iterations: 2000,
            initial_step: 0.005,
        };
        let result = steepest_descent(&mut sys, opts, &mut force).unwrap();
        // The bond length should end near equilibrium.
        let len = (sys.positions[0] - sys.positions[1]).norm();
        assert!((len - 0.15).abs() < 0.01, "bond ended at {len}");
        assert!(result.final_energy < 1.0);
    }

    #[test]
    fn conjugate_gradient_relaxes_a_strained_bond() {
        let (mut sys, term) = strained_bond(0.22, 0.15);
        let mut force = force_closure(&term);
        let opts = MinimizeOptions {
            force_tolerance: 1.0,
            max_iterations: 2000,
            initial_step: 0.005,
        };
        let result = conjugate_gradient(&mut sys, opts, &mut force).unwrap();
        let len = (sys.positions[0] - sys.positions[1]).norm();
        assert!((len - 0.15).abs() < 0.01, "CG bond ended at {len}");
        assert!(result.iterations >= 1);
    }

    #[test]
    fn lbfgs_relaxes_a_strained_bond() {
        let (mut sys, term) = strained_bond(0.09, 0.15); // compressed
        let mut force = force_closure(&term);
        let opts = MinimizeOptions {
            force_tolerance: 1.0,
            max_iterations: 2000,
            initial_step: 0.005,
        };
        let result = lbfgs(&mut sys, opts, 6, &mut force).unwrap();
        let len = (sys.positions[0] - sys.positions[1]).norm();
        assert!((len - 0.15).abs() < 0.01, "L-BFGS bond ended at {len}");
        assert!(result.final_max_force < 5.0);
    }

    #[test]
    fn minimisers_lower_the_energy() {
        for which in 0..3 {
            let (mut sys, term) = strained_bond(0.25, 0.15);
            let mut force = force_closure(&term);
            let e_before = {
                let mut ef = EnergyForce::zeros(sys.len());
                term.accumulate(&sys, &mut ef).unwrap();
                ef.energy
            };
            let opts = MinimizeOptions::default();
            let result = match which {
                0 => steepest_descent(&mut sys, opts, &mut force).unwrap(),
                1 => conjugate_gradient(&mut sys, opts, &mut force).unwrap(),
                _ => lbfgs(&mut sys, opts, 5, &mut force).unwrap(),
            };
            assert!(
                result.final_energy < e_before,
                "minimiser {which}: {} !< {}",
                result.final_energy,
                e_before
            );
        }
    }

    #[test]
    fn rejects_bad_options() {
        let (mut sys, term) = strained_bond(0.2, 0.15);
        let mut force = force_closure(&term);
        let bad = MinimizeOptions {
            force_tolerance: -1.0,
            ..MinimizeOptions::default()
        };
        assert!(steepest_descent(&mut sys, bad, &mut force).is_err());
        assert!(lbfgs(&mut sys, MinimizeOptions::default(), 0, &mut force).is_err());
    }
}
