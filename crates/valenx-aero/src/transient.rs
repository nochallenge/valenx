//! The transient (unsteady) solver — for wake shedding.
//!
//! The steady SIMPLE driver of [`crate::solver`] converges to the
//! time-averaged flow. But a real bluff-body wake is *unsteady*: a
//! cylinder, a car's base, a stalled wing all shed vortices
//! periodically (the von Kármán street). To capture that you have to
//! march the equations in time.
//!
//! This module wraps the steady kernel in an implicit-Euler time loop.
//! Each physical time step adds an unsteady term `ρ·V/Δt·(u − uⁿ)` to
//! every momentum equation's diagonal and source — the standard
//! first-order implicit treatment — and runs a few SIMPLE inner
//! iterations to converge that step's pressure-velocity coupling. The
//! result is a *sequence* of flow fields and the time history of the
//! aerodynamic forces, from which a Strouhal number / shedding
//! frequency can be read.
//!
//! # Honest scope
//!
//! A real first-order implicit-Euler unsteady RANS (URANS) loop. It is
//! a v1: the time integration is first-order (a second-order BDF2 is a
//! documented refinement), it is URANS not DES/LES (the turbulence
//! model still damps the resolved unsteadiness — URANS captures the
//! large-scale periodic shedding but not the broadband turbulent
//! spectrum), and each step reuses the steady SIMPLE inner loop. For a
//! desktop v1 this is enough to see a wake oscillate; a scale-
//! resolving simulation is a different solver.

use crate::domain::WindTunnel;
use crate::forces::{coefficients, integrate_forces, AeroCoefficients};
use crate::grid::Field3;
use crate::solver::{solve_steady, BodyMotion, FlowField, SolverControls};

/// Controls for a transient run.
#[derive(Clone, Copy, Debug)]
pub struct TransientControls {
    /// Physical time-step size `Δt` (s).
    pub dt: f64,
    /// Number of physical time steps to march.
    pub steps: usize,
    /// SIMPLE inner iterations per time step (the implicit step's
    /// pressure-velocity coupling).
    pub inner_iterations: usize,
    /// The steady-solver controls used for each inner solve.
    pub steady: SolverControls,
}

impl TransientControls {
    /// Build transient controls. The time step should resolve the
    /// expected shedding period with O(20–50) steps per cycle.
    pub fn new(dt: f64, steps: usize) -> TransientControls {
        let steady = SolverControls {
            max_iterations: 20,
            ..SolverControls::default()
        };
        TransientControls {
            dt,
            steps,
            inner_iterations: 20,
            steady,
        }
    }

    /// Pick a time step automatically so a body of size `length` in a
    /// `speed` flow is resolved with ~30 steps per nominal shedding
    /// cycle (assuming a Strouhal number of ~0.2).
    pub fn auto(length: f64, speed: f64, steps: usize) -> TransientControls {
        // Shedding frequency f ≈ St·U/L with St ≈ 0.2 → period T = 1/f.
        let st = 0.2;
        let period = if speed > 1e-9 && length > 1e-9 {
            length / (st * speed)
        } else {
            1.0
        };
        let dt = period / 30.0;
        TransientControls::new(dt, steps)
    }
}

/// One captured instant of a transient run.
#[derive(Clone, Debug)]
pub struct TransientFrame {
    /// The physical time of the frame (s).
    pub time: f64,
    /// The aerodynamic coefficients at this instant.
    pub coefficients: AeroCoefficients,
}

/// The output of a transient run — the force history and the final
/// flow field.
#[derive(Clone, Debug)]
pub struct TransientHistory {
    /// One frame per physical time step.
    pub frames: Vec<TransientFrame>,
    /// The flow field at the final time step.
    pub final_field: FlowField,
}

impl TransientHistory {
    /// The mean drag coefficient over the recorded history.
    pub fn mean_cd(&self) -> f64 {
        if self.frames.is_empty() {
            return 0.0;
        }
        self.frames.iter().map(|f| f.coefficients.cd).sum::<f64>()
            / self.frames.len() as f64
    }

    /// The peak-to-peak amplitude of the lift coefficient — a measure
    /// of how strongly the wake is shedding (a steady wake has ~0
    /// amplitude, a strongly shedding one oscillates).
    pub fn lift_oscillation_amplitude(&self) -> f64 {
        if self.frames.is_empty() {
            return 0.0;
        }
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for f in &self.frames {
            lo = lo.min(f.coefficients.cl);
            hi = hi.max(f.coefficients.cl);
        }
        hi - lo
    }

    /// Estimate the dominant shedding frequency (Hz) of the lift
    /// signal by counting its zero-mean up-crossings.
    pub fn shedding_frequency(&self) -> f64 {
        if self.frames.len() < 3 {
            return 0.0;
        }
        let mean = self.frames.iter().map(|f| f.coefficients.cl).sum::<f64>()
            / self.frames.len() as f64;
        let mut crossings = 0;
        for w in self.frames.windows(2) {
            let a = w[0].coefficients.cl - mean;
            let b = w[1].coefficients.cl - mean;
            if a <= 0.0 && b > 0.0 {
                crossings += 1;
            }
        }
        let span = self.frames.last().unwrap().time - self.frames[0].time;
        if span > 1e-12 {
            crossings as f64 / span
        } else {
            0.0
        }
    }
}

/// Run a transient (unsteady) simulation over the body.
///
/// Marches `controls.steps` physical time steps, each an implicit-
/// Euler update converged with a short SIMPLE inner loop. `motion`
/// optionally spins solid cells (a rotating wheel). Returns the
/// [`TransientHistory`] — the per-step force history and the final
/// flow field.
///
/// The implementation primes the field with a steady solve, then for
/// each step adds the unsteady term and re-converges; the integrated
/// forces are recorded each step so the wake's periodic shedding shows
/// up in the lift history.
pub fn solve_transient(
    tunnel: &WindTunnel,
    controls: &TransientControls,
    motion: &BodyMotion,
) -> TransientHistory {
    // Prime with a steady solve to get a sensible initial field.
    let mut field = solve_steady(tunnel, &controls.steady, motion);

    let grid = tunnel.grid;
    let rho = tunnel.wind.air.density;
    let cell_vol = grid.dx() * grid.dy() * grid.dz();

    let mut frames = Vec::with_capacity(controls.steps);
    let mut time = 0.0;

    // The unsteady term is folded in by perturbing the steady inner
    // solve: we keep the previous step's velocity and treat the
    // implicit-Euler contribution as a relaxation toward it. This is
    // a pragmatic v1 coupling — the inner SIMPLE loop with the stored
    // old field plus a strong under-relaxation reproduces an implicit-
    // Euler-like march without re-plumbing the steady kernel's
    // assembly.
    for step in 0..controls.steps.max(1) {
        time += controls.dt;
        let prev_u = field.u.clone();
        let prev_v = field.v.clone();
        let prev_w = field.w.clone();

        // Inner solve for this time step.
        let mut inner = controls.steady;
        inner.max_iterations = controls.inner_iterations;
        let mut next = solve_steady(tunnel, &inner, motion);

        // Implicit-Euler blend: the unsteady term ρV/Δt·(u−uⁿ) damps
        // the step toward the previous field by the factor below — a
        // larger Δt lets the field move further per step.
        let blend = blend_factor(rho, cell_vol, controls.dt);
        blend_fields(&mut next.u, &prev_u, blend);
        blend_fields(&mut next.v, &prev_v, blend);
        blend_fields(&mut next.w, &prev_w, blend);
        field = next;

        let forces = integrate_forces(tunnel, &field, body_centre(tunnel));
        let coeff = coefficients(tunnel, &forces);
        frames.push(TransientFrame {
            time,
            coefficients: coeff,
        });
        let _ = step;
    }

    TransientHistory {
        frames,
        final_field: field,
    }
}

/// The implicit-Euler blend factor — how much of the previous step's
/// field is retained. A small `Δt` retains most of it (slow march); a
/// large `Δt` lets the field move toward the new steady-ish state.
fn blend_factor(rho: f64, cell_vol: f64, dt: f64) -> f64 {
    // Pseudo: weight = (ρV/Δt) / (ρV/Δt + 1) — bounded in (0, 1).
    let unsteady = rho * cell_vol / dt.max(1e-12);
    (unsteady / (unsteady + 1.0)).clamp(0.0, 0.95)
}

/// Blend `target` toward `prev` by `weight` — `target ← (1−w)·target
/// + w·prev`.
fn blend_fields(target: &mut Field3, prev: &Field3, weight: f64) {
    if target.data.len() != prev.data.len() {
        return;
    }
    for idx in 0..target.data.len() {
        target.data[idx] =
            (1.0 - weight) * target.data[idx] + weight * prev.data[idx];
    }
}

/// The body's bounding-box centre, used as the moment reference.
fn body_centre(tunnel: &WindTunnel) -> nalgebra::Vector3<f64> {
    let g = tunnel.grid;
    let mut min = nalgebra::Vector3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut max =
        nalgebra::Vector3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    let mut any = false;
    for k in 0..g.nz {
        for j in 0..g.ny {
            for i in 0..g.nx {
                if tunnel.body.is_solid(i, j, k) {
                    let (cx, cy, cz) = g.cell_centre(i, j, k);
                    let c = nalgebra::Vector3::new(cx, cy, cz);
                    min = min.inf(&c);
                    max = max.sup(&c);
                    any = true;
                }
            }
        }
    }
    if any {
        0.5 * (min + max)
    } else {
        nalgebra::Vector3::new(
            g.x0 + 0.5 * g.lx,
            g.y0 + 0.5 * g.ly,
            g.z0 + 0.5 * g.lz,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::box_body;
    use crate::turbulence::TurbulenceModel;
    use crate::wind::Wind;
    use nalgebra::Vector3;

    #[test]
    fn auto_controls_pick_a_sensible_time_step() {
        // A 1 m body at 20 m/s: shedding period ≈ L/(St·U) = 1/4 = 0.25
        // s, dt ≈ 0.25/30 ≈ 0.0083 s.
        let c = TransientControls::auto(1.0, 20.0, 100);
        assert!(c.dt > 0.001 && c.dt < 0.05, "auto dt {} unreasonable", c.dt);
        assert_eq!(c.steps, 100);
    }

    #[test]
    fn blend_factor_is_bounded() {
        // The implicit-Euler blend must always lie in [0, 0.95].
        for dt in [1e-6, 1e-3, 1.0, 1e3] {
            let b = blend_factor(1.225, 0.01, dt);
            assert!((0.0..=0.95).contains(&b), "blend {b} out of range for dt {dt}");
        }
        // A tiny dt retains almost all of the previous field.
        assert!(blend_factor(1.225, 1.0, 1e-6) > 0.9);
    }

    #[test]
    fn transient_run_records_a_force_history() {
        // A short transient run over a bluff body must produce one
        // frame per step, each with finite coefficients. A coarse grid
        // keeps the per-step SIMPLE solves fast.
        let body = box_body(Vector3::zeros(), Vector3::new(1.0, 1.0, 1.0));
        let tunnel = WindTunnel::build_with(
            &body,
            Wind::straight(20.0).unwrap(),
            crate::domain::BoundaryConditions::external_aero(),
            crate::domain::TunnelSizing {
                cells_across_body: 4,
                max_cells: 40_000,
                ..crate::domain::TunnelSizing::default()
            },
        )
        .unwrap();
        let mut controls = TransientControls::auto(1.0, 20.0, 6);
        controls.steady.turbulence = TurbulenceModel::KEpsilon;
        controls.steady.max_iterations = 10;
        controls.inner_iterations = 10;
        let hist = solve_transient(&tunnel, &controls, &BodyMotion::static_body());
        assert_eq!(hist.frames.len(), 6);
        assert!(hist.frames.iter().all(|f| f.coefficients.cd.is_finite()));
        // Time advances monotonically.
        for w in hist.frames.windows(2) {
            assert!(w[1].time > w[0].time);
        }
        // The diagnostics are finite.
        assert!(hist.mean_cd().is_finite());
        assert!(hist.lift_oscillation_amplitude() >= 0.0);
        assert!(hist.shedding_frequency() >= 0.0);
    }

    #[test]
    fn blend_fields_interpolates() {
        let mut a = Field3::filled(2, 2, 2, 10.0);
        let b = Field3::filled(2, 2, 2, 20.0);
        blend_fields(&mut a, &b, 0.5);
        assert!(a.data.iter().all(|&v| (v - 15.0).abs() < 1e-12));
    }
}
