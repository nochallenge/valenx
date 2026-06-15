//! Fixed-step RK4 integration of the target-cell-limited model.
//!
//! [`simulate`] advances an initial [`State`] under a [`Parameters`] set
//! from `t = 0` to `t_end` with the classic explicit fourth-order
//! Runge-Kutta scheme and a fixed step `dt`, returning a [`Trajectory`]
//! of paired time / state samples.
//!
//! Because the populations of this model are physical *amounts*, each
//! accepted RK4 step is clamped onto the non-negative orthant
//! (`max(0, ·)`). The continuous flow keeps the state in the positive
//! cone — `dT/dt = -beta*T*V` vanishes as `T -> 0`, and likewise for the
//! other compartments — but a finite step can carry a fast-decaying
//! component a hair below zero; the clamp keeps the discrete trajectory
//! inside the model's valid domain without changing the dynamics
//! anywhere the continuous solution is positive. The excursion clamped
//! away is within the method's truncation error.

use serde::{Deserialize, Serialize};

use crate::error::{ImmunoError, Result};
use crate::model::{Parameters, State};

/// A solved trajectory: parallel time and state samples.
///
/// `times[k]` is the simulation time of `states[k]`. The first sample is
/// always the initial condition at `t = 0`; the last is at (or just
/// past) `t_end`.
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub struct Trajectory {
    /// Sample times, strictly increasing.
    pub times: Vec<f64>,
    /// State vectors, one per time (`states.len() == times.len()`).
    pub states: Vec<State>,
}

impl Trajectory {
    /// Number of stored samples.
    pub fn len(&self) -> usize {
        self.times.len()
    }

    /// Whether the trajectory is empty.
    pub fn is_empty(&self) -> bool {
        self.times.is_empty()
    }

    /// The final state, if any.
    pub fn final_state(&self) -> Option<&State> {
        self.states.last()
    }

    /// The viral-load (`V`) time series across every sample.
    pub fn viral_load(&self) -> Vec<f64> {
        self.states.iter().map(|s| s.virus).collect()
    }

    /// The peak free-virus load and the sample at which it occurs.
    ///
    /// Returns `(index, time, peak_v)` for the sample with the largest
    /// `V`, or `None` for an empty trajectory. The first maximal sample
    /// wins on ties.
    pub fn peak_viral_load(&self) -> Option<(usize, f64, f64)> {
        let mut best: Option<(usize, f64, f64)> = None;
        for (i, s) in self.states.iter().enumerate() {
            match best {
                Some((_, _, bv)) if s.virus <= bv => {}
                _ => best = Some((i, self.times[i], s.virus)),
            }
        }
        best
    }

    /// Whether every sampled state is non-negative within `tol`.
    pub fn all_non_negative(&self, tol: f64) -> bool {
        self.states.iter().all(|s| s.is_non_negative(tol))
    }
}

/// One explicit RK4 step of size `h` from state `y` under `params`.
///
/// Returns the (un-clamped) `y(t + h)`; the model is autonomous, so no
/// explicit time argument is needed. The four stage derivatives are the
/// standard `k1..k4` of the classic Runge-Kutta method.
pub fn rk4_step(params: &Parameters, y: &State, h: f64) -> State {
    let k1 = y.derivative(params);
    let y2 = axpy(y, &k1, 0.5 * h);
    let k2 = y2.derivative(params);
    let y3 = axpy(y, &k2, 0.5 * h);
    let k3 = y3.derivative(params);
    let y4 = axpy(y, &k3, h);
    let k4 = y4.derivative(params);
    State {
        target: y.target + (h / 6.0) * (k1.target + 2.0 * k2.target + 2.0 * k3.target + k4.target),
        infected: y.infected
            + (h / 6.0) * (k1.infected + 2.0 * k2.infected + 2.0 * k3.infected + k4.infected),
        virus: y.virus + (h / 6.0) * (k1.virus + 2.0 * k2.virus + 2.0 * k3.virus + k4.virus),
    }
}

/// `base + scale * delta`, component-wise (an RK4 stage probe). The
/// result is an intermediate, so it is deliberately *not* validated for
/// non-negativity.
fn axpy(base: &State, delta: &State, scale: f64) -> State {
    State {
        target: base.target + scale * delta.target,
        infected: base.infected + scale * delta.infected,
        virus: base.virus + scale * delta.virus,
    }
}

/// Clamp a state onto the non-negative orthant.
fn clamp_non_negative(s: &State) -> State {
    State {
        target: s.target.max(0.0),
        infected: s.infected.max(0.0),
        virus: s.virus.max(0.0),
    }
}

/// Integrate the model from `t = 0` to `t_end` with fixed step `dt`.
///
/// A sample is stored every `n_out`-th step (and always the start point
/// and the final point). `n_out` is treated as at least `1`.
///
/// # Errors
///
/// - [`ImmunoError::Invalid`] if `dt <= 0`, `t_end <= 0`, or either is
///   non-finite.
/// - [`ImmunoError::NotFinite`] if the explicit scheme diverges to a
///   non-finite state (which happens only for a wildly over-large `dt`);
///   the error carries the offending step index and time.
pub fn simulate(
    params: &Parameters,
    y0: &State,
    t_end: f64,
    dt: f64,
    n_out: usize,
) -> Result<Trajectory> {
    if !dt.is_finite() || dt <= 0.0 {
        return Err(ImmunoError::invalid(
            "dt",
            "step must be positive and finite",
        ));
    }
    if !t_end.is_finite() || t_end <= 0.0 {
        return Err(ImmunoError::invalid(
            "t_end",
            "end time must be positive and finite",
        ));
    }
    let n_out = n_out.max(1);
    let steps = (t_end / dt).ceil() as usize;

    let mut traj = Trajectory {
        times: vec![0.0],
        states: vec![*y0],
    };
    let mut t = 0.0;
    let mut y = *y0;

    for step in 0..steps {
        let h = dt.min(t_end - t);
        let next = clamp_non_negative(&rk4_step(params, &y, h));
        if !next.is_finite() {
            return Err(ImmunoError::not_finite(step + 1, t + h));
        }
        y = next;
        t += h;
        if (step + 1) % n_out == 0 || step + 1 == steps {
            traj.times.push(t);
            traj.states.push(y);
        }
    }
    Ok(traj)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A representative acute-infection parameter set with `R0 > 1`,
    /// scaled so the populations stay in a numerically comfortable
    /// range. T0 = 1e5 target cells, one infected cell to seed.
    fn acute(c: f64) -> (Parameters, State) {
        // beta=2e-5, delta=1.0, p=50, with T0=1e5:
        // R0 = beta*T0*p/(delta*c) = (2e-5 * 1e5 * 50)/(1*c) = 100/c.
        let params = Parameters::new(2e-5, 1.0, 50.0, c).unwrap();
        let y0 = State::new(1e5, 1.0, 0.0).unwrap();
        (params, y0)
    }

    #[test]
    fn pure_virion_clearance_matches_analytic() {
        // With no infected cells and no infection (I0=0, T0=0 so no new
        // infection), dV/dt = -c V, an exact exponential decay.
        // V(t) = V0 exp(-c t).
        let params = Parameters::new(0.0, 1.0, 50.0, 3.0).unwrap();
        let y0 = State::new(0.0, 0.0, 1000.0).unwrap();
        let traj = simulate(&params, &y0, 2.0, 1e-3, 1).unwrap();
        let v = traj.final_state().unwrap().virus;
        let expect = 1000.0 * (-3.0_f64 * 2.0).exp();
        assert!((v - expect).abs() < 1e-3, "got {v}, want {expect}");
    }

    #[test]
    fn infected_cell_decay_matches_analytic() {
        // I0>0 but V0=0 and p=0 so V stays 0 and no infection occurs;
        // dI/dt = -delta I, exact exponential decay.
        // I(t) = I0 exp(-delta t).
        let params = Parameters::new(0.0, 0.7, 0.0, 1.0).unwrap();
        let y0 = State::new(0.0, 500.0, 0.0).unwrap();
        let traj = simulate(&params, &y0, 3.0, 1e-3, 1).unwrap();
        let i = traj.final_state().unwrap().infected;
        let expect = 500.0 * (-0.7_f64 * 3.0).exp();
        assert!((i - expect).abs() < 1e-3, "got {i}, want {expect}");
    }

    #[test]
    fn populations_stay_non_negative() {
        let (params, y0) = acute(5.0);
        let traj = simulate(&params, &y0, 20.0, 1e-3, 10).unwrap();
        // Every sampled population must be non-negative.
        assert!(traj.all_non_negative(0.0), "a population went negative");
    }

    #[test]
    fn target_cells_monotonically_deplete() {
        // dT/dt = -beta*T*V <= 0 always, so the target-cell series must
        // be non-increasing across the whole trajectory.
        let (params, y0) = acute(5.0);
        let traj = simulate(&params, &y0, 25.0, 1e-3, 1).unwrap();
        let mut prev = traj.states[0].target;
        for s in &traj.states {
            assert!(
                s.target <= prev + 1e-6,
                "target cells increased: {prev} -> {}",
                s.target
            );
            prev = s.target;
        }
    }

    #[test]
    fn viral_load_peaks_then_declines() {
        // The hallmark of the target-cell-limited model: V rises to a
        // peak, then falls as target cells deplete.
        let (params, y0) = acute(5.0);
        let traj = simulate(&params, &y0, 30.0, 1e-3, 5).unwrap();
        let (peak_idx, _peak_t, peak_v) = traj.peak_viral_load().unwrap();

        // The peak is interior (not at the first or last sample): it
        // genuinely rose and then fell.
        assert!(peak_idx > 0, "peak at the very first sample (no rise)");
        assert!(
            peak_idx < traj.len() - 1,
            "peak at the last sample (V never declined)"
        );

        // The peak is well above the seed value (a real outbreak).
        assert!(peak_v > 1.0, "viral load never grew: peak = {peak_v}");

        // The final load is below the peak (genuine decline).
        let final_v = traj.final_state().unwrap().virus;
        assert!(
            final_v < peak_v,
            "viral load did not decline: final {final_v} >= peak {peak_v}"
        );
    }

    #[test]
    fn target_cells_deplete_by_the_peak() {
        // By the time V peaks, a substantial fraction of target cells
        // should be consumed — depletion is what turns the epidemic
        // over. Confirm T at the peak is well below T0.
        let (params, y0) = acute(5.0);
        let traj = simulate(&params, &y0, 30.0, 1e-3, 1).unwrap();
        let (peak_idx, _, _) = traj.peak_viral_load().unwrap();
        let t_at_peak = traj.states[peak_idx].target;
        assert!(
            t_at_peak < 0.9 * y0.target,
            "target cells not depleted at the peak: {t_at_peak} vs T0 {}",
            y0.target
        );
    }

    #[test]
    fn higher_clearance_lowers_peak_viral_load() {
        // Increasing the virion clearance rate c lowers R0 and therefore
        // the peak viral load. Compare two runs differing only in c.
        let (p_low_c, y0) = acute(5.0); // c = 5  -> R0 = 20
        let (p_high_c, _) = acute(20.0); // c = 20 -> R0 = 5

        let traj_low = simulate(&p_low_c, &y0, 30.0, 1e-3, 1).unwrap();
        let traj_high = simulate(&p_high_c, &y0, 30.0, 1e-3, 1).unwrap();

        let (_, _, peak_low) = traj_low.peak_viral_load().unwrap();
        let (_, _, peak_high) = traj_high.peak_viral_load().unwrap();

        assert!(
            peak_high < peak_low,
            "higher clearance did not lower the peak: c=20 peak {peak_high} >= c=5 peak {peak_low}"
        );
    }

    #[test]
    fn subthreshold_r0_gives_no_outbreak() {
        // With R0 < 1 (very high clearance) the seeded infection fades
        // without an interior peak — V is largest at the start and only
        // declines. R0 = 100/c; pick c = 200 -> R0 = 0.5.
        let (params, y0) = acute(200.0);
        assert!(params.r0(y0.target).unwrap() < 1.0);
        let traj = simulate(&params, &y0, 30.0, 1e-3, 5).unwrap();
        // Free virus starts at 0, gets a tiny bump from the single
        // infected cell, but never mounts a real outbreak: the peak load
        // stays tiny compared with an R0=20 run.
        let (_, _, peak_sub) = traj.peak_viral_load().unwrap();
        let (params_hi, y0_hi) = acute(5.0);
        let traj_hi = simulate(&params_hi, &y0_hi, 30.0, 1e-3, 5).unwrap();
        let (_, _, peak_hi) = traj_hi.peak_viral_load().unwrap();
        assert!(
            peak_sub < peak_hi,
            "sub-threshold peak {peak_sub} not below outbreak peak {peak_hi}"
        );
    }

    #[test]
    fn rk4_step_reduces_to_decay() {
        // A single RK4 step of pure virion decay (beta=p=0) should match
        // exp(-c h) closely for a small h.
        let params = Parameters::new(0.0, 1.0, 0.0, 2.0).unwrap();
        let y = State::new(0.0, 0.0, 100.0).unwrap();
        let stepped = rk4_step(&params, &y, 0.01);
        let expect = 100.0 * (-2.0_f64 * 0.01).exp();
        // RK4 local error for this stiff-free decay is ~h^5; well under 1e-6.
        assert!(
            (stepped.virus - expect).abs() < 1e-6,
            "got {}, want {expect}",
            stepped.virus
        );
    }

    #[test]
    fn rejects_bad_step_and_horizon() {
        let (params, y0) = acute(5.0);
        assert!(simulate(&params, &y0, 10.0, 0.0, 1).is_err());
        assert!(simulate(&params, &y0, 10.0, -1.0, 1).is_err());
        assert!(simulate(&params, &y0, 0.0, 0.1, 1).is_err());
        assert!(simulate(&params, &y0, -5.0, 0.1, 1).is_err());
        assert!(simulate(&params, &y0, f64::NAN, 0.1, 1).is_err());
    }

    #[test]
    fn trajectory_starts_at_initial_condition() {
        let (params, y0) = acute(5.0);
        let traj = simulate(&params, &y0, 5.0, 0.01, 7).unwrap();
        assert!(!traj.is_empty());
        let first = &traj.states[0];
        assert!((first.target - y0.target).abs() < 1e-12);
        assert!((first.infected - y0.infected).abs() < 1e-12);
        assert!((first.virus - y0.virus).abs() < 1e-12);
        assert!((traj.times[0]).abs() < 1e-12);
        // Final time reaches the horizon.
        assert!((traj.times.last().unwrap() - 5.0).abs() < 1e-9);
    }

    #[test]
    fn peak_viral_load_picks_the_maximum() {
        // Construct a hand-made trajectory and confirm the picker.
        let traj = Trajectory {
            times: vec![0.0, 1.0, 2.0, 3.0],
            states: vec![
                State::new(0.0, 0.0, 1.0).unwrap(),
                State::new(0.0, 0.0, 5.0).unwrap(),
                State::new(0.0, 0.0, 9.0).unwrap(),
                State::new(0.0, 0.0, 4.0).unwrap(),
            ],
        };
        let (idx, t, v) = traj.peak_viral_load().unwrap();
        assert_eq!(idx, 2);
        assert!((t - 2.0).abs() < 1e-12);
        assert!((v - 9.0).abs() < 1e-12);
        assert!(Trajectory::default().peak_viral_load().is_none());
    }
}
