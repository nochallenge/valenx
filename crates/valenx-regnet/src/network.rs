//! The gene-regulatory-network model and its RK4 integrator.
//!
//! A [`GeneRegulatoryNetwork`] is `N` genes, each described by a [`Gene`]:
//! a basal **production** rate, a first-order **degradation** rate, a
//! **basal** (regulator-independent) fraction, and a list of
//! [`Regulator`] inputs. The state is the vector of expression levels
//! `x = (x_0, …, x_{N-1})`, stored as a plain `Vec<f64>`.
//!
//! ## Governing equation
//!
//! For each gene `i`,
//!
//! ```text
//! dx_i/dt = production_i * regulation_i(x) - degradation_i * x_i
//! ```
//!
//! The dimensionless **regulation term** `regulation_i(x) ∈ [0, 1]` is
//! assembled from the gene's regulators. Each regulator contributes a
//! Hill factor — [`hill_activate`] for an
//! activator, [`hill_repress`] for a
//! repressor, evaluated at the regulating gene's current level — and the
//! factors are **combined multiplicatively** (AND-like logic: every
//! activator must be present and every repressor absent for full
//! transcription). A gene with no regulators uses its `basal` fraction
//! directly (constitutive expression). When a gene *has* regulators, the
//! `basal` term is added as a leak floor and the result is capped at `1`,
//! so a fully repressed promoter still produces at the basal rate.
//!
//! ## Integration
//!
//! [`GeneRegulatoryNetwork::rk4_step`] advances the state by one fixed
//! step `dt` with the classical 4th-order Runge-Kutta method, and
//! [`GeneRegulatoryNetwork::simulate`] runs `steps` of it, returning the
//! full [`Trajectory`] (the initial state plus one row per step).

use crate::error::{RegnetError, Result};
use crate::hill::{hill_activate, hill_repress};
use serde::{Deserialize, Serialize};

/// The sign of a regulatory interaction.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum RegulatorKind {
    /// The regulator is an **activator**: its presence raises the
    /// target's transcription via [`hill_activate`].
    Activate,
    /// The regulator is a **repressor**: its presence lowers the
    /// target's transcription via [`hill_repress`].
    Repress,
}

/// One regulatory input to a gene: a directed, signed Hill interaction
/// from a source gene (`regulator`) onto the gene that owns this struct.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Regulator {
    /// Index of the **source** gene whose expression level drives this
    /// interaction. Must be a valid gene index of the owning network.
    pub regulator: usize,
    /// Whether the source activates or represses the target.
    pub kind: RegulatorKind,
    /// Half-saturation threshold `k` of the Hill factor (`k > 0`).
    pub k: f64,
    /// Hill (cooperativity) coefficient `n` (`n > 0`).
    pub n: f64,
}

impl Regulator {
    /// Construct an **activating** regulator from gene `regulator` with
    /// threshold `k` and Hill coefficient `n`.
    #[must_use]
    pub fn activate(regulator: usize, k: f64, n: f64) -> Self {
        Self {
            regulator,
            kind: RegulatorKind::Activate,
            k,
            n,
        }
    }

    /// Construct a **repressing** regulator from gene `regulator` with
    /// threshold `k` and Hill coefficient `n`.
    #[must_use]
    pub fn repress(regulator: usize, k: f64, n: f64) -> Self {
        Self {
            regulator,
            kind: RegulatorKind::Repress,
            k,
            n,
        }
    }

    /// Evaluate this regulator's Hill factor at network state `state`.
    /// The factor lies in `[0, 1]`.
    #[must_use]
    fn factor(&self, state: &[f64]) -> f64 {
        let x = state[self.regulator];
        match self.kind {
            RegulatorKind::Activate => hill_activate(x, self.k, self.n),
            RegulatorKind::Repress => hill_repress(x, self.k, self.n),
        }
    }
}

/// A single gene: its kinetic rates plus its incoming regulatory inputs.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Gene {
    /// Maximal production (transcription) rate when fully active.
    pub production: f64,
    /// First-order degradation (decay) rate constant.
    pub degradation: f64,
    /// Basal expression fraction in `[0, 1]`: the leak floor of the
    /// regulation term. A gene with **no** regulators is expressed at
    /// exactly this fraction (constitutive); a regulated gene adds this
    /// as a floor so a fully repressed promoter is not perfectly silent.
    pub basal: f64,
    /// Incoming regulatory inputs (may be empty for a constitutive gene).
    pub regulators: Vec<Regulator>,
}

impl Gene {
    /// Build a constitutive gene (no regulators) expressed at the given
    /// `basal` fraction.
    #[must_use]
    pub fn constitutive(production: f64, degradation: f64, basal: f64) -> Self {
        Self {
            production,
            degradation,
            basal,
            regulators: Vec::new(),
        }
    }

    /// Build a regulated gene from a production rate, degradation rate,
    /// basal leak and a set of regulators.
    #[must_use]
    pub fn regulated(
        production: f64,
        degradation: f64,
        basal: f64,
        regulators: Vec<Regulator>,
    ) -> Self {
        Self {
            production,
            degradation,
            basal,
            regulators,
        }
    }

    /// The dimensionless regulation term `regulation_i(x) ∈ [0, 1]` for
    /// this gene at network state `state`.
    ///
    /// With no regulators this is the constitutive `basal` fraction. With
    /// regulators it is the product of their Hill factors, lifted by the
    /// `basal` leak floor and capped at `1`.
    #[must_use]
    fn regulation(&self, state: &[f64]) -> f64 {
        if self.regulators.is_empty() {
            return self.basal.clamp(0.0, 1.0);
        }
        let mut product = 1.0;
        for r in &self.regulators {
            product *= r.factor(state);
        }
        (self.basal + product).clamp(0.0, 1.0)
    }
}

/// A gene-regulatory network: an ordered collection of [`Gene`]s.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GeneRegulatoryNetwork {
    /// The genes, indexed `0..len`. A [`Regulator::regulator`] field
    /// indexes into this vector.
    pub genes: Vec<Gene>,
}

impl GeneRegulatoryNetwork {
    /// Build a network from a vector of genes, validating every kinetic
    /// rate, Hill parameter and regulator index.
    ///
    /// # Errors
    ///
    /// - [`RegnetError::InvalidRate`] if any production or degradation
    ///   rate is negative (or non-finite).
    /// - [`RegnetError::InvalidHill`] if any regulator's `k` or `n` is not
    ///   strictly positive.
    /// - [`RegnetError::GeneIndexOutOfRange`] if any regulator references
    ///   a gene index outside `0..genes.len()`.
    pub fn new(genes: Vec<Gene>) -> Result<Self> {
        let count = genes.len();
        for (i, g) in genes.iter().enumerate() {
            if !(g.production.is_finite() && g.production >= 0.0) {
                return Err(RegnetError::InvalidRate {
                    what: "production",
                    gene: i,
                    value: g.production,
                });
            }
            if !(g.degradation.is_finite() && g.degradation >= 0.0) {
                return Err(RegnetError::InvalidRate {
                    what: "degradation",
                    gene: i,
                    value: g.degradation,
                });
            }
            for r in &g.regulators {
                if !(r.k.is_finite() && r.k > 0.0) {
                    return Err(RegnetError::InvalidHill {
                        what: "threshold k",
                        value: r.k,
                    });
                }
                if !(r.n.is_finite() && r.n > 0.0) {
                    return Err(RegnetError::InvalidHill {
                        what: "coefficient n",
                        value: r.n,
                    });
                }
                if r.regulator >= count {
                    return Err(RegnetError::GeneIndexOutOfRange {
                        index: r.regulator,
                        count,
                    });
                }
            }
        }
        Ok(Self { genes })
    }

    /// The number of genes `N` in the network.
    #[must_use]
    pub fn len(&self) -> usize {
        self.genes.len()
    }

    /// Whether the network has no genes.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.genes.is_empty()
    }

    /// Evaluate the right-hand side `dx/dt` of the ODE system at `state`,
    /// writing the `N` derivatives into `deriv`.
    ///
    /// `deriv` must already have length `N` (the caller owns its
    /// allocation so the hot RK4 loop does not allocate).
    fn rhs_into(&self, state: &[f64], deriv: &mut [f64]) {
        for (i, g) in self.genes.iter().enumerate() {
            let reg = g.regulation(state);
            deriv[i] = g.production * reg - g.degradation * state[i];
        }
    }

    /// Evaluate the right-hand side `dx/dt` at `state` and return it as a
    /// fresh `Vec<f64>` of length `N`.
    ///
    /// # Errors
    ///
    /// Returns [`RegnetError::DimensionMismatch`] if `state.len() != N`.
    pub fn derivative(&self, state: &[f64]) -> Result<Vec<f64>> {
        self.check_state(state)?;
        let mut d = vec![0.0; self.genes.len()];
        self.rhs_into(state, &mut d);
        Ok(d)
    }

    /// Validate that `state` has exactly `N` entries.
    fn check_state(&self, state: &[f64]) -> Result<()> {
        if state.len() != self.genes.len() {
            return Err(RegnetError::DimensionMismatch {
                what: "state",
                expected: self.genes.len(),
                found: state.len(),
            });
        }
        Ok(())
    }

    /// Advance `state` by one classical RK4 step of size `dt`, returning
    /// the new state.
    ///
    /// The classical fourth-order Runge-Kutta update for `dx/dt = f(x)`:
    ///
    /// ```text
    /// k1 = f(x)
    /// k2 = f(x + dt/2 * k1)
    /// k3 = f(x + dt/2 * k2)
    /// k4 = f(x + dt   * k3)
    /// x' = x + dt/6 * (k1 + 2 k2 + 2 k3 + k4)
    /// ```
    ///
    /// # Errors
    ///
    /// - [`RegnetError::DimensionMismatch`] if `state.len() != N`.
    /// - [`RegnetError::InvalidStep`] if `dt` is not strictly positive
    ///   and finite.
    pub fn rk4_step(&self, state: &[f64], dt: f64) -> Result<Vec<f64>> {
        self.check_state(state)?;
        if !(dt.is_finite() && dt > 0.0) {
            return Err(RegnetError::invalid_step(format!(
                "dt must be a finite value > 0, got {dt}"
            )));
        }
        Ok(self.rk4_step_unchecked(state, dt))
    }

    /// One RK4 step without revalidating `state` length or `dt`. Used by
    /// the inner [`simulate`](Self::simulate) loop after a single
    /// up-front validation.
    fn rk4_step_unchecked(&self, state: &[f64], dt: f64) -> Vec<f64> {
        let n = state.len();
        let mut k1 = vec![0.0; n];
        let mut k2 = vec![0.0; n];
        let mut k3 = vec![0.0; n];
        let mut k4 = vec![0.0; n];
        let mut tmp = vec![0.0; n];

        self.rhs_into(state, &mut k1);

        for i in 0..n {
            tmp[i] = state[i] + 0.5 * dt * k1[i];
        }
        self.rhs_into(&tmp, &mut k2);

        for i in 0..n {
            tmp[i] = state[i] + 0.5 * dt * k2[i];
        }
        self.rhs_into(&tmp, &mut k3);

        for i in 0..n {
            tmp[i] = state[i] + dt * k3[i];
        }
        self.rhs_into(&tmp, &mut k4);

        let mut out = vec![0.0; n];
        for i in 0..n {
            out[i] = state[i] + (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
        out
    }

    /// Integrate the system from `initial` for `steps` fixed RK4 steps of
    /// size `dt`, returning the full [`Trajectory`].
    ///
    /// The returned trajectory has `steps + 1` rows: row `0` is `initial`
    /// and row `j` is the state after `j` steps (at time `j * dt`).
    ///
    /// # Errors
    ///
    /// - [`RegnetError::DimensionMismatch`] if `initial.len() != N`.
    /// - [`RegnetError::InvalidStep`] if `dt <= 0` (or non-finite) or
    ///   `steps == 0`.
    pub fn simulate(&self, initial: &[f64], dt: f64, steps: usize) -> Result<Trajectory> {
        self.check_state(initial)?;
        if !(dt.is_finite() && dt > 0.0) {
            return Err(RegnetError::invalid_step(format!(
                "dt must be a finite value > 0, got {dt}"
            )));
        }
        if steps == 0 {
            return Err(RegnetError::invalid_step("steps must be >= 1, got 0"));
        }

        let mut states = Vec::with_capacity(steps + 1);
        states.push(initial.to_vec());
        let mut current = initial.to_vec();
        for _ in 0..steps {
            current = self.rk4_step_unchecked(&current, dt);
            states.push(current.clone());
        }
        Ok(Trajectory { dt, states })
    }
}

/// The output of [`GeneRegulatoryNetwork::simulate`]: a fixed-step time
/// series of network states.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Trajectory {
    /// The fixed timestep used to produce this trajectory.
    pub dt: f64,
    /// One row per recorded time point. `states[0]` is the initial
    /// condition; `states[j]` is the state at time `j * dt`. Every row has
    /// the same length `N` (the gene count).
    pub states: Vec<Vec<f64>>,
}

impl Trajectory {
    /// The number of recorded time points (`steps + 1`).
    #[must_use]
    pub fn len(&self) -> usize {
        self.states.len()
    }

    /// Whether the trajectory has no recorded states.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.states.is_empty()
    }

    /// The simulated time of row `index`, i.e. `index * dt`.
    #[must_use]
    pub fn time_at(&self, index: usize) -> f64 {
        index as f64 * self.dt
    }

    /// The final recorded state, or `None` if the trajectory is empty.
    #[must_use]
    pub fn final_state(&self) -> Option<&[f64]> {
        self.states.last().map(Vec::as_slice)
    }

    /// Extract the time series of a single gene `gene` across all recorded
    /// time points.
    ///
    /// # Errors
    ///
    /// Returns [`RegnetError::GeneIndexOutOfRange`] if `gene` is not a
    /// valid column index of this trajectory.
    pub fn series(&self, gene: usize) -> Result<Vec<f64>> {
        let width = self.states.first().map_or(0, Vec::len);
        if gene >= width {
            return Err(RegnetError::GeneIndexOutOfRange {
                index: gene,
                count: width,
            });
        }
        Ok(self.states.iter().map(|row| row[gene]).collect())
    }

    /// Count the **strict interior local maxima** of gene `gene`'s time
    /// series: indices `j` (with `0 < j < len-1`) where the value is
    /// strictly greater than both neighbours.
    ///
    /// This is the oscillation diagnostic used by the repressilator test —
    /// a sustained oscillator produces several such peaks over a long run,
    /// while a monotone approach to steady state produces none.
    ///
    /// # Errors
    ///
    /// Returns [`RegnetError::GeneIndexOutOfRange`] if `gene` is not a
    /// valid column index.
    pub fn count_local_maxima(&self, gene: usize) -> Result<usize> {
        let s = self.series(gene)?;
        if s.len() < 3 {
            return Ok(0);
        }
        let mut peaks = 0;
        for j in 1..s.len() - 1 {
            if s[j] > s[j - 1] && s[j] > s[j + 1] {
                peaks += 1;
            }
        }
        Ok(peaks)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    /// A one-gene "network" implementing pure decay `dx/dt = -decay * x`
    /// (zero production), used to check RK4 against the analytic
    /// `x(t) = x0 * exp(-decay * t)`.
    fn decay_net(decay: f64) -> GeneRegulatoryNetwork {
        GeneRegulatoryNetwork::new(vec![Gene::constitutive(0.0, decay, 0.0)]).unwrap()
    }

    #[test]
    fn rk4_matches_exponential_decay() {
        // dx/dt = -x, x0 = 1  =>  x(t) = exp(-t).
        let net = decay_net(1.0);
        let dt = 0.01;
        let steps = 500; // t_final = 5.0
        let traj = net.simulate(&[1.0], dt, steps).unwrap();
        for j in 0..=steps {
            let t = j as f64 * dt;
            let got = traj.states[j][0];
            let want = (-t).exp();
            assert!(
                (got - want).abs() < 1e-3,
                "t={t}: got {got}, want {want}, |err|={}",
                (got - want).abs()
            );
        }
    }

    #[test]
    fn rk4_decay_rate_two() {
        // dx/dt = -2x, x0 = 3  =>  x(t) = 3 exp(-2t).
        let net = decay_net(2.0);
        let dt = 0.005;
        let steps = 400; // t_final = 2.0
        let traj = net.simulate(&[3.0], dt, steps).unwrap();
        let t: f64 = 2.0;
        let got = *traj.final_state().unwrap().first().unwrap();
        let want = 3.0 * (-2.0 * t).exp();
        assert!((got - want).abs() < 1e-3, "got {got}, want {want}");
    }

    #[test]
    fn trajectory_shape_and_times() {
        let net = decay_net(1.0);
        let traj = net.simulate(&[1.0], 0.1, 10).unwrap();
        assert_eq!(traj.len(), 11); // steps + 1
        assert!((traj.time_at(0) - 0.0).abs() < EPS);
        assert!((traj.time_at(10) - 1.0).abs() < EPS);
        assert_eq!(traj.states[0], vec![1.0]);
    }

    #[test]
    fn derivative_law_is_production_times_regulation_minus_decay() {
        // Single constitutive gene: dx/dt = prod*basal - deg*x.
        // prod=2, basal=1 (so regulation=1), deg=0.5, x=4 => 2 - 2 = 0.
        let net = GeneRegulatoryNetwork::new(vec![Gene::constitutive(2.0, 0.5, 1.0)]).unwrap();
        let d = net.derivative(&[4.0]).unwrap();
        assert!((d[0] - 0.0).abs() < EPS, "got {}", d[0]);
        // At x=0: dx/dt = 2*1 - 0 = 2.
        let d0 = net.derivative(&[0.0]).unwrap();
        assert!((d0[0] - 2.0).abs() < EPS, "got {}", d0[0]);
    }

    #[test]
    fn constitutive_gene_reaches_steady_state_prod_basal_over_deg() {
        // Steady state of dx/dt = p*b - d*x is x* = p*b/d.
        // p=5, b=1, d=2 => x* = 2.5.
        let net = GeneRegulatoryNetwork::new(vec![Gene::constitutive(5.0, 2.0, 1.0)]).unwrap();
        let traj = net.simulate(&[0.0], 0.01, 2000).unwrap();
        let x = *traj.final_state().unwrap().first().unwrap();
        assert!((x - 2.5).abs() < 1e-3, "got {x}");
    }

    #[test]
    fn new_rejects_negative_rates() {
        assert!(GeneRegulatoryNetwork::new(vec![Gene::constitutive(-1.0, 1.0, 0.5)]).is_err());
        assert!(GeneRegulatoryNetwork::new(vec![Gene::constitutive(1.0, -1.0, 0.5)]).is_err());
    }

    #[test]
    fn new_rejects_out_of_range_regulator() {
        // Gene 0 references gene 5 in a 1-gene network.
        let g = Gene::regulated(1.0, 1.0, 0.0, vec![Regulator::repress(5, 1.0, 2.0)]);
        let err = GeneRegulatoryNetwork::new(vec![g]).unwrap_err();
        assert!(matches!(err, RegnetError::GeneIndexOutOfRange { .. }));
    }

    #[test]
    fn new_rejects_bad_hill_params_on_regulator() {
        let g = Gene::regulated(1.0, 1.0, 0.0, vec![Regulator::repress(0, 0.0, 2.0)]);
        assert!(GeneRegulatoryNetwork::new(vec![g]).is_err());
        let g2 = Gene::regulated(1.0, 1.0, 0.0, vec![Regulator::activate(0, 1.0, 0.0)]);
        assert!(GeneRegulatoryNetwork::new(vec![g2]).is_err());
    }

    #[test]
    fn rk4_step_rejects_bad_dt_and_dim() {
        let net = decay_net(1.0);
        assert!(net.rk4_step(&[1.0], 0.0).is_err());
        assert!(net.rk4_step(&[1.0], -0.1).is_err());
        assert!(net.rk4_step(&[1.0, 2.0], 0.1).is_err()); // wrong length
    }

    #[test]
    fn simulate_rejects_zero_steps_and_bad_dim() {
        let net = decay_net(1.0);
        assert!(net.simulate(&[1.0], 0.1, 0).is_err());
        assert!(net.simulate(&[1.0, 2.0], 0.1, 5).is_err());
    }

    #[test]
    fn series_and_local_maxima_on_a_known_signal() {
        // Build a trajectory by hand: a clear up-down-up-down sawtooth on
        // gene 0 has two interior local maxima.
        let traj = Trajectory {
            dt: 1.0,
            states: vec![
                vec![0.0],
                vec![2.0], // peak (>0 and >1)
                vec![1.0],
                vec![3.0], // peak (>1 and >0)
                vec![0.0],
            ],
        };
        assert_eq!(traj.count_local_maxima(0).unwrap(), 2);
        assert_eq!(traj.series(0).unwrap(), vec![0.0, 2.0, 1.0, 3.0, 0.0]);
        assert!(traj.count_local_maxima(1).is_err()); // no such column
    }

    #[test]
    fn monotone_decay_has_no_interior_maxima() {
        let net = decay_net(1.0);
        let traj = net.simulate(&[1.0], 0.05, 100).unwrap();
        assert_eq!(traj.count_local_maxima(0).unwrap(), 0);
    }

    #[test]
    fn regulation_term_is_complement_pair_for_repressor() {
        // A gene repressed by gene 0 with basal 0: regulation = hill_repress(x0).
        // At x0 = k, that's 0.5; production*0.5 is the instantaneous max rate term.
        let net = GeneRegulatoryNetwork::new(vec![
            Gene::constitutive(0.0, 0.0, 0.0), // gene 0: inert source we set by state
            Gene::regulated(4.0, 1.0, 0.0, vec![Regulator::repress(0, 2.0, 2.0)]),
        ])
        .unwrap();
        // state: gene0 = 2.0 (== k) -> hill_repress = 0.5; gene1 = 0.
        let d = net.derivative(&[2.0, 0.0]).unwrap();
        // dx1/dt = 4 * 0.5 - 1*0 = 2.0.
        assert!((d[1] - 2.0).abs() < EPS, "got {}", d[1]);
    }
}
