//! The M/M/1 single-server queue.
//!
//! In Kendall's notation `A/S/c`, **M/M/1** is the queue with
//! **M**arkovian (Poisson) arrivals at rate `lambda`, **M**arkovian
//! (exponential) service at rate `mu`, and **1** server, an infinite
//! waiting room, and first-come-first-served discipline. It is the
//! birth-death chain with constant birth rate `lambda` and constant
//! death rate `mu`.
//!
//! # Steady-state results
//!
//! Writing the **traffic intensity** (offered load / server
//! utilization) as `rho = lambda / mu`, the queue has a stationary
//! distribution iff `rho < 1`, and then:
//!
//! - stationary state probability `P(n) = (1 - rho) rho^n` for
//!   `n = 0, 1, 2, ...` (a geometric law);
//! - idle probability `P(0) = 1 - rho`;
//! - mean number in system `L = rho / (1 - rho)`;
//! - mean number waiting in queue `Lq = rho^2 / (1 - rho) = L - rho`;
//! - mean time in system `W = 1 / (mu - lambda)`;
//! - mean waiting time in queue `Wq = rho / (mu - lambda) = W - 1 / mu`.
//!
//! These satisfy **Little's law** `L = lambda W` (and `Lq = lambda Wq`)
//! exactly; [`Mm1::little_residual`] returns the numeric residual as a
//! cross-check. As `rho -> 1` from below, every mean
//! (`L`, `Lq`, `W`, `Wq`) diverges to `+inf`.
//!
//! # Honest scope
//!
//! This is the textbook closed-form steady-state M/M/1 model: a single
//! server, Poisson arrivals, exponential service, infinite queue,
//! work-conserving FCFS, and the long-run *average* behaviour. It is
//! **not** a discrete-event simulator, says nothing about transient or
//! finite-horizon behaviour, and is research/educational grade — not a
//! production capacity-planning tool.

use crate::error::{QueueingError, Result};
use serde::{Deserialize, Serialize};

/// A validated M/M/1 queue: Poisson arrivals at rate `lambda`,
/// exponential service at rate `mu`, one server.
///
/// Construct with [`Mm1::new`], which rejects non-positive or
/// non-finite rates and unstable loads (`rho = lambda / mu >= 1`). Once
/// built, every steady-state accessor is total (returns a finite
/// number) because stability is an invariant of the type.
///
/// Rates are dimensionless *per unit time* and must share the same time
/// unit (customers per second, jobs per hour, …); the time-based
/// outputs ([`W`](Self::w), [`Wq`](Self::wq)) come back in that same
/// unit.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Mm1 {
    /// Mean arrival rate `lambda` (customers per unit time), `> 0`.
    lambda: f64,
    /// Mean service rate `mu` (customers per unit time), `> 0` and
    /// strictly greater than `lambda`.
    mu: f64,
}

impl Mm1 {
    /// Build an M/M/1 queue from an arrival rate `lambda` and a service
    /// rate `mu` (both in the same per-unit-time units).
    ///
    /// # Errors
    ///
    /// - [`QueueingError::Invalid`] if either rate is non-finite or not
    ///   strictly positive.
    /// - [`QueueingError::Unstable`] if `rho = lambda / mu >= 1` (i.e.
    ///   `lambda >= mu`), since no finite steady state then exists.
    pub fn new(lambda: f64, mu: f64) -> Result<Self> {
        if !lambda.is_finite() {
            return Err(QueueingError::invalid(
                "lambda",
                format!("arrival rate must be finite, got {lambda}"),
            ));
        }
        if !mu.is_finite() {
            return Err(QueueingError::invalid(
                "mu",
                format!("service rate must be finite, got {mu}"),
            ));
        }
        if lambda <= 0.0 {
            return Err(QueueingError::invalid(
                "lambda",
                format!("arrival rate must be > 0, got {lambda}"),
            ));
        }
        if mu <= 0.0 {
            return Err(QueueingError::invalid(
                "mu",
                format!("service rate must be > 0, got {mu}"),
            ));
        }
        let rho = lambda / mu;
        if rho >= 1.0 {
            return Err(QueueingError::unstable(lambda, mu, rho));
        }
        Ok(Mm1 { lambda, mu })
    }

    /// The service rate `mu` a single server must sustain so the mean time
    /// in system (response time) [`W`](Self::w) equals `target_w` at an
    /// arrival rate `arrival_rate` — the capacity-sizing inverse of
    /// [`Mm1::w`].
    ///
    /// Inverting `W = 1 / (mu - lambda)` gives `mu = lambda + 1 / target_w`.
    /// Because `target_w > 0` the result always exceeds `lambda`, so a queue
    /// built with it is always stable and its [`w`](Self::w) reproduces
    /// `target_w`. A tighter (smaller) target response time demands a higher
    /// service rate.
    ///
    /// # Errors
    ///
    /// [`QueueingError::Invalid`] if `arrival_rate` or `target_w` is not
    /// finite and strictly positive.
    pub fn service_rate_for_mean_response_time(arrival_rate: f64, target_w: f64) -> Result<f64> {
        if !arrival_rate.is_finite() || arrival_rate <= 0.0 {
            return Err(QueueingError::invalid(
                "arrival_rate",
                format!("arrival rate must be finite and > 0, got {arrival_rate}"),
            ));
        }
        if !target_w.is_finite() || target_w <= 0.0 {
            return Err(QueueingError::invalid(
                "target_w",
                format!("target response time must be finite and > 0, got {target_w}"),
            ));
        }
        Ok(arrival_rate + 1.0 / target_w)
    }

    /// The arrival rate `lambda` a single server of rate `service_rate`
    /// can accept while holding the mean time in system
    /// [`W`](Self::w) at `target_w` — the load-sizing inverse of
    /// [`Mm1::w`], complementary to
    /// [`service_rate_for_mean_response_time`](Self::service_rate_for_mean_response_time).
    ///
    /// Inverting `W = 1 / (mu - lambda)` for the arrival rate gives
    /// `lambda = mu - 1 / target_w`. A queue built with it always
    /// reproduces `target_w` and is stable (`lambda < mu`). A more
    /// tolerant (larger) `target_w` admits a higher arrival rate.
    ///
    /// # Errors
    ///
    /// [`QueueingError::Invalid`] if `service_rate` or `target_w` is not
    /// finite and strictly positive, or if `target_w` is at or below the
    /// no-wait service time `1 / mu` (which would demand a non-positive
    /// arrival rate — no load can be served faster than the bare service
    /// time).
    pub fn arrival_rate_for_mean_response_time(service_rate: f64, target_w: f64) -> Result<f64> {
        if !service_rate.is_finite() || service_rate <= 0.0 {
            return Err(QueueingError::invalid(
                "service_rate",
                format!("service rate must be finite and > 0, got {service_rate}"),
            ));
        }
        if !target_w.is_finite() || target_w <= 0.0 {
            return Err(QueueingError::invalid(
                "target_w",
                format!("target response time must be finite and > 0, got {target_w}"),
            ));
        }
        let lambda = service_rate - 1.0 / target_w;
        if lambda <= 0.0 {
            return Err(QueueingError::invalid(
                "target_w",
                format!(
                    "target response time {target_w} is unreachable: it must exceed the \
                     no-wait service time 1/mu = {}",
                    1.0 / service_rate
                ),
            ));
        }
        Ok(lambda)
    }

    /// The arrival rate `lambda` this queue was built with.
    pub fn lambda(&self) -> f64 {
        self.lambda
    }

    /// The service rate `mu` this queue was built with.
    pub fn mu(&self) -> f64 {
        self.mu
    }

    /// Traffic intensity (offered load / server utilization)
    /// `rho = lambda / mu`. Always in `(0, 1)` for a constructed queue.
    ///
    /// `rho` is also the long-run fraction of time the single server is
    /// busy, so it doubles as the server utilization.
    pub fn rho(&self) -> f64 {
        self.lambda / self.mu
    }

    /// Server utilization — the long-run fraction of time the server is
    /// busy. For M/M/1 this equals [`rho`](Self::rho); provided as a
    /// named alias for callers thinking in utilization terms.
    pub fn utilization(&self) -> f64 {
        self.rho()
    }

    /// Idle probability `P(0) = 1 - rho` — the long-run fraction of time
    /// the system is empty (zero customers present).
    pub fn p0(&self) -> f64 {
        1.0 - self.rho()
    }

    /// Mean number of customers **in the system**
    /// `L = rho / (1 - rho)` (those waiting plus the one in service).
    pub fn l(&self) -> f64 {
        let rho = self.rho();
        rho / (1.0 - rho)
    }

    /// Mean number of customers **waiting in the queue**
    /// `Lq = rho^2 / (1 - rho)`, excluding any customer currently in
    /// service. Equivalently `Lq = L - rho`.
    pub fn lq(&self) -> f64 {
        let rho = self.rho();
        rho * rho / (1.0 - rho)
    }

    /// Mean total time a customer spends **in the system**
    /// `W = 1 / (mu - lambda)` (wait plus service), in the queue's time
    /// unit.
    pub fn w(&self) -> f64 {
        1.0 / (self.mu - self.lambda)
    }

    /// Mean time a customer spends **waiting in the queue**
    /// `Wq = rho / (mu - lambda)` before service begins. Equivalently
    /// `Wq = W - 1 / mu`.
    pub fn wq(&self) -> f64 {
        self.rho() / (self.mu - self.lambda)
    }

    /// Stationary probability of exactly `n` customers in the system,
    /// `P(n) = (1 - rho) rho^n` (a geometric distribution).
    ///
    /// # Errors
    ///
    /// This never fails for a `usize` index — every non-negative `n` is
    /// in domain — but returns [`Result`] for symmetry with the
    /// signed-index helper [`state_probability_i`](Self::state_probability_i)
    /// and to leave room for future variants. The probability decays
    /// geometrically and underflows smoothly to `0.0` for large `n`.
    pub fn state_probability(&self, n: usize) -> Result<f64> {
        let rho = self.rho();
        Ok(self.p0() * rho.powi(n as i32))
    }

    /// Stationary probability `P(n)` accepting a signed index.
    ///
    /// # Errors
    ///
    /// [`QueueingError::Domain`] if `n < 0`, since a customer count
    /// cannot be negative.
    pub fn state_probability_i(&self, n: i64) -> Result<f64> {
        if n < 0 {
            return Err(QueueingError::domain(
                "state_probability",
                format!("customer count n must be >= 0, got {n}"),
            ));
        }
        let rho = self.rho();
        Ok(self.p0() * rho.powi(n as i32))
    }

    /// Cumulative probability of **at most** `n` customers in the
    /// system, `P(N <= n) = 1 - rho^(n + 1)` (the geometric CDF).
    pub fn cdf(&self, n: usize) -> f64 {
        let rho = self.rho();
        1.0 - rho.powi(n as i32 + 1)
    }

    /// Little's-law residual `L - lambda * W`.
    ///
    /// Little's law guarantees `L = lambda W` for any stable queue, so
    /// this is exactly `0` in theory; the returned value is the
    /// floating-point residual and is useful as a self-consistency
    /// check (its magnitude is at the rounding-error scale).
    pub fn little_residual(&self) -> f64 {
        self.l() - self.lambda * self.w()
    }

    /// Compute every steady-state metric at once and bundle them into a
    /// [`Metrics`] record (handy for serialization or display).
    pub fn metrics(&self) -> Metrics {
        Metrics {
            lambda: self.lambda,
            mu: self.mu,
            rho: self.rho(),
            p0: self.p0(),
            l: self.l(),
            lq: self.lq(),
            w: self.w(),
            wq: self.wq(),
        }
    }
}

/// A snapshot of all M/M/1 steady-state metrics, produced by
/// [`Mm1::metrics`].
///
/// Every field is finite for a queue built through [`Mm1::new`]. The
/// record is `Serialize`/`Deserialize` so a computed result can be
/// stored or sent over the wire.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Metrics {
    /// Arrival rate `lambda`.
    pub lambda: f64,
    /// Service rate `mu`.
    pub mu: f64,
    /// Traffic intensity / utilization `rho = lambda / mu`.
    pub rho: f64,
    /// Idle probability `P(0) = 1 - rho`.
    pub p0: f64,
    /// Mean number in system `L`.
    pub l: f64,
    /// Mean number in queue `Lq`.
    pub lq: f64,
    /// Mean time in system `W`.
    pub w: f64,
    /// Mean wait in queue `Wq`.
    pub wq: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute-difference assertion helper — never compares floats with
    /// `==`.
    fn close(a: f64, b: f64, eps: f64) {
        assert!(
            (a - b).abs() < eps,
            "expected {a} approx {b} (|diff| = {diff} >= eps = {eps})",
            diff = (a - b).abs()
        );
    }

    const EPS: f64 = 1e-12;

    // --- Ground-truth worked example -----------------------------------
    //
    // lambda = 2, mu = 3 (so rho = 2/3). Hand-computed exact values:
    //   rho = 2/3                = 0.6666666666...
    //   P0  = 1 - rho = 1/3      = 0.3333333333...
    //   L   = rho/(1-rho)
    //       = (2/3)/(1/3)        = 2
    //   Lq  = rho^2/(1-rho)
    //       = (4/9)/(1/3)        = 4/3 = 1.3333333333...
    //   W   = 1/(mu-lambda)
    //       = 1/1                = 1
    //   Wq  = rho/(mu-lambda)
    //       = (2/3)/1            = 2/3 = 0.6666666666...

    #[test]
    fn worked_example_lambda2_mu3() {
        let q = Mm1::new(2.0, 3.0).unwrap();
        close(q.rho(), 2.0 / 3.0, EPS);
        close(q.utilization(), 2.0 / 3.0, EPS);
        close(q.p0(), 1.0 / 3.0, EPS);
        close(q.l(), 2.0, EPS);
        close(q.lq(), 4.0 / 3.0, EPS);
        close(q.w(), 1.0, EPS);
        close(q.wq(), 2.0 / 3.0, EPS);
    }

    /// Second independent ground-truth point: lambda = 1, mu = 2,
    /// rho = 1/2. Then L = 1, Lq = 1/2, W = 1, Wq = 1/2, P0 = 1/2.
    #[test]
    fn worked_example_lambda1_mu2() {
        let q = Mm1::new(1.0, 2.0).unwrap();
        close(q.rho(), 0.5, EPS);
        close(q.p0(), 0.5, EPS);
        close(q.l(), 1.0, EPS);
        close(q.lq(), 0.5, EPS);
        close(q.w(), 1.0, EPS);
        close(q.wq(), 0.5, EPS);
    }

    // --- Identities that must hold for ALL stable (lambda, mu) ----------

    #[test]
    fn rho_is_lambda_over_mu() {
        for &(lambda, mu) in &[(2.0, 3.0), (1.0, 10.0), (9.9, 10.0), (0.5, 0.6)] {
            let q = Mm1::new(lambda, mu).unwrap();
            close(q.rho(), lambda / mu, EPS);
        }
    }

    #[test]
    fn lq_equals_l_minus_rho() {
        // Lq = L - rho is the standard reduction of rho^2/(1-rho).
        for &(lambda, mu) in &[(2.0, 3.0), (1.0, 4.0), (7.0, 8.0), (3.0, 100.0)] {
            let q = Mm1::new(lambda, mu).unwrap();
            close(q.lq(), q.l() - q.rho(), EPS);
        }
    }

    #[test]
    fn w_is_reciprocal_of_mu_minus_lambda() {
        for &(lambda, mu) in &[(2.0, 3.0), (5.0, 9.0), (1.0, 1.5)] {
            let q = Mm1::new(lambda, mu).unwrap();
            close(q.w(), 1.0 / (mu - lambda), EPS);
        }
    }

    #[test]
    fn wq_equals_w_minus_service_time() {
        // Wq = W - 1/mu: total time minus the mean service time.
        for &(lambda, mu) in &[(2.0, 3.0), (4.0, 7.0), (1.0, 100.0)] {
            let q = Mm1::new(lambda, mu).unwrap();
            close(q.wq(), q.w() - 1.0 / mu, EPS);
        }
    }

    #[test]
    fn littles_law_l_equals_lambda_w() {
        // L = lambda * W exactly; residual at rounding scale.
        for &(lambda, mu) in &[(2.0, 3.0), (1.0, 4.0), (9.0, 10.0), (0.25, 0.30)] {
            let q = Mm1::new(lambda, mu).unwrap();
            close(q.l(), lambda * q.w(), EPS);
            close(q.little_residual(), 0.0, EPS);
        }
    }

    #[test]
    fn littles_law_lq_equals_lambda_wq() {
        // Little's law also holds for the in-queue quantities.
        for &(lambda, mu) in &[(2.0, 3.0), (1.0, 4.0), (9.0, 10.0)] {
            let q = Mm1::new(lambda, mu).unwrap();
            close(q.lq(), lambda * q.wq(), EPS);
        }
    }

    // --- Stationary distribution ----------------------------------------

    #[test]
    fn state_probabilities_are_geometric_and_sum_to_one() {
        let q = Mm1::new(2.0, 3.0).unwrap();
        // P(0) = 1 - rho.
        close(q.state_probability(0).unwrap(), q.p0(), EPS);
        // P(n) = P(n-1) * rho — geometric ratio.
        let rho = q.rho();
        let mut prev = q.state_probability(0).unwrap();
        let mut sum = prev;
        for n in 1..=200 {
            let p = q.state_probability(n).unwrap();
            close(p, prev * rho, EPS);
            sum += p;
            prev = p;
        }
        // Geometric series sums to 1 (tail beyond n = 200 is tiny).
        close(sum, 1.0, 1e-9);
    }

    #[test]
    fn mean_of_state_distribution_recovers_l() {
        // E[N] = sum_n n P(n) must equal L = rho/(1-rho).
        let q = Mm1::new(0.75, 1.0).unwrap(); // rho = 0.75, L = 3.
        let mut mean = 0.0;
        for n in 0..=2000 {
            mean += (n as f64) * q.state_probability(n).unwrap();
        }
        close(mean, q.l(), 1e-6);
        close(q.l(), 3.0, EPS);
    }

    #[test]
    fn cdf_matches_summed_pmf() {
        let q = Mm1::new(3.0, 5.0).unwrap();
        for n in 0..=10 {
            let mut summed = 0.0;
            for k in 0..=n {
                summed += q.state_probability(k).unwrap();
            }
            close(q.cdf(n), summed, EPS);
        }
        // CDF closed form: 1 - rho^(n+1).
        close(q.cdf(0), 1.0 - q.rho(), EPS);
    }

    #[test]
    fn signed_state_probability_matches_unsigned() {
        let q = Mm1::new(1.0, 4.0).unwrap();
        for n in 0..=5 {
            close(
                q.state_probability_i(n as i64).unwrap(),
                q.state_probability(n).unwrap(),
                EPS,
            );
        }
    }

    // --- Limiting behaviour: rho -> 1 means everything diverges ---------

    #[test]
    fn metrics_diverge_as_rho_approaches_one() {
        // As lambda -> mu from below, L, Lq, W, Wq all grow without
        // bound. Check monotone growth toward +inf across a ramp.
        let mu = 1.0;
        let lambdas = [0.5, 0.9, 0.99, 0.999, 0.9999];
        let mut last_l = f64::NEG_INFINITY;
        let mut last_w = f64::NEG_INFINITY;
        for &lambda in &lambdas {
            let q = Mm1::new(lambda, mu).unwrap();
            assert!(q.l() > last_l, "L should increase toward rho=1");
            assert!(q.w() > last_w, "W should increase toward rho=1");
            assert!(q.l().is_finite());
            last_l = q.l();
            last_w = q.w();
        }
        // Very close to saturation, L is huge: rho = 0.9999 -> L ~ 9999.
        let q = Mm1::new(0.9999, 1.0).unwrap();
        close(q.l(), 9999.0, 1e-3);
    }

    // --- Metrics bundle mirrors the accessors ---------------------------

    #[test]
    fn metrics_bundle_matches_accessors() {
        let q = Mm1::new(2.0, 3.0).unwrap();
        let m = q.metrics();
        close(m.lambda, q.lambda(), EPS);
        close(m.mu, q.mu(), EPS);
        close(m.rho, q.rho(), EPS);
        close(m.p0, q.p0(), EPS);
        close(m.l, q.l(), EPS);
        close(m.lq, q.lq(), EPS);
        close(m.w, q.w(), EPS);
        close(m.wq, q.wq(), EPS);
    }

    #[test]
    fn metrics_round_trip_through_json() {
        let q = Mm1::new(2.0, 5.0).unwrap();
        let m = q.metrics();
        let json = serde_json::to_string(&m).unwrap();
        let back: Metrics = serde_json::from_str(&json).unwrap();
        close(back.l, m.l, EPS);
        close(back.wq, m.wq, EPS);
        // The queue itself round-trips too.
        let qj = serde_json::to_string(&q).unwrap();
        let q2: Mm1 = serde_json::from_str(&qj).unwrap();
        close(q2.lambda(), q.lambda(), EPS);
        close(q2.mu(), q.mu(), EPS);
    }

    // --- Validation / error paths ---------------------------------------

    #[test]
    fn rejects_unstable_load() {
        // rho = 1 exactly (lambda = mu): unstable.
        let err = Mm1::new(3.0, 3.0).unwrap_err();
        assert_eq!(err.code(), "queueing.unstable");
        // rho > 1.
        let err = Mm1::new(5.0, 4.0).unwrap_err();
        assert_eq!(err.code(), "queueing.unstable");
        match err {
            QueueingError::Unstable { lambda, mu, rho } => {
                close(lambda, 5.0, EPS);
                close(mu, 4.0, EPS);
                close(rho, 1.25, EPS);
            }
            other => panic!("expected Unstable, got {other:?}"),
        }
    }

    #[test]
    fn rejects_non_positive_rates() {
        assert_eq!(Mm1::new(0.0, 3.0).unwrap_err().code(), "queueing.invalid");
        assert_eq!(Mm1::new(-1.0, 3.0).unwrap_err().code(), "queueing.invalid");
        assert_eq!(Mm1::new(2.0, 0.0).unwrap_err().code(), "queueing.invalid");
        assert_eq!(Mm1::new(2.0, -3.0).unwrap_err().code(), "queueing.invalid");
    }

    #[test]
    fn rejects_non_finite_rates() {
        assert_eq!(
            Mm1::new(f64::NAN, 3.0).unwrap_err().code(),
            "queueing.invalid"
        );
        assert_eq!(
            Mm1::new(2.0, f64::INFINITY).unwrap_err().code(),
            "queueing.invalid"
        );
        // An infinite arrival rate is also rejected (non-finite check
        // fires before the stability check).
        assert_eq!(
            Mm1::new(f64::INFINITY, 3.0).unwrap_err().code(),
            "queueing.invalid"
        );
    }

    #[test]
    fn rejects_negative_state_index() {
        let q = Mm1::new(2.0, 3.0).unwrap();
        let err = q.state_probability_i(-1).unwrap_err();
        assert_eq!(err.code(), "queueing.domain");
    }

    // --- Capacity sizing: service rate for a target response time -------

    #[test]
    fn service_rate_for_response_time_inverts_w() {
        // Inverse of the headline example: lambda = 2, W = 1 -> mu = 3.
        let mu = Mm1::service_rate_for_mean_response_time(2.0, 1.0).unwrap();
        close(mu, 3.0, EPS);
        // Round-trip over a sweep: build the queue and recover the target W.
        for &(lambda, w) in &[(2.0_f64, 1.0), (5.0, 0.25), (1.0, 4.0), (10.0, 0.05)] {
            let mu = Mm1::service_rate_for_mean_response_time(lambda, w).unwrap();
            let q = Mm1::new(lambda, mu).unwrap();
            close(q.w(), w, 1e-12 * w.max(1.0));
        }
    }

    #[test]
    fn service_rate_for_response_time_is_always_stable_and_monotone() {
        // mu = lambda + 1/W > lambda, so the queue is always constructible.
        let mu = Mm1::service_rate_for_mean_response_time(100.0, 1e-3).unwrap();
        assert!(mu > 100.0);
        assert!(Mm1::new(100.0, mu).is_ok());
        // A tighter (smaller) target response time demands a higher rate.
        let slow = Mm1::service_rate_for_mean_response_time(5.0, 2.0).unwrap();
        let fast = Mm1::service_rate_for_mean_response_time(5.0, 0.5).unwrap();
        assert!(fast > slow, "fast {fast} should exceed slow {slow}");
    }

    #[test]
    fn service_rate_for_response_time_rejects_bad_inputs() {
        assert!(Mm1::service_rate_for_mean_response_time(0.0, 1.0).is_err()); // lambda
        assert!(Mm1::service_rate_for_mean_response_time(2.0, 0.0).is_err()); // W
        assert!(Mm1::service_rate_for_mean_response_time(-1.0, 1.0).is_err());
        assert!(Mm1::service_rate_for_mean_response_time(2.0, f64::NAN).is_err());
    }

    // --- Load sizing: arrival rate for a target response time ----------

    #[test]
    fn arrival_rate_for_response_time_inverts_w() {
        // Inverse of the headline example: mu = 3, W = 1 -> lambda = 2.
        let lambda = Mm1::arrival_rate_for_mean_response_time(3.0, 1.0).unwrap();
        close(lambda, 2.0, EPS);
        // Round-trip over a sweep: build the queue and recover the target W.
        for &(mu, w) in &[(3.0_f64, 1.0), (10.0, 0.2), (2.0, 5.0), (100.0, 0.02)] {
            let lambda = Mm1::arrival_rate_for_mean_response_time(mu, w).unwrap();
            let q = Mm1::new(lambda, mu).unwrap();
            close(q.w(), w, 1e-12 * w.max(1.0));
        }
    }

    #[test]
    fn arrival_and_service_inverses_are_consistent() {
        // The two inverses of W = 1/(mu - lambda) compose: sizing mu for
        // (lambda, W) then the arrival rate for (mu, W) recovers lambda.
        let (lambda, w) = (5.0, 0.25);
        let mu = Mm1::service_rate_for_mean_response_time(lambda, w).unwrap();
        let lambda_back = Mm1::arrival_rate_for_mean_response_time(mu, w).unwrap();
        close(lambda_back, lambda, EPS);
    }

    #[test]
    fn arrival_rate_for_response_time_more_tolerant_admits_more_load() {
        // A larger (more tolerant) target W admits a higher arrival rate,
        // but it always stays below mu for stability.
        let mu = 10.0;
        let tight = Mm1::arrival_rate_for_mean_response_time(mu, 0.2).unwrap();
        let loose = Mm1::arrival_rate_for_mean_response_time(mu, 1.0).unwrap();
        assert!(loose > tight, "loose {loose} should exceed tight {tight}");
        assert!(loose < mu, "arrival rate must stay below mu for stability");
    }

    #[test]
    fn arrival_rate_for_response_time_rejects_unreachable_and_bad() {
        // W at or below 1/mu is unreachable (would need lambda <= 0).
        assert_eq!(
            Mm1::arrival_rate_for_mean_response_time(3.0, 1.0 / 3.0)
                .unwrap_err()
                .code(),
            "queueing.invalid"
        );
        assert!(Mm1::arrival_rate_for_mean_response_time(3.0, 0.1).is_err()); // < 1/mu
        assert!(Mm1::arrival_rate_for_mean_response_time(0.0, 1.0).is_err()); // mu
        assert!(Mm1::arrival_rate_for_mean_response_time(3.0, 0.0).is_err()); // W
        assert!(Mm1::arrival_rate_for_mean_response_time(3.0, f64::NAN).is_err());
    }
}
