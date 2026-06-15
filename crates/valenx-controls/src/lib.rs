//! # valenx-controls — linear control-systems toolkit
//!
//! A small, dependency-light library of the closed-form and discrete
//! building blocks that every introductory controls course leans on:
//! second-order step-response metrics, a discrete PID controller, and the
//! Routh-Hurwitz stability test for a quadratic characteristic equation.
//!
//! ## What
//!
//! - **Second-order step metrics** ([`second_order`]) — wrap a system in
//!   a [`SecondOrder`] `{ wn, zeta }` and read back its unit-step
//!   performance: fractional/percent **overshoot**, **peak time**,
//!   2% / 5% **settling time**, **rise time** and **damped frequency**,
//!   each a textbook closed form.
//! - **Discrete PID** ([`pid`]) — a stateful fixed-step
//!   proportional-integral-derivative controller
//!   (`u = Kp*e + Ki*I + Kd*D`) with output clamping and integral
//!   anti-windup, stepped once per sample by [`Pid::update`].
//! - **Stability** ([`stability`]) — a validated quadratic characteristic
//!   polynomial [`QuadraticChar`] with a Routh-Hurwitz
//!   [`is_stable`](QuadraticChar::is_stable) test and closed-form
//!   [`roots`](QuadraticChar::roots) so callers can inspect the pole real
//!   parts directly.
//!
//! ```
//! use valenx_controls::{SecondOrder, QuadraticChar};
//!
//! // A lightly-damped system overshoots and is stable.
//! let sys = SecondOrder::new(4.0, 0.5).expect("valid wn, zeta");
//! let m = sys.step_metrics().expect("underdamped");
//! assert!((m.overshoot - 0.163).abs() < 1e-3); // ~16% overshoot
//!
//! // Its characteristic polynomial s^2 + 4s + 16 is stable.
//! let poly = QuadraticChar::from_wn_zeta(4.0, 0.5).expect("finite");
//! assert!(poly.is_stable());
//! ```
//!
//! ## Model
//!
//! A canonical second-order system has the transfer function
//! `wn^2 / (s^2 + 2*zeta*wn*s + wn^2)`. For the underdamped case
//! `0 <= zeta < 1` the damped frequency is `wd = wn*sqrt(1 - zeta^2)` and
//! the step-response figures are
//!
//! - overshoot `Mp = exp(-pi*zeta/sqrt(1 - zeta^2))`,
//! - peak time `tp = pi/wd`,
//! - settling time `ts ~ 4/(zeta*wn)` (2% band; `3/(zeta*wn)` for 5%),
//! - rise time `tr = (pi - acos(zeta))/wd`.
//!
//! The PID law is integrated by the rectangle rule (`I += e*dt`) and
//! differentiated by a backward difference (`D = (e - e_prev)/dt`). The
//! Routh-Hurwitz test for `a*s^2 + b*s + c` reduces to "stable iff all
//! coefficients share a sign", i.e. `b > 0` and `c > 0` after normalising
//! `a > 0`; equivalently both roots have negative real part.
//!
//! ## Honest scope
//!
//! Research/educational grade. These are well-established **textbook
//! closed-form and discrete numerical models** — the formulae found in
//! Ogata or Franklin/Powell/Emami-Naeini and in countless embedded loops.
//! They are exact for the idealised models they describe, but the library
//! is deliberately narrow: it covers only the standard second-order
//! prototype, a single-loop fixed-step PID, and quadratic stability. It
//! does **not** model higher-order or MIMO plants, frequency-domain
//! design (Bode / Nyquist / root-locus), state-space / pole-placement /
//! LQR/LQG, discretisation of arbitrary transfer functions, derivative
//! filtering, or robustness margins. It is **NOT** a clinical, medical,
//! safety-certified, or production control-engineering tool, and must not
//! be used to design or certify real-world safety-critical control loops.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod pid;
pub mod second_order;
pub mod stability;

// --- Convenience re-exports of the most-used types --------------------

pub use error::{ControlsError, ErrorCategory, Result};
pub use pid::{Pid, PidConfig};
pub use second_order::{DampingRegime, SecondOrder, StepMetrics};
pub use stability::{QuadraticChar, Root};

#[cfg(test)]
mod tests {
    use super::*;

    /// End-to-end: a lightly-damped second-order system overshoots, has a
    /// finite settling time, and its characteristic polynomial is stable
    /// — and the two stability views (the `SecondOrder` damping sign and
    /// the Routh-Hurwitz test) agree.
    #[test]
    fn second_order_metrics_and_stability_agree() {
        let sys = SecondOrder::new(5.0, 0.3).unwrap();
        let m = sys.step_metrics().unwrap();

        // Underdamped -> positive overshoot, positive finite times.
        assert!(m.overshoot > 0.0 && m.overshoot < 1.0);
        assert!(m.peak_time > 0.0);
        assert!(m.settling_time > 0.0);
        assert!(m.rise_time > 0.0);

        // The matching characteristic polynomial is stable (zeta > 0).
        let poly = QuadraticChar::from_wn_zeta(5.0, 0.3).unwrap();
        assert!(poly.is_stable());
        assert!(poly.max_real_part() < 0.0);

        // Drop the damping to zero: marginal/unstable boundary, and the
        // oscillatory metrics become a domain error.
        let undamped = SecondOrder::new(5.0, 0.0).unwrap();
        assert!(undamped.settling_time_2pct().is_err());
        assert!(!QuadraticChar::from_wn_zeta(5.0, 0.0).unwrap().is_stable());
    }

    /// End-to-end: a PID controller drives a trivial integrator plant
    /// `x[k+1] = x[k] + u*dt` from 0 toward a unit setpoint, and the
    /// tracking error shrinks over time.
    #[test]
    fn pid_drives_a_simple_plant_toward_setpoint() {
        let cfg = PidConfig::new(2.0, 0.5, 0.0, 0.05)
            .unwrap()
            .with_limits(-10.0, 10.0)
            .unwrap();
        let mut pid = Pid::new(cfg);

        let setpoint = 1.0_f64;
        let mut x = 0.0_f64; // plant state
        let dt = 0.05;
        let e0 = (setpoint - x).abs();
        for _ in 0..400 {
            let u = pid.update(setpoint, x);
            x += u * dt; // integrator plant
        }
        let e_final = (setpoint - x).abs();
        assert!(
            e_final < e0,
            "error should shrink: started {e0}, ended {e_final}"
        );
        assert!(e_final < 0.05, "did not converge: final error {e_final}");
    }
}
