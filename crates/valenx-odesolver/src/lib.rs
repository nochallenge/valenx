//! # valenx-odesolver
//!
//! Fixed-step explicit integrators for first-order ordinary differential
//! equations (ODEs).
//!
//! ## What
//!
//! Three classical, explicit, single-step Runge-Kutta schemes for the
//! initial-value problem `dy/dt = f(t, y)`, `y(t0) = y0`:
//!
//! - explicit (forward) Euler, classical order 1;
//! - Heun's method (improved Euler / trapezoidal RK2), classical order 2;
//! - classical 4-stage Runge-Kutta (RK4), classical order 4.
//!
//! Each scheme is offered in two forms:
//!
//! - a [`scalar`] form where the state is a single `f64`, and
//! - a [`system`] form where the state is a `Vec<f64>` of fixed length (used
//!   for higher-order ODEs reduced to first order, such as the harmonic
//!   oscillator `[x, v]`).
//!
//! The step size `dt` is fixed for an entire run; there is no adaptive
//! step-size control. Pick the method with [`Method`].
//!
//! ## Model
//!
//! Writing `h` for the step size, one step from `(t, y)` to `(t + h, y_next)`:
//!
//! - Euler: `y_next = y + h * f(t, y)`.
//! - Heun: `k1 = f(t, y)`; `k2 = f(t + h, y + h*k1)`;
//!   `y_next = y + (h/2) * (k1 + k2)`.
//! - RK4: `k1 = f(t, y)`; `k2 = f(t + h/2, y + (h/2)*k1)`;
//!   `k3 = f(t + h/2, y + (h/2)*k2)`; `k4 = f(t + h, y + h*k3)`;
//!   `y_next = y + (h/6) * (k1 + 2*k2 + 2*k3 + k4)`.
//!
//! The global (accumulated) error of an order-`p` method scales like `O(h^p)`,
//! so halving `h` shrinks the error by roughly `2^p`. The crate's test suite
//! verifies this against analytic ground truth: RK4 reproducing `exp(t)`
//! (`dy/dt = y`), the `~16x` error reduction under step-halving for RK4 (`~2x`
//! for Euler), the decay `dy/dt = -y` to `exp(-t)`, and approximate energy
//! conservation of the undamped harmonic oscillator under RK4.
//!
//! ## Honest scope
//!
//! Research/educational grade. These are textbook closed-form / numerical
//! models implemented for clarity and pedagogy, validated against analytic
//! solutions; they are NOT a clinical/medical/production engineering tool. In
//! particular: the methods are explicit and unconditionally suited only to
//! non-stiff problems, there is no adaptive error control, no implicit /
//! symplectic / stiff solver, and no dense output or event detection. For
//! stiff systems, long-horizon Hamiltonian integration, or any
//! safety-critical use, reach for a purpose-built, validated solver instead.
//!
//! ## Example
//!
//! ```
//! use valenx_odesolver::{scalar, Method};
//!
//! // Solve dy/dt = y, y(0) = 1 to t = 1; the answer is e.
//! let yf = scalar::integrate_final(Method::Rk4, |_t, y| y, 0.0, 1.0, 0.01, 100)
//!     .expect("valid inputs");
//! assert!((yf - std::f64::consts::E).abs() < 1e-8);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod method;
pub mod scalar;
pub mod system;

pub use error::OdeError;
pub use method::Method;

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    #[test]
    fn method_metadata_is_stable() {
        assert_eq!(Method::Euler.order(), 1);
        assert_eq!(Method::Heun.order(), 2);
        assert_eq!(Method::Rk4.order(), 4);
        assert_eq!(Method::Euler.name(), "euler");
        assert_eq!(Method::Heun.name(), "heun");
        assert_eq!(Method::Rk4.name(), "rk4");
    }

    #[test]
    fn error_constructors_validate() {
        // Steps.
        assert!(OdeError::bad_step(0.0).is_some());
        assert!(OdeError::bad_step(-1.0).is_some());
        assert!(OdeError::bad_step(f64::NAN).is_some());
        assert!(OdeError::bad_step(f64::INFINITY).is_some());
        assert!(OdeError::bad_step(1e-6).is_none());
        // Counts.
        assert!(OdeError::bad_step_count(0).is_some());
        assert!(OdeError::bad_step_count(1).is_none());
        // Finiteness.
        assert!(OdeError::non_finite("y", f64::NAN).is_some());
        assert!(OdeError::non_finite("y", 1.0).is_none());
        // Dimensions.
        assert!(OdeError::dimension_mismatch(2, 3).is_some());
        assert!(OdeError::dimension_mismatch(2, 2).is_none());
    }

    #[test]
    fn error_codes_are_distinct() {
        let codes = [
            OdeError::BadStep {
                dt: 0.0,
                reason: "x",
            }
            .code(),
            OdeError::BadStepCount.code(),
            OdeError::NonFinite {
                name: "y",
                value: f64::NAN,
            }
            .code(),
            OdeError::DimensionMismatch {
                expected: 1,
                actual: 2,
            }
            .code(),
        ];
        // All four codes are unique.
        for (i, a) in codes.iter().enumerate() {
            for b in &codes[i + 1..] {
                assert_ne!(a, b, "duplicate error code {a}");
            }
        }
    }

    #[test]
    fn scalar_and_system_agree_on_the_same_scalar_problem() {
        // dy/dt = y solved as a scalar and as a length-1 system must match.
        let s = scalar::integrate_final(Method::Rk4, |_t, y| y, 0.0, 1.0, 0.01, 100).unwrap();
        let v = system::integrate_final(
            Method::Rk4,
            |_t, y: &[f64]| vec![y[0]],
            0.0,
            &[1.0],
            0.01,
            100,
        )
        .unwrap();
        assert!(close(s, v[0], 1e-12), "scalar {s} vs system {}", v[0]);
    }
}
