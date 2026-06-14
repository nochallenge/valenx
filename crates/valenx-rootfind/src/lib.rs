//! # valenx-rootfind — scalar root finding
//!
//! ## What
//!
//! Solve `f(x) = 0` for a single real unknown, given any closure
//! `f: Fn(f64) -> f64`. The crate provides the three classical
//! one-dimensional methods, each returning a [`Result`] that carries the
//! root and its diagnostics on success, or a typed
//! [`RootError`] on failure:
//!
//! - [`bisection`] — the bracketing method. Needs a sign-change bracket
//!   `[a, b]` (one where `f(a)` and `f(b)` have opposite signs) and is
//!   the only method here whose convergence is *guaranteed*.
//! - [`newton`] — the Newton-Raphson tangent method. Needs the analytic
//!   derivative `f'`; converges quadratically near a simple root.
//! - [`secant`] — the chord method. Needs two seed points but **no
//!   derivative**; converges super-linearly (order ≈ 1.618).
//!
//! All three share a [`Settings`] (tolerance + iteration cap) and a
//! [`Root`] result (root, residual, iteration count, method tag).
//!
//! ```
//! use valenx_rootfind::{bisection, newton, secant, Settings};
//!
//! // f(x) = x^2 - 2  ⇒  positive root is sqrt(2).
//! let f = |x: f64| x * x - 2.0;
//! let df = |x: f64| 2.0 * x;
//! let cfg = Settings::default();
//!
//! let by_bisection = bisection(f, 0.0, 2.0, cfg).unwrap();
//! let by_newton = newton(f, df, 1.0, cfg).unwrap();
//! let by_secant = secant(f, 1.0, 2.0, cfg).unwrap();
//!
//! let truth = 2.0_f64.sqrt();
//! assert!((by_bisection.root - truth).abs() < 1e-10);
//! assert!((by_newton.root - truth).abs() < 1e-10);
//! assert!((by_secant.root - truth).abs() < 1e-10);
//! ```
//!
//! ## Model
//!
//! Every method is the textbook iteration, stopped on a residual /
//! width / step tolerance:
//!
//! - **Bisection** repeatedly halves the bracket, keeping the half that
//!   still straddles the sign change. After `n` halvings the root is
//!   localised to `(b - a) / 2^n`, so it converges linearly and the step
//!   count is bounded by `log2((b - a) / tol)`.
//! - **Newton-Raphson** follows the tangent line at the current iterate
//!   to its `x`-intercept: `x_{n+1} = x_n - f(x_n) / f'(x_n)`. Near a
//!   simple root the error squares each step (quadratic convergence).
//! - **Secant** replaces the analytic derivative with the slope through
//!   the two latest iterates:
//!   `x_{n+1} = x_n - f(x_n) (x_n - x_{n-1}) / (f(x_n) - f(x_{n-1}))`.
//!   The convergence order is the golden ratio `φ ≈ 1.618`.
//!
//! Each solver also short-circuits when an input point is already a root
//! and validates its arguments up front (positive finite tolerance,
//! non-zero budget, finite seeds, distinct secant seeds).
//!
//! ## Honest scope
//!
//! Research / educational grade: these are textbook closed-form /
//! numerical models, **not** a clinical / medical / production
//! engineering tool. They are the methods exactly as taught in a first
//! numerical-analysis course and make no claim to the robustness of a
//! production hybrid solver (e.g. Brent's method, which safeguards the
//! open steps inside a maintained bracket):
//!
//! - Newton and the secant method are *unsafeguarded* — from a poor
//!   start they can diverge, oscillate, or wander off to a different
//!   root; only bisection is guaranteed to converge.
//! - A vanishing slope (a flat Newton tangent or a horizontal secant
//!   chord) is reported as [`RootError::ZeroDerivative`] rather than
//!   recovered from.
//! - Everything is scalar `f64`; there is no complex-root, no
//!   multiple-root acceleration, and no systems-of-equations support.
//! - Convergence near a *multiple* root degrades from the quoted orders
//!   (Newton drops to linear) — this is a property of the classical
//!   methods, not a bug.
//!
//! All numerics are pure and deterministic; there are no external
//! processes and the crate contains no `unsafe` code.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod solver;

// --- Convenience re-exports of the most-used items --------------------

pub use error::{ErrorCategory, Result, RootError};
pub use solver::{
    bisection, newton, secant, Method, Root, Settings, DEFAULT_MAX_ITER, DEFAULT_TOL, SLOPE_GUARD,
};

#[cfg(test)]
mod tests {
    use super::*;

    /// End-to-end: the three methods all converge to the same root of a
    /// non-trivial transcendental equation, `x e^x = 1` (whose root is
    /// the Omega constant `≈ 0.5671432904`).
    #[test]
    fn three_methods_agree_on_omega_constant() {
        let f = |x: f64| x * x.exp() - 1.0;
        let df = |x: f64| (x + 1.0) * x.exp();
        let omega = 0.567_143_290_409_784;
        let eps = 1e-10;

        let b = bisection(f, 0.0, 1.0, Settings::default()).unwrap();
        let n = newton(f, df, 0.5, Settings::default()).unwrap();
        let s = secant(f, 0.0, 1.0, Settings::default()).unwrap();

        assert!((b.root - omega).abs() < eps, "bisection {}", b.root);
        assert!((n.root - omega).abs() < eps, "newton {}", n.root);
        assert!((s.root - omega).abs() < eps, "secant {}", s.root);

        // Newton should not need more steps than bisection here.
        assert!(n.iterations <= b.iterations);
    }

    #[test]
    fn re_exports_are_reachable() {
        // The public surface is wired up.
        let _ = (DEFAULT_TOL, DEFAULT_MAX_ITER, SLOPE_GUARD);
        let cfg = Settings::default();
        let r = newton(|x| x - 3.0, |_| 1.0, 0.0, cfg).unwrap();
        assert!((r.root - 3.0).abs() < 1e-12);
        assert_eq!(r.method, Method::Newton);
    }
}
