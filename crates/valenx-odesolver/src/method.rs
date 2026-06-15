//! The integration-method selector shared by the scalar and system solvers.
//!
//! ## Model
//!
//! All three methods are explicit, fixed-step, single-step Runge-Kutta
//! schemes for the first-order initial-value problem
//!
//! ```text
//!   dy/dt = f(t, y),   y(t0) = y0.
//! ```
//!
//! Writing `h` for the step size, one step from `(t, y)` to `(t + h, y_next)`
//! is:
//!
//! - [`Method::Euler`] (explicit / forward Euler, order 1):
//!   `y_next = y + h * f(t, y)`.
//! - [`Method::Heun`] (improved Euler / trapezoidal RK2, order 2):
//!   `k1 = f(t, y)`, `k2 = f(t + h, y + h*k1)`,
//!   `y_next = y + (h/2) * (k1 + k2)`.
//! - [`Method::Rk4`] (classical Runge-Kutta, order 4):
//!   `k1 = f(t, y)`,
//!   `k2 = f(t + h/2, y + (h/2)*k1)`,
//!   `k3 = f(t + h/2, y + (h/2)*k2)`,
//!   `k4 = f(t + h, y + h*k3)`,
//!   `y_next = y + (h/6) * (k1 + 2*k2 + 2*k3 + k4)`.
//!
//! The global error of a fixed-step method of order `p` scales like `O(h^p)`,
//! so halving `h` shrinks the error by roughly `2^p` (about 2x for Euler, 4x
//! for Heun, 16x for RK4). The crate's tests pin these factors against
//! analytic solutions.

use serde::{Deserialize, Serialize};

/// Which fixed-step explicit Runge-Kutta scheme to use.
///
/// The variants are ordered by classical order of accuracy; see the
/// [module documentation](self) for the per-step formulas.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub enum Method {
    /// Explicit (forward) Euler. Classical order 1.
    Euler,
    /// Heun's method (improved Euler, trapezoidal RK2). Classical order 2.
    Heun,
    /// Classical 4-stage Runge-Kutta. Classical order 4.
    Rk4,
}

impl Method {
    /// The classical order of accuracy `p` of the global truncation error,
    /// i.e. global error scales like `O(h^p)`.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_odesolver::Method;
    /// assert_eq!(Method::Euler.order(), 1);
    /// assert_eq!(Method::Heun.order(), 2);
    /// assert_eq!(Method::Rk4.order(), 4);
    /// ```
    #[must_use]
    pub fn order(self) -> u32 {
        match self {
            Method::Euler => 1,
            Method::Heun => 2,
            Method::Rk4 => 4,
        }
    }

    /// Short, stable, lower-case name for logging and display.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_odesolver::Method;
    /// assert_eq!(Method::Rk4.name(), "rk4");
    /// ```
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Method::Euler => "euler",
            Method::Heun => "heun",
            Method::Rk4 => "rk4",
        }
    }
}
