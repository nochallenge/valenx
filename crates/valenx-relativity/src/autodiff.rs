//! Forward-mode automatic differentiation.
//!
//! The curvature engine needs the metric's first derivatives (for Christoffel
//! symbols) and second derivatives (for the Riemann tensor). Rather than
//! finite-differencing, we evaluate the *same* generic metric function over
//! special number types that carry derivative information exactly:
//!
//! * [`Dual`] — value plus one first derivative.
//! * [`HyperDual`] — value, two first derivatives, and the mixed second
//!   derivative `∂²/(∂dir1 ∂dir2)`.
//!
//! Seeding a coordinate's derivative slot to `1` and evaluating the metric then
//! reads the partial derivative straight out of the result, to machine
//! precision and with no step-size to tune.
//!
//! All three of [`f64`], [`Dual`] and [`HyperDual`] implement [`Scalar`], so a
//! metric written once as `fn metric<T: Scalar>(..)` is reused unchanged for
//! plain evaluation and for both derivative passes.

use core::ops::{Add, Div, Mul, Neg, Sub};

/// A real-like scalar usable inside metric definitions: field arithmetic plus
/// the handful of transcendental functions the built-in metrics need.
pub trait Scalar:
    Copy
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
{
    /// Lift a real constant into the scalar type (all derivatives zero).
    fn from_f64(x: f64) -> Self;
    /// The real (value) component.
    fn re(self) -> f64;
    /// Sine.
    fn sin(self) -> Self;
    /// Cosine.
    fn cos(self) -> Self;
    /// Square root.
    fn sqrt(self) -> Self;
    /// Reciprocal `1 / self`.
    fn recip(self) -> Self;
    /// Square, `self * self`.
    fn sq(self) -> Self {
        self * self
    }
}

impl Scalar for f64 {
    fn from_f64(x: f64) -> Self {
        x
    }
    fn re(self) -> f64 {
        self
    }
    fn sin(self) -> Self {
        f64::sin(self)
    }
    fn cos(self) -> Self {
        f64::cos(self)
    }
    fn sqrt(self) -> Self {
        f64::sqrt(self)
    }
    fn recip(self) -> Self {
        1.0 / self
    }
}

/// A first-order dual number `v + d·ε` with `ε² = 0`; tracks one derivative.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Dual {
    /// Value component.
    pub v: f64,
    /// First-derivative component.
    pub d: f64,
}

impl Dual {
    /// A constant: value `v`, zero derivative.
    pub fn constant(v: f64) -> Self {
        Self { v, d: 0.0 }
    }
    /// A seed variable: value `v`, unit derivative.
    pub fn variable(v: f64) -> Self {
        Self { v, d: 1.0 }
    }
}

impl Add for Dual {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self {
            v: self.v + o.v,
            d: self.d + o.d,
        }
    }
}
impl Sub for Dual {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self {
            v: self.v - o.v,
            d: self.d - o.d,
        }
    }
}
impl Mul for Dual {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        Self {
            v: self.v * o.v,
            d: self.v * o.d + self.d * o.v,
        }
    }
}
impl Neg for Dual {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            v: -self.v,
            d: -self.d,
        }
    }
}
// `div` is multiplication by the reciprocal, so the differentiation rule lives
// in one place (`recip`). clippy::suspicious_arithmetic_impl flags the `*`,
// which is exactly what we intend here; correctness is covered by tests.
#[allow(clippy::suspicious_arithmetic_impl)]
impl Div for Dual {
    type Output = Self;
    fn div(self, o: Self) -> Self {
        self * o.recip()
    }
}

impl Scalar for Dual {
    fn from_f64(x: f64) -> Self {
        Self::constant(x)
    }
    fn re(self) -> f64 {
        self.v
    }
    fn sin(self) -> Self {
        Self {
            v: self.v.sin(),
            d: self.d * self.v.cos(),
        }
    }
    fn cos(self) -> Self {
        Self {
            v: self.v.cos(),
            d: -self.d * self.v.sin(),
        }
    }
    fn sqrt(self) -> Self {
        let s = self.v.sqrt();
        Self {
            v: s,
            d: self.d / (2.0 * s),
        }
    }
    fn recip(self) -> Self {
        Self {
            v: 1.0 / self.v,
            d: -self.d / (self.v * self.v),
        }
    }
}

/// A hyper-dual number carrying value, two first derivatives (along independent
/// seed directions 1 and 2), and the mixed second derivative `∂²/(∂1 ∂2)`.
///
/// Seeding direction 1 on coordinate `c` and direction 2 on coordinate `d`,
/// then evaluating the metric, yields in one pass: `∂_c g` (`d1`), `∂_d g`
/// (`d2`) and `∂_c∂_d g` (`d12`). Seeding both directions on the same
/// coordinate gives the pure second derivative.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HyperDual {
    /// Value component.
    pub v: f64,
    /// First derivative along seed direction 1.
    pub d1: f64,
    /// First derivative along seed direction 2.
    pub d2: f64,
    /// Mixed second derivative `∂²/(∂1 ∂2)`.
    pub d12: f64,
}

impl HyperDual {
    /// A constant: value `v`, all derivatives zero.
    pub fn constant(v: f64) -> Self {
        Self {
            v,
            d1: 0.0,
            d2: 0.0,
            d12: 0.0,
        }
    }
    /// A seed at value `v`, with unit first derivative along the chosen
    /// direction(s). Set `dir1`/`dir2` to mark this variable as the one being
    /// differentiated along seed direction 1 / direction 2.
    pub fn seed(v: f64, dir1: bool, dir2: bool) -> Self {
        Self {
            v,
            d1: if dir1 { 1.0 } else { 0.0 },
            d2: if dir2 { 1.0 } else { 0.0 },
            d12: 0.0,
        }
    }
}

impl Add for HyperDual {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self {
            v: self.v + o.v,
            d1: self.d1 + o.d1,
            d2: self.d2 + o.d2,
            d12: self.d12 + o.d12,
        }
    }
}
impl Sub for HyperDual {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self {
            v: self.v - o.v,
            d1: self.d1 - o.d1,
            d2: self.d2 - o.d2,
            d12: self.d12 - o.d12,
        }
    }
}
impl Mul for HyperDual {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        Self {
            v: self.v * o.v,
            d1: self.v * o.d1 + self.d1 * o.v,
            d2: self.v * o.d2 + self.d2 * o.v,
            d12: self.v * o.d12 + self.d1 * o.d2 + self.d2 * o.d1 + self.d12 * o.v,
        }
    }
}
impl Neg for HyperDual {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            v: -self.v,
            d1: -self.d1,
            d2: -self.d2,
            d12: -self.d12,
        }
    }
}
#[allow(clippy::suspicious_arithmetic_impl)] // see the note on Dual's Div impl
impl Div for HyperDual {
    type Output = Self;
    fn div(self, o: Self) -> Self {
        self * o.recip()
    }
}

impl HyperDual {
    /// Apply a smooth real function via the chain rule, given `f`, `f'` and
    /// `f''` evaluated at the value component. This is the single place the
    /// second-order propagation rule lives.
    fn chain(self, f: f64, df: f64, ddf: f64) -> Self {
        Self {
            v: f,
            d1: df * self.d1,
            d2: df * self.d2,
            d12: df * self.d12 + ddf * self.d1 * self.d2,
        }
    }
}

impl Scalar for HyperDual {
    fn from_f64(x: f64) -> Self {
        Self::constant(x)
    }
    fn re(self) -> f64 {
        self.v
    }
    fn sin(self) -> Self {
        self.chain(self.v.sin(), self.v.cos(), -self.v.sin())
    }
    fn cos(self) -> Self {
        self.chain(self.v.cos(), -self.v.sin(), -self.v.cos())
    }
    fn sqrt(self) -> Self {
        let s = self.v.sqrt();
        self.chain(s, 0.5 / s, -0.25 / (self.v * s))
    }
    fn recip(self) -> Self {
        let inv = 1.0 / self.v;
        self.chain(inv, -inv * inv, 2.0 * inv * inv * inv)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) {
        assert!((a - b).abs() <= tol, "{a} vs {b} (tol {tol})");
    }

    #[test]
    fn dual_polynomial_derivative() {
        // f(x) = x^3 - 2x ; f'(x) = 3x^2 - 2. At x = 4: f=56, f'=46.
        let x = Dual::variable(4.0);
        let f = x * x * x - Dual::constant(2.0) * x;
        approx(f.v, 56.0, 1e-12);
        approx(f.d, 46.0, 1e-12);
    }

    #[test]
    fn dual_trig_and_div() {
        // f(x) = sin(x)/x at x = 1: value sin1, derivative (cos1 - sin1)/1.
        let x = Dual::variable(1.0);
        let f = x.sin() / x;
        approx(f.v, 1.0_f64.sin(), 1e-12);
        approx(f.d, 1.0_f64.cos() - 1.0_f64.sin(), 1e-12);
    }

    #[test]
    fn hyperdual_second_derivative_pure() {
        // f(x) = x^4 ; f'=4x^3, f''=12x^2. At x=2: 16,32,48.
        let x = HyperDual::seed(2.0, true, true);
        let f = x * x * x * x;
        approx(f.v, 16.0, 1e-10);
        approx(f.d1, 32.0, 1e-10);
        approx(f.d2, 32.0, 1e-10);
        approx(f.d12, 48.0, 1e-10);
    }

    #[test]
    fn hyperdual_mixed_second_derivative() {
        // f(x,y) = x^2 * y^3. ∂x = 2x y^3, ∂y = 3 x^2 y^2, ∂x∂y = 6 x y^2.
        // At (x,y) = (3, 2): seed x along dir1, y along dir2.
        let x = HyperDual::seed(3.0, true, false);
        let y = HyperDual::seed(2.0, false, true);
        let f = x * x * y * y * y;
        approx(f.v, 9.0 * 8.0, 1e-10); // 72
        approx(f.d1, 2.0 * 3.0 * 8.0, 1e-10); // 48
        approx(f.d2, 3.0 * 9.0 * 4.0, 1e-10); // 108
        approx(f.d12, 6.0 * 3.0 * 4.0, 1e-10); // 72
    }

    #[test]
    fn hyperdual_recip_and_sin_second_derivative() {
        // f(x) = sin(x)/x; check value and that d1==d2 (same seed) and the
        // second derivative matches the analytic f''(x).
        let x = HyperDual::seed(1.3, true, true);
        let f = x.sin() / x;
        let (s, c) = (1.3_f64.sin(), 1.3_f64.cos());
        // f = sin/x. f' = (x cos - sin)/x^2. f'' = (-x^2 sin - 2x cos + 2 sin)/x^3
        let x0 = 1.3;
        let fp = (x0 * c - s) / (x0 * x0);
        let fpp = (-x0 * x0 * s - 2.0 * x0 * c + 2.0 * s) / (x0 * x0 * x0);
        approx(f.v, s / x0, 1e-12);
        approx(f.d1, fp, 1e-10);
        approx(f.d2, fp, 1e-10);
        approx(f.d12, fpp, 1e-9);
    }

    #[test]
    fn sqrt_second_derivative() {
        // f(x)=sqrt(x); f'=1/(2 sqrt x); f''=-1/(4 x^{3/2}). At x=4: 2, .25, -1/32.
        let x = HyperDual::seed(4.0, true, true);
        let f = x.sqrt();
        approx(f.v, 2.0, 1e-12);
        approx(f.d1, 0.25, 1e-12);
        approx(f.d12, -1.0 / 32.0, 1e-12);
    }
}
