//! A minimal `f64` complex number, `(re, im)`.
//!
//! Deliberately tiny and dependency-free: the crate's promise is "no
//! external complex crate". Only the handful of operations the DFT
//! needs are provided, exposed through the standard [`std::ops`]
//! arithmetic traits (`Add`, `Mul`, plus scalar `Mul<f64>`) so call
//! sites read like ordinary algebra.

use serde::{Deserialize, Serialize};
use std::ops::{Add, AddAssign, Mul};

/// A complex number stored as a real / imaginary `f64` pair,
/// `value = re + im i`.
#[derive(Copy, Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Complex {
    /// Real part.
    pub re: f64,
    /// Imaginary part.
    pub im: f64,
}

impl Complex {
    /// The additive identity, `0 + 0 i`.
    pub const ZERO: Complex = Complex { re: 0.0, im: 0.0 };

    /// Construct from real and imaginary parts.
    #[inline]
    pub const fn new(re: f64, im: f64) -> Self {
        Complex { re, im }
    }

    /// A purely real value, `re + 0 i`.
    #[inline]
    pub const fn real(re: f64) -> Self {
        Complex { re, im: 0.0 }
    }

    /// `e^{i theta} = cos(theta) + i sin(theta)`, the point on the unit
    /// circle at angle `theta` radians.
    #[inline]
    pub fn expi(theta: f64) -> Self {
        Complex {
            re: theta.cos(),
            im: theta.sin(),
        }
    }

    /// Complex conjugate, `re - im i`.
    #[inline]
    pub fn conj(self) -> Complex {
        Complex {
            re: self.re,
            im: -self.im,
        }
    }

    /// Squared magnitude, `re^2 + im^2`. Cheaper than [`Complex::norm`]
    /// (no square root) and handy for energy / Parseval sums.
    #[inline]
    pub fn norm_sqr(self) -> f64 {
        self.re * self.re + self.im * self.im
    }

    /// Magnitude (modulus), `sqrt(re^2 + im^2)`.
    #[inline]
    pub fn norm(self) -> f64 {
        self.norm_sqr().sqrt()
    }
}

impl Add for Complex {
    type Output = Complex;

    /// Complex addition, `(a + b)`.
    #[inline]
    fn add(self, other: Complex) -> Complex {
        Complex {
            re: self.re + other.re,
            im: self.im + other.im,
        }
    }
}

impl AddAssign for Complex {
    /// In-place complex addition, `a += b`.
    #[inline]
    fn add_assign(&mut self, other: Complex) {
        self.re += other.re;
        self.im += other.im;
    }
}

impl Mul for Complex {
    type Output = Complex;

    /// Complex multiplication, `(a * b)`.
    #[inline]
    fn mul(self, other: Complex) -> Complex {
        Complex {
            re: self.re * other.re - self.im * other.im,
            im: self.re * other.im + self.im * other.re,
        }
    }
}

impl Mul<f64> for Complex {
    type Output = Complex;

    /// Scale by a real scalar, `(a * s)`.
    #[inline]
    fn mul(self, s: f64) -> Complex {
        Complex {
            re: self.re * s,
            im: self.im * s,
        }
    }
}

impl From<f64> for Complex {
    /// Lift a real `f64` into the complex plane with zero imaginary
    /// part.
    #[inline]
    fn from(re: f64) -> Self {
        Complex::real(re)
    }
}
