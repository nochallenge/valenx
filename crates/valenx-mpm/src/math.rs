//! Minimal 2D linear-algebra primitives for the MPM solver.
//!
//! Kept dependency-free and inline so the solver core is fully deterministic
//! and self-contained. Only the operations the Material Point Method needs are
//! provided: 2-vectors and 2x2 matrices with the products used in the
//! particle-to-grid / grid-to-particle transfers and the elastic stress.

/// A 2D column vector `(x, y)` of `f64` components.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Vec2 {
    /// First (x) component.
    pub x: f64,
    /// Second (y) component.
    pub y: f64,
}

impl Vec2 {
    /// Constructs a vector from its components.
    #[inline]
    #[must_use]
    pub const fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    /// The zero vector `(0, 0)`.
    #[inline]
    #[must_use]
    pub const fn zero() -> Self {
        Self { x: 0.0, y: 0.0 }
    }

    /// Returns `true` iff both components are finite (no `NaN`/`inf`).
    #[inline]
    #[must_use]
    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite()
    }

    /// Euclidean dot product `self · rhs`.
    #[inline]
    #[must_use]
    pub fn dot(self, rhs: Self) -> f64 {
        self.x * rhs.x + self.y * rhs.y
    }

    /// Squared Euclidean length `‖self‖²`.
    #[inline]
    #[must_use]
    pub fn length_squared(self) -> f64 {
        self.dot(self)
    }

    /// Outer product `self ⊗ rhs`, the 2x2 matrix `self · rhsᵀ`.
    #[inline]
    #[must_use]
    pub fn outer(self, rhs: Self) -> Mat2 {
        Mat2::new(
            self.x * rhs.x,
            self.x * rhs.y,
            self.y * rhs.x,
            self.y * rhs.y,
        )
    }
}

impl core::ops::Add for Vec2 {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self::new(self.x + rhs.x, self.y + rhs.y)
    }
}

impl core::ops::Sub for Vec2 {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.x - rhs.x, self.y - rhs.y)
    }
}

impl core::ops::Mul<f64> for Vec2 {
    type Output = Self;
    #[inline]
    fn mul(self, s: f64) -> Self {
        Self::new(self.x * s, self.y * s)
    }
}

impl core::ops::AddAssign for Vec2 {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.x += rhs.x;
        self.y += rhs.y;
    }
}

/// A 2x2 matrix stored row-major as `[[m00, m01], [m10, m11]]`.
///
/// Used for the affine velocity field `C` (APIC), the deformation gradient
/// `F`, and the Cauchy/PK1 elastic stress tensors.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Mat2 {
    /// Row 0, column 0.
    pub m00: f64,
    /// Row 0, column 1.
    pub m01: f64,
    /// Row 1, column 0.
    pub m10: f64,
    /// Row 1, column 1.
    pub m11: f64,
}

impl Mat2 {
    /// Constructs a matrix from its four row-major entries.
    #[inline]
    #[must_use]
    pub const fn new(m00: f64, m01: f64, m10: f64, m11: f64) -> Self {
        Self { m00, m01, m10, m11 }
    }

    /// The 2x2 zero matrix.
    #[inline]
    #[must_use]
    pub const fn zero() -> Self {
        Self::new(0.0, 0.0, 0.0, 0.0)
    }

    /// The 2x2 identity matrix.
    #[inline]
    #[must_use]
    pub const fn identity() -> Self {
        Self::new(1.0, 0.0, 0.0, 1.0)
    }

    /// Returns `true` iff all entries are finite.
    #[inline]
    #[must_use]
    pub fn is_finite(self) -> bool {
        self.m00.is_finite() && self.m01.is_finite() && self.m10.is_finite() && self.m11.is_finite()
    }

    /// Determinant `det(self)`.
    #[inline]
    #[must_use]
    pub fn determinant(self) -> f64 {
        self.m00 * self.m11 - self.m01 * self.m10
    }

    /// Transpose `selfᵀ`.
    #[inline]
    #[must_use]
    pub fn transpose(self) -> Self {
        Self::new(self.m00, self.m10, self.m01, self.m11)
    }

    /// Matrix-vector product `self · v`.
    #[inline]
    #[must_use]
    pub fn mul_vec(self, v: Vec2) -> Vec2 {
        Vec2::new(
            self.m00 * v.x + self.m01 * v.y,
            self.m10 * v.x + self.m11 * v.y,
        )
    }

    /// Matrix-matrix product `self · rhs`.
    #[inline]
    #[must_use]
    pub fn mul_mat(self, rhs: Self) -> Self {
        Self::new(
            self.m00 * rhs.m00 + self.m01 * rhs.m10,
            self.m00 * rhs.m01 + self.m01 * rhs.m11,
            self.m10 * rhs.m00 + self.m11 * rhs.m10,
            self.m10 * rhs.m01 + self.m11 * rhs.m11,
        )
    }

    /// Scales every entry by `s`.
    #[inline]
    #[must_use]
    pub fn scale(self, s: f64) -> Self {
        Self::new(self.m00 * s, self.m01 * s, self.m10 * s, self.m11 * s)
    }

    /// Entry-wise sum `self + rhs`.
    #[inline]
    #[must_use]
    pub fn plus(self, rhs: Self) -> Self {
        Self::new(
            self.m00 + rhs.m00,
            self.m01 + rhs.m01,
            self.m10 + rhs.m10,
            self.m11 + rhs.m11,
        )
    }

    /// Entry-wise difference `self - rhs`.
    #[inline]
    #[must_use]
    pub fn minus(self, rhs: Self) -> Self {
        Self::new(
            self.m00 - rhs.m00,
            self.m01 - rhs.m01,
            self.m10 - rhs.m10,
            self.m11 - rhs.m11,
        )
    }

    /// Polar decomposition `self = R · S`, returning the rotation `R`.
    ///
    /// For a 2x2 matrix the rotation factor has a closed form: with
    /// `c = m00 + m11`, `s = m10 - m01`, and `d = √(c² + s²)`, the rotation is
    /// `R = [[c, -s], [s, c]] / d`. When `self` is (numerically) singular the
    /// identity is returned. Used by the fixed-corotated elastic model.
    #[inline]
    #[must_use]
    pub fn polar_rotation(self) -> Self {
        let c = self.m00 + self.m11;
        let s = self.m10 - self.m01;
        let d = (c * c + s * s).sqrt();
        if d <= f64::EPSILON {
            return Self::identity();
        }
        let inv = 1.0 / d;
        Self::new(c * inv, -s * inv, s * inv, c * inv)
    }
}
