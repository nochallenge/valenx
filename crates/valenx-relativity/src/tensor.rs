//! Small fixed-size 4-D tensor helpers used by the curvature engine.
//!
//! Spacetime is four-dimensional, so every tensor here is backed by a plain
//! fixed-size array (`[f64; 4]`, `[[f64; 4]; 4]`, …) rather than a dynamic
//! matrix type. The only non-trivial linear-algebra operation we need is the
//! inverse of the metric, which we delegate to `nalgebra`.

use nalgebra::Matrix4;

/// A contravariant 4-vector or a coordinate point `(t, r, θ, φ)`.
pub type Vec4 = [f64; 4];
/// A rank-2 tensor / matrix in the coordinate basis (e.g. the metric `g_μν`).
pub type Mat4 = [[f64; 4]; 4];

/// Invert a 4×4 matrix, returning `None` if it is (numerically) singular.
///
/// Used to obtain the inverse metric `g^{μν}` from `g_{μν}`; a `None` result
/// signals a coordinate singularity (e.g. on a horizon) and is surfaced as an
/// error rather than producing NaNs downstream.
pub fn inverse(m: &Mat4) -> Option<Mat4> {
    let nm = Matrix4::from_fn(|i, j| m[i][j]);
    nm.try_inverse().map(|inv| {
        let mut out = [[0.0_f64; 4]; 4];
        for (i, row) in out.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                *cell = inv[(i, j)];
            }
        }
        out
    })
}

/// True iff every component of the matrix is finite (no NaN/∞).
pub fn all_finite(m: &Mat4) -> bool {
    m.iter().all(|row| row.iter().all(|x| x.is_finite()))
}
