//! Steady-state 1-D diffusion between two fixed-concentration walls.
//!
//! ## Model
//!
//! At steady state Fick's second law reduces to Laplace's equation in
//! one dimension,
//!
//! ```text
//!   d2C/dx2 = 0,
//! ```
//!
//! whose only solution between two Dirichlet walls is the straight line
//! joining them. For a slab of length `L` with `C(0) = c_left` and
//! `C(L) = c_right`,
//!
//! ```text
//!   C(x) = c_left + (c_right - c_left) * x / L,
//! ```
//!
//! a constant gradient `(c_right - c_left) / L` and hence (by Fick's
//! first law) a constant flux `J = -D (c_right - c_left) / L` everywhere
//! in the slab. [`steady_profile`] samples the line on a grid,
//! [`steady_gradient`] returns the gradient, and [`steady_flux`] returns
//! the through-slab flux.

use crate::error::{DiffusionError, Result};
use crate::explicit::Grid;

/// The steady-state concentration at a point a fraction `frac` of the
/// way from the left wall (`frac = 0`) to the right wall (`frac = 1`).
///
/// Returns the linear interpolant `c_left + (c_right - c_left) * frac`.
/// `frac` may lie outside `[0, 1]` (linear extrapolation) as long as it
/// is finite.
///
/// # Errors
///
/// [`DiffusionError::BadParameter`] if any argument is non-finite.
///
/// # Examples
///
/// ```
/// use valenx_diffusion::steady_value;
///
/// // Halfway between walls at 0 and 10 sits at 5.
/// assert!((steady_value(0.0, 10.0, 0.5).unwrap() - 5.0).abs() < 1e-12);
/// ```
pub fn steady_value(c_left: f64, c_right: f64, frac: f64) -> Result<f64> {
    if !c_left.is_finite() || !c_right.is_finite() || !frac.is_finite() {
        return Err(DiffusionError::bad_parameter(
            "steady",
            "boundary values and fraction must be finite",
        ));
    }
    Ok(c_left + (c_right - c_left) * frac)
}

/// Sample the steady-state linear profile on `grid`, pinning the first
/// node to `c_left` and the last node to `c_right`.
///
/// The returned vector has length `grid.len()` and its values rise (or
/// fall) linearly across the domain.
///
/// # Errors
///
/// [`DiffusionError::BadParameter`] if either boundary value is
/// non-finite.
///
/// # Examples
///
/// ```
/// use valenx_diffusion::{steady_profile, Grid};
///
/// let g = Grid::new(5, 1.0).unwrap();
/// let c = steady_profile(&g, 0.0, 8.0).unwrap();
/// // Evenly spaced 0, 2, 4, 6, 8 across the 4-unit span.
/// assert!((c[2] - 4.0).abs() < 1e-12);
/// ```
pub fn steady_profile(grid: &Grid, c_left: f64, c_right: f64) -> Result<Vec<f64>> {
    if !c_left.is_finite() || !c_right.is_finite() {
        return Err(DiffusionError::bad_parameter(
            "boundary",
            "fixed-wall concentrations must be finite",
        ));
    }
    let last = grid.len() - 1;
    let denom = last as f64;
    let mut c = Vec::with_capacity(grid.len());
    for i in 0..grid.len() {
        let frac = i as f64 / denom;
        c.push(c_left + (c_right - c_left) * frac);
    }
    Ok(c)
}

/// The constant concentration gradient of the steady profile across a
/// slab of length `length`: `(c_right - c_left) / length`.
///
/// # Errors
///
/// [`DiffusionError::BadParameter`] if `length` is not finite and
/// strictly positive, or if either boundary value is non-finite.
///
/// # Examples
///
/// ```
/// use valenx_diffusion::steady_gradient;
///
/// // From 0 to 10 over a length of 5: gradient 2.
/// assert!((steady_gradient(0.0, 10.0, 5.0).unwrap() - 2.0).abs() < 1e-12);
/// ```
pub fn steady_gradient(c_left: f64, c_right: f64, length: f64) -> Result<f64> {
    if !c_left.is_finite() || !c_right.is_finite() {
        return Err(DiffusionError::bad_parameter(
            "boundary",
            "fixed-wall concentrations must be finite",
        ));
    }
    if !length.is_finite() || length <= 0.0 {
        return Err(DiffusionError::bad_parameter(
            "length",
            "slab length must be finite and strictly positive",
        ));
    }
    Ok((c_right - c_left) / length)
}

/// The constant diffusive flux through a steady slab of length `length`
/// and diffusion coefficient `d`: `J = -D (c_right - c_left) / length`.
///
/// Combines Fick's first law with the steady linear gradient. The flux
/// is uniform across the slab and points from the high-concentration
/// wall to the low one.
///
/// # Errors
///
/// [`DiffusionError::BadParameter`] if `d` is not finite and strictly
/// positive, if `length` is not finite and strictly positive, or if
/// either boundary value is non-finite.
///
/// # Examples
///
/// ```
/// use valenx_diffusion::steady_flux;
///
/// // High at the left wall drives flux to the right (positive x): J > 0.
/// let j = steady_flux(2.0, 10.0, 0.0, 5.0).unwrap();
/// assert!(j > 0.0);
/// assert!((j - 4.0).abs() < 1e-12);
/// ```
pub fn steady_flux(d: f64, c_left: f64, c_right: f64, length: f64) -> Result<f64> {
    if !d.is_finite() || d <= 0.0 {
        return Err(DiffusionError::bad_parameter(
            "D",
            "diffusion coefficient must be finite and strictly positive",
        ));
    }
    let grad = steady_gradient(c_left, c_right, length)?;
    Ok(-d * grad)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    #[test]
    fn profile_is_linear_and_pins_ends() {
        let g = Grid::new(11, 1.0).unwrap();
        let c = steady_profile(&g, 0.0, 10.0).unwrap();
        assert_eq!(c.len(), 11);
        // Ends pinned exactly.
        assert!((c[0] - 0.0).abs() < EPS);
        assert!((c[10] - 10.0).abs() < EPS);
        // Each node equals its index (0..10 over a 10-unit span).
        for (i, &v) in c.iter().enumerate() {
            assert!((v - i as f64).abs() < EPS, "node {i} = {v}");
        }
        // Constant spacing between consecutive nodes (second difference 0).
        for i in 1..10 {
            let second = c[i + 1] - 2.0 * c[i] + c[i - 1];
            assert!(second.abs() < EPS, "curvature at {i} = {second}");
        }
    }

    #[test]
    fn descending_profile() {
        let g = Grid::new(5, 0.5).unwrap();
        let c = steady_profile(&g, 8.0, 0.0).unwrap();
        assert!((c[0] - 8.0).abs() < EPS);
        assert!((c[4] - 0.0).abs() < EPS);
        assert!((c[2] - 4.0).abs() < EPS); // midpoint
    }

    #[test]
    fn gradient_and_value_are_consistent() {
        // value at fraction f == c_left + grad * (f * length)
        let (cl, cr, len) = (1.0, 7.0, 3.0);
        let grad = steady_gradient(cl, cr, len).unwrap();
        assert!((grad - 2.0).abs() < EPS, "grad {grad}");
        let mid = steady_value(cl, cr, 0.5).unwrap();
        assert!((mid - (cl + grad * 0.5 * len)).abs() < EPS, "mid {mid}");
    }

    #[test]
    fn flux_is_uniform_and_down_gradient() {
        // Left wall high -> matter flows toward +x -> J > 0.
        let j = steady_flux(2.0, 10.0, 0.0, 5.0).unwrap();
        assert!(j > 0.0);
        assert!((j - 4.0).abs() < EPS, "j {j}");

        // Right wall high -> matter flows toward -x -> J < 0.
        let j = steady_flux(2.0, 0.0, 10.0, 5.0).unwrap();
        assert!(j < 0.0);
        assert!((j - (-4.0)).abs() < EPS, "j {j}");

        // Equal walls -> no gradient -> no flux.
        let j = steady_flux(2.0, 3.0, 3.0, 5.0).unwrap();
        assert!(j.abs() < EPS, "j {j}");
    }

    #[test]
    fn steady_flux_matches_grid_length() {
        // The slab flux computed from the grid length agrees with the
        // explicit-length form.
        let g = Grid::new(6, 0.4).unwrap(); // length = 5 * 0.4 = 2.0
        let len = g.length();
        let j_grid = steady_flux(1.5, 4.0, 0.0, len).unwrap();
        let j_direct = steady_flux(1.5, 4.0, 0.0, 2.0).unwrap();
        assert!((j_grid - j_direct).abs() < EPS, "{j_grid} {j_direct}");
    }

    #[test]
    fn rejects_bad_inputs() {
        assert!(steady_value(f64::NAN, 1.0, 0.5).is_err());
        let g = Grid::new(4, 1.0).unwrap();
        assert!(steady_profile(&g, f64::INFINITY, 1.0).is_err());
        assert!(steady_gradient(0.0, 1.0, 0.0).is_err());
        assert!(steady_gradient(0.0, 1.0, -1.0).is_err());
        assert!(steady_flux(0.0, 0.0, 1.0, 1.0).is_err());
        assert!(steady_flux(1.0, 0.0, 1.0, 0.0).is_err());
    }
}
