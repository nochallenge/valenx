//! Explicit forward-time centred-space (FTCS) solver for Fick's second
//! law in one dimension.
//!
//! ## Model
//!
//! Fick's second law,
//!
//! ```text
//!   dC/dt = D d2C/dx2,
//! ```
//!
//! is discretised on a uniform grid of `n` nodes spaced `dx` apart. The
//! Laplacian uses the three-point centred stencil and time advances by
//! an explicit Euler step:
//!
//! ```text
//!   C_i^{k+1} = C_i^k + r ( C_{i+1}^k - 2 C_i^k + C_{i-1}^k ),
//!   r = D dt / dx^2.
//! ```
//!
//! ## Stability
//!
//! The scheme is conditionally stable: the dimensionless diffusion
//! number `r` must satisfy `r <= 1/2`, i.e.
//!
//! ```text
//!   dt <= dx^2 / (2 D).
//! ```
//!
//! [`max_stable_dt`] returns the limit, [`is_stable_dt`] tests a
//! candidate, and [`Grid::stable_dt`] gives a safe step at a chosen
//! fraction of the limit. The stepping routines reject an over-large
//! step with [`DiffusionError::Unstable`] rather than silently
//! producing a blowing-up field.
//!
//! ## Boundaries
//!
//! Each end is one of:
//!
//! - [`Boundary::Dirichlet`] — the boundary node is pinned to a fixed
//!   concentration every step;
//! - [`Boundary::NoFlux`] — a zero-flux (reflecting, homogeneous
//!   Neumann) wall: the wall face simply carries no diffusive flux, so
//!   `dC/dx = 0` there.
//!
//! The step is written in **conservative flux form** — each interior
//! face's flux is added to one neighbour and subtracted from the other —
//! so a domain closed by [`Boundary::NoFlux`] on both ends conserves
//! total mass exactly (to round-off) under this scheme; see
//! [`Field::total_mass`] and the crate tests. (On a uniform interior the
//! flux form is algebraically identical to the centred three-point
//! Laplacian stencil above; the two differ only in how cleanly they
//! close the no-flux wall.)

use crate::error::{DiffusionError, Result};
use serde::{Deserialize, Serialize};

/// A uniform 1-D grid: `n` equally spaced nodes a distance `dx` apart.
///
/// Node `i` sits at physical position `x_0 + i * dx`. The grid carries
/// no field data itself; it is the geometry a [`Field`] is defined on
/// and the object that knows the explicit-scheme stability limit.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Grid {
    /// Number of nodes (at least two).
    n: usize,
    /// Spacing between adjacent nodes (strictly positive).
    dx: f64,
    /// Physical coordinate of node 0.
    x0: f64,
}

impl Grid {
    /// Build a grid of `n` nodes spaced `dx` apart with node 0 at the
    /// origin.
    ///
    /// # Errors
    ///
    /// [`DiffusionError::BadParameter`] if `n < 2` or `dx` is not finite
    /// and strictly positive.
    pub fn new(n: usize, dx: f64) -> Result<Self> {
        Self::with_origin(n, dx, 0.0)
    }

    /// Build a grid of `n` nodes spaced `dx` apart with node 0 at `x0`.
    ///
    /// # Errors
    ///
    /// [`DiffusionError::BadParameter`] if `n < 2`, if `dx` is not finite
    /// and strictly positive, or if `x0` is not finite.
    pub fn with_origin(n: usize, dx: f64, x0: f64) -> Result<Self> {
        if n < 2 {
            return Err(DiffusionError::bad_parameter(
                "n",
                "grid needs at least two nodes",
            ));
        }
        if !dx.is_finite() || dx <= 0.0 {
            return Err(DiffusionError::bad_parameter(
                "dx",
                "cell spacing must be finite and strictly positive",
            ));
        }
        if !x0.is_finite() {
            return Err(DiffusionError::bad_parameter("x0", "must be finite"));
        }
        Ok(Grid { n, dx, x0 })
    }

    /// Number of nodes.
    pub fn len(&self) -> usize {
        self.n
    }

    /// Always `false` — a [`Grid`] has at least two nodes by
    /// construction. Provided so Clippy's `len_without_is_empty` lint is
    /// satisfied and callers reading `.len()` have the companion method.
    pub fn is_empty(&self) -> bool {
        false
    }

    /// Node spacing.
    pub fn dx(&self) -> f64 {
        self.dx
    }

    /// Coordinate of node 0.
    pub fn x0(&self) -> f64 {
        self.x0
    }

    /// Total span of the domain, `(n - 1) * dx`.
    pub fn length(&self) -> f64 {
        (self.n as f64 - 1.0) * self.dx
    }

    /// Physical coordinate of node `i`.
    ///
    /// # Panics
    ///
    /// Panics if `i >= self.len()`.
    pub fn x(&self, i: usize) -> f64 {
        assert!(i < self.n, "node index {i} out of range 0..{}", self.n);
        self.x0 + self.dx * i as f64
    }

    /// The largest explicit time step that keeps the FTCS scheme stable
    /// for diffusion coefficient `d`: `dx^2 / (2 D)`.
    ///
    /// # Errors
    ///
    /// [`DiffusionError::BadParameter`] if `d` is not finite and
    /// strictly positive.
    pub fn max_stable_dt(&self, d: f64) -> Result<f64> {
        max_stable_dt(d, self.dx)
    }

    /// A safe time step: `safety * dx^2 / (2 D)`.
    ///
    /// `safety` must lie in `(0, 1]`; a value below `1` keeps the run
    /// comfortably inside the stability boundary.
    ///
    /// # Errors
    ///
    /// [`DiffusionError::BadParameter`] if `d` is not positive or
    /// `safety` is not in `(0, 1]`.
    pub fn stable_dt(&self, d: f64, safety: f64) -> Result<f64> {
        if !safety.is_finite() || safety <= 0.0 || safety > 1.0 {
            return Err(DiffusionError::bad_parameter(
                "safety",
                "must lie in (0, 1]",
            ));
        }
        Ok(safety * self.max_stable_dt(d)?)
    }
}

/// Which boundary condition closes one end of the domain.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum Boundary {
    /// A fixed-concentration (Dirichlet) wall: the boundary node is
    /// held at this value every step.
    Dirichlet(f64),
    /// A zero-flux (reflecting / homogeneous-Neumann) wall: no matter
    /// crosses, so `dC/dx = 0` there.
    NoFlux,
}

/// A discrete concentration field on a [`Grid`], plus the boundary
/// condition at each end.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Field {
    grid: Grid,
    /// Per-node concentration, length `grid.len()`.
    c: Vec<f64>,
    /// Condition at node 0.
    left: Boundary,
    /// Condition at node `n - 1`.
    right: Boundary,
}

impl Field {
    /// Build a field from an explicit per-node concentration vector and
    /// the two boundary conditions.
    ///
    /// If a boundary is [`Boundary::Dirichlet`], the corresponding end
    /// node is overwritten with the pinned value so the initial state is
    /// already consistent with the boundary.
    ///
    /// # Errors
    ///
    /// [`DiffusionError::DimensionMismatch`] if `c.len() != grid.len()`;
    /// [`DiffusionError::BadParameter`] if any sample (or Dirichlet
    /// value) is non-finite.
    pub fn new(grid: Grid, c: Vec<f64>, left: Boundary, right: Boundary) -> Result<Self> {
        if c.len() != grid.len() {
            return Err(DiffusionError::dimension(format!(
                "initial condition has {} values but the grid has {} nodes",
                c.len(),
                grid.len()
            )));
        }
        if c.iter().any(|v| !v.is_finite()) {
            return Err(DiffusionError::bad_parameter(
                "concentration",
                "every initial sample must be finite",
            ));
        }
        check_boundary_finite(left)?;
        check_boundary_finite(right)?;
        let mut field = Field {
            grid,
            c,
            left,
            right,
        };
        field.apply_dirichlet();
        Ok(field)
    }

    /// Build a field that is uniformly zero except for a single unit of
    /// mass deposited at node `spike` — a discrete approximation to an
    /// instantaneous point source. The spike value is `1 / dx` so the
    /// discrete integral (sum times `dx`) is exactly one.
    ///
    /// # Errors
    ///
    /// [`DiffusionError::BadParameter`] if `spike >= grid.len()`.
    pub fn point_source(grid: Grid, spike: usize, left: Boundary, right: Boundary) -> Result<Self> {
        if spike >= grid.len() {
            return Err(DiffusionError::bad_parameter(
                "spike",
                "spike index is outside the grid",
            ));
        }
        let mut c = vec![0.0; grid.len()];
        c[spike] = 1.0 / grid.dx();
        Field::new(grid, c, left, right)
    }

    /// The grid this field is defined on.
    pub fn grid(&self) -> &Grid {
        &self.grid
    }

    /// Read-only view of the per-node concentrations.
    pub fn values(&self) -> &[f64] {
        &self.c
    }

    /// The left-end boundary condition.
    pub fn left(&self) -> Boundary {
        self.left
    }

    /// The right-end boundary condition.
    pub fn right(&self) -> Boundary {
        self.right
    }

    /// Total mass on the grid, the midpoint/rectangle quadrature
    /// `dx * sum_i C_i`.
    ///
    /// Under two [`Boundary::NoFlux`] ends this quantity is invariant
    /// (to round-off) across [`step`](Self::step) / [`advance`](Self::advance).
    pub fn total_mass(&self) -> f64 {
        self.grid.dx() * self.c.iter().sum::<f64>()
    }

    /// The smallest and largest node concentrations, as `(min, max)`.
    ///
    /// Useful for the discrete maximum principle: under a stable step a
    /// no-source field never develops a new interior extremum, so the
    /// global range can only shrink (or hold) over time.
    pub fn min_max(&self) -> (f64, f64) {
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for &v in &self.c {
            if v < lo {
                lo = v;
            }
            if v > hi {
                hi = v;
            }
        }
        (lo, hi)
    }

    /// Advance the field by one explicit FTCS step of size `dt` for a
    /// medium of diffusion coefficient `d`.
    ///
    /// # Errors
    ///
    /// [`DiffusionError::BadParameter`] if `d` or `dt` is not finite and
    /// strictly positive; [`DiffusionError::Unstable`] if
    /// `dt > dx^2 / (2 D)`.
    pub fn step(&mut self, d: f64, dt: f64) -> Result<()> {
        if !d.is_finite() || d <= 0.0 {
            return Err(DiffusionError::bad_parameter(
                "D",
                "diffusion coefficient must be finite and strictly positive",
            ));
        }
        if !dt.is_finite() || dt <= 0.0 {
            return Err(DiffusionError::bad_parameter(
                "dt",
                "time step must be finite and strictly positive",
            ));
        }
        let dx = self.grid.dx();
        let dt_max = max_stable_dt(d, dx)?;
        // Allow a hair of slack for the dt == dt_max case under
        // floating-point rounding of dx*dx.
        if dt > dt_max * (1.0 + 1e-12) {
            return Err(DiffusionError::Unstable { dt, dt_max, dx, d });
        }

        let r = d * dt / (dx * dx);
        let n = self.c.len();
        let old = self.c.clone();

        // Conservative flux form: accumulate the dimensionless diffusive
        // flux `r * (C_{i+1} - C_i)` across each *interior* face once,
        // adding it to the lower node and subtracting it from the upper.
        // Because every interior face contributes equal-and-opposite
        // amounts to its two neighbours, the total of all increments is
        // exactly zero — so a domain with no flux through either wall
        // conserves `dx * sum(C)` to round-off. A `NoFlux` wall carries
        // no wall-face flux (it is simply omitted); a `Dirichlet` end
        // node is re-pinned afterwards, so its own increment is moot, but
        // it still participates correctly as the fixed neighbour of the
        // first interior node through the shared interior face.
        let mut delta = vec![0.0_f64; n];
        for i in 0..n - 1 {
            let flux = r * (old[i + 1] - old[i]);
            delta[i] += flux;
            delta[i + 1] -= flux;
        }
        for i in 0..n {
            self.c[i] = old[i] + delta[i];
        }

        self.apply_dirichlet();
        Ok(())
    }

    /// Advance by `steps` explicit steps of size `dt`.
    ///
    /// Equivalent to calling [`step`](Self::step) `steps` times; returns
    /// the total elapsed time `steps * dt`.
    ///
    /// # Errors
    ///
    /// Same conditions as [`step`](Self::step) (checked on the first
    /// step; the parameters do not change between steps).
    pub fn advance(&mut self, d: f64, dt: f64, steps: usize) -> Result<f64> {
        for _ in 0..steps {
            self.step(d, dt)?;
        }
        Ok(dt * steps as f64)
    }

    /// Re-pin any Dirichlet ends to their fixed value.
    fn apply_dirichlet(&mut self) {
        let n = self.c.len();
        if let Boundary::Dirichlet(v) = self.left {
            self.c[0] = v;
        }
        if let Boundary::Dirichlet(v) = self.right {
            self.c[n - 1] = v;
        }
    }
}

/// The largest explicit time step that keeps the FTCS scheme stable for
/// diffusion coefficient `d` on a grid of spacing `dx`: `dx^2 / (2 D)`.
///
/// # Errors
///
/// [`DiffusionError::BadParameter`] if `d` or `dx` is not finite and
/// strictly positive.
///
/// # Examples
///
/// ```
/// use valenx_diffusion::max_stable_dt;
///
/// // dx = 0.1, D = 1 -> dt_max = 0.01 / 2 = 0.005.
/// assert!((max_stable_dt(1.0, 0.1).unwrap() - 0.005).abs() < 1e-15);
/// ```
pub fn max_stable_dt(d: f64, dx: f64) -> Result<f64> {
    if !d.is_finite() || d <= 0.0 {
        return Err(DiffusionError::bad_parameter(
            "D",
            "diffusion coefficient must be finite and strictly positive",
        ));
    }
    if !dx.is_finite() || dx <= 0.0 {
        return Err(DiffusionError::bad_parameter(
            "dx",
            "cell spacing must be finite and strictly positive",
        ));
    }
    Ok(dx * dx / (2.0 * d))
}

/// Whether the candidate step `dt` is within the FTCS stability limit
/// for `(d, dx)` (i.e. `dt <= dx^2 / (2 D)`).
///
/// # Errors
///
/// [`DiffusionError::BadParameter`] if `d` or `dx` is invalid, or if
/// `dt` is non-finite.
///
/// # Examples
///
/// ```
/// use valenx_diffusion::is_stable_dt;
///
/// assert!(is_stable_dt(0.004, 1.0, 0.1).unwrap());   // below 0.005
/// assert!(!is_stable_dt(0.006, 1.0, 0.1).unwrap());  // above 0.005
/// ```
pub fn is_stable_dt(dt: f64, d: f64, dx: f64) -> Result<bool> {
    if !dt.is_finite() {
        return Err(DiffusionError::bad_parameter("dt", "must be finite"));
    }
    let dt_max = max_stable_dt(d, dx)?;
    Ok(dt <= dt_max * (1.0 + 1e-12))
}

/// Reject a Dirichlet boundary that carries a non-finite value.
fn check_boundary_finite(b: Boundary) -> Result<()> {
    if let Boundary::Dirichlet(v) = b {
        if !v.is_finite() {
            return Err(DiffusionError::bad_parameter(
                "boundary",
                "Dirichlet value must be finite",
            ));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    #[test]
    fn grid_geometry() {
        let g = Grid::with_origin(5, 0.25, 1.0).unwrap();
        assert_eq!(g.len(), 5);
        assert!(!g.is_empty());
        assert!((g.dx() - 0.25).abs() < EPS);
        assert!((g.x0() - 1.0).abs() < EPS);
        assert!((g.x(0) - 1.0).abs() < EPS);
        assert!((g.x(4) - 2.0).abs() < EPS);
        assert!((g.length() - 1.0).abs() < EPS);
    }

    #[test]
    fn grid_rejects_bad_inputs() {
        assert!(Grid::new(1, 1.0).is_err());
        assert!(Grid::new(0, 1.0).is_err());
        assert!(Grid::new(3, 0.0).is_err());
        assert!(Grid::new(3, -1.0).is_err());
        assert!(Grid::with_origin(3, 1.0, f64::NAN).is_err());
    }

    #[test]
    fn stability_limit_value() {
        // dt_max = dx^2 / (2 D).
        assert!((max_stable_dt(2.0, 0.5).unwrap() - 0.0625).abs() < EPS);
        let g = Grid::new(10, 0.5).unwrap();
        assert!((g.max_stable_dt(2.0).unwrap() - 0.0625).abs() < EPS);
        // stable_dt scales by the safety factor.
        assert!((g.stable_dt(2.0, 0.5).unwrap() - 0.03125).abs() < EPS);
        assert!(g.stable_dt(2.0, 0.0).is_err());
        assert!(g.stable_dt(2.0, 1.5).is_err());
    }

    #[test]
    fn is_stable_dt_boundary() {
        let (d, dx) = (1.0, 0.1);
        let dt_max = max_stable_dt(d, dx).unwrap();
        assert!(is_stable_dt(dt_max, d, dx).unwrap());
        assert!(is_stable_dt(0.5 * dt_max, d, dx).unwrap());
        assert!(!is_stable_dt(1.5 * dt_max, d, dx).unwrap());
    }

    #[test]
    fn step_rejects_unstable_dt() {
        let g = Grid::new(6, 0.1).unwrap();
        let mut f = Field::new(g, vec![0.0; 6], Boundary::NoFlux, Boundary::NoFlux).unwrap();
        let dt_max = g.max_stable_dt(1.0).unwrap();
        // Just over the limit -> Unstable.
        let err = f.step(1.0, dt_max * 2.0).unwrap_err();
        assert_eq!(err.code(), "diffusion.unstable");
        // Exactly at the limit -> fine.
        assert!(f.step(1.0, dt_max).is_ok());
    }

    #[test]
    fn dirichlet_pins_ends_on_construction() {
        let g = Grid::new(4, 1.0).unwrap();
        let f = Field::new(
            g,
            vec![9.0, 9.0, 9.0, 9.0],
            Boundary::Dirichlet(0.0),
            Boundary::Dirichlet(1.0),
        )
        .unwrap();
        assert!((f.values()[0] - 0.0).abs() < EPS);
        assert!((f.values()[3] - 1.0).abs() < EPS);
    }

    #[test]
    fn no_flux_conserves_mass() {
        // A lopsided initial profile in a closed (no-flux) box: the
        // total mass must be invariant to round-off across many steps.
        let g = Grid::new(21, 0.1).unwrap();
        let mut c = vec![0.0; 21];
        c[5] = 3.0;
        c[6] = 5.0;
        c[7] = 2.0;
        let mut f = Field::new(g, c, Boundary::NoFlux, Boundary::NoFlux).unwrap();
        let m0 = f.total_mass();
        let dt = g.stable_dt(1.0, 0.4).unwrap();
        for _ in 0..2_000 {
            f.step(1.0, dt).unwrap();
        }
        let m1 = f.total_mass();
        assert!((m1 - m0).abs() < 1e-9, "mass drifted: {m0} -> {m1}");
        // And it should have flattened toward a uniform value.
        let (lo, hi) = f.min_max();
        assert!(hi - lo < 0.05, "did not flatten: range {}", hi - lo);
        // The conservative scheme holds `sum(C)` fixed, so the uniform
        // equilibrium every node approaches is the node-arithmetic mean
        // `sum(C) / n = total_mass / (dx * n)`.
        let node_mean = m0 / (g.dx() * g.len() as f64);
        assert!(
            (f.values()[10] - node_mean).abs() < 0.02,
            "interior node {} did not reach the node mean {node_mean}",
            f.values()[10]
        );
    }

    #[test]
    fn no_flux_relaxes_to_flat_mean() {
        // Closed box, two-level start: equilibrium is the flat average.
        let g = Grid::new(11, 1.0).unwrap();
        let mut c = vec![0.0; 11];
        for ci in c.iter_mut().take(6) {
            *ci = 2.0;
        }
        let avg = c.iter().sum::<f64>() / c.len() as f64;
        let mut f = Field::new(g, c, Boundary::NoFlux, Boundary::NoFlux).unwrap();
        let dt = g.stable_dt(0.5, 0.5).unwrap();
        f.advance(0.5, dt, 20_000).unwrap();
        for (i, &v) in f.values().iter().enumerate() {
            assert!((v - avg).abs() < 1e-3, "node {i} = {v}, mean {avg}");
        }
    }

    #[test]
    fn maximum_principle_holds() {
        // Under a stable step a no-source field cannot exceed its
        // initial range (discrete maximum principle).
        let g = Grid::new(31, 0.1).unwrap();
        let mut c = vec![0.0; 31];
        c[15] = 1.0;
        let mut f = Field::new(g, c, Boundary::NoFlux, Boundary::NoFlux).unwrap();
        let (lo0, hi0) = f.min_max();
        let dt = g.stable_dt(1.0, 0.5).unwrap();
        for _ in 0..500 {
            f.step(1.0, dt).unwrap();
            let (lo, hi) = f.min_max();
            assert!(hi <= hi0 + 1e-12, "max grew to {hi}");
            assert!(lo >= lo0 - 1e-12, "min dropped to {lo}");
        }
    }

    #[test]
    fn dirichlet_step_reaches_linear_steady_state() {
        // Fixed 0 at the left, 10 at the right: the long-time solution
        // of dC/dt = D d2C/dx2 is the linear profile.
        let g = Grid::new(11, 1.0).unwrap();
        let mut f = Field::new(
            g,
            vec![0.0; 11],
            Boundary::Dirichlet(0.0),
            Boundary::Dirichlet(10.0),
        )
        .unwrap();
        let dt = g.stable_dt(1.0, 0.5).unwrap();
        f.advance(1.0, dt, 40_000).unwrap();
        for (i, &v) in f.values().iter().enumerate() {
            let expected = i as f64; // 0,1,...,10 across the 10-unit span
            assert!((v - expected).abs() < 1e-3, "node {i} = {v}, want {expected}");
        }
    }

    #[test]
    fn point_source_constructor_unit_mass() {
        let g = Grid::new(11, 0.2).unwrap();
        let f = Field::point_source(g, 5, Boundary::NoFlux, Boundary::NoFlux).unwrap();
        // Discrete integral dx * sum = 1.
        assert!((f.total_mass() - 1.0).abs() < EPS);
        assert!(Field::point_source(g, 99, Boundary::NoFlux, Boundary::NoFlux).is_err());
    }

    #[test]
    fn discrete_spreading_tracks_two_d_t() {
        // A unit point source in a wide closed box. While the cloud
        // stays away from the walls the discrete variance should grow
        // close to 2 D t (the analytic law) plus the initial-spike
        // discretisation offset; we check the *slope* dVar/dt ~ 2 D.
        let g = Grid::new(401, 0.05).unwrap();
        let mid = 200;
        let mut f = Field::point_source(g, mid, Boundary::NoFlux, Boundary::NoFlux).unwrap();
        let d = 1.0;
        let dt = g.stable_dt(d, 0.4).unwrap();

        // Second central moment of the discrete profile: the variance
        // of position weighted by concentration. The cell width dx is a
        // common factor of every term in both the moment sums and the
        // normalising mass, so it cancels and need not appear here.
        let var_of = |f: &Field| {
            let total: f64 = f.values().iter().sum::<f64>();
            let mean: f64 = f
                .values()
                .iter()
                .enumerate()
                .map(|(i, &c)| f.grid().x(i) * c)
                .sum::<f64>()
                / total;
            f.values()
                .iter()
                .enumerate()
                .map(|(i, &c)| {
                    let dxi = f.grid().x(i) - mean;
                    dxi * dxi * c
                })
                .sum::<f64>()
                / total
        };

        f.advance(d, dt, 200).unwrap();
        let t1 = dt * 200.0;
        let v1 = var_of(&f);
        f.advance(d, dt, 200).unwrap();
        let t2 = dt * 400.0;
        let v2 = var_of(&f);

        let slope = (v2 - v1) / (t2 - t1);
        // dVar/dt = 2 D exactly for the continuous heat equation; the
        // FTCS scheme reproduces it to a few percent here.
        assert!(
            (slope - 2.0 * d).abs() < 0.05,
            "dVar/dt = {slope}, want {}",
            2.0 * d
        );
    }

    #[test]
    fn field_rejects_bad_inputs() {
        let g = Grid::new(4, 1.0).unwrap();
        // Wrong length.
        assert!(Field::new(g, vec![0.0; 3], Boundary::NoFlux, Boundary::NoFlux).is_err());
        // Non-finite sample.
        assert!(Field::new(
            g,
            vec![0.0, f64::NAN, 0.0, 0.0],
            Boundary::NoFlux,
            Boundary::NoFlux
        )
        .is_err());
        // Non-finite Dirichlet value.
        assert!(Field::new(
            g,
            vec![0.0; 4],
            Boundary::Dirichlet(f64::INFINITY),
            Boundary::NoFlux
        )
        .is_err());
        let mut f = Field::new(g, vec![0.0; 4], Boundary::NoFlux, Boundary::NoFlux).unwrap();
        assert!(f.step(0.0, 0.1).is_err());
        assert!(f.step(1.0, 0.0).is_err());
        assert!(f.step(1.0, f64::NAN).is_err());
    }
}
