//! Cross-validation of the explicit FTCS solver against the closed-form
//! analytic results.
//!
//! These exercise the whole crate end to end:
//!
//! 1. the transient FTCS solution of a point source converges to the
//!    analytic Gaussian heat kernel;
//! 2. the long-time Dirichlet solution converges to the analytic linear
//!    steady state, and its through-domain flux matches Fick's first law;
//! 3. a closed (no-flux) domain conserves mass exactly under the scheme.

use valenx_diffusion::{gaussian_point_source, steady_flux, steady_profile, Boundary, Field, Grid};

const EPS: f64 = 1e-12;

/// The FTCS evolution of a discrete point source matches the analytic
/// Gaussian `C(x, t) = M / sqrt(4 pi D t) exp(-x^2 / 4 D t)` away from
/// the walls.
#[test]
fn ftcs_point_source_matches_gaussian() {
    // Fine, wide grid so the cloud stays clear of the boundaries.
    let dx = 0.05;
    let n = 401;
    let mid = 200;
    let grid = Grid::new(n, dx).unwrap();
    let d = 1.0;

    // Unit-mass spike: value 1/dx at one node (discrete integral = 1).
    let mut field = Field::point_source(grid, mid, Boundary::NoFlux, Boundary::NoFlux).unwrap();

    // March to a time at which the spike has spread over many cells but
    // not reached the walls. sigma = sqrt(2 D t); pick t so sigma ~ 1.0.
    let dt = grid.stable_dt(d, 0.4).unwrap();
    let target_t = 0.5; // sigma = sqrt(2*1*0.5) = 1.0
    let steps = (target_t / dt).round() as usize;
    let elapsed = field.advance(d, dt, steps).unwrap();

    let x_center = grid.x(mid);
    let mass = 1.0;

    // Compare node-by-node over the central region (|x - x0| <= 3 sigma).
    let sigma = (2.0 * d * elapsed).sqrt();
    let mut max_abs_err = 0.0_f64;
    for i in 0..n {
        let x = grid.x(i);
        if (x - x_center).abs() > 3.0 * sigma {
            continue;
        }
        let analytic = gaussian_point_source(mass, d, x - x_center, elapsed).unwrap();
        let err = (field.values()[i] - analytic).abs();
        if err > max_abs_err {
            max_abs_err = err;
        }
    }
    // The peak analytic value is ~M/sqrt(4 pi D t) ~ 0.4; a few-percent
    // absolute agreement across the cloud confirms the scheme solves the
    // right equation.
    assert!(
        max_abs_err < 0.01,
        "FTCS vs Gaussian max abs error = {max_abs_err}"
    );
}

/// The long-time Dirichlet solution converges to the analytic linear
/// steady profile, and the steady flux follows Fick's first law.
#[test]
fn ftcs_converges_to_linear_steady_state() {
    let grid = Grid::new(21, 0.5).unwrap(); // length = 10
    let (c_lo, c_hi) = (2.0, 12.0);
    let mut field = Field::new(
        grid,
        vec![c_lo; 21],
        Boundary::Dirichlet(c_lo),
        Boundary::Dirichlet(c_hi),
    )
    .unwrap();

    let d = 1.0;
    let dt = grid.stable_dt(d, 0.5).unwrap();
    field.advance(d, dt, 60_000).unwrap();

    let analytic = steady_profile(&grid, c_lo, c_hi).unwrap();
    for (i, (&num, &exact)) in field.values().iter().zip(analytic.iter()).enumerate() {
        assert!(
            (num - exact).abs() < 1e-3,
            "node {i}: numeric {num} vs analytic {exact}"
        );
    }

    // The flux implied by the converged gradient equals the closed-form
    // steady flux -D (c_hi - c_lo) / L.
    let length = grid.length();
    let j_analytic = steady_flux(d, c_lo, c_hi, length).unwrap();
    // Estimate the numeric flux from a central difference at mid-domain.
    let i = 10;
    let grad_num = (field.values()[i + 1] - field.values()[i - 1]) / (2.0 * grid.dx());
    let j_num = -d * grad_num;
    assert!(
        (j_num - j_analytic).abs() < 1e-3,
        "flux numeric {j_num} vs analytic {j_analytic}"
    );
    // High wall on the right -> flux toward -x.
    assert!(j_analytic < 0.0);
}

/// A closed (no-flux) domain conserves the discrete mass exactly under
/// many FTCS steps, regardless of the initial profile.
#[test]
fn closed_domain_conserves_mass_exactly() {
    let grid = Grid::new(41, 0.25).unwrap();
    let mut c = vec![0.0; 41];
    // An asymmetric blob.
    for (i, ci) in c.iter_mut().enumerate() {
        *ci = ((i as f64 - 8.0).powi(2) / -10.0).exp() + 0.3 * ((i as f64 - 30.0) / 4.0).cos();
    }
    let mut field = Field::new(grid, c, Boundary::NoFlux, Boundary::NoFlux).unwrap();
    let m0 = field.total_mass();

    let d = 0.7;
    let dt = grid.stable_dt(d, 0.45).unwrap();
    for _ in 0..5_000 {
        field.step(d, dt).unwrap();
        // Mass holds at *every* step, not just at the end.
        assert!(
            (field.total_mass() - m0).abs() < 1e-8,
            "mass drift: {} vs {m0}",
            field.total_mass()
        );
    }
    assert!((field.total_mass() - m0).abs() < EPS.max(1e-9));
}
