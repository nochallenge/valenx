//! Orbital and density evaluation on a 3-D grid.
//!
//! For visualisation a molecular orbital or the electron density must
//! be sampled on a regular Cartesian grid — the same data a Gaussian
//! `.cube` file holds. This module evaluates a basis function, a
//! molecular orbital, or the total density at any point, and fills a
//! [`VolumetricGrid`] over a box.
//!
//! ## The value of a function at a point
//!
//! A Cartesian basis function at `r` is
//!
//! ```text
//! φ(r) = (x-Ax)^i (y-Ay)^j (z-Az)^k Σ_p c_p exp(-α_p |r-A|²)
//! ```
//!
//! A molecular orbital is `ψ_n(r) = Σ_μ C_{μn} φ_μ(r)`; the closed-shell
//! density is `ρ(r) = Σ_{μν} D_{μν} φ_μ(r) φ_ν(r)`.
//!
//! ## v1 note
//!
//! The grid fill is a direct triple loop — `O(n_grid · n_basis)` — with
//! no distance screening of negligible Gaussian tails. Fine for the
//! modest grids used to preview an orbital; a production cube routine
//! would screen and block.

use crate::basis::{BasisFunction, BasisSet};
use nalgebra::DMatrix;

/// Evaluate one contracted Cartesian Gaussian basis function at `point`.
pub fn basis_function_value(f: &BasisFunction, point: [f64; 3]) -> f64 {
    let dx = point[0] - f.centre[0];
    let dy = point[1] - f.centre[1];
    let dz = point[2] - f.centre[2];
    let r2 = dx * dx + dy * dy + dz * dz;
    let mut radial = 0.0;
    for p in &f.primitives {
        radial += p.coefficient * (-p.exponent * r2).exp();
    }
    let angular = dx.powi(f.cart.0 as i32) * dy.powi(f.cart.1 as i32) * dz.powi(f.cart.2 as i32);
    angular * radial
}

/// Evaluate molecular orbital `mo_index` at `point`:
/// `ψ(r) = Σ_μ C_{μ,mo_index} φ_μ(r)`.
pub fn orbital_value(
    basis: &BasisSet,
    coefficients: &DMatrix<f64>,
    mo_index: usize,
    point: [f64; 3],
) -> f64 {
    let mut psi = 0.0;
    for (mu, f) in basis.functions.iter().enumerate() {
        psi += coefficients[(mu, mo_index)] * basis_function_value(f, point);
    }
    psi
}

/// Evaluate the closed-shell electron density at `point`:
/// `ρ(r) = Σ_{μν} D_{μν} φ_μ(r) φ_ν(r)`.
pub fn density_value(basis: &BasisSet, density: &DMatrix<f64>, point: [f64; 3]) -> f64 {
    let phi: Vec<f64> = basis
        .functions
        .iter()
        .map(|f| basis_function_value(f, point))
        .collect();
    let n = phi.len();
    let mut rho = 0.0;
    for mu in 0..n {
        for nu in 0..n {
            rho += density[(mu, nu)] * phi[mu] * phi[nu];
        }
    }
    rho
}

/// A scalar field sampled on a regular Cartesian grid — `.cube`-style
/// volumetric data.
#[derive(Clone, Debug)]
pub struct VolumetricGrid {
    /// Grid origin (the lowest-corner point) in bohr.
    pub origin: [f64; 3],
    /// Number of points along each axis `[nx, ny, nz]`.
    pub counts: [usize; 3],
    /// Spacing between points along each axis (bohr).
    pub spacing: [f64; 3],
    /// Sampled values in `x`-fastest (then `y`, then `z`) order.
    pub values: Vec<f64>,
}

impl VolumetricGrid {
    /// Total number of grid points.
    pub fn n_points(&self) -> usize {
        self.counts[0] * self.counts[1] * self.counts[2]
    }

    /// The value at grid index `(ix, iy, iz)`.
    pub fn at(&self, ix: usize, iy: usize, iz: usize) -> f64 {
        self.values[(iz * self.counts[1] + iy) * self.counts[0] + ix]
    }

    /// The Cartesian coordinates of grid index `(ix, iy, iz)`.
    pub fn coordinate(&self, ix: usize, iy: usize, iz: usize) -> [f64; 3] {
        [
            self.origin[0] + ix as f64 * self.spacing[0],
            self.origin[1] + iy as f64 * self.spacing[1],
            self.origin[2] + iz as f64 * self.spacing[2],
        ]
    }
}

/// Sample an arbitrary scalar function over a regular grid.
fn fill_grid<F>(origin: [f64; 3], counts: [usize; 3], spacing: [f64; 3], f: F) -> VolumetricGrid
where
    F: Fn([f64; 3]) -> f64,
{
    let mut values = Vec::with_capacity(counts[0] * counts[1] * counts[2]);
    for iz in 0..counts[2] {
        for iy in 0..counts[1] {
            for ix in 0..counts[0] {
                let p = [
                    origin[0] + ix as f64 * spacing[0],
                    origin[1] + iy as f64 * spacing[1],
                    origin[2] + iz as f64 * spacing[2],
                ];
                values.push(f(p));
            }
        }
    }
    VolumetricGrid {
        origin,
        counts,
        spacing,
        values,
    }
}

/// Sample molecular orbital `mo_index` over a regular grid.
pub fn orbital_grid(
    basis: &BasisSet,
    coefficients: &DMatrix<f64>,
    mo_index: usize,
    origin: [f64; 3],
    counts: [usize; 3],
    spacing: [f64; 3],
) -> VolumetricGrid {
    fill_grid(origin, counts, spacing, |p| {
        orbital_value(basis, coefficients, mo_index, p)
    })
}

/// Sample the electron density over a regular grid.
pub fn density_grid(
    basis: &BasisSet,
    density: &DMatrix<f64>,
    origin: [f64; 3],
    counts: [usize; 3],
    spacing: [f64; 3],
) -> VolumetricGrid {
    fill_grid(origin, counts, spacing, |p| {
        density_value(basis, density, p)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{Atom, MolecularGeometry};
    use crate::integrals::IntegralSet;
    use crate::scf::rhf::{run_rhf_scf, ScfSettings};

    fn h2() -> (MolecularGeometry, BasisSet) {
        let geom = MolecularGeometry::new(vec![
            Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.0]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.7414]).unwrap(),
        ]);
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        (geom, basis)
    }

    #[test]
    fn basis_function_peaks_at_its_centre() {
        let (_, basis) = h2();
        let f = &basis.functions[0];
        let at_centre = basis_function_value(f, f.centre);
        let away = basis_function_value(f, [0.0, 0.0, 5.0]);
        assert!(at_centre > away.abs());
        assert!(at_centre > 0.0);
    }

    #[test]
    fn h2_bonding_orbital_is_nonzero_between_atoms() {
        let (geom, basis) = h2();
        let ints = IntegralSet::compute(&geom, &basis);
        let res = run_rhf_scf(&ints, 2, ScfSettings::default()).unwrap();
        // Midpoint of the bond, in bohr.
        let mid = 0.7414 * crate::geometry::BOHR_PER_ANGSTROM / 2.0;
        let psi = orbital_value(&basis, &res.orbital_coefficients, 0, [0.0, 0.0, mid]);
        assert!(psi.abs() > 0.1, "bonding MO at midpoint = {psi}");
    }

    #[test]
    fn density_integrates_roughly_to_electron_count() {
        // A coarse Riemann sum of ρ over a generous box should land
        // near the 2 electrons of H2.
        let (geom, basis) = h2();
        let ints = IntegralSet::compute(&geom, &basis);
        let res = run_rhf_scf(&ints, 2, ScfSettings::default()).unwrap();
        let n = 40usize;
        let span = 12.0;
        let h = span / n as f64;
        let origin = [-span / 2.0, -span / 2.0, -span / 2.0 + 0.7];
        let grid = density_grid(&basis, &res.density, origin, [n, n, n], [h, h, h]);
        let integral: f64 = grid.values.iter().sum::<f64>() * h * h * h;
        assert!(
            (integral - 2.0).abs() < 0.2,
            "∫ρ ≈ {integral} (expected ~2)"
        );
    }

    #[test]
    fn grid_indexing_round_trips() {
        let (geom, basis) = h2();
        let ints = IntegralSet::compute(&geom, &basis);
        let res = run_rhf_scf(&ints, 2, ScfSettings::default()).unwrap();
        let grid = orbital_grid(
            &basis,
            &res.orbital_coefficients,
            0,
            [0.0, 0.0, 0.0],
            [3, 4, 5],
            [0.5, 0.5, 0.5],
        );
        assert_eq!(grid.n_points(), 60);
        let c = grid.coordinate(2, 3, 4);
        assert!((c[0] - 1.0).abs() < 1.0e-12);
        assert!((c[2] - 2.0).abs() < 1.0e-12);
    }
}
