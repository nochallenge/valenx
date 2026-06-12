//! The atom-centred molecular integration grid.
//!
//! This module assembles the full molecular quadrature used by the
//! exchange-correlation integration. The construction follows the
//! standard recipe:
//!
//! 1. **Per atom** — take the tensor product of a
//!    [`radial`](crate::dft::radial) quadrature (Treutler-Ahlrichs) and
//!    a [`lebedev`](crate::dft::lebedev) angular quadrature. The radial
//!    weight already carries the `r²` volume element, the Lebedev
//!    weight carries the `4π`, so a point's *atomic* weight is
//!    `4π · w_radial · w_lebedev`.
//! 2. **Becke partition** — multiply every point's atomic weight by the
//!    [`becke`](crate::dft::becke) fuzzy-cell weight so the atomic
//!    grids tile space without double-counting.
//!
//! The result is a flat list of [`GridPoint`]s — a position and a
//! weight — such that `∫ f(r) dr ≈ Σ_k weight_k · f(r_k)` over all of
//! space.
//!
//! ## Coarseness
//!
//! [`GridQuality`] selects the radial-point count and the Lebedev grid
//! order. Three levels — `Coarse`, `Medium`, `Fine` — span quick
//! checks to a converged exchange-correlation integration; `Medium` is
//! the default.
//!
//! ## On-grid quantities
//!
//! [`GridDensity`] evaluates, on every grid point, the electron density
//! `ρ` and its gradient `∇ρ` from a density matrix and the basis. The
//! gradient is what a GGA functional (PBE, B88) needs; an LDA uses only
//! `ρ`. Basis-function values and gradients are cached per point so the
//! XC build and the `V_xc`-matrix build can reuse them.

use crate::basis::{BasisFunction, BasisSet};
use crate::dft::becke::becke_weight;
use crate::dft::lebedev::{lebedev_110, lebedev_26, lebedev_6};
use crate::dft::radial::{atomic_radial_scale, treutler_ahlrichs};
use crate::geometry::MolecularGeometry;
use nalgebra::DMatrix;

/// Coarseness of the molecular integration grid.
///
/// Trades accuracy of the XC integration against cost. The electron
/// count integrates to within a few parts in `10²` on `Coarse`, a few
/// in `10³` on `Medium` and a few in `10⁴` on `Fine` for the small
/// molecules this crate targets.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub enum GridQuality {
    /// Smallest grid — ~20 radial × 6 angular points per atom. Fast
    /// sanity checks.
    Coarse,
    /// Default grid — ~50 radial × 26 angular points per atom.
    #[default]
    Medium,
    /// Largest grid — ~85 radial × 110 angular points per atom
    /// (Lebedev order 17). Converged XC integration for small
    /// molecules.
    Fine,
}

impl GridQuality {
    /// Number of radial quadrature points per atom.
    fn radial_points(self) -> usize {
        match self {
            GridQuality::Coarse => 20,
            GridQuality::Medium => 50,
            GridQuality::Fine => 85,
        }
    }

    /// The Lebedev angular grid for this quality.
    fn angular(self) -> Vec<crate::dft::lebedev::LebedevPoint> {
        match self {
            GridQuality::Coarse => lebedev_6(),
            GridQuality::Medium => lebedev_26(),
            GridQuality::Fine => lebedev_110(),
        }
    }
}

/// One point of the molecular integration grid.
#[derive(Copy, Clone, Debug)]
pub struct GridPoint {
    /// Cartesian position in bohr.
    pub position: [f64; 3],
    /// Full quadrature weight — `∫ f dr ≈ Σ_k weight_k f(r_k)`.
    pub weight: f64,
}

/// The assembled atom-centred molecular grid.
#[derive(Clone, Debug)]
pub struct MolecularGrid {
    /// Every grid point with its Becke-partitioned weight.
    pub points: Vec<GridPoint>,
    /// The quality level the grid was built at.
    pub quality: GridQuality,
}

impl MolecularGrid {
    /// Build the molecular grid for `geometry` at the requested
    /// `quality`.
    ///
    /// Each atom contributes a Treutler-Ahlrichs × Lebedev product
    /// grid; the Becke scheme then partitions the overlapping atomic
    /// grids into a single space-filling quadrature.
    pub fn build(geometry: &MolecularGeometry, quality: GridQuality) -> Self {
        let positions: Vec<[f64; 3]> = geometry.atoms.iter().map(|a| a.position).collect();
        let atom_z: Vec<u8> = geometry
            .atoms
            .iter()
            .map(|a| a.element.atomic_number())
            .collect();
        let angular = quality.angular();
        let n_radial = quality.radial_points();

        let mut points = Vec::new();
        for (atom_idx, atom) in geometry.atoms.iter().enumerate() {
            let xi = atomic_radial_scale(atom.element.atomic_number());
            let radial = treutler_ahlrichs(n_radial, xi);
            let centre = atom.position;
            for rp in &radial {
                for ap in &angular {
                    let position = [
                        centre[0] + rp.radius * ap.dir[0],
                        centre[1] + rp.radius * ap.dir[1],
                        centre[2] + rp.radius * ap.dir[2],
                    ];
                    // Atomic weight: 4π folds in here (Lebedev weights
                    // sum to 1, radial weight carries r²).
                    let atomic_w = 4.0 * std::f64::consts::PI * rp.weight * ap.weight;
                    let becke = becke_weight(position, atom_idx, &positions, &atom_z);
                    let weight = atomic_w * becke;
                    // Drop negligible-weight points to keep the grid
                    // lean; they contribute nothing to any integral.
                    if weight > 1.0e-14 {
                        points.push(GridPoint { position, weight });
                    }
                }
            }
        }
        MolecularGrid { points, quality }
    }

    /// Number of grid points.
    #[inline]
    pub fn n_points(&self) -> usize {
        self.points.len()
    }

    /// Sum of every grid weight — should approximate the volume
    /// measure; mostly a diagnostic.
    pub fn total_weight(&self) -> f64 {
        self.points.iter().map(|p| p.weight).sum()
    }
}

/// Value and gradient of one contracted Cartesian Gaussian at a point.
///
/// Returns `(φ, ∂φ/∂x, ∂φ/∂y, ∂φ/∂z)`.
pub fn basis_value_grad(f: &BasisFunction, point: [f64; 3]) -> (f64, [f64; 3]) {
    let dx = point[0] - f.centre[0];
    let dy = point[1] - f.centre[1];
    let dz = point[2] - f.centre[2];
    let r2 = dx * dx + dy * dy + dz * dz;

    // Radial part R = Σ c_p e^{-α_p r²} and its r-derivative factor
    // dR = Σ c_p (−2α_p) e^{-α_p r²}.
    let mut radial = 0.0;
    let mut radial_d = 0.0;
    for p in &f.primitives {
        let e = (-p.exponent * r2).exp();
        radial += p.coefficient * e;
        radial_d += p.coefficient * (-2.0 * p.exponent) * e;
    }

    let (i, j, k) = (f.cart.0 as i32, f.cart.1 as i32, f.cart.2 as i32);
    let xi = dx.powi(i);
    let yj = dy.powi(j);
    let zk = dz.powi(k);
    let angular = xi * yj * zk;
    let value = angular * radial;

    // ∂φ/∂x = [ i x^{i-1} y^j z^k ] R + x^i y^j z^k · R' · x
    // (the chain rule on r² gives the extra dx factor on R').
    let dxi = if i > 0 {
        i as f64 * dx.powi(i - 1)
    } else {
        0.0
    };
    let dyj = if j > 0 {
        j as f64 * dy.powi(j - 1)
    } else {
        0.0
    };
    let dzk = if k > 0 {
        k as f64 * dz.powi(k - 1)
    } else {
        0.0
    };
    let grad = [
        dxi * yj * zk * radial + angular * radial_d * dx,
        xi * dyj * zk * radial + angular * radial_d * dy,
        xi * yj * dzk * radial + angular * radial_d * dz,
    ];
    (value, grad)
}

/// Per-grid-point electron density and density gradient.
///
/// Built once from a density matrix and the basis; the cached
/// basis-function values and gradients let the XC energy build and the
/// `V_xc`-matrix build share the work.
#[derive(Clone, Debug)]
pub struct GridDensity {
    /// Electron density `ρ` at each grid point.
    pub rho: Vec<f64>,
    /// Density gradient `∇ρ = (∂ρ/∂x, ∂ρ/∂y, ∂ρ/∂z)` at each point.
    pub grad: Vec<[f64; 3]>,
    /// Cached basis-function values: `phi[point][μ]`.
    pub phi: Vec<Vec<f64>>,
    /// Cached basis-function gradients: `dphi[point][μ] = ∇φ_μ`.
    pub dphi: Vec<Vec<[f64; 3]>>,
}

impl GridDensity {
    /// Evaluate the density and its gradient on every grid point.
    ///
    /// `density` is the closed-shell (or total) AO density matrix `D`;
    /// the density is `ρ = Σ_{μν} D_{μν} φ_μ φ_ν` and the gradient is
    /// `∇ρ = 2 Σ_{μν} D_{μν} φ_μ ∇φ_ν` (using `D` symmetric).
    pub fn evaluate(grid: &MolecularGrid, basis: &BasisSet, density: &DMatrix<f64>) -> GridDensity {
        let n = basis.n_functions();
        let np = grid.points.len();
        let mut rho = vec![0.0; np];
        let mut grad = vec![[0.0; 3]; np];
        let mut phi_all = Vec::with_capacity(np);
        let mut dphi_all = Vec::with_capacity(np);

        for (pi, gp) in grid.points.iter().enumerate() {
            let mut phi = vec![0.0; n];
            let mut dphi = vec![[0.0; 3]; n];
            for (mu, bf) in basis.functions.iter().enumerate() {
                let (v, g) = basis_value_grad(bf, gp.position);
                phi[mu] = v;
                dphi[mu] = g;
            }
            // ρ = Σ_{μν} D_{μν} φ_μ φ_ν ; first form t_ν = Σ_μ D_{μν} φ_μ.
            let mut t = vec![0.0; n];
            for nu in 0..n {
                let mut acc = 0.0;
                for mu in 0..n {
                    acc += density[(mu, nu)] * phi[mu];
                }
                t[nu] = acc;
            }
            let mut r = 0.0;
            let mut gx = 0.0;
            let mut gy = 0.0;
            let mut gz = 0.0;
            for nu in 0..n {
                r += t[nu] * phi[nu];
                // ∇ρ = 2 Σ t_ν ∇φ_ν.
                gx += 2.0 * t[nu] * dphi[nu][0];
                gy += 2.0 * t[nu] * dphi[nu][1];
                gz += 2.0 * t[nu] * dphi[nu][2];
            }
            rho[pi] = r;
            grad[pi] = [gx, gy, gz];
            phi_all.push(phi);
            dphi_all.push(dphi);
        }
        GridDensity {
            rho,
            grad,
            phi: phi_all,
            dphi: dphi_all,
        }
    }

    /// Integrate the electron density over the grid — `∫ ρ dr`. For a
    /// converged density this equals the electron count.
    pub fn integrate_electrons(&self, grid: &MolecularGrid) -> f64 {
        self.rho
            .iter()
            .zip(&grid.points)
            .map(|(&r, gp)| r * gp.weight)
            .sum()
    }

    /// The reduced density gradient norm `|∇ρ|` at grid point `i`.
    pub fn grad_norm(&self, i: usize) -> f64 {
        let g = self.grad[i];
        (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Atom;
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
    fn grid_grows_with_quality() {
        let (geom, _) = h2();
        let coarse = MolecularGrid::build(&geom, GridQuality::Coarse);
        let medium = MolecularGrid::build(&geom, GridQuality::Medium);
        let fine = MolecularGrid::build(&geom, GridQuality::Fine);
        assert!(coarse.n_points() < medium.n_points());
        assert!(medium.n_points() < fine.n_points());
    }

    #[test]
    fn grid_weights_are_positive() {
        let (geom, _) = h2();
        let grid = MolecularGrid::build(&geom, GridQuality::Medium);
        for p in &grid.points {
            assert!(p.weight > 0.0, "weight {}", p.weight);
        }
    }

    /// The grid must integrate a normalised Gaussian to 1. Place a
    /// unit-normalised `e^{-r²}` Gaussian on the first atom and sum.
    #[test]
    fn grid_integrates_a_gaussian_to_unity() {
        let (geom, _) = h2();
        let grid = MolecularGrid::build(&geom, GridQuality::Fine);
        // ∫ (α/π)^{3/2} e^{-α r²} dr = 1.
        let alpha = 1.3;
        let norm = (alpha / std::f64::consts::PI).powf(1.5);
        let centre = geom.atoms[0].position;
        let integral: f64 = grid
            .points
            .iter()
            .map(|p| {
                let dx = p.position[0] - centre[0];
                let dy = p.position[1] - centre[1];
                let dz = p.position[2] - centre[2];
                let r2 = dx * dx + dy * dy + dz * dz;
                p.weight * norm * (-alpha * r2).exp()
            })
            .sum();
        assert!((integral - 1.0).abs() < 1.0e-3, "∫ Gaussian = {integral}");
    }

    /// The grid must integrate the SCF electron density to the electron
    /// count. For H₂/STO-3G that is 2 electrons.
    #[test]
    fn grid_integrates_density_to_electron_count() {
        let (geom, basis) = h2();
        let ints = IntegralSet::compute(&geom, &basis);
        let rhf = run_rhf_scf(&ints, 2, ScfSettings::default()).unwrap();
        let grid = MolecularGrid::build(&geom, GridQuality::Fine);
        let gd = GridDensity::evaluate(&grid, &basis, &rhf.density);
        let n = gd.integrate_electrons(&grid);
        assert!((n - 2.0).abs() < 1.0e-3, "∫ρ = {n} (expected 2)");
    }

    /// The analytic basis-function gradient must match a central
    /// finite difference.
    #[test]
    fn basis_gradient_matches_finite_difference() {
        let (_, basis) = h2();
        let f = &basis.functions[0];
        let p = [0.3, -0.2, 0.5];
        let (_, grad) = basis_value_grad(f, p);
        let h = 1.0e-5;
        for axis in 0..3 {
            let mut pp = p;
            let mut pm = p;
            pp[axis] += h;
            pm[axis] -= h;
            let (vp, _) = basis_value_grad(f, pp);
            let (vm, _) = basis_value_grad(f, pm);
            let fd = (vp - vm) / (2.0 * h);
            assert!(
                (grad[axis] - fd).abs() < 1.0e-6,
                "axis {axis}: analytic {} vs FD {}",
                grad[axis],
                fd
            );
        }
    }

    /// The density gradient must match a central finite difference of
    /// the density.
    #[test]
    fn density_gradient_matches_finite_difference() {
        let (geom, basis) = h2();
        let ints = IntegralSet::compute(&geom, &basis);
        let rhf = run_rhf_scf(&ints, 2, ScfSettings::default()).unwrap();
        let d = &rhf.density;
        let point = [0.1, 0.2, 0.7];
        // Analytic gradient via a one-point grid.
        let one = MolecularGrid {
            points: vec![GridPoint {
                position: point,
                weight: 1.0,
            }],
            quality: GridQuality::Coarse,
        };
        let gd = GridDensity::evaluate(&one, &basis, d);
        let h = 1.0e-5;
        for axis in 0..3 {
            let mut pp = point;
            let mut pm = point;
            pp[axis] += h;
            pm[axis] -= h;
            let rho_p = crate::properties::grid::density_value(&basis, d, pp);
            let rho_m = crate::properties::grid::density_value(&basis, d, pm);
            let fd = (rho_p - rho_m) / (2.0 * h);
            assert!(
                (gd.grad[0][axis] - fd).abs() < 1.0e-6,
                "axis {axis}: analytic {} vs FD {}",
                gd.grad[0][axis],
                fd
            );
        }
    }
}
