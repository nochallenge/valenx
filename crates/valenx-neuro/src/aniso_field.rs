//! Anisotropic / heterogeneous extracellular field by direct FEM.
//!
//! The v1 [`crate::field`] module reuses valenx-fem's *scalar* steady-conduction
//! solver, so conductivity is a single number. Real nervous tissue is
//! **anisotropic** — current flows several times more easily along a
//! white-matter tract than across it — and **heterogeneous** (grey vs white
//! matter differ). This module solves `−∇·(σ∇φ)=I` with a full per-element
//! 3×3 conductivity tensor `σ`, assembling the scalar stiffness on linear
//! tetrahedra (`Kₑ[i][j] = Vₑ · ∇Nᵢᵀ σ ∇Nⱼ`) and solving the symmetric
//! positive-definite system by conjugate gradient.
//!
//! Validated against the closed-form infinite-medium anisotropic point source
//! `φ = I / (4π·√(det σ)·√(rᵀ σ⁻¹ r))`, which reduces to `I/(4πσr)` when σ is
//! isotropic.
//!
//! Units: node coordinates metres, conductivity S/m, current amperes at the
//! FEM boundary; potentials returned in mV.

use nalgebra::{Matrix3, Vector3};
use std::collections::HashMap;

/// A symmetric positive-definite 3×3 conductivity tensor (S/m).
pub type Conductivity = Matrix3<f64>;

/// An anisotropic tissue block: a structured cubic tetrahedral mesh (Kuhn
/// 6-tet split, centred at the origin) carrying a per-element conductivity
/// tensor.
pub struct AnisoTissue {
    nodes: Vec<Vector3<f64>>,
    tets: Vec<[u32; 4]>,
    sigma: Vec<Conductivity>,
    n: usize,
    spacing_m: f64,
}

/// Constant (linear-tet) global shape-function gradients and the element
/// volume, from the four vertex positions.
fn tet_grads(p: &[Vector3<f64>; 4]) -> ([Vector3<f64>; 4], f64) {
    let j = Matrix3::from_columns(&[p[1] - p[0], p[2] - p[0], p[3] - p[0]]);
    let vol = j.determinant().abs() / 6.0;
    let jit = j.try_inverse().expect("non-degenerate tetrahedron").transpose();
    let gl = [
        Vector3::new(-1.0, -1.0, -1.0),
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
    ];
    ([jit * gl[0], jit * gl[1], jit * gl[2], jit * gl[3]], vol)
}

impl AnisoTissue {
    /// Build a `side_mm` cube with `n` nodes per edge (`n ≥ 3`), assigning each
    /// tetrahedron the conductivity tensor `sigma_at(centroid_m)`.
    pub fn cube(side_mm: f64, n: usize, sigma_at: impl Fn(Vector3<f64>) -> Conductivity) -> Self {
        assert!(n >= 3, "need at least 3 nodes per edge for an interior source");
        let side_m = side_mm * 1.0e-3;
        let spacing = side_m / (n as f64 - 1.0);
        let half = side_m / 2.0;

        let mut nodes = Vec::with_capacity(n * n * n);
        for k in 0..n {
            for j in 0..n {
                for i in 0..n {
                    nodes.push(Vector3::new(
                        i as f64 * spacing - half,
                        j as f64 * spacing - half,
                        k as f64 * spacing - half,
                    ));
                }
            }
        }
        let idx = |i: usize, j: usize, k: usize| (i + n * j + n * n * k) as u32;

        let mut tets = Vec::new();
        let mut sigma = Vec::new();
        for k in 0..n - 1 {
            for j in 0..n - 1 {
                for i in 0..n - 1 {
                    let c = |di: usize, dj: usize, dk: usize| idx(i + di, j + dj, k + dk);
                    let cell = [
                        [c(0, 0, 0), c(1, 0, 0), c(1, 1, 0), c(1, 1, 1)],
                        [c(0, 0, 0), c(1, 0, 0), c(1, 0, 1), c(1, 1, 1)],
                        [c(0, 0, 0), c(0, 1, 0), c(1, 1, 0), c(1, 1, 1)],
                        [c(0, 0, 0), c(0, 1, 0), c(0, 1, 1), c(1, 1, 1)],
                        [c(0, 0, 0), c(0, 0, 1), c(1, 0, 1), c(1, 1, 1)],
                        [c(0, 0, 0), c(0, 0, 1), c(0, 1, 1), c(1, 1, 1)],
                    ];
                    for t in cell {
                        let centroid = (nodes[t[0] as usize]
                            + nodes[t[1] as usize]
                            + nodes[t[2] as usize]
                            + nodes[t[3] as usize])
                            / 4.0;
                        sigma.push(sigma_at(centroid));
                        tets.push(t);
                    }
                }
            }
        }
        Self { nodes, tets, sigma, n, spacing_m: spacing }
    }

    /// A homogeneous anisotropic block — one conductivity tensor everywhere.
    pub fn homogeneous(side_mm: f64, n: usize, sigma: Conductivity) -> Self {
        Self::cube(side_mm, n, move |_| sigma)
    }

    fn node_index(&self, i: usize, j: usize, k: usize) -> usize {
        i + self.n * j + self.n * self.n * k
    }

    fn is_boundary(&self, node: usize) -> bool {
        let n = self.n;
        let x = node % n;
        let y = (node / n) % n;
        let z = node / (n * n);
        x == 0 || x == n - 1 || y == 0 || y == n - 1 || z == 0 || z == n - 1
    }

    /// Inject `current_ua` (µA) at the central node, ground the outer faces
    /// (φ = 0), and solve the extracellular potential field by FEM + CG.
    pub fn solve_point_source(&self, current_ua: f64) -> SolvedField {
        self.solve_with_boundary(current_ua, |_| 0.0)
    }

    /// As [`Self::solve_point_source`], but with the outer boundary held at
    /// `boundary_mv(pos)` (mV) rather than ground. Imposing the analytic far
    /// field on the boundary removes finite-domain truncation, so the interior
    /// solution can be checked against a closed form to mesh accuracy.
    pub fn solve_with_boundary(
        &self,
        current_ua: f64,
        boundary_mv: impl Fn(Vector3<f64>) -> f64,
    ) -> SolvedField {
        let center = self.node_index(self.n / 2, self.n / 2, self.n / 2);
        self.solve_currents(&[(center, current_ua)], boundary_mv)
    }

    /// Grid node index at `(di, dj, dk)` steps from the centre — for placing
    /// electrode contacts.
    pub fn node_at_steps(&self, di: i32, dj: i32, dk: i32) -> usize {
        let c = (self.n / 2) as i32;
        self.node_index((c + di) as usize, (c + dj) as usize, (c + dk) as usize)
    }

    /// Solve with arbitrary nodal current `loads` (interior node index, µA) and
    /// a Dirichlet `boundary_mv(pos)` (mV). Several contacts superpose because
    /// the operator is linear — the basis for current steering.
    pub fn solve_currents(
        &self,
        loads: &[(usize, f64)],
        boundary_mv: impl Fn(Vector3<f64>) -> f64,
    ) -> SolvedField {
        let n_nodes = self.nodes.len();
        // Free DOFs = interior nodes; boundary nodes are fixed (Dirichlet), so
        // their columns move to the right-hand side of K_ff φ_f = b_f.
        let mut reduced = vec![usize::MAX; n_nodes];
        let mut n_free = 0;
        for (node, slot) in reduced.iter_mut().enumerate() {
            if !self.is_boundary(node) {
                *slot = n_free;
                n_free += 1;
            }
        }
        let mut phi_fixed = vec![0.0; n_nodes]; // volts
        for (node, slot) in phi_fixed.iter_mut().enumerate() {
            if self.is_boundary(node) {
                *slot = boundary_mv(self.nodes[node]) * 1.0e-3;
            }
        }

        let mut rows: Vec<HashMap<usize, f64>> = vec![HashMap::new(); n_free];
        let mut b = vec![0.0; n_free];
        for (e, t) in self.tets.iter().enumerate() {
            let p = [
                self.nodes[t[0] as usize],
                self.nodes[t[1] as usize],
                self.nodes[t[2] as usize],
                self.nodes[t[3] as usize],
            ];
            let (g, vol) = tet_grads(&p);
            let s = self.sigma[e];
            for a in 0..4 {
                let ra = reduced[t[a] as usize];
                if ra == usize::MAX {
                    continue;
                }
                for b_idx in 0..4 {
                    let gb = t[b_idx] as usize;
                    let kab = vol * g[a].dot(&(s * g[b_idx]));
                    if reduced[gb] == usize::MAX {
                        b[ra] -= kab * phi_fixed[gb]; // known boundary value → RHS
                    } else {
                        *rows[ra].entry(reduced[gb]).or_insert(0.0) += kab;
                    }
                }
            }
        }

        for &(node, ua) in loads {
            if reduced[node] != usize::MAX {
                b[reduced[node]] += ua * 1.0e-6; // amperes
            }
        }

        let (indptr, indices, data) = to_csr(&rows);
        let phi = conjugate_gradient(&indptr, &indices, &data, &b);

        let mut phi_mv = vec![0.0; n_nodes];
        for node in 0..n_nodes {
            phi_mv[node] = if reduced[node] != usize::MAX {
                phi[reduced[node]] * 1.0e3
            } else {
                phi_fixed[node] * 1.0e3
            };
        }
        SolvedField { phi_mv, n: self.n, spacing_m: self.spacing_m }
    }
}

/// Flatten per-row sparse maps into compressed-sparse-row arrays.
fn to_csr(rows: &[HashMap<usize, f64>]) -> (Vec<usize>, Vec<usize>, Vec<f64>) {
    let mut indptr = vec![0usize; rows.len() + 1];
    for (i, r) in rows.iter().enumerate() {
        indptr[i + 1] = indptr[i] + r.len();
    }
    let nnz = indptr[rows.len()];
    let mut indices = vec![0usize; nnz];
    let mut data = vec![0.0; nnz];
    for (i, r) in rows.iter().enumerate() {
        let mut off = indptr[i];
        for (&col, &val) in r {
            indices[off] = col;
            data[off] = val;
            off += 1;
        }
    }
    (indptr, indices, data)
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Solve the SPD system `A x = b` (A in CSR) by conjugate gradient.
fn conjugate_gradient(indptr: &[usize], indices: &[usize], data: &[f64], b: &[f64]) -> Vec<f64> {
    let n = b.len();
    let matvec = |x: &[f64], y: &mut [f64]| {
        for i in 0..n {
            let mut s = 0.0;
            for k in indptr[i]..indptr[i + 1] {
                s += data[k] * x[indices[k]];
            }
            y[i] = s;
        }
    };
    let mut x = vec![0.0; n];
    let mut r = b.to_vec();
    let mut p = r.clone();
    let mut ap = vec![0.0; n];
    let mut rsold = dot(&r, &r);
    let tol2 = 1.0e-20 * rsold.max(1.0e-30); // relative ‖r‖² target
    for _ in 0..(10 * n + 50) {
        matvec(&p, &mut ap);
        let pap = dot(&p, &ap);
        if pap <= 0.0 {
            break;
        }
        let alpha = rsold / pap;
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        let rsnew = dot(&r, &r);
        if rsnew < tol2 {
            break;
        }
        let beta = rsnew / rsold;
        for i in 0..n {
            p[i] = r[i] + beta * p[i];
        }
        rsold = rsnew;
    }
    x
}

/// A solved potential field over an [`AnisoTissue`].
pub struct SolvedField {
    phi_mv: Vec<f64>,
    n: usize,
    spacing_m: f64,
}

impl SolvedField {
    fn node_index(&self, i: usize, j: usize, k: usize) -> usize {
        i + self.n * j + self.n * self.n * k
    }

    /// Potential (mV) at the grid node `(di, dj, dk)` steps from the centre.
    pub fn phi_at_steps(&self, di: i32, dj: i32, dk: i32) -> f64 {
        let c = (self.n / 2) as i32;
        self.phi_mv[self.node_index((c + di) as usize, (c + dj) as usize, (c + dk) as usize)]
    }

    /// Grid node spacing (metres).
    pub fn spacing_m(&self) -> f64 {
        self.spacing_m
    }

    /// Per-node potential (mV), indexed `i + n·j + n²·k`.
    pub fn phi_mv(&self) -> &[f64] {
        &self.phi_mv
    }

    /// The x-coordinate (metres) of the positive-potential-weighted field
    /// centroid — a simple scalar measure of where stimulation is focused,
    /// used to quantify current steering.
    pub fn focus_x_m(&self) -> f64 {
        let c = (self.n as f64 - 1.0) / 2.0;
        let mut num = 0.0;
        let mut den = 0.0;
        for (idx, &phi) in self.phi_mv.iter().enumerate() {
            if phi > 0.0 {
                let i = (idx % self.n) as f64;
                num += phi * (i - c) * self.spacing_m;
                den += phi;
            }
        }
        if den > 0.0 {
            num / den
        } else {
            0.0
        }
    }
}

/// Analytic potential (mV) of a point source `current_ua` (µA) in an infinite
/// homogeneous **anisotropic** medium of conductivity tensor `sigma` (S/m), at
/// displacement `r_m` (metres): `φ = I / (4π·√(det σ)·√(rᵀ σ⁻¹ r))`. Reduces
/// to `I/(4πσr)` for isotropic σ.
pub fn analytic_aniso_point_source_mv(
    current_ua: f64,
    sigma: &Conductivity,
    r_m: Vector3<f64>,
) -> f64 {
    let i = current_ua * 1.0e-6;
    let det = sigma.determinant();
    let sinv = sigma.try_inverse().expect("conductivity tensor invertible");
    let rsr = r_m.dot(&(sinv * r_m));
    (i / (4.0 * std::f64::consts::PI * det.sqrt() * rsr.sqrt())) * 1.0e3
}

#[cfg(test)]
mod tests {
    use super::*;

    fn diag(sx: f64, sy: f64, sz: f64) -> Conductivity {
        Matrix3::from_diagonal(&Vector3::new(sx, sy, sz))
    }

    #[test]
    fn isotropic_recovers_coulomb_point_source() {
        // Isotropic σ: the anisotropic closed form reduces to I/(4πσr). With the
        // analytic far field on the boundary, the interior FEM matches it.
        let s = diag(0.3, 0.3, 0.3);
        let f = AnisoTissue::homogeneous(40.0, 21, s)
            .solve_with_boundary(100.0, |pos| analytic_aniso_point_source_mv(100.0, &s, pos));
        let dx = f.spacing_m();
        for m in [3, 4, 5] {
            let fem = f.phi_at_steps(m, 0, 0);
            let ana =
                analytic_aniso_point_source_mv(100.0, &s, Vector3::new(m as f64 * dx, 0.0, 0.0));
            assert!(
                (fem / ana - 1.0).abs() < 0.12,
                "isotropic FEM ≈ Coulomb at m={m}: {fem:.2} vs {ana:.2}"
            );
        }
    }

    #[test]
    fn anisotropic_point_source_matches_closed_form() {
        // White-matter-like anisotropy (longitudinal 3× transverse).
        let s = diag(0.6, 0.2, 0.2);
        let f = AnisoTissue::homogeneous(40.0, 21, s)
            .solve_with_boundary(100.0, |pos| analytic_aniso_point_source_mv(100.0, &s, pos));
        let dx = f.spacing_m();
        for m in [3, 4, 5] {
            let r = m as f64 * dx;
            let fem_x = f.phi_at_steps(m, 0, 0);
            let ana_x = analytic_aniso_point_source_mv(100.0, &s, Vector3::new(r, 0.0, 0.0));
            let fem_y = f.phi_at_steps(0, m, 0);
            let ana_y = analytic_aniso_point_source_mv(100.0, &s, Vector3::new(0.0, r, 0.0));
            assert!(
                (fem_x / ana_x - 1.0).abs() < 0.15,
                "FEM ≈ analytic along x at m={m}: {fem_x:.2} vs {ana_x:.2}"
            );
            assert!(
                (fem_y / ana_y - 1.0).abs() < 0.15,
                "FEM ≈ analytic along y at m={m}: {fem_y:.2} vs {ana_y:.2}"
            );
        }
    }

    #[test]
    fn anisotropy_extends_field_along_high_conductivity_axis() {
        // Grounded box, σx = 3·σy: at equal distance the potential is larger
        // along the high-conductivity x-axis than across it (the y-axis).
        let s = diag(0.6, 0.2, 0.2);
        let f = AnisoTissue::homogeneous(40.0, 21, s).solve_point_source(100.0);
        let px = f.phi_at_steps(4, 0, 0);
        let py = f.phi_at_steps(0, 4, 0);
        assert!(px > py * 1.4, "field extends along high-σ axis: φx={px:.2} φy={py:.2}");
    }

    #[test]
    fn heterogeneous_field_solves_and_is_finite() {
        // Two slabs: lower-σ grey matter (z<0), higher-σ white matter (z>0).
        let f = AnisoTissue::cube(40.0, 15, |pos| {
            if pos.z < 0.0 {
                diag(0.2, 0.2, 0.2)
            } else {
                diag(0.6, 0.6, 0.6)
            }
        })
        .solve_point_source(100.0);
        assert!(f.phi_mv().iter().all(|v| v.is_finite()), "all potentials finite");
        let center = f.phi_at_steps(0, 0, 0);
        assert!(
            center > f.phi_at_steps(0, 0, 3) && center > f.phi_at_steps(0, 0, -3),
            "source node is the peak; field falls off into both slabs"
        );
    }
}
