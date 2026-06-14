//! Native in-process steady-state thermal FEA solver (Phase 24.7).
//!
//! ## What this is
//!
//! A genuine, self-contained **steady-state heat-conduction finite-
//! element solver** for 4-node linear tetrahedra (`C3D4` /
//! [`ElementType::Tet4`]). It solves the steady heat equation
//!
//! ```text
//!   −∇·(k ∇T) = q     in the body
//!            T = T̄    on the Dirichlet boundary  (fixed temperature)
//!     −k ∂T/∂n = q̄    on the Neumann boundary    (prescribed flux)
//! ```
//!
//! discretised into the linear system `Kₜ · T = Q`, where `Kₜ` is the
//! global **conductivity matrix** (the thermal analogue of the elastic
//! stiffness matrix), `T` is the unknown nodal temperature field, and
//! `Q` is the heat-load vector assembled from internal heat generation
//! and boundary heat fluxes.
//!
//! This is the in-process alternative to dispatching a `*HEAT TRANSFER,
//! STEADY STATE` step to the CalculiX / Elmer subprocess adapters.
//!
//! ## Method
//!
//! 1. **Element conductivity.** Heat conduction in a linear tet is a
//!    *scalar* field problem — one temperature DOF per node, not three
//!    displacement DOFs. The four shape functions are linear, so their
//!    gradients `∇Nₐ` are constant over the element (exactly as in the
//!    elastic solver). The element conductivity matrix is the closed
//!    form
//!
//!    ```text
//!      Kₜₑ = k · Vₑ · Gᵀ · G
//!    ```
//!
//!    where `G` is the 3×4 gradient matrix whose columns are the four
//!    `∇Nₐ`, `k` is the (isotropic) thermal conductivity, and `Vₑ` the
//!    tet volume. The integrand `∇Nₐ·∇Nᵦ` is constant, so a single
//!    point integrates it exactly.
//! 2. **Assembly.** Each 4×4 `Kₜₑ` is scatter-added into a global
//!    sparse matrix of size `n_nodes`.
//! 3. **Boundary conditions.**
//!    - **Neumann (heat flux / heat source).** A prescribed nodal heat
//!      input goes straight into the load vector `Q` — these are the
//!      natural BCs of the variational form, so no matrix surgery is
//!      needed. Distributed surface fluxes are pre-integrated by the
//!      caller into equivalent nodal heat inputs.
//!    - **Dirichlet (fixed temperature).** Imposed by the large-penalty
//!      method: a stiff term `k·penalty` is added to the constrained
//!      node's diagonal and `k·penalty·T̄` to its load entry, so the
//!      row reads ≈ `penalty·T = penalty·T̄`. This keeps the matrix
//!      symmetric positive-definite, so the Cholesky path stays valid.
//! 4. **Solve.** The assembled system is converted to CSC and
//!    factorised with [`nalgebra_sparse::factorization::CscCholesky`];
//!    the temperature field is the back-substitution.
//! 5. **Flux recovery.** The element heat-flux vector
//!    `qₑ = −k·G·Tₑ` is constant over a linear tet; element fluxes are
//!    averaged into their nodes for a smooth nodal flux field.
//!
//! ## Honest scope — what it is and is not
//!
//! **It is a real thermal solver.** For a linear, isotropic, steady
//! conduction problem meshed with linear tets, the temperature field
//! it returns is the correct FE solution for that mesh. A bar with one
//! end held hot and the other held cold reproduces the analytic linear
//! temperature gradient to solver precision.
//!
//! **It is deliberately the 90 % case:**
//!
//! - **Steady state only.** No transient time-stepping (that needs the
//!   thermal *capacitance* matrix and a time integrator), no
//!   phase change, no radiation (radiation is non-linear in `T⁴`).
//! - **Linear conduction.** Conductivity is a single constant `k`; no
//!   temperature-dependent `k(T)`, no anisotropy.
//! - **Convection (Robin) BCs are supported** via
//!   [`solve_steady_thermal_with_convection`]: a film boundary
//!   `q = h·A·(T − T∞)` folds the conductance `h·A` onto the node
//!   diagonal and `h·A·T∞` into the load. Radiation (non-linear in
//!   `T⁴`) is still out of scope.
//! - **Tet4 only**; needs a volume mesh (see
//!   [`crate::native_solver::structured_box_mesh`]).
//!
//! Within that scope the numbers are real.

use nalgebra::{DMatrix, DVector, Matrix3, Matrix3x4, Vector3};
use nalgebra_sparse::factorization::CscCholesky;
use nalgebra_sparse::{CooMatrix, CscMatrix};
use thiserror::Error;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

/// Errors from the native thermal solver.
#[derive(Debug, Error)]
pub enum ThermalSolverError {
    /// The supplied mesh has no [`ElementType::Tet4`] block.
    #[error("mesh has no Tet4 element block — the thermal solver needs a tetrahedral volume mesh")]
    NoTetBlock,
    /// The mesh has zero nodes.
    #[error("mesh has no nodes")]
    EmptyMesh,
    /// A tetrahedron is degenerate (zero or negative signed volume).
    #[error("tetrahedron {0} is degenerate (|volume| < tolerance)")]
    DegenerateElement(usize),
    /// A connectivity index points past the end of the node array.
    #[error("element {elem} references node {node} but the mesh has only {n_nodes} nodes")]
    BadConnectivity {
        /// 0-based element index (or `usize::MAX` for a BC reference).
        elem: usize,
        /// Out-of-range node index.
        node: usize,
        /// Node-array length.
        n_nodes: usize,
    },
    /// The thermal conductivity is non-physical (`k ≤ 0`).
    #[error("invalid conductivity: must be finite and positive, got {0}")]
    BadConductivity(f64),
    /// A load or boundary-condition input carried a non-finite value
    /// (`NaN` or `±∞`). Round-1 H4: the conductivity and node coordinates
    /// were finiteness-checked, but a [`HeatLoad`] power and a
    /// [`FixedTemperature`] (prescribed Dirichlet BC) were pushed straight
    /// into the right-hand side `q` (the latter as `penalty·temperature`).
    /// The Cholesky back-substitution is pure arithmetic, so a non-finite
    /// RHS yields a non-finite temperature returned as `Ok(..)` — a single
    /// `power: NaN` silently corrupts the whole field. Validating the
    /// input up front lets the error name the cause.
    #[error("non-finite {kind} at node {node} (NaN or infinity is not a valid load/BC)")]
    InvalidLoad {
        /// 0-based node index carrying the bad value.
        node: usize,
        /// Which input was non-finite — e.g. `"heat power"`, `"prescribed
        /// temperature"`, or a convection coefficient / ambient temperature.
        kind: &'static str,
    },
    /// No Dirichlet (fixed-temperature) boundary condition was supplied.
    /// With only Neumann BCs the conductivity matrix is singular (the
    /// temperature is determined only up to an additive constant) and
    /// the solve has no unique answer.
    #[error("no fixed-temperature boundary — the steady problem is singular (temperature floats)")]
    NoDirichletBc,
    /// The Cholesky factorisation failed — the assembled system was not
    /// positive-definite.
    #[error("linear solve failed: conductivity matrix is not positive-definite")]
    SolveFailed,
}

/// A single-node Dirichlet (fixed-temperature) boundary condition.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct FixedTemperature {
    /// 0-based node index.
    pub node: usize,
    /// Prescribed temperature at the node (kelvin, or any consistent
    /// unit — the solver is unit-agnostic as long as `k`, the loads,
    /// and the geometry agree).
    pub temperature: f64,
}

/// A single-node Neumann boundary condition: a concentrated heat input
/// at a node (watts). Positive = heat flowing *into* the body.
///
/// A distributed surface flux `q̄` (W/m²) is converted to equivalent
/// nodal inputs by the caller — for a linear tet face of area `A` a
/// uniform flux contributes `q̄·A/3` to each of the face's three nodes.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct HeatLoad {
    /// 0-based node index.
    pub node: usize,
    /// Heat input at the node in watts (positive = into the body).
    pub power: f64,
}

/// A single-node **convection (Robin)** boundary condition modelling a film
/// boundary `q = h·A·(T − T∞)`: heat leaves the node in proportion to how far
/// it sits above the ambient temperature `T∞`.
///
/// `h_area` is the lumped conductance `h·A` (W/K) for the node's share of the
/// convecting surface — a distributed film coefficient `h` (W/(m²·K)) over a
/// linear-tet face of area `A` contributes `h·A/3` to each of the face's three
/// nodes, exactly like [`HeatLoad`]. The term folds `h·A` onto the node's
/// matrix diagonal and `h·A·T∞` into the load.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct ConvectionBc {
    /// 0-based node index.
    pub node: usize,
    /// Lumped convective conductance `h·A` at the node (W/K, non-negative).
    pub h_area: f64,
    /// Ambient (free-stream) temperature `T∞` the node convects to.
    pub ambient: f64,
}

/// Result of a steady-state thermal solve.
#[derive(Clone, Debug)]
pub struct ThermalSolution {
    /// Per-node temperature field, indexed by node.
    pub temperature: Vec<f64>,
    /// Per-node heat-flux vector `[qx, qy, qz]` in W/m², computed per
    /// element (constant over a linear tet) then averaged to the nodes.
    pub flux: Vec<[f64; 3]>,
}

impl ThermalSolution {
    /// Highest nodal temperature in the field.
    pub fn max_temperature(&self) -> f64 {
        self.temperature.iter().copied().fold(f64::MIN, f64::max)
    }

    /// Lowest nodal temperature in the field.
    pub fn min_temperature(&self) -> f64 {
        self.temperature.iter().copied().fold(f64::MAX, f64::min)
    }

    /// Largest nodal heat-flux magnitude.
    pub fn max_flux_magnitude(&self) -> f64 {
        self.flux
            .iter()
            .map(|q| (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]).sqrt())
            .fold(0.0, f64::max)
    }
}

/// Solve a steady-state heat-conduction problem on a tetrahedral mesh.
///
/// `mesh` must carry at least one [`ElementType::Tet4`] block.
/// `conductivity` is the isotropic thermal conductivity `k` in
/// W/(m·K). `fixed` supplies the Dirichlet (fixed-temperature) BCs (at
/// least one is required, or the steady problem is singular).
/// `loads` supplies the Neumann (heat-input) BCs.
///
/// Returns the per-node temperature field and the recovered nodal heat
/// flux.
///
/// # Errors
///
/// See [`ThermalSolverError`].
pub fn solve_steady_thermal(
    mesh: &Mesh,
    conductivity: f64,
    fixed: &[FixedTemperature],
    loads: &[HeatLoad],
) -> Result<ThermalSolution, ThermalSolverError> {
    solve_steady_thermal_with_convection(mesh, conductivity, fixed, loads, &[])
}

/// Solve a steady-state heat-conduction problem with **convection (Robin)**
/// boundary conditions, in addition to the Dirichlet/Neumann BCs that
/// [`solve_steady_thermal`] handles.
///
/// Each [`ConvectionBc`] models a film boundary `q = h·A·(T − T∞)` by folding
/// the conductance `h·A` onto the node's matrix diagonal and `h·A·T∞` into the
/// load vector — the natural way a Robin term enters the FE system. A
/// convection BC alone makes the conductivity matrix non-singular, so (unlike
/// pure conduction) **no Dirichlet BC is required** when at least one
/// convection BC is present.
///
/// # Errors
///
/// See [`ThermalSolverError`].
pub fn solve_steady_thermal_with_convection(
    mesh: &Mesh,
    conductivity: f64,
    fixed: &[FixedTemperature],
    loads: &[HeatLoad],
    convection: &[ConvectionBc],
) -> Result<ThermalSolution, ThermalSolverError> {
    let n_nodes = mesh.nodes.len();
    if n_nodes == 0 {
        return Err(ThermalSolverError::EmptyMesh);
    }
    if !conductivity.is_finite() || conductivity <= 0.0 {
        return Err(ThermalSolverError::BadConductivity(conductivity));
    }
    if fixed.is_empty() && convection.is_empty() {
        return Err(ThermalSolverError::NoDirichletBc);
    }
    let tets = tet_block(mesh).ok_or(ThermalSolverError::NoTetBlock)?;

    // --- assemble the global conductivity matrix (n_nodes square) ---
    let mut coo = CooMatrix::<f64>::new(n_nodes, n_nodes);
    let conn = &tets.connectivity;
    let n_elem = conn.len() / 4;
    for e in 0..n_elem {
        let nodes = [
            conn[4 * e] as usize,
            conn[4 * e + 1] as usize,
            conn[4 * e + 2] as usize,
            conn[4 * e + 3] as usize,
        ];
        for &nd in &nodes {
            if nd >= n_nodes {
                return Err(ThermalSolverError::BadConnectivity {
                    elem: e,
                    node: nd,
                    n_nodes,
                });
            }
        }
        let coords = [
            mesh.nodes[nodes[0]],
            mesh.nodes[nodes[1]],
            mesh.nodes[nodes[2]],
            mesh.nodes[nodes[3]],
        ];
        let kte = element_conductivity(&coords, conductivity)
            .ok_or(ThermalSolverError::DegenerateElement(e))?;
        // Scatter-add the 4×4 element matrix.
        for a in 0..4 {
            for b in 0..4 {
                let v = kte[(a, b)];
                if v != 0.0 {
                    coo.push(nodes[a], nodes[b], v);
                }
            }
        }
    }
    let mut csc = CscMatrix::from(&coo);

    // Find the peak diagonal for penalty scaling.
    let mut max_diag = 0.0_f64;
    for (r, c, v) in csc.triplet_iter() {
        if r == c {
            max_diag = max_diag.max(v.abs());
        }
    }
    if max_diag <= 0.0 {
        return Err(ThermalSolverError::SolveFailed);
    }
    let penalty = max_diag * 1.0e8;

    // --- heat-load vector ---
    let mut q = DVector::<f64>::zeros(n_nodes);
    for load in loads {
        if load.node >= n_nodes {
            return Err(ThermalSolverError::BadConnectivity {
                elem: usize::MAX,
                node: load.node,
                n_nodes,
            });
        }
        // Round-1 H4: a non-finite heat power flows straight into the RHS
        // q; the Cholesky solve would then return Ok(..) with a silently
        // non-finite temperature field. Reject it before it reaches q.
        if !load.power.is_finite() {
            return Err(ThermalSolverError::InvalidLoad {
                node: load.node,
                kind: "heat power",
            });
        }
        q[load.node] += load.power;
    }

    // --- Dirichlet BCs via the penalty method ---
    let mut penalty_diag = vec![0.0_f64; n_nodes];
    for fb in fixed {
        if fb.node >= n_nodes {
            return Err(ThermalSolverError::BadConnectivity {
                elem: usize::MAX,
                node: fb.node,
                n_nodes,
            });
        }
        // Round-1 H4: a prescribed temperature is folded into the RHS as
        // `penalty·temperature`; a non-finite value corrupts the solve
        // just like a non-finite power. Reject it before it reaches q.
        if !fb.temperature.is_finite() {
            return Err(ThermalSolverError::InvalidLoad {
                node: fb.node,
                kind: "prescribed temperature",
            });
        }
        penalty_diag[fb.node] += penalty;
        q[fb.node] += penalty * fb.temperature;
    }

    // --- convection (Robin) BCs: q = h·A·(T − T∞) folds h·A onto the
    // diagonal and h·A·T∞ into the load (reusing the penalty-diagonal vector
    // as a general "extra diagonal"). ---
    for cb in convection {
        if cb.node >= n_nodes {
            return Err(ThermalSolverError::BadConnectivity {
                elem: usize::MAX,
                node: cb.node,
                n_nodes,
            });
        }
        if !cb.h_area.is_finite() || cb.h_area < 0.0 {
            return Err(ThermalSolverError::InvalidLoad {
                node: cb.node,
                kind: "convection conductance",
            });
        }
        if !cb.ambient.is_finite() {
            return Err(ThermalSolverError::InvalidLoad {
                node: cb.node,
                kind: "convection ambient temperature",
            });
        }
        penalty_diag[cb.node] += cb.h_area;
        q[cb.node] += cb.h_area * cb.ambient;
    }

    // Fold the (Dirichlet penalty + convection) diagonal into the matrix.
    csc = add_diagonal(&csc, &penalty_diag);

    // --- solve ---
    let chol = CscCholesky::factor(&csc).map_err(|_| ThermalSolverError::SolveFailed)?;
    let t: DMatrix<f64> = chol.solve(&q);
    let t = t.column(0);

    // --- gather the nodal temperature field ---
    let temperature: Vec<f64> = (0..n_nodes).map(|n| t[n]).collect();

    // --- flux recovery: qₑ = −k·G·Tₑ per element, averaged to nodes ---
    let mut nodal_flux = vec![[0.0_f64; 3]; n_nodes];
    let mut nodal_count = vec![0_u32; n_nodes];
    for e in 0..n_elem {
        let nodes = [
            conn[4 * e] as usize,
            conn[4 * e + 1] as usize,
            conn[4 * e + 2] as usize,
            conn[4 * e + 3] as usize,
        ];
        let coords = [
            mesh.nodes[nodes[0]],
            mesh.nodes[nodes[1]],
            mesh.nodes[nodes[2]],
            mesh.nodes[nodes[3]],
        ];
        let Some(g) = gradient_matrix(&coords) else {
            continue;
        };
        // Element temperature vector (one DOF per node).
        let t_elem = nalgebra::Vector4::new(t[nodes[0]], t[nodes[1]], t[nodes[2]], t[nodes[3]]);
        // ∇T = G·Tₑ  (3×1);  q = −k·∇T.
        let grad_t = g * t_elem;
        let flux = grad_t * (-conductivity);
        for &nd in &nodes {
            nodal_flux[nd][0] += flux.x;
            nodal_flux[nd][1] += flux.y;
            nodal_flux[nd][2] += flux.z;
            nodal_count[nd] += 1;
        }
    }
    let mut flux = vec![[0.0_f64; 3]; n_nodes];
    for n in 0..n_nodes {
        let count = nodal_count[n].max(1) as f64;
        flux[n][0] = nodal_flux[n][0] / count;
        flux[n][1] = nodal_flux[n][1] / count;
        flux[n][2] = nodal_flux[n][2] / count;
    }

    Ok(ThermalSolution { temperature, flux })
}

/// First [`ElementType::Tet4`] block in the mesh, if any.
pub(crate) fn tet_block(mesh: &Mesh) -> Option<&ElementBlock> {
    mesh.element_blocks
        .iter()
        .find(|b| b.element_type == ElementType::Tet4)
        .filter(|b| b.connectivity.len() >= 4)
}

/// 3×4 gradient matrix `G` of a linear tet — column `a` is the
/// constant gradient `∇Nₐ` of shape function `a`. `None` for a
/// degenerate tet.
pub(crate) fn gradient_matrix(coords: &[Vector3<f64>; 4]) -> Option<Matrix3x4<f64>> {
    let p0 = coords[0];
    let e1 = coords[1] - p0;
    let e2 = coords[2] - p0;
    let e3 = coords[3] - p0;
    let jac = Matrix3::from_columns(&[e1, e2, e3]);
    let det = jac.determinant();
    if det.abs() < 1.0e-18 {
        return None;
    }
    let inv = jac.try_inverse()?;
    let invt = inv.transpose();
    let g1 = invt * Vector3::new(1.0, 0.0, 0.0);
    let g2 = invt * Vector3::new(0.0, 1.0, 0.0);
    let g3 = invt * Vector3::new(0.0, 0.0, 1.0);
    let g0 = -(g1 + g2 + g3);
    Some(Matrix3x4::from_columns(&[g0, g1, g2, g3]))
}

/// Signed volume of a tet.
pub(crate) fn tet_volume(coords: &[Vector3<f64>; 4]) -> f64 {
    let e1 = coords[1] - coords[0];
    let e2 = coords[2] - coords[0];
    let e3 = coords[3] - coords[0];
    e1.cross(&e2).dot(&e3) / 6.0
}

/// 4×4 element conductivity matrix `Kₜₑ = k·Vₑ·Gᵀ·G` for a linear tet.
/// `None` for a degenerate tet.
pub(crate) fn element_conductivity(coords: &[Vector3<f64>; 4], k: f64) -> Option<DMatrix<f64>> {
    let g = gradient_matrix(coords)?;
    let vol = tet_volume(coords).abs();
    if vol < 1.0e-18 {
        return None;
    }
    // Gᵀ·G is 4×4; scale by k·V.
    let gt_g = g.transpose() * g;
    let scaled = gt_g * (k * vol);
    // Convert the fixed-size 4×4 into a DMatrix for a uniform API with
    // the assembly loop.
    Some(DMatrix::from_iterator(4, 4, scaled.iter().copied()))
}

/// Return a copy of `csc` with `diag[i]` added to entry `(i, i)`.
/// Diagonal entries already present are incremented; missing ones are
/// inserted. Used to fold the Dirichlet penalty into the conductivity
/// matrix.
pub(crate) fn add_diagonal(csc: &CscMatrix<f64>, diag: &[f64]) -> CscMatrix<f64> {
    let n = csc.nrows();
    let mut coo = CooMatrix::<f64>::new(n, n);
    for (r, c, v) in csc.triplet_iter() {
        coo.push(r, c, *v);
    }
    for (i, &d) in diag.iter().enumerate() {
        if d != 0.0 {
            coo.push(i, i, d);
        }
    }
    CscMatrix::from(&coo)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::native_solver::structured_box_mesh;

    /// Node id of grid point `(i, j, k)` in a `structured_box_mesh`.
    fn nid(i: usize, j: usize, k: usize, nx: usize, ny: usize) -> usize {
        i + (nx + 1) * j + (nx + 1) * (ny + 1) * k
    }

    #[test]
    fn element_conductivity_is_symmetric_with_zero_row_sums() {
        // Gᵀ·G has zero row sums because the four shape-function
        // gradients of a linear tet sum to zero (Σ Nₐ ≡ 1 ⇒ Σ ∇Nₐ ≡ 0).
        // So each row of the element conductivity matrix sums to zero —
        // a constant temperature field carries no heat flux.
        let coords = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let kte = element_conductivity(&coords, 50.0).unwrap();
        for i in 0..4 {
            let row_sum: f64 = (0..4).map(|j| kte[(i, j)]).sum();
            assert!(row_sum.abs() < 1e-9, "row {i} sum {row_sum} ≠ 0");
            for j in 0..4 {
                assert!(
                    (kte[(i, j)] - kte[(j, i)]).abs() < 1e-9,
                    "Kₜₑ not symmetric"
                );
            }
        }
    }

    #[test]
    fn degenerate_tet_has_no_conductivity() {
        // Four coplanar points → zero volume → None.
        let coords = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
        ];
        assert!(element_conductivity(&coords, 1.0).is_none());
    }

    #[test]
    fn solve_rejects_bad_conductivity() {
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box params");
        let err = solve_steady_thermal(
            &mesh,
            -1.0,
            &[FixedTemperature {
                node: 0,
                temperature: 300.0,
            }],
            &[],
        )
        .unwrap_err();
        assert!(matches!(err, ThermalSolverError::BadConductivity(_)));
    }

    #[test]
    fn solve_rejects_missing_dirichlet_bc() {
        // Only Neumann loads → singular steady problem.
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box params");
        let err = solve_steady_thermal(&mesh, 50.0, &[], &[]).unwrap_err();
        assert!(matches!(err, ThermalSolverError::NoDirichletBc));
    }

    #[test]
    fn solve_rejects_mesh_without_tets() {
        let mut mesh = Mesh::new("surface-only");
        mesh.nodes.push(Vector3::new(0.0, 0.0, 0.0));
        mesh.element_blocks.push(ElementBlock {
            element_type: ElementType::Tri3,
            connectivity: vec![0, 0, 0],
        });
        let err = solve_steady_thermal(
            &mesh,
            50.0,
            &[FixedTemperature {
                node: 0,
                temperature: 300.0,
            }],
            &[],
        )
        .unwrap_err();
        assert!(matches!(err, ThermalSolverError::NoTetBlock));
    }

    #[test]
    fn uniform_dirichlet_gives_uniform_temperature() {
        // Hold every node at the same temperature → the whole body is
        // at that temperature and the flux is zero everywhere.
        let (nx, ny, nz) = (2, 2, 2);
        let mesh = structured_box_mesh(2.0, 2.0, 2.0, nx, ny, nz).expect("valid box params");
        let fixed: Vec<FixedTemperature> = (0..mesh.nodes.len())
            .map(|n| FixedTemperature {
                node: n,
                temperature: 350.0,
            })
            .collect();
        let sol = solve_steady_thermal(&mesh, 50.0, &fixed, &[]).unwrap();
        for &t in &sol.temperature {
            assert!((t - 350.0).abs() < 1e-3, "temperature {t} ≠ 350");
        }
        assert!(
            sol.max_flux_magnitude() < 1e-3,
            "uniform field should carry no flux"
        );
    }

    #[test]
    fn one_dimensional_conduction_reproduces_linear_gradient() {
        // A bar along X: hold the x=0 face at T_cold, the x=L face at
        // T_hot, leave the sides adiabatic (no BC = natural zero-flux).
        // The steady solution is the analytic *linear* temperature
        // profile  T(x) = T_cold + (T_hot − T_cold)·x/L,  and the heat
        // flux is the constant  q = −k·(T_hot − T_cold)/L.
        let (lx, ly, lz) = (4.0, 1.0, 1.0);
        let (nx, ny, nz) = (4, 1, 1);
        let mesh = structured_box_mesh(lx, ly, lz, nx, ny, nz).expect("valid box params");
        let k = 50.0;
        let t_cold = 300.0;
        let t_hot = 500.0;

        let mut fixed = Vec::new();
        for j in 0..=ny {
            for kk in 0..=nz {
                fixed.push(FixedTemperature {
                    node: nid(0, j, kk, nx, ny),
                    temperature: t_cold,
                });
                fixed.push(FixedTemperature {
                    node: nid(nx, j, kk, nx, ny),
                    temperature: t_hot,
                });
            }
        }
        let sol = solve_steady_thermal(&mesh, k, &fixed, &[]).unwrap();

        // Check the analytic linear profile at every node.
        for i in 0..=nx {
            for j in 0..=ny {
                for kk in 0..=nz {
                    let n = nid(i, j, kk, nx, ny);
                    let x = mesh.nodes[n].x;
                    let expected = t_cold + (t_hot - t_cold) * x / lx;
                    let got = sol.temperature[n];
                    assert!(
                        (got - expected).abs() < 1e-2,
                        "node {n} at x={x}: T {got} vs analytic {expected}"
                    );
                }
            }
        }
        // The flux must be (close to) the constant −k·ΔT/L along X.
        let expected_flux = -k * (t_hot - t_cold) / lx;
        let interior = nid(2, 0, 0, nx, ny);
        let qx = sol.flux[interior][0];
        assert!(
            (qx - expected_flux).abs() < 1e-2 * expected_flux.abs(),
            "interior flux qx {qx} vs analytic {expected_flux}"
        );
        // Heat flows from hot to cold: qx is negative (toward −X here
        // because hot is at +X).
        assert!(qx < 0.0, "heat should flow from the hot (+X) end");
    }

    #[test]
    fn heat_source_raises_temperature_above_the_fixed_boundary() {
        // Hold the whole x=0 face cold, inject heat at the opposite
        // face. With heat flowing in and only one cold sink, the
        // interior must be hotter than the fixed boundary.
        let (lx, ly, lz) = (4.0, 1.0, 1.0);
        let (nx, ny, nz) = (4, 1, 1);
        let mesh = structured_box_mesh(lx, ly, lz, nx, ny, nz).expect("valid box params");
        let k = 50.0;
        let t_cold = 300.0;

        let mut fixed = Vec::new();
        for j in 0..=ny {
            for kk in 0..=nz {
                fixed.push(FixedTemperature {
                    node: nid(0, j, kk, nx, ny),
                    temperature: t_cold,
                });
            }
        }
        // Inject 1 kW total, split over the x=L face's 4 nodes.
        let face: Vec<usize> = (0..=ny)
            .flat_map(|j| (0..=nz).map(move |kk| (j, kk)))
            .map(|(j, kk)| nid(nx, j, kk, nx, ny))
            .collect();
        let per = 1000.0 / face.len() as f64;
        let loads: Vec<HeatLoad> = face
            .iter()
            .map(|&n| HeatLoad {
                node: n,
                power: per,
            })
            .collect();

        let sol = solve_steady_thermal(&mesh, k, &fixed, &loads).unwrap();
        // The hot face must end up above the cold boundary.
        let hot_t = face
            .iter()
            .map(|&n| sol.temperature[n])
            .fold(f64::MIN, f64::max);
        assert!(
            hot_t > t_cold,
            "heated face {hot_t} should be hotter than the cold boundary {t_cold}"
        );
        // The cold boundary stayed put.
        let cold_n = nid(0, 0, 0, nx, ny);
        assert!(
            (sol.temperature[cold_n] - t_cold).abs() < 1e-2,
            "cold boundary moved: {}",
            sol.temperature[cold_n]
        );
        // Conservation sanity: the steady total heat in = heat out, so
        // the temperature rise ΔT relates to the conducted power.
        // ΔT_analytic = P·L/(k·A) for 1-D conduction with one cold end.
        let area = ly * lz;
        let dt_analytic = 1000.0 * lx / (k * area);
        let dt_fe = hot_t - t_cold;
        assert!(
            (dt_fe - dt_analytic).abs() < 0.1 * dt_analytic,
            "ΔT {dt_fe} vs analytic {dt_analytic}"
        );
    }

    #[test]
    fn solve_rejects_non_finite_heat_load() {
        // R3: a HeatLoad power is folded straight into the RHS vector q;
        // a non-finite power makes q non-finite, and the (pure-arithmetic)
        // Cholesky back-substitution then returns Ok(..) with a silently
        // NaN/Inf temperature field. Reject it up front instead.
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box params");
        let fixed = [FixedTemperature {
            node: 0,
            temperature: 300.0,
        }];
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let err = solve_steady_thermal(
                &mesh,
                50.0,
                &fixed,
                &[HeatLoad {
                    node: 1,
                    power: bad,
                }],
            )
            .unwrap_err();
            assert!(
                matches!(
                    err,
                    ThermalSolverError::InvalidLoad {
                        node: 1,
                        kind: "heat power"
                    }
                ),
                "expected InvalidLoad(heat power) for power {bad}, got {err:?}"
            );
        }
    }

    #[test]
    fn solve_rejects_non_finite_fixed_temperature() {
        // R3: a FixedTemperature is imposed by the penalty method, adding
        // `penalty·temperature` to the RHS entry; a non-finite prescribed
        // temperature corrupts the solve into a silently-NaN Ok(..) field.
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box params");
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let err = solve_steady_thermal(
                &mesh,
                50.0,
                &[FixedTemperature {
                    node: 0,
                    temperature: bad,
                }],
                &[],
            )
            .unwrap_err();
            assert!(
                matches!(
                    err,
                    ThermalSolverError::InvalidLoad {
                        node: 0,
                        kind: "prescribed temperature"
                    }
                ),
                "expected InvalidLoad(prescribed temperature) for T {bad}, got {err:?}"
            );
        }
    }

    #[test]
    fn solution_extrema_helpers_report_min_max() {
        let sol = ThermalSolution {
            temperature: vec![300.0, 450.0, 375.0],
            flux: vec![[3.0, 4.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        };
        assert!((sol.max_temperature() - 450.0).abs() < 1e-9);
        assert!((sol.min_temperature() - 300.0).abs() < 1e-9);
        assert!((sol.max_flux_magnitude() - 5.0).abs() < 1e-9);
    }

    #[test]
    fn convection_bar_matches_series_resistance() {
        // A bar: x=0 face fixed hot, x=L face convecting to ambient. The steady
        // heat flow is the series conduction+convection resistance
        //   q = (T_hot − T∞) / (L/(kA) + 1/(hA)),
        // and the convecting face sits at T∞ + q/(hA).
        let (lx, ly, lz) = (4.0, 1.0, 1.0);
        let (nx, ny, nz) = (4, 1, 1);
        let mesh = structured_box_mesh(lx, ly, lz, nx, ny, nz).expect("valid box");
        let k = 50.0;
        let area = ly * lz;
        let (t_hot, t_inf) = (500.0, 300.0);
        let h_area_total = 25.0; // h·A (W/K)

        let mut fixed = Vec::new();
        for j in 0..=ny {
            for kk in 0..=nz {
                fixed.push(FixedTemperature {
                    node: nid(0, j, kk, nx, ny),
                    temperature: t_hot,
                });
            }
        }
        let face: Vec<usize> = (0..=ny)
            .flat_map(|j| (0..=nz).map(move |kk| nid(nx, j, kk, nx, ny)))
            .collect();
        let conv: Vec<ConvectionBc> = face
            .iter()
            .map(|&n| ConvectionBc {
                node: n,
                h_area: h_area_total / face.len() as f64,
                ambient: t_inf,
            })
            .collect();
        let sol = solve_steady_thermal_with_convection(&mesh, k, &fixed, &[], &conv).unwrap();

        let r_cond = lx / (k * area);
        let r_conv = 1.0 / h_area_total;
        let q = (t_hot - t_inf) / (r_cond + r_conv);
        let t_face_analytic = t_inf + q * r_conv;
        let t_face_fe = face.iter().map(|&n| sol.temperature[n]).sum::<f64>() / face.len() as f64;
        assert!(
            (t_face_fe - t_face_analytic).abs() < 0.5,
            "convecting-face T {t_face_fe} vs series-resistance {t_face_analytic}"
        );
        assert!(
            (sol.temperature[nid(0, 0, 0, nx, ny)] - t_hot).abs() < 1e-2,
            "fixed face moved"
        );
    }

    #[test]
    fn large_convection_coefficient_pins_face_to_ambient() {
        // h → ∞ is the Dirichlet limit: the convecting face is pinned to T∞.
        let (nx, ny, nz) = (4, 1, 1);
        let mesh = structured_box_mesh(4.0, 1.0, 1.0, nx, ny, nz).expect("valid box");
        let (t_hot, t_inf) = (500.0, 300.0);
        let mut fixed = Vec::new();
        for j in 0..=ny {
            for kk in 0..=nz {
                fixed.push(FixedTemperature {
                    node: nid(0, j, kk, nx, ny),
                    temperature: t_hot,
                });
            }
        }
        let face: Vec<usize> = (0..=ny)
            .flat_map(|j| (0..=nz).map(move |kk| nid(nx, j, kk, nx, ny)))
            .collect();
        let conv: Vec<ConvectionBc> = face
            .iter()
            .map(|&n| ConvectionBc {
                node: n,
                h_area: 1.0e7,
                ambient: t_inf,
            })
            .collect();
        let sol = solve_steady_thermal_with_convection(&mesh, 50.0, &fixed, &[], &conv).unwrap();
        for &n in &face {
            assert!(
                (sol.temperature[n] - t_inf).abs() < 0.5,
                "node {n} T {} should be pinned to ambient {t_inf}",
                sol.temperature[n]
            );
        }
    }

    #[test]
    fn convection_alone_is_non_singular_and_conserves_heat() {
        // No Dirichlet BC, but a convection BC regularises the system: a heated
        // body convecting to ambient reaches a finite steady field, and at
        // steady state all injected heat leaves by convection (heat balance).
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 2, 2, 2).expect("valid box");
        let n = mesh.nodes.len();
        let h_area = 10.0;
        let t_inf = 300.0;
        let conv: Vec<ConvectionBc> = (0..n)
            .map(|i| ConvectionBc {
                node: i,
                h_area,
                ambient: t_inf,
            })
            .collect();
        let loads = [HeatLoad {
            node: 0,
            power: 1000.0,
        }];
        // Would return NoDirichletBc without convection; here it must solve.
        let sol = solve_steady_thermal_with_convection(&mesh, 50.0, &[], &loads, &conv).unwrap();
        for &t in &sol.temperature {
            assert!(t.is_finite() && t >= t_inf - 1e-6, "temperature {t}");
        }
        let convective_loss: f64 = (0..n).map(|i| h_area * (sol.temperature[i] - t_inf)).sum();
        assert!(
            (convective_loss - 1000.0).abs() < 1.0,
            "convective loss {convective_loss} W vs injected 1000 W (heat balance)"
        );
    }
}
