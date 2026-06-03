//! Native in-process **linear (eigenvalue) buckling** FEA solver.
//!
//! ## What this is
//!
//! A genuine, self-contained **linear buckling** solver for 4-node
//! linear tetrahedra. A linear-static solve ([`crate::native_solver`])
//! predicts the deflection under a load and is happy to return a
//! finite, "stable" answer for a slender column under *any* axial
//! load — but a real slender column does not just compress: past a
//! critical load it **buckles**, snapping sideways into a bowed shape.
//! Linear buckling analysis finds that critical load.
//!
//! The key physical idea is the **geometric (stress) stiffness**. When
//! a structure carries an internal stress, that stress changes its
//! *apparent* stiffness against further deformation: a member in
//! tension is stiffened (think of a guitar string), a member in
//! compression is **softened**. Buckling is the load at which a
//! compressive stress state has softened the structure so much that
//! its effective stiffness vanishes along some deformation pattern —
//! at that load the structure can deflect into that pattern with no
//! extra force, which is the onset of instability.
//!
//! ## The eigenvalue problem
//!
//! Linear buckling formalises that as a **generalised eigenvalue
//! problem**. From a reference linear-static solve under a reference
//! load `f_ref`, recover the reference stress state and assemble the
//! **geometric stiffness matrix** `K_g` it produces ([`geometric_stiffness`]).
//! The structure's tangent stiffness under `λ·f_ref` is then
//! `K + λ·K_g`. Buckling is the load factor `λ` at which that tangent
//! becomes singular along some mode `φ`:
//!
//! ```text
//!   (K + λ·K_g)·φ = 0
//! ```
//!
//! The smallest positive eigenvalue `λ_cr` is the **critical load
//! factor** — the buckling load is `λ_cr · f_ref` — and the companion
//! `φ` is the **buckling mode shape** (the bowed form the structure
//! collapses into).
//!
//! ## Method
//!
//! 1. Run a reference linear-static solve under `f_ref`
//!    ([`crate::native_solver::solve_linear_static`]); recover the
//!    per-element reference stress.
//! 2. Assemble `K` (the elastic stiffness) and `K_g` (the geometric
//!    stiffness from the reference stress); apply the displacement
//!    constraints by **DOF elimination** — the correct treatment for an
//!    eigenproblem.
//! 3. `K` is SPD; Cholesky-factor it, `K = L·Lᵀ`. The substitution
//!    `φ = L⁻ᵀ·ψ` turns `(K + λ·K_g)·φ = 0` into the **symmetric
//!    standard** eigenproblem `C·ψ = θ·ψ` with `C = L⁻¹·K_g·L⁻ᵀ` and
//!    `θ = −1/λ`.
//! 4. Solve `C` with the same symmetric QR eigensolver the modal
//!    solver uses ([`nalgebra::SymmetricEigen`]); convert each
//!    eigenvalue back, `λ = −1/θ`, and keep the smallest positive ones
//!    — those are the buckling load factors.
//!
//! This is the same generalised-eigenvalue machinery as
//! [`crate::modal_solver`], with the consistent mass `M` swapped for
//! the geometric stiffness `K_g`.
//!
//! ## Honest scope
//!
//! This is a **real linear-buckling v1** — the [tests](self) verify
//! that a slender column's lowest buckling load approaches the analytic
//! **Euler critical load** `P_cr = π²·E·I/(K·L)²`. It is deliberately a
//! v1:
//!
//! - **Linear (eigenvalue) buckling.** It finds the *bifurcation* load
//!   of the perfect structure from one reference stress state. It does
//!   **not** trace the post-buckling path, handle imperfection
//!   sensitivity, or do a nonlinear collapse analysis — those need an
//!   arc-length nonlinear solver.
//! - **Tet4, isotropic, small pre-buckling deformation** — the
//!   reference stress is taken from a *linear* static solve. Linear
//!   tets stiffen in bending (mild shear locking), so the FE buckling
//!   load sits *above* the slender-beam Euler value by a
//!   mesh-dependent factor — the verification test bounds it within a
//!   sane range rather than asserting an exact match.
//! - **Consistent geometric stiffness** for the constant-strain tet.
//! - End conditions are imposed with whole-face displacement
//!   constraints; the verification uses the **fixed-free** ("cantilever
//!   strut", `K = 2`) Euler case, which is exactly representable that
//!   way. A genuine *pinned* end (the face free to rotate) needs a
//!   multi-point coupling — a documented follow-up.
//!
//! Each of those is a documented, well-understood extension; none
//! affects the correctness of the bifurcation load this module
//! computes. It is *not* a CalculiX replacement.

use nalgebra::{DMatrix, Matrix3, Vector3};
use thiserror::Error;

use valenx_mesh::element::ElementType;
use valenx_mesh::Mesh;

use crate::material::FemMaterial;
use crate::native_solver::{
    check_dense_dofs, elasticity_matrix, element_stiffness, shape_gradients,
    solve_linear_static, tet_volume, NativeSolverError, NodalConstraint, NodalForce,
};

/// Errors from the native buckling solver.
#[derive(Debug, Error)]
pub enum BucklingSolverError {
    /// The underlying element assembly or the reference linear-static
    /// solve failed. Wraps the [`NativeSolverError`] raised.
    #[error("buckling pre-solve failed: {0}")]
    Native(#[from] NativeSolverError),
    /// Fewer modes were requested than 1, or more than the constrained
    /// system can supply.
    #[error("requested {requested} buckling modes but the constrained system has only {available} DOFs")]
    TooManyModes {
        /// Modes the caller asked for.
        requested: usize,
        /// Free DOFs left after the constraints.
        available: usize,
    },
    /// Every DOF was constrained — nothing is free to buckle.
    #[error("all degrees of freedom are constrained — nothing can buckle")]
    FullyConstrained,
    /// The elastic stiffness was not positive-definite, so the Cholesky
    /// transform to standard form failed (an under-constrained model).
    #[error("elastic stiffness is not positive-definite (model is under-constrained)")]
    StiffnessNotPositiveDefinite,
    /// The symmetric eigensolver did not converge.
    #[error("symmetric eigensolver failed to converge")]
    EigenFailed,
    /// No positive buckling load factor was found — the reference load
    /// does not destabilise the structure in any extracted mode (e.g.
    /// the reference load is tensile, which only stiffens).
    #[error("no positive buckling load factor found (the reference load may be stabilising)")]
    NoPositiveBucklingLoad,
}

/// One linear-buckling mode.
#[derive(Clone, Debug)]
pub struct BucklingMode {
    /// The **critical load factor** `λ` — the buckling load is
    /// `λ · f_ref`, the reference load scaled by this factor.
    pub load_factor: f64,
    /// The buckling **mode shape**: per-node relative displacement
    /// `[ux, uy, uz]`, indexed by node. The bowed form the structure
    /// snaps into; its absolute scale is arbitrary (normalised so the
    /// largest nodal magnitude is 1).
    pub shape: Vec<[f64; 3]>,
}

impl BucklingMode {
    /// Largest displacement magnitude across the mode shape.
    pub fn max_amplitude(&self) -> f64 {
        self.shape
            .iter()
            .map(|d| (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt())
            .fold(0.0, f64::max)
    }
}

/// Result of a linear-buckling solve — the lowest critical load
/// factors, ascending.
#[derive(Clone, Debug)]
pub struct BucklingSolution {
    /// The buckling modes, sorted by ascending critical load factor.
    /// `modes[0]` is the first (lowest-load) buckling mode.
    pub modes: Vec<BucklingMode>,
}

impl BucklingSolution {
    /// The **critical load factor** of the first buckling mode — the
    /// reference load multiplier at which the structure first buckles.
    /// `None` if no modes were extracted.
    pub fn critical_load_factor(&self) -> Option<f64> {
        self.modes.first().map(|m| m.load_factor)
    }
}

/// Solve a linear (eigenvalue) buckling problem on a tetrahedral mesh.
///
/// `mesh` carries the Tet4 body; `material` its isotropic elastic
/// constants; `constraints` the displacement boundary conditions (a
/// pin-ended column is pinned at both ends). `reference_forces` is the
/// **reference load** — the buckling load comes out as a multiple of
/// it. `n_modes` is how many of the lowest buckling modes to return.
///
/// Returns the lowest `n_modes` [`BucklingMode`]s, ascending in
/// critical load factor.
///
/// # Method
///
/// See the [module documentation](self): a reference linear-static
/// solve recovers the stress state, the geometric stiffness `K_g` is
/// assembled from it, and the generalised eigenproblem
/// `(K + λ·K_g)·φ = 0` is solved by transforming to the symmetric
/// standard form via the Cholesky factor of `K`.
///
/// # Errors
///
/// See [`BucklingSolverError`].
pub fn solve_buckling(
    mesh: &Mesh,
    material: &FemMaterial,
    constraints: &[NodalConstraint],
    reference_forces: &[NodalForce],
    n_modes: usize,
) -> Result<BucklingSolution, BucklingSolverError> {
    if n_modes == 0 {
        return Err(BucklingSolverError::TooManyModes {
            requested: 0,
            available: 0,
        });
    }
    let n_nodes = mesh.nodes.len();
    if n_nodes == 0 {
        return Err(BucklingSolverError::Native(
            NativeSolverError::EmptyMesh,
        ));
    }
    // Round-1 H1–H3: buckling builds TWO dense `n_dof × n_dof` matrices
    // (K and K_g) with NO prior cap. Reject an oversized mesh BEFORE the
    // allocations below — and before the (cheaper, sparse) reference
    // solve does pointless work for a model that can never be assembled
    // densely. `check_dense_dofs` is the shared, tested guard.
    let n_dof = check_dense_dofs(n_nodes)?;

    // --- (1) reference linear-static solve → reference stress state ---
    let static_sol =
        solve_linear_static(mesh, material, constraints, reference_forces)?;

    // --- (2) assemble the elastic K and the geometric K_g ---
    let d_matrix = elasticity_matrix(material)?;
    let tets = mesh
        .element_blocks
        .iter()
        .find(|b| b.element_type == ElementType::Tet4)
        .filter(|b| b.connectivity.len() >= 4)
        .ok_or(BucklingSolverError::Native(NativeSolverError::NoTetBlock))?;
    let conn = &tets.connectivity;
    let n_elem = conn.len() / 4;

    let mut k = DMatrix::<f64>::zeros(n_dof, n_dof);
    let mut kg = DMatrix::<f64>::zeros(n_dof, n_dof);
    for e in 0..n_elem {
        let nodes = [
            conn[4 * e] as usize,
            conn[4 * e + 1] as usize,
            conn[4 * e + 2] as usize,
            conn[4 * e + 3] as usize,
        ];
        for &nd in &nodes {
            if nd >= n_nodes {
                return Err(BucklingSolverError::Native(
                    NativeSolverError::BadConnectivity {
                        elem: e,
                        node: nd,
                        n_nodes,
                    },
                ));
            }
        }
        let coords = [
            mesh.nodes[nodes[0]],
            mesh.nodes[nodes[1]],
            mesh.nodes[nodes[2]],
            mesh.nodes[nodes[3]],
        ];
        // Element elastic stiffness.
        let ke = element_stiffness(&coords, &d_matrix).ok_or(
            BucklingSolverError::Native(NativeSolverError::DegenerateElement(
                e,
            )),
        )?;
        // The element's reference stress — the average of its four
        // nodal stresses (the linear-static solve recovered nodal
        // stress; for a constant-strain tet the element stress is
        // recovered the same way and the nodal average reproduces it).
        let mut elem_stress = [0.0_f64; 6];
        for &nd in &nodes {
            for (k_i, es) in elem_stress.iter_mut().enumerate() {
                *es += static_sol.stress[nd][k_i] / 4.0;
            }
        }
        // Element geometric stiffness from that reference stress.
        let kge = geometric_stiffness(&coords, &elem_stress).ok_or(
            BucklingSolverError::Native(NativeSolverError::DegenerateElement(
                e,
            )),
        )?;
        // Scatter both 12×12 matrices into the global system.
        for a in 0..4 {
            for i in 0..3 {
                let gi = 3 * nodes[a] + i;
                for b in 0..4 {
                    for j in 0..3 {
                        let gj = 3 * nodes[b] + j;
                        k[(gi, gj)] += ke[(3 * a + i, 3 * b + j)];
                        kg[(gi, gj)] += kge[(3 * a + i, 3 * b + j)];
                    }
                }
            }
        }
    }

    // --- which DOFs are free? constrained DOFs eliminated ---
    let mut constrained = vec![false; n_dof];
    for c in constraints {
        if c.node >= n_nodes {
            return Err(BucklingSolverError::Native(
                NativeSolverError::BadConnectivity {
                    elem: usize::MAX,
                    node: c.node,
                    n_nodes,
                },
            ));
        }
        for (i, fixed) in c.fixed.iter().enumerate() {
            if fixed.is_some() {
                constrained[3 * c.node + i] = true;
            }
        }
    }
    let free: Vec<usize> = (0..n_dof).filter(|&d| !constrained[d]).collect();
    if free.is_empty() {
        return Err(BucklingSolverError::FullyConstrained);
    }
    let n_free = free.len();
    if n_modes > n_free {
        return Err(BucklingSolverError::TooManyModes {
            requested: n_modes,
            available: n_free,
        });
    }

    // --- reduce K and K_g to the free DOFs ---
    let mut k_ff = DMatrix::<f64>::zeros(n_free, n_free);
    let mut kg_ff = DMatrix::<f64>::zeros(n_free, n_free);
    for (ri, &gr) in free.iter().enumerate() {
        for (ci, &gc) in free.iter().enumerate() {
            k_ff[(ri, ci)] = k[(gr, gc)];
            kg_ff[(ri, ci)] = kg[(gr, gc)];
        }
    }

    // --- generalised → standard via the Cholesky factor of K_ff ---
    // K_ff = L·Lᵀ. The substitution φ = L⁻ᵀ·ψ turns
    // (K + λ·K_g)·φ = 0  →  K φ = −λ K_g φ  into the symmetric
    // standard problem C·ψ = θ·ψ with C = L⁻¹·K_g·L⁻ᵀ and θ = −1/λ.
    let chol = k_ff
        .clone()
        .cholesky()
        .ok_or(BucklingSolverError::StiffnessNotPositiveDefinite)?;
    let l = chol.l();
    let l_inv = invert_lower_triangular(&l)
        .ok_or(BucklingSolverError::StiffnessNotPositiveDefinite)?;
    let l_inv_t = l_inv.transpose();
    let mut c = &l_inv * &kg_ff * &l_inv_t;
    // Symmetrise to kill floating-point asymmetry.
    let c_t = c.transpose();
    c = (&c + &c_t) * 0.5;

    // --- symmetric eigensolve of C ---
    let eigen = nalgebra::SymmetricEigen::try_new(c, 1.0e-12, 0)
        .ok_or(BucklingSolverError::EigenFailed)?;
    let eigvals = &eigen.eigenvalues;
    let eigvecs = &eigen.eigenvectors;

    // --- convert θ → λ and keep the positive load factors ---
    // θ = −1/λ  ⇒  λ = −1/θ. A buckling load factor must be positive;
    // a non-positive λ means that mode is stabilised (or unaffected)
    // by the reference load — discard it. (θ = 0 ⇒ λ = ∞, never
    // buckles; θ > 0 ⇒ λ < 0, the load would have to reverse.)
    let mut candidates: Vec<(f64, usize)> = Vec::new();
    for idx in 0..n_free {
        let theta = eigvals[idx];
        if theta < -1e-14 {
            // θ < 0 ⇒ λ = −1/θ > 0 — a genuine buckling factor.
            let lambda = -1.0 / theta;
            if lambda.is_finite() && lambda > 0.0 {
                candidates.push((lambda, idx));
            }
        }
    }
    if candidates.is_empty() {
        return Err(BucklingSolverError::NoPositiveBucklingLoad);
    }
    // Sort ascending by load factor — the lowest buckles first.
    candidates.sort_by(|a, b| {
        a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal)
    });

    // --- build the lowest n_modes BucklingModes ---
    let take = n_modes.min(candidates.len());
    let mut modes = Vec::with_capacity(take);
    for &(lambda, idx) in candidates.iter().take(take) {
        // Recover the buckling mode φ = L⁻ᵀ·ψ on the free DOFs.
        let psi = eigvecs.column(idx).into_owned();
        let phi_free = &l_inv_t * &psi;
        // Scatter into a full 3·n_nodes shape; normalise so the largest
        // nodal magnitude is 1 (the absolute scale of a buckling mode
        // is arbitrary).
        let mut shape = vec![[0.0_f64; 3]; n_nodes];
        for (fi, &gd) in free.iter().enumerate() {
            shape[gd / 3][gd % 3] = phi_free[fi];
        }
        let max_mag = shape
            .iter()
            .map(|d| (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt())
            .fold(0.0_f64, f64::max);
        if max_mag > 0.0 {
            for d in shape.iter_mut() {
                d[0] /= max_mag;
                d[1] /= max_mag;
                d[2] /= max_mag;
            }
        }
        modes.push(BucklingMode {
            load_factor: lambda,
            shape,
        });
    }

    Ok(BucklingSolution { modes })
}

/// The 12×12 **geometric (stress) stiffness matrix** `K_g` of a linear
/// tetrahedron carrying a constant reference stress.
///
/// The geometric stiffness captures how an existing internal stress
/// changes the element's stiffness against further deformation — the
/// term that makes a compressed member soften and so makes buckling
/// possible. For a constant-strain tet with constant shape-function
/// gradients `∇Nₐ` and a constant Cauchy stress `σ` (passed in Voigt
/// order `[σxx σyy σzz σxy σyz σzx]`), it is the closed form
///
/// ```text
///   K_g[3a+k, 3b+k] = Vₑ · (∇Nₐ)ᵀ · σ · (∇Nᵦ)
/// ```
///
/// for each spatial component `k` (the geometric stiffness couples the
/// *same* component of two nodes — the `δ_kl` structure), with `σ` the
/// 3×3 stress tensor and `Vₑ` the element volume. Returns `None` for a
/// degenerate (zero-volume) tet.
pub fn geometric_stiffness(
    coords: &[Vector3<f64>; 4],
    stress_voigt: &[f64; 6],
) -> Option<DMatrix<f64>> {
    let grads = shape_gradients(coords)?;
    let vol = tet_volume(coords).abs();
    if vol < 1.0e-18 {
        return None;
    }
    // The 3×3 Cauchy stress tensor from the Voigt vector.
    let s = Matrix3::new(
        stress_voigt[0],
        stress_voigt[3],
        stress_voigt[5],
        stress_voigt[3],
        stress_voigt[1],
        stress_voigt[4],
        stress_voigt[5],
        stress_voigt[4],
        stress_voigt[2],
    );
    let mut kg = DMatrix::<f64>::zeros(12, 12);
    for a in 0..4 {
        for b in 0..4 {
            // Scalar coupling g_ab = Vₑ · (∇Nₐ)ᵀ·σ·(∇Nᵦ).
            let g_ab = vol * grads[a].dot(&(s * grads[b]));
            // It loads the diagonal of every spatial component pair.
            for k in 0..3 {
                kg[(3 * a + k, 3 * b + k)] += g_ab;
            }
        }
    }
    Some(kg)
}

/// Invert a lower-triangular matrix `L` by forward substitution —
/// the same routine [`crate::modal_solver`] uses for its
/// generalised → standard transform, reproduced here so the buckling
/// module is self-contained. Returns `None` for a singular factor.
fn invert_lower_triangular(l: &DMatrix<f64>) -> Option<DMatrix<f64>> {
    let n = l.nrows();
    debug_assert_eq!(n, l.ncols(), "L must be square");
    let mut inv = DMatrix::<f64>::zeros(n, n);
    for col in 0..n {
        for row in 0..n {
            let mut sum = if row == col { 1.0 } else { 0.0 };
            for k in 0..row {
                sum -= l[(row, k)] * inv[(k, col)];
            }
            let diag = l[(row, row)];
            if diag.abs() < 1.0e-300 {
                return None;
            }
            inv[(row, col)] = sum / diag;
        }
    }
    Some(inv)
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
    fn geometric_stiffness_is_symmetric() {
        // K_g is a symmetric matrix.
        let coords = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        // An arbitrary stress state.
        let stress = [1e6, -2e6, 3e6, 0.5e6, -1e6, 0.2e6];
        let kg = geometric_stiffness(&coords, &stress).unwrap();
        assert_eq!(kg.nrows(), 12);
        for i in 0..12 {
            for j in 0..12 {
                assert!(
                    (kg[(i, j)] - kg[(j, i)]).abs()
                        < 1e-3 * kg[(i, i)].abs().max(1.0),
                    "K_g must be symmetric"
                );
            }
        }
    }

    #[test]
    fn geometric_stiffness_vanishes_with_no_stress() {
        // No reference stress → no geometric stiffness.
        let coords = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(2.0, 0.0, 0.0),
            Vector3::new(0.0, 1.5, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let kg = geometric_stiffness(&coords, &[0.0; 6]).unwrap();
        for v in kg.iter() {
            assert!(v.abs() < 1e-9, "K_g must vanish for a stress-free element");
        }
    }

    #[test]
    fn geometric_stiffness_flips_sign_with_the_stress() {
        // K_g is linear in the stress — reversing the stress (tension ↔
        // compression) negates the geometric stiffness. This is the
        // physical asymmetry that makes a compressed member soften and
        // a tensioned one stiffen.
        let coords = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let tension = [5e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        let compression = [-5e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        let kg_t = geometric_stiffness(&coords, &tension).unwrap();
        let kg_c = geometric_stiffness(&coords, &compression).unwrap();
        for i in 0..12 {
            for j in 0..12 {
                assert!(
                    (kg_t[(i, j)] + kg_c[(i, j)]).abs()
                        < 1e-3 * kg_t[(i, i)].abs().max(1.0),
                    "K_g(−σ) must equal −K_g(σ)"
                );
            }
        }
    }

    #[test]
    fn degenerate_tet_has_no_geometric_stiffness() {
        // Four coplanar points → zero volume → None.
        let coords = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
        ];
        assert!(geometric_stiffness(&coords, &[1e6; 6]).is_none());
    }

    #[test]
    fn rejects_zero_modes() {
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box params");
        let err = solve_buckling(
            &mesh,
            &FemMaterial::default(),
            &[NodalConstraint::fixed(0)],
            &[],
            0,
        )
        .unwrap_err();
        assert!(matches!(err, BucklingSolverError::TooManyModes { .. }));
    }

    #[test]
    fn rejects_too_many_modes() {
        // A single hex (8 nodes, 24 DOF). The x=0 face (nodes 0,2,4,6 in
        // this mesh's ordering — every node with i=0) is fully clamped
        // so the buckling solver's reference linear-static pre-solve is
        // well-posed; the ONLY thing wrong here is the absurd mode count.
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box params");
        let mut constraints = Vec::new();
        for k in 0..=1 {
            for j in 0..=1 {
                constraints.push(NodalConstraint::fixed(nid(0, j, k, 1, 1)));
            }
        }
        let err = solve_buckling(
            &mesh,
            &FemMaterial::default(),
            &constraints,
            &[NodalForce {
                node: nid(1, 1, 1, 1, 1),
                force: [-1.0, 0.0, 0.0],
            }],
            100,
        )
        .unwrap_err();
        assert!(
            matches!(
                err,
                BucklingSolverError::TooManyModes { requested: 100, .. }
            ),
            "expected TooManyModes, got {err:?}"
        );
    }

    #[test]
    fn column_buckling_load_approaches_the_euler_critical_load() {
        // The headline buckling verification. A slender column, axially
        // compressed, buckles at the **Euler critical load**
        //   P_cr = π²·E·I / (K·L)²
        // with `K` the effective-length factor of the end conditions.
        //
        // We model the unambiguous **fixed-free** ("cantilever strut")
        // case — the one end condition a solid-element column can be
        // restrained for *exactly* with whole-face constraints: the
        // bottom face is fully clamped, the top is entirely free except
        // for the axial compressive load. Its effective length factor
        // is K = 2, so the Euler load is
        //   P_cr = π²·E·I / (2L)².
        // (A genuine *pinned* end would have to let the end face rotate,
        // which needs a multi-point coupling this v1 does not carry —
        // the cantilever strut is the honest, exactly-modelable Euler
        // case. Both are Euler buckling; only `K` differs.)
        //
        // The column runs along Z, square in section (b × b), slender
        // (L ≫ b). The linear-buckling eigenvalue × the reference load
        // must approach the analytic Euler value.
        let b = 1.0;
        let l = 24.0; // slender: L/b = 24
        // Enough elements along the length to resolve the smooth
        // quarter-sine buckling mode.
        let (nx, ny, nz) = (2, 2, 24);
        let mesh = structured_box_mesh(b, b, l, nx, ny, nz).expect("valid box params");
        let e_mod = 200.0e9;
        let mat = FemMaterial {
            youngs_modulus: e_mod,
            poisson_ratio: 0.3,
            ..FemMaterial::default()
        };

        // Bottom face (z=0): fully clamped.
        let mut constraints = Vec::new();
        for j in 0..=ny {
            for i in 0..=nx {
                constraints.push(NodalConstraint::fixed(nid(i, j, 0, nx, ny)));
            }
        }
        // Top face (z=nz): completely free — only loaded, not
        // restrained. (Fixed-free / cantilever strut.)
        // Reference load: a total unit compressive force on the top
        // face, pushing in −Z.
        let top_nodes: Vec<usize> = (0..=ny)
            .flat_map(|j| (0..=nx).map(move |i| (i, j)))
            .map(|(i, j)| nid(i, j, nz, nx, ny))
            .collect();
        let total_ref = 1.0; // unit reference load
        let per = total_ref / top_nodes.len() as f64;
        let forces: Vec<NodalForce> = top_nodes
            .iter()
            .map(|&n| NodalForce {
                node: n,
                force: [0.0, 0.0, -per],
            })
            .collect();

        let sol =
            solve_buckling(&mesh, &mat, &constraints, &forces, 3).unwrap();
        let lambda = sol.critical_load_factor().unwrap();
        assert!(lambda > 0.0, "the critical load factor must be positive");
        // Since the reference load is unit, the buckling load IS the
        // load factor.
        let p_cr_fe = lambda * total_ref;

        // Analytic Euler load for the fixed-free column (K = 2):
        //   P_cr = π²·E·I / (2L)²,  I = b·b³/12 for the square section.
        let i_section = b * b.powi(3) / 12.0;
        let effective_length = 2.0 * l;
        let p_cr_euler = std::f64::consts::PI.powi(2) * e_mod * i_section
            / effective_length.powi(2);

        assert!(p_cr_euler > 0.0);
        // A linear-tet column is stiffer in bending than slender-beam
        // theory (mild shear locking — worse for a coarse cross
        // section), so the FE buckling load sits *above* the Euler
        // value. On this mesh it must land within a sane factor —
        // confirming the geometric stiffness + the eigensolve are
        // correct. The upper bound is deliberately generous: the exact
        // factor is mesh-dependent for linear tets (a documented
        // caveat), the load-direction sense and order of magnitude are
        // what this test pins.
        let ratio = p_cr_fe / p_cr_euler;
        assert!(
            ratio > 0.8,
            "FE buckling load {p_cr_fe} far below Euler {p_cr_euler} (ratio {ratio})"
        );
        assert!(
            ratio < 6.0,
            "FE buckling load {p_cr_fe} far above Euler {p_cr_euler} (ratio {ratio})"
        );

        // The buckling modes are sorted ascending in load factor.
        for w in sol.modes.windows(2) {
            assert!(
                w[1].load_factor >= w[0].load_factor - 1e-6,
                "buckling modes must be ascending"
            );
        }
        // The first buckling mode is a genuine non-trivial sideways
        // bow — it has lateral (X or Y) motion, and the clamped bottom
        // does not move.
        let m0 = &sol.modes[0];
        assert!(m0.max_amplitude() > 0.0, "buckling mode is all zero");
        let bottom = nid(nx / 2, ny / 2, 0, nx, ny);
        let bottom_amp = (m0.shape[bottom][0].powi(2)
            + m0.shape[bottom][1].powi(2)
            + m0.shape[bottom][2].powi(2))
        .sqrt();
        assert!(
            bottom_amp < 1e-6,
            "the clamped bottom must not move in the buckling mode: {bottom_amp}"
        );
        // The free tip of the cantilever strut bows sideways the most —
        // its lateral displacement must be well away from zero.
        let tip = nid(nx / 2, ny / 2, nz, nx, ny);
        let tip_lateral =
            (m0.shape[tip][0].powi(2) + m0.shape[tip][1].powi(2)).sqrt();
        assert!(
            tip_lateral > 0.3,
            "the cantilever tip should bow sideways, lateral motion {tip_lateral}"
        );
    }

    #[test]
    fn stiffer_column_buckles_at_a_higher_load() {
        // A physical trend check: a column of a stiffer material (larger
        // E) buckles at a proportionally higher load — P_cr ∝ E. We
        // check the monotone trend.
        let b = 1.0;
        let l = 16.0;
        let (nx, ny, nz) = (2, 2, 16);
        let mesh = structured_box_mesh(b, b, l, nx, ny, nz).expect("valid box params");

        let buckle_load = |e_mod: f64| -> f64 {
            let mat = FemMaterial {
                youngs_modulus: e_mod,
                poisson_ratio: 0.3,
                ..FemMaterial::default()
            };
            // Fixed-free column — clamp the bottom face, leave the top
            // free except for the compressive load.
            let mut constraints = Vec::new();
            for j in 0..=ny {
                for i in 0..=nx {
                    constraints
                        .push(NodalConstraint::fixed(nid(i, j, 0, nx, ny)));
                }
            }
            let top_nodes: Vec<usize> = (0..=ny)
                .flat_map(|j| (0..=nx).map(move |i| (i, j)))
                .map(|(i, j)| nid(i, j, nz, nx, ny))
                .collect();
            let per = 1.0 / top_nodes.len() as f64;
            let forces: Vec<NodalForce> = top_nodes
                .iter()
                .map(|&n| NodalForce {
                    node: n,
                    force: [0.0, 0.0, -per],
                })
                .collect();
            let sol =
                solve_buckling(&mesh, &mat, &constraints, &forces, 1).unwrap();
            sol.critical_load_factor().unwrap()
        };

        let soft = buckle_load(70.0e9); // aluminium-ish
        let stiff = buckle_load(200.0e9); // steel-ish
        assert!(
            stiff > soft,
            "a stiffer column should buckle at a higher load: {soft} → {stiff}"
        );
        // P_cr ∝ E, so the ratio should track the E ratio (~2.86).
        let ratio = stiff / soft;
        assert!(
            ratio > 2.0 && ratio < 4.0,
            "buckling load should scale with E: ratio {ratio} vs E-ratio 2.86"
        );
    }

    #[test]
    fn rejects_a_mesh_without_tets() {
        use valenx_mesh::element::ElementBlock;
        let mut mesh = Mesh::new("surface");
        mesh.nodes.push(Vector3::new(0.0, 0.0, 0.0));
        mesh.element_blocks.push(ElementBlock {
            element_type: ElementType::Tri3,
            connectivity: vec![0, 0, 0],
        });
        let err = solve_buckling(
            &mesh,
            &FemMaterial::default(),
            &[NodalConstraint::fixed(0)],
            &[NodalForce {
                node: 0,
                force: [-1.0, 0.0, 0.0],
            }],
            1,
        )
        .unwrap_err();
        // The reference linear-static solve fails first on a tet-less
        // mesh — surfaced as a wrapped NativeSolverError.
        assert!(matches!(
            err,
            BucklingSolverError::Native(NativeSolverError::NoTetBlock)
        ));
    }
}
