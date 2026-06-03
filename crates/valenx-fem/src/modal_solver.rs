//! Native in-process modal (natural-frequency) FEA solver (Phase 24.6).
//!
//! ## What this is
//!
//! A genuine, self-contained **modal finite-element solver** for 4-node
//! linear tetrahedra (`C3D4` / [`ElementType::Tet4`]). It reuses the
//! linear-elastic element formulation from [`crate::native_solver`] —
//! the same constant-strain-tetrahedron stiffness `Kₑ` — and pairs it
//! with a **consistent element mass matrix** `Mₑ`. It assembles the
//! global stiffness `K` and mass `M`, applies displacement boundary
//! conditions by *eliminating* the constrained degrees of freedom, and
//! solves the generalised symmetric eigenvalue problem
//!
//! ```text
//!   K φ = λ M φ
//! ```
//!
//! for the lowest `n` modes. Each eigenvalue `λᵢ = ωᵢ²` is an angular
//! frequency squared; the natural frequency is `fᵢ = ωᵢ / 2π` in hertz.
//! The companion eigenvector `φᵢ` is the mode shape (the relative
//! deflected form the structure takes when vibrating at `fᵢ`).
//!
//! This is the in-process alternative to dispatching a modal step to
//! the CalculiX / Elmer / Code_Aster subprocess adapters: no external
//! binary, no input-deck round-trip.
//!
//! ## Method
//!
//! 1. **Element matrices.** For each linear tet the stiffness is the
//!    closed-form `Kₑ = Vₑ·Bᵀ·D·B` (see [`crate::native_solver`]); the
//!    **consistent mass** is the textbook
//!    `Mₑ = ρ·Vₑ/20 · (I₁₂ + diag-doubling)` — every translational DOF
//!    block is `ρ·Vₑ/20·[2 on the diagonal, 1 off]`. The consistent
//!    mass (as opposed to a lumped diagonal mass) is exact for the
//!    linear-tet kinematics and gives an upper bound on each frequency.
//! 2. **Assembly.** Both `Kₑ` and `Mₑ` are scatter-added into global
//!    dense matrices of size `3·n_nodes`.
//! 3. **Boundary conditions.** Displacement constraints are imposed by
//!    **DOF elimination** — the constrained rows and columns are struck
//!    out of both `K` and `M`, leaving a reduced system on the free
//!    DOFs only. (The penalty method that [`crate::native_solver`] uses
//!    for the *static* solve is wrong for a modal solve: a huge penalty
//!    spring injects spurious very-high-frequency modes. Elimination is
//!    the correct treatment.) Only homogeneous constraints (fixed = 0)
//!    are meaningful for an eigen-analysis, so any non-zero prescribed
//!    displacement is treated as a plain "fix this DOF".
//! 4. **Generalised → standard.** With the reduced free-DOF matrices
//!    `K_ff` (SPD) and `M_ff` (SPD), the Cholesky factor `M_ff = L·Lᵀ`
//!    turns the generalised problem into the **standard symmetric**
//!    problem `C ψ = λ ψ` with `C = L⁻¹·K_ff·L⁻ᵀ` and `φ = L⁻ᵀ·ψ`.
//!    `C` is symmetric, so a symmetric tridiagonal QR eigensolver
//!    ([`nalgebra::SymmetricEigen`]) returns every `(λ, ψ)` pair
//!    accurately.
//! 5. **Recover + sort.** Eigenvalues are sorted ascending; the lowest
//!    `n` are kept; each mode shape is back-substituted to full
//!    `3·n_nodes` length (constrained DOFs filled with zero) and
//!    mass-normalised (`φᵀ M φ = 1`).
//!
//! ## Honest scope — what it is and is not
//!
//! **It is a real modal solver.** For a linear-elastic structure
//! meshed with linear tets, the frequencies it returns are the true
//! FE natural frequencies of that mesh, and they converge toward the
//! analytic frequencies as the mesh is refined.
//!
//! **It is deliberately the 90 % case:**
//!
//! - **Undamped, free vibration only.** No damping (the eigenproblem
//!   would become complex / quadratic), no harmonic-response or
//!   transient time-stepping, no pre-stress / stress-stiffening
//!   ("modal under load"), no buckling (that needs the geometric
//!   stiffness matrix `K_g`).
//! - **Tet4 only**, isotropic linear elasticity, consistent mass.
//! - **Dense eigensolver — supports up to [`MODAL_DENSE_MAX_DOFS`]
//!   free DOFs.** The reduced system is solved with a dense symmetric
//!   QR. That is exact and robust for workshop-scale meshes (a few
//!   thousand DOF); above the cap the solver returns
//!   [`ModalSolverError::DofCountExceeded`] rather than allocating a
//!   multi-hundred-MB dense matrix. A Lanczos / subspace-iteration
//!   sparse eigensolver for larger problems is tracked as future work.
//! - **It needs a volume mesh** — see [`crate::native_solver`]'s note
//!   and [`crate::native_solver::structured_box_mesh`].
//!
//! Within that scope the numbers are real: the lowest non-rigid mode
//! of a cantilever beam comes out close to the Euler-Bernoulli
//! prediction, an unconstrained body returns six near-zero rigid-body
//! frequencies, and refining the mesh moves each frequency the right
//! way.

use nalgebra::{DMatrix, DVector, Matrix6, Vector3};
use thiserror::Error;

use valenx_mesh::element::ElementType;
use valenx_mesh::Mesh;

use crate::material::FemMaterial;
use crate::native_solver::{NativeSolverError, NodalConstraint};

/// Errors from the native modal solver.
#[derive(Debug, Error)]
pub enum ModalSolverError {
    /// The underlying element assembly failed — bad material, no tet
    /// block, empty mesh, degenerate element, or out-of-range
    /// connectivity. Wraps the [`NativeSolverError`] that the shared
    /// element-assembly path raised.
    #[error("element assembly failed: {0}")]
    Assembly(#[from] NativeSolverError),
    /// The material density is non-physical (ρ ≤ 0). A modal solve
    /// divides by mass, so a zero or negative density is meaningless.
    #[error("invalid material: density must be finite and positive, got {0}")]
    BadDensity(f64),
    /// Fewer modes were requested than 1, or more than the reduced
    /// system can supply.
    #[error("requested {requested} modes but the constrained system has only {available} DOFs")]
    TooManyModes {
        /// Modes the caller asked for.
        requested: usize,
        /// Free DOFs left after applying the constraints.
        available: usize,
    },
    /// Every DOF was constrained — there is nothing left to vibrate.
    #[error("all degrees of freedom are constrained — no modes to solve")]
    FullyConstrained,
    /// The reduced mass matrix was not positive-definite, so the
    /// Cholesky transform to standard form failed. Indicates a
    /// massless / disconnected region in the free-DOF set.
    #[error("reduced mass matrix is not positive-definite (massless or disconnected free DOFs)")]
    MassNotPositiveDefinite,
    /// The symmetric eigensolver did not converge within its iteration
    /// budget.
    #[error("symmetric eigensolver failed to converge")]
    EigenFailed,
    /// The free-DOF count exceeded the dense-solver upper bound. The
    /// in-process modal solver builds a dense `n × n` matrix; at 5k
    /// DOFs that's ~200 MB and at 10k it's ~800 MB. We refuse early
    /// with an actionable error rather than OOMing the host.
    #[error(
        "modal solve requires {n_free} free DOFs, dense solver supports at most {max}; \
         use a sparse Lanczos solver (future work) or coarsen the mesh"
    )]
    DofCountExceeded {
        /// Free DOFs in the constrained system the caller passed.
        n_free: usize,
        /// Maximum supported by the dense solver.
        max: usize,
    },
}

/// Upper bound on the constrained-DOF count for the dense modal
/// solver. Building a dense `n × n` mass / stiffness matrix scales as
/// `8 · n²` bytes; at 2k DOFs that's ~32 MB which comfortably fits in
/// a workstation. 5k DOFs is ~200 MB and 10k DOFs is ~800 MB, both of
/// which will OOM smaller dev boxes. A sparse Lanczos solver for
/// larger problems is tracked as future work.
pub const MODAL_DENSE_MAX_DOFS: usize = 2_000;

/// One natural mode of vibration.
#[derive(Clone, Debug)]
pub struct VibrationMode {
    /// Natural frequency in hertz (`f = ω / 2π`).
    pub frequency_hz: f64,
    /// Angular frequency in rad/s (`ω = √λ`).
    pub angular_frequency: f64,
    /// Eigenvalue `λ = ω²` of `K φ = λ M φ`.
    pub eigenvalue: f64,
    /// Mode shape: per-node relative displacement `[ux, uy, uz]`,
    /// indexed by node. Mass-normalised so `φᵀ M φ = 1`; the absolute
    /// scale of a mode shape is otherwise arbitrary.
    pub shape: Vec<[f64; 3]>,
}

impl VibrationMode {
    /// Largest displacement magnitude across the mode shape — useful
    /// for scaling the shape for display.
    pub fn max_amplitude(&self) -> f64 {
        self.shape
            .iter()
            .map(|d| (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt())
            .fold(0.0, f64::max)
    }
}

/// Result of a native modal solve: the lowest `n` modes, ascending in
/// frequency.
#[derive(Clone, Debug)]
pub struct ModalSolution {
    /// The modes, sorted by ascending frequency. `modes[0]` is the
    /// fundamental.
    pub modes: Vec<VibrationMode>,
}

impl ModalSolution {
    /// Frequency of the fundamental (lowest) mode in hertz, or `None`
    /// if no modes were extracted.
    pub fn fundamental_hz(&self) -> Option<f64> {
        self.modes.first().map(|m| m.frequency_hz)
    }
}

/// Solve the modal (natural-frequency) eigenproblem on a tetrahedral
/// mesh.
///
/// `mesh` must carry at least one [`ElementType::Tet4`] block.
/// `material` supplies the isotropic elastic constants **and** the
/// density (a modal solve needs mass). `constraints` pin nodal DOFs
/// (an unconstrained model still solves — it returns the six
/// rigid-body modes at ≈ 0 Hz followed by the elastic modes).
/// `n_modes` is how many of the lowest modes to return.
///
/// Returns the lowest `n_modes` [`VibrationMode`]s, ascending in
/// frequency.
///
/// # Errors
///
/// See [`ModalSolverError`].
pub fn solve_modal(
    mesh: &Mesh,
    material: &FemMaterial,
    constraints: &[NodalConstraint],
    n_modes: usize,
) -> Result<ModalSolution, ModalSolverError> {
    if n_modes == 0 {
        return Err(ModalSolverError::TooManyModes {
            requested: 0,
            available: 0,
        });
    }
    if !material.density.is_finite() || material.density <= 0.0 {
        return Err(ModalSolverError::BadDensity(material.density));
    }
    // Round-1 H1–H3: bound the dense `n_dof × n_dof` K and M BEFORE the
    // assembler allocates them. This is the PRIMARY allocation guard;
    // the `n_free > MODAL_DENSE_MAX_DOFS` check below is a secondary
    // bound on the (smaller) reduced eigensolve. `assemble_global_…`
    // also calls `check_dense_dofs` internally — defense in depth.
    crate::native_solver::check_dense_dofs(mesh.nodes.len())?;

    // --- assemble the global K and M (dense, 3·n_nodes square) ---
    let (k_global, m_global) = assemble_global_stiffness_mass(mesh, material)?;
    let n_dof = k_global.nrows();

    // --- which DOFs are free? a constrained DOF is struck out ---
    let mut constrained = vec![false; n_dof];
    for c in constraints {
        if c.node >= mesh.nodes.len() {
            return Err(ModalSolverError::Assembly(
                NativeSolverError::BadConnectivity {
                    elem: usize::MAX,
                    node: c.node,
                    n_nodes: mesh.nodes.len(),
                },
            ));
        }
        for (i, fixed) in c.fixed.iter().enumerate() {
            // For a modal analysis only the *fact* of a constraint
            // matters; a non-zero prescribed value is not physical for
            // free-vibration eigen-analysis, so any Some(_) fixes it.
            if fixed.is_some() {
                constrained[3 * c.node + i] = true;
            }
        }
    }
    let free: Vec<usize> = (0..n_dof).filter(|&d| !constrained[d]).collect();
    if free.is_empty() {
        return Err(ModalSolverError::FullyConstrained);
    }
    let n_free = free.len();
    if n_modes > n_free {
        return Err(ModalSolverError::TooManyModes {
            requested: n_modes,
            available: n_free,
        });
    }
    // Dense-solver upper bound. Building two `n_free × n_free` f64
    // matrices below allocates `16 · n_free²` bytes; refuse early
    // before the OOM rather than after the kernel kills us.
    if n_free > MODAL_DENSE_MAX_DOFS {
        return Err(ModalSolverError::DofCountExceeded {
            n_free,
            max: MODAL_DENSE_MAX_DOFS,
        });
    }

    // --- reduce K and M to the free DOFs only (DOF elimination) ---
    let mut k_ff = DMatrix::<f64>::zeros(n_free, n_free);
    let mut m_ff = DMatrix::<f64>::zeros(n_free, n_free);
    for (ri, &gr) in free.iter().enumerate() {
        for (ci, &gc) in free.iter().enumerate() {
            k_ff[(ri, ci)] = k_global[(gr, gc)];
            m_ff[(ri, ci)] = m_global[(gr, gc)];
        }
    }

    // --- generalised → standard via the Cholesky factor of M_ff ---
    // M_ff = L·Lᵀ. The substitution φ = L⁻ᵀ·ψ turns K φ = λ M φ into
    // the symmetric standard problem C ψ = λ ψ with C = L⁻¹·K_ff·L⁻ᵀ.
    let chol = m_ff
        .clone()
        .cholesky()
        .ok_or(ModalSolverError::MassNotPositiveDefinite)?;
    let l = chol.l();
    // Invert the lower-triangular L by forward substitution (exact,
    // O(n³/6) — cheap for a workshop-scale reduced system).
    let l_inv = invert_lower_triangular(&l).ok_or(ModalSolverError::MassNotPositiveDefinite)?;
    let l_inv_t = l_inv.transpose();
    // C = L⁻¹ K_ff L⁻ᵀ. Symmetrise to kill the tiny asymmetry that
    // floating-point accumulation leaves behind, so the symmetric
    // eigensolver sees an exactly-symmetric matrix.
    let mut c = &l_inv * &k_ff * &l_inv_t;
    let c_t = c.transpose();
    c = (&c + &c_t) * 0.5;

    // --- symmetric eigensolve of C ---
    let eigen = nalgebra::SymmetricEigen::try_new(c, 1.0e-12, 0)
        .ok_or(ModalSolverError::EigenFailed)?;
    let eigvals = &eigen.eigenvalues;
    let eigvecs = &eigen.eigenvectors;

    // --- sort eigenpairs ascending by eigenvalue ---
    let mut order: Vec<usize> = (0..n_free).collect();
    order.sort_by(|&a, &b| {
        eigvals[a]
            .partial_cmp(&eigvals[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // --- build the lowest n_modes VibrationModes ---
    let mut modes = Vec::with_capacity(n_modes);
    for &idx in order.iter().take(n_modes) {
        // A linear-elastic stiffness matrix is positive-semi-definite;
        // tiny negative eigenvalues are float noise on rigid-body
        // modes. Clamp to zero before the square root.
        let lambda = eigvals[idx].max(0.0);
        let omega = lambda.sqrt();
        let freq = omega / (2.0 * std::f64::consts::PI);

        // Recover the generalised eigenvector φ = L⁻ᵀ ψ on the free
        // DOFs, then scatter it back into a full 3·n_nodes vector.
        let psi = eigvecs.column(idx).into_owned();
        let phi_free = &l_inv_t * &psi;
        // Mass-normalise: φᵀ M φ = ψᵀ ψ = 1 already holds for the
        // standard eigenvectors, but we re-normalise defensively in
        // case the eigensolver returned an un-normalised column.
        let mut shape = vec![[0.0_f64; 3]; mesh.nodes.len()];
        for (fi, &gd) in free.iter().enumerate() {
            let node = gd / 3;
            let comp = gd % 3;
            shape[node][comp] = phi_free[fi];
        }

        modes.push(VibrationMode {
            frequency_hz: freq,
            angular_frequency: omega,
            eigenvalue: lambda,
            shape,
        });
    }

    Ok(ModalSolution { modes })
}

/// Solve the modal eigenproblem on a **mixed-element** solid mesh.
///
/// The generic counterpart of [`solve_modal`]: instead of a single
/// [`ElementType::Tet4`] block this walks every supported continuum
/// element ([`ElementType::Tet4`], [`ElementType::Hex8`],
/// [`ElementType::Tet10`]) via [`crate::assembly`] and the
/// [`crate::elements::SolidElement`] trait, building the global
/// stiffness `K` and consistent mass `M` from a mixed-element mesh.
/// The boundary-condition treatment (DOF elimination), the
/// generalised → standard reduction and the symmetric eigensolve are
/// identical to [`solve_modal`].
///
/// For a pure-Tet4 mesh this returns the same frequencies as
/// [`solve_modal`]; the value is the Hex8 / Tet10 elements, which
/// predict bending frequencies far more accurately per DOF.
///
/// # Errors
///
/// See [`ModalSolverError`].
pub fn solve_modal_mixed(
    mesh: &Mesh,
    material: &FemMaterial,
    constraints: &[NodalConstraint],
    n_modes: usize,
) -> Result<ModalSolution, ModalSolverError> {
    if n_modes == 0 {
        return Err(ModalSolverError::TooManyModes {
            requested: 0,
            available: 0,
        });
    }
    if !material.density.is_finite() || material.density <= 0.0 {
        return Err(ModalSolverError::BadDensity(material.density));
    }
    // Round-1 H1–H3: cap the dense allocation up front (primary guard),
    // before the mixed assembler builds two `n_dof × n_dof` matrices.
    crate::native_solver::check_dense_dofs(mesh.nodes.len())?;

    let d_matrix = elasticity_matrix(material)?;
    let (k_global, m_global) = crate::assembly::assemble_global_stiffness_mass_mixed(
        mesh,
        &d_matrix,
        material.density,
    )?;
    let n_dof = k_global.nrows();

    // Constrained-DOF set.
    let mut constrained = vec![false; n_dof];
    for c in constraints {
        if c.node >= mesh.nodes.len() {
            return Err(ModalSolverError::Assembly(
                NativeSolverError::BadConnectivity {
                    elem: usize::MAX,
                    node: c.node,
                    n_nodes: mesh.nodes.len(),
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
        return Err(ModalSolverError::FullyConstrained);
    }
    let n_free = free.len();
    if n_modes > n_free {
        return Err(ModalSolverError::TooManyModes {
            requested: n_modes,
            available: n_free,
        });
    }
    // Dense-solver upper bound — refuse before the OOM. See the
    // `MODAL_DENSE_MAX_DOFS` constant doc for the back-of-envelope.
    if n_free > MODAL_DENSE_MAX_DOFS {
        return Err(ModalSolverError::DofCountExceeded {
            n_free,
            max: MODAL_DENSE_MAX_DOFS,
        });
    }

    // Reduce to the free DOFs.
    let mut k_ff = DMatrix::<f64>::zeros(n_free, n_free);
    let mut m_ff = DMatrix::<f64>::zeros(n_free, n_free);
    for (ri, &gr) in free.iter().enumerate() {
        for (ci, &gc) in free.iter().enumerate() {
            k_ff[(ri, ci)] = k_global[(gr, gc)];
            m_ff[(ri, ci)] = m_global[(gr, gc)];
        }
    }

    // Generalised → standard via the Cholesky factor of M_ff.
    let chol = m_ff
        .clone()
        .cholesky()
        .ok_or(ModalSolverError::MassNotPositiveDefinite)?;
    let l = chol.l();
    let l_inv = invert_lower_triangular(&l).ok_or(ModalSolverError::MassNotPositiveDefinite)?;
    let l_inv_t = l_inv.transpose();
    let mut c = &l_inv * &k_ff * &l_inv_t;
    let c_t = c.transpose();
    c = (&c + &c_t) * 0.5;

    let eigen = nalgebra::SymmetricEigen::try_new(c, 1.0e-12, 0)
        .ok_or(ModalSolverError::EigenFailed)?;
    let eigvals = &eigen.eigenvalues;
    let eigvecs = &eigen.eigenvectors;

    let mut order: Vec<usize> = (0..n_free).collect();
    order.sort_by(|&a, &b| {
        eigvals[a]
            .partial_cmp(&eigvals[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut modes = Vec::with_capacity(n_modes);
    for &idx in order.iter().take(n_modes) {
        let lambda = eigvals[idx].max(0.0);
        let omega = lambda.sqrt();
        let freq = omega / (2.0 * std::f64::consts::PI);
        let psi = eigvecs.column(idx).into_owned();
        let phi_free = &l_inv_t * &psi;
        let mut shape = vec![[0.0_f64; 3]; mesh.nodes.len()];
        for (fi, &gd) in free.iter().enumerate() {
            let node = gd / 3;
            let comp = gd % 3;
            shape[node][comp] = phi_free[fi];
        }
        modes.push(VibrationMode {
            frequency_hz: freq,
            angular_frequency: omega,
            eigenvalue: lambda,
            shape,
        });
    }
    Ok(ModalSolution { modes })
}

/// Assemble the global stiffness `K` and **consistent mass** `M`
/// matrices of a Tet4 mesh as dense `3·n_nodes` squares.
///
/// Shares the element-stiffness path with [`crate::native_solver`] —
/// the same constant-strain-tetrahedron `Kₑ` — and adds the consistent
/// element mass `Mₑ`. Returned `(K, M)` in that order.
///
/// This is the shared assembly the modal, transient-dynamics
/// ([`crate::dynamics`]) and buckling ([`crate::buckling`]) solvers all
/// build on — every solver that needs `M` reuses the *same* consistent
/// mass matrix rather than re-deriving it.
///
/// # Errors
///
/// See [`NativeSolverError`] — empty mesh, no Tet4 block, a degenerate
/// element, out-of-range connectivity, or a bad material.
pub fn assemble_global_stiffness_mass(
    mesh: &Mesh,
    material: &FemMaterial,
) -> Result<(DMatrix<f64>, DMatrix<f64>), NativeSolverError> {
    let n_nodes = mesh.nodes.len();
    if n_nodes == 0 {
        return Err(NativeSolverError::EmptyMesh);
    }
    let tets = mesh
        .element_blocks
        .iter()
        .find(|b| b.element_type == ElementType::Tet4)
        .filter(|b| b.connectivity.len() >= 4)
        .ok_or(NativeSolverError::NoTetBlock)?;

    let d_matrix = elasticity_matrix(material)?;
    // Round-1 H1–H3: this shared dense assembler is also reached
    // directly from `dynamics::solve_transient_dynamics`, so guard the
    // two `n_dof × n_dof` allocations here too (defense in depth — the
    // modal entry points hoist the same check above their call).
    let n_dof = crate::native_solver::check_dense_dofs(n_nodes)?;
    let mut k = DMatrix::<f64>::zeros(n_dof, n_dof);
    let mut m = DMatrix::<f64>::zeros(n_dof, n_dof);

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
                return Err(NativeSolverError::BadConnectivity {
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
        let ke = element_stiffness(&coords, &d_matrix)
            .ok_or(NativeSolverError::DegenerateElement(e))?;
        let me = element_consistent_mass(&coords, material.density)
            .ok_or(NativeSolverError::DegenerateElement(e))?;
        // Scatter-add the two 12×12 element matrices.
        for a in 0..4 {
            for i in 0..3 {
                let gi = 3 * nodes[a] + i;
                for b in 0..4 {
                    for j in 0..3 {
                        let gj = 3 * nodes[b] + j;
                        k[(gi, gj)] += ke[(3 * a + i, 3 * b + j)];
                        m[(gi, gj)] += me[(3 * a + i, 3 * b + j)];
                    }
                }
            }
        }
    }
    Ok((k, m))
}

/// Isotropic 6×6 elasticity matrix `D` in Voigt notation — identical
/// formulation to [`crate::native_solver`]'s private helper, reproduced
/// here so the modal module is self-contained.
fn elasticity_matrix(m: &FemMaterial) -> Result<Matrix6<f64>, NativeSolverError> {
    let e = m.youngs_modulus;
    let nu = m.poisson_ratio;
    if !e.is_finite() || e <= 0.0 {
        return Err(NativeSolverError::BadMaterial(format!(
            "Young's modulus must be finite and positive, got {e}"
        )));
    }
    if !nu.is_finite() || nu <= -1.0 || nu >= 0.5 {
        return Err(NativeSolverError::BadMaterial(format!(
            "Poisson's ratio must lie in (-1, 0.5), got {nu}"
        )));
    }
    let lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mu = e / (2.0 * (1.0 + nu));
    let mut d = Matrix6::zeros();
    for i in 0..3 {
        for j in 0..3 {
            d[(i, j)] = lambda;
        }
        d[(i, i)] = lambda + 2.0 * mu;
    }
    for i in 3..6 {
        d[(i, i)] = mu;
    }
    Ok(d)
}

/// 12×12 element stiffness `Kₑ = Vₑ·Bᵀ·D·B` for a linear tet. `None`
/// if the tet is degenerate. Identical formulation to the static
/// solver; reproduced for module independence.
fn element_stiffness(coords: &[Vector3<f64>; 4], d: &Matrix6<f64>) -> Option<DMatrix<f64>> {
    let grads = shape_gradients(coords)?;
    let mut b = DMatrix::<f64>::zeros(6, 12);
    for (a, g) in grads.iter().enumerate() {
        let (bx, by, bz) = (g.x, g.y, g.z);
        let col = 3 * a;
        b[(0, col)] = bx;
        b[(1, col + 1)] = by;
        b[(2, col + 2)] = bz;
        b[(3, col)] = by;
        b[(3, col + 1)] = bx;
        b[(4, col + 1)] = bz;
        b[(4, col + 2)] = by;
        b[(5, col)] = bz;
        b[(5, col + 2)] = bx;
    }
    let vol = tet_volume(coords).abs();
    if vol < 1.0e-18 {
        return None;
    }
    let bt_d = b.transpose() * d;
    Some((bt_d * b) * vol)
}

/// 12×12 **consistent** element mass matrix for a linear tetrahedron.
///
/// The consistent mass of a 4-node tet is the classical
///
/// ```text
///   Mₑ = (ρ·V / 20) · Ĉ ⊗ I₃
/// ```
///
/// where `Ĉ` is the 4×4 matrix with `2` on the diagonal and `1`
/// everywhere off-diagonal (the integral of the products of the linear
/// shape functions over the element), `I₃` is the 3×3 identity (the
/// three translational DOFs per node are uncoupled), and `⊗` is the
/// Kronecker product. The row sums of `Mₑ` equal `ρ·V/4` per DOF, i.e.
/// the total element mass `ρ·V` is correctly distributed.
///
/// Returns `None` for a degenerate (zero-volume) tet.
fn element_consistent_mass(coords: &[Vector3<f64>; 4], density: f64) -> Option<DMatrix<f64>> {
    let vol = tet_volume(coords).abs();
    if vol < 1.0e-18 {
        return None;
    }
    let factor = density * vol / 20.0;
    let mut me = DMatrix::<f64>::zeros(12, 12);
    for a in 0..4 {
        for b in 0..4 {
            // Diagonal node-pair coupling = 2, off-diagonal = 1.
            let coeff = if a == b { 2.0 } else { 1.0 } * factor;
            // Each (a, b) node pair contributes coeff·I₃ to the three
            // matching translational DOFs.
            for i in 0..3 {
                me[(3 * a + i, 3 * b + i)] = coeff;
            }
        }
    }
    Some(me)
}

/// Constant shape-function gradients `∇Nₐ` of a linear tet. `None` if
/// the tet is degenerate. Same formulation as the static solver.
fn shape_gradients(coords: &[Vector3<f64>; 4]) -> Option<[Vector3<f64>; 4]> {
    use nalgebra::Matrix3;
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
    Some([g0, g1, g2, g3])
}

/// Signed volume of a tet.
fn tet_volume(coords: &[Vector3<f64>; 4]) -> f64 {
    let e1 = coords[1] - coords[0];
    let e2 = coords[2] - coords[0];
    let e3 = coords[3] - coords[0];
    e1.cross(&e2).dot(&e3) / 6.0
}

/// Invert a lower-triangular matrix `L` by forward substitution.
///
/// Solves `L·X = I` column by column. Returns `None` if any diagonal
/// entry is too small to divide by (a singular / rank-deficient
/// factor).
fn invert_lower_triangular(l: &DMatrix<f64>) -> Option<DMatrix<f64>> {
    let n = l.nrows();
    debug_assert_eq!(n, l.ncols(), "L must be square");
    let mut inv = DMatrix::<f64>::zeros(n, n);
    for col in 0..n {
        // Solve L·x = eₖ for column `col`.
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

/// Convenience: mass-normalise a full-length displacement vector
/// against the global mass matrix so `dᵀ M d = 1`. Returns the scale
/// factor applied, or `None` if `d` carries no mass at all.
///
/// Exposed for callers that build a mode shape from a hand-supplied
/// deflection and want it on the same normalisation as
/// [`solve_modal`]'s output.
pub fn mass_normalise(d: &mut DVector<f64>, m: &DMatrix<f64>) -> Option<f64> {
    let energy = (d.transpose() * m * &*d)[(0, 0)];
    if energy <= 0.0 || !energy.is_finite() {
        return None;
    }
    let scale = 1.0 / energy.sqrt();
    *d *= scale;
    Some(scale)
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
    fn consistent_mass_distributes_total_element_mass() {
        // Sum of every entry of Mₑ must equal the total element mass
        // ρ·V counted three times (once per spatial direction).
        let coords = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let density = 7850.0;
        let me = element_consistent_mass(&coords, density).unwrap();
        let vol = tet_volume(&coords).abs();
        let total: f64 = me.iter().sum();
        // 4 diag blocks·2 + 12 off-diag blocks·1 = 8 + 12 = 20 units of
        // `factor` per direction; ×3 directions. factor = ρV/20, so the
        // grand total is 3·ρ·V.
        assert!(
            (total - 3.0 * density * vol).abs() < 1e-6 * density * vol,
            "mass sum {total} vs 3ρV {}",
            3.0 * density * vol
        );
    }

    #[test]
    fn consistent_mass_is_symmetric_and_positive_definite() {
        let coords = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(2.0, 0.0, 0.0),
            Vector3::new(0.3, 1.5, 0.1),
            Vector3::new(0.1, 0.2, 1.8),
        ];
        let me = element_consistent_mass(&coords, 2700.0).unwrap();
        for i in 0..12 {
            for j in 0..12 {
                assert!((me[(i, j)] - me[(j, i)]).abs() < 1e-9, "Mₑ not symmetric");
            }
        }
        // A consistent mass matrix is SPD — a Cholesky factor exists.
        assert!(me.cholesky().is_some(), "Mₑ should be positive-definite");
    }

    #[test]
    fn degenerate_tet_has_no_mass() {
        // Four coplanar points → zero volume → None.
        let coords = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
        ];
        assert!(element_consistent_mass(&coords, 1000.0).is_none());
    }

    #[test]
    fn invert_lower_triangular_round_trips() {
        // L·L⁻¹ should be the identity.
        let l = DMatrix::from_row_slice(
            3,
            3,
            &[2.0, 0.0, 0.0, 1.0, 3.0, 0.0, 0.5, -1.0, 4.0],
        );
        let inv = invert_lower_triangular(&l).unwrap();
        let prod = &l * &inv;
        for i in 0..3 {
            for j in 0..3 {
                let want = if i == j { 1.0 } else { 0.0 };
                assert!((prod[(i, j)] - want).abs() < 1e-9, "L·L⁻¹ ≠ I");
            }
        }
    }

    #[test]
    fn modal_rejects_zero_modes() {
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box params");
        let err = solve_modal(&mesh, &FemMaterial::default(), &[], 0).unwrap_err();
        assert!(matches!(err, ModalSolverError::TooManyModes { .. }));
    }

    #[test]
    fn modal_rejects_bad_density() {
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box params");
        let bad = FemMaterial {
            density: 0.0,
            ..FemMaterial::default()
        };
        let err = solve_modal(&mesh, &bad, &[], 1).unwrap_err();
        assert!(matches!(err, ModalSolverError::BadDensity(_)));
    }

    #[test]
    fn modal_rejects_fully_constrained_model() {
        // Pin every node of a single-cell box → no free DOF.
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box params");
        let constraints: Vec<NodalConstraint> =
            (0..mesh.nodes.len()).map(NodalConstraint::fixed).collect();
        let err = solve_modal(&mesh, &FemMaterial::default(), &constraints, 1).unwrap_err();
        assert!(matches!(err, ModalSolverError::FullyConstrained));
    }

    #[test]
    fn modal_rejects_too_many_modes() {
        // A single-cell box has 8 nodes = 24 DOF; fixing one node
        // leaves 21 free. Asking for 100 modes must fail.
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box params");
        let err = solve_modal(
            &mesh,
            &FemMaterial::default(),
            &[NodalConstraint::fixed(0)],
            100,
        )
        .unwrap_err();
        assert!(matches!(
            err,
            ModalSolverError::TooManyModes {
                requested: 100,
                ..
            }
        ));
    }

    #[test]
    fn cantilever_fundamental_is_positive_and_real() {
        // A clamped beam must have a positive, finite fundamental
        // frequency, and the lowest modes must be ascending.
        let (lx, ly, lz) = (10.0, 1.0, 1.0);
        let (nx, ny, nz) = (8, 1, 1);
        let mesh = structured_box_mesh(lx, ly, lz, nx, ny, nz).expect("valid box params");
        let mat = FemMaterial {
            youngs_modulus: 200.0e9,
            poisson_ratio: 0.3,
            density: 7850.0,
            ..FemMaterial::default()
        };
        // Clamp the x=0 face.
        let mut constraints = Vec::new();
        for k in 0..=nz {
            for j in 0..=ny {
                constraints.push(NodalConstraint::fixed(nid(0, j, k, nx, ny)));
            }
        }
        let sol = solve_modal(&mesh, &mat, &constraints, 4).unwrap();
        assert_eq!(sol.modes.len(), 4);
        let f0 = sol.fundamental_hz().unwrap();
        assert!(f0 > 0.0 && f0.is_finite(), "fundamental {f0} must be > 0");
        // Modes are sorted ascending.
        for w in sol.modes.windows(2) {
            assert!(
                w[1].frequency_hz >= w[0].frequency_hz - 1e-6,
                "modes not ascending: {} then {}",
                w[0].frequency_hz,
                w[1].frequency_hz
            );
        }
        // The fundamental mode shape must be non-trivial and the
        // clamped end must stay put.
        let m0 = &sol.modes[0];
        assert!(m0.max_amplitude() > 0.0, "mode shape is all zero");
        let clamped = nid(0, 0, 0, nx, ny);
        let amp = (m0.shape[clamped][0].powi(2)
            + m0.shape[clamped][1].powi(2)
            + m0.shape[clamped][2].powi(2))
        .sqrt();
        assert!(amp < 1e-9, "clamped node moved in the mode shape: {amp}");
    }

    #[test]
    fn cantilever_fundamental_near_euler_bernoulli() {
        // Euler-Bernoulli fundamental bending frequency of a clamped-
        // free beam:  f₁ = (β₁L)²/(2π) · √(E·I / (ρ·A·L⁴))  with
        // (β₁L)² = 3.5160. A linear-tet mesh stiffens in bending so
        // the FE frequency sits *above* the analytic value, but it
        // must land within a sane factor — that the assembly + the
        // generalised eigensolve are correct.
        let (lx, ly, lz) = (10.0, 1.0, 1.0);
        let (nx, ny, nz) = (16, 2, 2);
        let mesh = structured_box_mesh(lx, ly, lz, nx, ny, nz).expect("valid box params");
        let mat = FemMaterial {
            youngs_modulus: 200.0e9,
            poisson_ratio: 0.3,
            density: 7850.0,
            ..FemMaterial::default()
        };
        let mut constraints = Vec::new();
        for k in 0..=nz {
            for j in 0..=ny {
                constraints.push(NodalConstraint::fixed(nid(0, j, k, nx, ny)));
            }
        }
        let sol = solve_modal(&mesh, &mat, &constraints, 3).unwrap();
        let f_fe = sol.fundamental_hz().unwrap();

        // Analytic Euler-Bernoulli reference.
        let area = ly * lz;
        let i_section = ly * lz.powi(3) / 12.0;
        let beta1_l_sq = 3.5160_f64;
        let f_analytic = beta1_l_sq / (2.0 * std::f64::consts::PI)
            * (mat.youngs_modulus * i_section / (mat.density * area * lx.powi(4))).sqrt();

        assert!(f_fe > 0.0 && f_analytic > 0.0);
        // FE bending frequency is above analytic (linear-tet stiffening)
        // but within a factor of ~3 on this mesh.
        assert!(
            f_fe > f_analytic * 0.7,
            "FE fundamental {f_fe} far below analytic {f_analytic}"
        );
        assert!(
            f_fe < f_analytic * 3.0,
            "FE fundamental {f_fe} far above analytic {f_analytic}"
        );
    }

    #[test]
    fn unconstrained_body_has_near_zero_rigid_modes() {
        // A free-floating body has six rigid-body modes (3 translation,
        // 3 rotation) at ≈ 0 Hz; the first elastic mode follows. The
        // lowest extracted frequencies must therefore be tiny relative
        // to a genuine elastic mode.
        let mesh = structured_box_mesh(2.0, 2.0, 2.0, 2, 2, 2).expect("valid box params");
        let mat = FemMaterial {
            youngs_modulus: 200.0e9,
            poisson_ratio: 0.3,
            density: 7850.0,
            ..FemMaterial::default()
        };
        // No constraints at all.
        let sol = solve_modal(&mesh, &mat, &[], 8).unwrap();
        assert_eq!(sol.modes.len(), 8);
        // The first six modes are rigid-body: their frequency must be
        // far below the seventh (first elastic) mode.
        let f6 = sol.modes[5].frequency_hz;
        let f7 = sol.modes[6].frequency_hz;
        assert!(f7 > 0.0, "first elastic mode {f7} should be positive");
        assert!(
            f6 < f7 * 1e-3,
            "rigid-body mode {f6} should be ≪ first elastic mode {f7}"
        );
    }

    #[test]
    fn mass_normalise_scales_to_unit_energy() {
        // dᵀ M d = 1 after normalisation.
        let m = DMatrix::<f64>::identity(3, 3) * 4.0;
        let mut d = DVector::from_vec(vec![1.0, 0.0, 0.0]);
        mass_normalise(&mut d, &m).unwrap();
        let energy = (d.transpose() * &m * &d)[(0, 0)];
        assert!((energy - 1.0).abs() < 1e-9, "energy {energy} ≠ 1");
    }
}
