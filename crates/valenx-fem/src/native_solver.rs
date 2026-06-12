//! Native in-process linear-static FEA solver (Phase 24.5).
//!
//! ## What this is
//!
//! A genuine, self-contained **linear elastostatic finite-element
//! solver** for 4-node linear tetrahedra (`C3D4` / [`ElementType::Tet4`]).
//! It assembles the global stiffness matrix, applies displacement and
//! force boundary conditions, solves the sparse symmetric-positive-
//! definite system with a Cholesky factorisation
//! ([`nalgebra_sparse::factorization::CscCholesky`]), and recovers the
//! per-node displacement field plus per-node von Mises stress.
//!
//! This is the in-process alternative to dispatching to the CalculiX /
//! Elmer / Code_Aster subprocess adapters ([`crate::FemSolverChoice`]):
//! no external binary, no input-deck round-trip, runs anywhere the
//! crate compiles.
//!
//! ## Honest scope — what it is and is not
//!
//! **It is a real solver.** The element formulation is the textbook
//! constant-strain tetrahedron (CST-3D): the shape-function gradients
//! are exact for a linear tet, the element stiffness
//! `Kₑ = Vₑ · Bᵀ · D · B` is closed-form (no numerical quadrature
//! needed — the integrand is constant over the element), the assembly
//! is the standard scatter-add, and the linear solve is an exact
//! direct factorisation. For a linear-elastic small-deformation
//! problem meshed with linear tets, the displacements it returns are
//! the correct FE solution for that mesh.
//!
//! **It is not CalculiX.** It deliberately ships only the 90 % case:
//!
//! - **Linear static only.** No modal / buckling / transient /
//!   nonlinear / contact / thermal coupling. Those are whole solvers.
//! - **Tet4 only.** Linear tets stiffen in bending ("shear locking"
//!   is mild but present); a quadratic Tet10 element would be far more
//!   accurate per DOF. Tet10 is a bounded follow-up (same assembly,
//!   10-node shape functions, 4-point quadrature).
//! - **It needs a volume mesh.** A linear-static solve runs on a
//!   tetrahedral *volume* mesh. Tetrahedralising an arbitrary B-rep
//!   solid (constrained Delaunay) is itself a major project and stays
//!   an external-adapter dependency (gmsh / Netgen). To make the
//!   solver testable and usable end-to-end without a mesher this
//!   module ships [`structured_box_mesh`], a structured 6-tets-per-
//!   cell tetrahedralisation of an axis-aligned box.
//! - **Isotropic linear elasticity only** — two constants (E, ν).
//!
//! Within that scope the numbers are real: a cantilever beam under an
//! end load deflects in the correct direction by an amount that
//! converges toward Euler-Bernoulli theory as the mesh refines, and a
//! bar in uniaxial tension reproduces `σ = E·ε` to solver precision.

use nalgebra::{DMatrix, DVector, Matrix3, Matrix6, Vector3};
use nalgebra_sparse::factorization::CscCholesky;
use nalgebra_sparse::{CooMatrix, CscMatrix};
use thiserror::Error;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::material::FemMaterial;

/// Errors from the native FEA solver.
#[derive(Debug, Error)]
pub enum NativeSolverError {
    /// The supplied mesh has no [`ElementType::Tet4`] block, so there
    /// is nothing volumetric to solve.
    #[error("mesh has no Tet4 element block — the native solver needs a tetrahedral volume mesh")]
    NoTetBlock,
    /// The mesh has zero nodes.
    #[error("mesh has no nodes")]
    EmptyMesh,
    /// A tetrahedron is degenerate (zero or negative signed volume) —
    /// the mesh is invalid or inside-out.
    #[error("tetrahedron {0} is degenerate (|volume| < tolerance)")]
    DegenerateElement(usize),
    /// A connectivity index points past the end of the node array.
    #[error("element {elem} references node {node} but the mesh has only {n_nodes} nodes")]
    BadConnectivity {
        /// 0-based element index.
        elem: usize,
        /// Out-of-range node index.
        node: usize,
        /// Node-array length.
        n_nodes: usize,
    },
    /// No displacement boundary condition was supplied, so the
    /// stiffness matrix is singular (the body can float / rotate
    /// freely) and the solve has no unique answer.
    #[error("no displacement constraint — the model is unrestrained (rigid-body singular)")]
    Unconstrained,
    /// A material constant is non-physical (E ≤ 0, or ν outside the
    /// open interval where the elasticity tensor stays positive-definite).
    #[error("invalid material: {0}")]
    BadMaterial(String),
    /// The Cholesky factorisation failed — the assembled system was
    /// not positive-definite (typically an under-constrained model
    /// the cheap singularity check did not catch).
    #[error("linear solve failed: stiffness matrix is not positive-definite")]
    SolveFailed,
    /// A mesh-builder or input parameter is invalid (zero subdivisions,
    /// negative dimensions, overflow). Round-3 fix: replaces the old
    /// `assert!()` calls in [`structured_box_mesh`] with a structured
    /// error so callers can recover gracefully.
    #[error("invalid params: {reason}")]
    InvalidParams {
        /// Human-readable explanation.
        reason: String,
    },
    /// The model is too large for the **dense** assembly path. A dense
    /// `n_dof × n_dof` `f64` matrix needs `8·n_dof²` bytes (and the
    /// modal / buckling / dynamics paths build two of them), so a large
    /// mesh would OOM the host. Round-1 fix: [`check_dense_dofs`] rejects
    /// it up front — before any [`nalgebra::DMatrix::zeros`] allocation —
    /// rather than letting the kernel kill the process. The sparse
    /// static / thermal / contact / nonlinear paths are unaffected.
    #[error(
        "dense solve needs {dofs} DOFs, the dense path supports at most {max}; \
         coarsen the mesh or use a sparse solver"
    )]
    TooLarge {
        /// The model's degree-of-freedom count `3·n_nodes` (or the
        /// overflow sentinel `usize::MAX` if `3·n_nodes` itself
        /// overflowed).
        dofs: usize,
        /// The dense-path upper bound, [`MAX_DENSE_DOFS`].
        max: usize,
    },
    /// A load or boundary-condition input carried a non-finite value
    /// (`NaN` or `±∞`). Round-1 fix: materials and node coordinates were
    /// finiteness-checked, but nodal forces / prescribed displacements /
    /// initial state were pushed straight into the right-hand side. The
    /// Cholesky back-substitution is pure arithmetic, so a non-finite
    /// RHS yields a non-finite displacement returned as `Ok(..)` — a
    /// single `force: [NaN, 0, 0]` silently corrupts the whole solution.
    /// Validating the input up front lets the error name the cause.
    #[error("non-finite {kind} at node {node} (NaN or infinity is not a valid load/BC)")]
    InvalidLoad {
        /// 0-based node index carrying the bad value.
        node: usize,
        /// Which input was non-finite — e.g. `"force"`, `"prescribed
        /// displacement"`, `"initial displacement"`, `"initial
        /// velocity"`.
        kind: &'static str,
    },
}

/// Reject a non-finite (`NaN` / `±∞`) component in a 3-vector load or
/// boundary-condition input.
///
/// Returns [`NativeSolverError::InvalidLoad`] naming `node` and `kind`
/// at the first non-finite component, else `Ok(())`. Used to validate
/// nodal forces, initial displacement / velocity, etc. *before* they
/// reach the right-hand side — a non-finite RHS would otherwise produce
/// a silently-non-finite solution from the (pure-substitution) Cholesky
/// solve.
pub(crate) fn ensure_finite3(
    values: &[f64; 3],
    node: usize,
    kind: &'static str,
) -> Result<(), NativeSolverError> {
    if values.iter().any(|v| !v.is_finite()) {
        return Err(NativeSolverError::InvalidLoad { node, kind });
    }
    Ok(())
}

/// Upper bound on the **global** degree-of-freedom count `3·n_nodes`
/// the *dense* assembly path accepts.
///
/// The dense buckling / modal / dynamics solvers — and the dense
/// [`crate::assembly::assemble_global_stiffness`] /
/// [`crate::assembly::assemble_global_stiffness_mass_mixed`] entry
/// points — allocate one or two `n_dof × n_dof` `f64` matrices. At
/// `8·n_dof²` bytes each, a 16 000-DOF matrix is ≈ 2 GB, and the
/// two-matrix paths (K + M, or K + K_g) reach ≈ 4 GB at the cap — the
/// outer edge of what a workstation tolerates. Beyond it the helper
/// [`check_dense_dofs`] returns [`NativeSolverError::TooLarge`] *before*
/// any allocation rather than letting the host OOM.
///
/// This caps the **full** global system; the modal/buckling solvers
/// additionally cap the post-elimination *free*-DOF count at the much
/// smaller [`crate::modal_solver::MODAL_DENSE_MAX_DOFS`] (the size of
/// the matrix actually handed to the dense eigensolver). The two checks
/// are complementary: this one bounds the allocation, that one bounds
/// the `O(n³)` eigensolve.
///
/// The sparse static / thermal / contact / nonlinear paths build no
/// dense global matrix and are not subject to this cap.
pub const MAX_DENSE_DOFS: usize = 16_000;

/// Reject an **already-computed** dense DOF count too large for the
/// dense assembly path, **before** any `n_dof × n_dof` allocation.
///
/// Returns `Ok(n_dof)` when `n_dof ≤ `[`MAX_DENSE_DOFS`] *and* `n_dof²`
/// fits in `usize`; otherwise [`NativeSolverError::TooLarge`]. This is
/// the single capping path: callers that know their per-node DOF
/// multiplier compute `n_dof` themselves (via [`usize::checked_mul`])
/// and pass it here. The 3-DOF/node volumetric path uses the
/// [`check_dense_dofs`] convenience wrapper; the 6-DOF/node beam path
/// ([`crate::beam`]) computes `6·n_nodes` and routes the result through
/// this helper — so a beam mesh is no longer under-counted 2× by the
/// hardcoded `3·n_nodes` wrapper.
///
/// The check is pure arithmetic — `O(1)`, no allocation — so it is safe
/// to call at the very top of a solver on an untrusted DOF count.
pub fn check_dense_dof_count(n_dof: usize) -> Result<usize, NativeSolverError> {
    if n_dof > MAX_DENSE_DOFS {
        return Err(NativeSolverError::TooLarge {
            dofs: n_dof,
            max: MAX_DENSE_DOFS,
        });
    }
    // The dense matrix is `n_dof²` elements; even under the DOF cap,
    // guard the element-count multiplication so the contract ("n_dof²
    // is representable") holds for any future, larger MAX_DENSE_DOFS.
    if n_dof.checked_mul(n_dof).is_none() {
        return Err(NativeSolverError::TooLarge {
            dofs: n_dof,
            max: MAX_DENSE_DOFS,
        });
    }
    Ok(n_dof)
}

/// Compute the dense global DOF count `n_dof = 3·n_nodes` and reject a
/// model too large for the dense assembly path **before** any
/// `n_dof × n_dof` allocation is attempted.
///
/// Returns `Ok(n_dof)` when `n_dof ≤ `[`MAX_DENSE_DOFS`] *and* `n_dof²`
/// fits in `usize`; otherwise [`NativeSolverError::TooLarge`]. Every
/// arithmetic step uses [`usize::checked_mul`] so a pathological
/// `n_nodes` surfaces as a structured error rather than a silent wrap
/// (this also closes the Low-1 `3·n_nodes` overflow concern). The
/// per-element count and `n_dof²` bound are delegated to the single
/// capping path [`check_dense_dof_count`]. The check is pure arithmetic
/// — `O(1)`, no allocation — so it is safe to call at the very top of a
/// solver on an untrusted node count.
pub fn check_dense_dofs(n_nodes: usize) -> Result<usize, NativeSolverError> {
    // `3·n_nodes` — the global DOF count for the 3-DOF/node volumetric
    // path. A wrap here would understate the size and defeat the cap,
    // so surface the overflow.
    let n_dof = n_nodes.checked_mul(3).ok_or(NativeSolverError::TooLarge {
        dofs: usize::MAX,
        max: MAX_DENSE_DOFS,
    })?;
    check_dense_dof_count(n_dof)
}

/// A single-node Dirichlet (displacement) boundary condition for the
/// native solver. Pins one, two, or three translational DOFs of one
/// node to prescribed values.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct NodalConstraint {
    /// 0-based node index.
    pub node: usize,
    /// Per-axis pin: `Some(v)` fixes that DOF to displacement `v`
    /// (metres), `None` leaves it free. Order is `[x, y, z]`.
    pub fixed: [Option<f64>; 3],
}

impl NodalConstraint {
    /// Fully fix a node (encastré) — all three translations zero.
    pub fn fixed(node: usize) -> Self {
        Self {
            node,
            fixed: [Some(0.0), Some(0.0), Some(0.0)],
        }
    }

    /// Pin a node to a prescribed displacement vector.
    pub fn displaced(node: usize, value: [f64; 3]) -> Self {
        Self {
            node,
            fixed: [Some(value[0]), Some(value[1]), Some(value[2])],
        }
    }
}

/// A single-node concentrated force for the native solver.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct NodalForce {
    /// 0-based node index.
    pub node: usize,
    /// Force vector `[Fx, Fy, Fz]` in newtons.
    pub force: [f64; 3],
}

/// Result of a native linear-static solve.
#[derive(Clone, Debug)]
pub struct NativeSolution {
    /// Per-node displacement `[ux, uy, uz]` in metres. Indexed by node.
    pub displacement: Vec<[f64; 3]>,
    /// Per-node von Mises stress in pascals. Computed per element
    /// (constant over a linear tet) then averaged to the nodes.
    pub von_mises: Vec<f64>,
    /// Per-node full Cauchy stress tensor in Voigt order
    /// `[σxx, σyy, σzz, σxy, σyz, σzx]`, averaged to the nodes.
    pub stress: Vec<[f64; 6]>,
}

impl NativeSolution {
    /// Largest displacement magnitude across all nodes.
    pub fn max_displacement(&self) -> f64 {
        self.displacement
            .iter()
            .map(|d| (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt())
            .fold(0.0, f64::max)
    }

    /// Largest nodal von Mises stress.
    pub fn max_von_mises(&self) -> f64 {
        self.von_mises.iter().copied().fold(0.0, f64::max)
    }

    /// Mean nodal von Mises stress over the field — the typical stress level, as
    /// opposed to the peak that [`NativeSolution::max_von_mises`] reveals. Their
    /// ratio (peak / mean) is the stress-concentration factor `Kt`: ~1 for a
    /// uniform field, larger where stress crowds into a notch or fixed root.
    /// Returns `0` for an empty field.
    pub fn mean_von_mises(&self) -> f64 {
        if self.von_mises.is_empty() {
            0.0
        } else {
            self.von_mises.iter().sum::<f64>() / self.von_mises.len() as f64
        }
    }

    /// Largest principal (maximum normal) stress over all nodes (Pa) — the
    /// Rankine / maximum-normal-stress measure for brittle failure, the tensile
    /// counterpart to the von Mises (ductile) stress. Returns
    /// `f64::NEG_INFINITY` for an empty field.
    pub fn max_principal_stress(&self) -> f64 {
        self.stress
            .iter()
            .map(|&s| max_principal_stress_voigt(s))
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Smallest principal (most compressive) stress over all nodes (Pa) — the
    /// compressive counterpart to [`NativeSolution::max_principal_stress`].
    /// Returns `f64::INFINITY` for an empty field.
    pub fn min_principal_stress(&self) -> f64 {
        self.stress
            .iter()
            .map(|&s| min_principal_stress_voigt(s))
            .fold(f64::INFINITY, f64::min)
    }

    /// Maximum shear stress over all nodes (Pa) — `(σ₁ − σ₃)/2` evaluated per
    /// node, the Tresca / maximum-shear-stress measure (the conservative
    /// ductile-yield criterion alongside von Mises). Computed per node so the
    /// peak is a true point value, not the global principal range. Returns `0`
    /// for an empty field.
    pub fn max_shear_stress(&self) -> f64 {
        self.stress
            .iter()
            .map(|&s| {
                let (lo, hi) = principal_extremes_voigt(s);
                0.5 * (hi - lo)
            })
            .fold(0.0, f64::max)
    }

    /// The maximum (most tensile) hydrostatic / mean stress over all nodes (Pa)
    /// — `σ_m = (σxx + σyy + σzz)/3`, the volumetric (pressure-like) part of the
    /// stress, complementary to the deviatoric von Mises measure. A high tensile
    /// hydrostatic stress promotes void growth and cavitation. Returns
    /// `f64::NEG_INFINITY` for an empty field.
    pub fn max_hydrostatic_stress(&self) -> f64 {
        self.stress
            .iter()
            .map(|s| (s[0] + s[1] + s[2]) / 3.0)
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// The minimum (most compressive) hydrostatic / mean stress over all nodes
    /// (Pa) — the compressive counterpart to
    /// [`NativeSolution::max_hydrostatic_stress`]. Returns `f64::INFINITY` for an
    /// empty field.
    pub fn min_hydrostatic_stress(&self) -> f64 {
        self.stress
            .iter()
            .map(|s| (s[0] + s[1] + s[2]) / 3.0)
            .fold(f64::INFINITY, f64::min)
    }

    /// Stress triaxiality `T = σ_m / σ_vm` (mean stress over von Mises stress)
    /// at the node of peak von Mises stress — the site where ductile failure
    /// initiates. Triaxiality sets the ductile-fracture mode: `1/3` in uniaxial
    /// tension, `2/3` in equibiaxial tension, larger under hydrostatic tension
    /// (void growth, brittle-like), negative under net compression. Returns `0`
    /// for an empty solution or an unstressed body (zero peak von Mises).
    pub fn peak_triaxiality(&self) -> f64 {
        let Some((idx, &peak)) = self
            .von_mises
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
        else {
            return 0.0;
        };
        if peak <= 1e-12 {
            return 0.0;
        }
        match self.stress.get(idx) {
            // Mean (hydrostatic) stress = trace/3, divided by the von Mises peak.
            Some(s) => (s[0] + s[1] + s[2]) / 3.0 / peak,
            None => 0.0,
        }
    }

    /// The **Lode parameter** `μ_L = (2σ₂ − σ₁ − σ₃)/(σ₁ − σ₃)` at the node of
    /// peak von Mises stress — the *deviatoric* descriptor of the stress state,
    /// the companion to [`NativeSolution::peak_triaxiality`]. Where triaxiality
    /// fixes the hydrostatic axis, the Lode parameter places the state on the
    /// deviatoric (π-)plane: `−1` in axisymmetric / uniaxial **tension**, `0` in
    /// pure shear or plane strain, `+1` in axisymmetric / uniaxial
    /// **compression**. Modern ductile-fracture loci depend on *both* the
    /// triaxiality and the Lode parameter, so together they characterise the
    /// failure-driving state at the critical point. Returns `0` for an empty
    /// solution, an unstressed body, or the degenerate hydrostatic state
    /// `σ₁ = σ₃` where the parameter is undefined.
    pub fn peak_lode_parameter(&self) -> f64 {
        let Some((idx, &peak)) = self
            .von_mises
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
        else {
            return 0.0;
        };
        if peak <= 1e-12 {
            return 0.0;
        }
        match self.stress.get(idx) {
            Some(&s) => {
                let [s1, s2, s3] = principal_triple_voigt(s);
                let range = s1 - s3;
                if range <= 1e-12 {
                    0.0 // σ₁ = σ₃ (hydrostatic) → undefined
                } else {
                    (2.0 * s2 - s1 - s3) / range
                }
            }
            None => 0.0,
        }
    }

    /// The index of the node carrying the peak von Mises stress — the critical
    /// location where ductile failure initiates (for a tip-loaded cantilever,
    /// the fixed root). Pairs with the mesh node list to recover its
    /// coordinates. `None` for an empty field.
    pub fn peak_von_mises_index(&self) -> Option<usize> {
        self.von_mises
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| i)
    }
}

/// The (minimum, maximum) principal stresses (eigenvalues) of a symmetric
/// Cauchy stress tensor in Voigt order `[σxx, σyy, σzz, σxy, σyz, σzx]`, via the
/// closed-form trigonometric solution for a symmetric 3×3 matrix (Smith 1961).
fn principal_extremes_voigt(s: [f64; 6]) -> (f64, f64) {
    let (sxx, syy, szz, sxy, syz, szx) = (s[0], s[1], s[2], s[3], s[4], s[5]);
    let p1 = sxy * sxy + syz * syz + szx * szx;
    if p1 <= 0.0 {
        // Already diagonal: the principals are the diagonal entries.
        let lo = sxx.min(syy).min(szz);
        let hi = sxx.max(syy).max(szz);
        return (lo, hi);
    }
    let q = (sxx + syy + szz) / 3.0;
    let a = sxx - q;
    let b = syy - q;
    let c = szz - q;
    let p2 = a * a + b * b + c * c + 2.0 * p1;
    let p = (p2 / 6.0).sqrt();
    // det(A − qI), then r = det((A − qI)/p) / 2 clamped to the acos domain.
    let det = a * (b * c - syz * syz) - sxy * (sxy * c - syz * szx) + szx * (sxy * syz - b * szx);
    let r = (det / (p * p * p) / 2.0).clamp(-1.0, 1.0);
    let phi = r.acos() / 3.0;
    let two_thirds_pi = 2.0 * std::f64::consts::PI / 3.0;
    // The cos(φ) branch is the largest eigenvalue; cos(φ + 2π/3) the smallest.
    let eig_max = q + 2.0 * p * phi.cos();
    let eig_min = q + 2.0 * p * (phi + two_thirds_pi).cos();
    (eig_min, eig_max)
}

/// The three principal stresses (eigenvalues) of a symmetric Cauchy stress
/// tensor in Voigt order, sorted descending `[σ₁, σ₂, σ₃]` (`σ₁ ≥ σ₂ ≥ σ₃`).
/// The same closed-form symmetric-3×3 solution as [`principal_extremes_voigt`],
/// recovering the middle eigenvalue from the trace (`σ₂ = tr − σ₁ − σ₃`).
fn principal_triple_voigt(s: [f64; 6]) -> [f64; 3] {
    let (sxx, syy, szz, sxy, syz, szx) = (s[0], s[1], s[2], s[3], s[4], s[5]);
    let p1 = sxy * sxy + syz * syz + szx * szx;
    if p1 <= 0.0 {
        // Already diagonal: sort the diagonal entries descending.
        let mut e = [sxx, syy, szz];
        e.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
        return e;
    }
    let q = (sxx + syy + szz) / 3.0;
    let a = sxx - q;
    let b = syy - q;
    let c = szz - q;
    let p2 = a * a + b * b + c * c + 2.0 * p1;
    let p = (p2 / 6.0).sqrt();
    let det = a * (b * c - syz * syz) - sxy * (sxy * c - syz * szx) + szx * (sxy * syz - b * szx);
    let r = (det / (p * p * p) / 2.0).clamp(-1.0, 1.0);
    let phi = r.acos() / 3.0;
    let two_thirds_pi = 2.0 * std::f64::consts::PI / 3.0;
    let s1 = q + 2.0 * p * phi.cos(); // largest (cos φ branch)
    let s3 = q + 2.0 * p * (phi + two_thirds_pi).cos(); // smallest
    let s2 = 3.0 * q - s1 - s3; // middle: the three eigenvalues sum to the trace 3q
    [s1, s2, s3]
}

/// Largest principal stress (the maximum normal stress, Rankine criterion).
fn max_principal_stress_voigt(s: [f64; 6]) -> f64 {
    principal_extremes_voigt(s).1
}

/// Smallest principal stress (the most compressive normal stress).
fn min_principal_stress_voigt(s: [f64; 6]) -> f64 {
    principal_extremes_voigt(s).0
}

/// Solve a linear-static elasticity problem on a tetrahedral mesh.
///
/// `mesh` must carry at least one [`ElementType::Tet4`] block (other
/// block types are ignored). `material` supplies the isotropic
/// elastic constants. `constraints` pin nodal DOFs (at least one is
/// required, or the system is rigid-body singular). `forces` apply
/// concentrated nodal loads.
///
/// Returns the per-node displacement field and the recovered nodal
/// stress.
///
/// # Method
///
/// 1. **Element stiffness.** For each linear tet, the four shape
///    functions are linear in `(x, y, z)`; their gradients are
///    constant and obtained from the inverse of the tet's
///    coordinate-difference matrix. The strain-displacement matrix
///    `B` (6×12) is assembled from those gradients; the element
///    stiffness is `Kₑ = Vₑ · Bᵀ · D · B` where `D` is the 6×6
///    isotropic elasticity matrix and `Vₑ` the tet volume.
/// 2. **Assembly.** Each `Kₑ` is scatter-added into a global sparse
///    matrix of size `3·n_nodes`.
/// 3. **Boundary conditions.** Concentrated forces go straight into
///    the load vector. Displacement constraints are imposed by the
///    large-penalty method — a stiff spring `k·penalty` is added to
///    the constrained DOF's diagonal and `k·penalty·value` to its
///    load entry. This keeps the matrix symmetric positive-definite
///    so the Cholesky path stays valid.
/// 4. **Solve.** The system is converted to CSC and factorised with
///    [`CscCholesky`]; the displacement vector is the back-substitution.
/// 5. **Stress recovery.** `σ = D·B·uₑ` per element (constant over a
///    linear tet); element stresses are averaged into their nodes.
///
/// # Errors
///
/// See [`NativeSolverError`] — bad material, no tet block, empty
/// mesh, degenerate element, out-of-range connectivity, an
/// unconstrained (singular) model, or a failed factorisation.
pub fn solve_linear_static(
    mesh: &Mesh,
    material: &FemMaterial,
    constraints: &[NodalConstraint],
    forces: &[NodalForce],
) -> Result<NativeSolution, NativeSolverError> {
    let n_nodes = mesh.nodes.len();
    if n_nodes == 0 {
        return Err(NativeSolverError::EmptyMesh);
    }
    let tets = tet_block(mesh).ok_or(NativeSolverError::NoTetBlock)?;
    if constraints.is_empty() {
        return Err(NativeSolverError::Unconstrained);
    }

    let d_matrix = elasticity_matrix(material)?;
    let n_dof = 3 * n_nodes;

    // --- assemble the global stiffness matrix ---
    let mut coo = CooMatrix::<f64>::new(n_dof, n_dof);
    // Track the largest diagonal magnitude so the BC penalty scales
    // with the problem (a fixed penalty is either too weak — leaks —
    // or too strong — wrecks conditioning).
    let mut max_diag = 0.0_f64;

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
        let ke =
            element_stiffness(&coords, &d_matrix).ok_or(NativeSolverError::DegenerateElement(e))?;
        // Scatter-add the 12×12 element matrix into the global system.
        // Local DOF `3a+i` maps to global DOF `3·nodes[a]+i`.
        for a in 0..4 {
            for i in 0..3 {
                let gi = 3 * nodes[a] + i;
                for b in 0..4 {
                    for j in 0..3 {
                        let gj = 3 * nodes[b] + j;
                        let v = ke[(3 * a + i, 3 * b + j)];
                        if v != 0.0 {
                            coo.push(gi, gj, v);
                        }
                    }
                }
            }
        }
    }

    // The CSC conversion sums duplicate (i, j) entries, so the
    // scatter-add above is correct even though each global DOF is
    // touched by several elements.
    let mut csc = CscMatrix::from(&coo);
    // Find the peak diagonal for penalty scaling.
    for (r, c, v) in csc.triplet_iter() {
        if r == c {
            max_diag = max_diag.max(v.abs());
        }
    }
    if max_diag <= 0.0 {
        // No stiffness at all — empty / disconnected mesh.
        return Err(NativeSolverError::SolveFailed);
    }
    // Penalty: ~1e8× the stiffest diagonal pins a DOF to ~8 digits
    // without destroying the Cholesky conditioning.
    let penalty = max_diag * 1.0e8;

    // --- load vector ---
    let mut f = DVector::<f64>::zeros(n_dof);
    for force in forces {
        if force.node >= n_nodes {
            return Err(NativeSolverError::BadConnectivity {
                elem: usize::MAX,
                node: force.node,
                n_nodes,
            });
        }
        // Round-1 H4: reject a non-finite force before it reaches the
        // RHS — the Cholesky solve would otherwise back-substitute it
        // into a silently-NaN displacement returned as Ok(..).
        ensure_finite3(&force.force, force.node, "force")?;
        for i in 0..3 {
            f[3 * force.node + i] += force.force[i];
        }
    }

    // --- displacement BCs via the penalty method ---
    // Add `penalty` to the constrained diagonal and `penalty·value`
    // to the load: the row then reads ≈ `penalty·u = penalty·value`,
    // i.e. `u ≈ value`, while the matrix stays SPD.
    let mut penalty_diag = vec![0.0_f64; n_dof];
    let mut any_constraint = false;
    for c in constraints {
        if c.node >= n_nodes {
            return Err(NativeSolverError::BadConnectivity {
                elem: usize::MAX,
                node: c.node,
                n_nodes,
            });
        }
        for i in 0..3 {
            if let Some(value) = c.fixed[i] {
                // Round-1 H4: a prescribed displacement is folded into
                // the RHS as `penalty·value`; a non-finite value
                // corrupts the solve just like a non-finite force.
                if !value.is_finite() {
                    return Err(NativeSolverError::InvalidLoad {
                        node: c.node,
                        kind: "prescribed displacement",
                    });
                }
                let dof = 3 * c.node + i;
                penalty_diag[dof] += penalty;
                f[dof] += penalty * value;
                any_constraint = true;
            }
        }
    }
    if !any_constraint {
        return Err(NativeSolverError::Unconstrained);
    }

    // Fold the penalty diagonal back into the CSC matrix. Rebuild via
    // a COO so we can add to existing diagonal entries cleanly.
    let stiffened = add_diagonal(&csc, &penalty_diag);
    csc = stiffened;

    // --- solve ---
    let chol = CscCholesky::factor(&csc).map_err(|_| NativeSolverError::SolveFailed)?;
    let u: DMatrix<f64> = chol.solve(&f);
    let u = u.column(0);

    // --- gather the nodal displacement field ---
    let mut displacement = vec![[0.0_f64; 3]; n_nodes];
    for (n, d) in displacement.iter_mut().enumerate() {
        d[0] = u[3 * n];
        d[1] = u[3 * n + 1];
        d[2] = u[3 * n + 2];
    }

    // --- stress recovery: σ = D·B·uₑ per element, averaged to nodes ---
    let mut nodal_stress = vec![[0.0_f64; 6]; n_nodes];
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
        let Some(b) = strain_displacement(&coords) else {
            continue;
        };
        // Element displacement vector (12×1).
        let mut ue = DVector::<f64>::zeros(12);
        for a in 0..4 {
            for i in 0..3 {
                ue[3 * a + i] = u[3 * nodes[a] + i];
            }
        }
        // ε = B·uₑ  (6×1), σ = D·ε  (6×1).
        let strain = &b * &ue;
        let sigma = d_matrix * &strain;
        for &nd in &nodes {
            for k in 0..6 {
                nodal_stress[nd][k] += sigma[k];
            }
            nodal_count[nd] += 1;
        }
    }
    let mut stress = vec![[0.0_f64; 6]; n_nodes];
    let mut von_mises = vec![0.0_f64; n_nodes];
    for n in 0..n_nodes {
        let count = nodal_count[n].max(1) as f64;
        for k in 0..6 {
            stress[n][k] = nodal_stress[n][k] / count;
        }
        von_mises[n] = von_mises_from_voigt(&stress[n]);
    }

    Ok(NativeSolution {
        displacement,
        von_mises,
        stress,
    })
}

/// Solve a linear-static elasticity problem on a **mixed-element**
/// solid mesh.
///
/// The generic counterpart of [`solve_linear_static`]: instead of
/// requiring a single [`ElementType::Tet4`] block, this walks every
/// supported continuum-element block — [`ElementType::Tet4`],
/// [`ElementType::Hex8`] and [`ElementType::Tet10`], in any mix — and
/// assembles them into one global system through the
/// [`crate::elements::SolidElement`] trait and
/// [`crate::assembly`]. A mesh that carries a Hex8 block *and* a Tet10
/// block solves as a single coupled model.
///
/// `material`, `constraints` and `forces` have the same meaning as in
/// [`solve_linear_static`]. The method is identical — global assembly,
/// the large-penalty treatment of displacement BCs, a sparse Cholesky
/// solve, then `σ = D·B·uₑ` stress recovery at each element centroid
/// averaged to the nodes.
///
/// For a pure-Tet4 mesh this returns the same answer as
/// [`solve_linear_static`] (the Tet4 element is the same constant-
/// strain tetrahedron); the value of this entry point is the Hex8 and
/// Tet10 elements, which are far more accurate in bending.
///
/// # Errors
///
/// See [`NativeSolverError`].
pub fn solve_linear_static_mixed(
    mesh: &Mesh,
    material: &FemMaterial,
    constraints: &[NodalConstraint],
    forces: &[NodalForce],
) -> Result<NativeSolution, NativeSolverError> {
    use crate::assembly::{
        assemble_global_stiffness_sparse, has_solid_elements, recover_nodal_stress_mixed,
    };

    let n_nodes = mesh.nodes.len();
    if n_nodes == 0 {
        return Err(NativeSolverError::EmptyMesh);
    }
    if !has_solid_elements(mesh) {
        return Err(NativeSolverError::NoTetBlock);
    }
    if constraints.is_empty() {
        return Err(NativeSolverError::Unconstrained);
    }
    let d_matrix = elasticity_matrix(material)?;
    let n_dof = 3 * n_nodes;

    // --- assemble the global stiffness (sparse, mixed-element) ---
    let k_sparse = assemble_global_stiffness_sparse(mesh, &d_matrix)?;

    // Peak diagonal → penalty scale.
    let mut max_diag = 0.0_f64;
    for (r, c, v) in k_sparse.triplet_iter() {
        if r == c {
            max_diag = max_diag.max(v.abs());
        }
    }
    if max_diag <= 0.0 {
        return Err(NativeSolverError::SolveFailed);
    }
    let penalty = max_diag * 1.0e8;

    // --- load vector ---
    let mut f = DVector::<f64>::zeros(n_dof);
    for force in forces {
        if force.node >= n_nodes {
            return Err(NativeSolverError::BadConnectivity {
                elem: usize::MAX,
                node: force.node,
                n_nodes,
            });
        }
        // Round-1 H4: reject a non-finite force up front.
        ensure_finite3(&force.force, force.node, "force")?;
        for i in 0..3 {
            f[3 * force.node + i] += force.force[i];
        }
    }

    // --- penalty displacement BCs ---
    let mut penalty_diag = vec![0.0_f64; n_dof];
    let mut any_constraint = false;
    for c in constraints {
        if c.node >= n_nodes {
            return Err(NativeSolverError::BadConnectivity {
                elem: usize::MAX,
                node: c.node,
                n_nodes,
            });
        }
        for i in 0..3 {
            if let Some(value) = c.fixed[i] {
                // Round-1 H4: reject a non-finite prescribed displacement.
                if !value.is_finite() {
                    return Err(NativeSolverError::InvalidLoad {
                        node: c.node,
                        kind: "prescribed displacement",
                    });
                }
                let dof = 3 * c.node + i;
                penalty_diag[dof] += penalty;
                f[dof] += penalty * value;
                any_constraint = true;
            }
        }
    }
    if !any_constraint {
        return Err(NativeSolverError::Unconstrained);
    }

    // --- fold the penalty diagonal in ---
    let csc = add_diagonal(&k_sparse, &penalty_diag);

    // --- fill-reducing reorder, factorise, solve, unpermute ---
    //
    // CscCholesky does no fill reduction, so a mesh whose node
    // numbering scatters coupled DOFs (notably the Tet10 mesh, whose
    // mid-edge nodes follow all corner nodes) would produce a near-
    // dense factor. A Reverse Cuthill-McKee node permutation shrinks
    // the bandwidth so the sparse factor stays sparse.
    use crate::ordering::{
        node_adjacency, node_perm_to_dof_perm, permute_csc, permute_vector, reverse_cuthill_mckee,
        unpermute_vector,
    };
    let elem_nodes = crate::assembly::solid_element_node_lists(mesh);
    let adjacency = node_adjacency(n_nodes, &elem_nodes);
    let node_perm = reverse_cuthill_mckee(&adjacency);
    let dof_perm = node_perm_to_dof_perm(&node_perm, 3);

    let csc_p = permute_csc(&csc, &dof_perm);
    let f_p = permute_vector(f.as_slice(), &dof_perm);
    let chol = CscCholesky::factor(&csc_p).map_err(|_| NativeSolverError::SolveFailed)?;
    let u_p: DMatrix<f64> = chol.solve(&DVector::from_vec(f_p));
    let u = DVector::from_vec(unpermute_vector(u_p.column(0).as_slice(), &dof_perm));

    let mut displacement = vec![[0.0_f64; 3]; n_nodes];
    for (n, d) in displacement.iter_mut().enumerate() {
        d[0] = u[3 * n];
        d[1] = u[3 * n + 1];
        d[2] = u[3 * n + 2];
    }

    // --- stress recovery (mixed-element, centroid → nodes) ---
    let stress = recover_nodal_stress_mixed(mesh, &d_matrix, &u);
    let von_mises: Vec<f64> = stress.iter().map(von_mises_from_voigt).collect();

    Ok(NativeSolution {
        displacement,
        von_mises,
        stress,
    })
}

/// First [`ElementType::Tet4`] block in the mesh, if any.
pub(crate) fn tet_block(mesh: &Mesh) -> Option<&ElementBlock> {
    mesh.element_blocks
        .iter()
        .find(|b| b.element_type == ElementType::Tet4)
        .filter(|b| b.connectivity.len() >= 4)
}

/// Isotropic 6×6 elasticity matrix `D` in Voigt notation
/// `[σxx σyy σzz σxy σyz σzx] = D · [εxx εyy εzz γxy γyz γzx]`.
///
/// Uses engineering shear strain (`γ = 2ε`); the lower-right block is
/// therefore `G` (not `2G`).
pub(crate) fn elasticity_matrix(m: &FemMaterial) -> Result<Matrix6<f64>, NativeSolverError> {
    let e = m.youngs_modulus;
    let nu = m.poisson_ratio;
    if !(e.is_finite()) || e <= 0.0 {
        return Err(NativeSolverError::BadMaterial(format!(
            "Young's modulus must be finite and positive, got {e}"
        )));
    }
    // D stays positive-definite only for -1 < ν < 0.5.
    if !nu.is_finite() || nu <= -1.0 || nu >= 0.5 {
        return Err(NativeSolverError::BadMaterial(format!(
            "Poisson's ratio must lie in (-1, 0.5), got {nu}"
        )));
    }
    let lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mu = e / (2.0 * (1.0 + nu)); // shear modulus G
    let mut d = Matrix6::zeros();
    // Normal-stress block.
    for i in 0..3 {
        for j in 0..3 {
            d[(i, j)] = lambda;
        }
        d[(i, i)] = lambda + 2.0 * mu;
    }
    // Shear block (engineering strain).
    for i in 3..6 {
        d[(i, i)] = mu;
    }
    Ok(d)
}

/// Shape-function gradients of a linear tet, packed as the 6×12
/// strain-displacement matrix `B`. Returns `None` if the tet is
/// degenerate.
///
/// The four linear shape functions satisfy `Nₐ(xᵦ) = δₐᵦ`. Their
/// gradients `∇Nₐ = (bₐ, cₐ, dₐ)` are constant and come from the
/// inverse of the 4×4 matrix `[1 x y z]`. `B` maps the 12 nodal
/// displacement components to the 6 Voigt strain components.
pub(crate) fn strain_displacement(coords: &[Vector3<f64>; 4]) -> Option<DMatrix<f64>> {
    let grads = shape_gradients(coords)?;
    let mut b = DMatrix::<f64>::zeros(6, 12);
    for (a, g) in grads.iter().enumerate() {
        let (bx, by, bz) = (g.x, g.y, g.z);
        let col = 3 * a;
        // εxx = ∂u/∂x
        b[(0, col)] = bx;
        // εyy = ∂v/∂y
        b[(1, col + 1)] = by;
        // εzz = ∂w/∂z
        b[(2, col + 2)] = bz;
        // γxy = ∂u/∂y + ∂v/∂x
        b[(3, col)] = by;
        b[(3, col + 1)] = bx;
        // γyz = ∂v/∂z + ∂w/∂y
        b[(4, col + 1)] = bz;
        b[(4, col + 2)] = by;
        // γzx = ∂w/∂x + ∂u/∂z
        b[(5, col)] = bz;
        b[(5, col + 2)] = bx;
    }
    Some(b)
}

/// Constant shape-function gradients `∇Nₐ` of a linear tet (one
/// [`Vector3`] per node). `None` if the tet is degenerate.
pub(crate) fn shape_gradients(coords: &[Vector3<f64>; 4]) -> Option<[Vector3<f64>; 4]> {
    // The gradient of shape function `a` is the `a`-th row of the
    // inverse of `J = [[x1-x0, x2-x0, x3-x0] ...]ᵀ` applied to the
    // reference gradients. Equivalent and numerically simple: build
    // the 3×3 edge matrix and invert it.
    let p0 = coords[0];
    let e1 = coords[1] - p0;
    let e2 = coords[2] - p0;
    let e3 = coords[3] - p0;
    // Jacobian columns are the edge vectors.
    let jac = Matrix3::from_columns(&[e1, e2, e3]);
    let det = jac.determinant();
    // Tet volume = |det| / 6; reject near-zero.
    if det.abs() < 1.0e-18 {
        return None;
    }
    let inv = jac.try_inverse()?;
    // Reference gradients of N1,N2,N3 are the rows of the identity in
    // reference space; N0 = 1 - N1 - N2 - N3.
    // ∇Nₐ (physical) = invᵀ · (reference gradient).
    let invt = inv.transpose();
    let g1 = invt * Vector3::new(1.0, 0.0, 0.0);
    let g2 = invt * Vector3::new(0.0, 1.0, 0.0);
    let g3 = invt * Vector3::new(0.0, 0.0, 1.0);
    let g0 = -(g1 + g2 + g3);
    Some([g0, g1, g2, g3])
}

/// Signed volume of a tet (positive for a right-handed node order).
pub(crate) fn tet_volume(coords: &[Vector3<f64>; 4]) -> f64 {
    let e1 = coords[1] - coords[0];
    let e2 = coords[2] - coords[0];
    let e3 = coords[3] - coords[0];
    e1.cross(&e2).dot(&e3) / 6.0
}

/// 12×12 element stiffness `Kₑ = Vₑ · Bᵀ · D · B`. `None` if the tet
/// is degenerate.
pub(crate) fn element_stiffness(
    coords: &[Vector3<f64>; 4],
    d_matrix: &Matrix6<f64>,
) -> Option<DMatrix<f64>> {
    let b = strain_displacement(coords)?;
    let vol = tet_volume(coords).abs();
    if vol < 1.0e-18 {
        return None;
    }
    // Kₑ = V · Bᵀ D B. The integrand is constant over a linear tet,
    // so a single point (the element volume) integrates it exactly.
    let bt_d = b.transpose() * d_matrix;
    let ke = (bt_d * b) * vol;
    Some(ke)
}

/// Return a copy of `csc` with `diag[i]` added to entry `(i, i)`.
///
/// Diagonal entries already present are incremented; missing ones are
/// inserted. Used to fold the boundary-condition penalty into the
/// assembled stiffness matrix.
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

/// von Mises equivalent stress from a Voigt stress vector
/// `[σxx σyy σzz σxy σyz σzx]`.
pub(crate) fn von_mises_from_voigt(s: &[f64; 6]) -> f64 {
    let (sx, sy, sz, sxy, syz, szx) = (s[0], s[1], s[2], s[3], s[4], s[5]);
    (0.5 * ((sx - sy).powi(2) + (sy - sz).powi(2) + (sz - sx).powi(2))
        + 3.0 * (sxy * sxy + syz * syz + szx * szx))
        .sqrt()
}

/// Build a structured tetrahedral mesh of an axis-aligned box.
///
/// The box `[0, lx] × [0, ly] × [0, lz]` is divided into
/// `nx × ny × nz` hexahedral cells; each cell is split into **6
/// tetrahedra** by the standard Kuhn / diagonal triangulation (every
/// tet shares the cell's main diagonal). The result is a conforming
/// [`ElementType::Tet4`] volume mesh with `(nx+1)(ny+1)(nz+1)` nodes —
/// exactly the input the native solver wants for testing and for
/// box-domain problems without an external mesher.
///
/// Nodes are laid out in `x`-fastest, then `y`, then `z` order, so the
/// node index of grid point `(i, j, k)` is
/// `i + (nx+1)·j + (nx+1)(ny+1)·k`. Callers use that formula to find
/// the node ids of the faces they want to constrain or load.
///
/// # Panics
///
/// # Panics
///
/// Construct a structured hex-mesh of a unit box with `nx`/`ny`/`nz`
/// subdivisions. Each hex is split into 6 Kuhn tets so the output
/// mesh is a single [`ElementBlock`] of [`ElementType::Tet4`].
///
/// # Errors
///
/// Returns [`NativeSolverError::InvalidParams`] when any of `nx`, `ny`, `nz` is
/// zero, when any length is non-positive / non-finite, or when the
/// implied node-count `(nx+1)*(ny+1)*(nz+1)` would overflow `usize`
/// (32-bit platforms can hit this with surprisingly small inputs —
/// `1000^3` already overflows `usize::MAX` on Wasm32, and the cap
/// rejects it before any allocation is attempted).
///
/// Round-5 cleanup: this used to be `try_structured_box_mesh`, with a
/// panicking `structured_box_mesh` wrapper that was bound to UI knobs
/// and crashed the process on `nx=0`. The panicking variant has been
/// deleted and all 45 internal callers updated to handle the Result.
pub fn structured_box_mesh(
    lx: f64,
    ly: f64,
    lz: f64,
    nx: usize,
    ny: usize,
    nz: usize,
) -> Result<Mesh, NativeSolverError> {
    if nx == 0 || ny == 0 || nz == 0 {
        return Err(NativeSolverError::InvalidParams {
            reason: format!("subdivisions must be positive (got nx={nx}, ny={ny}, nz={nz})"),
        });
    }
    if !(lx.is_finite() && ly.is_finite() && lz.is_finite()) {
        return Err(NativeSolverError::InvalidParams {
            reason: format!("box dimensions must be finite (got lx={lx}, ly={ly}, lz={lz})"),
        });
    }
    if !(lx > 0.0 && ly > 0.0 && lz > 0.0) {
        return Err(NativeSolverError::InvalidParams {
            reason: format!("box dimensions must be positive (got lx={lx}, ly={ly}, lz={lz})"),
        });
    }
    // Round-3 fix: bound the implied allocations BEFORE we walk the
    // `for k in 0..=nz` loops. `(nx+1) * (ny+1) * (nz+1)` is the
    // node count; raw `*` would silently wrap on 32-bit targets and
    // produce a tiny mesh whose later indexing then panics. checked_mul
    // surfaces the overflow as a structured error.
    let nx1 = nx
        .checked_add(1)
        .ok_or_else(|| NativeSolverError::InvalidParams {
            reason: format!("nx+1 overflowed usize (nx={nx})"),
        })?;
    let ny1 = ny
        .checked_add(1)
        .ok_or_else(|| NativeSolverError::InvalidParams {
            reason: format!("ny+1 overflowed usize (ny={ny})"),
        })?;
    let nz1 = nz
        .checked_add(1)
        .ok_or_else(|| NativeSolverError::InvalidParams {
            reason: format!("nz+1 overflowed usize (nz={nz})"),
        })?;
    let node_count = nx1
        .checked_mul(ny1)
        .and_then(|v| v.checked_mul(nz1))
        .ok_or_else(|| NativeSolverError::InvalidParams {
            reason: format!(
                "node count (nx+1)*(ny+1)*(nz+1) overflowed usize \
                 (nx={nx}, ny={ny}, nz={nz})"
            ),
        })?;
    // Element count is nx*ny*nz hex cells × 6 tets × 4 nodes = 24x
    // entries in the connectivity vec. Check the multiplication.
    let cell_count = nx
        .checked_mul(ny)
        .and_then(|v| v.checked_mul(nz))
        .ok_or_else(|| NativeSolverError::InvalidParams {
            reason: "cell count nx*ny*nz overflowed usize".to_string(),
        })?;
    let _conn_capacity =
        cell_count
            .checked_mul(24)
            .ok_or_else(|| NativeSolverError::InvalidParams {
                reason: "connectivity capacity nx*ny*nz*24 overflowed usize".to_string(),
            })?;
    // Sanity cap on total node count — beyond this the implied
    // allocation is in the multi-gigabyte range and definitely a
    // programmer error rather than a real mesh. Mirrors MAX_MODAL_SIZE.
    const MAX_NODES: usize = 200_000_000;
    if node_count > MAX_NODES {
        return Err(NativeSolverError::InvalidParams {
            reason: format!(
                "structured box mesh would have {node_count} nodes, \
                 cap is {MAX_NODES}"
            ),
        });
    }
    let mut mesh = Mesh::new("native-fem-box");
    let stride_y = nx + 1;
    let stride_z = (nx + 1) * (ny + 1);
    let node_id =
        |i: usize, j: usize, k: usize| -> u32 { (i + stride_y * j + stride_z * k) as u32 };
    // Nodes.
    for k in 0..=nz {
        for j in 0..=ny {
            for i in 0..=nx {
                mesh.nodes.push(Vector3::new(
                    lx * i as f64 / nx as f64,
                    ly * j as f64 / ny as f64,
                    lz * k as f64 / nz as f64,
                ));
            }
        }
    }
    // Six tets per hex cell. The eight corners of cell (i, j, k):
    //   n0=(i,j,k)     n1=(i+1,j,k)   n2=(i+1,j+1,k)   n3=(i,j+1,k)
    //   n4=(i,j,k+1)   n5=(i+1,j,k+1) n6=(i+1,j+1,k+1) n7=(i,j+1,k+1)
    // Kuhn triangulation around the n0-n6 main diagonal. Each tet
    // below is wound so its signed volume is positive.
    let mut conn: Vec<u32> = Vec::with_capacity(nx * ny * nz * 24);
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let n0 = node_id(i, j, k);
                let n1 = node_id(i + 1, j, k);
                let n2 = node_id(i + 1, j + 1, k);
                let n3 = node_id(i, j + 1, k);
                let n4 = node_id(i, j, k + 1);
                let n5 = node_id(i + 1, j, k + 1);
                let n6 = node_id(i + 1, j + 1, k + 1);
                let n7 = node_id(i, j + 1, k + 1);
                // Kuhn's six tets, all sharing edge n0-n6.
                for tet in [
                    [n0, n1, n2, n6],
                    [n0, n2, n3, n6],
                    [n0, n3, n7, n6],
                    [n0, n7, n4, n6],
                    [n0, n4, n5, n6],
                    [n0, n5, n1, n6],
                ] {
                    conn.extend_from_slice(&tet);
                }
            }
        }
    }
    mesh.element_blocks.push(ElementBlock {
        element_type: ElementType::Tet4,
        connectivity: conn,
    });
    mesh.recompute_stats();
    Ok(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Node id of grid point `(i, j, k)` in a `structured_box_mesh`.
    fn nid(i: usize, j: usize, k: usize, nx: usize, ny: usize) -> usize {
        i + (nx + 1) * j + (nx + 1) * (ny + 1) * k
    }

    #[test]
    fn max_principal_stress_of_known_tensors() {
        // Diagonal: principals are the diagonal entries → max is the largest.
        assert!((max_principal_stress_voigt([10.0, 5.0, 2.0, 0.0, 0.0, 0.0]) - 10.0).abs() < 1e-9);
        // Uniaxial tension σxx = 7 → max principal = 7.
        assert!((max_principal_stress_voigt([7.0, 0.0, 0.0, 0.0, 0.0, 0.0]) - 7.0).abs() < 1e-9);
        // Hydrostatic: every principal equals the pressure.
        assert!((max_principal_stress_voigt([4.0, 4.0, 4.0, 0.0, 0.0, 0.0]) - 4.0).abs() < 1e-9);
        // Pure shear σxy = τ → principals ±τ, 0 → max = τ.
        assert!((max_principal_stress_voigt([0.0, 0.0, 0.0, 3.0, 0.0, 0.0]) - 3.0).abs() < 1e-9);
        // A 2-D state (σxx=3, σyy=1, σxy=2): Mohr's circle gives 2 ± √5, so the
        // maximum principal stress is 2 + √5 ≈ 4.236.
        let mohr = max_principal_stress_voigt([3.0, 1.0, 0.0, 2.0, 0.0, 0.0]);
        assert!(
            (mohr - (2.0 + 5.0_f64.sqrt())).abs() < 1e-9,
            "Mohr max principal {mohr}"
        );
    }

    #[test]
    fn min_principal_stress_of_known_tensors() {
        // Diagonal: min principal is the smallest diagonal entry.
        assert!((min_principal_stress_voigt([10.0, 5.0, 2.0, 0.0, 0.0, 0.0]) - 2.0).abs() < 1e-9);
        // Uniaxial tension σxx = 7 (σyy = σzz = 0): min principal = 0.
        assert!(min_principal_stress_voigt([7.0, 0.0, 0.0, 0.0, 0.0, 0.0]).abs() < 1e-9);
        // Hydrostatic: every principal equals the pressure.
        assert!((min_principal_stress_voigt([4.0, 4.0, 4.0, 0.0, 0.0, 0.0]) - 4.0).abs() < 1e-9);
        // Pure shear σxy = τ → principals τ, 0, −τ → min = −τ.
        assert!((min_principal_stress_voigt([0.0, 0.0, 0.0, 3.0, 0.0, 0.0]) - (-3.0)).abs() < 1e-9);
        // Mohr's circle (σxx=3, σyy=1, σxy=2): principals 2 ± √5 and 0, so the
        // minimum principal stress is 2 − √5 ≈ −0.236.
        let mohr = min_principal_stress_voigt([3.0, 1.0, 0.0, 2.0, 0.0, 0.0]);
        assert!(
            (mohr - (2.0 - 5.0_f64.sqrt())).abs() < 1e-9,
            "Mohr min principal {mohr}"
        );
    }

    #[test]
    fn max_shear_stress_is_the_peak_half_principal_range() {
        let sol = NativeSolution {
            displacement: vec![],
            von_mises: vec![],
            stress: vec![
                [10.0, 0.0, 0.0, 0.0, 0.0, 0.0], // uniaxial 10 → max shear 5
                [0.0, 0.0, 0.0, 4.0, 0.0, 0.0],  // pure shear 4 → max shear 4
            ],
        };
        // The peak per-node (σ₁ − σ₃)/2 is 5 (from the uniaxial node).
        assert!(
            (sol.max_shear_stress() - 5.0).abs() < 1e-9,
            "{}",
            sol.max_shear_stress()
        );
        // An empty field yields zero.
        let empty = NativeSolution {
            displacement: vec![],
            von_mises: vec![],
            stress: vec![],
        };
        assert_eq!(empty.max_shear_stress(), 0.0);
    }

    #[test]
    fn peak_lode_parameter_classifies_the_deviatoric_state() {
        // Single-node solutions so the peak-von-Mises node is unambiguous. The
        // Lode parameter classifies the deviatoric state at that critical point:
        // −1 for axisymmetric/uniaxial tension, 0 for pure shear, +1 for
        // axisymmetric/uniaxial compression.
        let mk = |stress: [f64; 6], vm: f64| NativeSolution {
            displacement: vec![[0.0, 0.0, 0.0]],
            von_mises: vec![vm],
            stress: vec![stress],
        };
        let sigma = 100.0e6;
        let tension = mk([sigma, 0.0, 0.0, 0.0, 0.0, 0.0], sigma);
        assert!(
            (tension.peak_lode_parameter() + 1.0).abs() < 1e-9,
            "tension {}",
            tension.peak_lode_parameter()
        );
        let compression = mk([-sigma, 0.0, 0.0, 0.0, 0.0, 0.0], sigma);
        assert!(
            (compression.peak_lode_parameter() - 1.0).abs() < 1e-9,
            "compression {}",
            compression.peak_lode_parameter()
        );
        let shear = mk([0.0, 0.0, 0.0, sigma, 0.0, 0.0], sigma);
        assert!(
            shear.peak_lode_parameter().abs() < 1e-9,
            "shear {}",
            shear.peak_lode_parameter()
        );
        // Always within [−1, 1].
        for s in [&tension, &compression, &shear] {
            assert!((-1.0..=1.0).contains(&s.peak_lode_parameter()));
        }
        // Empty solution → 0; an unstressed body (zero peak von Mises) → 0.
        let empty = NativeSolution {
            displacement: vec![],
            von_mises: vec![],
            stress: vec![],
        };
        assert_eq!(empty.peak_lode_parameter(), 0.0);
        let unstressed = mk([0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 0.0);
        assert_eq!(unstressed.peak_lode_parameter(), 0.0);
    }

    #[test]
    fn hydrostatic_stress_extremes_are_the_mean_normal_stress() {
        let sol = NativeSolution {
            displacement: vec![],
            von_mises: vec![],
            stress: vec![
                [9.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // uniaxial tension → σ_m = 3
                [-3.0, -3.0, -3.0, 0.0, 0.0, 0.0], // hydrostatic compression → σ_m = −3
                [0.0, 0.0, 0.0, 5.0, 0.0, 0.0],    // pure shear → σ_m = 0
            ],
        };
        // σ_m = trace/3 per node; the field extremes are +3 and −3.
        assert!(
            (sol.max_hydrostatic_stress() - 3.0).abs() < 1e-12,
            "{}",
            sol.max_hydrostatic_stress()
        );
        assert!(
            (sol.min_hydrostatic_stress() + 3.0).abs() < 1e-12,
            "{}",
            sol.min_hydrostatic_stress()
        );
    }

    #[test]
    fn peak_triaxiality_matches_textbook_uniaxial_and_biaxial_values() {
        let sigma = 200.0e6;
        // Uniaxial tension: σ_vm = σ, σ_m = σ/3 → T = 1/3.
        let uni = NativeSolution {
            displacement: vec![],
            von_mises: vec![sigma],
            stress: vec![[sigma, 0.0, 0.0, 0.0, 0.0, 0.0]],
        };
        assert!(
            (uni.peak_triaxiality() - 1.0 / 3.0).abs() < 1e-9,
            "uniaxial {}",
            uni.peak_triaxiality()
        );
        // Equibiaxial tension (σxx = σyy = σ): σ_vm = σ, σ_m = 2σ/3 → T = 2/3.
        let bi = NativeSolution {
            displacement: vec![],
            von_mises: vec![sigma],
            stress: vec![[sigma, sigma, 0.0, 0.0, 0.0, 0.0]],
        };
        assert!(
            (bi.peak_triaxiality() - 2.0 / 3.0).abs() < 1e-9,
            "equibiaxial {}",
            bi.peak_triaxiality()
        );
        // Evaluated at the PEAK von Mises node: a low-stress uniaxial node is
        // ignored in favour of the high-stress equibiaxial node.
        let mixed = NativeSolution {
            displacement: vec![],
            von_mises: vec![0.1 * sigma, sigma],
            stress: vec![
                [0.1 * sigma, 0.0, 0.0, 0.0, 0.0, 0.0],
                [sigma, sigma, 0.0, 0.0, 0.0, 0.0],
            ],
        };
        assert!(
            (mixed.peak_triaxiality() - 2.0 / 3.0).abs() < 1e-9,
            "peak-node"
        );
        // Empty / unstressed → 0, never a divide-by-zero.
        let empty = NativeSolution {
            displacement: vec![],
            von_mises: vec![],
            stress: vec![],
        };
        assert_eq!(empty.peak_triaxiality(), 0.0);
    }

    #[test]
    fn peak_von_mises_index_finds_the_most_stressed_node() {
        let sol = NativeSolution {
            displacement: vec![],
            von_mises: vec![1.0e6, 5.0e6, 3.0e6],
            stress: vec![],
        };
        assert_eq!(sol.peak_von_mises_index(), Some(1));
        // An empty field has no peak node.
        let empty = NativeSolution {
            displacement: vec![],
            von_mises: vec![],
            stress: vec![],
        };
        assert_eq!(empty.peak_von_mises_index(), None);
    }

    #[test]
    fn mean_von_mises_is_the_nodal_average_and_gives_a_concentration_factor() {
        let sol = NativeSolution {
            displacement: vec![],
            von_mises: vec![2.0e6, 4.0e6, 6.0e6],
            stress: vec![],
        };
        assert!(
            (sol.mean_von_mises() - 4.0e6).abs() < 1e-3,
            "{}",
            sol.mean_von_mises()
        );
        // Stress-concentration factor Kt = peak/mean = 6/4 = 1.5 here.
        assert!((sol.max_von_mises() / sol.mean_von_mises() - 1.5).abs() < 1e-9);
        // Empty field → 0 (no divide-by-zero).
        let empty = NativeSolution {
            displacement: vec![],
            von_mises: vec![],
            stress: vec![],
        };
        assert_eq!(empty.mean_von_mises(), 0.0);
    }

    #[test]
    fn single_cell_box_has_six_tets() {
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid params");
        assert_eq!(mesh.nodes.len(), 8);
        assert_eq!(mesh.element_blocks.len(), 1);
        assert_eq!(mesh.element_blocks[0].element_type, ElementType::Tet4);
        // 6 tets × 4 nodes.
        assert_eq!(mesh.element_blocks[0].connectivity.len(), 24);
    }

    /// Round-3 fix (now round-5 GREEN): `structured_box_mesh(0, 1, 1, ..)`
    /// used to panic via `assert!()`. The fallible variant returns a
    /// structured error so callers handling user input recover. Round-5:
    /// the panicking wrapper has been deleted entirely — there is now
    /// only one `structured_box_mesh` and it is fallible.
    #[test]
    fn structured_box_mesh_rejects_zero_subdivisions() {
        let err = structured_box_mesh(1.0, 1.0, 1.0, 0, 1, 1).expect_err("nx=0 must be rejected");
        let msg = format!("{err}");
        assert!(msg.contains("positive"), "msg: {msg}");
    }

    #[test]
    fn structured_box_mesh_rejects_zero_dimension() {
        let err = structured_box_mesh(0.0, 1.0, 1.0, 1, 1, 1).expect_err("lx=0 must be rejected");
        let msg = format!("{err}");
        assert!(msg.contains("positive"), "msg: {msg}");
    }

    #[test]
    fn structured_box_mesh_rejects_negative_dimension() {
        let err =
            structured_box_mesh(-1.0, 1.0, 1.0, 1, 1, 1).expect_err("negative lx must be rejected");
        let msg = format!("{err}");
        assert!(msg.contains("positive"), "msg: {msg}");
    }

    #[test]
    fn structured_box_mesh_rejects_non_finite_dimension() {
        let err =
            structured_box_mesh(f64::NAN, 1.0, 1.0, 1, 1, 1).expect_err("NaN lx must be rejected");
        let msg = format!("{err}");
        assert!(msg.contains("finite"), "msg: {msg}");
    }

    /// Round-3 fix: the cap also defends against overflow of the
    /// implied node count `(nx+1)*(ny+1)*(nz+1)`. We can't reach
    /// `usize::MAX` portably without running the test on Wasm32, but
    /// the 200M-node cap catches the asymptotic case on 64-bit too.
    #[test]
    fn structured_box_mesh_rejects_oversized_grid() {
        // 1000^3 = 10^9 nodes — over the 200M cap.
        let err = structured_box_mesh(1.0, 1.0, 1.0, 1000, 1000, 1000)
            .expect_err("1000^3 nodes must be over the cap");
        let msg = format!("{err}");
        assert!(msg.contains("cap") || msg.contains("nodes"), "msg: {msg}");
    }

    #[test]
    fn box_mesh_tets_all_have_positive_volume() {
        // Every tet from the Kuhn triangulation must be right-handed.
        let mesh = structured_box_mesh(2.0, 3.0, 1.5, 2, 2, 2).expect("valid box params");
        let conn = &mesh.element_blocks[0].connectivity;
        for e in 0..conn.len() / 4 {
            let c = [
                mesh.nodes[conn[4 * e] as usize],
                mesh.nodes[conn[4 * e + 1] as usize],
                mesh.nodes[conn[4 * e + 2] as usize],
                mesh.nodes[conn[4 * e + 3] as usize],
            ];
            let v = tet_volume(&c);
            assert!(v > 0.0, "tet {e} has non-positive volume {v}");
        }
        // The six tets of a unit cell partition its volume exactly.
        let total: f64 = (0..conn.len() / 4)
            .map(|e| {
                tet_volume(&[
                    mesh.nodes[conn[4 * e] as usize],
                    mesh.nodes[conn[4 * e + 1] as usize],
                    mesh.nodes[conn[4 * e + 2] as usize],
                    mesh.nodes[conn[4 * e + 3] as usize],
                ])
            })
            .sum();
        assert!((total - 2.0 * 3.0 * 1.5).abs() < 1e-9, "volume sum {total}");
    }

    #[test]
    fn elasticity_matrix_is_symmetric_and_positive_diagonal() {
        let d = elasticity_matrix(&FemMaterial::default()).unwrap();
        for i in 0..6 {
            assert!(d[(i, i)] > 0.0);
            for j in 0..6 {
                assert!((d[(i, j)] - d[(j, i)]).abs() < 1e-3 * d[(i, i)].abs().max(1.0));
            }
        }
    }

    #[test]
    fn elasticity_matrix_rejects_bad_poisson() {
        // Incompressible (ν = 0.5) — the elasticity tensor is singular.
        let incompressible = FemMaterial {
            poisson_ratio: 0.5,
            ..FemMaterial::default()
        };
        assert!(matches!(
            elasticity_matrix(&incompressible),
            Err(NativeSolverError::BadMaterial(_))
        ));
        // Non-physical negative modulus.
        let negative_e = FemMaterial {
            poisson_ratio: 0.3,
            youngs_modulus: -1.0,
            ..FemMaterial::default()
        };
        assert!(matches!(
            elasticity_matrix(&negative_e),
            Err(NativeSolverError::BadMaterial(_))
        ));
    }

    #[test]
    fn element_stiffness_is_symmetric() {
        // A reference unit tet.
        let coords = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let d = elasticity_matrix(&FemMaterial::default()).unwrap();
        let ke = element_stiffness(&coords, &d).unwrap();
        assert_eq!(ke.nrows(), 12);
        for i in 0..12 {
            for j in 0..12 {
                let asym = (ke[(i, j)] - ke[(j, i)]).abs();
                assert!(asym < 1e-3 * ke[(i, i)].abs().max(1.0), "Ke not symmetric");
            }
        }
    }

    #[test]
    fn element_stiffness_has_zero_energy_rigid_modes() {
        // A rigid-body translation produces zero strain → zero force.
        let coords = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let d = elasticity_matrix(&FemMaterial::default()).unwrap();
        let ke = element_stiffness(&coords, &d).unwrap();
        // Uniform +X translation of all four nodes.
        let mut t = DVector::<f64>::zeros(12);
        for a in 0..4 {
            t[3 * a] = 1.0;
        }
        let f = &ke * &t;
        assert!(f.norm() < 1e-3, "rigid translation should yield no force");
    }

    #[test]
    fn degenerate_tet_is_rejected() {
        // Four coplanar points — zero volume.
        let coords = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
        ];
        let d = elasticity_matrix(&FemMaterial::default()).unwrap();
        assert!(element_stiffness(&coords, &d).is_none());
    }

    #[test]
    fn solve_rejects_unconstrained_model() {
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box params");
        let err = solve_linear_static(&mesh, &FemMaterial::default(), &[], &[]).unwrap_err();
        assert!(matches!(err, NativeSolverError::Unconstrained));
    }

    #[test]
    fn solve_rejects_mesh_without_tets() {
        let mut mesh = Mesh::new("surface-only");
        mesh.nodes.push(Vector3::new(0.0, 0.0, 0.0));
        mesh.element_blocks.push(ElementBlock {
            element_type: ElementType::Tri3,
            connectivity: vec![0, 0, 0],
        });
        let err = solve_linear_static(
            &mesh,
            &FemMaterial::default(),
            &[NodalConstraint::fixed(0)],
            &[],
        )
        .unwrap_err();
        assert!(matches!(err, NativeSolverError::NoTetBlock));
    }

    #[test]
    fn uniaxial_tension_reproduces_hookes_law() {
        // Hooke's-law PATCH TEST for the linear tet. A bar of length L
        // along X, square 1×1 cross section. The x=0 face is fixed in X
        // and the x=L face is given a prescribed elongation ΔL; the
        // recovered field must be the exact uniform-strain state
        // ε = ΔL/L, σ = E·ε, and a constant-strain field is reproduced
        // exactly by linear tetrahedra.
        //
        // NB: this uses a *displacement-driven* load on purpose. A
        // force-driven test that puts equal point loads on the four
        // face nodes is wrong for a Kuhn-split tet mesh — the consistent
        // nodal forces of a uniform face traction are unequal on an
        // asymmetric triangulation, so equal point loads give a
        // (correct) NON-uniform response and would not match ΔL=FL/EA.
        let (lx, ly, lz) = (4.0, 1.0, 1.0);
        let (nx, ny, nz) = (4, 1, 1);
        let mesh = structured_box_mesh(lx, ly, lz, nx, ny, nz).expect("valid box params");
        // ν = 0 → pure uniaxial, no lateral (Poisson) coupling.
        let mat = FemMaterial {
            youngs_modulus: 200.0e9,
            poisson_ratio: 0.0,
            ..FemMaterial::default()
        };

        // Analytic state for an applied total force F = 1e6 N over
        // A = ly·lz: ΔL = F·L/(E·A), ε = ΔL/L, σ = E·ε = F/A.
        let area = ly * lz;
        let total_f = 1.0e6;
        let expected_dl = total_f * lx / (mat.youngs_modulus * area);

        // x=0 face fixed in X (one node fully pinned to kill rigid-body
        // modes); x=L face given the prescribed elongation in X.
        let mut constraints = Vec::new();
        for k in 0..=nz {
            for j in 0..=ny {
                let n0 = nid(0, j, k, nx, ny);
                let c0 = if j == 0 && k == 0 {
                    NodalConstraint::fixed(n0)
                } else {
                    NodalConstraint {
                        node: n0,
                        fixed: [Some(0.0), None, None],
                    }
                };
                constraints.push(c0);
                constraints.push(NodalConstraint {
                    node: nid(nx, j, k, nx, ny),
                    fixed: [Some(expected_dl), None, None],
                });
            }
        }

        let sol = solve_linear_static(&mesh, &mat, &constraints, &[]).unwrap();

        // Every node's X displacement must be exactly linear in x — the
        // patch-test pass criterion (a constant-strain field). The
        // tolerance allows for the penalty method's ~8-digit pinning of
        // the prescribed DOFs.
        for (n, p) in mesh.nodes.iter().enumerate() {
            let expect = expected_dl * p[0] / lx;
            let ux = sol.displacement[n][0];
            assert!(
                (ux - expect).abs() < 1e-6 * expected_dl,
                "node {n} u_x {ux} vs linear field {expect}"
            );
        }

        // Normal stress σxx = F/A everywhere (constant over the bar).
        let interior = nid(2, 0, 0, nx, ny);
        let sigma_xx = sol.stress[interior][0];
        let expected_sigma = total_f / area;
        let rel_s = (sigma_xx - expected_sigma).abs() / expected_sigma;
        assert!(
            rel_s < 1e-6,
            "σxx {sigma_xx} vs analytic {expected_sigma} (rel {rel_s})"
        );
    }

    #[test]
    fn cantilever_beam_deflects_in_the_load_direction() {
        // A slender beam along X, clamped at x=0, loaded with a
        // downward (-Z) tip force. The free end must deflect downward,
        // and the deflection must be far larger than any spurious
        // sideways (Y) motion. We do not assert the exact Euler-
        // Bernoulli value (linear tets are too stiff in bending on a
        // coarse mesh) — only the qualitative, sign-correct behaviour.
        let (lx, ly, lz) = (10.0, 1.0, 1.0);
        let (nx, ny, nz) = (10, 1, 1);
        let mesh = structured_box_mesh(lx, ly, lz, nx, ny, nz).expect("valid box params");
        let mat = FemMaterial::default();

        let mut constraints = Vec::new();
        for k in 0..=nz {
            for j in 0..=ny {
                constraints.push(NodalConstraint::fixed(nid(0, j, k, nx, ny)));
            }
        }
        let tip_nodes: Vec<usize> = (0..=nz)
            .flat_map(|k| (0..=ny).map(move |j| (j, k)))
            .map(|(j, k)| nid(nx, j, k, nx, ny))
            .collect();
        let forces: Vec<NodalForce> = tip_nodes
            .iter()
            .map(|&n| NodalForce {
                node: n,
                force: [0.0, 0.0, -2.5e4],
            })
            .collect();

        let sol = solve_linear_static(&mesh, &mat, &constraints, &forces).unwrap();
        let tip_dz = tip_nodes
            .iter()
            .map(|&n| sol.displacement[n][2])
            .sum::<f64>()
            / tip_nodes.len() as f64;
        assert!(tip_dz < 0.0, "tip should deflect downward, got {tip_dz}");
        // The clamped end must not move.
        let clamped_dz = sol.displacement[nid(0, 0, 0, nx, ny)][2];
        assert!(clamped_dz.abs() < 1e-9, "clamped end moved: {clamped_dz}");
        // A finer beam deflects more (a longer lever / more elements).
        assert!(sol.max_displacement() > tip_dz.abs() * 0.5);
    }

    #[test]
    fn cantilever_refinement_converges_toward_beam_theory() {
        // Euler-Bernoulli tip deflection of an end-loaded cantilever:
        //   δ = P·L³ / (3·E·I),  I = b·h³/12 for a rectangular section.
        // A linear-tet mesh under-predicts δ (shear locking) but the
        // prediction must *increase toward* the analytic value as the
        // mesh is refined — that monotone approach is the real test
        // that the assembly + solve are correct.
        let (lx, ly, lz) = (10.0, 1.0, 1.0);
        let mat = FemMaterial {
            youngs_modulus: 200.0e9,
            poisson_ratio: 0.3,
            ..FemMaterial::default()
        };
        let load = 1.0e4;

        let solve_tip = |nx: usize, ny: usize, nz: usize| -> f64 {
            let mesh = structured_box_mesh(lx, ly, lz, nx, ny, nz).expect("valid box params");
            let mut constraints = Vec::new();
            for k in 0..=nz {
                for j in 0..=ny {
                    constraints.push(NodalConstraint::fixed(nid(0, j, k, nx, ny)));
                }
            }
            let tip_nodes: Vec<usize> = (0..=nz)
                .flat_map(|k| (0..=ny).map(move |j| (j, k)))
                .map(|(j, k)| nid(nx, j, k, nx, ny))
                .collect();
            let per = -load / tip_nodes.len() as f64;
            let forces: Vec<NodalForce> = tip_nodes
                .iter()
                .map(|&n| NodalForce {
                    node: n,
                    force: [0.0, 0.0, per],
                })
                .collect();
            let sol = solve_linear_static(&mesh, &mat, &constraints, &forces).unwrap();
            tip_nodes
                .iter()
                .map(|&n| sol.displacement[n][2].abs())
                .sum::<f64>()
                / tip_nodes.len() as f64
        };

        let coarse = solve_tip(4, 1, 1);
        let fine = solve_tip(16, 2, 2);
        // Analytic Euler-Bernoulli reference.
        let i_section = ly * lz.powi(3) / 12.0;
        let analytic = load * lx.powi(3) / (3.0 * mat.youngs_modulus * i_section);

        // Both are positive, finite, and the refined mesh is closer to
        // (and still below) the analytic beam-theory deflection.
        assert!(coarse > 0.0 && fine > 0.0);
        assert!(
            fine > coarse,
            "refinement should increase δ: {coarse} → {fine}"
        );
        assert!(
            fine < analytic * 1.2,
            "refined δ {fine} should not overshoot analytic {analytic}"
        );
        assert!(
            fine > analytic * 0.4,
            "refined δ {fine} should be within a factor of analytic {analytic}"
        );
    }

    #[test]
    fn von_mises_of_pure_hydrostatic_stress_is_zero() {
        // Equal normal stresses, no shear → von Mises = 0.
        let s = [5.0e6, 5.0e6, 5.0e6, 0.0, 0.0, 0.0];
        assert!(von_mises_from_voigt(&s) < 1e-3);
        // Uniaxial stress σ → von Mises = σ.
        let u = [3.0e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        assert!((von_mises_from_voigt(&u) - 3.0e6).abs() < 1.0);
    }

    #[test]
    fn solution_helpers_report_extrema() {
        let sol = NativeSolution {
            displacement: vec![[3.0, 4.0, 0.0], [0.0, 0.0, 0.0]],
            von_mises: vec![10.0, 7.0],
            stress: vec![[0.0; 6]; 2],
        };
        assert!((sol.max_displacement() - 5.0).abs() < 1e-9);
        assert!((sol.max_von_mises() - 10.0).abs() < 1e-9);
    }

    // ----- Round-1 H1–H3: dense-allocation cap --------------------------

    #[test]
    fn check_dense_dofs_accepts_a_small_model() {
        // A handful of nodes is comfortably under the cap and returns the
        // computed DOF count `3·n_nodes` untouched.
        assert_eq!(check_dense_dofs(8).unwrap(), 24);
        assert_eq!(check_dense_dofs(1).unwrap(), 3);
        // Exactly at the cap is still allowed.
        let at = MAX_DENSE_DOFS / 3;
        assert!(check_dense_dofs(at).is_ok());
    }

    #[test]
    fn check_dense_dofs_rejects_an_oversized_model() {
        // A node count whose `3·n_nodes` exceeds MAX_DENSE_DOFS must be
        // refused WITHOUT allocating any matrix — the helper does pure
        // arithmetic, so this is instant even for a huge `n_nodes`.
        let too_many = MAX_DENSE_DOFS / 3 + 1;
        let err = check_dense_dofs(too_many).unwrap_err();
        match err {
            NativeSolverError::TooLarge { dofs, max } => {
                assert_eq!(dofs, too_many * 3);
                assert_eq!(max, MAX_DENSE_DOFS);
            }
            other => panic!("expected TooLarge, got {other:?}"),
        }
    }

    #[test]
    fn check_dense_dofs_rejects_overflowing_node_count() {
        // `3·n_nodes` overflows usize → a structured error, never a wrap.
        let err = check_dense_dofs(usize::MAX).unwrap_err();
        assert!(matches!(err, NativeSolverError::TooLarge { .. }));
        // A node count that does not overflow `3·n` but whose `n_dof²`
        // would overflow usize is likewise rejected.
        let big = (MAX_DENSE_DOFS / 3).max(usize::MAX / 4);
        assert!(check_dense_dofs(big).is_err());
    }

    // ----- Round-2 F1: dense cap on an already-computed DOF count -------

    #[test]
    fn check_dense_dof_count_accepts_at_cap_rejects_above() {
        // The already-computed-count entry point returns its argument
        // unchanged up to and including the cap, then rejects.
        assert_eq!(check_dense_dof_count(0).unwrap(), 0);
        assert_eq!(check_dense_dof_count(48).unwrap(), 48);
        assert_eq!(
            check_dense_dof_count(MAX_DENSE_DOFS).unwrap(),
            MAX_DENSE_DOFS
        );
        let err = check_dense_dof_count(MAX_DENSE_DOFS + 1).unwrap_err();
        match err {
            NativeSolverError::TooLarge { dofs, max } => {
                assert_eq!(dofs, MAX_DENSE_DOFS + 1);
                assert_eq!(max, MAX_DENSE_DOFS);
            }
            other => panic!("expected TooLarge, got {other:?}"),
        }
    }

    #[test]
    fn check_dense_dof_count_catches_beams_6dof_undercount() {
        // The regression motivating round-2: a node count whose 3·n is
        // under the cap but whose 6·n (beam: 6 DOF/node) is over it must
        // be rejected when the caller passes the *beam* DOF count. The
        // old 3·n-hardcoded check_dense_dofs would wave it through.
        let n_nodes = MAX_DENSE_DOFS / 6 + 1; // 6·n > cap, 3·n < cap
        assert!(
            check_dense_dofs(n_nodes).is_ok(),
            "3·n is still under the cap — proves the 2× undercount"
        );
        let n_dof = n_nodes * 6;
        let err = check_dense_dof_count(n_dof).unwrap_err();
        assert!(
            matches!(err, NativeSolverError::TooLarge { .. }),
            "beam's 6·n DOF count must be capped, got {err:?}"
        );
    }

    // ----- Round-1 H4: non-finite load / BC rejection ------------------

    #[test]
    fn solve_rejects_nan_force() {
        // A single NodalForce carrying NaN must be rejected up front —
        // not pushed into the RHS where CscCholesky::solve would
        // back-substitute it into a silently-NaN displacement returned
        // as Ok(..). The error names the cause.
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box params");
        let forces = vec![NodalForce {
            node: 7,
            force: [f64::NAN, 0.0, 0.0],
        }];
        let err = solve_linear_static(
            &mesh,
            &FemMaterial::default(),
            &[NodalConstraint::fixed(0)],
            &forces,
        )
        .unwrap_err();
        assert!(
            matches!(err, NativeSolverError::InvalidLoad { .. }),
            "expected InvalidLoad, got {err:?}"
        );
    }

    #[test]
    fn solve_rejects_infinite_force() {
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box params");
        let forces = vec![NodalForce {
            node: 7,
            force: [0.0, f64::INFINITY, 0.0],
        }];
        let err = solve_linear_static(
            &mesh,
            &FemMaterial::default(),
            &[NodalConstraint::fixed(0)],
            &forces,
        )
        .unwrap_err();
        assert!(matches!(err, NativeSolverError::InvalidLoad { .. }));
    }

    #[test]
    fn solve_rejects_non_finite_prescribed_displacement() {
        // A prescribed (Dirichlet) displacement value is folded into the
        // RHS as `penalty·value`; a NaN there corrupts the solve exactly
        // like a NaN force. Reject it too.
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box params");
        let constraints = vec![
            NodalConstraint::fixed(0),
            NodalConstraint {
                node: 7,
                fixed: [Some(f64::NAN), None, None],
            },
        ];
        let err =
            solve_linear_static(&mesh, &FemMaterial::default(), &constraints, &[]).unwrap_err();
        assert!(matches!(err, NativeSolverError::InvalidLoad { .. }));
    }

    #[test]
    fn solve_mixed_rejects_nan_force() {
        // The mixed-element static path validates loads too.
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box params");
        let forces = vec![NodalForce {
            node: 7,
            force: [0.0, 0.0, f64::NAN],
        }];
        let err = solve_linear_static_mixed(
            &mesh,
            &FemMaterial::default(),
            &[NodalConstraint::fixed(0)],
            &forces,
        )
        .unwrap_err();
        assert!(matches!(err, NativeSolverError::InvalidLoad { .. }));
    }

    #[test]
    fn solve_still_accepts_finite_loads() {
        // The validation must not reject a normal finite load — the
        // small model still solves with a finite result. Clamp the x=0
        // face (a single pinned node would be rigid-body singular) and
        // pull the x=L face in +X.
        let (nx, ny, nz) = (1, 1, 1);
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, nx, ny, nz).expect("valid box params");
        let mut constraints = Vec::new();
        for k in 0..=nz {
            for j in 0..=ny {
                constraints.push(NodalConstraint::fixed(nid(0, j, k, nx, ny)));
            }
        }
        let forces = vec![NodalForce {
            node: nid(nx, 0, 0, nx, ny),
            force: [1.0e5, 0.0, 0.0],
        }];
        let sol =
            solve_linear_static(&mesh, &FemMaterial::default(), &constraints, &forces).unwrap();
        assert!(sol.max_displacement().is_finite());
        assert!(sol
            .displacement
            .iter()
            .all(|d| d.iter().all(|c| c.is_finite())));
        assert!(
            sol.max_displacement() > 0.0,
            "a real load should move the body"
        );
    }

    #[test]
    fn assemble_global_stiffness_rejects_oversized_mesh() {
        // The dense assembler must check the cap BEFORE the
        // `DMatrix::zeros(n_dof, n_dof)` allocation. We can't build a
        // real mesh that large, so we forge a mesh whose node array is
        // huge but whose element block is empty-ish — the cap fires on
        // `mesh.nodes.len()` before any element work or allocation.
        //
        // Building 16k+ real nodes is cheap and well-defined, so use a
        // node count one past the cap with a single (valid) Tet4 element
        // referencing the first four nodes.
        let n_nodes = MAX_DENSE_DOFS / 3 + 1;
        let mut mesh = Mesh::new("oversized");
        mesh.nodes = vec![Vector3::zeros(); n_nodes];
        mesh.element_blocks.push(ElementBlock {
            element_type: ElementType::Tet4,
            connectivity: vec![0, 1, 2, 3],
        });
        let d = elasticity_matrix(&FemMaterial::default()).unwrap();
        let err = crate::assembly::assemble_global_stiffness(&mesh, &d).unwrap_err();
        assert!(
            matches!(err, NativeSolverError::TooLarge { .. }),
            "expected TooLarge, got {err:?}"
        );
    }
}
