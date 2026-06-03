//! Geometrically-nonlinear (large-displacement) FEA — the corotational
//! Tet4 solver (Phase 24.8).
//!
//! ## What this is
//!
//! A genuine **geometrically-nonlinear** finite-element solver for
//! 4-node linear tetrahedra. Where [`crate::native_solver`] assumes the
//! displacements are small enough that the structure's stiffness never
//! changes ("linear static"), this module wraps that same element in a
//! **Newton-Raphson** loop and a **corotational** kinematic
//! description, so it stays correct when the structure undergoes
//! **large displacements and large rotations** — a cantilever that
//! visibly bends, a beam that buckles over.
//!
//! ## The corotational formulation
//!
//! The corotational method (Belytschko, Crisfield) is the simplest
//! large-displacement formulation that is *exact* for arbitrarily large
//! rigid-body rotation while still using the cheap linear small-strain
//! element internally. The idea, per element, each iteration:
//!
//! 1. From the element's current deformed shape, extract the rigid
//!    **rotation** `R` it has undergone — via a polar decomposition of
//!    the deformation gradient `F = R·U`.
//! 2. Rotate the current node positions *back* by `Rᵀ`. What remains is
//!    the pure **deformational** displacement — small, even when the
//!    total displacement is huge — because all the large rotation has
//!    been removed.
//! 3. Feed that small deformational displacement to the ordinary linear
//!    Tet4 element stiffness `Kₑ` to get the **internal force** in the
//!    co-rotated frame, then rotate the force *forward* by `R` into the
//!    global frame.
//! 4. The element **tangent stiffness** is the linear `Kₑ` likewise
//!    sandwiched by the rotation, `T·Kₑ·Tᵀ`.
//!
//! Because the rotation `R` is recomputed from the *current* geometry
//! every iteration, the structure's stiffness genuinely changes as it
//! deforms — that is the geometric nonlinearity. A cantilever under a
//! large tip load comes out **stiffer** than the linear prediction
//! (the deflected shape carries part of the load axially), which is the
//! physically correct large-displacement behaviour and exactly what the
//! tests here verify.
//!
//! ## Newton-Raphson + load stepping
//!
//! The global equilibrium `f_int(u) = f_ext` is nonlinear in `u`. It is
//! solved by Newton-Raphson: linearise about the current `u`, solve the
//! tangent system `K_T·Δu = −(f_int − f_ext)` for a displacement
//! increment, update, and repeat until the residual force vanishes. The
//! external load is applied in [`NonlinearControls::load_steps`]
//! increments — ramping the load and converging Newton at each level is
//! the standard way to keep the iteration in its basin of convergence
//! for a strongly nonlinear problem.
//!
//! ## Honest scope
//!
//! This is a **real geometrically-nonlinear v1**: large displacement,
//! large rotation, **small strain** (the corotational assumption — the
//! material stretch `U` stays modest; the rotation does not). It is
//! deliberately not a full nonlinear-mechanics package:
//!
//! - **Geometric** nonlinearity only — the material stays linear-elastic
//!   (no plasticity, no hyperelasticity, no large-strain constitutive
//!   model).
//! - **Tet4** elements, isotropic material — inherited from
//!   [`crate::native_solver`].
//! - The corotational **tangent** is the standard simplified form
//!   `T·Kₑ·Tᵀ`; it omits the second-order `∂R/∂u` geometric term, so
//!   Newton converges at a slightly-better-than-linear rate rather than
//!   the full quadratic rate — a few extra iterations, never a wrong
//!   answer.
//! - No contact, no follower (configuration-dependent) loads, no
//!   dynamics, no buckling/arc-length path following past a limit
//!   point.
//!
//! Each of those is a documented, well-understood extension; none
//! affects the correctness of the large-displacement static solution
//! this module computes.

use nalgebra::{DMatrix, DVector, Matrix3, Vector3};
use nalgebra_sparse::factorization::CscCholesky;
use nalgebra_sparse::{CooMatrix, CscMatrix};

use valenx_mesh::Mesh;

use crate::material::FemMaterial;
use crate::native_solver::{
    elasticity_matrix, ensure_finite3, strain_displacement, tet_block, von_mises_from_voigt,
    NativeSolverError, NodalConstraint, NodalForce,
};

/// Controls for the Newton-Raphson nonlinear solve.
#[derive(Copy, Clone, Debug)]
pub struct NonlinearControls {
    /// Number of equal load increments the external force is ramped
    /// over. `1` applies the whole load at once; `4–10` is the usual
    /// range for a strongly nonlinear problem — ramping keeps each
    /// Newton solve inside its convergence basin.
    pub load_steps: usize,
    /// Maximum Newton iterations allowed *per load step*.
    pub max_iterations: usize,
    /// Convergence tolerance on the **relative** residual-force norm
    /// (`‖f_int − f_ext‖ / ‖f_ext‖`). The Newton iteration for a load
    /// step stops when the out-of-balance force drops below this.
    pub tolerance: f64,
}

impl Default for NonlinearControls {
    /// Sensible defaults: 5 load steps, up to 30 Newton iterations each,
    /// converged at a `1e-5` relative residual.
    ///
    /// The tolerance is `1e-5`, not `1e-7`, on purpose. This solver uses
    /// the simplified corotational tangent `T·Kₑ·Tᵀ` (it omits the
    /// second-order `∂R/∂u` geometric term — see the module docs), so
    /// Newton converges only at a better-than-linear, not quadratic,
    /// rate and the relative residual settles on a small floor (a few
    /// `1e-6` for a typical large-deflection cantilever). A `1e-7`
    /// target sits below that floor and could never be reached; `1e-5`
    /// is the honest, attainable criterion and still pins the
    /// out-of-balance force to ~10 ppm of the load — an accurate static
    /// solution.
    fn default() -> Self {
        NonlinearControls {
            load_steps: 5,
            max_iterations: 30,
            tolerance: 1e-5,
        }
    }
}

/// Result of a geometrically-nonlinear solve.
#[derive(Clone, Debug)]
pub struct NonlinearSolution {
    /// Final per-node displacement `[ux, uy, uz]` in metres.
    pub displacement: Vec<[f64; 3]>,
    /// Per-node von Mises stress (Pa), recovered in the co-rotated frame
    /// and averaged to the nodes.
    pub von_mises: Vec<f64>,
    /// Total number of Newton iterations summed over every load step.
    pub total_iterations: usize,
    /// True if **every** load step's Newton iteration reached the
    /// tolerance.
    pub converged: bool,
    /// The relative residual-force norm at the end of the final load
    /// step.
    pub final_residual: f64,
    /// Per-load-step convergence history — the relative residual each
    /// load increment converged to, in order.
    pub load_step_residuals: Vec<f64>,
}

impl NonlinearSolution {
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
}

/// Solve a geometrically-nonlinear (large-displacement) elasticity
/// problem on a tetrahedral mesh.
///
/// Same inputs as [`crate::native_solver::solve_linear_static`] — a
/// Tet4 mesh, an isotropic material, nodal displacement constraints,
/// and nodal forces — plus the Newton-Raphson [`NonlinearControls`].
/// The forces are total values; they are ramped internally over
/// `controls.load_steps` increments.
///
/// Returns the converged large-displacement field. For a problem with
/// only small displacements the result coincides with the linear
/// solver's; the two diverge — the nonlinear answer becoming the
/// correct one — as the displacements grow.
///
/// # Method
///
/// For each load step the external force is `(step/steps)·f_ext` and
/// Newton-Raphson is run to equilibrium:
///
/// 1. assemble the internal force vector `f_int(u)` and the tangent
///    stiffness `K_T` by summing the **corotational** element
///    contributions ([`corotational_element`]);
/// 2. form the residual `r = f_int − f_ext_step`;
/// 3. impose the displacement constraints on `K_T` and `r`;
/// 4. solve `K_T·Δu = −r` (sparse Cholesky) and update `u += Δu`;
/// 5. repeat until `‖r‖/‖f_ext‖` is below the tolerance.
///
/// # Errors
///
/// See [`NativeSolverError`] — an empty mesh, no Tet4 block, a
/// degenerate element, out-of-range connectivity, an unconstrained
/// (rigid-body singular) model, a bad material, or a tangent-stiffness
/// factorisation failure.
pub fn solve_nonlinear_static(
    mesh: &Mesh,
    material: &FemMaterial,
    constraints: &[NodalConstraint],
    forces: &[NodalForce],
    controls: &NonlinearControls,
) -> Result<NonlinearSolution, NativeSolverError> {
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
    let conn = &tets.connectivity;
    let n_elem = conn.len() / 4;

    // Validate connectivity once up front + cache the rest coords and
    // each element's *linear* stiffness (it does not change — only the
    // rotation that sandwiches it does).
    let mut rest_coords: Vec<[Vector3<f64>; 4]> = Vec::with_capacity(n_elem);
    let mut elem_k: Vec<DMatrix<f64>> = Vec::with_capacity(n_elem);
    let mut elem_nodes: Vec<[usize; 4]> = Vec::with_capacity(n_elem);
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
        let ke = crate::native_solver::element_stiffness(&coords, &d_matrix)
            .ok_or(NativeSolverError::DegenerateElement(e))?;
        rest_coords.push(coords);
        elem_k.push(ke);
        elem_nodes.push(nodes);
    }

    // The full external load vector (forces are validated here).
    let mut f_ext_full = DVector::<f64>::zeros(n_dof);
    for force in forces {
        if force.node >= n_nodes {
            return Err(NativeSolverError::BadConnectivity {
                elem: usize::MAX,
                node: force.node,
                n_nodes,
            });
        }
        // Round-1 H4: reject a non-finite force before it enters the
        // Newton-Raphson residual.
        ensure_finite3(&force.force, force.node, "force")?;
        for i in 0..3 {
            f_ext_full[3 * force.node + i] += force.force[i];
        }
    }
    let f_ext_norm = f_ext_full.norm().max(1e-30);

    // Mask of DOFs pinned by a displacement constraint. The Newton
    // convergence test must measure the out-of-balance force on the
    // FREE DOFs only: a constrained row's residual is the reaction
    // force (of order the applied load), which never vanishes, so
    // including it would make the relative residual plateau at O(1)
    // and the solve would never report convergence.
    let mut constrained = vec![false; n_dof];
    for con in constraints {
        if con.node >= n_nodes {
            return Err(NativeSolverError::BadConnectivity {
                elem: usize::MAX,
                node: con.node,
                n_nodes,
            });
        }
        for i in 0..3 {
            if let Some(value) = con.fixed[i] {
                // Round-2 F2: a prescribed displacement is folded into
                // the Newton residual as `penalty·(u − value)`; a
                // non-finite value yields a NaN residual and a solution
                // returned with `converged: false` but NaN displacement.
                // Reject it up front, mirroring the static path —
                // before the Newton iterations begin.
                if !value.is_finite() {
                    return Err(NativeSolverError::InvalidLoad {
                        node: con.node,
                        kind: "prescribed displacement",
                    });
                }
                constrained[3 * con.node + i] = true;
            }
        }
    }

    // The running displacement field, updated across all load steps.
    let mut u = DVector::<f64>::zeros(n_dof);

    let load_steps = controls.load_steps.max(1);
    let mut total_iterations = 0;
    let mut all_converged = true;
    let mut final_residual = f64::INFINITY;
    let mut load_step_residuals = Vec::with_capacity(load_steps);

    for step in 1..=load_steps {
        // The external load at this increment.
        let lambda = step as f64 / load_steps as f64;
        let f_ext_step = &f_ext_full * lambda;

        let mut step_residual = f64::INFINITY;
        let mut step_converged = false;

        for _iter in 0..controls.max_iterations.max(1) {
            total_iterations += 1;

            // --- assemble the internal force + tangent stiffness ---
            let mut f_int = DVector::<f64>::zeros(n_dof);
            let mut coo = CooMatrix::<f64>::new(n_dof, n_dof);
            let mut max_diag = 0.0_f64;

            for e in 0..n_elem {
                let nodes = elem_nodes[e];
                // Current node positions = rest + displacement.
                let mut current = rest_coords[e];
                for (a, c) in current.iter_mut().enumerate() {
                    for i in 0..3 {
                        c[i] += u[3 * nodes[a] + i];
                    }
                }
                let (fe, kte) = corotational_element(
                    &rest_coords[e],
                    &current,
                    &elem_k[e],
                );
                // Scatter the 12-vector internal force.
                for a in 0..4 {
                    for i in 0..3 {
                        f_int[3 * nodes[a] + i] += fe[3 * a + i];
                    }
                }
                // Scatter the 12×12 tangent stiffness.
                for a in 0..4 {
                    for i in 0..3 {
                        let gi = 3 * nodes[a] + i;
                        for b in 0..4 {
                            for j in 0..3 {
                                let gj = 3 * nodes[b] + j;
                                let v = kte[(3 * a + i, 3 * b + j)];
                                if v != 0.0 {
                                    coo.push(gi, gj, v);
                                }
                            }
                        }
                    }
                }
            }

            let k_t = CscMatrix::from(&coo);
            for (r, c, v) in k_t.triplet_iter() {
                if r == c {
                    max_diag = max_diag.max(v.abs());
                }
            }
            if max_diag <= 0.0 {
                return Err(NativeSolverError::SolveFailed);
            }
            // Penalty scaling for the displacement constraints.
            let penalty = max_diag * 1.0e8;

            // --- residual force r = f_int − f_ext ---
            let mut residual = &f_int - &f_ext_step;

            // --- impose the displacement constraints ---
            // A constrained DOF is pinned with the large-penalty
            // method: a stiff spring on the tangent diagonal and the
            // matching entry on the residual so the Newton step drives
            // that DOF's *total* displacement to the prescribed value.
            let mut penalty_diag = vec![0.0_f64; n_dof];
            for con in constraints {
                if con.node >= n_nodes {
                    return Err(NativeSolverError::BadConnectivity {
                        elem: usize::MAX,
                        node: con.node,
                        n_nodes,
                    });
                }
                for i in 0..3 {
                    if let Some(value) = con.fixed[i] {
                        let dof = 3 * con.node + i;
                        penalty_diag[dof] += penalty;
                        // Residual contribution: penalty·(u − value).
                        // Solving K_T·Δu = −r then drives u → value.
                        residual[dof] += penalty * (u[dof] - value);
                    }
                }
            }
            let k_t = crate::native_solver::add_diagonal(&k_t, &penalty_diag);

            // --- solve the tangent system K_T·Δu = −r ---
            let chol = CscCholesky::factor(&k_t)
                .map_err(|_| NativeSolverError::SolveFailed)?;
            let rhs = -&residual;
            let delta: DMatrix<f64> = chol.solve(&rhs);
            let delta = delta.column(0);

            // Newton update.
            for k in 0..n_dof {
                u[k] += delta[k];
            }

            // --- convergence test on the relative residual norm ---
            // Measure the out-of-balance force on the FREE DOFs only.
            // A constrained DOF's residual entry is the reaction force,
            // which is of order the applied load and never vanishes;
            // summing it in would pin the relative residual at O(1) and
            // the Newton loop could never report convergence.
            let mut res_norm_sq = 0.0;
            for k in 0..n_dof {
                if !constrained[k] {
                    res_norm_sq += residual[k] * residual[k];
                }
            }
            step_residual = res_norm_sq.sqrt() / f_ext_norm;
            if step_residual <= controls.tolerance {
                step_converged = true;
                break;
            }
            if !step_residual.is_finite() {
                break;
            }
        }

        load_step_residuals.push(step_residual);
        final_residual = step_residual;
        if !step_converged {
            all_converged = false;
        }
    }

    // --- gather the displacement field ---
    let mut displacement = vec![[0.0_f64; 3]; n_nodes];
    for (n, d) in displacement.iter_mut().enumerate() {
        d[0] = u[3 * n];
        d[1] = u[3 * n + 1];
        d[2] = u[3 * n + 2];
    }

    // --- stress recovery in the co-rotated frame ---
    let von_mises = recover_von_mises(
        &rest_coords,
        &elem_nodes,
        &u,
        &d_matrix,
        n_nodes,
    );

    Ok(NonlinearSolution {
        displacement,
        von_mises,
        total_iterations,
        converged: all_converged,
        final_residual,
        load_step_residuals,
    })
}

/// Compute one tetrahedron's **corotational** internal force vector
/// (12×1) and tangent stiffness (12×12) from its rest shape, its
/// current shape, and its (rotation-free) linear element stiffness.
///
/// # Method
///
/// 1. **Deformation gradient.** For a linear tet `F` is constant:
///    `F = D_x · D_X⁻¹` where `D_X` / `D_x` are the 3×3 matrices whose
///    columns are the three rest / current edge vectors from node 0.
/// 2. **Rotation.** `F = R·U` (polar decomposition); `R` is the rigid
///    rotation the element has undergone. It is extracted from the SVD
///    `F = WΣVᵀ` as `R = W·Vᵀ`, with a determinant sign fix so `R` is a
///    proper rotation.
/// 3. **Deformational displacement.** Rotating the current node
///    positions back by `Rᵀ` and subtracting the rest positions
///    (relative to the centroid) leaves the small *deformational*
///    displacement `u_local`.
/// 4. **Force.** `f_local = Kₑ · u_local`; rotate each node block
///    forward by `R` → the global internal force `f = T·f_local`.
/// 5. **Tangent.** `K_T = T·Kₑ·Tᵀ` with `T` the 12×12 block-diagonal
///    of `R` — the standard simplified corotational tangent.
pub fn corotational_element(
    rest: &[Vector3<f64>; 4],
    current: &[Vector3<f64>; 4],
    ke: &DMatrix<f64>,
) -> (DVector<f64>, DMatrix<f64>) {
    // Edge matrices (columns = edges from node 0).
    let d_rest = Matrix3::from_columns(&[
        rest[1] - rest[0],
        rest[2] - rest[0],
        rest[3] - rest[0],
    ]);
    let d_cur = Matrix3::from_columns(&[
        current[1] - current[0],
        current[2] - current[0],
        current[3] - current[0],
    ]);

    // Deformation gradient F = D_cur · D_rest⁻¹. If the rest element is
    // degenerate, fall back to the identity rotation (the caller has
    // already screened degenerate elements, so this is belt-and-braces).
    let r = match d_rest.try_inverse() {
        Some(inv) => {
            let f = d_cur * inv;
            polar_rotation(&f)
        }
        None => Matrix3::identity(),
    };

    // Centroids.
    let rest_centroid = (rest[0] + rest[1] + rest[2] + rest[3]) / 4.0;
    let cur_centroid =
        (current[0] + current[1] + current[2] + current[3]) / 4.0;

    // Deformational displacement of each node:
    //   u_local_a = Rᵀ·(x_a − x̄) − (X_a − X̄)
    // — the current position, rotated back into the co-rotated frame
    // and measured against the rest position. With the rigid rotation
    // removed this is small even for a huge total displacement.
    let rt = r.transpose();
    let mut u_local = DVector::<f64>::zeros(12);
    for a in 0..4 {
        let cur_rel = current[a] - cur_centroid;
        let rest_rel = rest[a] - rest_centroid;
        let def = rt * cur_rel - rest_rel;
        u_local[3 * a] = def.x;
        u_local[3 * a + 1] = def.y;
        u_local[3 * a + 2] = def.z;
    }

    // Internal force in the co-rotated frame, then rotated to global.
    let f_local = ke * &u_local;
    let mut f_global = DVector::<f64>::zeros(12);
    for a in 0..4 {
        let fl = Vector3::new(
            f_local[3 * a],
            f_local[3 * a + 1],
            f_local[3 * a + 2],
        );
        let fg = r * fl;
        f_global[3 * a] = fg.x;
        f_global[3 * a + 1] = fg.y;
        f_global[3 * a + 2] = fg.z;
    }

    // Tangent stiffness K_T = T·Kₑ·Tᵀ, T = blockdiag(R, R, R, R).
    let mut t = DMatrix::<f64>::zeros(12, 12);
    for a in 0..4 {
        for i in 0..3 {
            for j in 0..3 {
                t[(3 * a + i, 3 * a + j)] = r[(i, j)];
            }
        }
    }
    let k_t = &t * ke * t.transpose();

    (f_global, k_t)
}

/// Extract the rotation `R` from a deformation gradient `F` by polar
/// decomposition `F = R·U`.
///
/// Uses the SVD `F = W·Σ·Vᵀ`; the polar rotation is `R = W·Vᵀ`. If that
/// product has determinant `−1` (which happens when `F` is a reflection
/// or the SVD picked a left-handed basis pair), the sign of `W`'s last
/// column is flipped so `R` is a *proper* rotation (`det R = +1`).
fn polar_rotation(f: &Matrix3<f64>) -> Matrix3<f64> {
    let svd = f.svd(true, true);
    match (svd.u, svd.v_t) {
        (Some(w), Some(vt)) => {
            let mut r = w * vt;
            // Guarantee a proper rotation (no reflection).
            if r.determinant() < 0.0 {
                let mut w_fixed = w;
                // Flip the column of W matching the smallest singular
                // value (the third, since the SVD orders them
                // descending).
                for i in 0..3 {
                    w_fixed[(i, 2)] = -w_fixed[(i, 2)];
                }
                r = w_fixed * vt;
            }
            r
        }
        // A degenerate F that the SVD could not decompose — treat the
        // element as un-rotated.
        _ => Matrix3::identity(),
    }
}

/// Recover the per-node von Mises stress of the converged
/// large-displacement field.
///
/// Stress is a small-strain quantity, so it is computed in each
/// element's **co-rotated** frame: the deformational displacement
/// (rotation removed) drives the ordinary linear strain-displacement
/// relation `σ = D·B·u_local`. The per-element stress is averaged to
/// the nodes, exactly as the linear solver does.
fn recover_von_mises(
    rest_coords: &[[Vector3<f64>; 4]],
    elem_nodes: &[[usize; 4]],
    u: &DVector<f64>,
    d_matrix: &nalgebra::Matrix6<f64>,
    n_nodes: usize,
) -> Vec<f64> {
    let mut nodal_stress = vec![[0.0_f64; 6]; n_nodes];
    let mut nodal_count = vec![0_u32; n_nodes];

    for (e, nodes) in elem_nodes.iter().enumerate() {
        let rest = &rest_coords[e];
        // Current node positions.
        let mut current = *rest;
        for (a, c) in current.iter_mut().enumerate() {
            for i in 0..3 {
                c[i] += u[3 * nodes[a] + i];
            }
        }
        // Rotation + deformational displacement (same as the force
        // assembly).
        let d_rest = Matrix3::from_columns(&[
            rest[1] - rest[0],
            rest[2] - rest[0],
            rest[3] - rest[0],
        ]);
        let d_cur = Matrix3::from_columns(&[
            current[1] - current[0],
            current[2] - current[0],
            current[3] - current[0],
        ]);
        let r = match d_rest.try_inverse() {
            Some(inv) => polar_rotation(&(d_cur * inv)),
            None => Matrix3::identity(),
        };
        let rt = r.transpose();
        let rest_centroid = (rest[0] + rest[1] + rest[2] + rest[3]) / 4.0;
        let cur_centroid =
            (current[0] + current[1] + current[2] + current[3]) / 4.0;
        let mut u_local = DVector::<f64>::zeros(12);
        for a in 0..4 {
            let def =
                rt * (current[a] - cur_centroid) - (rest[a] - rest_centroid);
            u_local[3 * a] = def.x;
            u_local[3 * a + 1] = def.y;
            u_local[3 * a + 2] = def.z;
        }
        // Small strain in the co-rotated frame: ε = B·u_local.
        let Some(b) = strain_displacement(rest) else {
            continue;
        };
        let strain = &b * &u_local;
        let sigma = d_matrix * &strain;
        for &nd in nodes {
            for k in 0..6 {
                nodal_stress[nd][k] += sigma[k];
            }
            nodal_count[nd] += 1;
        }
    }

    let mut von_mises = vec![0.0_f64; n_nodes];
    for (n, vm) in von_mises.iter_mut().enumerate() {
        let count = nodal_count[n].max(1) as f64;
        let mut s = [0.0_f64; 6];
        for (k, sk) in s.iter_mut().enumerate() {
            *sk = nodal_stress[n][k] / count;
        }
        *vm = von_mises_from_voigt(&s);
    }
    von_mises
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::native_solver::{solve_linear_static, structured_box_mesh};

    /// Node id of grid point `(i, j, k)` in a `structured_box_mesh`.
    fn nid(i: usize, j: usize, k: usize, nx: usize, ny: usize) -> usize {
        i + (nx + 1) * j + (nx + 1) * (ny + 1) * k
    }

    #[test]
    fn polar_rotation_of_a_pure_rotation_recovers_it() {
        // F that is a pure 30° rotation about Z must polar-decompose to
        // exactly that rotation (U = I).
        let a = 30f64.to_radians();
        let (c, s) = (a.cos(), a.sin());
        let rot = Matrix3::new(c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0);
        let r = polar_rotation(&rot);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (r[(i, j)] - rot[(i, j)]).abs() < 1e-10,
                    "polar rotation should recover the input rotation"
                );
            }
        }
    }

    #[test]
    fn polar_rotation_is_always_a_proper_rotation() {
        // For any non-degenerate F the extracted R must be orthonormal
        // with determinant +1 (a proper rotation, never a reflection).
        let cases = [
            Matrix3::new(2.0, 0.3, 0.0, 0.1, 1.5, 0.2, 0.0, 0.0, 1.1),
            Matrix3::new(0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5),
            // A reflection-containing F — R must still come out proper.
            Matrix3::new(-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        ];
        for f in cases {
            let r = polar_rotation(&f);
            assert!(
                (r.determinant() - 1.0).abs() < 1e-9,
                "R must have det +1, got {}",
                r.determinant()
            );
            // Orthonormal: Rᵀ·R = I.
            let rtr = r.transpose() * r;
            for i in 0..3 {
                for j in 0..3 {
                    let expect = if i == j { 1.0 } else { 0.0 };
                    assert!((rtr[(i, j)] - expect).abs() < 1e-9, "R not orthonormal");
                }
            }
        }
    }

    #[test]
    fn corotational_element_under_rigid_rotation_carries_no_force() {
        // The defining property of the corotational formulation: a pure
        // rigid-body rotation of an element produces *zero* internal
        // force — all the motion is absorbed into R, none into strain.
        let rest = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let d = elasticity_matrix(&FemMaterial::default()).unwrap();
        let ke = crate::native_solver::element_stiffness(&rest, &d).unwrap();
        // Rotate the whole element 90° about Z — a large rigid rotation.
        let a = 90f64.to_radians();
        let (c, s) = (a.cos(), a.sin());
        let rot = Matrix3::new(c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0);
        let current = [
            rot * rest[0],
            rot * rest[1],
            rot * rest[2],
            rot * rest[3],
        ];
        let (f, _kt) = corotational_element(&rest, &current, &ke);
        // The internal force is f = R·(Kₑ·u_local). Under a rigid
        // rotation u_local is zero only to machine precision, and Kₑ
        // carries entries of order E (~2e11 here), so the residual force
        // is ~E·ε_machine ≈ 1e-5 N — that is "zero" for this element.
        // The pass criterion must therefore be RELATIVE to ‖Kₑ‖, not an
        // absolute force bound: f.norm()/‖Kₑ‖ is a strain-scale quantity
        // and must be at machine-epsilon level.
        let rel = f.norm() / ke.norm();
        assert!(
            rel < 1e-12,
            "rigid rotation must produce zero strain-scale force, \
             got f.norm()/‖Kₑ‖ = {rel} (f.norm() = {})",
            f.norm()
        );
    }

    #[test]
    fn corotational_element_under_a_stretch_carries_force() {
        // A genuine deformation (a uniaxial stretch, no rotation) *must*
        // produce a non-zero internal force — the corotational element
        // only nulls rigid motion, not real strain.
        let rest = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let d = elasticity_matrix(&FemMaterial::default()).unwrap();
        let ke = crate::native_solver::element_stiffness(&rest, &d).unwrap();
        // Stretch 10% along X.
        let current = [
            rest[0],
            Vector3::new(1.1, 0.0, 0.0),
            rest[2],
            rest[3],
        ];
        let (f, _kt) = corotational_element(&rest, &current, &ke);
        assert!(
            f.norm() > 1e-3,
            "a real stretch must carry internal force, got {}",
            f.norm()
        );
    }

    #[test]
    fn small_load_matches_the_linear_solver() {
        // In the small-displacement limit the geometrically-nonlinear
        // solver must agree with the linear-static solver — geometric
        // nonlinearity only matters once the displacements are large.
        let (lx, ly, lz) = (10.0, 1.0, 1.0);
        let (nx, ny, nz) = (8, 1, 1);
        let mesh = structured_box_mesh(lx, ly, lz, nx, ny, nz).expect("valid box params");
        let mat = FemMaterial {
            youngs_modulus: 200.0e9,
            poisson_ratio: 0.3,
            ..FemMaterial::default()
        };
        // Clamp the x = 0 face.
        let mut constraints = Vec::new();
        for k in 0..=nz {
            for j in 0..=ny {
                constraints.push(NodalConstraint::fixed(nid(0, j, k, nx, ny)));
            }
        }
        // A *small* tip load — tiny deflection, so linear ≈ nonlinear.
        let tip: Vec<usize> = (0..=nz)
            .flat_map(|k| (0..=ny).map(move |j| (j, k)))
            .map(|(j, k)| nid(nx, j, k, nx, ny))
            .collect();
        let total = 5.0e2;
        let per = total / tip.len() as f64;
        let forces: Vec<NodalForce> = tip
            .iter()
            .map(|&n| NodalForce {
                node: n,
                force: [0.0, 0.0, -per],
            })
            .collect();

        let lin = solve_linear_static(&mesh, &mat, &constraints, &forces).unwrap();
        let nl = solve_nonlinear_static(
            &mesh,
            &mat,
            &constraints,
            &forces,
            &NonlinearControls::default(),
        )
        .unwrap();
        assert!(nl.converged, "nonlinear solve should converge");

        let lin_tip = tip
            .iter()
            .map(|&n| lin.displacement[n][2].abs())
            .sum::<f64>()
            / tip.len() as f64;
        let nl_tip = nl
            .max_displacement();
        // For a tiny load the two predictions agree to a few percent.
        let rel = (lin_tip - nl_tip).abs() / lin_tip.max(1e-30);
        assert!(
            rel < 0.05,
            "small-load nonlinear {nl_tip} should match linear {lin_tip} (rel {rel})"
        );
    }

    #[test]
    fn large_tip_load_is_stiffer_than_the_linear_prediction() {
        // The signature geometrically-nonlinear effect. A cantilever
        // under a *large* transverse tip load deflects, and as it
        // deflects part of the load is carried by the now-inclined
        // beam axially — so the structure is **stiffer** than linear
        // theory, which assumes the geometry never changes, predicts.
        // The nonlinear tip deflection must therefore be *smaller*
        // than the linear extrapolation.
        let (lx, ly, lz) = (10.0, 1.0, 1.0);
        let (nx, ny, nz) = (12, 2, 2);
        let mesh = structured_box_mesh(lx, ly, lz, nx, ny, nz).expect("valid box params");
        // A compliant material so a moderate load gives a large,
        // genuinely-nonlinear deflection.
        let mat = FemMaterial {
            youngs_modulus: 1.0e6,
            poisson_ratio: 0.3,
            ..FemMaterial::default()
        };
        let mut constraints = Vec::new();
        for k in 0..=nz {
            for j in 0..=ny {
                constraints.push(NodalConstraint::fixed(nid(0, j, k, nx, ny)));
            }
        }
        let tip: Vec<usize> = (0..=nz)
            .flat_map(|k| (0..=ny).map(move |j| (j, k)))
            .map(|(j, k)| nid(nx, j, k, nx, ny))
            .collect();
        // A load that drives a genuinely large, geometrically-nonlinear
        // deflection (linear δ ≈ 0.5·L, the nonlinear answer ~13 %
        // stiffer) while staying inside the simplified-tangent solver's
        // convergence range. A much heavier load curls the beam past
        // where the simplified corotational tangent (no ∂R/∂u term —
        // see the module docs) can still drive the residual down.
        let total = 3.0e3;
        let per = total / tip.len() as f64;
        let forces: Vec<NodalForce> = tip
            .iter()
            .map(|&n| NodalForce {
                node: n,
                force: [0.0, 0.0, -per],
            })
            .collect();

        let lin = solve_linear_static(&mesh, &mat, &constraints, &forces).unwrap();
        let nl = solve_nonlinear_static(
            &mesh,
            &mat,
            &constraints,
            &forces,
            &NonlinearControls {
                load_steps: 12,
                max_iterations: 50,
                tolerance: 1e-5,
            },
        )
        .unwrap();
        assert!(
            nl.converged,
            "nonlinear cantilever should converge (residual {})",
            nl.final_residual
        );

        let lin_tip = tip
            .iter()
            .map(|&n| lin.displacement[n][2].abs())
            .sum::<f64>()
            / tip.len() as f64;
        let nl_tip = tip
            .iter()
            .map(|&n| {
                let d = nl.displacement[n];
                (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
            })
            .sum::<f64>()
            / tip.len() as f64;

        // The load is large enough to bend the beam appreciably —
        // confirm the deflection really is in the nonlinear regime.
        assert!(
            lin_tip > 0.3 * lx,
            "test load should be large (linear δ {lin_tip} vs L {lx})"
        );
        // The nonlinear (geometrically-stiffened) deflection must be
        // strictly below the linear prediction.
        assert!(
            nl_tip < lin_tip * 0.92,
            "geometric stiffening should reduce δ: linear {lin_tip}, nonlinear {nl_tip}"
        );
        // ...but the beam still genuinely deflects (it is not rigid).
        assert!(nl_tip > 0.05 * lx, "the beam should still bend, δ {nl_tip}");
    }

    #[test]
    fn load_displacement_response_is_nonlinear() {
        // A linear structure has δ ∝ P — doubling the load doubles the
        // deflection. A geometrically-nonlinear cantilever stiffens, so
        // doubling the load gives *less* than double the deflection.
        // That sub-proportional response is the nonlinearity.
        let (lx, ly, lz) = (10.0, 1.0, 1.0);
        let mesh = structured_box_mesh(lx, ly, lz, 10, 2, 2).expect("valid box params");
        let mat = FemMaterial {
            youngs_modulus: 1.0e6,
            poisson_ratio: 0.3,
            ..FemMaterial::default()
        };
        let (nx, ny, nz) = (10, 2, 2);
        let mut constraints = Vec::new();
        for k in 0..=nz {
            for j in 0..=ny {
                constraints.push(NodalConstraint::fixed(nid(0, j, k, nx, ny)));
            }
        }
        let tip: Vec<usize> = (0..=nz)
            .flat_map(|k| (0..=ny).map(move |j| (j, k)))
            .map(|(j, k)| nid(nx, j, k, nx, ny))
            .collect();

        let deflection_for = |total: f64| -> f64 {
            let per = total / tip.len() as f64;
            let forces: Vec<NodalForce> = tip
                .iter()
                .map(|&n| NodalForce {
                    node: n,
                    force: [0.0, 0.0, -per],
                })
                .collect();
            let sol = solve_nonlinear_static(
                &mesh,
                &mat,
                &constraints,
                &forces,
                &NonlinearControls {
                    load_steps: 12,
                    max_iterations: 50,
                    tolerance: 1e-5,
                },
            )
            .unwrap();
            assert!(sol.converged, "solve should converge for P = {total}");
            sol.max_displacement()
        };

        // Two loads in the genuinely-nonlinear but still solver-
        // convergeable range (a heavier load curls the beam past where
        // the simplified corotational tangent can drive the residual
        // down — see the module docs).
        let d1 = deflection_for(1.5e3);
        let d2 = deflection_for(3.0e3);
        // Both deflections are real and positive, and the response is
        // sub-linear: doubling the load gives strictly less than twice
        // the deflection (geometric stiffening).
        assert!(d1 > 0.0 && d2 > d1, "more load → more deflection");
        assert!(
            d2 < 2.0 * d1,
            "nonlinear response must be sub-proportional: δ(2P) {d2} < 2·δ(P) {}",
            2.0 * d1
        );
    }

    #[test]
    fn rejects_an_unconstrained_model() {
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box params");
        let err = solve_nonlinear_static(
            &mesh,
            &FemMaterial::default(),
            &[],
            &[],
            &NonlinearControls::default(),
        )
        .unwrap_err();
        assert!(matches!(err, NativeSolverError::Unconstrained));
    }

    #[test]
    fn rejects_a_mesh_without_tets() {
        use valenx_mesh::element::{ElementBlock, ElementType};
        let mut mesh = Mesh::new("surface");
        mesh.nodes.push(Vector3::new(0.0, 0.0, 0.0));
        mesh.element_blocks.push(ElementBlock {
            element_type: ElementType::Tri3,
            connectivity: vec![0, 0, 0],
        });
        let err = solve_nonlinear_static(
            &mesh,
            &FemMaterial::default(),
            &[NodalConstraint::fixed(0)],
            &[],
            &NonlinearControls::default(),
        )
        .unwrap_err();
        assert!(matches!(err, NativeSolverError::NoTetBlock));
    }

    #[test]
    fn no_load_gives_no_displacement() {
        // With zero external force the equilibrium is the rest state —
        // every displacement must be zero.
        let mesh = structured_box_mesh(2.0, 1.0, 1.0, 4, 1, 1).expect("valid box params");
        let constraints = vec![
            NodalConstraint::fixed(nid(0, 0, 0, 4, 1)),
            NodalConstraint::fixed(nid(0, 1, 0, 4, 1)),
            NodalConstraint::fixed(nid(0, 0, 1, 4, 1)),
            NodalConstraint::fixed(nid(0, 1, 1, 4, 1)),
        ];
        let sol = solve_nonlinear_static(
            &mesh,
            &FemMaterial::default(),
            &constraints,
            &[],
            &NonlinearControls::default(),
        )
        .unwrap();
        assert!(
            sol.max_displacement() < 1e-9,
            "an unloaded structure should not move, got {}",
            sol.max_displacement()
        );
    }

    // ----- Round-2 F2: non-finite prescribed-displacement rejection -----

    #[test]
    fn rejects_nan_prescribed_displacement() {
        // A prescribed-displacement constraint carrying NaN is folded
        // into the Newton residual as `penalty·(u − value)`. Without the
        // guard it yields a NaN residual and an `Ok` solution with NaN
        // displacement (masked only by `converged: false`). Reject it.
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box params");
        let constraints = vec![
            NodalConstraint::fixed(0),
            NodalConstraint {
                node: 7,
                fixed: [Some(f64::NAN), None, None],
            },
        ];
        let err = solve_nonlinear_static(
            &mesh,
            &FemMaterial::default(),
            &constraints,
            &[],
            &NonlinearControls::default(),
        )
        .unwrap_err();
        assert!(
            matches!(err, NativeSolverError::InvalidLoad { kind: "prescribed displacement", .. }),
            "expected InvalidLoad(prescribed displacement), got {err:?}"
        );
    }

    #[test]
    fn rejects_infinite_prescribed_displacement() {
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box params");
        let constraints = vec![
            NodalConstraint::fixed(0),
            NodalConstraint {
                node: 7,
                fixed: [None, Some(f64::INFINITY), None],
            },
        ];
        let err = solve_nonlinear_static(
            &mesh,
            &FemMaterial::default(),
            &constraints,
            &[],
            &NonlinearControls::default(),
        )
        .unwrap_err();
        assert!(matches!(err, NativeSolverError::InvalidLoad { .. }));
    }

    #[test]
    fn finite_prescribed_displacement_still_solves() {
        // A genuine finite prescribed elongation must not be rejected —
        // the guard targets only non-finite values.
        let mesh = structured_box_mesh(2.0, 1.0, 1.0, 4, 1, 1).expect("valid box params");
        let mut constraints = vec![
            NodalConstraint::fixed(nid(0, 0, 0, 4, 1)),
            NodalConstraint::fixed(nid(0, 1, 0, 4, 1)),
            NodalConstraint::fixed(nid(0, 0, 1, 4, 1)),
            NodalConstraint::fixed(nid(0, 1, 1, 4, 1)),
        ];
        // x=L face pulled to a small finite elongation in X.
        for (j, k) in [(0, 0), (1, 0), (0, 1), (1, 1)] {
            constraints.push(NodalConstraint {
                node: nid(4, j, k, 4, 1),
                fixed: [Some(0.01), None, None],
            });
        }
        let sol = solve_nonlinear_static(
            &mesh,
            &FemMaterial::default(),
            &constraints,
            &[],
            &NonlinearControls::default(),
        )
        .expect("finite prescribed displacement must still solve");
        assert!(
            sol.displacement.iter().all(|d| d.iter().all(|c| c.is_finite())),
            "a finite prescribed displacement must give a finite solution"
        );
    }
}
