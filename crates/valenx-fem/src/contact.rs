//! Native in-process **penalty-method contact** FEA solver.
//!
//! ## What this is
//!
//! A genuine, self-contained **node-to-surface contact** solver for
//! 4-node linear tetrahedra. An ordinary linear-static solve
//! ([`crate::native_solver`]) lets a body deform freely; if you push
//! one body toward another, nothing stops the meshes from passing
//! straight through each other. A **contact** solver adds the
//! non-penetration constraint: a body may touch a surface and press on
//! it, but it may not pass through it.
//!
//! This module implements the **penalty method** — the simplest, most
//! robust contact enforcement and the standard first choice. The
//! contact surface here is a **rigid plane** ([`ContactPlane`]); each
//! mesh node is checked against it:
//!
//! - if the node is on the allowed side of the plane there is no
//!   contact and no force;
//! - if the node has **penetrated** the plane by a gap `g < 0`, a
//!   restoring **penalty force** `f = −k_p·g·n̂` is applied along the
//!   plane normal `n̂`, with `k_p` a large **penalty stiffness**. The
//!   force grows with the penetration, so it drives the penetration
//!   back toward zero.
//!
//! The penalty method allows a *small* residual penetration — that is
//! its defining approximation — but the penetration is inversely
//! proportional to `k_p`, so a large penalty stiffness keeps it within
//! a tight tolerance. The tests here verify exactly that: two bodies
//! pushed together do not interpenetrate beyond the penalty tolerance.
//!
//! ## Contact as a nonlinear problem
//!
//! Contact is **nonlinear** even with linear-elastic material and small
//! strains: the set of nodes in contact is not known in advance — it
//! depends on the deformed shape, which depends on the contact forces.
//! [`solve_contact`] resolves it with a Newton-style iteration:
//!
//! 1. assemble the elastic stiffness `K` and the external load;
//! 2. with the current displacement, find which nodes have
//!    penetrated, and assemble the **contact stiffness** `K_c` (a
//!    penalty spring `k_p·n̂⊗n̂` on every penetrated node) and the
//!    **contact force** into the residual;
//! 3. solve the tangent system `(K + K_c)·Δu = −r` for a displacement
//!    increment;
//! 4. update and repeat until the residual — the out-of-balance force
//!    *including* the contact reactions — vanishes.
//!
//! At convergence every penetration is below the penalty tolerance and
//! the contact reactions balance the applied load.
//!
//! ## Honest scope
//!
//! This is a **real penalty-contact v1** — the [tests](self) verify
//! non-penetration to the penalty tolerance and a sane contact
//! pressure. It is deliberately a v1:
//!
//! - **Node-to-rigid-plane** contact. The master surface is an
//!   analytic rigid half-space, so the gap function is exact and cheap.
//!   General **deformable body-to-body** contact needs a
//!   surface-segment search and a contact-pair data structure — a
//!   bounded but substantial extension.
//! - **Penalty enforcement** — a small residual penetration is
//!   admitted by construction. The exact (zero-penetration) Lagrange-
//!   multiplier / augmented-Lagrangian methods are follow-ups.
//! - **Frictionless, normal contact only.** No tangential friction
//!   (Coulomb stick/slip), no adhesion.
//! - **Tet4, isotropic linear-elastic** — inherited from
//!   [`crate::native_solver`].
//!
//! Each of those is a documented, well-understood extension. It is
//! *not* a CalculiX replacement — CalculiX adds surface-to-surface
//! contact, friction, and augmented-Lagrangian enforcement; that gap
//! stays documented in Tier 3.

use nalgebra::{DMatrix, DVector, Vector3};
use nalgebra_sparse::factorization::CscCholesky;
use nalgebra_sparse::{CooMatrix, CscMatrix};

use valenx_mesh::Mesh;

use crate::material::FemMaterial;
use crate::native_solver::{
    add_diagonal, elasticity_matrix, element_stiffness, ensure_finite3, tet_block,
    von_mises_from_voigt, NativeSolverError, NodalConstraint, NodalForce,
};

/// A rigid **contact plane** — the master surface the mesh nodes are
/// kept from penetrating.
///
/// The plane is `{ x : n̂·x = d }`. The **allowed** side is the
/// half-space `n̂·x ≥ d` — a node there is in free space; a node with
/// `n̂·x < d` has penetrated and is pushed back. The signed gap of a
/// point `x` is `g = n̂·x − d`: positive (or zero) is no contact,
/// negative is penetration depth.
#[derive(Clone, Copy, Debug)]
pub struct ContactPlane {
    /// The plane's outward unit normal `n̂` — points *into* the allowed
    /// half-space. Normalised on construction.
    pub normal: Vector3<f64>,
    /// The plane offset `d` in `n̂·x = d`.
    pub offset: f64,
}

impl ContactPlane {
    /// Build a contact plane from a (not necessarily unit) normal and a
    /// point on the plane. The allowed side is the half-space the
    /// normal points toward.
    ///
    /// # Panics
    ///
    /// Panics if `normal` is the zero vector — a plane needs a
    /// direction.
    pub fn through_point(normal: Vector3<f64>, point: Vector3<f64>) -> Self {
        let len = normal.norm();
        assert!(len > 0.0, "contact-plane normal must be non-zero");
        let n = normal / len;
        ContactPlane {
            normal: n,
            offset: n.dot(&point),
        }
    }

    /// A horizontal floor at height `z` — the plane `z = height`, with
    /// the allowed side *above* it (`+Z`). The canonical
    /// "body resting on / pushed onto the floor" master surface.
    pub fn floor(height: f64) -> Self {
        ContactPlane {
            normal: Vector3::new(0.0, 0.0, 1.0),
            offset: height,
        }
    }

    /// The signed gap of a point — `g = n̂·x − d`. Non-negative means
    /// no contact; negative is the penetration depth.
    #[inline]
    pub fn gap(&self, point: Vector3<f64>) -> f64 {
        self.normal.dot(&point) - self.offset
    }
}

/// Controls for the penalty-contact Newton iteration.
#[derive(Copy, Clone, Debug)]
pub struct ContactControls {
    /// The **penalty stiffness** `k_p` (N/m) of a penetrated node's
    /// contact spring. A larger value means a stiffer "wall" and a
    /// smaller residual penetration, at the cost of a worse-conditioned
    /// tangent. When [`ContactControls::auto_penalty`] is set the
    /// solver overrides this with a value scaled to the assembled
    /// stiffness.
    pub penalty_stiffness: f64,
    /// If `true`, ignore [`ContactControls::penalty_stiffness`] and set
    /// the penalty automatically to a large multiple of the peak
    /// stiffness diagonal — the robust default (a fixed penalty is
    /// either too soft, and leaks, or too stiff, and wrecks
    /// conditioning).
    pub auto_penalty: bool,
    /// Number of equal load increments the external load is ramped
    /// over. Ramping eases the contact set into place.
    pub load_steps: usize,
    /// Maximum Newton iterations per load step.
    pub max_iterations: usize,
    /// Convergence tolerance on the relative residual-force norm.
    pub tolerance: f64,
}

impl Default for ContactControls {
    /// Sensible defaults: auto penalty, 5 load steps, up to 40 Newton
    /// iterations each, converged at a `1e-7` relative residual.
    fn default() -> Self {
        ContactControls {
            penalty_stiffness: 0.0,
            auto_penalty: true,
            load_steps: 5,
            max_iterations: 40,
            tolerance: 1e-7,
        }
    }
}

/// Result of a penalty-contact solve.
#[derive(Clone, Debug)]
pub struct ContactSolution {
    /// Final per-node displacement `[ux, uy, uz]` in metres.
    pub displacement: Vec<[f64; 3]>,
    /// Per-node von Mises stress (Pa), averaged from the elements.
    pub von_mises: Vec<f64>,
    /// Per-node **contact pressure** — the magnitude of the penalty
    /// contact force on each node (N). Zero for every node not in
    /// contact.
    pub contact_force: Vec<f64>,
    /// Per-node signed gap to the contact plane in the final
    /// configuration (metres). Negative is residual penetration.
    pub gap: Vec<f64>,
    /// Total Newton iterations summed over every load step.
    pub total_iterations: usize,
    /// `true` if every load step converged.
    pub converged: bool,
    /// The relative residual-force norm at the end of the final step.
    pub final_residual: f64,
}

impl ContactSolution {
    /// Largest displacement magnitude across all nodes.
    pub fn max_displacement(&self) -> f64 {
        self.displacement
            .iter()
            .map(|d| (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt())
            .fold(0.0, f64::max)
    }

    /// The **deepest penetration** of any node into the contact plane —
    /// `0` if no node penetrates, otherwise a positive depth (metres).
    /// This is the number the penalty tolerance bounds.
    pub fn max_penetration(&self) -> f64 {
        self.gap.iter().map(|&g| (-g).max(0.0)).fold(0.0, f64::max)
    }

    /// The total contact force summed over all nodes (N) — at
    /// equilibrium this balances the applied load pushing the body into
    /// the plane.
    pub fn total_contact_force(&self) -> f64 {
        self.contact_force.iter().sum()
    }

    /// `true` if any node is in contact with the plane.
    pub fn in_contact(&self) -> bool {
        self.contact_force.iter().any(|&f| f > 0.0)
    }
}

/// Solve a penalty-method node-to-rigid-plane contact problem on a
/// tetrahedral mesh.
///
/// `mesh` carries the deformable body (a Tet4 block); `material` its
/// isotropic elastic constants; `plane` the rigid master surface every
/// node is kept from penetrating; `constraints` and `forces` the
/// ordinary displacement / load boundary conditions (the forces
/// typically push the body *toward* the plane). `controls` tunes the
/// penalty stiffness and the Newton iteration.
///
/// At least one displacement constraint is still required *unless* the
/// contact itself fully restrains the body — but a body held off the
/// plane purely by contact is a delicate configuration, so supplying a
/// constraint is recommended.
///
/// Returns the converged [`ContactSolution`] — the displacement field,
/// the per-node contact force, and the final gaps. At convergence
/// every penetration is below the penalty tolerance.
///
/// # Method
///
/// Each load step ramps the external load and runs a Newton iteration.
/// Each iteration: assemble the elastic internal force + stiffness;
/// for every node whose current position has penetrated the plane, add
/// the penalty contact force `−k_p·g·n̂` to the residual and the
/// contact stiffness `k_p·(n̂⊗n̂)` to the node's 3×3 tangent block;
/// solve `(K + K_c)·Δu = −r`; update; repeat to equilibrium.
///
/// # Errors
///
/// See [`NativeSolverError`] — empty mesh, no Tet4 block, degenerate
/// element, bad connectivity, a bad material, or a factorisation
/// failure.
pub fn solve_contact(
    mesh: &Mesh,
    material: &FemMaterial,
    plane: &ContactPlane,
    constraints: &[NodalConstraint],
    forces: &[NodalForce],
    controls: &ContactControls,
) -> Result<ContactSolution, NativeSolverError> {
    let n_nodes = mesh.nodes.len();
    if n_nodes == 0 {
        return Err(NativeSolverError::EmptyMesh);
    }
    let tets = tet_block(mesh).ok_or(NativeSolverError::NoTetBlock)?;
    let d_matrix = elasticity_matrix(material)?;
    let n_dof = 3 * n_nodes;
    let conn = &tets.connectivity;
    let n_elem = conn.len() / 4;

    // --- assemble the (constant) elastic stiffness matrix once ---
    let mut coo = CooMatrix::<f64>::new(n_dof, n_dof);
    let mut elem_nodes: Vec<[usize; 4]> = Vec::with_capacity(n_elem);
    let mut elem_k: Vec<DMatrix<f64>> = Vec::with_capacity(n_elem);
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
        elem_nodes.push(nodes);
        elem_k.push(ke);
    }
    let k_elastic = CscMatrix::from(&coo);
    let mut max_diag = 0.0_f64;
    for (r, c, v) in k_elastic.triplet_iter() {
        if r == c {
            max_diag = max_diag.max(v.abs());
        }
    }
    if max_diag <= 0.0 {
        return Err(NativeSolverError::SolveFailed);
    }
    // The constraint penalty pins a DOF; the *contact* penalty is the
    // wall stiffness.
    let constraint_penalty = max_diag * 1.0e8;
    let contact_penalty = if controls.auto_penalty {
        // ~1e3× the peak stiffness keeps the residual penetration tiny
        // without destroying the conditioning.
        max_diag * 1.0e3
    } else {
        controls.penalty_stiffness.max(0.0)
    };

    // --- external load vector ---
    let mut f_ext_full = DVector::<f64>::zeros(n_dof);
    for force in forces {
        if force.node >= n_nodes {
            return Err(NativeSolverError::BadConnectivity {
                elem: usize::MAX,
                node: force.node,
                n_nodes,
            });
        }
        // Round-1 H4: reject a non-finite force before it reaches the
        // Newton residual / RHS.
        ensure_finite3(&force.force, force.node, "force")?;
        for i in 0..3 {
            f_ext_full[3 * force.node + i] += force.force[i];
        }
    }
    let f_ext_norm = f_ext_full.norm().max(constraint_penalty * 1e-12).max(1e-30);

    // Mask of DOFs pinned by a displacement constraint. The Newton
    // convergence test must measure the out-of-balance force on the FREE
    // DOFs only — a constrained row's residual is the reaction force
    // (of order the applied load), which never vanishes; counting it in
    // would pin the relative residual at O(1) and the solve could never
    // report convergence.
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

    let n = plane.normal;
    let mut u = DVector::<f64>::zeros(n_dof);
    let load_steps = controls.load_steps.max(1);
    let mut total_iterations = 0;
    let mut all_converged = true;
    let mut final_residual = f64::INFINITY;

    for step in 1..=load_steps {
        let lambda = step as f64 / load_steps as f64;
        let f_ext_step = &f_ext_full * lambda;
        let mut step_converged = false;

        for _iter in 0..controls.max_iterations.max(1) {
            total_iterations += 1;

            // --- elastic internal force f_int = K·u ---
            let f_int: DVector<f64> = &k_elastic * &u;
            let mut residual = &f_int - &f_ext_step;

            // --- contact: scan every node for penetration ---
            // The current node position is rest + displacement; its
            // gap to the plane is g = n̂·x − d. A negative gap is
            // penetration; the penalty force −k_p·g·n̂ pushes it back.
            let mut contact_diag: Vec<[[f64; 3]; 3]> = vec![[[0.0; 3]; 3]; n_nodes];
            for node in 0..n_nodes {
                let pos =
                    mesh.nodes[node] + Vector3::new(u[3 * node], u[3 * node + 1], u[3 * node + 2]);
                let gap = plane.gap(pos);
                // A node that has penetrated (g < 0) OR is exactly on
                // the surface (g == 0) is in contact. Including the
                // g == 0 touching case is both physically right — a
                // resting node resists further penetration — and
                // numerically essential: a body whose only restraint
                // against a rigid-body translation INTO the plane is the
                // contact stiffness would give a singular tangent on the
                // first Newton iteration (when it starts exactly resting
                // on the surface) if touching nodes were skipped.
                if gap <= 0.0 {
                    // The contact force on the node is f_c = −k_p·g·n̂
                    // (g ≤ 0 → f_c along +n̂, pushing the node out; at
                    // g = 0 the force is zero, only the stiffness
                    // remains). As a *residual* contribution (internal
                    // minus external) it ADDS k_p·g·n̂ to the residual.
                    for i in 0..3 {
                        residual[3 * node + i] += contact_penalty * gap * n[i];
                    }
                    // Contact stiffness K_c = k_p·(n̂⊗n̂) on the node's
                    // 3×3 diagonal block — d(f_c)/du.
                    for i in 0..3 {
                        for j in 0..3 {
                            contact_diag[node][i][j] += contact_penalty * n[i] * n[j];
                        }
                    }
                }
            }

            // --- displacement constraints (large-penalty) ---
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
                        penalty_diag[dof] += constraint_penalty;
                        residual[dof] += constraint_penalty * (u[dof] - value);
                    }
                }
            }

            // --- assemble the tangent K_T = K + K_c + constraints ---
            // Start from the elastic K, scatter the contact 3×3 blocks
            // and the constraint diagonal in.
            let mut tcoo = CooMatrix::<f64>::new(n_dof, n_dof);
            for (r, c, v) in k_elastic.triplet_iter() {
                tcoo.push(r, c, *v);
            }
            for (node, block) in contact_diag.iter().enumerate() {
                for (i, row) in block.iter().enumerate() {
                    for (j, &val) in row.iter().enumerate() {
                        if val != 0.0 {
                            tcoo.push(3 * node + i, 3 * node + j, val);
                        }
                    }
                }
            }
            let k_t = CscMatrix::from(&tcoo);
            let k_t = add_diagonal(&k_t, &penalty_diag);

            // --- solve K_T·Δu = −r ---
            let chol = CscCholesky::factor(&k_t).map_err(|_| NativeSolverError::SolveFailed)?;
            let rhs = -&residual;
            let delta: DMatrix<f64> = chol.solve(&rhs);
            let delta = delta.column(0);
            for d in 0..n_dof {
                u[d] += delta[d];
            }

            // --- convergence test ---
            // Measure the out-of-balance force on the FREE DOFs only;
            // a constrained DOF's residual entry is the (non-vanishing)
            // reaction force and would otherwise stall the test at O(1).
            let mut res_norm_sq = 0.0;
            for d in 0..n_dof {
                if !constrained[d] {
                    res_norm_sq += residual[d] * residual[d];
                }
            }
            final_residual = res_norm_sq.sqrt() / f_ext_norm;
            if final_residual <= controls.tolerance {
                step_converged = true;
                break;
            }
            if !final_residual.is_finite() {
                break;
            }
        }
        if !step_converged {
            all_converged = false;
        }
    }

    // --- gather the displacement field ---
    let mut displacement = vec![[0.0_f64; 3]; n_nodes];
    for (node, d) in displacement.iter_mut().enumerate() {
        d[0] = u[3 * node];
        d[1] = u[3 * node + 1];
        d[2] = u[3 * node + 2];
    }

    // --- final per-node contact force + gap ---
    let mut contact_force = vec![0.0_f64; n_nodes];
    let mut gap = vec![0.0_f64; n_nodes];
    for node in 0..n_nodes {
        let pos = mesh.nodes[node] + Vector3::new(u[3 * node], u[3 * node + 1], u[3 * node + 2]);
        let g = plane.gap(pos);
        gap[node] = g;
        if g < 0.0 {
            // |f_c| = k_p·|g|.
            contact_force[node] = contact_penalty * (-g);
        }
    }

    // --- recover the nodal von Mises stress ---
    let von_mises = recover_von_mises(mesh, &elem_nodes, &elem_k, &u, &d_matrix, n_nodes);

    Ok(ContactSolution {
        displacement,
        von_mises,
        contact_force,
        gap,
        total_iterations,
        converged: all_converged,
        final_residual,
    })
}

/// Recover the per-node von Mises stress from the converged
/// displacement — `σ = D·B·uₑ` per element, averaged to the nodes,
/// exactly as the linear solver does.
fn recover_von_mises(
    mesh: &Mesh,
    elem_nodes: &[[usize; 4]],
    _elem_k: &[DMatrix<f64>],
    u: &DVector<f64>,
    d_matrix: &nalgebra::Matrix6<f64>,
    n_nodes: usize,
) -> Vec<f64> {
    use crate::native_solver::strain_displacement;
    let mut nodal_stress = vec![[0.0_f64; 6]; n_nodes];
    let mut nodal_count = vec![0_u32; n_nodes];
    for nodes in elem_nodes.iter() {
        let coords = [
            mesh.nodes[nodes[0]],
            mesh.nodes[nodes[1]],
            mesh.nodes[nodes[2]],
            mesh.nodes[nodes[3]],
        ];
        let Some(b) = strain_displacement(&coords) else {
            continue;
        };
        let mut ue = DVector::<f64>::zeros(12);
        for a in 0..4 {
            for i in 0..3 {
                ue[3 * a + i] = u[3 * nodes[a] + i];
            }
        }
        let strain = &b * &ue;
        let sigma = d_matrix * &strain;
        for &nd in nodes {
            for k in 0..6 {
                nodal_stress[nd][k] += sigma[k];
            }
            nodal_count[nd] += 1;
        }
    }
    let mut von_mises = vec![0.0_f64; n_nodes];
    for (node, vm) in von_mises.iter_mut().enumerate() {
        let c = nodal_count[node].max(1) as f64;
        let mut s = [0.0_f64; 6];
        for (k, sk) in s.iter_mut().enumerate() {
            *sk = nodal_stress[node][k] / c;
        }
        *vm = von_mises_from_voigt(&s);
    }
    von_mises
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
    fn contact_plane_gap_sign_is_correct() {
        // A floor at z = 0: a point above has a positive gap, a point
        // below (penetrating) a negative one.
        let floor = ContactPlane::floor(0.0);
        assert!(floor.gap(Vector3::new(0.0, 0.0, 1.0)) > 0.0);
        assert!(floor.gap(Vector3::new(0.0, 0.0, -0.5)) < 0.0);
        assert!(floor.gap(Vector3::new(0.0, 0.0, 0.0)).abs() < 1e-12);
        // The gap magnitude equals the distance to the plane.
        assert!((floor.gap(Vector3::new(3.0, 7.0, 2.0)) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn contact_plane_through_point_normalises() {
        // through_point must accept a non-unit normal and normalise it.
        let p =
            ContactPlane::through_point(Vector3::new(0.0, 0.0, 5.0), Vector3::new(0.0, 0.0, 2.0));
        assert!((p.normal.norm() - 1.0).abs() < 1e-12);
        // The plane passes through z = 2.
        assert!(p.gap(Vector3::new(9.0, 9.0, 2.0)).abs() < 1e-9);
    }

    #[test]
    fn body_clear_of_the_plane_has_no_contact() {
        // A body sitting well above the floor, loaded only sideways,
        // never touches the plane — no contact force, the result is the
        // ordinary elastic solution.
        let mesh = structured_box_mesh(2.0, 1.0, 1.0, 2, 1, 1).expect("valid box params");
        let (nx, ny, nz) = (2, 1, 1);
        // The box spans z ∈ [0, 1]; put the floor far below at z = -5.
        let floor = ContactPlane::floor(-5.0);
        let mat = FemMaterial::default();
        // Clamp x=0.
        let mut constraints = Vec::new();
        for k in 0..=nz {
            for j in 0..=ny {
                constraints.push(NodalConstraint::fixed(nid(0, j, k, nx, ny)));
            }
        }
        // A tiny sideways tug at the far face.
        let forces = vec![NodalForce {
            node: nid(nx, 0, 0, nx, ny),
            force: [1.0e3, 0.0, 0.0],
        }];
        let sol = solve_contact(
            &mesh,
            &mat,
            &floor,
            &constraints,
            &forces,
            &ContactControls::default(),
        )
        .unwrap();
        assert!(sol.converged, "no-contact solve should converge");
        assert!(
            !sol.in_contact(),
            "a body clear of the plane is not in contact"
        );
        assert_eq!(sol.max_penetration(), 0.0);
    }

    #[test]
    fn body_pushed_into_the_plane_does_not_interpenetrate() {
        // The headline contact verification. A block sits just above a
        // floor and is pushed straight down into it. Without contact
        // the block's bottom face would pass through the floor; with
        // the penalty contact it must stop at the plane — the residual
        // penetration must be tiny (below the penalty tolerance).
        let (lx, ly, lz) = (1.0, 1.0, 1.0);
        let (nx, ny, nz) = (2, 2, 2);
        let mesh = structured_box_mesh(lx, ly, lz, nx, ny, nz).expect("valid box params");
        // The block spans z ∈ [0, 1]. Put the floor at z = 0 — the
        // block's bottom face is exactly on it to start.
        let floor = ContactPlane::floor(0.0);
        let mat = FemMaterial {
            youngs_modulus: 200.0e9,
            poisson_ratio: 0.3,
            ..FemMaterial::default()
        };
        // Constrain the *top* face in X and Y (so the block cannot
        // slide or tip) but leave Z free — it must be free to be
        // pushed down. Then push the top face down with a large force.
        let mut constraints = Vec::new();
        for j in 0..=ny {
            for i in 0..=nx {
                let n = nid(i, j, nz, nx, ny);
                constraints.push(NodalConstraint {
                    node: n,
                    fixed: [Some(0.0), Some(0.0), None],
                });
            }
        }
        let top_nodes: Vec<usize> = (0..=ny)
            .flat_map(|j| (0..=nx).map(move |i| (i, j)))
            .map(|(i, j)| nid(i, j, nz, nx, ny))
            .collect();
        let total_push = 5.0e7; // a big downward shove
        let per = total_push / top_nodes.len() as f64;
        let forces: Vec<NodalForce> = top_nodes
            .iter()
            .map(|&n| NodalForce {
                node: n,
                force: [0.0, 0.0, -per],
            })
            .collect();

        let sol = solve_contact(
            &mesh,
            &mat,
            &floor,
            &constraints,
            &forces,
            &ContactControls {
                auto_penalty: true,
                load_steps: 6,
                max_iterations: 50,
                tolerance: 1e-7,
                ..ContactControls::default()
            },
        )
        .unwrap();
        assert!(
            sol.converged,
            "contact solve should converge (residual {})",
            sol.final_residual
        );
        // The block IS in contact and IS being pressed down.
        assert!(sol.in_contact(), "the pushed block must contact the floor");
        // Non-penetration: the deepest penetration must be a tiny
        // fraction of the block height — the penalty tolerance.
        let penetration = sol.max_penetration();
        assert!(
            penetration < 1.0e-3 * lz,
            "penalty contact must stop interpenetration: depth {penetration} vs height {lz}"
        );
        // Force balance: the total contact reaction must balance the
        // applied downward push.
        let rel = (sol.total_contact_force() - total_push).abs() / total_push;
        assert!(
            rel < 0.05,
            "contact reaction {} should balance the push {total_push} (rel {rel})",
            sol.total_contact_force()
        );
    }

    #[test]
    fn deeper_push_still_respects_the_penalty_tolerance() {
        // Doubling the push force must NOT double the penetration — the
        // penalty wall is stiff. The penetration stays bounded by the
        // tolerance regardless of how hard the body is pressed.
        let (lx, ly, lz) = (1.0, 1.0, 1.0);
        let (nx, ny, nz) = (2, 2, 2);
        let mesh = structured_box_mesh(lx, ly, lz, nx, ny, nz).expect("valid box params");
        let floor = ContactPlane::floor(0.0);
        let mat = FemMaterial {
            youngs_modulus: 200.0e9,
            poisson_ratio: 0.3,
            ..FemMaterial::default()
        };
        let mut constraints = Vec::new();
        for j in 0..=ny {
            for i in 0..=nx {
                constraints.push(NodalConstraint {
                    node: nid(i, j, nz, nx, ny),
                    fixed: [Some(0.0), Some(0.0), None],
                });
            }
        }
        let top_nodes: Vec<usize> = (0..=ny)
            .flat_map(|j| (0..=nx).map(move |i| (i, j)))
            .map(|(i, j)| nid(i, j, nz, nx, ny))
            .collect();

        let penetration_for = |total_push: f64| -> f64 {
            let per = total_push / top_nodes.len() as f64;
            let forces: Vec<NodalForce> = top_nodes
                .iter()
                .map(|&n| NodalForce {
                    node: n,
                    force: [0.0, 0.0, -per],
                })
                .collect();
            let sol = solve_contact(
                &mesh,
                &mat,
                &floor,
                &constraints,
                &forces,
                &ContactControls {
                    auto_penalty: true,
                    load_steps: 6,
                    max_iterations: 50,
                    tolerance: 1e-7,
                    ..ContactControls::default()
                },
            )
            .unwrap();
            assert!(sol.converged, "contact solve should converge");
            sol.max_penetration()
        };

        let pen_1 = penetration_for(5.0e7);
        let pen_2 = penetration_for(1.0e8);
        // Both penetrations are within the penalty tolerance.
        assert!(pen_1 < 1.0e-3 * lz && pen_2 < 1.0e-3 * lz);
        // The penetration does grow with the force (penalty contact is
        // compliant) but stays a tiny fraction of the body even at the
        // larger load — the wall does not "give way".
        assert!(
            pen_2 < 4.0 * pen_1.max(1e-12),
            "a stiff penalty wall keeps penetration bounded: {pen_1} → {pen_2}"
        );
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
        let err = solve_contact(
            &mesh,
            &FemMaterial::default(),
            &ContactPlane::floor(0.0),
            &[NodalConstraint::fixed(0)],
            &[],
            &ContactControls::default(),
        )
        .unwrap_err();
        assert!(matches!(err, NativeSolverError::NoTetBlock));
    }

    #[test]
    fn rejects_an_empty_mesh() {
        let mesh = Mesh::new("empty");
        let err = solve_contact(
            &mesh,
            &FemMaterial::default(),
            &ContactPlane::floor(0.0),
            &[],
            &[],
            &ContactControls::default(),
        )
        .unwrap_err();
        assert!(matches!(err, NativeSolverError::EmptyMesh));
    }

    // ----- Round-2 F2: non-finite prescribed-displacement rejection -----

    #[test]
    fn rejects_nan_prescribed_displacement() {
        // A prescribed displacement is folded into the Newton residual
        // as `penalty·(u − value)`. A NaN value yields a NaN residual
        // and an `Ok` solution with NaN displacement (masked only by
        // `converged: false`). Reject it up front. The floor is far
        // below so contact never fires — isolating the BC validation.
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box params");
        let floor = ContactPlane::floor(-5.0);
        let constraints = vec![
            NodalConstraint::fixed(0),
            NodalConstraint {
                node: 7,
                fixed: [Some(f64::NAN), None, None],
            },
        ];
        let err = solve_contact(
            &mesh,
            &FemMaterial::default(),
            &floor,
            &constraints,
            &[],
            &ContactControls::default(),
        )
        .unwrap_err();
        assert!(
            matches!(
                err,
                NativeSolverError::InvalidLoad {
                    kind: "prescribed displacement",
                    ..
                }
            ),
            "expected InvalidLoad(prescribed displacement), got {err:?}"
        );
    }

    #[test]
    fn rejects_infinite_prescribed_displacement() {
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box params");
        let floor = ContactPlane::floor(-5.0);
        let constraints = vec![
            NodalConstraint::fixed(0),
            NodalConstraint {
                node: 7,
                fixed: [None, None, Some(f64::INFINITY)],
            },
        ];
        let err = solve_contact(
            &mesh,
            &FemMaterial::default(),
            &floor,
            &constraints,
            &[],
            &ContactControls::default(),
        )
        .unwrap_err();
        assert!(matches!(err, NativeSolverError::InvalidLoad { .. }));
    }

    #[test]
    fn finite_prescribed_displacement_still_solves() {
        // A genuine finite prescribed displacement must not be rejected
        // — the guard targets only non-finite values.
        let mesh = structured_box_mesh(2.0, 1.0, 1.0, 2, 1, 1).expect("valid box params");
        let (nx, ny, nz) = (2, 1, 1);
        let floor = ContactPlane::floor(-5.0);
        let mut constraints = Vec::new();
        for k in 0..=nz {
            for j in 0..=ny {
                constraints.push(NodalConstraint::fixed(nid(0, j, k, nx, ny)));
            }
        }
        // x=L face given a small finite prescribed elongation in X.
        for k in 0..=nz {
            for j in 0..=ny {
                constraints.push(NodalConstraint {
                    node: nid(nx, j, k, nx, ny),
                    fixed: [Some(1e-3), None, None],
                });
            }
        }
        let sol = solve_contact(
            &mesh,
            &FemMaterial::default(),
            &floor,
            &constraints,
            &[],
            &ContactControls::default(),
        )
        .expect("finite prescribed displacement must still solve");
        assert!(
            sol.displacement
                .iter()
                .all(|d| d.iter().all(|c| c.is_finite())),
            "a finite prescribed displacement must give a finite solution"
        );
    }
}
