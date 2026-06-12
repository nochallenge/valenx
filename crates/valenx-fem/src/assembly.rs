//! **Mixed-element global assembly** for the native FEA solvers
//! (Phase 24.8).
//!
//! ## What this is
//!
//! The original solvers assembled the global stiffness / mass matrix by
//! walking a *single* [`ElementType::Tet4`] block. This module
//! generalises that: it walks **every** continuum-element block in a
//! [`Mesh`] — Tet4, [`ElementType::Hex8`], [`ElementType::Tet10`], in
//! any mix — builds each element's `Kₑ` / `Mₑ` through the
//! [`crate::elements::SolidElement`] trait, and scatter-adds them into
//! one global system. A mesh that carries a Hex8 block *and* a Tet10
//! block assembles into a single coupled `K` (the two element families
//! simply share the global node numbering).
//!
//! It is the substrate the mixed-element linear-static
//! ([`crate::native_solver::solve_linear_static_mixed`]) and modal
//! ([`crate::modal_solver::solve_modal_mixed`]) solvers build on.
//!
//! ## What it does *not* touch
//!
//! The structural **beam** element ([`crate::beam`]) is assembled
//! separately — it carries 6 DOFs per node (3 translations + 3
//! rotations), so it cannot share the 3-DOF-per-node continuum global
//! numbering this module assumes.

use nalgebra::{DMatrix, DVector, Matrix6};
use nalgebra_sparse::{CooMatrix, CscMatrix};

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::elements::{Hex8, SolidElement, Tet10, Tet4};
use crate::native_solver::NativeSolverError;

/// The continuum-element types this assembly understands.
///
/// `Pyr5` / `Prism6` / `Hex20` exist in the mesh [`ElementType`] enum
/// but have no element implementation yet, so they are skipped (with
/// the rest of the mesh still assembled). 1D / 2D blocks (`Line2`,
/// `Tri3`, …) are likewise skipped — this is a *solid* assembly.
pub fn is_supported_solid(t: ElementType) -> bool {
    matches!(
        t,
        ElementType::Tet4 | ElementType::Hex8 | ElementType::Tet10
    )
}

/// Does this mesh carry at least one supported solid-element block with
/// a non-empty connectivity?
pub fn has_solid_elements(mesh: &Mesh) -> bool {
    mesh.element_blocks.iter().any(|b| {
        is_supported_solid(b.element_type)
            && b.connectivity.len() >= b.element_type.nodes_per_element()
    })
}

/// Collect the node-index list of every supported solid element in the
/// mesh — one `Vec<usize>` per element.
///
/// Used to build the node adjacency graph for the fill-reducing
/// reordering ([`crate::ordering`]).
pub fn solid_element_node_lists(mesh: &Mesh) -> Vec<Vec<usize>> {
    let mut lists = Vec::new();
    for block in &mesh.element_blocks {
        if !is_supported_solid(block.element_type) {
            continue;
        }
        let npe = block.element_type.nodes_per_element();
        if block.connectivity.len() < npe {
            continue;
        }
        let n_elem = block.connectivity.len() / npe;
        for e in 0..n_elem {
            lists.push(
                (0..npe)
                    .map(|k| block.connectivity[npe * e + k] as usize)
                    .collect(),
            );
        }
    }
    lists
}

/// Gather one element's physical node coordinates from a connectivity
/// block, returning `(global node indices, coordinates)`.
///
/// `None` (with the offending element / node recorded) if any
/// connectivity index is out of range.
fn gather_element(
    mesh: &Mesh,
    block: &ElementBlock,
    e: usize,
) -> Result<(Vec<usize>, Vec<nalgebra::Vector3<f64>>), NativeSolverError> {
    let npe = block.element_type.nodes_per_element();
    let mut nodes = Vec::with_capacity(npe);
    let mut coords = Vec::with_capacity(npe);
    for k in 0..npe {
        let nd = block.connectivity[npe * e + k] as usize;
        if nd >= mesh.nodes.len() {
            return Err(NativeSolverError::BadConnectivity {
                elem: e,
                node: nd,
                n_nodes: mesh.nodes.len(),
            });
        }
        nodes.push(nd);
        coords.push(mesh.nodes[nd]);
    }
    Ok((nodes, coords))
}

/// Build the concrete [`SolidElement`] for one element of a block.
///
/// Returns `None` if the block's element type is not a supported solid
/// element (the caller skips it).
fn make_element(
    element_type: ElementType,
    coords: &[nalgebra::Vector3<f64>],
) -> Option<Box<dyn SolidElement>> {
    match element_type {
        ElementType::Tet4 => {
            let c: [_; 4] = coords.try_into().ok()?;
            Some(Box::new(Tet4::new(c)))
        }
        ElementType::Hex8 => {
            let c: [_; 8] = coords.try_into().ok()?;
            Some(Box::new(Hex8::new(c)))
        }
        ElementType::Tet10 => {
            let c: [_; 10] = coords.try_into().ok()?;
            Some(Box::new(Tet10::new(c)))
        }
        _ => None,
    }
}

/// Assemble the global stiffness matrix `K` of a **mixed-element**
/// solid mesh as a dense `3·n_nodes` square.
///
/// Walks every supported solid-element block (Tet4 / Hex8 / Tet10),
/// builds each `Kₑ` through the [`SolidElement`] trait, and scatter-adds
/// it. Local DOF `3a+i` of element node `a` maps to global DOF
/// `3·node[a]+i`.
///
/// # Errors
///
/// [`NativeSolverError`] — empty mesh, no supported solid block, a
/// degenerate element, or out-of-range connectivity.
pub fn assemble_global_stiffness(
    mesh: &Mesh,
    d_matrix: &Matrix6<f64>,
) -> Result<DMatrix<f64>, NativeSolverError> {
    let n_nodes = mesh.nodes.len();
    if n_nodes == 0 {
        return Err(NativeSolverError::EmptyMesh);
    }
    if !has_solid_elements(mesh) {
        return Err(NativeSolverError::NoTetBlock);
    }
    // Round-1 H1–H3: bound the dense `n_dof × n_dof` allocation BEFORE
    // `DMatrix::zeros` — `check_dense_dofs` rejects an oversized mesh
    // (and a `3·n_nodes` overflow) up front rather than OOMing the host.
    let n_dof = crate::native_solver::check_dense_dofs(n_nodes)?;
    let mut k = DMatrix::<f64>::zeros(n_dof, n_dof);

    for block in &mesh.element_blocks {
        if !is_supported_solid(block.element_type) {
            continue;
        }
        let npe = block.element_type.nodes_per_element();
        if block.connectivity.len() < npe {
            continue;
        }
        let n_elem = block.connectivity.len() / npe;
        for e in 0..n_elem {
            let (nodes, coords) = gather_element(mesh, block, e)?;
            let elem = make_element(block.element_type, &coords)
                .ok_or(NativeSolverError::DegenerateElement(e))?;
            let ke = elem
                .stiffness(d_matrix)
                .ok_or(NativeSolverError::DegenerateElement(e))?;
            scatter_add(&mut k, &ke, &nodes);
        }
    }
    Ok(k)
}

/// Assemble the global stiffness matrix `K` of a mixed-element solid
/// mesh as a **sparse** [`CscMatrix`].
///
/// The sparse counterpart of [`assemble_global_stiffness`]. A finite-
/// element stiffness matrix is overwhelmingly sparse (a DOF couples
/// only to the DOFs of the elements that touch its node), so for any
/// non-trivial mesh the sparse path is dramatically faster and lighter
/// than the dense one — it is what the mixed-element static solver
/// uses. The duplicate `(i,j)` triplets the element scatter-add
/// produces are summed by the COO → CSC conversion.
///
/// # Errors
///
/// See [`NativeSolverError`].
pub fn assemble_global_stiffness_sparse(
    mesh: &Mesh,
    d_matrix: &Matrix6<f64>,
) -> Result<CscMatrix<f64>, NativeSolverError> {
    let n_nodes = mesh.nodes.len();
    if n_nodes == 0 {
        return Err(NativeSolverError::EmptyMesh);
    }
    if !has_solid_elements(mesh) {
        return Err(NativeSolverError::NoTetBlock);
    }
    let n_dof = 3 * n_nodes;
    let mut coo = CooMatrix::<f64>::new(n_dof, n_dof);

    for block in &mesh.element_blocks {
        if !is_supported_solid(block.element_type) {
            continue;
        }
        let npe = block.element_type.nodes_per_element();
        if block.connectivity.len() < npe {
            continue;
        }
        let n_elem = block.connectivity.len() / npe;
        for e in 0..n_elem {
            let (nodes, coords) = gather_element(mesh, block, e)?;
            let elem = make_element(block.element_type, &coords)
                .ok_or(NativeSolverError::DegenerateElement(e))?;
            let ke = elem
                .stiffness(d_matrix)
                .ok_or(NativeSolverError::DegenerateElement(e))?;
            scatter_add_sparse(&mut coo, &ke, &nodes);
        }
    }
    Ok(CscMatrix::from(&coo))
}

/// Assemble both the global stiffness `K` and **consistent mass** `M`
/// of a mixed-element solid mesh as dense `3·n_nodes` squares.
///
/// The mixed-element analogue of
/// [`crate::modal_solver::assemble_global_stiffness_mass`]. Returned
/// `(K, M)` in that order.
///
/// # Errors
///
/// See [`NativeSolverError`].
pub fn assemble_global_stiffness_mass_mixed(
    mesh: &Mesh,
    d_matrix: &Matrix6<f64>,
    density: f64,
) -> Result<(DMatrix<f64>, DMatrix<f64>), NativeSolverError> {
    let n_nodes = mesh.nodes.len();
    if n_nodes == 0 {
        return Err(NativeSolverError::EmptyMesh);
    }
    if !has_solid_elements(mesh) {
        return Err(NativeSolverError::NoTetBlock);
    }
    // Round-1 H1–H3: this path allocates TWO dense `n_dof × n_dof`
    // matrices, so the cap is doubly important. Check before either
    // `DMatrix::zeros`.
    let n_dof = crate::native_solver::check_dense_dofs(n_nodes)?;
    let mut k = DMatrix::<f64>::zeros(n_dof, n_dof);
    let mut m = DMatrix::<f64>::zeros(n_dof, n_dof);

    for block in &mesh.element_blocks {
        if !is_supported_solid(block.element_type) {
            continue;
        }
        let npe = block.element_type.nodes_per_element();
        if block.connectivity.len() < npe {
            continue;
        }
        let n_elem = block.connectivity.len() / npe;
        for e in 0..n_elem {
            let (nodes, coords) = gather_element(mesh, block, e)?;
            let elem = make_element(block.element_type, &coords)
                .ok_or(NativeSolverError::DegenerateElement(e))?;
            let ke = elem
                .stiffness(d_matrix)
                .ok_or(NativeSolverError::DegenerateElement(e))?;
            let me = elem
                .consistent_mass(density)
                .ok_or(NativeSolverError::DegenerateElement(e))?;
            scatter_add(&mut k, &ke, &nodes);
            scatter_add(&mut m, &me, &nodes);
        }
    }
    Ok((k, m))
}

/// Scatter-add a `3n × 3n` element matrix into the global dense matrix.
/// Local DOF `3a+i` → global DOF `3·nodes[a]+i`.
fn scatter_add(global: &mut DMatrix<f64>, elem: &DMatrix<f64>, nodes: &[usize]) {
    let n = nodes.len();
    for a in 0..n {
        for i in 0..3 {
            let gi = 3 * nodes[a] + i;
            for b in 0..n {
                for j in 0..3 {
                    let gj = 3 * nodes[b] + j;
                    global[(gi, gj)] += elem[(3 * a + i, 3 * b + j)];
                }
            }
        }
    }
}

/// Scatter-add a `3n × 3n` element matrix into a global sparse COO.
/// Duplicate `(i,j)` triplets are summed by the later CSC conversion.
fn scatter_add_sparse(coo: &mut CooMatrix<f64>, elem: &DMatrix<f64>, nodes: &[usize]) {
    let n = nodes.len();
    for a in 0..n {
        for i in 0..3 {
            let gi = 3 * nodes[a] + i;
            for b in 0..n {
                for j in 0..3 {
                    let v = elem[(3 * a + i, 3 * b + j)];
                    if v != 0.0 {
                        coo.push(gi, 3 * nodes[b] + j, v);
                    }
                }
            }
        }
    }
}

/// Recover per-node averaged Voigt stress from a solved displacement
/// field on a mixed-element mesh.
///
/// Each element's centroid stress `σ = D·B·uₑ` is computed through the
/// [`SolidElement`] trait and accumulated to its nodes; the per-node
/// value is the average over the elements that touch it.
pub fn recover_nodal_stress_mixed(
    mesh: &Mesh,
    d_matrix: &Matrix6<f64>,
    u: &DVector<f64>,
) -> Vec<[f64; 6]> {
    let n_nodes = mesh.nodes.len();
    let mut acc = vec![[0.0_f64; 6]; n_nodes];
    let mut count = vec![0u32; n_nodes];

    for block in &mesh.element_blocks {
        if !is_supported_solid(block.element_type) {
            continue;
        }
        let npe = block.element_type.nodes_per_element();
        if block.connectivity.len() < npe {
            continue;
        }
        let n_elem = block.connectivity.len() / npe;
        for e in 0..n_elem {
            let Ok((nodes, coords)) = gather_element(mesh, block, e) else {
                continue;
            };
            let Some(elem) = make_element(block.element_type, &coords) else {
                continue;
            };
            // Element displacement vector.
            let mut ue = DVector::zeros(3 * npe);
            for (a, &nd) in nodes.iter().enumerate() {
                for i in 0..3 {
                    ue[3 * a + i] = u[3 * nd + i];
                }
            }
            let Some(sigma) = elem.recover_stress(d_matrix, &ue) else {
                continue;
            };
            for &nd in &nodes {
                for (acc_k, &sig_k) in acc[nd].iter_mut().zip(sigma.iter()) {
                    *acc_k += sig_k;
                }
                count[nd] += 1;
            }
        }
    }
    for n in 0..n_nodes {
        let c = count[n].max(1) as f64;
        for k in 0..6 {
            acc[n][k] /= c;
        }
    }
    acc
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::meshgen::{structured_hex_mesh, structured_tet10_mesh};
    use crate::native_solver::{elasticity_matrix, structured_box_mesh};
    use crate::FemMaterial;

    #[test]
    fn assembly_rejects_empty_mesh() {
        let mesh = Mesh::new("empty");
        let d = elasticity_matrix(&FemMaterial::default()).unwrap();
        assert!(matches!(
            assemble_global_stiffness(&mesh, &d),
            Err(NativeSolverError::EmptyMesh)
        ));
    }

    #[test]
    fn assembly_rejects_mesh_without_solid_elements() {
        let mut mesh = Mesh::new("surface");
        mesh.nodes.push(nalgebra::Vector3::zeros());
        mesh.element_blocks.push(ElementBlock {
            element_type: ElementType::Tri3,
            connectivity: vec![0, 0, 0],
        });
        let d = elasticity_matrix(&FemMaterial::default()).unwrap();
        assert!(matches!(
            assemble_global_stiffness(&mesh, &d),
            Err(NativeSolverError::NoTetBlock)
        ));
    }

    #[test]
    fn global_stiffness_is_symmetric_for_each_element_family() {
        let d = elasticity_matrix(&FemMaterial::default()).unwrap();
        // Round-10: structured_hex_mesh / _tet10_mesh are now fallible.
        // The hard-coded 2×1×1 / 1×1×1 dimensions are well within the
        // allocation cap — `expect` is the honest "this can't fail
        // with these args" assertion.
        for mesh in [
            structured_box_mesh(2.0, 1.0, 1.0, 2, 1, 1).expect("valid box params"),
            structured_hex_mesh(2.0, 1.0, 1.0, 2, 1, 1).expect("valid hex params"),
            structured_tet10_mesh(2.0, 1.0, 1.0, 1, 1, 1).expect("valid tet10 params"),
        ] {
            let k = assemble_global_stiffness(&mesh, &d).unwrap();
            for i in 0..k.nrows() {
                for j in 0..k.ncols() {
                    let asym = (k[(i, j)] - k[(j, i)]).abs();
                    assert!(
                        asym < 1e-6 * k[(i, i)].abs().max(1.0),
                        "global K not symmetric"
                    );
                }
            }
        }
    }

    #[test]
    fn global_stiffness_annihilates_rigid_translation() {
        // A uniform translation of every node yields zero global force
        // — the assembly preserved the element rigid-body property.
        let d = elasticity_matrix(&FemMaterial::default()).unwrap();
        // Round-10: structured_hex_mesh / _tet10_mesh are now fallible.
        for mesh in [
            structured_hex_mesh(2.0, 1.0, 1.0, 2, 1, 1).expect("valid hex params"),
            structured_tet10_mesh(2.0, 1.0, 1.0, 1, 1, 1).expect("valid tet10 params"),
        ] {
            let k = assemble_global_stiffness(&mesh, &d).unwrap();
            let mut t = DVector::zeros(k.nrows());
            for n in 0..mesh.nodes.len() {
                t[3 * n] = 1.0; // uniform +X
            }
            let f = &k * &t;
            assert!(f.norm() < 1e-6 * k.norm(), "rigid translation gave force");
        }
    }

    #[test]
    fn dense_and_sparse_assembly_agree_entry_by_entry() {
        // Round-1 L2: the dense `assemble_global_stiffness` and the
        // sparse `assemble_global_stiffness_sparse` must produce the
        // identical global K. They share the element path but scatter it
        // into different containers (a DMatrix vs a COO→CSC); this test
        // pins them together so a future edit to one cannot silently
        // drift from the other.
        //
        // Use a genuinely MIXED-element mesh: a Hex8 cube on nodes 0–7
        // plus a separate, non-degenerate Tet4 on nodes 8–11, so the
        // assembly walks two element families.
        let mut mesh = Mesh::new("mixed");
        // Hex8 unit cube, standard node order (bottom CCW then top).
        mesh.nodes.extend_from_slice(&[
            nalgebra::Vector3::new(0.0, 0.0, 0.0),
            nalgebra::Vector3::new(1.0, 0.0, 0.0),
            nalgebra::Vector3::new(1.0, 1.0, 0.0),
            nalgebra::Vector3::new(0.0, 1.0, 0.0),
            nalgebra::Vector3::new(0.0, 0.0, 1.0),
            nalgebra::Vector3::new(1.0, 0.0, 1.0),
            nalgebra::Vector3::new(1.0, 1.0, 1.0),
            nalgebra::Vector3::new(0.0, 1.0, 1.0),
        ]);
        // A second, well-separated unit tet on fresh nodes 8–11.
        mesh.nodes.extend_from_slice(&[
            nalgebra::Vector3::new(5.0, 0.0, 0.0),
            nalgebra::Vector3::new(6.0, 0.0, 0.0),
            nalgebra::Vector3::new(5.0, 1.0, 0.0),
            nalgebra::Vector3::new(5.0, 0.0, 1.0),
        ]);
        mesh.element_blocks.push(ElementBlock {
            element_type: ElementType::Hex8,
            connectivity: vec![0, 1, 2, 3, 4, 5, 6, 7],
        });
        mesh.element_blocks.push(ElementBlock {
            element_type: ElementType::Tet4,
            connectivity: vec![8, 9, 10, 11],
        });

        let d = elasticity_matrix(&FemMaterial::default()).unwrap();
        let dense = assemble_global_stiffness(&mesh, &d).unwrap();
        let sparse = assemble_global_stiffness_sparse(&mesh, &d).unwrap();
        let sparse_as_dense = DMatrix::from(&sparse);

        assert_eq!(dense.shape(), sparse_as_dense.shape());
        // Tight tolerance — the two paths do the same f64 arithmetic, so
        // they should agree to round-off.
        let max_diff = (&dense - &sparse_as_dense).amax();
        assert!(
            max_diff < 1e-9,
            "dense and sparse global K differ by {max_diff} (must be < 1e-9)"
        );
    }

    #[test]
    fn consistent_mass_total_equals_three_rho_v() {
        // The assembled mass matrix must carry the whole mesh mass.
        let d = elasticity_matrix(&FemMaterial::default()).unwrap();
        let density = 7850.0;
        let (lx, ly, lz) = (2.0, 1.0, 3.0);
        // Round-10: structured_hex_mesh is now fallible.
        let mesh = structured_hex_mesh(lx, ly, lz, 2, 2, 2).expect("valid hex params");
        let (_k, m) = assemble_global_stiffness_mass_mixed(&mesh, &d, density).unwrap();
        let total: f64 = m.iter().sum();
        let want = 3.0 * density * lx * ly * lz;
        assert!(
            (total - want).abs() < 1e-6 * want,
            "mass sum {total} ≠ 3ρV {want}"
        );
    }
}
