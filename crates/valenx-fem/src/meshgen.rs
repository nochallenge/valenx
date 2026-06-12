//! Structured box-mesh generators for the new element types
//! (Phase 24.8).
//!
//! [`crate::native_solver::structured_box_mesh`] builds a Tet4 box so
//! the original solvers are testable without an external mesher. This
//! module is the same idea for the two new continuum elements:
//!
//! - [`structured_hex_mesh`] — an axis-aligned box meshed with one
//!   [`ElementType::Hex8`] per cell (the natural, distortion-free hex
//!   mesh of a box).
//! - [`structured_tet10_mesh`] — an axis-aligned box meshed with
//!   10-node [`ElementType::Tet10`] tetrahedra (six per cell, the same
//!   Kuhn split as the Tet4 mesh, with the mid-edge nodes added).
//!
//! Both lay their **corner** nodes out `x`-fastest then `y` then `z`,
//! so corner `(i,j,k)` is corner-node index
//! `i + (nx+1)·j + (nx+1)(ny+1)·k` — identical to the Tet4 generator,
//! so the same face-selection arithmetic picks constraint / load nodes.
//! The Tet10 generator appends its mid-edge nodes *after* all corner
//! nodes, so corner indexing is unchanged.
//!
//! Round-10 cleanup: both generators used to `assert!` on `nx > 0`
//! and `Vec::with_capacity(nx*ny*nz*8)` unchecked. They now return
//! [`Result<Mesh, NativeSolverError>`] so a hostile / mistyped
//! parameter set can't panic the host process; the cap matches the
//! round-5 fix on
//! [`crate::native_solver::structured_box_mesh`].

use std::collections::HashMap;

use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::native_solver::NativeSolverError;

/// Maximum node count either generator will allocate. Matches the
/// cap in [`crate::native_solver::structured_box_mesh`] —
/// `MAX_NODES = 200_000_000` corresponds to a multi-GiB mesh, well
/// beyond any realistic test or validation grid. Past the cap, the
/// implied allocation is in the same class as a hostile recipe trying
/// to OOM the process.
const MAX_NODES: usize = 200_000_000;

/// Bound-check the `(nx+1)·(ny+1)·(nz+1)` corner-node count *and*
/// the `nx·ny·nz·entries_per_cell` connectivity capacity before
/// either `Vec::with_capacity` runs. Returns the validated
/// `(corner_node_count, cell_count)` pair.
///
/// `entries_per_cell` is 8 for the Hex8 generator (one hex per cell)
/// and `60 = 6·10` for the Tet10 generator (six 10-node tets per
/// cell). Either way the connectivity vector is bounded by
/// `cell_count * entries_per_cell` and must not silently wrap.
fn validate_grid(
    nx: usize,
    ny: usize,
    nz: usize,
    lx: f64,
    ly: f64,
    lz: f64,
    entries_per_cell: usize,
) -> Result<(usize, usize), NativeSolverError> {
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
    let cell_count = nx
        .checked_mul(ny)
        .and_then(|v| v.checked_mul(nz))
        .ok_or_else(|| NativeSolverError::InvalidParams {
            reason: "cell count nx*ny*nz overflowed usize".to_string(),
        })?;
    // Check the connectivity capacity (cell_count * entries_per_cell)
    // for overflow too. Pre-fix `Vec::with_capacity(nx*ny*nz*8)` would
    // wrap on 32-bit targets / huge dims and back-fill with bogus data.
    let _conn_capacity = cell_count.checked_mul(entries_per_cell).ok_or_else(|| {
        NativeSolverError::InvalidParams {
            reason: format!("connectivity capacity nx*ny*nz*{entries_per_cell} overflowed usize"),
        }
    })?;
    if node_count > MAX_NODES {
        return Err(NativeSolverError::InvalidParams {
            reason: format!(
                "structured box mesh would have {node_count} nodes, \
                 cap is {MAX_NODES}"
            ),
        });
    }
    Ok((node_count, cell_count))
}

/// Build a structured **Hex8** mesh of the axis-aligned box
/// `[0,lx]×[0,ly]×[0,lz]`, `nx×ny×nz` cells, one 8-node hexahedron per
/// cell.
///
/// Corner-node layout is `x`-fastest then `y` then `z`: node
/// `(i,j,k) = i + (nx+1)·j + (nx+1)(ny+1)·k`. Each cell's eight nodes
/// are emitted in the standard hex order
/// [`crate::elements::Hex8`] expects — the `ζ=-1` face CCW
/// (`n0..n4`) then the `ζ=+1` face (`n4..n8`).
///
/// # Errors
///
/// Returns [`NativeSolverError::InvalidParams`] when any of `nx`, `ny`,
/// `nz` is zero, when any length is non-positive / non-finite, or when
/// the implied node-count `(nx+1)*(ny+1)*(nz+1)` overflows `usize` or
/// exceeds the `MAX_NODES` cap.
///
/// Round-10 fix: this used to `assert!` on the inputs and panic.
/// Migrated to `Result` so a hostile / mistyped configuration can't
/// crash the host process.
pub fn structured_hex_mesh(
    lx: f64,
    ly: f64,
    lz: f64,
    nx: usize,
    ny: usize,
    nz: usize,
) -> Result<Mesh, NativeSolverError> {
    // One Hex8 per cell → 8 connectivity entries per cell.
    let (_node_count, cell_count) = validate_grid(nx, ny, nz, lx, ly, lz, 8)?;
    let mut mesh = Mesh::new("native-fem-hex-box");
    let stride_y = nx + 1;
    let stride_z = (nx + 1) * (ny + 1);
    let node_id =
        |i: usize, j: usize, k: usize| -> u32 { (i + stride_y * j + stride_z * k) as u32 };
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
    let mut conn: Vec<u32> = Vec::with_capacity(cell_count * 8);
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                // Standard hex node order: bottom face CCW, then top.
                conn.extend_from_slice(&[
                    node_id(i, j, k),
                    node_id(i + 1, j, k),
                    node_id(i + 1, j + 1, k),
                    node_id(i, j + 1, k),
                    node_id(i, j, k + 1),
                    node_id(i + 1, j, k + 1),
                    node_id(i + 1, j + 1, k + 1),
                    node_id(i, j + 1, k + 1),
                ]);
            }
        }
    }
    mesh.element_blocks.push(ElementBlock {
        element_type: ElementType::Hex8,
        connectivity: conn,
    });
    mesh.recompute_stats();
    Ok(mesh)
}

/// Build a structured **Tet10** mesh of the axis-aligned box
/// `[0,lx]×[0,ly]×[0,lz]`, `nx×ny×nz` cells, six 10-node quadratic
/// tetrahedra per cell.
///
/// The corner nodes are exactly the corner grid of
/// [`crate::native_solver::structured_box_mesh`] — node `(i,j,k)` keeps
/// the index `i + (nx+1)·j + (nx+1)(ny+1)·k`. The six tets per cell use
/// the same Kuhn split; each tet then gains its six mid-edge nodes.
/// Mid-edge nodes are **deduplicated** across elements (an edge shared
/// by two tets gets one node) and appended after the corner block, so
/// the resulting mesh is conforming.
///
/// # Errors
///
/// Returns [`NativeSolverError::InvalidParams`] when any of `nx`, `ny`,
/// `nz` is zero, when any length is non-positive / non-finite, or when
/// the implied node-count overflows `usize` or exceeds the
/// `MAX_NODES` cap.
///
/// Round-10 fix: this used to `assert!` on the inputs and panic.
/// Migrated to `Result` so a hostile / mistyped configuration can't
/// crash the host process.
pub fn structured_tet10_mesh(
    lx: f64,
    ly: f64,
    lz: f64,
    nx: usize,
    ny: usize,
    nz: usize,
) -> Result<Mesh, NativeSolverError> {
    // Six tets per cell × 10 nodes per tet = 60 connectivity entries
    // per cell. Check the multiplication so the
    // `Vec::with_capacity(corner_conn.len() * 10)` below cannot wrap.
    let (_node_count, cell_count) = validate_grid(nx, ny, nz, lx, ly, lz, 60)?;
    let mut mesh = Mesh::new("native-fem-tet10-box");
    let stride_y = nx + 1;
    let stride_z = (nx + 1) * (ny + 1);
    let node_id =
        |i: usize, j: usize, k: usize| -> u32 { (i + stride_y * j + stride_z * k) as u32 };
    // Corner nodes — identical layout to the Tet4 generator.
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
    let n_corner = mesh.nodes.len();

    // First collect the corner connectivity of the six tets per cell.
    let mut corner_conn: Vec<[u32; 4]> = Vec::with_capacity(cell_count * 6);
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
                for tet in [
                    [n0, n1, n2, n6],
                    [n0, n2, n3, n6],
                    [n0, n3, n7, n6],
                    [n0, n7, n4, n6],
                    [n0, n4, n5, n6],
                    [n0, n5, n1, n6],
                ] {
                    corner_conn.push(tet);
                }
            }
        }
    }

    // Edge midpoint table: a sorted corner-index pair → its mid-edge
    // node index. An edge shared by several tets resolves to one node.
    let mut edge_node: HashMap<(u32, u32), u32> = HashMap::new();
    let mut connectivity: Vec<u32> = Vec::with_capacity(corner_conn.len() * 10);
    // Tet10 edge order — must match `crate::elements::TET10_EDGES`.
    const EDGES: [(usize, usize); 6] = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)];
    for tet in &corner_conn {
        // Corner nodes first.
        connectivity.extend_from_slice(tet);
        // Then the six mid-edge nodes.
        for &(p, q) in &EDGES {
            let (a, b) = (tet[p], tet[q]);
            let key = if a < b { (a, b) } else { (b, a) };
            let mid = *edge_node.entry(key).or_insert_with(|| {
                let coord = (mesh.nodes[a as usize] + mesh.nodes[b as usize]) * 0.5;
                mesh.nodes.push(coord);
                (mesh.nodes.len() - 1) as u32
            });
            connectivity.push(mid);
        }
    }
    debug_assert!(mesh.nodes.len() >= n_corner);

    mesh.element_blocks.push(ElementBlock {
        element_type: ElementType::Tet10,
        connectivity,
    });
    mesh.recompute_stats();
    Ok(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::elements::{Hex8, SolidElement, Tet10};

    #[test]
    fn hex_mesh_single_cell_has_one_brick_eight_nodes() {
        let mesh = structured_hex_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid params");
        assert_eq!(mesh.nodes.len(), 8);
        assert_eq!(mesh.element_blocks.len(), 1);
        assert_eq!(mesh.element_blocks[0].element_type, ElementType::Hex8);
        assert_eq!(mesh.element_blocks[0].connectivity.len(), 8);
    }

    #[test]
    fn hex_mesh_cells_have_positive_volume_summing_to_box() {
        let (lx, ly, lz) = (2.0, 3.0, 1.5);
        let (nx, ny, nz) = (2, 2, 2);
        let mesh = structured_hex_mesh(lx, ly, lz, nx, ny, nz).expect("valid params");
        let block = &mesh.element_blocks[0];
        let mut total = 0.0;
        for e in 0..block.connectivity.len() / 8 {
            let mut c = [Vector3::zeros(); 8];
            for (k, ci) in c.iter_mut().enumerate() {
                *ci = mesh.nodes[block.connectivity[8 * e + k] as usize];
            }
            let v = Hex8::new(c).volume().unwrap();
            assert!(v > 0.0, "hex {e} has non-positive volume {v}");
            total += v;
        }
        assert!(
            (total - lx * ly * lz).abs() < 1e-9,
            "hex volume sum {total} ≠ box {}",
            lx * ly * lz
        );
    }

    #[test]
    fn tet10_mesh_single_cell_has_six_tets() {
        let mesh = structured_tet10_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid params");
        assert_eq!(mesh.element_blocks.len(), 1);
        assert_eq!(mesh.element_blocks[0].element_type, ElementType::Tet10);
        // 6 tets × 10 nodes.
        assert_eq!(mesh.element_blocks[0].connectivity.len(), 60);
        // 8 corner + the unique edge nodes of a single cell.
        assert!(mesh.nodes.len() > 8, "mid-edge nodes were not added");
    }

    #[test]
    fn tet10_mesh_cells_have_positive_volume_summing_to_box() {
        let (lx, ly, lz) = (2.0, 1.0, 3.0);
        let (nx, ny, nz) = (2, 1, 2);
        let mesh = structured_tet10_mesh(lx, ly, lz, nx, ny, nz).expect("valid params");
        let block = &mesh.element_blocks[0];
        let mut total = 0.0;
        for e in 0..block.connectivity.len() / 10 {
            let mut c = [Vector3::zeros(); 10];
            for (k, ci) in c.iter_mut().enumerate() {
                *ci = mesh.nodes[block.connectivity[10 * e + k] as usize];
            }
            let v = Tet10::new(c).volume().unwrap();
            assert!(v > 0.0, "tet10 {e} has non-positive volume {v}");
            total += v;
        }
        assert!(
            (total - lx * ly * lz).abs() < 1e-9,
            "tet10 volume sum {total} ≠ box {}",
            lx * ly * lz
        );
    }

    #[test]
    fn tet10_mesh_mid_edge_nodes_are_deduplicated() {
        // A 2×1×1 box: the internal face between the two cells is
        // shared, so its edge nodes must not be duplicated. Count the
        // unique node coordinates and compare to the node array length.
        let mesh = structured_tet10_mesh(2.0, 1.0, 1.0, 2, 1, 1).expect("valid params");
        let mut seen = std::collections::HashSet::new();
        for p in &mesh.nodes {
            let key = (
                (p.x * 1e9).round() as i64,
                (p.y * 1e9).round() as i64,
                (p.z * 1e9).round() as i64,
            );
            assert!(seen.insert(key), "duplicate node coordinate {p}");
        }
        assert_eq!(seen.len(), mesh.nodes.len());
    }

    #[test]
    fn tet10_mid_edge_nodes_lie_on_the_edges() {
        // Every mid-edge node of every element must sit exactly at the
        // midpoint of its two corner nodes.
        use crate::elements::TET10_EDGES;
        let mesh = structured_tet10_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid params");
        let block = &mesh.element_blocks[0];
        for e in 0..block.connectivity.len() / 10 {
            let base = 10 * e;
            for (m, &(p, q)) in TET10_EDGES.iter().enumerate() {
                let cp = mesh.nodes[block.connectivity[base + p] as usize];
                let cq = mesh.nodes[block.connectivity[base + q] as usize];
                let cm = mesh.nodes[block.connectivity[base + 4 + m] as usize];
                let mid = (cp + cq) * 0.5;
                assert!((cm - mid).norm() < 1e-12, "mid-edge node off the edge");
            }
        }
    }

    // --- Round-10 M4 RED→GREEN tests: meshgen Result migration ---

    /// Round-10 M4: pre-fix `structured_hex_mesh(0, 1, 1, ...)`
    /// panicked via `assert!(nx > 0 ...)`. Migrated to Result so a
    /// hostile / mistyped configuration can't crash the host process.
    #[test]
    fn structured_hex_mesh_rejects_zero_subdivisions() {
        let err = structured_hex_mesh(1.0, 1.0, 1.0, 0, 1, 1).expect_err("nx=0 must be rejected");
        let msg = format!("{err}");
        assert!(msg.contains("positive"), "msg: {msg}");
    }

    #[test]
    fn structured_hex_mesh_rejects_non_finite_dimension() {
        let err = structured_hex_mesh(f64::NAN, 1.0, 1.0, 1, 1, 1)
            .expect_err("NaN length must be rejected");
        let msg = format!("{err}");
        assert!(msg.contains("finite"), "msg: {msg}");
    }

    #[test]
    fn structured_hex_mesh_rejects_oversized_grid() {
        // Past the cap but well inside usize.
        let err = structured_hex_mesh(1.0, 1.0, 1.0, 1000, 1000, 1000)
            .expect_err("oversized grid must be rejected");
        let msg = format!("{err}");
        assert!(msg.contains("cap"), "msg: {msg}");
    }

    #[test]
    fn structured_tet10_mesh_rejects_zero_subdivisions() {
        let err = structured_tet10_mesh(1.0, 1.0, 1.0, 1, 0, 1).expect_err("ny=0 must be rejected");
        let msg = format!("{err}");
        assert!(msg.contains("positive"), "msg: {msg}");
    }
}
